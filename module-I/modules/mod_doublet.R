mod_doublet_ui <- function(id) {
  ns <- NS(id)
  tabPanel("Doublet Removal",
           fluidRow(
             column(4, radioButtons(ns("dbl.assay"), "10X Genomics Assay", 
                                    choices = c("High Throughput v3.1", "Standard v3.1"), inline = TRUE)),
             column(4, actionButton(ns("dbl.assay.run"), "Calculate Default Parameters", class = "btn-success"))
           ),
           fluidRow(
             column(2, numericInput(ns("dbl.percent"), "Doublet Rate (%)", min = 0, max = 20, value = 0, step = 0.1)),
             column(2, numericInput(ns("dbl.pK"), "pK (Optimal)", value = 0.09, step = 0.01)),
             column(2, numericInput(ns("dbl.pN"), "pN (Default)", value = 0.25, min = 0, max = 1, step = 0.01)),
             column(2, actionButton(ns("dbl.run"), "Run Doublet Finder", class = "btn-success"))
           ),
           fluidRow(
             column(5, plotOutput(ns("plot.dbl.umap"), height = "500px")),
             column(4, plotOutput(ns("plot.dbl.scatter"), height = "500px")),
             column(3, 
                    h4("Classification Summary"),
                    tableOutput(ns("dbl.table")), # Summary table
                    hr(),
                    actionButton(ns("dbl.remove.run"), "Remove Doublets", class = "btn-danger", style="width: 100%"),
                    br(), br(),
                    wellPanel(textOutput(ns("dbl.cell.count"))),
                    downloadButton(ns("download.dbl"), "Download Seurat Object")
             )
           )
  )
}


mod_doublet_server <- function(id, seurat_obj_pca, pca_dims) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # --- 1. Instant Update for Doublet Rate ---
    # We use observeEvent on BOTH the assay choice and the data object
    observeEvent({
      input$dbl.assay
      seurat_obj_pca()
    }, {
      req(seurat_obj_pca())
      srt <- seurat_obj_pca()
      n_cells <- ncol(srt)
      
      rate_mult <- if(input$dbl.assay == "Standard v3.1") 0.8 else 0.4
      calc_percent <- (n_cells / 1000) * rate_mult
      final_percent <- min(max(calc_percent, 0.5), 15)
      
      updateNumericInput(session, "dbl.percent", value = round(final_percent, 2))
    })
    
    # --- 2. Prep Data (UMAP) ---
    # Triggered by EITHER the "Calculate" button OR the "Run" button if prep is missing
    data.dbl.prep <- eventReactive({
      input$dbl.assay.run
      input$dbl.run
    }, {
      req(seurat_obj_pca())
      dims_to_use <- if(is.reactive(pca_dims)) pca_dims() else 20
      
      withProgress(message = 'Preparing UMAP...', value = 0.5, {
        srt <- seurat_obj_pca()
        # Only run UMAP if it hasn't been run or if we explicitly clicked Calculate
        srt <- RunUMAP(srt, dims = 1:dims_to_use, verbose = FALSE)
      })
      return(srt)
    })
    
    # --- 3. pK Calculation ---
    observeEvent(input$dbl.assay.run, {
      req(data.dbl.prep())
      srt <- data.dbl.prep()
      dims_to_use <- if(is.reactive(pca_dims)) pca_dims() else 20
      
      withProgress(message = 'Calculating optimal pK...', {
        sweep.res.list <- paramSweep(srt, PCs = 1:dims_to_use, sct = FALSE)
        sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
        bcmvn <- find.pK(sweep.stats)
        optimal_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
        updateNumericInput(session, "dbl.pK", value = optimal_pK)
      })
    })
    
    # --- 4. Main Doublet Finder Logic ---
    data.dbl <- eventReactive(input$dbl.run, {
      req(data.dbl.prep())
      srt <- data.dbl.prep()
      dims_to_use <- if(is.reactive(pca_dims)) pca_dims() else 20
      
      withProgress(message = 'Identifying Doublets...', {
        # Use the latest input values
        pct <- as.numeric(input$dbl.percent)
        nExp_poi <- round((pct/100) * ncol(srt))
        
        srt <- doubletFinder(srt, PCs = 1:dims_to_use, 
                             pN = input$dbl.pN, 
                             pK = input$dbl.pK, 
                             nExp = nExp_poi, sct = FALSE)
        
        df_col <- tail(grep("DF.classifications", colnames(srt@meta.data), value = TRUE), 1)
        pANN_col <- tail(grep("pANN", colnames(srt@meta.data), value = TRUE), 1)
        
        srt$doublet.class <- srt@meta.data[[df_col]]
        srt$doublet.score <- srt@meta.data[[pANN_col]]
      })
      return(srt)
    })
    
    # --- 5. Outputs ---
    output$dbl.table <- renderTable({
      req(data.dbl())
      df <- as.data.frame(table(data.dbl()$doublet.class))
      colnames(df) <- c("Classification", "Count")
      df
    })
    
    output$plot.dbl.umap <- renderPlot({
      req(data.dbl())
      DimPlot(data.dbl(), reduction = "umap", group.by = "doublet.class", 
              cols = c("Doublet" = "#F8766D", "Singlet" = "darkgrey")) + 
        theme(legend.position = "bottom") + ggtitle("UMAP: Doublet Identification")
    })
    
    output$plot.dbl.scatter <- renderPlot({
      req(data.dbl())
      FeatureScatter(data.dbl(), feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
                     group.by = "doublet.class", cols = c("Doublet" = "#F8766D", "Singlet" = "darkgrey")) +
        theme(legend.position = "none")
    })
    
    # --- 6. Filtering & Return ---
    data.dbl.final <- eventReactive(input$dbl.remove.run, {
      req(data.dbl())
      subset(data.dbl(), doublet.class == "Singlet")
    })
    
    output$dbl.cell.count <- renderText({
      if(is.null(data.dbl.final())) {
        n <- if(!is.null(data.dbl())) ncol(data.dbl()) else 0
        return(paste0("Cells before filtering: ", n))
      }
      paste0("Singlets remaining: ", ncol(data.dbl.final()))
    })
    
    return(data.dbl.final)
  })
}