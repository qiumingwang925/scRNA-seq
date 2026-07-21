## ABOUTME: Shiny module for doublet detection and removal using DoubletFinder.
## ABOUTME: Estimates doublet rate by assay type, runs pK optimization, and visualizes results on UMAP.

mod.doublet.ui <- function(id) {
  ns <- NS(id)
  tabPanel("Doublet Removal", value = "tab.doublet",
           wellPanel(
             strong("Assay Settings"),
             fluidRow(
               column(4, radioButtons(ns("dbl.assay"), "10X Genomics Assay",
                                      choices = c("High Throughput v3.1", "Standard v3.1"), inline = TRUE)),
               column(4, actionButton(ns("dbl.assay.run"), "Calculate Default Parameters", class = "btn-success"))
             ),
             fluidRow(
               column(2, numericInput(ns("dbl.percent"), "Doublet Rate (%)", min = 0, max = 20, value = 0, step = 0.1)),
               column(2, numericInput(ns("dbl.pK"), "pK (Optimal)", value = 0.09, step = 0.01)),
               column(2, numericInput(ns("dbl.pN"), "pN (Default)", value = 0.25, min = 0, max = 1, step = 0.01)),
               column(2, actionButton(ns("dbl.run"), "Run Doublet Finder", class = "btn-success")),
               column(4, uiOutput(ns("dbl.run.hint")))
             )
           ),
           wellPanel(
             strong("Doublet Results"),
             fluidRow(
               column(5, plotOutput(ns("plot.dbl.umap"), height = "500px")),
               column(4, plotOutput(ns("plot.dbl.scatter"), height = "500px")),
               column(3,
                      h4("Classification Summary"),
                      tableOutput(ns("dbl.table")),
                      hr(),
                      actionButton(ns("dbl.remove.run"), "Remove Doublets", class = "btn-danger", style="width: 100%"),
                      br(), br(),
                      wellPanel(textOutput(ns("dbl.cell.count")))
               )
             )
           )
  )
}


mod.doublet.server <- function(id, seurat.obj.pca, pca.dims, upstream.completed = reactive(TRUE)) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    completed <- reactiveVal(FALSE)

    # Soft guard: warn if the PCA step isn't done, but let the user proceed
    dbl.run.trigger <- reactiveVal(0)
    observeEvent(input$dbl.run, {
      if (isTRUE(upstream.completed())) {
        dbl.run.trigger(dbl.run.trigger() + 1)
      } else {
        showModal(modalDialog(
          title = "Upstream step not completed",
          "The PCA step hasn't been completed yet. Proceed anyway?",
          footer = tagList(
            modalButton("Cancel"),
            actionButton(ns("dbl.run.proceed"), "Proceed anyway", class = "btn-warning")
          )
        ))
      }
    })
    observeEvent(input$dbl.run.proceed, {
      removeModal()
      dbl.run.trigger(dbl.run.trigger() + 1)
    })

    output$dbl.run.hint <- renderUI({
      if (!isTRUE(upstream.completed())) {
        tags$small(style = "color:#c0392b;", "PCA step not completed yet.")
      }
    })

    # Number of PCA dimensions to use for UMAP and doublet detection
    dims.to.use <- reactive({ pca.dims() })

    # --- 1. Instant Update for Doublet Rate ---
    # We use observeEvent on BOTH the assay choice and the data object
    observeEvent({
      input$dbl.assay
      seurat.obj.pca()
    }, {
      req(seurat.obj.pca())

      doublet.rate.percent <- if (input$dbl.assay == "Standard v3.1") 8 else 4

      updateNumericInput(session, "dbl.percent", value = doublet.rate.percent)
    })
    
    # --- 2. Prep Data (UMAP) ---
    # Triggered by EITHER the "Calculate" button OR the "Run" button if prep is missing
    data.dbl.prep <- eventReactive({
      input$dbl.assay.run
      dbl.run.trigger()
    }, {
      req(seurat.obj.pca())

      withProgress(message = 'Preparing UMAP...', value = 0.5, {
        srt <- seurat.obj.pca()
        srt <- RunUMAP(srt, dims = 1:dims.to.use(), verbose = FALSE)
      })
      return(srt)
    })
    
    # --- 3. pK Calculation ---
    observeEvent(input$dbl.assay.run, {
      req(data.dbl.prep())
      srt <- data.dbl.prep()

      withProgress(message = 'Calculating optimal pK...', {
        sweep.res.list <- paramSweep(srt, PCs = 1:dims.to.use(), sct = FALSE)
        sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
        bcmvn <- find.pK(sweep.stats)
        optimal.pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
        updateNumericInput(session, "dbl.pK", value = optimal.pK)
      })
    })
    
    # --- 4. Main Doublet Finder Logic ---
    data.dbl <- eventReactive(dbl.run.trigger(), {
      req(data.dbl.prep())
      srt <- data.dbl.prep()

      withProgress(message = 'Identifying Doublets...', {
        # Use the latest input values
        pct <- as.numeric(input$dbl.percent)
        nExp.poi <- round((pct/100) * ncol(srt))
        
        srt <- doubletFinder(srt, PCs = 1:dims.to.use(),
                             pN = input$dbl.pN, 
                             pK = input$dbl.pK, 
                             nExp = nExp.poi, sct = FALSE)
        
        df.col <- tail(grep("DF.classifications", colnames(srt@meta.data), value = TRUE), 1)
        pANN.col <- tail(grep("pANN", colnames(srt@meta.data), value = TRUE), 1)
        
        srt$doublet.class <- srt@meta.data[[df.col]]
        srt$doublet.score <- srt@meta.data[[pANN.col]]
      })
      return(srt)
    }, ignoreInit = TRUE)

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
      completed(TRUE)
      srt <- subset(data.dbl(), doublet.class == "Singlet")

      # Every surviving cell is a singlet, so the classification is constant
      # from here on; DoubletFinder's own run-parameterised columns are
      # superseded by doublet.score, which is kept because it still varies.
      # Dropped here rather than in data.dbl() so this tab's own plots keep
      # working when the user navigates back.
      run.cols <- grep("^(pANN|DF\\.classifications)", colnames(srt@meta.data), value = TRUE)
      for (col in c("doublet.class", run.cols)) srt[[col]] <- NULL
      srt
    })
    
    output$dbl.cell.count <- renderText({
      if(is.null(data.dbl.final())) {
        n <- if(!is.null(data.dbl())) ncol(data.dbl()) else 0
        return(paste0("Cells before filtering: ", n))
      }
      paste0("Singlets remaining: ", ncol(data.dbl.final()))
    })
    
    return(list(seurat.obj = data.dbl.final, completed = completed))
  })
}