mod_pca_ui <- function(id) {
  ns <- NS(id)
  
  tabPanel("PCA",
           fluidRow(
             # Added ns() to all IDs below
             column(5, radioButtons(ns("pca.norm"), "Normalization and Scaling", 
                                    c("Log Normalization" = "LogNorm", "SCTransform"="SCTransform"), inline = TRUE)),
             column(2, actionButton(ns("pca.run"), "Run PCA", class = "btn-success"))
           ), 
           fluidRow(
             column(7, radioButtons(ns("pca.plot.type"), "Plot Type", 
                                    c("PC Variance", "Feature Loading", "Heatmap", "2-D PCA"), inline = TRUE)),
             column(2, actionButton(ns("pca.plot.run"), "Plot", class = "btn-success"))
           ),
           fluidRow(
             column(7, plotOutput(ns("plot.pca"), height = "500px", width = "600px")),
             column(3, 
                    numericInput(ns("pca.number"), "PC Variance: Number of PCs", min = 1, max = 50, value = 20),
                    numericInput(ns("pca.load"), "Feature Loading: PC #", min = 1, max = 50, value = 1),
                    numericInput(ns("pca.heatmap"), "Heatmap: PC#", min = 1, max = 50, value = 1),
                    numericInput(ns("pca.2d.1"), "2-D PCA: PC # (x-axis)", min = 1, max = 50, value = 1),
                    numericInput(ns("pca.2d.2"), "2-D PCA: PC # (y-axis)", min = 1, max = 50, value = 2),
                    actionButton(ns("pca.filter.run"), "Confirm PCs for Next Step", class = "btn-success"))
           ),
           fluidRow(
             column(2, downloadButton(ns("download.pca"), "Download Seurat Object"))
           )
  )
}

mod_pca_server <- function(id, seurat_obj_qc) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    completed <- reactiveVal(FALSE)

    # Logic to handle the PCA computation
    data.pca <- eventReactive(input$pca.run, {
      req(seurat_obj_qc()) # Ensure data exists from previous module
      
      # Use the object passed from the previous module
      srt <- seurat_obj_qc()
      
      withProgress(message = 'Running PCA...', value = 0, {
        if (input$pca.norm == "SCTransform") {
          srt <- SCTransform(object = srt, method = "glmGamPoi", verbose = FALSE)
          srt <- RunPCA(srt, verbose = FALSE)
        } else {
          srt <- NormalizeData(object = srt, verbose = FALSE)
          srt <- FindVariableFeatures(object = srt, verbose = FALSE)
          srt <- ScaleData(object = srt, verbose = FALSE)
          srt <- RunPCA(srt, verbose = FALSE)
        }
      })
      return(srt)
    })
    
    # Plotting Logic
    plotInput.pca <- eventReactive(input$pca.plot.run, {
      req(data.pca())
      
      if (input$pca.plot.type == "PC Variance") {
        ElbowPlot(data.pca(), ndims = input$pca.number)
      } else if (input$pca.plot.type == "Feature Loading") {
        VizDimLoadings(data.pca(), dims = input$pca.load, reduction = "pca")
      } else if (input$pca.plot.type == "Heatmap") {
        DimHeatmap(data.pca(), dims = input$pca.heatmap, cells = 500, balanced = TRUE)
      } else {
        DimPlot(data.pca(), reduction = "pca", dims = c(input$pca.2d.1, input$pca.2d.2)) + NoLegend()
      }
    })
    
    output$plot.pca <- renderPlot({
      plotInput.pca()
    }, res = 96)
    
    observeEvent(input$pca.filter.run, {
      completed(TRUE)
    })

    return(list(seurat_obj = data.pca, completed = completed))
    
    # Download handler
    output$download.pca <- downloadHandler(
      filename = function() { paste0("seurat_pca_", Sys.Date(), ".rds") },
      content = function(file) { saveRDS(data.pca(), file) }
    )
  })
}