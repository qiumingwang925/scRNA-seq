## ABOUTME: Shiny module for normalization and principal component analysis.
## ABOUTME: Supports LogNormalize and SCTransform, with elbow/loading/heatmap/2D PCA plots.

mod.pca.ui <- function(id) {
  ns <- NS(id)
  
  tabPanel("PCA",
           wellPanel(
             strong("Normalization"),
             fluidRow(
               column(5, radioButtons(ns("pca.norm"), "Normalization and Scaling",
                                      c("Log Normalization" = "LogNorm", "SCTransform"="SCTransform"), inline = TRUE)),
               column(2, actionButton(ns("pca.run"), "Run PCA", class = "btn-success"))
             )
           ),
           wellPanel(
             fluidRow(
               column(3, strong("PCA Plot")),
               column(5, radioButtons(ns("pca.plot.type"), "Plot Type",
                                      c("PC Variance", "Feature Loading", "Heatmap", "2-D PCA"), inline = TRUE)),
               column(2, actionButton(ns("pca.plot.run"), "Plot", class = "btn-success"))
             ),
             plotOutput(ns("plot.pca"), height = "500px", width = "700px")
           ),
           wellPanel(
             strong("PC Selection"),
             fluidRow(
               column(2, numericInput(ns("pca.number"), "PC Variance: Number of PCs", min = 1, max = 50, value = 20)),
               column(2, numericInput(ns("pca.load"), "Feature Loading: PC #", min = 1, max = 50, value = 1)),
               column(2, numericInput(ns("pca.heatmap"), "Heatmap: PC#", min = 1, max = 50, value = 1)),
               column(2, numericInput(ns("pca.2d.1"), "2-D PCA: PC # (x-axis)", min = 1, max = 50, value = 1)),
               column(2, numericInput(ns("pca.2d.2"), "2-D PCA: PC # (y-axis)", min = 1, max = 50, value = 2))
             ),
             actionButton(ns("pca.filter.run"), "Confirm PCs for Next Step", class = "btn-success")
           )
  )
}

mod.pca.server <- function(id, seurat.obj.qc) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    completed <- reactiveVal(FALSE)
    pca.plotted <- reactiveVal(FALSE)

    shinyjs::disable("pca.filter.run")

    # Logic to handle the PCA computation
    data.pca <- eventReactive(input$pca.run, {
      req(seurat.obj.qc()) # Ensure data exists from previous module
      
      # Use the object passed from the previous module
      srt <- seurat.obj.qc()
      
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
    plot.input.pca <- eventReactive(input$pca.plot.run, {
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
      plot.input.pca()
    }, res = 96)

    # Invalidate plotted state when PCA is re-run
    observeEvent(input$pca.run, {
      pca.plotted(FALSE)
    })

    # Invalidate plotted state when any plot parameter changes
    observe({
      input$pca.number; input$pca.load; input$pca.heatmap
      input$pca.2d.1; input$pca.2d.2; input$pca.plot.type
      pca.plotted(FALSE)
    })

    # Mark as plotted after plot button click
    observeEvent(input$pca.plot.run, {
      pca.plotted(TRUE)
    })

    # Enable/disable confirm button based on plotted state
    observe({
      if (pca.plotted()) {
        shinyjs::enable("pca.filter.run")
      } else {
        shinyjs::disable("pca.filter.run")
      }
    })

    # Rename Plot button after first use
    observeEvent(input$pca.plot.run, {
      updateActionButton(session, "pca.plot.run", label = "Update Plot")
    }, once = TRUE)

    observeEvent(input$pca.filter.run, {
      completed(TRUE)
    })

    return(list(seurat.obj = data.pca, completed = completed))
  })
}