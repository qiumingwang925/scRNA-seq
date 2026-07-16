## ABOUTME: Shiny module for normalization and principal component analysis.
## ABOUTME: Supports LogNormalize and SCTransform, with elbow/loading/heatmap/2D PCA plots.

mod.pca.ui <- function(id) {
  ns <- NS(id)
  
  tabPanel("PCA", value = "tab.pca",
           wellPanel(
             strong("Normalization"),
             fluidRow(
               column(5, radioButtons(ns("pca.norm"), "Normalization and Scaling",
                                      c("Log Normalization" = "LogNorm", "SCTransform"="SCTransform"), inline = TRUE)),
               column(2, actionButton(ns("pca.run"), "Run PCA", class = "btn-success")),
               column(5, uiOutput(ns("pca.run.hint")))
             )
           ),
           fluidRow(class = "plot.params.row",
             column(8,
               wellPanel(
                 fluidRow(
                   column(6, strong("PCA Plot"),
                     radioButtons(ns("pca.plot.type"), NULL,
                                  c("PC Variance", "Feature Loading", "Heatmap", "2-D PCA"), inline = TRUE)),
                   column(3, offset = 3, actionButton(ns("pca.plot.run"), "Plot", class = "btn-success", style = "width: 100%"))
                 ),
                 tags$div(class = "square.plot",
                   plotOutput(ns("plot.pca"), width = "100%", height = "100%")
                 )
               )
             ),
             column(4,
               wellPanel(
                 strong("PC Selection"),
                 numericInput(ns("pca.number"), "PC Variance: Number of PCs", min = 1, max = 50, value = 20),
                 numericInput(ns("pca.load"), "Feature Loading: PC #", min = 1, max = 50, value = 1),
                 numericInput(ns("pca.heatmap"), "Heatmap: PC#", min = 1, max = 50, value = 1),
                 numericInput(ns("pca.2d.1"), "2-D PCA: PC # (x-axis)", min = 1, max = 50, value = 1),
                 numericInput(ns("pca.2d.2"), "2-D PCA: PC # (y-axis)", min = 1, max = 50, value = 2),
                 actionButton(ns("pca.filter.run"), "Confirm PCs for Next Step", class = "btn-success", style = "width: 100%")
               )
             )
           )
  )
}

mod.pca.server <- function(id, seurat.obj.qc, upstream.completed = reactive(TRUE)) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    completed <- reactiveVal(FALSE)

    # Soft guard: warn if the QC step isn't done, but let the user proceed
    pca.run.trigger <- reactiveVal(0)
    observeEvent(input$pca.run, {
      if (isTRUE(upstream.completed())) {
        pca.run.trigger(pca.run.trigger() + 1)
      } else {
        showModal(modalDialog(
          title = "Upstream step not completed",
          "The QC step hasn't been completed yet. Proceed anyway?",
          footer = tagList(
            modalButton("Cancel"),
            actionButton(ns("pca.run.proceed"), "Proceed anyway", class = "btn-warning")
          )
        ))
      }
    })
    observeEvent(input$pca.run.proceed, {
      removeModal()
      pca.run.trigger(pca.run.trigger() + 1)
    })

    output$pca.run.hint <- renderUI({
      if (!isTRUE(upstream.completed())) {
        tags$small(style = "color:#c0392b;", "QC step not completed yet.")
      }
    })

    # Logic to handle the PCA computation
    data.pca <- eventReactive(pca.run.trigger(), {
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
    }, ignoreInit = TRUE)

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

    # Rename Plot button after first use
    observeEvent(input$pca.plot.run, {
      updateActionButton(session, "pca.plot.run", label = "Update Plot")
    }, once = TRUE)

    observeEvent(input$pca.filter.run, {
      req(data.pca())
      completed(TRUE)
    })

    return(list(seurat.obj = data.pca, completed = completed, pca.dims = reactive(input$pca.number)))
  })
}