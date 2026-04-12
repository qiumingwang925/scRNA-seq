## ABOUTME: Shiny sub-module for automated cell type annotation using SingleR.
## ABOUTME: Runs SingleR prediction against a reference SCE and stores labels in the shared Seurat object.

mod.annotation.singler.ui <- function(id) {
  ns <- NS(id)

  tagList(
    wellPanel(
      strong("SingleR Annotation"),
      fluidRow(
        column(6, fileInput(ns("ref.upload"), "Upload Reference SCE (.rds)", accept = ".rds")),
        column(3, textInput(ns("label.name"), "Label Name", value = "labels")),
        column(3, actionButton(ns("annotation.run"), "Run SingleR", class = "btn-success", style = "margin-top: 25px"))
      )
    ),
    fluidRow(
      column(8,
        wellPanel(
          strong("UMAP by Annotation"),
          plotOutput(ns("plot.annotation"), height = "500px", width = "100%")
        )
      ),
      column(4,
        wellPanel(
          strong("Annotation Summary"),
          verbatimTextOutput(ns("annotation.summary"))
        )
      )
    ),
    wellPanel(
      mod.save.config.ui(ns("save"), label = "Download Seurat Object")
    )
  )
}

mod.annotation.singler.server <- function(id, current.obj) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    completed <- reactiveVal(FALSE)

    shinyjs::disable("annotation.run")

    # Load reference SCE from upload
    data.ref <- reactive({
      req(input$ref.upload)
      readRDS(input$ref.upload$datapath)
    })

    # Enable run button when reference is uploaded and data is available
    observe({
      if (!is.null(current.obj()) && !is.null(input$ref.upload)) {
        shinyjs::enable("annotation.run")
      } else {
        shinyjs::disable("annotation.run")
      }
    })

    # Track the column name for the most recent annotation
    annotation.col <- reactiveVal(NULL)

    # Run SingleR annotation
    observeEvent(input$annotation.run, {
      req(current.obj(), data.ref())

      withProgress(message = "Running SingleR annotation...", value = 0.5, {
        srt <- current.obj()
        ref <- data.ref()

        sce <- as.SingleCellExperiment(srt)
        ref.labels <- ref$celltype_level3

        pred <- SingleR(test = sce, ref = ref, labels = ref.labels)
        col.name <- paste0("SingleR.", input$label.name)
        srt[[col.name]] <- pred$labels

        current.obj(srt)
        annotation.col(col.name)
      })

      completed(TRUE)
      showNotification("SingleR annotation complete!", type = "message")
    })

    # UMAP colored by SingleR labels
    output$plot.annotation <- renderPlot({
      req(current.obj(), annotation.col())
      DimPlot(current.obj(), group.by = annotation.col(), reduction = "umap", label = TRUE) +
        theme(aspect.ratio = 1)
    }, res = 96)

    # Annotation summary
    output$annotation.summary <- renderPrint({
      req(current.obj(), annotation.col())
      table(current.obj()[[annotation.col()]])
    })

    # Export with save configuration
    mod.save.config.server("save", current.obj)

    return(list(completed = completed))
  })
}
