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
        column(3,
               actionButton(ns("annotation.run"), "Run SingleR", class = "btn-success", style = "margin-top: 25px"),
               uiOutput(ns("annotation.run.hint")))
      )
    ),
    fluidRow(
      column(12,
        wellPanel(
          strong("Annotation Summary"),
          verbatimTextOutput(ns("annotation.summary")),
          helpText("To view the predicted cell types on a UMAP, open the Manual Annotation tab and colour by the SingleR annotation column.")
        )
      )
    ),
    wellPanel(
      mod.save.config.ui(ns("save"), label = "Download Seurat Object")
    )
  )
}

mod.annotation.singler.server <- function(id, current.obj, upstream.completed = reactive(TRUE)) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    completed <- reactiveVal(FALSE)

    # Load reference SCE from upload
    data.ref <- reactive({
      req(input$ref.upload)
      readRDS(input$ref.upload$datapath)
    })

    # Soft guard: warn if the Cell Cycle step isn't done, but let the user proceed
    annotation.run.trigger <- reactiveVal(0)
    observeEvent(input$annotation.run, {
      if (isTRUE(upstream.completed())) {
        annotation.run.trigger(annotation.run.trigger() + 1)
      } else {
        showModal(modalDialog(
          title = "Upstream step not completed",
          "The Cell Cycle step hasn't been completed yet. Proceed anyway?",
          footer = tagList(
            modalButton("Cancel"),
            actionButton(ns("annotation.run.proceed"), "Proceed anyway", class = "btn-warning")
          )
        ))
      }
    })
    observeEvent(input$annotation.run.proceed, {
      removeModal()
      annotation.run.trigger(annotation.run.trigger() + 1)
    })

    output$annotation.run.hint <- renderUI({
      if (!isTRUE(upstream.completed())) {
        tags$small(style = "color:#c0392b;", "Cell Cycle step not completed yet.")
      }
    })

    # Track the column name for the most recent annotation
    annotation.col <- reactiveVal(NULL)

    # Run SingleR annotation
    observeEvent(annotation.run.trigger(), {
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
    }, ignoreInit = TRUE)

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
