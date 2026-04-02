## ABOUTME: Shiny module for cell type annotation using SingleR.
## ABOUTME: Accepts a reference SCE file, runs SingleR prediction, and stores labels in the Seurat object.

mod.annotation.ui <- function(id) {
  ns <- NS(id)

  tabPanel("Annotation",
           wellPanel(
             strong("Input Data"),
             fluidRow(
               column(6, fileInput(ns("data.upload"), "Upload Seurat Object (.rds)", accept = ".rds")),
               column(6, fileInput(ns("ref.upload"), "Upload Reference SCE (.rds)", accept = ".rds"))
             )
           ),
           wellPanel(
             strong("SingleR Annotation"),
             fluidRow(
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
             downloadButton(ns("download.annotation"), "Download Seurat Object")
           )
  )
}

mod.annotation.server <- function(id, seurat.obj.cellcycle) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    completed <- reactiveVal(FALSE)

    shinyjs::disable("annotation.run")

    # Load Seurat object from upload (temporary for testing)
    data.seurat <- reactive({
      req(input$data.upload)
      readRDS(input$data.upload$datapath)
    })

    # Load reference SCE from upload
    data.ref <- reactive({
      req(input$ref.upload)
      readRDS(input$ref.upload$datapath)
    })

    # Enable run button when both files are uploaded
    observe({
      if (!is.null(input$data.upload) && !is.null(input$ref.upload)) {
        shinyjs::enable("annotation.run")
      } else {
        shinyjs::disable("annotation.run")
      }
    })

    # Run SingleR annotation
    data.annotation <- eventReactive(input$annotation.run, {
      req(data.seurat(), data.ref())

      withProgress(message = "Running SingleR annotation...", value = 0.5, {
        srt <- data.seurat()
        ref <- data.ref()

        sce <- as.SingleCellExperiment(srt)
        ref.labels <- ref$celltype_level3

        pred <- SingleR(test = sce, ref = ref, labels = ref.labels)
        col.name <- paste0("SingleR.", input$label.name)
        srt[[col.name]] <- pred$labels
      })

      completed(TRUE)
      return(srt)
    })

    # Column name used for the most recent annotation run
    annotation.col <- reactiveVal(NULL)

    observeEvent(input$annotation.run, {
      annotation.col(paste0("SingleR.", input$label.name))
    })

    # UMAP colored by SingleR labels
    output$plot.annotation <- renderPlot({
      req(data.annotation(), annotation.col())
      DimPlot(data.annotation(), group.by = annotation.col(), reduction = "umap", label = TRUE) +
        theme(aspect.ratio = 1)
    }, res = 96)

    # Show annotation summary
    output$annotation.summary <- renderPrint({
      req(data.annotation(), annotation.col())
      table(data.annotation()[[annotation.col()]])
    })

    output$download.annotation <- downloadHandler(
      filename = function() {
        paste0("seurat_annotation_", Sys.Date(), ".rds")
      },
      content = function(file) {
        saveRDS(data.annotation(), file = file)
      }
    )

    return(list(seurat.obj = data.annotation, completed = completed))
  })
}
