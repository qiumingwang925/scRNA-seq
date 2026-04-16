## ABOUTME: Upload module for loading pre-processed Seurat .rds objects.
## ABOUTME: Validates object structure and reports available assays, reductions, and metadata.

mod.explore.upload.ui <- function(id) {
  ns <- NS(id)
  tabPanel("Upload Data",
    sidebarLayout(
      sidebarPanel(width = 4,
        fileInput(ns("file.in"), "Upload Seurat Object (.rds)", accept = ".rds"),
        hr(),
        helpText("Ensure object is slimmed down: SCT/RNA 'data' slot only,",
                 "PCA + UMAP reductions, matching graph slot,",
                 "and active.ident set to cell type labels.")
      ),
      mainPanel(width = 8,
        h4("Object Summary Report"),
        tableOutput(ns("report.table")),
        verbatimTextOutput(ns("validation.msg"))
      )
    )
  )
}

mod.explore.upload.server <- function(id) {
  moduleServer(id, function(input, output, session) {

    raw.obj <- reactive({
      req(input$file.in)
      withProgress(message = "Loading Seurat object...", value = 0.3, {
        tryCatch({
          obj <- readRDS(input$file.in$datapath)
          incProgress(0.7, detail = "Done")
          obj
        }, error = function(e) {
          message("Upload error: ", e$message)
          showNotification(paste("Error loading file:", e$message), type = "error")
          NULL
        })
      })
    })

    report.data <- reactive({
      obj <- raw.obj()
      req(obj)
      data.frame(
        Feature = c("Assays", "Metadata Columns", "Reductions", "Graphs", "Active Ident"),
        Status = c(
          paste(Assays(obj), collapse = ", "),
          paste(names(obj@meta.data), collapse = ", "),
          paste(Reductions(obj), collapse = ", "),
          paste(Graphs(obj), collapse = ", "),
          paste(unique(Idents(obj)), collapse = ", ")
        ),
        stringsAsFactors = FALSE
      )
    })

    output$report.table <- renderTable({ report.data() })

    processed.obj <- reactive({
      obj <- raw.obj()
      req(obj)

      main.assay <- ifelse("SCT" %in% Assays(obj), "SCT", "RNA")
      DefaultAssay(obj) <- main.assay

      output$validation.msg <- renderText({
        paste0("Default assay set to: ", main.assay,
               "\nCells: ", ncol(obj),
               "\nFeatures: ", nrow(obj))
      })

      obj
    })

    return(processed.obj)
  })
}
