## ABOUTME: Upload module for Module IV. Loads annotated Seurat object, joins layers,
## ABOUTME: validates the normalized data slot, and assigns a user-chosen metadata column as "group".

mod.interact.upload.ui <- function(id) {
  ns <- NS(id)
  tabPanel("Upload Data",
    sidebarLayout(
      sidebarPanel(width = 4,
        fileInput(ns("file.in"), "Upload Seurat Object (.rds)", accept = ".rds"),
        hr(),
        uiOutput(ns("group.selector.ui")),
        helpText("Upload an annotated Seurat object with Idents set to cell types.",
                 "Pick the metadata column that labels conditions / treatments.",
                 "At least two levels are required for cross-group comparison.")
      ),
      mainPanel(width = 8,
        h4("Object Summary"),
        tableOutput(ns("report.table")),
        h4("Validation"),
        verbatimTextOutput(ns("validation.msg"))
      )
    )
  )
}

mod.interact.upload.server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    raw.obj <- reactive({
      req(input$file.in)
      withProgress(message = "Loading Seurat object...", value = 0.3, {
        tryCatch({
          obj <- readRDS(input$file.in$datapath)
          incProgress(0.7, detail = "Done")
          obj
        }, error = function(e) {
          showNotification(paste("Error loading file:", e$message), type = "error")
          NULL
        })
      })
    })

    report.data <- reactive({
      obj <- raw.obj()
      req(obj)
      data.frame(
        Feature = c("Assays", "Metadata Columns", "Active Ident Levels", "Cells", "Features"),
        Value = c(
          paste(Assays(obj), collapse = ", "),
          paste(names(obj@meta.data), collapse = ", "),
          paste(levels(Idents(obj)), collapse = ", "),
          as.character(ncol(obj)),
          as.character(nrow(obj))
        ),
        stringsAsFactors = FALSE
      )
    })

    output$report.table <- renderTable({ report.data() })

    output$group.selector.ui <- renderUI({
      obj <- raw.obj()
      req(obj)
      cat.cols <- get.categorical.meta(obj)
      default <- if ("group" %in% cat.cols) "group" else cat.cols[1]
      tagList(
        selectInput(ns("group.col"), "Condition / Group column:",
                    choices = cat.cols, selected = default),
        actionButton(ns("confirm"), "Confirm & Prepare Object",
                     class = "btn-success", style = "width:100%")
      )
    })

    processed.obj <- eventReactive(input$confirm, {
      obj <- raw.obj()
      req(obj, input$group.col)

      withProgress(message = "Preparing object...", value = 0.2, {
        main.assay <- ifelse("SCT" %in% Assays(obj), "SCT", "RNA")
        DefaultAssay(obj) <- main.assay

        incProgress(0.2, detail = "Joining layers...")
        obj <- tryCatch({
          JoinLayers(obj, assay = "RNA")
        }, error = function(e) {
          message("JoinLayers skipped: ", e$message)
          obj
        })

        incProgress(0.2, detail = "Validating normalized data...")
        data.mat <- GetAssayData(obj, assay = main.assay, layer = "data")
        if (is.null(data.mat) || nrow(data.mat) == 0 || ncol(data.mat) == 0) {
          showNotification(
            paste0("Assay '", main.assay, "' has an empty 'data' layer. ",
                   "CellChat requires log-normalized expression data."),
            type = "error", duration = NULL
          )
          return(NULL)
        }

        incProgress(0.2, detail = "Assigning group column...")
        obj$group <- as.character(obj@meta.data[[input$group.col]])

        n.groups <- length(unique(obj$group))
        if (n.groups < 2) {
          showNotification(
            paste0("Column '", input$group.col, "' has only ", n.groups,
                   " unique value(s). At least 2 are required for CellChat comparison."),
            type = "error", duration = NULL
          )
          return(NULL)
        }

        incProgress(0.2, detail = "Done")
        obj
      })
    })

    output$validation.msg <- renderText({
      obj <- processed.obj()
      req(obj)
      paste0(
        "Default assay: ", DefaultAssay(obj), "\n",
        "Cells: ", ncol(obj), "\n",
        "Cell-type idents: ", length(levels(Idents(obj))), "\n",
        "Group column: ", input$group.col, " (copied to 'group')\n",
        "Group levels: ", paste(sort(unique(obj$group)), collapse = ", ")
      )
    })

    return(processed.obj)
  })
}
