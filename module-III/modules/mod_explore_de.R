## ABOUTME: Differential expression module comparing configurable cell populations.
## ABOUTME: Supports cell type and metadata-based subsetting for ident.1/ident.2, with DT table and download.

mod.explore.de.ui <- function(id) {
  ns <- NS(id)
  tabPanel("Differential Expression",
    sidebarLayout(
      sidebarPanel(width = 4,
        # ident.1 (target)
        wellPanel(
          h4("Target Group (ident.1)"),
          selectInput(ns("ident1.celltypes"), "Cell Type(s):",
                      choices = NULL, multiple = TRUE),
          selectInput(ns("ident1.meta.col"), "Metadata Column:",
                      choices = c("None"), selected = "None"),
          conditionalPanel(
            condition = sprintf("input['%s'] != 'None'", ns("ident1.meta.col")),
            selectInput(ns("ident1.meta.vals"), "Select Value(s):",
                        choices = NULL, multiple = TRUE)
          )
        ),
        # ident.2 (baseline)
        wellPanel(
          h4("Baseline Group (ident.2)"),
          selectInput(ns("ident2.celltypes"), "Cell Type(s):",
                      choices = NULL, multiple = TRUE),
          selectInput(ns("ident2.meta.col"), "Metadata Column:",
                      choices = c("None"), selected = "None"),
          conditionalPanel(
            condition = sprintf("input['%s'] != 'None'", ns("ident2.meta.col")),
            selectInput(ns("ident2.meta.vals"), "Select Value(s):",
                        choices = NULL, multiple = TRUE)
          )
        ),
        actionButton(ns("run.de"), "Run DE Analysis",
                     class = "btn-success", style = "width:100%")
      ),
      mainPanel(width = 8,
        radioButtons(ns("de.direction"), "Show:",
                     choices = c("Both", "Upregulated", "Downregulated"),
                     selected = "Both", inline = TRUE),
        DT::DTOutput(ns("de.table")),
        br(),
        downloadButton(ns("download.de"), "Download Table (.xlsx)", class = "btn-success")
      )
    )
  )
}

mod.explore.de.server <- function(id, shared.data) {
  moduleServer(id, function(input, output, session) {

    de.result <- reactiveVal(NULL)

    # Populate cell type and metadata choices
    observe({
      req(shared.data())
      obj <- shared.data()
      ident.levels <- levels(obj)
      cat.cols <- get.categorical.meta(obj)

      updateSelectInput(session, "ident1.celltypes",
                        choices = ident.levels, selected = ident.levels[1])
      updateSelectInput(session, "ident2.celltypes",
                        choices = ident.levels, selected = ident.levels[1])
      updateSelectInput(session, "ident1.meta.col",
                        choices = c("None", cat.cols), selected = "None")
      updateSelectInput(session, "ident2.meta.col",
                        choices = c("None", cat.cols), selected = "None")
    })

    # Cascade: ident1 metadata column -> values
    observeEvent(input$ident1.meta.col, {
      req(shared.data(), input$ident1.meta.col != "None")
      vals <- unique(as.character(shared.data()@meta.data[[input$ident1.meta.col]]))
      updateSelectInput(session, "ident1.meta.vals", choices = vals, selected = vals)
    })

    # Cascade: ident2 metadata column -> values
    observeEvent(input$ident2.meta.col, {
      req(shared.data(), input$ident2.meta.col != "None")
      vals <- unique(as.character(shared.data()@meta.data[[input$ident2.meta.col]]))
      updateSelectInput(session, "ident2.meta.vals", choices = vals, selected = vals)
    })

    # Helper: get cell barcodes matching selection criteria
    get.cells <- function(obj, celltypes, meta.col, meta.vals) {
      # Filter by cell type
      cells <- WhichCells(obj, idents = celltypes)
      # Further filter by metadata if specified
      if (meta.col != "None" && length(meta.vals) > 0) {
        meta.mask <- obj@meta.data[[meta.col]] %in% meta.vals
        meta.cells <- colnames(obj)[meta.mask]
        cells <- intersect(cells, meta.cells)
      }
      cells
    }

    observeEvent(input$run.de, {
      req(shared.data())
      obj <- shared.data()

      cells.1 <- get.cells(obj, input$ident1.celltypes,
                           input$ident1.meta.col, input$ident1.meta.vals)
      cells.2 <- get.cells(obj, input$ident2.celltypes,
                           input$ident2.meta.col, input$ident2.meta.vals)

      validate(need(length(cells.1) >= 3,
                    paste0("Target group has too few cells (", length(cells.1), "). Need at least 3.")))
      validate(need(length(cells.2) >= 3,
                    paste0("Baseline group has too few cells (", length(cells.2), "). Need at least 3.")))
      validate(need(length(intersect(cells.1, cells.2)) == 0,
                    "Target and baseline groups overlap. Please adjust selections."))

      withProgress(message = "Running DE analysis...", value = 0.2, {
        tryCatch({
          markers <- FindMarkers(obj, ident.1 = cells.1, ident.2 = cells.2)
          incProgress(0.7, detail = "Filtering results")

          markers$gene <- rownames(markers)
          markers <- markers %>%
            filter(p_val_adj <= 0.05) %>%
            arrange(p_val_adj)

          de.result(markers)
          incProgress(0.1, detail = "Done")
          showNotification(paste(nrow(markers), "significant DE genes found."), type = "message")
        }, error = function(e) {
          showNotification(paste("DE error:", e$message), type = "error")
        })
      })
    })

    # Filtered table by direction
    de.filtered <- reactive({
      req(de.result())
      df <- de.result()
      if (input$de.direction == "Upregulated") {
        df <- df %>% filter(avg_log2FC > 0)
      } else if (input$de.direction == "Downregulated") {
        df <- df %>% filter(avg_log2FC < 0)
      }
      df
    })

    output$de.table <- DT::renderDT({
      req(de.filtered())
      DT::datatable(de.filtered(), rownames = FALSE,
                    options = list(pageLength = 10, scrollX = TRUE))
    })

    output$download.de <- downloadHandler(
      filename = function() { "de_results.xlsx" },
      content = function(file) {
        openxlsx::write.xlsx(de.filtered(), file)
      }
    )

    return(de.result)
  })
}
