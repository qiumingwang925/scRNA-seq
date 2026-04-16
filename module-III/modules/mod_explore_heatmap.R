## ABOUTME: Heatmap module for scaled gene expression across cell types.
## ABOUTME: Supports HVG or custom gene lists, ScaleData with vars.to.regress, and cell sampling.

mod.explore.heatmap.ui <- function(id) {
  ns <- NS(id)
  tabPanel("Heatmap",
    sidebarLayout(
      sidebarPanel(width = 4,
        selectInput(ns("select.idents"), "Cell Type(s):",
                    choices = NULL, multiple = TRUE),
        fluidRow(
          column(6, actionButton(ns("btn.select.all"), "Select All",
                                 class = "btn-info", style = "width:100%")),
          column(6, actionButton(ns("btn.clear.all"), "Clear",
                                 class = "btn-default", style = "width:100%"))
        ),
        hr(),
        radioButtons(ns("gene.mode"), "Gene Selection:",
                     choices = c("Highly Variable Genes" = "hvg",
                                 "Custom Gene List" = "custom"),
                     selected = "hvg"),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'hvg'", ns("gene.mode")),
          numericInput(ns("top.n"), "Top N Variable Genes", value = 20, min = 5, max = 200)
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'custom'", ns("gene.mode")),
          textInput(ns("gene.input"), "Gene(s) (comma-separated)",
                    placeholder = "e.g. Cd68, Cx3cr1, Ccr2")
        ),
        hr(),
        selectInput(ns("vars.to.regress"), "Variables to Regress (optional):",
                    choices = NULL, multiple = TRUE),
        numericInput(ns("n.cells"), "Max Cells in Heatmap", value = 500, min = 50, max = 5000),
        actionButton(ns("run.heatmap"), "Generate Heatmap",
                     class = "btn-success", style = "width:100%"),
        hr(),
        h4("Download Figure"),
        numericInput(ns("fig.w"), "Width (inches)", value = 10, min = 2, max = 30),
        numericInput(ns("fig.h"), "Height (inches)", value = 8, min = 2, max = 30),
        downloadButton(ns("download.heatmap"), "Download Figure", class = "btn-success")
      ),
      mainPanel(width = 8,
        plotOutput(ns("plot.heatmap"), height = "800px")
      )
    )
  )
}

mod.explore.heatmap.server <- function(id, shared.data) {
  moduleServer(id, function(input, output, session) {

    observe({
      req(shared.data())
      obj <- shared.data()
      ident.levels <- levels(obj)
      updateSelectInput(session, "select.idents",
                        choices = ident.levels, selected = ident.levels)
      # Populate vars.to.regress with numeric metadata columns
      meta <- obj@meta.data
      num.cols <- names(meta)[sapply(meta, is.numeric)]
      updateSelectInput(session, "vars.to.regress", choices = num.cols)
    })

    observeEvent(input$btn.select.all, {
      req(shared.data())
      updateSelectInput(session, "select.idents",
                        choices = levels(shared.data()),
                        selected = levels(shared.data()))
    })

    observeEvent(input$btn.clear.all, {
      req(shared.data())
      updateSelectInput(session, "select.idents",
                        choices = levels(shared.data()),
                        selected = character(0))
    })

    plot.heatmap <- eventReactive(input$run.heatmap, {
      req(shared.data())
      obj <- shared.data()

      idents.selected <- input$select.idents
      validate(need(length(idents.selected) > 0, "Please select at least one cell type."))
      obj <- subset(obj, idents = idents.selected)

      withProgress(message = "Generating heatmap...", value = 0, {

        # Determine features
        if (input$gene.mode == "hvg") {
          incProgress(0.2, detail = "Finding variable features")
          # Use whichever layer is available (slim objects may only have "data", not "counts")
          available.layers <- Layers(obj)
          hvg.layer <- if ("counts" %in% available.layers) "counts" else "data"
          tryCatch({
            obj <- FindVariableFeatures(obj, layer = hvg.layer)
          }, error = function(e) {
            message("FindVariableFeatures error: ", e$message)
            # Fallback: try without specifying layer
            obj <<- FindVariableFeatures(obj)
          })
          validate(need(length(VariableFeatures(obj)) > 0,
                        "No variable features found. Try using Custom Gene List instead."))
          features <- VariableFeatures(obj)[1:min(input$top.n, length(VariableFeatures(obj)))]
        } else {
          features <- trimws(unlist(strsplit(input$gene.input, ",")))
          features <- features[nchar(features) > 0]
          validate(need(length(features) > 0, "Please enter at least one gene name."))
          missing <- features[!features %in% rownames(obj)]
          validate(need(length(missing) == 0,
                        paste0("Gene(s) not found: ", paste(missing, collapse = ", "))))
        }

        # Scale data
        incProgress(0.3, detail = "Scaling data")
        vars.regress <- if (length(input$vars.to.regress) == 0) NULL else input$vars.to.regress
        tryCatch({
          obj <- ScaleData(obj, features = features, vars.to.regress = vars.regress)
        }, error = function(e) {
          message("ScaleData error: ", e$message)
          showNotification(paste("ScaleData error:", e$message), type = "error")
          validate(need(FALSE, "ScaleData failed. Check vars.to.regress selection."))
        })

        # Sample cells
        n.cells <- min(input$n.cells, ncol(obj))
        sampled.cells <- sample(colnames(obj), n.cells)

        incProgress(0.4, detail = "Rendering heatmap")
        DoHeatmap(obj, features = features, cells = sampled.cells,
                  size = 4, angle = 90)
      })
    })

    output$plot.heatmap <- renderPlot({ plot.heatmap() }, res = 96)

    output$download.heatmap <- downloadHandler(
      filename = function() { "heatmap.pdf" },
      content = function(file) {
        ggsave(file, plot = plot.heatmap(),
               width = input$fig.w, height = input$fig.h)
      }
    )
  })
}
