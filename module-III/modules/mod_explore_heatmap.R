## ABOUTME: Heatmap module for scaled gene expression across cell types.
## ABOUTME: Supports HVG or custom gene lists, ScaleData with vars.to.regress, and cell sampling.

mod.explore.heatmap.ui <- function(id) {
  ns <- NS(id)
  tabPanel("Heatmap",
    sidebarLayout(
      sidebarPanel(width = 4,
        radioButtons(ns("heatmap.label"),  "Heatmap Label:",
                     choices = c("Cell Types" = "ct.label",
                                 "Other Metadata Column" = "other.label"),
                     selected = "ct.label"
                     ),
        selectInput(ns("select.idents"), "Cell Type(s):",
                    choices = NULL, multiple = TRUE),
        fluidRow(
          column(6, actionButton(ns("btn.select.all"), "Select All",
                                 class = "btn-info", style = "width:100%")),
          column(6, actionButton(ns("btn.clear.all"), "Clear",
                                 class = "btn-default", style = "width:100%"))
        ),
        selectInput(ns("subset.meta.col"), "Metadata Column:",
                    choices = c("None"), selected = "None"),
        conditionalPanel(
          condition = sprintf("input['%s'] != 'None'", ns("subset.meta.col")),
          selectInput(ns("subset.meta.vals"), "Select Value(s):",
                      choices = NULL, multiple = TRUE)
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
        # disabled max cells in heatmap
        #numericInput(ns("n.cells"), "Max Cells in Heatmap", value = 500, min = 50, max = 5000),
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
      updateSelectInput(session, "vars.to.regress",
                        choices = vars.to.regress.choices(obj))
      
      #updateSelectInput(session, "subset.meta.col",
                        #choices = c("None", split.by.choices(obj)), selected = "None")
      
      # ===== CHANGED: determine choices based on heatmap.label =====
      if (input$heatmap.label == "other.label") {
        subset.meta.choices <- split.by.choices(obj)       # No "None"
      } else {
        subset.meta.choices <- c("None", split.by.choices(obj))
      }
      
      updateSelectInput(
        session, "subset.meta.col",
        choices = subset.meta.choices,                     # CHANGED
        selected = subset.meta.choices[1]                  # CHANGED
      )
      # ===== END CHANGE =====
    })

    # Cascade: metadata column -> its values
    observeEvent(input$subset.meta.col, {
      req(shared.data(), 
          input$subset.meta.col, # added safety check
          input$subset.meta.col != "None")
      vals <- unique(as.character(shared.data()@meta.data[[input$subset.meta.col]]))
      updateSelectInput(session, "subset.meta.vals", choices = vals, selected = vals)
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
      

      if (input$subset.meta.col != "None" && length(input$subset.meta.vals) > 0) {
        meta.mask <- obj@meta.data[[input$subset.meta.col]] %in% input$subset.meta.vals
        cells.keep <- colnames(obj)[meta.mask]
        validate(need(length(cells.keep) > 0, "No cells match the selected metadata values."))
        obj <- subset(obj, cells = cells.keep)
      }

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

        # ------ disabled max cells in heatmap ---------------------
        # Sample cells
        #n.cells <- min(input$n.cells, ncol(obj))
        #sampled.cells <- sample(colnames(obj), n.cells)
        # ------ End -----------------------------------------------

        incProgress(0.4, detail = "Rendering heatmap")
        # ------ Using pheatmap instead of the seurat --------------
        #DoHeatmap(obj, features = features, cells = sampled.cells,
                  #size = 4, angle = 90)
        # get the scaled data
        mat <- GetAssayData(
          obj,
          assay = DefaultAssay(obj),
          layer = "scale.data"
        )
        
        # Order cells by group
        cell.label <- if (input$heatmap.label == "ct.label") {
          obj$manual_annotation
        }else{
          obj[[input$subset.meta.col]][,1]
        }
        
        cell.order <- order(cell.label)
        mat <- mat[, cell.order, drop = FALSE]
        
        
        # Column annotation
        annotation_col <- data.frame(Label = cell.label[cell.order])
        rownames(annotation_col) <- colnames(mat)
        
        pheatmap(
          mat = mat,
          scale = "none",                # Scale only selected genes for this plot
          cluster_rows = TRUE,           # Hierarchical clustering of genes (Y-axis)
          cluster_cols = FALSE,          # Preserve cell grouping on the X-axis
          annotation_col = annotation_col,
          show_colnames = FALSE,
          fontsize_row = 9,
          border_color = NA,
          color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
          breaks = seq(-2.5, 2.5, length.out = 101),
        )
        
        
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
