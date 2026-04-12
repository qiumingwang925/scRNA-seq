## ABOUTME: Shiny sub-module for interactive manual cell type annotation on UMAP plots.
## ABOUTME: Provides lasso selection, sub-clustering, biomarker expression, and DE analysis tools.

mod.annotation.manual.ui <- function(id) {
  ns <- NS(id)

  tagList(
    sidebarLayout(
      sidebarPanel(
        width = 4,

        radioButtons(ns("display.mode"), "Display UMAP from:",
                     choices = c("All cells" = "full",
                                 "Selected cells (subset)" = "subset"),
                     selected = "full"),
        numericInput(ns("subset.pcs"), "PCs for Subset UMAP:",
                     value = 10, min = 2, max = 50, step = 1),
        actionButton(ns("run.subset.umap"), "Run UMAP on Selection",
                     class = "btn-warning", style = "width:100%"),

        hr(),

        tabsetPanel(
          id = ns("viz.mode"),
          tabPanel("Metadata",
            br(),
            selectInput(ns("view.meta"), "Color by Metadata:", choices = NULL),
            hr(),
            numericInput(ns("cluster.res"), "Cluster Resolution:",
                         value = 0.1, min = 0, max = 2, step = 0.05),
            actionButton(ns("run.clustering"), "Run Clustering",
                         class = "btn-info", style = "width:100%")
          ),
          tabPanel("Expression",
            br(),
            fileInput(ns("marker.file"), "Upload Biomarker CSV", accept = ".csv"),
            selectInput(ns("gene.major"), "Select Major Cell Type:", choices = NULL),
            selectInput(ns("gene.fine"), "Select Fine Cell Type:", choices = NULL),
            selectInput(ns("gene.selected"), "Select Gene/Marker:", choices = NULL),
            hr(),
            textInput(ns("gene.entered"), "OR Type Gene Name:",
                      value = "", placeholder = "e.g., Cd68"),
            helpText("Typing a gene here will override the dropdown selection.")
          )
        ),

        hr(),

        h4("Differential Expression"),
        tags$small("Identify markers for the current selection vs. the rest."),
        actionButton(ns("run.de"), "Run DE Analysis",
                     class = "btn-success", style = "width:100%"),

        hr(),

        h4("Manual Annotation"),
        tags$small("Use Lasso/Box tool on UMAP to select cells."),
        textInput(ns("new.label"), "Enter New Label:",
                  placeholder = "e.g., Alveolar Macrophage"),
        actionButton(ns("apply.label"), "Apply Label to Selection",
                     class = "btn-primary", style = "width:100%"),
        br(), br(),
        mod.save.config.ui(ns("save"), label = "Export Annotated Object")
      ),

      mainPanel(
        width = 8,
        tabsetPanel(
          id = ns("main.tabs"),
          tabPanel("UMAP Visualization",
            br(),
            plotly::plotlyOutput(ns("umap.main"), height = "600px")
          ),
          tabPanel("DE Analysis",
            br(),
            DT::DTOutput(ns("de.table"))
          )
        )
      )
    )
  )
}

mod.annotation.manual.server <- function(id, current.obj) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    completed <- reactiveVal(FALSE)
    subset.obj <- reactiveVal(NULL)

    # --- Update metadata dropdown based on active object ---
    observe({
      obj <- if (input$display.mode == "subset") subset.obj() else current.obj()
      req(obj)

      meta.cols <- colnames(obj@meta.data)

      # Priority ordering: manual_annotation, SingleR columns, seurat_clusters, then rest
      singler.cols <- sort(meta.cols[grepl("^SingleR\\.", meta.cols)])
      priority <- c("manual_annotation", singler.cols)
      if ("seurat_clusters" %in% meta.cols) priority <- c(priority, "seurat_clusters")
      other <- setdiff(meta.cols, priority)
      ordered.choices <- c(intersect(priority, meta.cols), other)

      current.sel <- input$view.meta
      selected <- if (!is.null(current.sel) && current.sel %in% ordered.choices) {
        current.sel
      } else {
        ordered.choices[1]
      }

      updateSelectInput(session, "view.meta",
                        choices = ordered.choices, selected = selected)
    })

    # --- Biomarker CSV upload and cascading dropdowns ---
    marker.df <- reactive({
      req(input$marker.file)
      read.csv(input$marker.file$datapath)
    })

    observeEvent(marker.df(), {
      df <- marker.df()
      updateSelectInput(session, "gene.major",
                        choices = c("", unique(df$Label.main)))
    })

    observeEvent(input$gene.major, {
      req(input$gene.major != "", marker.df())
      df <- marker.df()
      sub.df <- df[df$Label.main == input$gene.major, ]
      updateSelectInput(session, "gene.fine",
                        choices = c("", unique(sub.df$Label.fine)))
    })

    observeEvent(input$gene.fine, {
      req(input$gene.fine != "", marker.df())
      df <- marker.df()
      sub.df <- df[df$Label.main == input$gene.major &
                     df$Label.fine == input$gene.fine, ]
      updateSelectInput(session, "gene.selected",
                        choices = c("", unique(sub.df$Markers)))
    })

    # --- Subset UMAP on selected cells ---
    observeEvent(input$run.subset.umap, {
      source.obj <- if (input$display.mode == "subset") subset.obj() else current.obj()
      req(source.obj)

      sel <- plotly::event_data("plotly_selected", source = ns("umap.select"))
      validate(need(!is.null(sel),
                    "Please lasso select cells on the UMAP first."))

      selected.barcodes <- unique(as.character(sel$customdata))
      req(length(selected.barcodes) >= 3)

      withProgress(message = "Calculating Subset UMAP...", value = 0, {
        tryCatch({
          sub <- subset(source.obj, cells = selected.barcodes)
          # Cap PCs to what the data supports: must be < min(cells, features)
          max.pcs <- min(ncol(sub), nrow(sub)) - 1
          n.pcs <- min(input$subset.pcs, max.pcs)
          sub <- RunPCA(sub, npcs = n.pcs)
          sub <- FindNeighbors(sub, dims = 1:n.pcs)
          sub <- RunUMAP(sub, dims = 1:n.pcs)

          subset.obj(sub)
          updateRadioButtons(session, "display.mode", selected = "subset")
          showNotification("Subset UMAP calculated!", type = "message")
        }, error = function(e) {
          showNotification(paste("Error:", e$message),
                           type = "error", duration = 10)
        })
      })
    })

    # --- Clustering ---
    observeEvent(input$run.clustering, {
      res <- input$cluster.res
      is.full <- input$display.mode == "full"
      obj <- if (is.full) current.obj() else subset.obj()
      req(obj)

      tryCatch({
        withProgress(message = "Running Clustering...", value = 0, {
          obj <- FindClusters(obj, resolution = res)

          # Keep only seurat_clusters, remove intermediate resolution columns
          extra.cols <- grepl("_res\\.", colnames(obj@meta.data))
          obj@meta.data <- obj@meta.data[, !extra.cols, drop = FALSE]

          if (is.full) current.obj(obj) else subset.obj(obj)
          showNotification("Clustering complete!", type = "message")
        })
      }, error = function(e) {
        showNotification(paste("Clustering error:", e$message),
                         type = "error", duration = 10)
      })
    })

    # --- UMAP plot ---
    output$umap.main <- plotly::renderPlotly({
      if (input$display.mode == "subset") {
        validate(need(!is.null(subset.obj()),
                      "No subset calculated. Select cells and click 'Run UMAP on Selection'."))
        obj <- subset.obj()
      } else {
        req(current.obj())
        obj <- current.obj()
      }

      req("umap" %in% names(obj@reductions))

      if (input$viz.mode == "Expression") {
        gene <- if (nzchar(input$gene.entered)) input$gene.entered else input$gene.selected
        req(gene)
        validate(need(gene %in% rownames(obj),
                      paste0("Gene '", gene, "' not found. Check spelling/case.")))
        p <- FeaturePlot(obj, features = gene, reduction = "umap")
      } else {
        req(input$view.meta)
        validate(need(input$view.meta %in% colnames(obj@meta.data),
                      "Metadata column not found."))
        meta.vals <- obj[[input$view.meta]][, 1]
        if (is.numeric(meta.vals)) {
          p <- FeaturePlot(obj, reduction = "umap", features = input$view.meta)
        } else {
          p <- DimPlot(obj, reduction = "umap", group.by = input$view.meta)
        }
      }

      # Attach barcodes for plotly lasso selection
      p$layers[[1]]$mapping$customdata <- ggplot2::aes(
        label = rownames(p$data))$label

      plotly::ggplotly(p + theme_minimal(), source = ns("umap.select")) %>%
        plotly::layout(dragmode = "lasso")
    })

    # --- Apply annotation label ---
    observeEvent(input$apply.label, {
      sel <- plotly::event_data("plotly_selected", source = ns("umap.select"))
      req(sel, input$new.label != "")

      selected.barcodes <- unique(as.character(sel$customdata))
      new.label <- input$new.label

      # Update master object
      master <- isolate(current.obj())
      if (!"manual_annotation" %in% colnames(master@meta.data)) {
        master$manual_annotation <- "Unlabeled"
      }
      master@meta.data[selected.barcodes, "manual_annotation"] <- new.label
      current.obj(master)

      # Sync to subset if it exists
      if (!is.null(subset.obj())) {
        sub <- isolate(subset.obj())
        if (!"manual_annotation" %in% colnames(sub@meta.data)) {
          sub$manual_annotation <- "Unlabeled"
        }
        cells.in.sub <- intersect(selected.barcodes, colnames(sub))
        if (length(cells.in.sub) > 0) {
          sub@meta.data[cells.in.sub, "manual_annotation"] <- new.label
          subset.obj(sub)
        }
      }

      completed(TRUE)
      showNotification(paste("Applied label:", new.label), type = "message")
    })

    # --- DE analysis ---
    de.results <- eventReactive(input$run.de, {
      sel <- plotly::event_data("plotly_selected", source = ns("umap.select"))
      validate(need(!is.null(sel),
                    "Please select cells on the UMAP first."))

      selected.cells <- unique(as.character(sel$customdata))

      universe <- if (input$display.mode == "full") current.obj() else subset.obj()
      req(universe)

      all.cells <- colnames(universe)
      selected.cells <- intersect(selected.cells, all.cells)
      rest.cells <- setdiff(all.cells, selected.cells)

      validate(
        need(length(selected.cells) > 3, "Select at least 4 cells for DE."),
        need(length(rest.cells) > 3, "Not enough background cells.")
      )

      withProgress(message = "Running DE Analysis...", value = 0, {
        tryCatch({
          results <- Seurat::FindMarkers(
            object = universe,
            ident.1 = selected.cells,
            ident.2 = rest.cells,
            only.pos = TRUE,
            min.pct = 0.1,
            logfc.threshold = 0.25
          )

          if (nrow(results) == 0) {
            showNotification("No significant markers found.", type = "warning")
            return(NULL)
          }

          showNotification("DE Analysis complete!", type = "message")
          cbind(Gene = rownames(results), results)
        }, error = function(e) {
          showNotification(paste("DE error:", e$message), type = "error")
          NULL
        })
      })
    })

    output$de.table <- DT::renderDT({
      req(de.results())
      DT::datatable(de.results(), rownames = FALSE,
                    options = list(pageLength = 10, scrollX = TRUE),
                    selection = "single")
    })

    # Switch to DE tab after running
    observeEvent(input$run.de, {
      updateTabsetPanel(session, "main.tabs", selected = "DE Analysis")
    })

    # --- Export ---
    mod.save.config.server("save", current.obj)

    return(list(completed = completed))
  })
}
