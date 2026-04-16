## ABOUTME: UMAP cell-type visualization with interactive selection, highlighting, and gene expression.
## ABOUTME: Three sub-tabs: interactive lasso/re-UMAP, high-res figure generation, and gene (co-)expression plots.

# Helper: detect best UMAP reduction and SNN graph for a given Seurat object
detect.slots.for <- function(obj) {
  red.name <- grep("umap", Reductions(obj), value = TRUE, ignore.case = TRUE)
  # Prefer umap_subset if it exists (from re-UMAP), otherwise first umap match
  if (length(red.name) == 0) {
    red.name <- "umap"
  } else if ("umap_subset" %in% red.name) {
    red.name <- "umap_subset"
  } else {
    red.name <- red.name[1]
  }
  graph.name <- grep("snn", Graphs(obj), value = TRUE, ignore.case = TRUE)[1]
  list(reduction = red.name, graph = graph.name)
}

# Helper: FeaturePlot with colored points on top and larger than grey points
reorder.feature.plot <- function(obj, features, reduction, cells = NULL,
                                 blend = FALSE, pt.size.bg = 0.5, pt.size.fg = 1.5) {
  if (blend) {
    p <- FeaturePlot(obj, features = features, reduction = reduction,
                     cells = cells, blend = TRUE, order = TRUE)
    # blend returns a list of plots; adjust point sizes on each
    if (is.list(p) && !inherits(p, "gg")) {
      p <- lapply(p, function(pp) {
        pp + theme_minimal()
      })
      return(patchwork::wrap_plots(p))
    }
    return(p + theme_minimal())
  }

  # Single gene: build plot with grey background layer + colored foreground layer
  p <- FeaturePlot(obj, features = features, reduction = reduction,
                   cells = cells, order = TRUE)
  # Reorder layers: grey cells as small background, expressing cells as larger foreground
  plot.data <- p$data
  expr.col <- features[1]
  if (!expr.col %in% colnames(plot.data)) {
    return(p + theme_minimal())
  }
  bg <- plot.data[plot.data[[expr.col]] == 0, ]
  fg <- plot.data[plot.data[[expr.col]] > 0, ]
  x.col <- colnames(plot.data)[1]
  y.col <- colnames(plot.data)[2]

  ggplot() +
    geom_point(data = bg, aes(x = .data[[x.col]], y = .data[[y.col]]),
               color = "lightgrey", size = pt.size.bg) +
    geom_point(data = fg, aes(x = .data[[x.col]], y = .data[[y.col]],
                               color = .data[[expr.col]]),
               size = pt.size.fg) +
    scale_color_gradientn(colors = c("lightgrey", "blue")) +
    labs(x = x.col, y = y.col, color = expr.col) +
    theme_minimal()
}

mod.explore.umap.ui <- function(id) {
  ns <- NS(id)
  tabPanel("UMAP",
    tabsetPanel(
      # Sub-tab 1: Interactive Selection
      tabPanel("Interactive Selection",
        sidebarLayout(
          sidebarPanel(width = 4,
            radioButtons(ns("display.mode"), "Display UMAP from:",
                         choices = c("All cells" = "full",
                                     "Selected cells (subset)" = "subset"),
                         selected = "full"),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'subset'", ns("display.mode")),
              helpText("1. Lasso select cells in 'All cells' mode."),
              helpText("2. Switch to 'Selected cells' mode."),
              numericInput(ns("subset.pcs"), "PCs for re-UMAP", value = 30, min = 5, max = 100),
              actionButton(ns("run.subset.umap"), "Run UMAP on Selection",
                           class = "btn-warning", style = "width:100%"),
              hr()
            ),
            actionButton(ns("reset.umap"), "Reset Everything",
                         class = "btn-danger", style = "width:100%"),
            br(), br(),
            verbatimTextOutput(ns("cell.stat.box"))
          ),
          mainPanel(width = 8,
            plotlyOutput(ns("umap.interactive"), height = "600px")
          )
        )
      ),
      # Sub-tab 2: Highlight View
      tabPanel("Highlight View",
        sidebarLayout(
          sidebarPanel(width = 4,
            selectInput(ns("select.idents"), "Select Cell Type(s) to Highlight:",
                        choices = NULL, multiple = TRUE),
            fluidRow(
              column(6, actionButton(ns("btn.select.all"), "Select All",
                                     class = "btn-info", style = "width:100%")),
              column(6, actionButton(ns("btn.clear.all"), "Clear",
                                     class = "btn-default", style = "width:100%"))
            ),
            hr(),
            h4("Download Figure"),
            numericInput(ns("fig.highlight.w"), "Width (inches)", value = 8, min = 2, max = 20),
            numericInput(ns("fig.highlight.h"), "Height (inches)", value = 6, min = 2, max = 20),
            downloadButton(ns("download.highlight"), "Download Figure", class = "btn-success")
          ),
          mainPanel(width = 8,
            plotOutput(ns("umap.static"), height = "600px")
          )
        )
      ),
      # Sub-tab 3: Gene Expression
      tabPanel("Gene Expression",
        sidebarLayout(
          sidebarPanel(width = 4,
            radioButtons(ns("expr.mode"), "Mode:",
                         choices = c("Expression" = "single",
                                     "Co-expression" = "blend"),
                         selected = "single"),
            # Single gene expression controls
            conditionalPanel(
              condition = sprintf("input['%s'] == 'single'", ns("expr.mode")),
              textInput(ns("gene.name"), "Gene Name", placeholder = "e.g. Cd68"),
              selectInput(ns("split.by"), "Split By (optional):",
                          choices = c("None"), selected = "None"),
              numericInput(ns("grid.ncol"), "Grid Columns", value = 2, min = 1, max = 6)
            ),
            # Co-expression controls
            conditionalPanel(
              condition = sprintf("input['%s'] == 'blend'", ns("expr.mode")),
              textInput(ns("gene1"), "Gene 1", placeholder = "e.g. Cd68"),
              textInput(ns("gene2"), "Gene 2", placeholder = "e.g. Cx3cr1"),
              selectInput(ns("coexpr.meta"), "Subset by Metadata Column:",
                          choices = c("None"), selected = "None"),
              conditionalPanel(
                condition = sprintf("input['%s'] != 'None'", ns("coexpr.meta")),
                selectInput(ns("coexpr.ident"), "Select Identity:",
                            choices = c("All"), selected = "All")
              )
            ),
            actionButton(ns("run.expr"), "Generate Plot",
                         class = "btn-success", style = "width:100%"),
            hr(),
            h4("Download Figure"),
            numericInput(ns("fig.expr.w"), "Width (inches)", value = 8, min = 2, max = 20),
            numericInput(ns("fig.expr.h"), "Height (inches)", value = 6, min = 2, max = 20),
            downloadButton(ns("download.expr"), "Download Figure", class = "btn-success")
          ),
          mainPanel(width = 8,
            plotOutput(ns("plot.expr"), height = "600px")
          )
        )
      )
    )
  )
}

mod.explore.umap.server <- function(id, shared.data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    active.obj <- reactiveVal(NULL)
    selected.cell.ids <- reactiveVal(NULL)

    # Build the default color palette from the full dataset (consistent across all sub-tabs)
    full.palette <- reactive({
      req(shared.data())
      all.idents <- levels(shared.data())
      pal <- scales::hue_pal()(length(all.idents))
      names(pal) <- all.idents
      pal
    })

    # Initialize when data is uploaded
    observe({
      req(shared.data())
      active.obj(shared.data())
      ident.levels <- levels(shared.data())
      updateSelectInput(session, "select.idents",
                        choices = ident.levels, selected = ident.levels)
      # Populate categorical metadata selectors
      cat.cols <- get.categorical.meta(shared.data())
      updateSelectInput(session, "split.by",
                        choices = c("None", cat.cols), selected = "None")
      updateSelectInput(session, "coexpr.meta",
                        choices = c("None", cat.cols), selected = "None")
    })

    # Reactive for detected slots on active object
    detect.slots <- reactive({
      req(active.obj())
      detect.slots.for(active.obj())
    })

    # ---- SUB-TAB 1: INTERACTIVE SELECTION ----

    # Capture lasso selection
    observeEvent(plotly::event_data("plotly_selected", source = "umap.lasso"), {
      sel <- plotly::event_data("plotly_selected", source = "umap.lasso")
      if (!is.null(sel) && !is.null(sel$customdata)) {
        ids <- unique(as.character(sel$customdata))
        selected.cell.ids(ids)
        showNotification(
          paste(length(ids), "cells selected. Switch to 'subset' mode to compute."),
          type = "message")
      }
    })

    # Run re-UMAP on selected cells
    observeEvent(input$run.subset.umap, {
      req(selected.cell.ids())

      withProgress(message = "Computing subset UMAP...", value = 0, {
        tryCatch({
          sub <- subset(shared.data(), cells = selected.cell.ids())
          slots <- detect.slots.for(sub)
          incProgress(0.3, detail = "Running UMAP")

          # Prefer PCA/integrated reduction; fall back to graph-based UMAP
          red.name <- grep("pca|integrated", Reductions(sub), value = TRUE, ignore.case = TRUE)[1]
          if (!is.na(red.name)) {
            n.dims <- min(input$subset.pcs, ncol(Embeddings(sub, red.name)))
            sub <- RunUMAP(sub, reduction = red.name, dims = 1:n.dims,
                           reduction.name = "umap_subset")
          } else if (!is.na(slots$graph)) {
            sub <- RunUMAP(sub, graph = slots$graph, reduction.name = "umap_subset")
          } else {
            stop("No PCA/integrated reduction or SNN graph found for re-UMAP.")
          }

          incProgress(0.7, detail = "Done")
          active.obj(sub)
          showNotification("Subset UMAP complete!", type = "default")
        }, error = function(e) {
          message("UMAP error: ", e$message)
          showNotification(paste("Error during UMAP:", e$message), type = "error")
        })
      })
    })

    # Reset to original data
    observeEvent(input$reset.umap, {
      active.obj(shared.data())
      selected.cell.ids(NULL)
      updateRadioButtons(session, "display.mode", selected = "full")
    })

    # Render interactive UMAP
    output$umap.interactive <- renderPlotly({
      req(active.obj())
      obj <- active.obj()
      slots <- detect.slots()

      p <- DimPlot(obj, reduction = slots$reduction, cols = full.palette())
      # Map cell barcodes to plotly customdata for lasso selection
      p$layers[[1]]$mapping$customdata <- ggplot2::aes(label = rownames(p$data))$label

      plotly::ggplotly(p + theme_minimal(), source = "umap.lasso") %>%
        plotly::layout(dragmode = "lasso")
    })

    # Cell stats
    output$cell.stat.box <- renderText({
      req(active.obj())
      paste0("Current Cell Count: ", ncol(active.obj()),
             "\nReduction: ", detect.slots()$reduction,
             if (!is.null(selected.cell.ids()))
               paste0("\nSelected: ", length(selected.cell.ids()), " cells")
             else "")
    })

    # ---- SUB-TAB 2: HIGHLIGHT VIEW ----

    # Select All button
    observeEvent(input$btn.select.all, {
      obj <- active.obj()
      req(obj)
      updateSelectInput(session, "select.idents",
                        choices = levels(obj), selected = levels(obj))
    })

    # Clear button
    observeEvent(input$btn.clear.all, {
      obj <- active.obj()
      req(obj)
      updateSelectInput(session, "select.idents",
                        choices = levels(obj), selected = character(0))
    })

    # Store the highlight plot for download
    plot.highlight <- reactive({
      req(active.obj())
      obj <- active.obj()
      slots <- detect.slots()
      idents.selected <- input$select.idents

      req(length(idents.selected) > 0)

      # Build color palette: selected types keep their original color, rest are grey
      all.idents <- levels(obj)
      pal <- full.palette()
      display.pal <- setNames(rep("lightgrey", length(all.idents)), all.idents)
      for (id in idents.selected) {
        if (id %in% names(pal)) display.pal[id] <- pal[id]
      }

      p <- DimPlot(obj, reduction = slots$reduction, cols = display.pal) + theme_minimal()

      # Keep only highlighted types in the legend
      if (length(idents.selected) < length(all.idents)) {
        p <- p + guides(color = guide_legend(override.aes = list(size = 3))) +
          scale_color_manual(values = display.pal,
                             breaks = idents.selected,
                             labels = idents.selected)
      }

      p
    })

    output$umap.static <- renderPlot({ plot.highlight() }, res = 96)

    output$download.highlight <- downloadHandler(
      filename = function() { "umap_highlight.pdf" },
      content = function(file) {
        ggsave(file, plot = plot.highlight(),
               width = input$fig.highlight.w, height = input$fig.highlight.h)
      }
    )

    # ---- SUB-TAB 3: GENE EXPRESSION ----

    # Update co-expression metadata identity choices
    observeEvent(input$coexpr.meta, {
      req(shared.data(), input$coexpr.meta != "None")
      vals <- unique(as.character(shared.data()@meta.data[[input$coexpr.meta]]))
      updateSelectInput(session, "coexpr.ident",
                        choices = c("All", vals), selected = "All")
    })

    plot.expr <- eventReactive(input$run.expr, {
      req(active.obj())
      obj <- active.obj()
      slots <- detect.slots()

      if (input$expr.mode == "single") {
        gene <- trimws(input$gene.name)
        validate(need(nchar(gene) > 0, "Please enter a gene name."))
        validate(need(gene %in% rownames(obj), paste0("Gene '", gene, "' not found.")))

        split.col <- if (input$split.by == "None") NULL else input$split.by

        if (is.null(split.col)) {
          reorder.feature.plot(obj, gene, slots$reduction)
        } else {
          # Compute shared ranges for consistent axes and color scale
          expr.vals <- GetAssayData(obj, layer = "data")[gene, ]
          expr.range <- range(expr.vals)
          umap.coords <- Embeddings(obj, reduction = slots$reduction)
          x.range <- range(umap.coords[, 1])
          y.range <- range(umap.coords[, 2])

          groups <- unique(as.character(obj@meta.data[[split.col]]))
          plot.list <- lapply(groups, function(grp) {
            cells <- colnames(obj)[obj@meta.data[[split.col]] == grp]
            reorder.feature.plot(obj, gene, slots$reduction, cells = cells) +
              scale_color_gradientn(colors = c("lightgrey", "blue"),
                                   limits = expr.range) +
              xlim(x.range) + ylim(y.range) +
              ggtitle(grp)
          })
          patchwork::wrap_plots(plot.list, ncol = input$grid.ncol) +
            patchwork::plot_layout(guides = "collect")
        }

      } else {
        gene1 <- trimws(input$gene1)
        gene2 <- trimws(input$gene2)
        validate(need(nchar(gene1) > 0 && nchar(gene2) > 0, "Please enter both gene names."))
        validate(need(gene1 %in% rownames(obj), paste0("Gene '", gene1, "' not found.")))
        validate(need(gene2 %in% rownames(obj), paste0("Gene '", gene2, "' not found.")))

        # Subset by metadata if selected
        plot.obj <- obj
        if (input$coexpr.meta != "None" && input$coexpr.ident != "All") {
          cells.keep <- colnames(obj)[obj@meta.data[[input$coexpr.meta]] == input$coexpr.ident]
          validate(need(length(cells.keep) > 0, "No cells match the selected identity."))
          plot.obj <- subset(obj, cells = cells.keep)
        }

        reorder.feature.plot(plot.obj, c(gene1, gene2), slots$reduction, blend = TRUE)
      }
    })

    output$plot.expr <- renderPlot({ plot.expr() }, res = 96)

    output$download.expr <- downloadHandler(
      filename = function() { "gene_expression.pdf" },
      content = function(file) {
        ggsave(file, plot = plot.expr(),
               width = input$fig.expr.w, height = input$fig.expr.h)
      }
    )
  })
}
