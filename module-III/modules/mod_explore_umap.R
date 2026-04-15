## ABOUTME: UMAP cell-type visualization with interactive selection, highlighting, and gene expression.
## ABOUTME: Three sub-tabs: interactive lasso/re-UMAP, high-res figure generation, and gene (co-)expression plots.

# Helper: detect best UMAP reduction and SNN graph for a given Seurat object
detect.slots.for <- function(obj) {
  red.name <- grep("umap", Reductions(obj), value = TRUE, ignore.case = TRUE)[1]
  graph.name <- grep("snn", Graphs(obj), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(red.name)) red.name <- "umap"
  list(reduction = red.name, graph = graph.name)
}

mod.explore.umap.ui <- function(id) {
  ns <- NS(id)
  tabPanel("UMAP Cell-Type",
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
            checkboxInput(ns("show.all"), "Select All", value = TRUE),
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
                          choices = c("None"), selected = "None")
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
    observeEvent(event_data("plotly_selected", source = "umap.lasso"), {
      sel <- event_data("plotly_selected", source = "umap.lasso")
      if (!is.null(sel) && !is.null(sel$customdata)) {
        selected.cell.ids(sel$customdata)
        showNotification(
          paste(length(sel$customdata), "cells selected. Switch to 'subset' mode to compute."),
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
          incProgress(0.3, detail = "Running UMAP on graph")

          if (!is.na(slots$graph)) {
            sub <- RunUMAP(sub, graph = slots$graph, reduction.name = "umap_subset")
          } else {
            red.name <- grep("pca|integrated", Reductions(sub), value = TRUE, ignore.case = TRUE)[1]
            if (is.na(red.name)) red.name <- "pca"
            sub <- RunUMAP(sub, reduction = red.name, dims = 1:input$subset.pcs,
                           reduction.name = "umap_subset")
          }

          incProgress(0.7, detail = "Done")
          active.obj(sub)
          showNotification("Subset UMAP complete!", type = "default")
        }, error = function(e) {
          showNotification(paste("Error during UMAP:", e$message), type = "error")
        })
      })
    })

    # Reset to original data
    observeEvent(input$reset.umap, {
      active.obj(shared.data())
      selected.cell.ids(NULL)
    })

    # Render interactive UMAP
    output$umap.interactive <- renderPlotly({
      req(active.obj())
      obj <- active.obj()
      slots <- detect.slots()

      p <- DimPlot(obj, reduction = slots$reduction)
      p$data$cell_id <- rownames(p$data)

      ggplotly(p, source = "umap.lasso") %>%
        layout(dragmode = "lasso") %>%
        event_register("plotly_selected")
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

    # Toggle select all
    observeEvent(input$show.all, {
      obj <- active.obj()
      req(obj)
      if (input$show.all) {
        updateSelectInput(session, "select.idents",
                          choices = levels(obj), selected = levels(obj))
      }
    })

    # Store the highlight plot for download
    plot.highlight <- reactive({
      req(active.obj())
      obj <- active.obj()
      slots <- detect.slots()

      idents.selected <- if (input$show.all) levels(obj) else input$select.idents
      req(length(idents.selected) > 0)

      if (length(idents.selected) < length(levels(obj))) {
        # Preserve original colors for selected types, grey for rest
        all.idents <- levels(obj)
        color.palette <- scales::hue_pal()(length(all.idents))
        names(color.palette) <- all.idents

        cells.to.highlight <- WhichCells(obj, idents = idents.selected)
        highlight.colors <- color.palette[idents.selected]

        DimPlot(obj,
                reduction = slots$reduction,
                cells.highlight = cells.to.highlight,
                cols.highlight = unname(highlight.colors),
                cols = "lightgrey",
                sizes.highlight = 1) +
          theme_minimal()
      } else {
        DimPlot(obj, reduction = slots$reduction) + theme_minimal()
      }
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

        split.by <- if (input$split.by == "None") NULL else input$split.by
        FeaturePlot(obj, features = gene, reduction = slots$reduction,
                    split.by = split.by) + theme_minimal()

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

        FeaturePlot(plot.obj, features = c(gene1, gene2), blend = TRUE,
                    reduction = slots$reduction) + theme_minimal()
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
