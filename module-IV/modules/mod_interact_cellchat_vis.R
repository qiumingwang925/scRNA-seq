## ABOUTME: CellChat visualization module with four subtabs: global network, pathway/LR
## ABOUTME: zoom-in, signaling role, and communication patterns. All plots render per-group.

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---------- Helpers ----------

check.python.umap <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return(list(ok = FALSE, msg = "R package 'reticulate' is not installed."))
  }
  ok <- tryCatch(reticulate::py_module_available("umap"), error = function(e) FALSE)
  if (!ok) {
    return(list(ok = FALSE,
                msg = paste("Python module 'umap-learn' not available.",
                            "Install via: reticulate::py_install('umap-learn').")))
  }
  list(ok = TRUE, msg = "")
}

# Arrange a list of ComplexHeatmap objects in a patchwork grid. Tolerates heatmaps
# with differing row/column dimensions (which horizontal `+` concatenation cannot).
draw.heatmap.grid <- function(ht.list, ncol, titles = NULL) {
  plots <- lapply(seq_along(ht.list), function(i) {
    tryCatch({
      patchwork::wrap_elements(grid::grid.grabExpr(ComplexHeatmap::draw(ht.list[[i]])))
    }, error = function(e) {
      label <- if (!is.null(titles)) titles[i] else paste("Panel", i)
      ggplot() + theme_void() +
        ggtitle(paste0(label, "\n[", conditionMessage(e), "]"))
    })
  })
  print(patchwork::wrap_plots(plots, ncol = ncol))
}

# Arrange a base-graphics render function (called per group) in a grid layout
plot.base.grid <- function(groups, ncol, fn, mar = c(2, 2, 3, 2)) {
  n <- length(groups)
  nrow <- ceiling(n / ncol)
  old <- par(mfrow = c(nrow, ncol), xpd = TRUE, mar = mar, no.readonly = TRUE)
  on.exit(par(old))
  for (g in groups) {
    tryCatch(fn(g), error = function(e) {
      plot.new()
      title(paste0(g, "\n[", conditionMessage(e), "]"))
    })
  }
}

# Wrap a thunk so renderPlot and download handler share identical logic
make.render.and.download <- function(output, session, plot.id, dl.id, thunk,
                                     width.in, height.in, filename.stem) {
  output[[plot.id]] <- renderPlot({ thunk() }, res = 96)
  output[[dl.id]] <- downloadHandler(
    filename = function() paste0(filename.stem(), "_", Sys.Date(), ".png"),
    content = function(file) {
      png(file, width = width.in() * 96, height = height.in() * 96, res = 96)
      tryCatch(thunk(), error = function(e) {
        plot.new()
        title(paste("Render failed:", conditionMessage(e)))
      })
      dev.off()
    }
  )
}

# Cell-type selector matching Module III: selectInput(multiple) + Select All / Clear
cell.type.selector.ui <- function(ns, id, label = "Cell Type(s):") {
  tagList(
    selectInput(ns(id), label, choices = NULL, multiple = TRUE),
    fluidRow(
      column(6, actionButton(ns(paste0(id, ".all")), "Select All",
                             class = "btn-info", style = "width:100%")),
      column(6, actionButton(ns(paste0(id, ".clear")), "Clear",
                             class = "btn-default", style = "width:100%"))
    )
  )
}

wire.cell.type.selector <- function(input, session, id, choices.reactive, selected.all = TRUE) {
  observe({
    choices <- choices.reactive()
    req(choices)
    updateSelectInput(session, id, choices = choices,
                      selected = if (selected.all) choices else character(0))
  })
  observeEvent(input[[paste0(id, ".all")]], {
    choices <- choices.reactive()
    req(choices)
    updateSelectInput(session, id, choices = choices, selected = choices)
  })
  observeEvent(input[[paste0(id, ".clear")]], {
    choices <- choices.reactive()
    req(choices)
    updateSelectInput(session, id, choices = choices, selected = character(0))
  })
}

# Common sidebar controls (width / height / download)
common.controls.ui <- function(ns, prefix, dl.label = "Download Figure") {
  tagList(
    hr(),
    h5("Figure controls"),
    numericInput(ns(paste0(prefix, ".width")), "Width (inches):",
                 value = 10, min = 3, max = 30, step = 1),
    numericInput(ns(paste0(prefix, ".height")), "Height (inches):",
                 value = 6, min = 3, max = 30, step = 1),
    downloadButton(ns(paste0(prefix, ".download")), dl.label, style = "width:100%")
  )
}

# ---------- UI ----------

mod.interact.cellchat.vis.ui <- function(id) {
  ns <- NS(id)
  tabPanel("Visualize CellChat",
    fluidRow(
      column(12,
        wellPanel(
          fileInput(ns("vis.file"),
                    "Upload processed CellChat result (.rds) — skip if computed above:",
                    accept = ".rds"),
          verbatimTextOutput(ns("input.summary")),
          hr(),
          radioButtons(ns("group.selected"), "Group to display:",
                       choices = character(0), inline = TRUE)
        )
      )
    ),

    tabsetPanel(id = ns("vis.tabs"),

      # ---------- Subtab 1: Global ----------
      tabPanel("Global Network",
        sidebarLayout(
          sidebarPanel(width = 3,
            h4("Global network"),
            radioButtons(ns("global.plot"), "Plot type:",
                         choices = c("Circle" = "circle", "Heatmap" = "heatmap"),
                         selected = "circle"),
            radioButtons(ns("global.measure"), "Measure:",
                         choices = c("Interaction count" = "count",
                                     "Interaction weight/strength" = "weight"),
                         selected = "weight"),
            cell.type.selector.ui(ns, "global.sources", "Sources:"),
            cell.type.selector.ui(ns, "global.targets", "Targets:"),
            common.controls.ui(ns, "global")
          ),
          mainPanel(width = 9,
            uiOutput(ns("global.plot.ui"))
          )
        )
      ),

      # ---------- Subtab 2: Zoom-in ----------
      tabPanel("Zoom-in (Pathway & LR)",
        sidebarLayout(
          sidebarPanel(width = 3,
            h4("Pathway / L-R zoom"),
            selectInput(ns("zoom.plot"), "Plot type:",
                        choices = c(
                          "Pathway — Hierarchy" = "pw.hier",
                          "Pathway — Circle" = "pw.circle",
                          "Pathway — Heatmap" = "pw.heat",
                          "Pathway — Chord" = "pw.chord",
                          "L-R pair — Hierarchy" = "lr.hier",
                          "L-R pair — Circle" = "lr.circle",
                          "L-R pair — Chord" = "lr.chord",
                          "Pathway/LR — Bubble" = "bubble",
                          "Pathway — Violin" = "violin"
                        )),
            selectInput(ns("zoom.pathway"), "Pathway:", choices = NULL),
            conditionalPanel(
              condition = sprintf("['lr.hier','lr.circle','lr.chord'].indexOf(input['%s']) >= 0",
                                  ns("zoom.plot")),
              selectInput(ns("zoom.lr"), "L-R pair:", choices = NULL)
            ),
            cell.type.selector.ui(ns, "zoom.sources", "Sources:"),
            cell.type.selector.ui(ns, "zoom.targets",
                                  "Targets (also used as hierarchy receivers):"),
            common.controls.ui(ns, "zoom")
          ),
          mainPanel(width = 9,
            uiOutput(ns("zoom.plot.ui"))
          )
        )
      ),

      # ---------- Subtab 3: Signaling role ----------
      tabPanel("Signaling-Focused",
        sidebarLayout(
          sidebarPanel(width = 3,
            h4("Signaling role"),
            radioButtons(ns("sig.plot"), "Plot type:",
                         choices = c("Score (per-pathway centrality)" = "score",
                                     "Scatter (in vs out strength)" = "scatter",
                                     "Heatmap (pathway × cell type)" = "heatmap"),
                         selected = "scatter"),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'score'", ns("sig.plot")),
              selectInput(ns("sig.pathway"), "Pathway (score only):", choices = NULL)
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'heatmap'", ns("sig.plot")),
              radioButtons(ns("sig.pattern"), "Direction (heatmap only):",
                           choices = c("Outgoing" = "outgoing",
                                       "Incoming" = "incoming",
                                       "All" = "all"),
                           selected = "all")
            ),
            common.controls.ui(ns, "sig")
          ),
          mainPanel(width = 9,
            uiOutput(ns("sig.plot.ui"))
          )
        )
      ),

      # ---------- Subtab 4: Communication patterns ----------
      tabPanel("Communication Patterns",
        sidebarLayout(
          sidebarPanel(width = 3,
            h4("Patterns"),
            radioButtons(ns("pat.plot"), "Plot type:",
                         choices = c("Pattern heatmap (select k)" = "heat",
                                     "River" = "river",
                                     "Dot" = "dot",
                                     "Manifold & Classification" = "manifold"),
                         selected = "heat"),
            conditionalPanel(
              condition = sprintf("input['%s'] != 'manifold'", ns("pat.plot")),
              radioButtons(ns("pat.direction"), "Direction:",
                           choices = c("Outgoing" = "outgoing", "Incoming" = "incoming"),
                           selected = "outgoing"),
              numericInput(ns("pat.k"), "Number of patterns (k):",
                           value = 3, min = 2, max = 10, step = 1)
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'manifold'", ns("pat.plot")),
              radioButtons(ns("pat.sim.type"), "Similarity type:",
                           choices = c("Functional" = "functional",
                                       "Structural" = "structural"),
                           selected = "functional"),
              helpText("Manifold requires Python 'umap-learn' via reticulate.")
            ),
            common.controls.ui(ns, "pat")
          ),
          mainPanel(width = 9,
            uiOutput(ns("pat.plot.ui"))
          )
        )
      )
    )
  )
}

# ---------- Server ----------

mod.interact.cellchat.vis.server <- function(id, cellchat.input) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Resolve cellchat result: upstream compute > uploaded file
    cellchat.data <- reactive({
      upstream <- NULL
      try(upstream <- cellchat.input(), silent = TRUE)
      if (!is.null(upstream) && is.list(upstream) &&
          all(c("cellchat.list", "cellchat.merged", "group.levels") %in% names(upstream))) {
        return(upstream)
      }

      req(input$vis.file)
      tryCatch({
        r <- readRDS(input$vis.file$datapath)
        if (!is.list(r) || !all(c("cellchat.list", "cellchat.merged", "group.levels") %in% names(r))) {
          showNotification("Uploaded file is not a valid CellChat result object.",
                           type = "error", duration = NULL)
          return(NULL)
        }
        r
      }, error = function(e) {
        showNotification(paste("Error loading file:", e$message),
                         type = "error", duration = NULL)
        NULL
      })
    })

    output$input.summary <- renderPrint({
      res <- cellchat.data()
      req(res)
      cat("CellChat result loaded.\n")
      cat("Groups: ", paste(res$group.levels, collapse = ", "), "\n")
      cat("Species:", res$species %||% "unknown", "\n")
      cat("Assay: ", res$assay %||% "unknown", "\n")
    })

    # Shared choice reactives
    pathway.choices <- reactive({
      res <- cellchat.data()
      req(res)
      sort(Reduce(union, lapply(res$cellchat.list, function(cc) cc@netP$pathways)))
    })

    ident.choices <- reactive({
      res <- cellchat.data()
      req(res)
      Reduce(union, lapply(res$cellchat.list, function(cc) levels(cc@idents)))
    })

    lr.choices <- reactive({
      res <- cellchat.data()
      req(res, input$zoom.pathway)
      Reduce(union, lapply(res$cellchat.list, function(cc) {
        tryCatch({
          df <- extractEnrichedLR(cc, signaling = input$zoom.pathway, geneLR.return = FALSE)
          if (is.null(df) || !"interaction_name" %in% colnames(df)) return(character())
          df$interaction_name
        }, error = function(e) character())
      }))
    })

    observe({
      pw <- pathway.choices()
      updateSelectInput(session, "zoom.pathway", choices = pw,
                        selected = if (length(pw)) pw[1] else NULL)
      updateSelectInput(session, "sig.pathway", choices = pw,
                        selected = if (length(pw)) pw[1] else NULL)
    })

    observe({
      lr <- lr.choices()
      updateSelectInput(session, "zoom.lr", choices = lr,
                        selected = if (length(lr)) lr[1] else NULL)
    })

    wire.cell.type.selector(input, session, "global.sources", ident.choices, selected.all = TRUE)
    wire.cell.type.selector(input, session, "global.targets", ident.choices, selected.all = TRUE)
    wire.cell.type.selector(input, session, "zoom.sources", ident.choices, selected.all = TRUE)
    wire.cell.type.selector(input, session, "zoom.targets", ident.choices, selected.all = TRUE)

    observe({
      res <- cellchat.data()
      req(res)
      updateRadioButtons(session, "group.selected",
                         choices = res$group.levels,
                         selected = res$group.levels[1], inline = TRUE)
    })

    # Auto-size plot output based on user width/height
    plot.output <- function(prefix) {
      renderUI({
        w <- input[[paste0(prefix, ".width")]]
        h <- input[[paste0(prefix, ".height")]]
        req(w, h)
        shinycssloaders::withSpinner(
          plotOutput(ns(paste0(prefix, ".plot")),
                     width = paste0(w * 96, "px"),
                     height = paste0(h * 96, "px")),
          type = 6
        )
      })
    }

    output$global.plot.ui <- plot.output("global")
    output$zoom.plot.ui <- plot.output("zoom")
    output$sig.plot.ui <- plot.output("sig")
    output$pat.plot.ui <- plot.output("pat")

    # =====================================================
    # Subtab 1: Global network
    # =====================================================

    # Resolve a selector value into: filter active? which values to use? is the
    # intersection with this group empty?
    resolve.sel <- function(cc, sel, all.idents) {
      if (is.null(sel) || length(sel) == 0 || setequal(sel, all.idents)) {
        return(list(active = FALSE, use = NULL, empty = FALSE))
      }
      kept <- intersect(sel, levels(cc@idents))
      list(active = TRUE, use = kept, empty = length(kept) == 0)
    }

    placeholder.base <- function(g, msg = "selection absent in this group") {
      plot.new()
      title(paste0(g, "\n(", msg, ")"))
    }
    placeholder.gg <- function(g, msg = "selection absent in this group") {
      ggplot() + theme_void() +
        ggtitle(paste0(g, "\n(", msg, ")"))
    }
    placeholder.ht <- function(g, msg = "selection absent") {
      ComplexHeatmap::Heatmap(
        matrix(0, 1, 1, dimnames = list("-", "-")),
        column_title = paste0(g, " (", msg, ")"),
        show_heatmap_legend = FALSE
      )
    }

    global.thunk <- function() {
      res <- cellchat.data(); req(res, input$global.plot, input$global.measure,
                                 input$group.selected)
      g <- input$group.selected
      cc <- res$cellchat.list[[g]]; req(cc)
      all.idents <- ident.choices()
      srcs <- input$global.sources; tgts <- input$global.targets
      s <- resolve.sel(cc, srcs, all.idents); t <- resolve.sel(cc, tgts, all.idents)
      validate(
        need(!s$empty, paste0("Selected source(s) not present in '", g, "'.")),
        need(!t$empty, paste0("Selected target(s) not present in '", g, "'."))
      )

      mat <- if (input$global.measure == "count") cc@net$count else cc@net$weight
      s.keep <- if (s$active) s$use else rownames(mat)
      t.keep <- if (t$active) t$use else colnames(mat)

      if (input$global.plot == "circle") {
        weight.max <- tryCatch(
          getMaxWeight(res$cellchat.list, attribute = c("idents", input$global.measure)),
          error = function(e) NULL
        )
        vertex.weight.max <- max(unlist(lapply(res$cellchat.list, function(cc) {
          as.numeric(table(cc@idents))
        })))
        nodes <- union(s.keep, t.keep)
        mat.sub <- mat[nodes, nodes, drop = FALSE]
        mat.sub[!(rownames(mat.sub) %in% s.keep), ] <- 0
        mat.sub[, !(colnames(mat.sub) %in% t.keep)] <- 0
        cell.counts <- as.numeric(table(cc@idents))
        names(cell.counts) <- levels(cc@idents)
        netVisual_circle(mat.sub,
                         vertex.weight = cell.counts[nodes],
                         vertex.weight.max = vertex.weight.max,
                         vertex.size.max = 15,
                         weight.scale = TRUE, label.edge = FALSE,
                         edge.weight.max = if (!is.null(weight.max)) weight.max[2] else NULL,
                         title.name = paste0(g, " — ", input$global.measure))
      } else {
        sub <- mat[s.keep, t.keep, drop = FALSE]
        vmax <- max(sub, na.rm = TRUE)
        col.fn <- circlize::colorRamp2(
          c(0, if (vmax > 0) vmax else 1),
          c("#FFF5F0", "#A50F15")
        )
        ht <- ComplexHeatmap::Heatmap(
          sub, name = input$global.measure, col = col.fn,
          column_title = g, row_title = "Sources",
          column_title_side = "top", row_names_side = "left",
          cluster_rows = FALSE, cluster_columns = FALSE
        )
        ComplexHeatmap::draw(ht)
      }
    }

    make.render.and.download(output, session, "global.plot", "global.download",
      global.thunk,
      width.in = reactive(input$global.width),
      height.in = reactive(input$global.height),
      filename.stem = reactive(paste0("global_", input$global.plot, "_", input$global.measure)))

    # =====================================================
    # Subtab 2: Zoom-in (pathway & LR)
    # =====================================================

    # Resolve source/target args (empty = NULL = all)
    nz <- function(x) if (!is.null(x) && length(x) > 0) x else NULL

    zoom.thunk <- function() {
      res <- cellchat.data(); req(res, input$zoom.plot, input$zoom.pathway,
                                 input$group.selected)
      g <- input$group.selected
      cc <- res$cellchat.list[[g]]; req(cc)
      pw <- input$zoom.pathway
      lr <- input$zoom.lr
      srcs <- nz(input$zoom.sources); tgts <- nz(input$zoom.targets)
      all.idents <- ident.choices()
      s <- resolve.sel(cc, srcs, all.idents); t <- resolve.sel(cc, tgts, all.idents)
      validate(
        need(!s$empty, paste0("Selected source(s) not present in '", g, "'.")),
        need(!t$empty, paste0("Selected target(s) not present in '", g, "'.")),
        need(pw %in% cc@netP$pathways,
             paste0("Pathway '", pw, "' not present in '", g, "'."))
      )

      pt <- input$zoom.plot
      is.lr <- startsWith(pt, "lr.")
      if (is.lr) {
        req(lr)
        enriched <- tryCatch(
          extractEnrichedLR(cc, signaling = pw, geneLR.return = FALSE)$interaction_name,
          error = function(e) character()
        )
        validate(need(lr %in% enriched,
                      paste0("L-R pair '", lr, "' not enriched in '", g, "'.")))
      }

      hier.receiver <- function() {
        lv <- levels(cc@idents)
        if (!is.null(tgts) && length(tgts) > 0) {
          idx <- which(lv %in% tgts)
          if (length(idx) > 0 && length(idx) < length(lv)) return(idx)
        }
        seq_len(ceiling(length(lv) / 2))
      }

      if (pt == "bubble") {
        return(print(netVisual_bubble(cc, signaling = pw,
                                      sources.use = if (s$active) s$use else NULL,
                                      targets.use = if (t$active) t$use else NULL,
                                      remove.isolate = FALSE) +
                       ggtitle(g)))
      }

      if (pt == "violin") {
        cc.sub <- cc
        keep <- union(if (s$active) s$use else NULL,
                      if (t$active) t$use else NULL)
        if (length(keep) > 0 && length(keep) < length(levels(cc@idents))) {
          cells <- names(cc@idents)[cc@idents %in% keep]
          cc.sub <- subsetCellChat(cc, cells.use = cells)
        }
        return(print(plotGeneExpression(cc.sub, signaling = pw, enriched.only = TRUE) +
                       patchwork::plot_annotation(title = g)))
      }

      if (pt == "pw.heat") {
        args <- list(object = cc, signaling = pw,
                     color.heatmap = "Reds", title.name = g)
        if (s$active) args$sources.use <- s$use
        if (t$active) args$targets.use <- t$use
        return(ComplexHeatmap::draw(do.call(netVisual_heatmap, args)))
      }

      layout.map <- c(pw.hier = "hierarchy", pw.circle = "circle", pw.chord = "chord",
                      lr.hier = "hierarchy", lr.circle = "circle", lr.chord = "chord")
      layout <- layout.map[[pt]]

      common <- list(object = cc, signaling = pw, layout = layout)
      if (layout == "hierarchy") {
        common$vertex.receiver <- hier.receiver()
      } else {
        if (s$active) common$sources.use <- s$use
        if (t$active) common$targets.use <- t$use
      }
      if (is.lr) {
        common$pairLR.use <- lr
        do.call(netVisual_individual, common)
      } else {
        do.call(netVisual_aggregate, common)
      }
    }

    make.render.and.download(output, session, "zoom.plot", "zoom.download",
      zoom.thunk,
      width.in = reactive(input$zoom.width),
      height.in = reactive(input$zoom.height),
      filename.stem = reactive(paste0("zoom_", input$zoom.plot, "_", input$zoom.pathway)))

    # =====================================================
    # Subtab 3: Signaling role
    # =====================================================

    compute.centrality <- function(cc) {
      # Ensure centrality scores are computed (idempotent)
      if (is.null(cc@netP$centr) || length(cc@netP$centr) == 0) {
        tryCatch(cc <- netAnalysis_computeCentrality(cc, slot.name = "netP"),
                 error = function(e) message("Centrality failed: ", e$message))
      }
      cc
    }

    sig.thunk <- function() {
      res <- cellchat.data(); req(res, input$sig.plot, input$group.selected)
      g <- input$group.selected
      cc <- compute.centrality(res$cellchat.list[[g]]); req(cc)

      if (input$sig.plot == "score") {
        req(input$sig.pathway)
        validate(need(input$sig.pathway %in% cc@netP$pathways,
                      paste0("Pathway '", input$sig.pathway, "' not in '", g, "'.")))
        netAnalysis_signalingRole_network(cc, signaling = input$sig.pathway,
                                          width = 8, height = 2.5, font.size = 10)
        return(invisible())
      }

      if (input$sig.plot == "scatter") {
        return(print(netAnalysis_signalingRole_scatter(cc) + ggtitle(g)))
      }

      if (input$sig.plot == "heatmap") {
        ComplexHeatmap::draw(netAnalysis_signalingRole_heatmap(
          cc, pattern = input$sig.pattern,
          width = 5, height = 8, title = g))
      }
    }

    make.render.and.download(output, session, "sig.plot", "sig.download",
      sig.thunk,
      width.in = reactive(input$sig.width),
      height.in = reactive(input$sig.height),
      filename.stem = reactive(paste0("signaling_", input$sig.plot)))

    # =====================================================
    # Subtab 4: Communication patterns
    # =====================================================

    pat.thunk <- function() {
      res <- cellchat.data(); req(res, input$pat.plot)

      if (input$pat.plot == "manifold") {
        check <- check.python.umap()
        if (!check$ok) {
          plot.new(); title(paste("Manifold unavailable:\n", check$msg))
          return(invisible())
        }
        req(input$pat.sim.type)
        merged <- res$cellchat.merged
        tryCatch({
          merged <- computeNetSimilarityPairwise(merged, type = input$pat.sim.type)
          merged <- netEmbedding(merged, type = input$pat.sim.type)
          merged <- netClustering(merged, type = input$pat.sim.type, do.parallel = FALSE)
          print(netVisual_embeddingPairwise(merged, type = input$pat.sim.type, label.size = 3.5))
        }, error = function(e) {
          plot.new(); title(paste("Manifold failed:", conditionMessage(e)))
        })
        return(invisible())
      }

      req(input$pat.direction, input$pat.k, input$group.selected)
      g <- input$group.selected
      direction <- input$pat.direction
      k <- input$pat.k
      cc <- compute.centrality(res$cellchat.list[[g]]); req(cc)
      cc <- identifyCommunicationPatterns(cc, pattern = direction, k = k,
                                          width = 5, height = 9)

      if (input$pat.plot == "heat") {
        ComplexHeatmap::draw(ComplexHeatmap::Heatmap(
          cc@netP$pattern[[direction]]$pattern$cell,
          name = paste0(g, " cell"), column_title = g))
        return(invisible())
      }
      if (input$pat.plot == "river") {
        print(netAnalysis_river(cc, pattern = direction) + ggtitle(g))
        return(invisible())
      }
      if (input$pat.plot == "dot") {
        print(netAnalysis_dot(cc, pattern = direction) + ggtitle(g))
      }
    }

    make.render.and.download(output, session, "pat.plot", "pat.download",
      pat.thunk,
      width.in = reactive(input$pat.width),
      height.in = reactive(input$pat.height),
      filename.stem = reactive(paste0("patterns_", input$pat.plot)))

  })
}
