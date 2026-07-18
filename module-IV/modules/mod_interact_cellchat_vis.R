## ABOUTME: CellChat visualization module with four subtabs: global network, pathway/LR
## ABOUTME: zoom-in, signaling role, and communication patterns. All plots render per-group.

# Shared helpers (%||%, plot.grid, make.render.and.download, cell.type.selector.*,
# common.controls.ui, resolve.sel, placeholder.*) come from mod_interact_vis_helpers.R
# which app.R sources before this file.

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
          verbatimTextOutput(ns("input.summary"))
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
            common.controls.ui(ns, "global", default.cols = 2)
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
            common.controls.ui(ns, "zoom", default.cols = 2)
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
            common.controls.ui(ns, "sig", default.cols = 2)
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
              selectInput(ns("pat.group"), "Group:", choices = NULL),
              numericInput(ns("pat.k"), "Number of patterns (k):",
                           value = 3, min = 2, max = 10, step = 1),
              helpText("Shows outgoing and incoming patterns side by side.")
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'manifold'", ns("pat.plot")),
              radioButtons(ns("pat.sim.type"), "Similarity type:",
                           choices = c("Functional" = "functional",
                                       "Structural" = "structural"),
                           selected = "functional")
            ),
            common.controls.ui(ns, "pat", default.cols = 2)
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
      updateSelectInput(session, "pat.group",
                        choices = res$group.levels,
                        selected = res$group.levels[1])
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

    global.thunk <- function() {
      res <- cellchat.data(); req(res, input$global.plot, input$global.measure)
      grps <- res$group.levels
      ncol <- input$global.cols %||% 2
      all.idents <- ident.choices()
      srcs <- input$global.sources; tgts <- input$global.targets
      measure <- input$global.measure

      weight.max <- tryCatch(
        getMaxWeight(res$cellchat.list, attribute = c("idents", measure)),
        error = function(e) NULL
      )
      vertex.weight.max <- max(unlist(lapply(res$cellchat.list, function(cc) {
        as.numeric(table(cc@idents))
      })))

      items <- lapply(grps, function(g) {
        force(g)
        cc <- res$cellchat.list[[g]]
        s <- resolve.sel(levels(cc@idents), srcs, all.idents)
        t <- resolve.sel(levels(cc@idents), tgts, all.idents)
        mat <- if (measure == "count") cc@net$count else cc@net$weight
        s.keep <- if (s$active) s$use else rownames(mat)
        t.keep <- if (t$active) t$use else colnames(mat)

        if (input$global.plot == "circle") {
          # Base graphics -> function (needs as.ggplot replay)
          function() {
            if (s$empty || t$empty) {
              plot.new(); title(paste0(g, "\n(selection absent)")); return(invisible())
            }
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
                             title.name = paste0(g, " — ", measure))
          }
        } else {
          # Heatmap -> return Heatmap object directly (fast)
          if (s$empty || t$empty) return(placeholder.ht(g))
          sub <- mat[s.keep, t.keep, drop = FALSE]
          vmax <- max(sub, na.rm = TRUE)
          col.fn <- circlize::colorRamp2(
            c(0, if (vmax > 0) vmax else 1),
            c("#FFF5F0", "#A50F15")
          )
          ComplexHeatmap::Heatmap(
            sub, name = measure, col = col.fn,
            column_title = g, row_title = "Sources",
            column_title_side = "top", row_names_side = "left",
            cluster_rows = FALSE, cluster_columns = FALSE
          )
        }
      })
      plot.grid(items, ncol, titles = grps)
    }

    make.render.and.download(output, session, "global.plot", "global.download",
      global.thunk,
      width.in = reactive(input$global.width),
      height.in = reactive(input$global.height),
      filename.stem = reactive(paste0("global_", input$global.plot, "_", input$global.measure)),
      trigger = reactive(input$global.plot.btn))

    # =====================================================
    # Subtab 2: Zoom-in (pathway & LR)
    # =====================================================

    zoom.thunk <- function() {
      res <- cellchat.data(); req(res, input$zoom.plot, input$zoom.pathway)
      grps <- res$group.levels
      ncol <- input$zoom.cols %||% 2
      pw <- input$zoom.pathway
      lr <- input$zoom.lr
      srcs <- nz(input$zoom.sources); tgts <- nz(input$zoom.targets)
      all.idents <- ident.choices()
      pt <- input$zoom.plot
      is.lr <- startsWith(pt, "lr.")
      layout.map <- c(pw.hier = "hierarchy", pw.circle = "circle", pw.chord = "chord",
                      lr.hier = "hierarchy", lr.circle = "circle", lr.chord = "chord")

      # Per-group status: the absence message (if any) for this group
      group.status <- function(cc) {
        s <- resolve.sel(levels(cc@idents), srcs, all.idents)
        t <- resolve.sel(levels(cc@idents), tgts, all.idents)
        msg <- NULL
        if (s$empty || t$empty) msg <- "selection absent"
        else if (!pw %in% cc@netP$pathways) msg <- "pathway absent"
        else if (is.lr) {
          enriched <- tryCatch(
            extractEnrichedLR(cc, signaling = pw, geneLR.return = FALSE)$interaction_name,
            error = function(e) character()
          )
          if (is.null(lr) || !(lr %in% enriched)) msg <- "L-R pair absent"
        }
        list(s = s, t = t, msg = msg)
      }

      items <- lapply(grps, function(g) {
        force(g)
        cc <- res$cellchat.list[[g]]
        st <- group.status(cc)
        s <- st$s; t <- st$t; missing.msg <- st$msg

        if (pt == "bubble") {
          if (!is.null(missing.msg))
            return(placeholder.gg(g, missing.msg))
          return(netVisual_bubble(cc, signaling = pw,
                                  sources.use = if (s$active) s$use else NULL,
                                  targets.use = if (t$active) t$use else NULL,
                                  remove.isolate = FALSE) + ggtitle(g))
        }

        if (pt == "violin") {
          if (!is.null(missing.msg))
            return(placeholder.gg(g, missing.msg))
          cc.sub <- cc
          keep <- union(if (s$active) s$use else NULL,
                        if (t$active) t$use else NULL)
          if (length(keep) > 0 && length(keep) < length(levels(cc@idents))) {
            cells <- names(cc@idents)[cc@idents %in% keep]
            cc.sub <- subsetCellChat(cc, cells.use = cells)
          }
          return(plotGeneExpression(cc.sub, signaling = pw, enriched.only = TRUE) +
                   patchwork::plot_annotation(title = g))
        }

        if (pt == "pw.heat") {
          if (!is.null(missing.msg))
            return(placeholder.ht(g, missing.msg))
          args <- list(object = cc, signaling = pw,
                       color.heatmap = "Reds", title.name = g)
          if (s$active) args$sources.use <- s$use
          if (t$active) args$targets.use <- t$use
          return(do.call(netVisual_heatmap, args))
        }

        # Base-graphics: circle / chord / hierarchy for pw or lr
        function() {
          if (!is.null(missing.msg)) {
            plot.new(); title(paste0(g, "\n(", missing.msg, ")")); return(invisible())
          }
          layout <- layout.map[[pt]]
          common <- list(object = cc, signaling = pw, layout = layout)
          if (layout == "hierarchy") {
            lv <- levels(cc@idents)
            idx <- if (!is.null(tgts) && length(tgts) > 0) which(lv %in% tgts) else integer(0)
            common$vertex.receiver <- if (length(idx) > 0 && length(idx) < length(lv))
              idx else seq_len(ceiling(length(lv) / 2))
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
      })
      plot.grid(items, ncol, titles = grps)
    }

    make.render.and.download(output, session, "zoom.plot", "zoom.download",
      zoom.thunk,
      width.in = reactive(input$zoom.width),
      height.in = reactive(input$zoom.height),
      filename.stem = reactive(paste0("zoom_", input$zoom.plot, "_", input$zoom.pathway)),
      trigger = reactive(input$zoom.plot.btn))

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
      res <- cellchat.data(); req(res, input$sig.plot)
      grps <- res$group.levels
      ncol <- input$sig.cols %||% 2

      items <- lapply(grps, function(g) {
        force(g)
        cc <- compute.centrality(res$cellchat.list[[g]])
        if (input$sig.plot == "scatter") {
          return(netAnalysis_signalingRole_scatter(cc) + ggtitle(g))
        }
        if (input$sig.plot == "heatmap") {
          return(netAnalysis_signalingRole_heatmap(
            cc, pattern = input$sig.pattern,
            width = 5, height = 8, title = g))
        }
        # score: base graphics via netAnalysis_signalingRole_network
        function() {
          req(input$sig.pathway)
          if (!input$sig.pathway %in% cc@netP$pathways) {
            plot.new(); title(paste0(g, "\n(pathway absent)")); return(invisible())
          }
          netAnalysis_signalingRole_network(cc, signaling = input$sig.pathway,
                                            width = 8, height = 2.5, font.size = 10)
        }
      })
      plot.grid(items, ncol, titles = grps)
    }

    make.render.and.download(output, session, "sig.plot", "sig.download",
      sig.thunk,
      width.in = reactive(input$sig.width),
      height.in = reactive(input$sig.height),
      filename.stem = reactive(paste0("signaling_", input$sig.plot)),
      trigger = reactive(input$sig.plot.btn))

    # =====================================================
    # Subtab 4: Communication patterns
    # =====================================================

    pat.thunk <- function() {
      res <- cellchat.data(); req(res, input$pat.plot)

      if (input$pat.plot == "manifold") {
        req(input$pat.sim.type)
        merged <- res$cellchat.merged
        tryCatch({
          merged <- computeNetSimilarityPairwise(merged, type = input$pat.sim.type)
          merged <- netEmbedding(merged, type = input$pat.sim.type, umap.method = "uwot")
          merged <- netClustering(merged, type = input$pat.sim.type, do.parallel = FALSE)
          print(netVisual_embeddingPairwise(merged, type = input$pat.sim.type, label.size = 3.5))
        }, error = function(e) {
          plot.new(); title(paste("Manifold failed:", conditionMessage(e)))
        })
        return(invisible())
      }

      req(input$pat.k, input$pat.group)
      g <- input$pat.group
      k <- input$pat.k
      ncol <- input$pat.cols %||% 2

      cc <- compute.centrality(res$cellchat.list[[g]]); req(cc)
      cc.out <- identifyCommunicationPatterns(cc, pattern = "outgoing", k = k,
                                              width = 5, height = 9)
      cc.in  <- identifyCommunicationPatterns(cc, pattern = "incoming", k = k,
                                              width = 5, height = 9)

      label <- c(outgoing = paste0(g, " — outgoing"),
                 incoming = paste0(g, " — incoming"))

      items <- if (input$pat.plot == "heat") {
        list(
          ComplexHeatmap::Heatmap(cc.out@netP$pattern$outgoing$pattern$cell,
            name = "outgoing", column_title = label[["outgoing"]]),
          ComplexHeatmap::Heatmap(cc.in@netP$pattern$incoming$pattern$cell,
            name = "incoming", column_title = label[["incoming"]])
        )
      } else if (input$pat.plot == "river") {
        list(
          netAnalysis_river(cc.out, pattern = "outgoing") + ggtitle(label[["outgoing"]]),
          netAnalysis_river(cc.in,  pattern = "incoming") + ggtitle(label[["incoming"]])
        )
      } else if (input$pat.plot == "dot") {
        list(
          netAnalysis_dot(cc.out, pattern = "outgoing") + ggtitle(label[["outgoing"]]),
          netAnalysis_dot(cc.in,  pattern = "incoming") + ggtitle(label[["incoming"]])
        )
      }
      plot.grid(items, ncol, titles = c("outgoing", "incoming"))
    }

    make.render.and.download(output, session, "pat.plot", "pat.download",
      pat.thunk,
      width.in = reactive(input$pat.width),
      height.in = reactive(input$pat.height),
      filename.stem = reactive(paste0("patterns_", input$pat.plot)),
      trigger = reactive(input$pat.plot.btn))

  })
}
