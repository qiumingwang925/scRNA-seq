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
          verbatimTextOutput(ns("input.summary")),
          hr(),
          h5("Cell types (applies to every tab below)"),
          cell.type.selector.ui(ns, "subset.idents", "Cell types to include:")
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
            common.controls.ui(ns, "global", default.cols = 1)
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
            conditionalPanel(
              condition = sprintf("input['%s'] != 'bubble'", ns("zoom.plot")),
              selectInput(ns("zoom.pathway"), "Pathway:", choices = NULL)
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'bubble'", ns("zoom.plot")),
              selectInput(ns("zoom.pathways"), "Pathways:", choices = NULL,
                          multiple = TRUE)
            ),
            conditionalPanel(
              condition = sprintf("['lr.hier','lr.circle','lr.chord'].indexOf(input['%s']) >= 0",
                                  ns("zoom.plot")),
              selectInput(ns("zoom.lr"), "L-R pair:", choices = NULL)
            ),
            conditionalPanel(
              condition = sprintf("['pw.hier','lr.hier'].indexOf(input['%s']) >= 0",
                                  ns("zoom.plot")),
              selectInput(ns("zoom.receivers"), "Vertex Receiver:",
                          choices = NULL, multiple = TRUE)
            ),
            conditionalPanel(
              condition = sprintf(
                "['pw.circle','pw.chord','lr.circle','lr.chord','bubble'].indexOf(input['%s']) >= 0",
                ns("zoom.plot")),
              cell.type.selector.ui(ns, "zoom.sources", "Sources:"),
              cell.type.selector.ui(ns, "zoom.targets", "Targets:"),
              checkboxInput(ns("zoom.remove.isolate"), "Remove isolate", value = FALSE)
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'violin'", ns("zoom.plot")),
              checkboxInput(ns("zoom.enriched.only"), "Enriched genes only", value = TRUE)
            ),
            common.controls.ui(ns, "zoom", default.cols = 1)
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
              condition = sprintf("['score','heatmap'].indexOf(input['%s']) >= 0",
                                  ns("sig.plot")),
              selectInput(ns("sig.pathway"), "Pathways:", choices = NULL, multiple = TRUE)
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'heatmap'", ns("sig.plot")),
              radioButtons(ns("sig.pattern"), "Direction (heatmap only):",
                           choices = c("Outgoing" = "outgoing",
                                       "Incoming" = "incoming",
                                       "All" = "all"),
                           selected = "all")
            ),
            common.controls.ui(ns, "sig", default.cols = 1)
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
            # "heat" (pattern heatmap) and "manifold" are temporarily hidden pending a
            # fix for deeper issues; their server branches below are kept for re-enabling.
            radioButtons(ns("pat.plot"), "Plot type:",
                         choices = c("River" = "river",
                                     "Dot" = "dot"),
                         selected = "river"),
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
            common.controls.ui(ns, "pat", default.cols = 1)
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

    # The cell types the user kept, in ident.choices() order rather than click
    # order, so "first half" defaults and index-based arguments stay stable.
    subset.universe <- reactive({
      all.idents <- ident.choices()
      req(all.idents)
      intersect(all.idents, input$subset.idents %||% all.idents)
    })

    # App-wide subset shared by all four subtabs. Lazy on purpose: subsetCellChat
    # re-derives pathway aggregation per group, so the cost lands on the first
    # "Generate Plot" after a selection change rather than on every click in the
    # selector. status[[group]] is NA when that group rendered fine, or the message
    # a panel should show instead of a plot.
    cellchat.subset <- reactive({
      res <- cellchat.data(); req(res)
      sel <- subset.universe()
      validate(need(length(sel) > 0,
                    "Select at least one cell type in the panel at the top of this tab."))
      grps <- res$group.levels
      status <- setNames(rep(NA_character_, length(grps)), grps)

      if (setequal(sel, ident.choices())) {
        return(list(cellchat.list = res$cellchat.list,
                    cellchat.merged = res$cellchat.merged,
                    status = status))
      }

      subs <- setNames(vector("list", length(grps)), grps)
      withProgress(message = "Subsetting CellChat objects", value = 0, {
        for (g in grps) {
          incProgress(1 / length(grps), detail = g)
          cc <- res$cellchat.list[[g]]
          present <- intersect(sel, levels(cc@idents))
          if (length(present) == 0) {
            status[[g]] <- "selected cell types absent in this group"
          } else if (length(present) < 2) {
            status[[g]] <- "select >= 2 cell types present in this group"
          } else if (setequal(present, levels(cc@idents))) {
            # Nothing to drop for this group; a no-op subsetCellChat recompute
            # errors on some objects, so skip it.
            subs[[g]] <- cc
          } else {
            out <- tryCatch({
              sub <- subsetCellChat(cc, idents.use = present)
              # subsetCellChat recomputes @netP$centr from the L-R-level prob,
              # leaving it indexed by L-R pair instead of pathway. The
              # signaling-role plots expect pathway-indexed centrality and
              # otherwise crash ("'from' must be a finite number").
              sub@netP$centr <- netAnalysis_computeCentrality(net = sub@netP$prob)
              sub
            }, error = function(e) paste("subset failed:", conditionMessage(e)))
            if (is.character(out)) status[[g]] <- out else subs[[g]] <- out
          }
        }
      })

      # Rebuild rather than subset the merged object: mergeCellChat is the
      # supported construction path, subsetCellChat is not written for merges.
      ok <- names(subs)[!vapply(subs, is.null, logical(1))]
      merged <- if (length(ok) >= 2)
        tryCatch(mergeCellChat(subs[ok], add.names = ok), error = function(e) NULL) else NULL

      list(cellchat.list = subs, cellchat.merged = merged, status = status)
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
      updateSelectInput(session, "zoom.pathways", choices = pw,
                        selected = if (length(pw)) pw[1] else NULL)
      # Multi-select shared by heatmap (empty = all pathways) and score (empty =
      # prompt for a pathway), so default to no selection rather than the first.
      updateSelectInput(session, "sig.pathway", choices = pw, selected = character(0))
    })

    observe({
      lr <- lr.choices()
      updateSelectInput(session, "zoom.lr", choices = lr,
                        selected = if (length(lr)) lr[1] else NULL)
    })

    wire.cell.type.selector(input, session, "subset.idents", ident.choices, selected.all = TRUE)
    wire.cell.type.selector(input, session, "global.sources", subset.universe, selected.all = TRUE)
    wire.cell.type.selector(input, session, "global.targets", subset.universe, selected.all = TRUE)
    wire.cell.type.selector(input, session, "zoom.sources", subset.universe, selected.all = TRUE)
    wire.cell.type.selector(input, session, "zoom.targets", subset.universe, selected.all = TRUE)

    # Receivers get no Select All button: selecting every cell type leaves the
    # hierarchy plot with no sender side. Default is the first half.
    observe({
      u <- subset.universe()
      req(length(u) > 0)
      updateSelectInput(session, "zoom.receivers", choices = u,
                        selected = u[seq_len(ceiling(length(u) / 2))])
    })

    observeEvent(input$zoom.plot, {
      updateCheckboxInput(session, "zoom.remove.isolate",
                          value = identical(input$zoom.plot, "bubble"))
    })

    observe({
      res <- cellchat.data()
      req(res)
      updateSelectInput(session, "pat.group",
                        choices = res$group.levels,
                        selected = res$group.levels[1])
    })


    # Number of panels each tab draws, so the canvas can grow with the grid rather
    # than cramming every panel into a fixed height (which squeezed titles at cols=1).
    n.panels.global <- reactive({ res <- cellchat.data(); req(res); length(res$group.levels) })
    n.panels.zoom   <- reactive({ res <- cellchat.data(); req(res); length(res$group.levels) })
    n.panels.sig <- reactive({
      res <- cellchat.data(); req(res)
      if (identical(input$sig.plot, "score"))
        length(res$group.levels) * max(1L, length(nz(input$sig.pathway)))
      else length(res$group.levels)
    })
    n.panels.pat <- reactive({ if (identical(input$pat.plot, "manifold")) 1L else 2L })

    grid.cols <- function(prefix) max(1L, input[[paste0(prefix, ".cols")]] %||% 1L)

    # Width/Height sliders are PER PANEL: total canvas = width x cols by height x rows,
    # so each panel keeps its full size no matter how the grid wraps.
    total.w <- function(prefix) reactive(input[[paste0(prefix, ".width")]] * grid.cols(prefix))
    total.h <- function(prefix, n.panels) reactive(
      input[[paste0(prefix, ".height")]] * ceiling(n.panels() / grid.cols(prefix)))

    plot.output <- function(prefix, n.panels) {
      renderUI({
        w <- input[[paste0(prefix, ".width")]]
        h <- input[[paste0(prefix, ".height")]]
        req(w, h)
        cols <- grid.cols(prefix)
        rows <- ceiling(n.panels() / cols)
        shinycssloaders::withSpinner(
          plotOutput(ns(paste0(prefix, ".plot")),
                     width = paste0(w * cols * 96, "px"),
                     height = paste0(h * rows * 96, "px")),
          type = 6
        )
      })
    }

    output$global.plot.ui <- plot.output("global", n.panels.global)
    output$zoom.plot.ui <- plot.output("zoom", n.panels.zoom)
    output$sig.plot.ui <- plot.output("sig", n.panels.sig)
    output$pat.plot.ui <- plot.output("pat", n.panels.pat)

    # =====================================================
    # Subtab 1: Global network
    # =====================================================

    global.thunk <- function() {
      res <- cellchat.data(); req(res, input$global.plot, input$global.measure)
      sub.res <- cellchat.subset(); req(sub.res)
      grps <- res$group.levels
      ncol <- input$global.cols %||% 1
      all.idents <- subset.universe()
      srcs <- input$global.sources; tgts <- input$global.targets
      measure <- input$global.measure

      cc.ok <- Filter(Negate(is.null), sub.res$cellchat.list)
      req(length(cc.ok) > 0)
      weight.max <- tryCatch(
        getMaxWeight(cc.ok, attribute = c("idents", measure)),
        error = function(e) NULL
      )
      vertex.weight.max <- max(unlist(lapply(cc.ok, function(cc) {
        as.numeric(table(cc@idents))
      })))

      items <- lapply(grps, function(g) {
        force(g)
        subset.msg <- sub.res$status[[g]]
        if (!is.na(subset.msg)) {
          return(if (input$global.plot == "circle")
            function() placeholder.base(g, subset.msg) else placeholder.ht(g, subset.msg))
        }
        cc <- sub.res$cellchat.list[[g]]
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
      width.in = total.w("global"),
      height.in = total.h("global", n.panels.global),
      filename.stem = reactive(paste0("global_", input$global.plot, "_", input$global.measure)),
      trigger = reactive(input$global.plot.btn))

    # =====================================================
    # Subtab 2: Zoom-in (pathway & LR)
    # =====================================================

    zoom.thunk <- function() {
      res <- cellchat.data(); req(res, input$zoom.plot)
      sub.res <- cellchat.subset(); req(sub.res)
      grps <- res$group.levels
      ncol <- input$zoom.cols %||% 1
      pt <- input$zoom.plot
      is.lr <- startsWith(pt, "lr.")
      is.hier <- pt %in% c("pw.hier", "lr.hier")
      # conditionalPanel only hides inputs, so a stale sources/targets selection is
      # still readable here — only honour it for the plot types that expose it.
      uses.st <- pt %in% c("pw.circle", "pw.chord", "lr.circle", "lr.chord", "bubble")

      pw <- if (pt == "bubble") nz(input$zoom.pathways) else input$zoom.pathway
      validate(need(length(pw) > 0, "Select at least one pathway."))
      lr <- input$zoom.lr
      srcs <- nz(input$zoom.sources); tgts <- nz(input$zoom.targets)
      rcvs <- nz(input$zoom.receivers)
      all.idents <- subset.universe()
      remove.isolate <- isTRUE(input$zoom.remove.isolate)
      layout.map <- c(pw.hier = "hierarchy", pw.circle = "circle", pw.chord = "chord",
                      lr.hier = "hierarchy", lr.circle = "circle", lr.chord = "chord")

      # Per-group status: the absence message (if any) for this group
      group.status <- function(cc) {
        lv <- levels(cc@idents)
        s <- resolve.sel(lv, srcs, all.idents)
        t <- resolve.sel(lv, tgts, all.idents)
        # vertex.receiver indexes THIS group's ident levels, so re-match per group
        # rather than reusing positions from the selector's own ordering.
        rcv.idx <- sort(match(intersect(rcvs, lv), lv))
        pw.use <- intersect(pw, cc@netP$pathways)
        msg <- NULL
        if (uses.st && (s$empty || t$empty)) msg <- "selection absent"
        else if (length(pw.use) == 0)
          msg <- if (length(pw) > 1) "pathways absent" else "pathway absent"
        else if (is.hier && length(rcv.idx) == 0) msg <- "vertex receiver absent"
        else if (is.hier && length(rcv.idx) == length(lv))
          msg <- "vertex receiver covers all cell types"
        else if (is.lr) {
          enriched <- tryCatch(
            extractEnrichedLR(cc, signaling = pw.use, geneLR.return = FALSE)$interaction_name,
            error = function(e) character()
          )
          if (is.null(lr) || !(lr %in% enriched)) msg <- "L-R pair absent"
        }
        list(s = s, t = t, rcv.idx = rcv.idx, pw.use = pw.use, msg = msg)
      }

      items <- lapply(grps, function(g) {
        force(g)
        subset.msg <- sub.res$status[[g]]
        cc <- sub.res$cellchat.list[[g]]
        st <- if (is.na(subset.msg)) group.status(cc) else
          list(s = NULL, t = NULL, rcv.idx = integer(0), pw.use = pw, msg = subset.msg)
        s <- st$s; t <- st$t; missing.msg <- st$msg

        if (pt == "bubble") {
          if (!is.null(missing.msg))
            return(placeholder.gg(g, missing.msg))
          return(netVisual_bubble(cc, signaling = st$pw.use,
                                  sources.use = if (s$active) s$use else NULL,
                                  targets.use = if (t$active) t$use else NULL,
                                  remove.isolate = remove.isolate) + ggtitle(g))
        }

        if (pt == "violin") {
          if (!is.null(missing.msg))
            return(placeholder.gg(g, missing.msg))
          return(plotGeneExpression(cc, signaling = st$pw.use, type = "violin",
                                    enriched.only = isTRUE(input$zoom.enriched.only)) +
                   patchwork::plot_annotation(title = g))
        }

        if (pt == "pw.heat") {
          if (!is.null(missing.msg))
            return(placeholder.ht(g, missing.msg))
          return(netVisual_heatmap(cc, signaling = st$pw.use,
                                   color.heatmap = "Reds", title.name = g))
        }

        # Base-graphics: circle / chord / hierarchy for pw or lr
        function() {
          if (!is.null(missing.msg)) {
            plot.new(); title(paste0(g, "\n(", missing.msg, ")")); return(invisible())
          }
          layout <- layout.map[[pt]]
          common <- list(object = cc, signaling = st$pw.use, layout = layout)
          if (layout == "hierarchy") {
            common$vertex.receiver <- st$rcv.idx
            # Hierarchy needs both sides drawn, so isolates always stay.
            common$remove.isolate <- FALSE
          } else {
            if (s$active) common$sources.use <- s$use
            if (t$active) common$targets.use <- t$use
            common$remove.isolate <- remove.isolate
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
      width.in = total.w("zoom"),
      height.in = total.h("zoom", n.panels.zoom),
      filename.stem = reactive(paste0(
        "zoom_", input$zoom.plot, "_",
        if (identical(input$zoom.plot, "bubble"))
          paste(nz(input$zoom.pathways), collapse = "-") else input$zoom.pathway)),
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
      sub.res <- cellchat.subset(); req(sub.res)
      grps <- res$group.levels
      ncol <- max(1, input$sig.cols %||% 1)
      sig.pw <- nz(input$sig.pathway)

      # Score is inherently per-pathway; without a pathway there is nothing to draw.
      if (input$sig.plot == "score" && is.null(sig.pw)) {
        return(plot.grid(list(placeholder.gg("Score", "select ≥ 1 pathway")), 1))
      }

      # CellChat's heatmap/network draw at absolute cm sizes, so scale them to the
      # per-panel figure size (each panel is width x height inches). The fractions
      # leave room for row/column labels, legend, and title.
      panel.w.cm <- (input$sig.width %||% 10) * 2.54
      panel.h.cm <- (input$sig.height %||% 6) * 2.54

      items <- list(); titles <- character()
      add <- function(item, title) {
        items[[length(items) + 1]] <<- item
        titles[[length(titles) + 1]] <<- title
      }

      for (g in grps) {
        subset.msg <- sub.res$status[[g]]
        if (!is.na(subset.msg)) {
          if (input$sig.plot == "score") {
            for (pw in sig.pw) add(placeholder.gg(g, subset.msg), paste0(g, " — ", pw))
          } else if (input$sig.plot == "heatmap") {
            add(placeholder.ht(g, subset.msg), g)
          } else {
            add(placeholder.gg(g, subset.msg), g)
          }
          next
        }
        cc <- compute.centrality(sub.res$cellchat.list[[g]])

        if (input$sig.plot == "scatter") {
          add(tryCatch(netAnalysis_signalingRole_scatter(cc) + ggtitle(g),
                       error = function(e) placeholder.gg(g, conditionMessage(e))), g)
        } else if (input$sig.plot == "heatmap") {
          pw.use <- if (is.null(sig.pw)) NULL else intersect(sig.pw, cc@netP$pathways)
          if (!is.null(sig.pw) && length(pw.use) == 0) {
            add(placeholder.ht(g, "selected pathways absent"), g)
          } else {
            add(tryCatch(netAnalysis_signalingRole_heatmap(
                  cc, signaling = pw.use, pattern = input$sig.pattern,
                  width = panel.w.cm * 0.5, height = panel.h.cm * 0.65, title = g),
                  error = function(e) placeholder.ht(g, conditionMessage(e))), g)
          }
        } else {  # score: one panel per selected pathway
          for (pw in sig.pw) local({
            pw. <- pw; cc. <- cc; g. <- g
            title <- paste0(g., " — ", pw.)
            if (!pw. %in% cc.@netP$pathways) {
              add(placeholder.gg(g., paste0(pw., " absent")), title)
            } else {
              add(function() {
                netAnalysis_signalingRole_network(
                  cc., signaling = pw.,
                  width = panel.w.cm * 0.6, height = panel.h.cm * 0.35, font.size = 10)
              }, title)
            }
          })
        }
      }
      plot.grid(items, ncol, titles = titles)
    }

    make.render.and.download(output, session, "sig.plot", "sig.download",
      sig.thunk,
      width.in = total.w("sig"),
      height.in = total.h("sig", n.panels.sig),
      filename.stem = reactive(paste0("signaling_", input$sig.plot)),
      trigger = reactive(input$sig.plot.btn))

    # =====================================================
    # Subtab 4: Communication patterns
    # =====================================================

    pat.thunk <- function() {
      res <- cellchat.data(); req(res, input$pat.plot)
      sub.res <- cellchat.subset(); req(sub.res)

      if (input$pat.plot == "manifold") {
        req(input$pat.sim.type)
        merged <- sub.res$cellchat.merged
        if (is.null(merged)) {
          plot.new()
          title("Manifold needs >= 2 groups surviving the cell-type selection")
          return(invisible())
        }
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
      ncol <- input$pat.cols %||% 1

      subset.msg <- sub.res$status[[g]]
      if (!is.na(subset.msg)) return(plot.grid(list(placeholder.gg(g, subset.msg)), 1))
      cc <- compute.centrality(sub.res$cellchat.list[[g]]); req(cc)

      n.types <- length(levels(cc@idents))
      if (k >= n.types) {
        return(plot.grid(list(placeholder.gg(g,
          paste0("k (", k, ") must be < number of cell types (", n.types, ")"))), 1))
      }

      cc.out <- tryCatch(identifyCommunicationPatterns(cc, pattern = "outgoing", k = k,
                                                       width = 5, height = 9),
                         error = function(e) NULL)
      cc.in  <- tryCatch(identifyCommunicationPatterns(cc, pattern = "incoming", k = k,
                                                       width = 5, height = 9),
                         error = function(e) NULL)
      if (is.null(cc.out) || is.null(cc.in)) {
        return(plot.grid(list(placeholder.gg(g,
          "pattern identification failed for this subset")), 1))
      }

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
      width.in = total.w("pat"),
      height.in = total.h("pat", n.panels.pat),
      filename.stem = reactive(paste0("patterns_", input$pat.plot)),
      trigger = reactive(input$pat.plot.btn))

  })
}
