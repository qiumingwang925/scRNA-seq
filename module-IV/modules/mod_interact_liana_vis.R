## ABOUTME: LIANA visualization module. A shared cell-type + method panel feeds four subtabs —
## ABOUTME: CCC Table (one raw method table, or liana_aggregate over a chosen subset), Dot Plot,
## ABOUTME: Freq Heatmap, Freq Chord Diagram. The table is the source of truth for the plots.

# liana_aggregate only has rank specs for these five. cytotalk and call_cellchat return a
# much smaller row support, so including them silently collapses the aggregate onto their
# interactions instead of raising — hence the hard block rather than a warning.
.liana.aggregatable <- c("connectome", "logfc", "natmi", "sca", "cellphonedb")

# Identity columns — never offered as a score, cutoff, or ranking metric.
.liana.id.cols <- c("source", "target",
                    "ligand", "ligand.complex",
                    "receptor", "receptor.complex")

# Columns where a *small* value means a stronger interaction. Drives the default cutoff
# direction, the top-N sort order, and the -log10 inversion in the dot plot.
.liana.is.pval.col <- function(cols) grepl("aggregate_rank|pvalue|padj|pval", cols)

.liana.score.cols <- function(tib) {
  numeric.cols <- names(tib)[vapply(tib, is.numeric, logical(1))]
  setdiff(numeric.cols, .liana.id.cols)
}

# First preference present in cols, else the fallback.index-th column.
.liana.pick.col <- function(cols, prefs, fallback.index = 1) {
  hit <- prefs[prefs %in% cols]
  if (length(hit) > 0) return(hit[1])
  if (length(cols) == 0) return(NULL)
  cols[min(fallback.index, length(cols))]
}

# liana_dotplot draws one row per ligand-receptor pair, so "top N" counts distinct pairs
# rather than table rows: rank the table, then keep every row belonging to the first N
# pairs so each pair keeps all of its source -> target combinations.
.liana.top.pairs <- function(tib, rank.col, n) {
  ordered <- tib[order(tib[[rank.col]],
                       decreasing = !.liana.is.pval.col(rank.col),
                       na.last = NA), ]
  pair.cols <- intersect(c("ligand.complex", "receptor.complex"), names(ordered))
  if (length(pair.cols) < 2) {
    pair.cols <- intersect(c("ligand", "receptor"), names(ordered))
  }
  if (length(pair.cols) < 2) return(utils::head(ordered, n))
  key <- do.call(paste, c(ordered[pair.cols], sep = "|"))
  ordered[key %in% utils::head(unique(key), n), ]
}

# How a method selection is turned into a table: one method reads its raw table, several
# aggregatable ones are combined. Returns list(ok, aggregate, msg); msg is the reason the
# selection cannot produce a table. Called on both the pending and the applied selection.
.liana.method.mode <- function(sel) {
  if (is.null(sel) || length(sel) == 0) {
    return(list(ok = FALSE, aggregate = FALSE,
                msg = "Select at least one method."))
  }
  if (length(sel) == 1) {
    return(list(ok = TRUE, aggregate = FALSE, msg = NULL))
  }
  blocked <- setdiff(sel, .liana.aggregatable)
  if (length(blocked) > 0) {
    return(list(ok = FALSE, aggregate = FALSE, msg = paste0(
      "Cannot aggregate ", paste(blocked, collapse = " and "), ". ",
      "liana_aggregate only supports ",
      paste(.liana.aggregatable, collapse = ", "),
      " — select ", paste(blocked, collapse = " or "),
      " on its own to see its raw table, or remove it from the selection.")))
  }
  list(ok = TRUE, aggregate = TRUE, msg = NULL)
}

# Returns the filtered table, or a character message explaining why it could not filter.
.liana.apply.cutoff <- function(tib, col, op, value) {
  if (is.null(col) || !nzchar(col)) return("choose a cutoff column")
  if (!col %in% names(tib)) {
    return(paste0("cutoff column '", col, "' absent in this group"))
  }
  if (is.null(value) || is.na(value)) return("enter a cutoff value")
  x <- tib[[col]]
  keep <- switch(op, gt = x > value, ge = x >= value,
                     lt = x < value, le = x <= value)
  if (is.null(keep)) return("unknown cutoff operator")
  tib[!is.na(keep) & keep, ]
}

liana.cutoff.controls.ui <- function(ns, prefix) {
  tagList(
    selectInput(ns(paste0(prefix, ".cutoff.col")), "Cutoff column:", choices = NULL),
    selectInput(ns(paste0(prefix, ".cutoff.op")), "Cutoff:",
                choices = c("greater than" = "gt",
                            "greater than or equal to" = "ge",
                            "less than" = "lt",
                            "less than or equal to" = "le"),
                selected = "le"),
    numericInput(ns(paste0(prefix, ".cutoff.val")), NULL, value = 0.01)
  )
}

mod.interact.liana.vis.ui <- function(id) {
  ns <- NS(id)
  tabPanel("Visualize LIANA",
    fluidRow(
      column(12,
        wellPanel(
          fileInput(ns("vis.file"),
                    "Upload processed LIANA result (.rds) — skip if computed above:",
                    accept = ".rds"),
          verbatimTextOutput(ns("input.summary")),
          hr(),
          h5("Applies to every tab below"),
          fluidRow(
            column(6,
              cell.type.selector.ui(ns, "subset.idents", "Cell types to include:")
            ),
            column(6,
              selectInput(ns("methods"), "Method(s):", choices = NULL, multiple = TRUE),
              helpText("One method shows its raw table. Two or more of ",
                       paste(.liana.aggregatable, collapse = ", "),
                       " are combined with liana_aggregate.",
                       "cytotalk and call_cellchat can only be viewed on their own."),
              uiOutput(ns("method.msg"))
            )
          ),
          hr(),
          fluidRow(
            column(3,
              actionButton(ns("apply.selection"), "Apply Selection",
                           class = "btn-success", style = "width:100%")
            ),
            column(9, uiOutput(ns("apply.msg")))
          )
        )
      )
    ),

    tabsetPanel(id = ns("liana.vis.tabs"),

      # ---------- Subtab 1: CCC Table ----------
      tabPanel("CCC Table",
        sidebarLayout(
          sidebarPanel(width = 3,
            h4("CCC Table"),
            selectInput(ns("table.group"), "Group:", choices = NULL),
            helpText("This group also drives the CCC Dot Plot. Rows selected here",
                     "can be plotted directly."),
            hr(),
            downloadButton(ns("table.download"), "Download Table (.csv)",
                           style = "width:100%")
          ),
          mainPanel(width = 9,
            verbatimTextOutput(ns("table.caption")),
            DT::DTOutput(ns("ccc.table"))
          )
        )
      ),

      # ---------- Subtab 2: CCC Dot Plot (the table's group) ----------
      tabPanel("CCC Dot Plot",
        sidebarLayout(
          sidebarPanel(width = 3,
            h4("CCC Dot Plot"),
            radioButtons(ns("dot.rows"), "L-R pairs:",
                         choices = c("Top N" = "top",
                                     "Selected in CCC Table" = "selected"),
                         selected = "top"),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'top'", ns("dot.rows")),
              numericInput(ns("dot.top"), "Top N L-R pairs:",
                           value = 15, min = 1, max = 200, step = 1),
              selectInput(ns("dot.rank"), "Rank by:", choices = NULL),
              cell.type.selector.ui(ns, "dot.sources", "Source populations:"),
              cell.type.selector.ui(ns, "dot.targets", "Target populations:")
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'selected'", ns("dot.rows")),
              helpText("Every source -> target combination present in the selected",
                       "rows is drawn; the population filters above do not apply.")
            ),
            selectInput(ns("dot.size"), "Dot size:", choices = NULL),
            selectInput(ns("dot.col"), "Dot colour:", choices = NULL),
            common.controls.ui(ns, "dot", default.cols = 1)
          ),
          mainPanel(width = 9,
            uiOutput(ns("dot.plot.ui"))
          )
        )
      ),

      # ---------- Subtab 3: CCC Freq Heatmap (N-group grid) ----------
      tabPanel("CCC Freq Heatmap",
        sidebarLayout(
          sidebarPanel(width = 3,
            h4("CCC Freq Heatmap"),
            liana.cutoff.controls.ui(ns, "heat"),
            common.controls.ui(ns, "heat", default.cols = 1)
          ),
          mainPanel(width = 9,
            uiOutput(ns("heat.plot.ui"))
          )
        )
      ),

      # ---------- Subtab 4: CCC Freq Chord Diagram (N-group grid) ----------
      tabPanel("CCC Freq Chord Diagram",
        sidebarLayout(
          sidebarPanel(width = 3,
            h4("CCC Freq Chord Diagram"),
            cell.type.selector.ui(ns, "chord.sources", "Source populations:"),
            cell.type.selector.ui(ns, "chord.targets", "Target populations:"),
            liana.cutoff.controls.ui(ns, "chord"),
            common.controls.ui(ns, "chord", default.cols = 1)
          ),
          mainPanel(width = 9,
            uiOutput(ns("chord.plot.ui"))
          )
        )
      )
    )
  )
}

mod.interact.liana.vis.server <- function(id, liana.input) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Resolve LIANA result: upstream compute first, then file upload
    liana.data <- reactive({
      upstream <- NULL
      try(upstream <- liana.input(), silent = TRUE)
      if (!is.null(upstream) && is.list(upstream) &&
          all(c("liana.list", "group.levels") %in% names(upstream))) {
        return(upstream)
      }

      req(input$vis.file)
      tryCatch({
        r <- readRDS(input$vis.file$datapath)
        if (!is.list(r) || !all(c("liana.list", "group.levels") %in% names(r))) {
          showNotification("Uploaded file is not a valid LIANA result.",
                           type = "error", duration = NULL)
          return(NULL)
        }
        if (inherits(r$liana.list[[1]], "data.frame")) {
          showNotification(
            paste("This file holds one pre-aggregated table per group, the format",
                  "produced before the per-method rewrite. Re-run the LIANA compute",
                  "tab to produce a result with one table per method."),
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
      res <- liana.data()
      req(res)
      cat("LIANA result loaded.\n")
      cat("Groups:   ", paste(res$group.levels, collapse = ", "), "\n")
      cat("Species:  ", res$species %||% "unknown", "\n")
      cat("Assay:    ", res$assay %||% "unknown", "\n")
      cat("Resource: ", res$resource %||% "unknown", "\n")
      cat("Methods per group:\n")
      for (g in res$group.levels) {
        per.method <- res$liana.list[[g]]
        cat("  ", g, ": ",
            paste0(names(per.method), " (",
                   vapply(per.method, nrow, integer(1)), " rows)",
                   collapse = ", "),
            "\n", sep = "")
      }
    })

    # ---------- Shared selections ----------

    ident.choices <- reactive({
      res <- liana.data(); req(res)
      all <- Reduce(union, lapply(res$liana.list, function(per.method) {
        Reduce(union, lapply(per.method, function(tib) {
          unique(c(tib$source, tib$target))
        }))
      }))
      sort(all[!is.na(all) & nzchar(all)])
    })

    # The shared selection is staged: the selectors are free to hold any interim state
    # while the user edits them, and only these two values — updated on Apply — reach
    # the table, the aggregation cache, and the plot thunks. Without the stage, every
    # click re-ran liana_aggregate, re-rendered the DT, and reset the per-tab
    # source/target selectors, none of which the user had asked for yet.
    applied.methods <- reactiveVal(NULL)
    applied.idents <- reactiveVal(NULL)

    # Kept cell types in ident.choices() order rather than click order.
    subset.universe <- reactive({
      all.idents <- ident.choices()
      req(all.idents)
      intersect(all.idents, applied.idents() %||% all.idents)
    })

    # Only methods that succeeded in *every* group are selectable: aggregating over
    # different method sets per group would make the resulting ranks incomparable
    # across conditions, which is the whole point of the N-panel grids.
    method.choices <- reactive({
      res <- liana.data(); req(res)
      Reduce(intersect, lapply(res$liana.list, names))
    })

    method.partial <- reactive({
      res <- liana.data(); req(res)
      setdiff(Reduce(union, lapply(res$liana.list, names)), method.choices())
    })

    # Seeding the applied values alongside the selectors means a freshly loaded result
    # renders straight away instead of waiting for a first Apply click.
    observe({
      methods <- method.choices()
      req(length(methods) > 0)
      default <- intersect(methods, .liana.aggregatable)
      if (length(default) < 2) default <- methods[1]
      updateSelectInput(session, "methods", choices = methods, selected = default)
      applied.methods(default)
    })

    observe({
      choices <- ident.choices()
      req(choices)
      applied.idents(choices)
    })

    observe({
      grps <- liana.data()$group.levels
      req(grps)
      updateSelectInput(session, "table.group", choices = grps, selected = grps[1])
    })

    wire.cell.type.selector(input, session, "subset.idents", ident.choices,
                            selected.all = TRUE)
    wire.cell.type.selector(input, session, "dot.sources", subset.universe,
                            selected.all = TRUE)
    wire.cell.type.selector(input, session, "dot.targets", subset.universe,
                            selected.all = TRUE)
    wire.cell.type.selector(input, session, "chord.sources", subset.universe,
                            selected.all = TRUE)
    wire.cell.type.selector(input, session, "chord.targets", subset.universe,
                            selected.all = TRUE)

    method.mode <- reactive(.liana.method.mode(applied.methods()))

    # Validation is judged on the *pending* selection so the user learns an interim
    # state is unusable while editing, rather than after committing to it.
    apply.blocked <- reactive({
      pending <- .liana.method.mode(input$methods)
      if (!pending$ok) return(pending$msg)
      if (length(input$subset.idents %||% character(0)) == 0) {
        return("Select at least one cell type.")
      }
      NULL
    })

    pending.changes <- reactive({
      !setequal(input$methods %||% character(0),
                applied.methods() %||% character(0)) ||
      !setequal(input$subset.idents %||% character(0),
                applied.idents() %||% character(0))
    })

    observe({
      shinyjs::toggleState(
        "apply.selection",
        condition = is.null(apply.blocked()) && isTRUE(pending.changes()))
    })

    observeEvent(input$apply.selection, {
      applied.methods(input$methods)
      applied.idents(input$subset.idents)
    })

    output$method.msg <- renderUI({
      partial <- method.partial()
      if (length(partial) == 0) return(NULL)
      tags$p(class = "text-warning", paste0(
        "Not selectable — succeeded in some groups only: ",
        paste(partial, collapse = ", "), "."))
    })

    output$apply.msg <- renderUI({
      blocked <- apply.blocked()
      if (!is.null(blocked)) return(tags$p(class = "text-danger", blocked))
      if (isTRUE(pending.changes())) {
        tags$p(class = "text-info",
               "Selection changed — click Apply Selection to update every tab.")
      } else {
        tags$p(class = "text-muted", "Every tab matches the current selection.")
      }
    })

    # ---------- Table construction ----------

    # liana_aggregate over the full group table is expensive enough to be worth
    # keeping, and cheap enough to run on demand. Keyed by group + method set.
    agg.cache <- new.env(parent = emptyenv())
    observeEvent(liana.data(), {
      rm(list = ls(envir = agg.cache), envir = agg.cache)
    }, ignoreNULL = TRUE)

    # Raw or aggregated table for one group, before the cell-type filter.
    # Returns a character message instead of a table when it cannot be built.
    group.table <- function(g) {
      res <- liana.data()
      mode <- method.mode()
      if (!mode$ok) return(mode$msg)
      sel <- applied.methods()
      per.method <- res$liana.list[[g]]
      if (is.null(per.method)) return(paste0("no result for group '", g, "'"))
      missing <- setdiff(sel, names(per.method))
      if (length(missing) > 0) {
        return(paste0("method(s) absent in this group: ",
                      paste(missing, collapse = ", ")))
      }
      if (!mode$aggregate) return(per.method[[sel]])

      key <- paste(c(g, sort(sel)), collapse = "|")
      if (!exists(key, envir = agg.cache)) {
        aggregated <- withProgress(
          message = paste("Aggregating", length(sel), "methods for", g),
          value = 0.5,
          tryCatch(liana::liana_aggregate(per.method[sort(sel)]),
                   error = function(e) paste("liana_aggregate failed:",
                                             conditionMessage(e)))
        )
        assign(key, aggregated, envir = agg.cache)
      }
      get(key, envir = agg.cache)
    }

    # Aggregation runs on the full group table so that ranks do not shift when the
    # cell-type selection changes; the filter is applied to the rows afterwards.
    table.for.group <- function(g) {
      tib <- group.table(g)
      if (is.character(tib)) return(tib)
      keep <- subset.universe()
      if (setequal(keep, ident.choices())) return(tib)
      tib[tib$source %in% keep & tib$target %in% keep, ]
    }

    active.table <- reactive({
      req(input$table.group)
      table.for.group(input$table.group)
    })

    selected.rows <- reactive({
      tib <- active.table()
      idx <- input$ccc.table_rows_selected
      if (is.character(tib) || is.null(idx) || length(idx) == 0) return(NULL)
      # Indices can outlive the table they were picked from for one reactive beat
      # after the cell-type filter changes and before DT reports the cleared selection.
      idx <- idx[idx <= nrow(tib)]
      if (length(idx) == 0) return(NULL)
      tib[idx, ]
    })

    # ---------- Score column selectors ----------

    score.choices <- reactive({
      tib <- active.table()
      if (is.character(tib)) return(character(0))
      .liana.score.cols(tib)
    })

    # Only push new choices when the column set actually changes, so that toggling a
    # cell type (which re-runs active.table) does not reset the user's metric picks.
    applied.score.cols <- reactiveVal(NULL)
    observe({
      cols <- score.choices()
      req(length(cols) > 0)
      if (identical(cols, applied.score.cols())) return()
      applied.score.cols(cols)

      size.default <- .liana.pick.col(cols, c("natmi.edge_specificity",
                                              "edge_specificity",
                                              "connectome.weight_sc", "weight_sc"), 1)
      colour.default <- .liana.pick.col(cols, c("sca.LRscore", "LRscore",
                                                "natmi.expr_prod", "expr_prod",
                                                "cellphonedb.lr.mean", "lr.mean"), 2)
      rank.default <- .liana.pick.col(cols, c("aggregate_rank", size.default), 1)

      updateSelectInput(session, "dot.size", choices = cols, selected = size.default)
      updateSelectInput(session, "dot.col", choices = cols, selected = colour.default)
      updateSelectInput(session, "dot.rank", choices = cols, selected = rank.default)
      updateSelectInput(session, "heat.cutoff.col", choices = cols,
                        selected = rank.default)
      updateSelectInput(session, "chord.cutoff.col", choices = cols,
                        selected = rank.default)
    })

    # Score columns span orders of magnitude between methods, so there is no sane
    # fixed default: seed the value from the column's own distribution and let the
    # user override it.
    seed.cutoff <- function(prefix) {
      observeEvent(input[[paste0(prefix, ".cutoff.col")]], {
        col <- input[[paste0(prefix, ".cutoff.col")]]
        tib <- active.table()
        req(!is.character(tib), col %in% names(tib))
        x <- tib[[col]]
        x <- x[is.finite(x)]
        req(length(x) > 0)
        if (.liana.is.pval.col(col)) {
          op <- "le"
          value <- if (col == "aggregate_rank") 0.01
                   else signif(stats::quantile(x, 0.1, names = FALSE), 3)
        } else {
          op <- "ge"
          value <- signif(stats::quantile(x, 0.9, names = FALSE), 3)
        }
        updateSelectInput(session, paste0(prefix, ".cutoff.op"), selected = op)
        updateNumericInput(session, paste0(prefix, ".cutoff.val"), value = value)
      })
    }
    seed.cutoff("heat")
    seed.cutoff("chord")

    # ---------- CCC Table ----------

    output$table.caption <- renderPrint({
      mode <- method.mode()
      validate(need(mode$ok, mode$msg))
      tib <- active.table()
      validate(need(!is.character(tib), tib))
      cat("Group:    ", input$table.group, "\n")
      cat("Table:    ",
          if (mode$aggregate) paste0("liana_aggregate of ",
                                     paste(sort(applied.methods()), collapse = ", "))
          else paste0(applied.methods(), " (raw)"), "\n")
      cat("Rows:     ", nrow(tib), " (after the cell-type filter)\n", sep = "")
      cat("Selected: ", length(input$ccc.table_rows_selected %||% integer(0)), "\n")
    })

    output$ccc.table <- DT::renderDT({
      mode <- method.mode()
      validate(need(mode$ok, mode$msg))
      tib <- active.table()
      validate(need(!is.character(tib), tib))
      validate(need(nrow(tib) > 0, "No rows left after the cell-type filter."))

      df <- as.data.frame(tib)
      numeric.cols <- names(df)[vapply(df, is.numeric, logical(1))]
      dt <- DT::datatable(
        df, rownames = FALSE, filter = "top", selection = "multiple",
        options = list(pageLength = 25, scrollX = TRUE, searchDelay = 500)
      )
      if (length(numeric.cols) > 0) {
        dt <- DT::formatSignif(dt, columns = numeric.cols, digits = 3)
      }
      dt
    })

    output$table.download <- downloadHandler(
      filename = function() {
        mode <- method.mode()
        tag <- if (isTRUE(mode$aggregate)) "aggregated"
               else paste(applied.methods(), collapse = "_")
        paste0("liana_table_",
               gsub("[^A-Za-z0-9]", "_", input$table.group %||% ""), "_",
               tag, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        tib <- active.table()
        validate(need(!is.character(tib), tib))
        write.csv(as.data.frame(tib), file, row.names = FALSE)
      }
    )

    # ---------- Plot outputs ----------

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

    output$dot.plot.ui   <- plot.output("dot")
    output$heat.plot.ui  <- plot.output("heat")
    output$chord.plot.ui <- plot.output("chord")

    dot.thunk <- function() {
      mode <- method.mode()
      validate(need(mode$ok, mode$msg))
      g <- input$table.group
      req(g)
      tib <- active.table()
      validate(need(!is.character(tib), tib))
      validate(need(nrow(tib) > 0, "No rows left after the cell-type filter."))
      validate(
        need(input$dot.size %in% names(tib), "Choose a dot size column."),
        need(input$dot.col %in% names(tib), "Choose a dot colour column.")
      )

      if (identical(input$dot.rows, "selected")) {
        tib.plot <- selected.rows()
        validate(need(!is.null(tib.plot) && nrow(tib.plot) > 0,
                      "No rows are selected in the CCC Table."))
        source.use <- NULL
        target.use <- NULL
      } else {
        validate(need(input$dot.rank %in% names(tib), "Choose a ranking column."))
        tib.idents <- unique(c(tib$source, tib$target))
        all.idents <- subset.universe()
        s <- resolve.sel(tib.idents, input$dot.sources, all.idents)
        t <- resolve.sel(tib.idents, input$dot.targets, all.idents)
        validate(
          need(!s$empty, paste0("Selected source(s) not present in '", g, "'.")),
          need(!t$empty, paste0("Selected target(s) not present in '", g, "'."))
        )
        source.use <- if (s$active) s$use else NULL
        target.use <- if (t$active) t$use else NULL
        # Rank within the selected populations, so "top 15" means the top pairs of
        # the interactions actually drawn rather than of the whole group.
        if (!is.null(source.use)) tib <- tib[tib$source %in% source.use, ]
        if (!is.null(target.use)) tib <- tib[tib$target %in% target.use, ]
        validate(need(nrow(tib) > 0,
                      "No interactions between the selected populations."))
        tib.plot <- .liana.top.pairs(tib, input$dot.rank, input$dot.top)
        validate(need(nrow(tib.plot) > 0,
                      paste0("No rows have a value for '", input$dot.rank, "'.")))
      }

      # call_cellchat's ligand/receptor already hold the complex identifier: CellChat
      # never decomplexifies, so unlike liana's own pipeline — where ligand is the
      # min-expressed subunit and ligand.complex the complex — there is no second
      # column pair for liana_dotplot to unite on.
      show.complex <- all(c("ligand.complex", "receptor.complex") %in% names(tib.plot))

      size.is.pval <- .liana.is.pval.col(input$dot.size)
      col.is.pval  <- .liana.is.pval.col(input$dot.col)
      size.label <- if (size.is.pval) paste0("-log10(", input$dot.size, ")")
                    else input$dot.size
      colour.label <- if (col.is.pval) paste0("-log10(", input$dot.col, ")")
                      else input$dot.col

      p <- tryCatch(
        liana::liana_dotplot(
          tib.plot,
          source_groups = source.use,
          target_groups = target.use,
          ntop = NULL,
          specificity = input$dot.size,
          magnitude = input$dot.col,
          y.label = "Interactions (Ligand -> Receptor)",
          size.label = size.label,
          colour.label = colour.label,
          show_complex = show.complex,
          size_range = c(2, 10),
          invert_specificity = size.is.pval,
          invert_magnitude = col.is.pval,
          invert_function = function(x) -log10(x + 1e-10)
        ) +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,
                                           color = "black", size = 13)) +
          ggtitle(g),
        error = function(e) {
          ggplot() + theme_void() +
            ggtitle(paste0(g, "\n[liana_dotplot: ", conditionMessage(e), "]"))
        }
      )
      print(p)
    }

    # Both grids apply the same method selection to every group, so the panels stay
    # comparable; only the cutoff and the cell-type filter vary with the controls.
    heat.thunk <- function() {
      mode <- method.mode()
      validate(need(mode$ok, mode$msg))
      res <- liana.data(); req(res)
      grps <- res$group.levels
      ncol <- input$heat.cols %||% 1

      items <- lapply(grps, function(g) {
        force(g)
        tib <- table.for.group(g)
        if (is.character(tib)) return(function() placeholder.base(g, tib))
        if (nrow(tib) == 0) {
          return(function() placeholder.base(g, "no rows after the cell-type filter"))
        }
        filtered <- .liana.apply.cutoff(tib, input$heat.cutoff.col,
                                        input$heat.cutoff.op, input$heat.cutoff.val)
        if (is.character(filtered)) {
          return(function() placeholder.base(g, filtered))
        }
        if (nrow(filtered) == 0) {
          return(function() placeholder.base(g, "no rows passed the cutoff"))
        }
        function() {
          ht <- tryCatch(liana::heat_freq(filtered), error = function(e) NULL)
          if (is.null(ht)) {
            placeholder.base(g, "heat_freq failed")
            return(invisible())
          }
          ComplexHeatmap::draw(ht)
        }
      })
      plot.grid(items, ncol, titles = grps)
    }

    chord.thunk <- function() {
      mode <- method.mode()
      validate(need(mode$ok, mode$msg))
      res <- liana.data(); req(res)
      grps <- res$group.levels
      ncol <- input$chord.cols %||% 1
      all.idents <- subset.universe()
      srcs <- input$chord.sources
      tgts <- input$chord.targets

      items <- lapply(grps, function(g) {
        force(g)
        tib <- table.for.group(g)
        if (is.character(tib)) return(function() placeholder.base(g, tib))
        if (nrow(tib) == 0) {
          return(function() placeholder.base(g, "no rows after the cell-type filter"))
        }
        s <- resolve.sel(unique(tib$source), srcs, all.idents)
        t <- resolve.sel(unique(tib$target), tgts, all.idents)
        filtered <- .liana.apply.cutoff(tib, input$chord.cutoff.col,
                                        input$chord.cutoff.op, input$chord.cutoff.val)

        function() {
          if (s$empty || t$empty) {
            placeholder.base(g, "selection absent in this group")
            return(invisible())
          }
          if (is.character(filtered)) {
            placeholder.base(g, filtered)
            return(invisible())
          }
          if (nrow(filtered) == 0) {
            placeholder.base(g, "no rows passed the cutoff")
            return(invisible())
          }
          tryCatch(
            liana::chord_freq(filtered,
                              source_groups = if (s$active) s$use else NULL,
                              target_groups = if (t$active) t$use else NULL),
            error = function(e) placeholder.base(g, conditionMessage(e))
          )
        }
      })
      plot.grid(items, ncol, titles = grps)
    }

    make.render.and.download(output, session, "dot.plot", "dot.download",
      dot.thunk,
      width.in = reactive(input$dot.width),
      height.in = reactive(input$dot.height),
      filename.stem = reactive(paste0("liana_dot_",
                                      gsub("[^A-Za-z0-9]", "_",
                                           input$table.group %||% ""))),
      trigger = reactive(input$dot.plot.btn))

    make.render.and.download(output, session, "heat.plot", "heat.download",
      heat.thunk,
      width.in = reactive(input$heat.width),
      height.in = reactive(input$heat.height),
      filename.stem = reactive("liana_freq_heatmap"),
      trigger = reactive(input$heat.plot.btn))

    make.render.and.download(output, session, "chord.plot", "chord.download",
      chord.thunk,
      width.in = reactive(input$chord.width),
      height.in = reactive(input$chord.height),
      filename.stem = reactive("liana_freq_chord"),
      trigger = reactive(input$chord.plot.btn))
  })
}
