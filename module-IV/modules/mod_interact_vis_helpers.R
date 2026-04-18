## ABOUTME: Shared UI + render helpers used by both CellChat and LIANA vis modules.
## ABOUTME: Source this file before either vis module.

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Treat length-0 or NULL selection as "no filter"; return NULL so callers can
# forward it straight into a function argument that defaults to "all".
nz <- function(x) if (!is.null(x) && length(x) > 0) x else NULL

# Cell-type selector matching Module III's pattern: selectInput(multiple) +
# Select All / Clear buttons.
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

wire.cell.type.selector <- function(input, session, id,
                                    choices.reactive, selected.all = TRUE) {
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

# Common sidebar controls: Generate Plot button, cols / width / height, Download
common.controls.ui <- function(ns, prefix, default.cols = 2,
                               dl.label = "Download Figure") {
  tagList(
    hr(),
    actionButton(ns(paste0(prefix, ".plot.btn")), "Generate Plot",
                 class = "btn-success", style = "width:100%"),
    hr(),
    h5("Figure controls"),
    numericInput(ns(paste0(prefix, ".cols")), "Columns per row:",
                 value = default.cols, min = 1, max = 6, step = 1),
    numericInput(ns(paste0(prefix, ".width")), "Width (inches):",
                 value = 10, min = 3, max = 30, step = 1),
    numericInput(ns(paste0(prefix, ".height")), "Height (inches):",
                 value = 6, min = 3, max = 30, step = 1),
    downloadButton(ns(paste0(prefix, ".download")), dl.label, style = "width:100%")
  )
}

# Bind a plot output to a Generate Plot button and attach a matching download
# handler. thunk produces the plot via side-effects (base graphics or print()).
make.render.and.download <- function(output, session, plot.id, dl.id, thunk,
                                     width.in, height.in, filename.stem, trigger) {
  output[[plot.id]] <- shiny::bindEvent(
    renderPlot({ thunk() }, res = 96),
    trigger(),
    ignoreInit = TRUE, ignoreNULL = TRUE
  )
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

# Arrange heterogeneous panel items in a patchwork grid. Each item may be:
#   - a ggplot / patchwork                  -> used directly (fast)
#   - a ComplexHeatmap Heatmap / HeatmapList -> captured via grid.grabExpr (fast)
#   - a function                              -> replayed via ggplotify::as.ggplot
#                                                (slow but required for base graphics
#                                                whose internal par(mfrow) otherwise
#                                                clashes with an outer par grid)
#   - a grob / gTree                          -> wrapped directly
# Function items get a bold centered title from `titles[i]`; other item types
# keep their own titles (CellChat's title.name, ggtitle(), column_title, ...).
plot.grid <- function(items, ncol, titles = NULL) {
  plots <- lapply(seq_along(items), function(i) {
    item <- items[[i]]
    title <- if (!is.null(titles) && length(titles) >= i) titles[i] else NULL
    tryCatch({
      if (inherits(item, c("ggplot", "patchwork", "gg"))) {
        item
      } else if (inherits(item, c("Heatmap", "HeatmapList"))) {
        patchwork::wrap_elements(
          full = grid::grid.grabExpr(ComplexHeatmap::draw(item)))
      } else if (is.function(item)) {
        p <- ggplotify::as.ggplot(item)
        if (!is.null(title)) p <- p + ggtitle(title) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        p
      } else if (inherits(item, c("grob", "gTree"))) {
        patchwork::wrap_elements(full = item)
      } else {
        ggplot() + theme_void() +
          ggtitle(paste("Unhandled plot type:",
                        paste(class(item), collapse = "/")))
      }
    }, error = function(e) {
      label <- title %||% paste("Panel", i)
      ggplot() + theme_void() +
        ggtitle(paste0(label, "\n[", conditionMessage(e), "]"))
    })
  })
  print(patchwork::wrap_plots(plots, ncol = ncol))
}

# Three-state resolution of a selection input against a context-specific
# set of available cell types. Shared across CellChat and LIANA vis.
#   available  : character vector of idents actually present in this context
#                (e.g. levels(cc@idents), or unique(c(tib$source, tib$target)))
#   sel        : the user's selection (may be NULL, character(0), or names)
#   all.idents : the full universe across all groups — used to decide "user
#                has not filtered" when sel equals the universe
# Returns list(active, use, empty):
#   active = FALSE -> user hasn't filtered; callers should ignore `use`
#   empty  = TRUE  -> user filtered but none of their picks exist here;
#                     callers should render an "absent" placeholder
resolve.sel <- function(available, sel, all.idents) {
  if (is.null(sel) || length(sel) == 0 || setequal(sel, all.idents)) {
    return(list(active = FALSE, use = NULL, empty = FALSE))
  }
  kept <- intersect(sel, available)
  list(active = TRUE, use = kept, empty = length(kept) == 0)
}

# Typed "absent selection / absent pathway" placeholders for the three
# render families: base graphics, ggplot, ComplexHeatmap.
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
