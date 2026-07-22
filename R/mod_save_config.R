## ABOUTME: Reusable Shiny sub-module for configurable Seurat object export.
## ABOUTME: Object-driven component picker (assays/layers/reductions/graphs); metadata is always kept.

# ── Export presets ───────────────────────────────────────────────────────────
# A preset pre-checks a subset of components so repeated or testing exports are
# lighter and faster to write. Each entry is a name -> function(obj) that
# returns which components to check; any name not present in the object is
# ignored, so a preset is safe to reuse across objects of different shapes.
# Metadata is always kept and is not part of the spec. "Full object"
# (everything checked) is always offered and needs no entry here.
#
# To define one, uncomment and edit — it then appears as a preset button:
#
# EXPORT.PRESETS[["Counts only (fast)"]] <- function(obj) list(
#   assays       = "RNA",            # assays to keep
#   layers       = "counts",         # layers to keep, applied across kept assays
#   reductions   = character(0),     # e.g. Reductions(obj) to keep all
#   graphs       = character(0),
#   active.ident = "manual_annotation"  # optional: metadata column to set as Idents
# )
EXPORT.PRESETS <- list()

mod.save.config.ui <- function(id, label = "Export Object") {
  ns <- NS(id)
  actionButton(ns("open"), label, style = "width:100%")
}

mod.save.config.server <- function(id, current.obj) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Breadcrumb logging for the export path: the last line printed before a
    # silent failure localizes the crash. Prefixed with the module namespace
    # so the three call sites (SingleR, Manual, Module II) are distinguishable.
    log.step <- function(...) message("[save ", session$ns(""), "] ", ...)

    # Assays() is qualified throughout: the apps load SummarizedExperiment
    # (via SingleCellExperiment in Module I, via batchelor in Module II) after
    # Seurat, and its Assays() masks SeuratObject's, returning an S4 assays
    # container instead of the character names the checkbox choices need.

    # Layers actually present across all assays, in canonical order.
    available.layers <- function(obj) {
      present <- unique(unlist(lapply(SeuratObject::Assays(obj), function(a) Layers(obj[[a]]))))
      intersect(c("counts", "data", "scale.data"), present)
    }

    # Only render a checkbox group when the object actually has that component.
    optional.group <- function(input.id, label, choices) {
      if (length(choices) == 0) return(NULL)
      checkboxGroupInput(ns(input.id), label, choices = choices, selected = choices)
    }

    # Metadata columns that make sense as a cell identity, intersected with what
    # the object actually has. SingleR* is a prefix match (SingleR.main, etc.);
    # the rest are exact. Kept in most- to least-specific order.
    identity.candidates <- function(obj) {
      cols <- colnames(obj@meta.data)
      singler <- grep("^SingleR", cols, value = TRUE)
      c(intersect("manual_annotation", cols), singler,
        intersect(c("condition", "group", "orig.ident"), cols))
    }

    observeEvent(input$open, {
      req(current.obj())
      obj <- current.obj()
      log.step("open dialog | size=", format(object.size(obj), units = "auto"),
               " | assays=", paste(SeuratObject::Assays(obj), collapse = ","),
               " | layers=", paste(available.layers(obj), collapse = ","),
               " | reductions=", paste(Reductions(obj), collapse = ","),
               " | graphs=", paste(Graphs(obj), collapse = ","))

      showModal(modalDialog(
        title = "Export Configuration",
        size = "m",
        p(strong("Current object size: "), format(object.size(obj), units = "auto")),
        radioButtons(ns("preset"), "Preset:",
                     choices = c("Full object", names(EXPORT.PRESETS)),
                     selected = "Full object"),
        hr(),
        optional.group("keep.assays", "Assays:", SeuratObject::Assays(obj)),
        optional.group("keep.layers", "Layers (applied across kept assays):",
                       available.layers(obj)),
        optional.group("keep.reductions", "Reductions:", Reductions(obj)),
        optional.group("keep.graphs", "Graphs:", Graphs(obj)),
        tags$label("Metadata columns:"),
        div(style = "max-height: 180px; overflow-y: auto; border: 1px solid #ddd; padding: 6px; border-radius: 4px;",
            optional.group("keep.meta", NULL, colnames(obj@meta.data))),
        if (length(identity.candidates(obj))) {
          selectInput(ns("active.ident.col"), "Set active identity to:",
                      choices = c("(leave unchanged)" = "", identity.candidates(obj)),
                      selected = "")
        },
        footer = tagList(
          downloadButton(ns("download"), "Download"),
          modalButton("Cancel")
        )
      ))
    })

    # Applying a preset just re-checks the boxes; the user can still adjust after.
    observeEvent(input$preset, {
      req(current.obj())
      obj <- current.obj()
      log.step("preset selected: ", input$preset)

      spec <- if (input$preset == "Full object") {
        list(assays = SeuratObject::Assays(obj), layers = available.layers(obj),
             reductions = Reductions(obj), graphs = Graphs(obj),
             metadata = colnames(obj@meta.data))
      } else {
        EXPORT.PRESETS[[input$preset]](obj)
      }
      pick <- function(key, available) {
        wanted <- spec[[key]]
        if (is.null(wanted)) wanted <- character(0)
        intersect(wanted, available)
      }

      updateCheckboxGroupInput(session, "keep.assays",     selected = pick("assays", SeuratObject::Assays(obj)))
      updateCheckboxGroupInput(session, "keep.layers",     selected = pick("layers", available.layers(obj)))
      updateCheckboxGroupInput(session, "keep.reductions", selected = pick("reductions", Reductions(obj)))
      updateCheckboxGroupInput(session, "keep.graphs",     selected = pick("graphs", Graphs(obj)))
      updateCheckboxGroupInput(session, "keep.meta",       selected = pick("metadata", colnames(obj@meta.data)))

      # A preset may name an active.ident column; otherwise leave the choice be.
      ident.wanted <- spec$active.ident
      if (!is.null(ident.wanted) && ident.wanted %in% identity.candidates(obj)) {
        updateSelectInput(session, "active.ident.col", selected = ident.wanted)
      }
    }, ignoreInit = TRUE)

    strip.object <- function(obj) {
      assays.keep <- input$keep.assays
      layers.keep <- input$keep.layers
      log.step("strip: assays=", paste(assays.keep, collapse = ","),
               " | layers=", paste(layers.keep, collapse = ","),
               " | reductions=", paste(input$keep.reductions, collapse = ","),
               " | graphs=", paste(input$keep.graphs, collapse = ","))
      if (length(assays.keep) == 0) stop("Select at least one assay to export.")
      if (length(layers.keep) == 0) stop("Select at least one layer to export.")

      # DietSeurat requires the active assay to survive; repoint it if the user
      # unchecked it, otherwise the call errors.
      if (!DefaultAssay(obj) %in% assays.keep) {
        log.step("repointing DefaultAssay ", DefaultAssay(obj), " -> ", assays.keep[1])
        DefaultAssay(obj) <- assays.keep[1]
      }

      dimreducs <- if (length(input$keep.reductions)) input$keep.reductions else NULL
      graphs    <- if (length(input$keep.graphs)) input$keep.graphs else NULL

      log.step("calling DietSeurat")
      out <- if (packageVersion("Seurat") >= "5.0.0") {
        DietSeurat(obj, layers = layers.keep, assays = assays.keep,
                   dimreducs = dimreducs, graphs = graphs)
      } else {
        DietSeurat(obj,
                   counts = "counts" %in% layers.keep,
                   data = "data" %in% layers.keep,
                   scale.data = "scale.data" %in% layers.keep,
                   assays = assays.keep, dimreducs = dimreducs, graphs = graphs)
      }

      # Set the active identity before dropping columns: Idents() copies the
      # values into the @active.ident slot, which is independent of meta.data,
      # so the identity survives even if that column is unchecked below.
      ident.col <- input$active.ident.col
      if (!is.null(ident.col) && nzchar(ident.col) && ident.col %in% colnames(out@meta.data)) {
        log.step("setting active.ident from column: ", ident.col)
        Idents(out) <- ident.col
      }

      # DietSeurat keeps every metadata column; drop the ones the user unchecked.
      meta.keep <- intersect(input$keep.meta, colnames(out@meta.data))
      dropped <- setdiff(colnames(out@meta.data), meta.keep)
      if (length(dropped)) {
        log.step("dropping metadata columns: ", paste(dropped, collapse = ","))
        out@meta.data <- out@meta.data[, meta.keep, drop = FALSE]
      }

      log.step("DietSeurat done | size=", format(object.size(out), units = "auto"))
      out
    }

    output$download <- downloadHandler(
      filename = function() {
        log.step("download requested")
        paste0("seurat_export_", Sys.Date(), ".rds")
      },
      content = function(file) {
        log.step("content start | writing to ", file)
        withProgress(message = "Preparing export...", value = 0.5, {
          stripped <- tryCatch(strip.object(current.obj()),
                               error = function(e) {
                                 log.step("ERROR in strip.object: ", conditionMessage(e))
                                 showNotification(conditionMessage(e), type = "error")
                                 NULL
                               })
          req(stripped)
          log.step("saveRDS start")
          tryCatch(
            saveRDS(stripped, file),
            error = function(e) {
              log.step("ERROR in saveRDS: ", conditionMessage(e))
              stop(e)
            }
          )
          log.step("saveRDS done")
        })
        removeModal()
        log.step("export complete")
      }
    )
  })
}
