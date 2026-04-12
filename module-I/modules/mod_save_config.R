## ABOUTME: Reusable Shiny sub-module for configurable Seurat object export.
## ABOUTME: Offers presets (full, analysis-ready, CSV) and custom component selection via DietSeurat.

mod.save.config.ui <- function(id, label = "Export Object") {
  ns <- NS(id)
  actionButton(ns("open"), label, style = "width:100%")
}

mod.save.config.server <- function(id, current.obj) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Open modal with export options
    observeEvent(input$open, {
      req(current.obj())

      obj.size <- format(object.size(current.obj()), units = "auto")

      showModal(modalDialog(
        title = "Export Configuration",
        size = "m",

        p(strong("Current object size: "), obj.size),

        radioButtons(ns("preset"), "Export preset:",
          choices = c(
            "Full object" = "full",
            "Analysis-ready (recommended)" = "analysis",
            "Custom" = "custom",
            "Annotations only (CSV)" = "csv"
          ),
          selected = "analysis"
        ),

        textOutput(ns("preset.description")),
        br(),
        uiOutput(ns("advanced.ui")),

        footer = tagList(
          downloadButton(ns("download"), "Download"),
          modalButton("Cancel")
        )
      ))
    })

    # Preset descriptions
    output$preset.description <- renderText({
      req(input$preset)
      switch(input$preset,
        "full" = "Keeps everything. Largest file size.",
        "analysis" = "Drops scaled data and neighbor graphs. Keeps counts, normalized data, PCA, UMAP, and all metadata.",
        "custom" = "Choose which components to include.",
        "csv" = "Exports metadata and UMAP coordinates as a CSV file. Smallest output."
      )
    })

    # Advanced checkboxes (only shown for Custom preset)
    output$advanced.ui <- renderUI({
      if (is.null(input$preset) || input$preset != "custom") return(NULL)

      obj <- current.obj()
      has.sct <- "SCT" %in% names(obj@assays)

      choices <- c(
        "Raw counts" = "counts",
        "Normalized data" = "data",
        "Scaled data" = "scale.data",
        "PCA embeddings" = "pca",
        "UMAP embeddings" = "umap",
        "Neighbor graphs" = "graphs"
      )
      if (has.sct) choices <- c(choices, c("SCT assay" = "sct"))

      wellPanel(
        checkboxGroupInput(ns("keep"), "Components to keep:",
          choices = choices,
          selected = c("counts", "data", "pca", "umap")
        ),
        helpText("Metadata (annotations, clusters) is always included.")
      )
    })

    # Strip a Seurat object based on selected components
    strip.object <- function(obj, keep) {
      keep.layers <- intersect(keep, c("counts", "data", "scale.data"))
      dimreducs <- c()
      if ("pca" %in% keep) dimreducs <- c(dimreducs, "pca")
      if ("umap" %in% keep) dimreducs <- c(dimreducs, "umap")
      if (length(dimreducs) == 0) dimreducs <- NULL

      graphs.keep <- if ("graphs" %in% keep) names(obj@graphs) else NULL

      assays <- DefaultAssay(obj)
      if ("sct" %in% keep && "SCT" %in% names(obj@assays)) {
        assays <- unique(c(assays, "SCT"))
      }

      if (packageVersion("Seurat") >= "5.0.0") {
        DietSeurat(obj, layers = keep.layers, dimreducs = dimreducs,
                   graphs = graphs.keep, assays = assays)
      } else {
        DietSeurat(obj,
                   counts = "counts" %in% keep.layers,
                   data = "data" %in% keep.layers,
                   scale.data = "scale.data" %in% keep.layers,
                   dimreducs = dimreducs,
                   graphs = graphs.keep,
                   assays = assays)
      }
    }

    # Download handler — stripping only runs on actual download
    output$download <- downloadHandler(
      filename = function() {
        if (input$preset == "csv") {
          paste0("annotations_", Sys.Date(), ".csv")
        } else {
          paste0("seurat_", input$preset, "_", Sys.Date(), ".rds")
        }
      },
      content = function(file) {
        withProgress(message = "Preparing export...", value = 0.5, {
          obj <- current.obj()

          if (input$preset == "csv") {
            umap.coords <- as.data.frame(Embeddings(obj, "umap"))
            meta <- obj@meta.data
            export.df <- cbind(barcode = rownames(meta), meta, umap.coords)
            write.csv(export.df, file, row.names = FALSE)
          } else if (input$preset == "full") {
            saveRDS(obj, file)
          } else {
            keep <- if (input$preset == "analysis") {
              c("counts", "data", "pca", "umap")
            } else {
              input$keep
            }
            stripped <- strip.object(obj, keep)
            saveRDS(stripped, file)
          }
        })
        removeModal()
      }
    )
  })
}
