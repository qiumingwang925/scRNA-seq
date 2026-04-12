## ABOUTME: Container module that wraps SingleR and manual annotation as sub-tabs.
## ABOUTME: Manages a shared Seurat object that both annotation methods can read and update.

mod.annotation.ui <- function(id) {
  ns <- NS(id)

  tabPanel("Annotation", value = "tab.annotation",
    tabsetPanel(id = ns("annotation.tabs"),
      tabPanel("SingleR Annotation",
        mod.annotation.singler.ui(ns("singler"))
      ),
      tabPanel("Manual Annotation",
        mod.annotation.manual.ui(ns("manual"))
      )
    )
  )
}

mod.annotation.server <- function(id, seurat.obj.cellcycle) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Shared Seurat object, initialized from cell cycle output
    current.obj <- reactiveVal(NULL)

    # Track the upstream object so we can detect re-runs
    prev.upstream <- reactiveVal(NULL)

    # Initialize from upstream pipeline, warn on re-run if annotations exist
    observe({
      req(seurat.obj.cellcycle())
      upstream <- seurat.obj.cellcycle()

      # isolate() prevents this observe from re-triggering when we set
      # current.obj or prev.upstream below — it should only fire when
      # seurat.obj.cellcycle() changes.
      if (is.null(isolate(current.obj()))) {
        # First time: initialize directly
        if (!"manual_annotation" %in% colnames(upstream@meta.data)) {
          upstream$manual_annotation <- "Unlabeled"
        }
        current.obj(upstream)
        prev.upstream(upstream)
      } else if (!identical(upstream, isolate(prev.upstream()))) {
        # Upstream changed: confirm before overwriting annotations
        showModal(modalDialog(
          title = "Upstream Data Changed",
          "An earlier pipeline step was re-run. Accepting the new data will discard any annotations you've made in this tab.",
          footer = tagList(
            actionButton(ns("accept.upstream"), "Accept New Data", class = "btn-danger"),
            modalButton("Keep Current Annotations")
          )
        ))
        prev.upstream(upstream)
      }
    })

    # Apply new upstream data when user confirms
    observeEvent(input$accept.upstream, {
      srt <- seurat.obj.cellcycle()
      if (!"manual_annotation" %in% colnames(srt@meta.data)) {
        srt$manual_annotation <- "Unlabeled"
      }
      current.obj(srt)
      completed(FALSE)
      removeModal()
    })

    # Wire sub-modules with shared object
    singler.result <- mod.annotation.singler.server("singler", current.obj)
    manual.result <- mod.annotation.manual.server("manual", current.obj)

    # Completed when either sub-module completes
    completed <- reactiveVal(FALSE)
    observe({
      if (singler.result$completed() || manual.result$completed()) {
        completed(TRUE)
      }
    })

    return(list(
      seurat.obj = reactive({ current.obj() }),
      completed = completed
    ))
  })
}
