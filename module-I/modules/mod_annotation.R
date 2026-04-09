## ABOUTME: Container module that wraps SingleR and manual annotation as sub-tabs.
## ABOUTME: Manages a shared Seurat object that both annotation methods can read and update.

mod.annotation.ui <- function(id) {
  ns <- NS(id)

  tabPanel("Annotation",
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

    # Initialize from upstream pipeline
    observe({
      req(seurat.obj.cellcycle())
      if (is.null(current.obj())) {
        srt <- seurat.obj.cellcycle()
        if (!"manual_annotation" %in% colnames(srt@meta.data)) {
          srt$manual_annotation <- "Unlabeled"
        }
        current.obj(srt)
      }
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
