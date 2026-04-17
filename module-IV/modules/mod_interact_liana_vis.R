## ABOUTME: LIANA visualization module (Phase 3 stub — subtabs added in Phase 5/6).
## ABOUTME: Resolves input from upstream comp or an uploaded .rds and reports a summary.

mod.interact.liana.vis.ui <- function(id) {
  ns <- NS(id)
  tabPanel("Visualize LIANA",
    fluidRow(
      column(12,
        wellPanel(
          fileInput(ns("vis.file"),
                    "Upload processed LIANA result (.rds) — skip if computed above:",
                    accept = ".rds"),
          verbatimTextOutput(ns("input.summary"))
        )
      )
    ),
    helpText("Visualization subtabs (CCC Dot Plot, Freq Heatmap, Freq Chord Diagram)",
             "will be added in Phase 5 and 6.")
  )
}

mod.interact.liana.vis.server <- function(id, liana.input) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

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
      cat("Methods:  ", paste(res$methods, collapse = ", "), "\n")
      cat("Rows per group:\n")
      for (g in res$group.levels) {
        cat("  ", g, ": ", nrow(res$liana.list[[g]]), "\n", sep = "")
      }
    })
  })
}

# Null-coalesce helper used by the input summary
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}
