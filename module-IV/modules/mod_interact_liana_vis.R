## ABOUTME: LIANA visualization module. Three subtabs — Dot Plot, Freq Heatmap,
## ABOUTME: Freq Chord Diagram. Dot Plot is single-group; the other two are N-group grids.

# LIANA method <-> score column map for the Dot Plot size / colour selectors.
# Keys are human-readable labels; values describe how to pull the value out of
# the aggregated LIANA tibble.
.liana.score.map <- list(
  "Interaction Specificity (NATMI)"  = list(col = "natmi.edge_specificity",
                                             label = "Interaction\nSpecificity"),
  "Interaction Weight (Connectome)"  = list(col = "connectome.weight_sc",
                                             label = "Interaction\nWeight"),
  "LogFC Mean (iTALK-style)"         = list(col = "logfc.logfc_comb",
                                             label = "LogFC\nMean"),
  "Expression Magnitude (SCSignalR)" = list(col = "sca.LRscore",
                                             label = "Expression\nMagnitude"),
  "P-value (CellPhoneDB)"            = list(col = "cellphonedb.pvalue",
                                             label = "-log10(p)")
)
.liana.score.labels <- names(.liana.score.map)

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

    tabsetPanel(id = ns("liana.vis.tabs"),

      # ---------- Subtab 1: CCC Dot Plot (single-group) ----------
      tabPanel("CCC Dot Plot",
        sidebarLayout(
          sidebarPanel(width = 3,
            h4("CCC Dot Plot"),
            selectInput(ns("dot.group"), "Group:", choices = NULL),
            cell.type.selector.ui(ns, "dot.sources", "Source populations:"),
            cell.type.selector.ui(ns, "dot.targets", "Target populations:"),
            numericInput(ns("dot.top"), "Top N L-R pairs:",
                         value = 15, min = 1, max = 200, step = 1),
            selectInput(ns("dot.size"), "Dot size:",
                        choices = .liana.score.labels,
                        selected = "Interaction Specificity (NATMI)"),
            selectInput(ns("dot.col"), "Dot colour:",
                        choices = .liana.score.labels,
                        selected = "Expression Magnitude (SCSignalR)"),
            common.controls.ui(ns, "dot", default.cols = 1)
          ),
          mainPanel(width = 9,
            uiOutput(ns("dot.plot.ui"))
          )
        )
      ),

      # ---------- Subtab 2: CCC Freq Heatmap (N-group grid) ----------
      tabPanel("CCC Freq Heatmap",
        sidebarLayout(
          sidebarPanel(width = 3,
            h4("CCC Freq Heatmap"),
            cell.type.selector.ui(ns, "heat.idents", "Populations (source & target):"),
            numericInput(ns("heat.alpha"), "aggregate_rank cutoff:",
                         value = 0.01, min = 0, max = 1, step = 0.005),
            common.controls.ui(ns, "heat", default.cols = 2)
          ),
          mainPanel(width = 9,
            uiOutput(ns("heat.plot.ui"))
          )
        )
      ),

      # ---------- Subtab 3: CCC Freq Chord Diagram (N-group grid) ----------
      tabPanel("CCC Freq Chord Diagram",
        sidebarLayout(
          sidebarPanel(width = 3,
            h4("CCC Freq Chord Diagram"),
            cell.type.selector.ui(ns, "chord.sources", "Source populations:"),
            cell.type.selector.ui(ns, "chord.targets", "Target populations:"),
            numericInput(ns("chord.alpha"), "aggregate_rank cutoff:",
                         value = 0.01, min = 0, max = 1, step = 0.005),
            common.controls.ui(ns, "chord", default.cols = 2)
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

    # Shared choice reactives — union of source + target cell types across groups
    ident.choices <- reactive({
      res <- liana.data(); req(res)
      all <- Reduce(union, lapply(res$liana.list, function(tib) {
        unique(c(tib$source, tib$target))
      }))
      sort(all[!is.na(all) & nzchar(all)])
    })

    group.choices <- reactive({
      res <- liana.data(); req(res)
      res$group.levels
    })

    # Populate dot-plot group dropdown
    observe({
      grps <- group.choices()
      req(grps)
      updateSelectInput(session, "dot.group",
                        choices = grps,
                        selected = grps[1])
    })

    wire.cell.type.selector(input, session, "dot.sources",  ident.choices, selected.all = TRUE)
    wire.cell.type.selector(input, session, "dot.targets",  ident.choices, selected.all = TRUE)
    wire.cell.type.selector(input, session, "heat.idents",  ident.choices, selected.all = TRUE)
    wire.cell.type.selector(input, session, "chord.sources", ident.choices, selected.all = TRUE)
    wire.cell.type.selector(input, session, "chord.targets", ident.choices, selected.all = TRUE)

    # Auto-size plot output for each subtab
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

    # ---------- Subtab thunks (Phase 5/6 will implement) ----------

    dot.thunk <- function() {
      res <- liana.data(); req(res, input$dot.group)
      plot.new()
      title("CCC Dot Plot — will be implemented in Phase 5")
    }

    heat.thunk <- function() {
      res <- liana.data(); req(res)
      plot.new()
      title("CCC Freq Heatmap — will be implemented in Phase 6")
    }

    chord.thunk <- function() {
      res <- liana.data(); req(res)
      plot.new()
      title("CCC Freq Chord Diagram — will be implemented in Phase 6")
    }

    make.render.and.download(output, session, "dot.plot", "dot.download",
      dot.thunk,
      width.in = reactive(input$dot.width),
      height.in = reactive(input$dot.height),
      filename.stem = reactive(paste0("liana_dot_",
                                      gsub("[^A-Za-z0-9]", "_",
                                           input$dot.group %||% ""))),
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
