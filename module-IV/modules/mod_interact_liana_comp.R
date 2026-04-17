## ABOUTME: LIANA computation module. Splits the uploaded Seurat object by group,
## ABOUTME: runs LIANA per group across selected methods, and aggregates the tibble.

mod.interact.liana.comp.ui <- function(id) {
  ns <- NS(id)
  tabPanel("Run LIANA",
    sidebarLayout(
      sidebarPanel(width = 4,
        radioButtons(ns("species"), "Species:",
                     choices = c("Mouse" = "mouse", "Human" = "human"),
                     selected = "mouse"),
        radioButtons(ns("assay"), "Assay:",
                     choices = c("RNA" = "RNA", "SCT" = "SCT"),
                     selected = "RNA"),
        selectInput(ns("resource"), "CCC resource:", choices = NULL),
        checkboxGroupInput(ns("methods"), "CCC methods:",
                           choices = c("natmi", "connectome", "logfc", "sca",
                                       "cellphonedb", "cytotalk", "call_cellchat"),
                           selected = c("natmi", "connectome", "sca", "cellphonedb")),
        numericInput(ns("workers"), "Parallel workers:",
                     value = 1, min = 1, max = 16, step = 1),
        numericInput(ns("min.cells"), "Min cells per group:",
                     value = 10, min = 1, max = 10000, step = 1),
        hr(),
        actionButton(ns("run"), "Run LIANA Analysis",
                     class = "btn-success", style = "width:100%"),
        hr(),
        downloadButton(ns("download"), "Download Result (.rds)",
                       style = "width:100%"),
        hr(),
        helpText("MouseConsensus is mouse-native and skips gene conversion.",
                 "Other resources trigger a mouse-human biomaRt mapping",
                 "when species is Mouse; results are mapped back to mouse",
                 "symbols before returning.")
      ),
      mainPanel(width = 8,
        h4("Pipeline Status"),
        verbatimTextOutput(ns("status.log")),
        hr(),
        h4("Per-Group Summary"),
        tableOutput(ns("summary.table"))
      )
    )
  )
}

mod.interact.liana.comp.server <- function(id, shared.data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Session-scoped gene conversion cache
    gene.conv.cache <- new.env(parent = emptyenv())

    # Populate resource dropdown from liana once
    observe({
      resources <- tryCatch(liana::show_resources(), error = function(e) character())
      resources <- setdiff(resources, "Default")
      if (length(resources) == 0) resources <- c("MouseConsensus", "Consensus")
      updateSelectInput(session, "resource",
                        choices = resources,
                        selected = if ("MouseConsensus" %in% resources)
                          "MouseConsensus" else resources[1])
    })

    # Track default assay from uploaded object
    observeEvent(shared.data(), {
      obj <- shared.data()
      req(obj)
      updateRadioButtons(session, "assay", selected = DefaultAssay(obj))
    }, ignoreInit = TRUE)

    liana.result <- eventReactive(input$run, {
      obj <- shared.data()
      req(obj)
      validate(need(length(input$methods) > 0, "Select at least one method."))

      # Copy Idents into an explicit column so liana_wrap finds cell types regardless
      # of whether the active Idents chain is preserved through subset().
      obj$liana_idents <- as.character(Idents(obj))

      groups <- sort(unique(obj$group))
      species <- input$species
      assay.name <- input$assay
      resource <- input$resource
      methods <- input$methods
      workers <- input$workers
      min.cells <- input$min.cells

      needs.conversion <- species == "mouse" && resource != "MouseConsensus"

      future::plan("multisession", workers = workers)
      options(future.globals.maxSize = 10 * 1024^3)

      withProgress(message = "Running LIANA...", value = 0, {
        n <- length(groups)
        step <- 1 / n
        liana.list <- list()
        conv.flag <- list()
        method.succeeded <- list()

        for (i in seq_along(groups)) {
          g <- groups[i]
          incProgress(0, detail = paste0("[", g, "] subsetting"))
          obj.sub <- subset(obj, subset = group == g)
          Idents(obj.sub) <- droplevels(Idents(obj.sub))

          if (ncol(obj.sub) < min.cells) {
            showNotification(
              paste0("Group '", g, "' has ", ncol(obj.sub), " cells (< min ",
                     min.cells, "). Skipping."),
              type = "warning", duration = NULL
            )
            next
          }
          if (length(levels(Idents(obj.sub))) < 2) {
            showNotification(
              paste0("Group '", g, "' has fewer than 2 cell-type idents. Skipping."),
              type = "warning", duration = NULL
            )
            next
          }

          if (needs.conversion) {
            incProgress(step * 0.2, detail = paste0("[", g, "] mouse->human gene map"))
            converted <- tryCatch(
              convert.mouse.to.human.rownames(obj.sub, gene.conv.cache),
              error = function(e) {
                showNotification(paste("biomaRt conversion failed for", g, ":",
                                       e$message),
                                 type = "error", duration = NULL)
                NULL
              }
            )
            if (is.null(converted)) next
            obj.sub <- converted$obj
          }

          # Build a SingleCellExperiment manually using the Seurat v5 `layer` API.
          # liana_wrap on a Seurat object calls GetAssayData(slot=...), which is
          # defunct under SeuratObject >= 5.0, so we sidestep that code path.
          incProgress(step * 0.1, detail = paste0("[", g, "] Seurat -> SCE"))
          sce <- tryCatch({
            data.mat <- GetAssayData(obj.sub, assay = assay.name, layer = "data")
            counts.mat <- tryCatch(
              GetAssayData(obj.sub, assay = assay.name, layer = "counts"),
              error = function(e) data.mat
            )
            sce.obj <- SingleCellExperiment::SingleCellExperiment(
              assays = list(counts = counts.mat, logcounts = data.mat),
              colData = obj.sub@meta.data
            )
            SingleCellExperiment::colLabels(sce.obj) <-
              as.character(obj.sub$liana_idents)
            sce.obj
          }, error = function(e) {
            showNotification(paste("[", g, "] SCE build failed:", e$message),
                             type = "error", duration = NULL)
            NULL
          })
          if (is.null(sce)) next

          # Run each method separately so one bad method doesn't collapse the group
          incProgress(step * 0.3, detail = paste0(
            "[", g, "] running LIANA (", length(methods), " methods)"))
          per.method <- list()
          for (m in methods) {
            tib <- tryCatch({
              result <- liana::liana_wrap(
                sce = sce,
                method = m,
                resource = resource,
                assay.type = "logcounts"
              )
              if (is.list(result) && !inherits(result, "data.frame") &&
                  !inherits(result, "tbl")) {
                result <- result[[1]]
              }
              result
            }, error = function(e) {
              message("[", g, "] method '", m, "' failed: ", e$message)
              NULL
            })
            if (!is.null(tib)) per.method[[m]] <- tib
          }

          if (length(per.method) == 0) {
            showNotification(
              paste0("No methods succeeded for group '", g, "'. Skipping."),
              type = "error", duration = NULL
            )
            next
          }

          incProgress(step * 0.3, detail = paste0("[", g, "] aggregating"))
          aggregated <- tryCatch(
            liana::liana_aggregate(per.method),
            error = function(e) {
              message("[", g, "] liana_aggregate failed: ", e$message)
              per.method[[1]]
            }
          )

          if (needs.conversion) {
            incProgress(step * 0.2, detail = paste0("[", g, "] human->mouse gene remap"))
            aggregated <- tryCatch(
              convert.human.to.mouse.lr(aggregated, gene.conv.cache),
              error = function(e) {
                showNotification(paste("Reverse conversion failed for", g, ":",
                                       e$message),
                                 type = "warning")
                aggregated
              }
            )
          }

          liana.list[[g]] <- aggregated
          conv.flag[[g]] <- needs.conversion
          method.succeeded[[g]] <- names(per.method)
        }

        if (length(liana.list) == 0) {
          showNotification("LIANA produced no results for any group.",
                           type = "error", duration = NULL)
          return(NULL)
        }

        list(
          liana.list = liana.list,
          group.levels = names(liana.list),
          species = species,
          assay = assay.name,
          resource = resource,
          methods = methods,
          method.succeeded = method.succeeded,
          converted = conv.flag
        )
      })
    })

    output$status.log <- renderPrint({
      res <- liana.result()
      req(res)
      cat("--- LIANA run complete ---\n")
      cat("Species:           ", res$species, "\n")
      cat("Assay:             ", res$assay, "\n")
      cat("Resource:          ", res$resource, "\n")
      cat("Methods requested: ", paste(res$methods, collapse = ", "), "\n")
      cat("Groups completed:  ", paste(res$group.levels, collapse = ", "), "\n")
    })

    output$summary.table <- renderTable({
      res <- liana.result()
      req(res)
      do.call(rbind, lapply(res$group.levels, function(g) {
        tib <- res$liana.list[[g]]
        data.frame(
          Group = g,
          CellTypes = length(unique(c(tib$source, tib$target))),
          Methods_Succeeded = length(res$method.succeeded[[g]]),
          LR_Rows = nrow(tib),
          Converted = if (isTRUE(res$converted[[g]])) "yes" else "no",
          stringsAsFactors = FALSE
        )
      }))
    })

    output$download <- downloadHandler(
      filename = function() paste0("liana_result_", Sys.Date(), ".rds"),
      content = function(file) {
        res <- liana.result()
        req(res)
        saveRDS(res, file)
      }
    )

    return(liana.result)
  })
}
