## ABOUTME: LIANA computation module. Splits the uploaded Seurat object by group and runs
## ABOUTME: LIANA per group, storing one raw table per method. Nothing is aggregated here —
## ABOUTME: liana_aggregate runs in the vis module over a user-chosen subset of methods.

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
      min.cells <- input$min.cells

      needs.conversion <- species == "mouse" && resource != "MouseConsensus"

      # CellChat picks its geneInfo/cofactor tables from this, so it must describe the
      # symbols in the matrix it receives -- which are human whenever we converted.
      cellchat.organism <- if (needs.conversion) "human" else species

      future::plan("sequential")
      options(future.globals.maxSize = 10 * 1024^3)

      withProgress(message = "Running LIANA...", value = 0, {
        n <- length(groups)
        step <- 1 / n
        liana.list <- list()
        conv.flag <- list()
        # Per requested method: NA when it produced a table, else the failure reason.
        method.status <- list()
        # Groups that never reached the method loop, keyed by group -> reason.
        group.skipped <- list()

        for (i in seq_along(groups)) {
          g <- groups[i]
          incProgress(0, detail = paste0("[", g, "] subsetting"))
          obj.sub <- subset(obj, subset = group == g)
          Idents(obj.sub) <- droplevels(Idents(obj.sub))

          if (ncol(obj.sub) < min.cells) {
            reason <- paste0(ncol(obj.sub), " cells (< min ", min.cells, ")")
            group.skipped[[g]] <- reason
            showNotification(paste0("Group '", g, "' has ", reason, ". Skipping."),
                             type = "warning", duration = NULL)
            next
          }
          if (length(levels(Idents(obj.sub))) < 2) {
            group.skipped[[g]] <- "fewer than 2 cell-type idents"
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
                msg <- paste0("biomaRt conversion failed for ", g, ": ",
                              e$message,
                              "\nEnsembl may be unreachable. Try again shortly, ",
                              "or switch Resource to 'MouseConsensus' to skip ",
                              "conversion entirely.")
                message(msg)
                showNotification(msg, type = "error", duration = NULL)
                NULL
              }
            )
            if (is.null(converted)) {
              group.skipped[[g]] <- "mouse->human gene conversion failed"
              next
            }
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
              error = function(e) NULL
            )
            # Slimmed-down objects return an empty 0x0 matrix rather than
            # erroring; fall back to data in either case.
            if (is.null(counts.mat) ||
                nrow(counts.mat) != nrow(data.mat) ||
                ncol(counts.mat) != ncol(data.mat)) {
              counts.mat <- data.mat
            }
            sce.obj <- SingleCellExperiment::SingleCellExperiment(
              assays = list(counts = counts.mat, logcounts = data.mat),
              colData = obj.sub@meta.data
            )
            SingleCellExperiment::colLabels(sce.obj) <-
              as.character(obj.sub$liana_idents)
            sce.obj
          }, error = function(e) {
            msg <- paste0("[", g, "] SCE build failed: ", e$message)
            message(msg)
            showNotification(msg, type = "error", duration = NULL)
            NULL
          })
          if (is.null(sce)) {
            group.skipped[[g]] <- "Seurat -> SCE conversion failed"
            next
          }

          # Run each method separately so one bad method doesn't collapse the group
          incProgress(step * 0.3, detail = paste0(
            "[", g, "] running LIANA (", length(methods), " methods)"))
          per.method <- list()
          status <- setNames(rep(NA_character_, length(methods)), methods)
          for (m in methods) {
            tib <- tryCatch({
              # liana_defaults() builds every method's parameter list unconditionally,
              # so cellchat.params is inert for the other methods.
              result <- liana::liana_wrap(
                sce = sce,
                method = m,
                resource = resource,
                assay.type = "logcounts",
                cellchat.params = list(organism = cellchat.organism)
              )
              if (is.list(result) && !inherits(result, "data.frame") &&
                  !inherits(result, "tbl")) {
                result <- result[[1]]
              }
              result
            }, error = function(e) {
              message("[", g, "] method '", m, "' failed: ", e$message)
              status[[m]] <<- conditionMessage(e)
              NULL
            })
            if (!is.null(tib)) per.method[[m]] <- tib
          }

          if (length(per.method) == 0) {
            group.skipped[[g]] <- "no method succeeded"
            showNotification(
              paste0("No methods succeeded for group '", g, "'. Skipping."),
              type = "error", duration = NULL
            )
            next
          }

          if (needs.conversion) {
            incProgress(step * 0.4, detail = paste0("[", g, "] human->mouse gene remap"))
            per.method <- setNames(lapply(names(per.method), function(m) {
              tryCatch(
                convert.human.to.mouse.lr(per.method[[m]], gene.conv.cache),
                error = function(e) {
                  showNotification(paste0("Reverse conversion failed for ", g,
                                          " / ", m, ": ", e$message),
                                   type = "warning")
                  per.method[[m]]
                }
              )
            }), names(per.method))
          }

          liana.list[[g]] <- per.method
          conv.flag[[g]] <- needs.conversion
          method.status[[g]] <- status
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
          method.status = method.status,
          group.skipped = group.skipped,
          converted = conv.flag
        )
      })
    })

    output$status.log <- renderPrint({
      res <- liana.result()
      req(res)
      shared <- Reduce(intersect, lapply(res$liana.list, names))
      partial <- setdiff(Reduce(union, lapply(res$liana.list, names)), shared)

      cat("--- LIANA run complete ---\n")
      cat("Species:           ", res$species, "\n")
      cat("Assay:             ", res$assay, "\n")
      cat("Resource:          ", res$resource, "\n")
      cat("Methods requested: ", paste(res$methods, collapse = ", "), "\n")
      cat("Groups completed:  ", paste(res$group.levels, collapse = ", "), "\n")
      cat("Methods in every group:", paste(shared, collapse = ", "), "\n")
      if (length(partial) > 0) {
        cat("Methods missing from at least one group (not selectable downstream):",
            paste(partial, collapse = ", "), "\n")
      }
      if (length(res$group.skipped) > 0) {
        cat("\nGroups skipped:\n")
        for (g in names(res$group.skipped)) {
          cat("  ", g, ": ", res$group.skipped[[g]], "\n", sep = "")
        }
      }
      failures <- unlist(lapply(res$group.levels, function(g) {
        st <- res$method.status[[g]]
        failed <- names(st)[!is.na(st)]
        if (length(failed) == 0) return(NULL)
        paste0("  ", g, " / ", failed, ": ", st[failed])
      }))
      if (length(failures) > 0) {
        cat("\nMethod failures:\n")
        cat(paste(failures, collapse = "\n"), "\n")
      }
    })

    output$summary.table <- renderTable({
      res <- liana.result()
      req(res)
      do.call(rbind, lapply(res$group.levels, function(g) {
        st <- res$method.status[[g]]
        do.call(rbind, lapply(names(st), function(m) {
          tib <- res$liana.list[[g]][[m]]
          data.frame(
            Group = g,
            Method = m,
            Status = if (is.null(tib)) "failed" else "ok",
            CellTypes = if (is.null(tib)) NA_integer_
                        else length(unique(c(tib$source, tib$target))),
            LR_Rows = if (is.null(tib)) NA_integer_ else nrow(tib),
            Converted = if (isTRUE(res$converted[[g]])) "yes" else "no",
            stringsAsFactors = FALSE
          )
        }))
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
