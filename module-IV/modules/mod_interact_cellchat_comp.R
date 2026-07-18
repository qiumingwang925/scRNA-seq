## ABOUTME: CellChat computation module. Splits the Seurat object by group, runs the
## ABOUTME: CellChat pipeline per group, and merges results for comparison visualization.

mod.interact.cellchat.comp.ui <- function(id) {
  ns <- NS(id)
  tabPanel("Run CellChat",
    sidebarLayout(
      sidebarPanel(width = 4,
        radioButtons(ns("species"), "Species:",
                     choices = c("Mouse" = "mouse", "Human" = "human"),
                     selected = "mouse"),
        radioButtons(ns("assay"), "Assay:",
                     choices = c("RNA" = "RNA", "SCT" = "SCT"),
                     selected = "RNA"),
        hr(),
        actionButton(ns("run"), "Run CellChat Analysis",
                     class = "btn-success", style = "width:100%"),
        hr(),
        downloadButton(ns("download"), "Download Result (.rds)",
                       style = "width:100%"),
        hr(),
        helpText("Signaling database is fixed to 'Secreted Signaling'.",
                 "CellChat runs once per group level, then merges.",
                 "Large datasets (>50k cells) may take 30+ minutes per group.")
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

mod.interact.cellchat.comp.server <- function(id, shared.data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observeEvent(shared.data(), {
      obj <- shared.data()
      req(obj)
      updateRadioButtons(session, "assay", selected = DefaultAssay(obj))
    }, ignoreInit = TRUE)

    cellchat.result <- eventReactive(input$run, {
      obj <- shared.data()
      req(obj)

      groups <- sort(unique(obj$group))
      species <- input$species
      assay.name <- input$assay

      cellchat.db <- if (species == "mouse") CellChatDB.mouse else CellChatDB.human
      cellchat.db <- subsetDB(cellchat.db, search = "Secreted Signaling", key = "annotation")

      future::plan("sequential")
      options(future.globals.maxSize = 10 * 1024^3)

      withProgress(message = "Running CellChat...", value = 0, {
        n <- length(groups)
        step <- 1 / n
        cellchat.list <- list()

        for (i in seq_along(groups)) {
          g <- groups[i]
          base <- (i - 1) * step

          incProgress(0, detail = paste0("[", g, "] subsetting"))
          obj.sub <- subset(obj, subset = group == g)
          Idents(obj.sub) <- droplevels(Idents(obj.sub))

          if (length(levels(Idents(obj.sub))) < 2) {
            showNotification(
              paste0("Group '", g, "' has fewer than 2 cell-type idents after subsetting. ",
                     "CellChat requires at least 2 cell types per group."),
              type = "error", duration = NULL
            )
            return(NULL)
          }

          cc <- tryCatch({
            incProgress(step * 0.1, detail = paste0("[", g, "] creating CellChat object"))
            x <- createCellChat(object = obj.sub, group.by = "ident", assay = assay.name)
            x@DB <- cellchat.db
            x <- subsetData(x)

            incProgress(step * 0.2, detail = paste0("[", g, "] over-expression analysis"))
            x <- identifyOverExpressedGenes(x)
            x <- identifyOverExpressedInteractions(x)

            incProgress(step * 0.4, detail = paste0("[", g, "] computing communication probabilities"))
            x <- computeCommunProb(x, type = "triMean")
            x <- computeCommunProbPathway(x)

            incProgress(step * 0.2, detail = paste0("[", g, "] aggregating network"))
            x <- aggregateNet(x)

            incProgress(step * 0.1, detail = paste0("[", g, "] done"))
            x
          }, error = function(e) {
            showNotification(
              paste0("CellChat failed for group '", g, "': ", e$message),
              type = "error", duration = NULL
            )
            message("CellChat error for group ", g, ": ", e$message)
            NULL
          })

          if (is.null(cc)) return(NULL)
          cellchat.list[[g]] <- cc
        }

        cellchat.merged <- tryCatch({
          mergeCellChat(cellchat.list, add.names = names(cellchat.list))
        }, error = function(e) {
          showNotification(paste("mergeCellChat failed:", e$message),
                           type = "error", duration = NULL)
          NULL
        })

        list(
          cellchat.list = cellchat.list,
          cellchat.merged = cellchat.merged,
          group.levels = groups,
          species = species,
          assay = assay.name
        )
      })
    })

    output$status.log <- renderPrint({
      res <- cellchat.result()
      req(res)
      cat("--- CellChat run complete ---\n")
      cat("Species:       ", res$species, "\n")
      cat("Assay:         ", res$assay, "\n")
      cat("Database:      ", "Secreted Signaling\n")
      cat("Groups:        ", paste(res$group.levels, collapse = ", "), "\n")
      cat("Merged object: ", if (is.null(res$cellchat.merged)) "FAILED" else "OK", "\n")
    })

    output$summary.table <- renderTable({
      res <- cellchat.result()
      req(res)
      do.call(rbind, lapply(names(res$cellchat.list), function(g) {
        cc <- res$cellchat.list[[g]]
        data.frame(
          Group = g,
          Cells = length(cc@idents),
          CellTypes = length(levels(cc@idents)),
          Pathways = length(cc@netP$pathways),
          LR_Pairs = if (!is.null(cc@LR$LRsig)) nrow(cc@LR$LRsig) else 0,
          stringsAsFactors = FALSE
        )
      }))
    })

    output$download <- downloadHandler(
      filename = function() paste0("cellchat_result_", Sys.Date(), ".rds"),
      content = function(file) {
        res <- cellchat.result()
        req(res)
        saveRDS(res, file)
      }
    )

    return(cellchat.result)
  })
}
