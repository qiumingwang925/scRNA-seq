## ABOUTME: Enrichment analysis module using enrichR on DE gene results.
## ABOUTME: Two sub-tabs: enrichment table with row selection, and bar plot visualization.

mod.explore.enrich.ui <- function(id) {
  ns <- NS(id)
  tabPanel("Enrichment Analysis",
    tabsetPanel(
      # Sub-tab 1: Enrichment Table
      tabPanel("Enrichment Table",
        sidebarLayout(
          sidebarPanel(width = 4,
            selectInput(ns("database"), "Enrichr Database:",
                        choices = c("GO_Biological_Process_2023",
                                    "GO_Molecular_Function_2023",
                                    "GO_Cellular_Component_2023",
                                    "KEGG_2021_Human",
                                    "MSigDB_Hallmark_2020",
                                    "Reactome_2022",
                                    "WikiPathway_2023_Human"),
                        selected = "GO_Biological_Process_2023"),
            actionButton(ns("run.ea"), "Run Enrichment Analysis",
                         class = "btn-success", style = "width:100%"),
            hr(),
            verbatimTextOutput(ns("selected.ea")),
            hr(),
            downloadButton(ns("download.ea.table"), "Download Table (.xlsx)",
                           class = "btn-success")
          ),
          mainPanel(width = 8,
            DT::DTOutput(ns("ea.table"))
          )
        )
      ),
      # Sub-tab 2: Bar Plot
      tabPanel("Bar Plot",
        sidebarLayout(
          sidebarPanel(width = 4,
            radioButtons(ns("bar.source"), "Pathways to Plot:",
                         choices = c("Top 10" = "top10", "Selected Rows" = "selected"),
                         inline = TRUE),
            selectInput(ns("bar.color"), "Bar Color:",
                        choices = c("red", "green", "blue", "purple", "orange"),
                        selected = "red"),
            actionButton(ns("run.bar"), "Generate Bar Plot",
                         class = "btn-success", style = "width:100%"),
            hr(),
            h4("Download Figure"),
            numericInput(ns("fig.w"), "Width (inches)", value = 10, min = 2, max = 30),
            numericInput(ns("fig.h"), "Height (inches)", value = 6, min = 2, max = 20),
            downloadButton(ns("download.bar"), "Download Figure", class = "btn-success")
          ),
          mainPanel(width = 8,
            plotOutput(ns("plot.bar"), height = "600px")
          )
        )
      )
    )
  )
}

mod.explore.enrich.server <- function(id, de.result) {
  moduleServer(id, function(input, output, session) {

    ea.table <- eventReactive(input$run.ea, {
      req(de.result())
      validate(need(nrow(de.result()) > 0, "No DE genes available. Run DE analysis first."))

      withProgress(message = "Running enrichment analysis...", value = 0.3, {
        tryCatch({
          genes <- de.result()$gene
          result <- enrichR::enrichr(genes, input$database)[[1]]
          incProgress(0.5, detail = "Filtering results")

          result <- result %>% filter(P.value < 0.05)
          validate(need(nrow(result) > 0, "No significant enrichment results (P < 0.05)."))

          result$P.value <- as.numeric(format(result$P.value, scientific = TRUE, digits = 3))
          result$Adjusted.P.value <- as.numeric(format(result$Adjusted.P.value, scientific = TRUE, digits = 3))
          result$Odds.Ratio <- as.numeric(format(result$Odds.Ratio, scientific = FALSE, digits = 3))
          result$Combined.Score <- as.numeric(format(result$Combined.Score, scientific = FALSE, digits = 4))

          incProgress(0.2, detail = "Done")
          result
        }, error = function(e) {
          message("Enrichment error: ", e$message)
          showNotification(paste("Enrichment error:", e$message), type = "error")
          NULL
        })
      })
    })

    output$ea.table <- DT::renderDT({
      req(ea.table())
      DT::datatable(ea.table(), rownames = FALSE,
                    selection = "multiple",
                    options = list(pageLength = 10, scrollX = TRUE))
    })

    # Display selected pathway names
    output$selected.ea <- renderPrint({
      s <- input$ea.table_rows_selected
      if (length(s)) {
        cat("Selected Pathways:\n\n")
        cat(ea.table()[s, ] %>% pull(Term), sep = "\n")
      } else {
        cat("Click rows in the table to select pathways.")
      }
    })

    output$download.ea.table <- downloadHandler(
      filename = function() { paste0("enrichment_", input$database, ".xlsx") },
      content = function(file) {
        openxlsx::write.xlsx(ea.table(), file)
      }
    )

    # Bar plot
    plot.bar <- eventReactive(input$run.bar, {
      req(ea.table())
      df.bar <- ea.table()
      df.bar$neg.log10.padj <- -log10(df.bar$Adjusted.P.value)

      if (input$bar.source == "top10") {
        df.bar <- head(df.bar, 10)
      } else {
        s <- input$ea.table_rows_selected
        validate(need(length(s) > 0, "Please select rows in the Enrichment Table first."))
        df.bar <- df.bar[s, ]
      }

      ggplot(df.bar, aes(x = neg.log10.padj,
                         y = reorder(Term, neg.log10.padj))) +
        geom_bar(stat = "identity", aes(fill = Odds.Ratio)) +
        scale_fill_gradient2(high = input$bar.color) +
        theme_classic() +
        labs(x = "-Log10(Adjusted P-value)", y = "")
    })

    output$plot.bar <- renderPlot({ plot.bar() }, res = 96)

    output$download.bar <- downloadHandler(
      filename = function() { paste0("enrichment_barplot_", input$database, ".pdf") },
      content = function(file) {
        ggsave(file, plot = plot.bar(),
               width = input$fig.w, height = input$fig.h)
      }
    )
  })
}
