## ABOUTME: Violin plot module for gene expression across cell types and metadata groups.
## ABOUTME: Supports single/stacked genes, split.by metadata, and figure download.

mod.explore.violin.ui <- function(id) {
  ns <- NS(id)
  tabPanel("Violin Plot",
    sidebarLayout(
      sidebarPanel(width = 4,
        selectInput(ns("select.idents"), "Cell Type(s):",
                    choices = NULL, multiple = TRUE),
        fluidRow(
          column(6, actionButton(ns("btn.select.all"), "Select All",
                                 class = "btn-info", style = "width:100%")),
          column(6, actionButton(ns("btn.clear.all"), "Clear",
                                 class = "btn-default", style = "width:100%"))
        ),
        hr(),
        textInput(ns("gene.input"), "Gene(s) (comma-separated)",
                  placeholder = "e.g. Cd68, Cx3cr1, Ccr2"),
        checkboxInput(ns("flip"), "Flip Axes", value = FALSE),
        selectInput(ns("split.by"), "Split By (optional):",
                    choices = c("None"), selected = "None"),
        actionButton(ns("run.vln"), "Generate Violin Plot",
                     class = "btn-success", style = "width:100%"),
        hr(),
        h4("Download Figure"),
        numericInput(ns("fig.w"), "Width (inches)", value = 8, min = 2, max = 20),
        numericInput(ns("fig.h"), "Height (inches)", value = 6, min = 2, max = 20),
        downloadButton(ns("download.vln"), "Download Figure", class = "btn-success")
      ),
      mainPanel(width = 8,
        plotOutput(ns("plot.vln"), height = "600px")
      )
    )
  )
}

mod.explore.violin.server <- function(id, shared.data) {
  moduleServer(id, function(input, output, session) {

    # Populate cell types and metadata on data load
    observe({
      req(shared.data())
      obj <- shared.data()
      ident.levels <- levels(obj)
      updateSelectInput(session, "select.idents",
                        choices = ident.levels, selected = ident.levels)
      cat.cols <- get.categorical.meta(obj)
      updateSelectInput(session, "split.by",
                        choices = c("None", cat.cols), selected = "None")
    })

    observeEvent(input$btn.select.all, {
      req(shared.data())
      updateSelectInput(session, "select.idents",
                        choices = levels(shared.data()),
                        selected = levels(shared.data()))
    })

    observeEvent(input$btn.clear.all, {
      req(shared.data())
      updateSelectInput(session, "select.idents",
                        choices = levels(shared.data()),
                        selected = character(0))
    })

    plot.vln <- eventReactive(input$run.vln, {
      req(shared.data())
      obj <- shared.data()

      # Parse genes
      genes <- trimws(unlist(strsplit(input$gene.input, ",")))
      genes <- genes[nchar(genes) > 0]
      validate(need(length(genes) > 0, "Please enter at least one gene name."))

      # Validate genes exist
      missing <- genes[!genes %in% rownames(obj)]
      validate(need(length(missing) == 0,
                    paste0("Gene(s) not found: ", paste(missing, collapse = ", "))))

      # Subset by selected cell types
      idents.selected <- input$select.idents
      validate(need(length(idents.selected) > 0, "Please select at least one cell type."))
      obj <- subset(obj, idents = idents.selected)

      split.by <- if (input$split.by == "None") NULL else input$split.by

      withProgress(message = "Generating violin plot...", value = 0.5, {
        if (length(genes) == 1) {
          p <- VlnPlot(obj, features = genes, split.by = split.by)
        } else {
          p <- VlnPlot(obj, features = genes, stack = TRUE,
                       flip = input$flip, split.by = split.by)
        }
        incProgress(0.5, detail = "Done")
        p
      })
    })

    output$plot.vln <- renderPlot({ plot.vln() }, res = 96)

    output$download.vln <- downloadHandler(
      filename = function() { "violin_plot.pdf" },
      content = function(file) {
        ggsave(file, plot = plot.vln(),
               width = input$fig.w, height = input$fig.h)
      }
    )
  })
}
