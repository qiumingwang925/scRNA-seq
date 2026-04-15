## ABOUTME: Dot plot module for gene expression across cell types with optional split.by.
## ABOUTME: Supports Seurat DotPlot or manual per-identity color scaling when split.by is used.

mod.explore.dot.ui <- function(id) {
  ns <- NS(id)
  tabPanel("Dot Plot",
    sidebarLayout(
      sidebarPanel(width = 4,
        selectInput(ns("select.idents"), "Cell Type(s):",
                    choices = NULL, multiple = TRUE),
        checkboxInput(ns("show.all.idents"), "Select All", value = TRUE),
        hr(),
        textInput(ns("gene.input"), "Gene(s) (comma-separated)",
                  placeholder = "e.g. Cd68, Cx3cr1, Ccr2, Ly6c2"),
        selectInput(ns("split.by"), "Split By (optional):",
                    choices = c("None"), selected = "None"),
        actionButton(ns("run.dot"), "Generate Dot Plot",
                     class = "btn-success", style = "width:100%"),
        hr(),
        h4("Download Figure"),
        numericInput(ns("fig.w"), "Width (inches)", value = 10, min = 2, max = 30),
        numericInput(ns("fig.h"), "Height (inches)", value = 6, min = 2, max = 20),
        downloadButton(ns("download.dot"), "Download Figure", class = "btn-success")
      ),
      mainPanel(width = 8,
        plotOutput(ns("plot.dot"), height = "600px")
      )
    )
  )
}

mod.explore.dot.server <- function(id, shared.data) {
  moduleServer(id, function(input, output, session) {

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

    observeEvent(input$show.all.idents, {
      req(shared.data())
      if (input$show.all.idents) {
        updateSelectInput(session, "select.idents",
                          choices = levels(shared.data()),
                          selected = levels(shared.data()))
      }
    })

    plot.dot <- eventReactive(input$run.dot, {
      req(shared.data())
      obj <- shared.data()

      genes <- trimws(unlist(strsplit(input$gene.input, ",")))
      genes <- genes[nchar(genes) > 0]
      validate(need(length(genes) > 0, "Please enter at least one gene name."))

      missing <- genes[!genes %in% rownames(obj)]
      validate(need(length(missing) == 0,
                    paste0("Gene(s) not found: ", paste(missing, collapse = ", "))))

      idents.selected <- if (input$show.all.idents) levels(obj) else input$select.idents
      validate(need(length(idents.selected) > 0, "Please select at least one cell type."))
      obj <- subset(obj, idents = idents.selected)

      split.by <- if (input$split.by == "None") NULL else input$split.by

      withProgress(message = "Generating dot plot...", value = 0.3, {
        if (is.null(split.by)) {
          p <- DotPlot(obj, features = genes) + RotatedAxis()
          incProgress(0.7, detail = "Done")
          p
        } else {
          # Manual computation for per-identity color scaling
          incProgress(0.2, detail = "Computing per-identity statistics")
          p <- build.split.dot.plot(obj, genes, split.by)
          incProgress(0.5, detail = "Done")
          p
        }
      })
    })

    output$plot.dot <- renderPlot({ plot.dot() }, res = 96)

    output$download.dot <- downloadHandler(
      filename = function() { "dot_plot.pdf" },
      content = function(file) {
        ggsave(file, plot = plot.dot(),
               width = input$fig.w, height = input$fig.h)
      }
    )
  })
}

# Compute average expression and percent expressed per cell type per split identity,
# then render a faceted dot plot with independent color scales per identity.
build.split.dot.plot <- function(obj, genes, split.by) {
  meta <- obj@meta.data
  meta$cell.type <- Idents(obj)
  meta$split.id <- as.character(meta[[split.by]])

  # Get expression data
  assay.data <- GetAssayData(obj, slot = "data")

  # Build stats table
  stats.list <- list()
  for (gene in genes) {
    expr.vals <- as.numeric(assay.data[gene, ])
    for (ct in unique(meta$cell.type)) {
      for (sid in unique(meta$split.id)) {
        idx <- which(meta$cell.type == ct & meta$split.id == sid)
        if (length(idx) == 0) next
        vals <- expr.vals[idx]
        stats.list[[length(stats.list) + 1]] <- data.frame(
          gene = gene,
          cell.type = ct,
          split.id = sid,
          avg.exp = mean(vals),
          pct.exp = sum(vals > 0) / length(vals) * 100,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  stats.df <- do.call(rbind, stats.list)

  # Assign color palette per split identity
  split.levels <- unique(stats.df$split.id)
  identity.colors <- scales::hue_pal()(length(split.levels))
  names(identity.colors) <- split.levels

  # Scale avg.exp within each split identity to [0, 1] for color mapping
  stats.df <- stats.df %>%
    group_by(split.id) %>%
    mutate(avg.exp.scaled = if (max(avg.exp) == min(avg.exp)) 0.5
           else (avg.exp - min(avg.exp)) / (max(avg.exp) - min(avg.exp))) %>%
    ungroup()

  # Build per-identity color columns (light grey to identity color)
  stats.df$fill.color <- mapply(function(scaled, sid) {
    grDevices::colorRamp(c("lightgrey", identity.colors[sid]))(scaled) %>%
      { grDevices::rgb(.[1], .[2], .[3], maxColorValue = 255) }
  }, stats.df$avg.exp.scaled, stats.df$split.id)

  stats.df$gene <- factor(stats.df$gene, levels = genes)
  stats.df$cell.type <- factor(stats.df$cell.type, levels = levels(Idents(obj)))

  ggplot(stats.df, aes(x = gene, y = cell.type, size = pct.exp, color = fill.color)) +
    geom_point() +
    scale_color_identity() +
    scale_size_continuous(range = c(1, 8), name = "% Expressed") +
    facet_wrap(~ split.id) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Gene", y = "Cell Type")
}
