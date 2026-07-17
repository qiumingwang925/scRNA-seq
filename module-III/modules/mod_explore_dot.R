## ABOUTME: Dot plot module for gene expression across cell types with optional split.by.
## ABOUTME: Supports Seurat DotPlot or manual per-identity color scaling when split.by is used.

mod.explore.dot.ui <- function(id) {
  ns <- NS(id)
  tabPanel("Dot Plot",
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
      updateSelectInput(session, "split.by",
                        choices = c("None", split.by.choices(obj)), selected = "None")
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

    plot.dot <- eventReactive(input$run.dot, {
      req(shared.data())
      obj <- shared.data()

      genes <- trimws(unlist(strsplit(input$gene.input, ",")))
      genes <- genes[nchar(genes) > 0]
      validate(need(length(genes) > 0, "Please enter at least one gene name."))

      missing <- genes[!genes %in% rownames(obj)]
      validate(need(length(missing) == 0,
                    paste0("Gene(s) not found: ", paste(missing, collapse = ", "))))

      idents.selected <- input$select.idents
      validate(need(length(idents.selected) > 0, "Please select at least one cell type."))
      obj <- subset(obj, idents = idents.selected)

      split.by <- if (input$split.by == "None") NULL else input$split.by

      withProgress(message = "Generating dot plot...", value = 0.3, {
        if (is.null(split.by)) {
          p <- DotPlot(obj, features = genes) + RotatedAxis() + theme_classic()
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
  assay.data <- GetAssayData(obj, layer = "data")

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

  stats.df$gene <- factor(stats.df$gene, levels = genes)

  # Combine cell type + split identity into a single y-axis label
  # Interleaved order (bottom to top): IM group1, IM group2, cMono group1, cMono group2, ...
  ct.levels <- levels(Idents(obj))
  split.levels <- unique(stats.df$split.id)
  y.levels <- unlist(lapply(rev(ct.levels), function(ct) {
    paste(ct, split.levels, sep = " ")
  }))

  print(y.levels)
  stats.df$y.label <- factor(paste(stats.df$cell.type, stats.df$split.id, sep = " "),
                              levels = y.levels)

  # Assign a distinct color per split identity, grey-scale color bar within each
  identity.colors <- scales::hue_pal()(length(split.levels))
  names(identity.colors) <- split.levels

  # Build plot layer by layer: one geom_point + color scale per split identity
  # Only show color bar legend for the first group
  p <- ggplot() +
    scale_size_continuous(range = c(1, 8), name = "% Expressed")

  for (i in seq_along(split.levels)) {
    sid <- split.levels[i]
    sub.df <- stats.df[stats.df$split.id == sid, ]
    show.legend <- TRUE 
    p <- p +
      geom_point(data = sub.df,
                 aes(x = gene, y = y.label, size = pct.exp, color = avg.exp),
                 show.legend = show.legend) +
      scale_color_gradient(low = "lightgrey", high = identity.colors[sid],
                           name = "Avg Expression",
                           guide = if (show.legend) "colorbar" else "none") +
      ggnewscale::new_scale_color()
  }

  p + scale_y_discrete(limits = y.levels) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Gene", y = "")
}
