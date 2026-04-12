## ABOUTME: Shiny module for quality control filtering of scRNA-seq data.
## ABOUTME: Provides interactive threshold sliders and visual pass/fail classification plots.

mod.qc.ui <- function(id) {
  ns <- NS(id)
  tabPanel("QC Removal", value = "tab.qc",
           wellPanel(
             strong("Plot Settings"),
             fluidRow(
               column(3, selectInput(ns("qc.metric.1"), "QC Matrix 1",
                                     choices = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp", "percent.hb"))),
               column(3, selectInput(ns("qc.metric.2"), "QC Matrix 2",
                                     choices = c("nCount_RNA","nFeature_RNA", "percent.mt", "percent.rp", "percent.hb"))),
               column(3, selectInput(ns("qc.plot.type"), "Plot Type",
                                     choices = c("Scatter","Violin", "Density")))
             )
           ),
           fluidRow(class = "plot.params.row",
             column(8,
               wellPanel(
                 fluidRow(
                   column(6, strong("QC Plot"), textOutput(ns("qc.selected.count"), inline = TRUE)),
                   column(3, offset = 3, actionButton(ns("qc.plot.run"), "Plot", class = "btn-success", style = "width: 100%"))
                 ),
                 tags$div(class = "square.plot",
                   plotOutput(ns("plot.qc"), width = "100%", height = "100%")
                 )
               )
             ),
             column(4,
               wellPanel(
                 strong("Filtering Controls"),
                 sliderInput(ns("nfeature"), "nFeature Range:", min = 0, max = 10000, value = c(0, 5000)),
                 sliderInput(ns("ncount"), "nCount Range:", min = 0, max = 100000, value = c(0, 60000)),
                 sliderInput(ns("mt"), "Mitochondrial %:", min = 0, max = 100, value = 10),
                 sliderInput(ns("rp"), "Ribosomal %:", min = 0, max = 100, value = 50),
                 sliderInput(ns("hb"), "Hemoglobin %:", min = 0, max = 100, value = 1),
                 actionButton(ns("qc.filter.run"), "Filter Low Quality Cells", class = "btn-success", style = "width: 100%")
               )
             )
           )
  )
}


mod.qc.server <- function(id, seurat.obj) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    completed <- reactiveVal(FALSE)

    # Filter button starts disabled until user visualizes the data
    shinyjs::disable("qc.filter.run")

    # Disable filter button when any slider changes (plot is stale)
    observe({
      input$nfeature; input$ncount; input$mt; input$rp; input$hb
      shinyjs::disable("qc.filter.run")
    })

    # Enable filter button after plot is rendered
    observeEvent(input$qc.plot.run, {
      shinyjs::enable("qc.filter.run")
    })

    # 1. Dynamic Slider Updates
    observe({
      req(seurat.obj())
      df <- seurat.obj()@meta.data
      
      updateSliderInput(session, "nfeature",
                        value = c(min(df$nFeature_RNA), max(df$nFeature_RNA)),
                        min = min(df$nFeature_RNA), max = max(df$nFeature_RNA))
      
      updateSliderInput(session, "ncount",
                        value = c(min(df$nCount_RNA), max(df$nCount_RNA)),
                        min = min(df$nCount_RNA), max = max(df$nCount_RNA))
      
      # Ensure max isn't 0 to avoid slider errors
      updateSliderInput(session, "mt", max = ceiling(max(df$percent.mt, na.rm = TRUE)))
      updateSliderInput(session, "rp", max = ceiling(max(df$percent.rp, na.rm = TRUE)))
      updateSliderInput(session, "hb", max = ceiling(max(df$percent.hb, na.rm = TRUE)))
    })
    
    # 2. QC Classification — triggered on Plot so user sees results before filtering
    data.qc <- eventReactive(input$qc.plot.run, {
      req(seurat.obj())
      srt <- seurat.obj()

      meta <- srt@meta.data %>%
        mutate(QC = case_when(
          percent.mt > input$mt ~ "To Be Filtered",
          percent.rp > input$rp ~ "To Be Filtered",
          percent.hb > input$hb ~ "To Be Filtered",
          nFeature_RNA < input$nfeature[1] | nFeature_RNA > input$nfeature[2] ~ "To Be Filtered",
          nCount_RNA < input$ncount[1] | nCount_RNA > input$ncount[2] ~ "To Be Filtered",
          TRUE ~ "Selected"
        ))
      
      srt@meta.data$QC <- meta$QC
      return(srt)
    })
    
    # Selected cell count next to QC Plot title
    output$qc.selected.count <- renderText({
      req(data.qc())
      n.selected <- sum(data.qc()@meta.data$QC == "Selected")
      paste0("(", n.selected, " samples selected)")
    })

    # 3. Plot rendering
    # Threshold lines update live with slider changes; point colors update on Plot click

    # Returns threshold value(s) for a given QC metric.
    # Range sliders (nFeature, nCount) return a 2-element vector — ggplot vectorizes
    # geom_vline/geom_hline over it, drawing one line per bound.
    # Percentage sliders (mt, rp, hb) return a single value.
    get.int <- function(m) {
      switch(m, "nFeature_RNA" = input$nfeature, "nCount_RNA" = input$ncount,
             "percent.mt" = input$mt, "percent.rp" = input$rp, "percent.hb" = input$hb)
    }

    output$plot.qc <- renderPlot({
      req(data.qc())
      srt <- data.qc()

      int1 <- get.int(input$qc.metric.1)
      int2 <- get.int(input$qc.metric.2)

      if(input$qc.plot.type == "Scatter"){
        # Draw "To Be Filtered" first so "Selected" points render on top
        df <- srt@meta.data
        df$QC <- factor(df$QC, levels = c("To Be Filtered", "Selected"))
        df <- df[order(df$QC), ]
        ggplot(df, aes(x = .data[[input$qc.metric.1]], y = .data[[input$qc.metric.2]], color = QC)) +
          geom_point(size = 1) +
          scale_color_manual(values = c("To Be Filtered" = "#d3d3d3", "Selected" = "black")) +
          geom_hline(yintercept = int2, linetype = "dashed", color = "blue") +
          geom_vline(xintercept = int1, linetype = "dashed", color = "blue") +
          theme_bw()

      } else if(input$qc.plot.type == "Violin"){
        p1 <- VlnPlot(srt, features = input$qc.metric.1, group.by = "QC") + geom_hline(yintercept = int1, linetype = "dashed")
        p2 <- VlnPlot(srt, features = input$qc.metric.2, group.by = "QC") + geom_hline(yintercept = int2, linetype = "dashed")
        ggpubr::ggarrange(p1, p2, ncol=2)

      } else {
        p1 <- ggplot(srt@meta.data, aes(x = .data[[input$qc.metric.1]], fill = QC)) +
          geom_density(alpha = 0.5) + theme_bw() + geom_vline(xintercept = int1, linetype = "dashed")
        p2 <- ggplot(srt@meta.data, aes(x = .data[[input$qc.metric.2]], fill = QC)) +
          geom_density(alpha = 0.5) + theme_bw() + geom_vline(xintercept = int2, linetype = "dashed")
        ggpubr::ggarrange(p1, p2, ncol=1)
      }
    }, res = 96)

    # Rename button to "Update Plot" after first plot render
    observeEvent(input$qc.plot.run, {
      updateActionButton(session, "qc.plot.run", label = "Update Selection")
    }, once = TRUE)
    
    # 4. Final Filtering
    observeEvent(input$qc.filter.run, {
      req(data.qc())
      completed(TRUE)
    })

    data.qc.filter <- eventReactive(input$qc.filter.run, {
      req(data.qc())
      subset(data.qc(), subset = QC == 'Selected')
    })
    
    return(list(seurat.obj = data.qc.filter, completed = completed))
  })
}