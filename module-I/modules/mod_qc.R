mod_qc_ui <- function(id) {
  ns <- NS(id)
  tabPanel("QC Removal",
           wellPanel(
             strong("Plot Settings"),
             fluidRow(
               column(3, selectInput(ns("qc.matric.1"), "QC Matrix 1",
                                     choices = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp", "percent.hb"))),
               column(3, selectInput(ns("qc.matric.2"), "QC Matrix 2",
                                     choices = c("nCount_RNA","nFeature_RNA", "percent.mt", "percent.rp", "percent.hb"))),
               column(3, selectInput(ns("qc.plot.type"), "Plot Type",
                                     choices = c("Scatter","Violin", "Density")))
             )
           ),
           wellPanel(
             strong("QC Plot"),
             plotOutput(ns("plot.qc"), height = "500px", width = "700px")
           ),
           wellPanel(
             strong("Filtering Controls"),
             fluidRow(
               column(4, sliderInput(ns("nfeature"), "nFeature Range:", min = 0, max = 10000, value = c(0, 5000))),
               column(4, sliderInput(ns("ncount"), "nCount Range:", min = 0, max = 100000, value = c(0, 60000))),
               column(4, sliderInput(ns("mt"), "Mitochondrial %:", min = 0, max = 100, value = 10))
             ),
             fluidRow(
               column(4, sliderInput(ns("rp"), "Ribosomal %:", min = 0, max = 100, value = 50)),
               column(4, sliderInput(ns("hb"), "Hemoglobin %:", min = 0, max = 100, value = 1))
             ),
             fluidRow(
               column(2, actionButton(ns("qc.plot.run"), "Plot", class = "btn-success", style = "width: 100%")),
               column(2, actionButton(ns("qc.filter.run"), "Filter Low Quality Cells", class = "btn-success", style = "width: 100%")),
               column(3, textOutput(ns("qc.cell.count")))
             )
           ),
           wellPanel(
             downloadButton(ns("download.qc"), "Download Seurat Object")
           )
  )
}


mod_qc_server <- function(id, seurat_obj) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    completed <- reactiveVal(FALSE)

    # 1. Dynamic Slider Updates
    observe({
      req(seurat_obj())
      df <- seurat_obj()@meta.data
      
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
    
    # 2. Refined QC Classification
    # We trigger this on Plot Run so the user can see the "Fail" red dots
    data.qc <- eventReactive(input$qc.plot.run, {
      req(seurat_obj())
      srt <- seurat_obj()
      
      # Logic: Start all as Pass, then flag Fails
      meta <- srt@meta.data %>%
        mutate(QC = case_when(
          percent.mt > input$mt ~ "Fail",
          percent.rp > input$rp ~ "Fail",
          percent.hb > input$hb ~ "Fail",
          nFeature_RNA < input$nfeature[1] | nFeature_RNA > input$nfeature[2] ~ "Fail",
          nCount_RNA < input$ncount[1] | nCount_RNA > input$ncount[2] ~ "Fail",
          TRUE ~ "Pass"
        ))
      
      srt@meta.data$QC <- meta$QC
      return(srt)
    })
    
    # 3. Plotting with Updated Mapping
    plotInput.qc <- eventReactive(input$qc.plot.run, {
      req(data.qc())
      srt <- data.qc()
      
      # Retrieve intercepts for current selection
      get_int <- function(m) {
        switch(m, "nFeature_RNA" = input$nfeature, "nCount_RNA" = input$ncount,
               "percent.mt" = input$mt, "percent.rp" = input$rp, "percent.hb" = input$hb)
      }
      
      int1 <- get_int(input$qc.matric.1)
      int2 <- get_int(input$qc.matric.2)
      
      if(input$qc.plot.type == "Scatter"){
        FeatureScatter(srt, feature1 = input$qc.matric.1, feature2 = input$qc.matric.2, group.by = "QC",
                       cols = c("Fail" = "#F8766D", "Pass" = "grey")) +
          geom_hline(yintercept = int2, linetype = "dashed", color = "blue") +
          geom_vline(xintercept = int1, linetype = "dashed", color = "blue")
        
      } else if(input$qc.plot.type == "Violin"){
        p1 <- VlnPlot(srt, features = input$qc.matric.1, group.by = "QC") + geom_hline(yintercept = int1, linetype = "dashed")
        p2 <- VlnPlot(srt, features = input$qc.matric.2, group.by = "QC") + geom_hline(yintercept = int2, linetype = "dashed")
        ggpubr::ggarrange(p1, p2, ncol=2)
        
      } else {
        # Density Plots using .data[[]] instead of aes_string
        p1 <- ggplot(srt@meta.data, aes(x = .data[[input$qc.matric.1]], fill = QC)) + 
          geom_density(alpha = 0.5) + theme_bw() + geom_vline(xintercept = int1, linetype = "dashed")
        p2 <- ggplot(srt@meta.data, aes(x = .data[[input$qc.matric.2]], fill = QC)) + 
          geom_density(alpha = 0.5) + theme_bw() + geom_vline(xintercept = int2, linetype = "dashed")
        ggpubr::ggarrange(p1, p2, ncol=1)
      }
    })
    
    output$plot.qc <- renderPlot({ plotInput.qc() }, res = 96)

    # Rename button to "Update Plot" after first plot render
    observeEvent(input$qc.plot.run, {
      updateActionButton(session, "qc.plot.run", label = "Update Plot")
    }, once = TRUE)
    
    # 4. Final Filtering
    data.qc.filter <- eventReactive(input$qc.filter.run, {
      req(data.qc())
      completed(TRUE)
      subset(data.qc(), subset = QC == 'Pass')
    })
    
    # Display count
    output$qc.cell.count <- renderText({
      req(data.qc.filter())
      paste0("Cells remaining: ", ncol(data.qc.filter()))
    })
    
    return(list(seurat_obj = data.qc.filter, completed = completed))
  })
}