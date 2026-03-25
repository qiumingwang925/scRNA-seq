mod_cellcycle_ui <- function(id) {
  ns <- NS(id)
  
  tabPanel("Cell Cycle",
           fluidRow(
             # Added ns() to inputs
             column(5, radioButtons(ns("cell.cycle.ref"), "Species", 
                                    choices = c("Mouse", "Human"), inline = TRUE)),
             column(2, actionButton(ns("cell.cycle.run"), "Run Cell Cycle", class = "btn-success"))
           ),
           hr(),
           fluidRow(
             # Added ns() to output
             column(12, plotOutput(ns("plot.cell.cycle"), height = "500px", width = "100%"))
           ),
           br(),
           fluidRow(
             column(2, downloadButton(ns("download.cellcycle"), "Download Seurat Object"))
           )
  )
}

mod_cellcycle_server <- function(id, seurat_obj_doublet) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    completed <- reactiveVal(FALSE)

    # 1. Process Cell Cycle Scoring
    data.cell.cycle <- eventReactive(input$cell.cycle.run, {
      # req() ensures the object from the previous module exists
      req(seurat_obj_doublet())
      
      withProgress(message = 'Calculating Cell Cycle Scores...', value = 0.5, {
        srt <- seurat_obj_doublet()
        
        # Determine gene lists based on species
        s.genes <- cc.genes.updated.2019$s.genes
        g2m.genes <- cc.genes.updated.2019$g2m.genes
        
        if (input$cell.cycle.ref == "Mouse") {
          # Convert HUMAN genes to Mouse format (Title Case)
          s.genes <- stringr::str_to_title(s.genes)
          g2m.genes <- stringr::str_to_title(g2m.genes)
        }
        
        # Run Seurat scoring
        srt <- CellCycleScoring(srt,
                                s.features = s.genes,
                                g2m.features = g2m.genes,
                                set.ident = TRUE)
      })
      completed(TRUE)
      return(srt)
    })
    
    # 2. Render Plot
    ### UMAP figure ####
    output$plot.cell.cycle <- renderPlot({
      # Ensure the data exists
      req(data.cell.cycle())
      
      # 1. Generate individual plots
      p1 <- FeaturePlot(data.cell.cycle(), 
                        features = "S.Score", 
                        reduction = "umap", 
                        raster = FALSE) + 
        ggtitle("S Phase Score") +
        theme(aspect.ratio = 1) # Keeps UMAPs square
      
      p2 <- FeaturePlot(data.cell.cycle(), 
                        features = "G2M.Score", 
                        reduction = "umap", 
                        raster = FALSE) + 
        ggtitle("G2M Phase Score") +
        theme(aspect.ratio = 1)
      
      # 2. Use patchwork to align horizontally
      # The '+' sign puts them side-by-side by default
      p_combined <- p1 | p2
      
      # 3. Print the combined plot
      p_combined
    }, res = 96)
    
    # 3. Download Handler
    output$download.cellcycle <- downloadHandler(
      filename = function() {
        paste0("seurat_cellcycle_", Sys.Date(), ".rds")
      },
      content = function(file) {
        saveRDS(data.cell.cycle(), file = file)
      }
    )
    
    return(list(seurat_obj = data.cell.cycle, completed = completed))
  })
}