## ABOUTME: Shiny module for cell cycle phase scoring using Seurat's CellCycleScoring.
## ABOUTME: Supports mouse and human gene lists, displays S and G2M phase scores on UMAP.

mod.cellcycle.ui <- function(id) {
  ns <- NS(id)
  
  tabPanel("Cell Cycle", value = "tab.cellcycle",
           wellPanel(
             strong("Cell Cycle Scoring"),
             fluidRow(
               column(5, radioButtons(ns("cell.cycle.ref"), "Species",
                                      choices = c("Mouse", "Human"), inline = TRUE)),
               column(2, actionButton(ns("cell.cycle.run"), "Run Cell Cycle", class = "btn-success"))
             )
           ),
           wellPanel(
             strong("Cell Cycle Plot"),
             fluidRow(
               column(12, plotOutput(ns("plot.cell.cycle"), height = "500px", width = "100%"))
             )
           ),
           wellPanel(
             mod.save.config.ui(ns("save"), label = "Download Seurat Object")
           )
  )
}

mod.cellcycle.server <- function(id, seurat.obj.doublet) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    completed <- reactiveVal(FALSE)

    # 1. Process Cell Cycle Scoring
    data.cell.cycle <- eventReactive(input$cell.cycle.run, {
      # req() ensures the object from the previous module exists
      req(seurat.obj.doublet())
      
      withProgress(message = 'Calculating Cell Cycle Scores...', value = 0.5, {
        srt <- seurat.obj.doublet()
        
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
      p.combined <- p1 | p2
      
      # 3. Print the combined plot
      p.combined
    }, res = 96)
    
    # 3. Export with save configuration
    mod.save.config.server("save", data.cell.cycle)

    return(list(seurat.obj = data.cell.cycle, completed = completed))
  })
}