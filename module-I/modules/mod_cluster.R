mod_pca_server <- function(id, seurat_obj){
  
  moduleServer(id, function(input, output, session){
    
    observeEvent(input$run,{
      
      srt <- seurat_obj()
      
      srt <- NormalizeData(srt)
      srt <- FindVariableFeatures(srt)
      srt <- ScaleData(srt)
      srt <- RunPCA(srt)
      
      seurat_obj(srt)
      
    })
    
    output$pca <- renderPlot({
      ElbowPlot(seurat_obj())
    })
    
  })
  
}