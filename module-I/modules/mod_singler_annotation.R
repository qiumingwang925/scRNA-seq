mod_annotation_singler_ui <- function(id){
  
  ns <- NS(id)
  
  tabPanel( "SingleR Annotation",
            
            fluidRow(
              column(5, selectInput("auto.annotation", "Reference", choices = c("BD1", "BD2", "DB3", "DB4"))),
              column(2, actionButton("auto.annotation.run", "Run SingleR Annotation", class ="btn-success"))
            ),
            
            fluidRow(
              column(10, plotOutput("plot.auto.annotation",  height = "450px", width = "500px"))
            ),
            
            fluidRow(
              column(2, downloadButton("download.annotation", "Download Seurat Object"))
            )
            
  )
}

mod_annotation_singler_server <- function(id, seurat_obj_cellcycle){
  
  moduleServer(id, function(input, output, session){
    
  })
}  