library(Seurat) # Seurat analysis and visualization
library(presto) # speed up differential expression
library(plotly) # UMAP plotly

options(shiny.maxRequestSize = 5 * 1024^3)

source("~/Documents/Research/UCDacis_Venosa/Shiny/updated_module_III/modules/mod_explore_upload.R")
source("~/Documents/Research/UCDacis_Venosa/Shiny/updated_module_III/modules/mod_explore_umap.R")
source("~/Documents/Research/UCDacis_Venosa/Shiny/updated_module_III/modules/mod_explore_violin.R")
source("~/Documents/Research/UCDacis_Venosa/Shiny/updated_module_III/modules/mod_explore_dot.R")
source("~/Documents/Research/UCDacis_Venosa/Shiny/updated_module_III/modules/mod_explore_heatmap.R")
source("~/Documents/Research/UCDacis_Venosa/Shiny/updated_module_III/modules/mod_explore_de.R")
source("~/Documents/Research/UCDacis_Venosa/Shiny/updated_module_III/modules/mod_explore_enrich.R")

ui <- fluidPage(
  titlePanel("scNexus-Explore"),
  tabsetPanel(
    mod_explore_upload_ui("upload"),
    mod_explore_umap_ui("umap"),
    #mod_explore_violin_ui("violin"),
    #mod_explore_dot_ui("dot"),
    #mod_explore_heatmap_ui("heatmap"),
    #mod_explore_de_ui("de"), 
    #mod_explore_enrich_ui("enrich")
  )
)

server <- function(input, output, session){
  uploaded_seurat <- mod_explore_upload_server("upload")
  mod_explore_umap_server("umap", shared_data = uploaded_seurat)
  #mod_explore_violin_server("violin", shared_data =uploaded_seurat )
  #mod_explore_dot_server("dot", shared_data = uploaded_seurat )
  #mod_explore_heatmap_server("heatmap", shared_data =uploaded_seurat )
  #mod_explore_de_server("de", shared_data =uploaded_seurat )
  #mod_explore_enrich_server("enrich", shared_data =uploaded_seurat )
}

shinyApp(ui, server)