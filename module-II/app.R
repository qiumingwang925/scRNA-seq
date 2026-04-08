library(shiny)
#library(shinyjs)
library(Seurat)
library(SeuratObject)
library(SeuratWrappers)
library(ggplot2)
library(tidyverse)
#library(openxlsx)
#library(shinyFiles)
library(plotly)
library(DT)
library(glmGamPoi) #SCTransform
library(ggpubr)
library(presto)

options(shiny.maxRequestSize = 10 * 1024^3)

source("~/Documents/Research/UCDacis_Venosa/Shiny/updated_module_II/modules/mod_upload_merge.R")
source("~/Documents/Research/UCDacis_Venosa/Shiny/updated_module_II/modules/mod_integrate.R")
source("~/Documents/Research/UCDacis_Venosa/Shiny/updated_module_II/modules/mod_benchmark.R")
source("~/Documents/Research/UCDacis_Venosa/Shiny/updated_module_II/modules/mod_annotation.R")
#source("~/Documents/Research/UCDacis_Venosa/Shiny/updated_module_II/modules/mod_download.R")

ui <- fluidPage(
  titlePanel("scNexus-Integrate"),
  tabsetPanel(
    mod_upload_merge_ui("upload_merge"),
    mod_integrate_ui("integrate"),
    mod_benchmark_ui("benchmark"),
    mod_annotation_ui("annotation"),
    #mod_download_ui("download")
  )
)

server <- function(input, output, session){
  uploaded_seurat <- mod_upload_merge_server("upload_merge")
  processed_seurat <- mod_integrate_server("integrate", shared_data = uploaded_seurat)
  mod_benchmark_server("benchmark", shared_data = processed_seurat)
  annotated_seurat <- mod_annotation_server("annotation", shared_data = processed_seurat )
  
}

shinyApp(ui, server)