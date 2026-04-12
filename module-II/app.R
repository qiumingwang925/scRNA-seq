## ABOUTME: Entry point for Module II (multi-sample integration and annotation).
## ABOUTME: Loads packages, defines UI layout, and wires module servers.

source("../R/utils.R")

load.or.install("shiny")
load.or.install("shinyjs")
load.or.install("Seurat")
load.or.install("SeuratObject")
load.or.install("SeuratWrappers", github.url = "satijalab/seurat-wrappers")
load.or.install("ggplot2")
load.or.install("tidyverse")
load.or.install("plotly")
load.or.install("DT")
load.or.install("glmGamPoi") #SCTransform
load.or.install("ggpubr")
load.or.install("presto", github.url = "immunogenomics/presto")
load.or.install("cluster")  #ASW
load.or.install("lisi", github.url = "immunogenomics/lisi") #LISI
load.or.install("Matrix")   #Graph LISI

options(shiny.maxRequestSize = 10 * 1024^3)

source("../R/mod_save_config.R")
source("modules/mod_upload_merge.R")
source("modules/mod_integrate.R")
source("modules/mod_benchmark.R")
source("modules/mod_annotation.R")

ui <- fluidPage(
  titlePanel("scNexus-Integrate"),
  tabsetPanel(
    mod_upload_merge_ui("upload_merge"),
    mod_integrate_ui("integrate"),
    mod_benchmark_ui("benchmark"),
    mod_annotation_ui("annotation")
  )
)

server <- function(input, output, session) {
  uploaded_seurat <- mod_upload_merge_server("upload_merge")
  processed_seurat <- mod_integrate_server("integrate", shared_data = uploaded_seurat)
  mod_benchmark_server("benchmark", shared_data = processed_seurat)
  annotated_seurat <- mod_annotation_server("annotation", shared_data = processed_seurat)
}

shinyApp(ui, server)
