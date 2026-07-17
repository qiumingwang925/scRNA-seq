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
load.or.install("harmony")
load.or.install("batchelor")  #FastMNN backend for SeuratWrappers
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
    mod.upload.merge.ui("upload.merge"),
    mod.integrate.ui("integrate"),
    mod.benchmark.ui("benchmark"),
    mod.annotation.ui("annotation")
  )
)

server <- function(input, output, session) {
  uploaded.seurat <- mod.upload.merge.server("upload.merge")
  processed.seurat <- mod.integrate.server("integrate", shared.data = uploaded.seurat)
  mod.benchmark.server("benchmark", shared.data = processed.seurat)
  annotated.seurat <- mod.annotation.server("annotation", shared.data = processed.seurat)
}

shinyApp(ui, server)
