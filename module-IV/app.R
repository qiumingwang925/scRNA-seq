## ABOUTME: Entry point for Module IV (cell-cell communication via CellChat).
## ABOUTME: Loads packages, defines UI layout, and wires module servers.

source("../R/utils.R")

load.or.install("shiny")
load.or.install("shinyjs")
load.or.install("Seurat")
load.or.install("SeuratObject")
load.or.install("ggplot2")
load.or.install("tidyverse")
load.or.install("patchwork")
load.or.install("future")
load.or.install("ComplexHeatmap")
load.or.install("CellChat", github.url = "jinworks/CellChat")

options(shiny.maxRequestSize = 5 * 1024^3)

source("modules/mod_interact_upload.R")
source("modules/mod_interact_cellchat_comp.R")

ui <- fluidPage(
  useShinyjs(),
  titlePanel("scNexus-Interact"),
  tabsetPanel(
    mod.interact.upload.ui("upload"),
    mod.interact.cellchat.comp.ui("comp")
  )
)

server <- function(input, output, session) {
  uploaded.seurat <- mod.interact.upload.server("upload")
  cellchat.result <- mod.interact.cellchat.comp.server("comp", shared.data = uploaded.seurat)
}

shinyApp(ui, server)
