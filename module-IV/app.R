## ABOUTME: Entry point for Module IV (cell-cell communication via CellChat).
## ABOUTME: Loads packages, defines UI layout, and wires module servers.

source("../R/utils.R")

load.or.install("shiny")
load.or.install("shinyjs")
load.or.install("shinycssloaders")
load.or.install("DT")
load.or.install("Seurat")
load.or.install("SeuratObject")
load.or.install("ggplot2")
load.or.install("tidyverse")
load.or.install("patchwork")
load.or.install("ggplotify")
load.or.install("future")
load.or.install("ComplexHeatmap")
load.or.install("NMF")
load.or.install("ggalluvial")
load.or.install("CellChat", github.url = "jinworks/CellChat")

# LIANA sub-module dependencies
load.or.install("digest")
load.or.install("entropy")
load.or.install("biomaRt")
load.or.install("OmnipathR")
load.or.install("liana", github.url = "saezlab/liana")

options(shiny.maxRequestSize = 5 * 1024^3)

source("modules/mod_interact_upload.R")
source("modules/mod_interact_vis_helpers.R")
source("modules/mod_interact_cellchat_comp.R")
source("modules/mod_interact_cellchat_vis.R")
source("modules/mod_interact_liana_comp.R")
source("modules/mod_interact_liana_vis.R")

ui <- fluidPage(
  useShinyjs(),
  titlePanel("scNexus-Interact"),
  tabsetPanel(
    mod.interact.upload.ui("upload"),
    mod.interact.cellchat.comp.ui("cc.comp"),
    mod.interact.cellchat.vis.ui("cc.vis"),
    mod.interact.liana.comp.ui("liana.comp"),
    mod.interact.liana.vis.ui("liana.vis")
  )
)

server <- function(input, output, session) {
  uploaded.seurat <- mod.interact.upload.server("upload")
  cellchat.result <- mod.interact.cellchat.comp.server("cc.comp", shared.data = uploaded.seurat)
  mod.interact.cellchat.vis.server("cc.vis", cellchat.input = cellchat.result)
  liana.result <- mod.interact.liana.comp.server("liana.comp", shared.data = uploaded.seurat)
  mod.interact.liana.vis.server("liana.vis", liana.input = liana.result)
}

shinyApp(ui, server)
