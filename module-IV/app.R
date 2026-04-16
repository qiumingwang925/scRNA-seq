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

options(shiny.maxRequestSize = 5 * 1024^3)

source("modules/mod_interact_upload.R")

ui <- fluidPage(
  useShinyjs(),
  titlePanel("scNexus-Interact"),
  tabsetPanel(
    mod.interact.upload.ui("upload")
  )
)

server <- function(input, output, session) {
  uploaded.seurat <- mod.interact.upload.server("upload")
}

shinyApp(ui, server)
