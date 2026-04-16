## ABOUTME: Entry point for Module III (post-annotation exploration and visualization).
## ABOUTME: Loads packages, defines UI layout, and wires module servers.

source("../R/utils.R")

load.or.install("shiny")
load.or.install("shinyjs")
load.or.install("Seurat")
load.or.install("SeuratObject")
load.or.install("ggplot2")
load.or.install("tidyverse")
load.or.install("plotly")
load.or.install("scales")
load.or.install("patchwork")
load.or.install("ggnewscale")
load.or.install("DT")
load.or.install("presto", github.url = "immunogenomics/presto")
load.or.install("enrichR")
load.or.install("openxlsx")

options(shiny.maxRequestSize = 5 * 1024^3)

source("modules/mod_explore_upload.R")
source("modules/mod_explore_umap.R")
source("modules/mod_explore_violin.R")
source("modules/mod_explore_dot.R")
source("modules/mod_explore_heatmap.R")
source("modules/mod_explore_de.R")
source("modules/mod_explore_enrich.R")

ui <- fluidPage(
  useShinyjs(),
  titlePanel("scNexus-Explore"),
  tabsetPanel(
    mod.explore.upload.ui("upload"),
    mod.explore.umap.ui("umap"),
    mod.explore.violin.ui("violin"),
    mod.explore.dot.ui("dot"),
    mod.explore.heatmap.ui("heatmap"),
    mod.explore.de.ui("de"),
    mod.explore.enrich.ui("enrich")
  )
)

server <- function(input, output, session) {
  uploaded.seurat <- mod.explore.upload.server("upload")
  mod.explore.umap.server("umap", shared.data = uploaded.seurat)
  mod.explore.violin.server("violin", shared.data = uploaded.seurat)
  mod.explore.dot.server("dot", shared.data = uploaded.seurat)
  mod.explore.heatmap.server("heatmap", shared.data = uploaded.seurat)
  de.result <- mod.explore.de.server("de", shared.data = uploaded.seurat)
  mod.explore.enrich.server("enrich", de.result = de.result)
}

shinyApp(ui, server)
