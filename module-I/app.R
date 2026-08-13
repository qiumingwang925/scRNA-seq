## ABOUTME: Defines the root Shiny app and loads the analysis modules for the workspace.
## ABOUTME: Ensures required R packages are available before creating the UI and server.

source("../R/utils.R")

load.or.install("shiny")
load.or.install("Seurat")
load.or.install("ggplot2")
load.or.install("tidyverse")
#load.or.install("openxlsx")
load.or.install("shinyFiles")
load.or.install("DoubletFinder", github.url = "chris-mcginnis-ucsf/DoubletFinder")
load.or.install("plotly")
load.or.install("glmGamPoi") #SCTransform
load.or.install("ggpubr")
load.or.install("shinycssloaders")
load.or.install("SingleR")
load.or.install("SingleCellExperiment")
load.or.install("DT")

options(shiny.maxRequestSize = 10 * 1024^3)
options(future.globals.maxSize = 10 * 1024^3) 

# biomarker database
#markers <- read.xlsx("data/biomarkers_mouse.xlsx")

# directory for rds files
#seurat.data.dir <- "Data/"

source("modules/mod_import.R")
source("modules/mod_qc.R")
source("modules/mod_pca.R")
source("modules/mod_doublet.R")
source("modules/mod_cellcycle.R")
source("../R/mod_save_config.R")
source("modules/mod_annotation_singler.R")
source("modules/mod_annotation_manual.R")
source("modules/mod_annotation.R")


ui <- fluidPage(
  tags$head(tags$style(HTML("
    .well {
      background-color: #ffffff !important;
    }
    .plot\\.params\\.row {
      display: flex;
      align-items: stretch;
    }
    .plot\\.params\\.row > [class*='col-'] {
      display: flex;
    }
    .plot\\.params\\.row .well {
      flex: 1;
      display: flex;
      flex-direction: column;
      width: 100%;
    }
    .square\\.plot {
      width: 100%;
      max-height: 70vh;
      display: flex;
      justify-content: center;
    }
    .square\\.plot .shiny-plot-output {
      aspect-ratio: 1 / 1;
      height: 70vh !important;
      width: auto !important;
    }
    @media (max-width: 1200px) {
      .plot\\.params\\.row {
        flex-direction: column;
      }
      .plot\\.params\\.row > [class*='col-'] {
        width: 100%;
      }
    }
  "))),

  titlePanel("Module I - scNexus-Process"),

  tabsetPanel(id = "main_tabs",
    
    mod.import.ui("import"),
    
    mod.qc.ui("qc"),
    
    mod.pca.ui("pca"),
    
    mod.doublet.ui("doublet"),
    
    mod.cellcycle.ui("cellcycle"),

    mod.annotation.ui("annotation")

  )
)

server <- function(input, output, session){
  
  seurat.obj <- reactiveVal(NULL)

  import.result <- mod.import.server("import")
  qc.result <- mod.qc.server("qc", import.result$seurat.obj, upstream.completed = import.result$completed)
  pca.result <- mod.pca.server("pca", qc.result$seurat.obj, upstream.completed = qc.result$completed)
  doublet.result <- mod.doublet.server("doublet", pca.result$seurat.obj, pca.result$pca.dims, upstream.completed = pca.result$completed)
  cellcycle.result <- mod.cellcycle.server("cellcycle", doublet.result$seurat.obj, upstream.completed = doublet.result$completed)
  annotation.result <- mod.annotation.server("annotation", cellcycle.result$seurat.obj, upstream.completed = cellcycle.result$completed)

}

shinyApp(ui, server)