## ABOUTME: Defines the root Shiny app and loads the analysis modules for the workspace.
## ABOUTME: Ensures required R packages are available before creating the UI and server.

cran.repo <- unname(getOption("repos")["CRAN"])
if (is.null(cran.repo) || length(cran.repo) == 0 || is.na(cran.repo) || cran.repo == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

load.or.install <- function(pkg, github.url = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!is.null(github.url)) {
      if (!requireNamespace("remotes", quietly = TRUE))
        install.packages("remotes")
      remotes::install_github(github.url)
    } else {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(pkg)
    }
  }
  library(pkg, character.only = TRUE)
}

load.or.install("shiny")
load.or.install("shinyjs")
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

# Set to TRUE to enable all tabs and buttons for UI testing
UI.TESTING <- FALSE

# biomarker database
#markers <- read.xlsx("data/biomarkers_mouse.xlsx")

# directory for rds files
#seurat.data.dir <- "Data/"

source("modules/mod_import.R")
source("modules/mod_qc.R")
source("modules/mod_pca.R")
source("modules/mod_doublet.R")
source("modules/mod_cellcycle.R")
source("modules/mod_save_config.R")
source("modules/mod_annotation_singler.R")
source("modules/mod_annotation_manual.R")
source("modules/mod_annotation.R")


ui <- fluidPage(
  useShinyjs(),
  tags$head(tags$style(HTML("
    .nav-tabs > li.disabled > a {
      color: #ccc !important;
      pointer-events: none !important;
      cursor: default !important;
    }
    .btn.disabled, .btn:disabled {
      background-color: #e0e0e0 !important;
      color: #aaa !important;
      border-color: #d0d0d0 !important;
    }
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

  titlePanel("scRNA-seq Analysis Platform"),

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

  # Pipeline tab values in order, used for progressive enable/disable
  pipeline.tabs <- c("tab.qc", "tab.pca", "tab.doublet", "tab.cellcycle", "tab.annotation")

  # Helper to enable/disable a tab by its value attribute
  enable.tab <- function(tab.value) {
    shinyjs::runjs(sprintf("$('#main_tabs > li > a[data-value=\"%s\"]').parent().removeClass('disabled');", tab.value))
  }
  disable.tab <- function(tab.value) {
    shinyjs::runjs(sprintf("$('#main_tabs > li > a[data-value=\"%s\"]').parent().addClass('disabled');", tab.value))
  }
  disable.tabs <- function(tab.values) { for (v in tab.values) disable.tab(v) }

  # Disable all downstream tabs on startup (unless UI testing)
  if (!UI.TESTING) disable.tabs(pipeline.tabs)

  import.result <- mod.import.server("import", ui.testing = UI.TESTING)
  qc.result <- mod.qc.server("qc", import.result$seurat.obj)
  pca.result <- mod.pca.server("pca", qc.result$seurat.obj)
  doublet.result <- mod.doublet.server("doublet", pca.result$seurat.obj, pca.result$pca.dims, ui.testing = UI.TESTING)
  cellcycle.result <- mod.cellcycle.server("cellcycle", doublet.result$seurat.obj)
  annotation.result <- mod.annotation.server("annotation", cellcycle.result$seurat.obj)

  # Progressive tab enabling: each step enables the next tab and disables all after it
  if (!UI.TESTING) {
    observe({
      if (import.result$completed()) {
        enable.tab("tab.qc")
      } else {
        disable.tabs(pipeline.tabs)
      }
    })

    observe({
      if (qc.result$completed()) { enable.tab("tab.pca") } else { disable.tabs(pipeline.tabs[2:5]) }
    })

    observe({
      if (pca.result$completed()) { enable.tab("tab.doublet") } else { disable.tabs(pipeline.tabs[3:5]) }
    })

    observe({
      if (doublet.result$completed()) { enable.tab("tab.cellcycle") } else { disable.tabs(pipeline.tabs[4:5]) }
    })

    observe({
      if (cellcycle.result$completed()) { enable.tab("tab.annotation") } else { disable.tab("tab.annotation") }
    })
  }
  
}

shinyApp(ui, server)