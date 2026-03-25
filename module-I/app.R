## ABOUTME: Defines the root Shiny app and loads the analysis modules for the workspace.
## ABOUTME: Ensures required R packages are available before creating the UI and server.

cran_repo <- unname(getOption("repos")["CRAN"])
if (is.null(cran_repo) || length(cran_repo) == 0 || is.na(cran_repo) || cran_repo == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

load_or_install <- function(pkg, github_url = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!is.null(github_url)) {
      if (!requireNamespace("remotes", quietly = TRUE))
        install.packages("remotes")
      remotes::install_github(github_url)
    } else {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(pkg)
    }
  }
  library(pkg, character.only = TRUE)
}

load_or_install("shiny")
load_or_install("shinyjs")
load_or_install("Seurat")
load_or_install("ggplot2")
load_or_install("tidyverse")
#load_or_install("openxlsx")
load_or_install("shinyFiles")
load_or_install("DoubletFinder", github_url = "chris-mcginnis-ucsf/DoubletFinder")
load_or_install("plotly")
load_or_install("glmGamPoi") #SCTransform
load_or_install("ggpubr")
load_or_install("shinycssloaders")

options(shiny.maxRequestSize = 10 * 1024^3)

# Set to TRUE to enable all tabs and buttons for UI testing
UI_TESTING <- FALSE

# biomarker database
#markers <- read.xlsx("data/biomarkers_mouse.xlsx")

# directory for rds files
#seurat_data_dir <- "Data/"

source("modules/mod_import.R")
source("modules/mod_qc.R")
source("modules/mod_pca.R")
source("modules/mod_doublet.R")
source("modules/mod_cellcycle.R")
#source("modules/mod_annotation_singler.R")
source("modules/mod_biomarker.R")
#source("modules/mod_annotation_manual.R")


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
  "))),

  titlePanel("scRNA-seq Analysis Platform"),

  tabsetPanel(id = "main_tabs",
    
    mod_import_ui("import"),
    
    mod_qc_ui("qc"),
    
    mod_pca_ui("pca"),
    
    mod_doublet_ui("doublet"),
    
    mod_cellcycle_ui("cellcycle"),
    
    #mod_annotation_singler_ui("annotation_singler")
    
    mod_biomarker_ui("biomarker")
    
    #mod_annotation_manual_ui("annotation_manual")
    
  )
)

server <- function(input, output, session){
  
  seurat_obj <- reactiveVal(NULL)

  # Helper to enable/disable a tab by index
  enable_tab <- function(i) {
    shinyjs::runjs(sprintf("$('#main_tabs > li:nth-child(%d)').removeClass('disabled');", i))
  }
  disable_tab <- function(i) {
    shinyjs::runjs(sprintf("$('#main_tabs > li:nth-child(%d)').addClass('disabled');", i))
  }
  disable_tabs <- function(indices) { for (i in indices) disable_tab(i) }

  # Disable all downstream tabs on startup (unless UI testing)
  if (!UI_TESTING) disable_tabs(2:6)

  import_result <- mod_import_server("import", ui_testing = UI_TESTING)
  qc_result <- mod_qc_server("qc", import_result$seurat_obj)
  pca_result <- mod_pca_server("pca", qc_result$seurat_obj)
  doublet_result <- mod_doublet_server("doublet", pca_result$seurat_obj, reactive({ 20 }), ui_testing = UI_TESTING)
  cellcycle_result <- mod_cellcycle_server("cellcycle", doublet_result$seurat_obj)

  #seurat_obj_annotation_singler <- mod_annotation_singler_server("annotation_singler", cellcycle_result$seurat_obj)
  # mod_biomarker_server("biomarker", cellcycle_result$seurat_obj)
  #mod_annotation_manual_server("annotation_manual", cellcycle_result$seurat_obj)

  # Progressive tab enabling: each step enables only the next tab
  if (!UI_TESTING) {
    # Import complete → enable QC (tab 2), disable 3-6
    observe({
      if (import_result$converted()) {
        enable_tab(2)
      } else {
        disable_tabs(2:6)
      }
    })

    # QC complete → enable PCA (tab 3)
    observe({
      if (qc_result$completed()) { enable_tab(3) } else { disable_tabs(3:6) }
    })

    # PCA complete → enable Doublet (tab 4)
    observe({
      if (pca_result$completed()) { enable_tab(4) } else { disable_tabs(4:6) }
    })

    # Doublet complete → enable Cell Cycle (tab 5)
    observe({
      if (doublet_result$completed()) { enable_tab(5) } else { disable_tabs(5:6) }
    })

    # Cell Cycle complete → enable Biomarker (tab 6)
    observe({
      if (cellcycle_result$completed()) { enable_tab(6) } else { disable_tab(6) }
    })
  }
  
}

shinyApp(ui, server)