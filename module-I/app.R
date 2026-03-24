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

  titlePanel("scRNA-seq Analysis Platform"),

  tabsetPanel(
    
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
  
  seurat_obj_import <- mod_import_server("import")
  
  seurat_obj_qc <- mod_qc_server("qc", seurat_obj_import)
  
  seurat_obj_pca <- mod_pca_server("pca", seurat_obj_qc)
  
  seurat_obj_doublet <- mod_doublet_server("doublet", seurat_obj_pca, reactive({ 20 }))
  
  seurat_obj_cellcycle <- mod_cellcycle_server("cellcycle", seurat_obj_doublet)
  
  #seurat_obj_annotation_singler <- mod_annotation_singler_server("annotation_singler", seurat_obj_cellcycle)
  
  # mod_biomarker_server("biomarker",seurat_obj_cellcycle)
  
  #mod_annotation_manual_server("annotation_manual", seurat_obj_cellcycle )
  
}

shinyApp(ui, server)