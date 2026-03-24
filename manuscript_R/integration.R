# 1. Load packages ####

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(tidyverse)
  library(Seurat) # toolkit for SC genomics
  library(SeuratWrappers) # fastMNN
  library(patchwork) #make plot composition
  library(ggplot2)
  library(SingleCellExperiment) # S4 data storage
  library(scater) # toolkit for scRNA seq visualization
  library(glmGamPoi) # fit a gamma poisson generalized linear model
  library(cowplot)# ggplot add-on for publication figures
  library(rstudioapi)
  library(ggpubr)
  #library(hciRdata) # mouse gene database ???
  library(presto)
  options(future.globals.maxSize = 1e20)
  set.seed(4321)
})

# 2. Merge samples (counts) ####
## samples number ####
samples <- c("22713X3","23957X4","21401X3", "20354X2", "23957X5", "22713X2", "20354X1", "24143X4")

## upload samples and make a list ####
datasets = list()

for(i in samples) {
  srt<- readRDS(paste0("Manual_Annotated_", i, "_edit.rds"))
  
  srt <- DietSeurat(
    srt,
    assays = "SCT",
    layers = "counts",
    features = NULL,
    dimreducs = NULL,
    graphs = NULL
  )
  srt@commands <- list()
  
  
  datasets <- append(datasets, srt)
  rm(i, srt)
}

## merge data ####
srt <- merge(x = datasets[[1]], y = c(datasets[-1]), add.cell.ids = samples)
rm(samples, datasets)
colnames(srt@meta.data)
srt@meta.data <- srt@meta.data[, c(1:7, 11, 13:22)]
## save merged data ####
saveRDS(srt,"Merged_raw.rds" )

# 3. Log_Transform Integration ####
# Load merged samples
srt <- readRDS("Merged_raw.rds")
# transform and regress data
srt <- NormalizeData(srt)
srt <- FindVariableFeatures(srt)
all.genes <- rownames(srt)
srt <- ScaleData(srt, features= all.genes, vars.to.regress = c("percent.mt", "nCount_RNA", "S.Score", "G2M.Score"))
srt <- ScaleData(srt, features= all.genes)

## unintegrated ####
# PCA, Cluster, and UMAP
srt <- RunPCA(srt)
srt <- FindNeighbors(srt, dims = 1:30)
srt <- RunUMAP(srt, dims = 1:30)

# Integration
## harmony ####
srt <- IntegrateLayers(
  object = srt, method = HarmonyIntegration,  
  orig.reduction = "pca", new.reduction = "integrated.harmony",
  verbose = FALSE
)
srt <- FindNeighbors(srt, reduction = "integrated.harmony", dims = 1:30)
srt <- RunUMAP(srt, reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.harmony")

## fastmnn ####
srt <- IntegrateLayers(
  object = srt, method = FastMNNIntegration, 
  new.reduction = "integrated.mnn",
  verbose = FALSE
)
srt <- FindNeighbors(srt, reduction = "integrated.mnn", dims = 1:30)
srt <- RunUMAP(srt, reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn")

## cca ####
srt <- IntegrateLayers(
  object = srt, method = CCAIntegration,  
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)
srt <- FindNeighbors(srt, reduction = "integrated.cca", dims = 1:30)
srt <- RunUMAP(srt, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")

## rpca ####
srt <- IntegrateLayers(
  object = srt, method = RPCAIntegration, 
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)
srt <- FindNeighbors(srt, reduction = "integrated.rpca", dims = 1:30)
srt <- RunUMAP(srt, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")

## save log_transform-based integration result ####
saveRDS(srt,"integrated.rds" )

# 4. SCTransform Integration ####
# Load merged samples
srt <- readRDS("Merged_raw.rds")

# transform and regress data
srt <-  SCTransform(srt,  vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), return.only.var.genes = FALSE)

## unintegrated ####
# PCA, Cluster and UMAP
srt <- RunPCA(srt)
srt <- FindNeighbors(srt, dims = 1:30)
srt <- RunUMAP(srt, dims = 1:30, verbose = F)

# Integration
## harmony ####
srt <- IntegrateLayers(
  object = srt, method = HarmonyIntegration,  
  normalization.method = "SCT",
  orig.reduction = "pca", new.reduction = "integrated.harmony",
  verbose = FALSE
)
srt <- FindNeighbors(srt, reduction = "integrated.harmony", dims = 1:30)
srt <- RunUMAP(srt, reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.harmony")

## cca ####
srt <- IntegrateLayers(
  object = srt, method = CCAIntegration,  
  normalization.method = "SCT",
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)
srt <- FindNeighbors(srt, reduction = "integrated.cca", dims = 1:30)
srt <- RunUMAP(srt, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")


## rpca ####
srt <- IntegrateLayers(
  object = srt, method = RPCAIntegration, 
  normalization.method = "SCT",
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)
srt <- FindNeighbors(srt, reduction = "integrated.rpca", dims = 1:30)
srt <- RunUMAP(srt, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")

## save SCTransform-based integration result ####
saveRDS(srt,"integrated_SCT.rds" )

# 5.scVI ####
# Optional: help R find the MKL libraries if needed
Sys.setenv(LD_LIBRARY_PATH = paste(
  Sys.getenv("LD_LIBRARY_PATH"),
  "r-miniconda/envs/scvi-env/lib",
  sep = ":"
))
library(reticulate)#miniconda for scVI
reticulate::use_condaenv("scvi-env", required = TRUE)
py_config()
reticulate::import("torch")
reticulate::import("scanpy")
reticulate::import("scvi")
# Load merged samples
srt <- readRDS("Merged_raw.rds")
srt <- IntegrateLayers(
  object = srt, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "scvi-env", verbose = FALSE
)
srt <- FindNeighbors(srt, reduction = "integrated.scvi", dims = 1:30)
srt <- RunUMAP(srt, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi")

## save scVI integration result ####
saveRDS(srt, "integrated_scvi.rds")




