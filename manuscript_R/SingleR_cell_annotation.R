# Load required libraries ####
library(SingleR)
library(SummarizedExperiment)
library(celldex)
library(Seurat)
library(ggplot2)
library(Polychrome) # color code
library(zellkonverter)
library(SingleCellExperiment)
library(biomaRt)
set.seed(4321)

# Reference Preparation####
## Tabula Muris ####
# Ref:TabulaMuris prepare
sce.ref <- readH5AD("Tabula_Muris_Senis_Aging_CZI.h5ad") # downloaded from TMS website
seurat_obj <- as.Seurat(sce.ref, counts = "X", data = "X")
ref <- as.SingleCellExperiment(seurat_obj)
# Get mapping
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  mart = mart
)
# Rename rownames in the reference
common_ids <- intersect(rownames(ref), gene_map$ensembl_gene_id)
gene_map <- gene_map[match(common_ids, gene_map$ensembl_gene_id), ]
# Filter and rename reference genes
ref <- ref[common_ids, ]
rownames(ref) <- gene_map$external_gene_name
rm(gene_map,mart, coomon_ids)
saveRDS(ref, "Tabula_Muris_Senis_Aging_CZI_SCE.rds")


## LungMap ####
## Ref: LungMap prepare
sce.ref <- readH5AD("LungMAP_MouseLungDevelopment_CellRef.v1.h5ad") #download from LungMap website
seurat_obj <- as.Seurat(sce.ref, counts = "X", data = "X")
saveRDS(seurat_obj, "LungMAP_MouseLungDevelopment_CellRef.v1.rds")
seurat_ref <- readRDS("LungMAP_MouseLungDevelopment_CellRef.v1.rds")
ref <- as.SingleCellExperiment(seurat_obj)
saveRDS(ref, "LungMAP_MouseLungDevelopment_CellRef_SCE.rds")


# Set samples ####
samples <- c("22713X3","23957X4","21401X3", "20354X2", "23957X5", "22713X2", "20354X1", "24143X4")

# SingleR 
for(i in samples){
  ## Load  Seurat_obj ####
  seurat_obj <- readRDS(paste0("seurat_object_",i ,".rds"))
  #Convert to SCE ####
  sce <- as.SingleCellExperiment(seurat_obj)
  
  # Tabula Muris Senis ####
  ### Load reference ####
  ref <- readRDS("Tabula_Muris_Senis_Aging_CZI_SCE.rds")
  
  ### Prepare reference ####
  # Get mapping
  gene_map <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    mart = mart
  )
  # Rename rownames in the reference
  common_ids <- intersect(rownames(ref), gene_map$ensembl_gene_id)
  gene_map <- gene_map[match(common_ids, gene_map$ensembl_gene_id), ]
  # Filter and rename reference genes
  ref <- ref[common_ids, ]
  rownames(ref) <- gene_map$external_gene_name
  rm(gene_map, common_ids)
  
  ### Label: free annotation ####
  ref_labels <- ref$free_annotation
  
  pred <- SingleR(
    test = sce,
    ref = ref,
    labels = ref_labels
  )
  
  seurat_obj$SingleR.TMS.free.annotation <- pred$labels
  rm(ref_labels, pred)
  
  ### Label: cell-type ####
  ref_labels <- ref$cell_type
  
  pred <- SingleR(
    test = sce,
    ref = ref,
    labels = ref_labels
  )
  
  seurat_obj$SingleR.TMS.cell_type <- pred$labels
  rm(ref_labels, pred)
  
  ### Remove reference ####
  rm(ref)
  
  #LungMap ####
  ## Load reference ####
  ref <-readRDS("LungMAP_MouseLungDevelopment_CellRef_SCE.rds")
  
  ## Assign label and prediction ####
  ref_labels <- ref$lineage_level1
  
  pred <- SingleR(
    test = sce,
    ref = ref,
    labels = ref_labels
  )
  
  ## Add SingleR annotations to Seurat object metadata ####
  seurat_obj$SingleR.LungMap.lineage_level1 <- pred$labels
  rm(ref_labels, pred)
  
  ## Assign label and prediction ####
  ref_labels <- ref$lineage_level2
  
  pred <- SingleR(
    test = sce,
    ref = ref,
    labels = ref_labels
  )
  
  ## Add SingleR annotations to Seurat object metadata ####
  seurat_obj$SingleR.LungMap.lineage_level2 <- pred$labels
  rm(ref_labels, pred)
  
  ## Assign label 1 and prediction ####
  ref_labels <- ref$celltype_level1
  
  pred <- SingleR(
    test = sce,
    ref = ref,
    labels = ref_labels
  )
  
  ## Add SingleR annotations to Seurat object metadata ####
  seurat_obj$SingleR.LungMap.celltype_level1 <- pred$labels
  rm(ref_labels, pred)
  
  ## Assign label and prediction ####
  ref_labels <- ref$celltype_level2
  
  pred <- SingleR(
    test = sce,
    ref = ref,
    labels = ref_labels
  )
  
  ## Add SingleR annotations to Seurat object metadata ####
  seurat_obj$SingleR.LungMap.celltype_level2 <- pred$labels
  rm(ref_labels, pred)
  
  ## Assign label and prediction ####
  ref_labels <- ref$celltype_level3
  
  pred <- SingleR(
    test = sce,
    ref = ref,
    labels = ref_labels
  )
  
  ## Add SingleR annotations to Seurat object metadata ####
  seurat_obj$SingleR.LungMap.celltype_level3 <- pred$labels
  rm(ref_labels, pred)
  
  ## remove reference 
  rm(ref)
  
  # Save annotated srt
  saveRDS(seurat_obj, paste0("SingleR_Annotated_", i, ".rds"))
  rm(seurat_obj,i, sce)
  
          
}

rm(mart, samples)

