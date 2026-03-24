library(Seurat)

# LISI (both) ####
#devtools::install_github("immunogenomics/lisi")
# load pacakge
library(lisi)

## log_tranform ####
# Load data
srt <- readRDS("integrated.rds")
srt <- subset(srt, manual_annotated_1 != "Others")
# assign tranformation
transform <- "Log"

### unintegrated ####
# Get PCA embedding for raw data
pca_data <- Embeddings(srt, reduction = "pca")[,1:30]
meta_data <- data.frame(Batch = srt$Batch, 
                        Cell_Type = srt$manual_annotated_1)
# Run lisi
result <- compute_lisi(pca_data, meta_data , 
                       label_colnames = c("Batch", "Cell_Type"))
rm(pca_data, meta_data)
# Add to lisi result list
lisi <- list("Log_unintegrated" = result)

# extract batch score
ilisi.score <- summary(result$Batch)
ilisi.df <- data.frame(matrix(ilisi.score, nrow = 1))
colnames(ilisi.df) <- names(ilisi.score)
ilisi.df$Integration <- "Log_unintegrated"
# add to summary table
ilisi.summary <- ilisi.df
#clean up
rm(ilisi.score, ilisi.df)

# extract celltype score
clisi.score <- summary(result$Cell_Type)
clisi.df <- data.frame(matrix(clisi.score, nrow = 1))
colnames(clisi.df) <- names(clisi.score)
clisi.df$Integration <- "Log_unintegrated"
# add to summary table
clisi.summary <- clisi.df
#clean up
rm(clisi.score, clisi.df)
rm(result)

### integrated ####
integration <- c("cca", "rpca", "harmony", "mnn")
for(i in integration){
  # asign methods
  method <- paste(transform, i, sep = "_")
  # Get PCA embedding for raw data
  pca_data <- Embeddings(srt, reduction = paste0("integrated.", i))[,1:30]
  meta_data <- data.frame(Batch = srt$Batch, 
                          Cell_Type = srt$manual_annotated_1)
  # Run lisi
  result <- compute_lisi(pca_data, meta_data , 
                         label_colnames = c("Batch", "Cell_Type"))
  rm(pca_data, meta_data)
  # Add to lisi result list
  lisi[[method]] <- result
  
  # extract batch score
  ilisi.score <- summary(result$Batch)
  ilisi.df <- data.frame(matrix(ilisi.score, nrow = 1))
  colnames(ilisi.df) <- names(ilisi.score)
  ilisi.df$Integration <- method
  # add to summary table
  ilisi.summary <- rbind(ilisi.summary, ilisi.df)
  #clean up
  rm(ilisi.score, ilisi.df)
  
  # extract celltype score
  clisi.score <- summary(result$Cell_Type)
  clisi.df <- data.frame(matrix(clisi.score, nrow = 1))
  colnames(clisi.df) <- names(clisi.score)
  clisi.df$Integration <- method
  # add to summary table
  clisi.summary <- rbind(clisi.summary, clisi.df)
  #clean up
  rm(clisi.score, clisi.df, result, i, method)
  
}
rm(transform, integration, srt)

## SCT ####
# Load data
srt <- readRDS("integrated_SCT.rds")
srt <- subset(srt, manual_annotated_1 != "Others")
# assign tranformation
transform <- "SCT"

### unintegrated ####
# Get PCA embedding for raw data
pca_data <- Embeddings(srt, reduction = "pca")[,1:30]
meta_data <- data.frame(Batch = srt$Batch, 
                        Cell_Type = srt$manual_annotated_1)
# Run lisi
result <- compute_lisi(pca_data, meta_data , 
                       label_colnames = c("Batch", "Cell_Type"))
rm(pca_data, meta_data)
# Add to lisi result list
lisi[["SCT_unintegrated"]] <- result

# extract batch score
ilisi.score <- summary(result$Batch)
ilisi.df <- data.frame(matrix(ilisi.score, nrow = 1))
colnames(ilisi.df) <- names(ilisi.score)
ilisi.df$Integration <- "SCT_unintegrated"
# add to summary table
ilisi.summary <- rbind(ilisi.summary, ilisi.df)
#clean up
rm(ilisi.score, ilisi.df)

# extract celltype score
clisi.score <- summary(result$Cell_Type)
clisi.df <- data.frame(matrix(clisi.score, nrow = 1))
colnames(clisi.df) <- names(clisi.score)
clisi.df$Integration <- "SCT_unintegrated"
# add to summary table
clisi.summary <- rbind(clisi.summary, clisi.df)
#clean up
rm(clisi.score, clisi.df)
rm(result)

### integrated ####
integration <- c("cca", "rpca", "harmony")
for(i in integration){
  # asign methods
  method <- paste(transform, i, sep = "_")
  # Get PCA embedding for raw data
  pca_data <- Embeddings(srt, reduction = paste0("integrated.", i))[,1:30]
  meta_data <- data.frame(Batch = srt$Batch, 
                          Cell_Type = srt$manual_annotated_1)
  # Run lisi
  result <- compute_lisi(pca_data, meta_data , 
                         label_colnames = c("Batch", "Cell_Type"))
  rm(pca_data, meta_data)
  # Add to lisi result list
  lisi[[method]] <- result
  
  # extract batch score
  ilisi.score <- summary(result$Batch)
  ilisi.df <- data.frame(matrix(ilisi.score, nrow = 1))
  colnames(ilisi.df) <- names(ilisi.score)
  ilisi.df$Integration <- method
  # add to summary table
  ilisi.summary <- rbind(ilisi.summary, ilisi.df)
  #clean up
  rm(ilisi.score, ilisi.df)
  
  # extract celltype score
  clisi.score <- summary(result$Cell_Type)
  clisi.df <- data.frame(matrix(clisi.score, nrow = 1))
  colnames(clisi.df) <- names(clisi.score)
  clisi.df$Integration <- method
  # add to summary table
  clisi.summary <- rbind(clisi.summary, clisi.df)
  #clean up
  rm(clisi.score, clisi.df, result, i, method)
  
}
rm(transform, integration, srt)

## scVI ####
# Load data
srt <- readRDS("integrated_scvi.rds")
srt <- subset(srt, manual_annotated_1 != "Others")
# Get PCA embedding for raw data
pca_data <- Embeddings(srt, reduction = "integrated.scvi")[,1:30]
meta_data <- data.frame(Batch = srt$Batch, 
                        Cell_Type = srt$manual_annotated_1)
# Run lisi
result <- compute_lisi(pca_data, meta_data , 
                       label_colnames = c("Batch", "Cell_Type"))
rm(pca_data, meta_data)

# Add to lisi result list
lisi[["scvi"]] <- result

# extract batch score
ilisi.score <- summary(result$Batch)
ilisi.df <- data.frame(matrix(ilisi.score, nrow = 1))
colnames(ilisi.df) <- names(ilisi.score)
ilisi.df$Integration <- "scvi"
# add to summary table
ilisi.summary <- rbind(ilisi.summary, ilisi.df)
#clean up
rm(ilisi.score, ilisi.df)

# extract celltype score
clisi.score <- summary(result$Cell_Type)
clisi.df <- data.frame(matrix(clisi.score, nrow = 1))
colnames(clisi.df) <- names(clisi.score)
clisi.df$Integration <- "scvi"
# add to summary table
clisi.summary <- rbind(clisi.summary, clisi.df)
#clean up
rm(clisi.score, clisi.df)
rm(result, srt)

## Save data ####
saveRDS(lisi, "lisi.rds")
saveRDS(ilisi.summary, "ilisi.summary.rds")
saveRDS(clisi.summary, "clisi.summary.rds")
rm(lisi,ilisi.summary,clisi.summary)



# ASW (batch)####
# load package
library(cluster)
library(tidyverse)
## log_tranform ####
# Load data
srt <- readRDS("integrated.rds")
srt <- subset(srt, manual_annotated_1 != "Others")

# assign tranformation
transform <- "Log"
### unintegrated ####
# Get PCA embedding for raw data
pca_data <- Embeddings(srt, reduction = "pca")[,1:30]
#  compute distance matrix
dists <- dist(pca_data)
# compute silhouette
result <- data.frame(silhouette(as.numeric(factor(srt$Batch)), dists))
rownames(result) <- colnames(srt)
# match cluter to batch
result$batch <- srt$Batch
# Add to basw result list
basw <- list("Log_unintegrated" = result)

# extract batch score
basw.score <- summary(result$sil_width)
basw.df <- data.frame(matrix(basw.score, nrow = 1))
colnames(basw.df) <- names(basw.score)
basw.df$Integration <- "Log_unintegrated"
# add to summary table
basw.summary <- basw.df

# clean up
rm(pca_data, result,basw.score, basw.df, dists)

### integrated ####
integration <- c("cca", "rpca", "harmony", "mnn")
for(i in integration){
  # asign methods
  method <- paste(transform, i, sep = "_")
  # Get PCA embedding for raw data
  pca_data <- Embeddings(srt, reduction = paste0("integrated.", i))[,1:30]
  #  compute distance matrix
  dists <- dist(pca_data)
  # compute silhouette
  result <- data.frame(silhouette(as.numeric(factor(srt$Batch)), dists))
  rownames(result) <- colnames(srt)
  # match cluter to batch
  result$batch <- srt$Batch
  # Add to basw result list
  basw[[method]] <- result
  
  # extract batch score
  basw.score <- summary(result$sil_width)
  basw.df <- data.frame(matrix(basw.score, nrow = 1))
  colnames(basw.df) <- names(basw.score)
  basw.df$Integration <- method
  # add to summary table
  basw.summary <- rbind(basw.summary,basw.df)
  # clean up
  rm(i, method, pca_data, dists, result, basw.score, basw.df)
}
rm(transform, integration, srt)
gc()
## SCT####
# Load data
srt <- readRDS("integrated_SCT.rds")
srt <- subset(srt, manual_annotated_1 != "Others")
# assign tranformation
transform <- "SCT"

### unintegrated ####
# Get PCA embedding for raw data
pca_data <- Embeddings(srt, reduction = "pca")[,1:30]
#  compute distance matrix
dists <- dist(pca_data)
# compute silhouette
result <- data.frame(silhouette(as.numeric(factor(srt$Batch)), dists))
rownames(result) <- colnames(srt)
# match cluter to batch
result$batch <- srt$Batch
# Add to basw result list
basw[["SCT_unintegrated"]] <- result

# extract batch score
basw.score <- summary(result$sil_width)
basw.df <- data.frame(matrix(basw.score, nrow = 1))
colnames(basw.df) <- names(basw.score)
basw.df$Integration <- "SCT_unintegrated"
# add to summary table
basw.summary <- rbind(basw.summary, basw.df)

# clean up
rm(pca_data, result,basw.score, basw.df, dists)

### integrated ####
integration <- c("cca", "rpca", "harmony")
for(i in integration){
  # asign methods
  method <- paste(transform, i, sep = "_")
  # Get PCA embedding for raw data
  pca_data <- Embeddings(srt, reduction = paste0("integrated.", i))[,1:30]
  #  compute distance matrix
  dists <- dist(pca_data)
  # compute silhouette
  result <- data.frame(silhouette(as.numeric(factor(srt$Batch)), dists))
  rownames(result) <- colnames(srt)
  # match cluter to batch
  result$batch <- srt$Batch
  # Add to basw result list
  basw[[method]] <- result
  
  # extract batch score
  basw.score <- summary(result$sil_width)
  basw.df <- data.frame(matrix(basw.score, nrow = 1))
  colnames(basw.df) <- names(basw.score)
  basw.df$Integration <- method
  # add to summary table
  basw.summary <- rbind(basw.summary,basw.df)
  # clean up
  rm(i, method, pca_data, dists, result, basw.score, basw.df)
}
rm(transform, integration, srt)



## scVI####
# Load data
srt <- readRDS("integrated_scvi.rds")
srt <- subset(srt, manual_annotated_1 != "Others")
# Get batch labels (from metadata)
batch_labels <- srt$Batch
# Get PCA embedding for raw data
pca_data <- Embeddings(srt, reduction = "integrated.scvi")[,1:30]
#  compute distance matrix
dists <- dist(pca_data)
# compute silhouette
result <- data.frame(silhouette(as.numeric(factor(srt$Batch)), dists))
rownames(result) <- colnames(srt)
# match cluter to batch
result$batch <- srt$Batch
# Add to basw result list
basw[["scvi"]] <- result
# extract batch score
basw.score <- summary(result$sil_width)
basw.df <- data.frame(matrix(basw.score, nrow = 1))
colnames(basw.df) <- names(basw.score)
basw.df$Integration <- "scvi"
# add to summary table
basw.summary <- rbind(basw.summary, basw.df)

# clean up
rm(pca_data, result, basw.score, basw.df, dists, batch_labels)

## Save data ####
saveRDS(basw, "basw.rds")
saveRDS(basw.summary, "basw.summary.rds")
rm(basw, basw.summary)

#ASW(cell_type)####
## log_tranform ####
# Load data
srt <- readRDS("integrated.rds")
srt <- subset(srt, manual_annotated_1 != "Others")

# assign tranformation
transform <- "Log"
### unintegrated ####
# Get PCA embedding for raw data
pca_data <- Embeddings(srt, reduction = "pca")[,1:30]
#  compute distance matrix
dists <- dist(pca_data)
# compute silhouette
result <- data.frame(silhouette(as.numeric(factor(srt$manual_annotated_1)), dists))
rownames(result) <- colnames(srt)
# match cluter to cell type
result$cell_type <- srt$manual_annotated_1
# Add to basw result list
casw <- list("Log_unintegrated" = result)

# extract cell type score
casw.score <- summary(result$sil_width)
casw.df <- data.frame(matrix(casw.score, nrow = 1))
colnames(casw.df) <- names(casw.score)
casw.df$Integration <- "Log_unintegrated"
# add to summary table
casw.summary <- casw.df

# clean up
rm(pca_data, result,casw.score, casw.df, dists)

### integrated ####
integration <- c("cca", "rpca", "harmony", "mnn")
for(i in integration){
  # asign methods
  method <- paste(transform, i, sep = "_")
  # Get PCA embedding for raw data
  pca_data <- Embeddings(srt, reduction = paste0("integrated.", i))[,1:30]
  #  compute distance matrix
  dists <- dist(pca_data)
  # compute silhouette
  result <- data.frame(silhouette(as.numeric(factor(srt$manual_annotated_1)), dists))
  rownames(result) <- colnames(srt)
  # match cluter to batch
  result$cell_type <- srt$manual_annotated_1
  # Add to basw result list
  casw[[method]] <- result
  
  # extract batch score
  casw.score <- summary(result$sil_width)
  casw.df <- data.frame(matrix(casw.score, nrow = 1))
  colnames(casw.df) <- names(casw.score)
  casw.df$Integration <- method
  # add to summary table
  casw.summary <- rbind(casw.summary, casw.df)
  # clean up
  rm(i, method, pca_data, dists, result, casw.score, casw.df)
}
rm(transform, integration, srt)
gc()
## SCT####
# Load data
srt <- readRDS("integrated_SCT.rds")
srt <- subset(srt, manual_annotated_1 != "Others")
# assign tranformation
transform <- "SCT"

### unintegrated ####
# Get PCA embedding for raw data
pca_data <- Embeddings(srt, reduction = "pca")[,1:30]
#  compute distance matrix
dists <- dist(pca_data)
# compute silhouette
result <- data.frame(silhouette(as.numeric(factor(srt$manual_annotated_1)), dists))
rownames(result) <- colnames(srt)
# match cluter to batch
result$cell_type <- srt$manual_annotated_1
# Add to basw result list
casw[["SCT_unintegrated"]] <- result

# extract batch score
casw.score <- summary(result$sil_width)
casw.df <- data.frame(matrix(casw.score, nrow = 1))
colnames(casw.df) <- names(casw.score)
casw.df$Integration <- "SCT_unintegrated"
# add to summary table
casw.summary <- rbind(casw.summary, casw.df)

# clean up
rm(pca_data, result,casw.score, casw.df, dists)

### integrated ####
integration <- c("cca", "rpca", "harmony")
for(i in integration){
  # asign methods
  method <- paste(transform, i, sep = "_")
  # Get PCA embedding for raw data
  pca_data <- Embeddings(srt, reduction = paste0("integrated.", i))[,1:30]
  #  compute distance matrix
  dists <- dist(pca_data)
  # compute silhouette
  result <- data.frame(silhouette(as.numeric(factor(srt$manual_annotated_1)), dists))
  rownames(result) <- colnames(srt)
  # match cluter to batch
  result$cell_type <- srt$manual_annotated_1
  # Add to basw result list
  casw[[method]] <- result
  
  # extract batch score
  casw.score <- summary(result$sil_width)
  casw.df <- data.frame(matrix(casw.score, nrow = 1))
  colnames(casw.df) <- names(casw.score)
  casw.df$Integration <- method
  # add to summary table
  casw.summary <- rbind(casw.summary,casw.df)
  # clean up
  rm(i, method, pca_data, dists, result, casw.score, casw.df)
}
rm(transform, integration, srt)

## scVI####
# Load data
srt <- readRDS("integrated_scvi.rds")
srt <- subset(srt, manual_annotated_1 != "Others")

# Get PCA embedding for raw data
pca_data <- Embeddings(srt, reduction = "integrated.scvi")[,1:30]
#  compute distance matrix
dists <- dist(pca_data)
# compute silhouette
result <- data.frame(silhouette(as.numeric(factor(srt$manual_annotated_1)), dists))
rownames(result) <- colnames(srt)
# match cluter to batch
result$cell_type <- srt$manual_annotated_1
# Add to basw result list
casw[["scvi"]] <- result

# extract batch score
casw.score <- summary(result$sil_width)
casw.df <- data.frame(matrix(casw.score, nrow = 1))
colnames(casw.df) <- names(casw.score)
casw.df$Integration <- "scvi"
# add to summary table
casw.summary <- rbind(casw.summary, casw.df)

# clean up
rm(pca_data, result,casw.score, casw.df, dists)

## Save data ####
saveRDS(casw, "casw.rds")
saveRDS(casw.summary, "casw.summary.rds")
rm(casw, casw.summary, srt)



# Graph LISI (both)####
library(Matrix)
# Function to get neighbors for each cell from graph
get_neighbors <- function(graph_matrix) {
  neighbors_list <- apply(graph_matrix, 1, function(row) {
    which(row > 0)
  })
  return(neighbors_list)
}
# Function to compute Inverse Simpson Index
inverse_simpson <- function(labels) {
  freq_table <- table(labels) / length(labels)
  isi <- 1 / sum(freq_table^2)
  return(isi)
}
## log_tranform ####
# Load data
srt <- readRDS("integrated.rds")
srt <- subset(srt, manual_annotated_1 != "Others")
# assign tranformation
transform <- "Log"
# create labels
batch_labels <- srt$Batch # or other batch column
cell_type_labels <- srt$manual_annotated_1 # Get cell type labels
### unintegrated ####
# compute snn
srt <- FindNeighbors(srt, reduction = "pca")
# Extract SNN graph from Seurat
graph <- as.matrix(srt@graphs$RNA_snn) # sparse matrix
# Get neighbors for each cell from graph
neighbors <- get_neighbors(graph)
rm(graph)
# Compute iLISI for each cell
ilisi_scores <- sapply(1:length(neighbors), function(i) {
  neighbor_ids <- neighbors[[i]]
  neighbor_batches <- batch_labels[neighbor_ids]
  inverse_simpson(neighbor_batches)
})
# graph iLISI summary
ilisi_score_summ <- summary(ilisi_scores)
ilisi_score_summ_df <- data.frame(matrix(ilisi_score_summ, nrow = 1))
colnames(ilisi_score_summ_df) <- names(ilisi_score_summ)
ilisi_score_summ_df$Integration <- "Log_unintegrated"
# add to summary table
gilisi.summary <- ilisi_score_summ_df
#clean up
rm(ilisi_score_summ, ilisi_score_summ_df)

# Compute cLISI using cell type labels in each cell's graph neighborhood
clisi_scores <- sapply(1:length(neighbors), function(i) {
  neighbor_ids <- neighbors[[i]]
  neighbor_types <- cell_type_labels[neighbor_ids]
  inverse_simpson(neighbor_types)
})
# graph iLISI summary
clisi_score_summ <- summary(clisi_scores)
clisi_score_summ_df <- data.frame(matrix(clisi_score_summ, nrow = 1))
colnames(clisi_score_summ_df) <- names(clisi_score_summ)
clisi_score_summ_df$Integration <- "Log_unintegrated"
# add to summary table
gclisi.summary <- clisi_score_summ_df
#clean up
rm(clisi_score_summ, clisi_score_summ_df)
# Creat a data frame
result <- data.frame(Batch = ilisi_scores, Cell_Type = clisi_scores)
rownames(result) <- colnames(srt)
glisi <- list("Log_unintegrated" = result)
rm(ilisi_scores, clisi_scores, result, neighbors)

### integrated ####
integration <- c("cca", "rpca", "harmony", "mnn")
for(i in integration){
  # asign methods
  method <- paste(transform, i, sep = "_")
  # compute snn
  srt <- FindNeighbors(srt, reduction = paste0("integrated.", i))
  # Extract SNN graph from Seurat
  graph <- as.matrix(srt@graphs$RNA_snn) # sparse matrix
  # Get neighbors for each cell from graph
  neighbors <- get_neighbors(graph)
  rm(graph)
  # Compute iLISI for each cell
  ilisi_scores <- sapply(1:length(neighbors), function(i) {
    neighbor_ids <- neighbors[[i]]
    neighbor_batches <- batch_labels[neighbor_ids]
    inverse_simpson(neighbor_batches)
  })
  # graph iLISI summary
  ilisi_score_summ <- summary(ilisi_scores)
  ilisi_score_summ_df <- data.frame(matrix(ilisi_score_summ, nrow = 1))
  colnames(ilisi_score_summ_df) <- names(ilisi_score_summ)
  ilisi_score_summ_df$Integration <- paste(transform, i, sep = "_")
  # add to summary table
  gilisi.summary <- rbind(gilisi.summary,ilisi_score_summ_df)
  #clean up
  rm(ilisi_score_summ, ilisi_score_summ_df)
  
  # Compute cLISI using cell type labels in each cell's graph neighborhood
  clisi_scores <- sapply(1:length(neighbors), function(i) {
    neighbor_ids <- neighbors[[i]]
    neighbor_types <- cell_type_labels[neighbor_ids]
    inverse_simpson(neighbor_types)
  })
  # graph cLISI summary
  clisi_score_summ <- summary(clisi_scores)
  clisi_score_summ_df <- data.frame(matrix(clisi_score_summ, nrow = 1))
  colnames(clisi_score_summ_df) <- names(clisi_score_summ)
  clisi_score_summ_df$Integration <- paste(transform, i, sep = "_")
  # add to summary table
  gclisi.summary <- rbind(gclisi.summary, clisi_score_summ_df)
  #clean up
  rm(clisi_score_summ, clisi_score_summ_df)
  
  # Creat a data frame
  result <- data.frame(Batch = ilisi_scores, Cell_Type = clisi_scores)
  rownames(result) <- colnames(srt)
  glisi[[method]] <-result
  rm(ilisi_scores, clisi_scores, result, neighbors, i, method)
  
}  
rm(srt, batch_labels, cell_type_labels, transform, integration)
gc()
## scVI ####
# Load data
srt <- readRDS("integrated_scvi.rds")
srt <- subset(srt, manual_annotated_1 != "Others")
# create labels
batch_labels <- srt$Batch # or other batch column
cell_type_labels <- srt$manual_annotated_1 # Get cell type labels
# compute snn
srt <- FindNeighbors(srt, reduction = "integrated.scvi")
# Extract SNN graph from Seurat
graph <- as.matrix(srt@graphs$RNA_snn) # sparse matrix
# Get neighbors for each cell from graph
neighbors <- get_neighbors(graph)
rm(graph)
# Compute iLISI for each cell
ilisi_scores <- sapply(1:length(neighbors), function(i) {
  neighbor_ids <- neighbors[[i]]
  neighbor_batches <- batch_labels[neighbor_ids]
  inverse_simpson(neighbor_batches)
})
# graph iLISI summary
ilisi_score_summ <- summary(ilisi_scores)
ilisi_score_summ_df <- data.frame(matrix(ilisi_score_summ, nrow = 1))
colnames(ilisi_score_summ_df) <- names(ilisi_score_summ)
ilisi_score_summ_df$Integration <- "scvi"
# add to summary table
gilisi.summary <- rbind(gilisi.summary, ilisi_score_summ_df)
#clean up
rm(ilisi_score_summ, ilisi_score_summ_df)

# Compute cLISI using cell type labels in each cell's graph neighborhood
clisi_scores <- sapply(1:length(neighbors), function(i) {
  neighbor_ids <- neighbors[[i]]
  neighbor_types <- cell_type_labels[neighbor_ids]
  inverse_simpson(neighbor_types)
})
# graph iLISI summary
clisi_score_summ <- summary(clisi_scores)
clisi_score_summ_df <- data.frame(matrix(clisi_score_summ, nrow = 1))
colnames(clisi_score_summ_df) <- names(clisi_score_summ)
clisi_score_summ_df$Integration <- "scvi"
# add to summary table
gclisi.summary <- rbind(gclisi.summary, clisi_score_summ_df)
#clean up
rm(clisi_score_summ, clisi_score_summ_df)
# Creat a data frame
result <- data.frame(Batch = ilisi_scores, Cell_Type = clisi_scores)
rownames(result) <- colnames(srt)
glisi[["scvi"]] <- result
rm(ilisi_scores, clisi_scores, result, neighbors)
rm(srt, batch_labels, cell_type_labels)
gc()
## SCT ####
# Load data
srt <- readRDS("integrated_SCT.rds")
srt <- subset(srt, manual_annotated_1 != "Others")
# assign tranformation
transform <- "SCT"
# create labels
batch_labels <- srt$Batch # or other batch column
cell_type_labels <- srt$manual_annotated_1 # Get cell type labels
### unintegrated ####
# compute snn
srt <- FindNeighbors(srt, reduction = "pca")
# Extract SNN graph from Seurat
graph <- as.matrix(srt@graphs$SCT_snn) # sparse matrix
# Get neighbors for each cell from graph
neighbors <- get_neighbors(graph)
rm(graph)
# Compute iLISI for each cell
ilisi_scores <- sapply(1:length(neighbors), function(i) {
  neighbor_ids <- neighbors[[i]]
  neighbor_batches <- batch_labels[neighbor_ids]
  inverse_simpson(neighbor_batches)
})
# graph iLISI summary
ilisi_score_summ <- summary(ilisi_scores)
ilisi_score_summ_df <- data.frame(matrix(ilisi_score_summ, nrow = 1))
colnames(ilisi_score_summ_df) <- names(ilisi_score_summ)
ilisi_score_summ_df$Integration <- "SCT_unintegrated"
# add to summary table
gilisi.summary <- rbind(gilisi.summary, ilisi_score_summ_df)
#clean up
rm(ilisi_score_summ, ilisi_score_summ_df)

# Compute cLISI using cell type labels in each cell's graph neighborhood
clisi_scores <- sapply(1:length(neighbors), function(i) {
  neighbor_ids <- neighbors[[i]]
  neighbor_types <- cell_type_labels[neighbor_ids]
  inverse_simpson(neighbor_types)
})
# graph iLISI summary
clisi_score_summ <- summary(clisi_scores)
clisi_score_summ_df <- data.frame(matrix(clisi_score_summ, nrow = 1))
colnames(clisi_score_summ_df) <- names(clisi_score_summ)
clisi_score_summ_df$Integration <- "SCT_unintegrated"
# add to summary table
gclisi.summary <- rbind(gclisi.summary, clisi_score_summ_df)
#clean up
rm(clisi_score_summ, clisi_score_summ_df)
# Creat a data frame
result <- data.frame(Batch = ilisi_scores, Cell_Type = clisi_scores)
rownames(result) <- colnames(srt)
glisi[["SCT_unintegrated"]] <- result
rm(ilisi_scores, clisi_scores, result, neighbors)
### integrated ####
integration <- c("cca", "rpca", "harmony")
for(i in integration){
  # asign methods
  method <- paste(transform, i, sep = "_")
  # compute snn
  srt <- FindNeighbors(srt, reduction = paste0("integrated.", i))
  # Extract SNN graph from Seurat
  graph <- as.matrix(srt@graphs$SCT_snn) # sparse matrix
  # Get neighbors for each cell from graph
  neighbors <- get_neighbors(graph)
  rm(graph)
  # Compute iLISI for each cell
  ilisi_scores <- sapply(1:length(neighbors), function(i) {
    neighbor_ids <- neighbors[[i]]
    neighbor_batches <- batch_labels[neighbor_ids]
    inverse_simpson(neighbor_batches)
  })
  # graph iLISI summary
  ilisi_score_summ <- summary(ilisi_scores)
  ilisi_score_summ_df <- data.frame(matrix(ilisi_score_summ, nrow = 1))
  colnames(ilisi_score_summ_df) <- names(ilisi_score_summ)
  ilisi_score_summ_df$Integration <- paste(transform, i, sep = "_")
  # add to summary table
  gilisi.summary <- rbind(gilisi.summary,ilisi_score_summ_df)
  #clean up
  rm(ilisi_score_summ, ilisi_score_summ_df)
  
  # Compute cLISI using cell type labels in each cell's graph neighborhood
  clisi_scores <- sapply(1:length(neighbors), function(i) {
    neighbor_ids <- neighbors[[i]]
    neighbor_types <- cell_type_labels[neighbor_ids]
    inverse_simpson(neighbor_types)
  })
  # graph cLISI summary
  clisi_score_summ <- summary(clisi_scores)
  clisi_score_summ_df <- data.frame(matrix(clisi_score_summ, nrow = 1))
  colnames(clisi_score_summ_df) <- names(clisi_score_summ)
  clisi_score_summ_df$Integration <- paste(transform, i, sep = "_")
  # add to summary table
  gclisi.summary <- rbind(gclisi.summary, clisi_score_summ_df)
  #clean up
  rm(clisi_score_summ, clisi_score_summ_df)
  
  # Creat a data frame
  result <- data.frame(Batch = ilisi_scores, Cell_Type = clisi_scores)
  rownames(result) <- colnames(srt)
  glisi[[method]] <-result
  rm(ilisi_scores, clisi_scores, result, neighbors, i, method)
  
}  
rm(srt, batch_labels, cell_type_labels, transform, integration)
rm(get_neighbors, inverse_simpson)
gc()
## Save data ####
saveRDS(glisi, "glisi.rds")
saveRDS(gilisi.summary, "gilisi.summary.rds")
saveRDS(gclisi.summary, "/gclisi.summary.rds")
rm(glisi,gilisi.summary,gclisi.summary)


