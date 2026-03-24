# load required packages
library(tidyverse)
library(magrittr)
library(liana)
library(Seurat)

# CCC resource
resources <- show_resources()[1:18]

# CCC methods
methods <- show_methods()[c(1,2,3,4,5,6, 8, 11)]

# Data input
srt <- readRDS("integrated_scvi.rds")
Idents(srt) <- srt$manual_annotated_1
srt <- JoinLayers(srt)
srt <- subset(srt, subset= manual_annotated_1 != "Others")
table(srt$manual_annotated_1)
gc()


#for cytotalk
#install.packages("entropy")
#for cellchat
#remotes::install_github("sqjin/CellChat")
#for italk
#remotes::install_github("Coolgenome/iTALK")
# Liana (resource= MouseConsensus, methods= all)
liana_result <- liana_wrap(
  sce = srt,
  assay = "RNA",
  slot = "data",
  method = methods ,  # etc.
  resource =  c("MouseConsensus") # ← specify the database here
)
saveRDS(liana_result,"liana_all_condition_MouseConsensus.rds")
rm(liana_result)
gc()

# Human Gene Symbol Conversion
#BiocManager::install("biomaRt")
library(biomaRt)
# Connect to Ensembl
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
# Extract mouse gene list from data set
mouse_genes <- rownames(srt)  
# Mouse gene match its human ortholog
mouse_to_human <- getLDS(
  attributes = c("mgi_symbol"), # Mouse gene symbolmapped_genes$MGI.symbol
  filters = "mgi_symbol",
  values = mouse_genes,
  mart = mouse,
  attributesL = c("hgnc_symbol"),# Human gene symbol
  martL = human,
  uniqueRows = TRUE
)
# Remove genes with no human ortholog
mapped_genes <- mouse_to_human
mapped_genes <- mapped_genes[!mapped_genes$HGNC.symbol == "", ]
mapped_genes <- mapped_genes[!duplicated(mapped_genes), ]
# Remove one mouse gene with multiple human ortholog
mapped_genes <- mapped_genes[!duplicated(mapped_genes$MGI.symbol), ]
mapped_genes$MGI.symbol[duplicated(mapped_genes$MGI.symbol)]
# Remove multiple mouse genes with one human ortholog
mapped_genes <- mapped_genes[!duplicated(mapped_genes$HGNC.symbol), ]
mapped_genes$HGNC.symbol[duplicated(mapped_genes$HGNC.symbol)]
# Filter to genes in your data
gene_map <- mapped_genes[match(rownames(srt), mapped_genes$MGI.symbol), ]
# Check how many genes have NA names
sum(is.na(gene_map$MGI.symbol))
sum(is.na(gene_map$HGNC.symbol))
# Replace rownames with human gene symbols
rownames(srt) <- gene_map$HGNC.symbol
# Check how many genes have NA names
sum(is.na(rownames(srt)))
# Remove genes with NA names
srt <- srt[!is.na(rownames(srt)), ]
# save new dataset
saveRDS(srt, "integrated_human_gene_name_dec2021.rds")
# Check how many gene removed may involved in CCC analysis
removed_gene_mouse <- setdiff(mouse_genes, mapped_genes$MGI.symbol)
df <- readRDS("liana_all_condition_MouseConsensus.rds")
# Store unique ligand and receptor genes for each method
results <- c()
# Loop over each method
for (method in methods) {
  # Access ligand and receptor columns dynamically
  ligand_genes <- df[[method]]$ligand
  receptor_genes <- df[[method]]$receptor
  
  # Combine and extract unique genes
  all_genes <- unique(c(ligand_genes, receptor_genes))
  
  # Store in results
  results<- unique(c(results,all_genes))
}

# Check overlap with removed genes
comm_genes_lost <- intersect(removed_gene_mouse, results)
comm_genes_lost
saveRDS(comm_genes_lost, "gene_lost_dec2021.rds")
