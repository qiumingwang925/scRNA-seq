# install cellchat package link: https://github.com/jinworks/CellChat
devtools::install_github("KlugerLab/ALRA")
devtools::install_github("JEFworks-Lab/MERINGUE")
devtools::install_github("jokergoo/ComplexHeatmap")
devtools::install_github("jinworks/CellChat")

library(Seurat)
library(CellChat)
library(patchwork)

options(stringsAsFactors = FALSE)

# -- mod_interact_upload --
# upload .rds seuart object from local.
##  check assays: print RNA or SCT
##  print column names in meta.data. 
##  Ask user to select one for condition/treatment. If it's not called "group", than save as colname "group"
##  print active.ident: levels(Ident(obj))
## join layers if more than one layer of "data" in RNA or SCT in assays obj <- JoinLayers(obj)

# -- mod_interact_cellchat_comp --
## use the .rds seurat object processed in mod_interact_upload as input
## Create a cellchat obj and specify assays RNA or SCT
cellchat <- createCellChat(object = obj, group.by = "ident", assay = "RNA")
cellchat <- createCellChat(object = obj, group.by = "ident", assay = "SCT")
## select if the study is mouse or human
CellChatDB <- CellChatDB.mouse
CellChatDB <- CellChatDB.human
## use a subset of CellChatDB for cell-cell communication analysis
CellChatDB <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
## subset the expression data of signaling genes for saving computation cost
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
## Compute the communication probability and infer cellular communication network
ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "triMean")
## Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
## Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
## add download button to save the cellchat result

# -- mod_interact_cellchat_vis --
## if no cellchat result file pass from mod_interact_cellchat_comp, upload a processed cellchat result
## subtab 1 global network analysis (
### cell-type-based Circle
### cell-type-based Heatmap
## subtab 2 Zooming in (pathway & pair level)
### Pathway-based Hierarchy
### Pathway-based Circle
### Pathway-based Heatmap
### Pathway-based Chord
### PairLR-based Hierarchy
### PairLR-based Circle
### Pathway&PairLR-based Chord
### Pathway&PairLR-based Bubble
### Pathway&PairLR-based Violin
## subtab 3 signaling-focused
## Signaling-based Score
## Signaling-based Scatter
## Signaling-based Heatmap
## subtab 4 Communication patterns
### Pattern-based Heatmap
### Pattern-based River
### Pattern-based Dot
### Manifold & Classification
## add a figure download button and numeric input for width and length adjustment for all subtabs
## reference tutorial https://github.com/jinworks/CellChat/blob/main/tutorial/CellChat-vignette.Rmd
## all the visualization and analyses need to stratify by group


