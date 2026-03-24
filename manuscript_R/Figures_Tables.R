library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)
library(ggrepel)
library(Polychrome)
library(presto)
library(tidyverse)
library(readxl)
library(ggpubr)
library(pheatmap)
#library(magick)
#library(pdftools)
library(grid)
library(gridExtra)
library(reshape2)
library(ComplexUpset)
#library(ComplexHeatmap)
#library(circlize)
library(openxlsx)

# 1. Load cell type key ####
cell_type_key <- read_excel("Documents/scRNAseq/Tool_paper/Data/cell_type_key.xlsx", sheet = 1)
# 2. Load SRT & DataPrep####
## Log ####
srt_log <- readRDS("~/Documents/scRNAseq/Tool_paper/Data/integrated.rds")
# remove cells with bad identifications
srt_log <- subset(srt_log, manual_annotated_1 != "Others")
# set "Condition" as factor
srt_log$Condition <- factor(srt_log$Condition, levels = c("Control", "Inflam", "Fib", "Fib_O3"))
# set manual_annotated_cell type names as factor
srt_log$manual_annotated_1 <- factor(srt_log$manual_annotated_1, 
                                     levels = c("AM", "MoDM", "IM", "ncMono", "cMono", 
                                                "MoDC", "cDC1", "cDC2", "Ccr7_DC", "pDC",
                                                "B", "Plasma",  "Neutrophil", "Eosinophil", "Mast",
                                                "CD4_T", "CD8_T","Treg", "DN_T", "NKT", "NK", "ILC",
                                                "Ciliated", "Club" , "Alv_Epi",
                                                "Artery", "Vein", "Lymphatics", "Cap_A", "Cap_G", "EPC",
                                                "Fibroblast", "SMC", 
                                                "Mesothelial", "Pericyte",
                                                "Megakaryocyte"))
### data summary ####
df <- data.frame(table(srt_log$manual_annotated_1))
colnames(df) <- c("Fibrosis_Label", "Number of Cells")
df_summ <- cell_type_key %>% dplyr::distinct(Category, Fibrosis_Label, Fibrosis_FullName) %>% inner_join(., df, by = "Fibrosis_Label")
total_cells <- sum(df_summ$`Number of Cells`)
df_summ$`Percent` <- round(df_summ$`Number of Cells` / total_cells *100, digits = 2)
write.xlsx(df_summ, "~/Documents/scRNAseq/Tool_paper/Data/Figure_1_Cell_Summary.xlsx")
rm(df_summ)
## SCT ####
srt_sct <- readRDS("~/Documents/scRNAseq/Tool_paper/Data/integrated_SCT.rds")
# remove cells with bad identifications
srt_sct <- subset(srt_sct, manual_annotated_1 != "Others")
# set "Condition" as factor
srt_sct$Condition <- factor(srt_sct$Condition, levels = c("Control", "Inflam", "Fib", "Fib_O3"))
# set manual_annotated_cell type names as factor
srt_sct$manual_annotated_1 <- factor(srt_sct$manual_annotated_1, 
                                     levels = c("AM", "MoDM", "IM", "ncMono", "cMono", 
                                                "MoDC", "cDC1", "cDC2", "Ccr7_DC", "pDC",
                                                "B", "Plasma",  "Neutrophil", "Eosinophil", "Mast",
                                                "CD4_T", "CD8_T","Treg", "DN_T", "NKT", "NK", "ILC",
                                                "Ciliated", "Club" , "Alv_Epi",
                                                "Artery", "Vein", "Lymphatics", "Cap_A", "Cap_G", "EPC",
                                                "Fibroblast", "SMC", 
                                                "Mesothelial", "Pericyte",
                                                "Megakaryocyte"))
## scVI ####
srt_scvi <- readRDS("~/Documents/scRNAseq/Tool_paper/Data/integrated_scvi.rds")
# remove cells with bad identifications
srt_scvi <- subset(srt_scvi, manual_annotated_1 != "Others")
# set "Condition" as factor
srt_scvi$Condition <- factor(srt_scvi$Condition, levels = c("Control", "Inflam", "Fib", "Fib_O3"))
# set manual_annotated_cell type names as factor
srt_scvi$manual_annotated_1 <- factor(srt_scvi$manual_annotated_1, 
                                     levels = c("AM", "MoDM", "IM", "ncMono", "cMono", 
                                                "MoDC", "cDC1", "cDC2", "Ccr7_DC", "pDC",
                                                "B", "Plasma",  "Neutrophil", "Eosinophil", "Mast",
                                                "CD4_T", "CD8_T","Treg", "DN_T", "NKT", "NK", "ILC",
                                                "Ciliated", "Club" , "Alv_Epi",
                                                "Artery", "Vein", "Lymphatics", "Cap_A", "Cap_G", "EPC",
                                                "Fibroblast", "SMC", 
                                                "Mesothelial", "Pericyte",
                                                "Megakaryocyte"))
                                     
# 3. UMAP label color code ####
# Check duplicates in labels 
duplicated(cell_type_key$TMS_Label[cell_type_key$TMS_Label != "N/A"])
duplicated(cell_type_key$LungMap_Label[cell_type_key$LungMap_Label != "N/A"])
duplicated(cell_type_key$Fibrosis_Label[cell_type_key$Fibrosis_Label != "N/A"])

# Create color code for umap (number of rows for my_color_summary + buplicated labels in Fibrosis label)
set.seed(168)
my_colors_backup <- createPalette(56, seedcolors = c("#F8000D", "#F9D216", "#228833", "#0DCFFE"))
my_colors_summary <- cell_type_key
my_colors_summary$my_colors <- my_colors_backup[1:53]
my_colors_summary$my_colors_tms <- my_colors_summary$my_colors
my_colors_summary$my_colors_lungmap <- my_colors_summary$my_colors
my_colors_summary <- my_colors_summary %>% select(TMS_Label,my_colors_tms, LungMap_Label,my_colors_lungmap ,Fibrosis_Label, my_colors)

# remove color code for missing labels in each categories
my_colors_summary$my_colors[my_colors_summary$Fibrosis_Label == "N/A"] <- NA
my_colors_summary$my_colors_tms[my_colors_summary$TMS_Label == "N/A"] <- NA
my_colors_summary$my_colors_lungmap[my_colors_summary$LungMap_Label == "N/A"] <- NA

# manual edit the duplicated label names in fibrosis color code
my_colors_summary[26,6] <- my_colors_backup[54]
my_colors_summary[27,6] <- NA
my_colors_summary[28,6] <- NA
my_colors_summary[41,6] <- my_colors_backup[55]
my_colors_summary[42,6] <- NA
my_colors_summary[43,6] <- NA
my_colors_summary[44,6] <- NA
my_colors_summary[45,6] <- NA
my_colors_summary[46,6] <- NA
my_colors_summary[47,6] <- my_colors_backup[56]
my_colors_summary[48,6] <- NA

# save the code
saveRDS(my_colors_summary, "~/Documents/scRNAseq/Tool_paper/Data/umap_colors_summary.rds")

# 4. UMAP_log ####
my_colors_summary <- readRDS("~/Documents/scRNAseq/Tool_paper/Data/umap_colors_summary.rds")
# fibrosis color code
my_color_fibrosis <- my_colors_summary %>% select(Fibrosis_Label, my_colors)%>% unique() %>%
  filter(!is.na(my_colors)) %>% column_to_rownames(var = "Fibrosis_Label")
my_color_fibrosis<- setNames(my_color_fibrosis[[1]], rownames(my_color_fibrosis))
## fastmnn (Fig2A&B) ####
p.umap.mnn.condition <- DimPlot(srt_log, reduction = "umap.mnn", 
        group.by = "manual_annotated_1", 
        split.by = "Condition", ncol =2, 
        raster = F, label = T, repel = T,
        label.size = 2.5)+ ggtitle(NULL)+
  scale_color_manual(values = my_color_fibrosis)

p.umap.mnn <- DimPlot(srt_log, reduction = "umap.mnn", 
        group.by = "manual_annotated_1", 
        raster = F, label = T, repel = T,
        label.size = 2.5)+ ggtitle(NULL)+
  scale_color_manual(values = my_color_fibrosis)
p.umap <- (p.umap.mnn + p.umap.mnn.condition) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")
  
## rpca (Supp) ####
p.umap.rpca.condition <- DimPlot(srt_log, reduction = "umap.rpca", 
                         group.by = "manual_annotated_1", 
                         split.by = "Condition", ncol =2, 
                         raster = F, label = T, repel = T,
                         label.size = 2.5)+ ggtitle(NULL)+ 
                         scale_color_manual(values = my_color_fibrosis) 
p.umap.rpca <- DimPlot(srt_log, reduction = "umap.rpca", 
                                group.by = "manual_annotated_1", 
                                raster = F, label = T, repel = T,
                                label.size = 2.5)+ ggtitle(NULL)+
                                scale_color_manual(values = my_color_fibrosis)

(p.umap.rpca + p.umap.rpca.condition) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
rm(p.umap.rpca, p.umap.rpca.condition)

## harmony (Supp) ####
p.umap.harmony.condition <- DimPlot(srt_log, reduction = "umap.harmony", 
        group.by = "manual_annotated_1", 
        split.by = "Condition", ncol =2, 
        raster = F, label = T, repel = T,
        label.size = 2.5)+ ggtitle(NULL)+
  scale_color_manual(values = my_color_fibrosis)
p.umap.harmony <- DimPlot(srt_log, reduction = "umap.harmony", 
                                    group.by = "manual_annotated_1", 
                                    raster = F, label = T, repel = T,
                                    label.size = 2.5)+ ggtitle(NULL)+
  scale_color_manual(values = my_color_fibrosis)

## scVI (Supp) ####
p.umap.scvi.condition <- DimPlot(srt_scvi, reduction = "umap.scvi", 
                                    group.by = "manual_annotated_1", 
                                    split.by = "Condition", ncol =2, 
                                    raster = F, label = T, repel = T,
                                    label.size = 2.5)+ ggtitle(NULL)+
  scale_color_manual(values = my_color_fibrosis)
p.umap.scvi <- DimPlot(srt_scvi, reduction = "umap.scvi", 
                          group.by = "manual_annotated_1", 
                          raster = F, label = T, repel = T,
                          label.size = 2.5)+ ggtitle(NULL)+
  scale_color_manual(values = my_color_fibrosis)
# 5.Bubble (Fig2D) ####
# DE
Idents(srt_log) <- srt_log$manual_annotated_1
srt_log <- JoinLayers(srt_log)
all_markers <- FindAllMarkers(
  srt_log,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)
all_markers_short <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

# plot 
features <- c("Atp6v0d2", "Marco", "Msr1","C1qa",
              "Cx3cr1", "Cd300e", "F13a1", "H2-Aa",
              "Flt3","Xcr1", "Cd209a","Fscn1", "Klk1",
              "Cd79a", "Jchain", "S100a9","Retnlg","Siglecf", 'Ms4a2',
              "Cd3e", "Cd4", "Cd8a", "Foxp3", "Il23r","Nkg7","Xcl1", "Ncr1", "Dach2",
              "Foxj1","Scgb3a2", "Ager", 
              "Cldn5","Bmx","Slc6a2", "Prox1","Ednrb", "Plvap", "Fgfr4",
              "Pdgfra", "Acta2", "Msln", "Trpc6", 
              "Ppbp"
              )

p.bubble <- DotPlot(srt_log, features = features)+ RotatedAxis()
rm(features)

# 6.Correlation heatmap (Fig2C) ####
# calculate means for each cell type
bulk <- AggregateExpression(srt_log, group.by = "manual_annotated_1", return.seurat = TRUE)
df <-bulk[["RNA"]]$data
# calculate correlation
res <- cor(df)
rm(bulk, df)
# convert the cell type names with "_" 
cell_types <- colnames(res)
cell_types <- gsub("-", "_", cell_types)
colnames(res) <- cell_types
rownames(res) <- cell_types
rm(cell_types)
# generate a annotation column for the heatmap label
categories <-data.frame(levels(srt_log$manual_annotated_1)) 
colnames(categories) <- "Fibrosis_Label"
categories <- cell_type_key %>% select(Category, Fibrosis_Label) %>% inner_join(categories, by = "Fibrosis_Label") %>% unique()
categories <- setNames(categories$Category, categories$Fibrosis_Label)
annotation_col <- data.frame(Category = factor(categories, levels = c("Immune", "Epithelial", "Endothelial", "Stromal", "Megakaryocyte")))
rm(categories)
# plot
p.cor.heatmap <- pheatmap(res,
         #annotation_col = annotation_col,
         annotation_row = annotation_col, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete")
rm(annotation_col, res)

# 7. UMAP TMS label (Fig3A) ####
my_color_tms <- my_colors_summary %>% select(TMS_Label, my_colors_tms)%>%unique() %>% 
  filter(!is.na(my_colors_tms)) %>% column_to_rownames(var = "TMS_Label")
my_color_tms<- setNames(my_color_tms[[1]], rownames(my_color_tms))

p.umap.tms <- DimPlot(srt_log, reduction = "umap.mnn", 
        group.by = "SingleR.TMS.cell_type", 
        raster = F, label = T, repel = T,
        label.size = 2.5)+ ggtitle("TMS")+
  scale_color_manual(values = my_color_tms)
rm(my_color_tms)
# 8. UMAP LungMap label (Fig3B) ####
my_color_lungmap <- my_colors_summary %>% select(LungMap_Label, my_colors_lungmap)%>%
  unique() %>% filter(!is.na(my_colors_lungmap)) %>%
  column_to_rownames(var = "LungMap_Label")
my_color_lungmap <- setNames(my_color_lungmap[[1]], rownames(my_color_lungmap))

p.umap.lungmap <-DimPlot(srt_log, reduction = "umap.mnn", 
        group.by = "SingleR.LungMap.celltype_level3", 
        raster = F, label = T, repel = T,
        label.size = 2.5)+ ggtitle("LungMap")+
  scale_color_manual(values = my_color_lungmap)
rm(my_color_lungmap)

# 9. Prediction eval ####
fib_key <- levels(srt_log$manual_annotated_1)
## TMS (Fig3C) ####
## create tms summary
df.tms <- data.frame(
  Cell_Type = character(),
  TP = numeric(),
  TN = numeric(),
  FP = numeric(),
  FN = numeric(),
  stringsAsFactors = FALSE
) 

## calculate for all cell type
for(i in fib_key){
  ## match label in TMS
  tms_label <- cell_type_key %>% filter(Fibrosis_Label == i) %>% pull(TMS_Label)
  if(length(tms_label)>1){ # if matched more than one label 
    if(sum(is.na(tms_label)) == length(tms_label)){ #if all of them are NA 
      tms_label= NA # then return NA
    }else{ #If not
      tms_label <- tms_label[!is.na(tms_label)] # remove na cell type in TMS
    }
  }else{tms_label = tms_label} # if matched just one label then return that label 
  
  if (length(tms_label) == 1) {
    if (is.na(tms_label)) {
      # No prediction for this label
      tp <- 0
      fp <- 0
      fn <- sum(srt_log$manual_annotated_1 == i)
      tn <- ncol(srt_log) - fn
    } else {
      pos_pred_cells <- WhichCells(srt_log, expression = SingleR.TMS.cell_type %in% tms_label)
      neg_pred_cells <- setdiff(Cells(srt_log), pos_pred_cells)
      
      tp <- sum(srt_log$manual_annotated_1[pos_pred_cells] == i, na.rm = T)
      fp <- length(pos_pred_cells) - tp
      
      tn <- sum(srt_log$manual_annotated_1[neg_pred_cells] != i)
      fn <- length(neg_pred_cells) - tn
      
      rm(pos_pred_cells, neg_pred_cells)
    }
  } else {
    pos_pred_cells <- WhichCells(srt_log, expression = SingleR.TMS.cell_type %in% tms_label)
    neg_pred_cells <- setdiff(Cells(srt_log), pos_pred_cells)
    
    tp <- sum(srt_log$manual_annotated_1[pos_pred_cells] == i)
    fp <- length(pos_pred_cells) - tp
    
    tn <- sum(srt_log$manual_annotated_1[neg_pred_cells] != i)
    fn <- length(neg_pred_cells) - tn
    
    rm(pos_pred_cells, neg_pred_cells)
  }
  
  # summary
  # Append row
  df.tms <- rbind(
    df.tms,
    data.frame(
      Cell_Type = i,
      TP = tp,
      TN = tn,
      FP = fp,
      FN = fn,
      stringsAsFactors = FALSE
    )
  )
  rm(fn, fp, tn, tp, tms_label, i)
  
}
df.tms$Total = df.tms$TP + df.tms$FP + df.tms$TN + df.tms$FN
## Precision
df.tms$Precision <- df.tms$TP / (df.tms$TP + df.tms$FP)
##Sensitivity
df.tms$Sensitivity = df.tms$TP / (df.tms$TP + df.tms$FN)
##Specificity
df.tms$Specificity = df.tms$TN / (df.tms$TN + df.tms$FP)
##Accuracy
df.tms$Accuracy = (df.tms$TP + df.tms$TN) / df.tms$Total

saveRDS(df.tms, "~/Documents/scRNAseq/Tool_paper/Data/tms_accuracy.csv")

df.tms <- readRDS("~/Documents/scRNAseq/Tool_paper/Data/tms_accuracy.csv")
df.tms.fig <- df.tms %>% 
  select(Cell_Type, Precision, Sensitivity) %>% 
  melt(id.vars = "Cell_Type",
       variable.name = "Performance_Metrics",
       value.name = "Value")

df.tms.fig$Cell_Type <- factor(df.tms.fig$Cell_Type, levels = fib_key)


p1 <- ggplot(df.tms.fig, aes(x = Cell_Type, y = Value, fill = Performance_Metrics)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("Precision" = "#4E79A7", "Sensitivity" = "#F28E2B"))+
  labs(title = "TMS") +
  theme_classic()+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, margin = margin(t = 5)), 
                    legend.position = "bottom")

rm(df.tms.fig)
write.xlsx(df.tms, "~/Documents/scRNAseq/Tool_paper/Data/tms_singleR_prediction.xlsx")
rm(df.tms)
## LungMap (Fig3D) ####
## create tms summary
df.lungmap <- data.frame(
  Cell_Type = character(),
  TP = numeric(),
  TN = numeric(),
  FP = numeric(),
  FN = numeric(),
  stringsAsFactors = FALSE
) 

## calculate for all cell type
for(i in fib_key){
  ## match label in TMS
  lungmap_label <- cell_type_key %>% filter(Fibrosis_Label == i) %>% pull(LungMap_Label)
  if(length(lungmap_label)>1){ # if matched more than one label 
    if(sum(is.na(lungmap_label)) == length(lungmap_label)){ #if all of them are NA 
      lungmap_label= NA # then return NA
    }else{ #If not
      lungmap_label <- lungmap_label[!is.na(lungmap_label)] # remove na cell type in TMS
    }
  }else{lungmap_label = lungmap_label} # if matched just one label then return that label 
  
  if (length(lungmap_label) == 1) {
    if (is.na(lungmap_label)) {
      # No prediction for this label
      tp <- 0
      fp <- 0
      fn <- sum(srt_log$manual_annotated_1 == i)
      tn <- ncol(srt_log) - fn
    } else {
      pos_pred_cells <- WhichCells(srt_log, expression = SingleR.LungMap.celltype_level3 %in% lungmap_label)
      neg_pred_cells <- setdiff(Cells(srt_log), pos_pred_cells)
      
      tp <- sum(srt_log$manual_annotated_1[pos_pred_cells] == i)
      fp <- length(pos_pred_cells) - tp
      
      tn <- sum(srt_log$manual_annotated_1[neg_pred_cells] != i)
      fn <- length(neg_pred_cells) - tn
      
      rm(pos_pred_cells, neg_pred_cells)
    }
  } else {
    pos_pred_cells <- WhichCells(srt_log, expression = SingleR.LungMap.celltype_level3 %in% lungmap_label)
    neg_pred_cells <- setdiff(Cells(srt_log), pos_pred_cells)
    
    tp <- sum(srt_log$manual_annotated_1[pos_pred_cells] == i)
    fp <- length(pos_pred_cells) - tp
    
    tn <- sum(srt_log$manual_annotated_1[neg_pred_cells] != i)
    fn <- length(neg_pred_cells) - tn
    
    rm(pos_pred_cells, neg_pred_cells)
  }
  
  # summary
  # Append row
  df.lungmap <- rbind(
    df.lungmap,
    data.frame(
      Cell_Type = i,
      TP = tp,
      TN = tn,
      FP = fp,
      FN = fn,
      stringsAsFactors = FALSE
    )
  )
  rm(fn, fp, tn, tp, lungmap_label, i)
  
}
df.lungmap$Total <- df.lungmap$TP + df.lungmap$FP + df.lungmap$TN + df.lungmap$FN
## Precision
df.lungmap$Precision <- df.lungmap$TP / (df.lungmap$TP + df.lungmap$FP)
##Sensitivity
df.lungmap$Sensitivity = df.lungmap$TP / (df.lungmap$TP + df.lungmap$FN)
##Specificity
df.lungmap$Specificity = df.lungmap$TN / (df.lungmap$TN + df.lungmap$FP)
##Accuracy
df.lungmap$Accuracy = (df.lungmap$TP + df.lungmap$TN) / df.lungmap$Total

# save data
saveRDS(df.lungmap, "~/Documents/scRNAseq/Tool_paper/Data/lungmap_accuracy.csv")

# Prepare data for figure
df.lungmap <- readRDS("~/Documents/scRNAseq/Tool_paper/Data/lungmap_accuracy.csv")
df.lungmap.fig <- df.lungmap %>% 
  select(Cell_Type, Precision, Sensitivity) %>% 
  melt(id.vars = "Cell_Type",
       variable.name = "Performance_Metrics",
       value.name = "Value")

df.lungmap.fig$Cell_Type <- factor(df.lungmap.fig$Cell_Type, levels = fib_key)

# plot
p2 <- ggplot(df.lungmap.fig, aes(x = Cell_Type, y = Value, fill = Performance_Metrics)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("Precision" = "#4E79A7", "Sensitivity" = "#F28E2B"))+
  labs(title = "LungMap") +
  theme_classic()+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, margin = margin(t = 5)),
        legend.position = "bottom")

rm(df.lungmap.fig)

write.xlsx(df.lungmap, "~/Documents/scRNAseq/Tool_paper/Data/lungmap_singleR_prediction.xlsx")
rm(df.lungmap)



# 10. UMAP Batch ####
## Top 6 (Fig4A-F) ####
# Log_integrated
p.umap.log.unintegrated <- DimPlot(srt_log, group.by = "Batch", raster = F )+ggtitle("")
# Log_rpca
p.umap.log.rpca <- DimPlot(srt_log, group.by = "Batch", reduction = "umap.rpca", raster = F)+theme(legend.position = "none")+ggtitle("")
# Log_harmony
p.umap.log.harmony <- DimPlot(srt_log, group.by = "Batch", reduction = "umap.harmony", raster = F)+theme(legend.position = "none")+ggtitle("")
# Log_fastmnn
p.umap.log.mnn <- DimPlot(srt_log, group.by = "Batch", reduction = "umap.mnn", raster = F)+theme(legend.position = "none")+ggtitle("")
# sct_unintegreted
p.umap.sct.unintegrated <-DimPlot(srt_sct, group.by = "Batch", raster = F)+theme(legend.position = "none")+ggtitle("")
# scVI 
p.umap.scvi <- DimPlot(srt_scvi, group.by = "Batch", reduction = "umap.scvi", raster = F)+theme(legend.position = "none")+ggtitle("")

## Rest 4 (Supp) ####
# Log_cca
p.umap.log.cca <- DimPlot(srt_log, group.by = "Batch", reduction = "umap.cca", raster = F)+ ggtitle("Log_CCA")
# SCT_cca
p.umap.sct.cca <- DimPlot(srt_sct, group.by = "Batch", reduction = "umap.cca", raster = F)+ ggtitle("SCT_CCA")+theme(legend.position = "none")
# SCT_rpca
p.umap.sct.rpca <- DimPlot(srt_sct, group.by = "Batch", reduction = "umap.rpca", raster = F)+ ggtitle("SCT_RPCA")+theme(legend.position = "none")
# SCT_harmony
p.umap.sct.harmony <- DimPlot(srt_sct, group.by = "Batch", reduction = "umap.harmony", raster = F)+ ggtitle("SCT_Harmony")+theme(legend.position = "none")

rm(srt_log, srt_sct, srt_scvi)

# 11. Integration Evaluation ####
## Eval score distribution (Supp) ####
# score calculation in a seperate file
### Summarize scores ####
#### lisi scores ####
lisi <- readRDS("~/Documents/scRNAseq/Tool_paper/Data/lisi.rds")
integration <- names(lisi)
result <- data.frame()
for(i in integration){
  df <- lisi[[i]] %>% dplyr::select(Batch) %>% rownames_to_column(., var = "Cell")
  names(df)[2] <- "lisi.batch"
  df$Integration <- i
  result <- rbind(result, df)
  rm(i, df)
}
result.summary <- result
rm(result)
# cell score
result <- data.frame()
for(i in integration){
  df <- lisi[[i]] %>% dplyr::select(Cell_Type) %>% rownames_to_column(., var = "Cell")
  names(df)[2] <- "lisi.cell"
  df$Integration <- i
  result <- rbind(result, df)
  rm(i, df)
}
result.summary <- inner_join(result.summary, result, by = c("Cell", "Integration"))
rm(result,lisi )

#### glisi scores ####
glisi <- readRDS("~/Documents/scRNAseq/Tool_paper/Data/glisi.rds")
# batch score
result <- data.frame()
for(i in integration){
  df <- glisi[[i]] %>% dplyr::select(Batch) %>% rownames_to_column(., var = "Cell")
  names(df)[2] <- "glisi.batch"
  df$Integration <- i
  result <- rbind(result, df)
  rm(i, df)
}
result.summary <- inner_join(result.summary, result, by = c("Cell", "Integration"))
rm(result)
# cell type score
result <- data.frame()
for(i in integration){
  df <- glisi[[i]] %>% dplyr::select(Cell_Type) %>% rownames_to_column(., var = "Cell")
  names(df)[2] <- "glisi.cell"
  df$Integration <- i
  result <- rbind(result, df)
  rm(i, df)
}
result.summary <- inner_join(result.summary, result, by = c("Cell", "Integration"))
rm(result, glisi)

#### asw scores ####
#batch score
basw <- readRDS("~/Documents/scRNAseq/Tool_paper/Data/basw.rds")
names(basw)[names(basw) == "Log_scvi"] <- "scvi" # type error from earlier step

result <- data.frame()
for(i in integration){
  df <- basw[[i]] %>% dplyr::select(sil_width) %>% rownames_to_column(., var = "Cell")
  names(df)[2] <- "asw.batch"
  df$Integration <- i
  result <- rbind(result, df)
  rm(i, df)
}
result.summary <- inner_join(result.summary, result, by = c("Cell", "Integration"))
rm(result, basw)
# cell type score
casw <-readRDS("~/Documents/scRNAseq/Tool_paper/Data/casw.rds")
result <- data.frame()
for(i in integration){
  df <- casw[[i]] %>% dplyr::select(sil_width) %>% rownames_to_column(., var = "Cell")
  names(df)[2] <- "asw.cell"
  df$Integration <- i
  result <- rbind(result, df)
  rm(i, df)
}
result.summary <- inner_join(result.summary, result, by = c("Cell", "Integration"))
rm(result, casw)
rm(integration)

result.summary$Integration[result.summary$Integration == "SCT_harmony"] <- "SCT_Harmony"
result.summary$Integration[result.summary$Integration == "scvi"] <- "scVI"
result.summary$Integration[result.summary$Integration == "SCT_unintegrated"] <- "SCT_Unintegrated"
result.summary$Integration[result.summary$Integration == "SCT_rpca"] <- "SCT_RPCA"
result.summary$Integration[result.summary$Integration == "SCT_cca"] <- "SCT_CCA"
result.summary$Integration[result.summary$Integration == "Log_unintegrated"] <- "Log_Unintegrated"
result.summary$Integration[result.summary$Integration == "Log_cca"] <- "Log_CCA"
result.summary$Integration[result.summary$Integration == "Log_harmony"] <- "Log_Harmony"
result.summary$Integration[result.summary$Integration == "Log_mnn"] <- "Log_MNN"
result.summary$Integration[result.summary$Integration == "Log_rpca"] <- "Log_RPCA"
result.summary$Integration <- factor(result.summary$Integration, 
                                     levels = c("Log_RPCA","Log_MNN","Log_Harmony","Log_CCA","Log_Unintegrated",
                                                "SCT_CCA","SCT_RPCA","SCT_Unintegrated","scVI","SCT_Harmony"))
#### plot####
cols_10 <- c("#1F77B4", "#FF7F0E", "#2CA02C","#D62728","#9467BD", "#8C564B",
  "#E377C2","#7F7F7F", "#BCBD22", "#17BECF")

p.basw <- ggplot(result.summary, aes(Integration,asw.batch, fill= Integration))+
  geom_violin(show.legend = F, trim = F)+ 
  geom_boxplot(width=0.1, show.legend = F, fill = "white")+
  scale_fill_manual(values = cols_10) +
  ylab("Scores")+ ggtitle("basw")+
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

p.ilisi <- ggplot(result.summary, aes(Integration,lisi.batch, fill = Integration))+
  geom_violin(show.legend = F, trim = F)+ 
  geom_boxplot(width=0.1, show.legend = F, fill = "white")+
  scale_fill_manual(values = cols_10) +
  ylab("Scores")+ ggtitle("ilisi")+
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

p.gilisi <- ggplot(result.summary, aes(Integration,glisi.batch, fill = Integration))+
  geom_violin(show.legend = F, trim = F)+ 
  geom_boxplot(width=0.1, show.legend = F, fill = "white")+
  scale_fill_manual(values = cols_10) +
  ylab("Scores")+ ggtitle("gilisi")+
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

p.casw <- ggplot(result.summary, aes(Integration,asw.cell, fill = Integration))+
  geom_violin(show.legend = F, trim = F)+ 
  geom_boxplot(width=0.1, show.legend = F, fill = "white")+
  scale_fill_manual(values = cols_10) +
  ylab("Scores")+ggtitle("casw")+
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

p.clisi <- ggplot(result.summary, aes(Integration,lisi.cell, fill = Integration))+
  geom_violin(show.legend = F, trim = F)+ 
  geom_boxplot(width=0.1, show.legend = F, fill = "white")+
  scale_fill_manual(values = cols_10) +
  ylab("Log2 (Scores)")+ggtitle("clisi") +scale_y_continuous(trans = "log2")+
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

p.gclisi <- ggplot(result.summary, aes(Integration,glisi.cell, fill = Integration))+
  geom_violin(show.legend = F, trim = F)+ 
  geom_boxplot(width=0.1, show.legend = F, fill = "white")+
  scale_fill_manual(values = cols_10) +
  ylab("Log2 (Scores)")+ ggtitle("gclisi")+scale_y_continuous(trans = "log2")+
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

rm(result.summary, cols_10)

## Eval Rank Bar (Fig4G) ####
ilisi.summary <-readRDS("~/Documents/scRNAseq/Tool_paper/Data/ilisi.summary.rds")
clisi.summary <- readRDS("~/Documents/scRNAseq/Tool_paper/Data/clisi.summary.rds")
gilisi.summary <- readRDS("~/Documents/scRNAseq/Tool_paper/Data/gilisi.summary.rds")
gclisi.summary <- readRDS("~/Documents/scRNAseq/Tool_paper/Data/gclisi.summary.rds")
basw.summary <-readRDS("~/Documents/scRNAseq/Tool_paper/Data/basw.summary.rds")
casw.summary <-readRDS("~/Documents/scRNAseq/Tool_paper/Data/casw.summary.rds")

eval.summary <- rbind(basw.summary, ilisi.summary, gilisi.summary, 
                      casw.summary, clisi.summary, gclisi.summary)
list <- c("basw", "ilisi", "gilisi", "casw", "clisi", "gclisi")
eval.summary$Method <- rep(list, each = 10)

# ranking scores for each method
rank.summary <- data.frame()
for(i in 1:length(list)){
  rank <- eval.summary %>% filter(Method == list[i]) %>% dplyr::select(Method, Integration, Median) 
  if(list[i] %in% c("ilisi", "gilisi", "casw")){
      rank <- rank %>% arrange(desc(Median))
    } else{
      rank <- rank %>% arrange(Median)
    }
  rank$Rank <- 10:1
  rank.summary <- rbind(rank.summary, rank)
  rm(rank,i)
}

eval.summary <- inner_join(eval.summary, rank.summary, by = c("Method", "Integration", "Median"))
rm(rank.summary, list)

saveRDS(eval.summary, "~/Documents/scRNAseq/Tool_paper/Data/integration_performance_rank.rds")
write.xlsx(eval.summary, "~/Documents/scRNAseq/Tool_paper/Data/integration_performance_rank.xlsx")
rm(basw.summary, casw.summary, clisi.summary, gclisi.summary, gilisi.summary, ilisi.summary)

eval.summary <- readRDS("~/Documents/scRNAseq/Tool_paper/Data/integration_performance_rank.rds")
df <- eval.summary[7:9]
df$Method_group <- ifelse(
  df$Method %in% c("basw", "ilisi", "gilisi"),
  "Batch Mixing",
  "Cell Type Separation"
)
df$Method_group <- factor(df$Method_group,
                          levels = c("Batch Mixing", "Cell Type Separation"))
df$Method <- factor(df$Method, levels= c("basw", "ilisi", "gilisi", 'casw', 'clisi', 'gclisi'))

df$Integration[df$Integration == "SCT_harmony"] <- "SCT_Harmony"
df$Integration[df$Integration == "scvi"] <- "scVI"
df$Integration[df$Integration == "SCT_unintegrated"] <- "SCT_Unintegrated"
df$Integration[df$Integration == "SCT_rpca"] <- "SCT_RPCA"
df$Integration[df$Integration == "SCT_cca"] <- "SCT_CCA"
df$Integration[df$Integration == "Log_unintegrated"] <- "Log_Unintegrated"
df$Integration[df$Integration == "Log_cca"] <- "Log_CCA"
df$Integration[df$Integration == "Log_harmony"] <- "Log_Harmony"
df$Integration[df$Integration == "Log_mnn"] <- "Log_MNN"
df$Integration[df$Integration == "Log_rpca"] <- "Log_RPCA"
df$Integration <- factor(df$Integration, levels = c("SCT_Harmony","scVI","SCT_Unintegrated","SCT_RPCA","SCT_CCA","Log_Unintegrated","Log_CCA","Log_Harmony","Log_MNN","Log_RPCA"))


cols_6 <- c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F","#EDC948")

p.eval.rank <- ggplot(df, aes(Rank, Integration, fill = Method)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.2) +
  facet_grid(. ~ Method_group) +
  scale_fill_manual(values = cols_6) +
  theme_classic() +
  labs(
    fill = "Methods",
    x = "Score Ranks",
    y = "Transformation and Integration"
  )
rm(df, eval.summary, cols_6)


# 12.CCC ####
df <- readRDS("~/Documents/scRNAseq/Tool_paper/Data/liana_all_condition_summary.edit.rds")
## Resource-based analysis ####
### Interaction ####
# Select unique interaction by resources
df_resource <- df %>%
  distinct(source, target, ligand, ligand.complex, receptor, receptor.complex, resource)

# create interaction column
df_venn <- df_resource %>%
  mutate(
    interaction = paste(
      source, target, ligand, ligand.complex,
      receptor, receptor.complex,
      sep = "_"
    )
  ) %>%
  select(resource, interaction)

# convert to wide logical format
df_upset <- df_venn %>%
  distinct() %>%   # important: avoid duplicates
  mutate(value = TRUE) %>%
  pivot_wider(
    names_from  = resource,
    values_from = value,
    values_fill = FALSE
  )

# Make the upset plot with ComplexUpset
resources <- setdiff(colnames(df_upset), "interaction")

p_upset <- upset(
  df_upset,
  intersect = resources,
  n_intersections = 20,
  sort_intersections_by = "cardinality",
  keep_empty_groups = FALSE,
  width_ratio = 0.25,
  base_annotations = list(
    "Intersection size" = intersection_size(
      text = list(size = 2)
    ) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  )
) 

rm(df_resource,df_upset, df_venn, resources)


### Ligand ####
df_venn <- df %>%
  distinct(ligand, resource) # ensure unique rows

# convert to wide logical format
df_upset <- df_venn %>%
  distinct() %>%   # important: avoid duplicates
  mutate(value = TRUE) %>%
  pivot_wider(
    names_from  = resource,
    values_from = value,
    values_fill = FALSE
  )
# Make the upset plot with ComplexUpset
resources <- setdiff(colnames(df_upset), "ligand")
p_upset_l <- upset(
  df_upset,
  intersect = resources,
  n_intersections = 20,
  sort_intersections_by = "cardinality",
  keep_empty_groups = FALSE,
  width_ratio = 0.25,
  base_annotations = list(
    "Intersection size" = intersection_size() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  )
) 

rm(df_upset, df_venn, resources)

### Receptor ####
df_venn <- df %>%
  distinct(receptor, resource) # ensure unique rows

# convert to wide logical format
df_upset <- df_venn %>%
  distinct() %>%   # important: avoid duplicates
  mutate(value = TRUE) %>%
  pivot_wider(
    names_from  = resource,
    values_from = value,
    values_fill = FALSE
  )
# Make the upset plot with ComplexUpset
resources <- setdiff(colnames(df_upset), "receptor")
p_upset_r <- upset(
  df_upset,
  intersect = resources,
  n_intersections = 20,
  sort_intersections_by = "cardinality",
  keep_empty_groups = FALSE,
  width_ratio = 0.25,
  base_annotations = list(
    "Intersection size" = intersection_size() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  )
) 

rm(df_upset, df_venn, resources)


###Ligand_Complex ####
df_venn <- df %>%
  distinct(ligand.complex, resource) %>%
  filter(grepl("_", ligand.complex))  # ensure unique rows

# convert to wide logical format
df_upset <- df_venn %>%
  distinct() %>%   # important: avoid duplicates
  mutate(value = TRUE) %>%
  pivot_wider(
    names_from  = resource,
    values_from = value,
    values_fill = FALSE
  )
# Make the upset plot with ComplexUpset
resources <- setdiff(colnames(df_upset), "ligand.complex")
p_upset_lc <- upset(
  df_upset,
  intersect = resources,
  n_intersections = 20,
  sort_intersections_by = "cardinality",
  keep_empty_groups = FALSE,
  width_ratio = 0.25,
  base_annotations = list(
    "Intersection size" = intersection_size() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  )
) 

rm(df_upset, df_venn, resources)

### Receptor_Complex ####
df_venn <- df %>%
  distinct(receptor.complex, resource) %>%
  filter(grepl("_", receptor.complex)) # ensure unique rows
# convert to wide logical format
df_upset <- df_venn %>%
  distinct() %>%   # important: avoid duplicates
  mutate(value = TRUE) %>%
  pivot_wider(
    names_from  = resource,
    values_from = value,
    values_fill = FALSE
  )
# Make the upset plot with ComplexUpset
resources <- setdiff(colnames(df_upset), "receptor.complex")
p_upset_rc <- upset(
  df_upset,
  intersect = resources,
  n_intersections = 20,
  sort_intersections_by = "cardinality",
  keep_empty_groups = FALSE,
  width_ratio = 0.25,
  base_annotations = list(
    "Intersection size" = intersection_size() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  )
) 

rm(df_upset, df_venn, resources)

## Method-based analysis ####
#### Interaction ####
df_method <- df %>% distinct(source, target, ligand, ligand.complex, receptor, receptor.complex, method)
# create interaction column
df_venn <- df_method %>%
  mutate(
    interaction = paste(
      source, target, ligand, ligand.complex,
      receptor, receptor.complex,
      sep = "_"
    )
  ) %>%
  select(method, interaction)

# convert to wide logical format
df_upset <- df_venn %>%
  distinct() %>%   # important: avoid duplicates
  mutate(value = TRUE) %>%
  pivot_wider(
    names_from  = method,
    values_from = value,
    values_fill = FALSE
  )

# Make the upset plot with ComplexUpset
methods <- setdiff(colnames(df_upset), "interaction")

p_upset_m <- upset(
  df_upset,
  intersect = methods,
  n_intersections = 20,
  sort_intersections_by = "cardinality",
  keep_empty_groups = FALSE,
  width_ratio = 0.25,
  base_annotations = list(
    "Intersection size" = intersection_size(
      text = list(size = 2)
    ) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  )
) 

rm(df_method,df_upset, df_venn, methods)


# 13.Figure Generation ####
## Figure 2 ####
## combine two umaps with shared legend 
p.umap <- (p.umap.mnn + p.umap.mnn.condition)+ 
  plot_layout(guides = "collect", ncol = 2) & 
  theme(legend.position = "right", legend.text = element_text(size = 9)) &
  plot_annotation(tag_levels = 'A')
p.umap
# convert patchwork to ggplot
#p.umap.grob <- ggplotify::as.ggplot(p.umap)
#p.umap.grob 

# convert pheatmap to ggplot
ph_grob<- p.cor.heatmap$gtable
# adjust bubble plot legend
p.bubble <- p.bubble+theme(legend.position = "bottom", legend.title = element_text(size = 11, face = "bold"))
# combine bubble plot and heatmap 
p.bubble.heatmap <- plot_grid(p.bubble, ph_grob, ncol = 2, rel_widths = c(1.3,1), labels = c("C", "D"))
# combine all
p.figure.2 <- plot_grid(p.umap, p.bubble.heatmap, nrow = 2, rel_heights = c(1.3,1))
# save plot 
save_plot("Documents/scRNAseq/Tool_paper/Figures/v3/Figure_2.pdf", 
          p.figure.2, 
          base_width = 19,     # Width in inches (7 inches is ~178mm)
          base_height = 16,    # Adjust height based on your panels
          device = cairo_pdf) # Ensures fonts are embedded

rm(p.umap.mnn, p.umap.mnn.condition, p.umap, p.cor.heatmap, ph_grob,p.bubble, p.bubble.heatmap, p.figure.2)

## Figure 3 ####
# umap
p.umap.tms.1 <- p.umap.tms + theme(legend.position = "bottom",legend.text = element_text(size = 9))
p.umap.lungmap.1 <- p.umap.lungmap + theme(legend.position = "right",legend.text = element_text(size = 9))
p.umap <- (p.umap.tms.1 + p.umap.lungmap.1)+ 
  plot_layout(ncol = 2)+
  plot_annotation(tag_levels = 'A')

p.umap

# prediction bar
p1.1 <- p1 + theme(legend.position = "none")
p2.1 <- p2 + theme(legend.position = "none")
legend <- get_legend(p1)
p.prediction <- plot_grid(
  plot_grid(p1.1, p2.1, ncol = 2, labels = c("C", "D")),
  legend, nrow =2,
  rel_heights = c(1, 0.15)
)

p.figure.3 <- plot_grid(p.umap, p.prediction, nrow = 2, rel_heights = c(1.6,1))
# save plot 
save_plot("Documents/scRNAseq/Tool_paper/Figures/v3/Figure_3.pdf", 
          p.figure.3, 
          base_width = 16,     # Width in inches (7 inches is ~178mm)
          base_height = 14,    # Adjust height based on your panels
          device = cairo_pdf) # Ensures fonts are embedded
rm(p.umap.tms, p.umap.tms.1, p.umap.lungmap, p.umap.lungmap.1, p.umap, p1, p1.1, p2, p2.1, p.prediction, p.figure.3, legend)

## SuppFig rpca ####
p.rpca <- (p.umap.rpca + p.umap.rpca.condition) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
save_plot("Documents/scRNAseq/Tool_paper/Figures/v3/SuppFig_UMAP_rpca.pdf", 
          p.rpca, 
          base_width = 15,     # Width in inches (7 inches is ~178mm)
          base_height = 10.5,    # Adjust height based on your panels
          device = cairo_pdf) # Ensures fonts are embedded
rm(p.umap.rpca, p.umap.rpca.condition, p.rpca)

## SuppFig harmony ####
p.harmony <- (p.umap.harmony + p.umap.harmony.condition) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
save_plot("Documents/scRNAseq/Tool_paper/Figures/v3/SuppFig_UMAP_harmony.pdf", 
          p.harmony, 
          base_width = 15,     # Width in inches (7 inches is ~178mm)
          base_height = 10.5,    # Adjust height based on your panels
          device = cairo_pdf)
rm(p.umap.harmony, p.umap.harmony.condition, p.harmony)

## SuppFig scvi ####
p.scvi <- (p.umap.scvi + p.umap.scvi.condition) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
save_plot("Documents/scRNAseq/Tool_paper/Figures/v3/SuppFig_UMAP_scvi.pdf", 
          p.scvi, 
          base_width = 15,     # Width in inches (7 inches is ~178mm)
          base_height = 10.5,    # Adjust height based on your panels
          device = cairo_pdf)
rm(p.umap.scvi, p.umap.scvi.condition, p.scvi)


## Figure 4 ####
# umap
p.umap <- (p.umap.log.unintegrated + p.umap.sct.unintegrated + p.umap.log.rpca + p.umap.log.harmony +
              p.umap.log.mnn + p.umap.scvi ) + 
  plot_layout(guides = "collect", ncol = 2) +
  plot_annotation(tag_levels = 'A')& 
  theme(legend.position = "bottom")

p.umap.supp <- (p.umap.log.cca + p.umap.sct.rpca + p.umap.sct.harmony +
             p.umap.sct.cca ) + 
  plot_layout(guides = "collect", ncol = 2)& 
  theme(legend.position = "bottom")
save_plot("Documents/scRNAseq/Tool_paper/Figures/v3/SuppFig_UMAP_Integration.pdf", 
          p.umap.supp, 
          base_width = 8,     # Width in inches (7 inches is ~178mm)
          base_height = 8,    # Adjust height based on your panels
          device = cairo_pdf)  

p.violin <- (p.basw+ p.casw + p.ilisi + p.clisi + p.gilisi + p.gclisi)+
  plot_layout(guides = "collect", ncol = 2)

p.integration <- plot_grid(p.umap, plot_grid(p.violin, p.eval.rank, nrow=2, rel_heights = c(1, 0.4), labels = c("G", "H")), 
          ncol = 2, rel_widths = c(1:1) )

save_plot("Documents/scRNAseq/Tool_paper/Figures/v3/Figure_4.pdf", 
          p.integration, 
          base_width = 16,     # Width in inches (7 inches is ~178mm)
          base_height = 13,    # Adjust height based on your panels
          device = cairo_pdf)  
rm(p.umap.log.unintegrated, p.umap.sct.unintegrated, p.umap.log.rpca, 
   p.umap.log.harmony, p.umap.log.mnn, p.umap.scvi,p.umap.log.cca,
   p.umap.sct.rpca, p.umap.sct.harmony, p.umap.sct.cca, p.umap.supp,
   p.violin, p.integration, p.umap)
rm(p.basw, p.ilisi, p.gilisi, p.casw, p.clisi, p.gclisi, p.eval.rank)

## Figure 5 ####
p1 <- wrap_elements(p_upset)
p2 <- wrap_elements(p_upset_m)

p.upset <-(p1 + p2) +
  plot_layout(ncol = 1, guides = "collect", heights = c(1.1, 0.7))+
  plot_annotation(tag_levels = 'A')
save_plot("Documents/scRNAseq/Tool_paper/Figures/v3/Figure_5.pdf", 
          p.upset, 
          base_width = 12,     # Width in inches (7 inches is ~178mm)
          base_height = 14,    # Adjust height based on your panels
          device = cairo_pdf)  
rm(p1, p2, p_upset, p_upset_m, p.upset)

p3 <- wrap_elements(p_upset_l)
p4 <- wrap_elements(p_upset_r)
p5 <- wrap_elements(p_upset_lc) 
p6 <- wrap_elements(p_upset_rc)

p.upset.s1 <-(p3 + p5)+
  plot_layout(guides = "collect", ncol = 1, nrow = 2, heights = c(1.3, 0.6))+
  plot_annotation(tag_levels = 'A')
save_plot("Documents/scRNAseq/Tool_paper/Figures/v3/SuppFig_UpSet_ligand.pdf", 
          p.upset.s1, 
          base_width = 10,     # Width in inches (7 inches is ~178mm)
          base_height = 16,    # Adjust height based on your panels
          device = cairo_pdf)  
rm(p3, p5, p_upset_l, p_upset_lc, p.upset.s1)

p.upset.s2<- (p4+ p6)+
  plot_layout(guides = "collect", ncol = 1, nrow = 2, heights = c(1.3, 0.6))+
  plot_annotation(tag_levels = 'A')

save_plot("Documents/scRNAseq/Tool_paper/Figures/v3/SuppFig_UpSet_receptor.pdf", 
          p.upset.s2, 
          base_width = 10,     # Width in inches (7 inches is ~178mm)
          base_height = 16,    # Adjust height based on your panels
          device = cairo_pdf)  
rm(p4, p6, p_upset_r, p_upset_rc, p.upset.s2)




################

install.packages("ggplotify")
ph.gg <- as_ggplot(p.cor.heatmap$gtable)

(p.bubble +ph.gg) + plot_layout()

plot_grid(p.bubble, ph.gg)
library(patchwork)
library(ggplotify)
library(cowplot)

plots_list <- patchwork::wrap_plots(p.umap.mnn)$plots
grobs <- lapply(plots_list, ggplotify::as.ggplot)

# Combine with cowplot
plot_grid(p.bubble, do.call(plot_grid, c(grobs, ncol = length(grobs))), ncol = 2)


# Read PDFs
p1 <- ggdraw() +draw_image(image_read_pdf("~/Downloads/rstudio-export-2/Fig2A&B.pdf"))
p2 <- ggdraw() +draw_image(image_read_pdf("~/Downloads/rstudio-export-2/Fig2C&D.pdf"))
final_plot <- plot_grid(p1, p2, 
                        #labels = "AUTO", # Automatically adds A, B
                        rel_heights = c(1, 1),
                        nrow = 2)
print(final_plot)
save_plot("~/Downloads/rstudio-export-2/Figure_2.pdf", final_plot)