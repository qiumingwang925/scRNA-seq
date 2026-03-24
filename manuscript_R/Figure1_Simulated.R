# Load library ####
library(tidyverse)
library(pheatmap)
library(igraph)
library(circlize)
library(ComplexHeatmap)
library(ggraph)
library(ggplot2)
library(dplyr)
library(ggalluvial)

#set.seed(123)
set.seed(42)

# 1. Integration ####
## Before integration ####
### Define centers for each batch ####
batch_A_centers <- matrix(c(0.3, -0.8,
                            1.5, 2.3), 
                          ncol = 2, byrow = TRUE)
batch_B_centers <- matrix(c(-2, 3.5,
                            -3, 0.5,
                            -2.5, -2), 
                          ncol = 2, byrow = TRUE)

### Function to generate cluster data ####
generate_cluster_data <- function(centers, batch_label, start_cluster) {
  data_list <- list()
  
  for (i in seq_len(nrow(centers))) {
    x <- rnorm(150, mean = centers[i, 1], sd = 0.5)
    y <- rnorm(150, mean = centers[i, 2], sd = 0.5)
    cluster <- start_cluster + i - 1
    data_list[[i]] <- tibble(x = x, y = y, batch = batch_label, cluster = cluster)
  }
  
  bind_rows(data_list)
}

### Generate Batch A and Batch B ####
data_A <- generate_cluster_data(batch_A_centers, "A", 1)
data_B <- generate_cluster_data(batch_B_centers, "B", 3)

### Combine into one dataframe ####
umap_data <- bind_rows(data_A, data_B, )

### Preview ####
head(umap_data)


### Visualize ####
ggplot(umap_data, aes(x, y, color = factor(cluster), shape = batch)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_classic() + xlab("UMAP1")+ylab("UMAP2")+ theme(legend.position = "none")

p1<- ggplot(umap_data, aes(x, y, color = batch)) +
  geom_point(size = 2, alpha = 0.8, ) +
  theme_classic() + xlab("UMAP1")+ylab("UMAP2")+ theme(legend.position = "none")

### Clean up ####
rm(batch_A_centers, batch_B_centers, data_A, data_B, generate_cluster_data, umap_data )

## After integration ####
### Define centers for each batch ####
batch_A_centers <- matrix(c(1.7, -1.2,
                            -1.3, -1.8), 
                          ncol = 2, byrow = TRUE)
batch_B_centers <- matrix(c(1.7, -1.2,
                            -1.3, -1.8,
                            0.1, 2), 
                          ncol = 2, byrow = TRUE)

### Function to generate cluster data ####
generate_cluster_data <- function(centers, batch_label, start_cluster) {
  data_list <- list()
  
  for (i in seq_len(nrow(centers))) {
    x <- rnorm(150, mean = centers[i, 1], sd = 0.5)
    y <- rnorm(150, mean = centers[i, 2], sd = 0.5)
    cluster <- start_cluster + i - 1
    data_list[[i]] <- tibble(x = x, y = y, batch = batch_label, cluster = cluster)
  }
  
  bind_rows(data_list)
}

### Generate Batch A and Batch B ####
data_A <- generate_cluster_data(batch_A_centers, "A", 1)
data_B <- generate_cluster_data(batch_B_centers, "B", 3)

### Combine into one dataframe ####
umap_data <- bind_rows(data_A, data_B)

### Preview ####
head(umap_data)

### Visualize ####
ggplot(umap_data, aes(x, y, color = batch)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_classic() + xlab("UMAP1")+ylab("UMAP2")+ theme(legend.position = "none")

p2 <- ggplot(umap_data, aes(x, y, color = batch)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_classic() + xlab("UMAP1")+ylab("UMAP2")+ theme(legend.position = "none")

### Clean up ####
rm(batch_A_centers, batch_B_centers, data_A, data_B, generate_cluster_data, umap_data )

## Combine before and after figures  ####
p1|p2

# 2. Integration eval ####
## Generate a dataframe ####
inte_eval <- data.frame("Methods" = c("Method 1", "Method 2", "Method 3", "Method 4", "Method 5"),
                        "Batch_Mixing" = 5:1,
                        "Cell_Type_Separation" = c(2,3,5,4,1 ))
## Convert to long format for ggplot ####
inte_eval_long <- inte_eval %>%
  pivot_longer(cols = c(Batch_Mixing, Cell_Type_Separation),
               names_to = "Metric",
               values_to = "Score")

## Visualize ####
ggplot(inte_eval_long, aes(x = Methods, y = Score, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_classic(base_size = 14) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        legend.position = "bottom") +
  labs(y = "Score") +
  scale_fill_brewer(palette = "Set2")

## Clean up ####
rm(inte_eval, inte_eval_long)

# 3. Bubble plot ####
## Simulate expression data ####
genes <- paste0("Gene", 1:5)
cell_types <- c("Cell Type 1", "Cell Type 2", "Cell Type 3", "Cell Type 4")

expr_data <- expand.grid(Gene = genes, CellType = cell_types) %>%
  mutate(Expression = runif(n(), 0, 10),
         PctExpr = runif(n(), 0.1, 1))

## Visualize ####
ggplot(expr_data, aes(x = CellType, y = Gene, size = PctExpr, color = Expression)) +
  geom_point(alpha = 0.8) +
  scale_color_viridis_c(option = "plasma") +
  scale_size(range = c(3, 10)) +
  theme_classic(base_size = 14) +
  labs(x = "Cell Type", y = "Gene") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

## Clean up ####
rm(genes, expr_data)

# 4. Violin plot ####
## Simulate single-cell expression values ####
violin_data <- data.frame(
  CellType = rep(cell_types, each = 100),
  Expression = c(
    rnorm(100, 5, 1),
    rnorm(100, 6, 1.5),
    rnorm(100, 3, 0.8),
    rnorm(100, 7, 1)
  )
)

## Visualize ####
ggplot(violin_data, aes(x = CellType, y = Expression, fill = CellType)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, color = "black", alpha = 0.7) +
  theme_classic(base_size = 14) +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Expression Level")

## Clean up ####
rm(violin_data)

# 5. Heatmap ####
## Create a matrix ####
heatmap_matrix <- expr_data %>%
  select(Gene, CellType, Expression) %>%
  pivot_wider(names_from = CellType, values_from = Expression) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

## Visualize ####
pheatmap(heatmap_matrix,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         cluster_rows = TRUE, cluster_cols = TRUE)

## Clean up ####
rm(heatmap_matrix)

# 6. Enrichment bar graph ####
## Enrichment data ####
enrich_df <- data.frame(
  Term = c(
    "Pathway1", "Pathway2", "Pathway3", "Pathway4", 
    "Pathway5", "Pathway6", "Pathway7", "Pathway8"
  ),
  pvalue = runif(8, 1e-6, 0.05),
  GeneRatio = runif(8, 0.1, 0.6)
)

# Add -log10(pvalue)
enrich_df <- enrich_df %>%
  mutate(logP = -log10(pvalue)) %>%
  arrange(logP)  # sort for plotting

## Visualize ####
ggplot(enrich_df, aes(x = logP, y = reorder(Term, logP), fill = GeneRatio)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_c(option = "plasma") +
  labs(
    x = expression(-log[10](pvalue)),
    y = NULL,
    fill = "Gene Ratio",
    title = "Enrichment Analysis"
  ) +
  theme_classic(base_size = 14) 

## Clean up ####
rm(enrich_df)

# 7. Cell-cell communication ####
## Simulate cell types ####
cell_types <- c("Cell Type 1", "Cell Type 2", "Cell Type 3", "Cell Type 4", "Cell Type 5")

## Simulate interaction strengths ####
interaction_matrix <- matrix(runif(25, 0, 1), nrow = 5)
rownames(interaction_matrix) <- cell_types
colnames(interaction_matrix) <- cell_types

interaction_df <- as.data.frame(as.table(interaction_matrix))
colnames(interaction_df) <- c("Sender", "Receiver", "Strength")
#devtools::install_github("thomasp85/ggraph")

## Visualize (Circle) ####
g <- graph_from_data_frame(interaction_df, directed = TRUE)

ggraph(g, layout = "circle") +
  geom_edge_link(aes(width = Strength, color = Strength), alpha = 0.8) +
  geom_node_point(size = 10, color = "grey90") +
  geom_node_text(aes(label = name), size = 4) +
  scale_edge_colour_gradient(low = "skyblue", high = "firebrick3") +
  theme_void() 

## Visualize (Heatmap) ####
Heatmap(
  interaction_matrix,
  name = "Interaction Strength",
  col = colorRampPalette(c("white", "pink", "firebrick3"))(50),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_title = "Receiver Cell Type",
  row_title = "Sender Cell Type"
)

# Prepare data
interaction_df$Sender <- factor(interaction_df$Sender, levels = cell_types)
interaction_df$Receiver <- factor(interaction_df$Receiver, levels = cell_types)

## Visualize (Chord) ####
circos.clear()
chordDiagram(
  x = interaction_df,
  grid.col = setNames(c("#E74C3C", "#3498DB", "#2ECC71", "#F1C40F", "#9B59B6"), cell_types),
  transparency = 0.4
)

## Clean up ####
rm(g, interaction_df, interaction_matrix, cell_types)

## Alluvial plots ####
set.seed(123)

### Simulate Cell group -> Pattern ####
left_data <- data.frame(
  cell_group = rep(c("C1", "C2", "C3", "C4", "C5"), each = 3),
  pattern = rep(c("P1", "P2", "P3"), times = 5),
  strength = runif(15, 1, 10)
)

left_data$strength <- left_data$strength / sum(left_data$strength)

### Simulate Pattern -> Pathway ####
right_data <- data.frame(
  pattern = rep(c("P1", "P2", "P3"), each = 8),
  pathway = rep(paste0("Pathway", 1:8), times = 3),
  strength = runif(24, 1, 10)
)

right_data$strength <- right_data$strength / sum(right_data$strength)

### Common color palette for patterns ####
pattern_colors <- c("P1" = "#66c2a5", "P2" = "#fc8d62", "P3" = "#8da0cb")

### Visualize (left) ####
p1 <- ggplot(left_data,
             aes(axis1 = cell_group, axis2 = pattern, y = strength)) +
  geom_alluvium(aes(fill = pattern), width = 1/12, alpha = 0.8) +
  geom_stratum(width = 1/12, fill = "grey90", color = "grey40") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_fill_manual(values = pattern_colors) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"   # 🔥 removes legend
  ) +
  scale_x_discrete(limits = c("Cell Group", "Pattern"), expand = c(.05, .05))

### Visualize (right) ####
p2 <- ggplot(right_data,
             aes(axis1 = pattern, axis2 = pathway, y = strength)) +
  geom_alluvium(aes(fill = pattern), width = 1/12, alpha = 0.8) +
  geom_stratum(width = 1/12, fill = "grey90", color = "grey40") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_fill_manual(values = pattern_colors) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"   # 🔥 removes legend
  ) +
  scale_x_discrete(limits = c("Pattern", "Pathway"), expand = c(.05, .05))

### Combine side-by-side ####
p1 + p2

### Clean up####
rm(left_data, right_data, pattern_colors, p1, p2)