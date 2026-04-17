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
##  Ask user to select one for group/treatment. If it's not called "group", than save as colname "group"
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

library(tidyverse)
library(ggplot2)
library(patchwork)
library(magrittr)
library(liana)
library(Seurat)
library(reshape2)
# -- mod_interact_liana_comp --
## use the .rds seurat object processed in mod_interact_upload as input
## subset the obj by group for LIANA cell-cell communication 
## CCC resource selection resources <- show_resources()[2:19] (single selection to avoid large memory use)
### MouseConsensus: no gene name conversion needed
#### liana_result <- liana_wrap(sce = obj_sub, assay = "RNA"/"SCT", slot = "data", method = methods, resource =  c("MouseConsensus"))
### other resources (human-based resources): gene name conversion needed
#### use following code as reference to convert mouse gene in the dataset to human gene before liana_wrap()
# Human Gene Symbol Conversion
#BiocManager::install("biomaRt")
library(biomaRt)
# Connect to Ensembl
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
# Extract mouse gene list from data set
mouse_genes <- rownames(obj_sub)  
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
gene_map <- mapped_genes[match(rownames(obj_sub), mapped_genes$MGI.symbol), ]
# Check how many genes have NA names
sum(is.na(gene_map$MGI.symbol))
sum(is.na(gene_map$HGNC.symbol))
# Replace rownames with human gene symbols
rownames(obj_sub) <- gene_map$HGNC.symbol
# Check how many genes have NA names
sum(is.na(rownames(obj_sub)))
# Remove genes with NA names
obj_sub <- obj_sub[!is.na(rownames(obj_sub)), ]
#### use following code  as refences to convert the human gene in liana result (ligand, Ligand complex, receptor, receptor complex) to mouse gene
# make a human gene list to convert
ligand.complex <- unique(df_human$ligand.complex) 
ligand.complex.list <- unlist(strsplit(ligand.complex, "_"))
receptor.complex <- unique(df_human$receptor.complex) 
receptor.complex.list <- unlist(strsplit(receptor.complex, "_"))
human_genes <- unique(c(df_human$ligand, df_human$receptor, ligand.complex.list, receptor.complex.list))
rm(ligand.complex.list, receptor.complex.list)

# Mouse gene match its human ortholog
human_to_mouse <- getLDS(
  attributes = c("hgnc_symbol"),       # FROM human
  filters = "hgnc_symbol",
  values = human_genes,
  mart = human,                        # FROM this mart
  attributesL = c("mgi_symbol"),       # TO mouse
  martL = mouse,                       # TO this mart
  uniqueRows = TRUE
)
#### Reference: ligand/receptor ####
# Remove genes with no human ortholog
mapped_genes <- human_to_mouse
mapped_genes <- mapped_genes[!mapped_genes$MGI.symbol == "", ]
mapped_genes <- mapped_genes[!duplicated(mapped_genes), ]
# Remove one mouse gene with multiple human ortholog
mapped_genes <- mapped_genes[!duplicated(mapped_genes$MGI.symbol), ]
mapped_genes$MGI.symbol[duplicated(mapped_genes$MGI.symbol)]
# Remove multiple mouse genes with one human ortholog
mapped_genes <- mapped_genes[!duplicated(mapped_genes$HGNC.symbol), ]
mapped_genes$HGNC.symbol[duplicated(mapped_genes$HGNC.symbol)]
# Filter to genes in your data
gene_map <- mapped_genes[match(human_genes, mapped_genes$HGNC.symbol), ]
# Check how many genes have NA names
sum(is.na(gene_map$MGI.symbol))
sum(is.na(gene_map$HGNC.symbol))
# remove na
gene_map <- na.omit(gene_map)
# check duplicated
sum(duplicated(gene_map$MGI.symbol))
sum(duplicated(gene_map$HGNC.symbol))

#### Reference: complex #### 
ligand.complex.df <- data.frame(complex = ligand.complex[grepl("_", ligand.complex)]) %>%
  mutate(ligand = complex) %>% 
  separate(ligand, into = c("gene1", "gene2", "gene3", "gene4"), sep = "_", fill = "right")
receptor.complex.df <- data.frame(complex =  receptor.complex[grepl("_",  receptor.complex)]) %>%
  mutate( receptor = complex) %>% 
  separate( receptor, into = c("gene1", "gene2", "gene3", "gene4"), sep = "_", fill = "right")
complex.df <- rbind(ligand.complex.df, receptor.complex.df)%>% distinct()
rm(ligand.complex.df, receptor.complex.df) 

complex.df <- melt(complex.df, id.vars = "complex", variable.name = "gene_order", value.name = "gene") %>%
  arrange(complex)
complex.df <- complex.df %>% 
  left_join(.,gene_map, by = c("gene" = "HGNC.symbol") ) %>%
  mutate(gene = ifelse(!is.na(MGI.symbol), MGI.symbol, gene)) %>% 
  dplyr::select(-MGI.symbol)
complex.df <-dcast(complex.df, complex ~ gene_order, fun.aggregate = function(x) paste(x, collapse = ", "))
complex.df[complex.df == "NA"] <- NA
complex.df$complex.edit <- apply(complex.df[, c("gene1", "gene2", "gene3", "gene4")], 1, function(x) {
  paste(na.omit(x), collapse = "_")
})
## CCC methods selction methods <- show_methods()[c(1,2,3,4,5,6, 8, 11)] (multiple selection)
### methods need extra packages
#### cytotalk: install.packages("entropy")
#### cellchat: remotes::install_github("sqjin/CellChat")
#### italk: remotes::install_github("Coolgenome/iTALK")
## aggregate the result from one resource and multiple methods
### show summary in tables
## add a download button to download the result as excel

# -- mod_interact_liana_vis --
## add a data upload
## reference UI 
### subtab 1 
tabPanel("CCC Dot Plot", 
         fluidRow(
           column(2, selectInput("ccc.dot.group", "group", choices = c(levels(srt$Treatment)))),
           column(3, selectInput("ccc.dot.source", "Source Populations (Select One or More)", choices = c(levels(Idents(srt))), multiple = T)),
           column(3, selectInput("ccc.dot.target", "Targted Populations (Select One or More)", choices = c(levels(Idents(srt))), multiple = T)),
           column(2, numericInput("ccc.dot.top", "Number of Top genes", value = 15)),
           column(2, actionButton("run.ccc.dot", "Run", class = "btn-primary" ))
         ), 
         fluidRow(
           column(4, selectInput("ccc.dot.col", "Dot Color Represention", choices = c("Interaction Specificity(NATMI)", 
                                                                                      "Interaction Weight(Connectome)", 
                                                                                      "LogFC Mean(iTALK)",
                                                                                      "Expression Magnitude(SingleCellSignalR)",
                                                                                      "P Value(CellPhoneDB)"),
                                 selected = "Expression Magnitude(SingleCellSignalR)")),
           column(4, selectInput("ccc.dot.size", "Dot Size Represention", choices = c("Interaction Specificity(NATMI)", 
                                                                                      "Interaction Weight(Connectome)", 
                                                                                      "LogFC Mean(iTALK)",
                                                                                      "Expression Magnituden/(SingleCellSignalR)",
                                                                                      "P Value(CellPhoneDB)"),
                                 selected = "Interaction Specificity(NATMI)"))
         ),
         fluidRow(
           column(12,  plotOutput("plot.ccc.dot", height = "600px", width = "1400px"))
         ),
         fluidRow(
           ### Download figure size
           column(3, numericInput("p.ccc.dot.w", "Width (Download Figure)", value = 8)),
           column(3, numericInput("p.ccc.dot.h", "Height (Download Figure)", value = 6))
         ),
         fluidRow(
           ### Download action button
           column(6, downloadButton("p.ccc.dot", "Download Figure", class="btn-success"))
         )
         
)
### subtab 2
tabPanel("CCC Freq Heatmap", 
         fluidRow(
           column(8, selectInput("ccc.freqheat", "Populations (Select One or More)", choices = c(levels(Idents(srt))), multiple = T))
         ), 
         fluidRow(
           column(6, selectInput("ccc.freqheat.group.1", "group(Left View)", choices = c(levels(srt$Treatment)))),
           column(4, selectInput("ccc.freqheat.group.2", "group(Right View)", choices = c(levels(srt$Treatment)))),
           column(2, actionButton("run.ccc.freqheat", "Run", class = "btn-primary" ))
         ),
         fluidRow(
           column(6,  plotOutput("plot.ccc.freqheat.1", height = "600px", width = "700px")),
           column(6,  plotOutput("plot.ccc.freqheat.2", height = "600px", width = "700px"))
         ),
         fluidRow(
           ### Download figure size
           column(3, numericInput("p.ccc.freqheat.1.w", "Width (Download Figure)", value = 7)),
           column(3, numericInput("p.ccc.freqheat.1.h", "Height (Download Figure)", value = 6)),
           column(3, numericInput("p.ccc.freqheat.2.w", "Width (Download Figure)", value = 7)),
           column(3, numericInput("p.ccc.freqheat.2.h", "Height (Download Figure)", value = 6))
         ),
         fluidRow(
           ### Download action button
           column(6, downloadButton("p.ccc.freqheat.1", "Download Figure", class="btn-success")),
           column(6, downloadButton("p.ccc.freqheat.2", "Download Figure", class="btn-success"))
         )
         
)

### subtab3
tabPanel("CCC Freq Chord Diagram", 
         fluidRow(
           column(3, selectInput("ccc.freqchord.source", "Source Populations (Select One or More)", choices = c(levels(Idents(srt))), multiple = T)),
           column(3, selectInput("ccc.freqchord.target", "Targted Populations (Select One or More)", choices = c(levels(Idents(srt))), multiple = T))
         ), 
         fluidRow(
           column(6, selectInput("ccc.freqchord.group.1", "group(Left View)", choices = c(levels(srt$Treatment)))),
           column(4, selectInput("ccc.freqchord.group.2", "group(Right View)", choices = c(levels(srt$Treatment)))),
           column(2, actionButton("run.ccc.freqchord", "Run", class = "btn-primary" ))
         ),
         fluidRow(
           column(6,  plotOutput("plot.ccc.freqchord.1", height = "400px", width = "400px")),
           column(6,  plotOutput("plot.ccc.freqchord.2", height = "400px", width = "400px"))
         ),
         fluidRow(
           ### Download figure size
           column(3, numericInput("p.ccc.freqchord.1.w", "Width (Download Figure)", value = 6)),
           column(3, numericInput("p.ccc.freqchord.1.h", "Height (Download Figure)", value = 6)),
           column(3, numericInput("p.ccc.freqchord.2.w", "Width (Download Figure)", value = 6)),
           column(3, numericInput("p.ccc.freqchord.2.h", "Height (Download Figure)", value = 6))
         ),
         fluidRow(
           ### Download action button
           column(6, downloadButton("p.ccc.freqchord.1", "Download Figure", class="btn-success")),
           column(6, downloadButton("p.ccc.freqchord.2", "Download Figure", class="btn-success"))
         )
         
)
## reference server
## CCC dot plot
plotInput.ccc.dot <- eventReactive(input$run.ccc.dot, {
  if(input$ccc.dot.size =="P Value(CellPhoneDB)"){
    size.rep <- "cellphonedb.pvalue"
    size.lab <- "-log10(p_value)"
    size.invert <- TRUE
  }else if(input$ccc.dot.size =="Interaction Specificity(NATMI)"){
    size.rep <- "natmi.edge_specificity"
    size.lab <- "Interaction\nSpecificity"
    size.invert <- FALSE
  }else if(input$ccc.dot.size == "Interaction Weight(Connectome)"){
    size.rep <- "connectome.weight_sc"
    size.lab <- "Interaction\nWeight"
    size.invert <- FALSE
  }else if(input$ccc.dot.size == "LogFC Mean(iTALK)"){
    size.rep <- "logfc.logfc_comb"
    size.lab <- "LogFC Mean"
    size.invert <- FALSE
  }else{
    size.rep <- "sca.LRscore"
    size.lab <- "Expression\nMagnitude"
    size.invert <- FALSE
  }
  
  if(input$ccc.dot.col =="P Value(CellPhoneDB)"){
    col.rep <- "cellphonedb.pvalue"
    col.lab <- "-log10(p_value)"
  }else if(input$ccc.dot.col =="Interaction Specificity(NATMI)"){
    col.rep <- "natmi.edge_specificity"
    col.lab <- "Interaction\nSpecificity"
  }else if(input$ccc.dot.col == "Interaction Weight(Connectome)"){
    col.rep <- "connectome.weight_sc"
    col.lab <- "Interaction\nWeight"
  }else if(input$ccc.dot.col == "LogFC Mean(iTALK)"){
    col.rep <- "logfc.logfc_comb"
    col.lab <- "LogFC Mean"
  }else{
    col.rep <- "sca.LRscore"
    col.lab <- "Expression\nMagnitude"
  }
  
  liana_dotplot(
    liana %>% filter(group == input$ccc.dot.group),
    source_groups = input$ccc.dot.source,
    target_groups = input$ccc.dot.target,
    ntop = input$ccc.dot.top,
    specificity = size.rep,
    magnitude = col.rep,
    y.label = "Interactions (Ligand -> Receptor)",
    size.label = size.lab,
    colour.label = col.lab,
    show_complex = TRUE,
    size_range = c(2, 10),
    invert_specificity = TRUE,
    invert_magnitude = size.invert,
    invert_function = function(x) -log10(x + 1e-10)
  )+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, color = "black", size = 15))
})

output$plot.ccc.dot <- renderPlot({
  print(plotInput.ccc.dot())
}, res = 96)

#download figure
output$p.ccc.dot <- downloadHandler(
  filename = function(){"Liana_CCC_dot.jpg"},
  content = function(file){
    ggsave(file, plot = plotInput.ccc.dot(), width = input$p.ccc.dot.w, height =  input$p.ccc.dot.h)
  }
)
# ccc freq heatmap
## left view
plotInput.ccc.freqheat.1 <- eventReactive(input$run.ccc.freqheat, {
  liana %>%
    filter(group == input$ccc.freqheat.group.1) %>% 
    # only keep interactions concordant between methods
    filter(aggregate_rank <= 0.01 & source %in% input$ccc.freqheat & target %in% input$ccc.freqheat) %>%  # note that these pvals are already corrected
    heat_freq()
})

output$plot.ccc.freqheat.1 <- renderPlot({
  print(plotInput.ccc.freqheat.1())
}, res = 96)

#download figure
output$p.ccc.freqheat.1 <- downloadHandler(
  filename = function(){paste0("Liana_CCC_freqheat_",input$ccc.freqheat.group.1,".jpg")},
  content = function(file){
    jpeg(file, width = input$p.ccc.freqheat.1.w, height =  input$p.ccc.freqheat.1.h, units = "in", res = 300)
    print(plotInput.ccc.freqheat.1())
    dev.off()
  }
)
## right view
plotInput.ccc.freqheat.2 <- eventReactive(input$run.ccc.freqheat, {
  liana %>%
    filter(group == input$ccc.freqheat.group.2) %>% 
    # only keep interactions concordant between methods
    filter(aggregate_rank <= 0.01 & source %in% input$ccc.freqheat & target %in% input$ccc.freqheat) %>%  # note that these pvals are already corrected
    heat_freq()
})

output$plot.ccc.freqheat.2 <- renderPlot({
  print(plotInput.ccc.freqheat.2())
}, res = 96)

#download figure
output$p.ccc.freqheat.2 <- downloadHandler(
  filename = function(){paste0("Liana_CCC_freqheat_",input$ccc.freqheat.group.2,".jpg")},
  content = function(file){
    jpeg(file, width = input$p.ccc.freqheat.2.w, height =  input$p.ccc.freqheat.2.h, units = "in", res = 300)
    print(plotInput.ccc.freqheat.2())
    dev.off()
  }
)


# ccc freqchord
## left view
plotInput.ccc.freqchord.1 <- eventReactive(input$run.ccc.freqchord, {
  liana %>%
    filter(group == input$ccc.freqchord.group.1) %>% 
    # only keep interactions concordant between methods
    filter(aggregate_rank <= 0.01 ) %>%  # note that these pvals are already corrected
    chord_freq(., source_groups = input$ccc.freqchord.source, target_groups = input$ccc.freqchord.target)
})

output$plot.ccc.freqchord.1 <- renderPlot({
  print(plotInput.ccc.freqchord.1())
}, res = 96)

#download figure
output$p.ccc.freqchord.1 <- downloadHandler(
  filename = function(){paste0("Liana_CCC_freqchord_",input$ccc.freqchord.group.1,".jpg")},
  content = function(file){
    jpeg(file, width = input$p.ccc.freqchord.1.w , height =  input$p.ccc.freqchord.1.h , units = "in", res = 300)
    print(plotInput.ccc.freqchord.1())
    dev.off()
  }
)
## right view
plotInput.ccc.freqchord.2 <- eventReactive(input$run.ccc.freqchord, {
  liana %>%
    filter(group == input$ccc.freqchord.group.2) %>%
    # only keep interactions concordant between methods
    filter(aggregate_rank <= 0.01 ) %>%  # note that these pvals are already corrected
    chord_freq(., source_groups = input$ccc.freqchord.source, target_groups = input$ccc.freqchord.target)
})

output$plot.ccc.freqchord.2 <- renderPlot({
  print(plotInput.ccc.freqchord.2())
}, res = 96)

#download figure
output$p.ccc.freqchord.2 <- downloadHandler(
  filename = function(){paste0("Liana_CCC_freqchord_",input$ccc.freqchord.group.2,".jpg")},
  content = function(file){
    jpeg(file, width = input$p.ccc.freqchord.2.w , height =  input$p.ccc.freqchord.2.h , units = "in", res = 300)
    print(plotInput.ccc.freqchord.2())
    dev.off()
  }
)

