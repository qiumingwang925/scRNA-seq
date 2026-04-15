library(Seurat)

# --- mod_explore_upload ---
# upload .rds seuart object from local. The obj should be a slim down version:
##  only "data" slot in RNA or SCT in assays
##  meta.data contains nCounts or more to scale data for heatmap 
##  treatment/condition should store in "group" of metadata
##  keep one pca and one umap from any interation or unintegrated data
##  matched graph slot should be exsit or disables re-umap for subset
##  active.ident needs to be cell type labels using for later analyses
obj <- readRDS("~/Documents/Research/UCDacis_Venosa/Shiny/data/Annotated_O3_fibrosis_data_slimdown.rds")

# --- mod_explore_umap_celltype ---
# UMAP cell type with two sub-tabs
## sub-tab 1: exploration purpose (simpler version module-I mod_annotation umap subtab)
### Using plotly
### full dataset and subset dataset (lasso or other tools for selection)
### full dataset uses stored umap coordinates and colored by active.ident 
### subset dataset needs a re-umap active button
### add a n.pca for the re-umap using stored pca or integrated.method (e.g. integrated.rpca)
### as soon as re-umap finished, show the selected cells in new umap coordinate 
## sub_tab 2: high-res fig generation purpose
### using regular ggplot (or DimPlot() from seurat)
### having matched full or subset view sub-tab 1 
### show all cell types or highlight one or more 
### cell type(s): (unique(Idents(ojb))) 
### umap with all cell types: legends show all cell types
### highlighted cell type(s) in umap have matched color(s) as the umap with all cell types
### not selected cell type(s) show as light grey
### umap with highlighted cell types: legends only show highlighted cell tyeps
### add a figure download button and numeric input for width and length adjustment
## sub_tab 3:gene (co-)expression purpose
### using regular ggplot (FeaturePlot from seurat)
### have matched full or subset view sub-tab 1 
### add a radio button for the selection of expression or co-expression 
### selection of expression: add a textInput for gene name input (e.g. FeaturePlot(obj, features = "Gene")
### selection of expression: add a selectInput to access the non-numerical columns in metadata for "split.by"(e.g. FeaturePlot(obj, features = "Gene", split.by = "group"))
### selection of co-expression: add a selectInput to access the non-numerical columns in metadata (single choice)
### selection of co-expression: followed by the selection of the non-numerical columns in metadata, show the unique identities and select all or just one to subset the data
### selection of co-expression: add two textInput for gene names (e.g. FeaturePlot(obj, features = c("Gene1", "Gene2"), blend = T)
### add a figure download button and numeric input for width and length adjustment

# --- mod_explore_violin ---
## show gene expression by different cell types and group/condition/treatment (metadata) in violin plots
## cell type selection (one, more, or all) based on active.ident
### only show selected cell types in the plot
## textinput for one or more gene name (feature)
### one gene: stack = FALSE
### multiple gene: stack = TRUE (only)
### multiple gene: flip = TRUE or FALSE (flip the x and y axis)
## split the data (split.by = "") by information in metadata (e.g. group/condition/treatment)
### add a selectInput to access the non-numerical columns in metadata (should be categorical, single choice)
## add a figure download button and numeric input for width and length adjustment
## examples
VlnPlot(obj, features = "Cx3cr1", flip = FALSE)
VlnPlot(obj, features = "Cx3cr1", flip = FALSE, split.by = "Group")
VlnPlot(obj, features = c("Cd68", "Cx3cr1","Ccr2","Ly6c2"), stack = TRUE, flip = TRUE)
VlnPlot(obj, features = c("Cd68", "Cx3cr1","Ccr2","Ly6c2"), stack = TRUE, flip = FALSE)
VlnPlot(obj, features = c("Cd68", "Cx3cr1","Ccr2","Ly6c2"), stack = TRUE, flip = TRUE, split.by = "Group")

# --- mod_explore_dot ---
## show gene expression by different cell types and group/condition/treatment (metadata) in dot plots
## cell type selection (one, more, or all) based on active.ident
### only show selected cell types in the plot
## textinput for one or more gene name (feature)
### use DotPlot(obj, features = c("Gene1", "Gene2","Gene3","Gene4")) 
### add + RotatedAxis() to rotate the gene names on x axis
## split the data by information in metadata (e.g. group/condition/treatment)
### add a selectInput to access the non-numerical columns in metadata (should be categorical, single choice)
### DotPlot(obj, features, split.by = "") is not able to make different color scales for split identities
### manual calculation of average expression (color scale) and percent expressed (dot size) for each split identities per cell type needed
### automatic define colors for each identity as the max average expression
### light grey will represent all identities' min average expression.
## add a figure download button and numeric input for width and length adjustment

# --- mod_explore_heatmap ---
## show scaled gene expression levels by different cell types and group/condition/treatment (metadata) in heatmap
## cell type selection (one, more, or all) based on active.ident
### subset the data if not all cell types selected
## choose High variable gene (HVG) or customized gene list
### HVG: run FindVariableFeatures() on subset_obj
### HVG: NumericInput to select how many top genes selected (n) e.g. features <- VariableFeatures(obj)[1:n]
### Customized gene list: textInput for more than one genes
## scale the HVG or customized gene list e.g. obj_scale <- ScaleData(obj, features = features)
### have the metadata column names for vars.to.regress e.g. obj_scale <- ScaleData(obj, features = features, vars.to.regress = "")
## numeric input for determine how many cells (m) are used in heatmap 
### e.g. DoHeatmap(obj_scale, features = features, cells = 1:m, size = 4, angle = 90)
## add a figure download button and numeric input for width and length adjustment

# --- mod_explore_de ---
## differential analysis between ident.1 and ident.2
## ident.1 (target)
### selection of one or more cell type(s) based on active.ident
### selection of information in metadata (e.g. group/condition/treatment) - single choice
### the following selection of identity (e.g if group contains group1, group2, group3, one, more, or all identities can be selected)
### subset the data
## ident.2 (baseline)
### selection of one or more cell type(s) based on active.ident
### selection of information in metadata (e.g. group/condition/treatment) - single choice
### the following selection of identity (e.g if group contains group1, group2, group3, one, more, or all identities can be selected)
### subset the data
## filter out adjusted.p.value > 0.05
## use DTOutput to show dynamic result table
## show table by upregulated, downregulated, or both
## add action button to trigger the analysis
## add a download button to download table as csv.or xlsx. file

# --- mod_explore_enrich ---
## edit based on the old module-III code
## Enrichment ui
tabPanel("Enrichment Analysis",
         tabsetPanel(
           # subtab 1
           tabPanel(
             fluidRow(
               column(3, selectInput("database", "Database", celldbs)),
               column(2, actionButton("run.ea", "Run", class = "btn-primary" ))
             ),
             
             fluidRow(
               column(6, DT::DTOutput("dynamic.ea")),
             ),
             fluidRow(
               column(12, verbatimTextOutput('select.ea'))
             ),
             fluidRow(
               column(12, downloadButton("ea.result.table", "Download Table", class="btn-success"))
             )
           ),
           # subtab 2
           tabPanel("EA Bar Plot",
                    fluidRow(
                      column(4, radioButtons("bar", "Populations (x-axis)", choices = c("EA Top10", "Selected"), inline=T)),
                      column(3, selectInput("bar.color", "Color", choices = c("red", "green","yellow","blue", "purple"))),
                      column(2, actionButton("run.bar", "Run", class = "btn-primary" ))
                    ),
                    fluidRow(
                      column(12,  plotOutput("plot.bar", height = "600px", width = "1400px"))
                    ),
                    fluidRow(
                      ### Download figure size
                      column(3, numericInput("p.bar.w", "Width (Download Figure)", value = 8)),
                      column(3, numericInput("p.bar.h", "Height (Download Figure)", value = 6))
                    ),
                    fluidRow(
                      ### Download action button
                      column(6, downloadButton("p.bar", "Download Figure", class="btn-success"))
                    )
           )
         )
)
## Enrichment server
# result
ea.table <- eventReactive(input$run.ea,{
  result <- de.table()
  result <- enrichr(result %>% pull(gene), input$database)[[1]] %>% filter(P.value < 0.05)
  result$P.value <- as.numeric(format(result$P.value, scientific = T, digits = 3))
  result$Adjusted.P.value <- as.numeric(format(result$Adjusted.P.value, scientific = T, digits = 3))
  result$Odds.Ratio <- as.numeric(format(result$Odds.Ratio, scientific = F, digits = 3))
  result$Combined.Score <- as.numeric(format(result$Combined.Score, scientific = F, digits = 4))
  result
}
)

# table
output$dynamic.ea <- DT::renderDT(ea.table(), rownames = FALSE, options = list(pageLength = 5))

# print selected row
output$select.ea = renderPrint({
  s = input$dynamic.ea_rows_selected
  if (length(s)) {
    cat('These Pathways were selected:\n\n')
    cat(ea.table()[s, ] %>% pull(Term), sep = ', ')
  }
})

# download table
output$ea.result.table <- downloadHandler(
  filename = function() {
    paste0("EA_", input$de.cluster, "_", input$de.condition.1, "_vs_", input$de.condition.2, "_", input$database ,".xlsx")
  },
  content = function(file) {
    write.xlsx(ea.table(), file)
  }
)
# figure
plotInput.bar <- eventReactive(input$run.bar,{
  df.bar<- ea.table()
  df.bar$`-Log10(Adjusted.P.value)` <- -log10(df.bar$Adjusted.P.value)
  if(input$bar == "EA Top10"){
    df.bar <- df.bar[1:10,]
  }else{
    {
      s = input$dynamic.ea_rows_selected
      if (length(s)){
        df.bar <- df.bar[s,]
      } 
    }
  }
  p <- ggplot(df.bar, aes(x= `-Log10(Adjusted.P.value)`, y = reorder(Term,`-Log10(Adjusted.P.value)`)))+
    geom_bar(stat = "identity", aes(fill = Odds.Ratio))+
    scale_fill_gradient2(high = input$bar.color)+
    theme_classic()+
    ylab("")
})
output$plot.bar <- renderPlot({
  print(plotInput.bar())
}, res = 96)

# download right figure
output$p.bar <- downloadHandler(
  filename = function(){paste0("EA_", input$de.cluster, "_", input$de.condition.1, "_vs_", input$de.condition.2, "_", input$database ,".jpg")},
  content = function(file){
    ggsave(file, plot = plotInput.bar(), width = input$p.bar.w, height =  input$p.bar.h)
  }
)

