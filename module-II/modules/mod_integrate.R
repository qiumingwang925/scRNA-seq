mod.integrate.ui <- function(id) {
  ns <- NS(id)
  tabPanel("Integation",
           sidebarLayout(
             sidebarPanel(
               h4("1. Normalization"),
               radioButtons(ns("norm.method"), "Method:",
                            choices = c("LogNorm" = "LogNormalize", "SCTransform" = "SCT")),

               selectInput(ns("vars.regress"), "Regress Factors:",
                           choices = c("nCount_RNA", "percent.mt", "percent.rp", "S.Score", "G2M.Score"),
                           selected = c("nCount_RNA", "percent.mt"),
                           multiple = TRUE),
               hr(),
               h4("2. Integration"),
               checkboxGroupInput(ns("int.methods"), "Methods to Run:",
                                  choices = c("CCA", "RPCA", "Harmony", "FastMNN"),
                                  selected = c("CCA", "Harmony")),
               
               actionButton(ns("run.flow"), "Run Pipeline", class = "btn-success btn-block")
             ),
             
             mainPanel(
               h4("Pipeline Status Log"),
               verbatimTextOutput(ns("status.log")),
               uiOutput(ns("finished.ui")) 
             )
           )
  )
}

mod.integrate.server <- function(id, shared.data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # --- UI LOGIC: Disable FastMNN if SCT is selected ---
    observeEvent(input$norm.method, {
      if (input$norm.method == "SCT") {
        # SCTransform already models sequencing depth, so nCount_RNA/percent.rp
        # are not offered as regressors here; FastMNN isn't supported on SCT.
        updateCheckboxGroupInput(session, "int.methods",
                                 choices = c("CCA", "RPCA", "Harmony"), # FastMNN removed
                                 selected = setdiff(input$int.methods, "FastMNN")
        )
        regress.choices <- c("percent.mt", "S.Score", "G2M.Score")
      } else {
        updateCheckboxGroupInput(session, "int.methods",
                                 choices = c("CCA", "RPCA", "Harmony", "FastMNN"),
                                 selected = input$int.methods
        )
        regress.choices <- c("nCount_RNA", "percent.mt", "percent.rp", "S.Score", "G2M.Score")
      }
      updateSelectInput(session, "vars.regress",
                        choices = regress.choices,
                        selected = intersect(input$vars.regress, regress.choices))
    })
    
    processed.obj <- eventReactive(input$run.flow, {
      req(shared.data())
      options(future.globals.maxSize = 10 * 1024^3)
      
      # Start with clean object
      raw.obj <- shared.data()
      
      withProgress(message = 'Running Pipeline...', value = 0, {
        
        # --- INITIALIZATION: SplitObject & Re-merge to create Layers ---
        incProgress(0.1, detail = "Initializing Batch Layers...")
        obj.list <- SplitObject(raw.obj, split.by = "batch")
        obj <- merge(obj.list[[1]], y = obj.list[2:length(obj.list)])
        
        # --- STEP 1: PREPROCESSING ---
        if (input$norm.method == "LogNormalize") {
          target.assay <- "RNA"
          incProgress(0.1, detail = "Log-Normalizing...")
          
          obj <- NormalizeData(obj)
          obj <- FindVariableFeatures(obj, nfeatures = 2000)
          
          # Sync features across layers manually to be safe for FastMNN/RPCA
          #all_layers <- Layers(obj[[target.assay]], search = "data")
          #common_genes <- Reduce(intersect, lapply(all_layers, function(l) rownames(GetAssayData(obj, layer = l))))
          #features_to_integrate <- intersect(VariableFeatures(obj), common_genes)
          #VariableFeatures(obj) <- features_to_integrate
          
          obj <- ScaleData(obj, vars.to.regress = input$vars.regress)
          
          obj <- RunPCA(obj, verbose = FALSE)
          obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, graph.name = "snn")
          obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, reduction.name = "umap", verbose = FALSE)
          
        } else {
          target.assay <- "SCT"
          incProgress(0.2, detail = "Running SCTransform...")
          
          # SCTransform handles the layered object automatically in v5. It picks
          # its own variable features (return.only.var.genes = TRUE), so no
          # separate FindVariableFeatures step is needed for the SCT path.
          obj <- SCTransform(obj, vars.to.regress = input$vars.regress,
                             variable.features.n = 3000,
                             return.only.var.genes = TRUE, verbose = FALSE)

          incProgress(0.1, detail = "Running PCA and UMAP on SCT...")
          obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
          obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, graph.name = "snn")
          obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, reduction.name = "umap", verbose = FALSE)
          
        }
        
        # --- STEP 2: INTEGRATION, CLUSTERING & UMAP ---
        methods <- input$int.methods
        method.map <- list(
          "CCA"     = Seurat::CCAIntegration, 
          "RPCA"    = Seurat::RPCAIntegration, 
          "Harmony" = Seurat::HarmonyIntegration, 
          "FastMNN" = SeuratWrappers::FastMNNIntegration
        )
        
        for (m in methods) {
          incProgress(0.1, detail = paste("Processing", m))
          red.name <- paste0("integrated.", tolower(m))
          umap.name <- paste0("umap.", tolower(m))
          
          obj <- tryCatch({
            if (input$norm.method == "SCT") {
              obj <- IntegrateLayers(
                object = obj,
                method = method.map[[m]],
                orig.reduction = "pca",
                new.reduction = red.name,
                assay = target.assay,
                normalization.method = "SCT",
                verbose = FALSE
              )
            } else {
              obj <- IntegrateLayers(
                object = obj,
                method = method.map[[m]],
                orig.reduction = "pca",
                new.reduction = red.name,
                assay = target.assay,
                verbose = FALSE
              )
            }
            obj
          }, error = function(e) {
            msg <- paste(m, "IntegrateLayers failed:", conditionMessage(e))
            message(msg, "\n", paste(capture.output(traceback()), collapse = "\n"))
            showNotification(msg, type = "error")
            return(obj)
          })

          # Only run downstream steps if integration produced the expected reduction
          if (red.name %in% names(obj@reductions)) {
            obj <- tryCatch({
              obj <- FindNeighbors(
                obj,
                reduction = red.name,
                dims = 1:30,
                graph.name = paste0("snn.", tolower(m))
              )
              obj <- RunUMAP(
                obj,
                reduction = red.name,
                dims = 1:30,
                reduction.name = umap.name,
                verbose = FALSE
              )
              obj
            }, error = function(e) {
              msg <- paste(m, "Neighbors/UMAP failed:", conditionMessage(e))
              message(msg, "\n", paste(capture.output(traceback()), collapse = "\n"))
              showNotification(msg, type = "error")
              return(obj)
            })
          }
        }
        
        
        
        # --- STEP 3: FINALIZATION ---
        
        incProgress(0.1, detail = "Finalizing Assay...")
        
        # 1. Identify the class of the target assay
        assay.class <- class(obj[[target.assay]])
        
        # 2. Logic-based Joining
        if (any(grepl("SCTAssay", assay.class))) {
          # IF SCT: We skip JoinLayers because SCTAssay doesn't support the method.
          # Most downstream v5 functions handle SCT layers automatically.
          message("SCTAssay detected. Skipping JoinLayers to prevent UseMethod error.")
          
        } else {
          # IF RNA (Assay5): JoinLayers is supported and recommended.
          obj <- tryCatch({
            JoinLayers(obj, assay = "RNA")
          }, error = function(e) {
            message("RNA JoinLayers skipped: ", e$message)
            return(obj)
          })
        }
        
        # 3. Set Default Assay
        DefaultAssay(obj) <- target.assay
        
        return(obj)
      })
    })
    
    output$status.log <- renderPrint({
      req(processed.obj())
      obj <- processed.obj()
      cat("--- Pipeline Complete ---\n")
      cat("Normalization: ", input$norm.method, "\n")
      cat("Reductions:    ", paste(names(obj@reductions), collapse = ", "), "\n")
    })
    
    return(processed.obj)
  })
}