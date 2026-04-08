mod_integrate_ui <- function(id) {
  ns <- NS(id)
  tabPanel("Integation",
           sidebarLayout(
             sidebarPanel(
               h4("1. Normalization"),
               radioButtons(ns("norm_method"), "Method:",
                            choices = c("Log-Normalization" = "LogNormalize", "SCTransform" = "SCT")),
               
               checkboxGroupInput(ns("vars_regress"), "Regress Factors:",
                                  choices = c("nCount_RNA", "percent.mt", "percent.ribo"),
                                  selected = "nCount_RNA"),
               hr(),
               h4("2. Integration"),
               checkboxGroupInput(ns("int_methods"), "Methods to Run:",
                                  choices = c("CCA", "RPCA", "Harmony", "FastMNN"),
                                  selected = c("CCA", "Harmony")),
               
               actionButton(ns("run_flow"), "Run Pipeline", class = "btn-success btn-block")
             ),
             
             mainPanel(
               h4("Pipeline Status Log"),
               verbatimTextOutput(ns("status_log")),
               uiOutput(ns("finished_ui")) 
             )
           )
  )
}

mod_integrate_server <- function(id, shared_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # --- UI LOGIC: Disable FastMNN if SCT is selected ---
    observeEvent(input$norm_method, {
      if (input$norm_method == "SCT") {
        updateCheckboxGroupInput(session, "int_methods",
                                 choices = c("CCA", "RPCA", "Harmony"), # FastMNN removed
                                 selected = setdiff(input$int_methods, "FastMNN")
        )
      } else {
        updateCheckboxGroupInput(session, "int_methods",
                                 choices = c("CCA", "RPCA", "Harmony", "FastMNN"),
                                 selected = input$int_methods
        )
      }
    })
    
    processed_obj <- eventReactive(input$run_flow, {
      req(shared_data())
      options(future.globals.maxSize = 10 * 1024^3)
      
      # Start with clean object
      raw_obj <- shared_data()
      
      withProgress(message = 'Running Pipeline...', value = 0, {
        
        # --- INITIALIZATION: SplitObject & Re-merge to create Layers ---
        incProgress(0.1, detail = "Initializing Batch Layers...")
        obj.list <- SplitObject(raw_obj, split.by = "batch")
        obj <- merge(obj.list[[1]], y = obj.list[2:length(obj.list)])
        
        # --- STEP 1: PREPROCESSING ---
        if (input$norm_method == "LogNormalize") {
          target_assay <- "RNA"
          incProgress(0.1, detail = "Log-Normalizing...")
          
          obj <- NormalizeData(obj)
          obj <- FindVariableFeatures(obj, nfeatures = 2000)
          
          # Sync features across layers manually to be safe for FastMNN/RPCA
          #all_layers <- Layers(obj[[target_assay]], search = "data")
          #common_genes <- Reduce(intersect, lapply(all_layers, function(l) rownames(GetAssayData(obj, layer = l))))
          #features_to_integrate <- intersect(VariableFeatures(obj), common_genes)
          #VariableFeatures(obj) <- features_to_integrate
          
          obj <- ScaleData(obj, vars.to.regress = input$vars_regress)
          
          obj <- RunPCA(obj, verbose = FALSE)
          obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, graph.name = "snn")
          obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, reduction.name = "umap", verbose = FALSE)
          
        } else {
          target_assay <- "SCT"
          incProgress(0.2, detail = "Running SCTransform...")
          
          # SCTransform handles the layered object automatically in v5
          obj <- SCTransform(obj, vars.to.regress = input$vars_regress, verbose = FALSE)
          
          # Sync Variable Features for SCT layers
          obj <- FindVariableFeatures(obj, assay = "SCT", nfeatures = 2000)
          #features_to_integrate <- VariableFeatures(obj, assay = "SCT")
          
          incProgress(0.1, detail = "Running PCA and UMAP on SCT...")
          obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
          obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, graph.name = "snn")
          obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, reduction.name = "umap", verbose = FALSE)
          
        }
        
        # --- STEP 2: INTEGRATION, CLUSTERING & UMAP ---
        methods <- input$int_methods
        method_map <- list(
          "CCA"     = Seurat::CCAIntegration, 
          "RPCA"    = Seurat::RPCAIntegration, 
          "Harmony" = Seurat::HarmonyIntegration, 
          "FastMNN" = SeuratWrappers::FastMNNIntegration
        )
        
        for (m in methods) {
          incProgress(0.1, detail = paste("Processing", m))
          red_name <- paste0("integrated.", tolower(m))
          umap_name <- paste0("umap.", tolower(m))
          
          obj <- tryCatch({
            # --- YOUR SPECIFIC INTEGRATION LOGIC ---
            if (input$norm_method == "SCT") {
              obj <- IntegrateLayers(
                object = obj, 
                method = method_map[[m]],
                orig.reduction = "pca", 
                new.reduction = red_name,
                assay = target_assay,
                normalization.method = "SCT",
                verbose = FALSE
              )
            } else {
              # LogNormalize path: normalization.method is omitted as per your discovery
              obj <- IntegrateLayers(
                object = obj, 
                method = method_map[[m]],
                orig.reduction = "pca", 
                new.reduction = red_name,
                assay = target_assay,
                verbose = FALSE
              )
            }
            
            # --- ADDING NEIGHBORS & UMAP ---
            # We run these on the newly created integrated reduction
            obj <- FindNeighbors(
              obj, 
              reduction = red_name, 
              dims = 1:30, 
              graph.name = paste0("snn.", tolower(m))
            )
            
            obj <- RunUMAP(
              obj, 
              reduction = red_name, 
              dims = 1:30, 
              reduction.name = umap_name, 
              verbose = FALSE
            )
            
            obj
          }, error = function(e) {
            showNotification(paste(m, "Integration/UMAP failed:", e$message), type = "error")
            return(obj)
          })
        }
        
        
        
        # --- STEP 3: FINALIZATION ---
        
        incProgress(0.1, detail = "Finalizing Assay...")
        
        # 1. Identify the class of the target assay
        assay_class <- class(obj[[target_assay]])
        
        # 2. Logic-based Joining
        if (any(grepl("SCTAssay", assay_class))) {
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
        DefaultAssay(obj) <- target_assay
        
        return(obj)
      })
    })
    
    output$status_log <- renderPrint({
      req(processed_obj())
      obj <- processed_obj()
      cat("--- Pipeline Complete ---\n")
      cat("Normalization: ", input$norm_method, "\n")
      cat("Reductions:    ", paste(names(obj@reductions), collapse = ", "), "\n")
    })
    
    return(processed_obj)
  })
}