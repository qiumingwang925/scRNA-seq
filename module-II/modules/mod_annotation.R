mod_annotation_ui <- function(id) {
  ns <- NS(id)
  tabPanel("Visual Inspection & Annotation",
           sidebarLayout(
             sidebarPanel(
               width = 4,
               
               # Integration selection
               selectInput(ns("view_red"), "Select Integration:", choices = NULL),
               # Display Mode Toggle: Controls what the UMAP shows 
               radioButtons(ns("display_mode"), "Display UMAP from:",
                            choices = c("All cells" = "full", "Selected cells (subset)" = "subset"),
                            selected = "full"),
               numericInput(ns("subset_pcs"), "PCs for Subset UMAP:",
                            value = 10, min = 2, max = 50, step = 1),
               actionButton(ns("run_subset_umap"), "📊 Run UMAP on Selection", class = "btn-warning", style="width:100%")
               
               hr(),
               tabsetPanel(
                 id = ns("viz_mode"),
                 tabPanel("Metadata",
                          br(),
                          # 2. FIXED UI ELEMENTS
                          selectInput(ns("view_meta"), "Color by Metadata:", choices = NULL),
                          #helpText("View categories like Batch or Cell Type."),
                          hr(),
                          numericInput(ns("cluster_res"), "Cluster Resolution:", value = 0.1, step = 0.05),
                          actionButton(ns("run_clustering"), "🔍 Run Clustering", class = "btn-info", style="width:100%"),
                 ),
                 tabPanel("Expression",
                          br(),
                          fileInput(ns("marker_file"), "Upload Biomarker CSV", accept = ".csv"),
                          
                          # Cascading Selectors
                          selectInput(ns("gene_major"), "Select Major Cell Type:", choices = NULL),
                          selectInput(ns("gene_fine"), "Select Fine Cell Type:", choices = NULL),
                          selectInput(ns("gene_selected"), "Select Gene/Marker:", choices = NULL),
                          
                          hr(),
                          # Manual Overwrite
                          textInput(ns("gene_entered"), "OR Type Gene Name:", value = "", placeholder = "e.g., Cd68"),
                          helpText("Typing a gene here will override the dropdown selection.")
                 ),
                 
               ),
               
               hr(),
               # --- NEW DE TOOL BLOCK ---
               h4("Differential Expression"),
               p(tags$small("Identify markers for the current selection vs. the rest.")),
               actionButton(ns("run_de"), "🧬 Run DE Analysis", 
                            class = "btn-success", style="width: 100%"),
               hr(),
               h4("Manual Annotation"),
               p(tags$small("Use Lasso/Box tool on UMAP to select cells.")),
               textInput(ns("new_label_name"), "Enter New Label:", placeholder = "e.g., Alveolar Macrophage"),
               actionButton(ns("update_label"), "Apply Label to Selection", 
                            class = "btn-primary", style="width: 100%"),
               br(), br(),
               mod.save.config.ui(ns("save"), label = "Export Annotated Object")
             ),
             
             mainPanel(
               width = 8,
               tabsetPanel(
                 id = ns("main_tabs"),
                 tabPanel("UMAP Visualization", 
                          br(),
                          plotly::plotlyOutput(ns("umap_main"), height = "600px")
                 ),
                 tabPanel("DE Analysis",
                          br(),
                          DT::DTOutput(ns("marker_table"))
                 )
               )
             )
           )
  )
}

mod_annotation_server <- function(id, shared_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    current_obj <- reactiveVal(NULL)
    subset_obj  <- reactiveVal(NULL) # Holds the newly calculated subset
    
    # Load initial data
    observe({
      req(shared_data())
      if (is.null(current_obj())) current_obj(shared_data())
    })
    
    # Update UI choices based on active object
    observeEvent(list(current_obj(), input$display_mode), {
      # Determine which object we are looking at
      obj <- if(input$display_mode == "subset") subset_obj() else current_obj()
      req(obj)
      
      # Update Metadata dropdown
      updateSelectInput(session, "view_meta", choices = colnames(obj@meta.data), selected = input$view_meta)
      
      # Update Reductions 
      reds <- names(obj@reductions)
      umap_reds <- reds[grepl("umap", reds, ignore.case = TRUE)]
      if (length(umap_reds) > 0) {
        choices_vec <- setNames(umap_reds, sapply(umap_reds, function(x) {
          if (tolower(x) == "umap") {
            return("Unintegrated")
          } else {
            # Remove "umap." prefix and capitalize the rest
            clean_name <- gsub("umap\\.", "", x, ignore.case = TRUE)
            return(toupper(clean_name))
          }
        }))
      updateSelectInput(session, "view_red", choices = choices_vec, selected = input$view_red)
      }
    })
    
    # --- Gene Expression ---
    ## upload biomaker_file
    marker_df <- reactive({
      req(input$marker_file)
      read.csv(input$marker_file$datapath)
    })
    
    ## When CSV is loaded, update Major Cell Types
    observeEvent(marker_df(), {
      df <- marker_df()
      updateSelectInput(session, "gene_major", choices = c("", unique(df$Label.main)))
    })
    
    ## When Major is selected, update Fine Cell Types
    observeEvent(input$gene_major, {
      req(input$gene_major != "", marker_df())
      df <- marker_df()
      sub_df <- df[df$Label.main == input$gene_major, ]
      updateSelectInput(session, "gene_fine", choices = c("", unique(sub_df$Label.fine)))
    })
    
    ## When Fine is selected, update Gene list
    observeEvent(input$gene_fine, {
      req(input$gene_fine != "", marker_df())
      df <- marker_df()
      sub_df <- df[df$Label.main == input$gene_major & df$Label.fine == input$gene_fine, ]
      updateSelectInput(session, "gene_selected", choices = c("", unique(sub_df$Markers)))
    })
    
    # --- Subset data logic---
    observeEvent(input$run_subset_umap, {
      # 1. Determine the source object based on current view
      source_obj <- if(input$display_mode == "subset") subset_obj() else current_obj()
      req(source_obj)
      
      # 2. Get the selection from Plotly
      sel <- plotly::event_data("plotly_selected", source = "umap_select")
      validate(need(!is.null(sel), "Please lasso select cells on the UMAP first."))
      
      # 3. Extract barcodes
      selected_barcodes <- unique(as.character(sel$customdata))
      req(length(selected_barcodes) > 3) # Need at least a few cells for PCA
      
      withProgress(message = 'Calculating Subset UMAP...', value = 0, {
        tryCatch({
          # 4. Subset from the SOURCE (allows subsetting a subset)
          sub <- subset(source_obj, cells = selected_barcodes)
          
          # 5. Determine reduction for processing
          # If we are already in a subset, we use the subset's PCA/UMAP
          if (input$view_red == "umap") {
            red_base <- "pca"
          } else {
            # This removes "umap." (case insensitive) to get "harmony", "rpca", etc.
            red_base <- gsub("umap\\.", "integrated.", input$view_red)
          }
          # safe check
          if (!red_base %in% Reductions(sub)) red_base <- "pca" 
          
          # Cap PCs to what the data supports: must be < min(cells, features)
          max.pcs <- min(ncol(sub), nrow(sub)) - 1
          n.pcs <- min(input$subset_pcs, max.pcs)
          sub <- RunPCA(sub, npcs = n.pcs)
          sub <- FindNeighbors(sub, dims = 1:n.pcs)
          sub <- RunUMAP(sub, dims = 1:n.pcs)
          # 7. Update the subset reactive
          subset_obj(sub)
          # 8. Force view to subset mode
          updateRadioButtons(session, "display_mode", selected = "subset")
          showNotification("New subset calculated!", type = "message")
          
        }, error = function(e) {
          # This catches the error and shows it as a red notification instead of crashing
          showNotification(paste("Error in UMAP calculation:", e$message), 
                           type = "error", duration = 10)
          print(paste("Detailed Error:", e)) # Logs to the R console for you
        })
      })
    })
    
    # --- CLUSTERING LOGIC --- 
    ##WITH METADATA CLEANING
    observeEvent(input$run_clustering, { 
      req(input$display_mode, input$view_red)
      
      # Normalize input once
      selected_view <- tolower(input$view_red)
      res <- input$cluster_res 
      
      tryCatch({
        withProgress(message = 'Running Clustering...', value = 0, {
          # Identify active object
          is_full <- input$display_mode == "full"
          obj <- if(is_full) current_obj() else subset_obj()
          req(obj)
          
          # Determine Graph Name
          target_graph <- if (selected_view == "umap") "snn" else gsub("umap\\.", "snn.", selected_view)
          
          # Validate graph exists
          if (!target_graph %in% names(obj@graphs)) {
            stop(paste("Graph", target_graph, "not found. Please run Neighbors first."))
          }
          
          # Run Seurat Clustering
          obj <- FindClusters(obj, resolution = res, graph.name = target_graph) 
          
          # Cleanup metadata (keeps only the final 'seurat_clusters')
          # This regex handles snn_res and snn.harmony_res
          pat <- paste0(target_graph, "_res\\.")
          obj@meta.data <- obj@meta.data[, !grepl(pat, colnames(obj@meta.data)), drop = FALSE] 
          
          # Update state
          if (is_full) current_obj(obj) else subset_obj(obj)
          
          showNotification(paste("Success: Clustered using", target_graph), type = "message")
        })
      }, error = function(e) {
        showNotification(paste("Clustering Error:", e$message), type = "error", duration = 10)
      })
    })
    
    # Render plot
    output$umap_main <- plotly::renderPlotly({
      # 1. Determine which object to use based on Display Mode
      if (input$display_mode == "subset") {
        # Check if subset exists; if not, show a message instead of an error
        validate(need(!is.null(subset_obj()), "No subset calculated yet. Please select cells and click 'Run UMAP'."))
        obj <- subset_obj()
      } else {
        req(current_obj())
        obj <- current_obj()
      }
      
      # 2. Basic requirements for plotting
      req(input$view_red, input$view_meta)
      # Ensure the selected reduction actually exists in the active object
      # (Important because subset objects only have the new 'umap' reduction)
      validate(need(input$view_red %in% names(obj@reductions), 
                    "The selected reduction is not available for this view."))
      
      # 3. Expression or Metadata
      if (input$viz_mode == "Expression") {
        # Determine which gene to use (Manual Entry takes priority)
        gene_to_plot <- if (nzchar(input$gene_entered)) {
          input$gene_entered 
        } else {
          input$gene_selected
        }
        req(gene_to_plot) 
        # CHECK IF GENE EXISTS: Provides a helpful message if there's a typo
        validate(
          need(gene_to_plot %in% rownames(obj), 
               paste0("Gene '", gene_to_plot, "' not found in this dataset. Check spelling/case."))
        )
        p <- FeaturePlot(obj, features = gene_to_plot, reduction = input$view_red)
      } else {
        req(input$view_meta)
        # Verify metadata column exists in the active object
        validate(need(input$view_meta %in% colnames(obj@meta.data), "Metadata column not found."))
        meta_data <- obj[[input$view_meta]][, 1]
        if(is.numeric(meta_data)){
          p <- FeaturePlot(obj, reduction = input$view_red, features = input$view_meta)
        }else{
          p <- DimPlot(obj, reduction = input$view_red, group.by = input$view_meta)
        }
      }
      # Carry barcodes for selection
      p$layers[[1]]$mapping$customdata <- aes(label = rownames(p$data))$label
      # plot
      plotly::ggplotly(p + theme_minimal(), source = "umap_select") %>% 
        layout(dragmode = "lasso")
    })
    
    
    
    # Updated Manual Annotation with Sync Support
    observeEvent(input$update_label, {
      select_data <- plotly::event_data("plotly_selected", source = "umap_select")
      req(select_data, input$new_label_name != "")
      
      selected_barcodes <- unique(as.character(select_data$customdata))
      new_label <- input$new_label_name
      
      # --- UPDATE MASTER OBJECT ---
      master <- isolate(current_obj())
      
      # Initialize column if missing
      if (!"manual_annotation" %in% colnames(master@meta.data)) {
        master$manual_annotation <- "Unlabeled"
      }
      
      master@meta.data[selected_barcodes, "manual_annotation"] <- new_label
      current_obj(master) # Push update to master
      
      # --- UPDATE SUBSET OBJECT (If it exists) ---
      if (!is.null(subset_obj())) {
        sub <- isolate(subset_obj())
        
        # Ensure column exists in subset too
        if (!"manual_annotation" %in% colnames(sub@meta.data)) {
          sub$manual_annotation <- "Unlabeled"
        }
        
        # Only update barcodes that are actually in this subset
        cells_to_update_in_sub <- intersect(selected_barcodes, colnames(sub))
        
        if (length(cells_to_update_in_sub) > 0) {
          sub@meta.data[cells_to_update_in_sub, "manual_annotation"] <- new_label
          subset_obj(sub) # Push update to subset so the plot refreshes
        }
      }
      
      showNotification(paste("Applied label:", new_label), type = "message")
    })
    
    # --- 3. DE Analysis ---
    marker_data <- eventReactive(input$run_de, {
      
      # 1. Capture Selection from Plotly
      # Note: We use the source "umap_select" defined in your ggplotly call
      selected_de <- plotly::event_data("plotly_selected", source = "umap_select")
      
      # Validation: Ensure something was selected
      validate(
        need(!is.null(selected_de), "Please use the Lasso or Box tool to select cells on the UMAP first.")
      )
      
      # Extract barcodes from the selection
      selected_cells <- unique(as.character(selected_de$customdata))
      
      # 2. Determine the Universe (Match your existing logic)
      # In your code, you use current_obj() and subset_obj()
      universe_data <- if (input$display_mode == "full") {
        current_obj() 
      } else {
        subset_obj()
      }
      
      req(universe_data) # Ensure object exists
      
      # 3. Define Groups
      all_cells_in_universe <- colnames(universe_data)
      
      # Ensure selected cells actually exist in the current universe 
      # (Prevents errors if switching modes with an old selection)
      selected_cells <- intersect(selected_cells, all_cells_in_universe)
      
      rest_of_cells <- setdiff(all_cells_in_universe, selected_cells)
      
      # 4. Validation: Check group sizes
      validate(
        need(length(selected_cells) > 3, "Please select at least 4 cells for DE."),
        need(length(rest_of_cells) > 3, "Not enough background cells for comparison.")
      )
      
      # 5. Run DE Analysis with Progress Bar
      withProgress(message = 'Calculating DE Markers...', value = 0, {
        
        incProgress(0.2, detail = "Preparing groups...")
        showNotification("Starting Differential Expression...", type = "message")
        
        # Run FindMarkers
        de_results <- tryCatch({
          incProgress(0.5, detail = "Running FindMarkers (Wilcoxon)...")
          Seurat::FindMarkers(
            object = universe_data,
            ident.1 = selected_cells,
            ident.2 = rest_of_cells,
            only.pos = TRUE,
            min.pct = 0.1,
            logfc.threshold = 0.25
          )
        }, error = function(e) {
          showNotification(paste("Error in DE:", e$message), type = "error")
          return(NULL)
        })
        
        incProgress(0.2, detail = "Finalizing table...")
        
        if (is.null(de_results) || nrow(de_results) == 0) {
          showNotification("No significant markers found.", type = "warning")
          return(NULL)
        }
        
        showNotification("DE Analysis Complete!", type = "default")
        
        # Return results with Gene column
        return(cbind(Gene = rownames(de_results), de_results))
      })
    })
    
    # Render the DT Table
    output$marker_table <- DT::renderDT({
      # Use req() so the table is empty until the button is pressed
      req(marker_data())
      
      DT::datatable(
        marker_data(),
        rownames = FALSE,
        options = list(pageLength = 10, scrollX = TRUE),
        selection = "single"
      )
    })
    
    # Switch to DE analysis sub-tab
    observeEvent(input$run_de, {
      updateTabsetPanel(session, "main_tabs", selected = "DE Analysis")
    })
    
    # --- EXPORT ---
    mod.save.config.server("save", current_obj)
    
    return(current_obj)
  })
}