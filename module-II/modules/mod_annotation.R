mod.annotation.ui <- function(id) {
  ns <- NS(id)
  tabPanel("Visual Inspection & Annotation",
           sidebarLayout(
             sidebarPanel(
               width = 4,
               
               # Integration selection
               selectInput(ns("view.red"), "Select Integration:", choices = NULL),
               # Display Mode Toggle: Controls what the UMAP shows 
               radioButtons(ns("display.mode"), "Display UMAP from:",
                            choices = c("All cells" = "full", "Selected cells (subset)" = "subset"),
                            selected = "full"),
               numericInput(ns("subset.pcs"), "PCs for Subset UMAP:",
                            value = 10, min = 2, max = 50, step = 1),
               actionButton(ns("run.subset.umap"), "📊 Run UMAP on Selection", class = "btn-warning", style="width:100%"),
               
               hr(),
               tabsetPanel(
                 id = ns("viz.mode"),
                 tabPanel("Metadata",
                          br(),
                          # 2. FIXED UI ELEMENTS
                          selectInput(ns("view.meta"), "Color by Metadata:", choices = NULL),
                          #helpText("View categories like Batch or Cell Type."),
                          hr(),
                          numericInput(ns("cluster.res"), "Cluster Resolution:", value = 0.1, step = 0.05),
                          actionButton(ns("run.clustering"), "🔍 Run Clustering", class = "btn-info", style="width:100%"),
                 ),
                 tabPanel("Expression",
                          br(),
                          fileInput(ns("marker.file"), "Upload Biomarker CSV", accept = ".csv"),
                          
                          # Cascading Selectors
                          selectInput(ns("gene.major"), "Select Major Cell Type:", choices = NULL),
                          selectInput(ns("gene.fine"), "Select Fine Cell Type:", choices = NULL),
                          selectInput(ns("gene.selected"), "Select Gene/Marker:", choices = NULL),
                          
                          hr(),
                          # Manual Overwrite
                          textInput(ns("gene.entered"), "OR Type Gene Name:", value = "", placeholder = "e.g., Cd68"),
                          helpText("Typing a gene here will override the dropdown selection.")
                 ),
                 
               ),
               
               hr(),
               # --- NEW DE TOOL BLOCK ---
               h4("Differential Expression"),
               p(tags$small("Identify markers for the current selection vs. the rest.")),
               actionButton(ns("run.de"), "🧬 Run DE Analysis", 
                            class = "btn-success", style="width: 100%"),
               hr(),
               h4("Manual Annotation"),
               p(tags$small("Use Lasso/Box tool on UMAP to select cells.")),
               textInput(ns("new.label.name"), "Enter New Label:", placeholder = "e.g., Alveolar Macrophage"),
               actionButton(ns("update.label"), "Apply Label to Selection", 
                            class = "btn-primary", style="width: 100%"),
               br(), br(),
               mod.save.config.ui(ns("save"), label = "Export Annotated Object")
             ),
             
             mainPanel(
               width = 8,
               tabsetPanel(
                 id = ns("main.tabs"),
                 tabPanel("UMAP Visualization", 
                          br(),
                          plotly::plotlyOutput(ns("umap.main"), height = "600px")
                 ),
                 tabPanel("DE Analysis",
                          br(),
                          DT::DTOutput(ns("marker.table"))
                 )
               )
             )
           )
  )
}

mod.annotation.server <- function(id, shared.data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    current.obj <- reactiveVal(NULL)
    subset.obj  <- reactiveVal(NULL) # Holds the newly calculated subset
    
    # Load initial data
    observe({
      req(shared.data())
      if (is.null(current.obj())) current.obj(shared.data())
    })
    
    # Update UI choices based on active object
    observeEvent(list(current.obj(), input$display.mode), {
      # Determine which object we are looking at
      obj <- if(input$display.mode == "subset") subset.obj() else current.obj()
      req(obj)
      
      # Update Metadata dropdown
      updateSelectInput(session, "view.meta", choices = colnames(obj@meta.data), selected = input$view.meta)
      
      # Update Reductions 
      reds <- names(obj@reductions)
      umap.reds <- reds[grepl("umap", reds, ignore.case = TRUE)]
      if (length(umap.reds) > 0) {
        choices.vec <- setNames(umap.reds, sapply(umap.reds, function(x) {
          if (tolower(x) == "umap") {
            return("Unintegrated")
          } else {
            # Remove "umap." prefix and capitalize the rest
            clean.name <- gsub("umap\\.", "", x, ignore.case = TRUE)
            return(toupper(clean.name))
          }
        }))
      updateSelectInput(session, "view.red", choices = choices.vec, selected = input$view.red)
      }
    })
    
    # --- Gene Expression ---
    ## upload biomaker_file
    marker.df <- reactive({
      req(input$marker.file)
      read.csv(input$marker.file$datapath)
    })
    
    ## When CSV is loaded, update Major Cell Types
    observeEvent(marker.df(), {
      df <- marker.df()
      updateSelectInput(session, "gene.major", choices = c("", unique(df$Label.main)))
    })
    
    ## When Major is selected, update Fine Cell Types
    observeEvent(input$gene.major, {
      req(input$gene.major != "", marker.df())
      df <- marker.df()
      sub.df <- df[df$Label.main == input$gene.major, ]
      updateSelectInput(session, "gene.fine", choices = c("", unique(sub.df$Label.fine)))
    })
    
    ## When Fine is selected, update Gene list
    observeEvent(input$gene.fine, {
      req(input$gene.fine != "", marker.df())
      df <- marker.df()
      sub.df <- df[df$Label.main == input$gene.major & df$Label.fine == input$gene.fine, ]
      updateSelectInput(session, "gene.selected", choices = c("", unique(sub.df$Markers)))
    })
    
    # --- Subset data logic---
    observeEvent(input$run.subset.umap, {
      # 1. Determine the source object based on current view
      source.obj <- if(input$display.mode == "subset") subset.obj() else current.obj()
      req(source.obj)
      
      # 2. Get the selection from Plotly
      sel <- plotly::event_data("plotly_selected", source = "umap_select")
      validate(need(!is.null(sel), "Please lasso select cells on the UMAP first."))
      
      # 3. Extract barcodes
      selected.barcodes <- unique(as.character(sel$customdata))
      req(length(selected.barcodes) > 3) # Need at least a few cells for PCA
      
      withProgress(message = 'Calculating Subset UMAP...', value = 0, {
        tryCatch({
          # 4. Subset from the SOURCE (allows subsetting a subset)
          sub <- subset(source.obj, cells = selected.barcodes)
          
          # 5. Determine reduction for processing
          # If we are already in a subset, we use the subset's PCA/UMAP
          if (input$view.red == "umap") {
            red.base <- "pca"
          } else {
            # This removes "umap." (case insensitive) to get "harmony", "rpca", etc.
            red.base <- gsub("umap\\.", "integrated.", input$view.red)
          }
          # safe check
          if (!red.base %in% Reductions(sub)) red.base <- "pca" 
          
          # Cap to available dimensions in the reduction
          max.pcs <- ncol(Embeddings(sub, red.base))
          n.pcs <- min(input$subset.pcs, max.pcs)
          print(paste("Number of PCs to use:", n.pcs))
          # sub <- RunPCA(sub, npcs = n.pcs)
          # sub <- FindNeighbors(sub, dims = 1:n.pcs)
          sub <- RunUMAP(sub, dims = 1:n.pcs, reduction = red.base)
          subset.obj(sub)
          updateRadioButtons(session, "display.mode", selected = "subset")
          # RunUMAP stores the result in "umap" — point the view there
          updateSelectInput(session, "view.red", selected = "umap")
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
    observeEvent(input$run.clustering, { 
      req(input$display.mode, input$view.red)
      
      # Normalize input once
      selected.view <- tolower(input$view.red)
      res <- input$cluster.res 
      
      tryCatch({
        withProgress(message = 'Running Clustering...', value = 0, {
          # Identify active object
          is.full <- input$display.mode == "full"
          obj <- if(is.full) current.obj() else subset.obj()
          req(obj)
          
          # Determine Graph Name
          target.graph <- if (selected.view == "umap") "snn" else gsub("umap\\.", "snn.", selected.view)
          
          # Validate graph exists
          if (!target.graph %in% names(obj@graphs)) {
            stop(paste("Graph", target.graph, "not found. Please run Neighbors first."))
          }
          
          # Run Seurat Clustering
          obj <- FindClusters(obj, resolution = res, graph.name = target.graph) 
          
          # Cleanup metadata (keeps only the final 'seurat_clusters')
          # This regex handles snn_res and snn.harmony_res
          pat <- paste0(target.graph, "_res\\.")
          obj@meta.data <- obj@meta.data[, !grepl(pat, colnames(obj@meta.data)), drop = FALSE] 
          
          # Update state
          if (is.full) current.obj(obj) else subset.obj(obj)
          
          showNotification(paste("Success: Clustered using", target.graph), type = "message")
        })
      }, error = function(e) {
        showNotification(paste("Clustering Error:", e$message), type = "error", duration = 10)
      })
    })
    
    # Render plot
    output$umap.main <- plotly::renderPlotly({
      # 1. Determine which object to use based on Display Mode
      if (input$display.mode == "subset") {
        # Check if subset exists; if not, show a message instead of an error
        validate(need(!is.null(subset.obj()), "No subset calculated yet. Please select cells and click 'Run UMAP'."))
        obj <- subset.obj()
      } else {
        req(current.obj())
        obj <- current.obj()
      }
      
      # 2. Basic requirements for plotting
      req(input$view.red, input$view.meta)
      # Ensure the selected reduction actually exists in the active object
      # (Important because subset objects only have the new 'umap' reduction)
      validate(need(input$view.red %in% names(obj@reductions), 
                    "The selected reduction is not available for this view."))
      
      # 3. Expression or Metadata
      if (input$viz.mode == "Expression") {
        # Determine which gene to use (Manual Entry takes priority)
        gene.to.plot <- if (nzchar(input$gene.entered)) {
          input$gene.entered 
        } else {
          input$gene.selected
        }
        req(gene.to.plot) 
        # CHECK IF GENE EXISTS: Provides a helpful message if there's a typo
        validate(
          need(gene.to.plot %in% rownames(obj), 
               paste0("Gene '", gene.to.plot, "' not found in this dataset. Check spelling/case."))
        )
        p <- FeaturePlot(obj, features = gene.to.plot, reduction = input$view.red)
      } else {
        req(input$view.meta)
        # Verify metadata column exists in the active object
        validate(need(input$view.meta %in% colnames(obj@meta.data), "Metadata column not found."))
        meta.val <- obj[[input$view.meta]][, 1]
        if(is.numeric(meta.val)){
          p <- FeaturePlot(obj, reduction = input$view.red, features = input$view.meta)
        }else{
          p <- DimPlot(obj, reduction = input$view.red, group.by = input$view.meta)
        }
      }
      # Carry barcodes for selection
      p$layers[[1]]$mapping$customdata <- aes(label = rownames(p$data))$label
      # plot
      plotly::ggplotly(p + theme_minimal(), source = "umap_select") %>% 
        layout(dragmode = "lasso")
    })
    
    
    
    # Updated Manual Annotation with Sync Support
    observeEvent(input$update.label, {
      select.data <- plotly::event_data("plotly_selected", source = "umap_select")
      req(select.data, input$new.label.name != "")
      
      selected.barcodes <- unique(as.character(select.data$customdata))
      new.label <- input$new.label.name
      
      # --- UPDATE MASTER OBJECT ---
      master <- isolate(current.obj())
      
      # Initialize column if missing
      if (!"manual_annotation" %in% colnames(master@meta.data)) {
        master$manual_annotation <- "Unlabeled"
      }
      
      master@meta.data[selected.barcodes, "manual_annotation"] <- new.label
      current.obj(master) # Push update to master
      
      # --- UPDATE SUBSET OBJECT (If it exists) ---
      if (!is.null(subset.obj())) {
        sub <- isolate(subset.obj())
        
        # Ensure column exists in subset too
        if (!"manual_annotation" %in% colnames(sub@meta.data)) {
          sub$manual_annotation <- "Unlabeled"
        }
        
        # Only update barcodes that are actually in this subset
        cells.in.sub <- intersect(selected.barcodes, colnames(sub))
        
        if (length(cells.in.sub) > 0) {
          sub@meta.data[cells.in.sub, "manual_annotation"] <- new.label
          subset.obj(sub) # Push update to subset so the plot refreshes
        }
      }
      
      showNotification(paste("Applied label:", new.label), type = "message")
    })
    
    # --- 3. DE Analysis ---
    marker.data <- eventReactive(input$run.de, {
      
      # 1. Capture Selection from Plotly
      # Note: We use the source "umap_select" defined in your ggplotly call
      selected.de <- plotly::event_data("plotly_selected", source = "umap_select")
      
      # Validation: Ensure something was selected
      validate(
        need(!is.null(selected.de), "Please use the Lasso or Box tool to select cells on the UMAP first.")
      )
      
      # Extract barcodes from the selection
      selected.cells <- unique(as.character(selected.de$customdata))
      
      # 2. Determine the Universe (Match your existing logic)
      # In your code, you use current.obj() and subset.obj()
      universe.data <- if (input$display.mode == "full") {
        current.obj() 
      } else {
        subset.obj()
      }
      
      req(universe.data) # Ensure object exists
      
      # 3. Define Groups
      all.cells <- colnames(universe.data)
      
      # Ensure selected cells actually exist in the current universe 
      # (Prevents errors if switching modes with an old selection)
      selected.cells <- intersect(selected.cells, all.cells)
      
      rest.cells <- setdiff(all.cells, selected.cells)
      
      # 4. Validation: Check group sizes
      validate(
        need(length(selected.cells) > 3, "Please select at least 4 cells for DE."),
        need(length(rest.cells) > 3, "Not enough background cells for comparison.")
      )
      
      # 5. Run DE Analysis with Progress Bar
      withProgress(message = 'Calculating DE Markers...', value = 0, {
        
        incProgress(0.2, detail = "Preparing groups...")
        showNotification("Starting Differential Expression...", type = "message")
        
        # Run FindMarkers
        de.results <- tryCatch({
          incProgress(0.5, detail = "Running FindMarkers (Wilcoxon)...")
          Seurat::FindMarkers(
            object = universe.data,
            ident.1 = selected.cells,
            ident.2 = rest.cells,
            only.pos = TRUE,
            min.pct = 0.1,
            logfc.threshold = 0.25
          )
        }, error = function(e) {
          showNotification(paste("Error in DE:", e$message), type = "error")
          return(NULL)
        })
        
        incProgress(0.2, detail = "Finalizing table...")
        
        if (is.null(de.results) || nrow(de.results) == 0) {
          showNotification("No significant markers found.", type = "warning")
          return(NULL)
        }
        
        showNotification("DE Analysis Complete!", type = "default")
        
        # Return results with Gene column
        return(cbind(Gene = rownames(de.results), de.results))
      })
    })
    
    # Render the DT Table
    output$marker.table <- DT::renderDT({
      # Use req() so the table is empty until the button is pressed
      req(marker.data())
      
      DT::datatable(
        marker.data(),
        rownames = FALSE,
        options = list(pageLength = 10, scrollX = TRUE),
        selection = "single"
      )
    })
    
    # Switch to DE analysis sub-tab
    observeEvent(input$run.de, {
      updateTabsetPanel(session, "main.tabs", selected = "DE Analysis")
    })
    
    # --- EXPORT ---
    mod.save.config.server("save", current.obj)
    
    return(current.obj)
  })
}