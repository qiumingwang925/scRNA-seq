mod_explore_umap_ui <- function(id) {
  ns <- NS(id)
  tabPanel("UMAP Cell-Type",
    tabsetPanel(
      # Sub-tab 1: Interactive Selection
      tabPanel("Interactive Selection",
               sidebarLayout(
                 sidebarPanel(
                   # Inside your sidebarPanel:
                   radioButtons(ns("display_mode"), "Display UMAP from:",
                                choices = c("All cells" = "full", "Selected cells (subset)" = "subset"),
                                selected = "full"),
                   
                   # Conditional Panel: Only show the Run button if 'subset' mode is chosen
                   conditionalPanel(
                     condition = sprintf("input['%s'] == 'subset'", ns("display_mode")),
                     helpText("1. Lasso select cells in 'All cells' mode."),
                     helpText("2. Switch to 'Selected cells' mode."),
                     actionButton(ns("run_subset_umap"), "📊 Run UMAP on Selection", 
                                  class = "btn-warning", style="width:100%"),
                     hr()
                   ),
                   
                   actionButton(ns("reset_umap"), "🔄 Reset Everything", 
                                class = "btn-danger", style="width:100%"),
                   br(),
                   verbatimTextOutput(ns("cell_stat_box"))
                 ),
                 mainPanel(
                   plotlyOutput(ns("umap_interactive"), height = "600px")
                 )
               )
      ),
      # Sub-tab 2: Highlight View
      tabPanel("Highlight View",
               sidebarLayout(
                 sidebarPanel(
                   selectInput(ns("select_idents"), "Select Cell Type(s) to Highlight:", 
                               choices = NULL, multiple = TRUE),
                   checkboxInput(ns("show_all"), "Select All", value = TRUE)
                 ),
                 mainPanel(
                   plotOutput(ns("umap_static"), height = "600px")
                 )
               )
      )
    )
  )
}

mod_explore_umap_server <- function(id, shared_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # 1. Reactive Value to store the 'Active' object
    active_obj <- reactiveVal(NULL)
    
    observe({
      req(shared_data())
      active_obj(shared_data())
      updateSelectInput(session, "idents_to_show", 
                        choices = levels(shared_data()), 
                        selected = levels(shared_data()))
    })
    
    # Helper: Detect the best Reduction and Graph
    # This looks for the first reduction containing 'umap' and first graph containing 'snn'
    detect_slots <- reactive({
      req(active_obj())
      obj <- active_obj()
      
      red_name <- grep("umap", Reductions(obj), value = TRUE, ignore.case = TRUE)[1]
      graph_name <- grep("snn", Graphs(obj), value = TRUE, ignore.case = TRUE)[1]
      
      # Fallback if specific names aren't found
      if(is.na(red_name)) red_name <- "umap" 
      
      list(reduction = red_name, graph = graph_name)
    })
    ######################
    # Store the subsetted object separately from the main object
    subset_obj <- reactiveVal(NULL)
    # Store cell IDs from the lasso
    selected_cell_ids <- reactiveVal(NULL)
    
    # 1. Capture the Lasso Selection
    observeEvent(event_data("plotly_selected", source = "umap_lasso"), {
      sel <- event_data("plotly_selected", source = "umap_lasso")
      if(!is.null(sel)) {
        selected_cell_ids(sel$customdata)
        showNotification(paste(length(sel$customdata), "cells selected. Switch to 'subset' mode to compute."), type = "message")
      }
    })
    
    # 2. Trigger the Re-UMAP Calculation
    observeEvent(input$run_subset_umap, {
      req(selected_cell_ids())
      
      # Visual feedback for heavy computation
      id_notif <- showNotification("Computing Subset UMAP via Graph...", duration = NULL, type = "message")
      
      tryCatch({
        # Generate subset from the original shared_data
        sub <- subset(shared_data(), cells = selected_cell_ids())
        
        # Detect Graph (e.g., snn.fastmnn or RNA_snn)
        slots <- detect_slots(sub) # Helper function to find 'snn'
        
        # Run UMAP on the graph (Fastest method)
        sub <- RunUMAP(sub, graph = slots$graph, reduction.name = "umap_subset")
        
        subset_obj(sub)
        removeNotification(id_notif)
        showNotification("Subset UMAP Complete!", type = "default")
        
      }, error = function(e) {
        removeNotification(id_notif)
        showNotification(paste("Error during UMAP:", e$message), type = "error")
      })
    })
    
    # --- SUB-TAB 1: INTERACTIVE ---
    
    output$umap_interactive <- renderPlotly({
      req(active_obj())
      obj <- active_obj()
      slots <- detect_slots()
      
      # Use the detected reduction
      p <- DimPlot(obj, reduction = slots$reduction)
      
      # Map cell names for lasso selection
      p$data$cell_id <- rownames(p$data)
      
      ggplotly(p, source = "umap_lasso") %>%
        layout(dragmode = "lasso") %>%
        event_register("plotly_selected")
    })
    
    observeEvent(event_data("plotly_selected", source = "umap_lasso"), {
      sel <- event_data("plotly_selected", source = "umap_lasso")
      req(sel)
      
      selected_cells <- sel$customdata
      slots <- detect_slots()
      
      if(length(selected_cells) > 10) {
        showNotification(paste("Re-computing UMAP using graph:", slots$graph), type = "message")
        
        # Subset from the ORIGINAL full dataset to maintain data integrity
        sub <- subset(shared_data(), cells = selected_cells)
        
        # Re-UMAP using the detected graph
        # We overwrite the default 'umap' slot or create a new one based on your preference
        sub <- RunUMAP(sub, 
                       graph = slots$graph, 
                       reduction.name = slots$reduction, 
                       reduction.key = "UMAP_")
        
        active_obj(sub)
      }
    })
    
    observeEvent(input$reset_umap, { active_obj(shared_data()) })
    
    # --- SUB-TAB 2: HIGHLIGHT ---
    
    output$umap_static <- renderPlot({
      req(active_obj())
      obj <- active_obj()
      slots <- detect_slots()
      
      idents <- if(input$select_all) levels(obj) else input$idents_to_show
      req(length(idents) > 0)
      
      if(length(idents) < length(levels(obj))) {
        DimPlot(obj, 
                reduction = slots$reduction,
                cells.highlight = WhichCells(obj, idents = idents),
                cols.highlight = "royalblue", 
                cols = "lightgrey", 
                sizes.highlight = 1) + 
          theme_minimal()
      } else {
        DimPlot(obj, reduction = slots$reduction) + theme_minimal()
      }
    })
    
    output$cell_count <- renderText({
      req(active_obj())
      paste0("Current Cell Count: ", ncol(active_obj()), 
             " | Reduction: ", detect_slots()$reduction)
    })
  })
}