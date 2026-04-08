# UI Function
mod_upload_merge_ui <- function(id) {
  ns <- NS(id)
  tabPanel("Upload & Merge",
           sidebarLayout(
             sidebarPanel(
               fileInput(ns("file_input"), "Upload Seurat Objects (.rds)", 
                         multiple = TRUE, accept = ".rds"),
               hr(),
               actionButton(ns("merge_btn"), "Merge Objects", class = "btn-primary"),
               #br(), br(),
               #verbatimTextOutput(ns("status"))
             ),
             mainPanel(
               h4("Metadata Assignment"),
               helpText("Double-click cells to edit Batch and Group names."),
               DT::DTOutput(ns("meta_table")),
               br(), br(),
               h4("Merge Status"),
               verbatimTextOutput(ns("status"))
             )
           )
  )
}

# Server Function
mod_upload_merge_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # This reactiveVal stores the user-editable metadata
    meta_df <- reactiveVal(NULL)
    
    # A. When files are uploaded, generate the initial table
    observeEvent(input$file_input, {
      df <- data.frame(
        File_Name = input$file_input$name,
        Batch = paste0("Batch_", 1:nrow(input$file_input)),
        Group = "Control",
        stringsAsFactors = FALSE
      )
      meta_df(df)
    })
    
    # B. Render the editable table
    output$meta_table <- DT::renderDT({
      req(meta_df())
      DT::datatable(meta_df(), editable = 'cell', rownames = FALSE, 
                    options = list(dom = 't', paging = FALSE))
    })
    
    # C. Critical: Save the edits the user makes in the UI
    observeEvent(input$meta_table_cell_edit, {
      info <- input$meta_table_cell_edit
      updated_df <- meta_df()
      # info$col + 1 because JS is 0-indexed and R is 1-indexed
      updated_df[info$row, info$col + 1] <- info$value
      meta_df(updated_df)
    })
    
    # D. The Merge Logic
    merge_data <- eventReactive(input$merge_btn, {
      req(input$file_input, meta_df())
      
      files <- input$file_input$datapath
      meta <- meta_df()
      
      # Load and annotate each object with the custom table values
      list_of_objects <- lapply(1:length(files), function(i) {
        obj <- readRDS(files[i])
        
        # Add the custom metadata from the table
        obj$batch <- meta$Batch[i]
        obj$group <- meta$Group[i]
        
        
        return(obj)
      })
      
      # Combine using Seurat's merge function
      if(length(list_of_objects) > 1) {
        merge_obj <- merge(
          x = list_of_objects[[1]], 
          y = list_of_objects[-1], 
          add.cell.ids = meta$Batch # Use your custom batch names as prefixes
        )
      } else {
        merge_obj <- list_of_objects[[1]]
      }
      
      # Essential for Seurat v5 to prevent layer errors later
      merge_obj <- JoinLayers(merge_obj)
      
      return(merge_obj)
    })
    
    output$status <- renderPrint({
      req(merge_data())
      obj <- merge_data()
      
      cat("--- Merge Successful ---\n")
      cat("Total Cells:   ", ncol(obj), "\n")
      cat("Total Features:", nrow(obj), "\n\n")
      
      # 1. Cells per Original Sample
      cat("Cells per Original Sample:\n")
      print(table(obj$orig.ident))
      cat("\n")
      
      # 2. Cells per Custom Batch
      cat("Cells per Custom Batch (Grouped):\n")
      print(table(obj$batch))
      cat("\n")
      
      # 3. Cells per Custom Group
      cat("Cells per Custom Group:\n")
      print(table(obj$group))
      
      cat("\n--- Full Seurat Object Summary ---\n")
      print(obj)
    })
    
    return(merge_data)
  })
}