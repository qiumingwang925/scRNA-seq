# UI Function
mod.upload.merge.ui <- function(id) {
  ns <- NS(id)
  tabPanel("Upload & Merge",
           sidebarLayout(
             sidebarPanel(
               fileInput(ns("file.input"), "Upload Seurat Objects (.rds)", 
                         multiple = TRUE, accept = ".rds"),
               hr(),
               actionButton(ns("merge.btn"), "Merge Objects", class = "btn-primary"),
               #br(), br(),
               #verbatimTextOutput(ns("status"))
             ),
             mainPanel(
               h4("Metadata Assignment"),
               helpText("Double-click cells to edit Batch and Group names."),
               DT::DTOutput(ns("meta.table")),
               br(), br(),
               h4("Merge Status"),
               verbatimTextOutput(ns("status"))
             )
           )
  )
}

# Server Function
mod.upload.merge.server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # This reactiveVal stores the user-editable metadata
    meta.df <- reactiveVal(NULL)
    
    # A. When files are uploaded, generate the initial table
    observeEvent(input$file.input, {
      df <- data.frame(
        File.Name = input$file.input$name,
        Batch = paste0("Batch.", 1:nrow(input$file.input)),
        Group = "Control",
        stringsAsFactors = FALSE
      )
      meta.df(df)
    })
    
    # B. Render the editable table
    output$meta.table <- DT::renderDT({
      req(meta.df())
      DT::datatable(meta.df(), editable = 'cell', rownames = FALSE, 
                    options = list(dom = 't', paging = FALSE))
    })
    
    # C. Critical: Save the edits the user makes in the UI
    observeEvent(input$meta.table_cell_edit, {
      info <- input$meta.table_cell_edit
      updated.df <- meta.df()
      # info$col + 1 because JS is 0-indexed and R is 1-indexed
      updated.df[info$row, info$col + 1] <- info$value
      meta.df(updated.df)
    })
    
    # D. The Merge Logic
    merge.data <- eventReactive(input$merge.btn, {
      req(input$file.input, meta.df())
      
      files <- input$file.input$datapath
      meta <- meta.df()
      
      # Load and annotate each object with the custom table values
      list.of.objects <- lapply(1:length(files), function(i) {
        obj <- readRDS(files[i])
        
        # Add the custom metadata from the table
        obj$batch <- meta$Batch[i]
        obj$group <- meta$Group[i]
        
        
        return(obj)
      })
      
      # Combine using Seurat's merge function
      if(length(list.of.objects) > 1) {
        merge.obj <- merge(
          x = list.of.objects[[1]], 
          y = list.of.objects[-1], 
          add.cell.ids = meta$Batch # Use your custom batch names as prefixes
        )
      } else {
        merge.obj <- list.of.objects[[1]]
      }
      
      # Essential for Seurat v5 to prevent layer errors later
      merge.obj <- JoinLayers(merge.obj)
      
      return(merge.obj)
    })
    
    output$status <- renderPrint({
      req(merge.data())
      obj <- merge.data()
      
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
    
    return(merge.data)
  })
}