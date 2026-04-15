mod_explore_upload_ui <- function(id) {
  ns <- NS(id)
  tabPanel("Upload Data",
           sidebarLayout(
             sidebarPanel(
               fileInput(ns("file_in"), "Upload Seurat Object (.rds)", accept = ".rds"),
               hr(),
               helpText("Ensure object is slimmed down (SCT/RNA 'data' slot only).")
             ),
             mainPanel(
               h4("Object Summary Report"),
               tableOutput(ns("report_table")),
               verbatimTextOutput(ns("validation_msg"))
             )
           )
  )
}

mod_explore_upload_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    # 1. Reactive for the uploaded object
    raw_obj <- reactive({
      req(input$file_in)
      readRDS(input$file_in$datapath)
    })
    
    # 2. Validation and Reporting Logic
    report_data <- reactive({
      obj <- raw_obj()
      
      # Check components
      res <- data.frame(
        Feature = c("Assays", "Metadata 'group'", "Reductions", "Graphs", "Active Ident"),
        Status = c(
          paste(Assays(obj), collapse = ", "),
          paste(names(obj@meta.data), collapse = ", "),
          paste(Reductions(obj), collapse = ", "),
          paste(Graphs(obj), collapse = ", "),
          paste(unique(Idents(obj)), collapse = ", ")
        ),
        stringsAsFactors = FALSE
      )
      return(res)
    })
    
    # Render the table
    output$report_table <- renderTable({ report_data() })
    
    # 3. Final Reactive Object (with safety checks)
    processed_obj <- reactive({
      obj <- raw_obj()
      
      # Ensure 'data' slot exists in SCT or RNA
      main_assay <- ifelse("SCT" %in% Assays(obj), "SCT", "RNA")
      DefaultAssay(obj) <- main_assay
      
      # Verify if 'scale.data' is missing (as requested for slim version)
      # Note: Heatmap will need 'data' or 'scale.data'
      
      return(obj)
    })
    
    return(processed_obj)
  })
}