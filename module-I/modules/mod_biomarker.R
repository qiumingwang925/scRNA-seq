## ABOUTME: Shiny module for visualizing gene expression from uploaded biomarker CSV lists.
## ABOUTME: Provides hierarchical gene selection (main type > sub type > gene) and UMAP feature plots.

mod.biomarker.ui <- function(id) {
  ns <- NS(id)
  
  tabPanel("Biomarker",
           fluidRow(
             column(4, strong("Upload Biomarkers (CSV)")),
             column(1, shinyFilesButton(ns("marker.upload"), "Select", title = "Select a biomarker CSV file", multiple = FALSE, class = "btn-success", style = "width: 100px;")),
             column(7, helpText("Expected columns: Label.main, Label.fine, Markers"))
           ),
           hr(),
           fluidRow(
             column(2, selectInput(ns("mainlabinput"), "Main cell type", choices = NULL)),
             column(2, selectInput(ns("finelabinput"), "Sub cell type", choices = NULL)),
             column(2, selectInput(ns("gene.input.1"), "Gene from list", choices = NULL)),
             column(2, textInput(ns("gene.input.2"), "OR type gene")),
             column(2, actionButton(ns("gene.input.run"), "Plot Expression", class = "btn-success", style="margin-top: 25px;"))
           ),
           fluidRow(
             column(12, plotOutput(ns("plot.gene.input"), height = "500px"))
           )
  )
}

mod.biomarker.server <- function(id, seurat.obj.cellcycle) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Setup shinyFiles for CSV selection
    volumes <- c("Project" = dirname(getwd()), shinyFiles::getVolumes()())
    shinyFileChoose(input, "marker.upload", roots = volumes, defaultRoot = "Project", filetypes = c("csv"))

    # Reactive storage for marker data
    # Defaults to the 'markers' object already in your environment
    current.markers <- reactiveVal(markers)

    # Handle file selection
    observeEvent(input$marker.upload, {
      req(input$marker.upload, !is.integer(input$marker.upload))
      csv.path <- shinyFiles::parseFilePaths(volumes, input$marker.upload)$datapath
      req(length(csv.path) > 0)
      df <- read.csv(csv.path, stringsAsFactors = FALSE)
      
      # Basic validation of columns
      required.cols <- c("Label.main", "Label.fine", "Markers")
      if(all(required.cols %in% colnames(df))) {
        current.markers(df)
        showNotification("Marker list updated successfully!", type = "message")
      } else {
        showNotification("Upload failed: CSV must have Label.main, Label.fine, and Markers columns.", type = "error")
      }
    })
    
    # 3. Update Dropdowns based on the current.markers()
    observe({
      df <- current.markers()
      updateSelectInput(session, "mainlabinput", choices = unique(df$Label.main))
    })
    
    observeEvent(input$mainlabinput, {
      df <- current.markers()
      sub.choices <- df$Label.fine[df$Label.main == input$mainlabinput]
      updateSelectInput(session, "finelabinput", choices = unique(sub.choices))
    })
    
    observeEvent(input$finelabinput, {
      df <- current.markers()
      gene.choices <- df$Markers[df$Label.fine == input$finelabinput]
      updateSelectInput(session, "gene.input.1", choices = unique(gene.choices))
    })
    
    # 4. Gene Input Logic
    geneinput <- eventReactive(input$gene.input.run, {
      if (input$gene.input.2 != "") return(input$gene.input.2)
      return(input$gene.input.1)
    })
    
    # 5. Render Plot
    output$plot.gene.input <- renderPlot({
      req(seurat.obj.cellcycle(), geneinput())
      srt <- seurat.obj.cellcycle()
      gene <- geneinput()
      
      if(!(gene %in% rownames(srt))) {
        showNotification(paste("Gene", gene, "not found."), type = "warning")
        return(NULL)
      }
      
      FeaturePlot(srt, features = gene, cols = c("lightgrey", "firebrick")) +
        theme(aspect.ratio = 1)
    }, res = 96)
  })
}