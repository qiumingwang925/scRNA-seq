mod_import_ui <- function(id) {
  ns <- NS(id)
  
  tabPanel(
    "Import raw data",
    fluidRow(
      column(10, strong("Upload a sample (Cell Ranger MEX Analysis Result Folder)"))
    ),
    fluidRow(
      column(1, shinyDirButton(ns("folder"), "Select", title ="Select a sample folder", multiple = FALSE, class = "btn-success")),
      column(2, verbatimTextOutput(ns("mex"), placeholder = FALSE)),
      column(3, textOutput(ns("cell.count")))
    ),
    fluidRow(
      column(10, strong("Convert to Seurat Object"))
    ),
    fluidRow(
      column(2, textInput(ns("project.name"), "Create a Sample ID", value = "")),
      column(2, numericInput(ns("min.cells"), "Cell Threshold", value = 3)),
      column(2, numericInput(ns("min.features"), "Feature Threshold", value = 100)),
      column(1, actionButton(ns("convert"), "Convert", class ="btn-success"))
    ),
    fluidRow(
      column(10, plotOutput(ns("plot.raw.vln"), height = "400px", width = "1000px"))
    ),
    fluidRow(
      column(10, plotOutput(ns("plot.raw.dst"), height = "150px", width = "1000px"))
    )
  )
}


mod_import_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    shinyjs::disable("convert")

    # Setup shinyFiles
    volumes <- c("Project" = dirname(getwd()), shinyFiles::getVolumes()())
    shinyDirChoose(input, 'folder', roots = volumes, defaultRoot = "Project")

    # Display selected folder name
    observe({
      req(input$folder, !is.integer(input$folder))
      folder_name <- tail(unlist(input$folder[1], use.names = FALSE), 1)
      output$mex <- renderText({ folder_name })
    })

    # Auto-update project name and enable Convert
    observeEvent(input$folder, {
      req(input$folder, !is.integer(input$folder))
      folder_name <- tail(unlist(input$folder[1], use.names = FALSE), 1)
      updateTextInput(session, "project.name", value = folder_name)
      shinyjs::enable("convert")
    })
    
    # Convert to Seurat object
    seurat_obj <- eventReactive(input$convert, {
      req(input$folder)
      path <- shinyFiles::parseDirPath(volumes, input$folder)
      srt <- CreateSeuratObject(
        counts = Read10X(data.dir = path),
        project = input$project.name,
        min.cells = input$min.cells,
        min.features = input$min.features
      )
      srt[["percent.mt"]] <- PercentageFeatureSet(srt, pattern = "^mt-")
      srt[["percent.rp"]] <- PercentageFeatureSet(srt, pattern = "^Rp[sl]")
      srt[["percent.hb"]] <- PercentageFeatureSet(srt, pattern = "^Hb[^(P)]")
      srt
    })
    
    # Violin plot
    plot_vln <- reactive({
      req(seurat_obj())
      VlnPlot(seurat_obj(), features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp", "percent.hb"), ncol = 5)
    })
    
    output$plot.raw.vln <- renderPlot({ plot_vln() }, res = 96)
    
    # Density plots
    plot_dst <- reactive({
      req(seurat_obj())
      df <- seurat_obj()@meta.data
      p1 <- ggplot(df, aes(x = nFeature_RNA)) + geom_density() + theme_bw()
      p2 <- ggplot(df, aes(x = nCount_RNA)) + geom_density() + theme_bw()
      p3 <- ggplot(df, aes(x = percent.mt)) + geom_density() + theme_bw()
      p4 <- ggplot(df, aes(x = percent.rp)) + geom_density() + theme_bw()
      p5 <- ggplot(df, aes(x = percent.hb)) + geom_density() + theme_bw()
      ggpubr::ggarrange(p1, p2, p3, p4, p5, ncol = 5, nrow = 1)
    })
    
    output$plot.raw.dst <- renderPlot({ plot_dst() }, res = 96)
    
    # Total cell count
    output$cell.count <- renderText({
      req(seurat_obj())
      paste0("Total cell counts: ", ncol(seurat_obj()))
    })
    
    return(seurat_obj)
  })
}