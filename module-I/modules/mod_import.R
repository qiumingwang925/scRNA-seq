## ABOUTME: Shiny module for importing Cell Ranger MEX data into a Seurat object.
## ABOUTME: Computes QC metrics (mito/ribo/hemoglobin %) and displays violin and density plots.

mod.import.ui <- function(id) {
  ns <- NS(id)
  
  tabPanel(
    "Import raw data", value = "tab.import",
    radioButtons(ns("import.mode"), "Load a sample using:", inline = TRUE,
      choices = c(
        "Option 1 — Cell Ranger MEX folder" = "mex",
        "Option 2 — Processed Seurat object (.rds)" = "rds"
      )
    ),
    conditionalPanel(
      condition = "input['import.mode'] == 'mex'", ns = ns,
      wellPanel(
        fluidRow(
          column(4, "Select a Cell Ranger MEX analysis result folder, then Convert"),
          column(1, shinyDirButton(ns("folder"), "Select", title = "Select a sample folder", multiple = FALSE, class = "btn-success", style = "width: 100px;"))
        ),
        fluidRow(
          column(4, verbatimTextOutput(ns("mex"), placeholder = FALSE))
        ),
        fluidRow(
          column(2, textInput(ns("project.name"), "Create a Sample ID", value = "")),
          column(2, numericInput(ns("min.cells"), "Cell Threshold", value = 3)),
          column(2, numericInput(ns("min.features"), "Feature Threshold", value = 100))
        ),
        fluidRow(
          column(2, actionButton(ns("convert"), "Convert", class = "btn-success", style = "width: 100px;"))
        )
      )
    ),
    conditionalPanel(
      condition = "input['import.mode'] == 'rds'", ns = ns,
      wellPanel(
        fluidRow(
          column(4, "Upload an already-processed Seurat object"),
          column(4, fileInput(ns("seurat_file"), NULL, accept = c(".rds")))
        )
      )
    ),
    wellPanel(
      strong("QC Overview"),
      textOutput(ns("cell.count")),
      shinycssloaders::withSpinner(plotOutput(ns("plot.raw.vln"), height = "400px", width = "1000px")),
      shinycssloaders::withSpinner(plotOutput(ns("plot.raw.dst"), height = "150px", width = "1000px"))
    )
  )
}


mod.import.server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Setup shinyFiles. In the Docker demo, SCNEXUS_DATA_ROOT points at the
    # bind-mounted datasets folder; restricting the browser to it keeps users
    # out of the rest of the container filesystem. Unset outside Docker, so
    # local use keeps the project dir + system volumes.
    data.root <- Sys.getenv("SCNEXUS_DATA_ROOT", unset = "")
    if (nzchar(data.root)) {
      volumes <- c("Datasets" = data.root)
      default.root <- "Datasets"
    } else {
      volumes <- c("Project" = dirname(getwd()), shinyFiles::getVolumes()())
      default.root <- "Project"
    }
    shinyDirChoose(input, 'folder', roots = volumes, defaultRoot = default.root)

    # Display selected folder name
    observe({
      req(input$folder, !is.integer(input$folder))
      folder.name <- tail(unlist(input$folder[1], use.names = FALSE), 1)
      output$mex <- renderText({ folder.name })
    })

    # Track whether import step is complete
    completed <- reactiveVal(FALSE)

    # Auto-update project name and enable Convert
    observeEvent(input$folder, {
      req(input$folder, !is.integer(input$folder))
      folder.name <- tail(unlist(input$folder[1], use.names = FALSE), 1)
      updateTextInput(session, "project.name", value = folder.name)
      completed(FALSE)
    })
    
    # Convert to Seurat object
    seurat.obj <- reactiveVal(NULL)
    observeEvent(input$convert, {
      req(input$folder)
      withProgress(message = "Converting to Seurat object...", value = 0, {
        path <- shinyFiles::parseDirPath(volumes, input$folder)
        setProgress(value = 0.2, detail = "Reading 10X data")
        counts <- Read10X(data.dir = path)
        setProgress(value = 0.5, detail = "Creating Seurat object")
        srt <- CreateSeuratObject(
          counts = counts,
          project = input$project.name,
          min.cells = input$min.cells,
          min.features = input$min.features
        )
        setProgress(value = 0.7, detail = "Computing QC metrics")
        srt[["percent.mt"]] <- PercentageFeatureSet(srt, pattern = "^mt-")
        srt[["percent.rp"]] <- PercentageFeatureSet(srt, pattern = "^Rp[sl]")
        srt[["percent.hb"]] <- PercentageFeatureSet(srt, pattern = "^Hb[^(P)]")
        setProgress(value = 1.0, detail = "Done")
      })
      seurat.obj(srt)
      completed(TRUE)
    })

    # Upload an already-processed Seurat object instead of converting a folder
    observeEvent(input$seurat_file, {
      req(input$seurat_file)
      withProgress(message = "Loading Seurat object...", value = 0, {
        srt <- readRDS(input$seurat_file$datapath)
        setProgress(value = 0.6, detail = "Computing missing QC metrics")
        # Mirror the mouse gene patterns used by the convert path; only fill in
        # metrics the uploaded object doesn't already carry.
        if (!"percent.mt" %in% colnames(srt@meta.data)) {
          srt[["percent.mt"]] <- PercentageFeatureSet(srt, pattern = "^mt-")
        }
        if (!"percent.rp" %in% colnames(srt@meta.data)) {
          srt[["percent.rp"]] <- PercentageFeatureSet(srt, pattern = "^Rp[sl]")
        }
        if (!"percent.hb" %in% colnames(srt@meta.data)) {
          srt[["percent.hb"]] <- PercentageFeatureSet(srt, pattern = "^Hb[^(P)]")
        }
        setProgress(value = 1.0, detail = "Done")
      })
      seurat.obj(srt)
      completed(TRUE)
    })

    # Violin plot
    output$plot.raw.vln <- renderPlot({
      req(completed())
      VlnPlot(seurat.obj(), features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp", "percent.hb"), ncol = 5)
    }, res = 96)

    # Density plots
    output$plot.raw.dst <- renderPlot({
      req(completed())
      df <- seurat.obj()@meta.data
      p1 <- ggplot(df, aes(x = nFeature_RNA)) + geom_density() + theme_bw()
      p2 <- ggplot(df, aes(x = nCount_RNA)) + geom_density() + theme_bw()
      p3 <- ggplot(df, aes(x = percent.mt)) + geom_density() + theme_bw()
      p4 <- ggplot(df, aes(x = percent.rp)) + geom_density() + theme_bw()
      p5 <- ggplot(df, aes(x = percent.hb)) + geom_density() + theme_bw()
      ggpubr::ggarrange(p1, p2, p3, p4, p5, ncol = 5, nrow = 1)
    }, res = 96)

    # Total cell count
    output$cell.count <- renderText({
      req(completed())
      paste0("Total cell counts: ", ncol(seurat.obj()))
    })
    
    return(list(seurat.obj = seurat.obj, completed = completed))
  })
}