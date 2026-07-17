# 1. Install the GNU Fortran Compiler
#The easiest way to get the exact libraries R is looking for is to download the official installer provided by the R-MacOS project.
#Go to the https://mac.r-project.org/tools/.
#Download the gfortran SDK that matches your macOS version (likely the gfortran-14.2-universal.pkg ).
#Run the .pkg installer. This will create the /opt/gfortran directory that the error message is complaining about.

#2. Update your Makevars (The "Secret Sauce")
#R uses a hidden file called .Makevars to know where compilers live. If the step above doesn't fix it immediately, you need to tell R where to find the libraries.
#Run these lines inside your R console:
# Create the directory if it doesn't exist
#if (!dir.exists("~/.R")) dir.create("~/.R")
# Add the correct paths to your Makevars file
#cat("
#FC = /opt/gfortran/bin/gfortran
#F77 = /opt/gfortran/bin/gfortran
#FLIBS = -L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/14.2.0 -L/opt/gfortran/lib -lgfortran -lquadmath
#", file = "~/.R/Makevars", append = TRUE)

#3. Install "lisi"
#devtools::install_github("immunogenomics/lisi")



mod.benchmark.ui <- function(id) {
  ns <- NS(id)
  tabPanel("Benchmarking",
           sidebarLayout(
             sidebarPanel(
               tags$h4("Evaluation Settings"),
               hr(),
               selectInput(ns("batch.label"), "Batch Label:", choices = NULL),
               selectInput(ns("celltype.label"), "Cell-Type Label:", choices = NULL),
               checkboxGroupInput(ns("eval.metrics"), "Metrics:",
                                  choices = c("ASW", "LISI", "GraphLISI"),
                                  selected = c("ASW", "LISI")),
               hr(),
               actionButton(ns("run.eval"), "Compute Scores",
                            class = "btn-primary btn-block", style = "color: white;")
             ),

             mainPanel(
               tabsetPanel(
                 tabPanel("Rank Summary",
                          plotOutput(ns("rank.barplot"), height = "500px"),
                          hr(),
                          tags$h4("Rank Summary"),
                          tableOutput(ns("rank.table"))
                 ),
                 tabPanel("Score Distributions",
                          plotOutput(ns("violin.plots"), height = "500px"),
                          hr(),
                          tags$h4("Median Score Summary"),
                          tableOutput(ns("score.table"))
                 )
               )
             )
           )
  )
}

mod.benchmark.server <- function(id, shared.data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Update selectors based on metadata
    observe({
      req(shared.data())
      meta <- shared.data()@meta.data
      cat.cols <- colnames(meta)[sapply(meta, function(x) is.character(x) || is.factor(x))]

      batch.choices <- cat.cols[grepl("batch", cat.cols, ignore.case = TRUE)]
      ct.keywords <- "TMS|LungMap|manual|annotation|celltype"
      celltype.choices <- cat.cols[grepl(ct.keywords, cat.cols, ignore.case = TRUE)]

      updateSelectInput(session, "batch.label", choices = c("None", batch.choices))
      updateSelectInput(session, "celltype.label", choices = c("None", celltype.choices))
    })

    # --- CORE COMPUTATION LOGIC ---
    eval.results <- eventReactive(input$run.eval, {
      req(shared.data(), input$batch.label != "None" | input$celltype.label != "None")
      obj <- shared.data()

      all.reds <- names(obj@reductions)
      target.reds <- all.reds[grepl("integrated|pca", all.reds) & !grepl("umap", all.reds)]

      detailed.scores <- data.frame()
      summary.stats <- data.frame()

      withProgress(message = 'Calculating Metrics...', value = 0, {
        for(red in target.reds) {
          #######################################
          #message("Currently evaluating: ", red)
          #emb <- Embeddings(obj, reduction = red)

          # Check for dimension mismatch!
          #if (nrow(emb) != nrow(obj@meta.data)) {
            #message("Warning: Skipping ", red, " due to cell count mismatch.")
            #next
          #}

          # Standardize to 30 dims (RPCA might have fewer)
          #dims_to_use <- min(30, ncol(emb))
          #emb <- emb[, 1:dims_to_use]
          ###################################
          incProgress(1/length(target.reds), detail = paste("Processing", red))

          # Clean Names for Plotting
          display.name <- if(red == "pca") "Unintegrated" else {
            clean <- gsub("integrated\\.", "", red)
            paste0(toupper(substring(clean, 1, 1)), substring(clean, 2))
          }
          # --- PCA embedding (ASW & LISI) ---
          emb <- Embeddings(obj, reduction = red)[, 1:30]

          # --- 1. ASW Block ---
          if("ASW" %in% input$eval.metrics) {

            dists <- dist(emb)
            if(input$batch.label != "None") {
              sw <- silhouette(as.numeric(factor(obj@meta.data[[input$batch.label]])), dists)[,3]
              detailed.scores <- rbind(detailed.scores, data.frame(Cell=rownames(emb), Score=sw, Integration=display.name, Metric="basw"))
              summary.stats <- rbind(summary.stats, data.frame(Integration=display.name, Method="basw", Median=median(sw)))
            }
            if(input$celltype.label != "None") {
              sw <- silhouette(as.numeric(factor(obj@meta.data[[input$celltype.label]])), dists)[,3]
              detailed.scores <- rbind(detailed.scores, data.frame(Cell=rownames(emb), Score=sw, Integration=display.name, Metric="casw"))
              summary.stats <- rbind(summary.stats, data.frame(Integration=display.name, Method="casw", Median=median(sw)))
            }
          }

          # --- 2. LISI Block ---
          if("LISI" %in% input$eval.metrics) {

            if(input$batch.label != "None") {
              lisi.res <- compute_lisi(emb, obj@meta.data, input$batch.label)
              detailed.scores <- rbind(detailed.scores, data.frame(Cell=rownames(emb), Score=lisi.res[[input$batch.label]], Integration=display.name, Metric="ilisi"))
              summary.stats <- rbind(summary.stats, data.frame(Integration=display.name, Method="ilisi", Median=median(lisi.res[[input$batch.label]])))
            }
            if(input$celltype.label != "None") {
              lisi.res <- compute_lisi(emb, obj@meta.data, input$celltype.label)
              detailed.scores <- rbind(detailed.scores, data.frame(Cell=rownames(emb), Score=lisi.res[[input$celltype.label]], Integration=display.name, Metric="clisi"))
              summary.stats <- rbind(summary.stats, data.frame(Integration=display.name, Method="clisi", Median=median(lisi.res[[input$celltype.label]])))
            }
          }

          # --- 3. GraphLISI Block (High Performance Sparse) ---
          if("GraphLISI" %in% input$eval.metrics) {

            target.suffix <- gsub("integrated\\.", "", red)
            # --- NAME MAPPING ---
            # If the reduction is 'pca', the graph is usually just 'snn'
            if(target.suffix == "pca") {
              graph.nm <- "snn"
            } else {
              graph.nm <- paste0("snn.", target.suffix)
            }

            snn.matrix <- obj@graphs[[graph.nm]]

            # Efficient neighbor extraction
            neighbors <- lapply(1:ncol(snn.matrix), function(i) {
              snn.matrix@i[(snn.matrix@p[i] + 1):snn.matrix@p[i + 1]] + 1
            })

            # Helper for local diversity
            get.isi <- function(labels.vec, n.list) {
              sapply(n.list, function(idx) {
                if(length(idx) == 0) return(NA)
                p <- table(labels.vec[idx]) / length(idx)
                1 / sum(p^2)
              })
            }

            if(input$batch.label != "None") {
              scores <- get.isi(obj[[input$batch.label, drop=TRUE]], neighbors)
              detailed.scores <- rbind(detailed.scores, data.frame(Cell=rownames(emb), Score=scores, Integration=display.name, Metric="gilisi"))
              summary.stats <- rbind(summary.stats, data.frame(Integration=display.name, Method="gilisi", Median=median(scores, na.rm=TRUE)))
            }
            if(input$celltype.label != "None") {
              scores <- get.isi(obj[[input$celltype.label, drop=TRUE]], neighbors)
              detailed.scores <- rbind(detailed.scores, data.frame(Cell=rownames(emb), Score=scores, Integration=display.name, Metric="gclisi"))
              summary.stats <- rbind(summary.stats, data.frame(Integration=display.name, Method="gclisi", Median=median(scores, na.rm=TRUE)))
            }
          }
        }
      })

      # --- RANKING LOGIC ---
      list.metrics <- unique(summary.stats$Method)
      rank.summary <- data.frame()

      for(m in list.metrics){
        rank.df <- summary.stats %>% filter(Method == m)
        # ilisi, gilisi, casw: Higher is Better. basw, clisi, gclisi: Lower is Better.
        # Orient so larger = better, then rank so the best gets the most points
        # and tied medians share the same rank (ties.method = "max").
        oriented <- if(m %in% c("ilisi", "gilisi", "casw")) rank.df$Median else -rank.df$Median
        rank.df$Rank <- rank(oriented, ties.method = "max")
        rank.summary <- rbind(rank.summary, rank.df)
      }

      return(list(details = detailed.scores, summary = rank.summary))
    })

    # --- VISUALIZATIONS ---

    output$violin.plots <- renderPlot({
      df <- eval.results()$details
      req(nrow(df) > 0)
      cols.10 <- c("#1F77B4", "#FF7F0E", "#2CA02C","#D62728","#9467BD", "#8C564B", "#E377C2","#7F7F7F", "#BCBD22", "#17BECF")

      # Lay out each metric family (ASW, LISI, GraphLISI) as its own column with
      # the batch metric on top and the cell-type metric below. facet_wrap fills
      # row-major, so order the batch row first, then the cell-type row, and set
      # ncol to the number of families present so subsets stay aligned.
      metric.family <- c(basw = "ASW", ilisi = "LISI", gilisi = "GraphLISI",
                         casw = "ASW", clisi = "LISI", gclisi = "GraphLISI")
      metric.order <- c("basw", "ilisi", "gilisi", "casw", "clisi", "gclisi")
      present <- metric.order[metric.order %in% df$Metric]
      df$Metric <- factor(df$Metric, levels = present)
      n.families <- length(unique(metric.family[present]))

      ggplot(df, aes(x = Integration, y = Score, fill = Integration)) +
        geom_violin(trim = FALSE, show.legend = FALSE) +
        geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, show.legend = FALSE) +
        facet_wrap(~Metric, scales = "free_y", ncol = n.families) +
        scale_fill_manual(values = cols.10) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Per-cell Score Distributions", y = "Metric Score")
    })

    output$rank.barplot <- renderPlot({
      df <- eval.results()$summary
      req(nrow(df) > 0)

      df$Method.group <- ifelse(df$Method %in% c("basw", "ilisi", "gilisi"),
                                "Batch Mixing", "Cell Type Separation")

      cols.6 <- c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F","#EDC948")

      ggplot(df, aes(x = Rank, y = Integration, fill = Method)) +
        geom_bar(stat = "identity", color = "white", linewidth = 0.2) +
        facet_grid(. ~ Method.group) +
        scale_fill_manual(values = cols.6) +
        theme_classic() +
        labs(title = "Integration Performance Ranking",
             subtitle = "Higher Rank Score = Better Performance",
             x = "Cumulative Rank Score", y = "Integration Method")
    })
    # --- SCORE TABLE ---
    output$score.table <- renderTable({
      df.details <- eval.results()$details
      req(nrow(df.details) > 0)

      df.details %>%
        group_by(Integration, Metric) %>%
        summarize(Median.Score = median(Score, na.rm = TRUE), .groups = 'drop') %>%
        pivot_wider(names_from = Metric, values_from = Median.Score) %>%
        dplyr::rename(Method = Integration) %>%
        dplyr::select(Method, dplyr::any_of(c("basw", "casw", "ilisi", "clisi", "gilisi", "gclisi"))) %>%
        arrange(Method)
    }, digits = 3)

   # --- RANK TABLE ###
    output$rank.table <- renderTable({
      # 1. Get the data
      df.rank <- eval.results()$summary
      req(nrow(df.rank) > 0)

      # 2. Pivot to Wide format
      # Each integration method gets one row
      rank.pivot <- df.rank %>%
        dplyr::select(Integration, Method, Rank) %>%
        tidyr::pivot_wider(
          names_from = Method,
          values_from = Rank
        )

      # 3. Calculate Total Score safely
      # 'across(where(is.numeric))' automatically ignores the 'Integration' text column
      # and works perfectly whether you have 1 metric or 6 metrics.
      rank.pivot <- rank.pivot %>%
        mutate(Total.Score = rowSums(across(where(is.numeric)), na.rm = TRUE)) %>%
        arrange(desc(Total.Score))

      # 4. Final Column Ordering
      # any_of() ensures that if a metric wasn't calculated, the code doesn't break
      rank.pivot <- rank.pivot %>%
        dplyr::select(
          Integration,
          any_of(c("basw", "ilisi", "gilisi", "casw", "clisi", "gclisi")),
          Total.Score
        ) %>%
        dplyr::rename(Method = Integration)
      return(rank.pivot)
    }, digits = 0)

  })
}
