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



mod_benchmark_ui <- function(id) {
  ns <- NS(id)
  tabPanel("Benchmarking", 
           sidebarLayout(
             sidebarPanel(
               tags$h4("Evaluation Settings"),
               hr(),
               selectInput(ns("batch_label"), "Batch Label:", choices = NULL),
               selectInput(ns("celltype_label"), "Cell-Type Label:", choices = NULL),
               checkboxGroupInput(ns("eval_metrics"), "Metrics:",
                                  choices = c("ASW", "LISI", "GraphLISI"),
                                  selected = c("ASW", "LISI")),
               hr(),
               actionButton(ns("run_eval"), "Compute Scores", 
                            class = "btn-primary btn-block", style = "color: white;")
             ),
             
             mainPanel(
               tabsetPanel(
                 tabPanel("Rank Summary", 
                          plotOutput(ns("rank_barplot"), height = "500px"),
                          hr(),
                          tags$h4("Rank Summary"),
                          tableOutput(ns("rank_table"))
                 ),
                 tabPanel("Score Distributions", 
                          plotOutput(ns("violin_plots"), height = "500px"),
                          hr(),
                          tags$h4("Median Score Summary"),
                          tableOutput(ns("score_table"))
                 )
               )
             )
           )
  )
}

mod_benchmark_server <- function(id, shared_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Update selectors based on metadata
    observe({
      req(shared_data())
      meta <- shared_data()@meta.data
      cat_cols <- colnames(meta)[sapply(meta, function(x) is.character(x) || is.factor(x))]
      
      batch_choices <- cat_cols[grepl("batch", cat_cols, ignore.case = TRUE)]
      ct_keywords <- "TMS|LungMap|manual|annotation|celltype"
      celltype_choices <- cat_cols[grepl(ct_keywords, cat_cols, ignore.case = TRUE)]
      
      updateSelectInput(session, "batch_label", choices = c("None", batch_choices))
      updateSelectInput(session, "celltype_label", choices = c("None", celltype_choices))
    })
    
    # --- CORE COMPUTATION LOGIC ---
    eval_results <- eventReactive(input$run_eval, {
      req(shared_data(), input$batch_label != "None" | input$celltype_label != "None")
      obj <- shared_data()
      
      all_reds <- names(obj@reductions)
      target_reds <- all_reds[grepl("integrated|pca", all_reds) & !grepl("umap", all_reds)]
      
      detailed_scores <- data.frame() 
      summary_stats <- data.frame()   
      
      withProgress(message = 'Calculating Metrics...', value = 0, {
        for(red in target_reds) {
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
          incProgress(1/length(target_reds), detail = paste("Processing", red))
          
          # Clean Names for Plotting
          display_name <- if(red == "pca") "Unintegrated" else {
            clean <- gsub("integrated\\.", "", red)
            paste0(toupper(substring(clean, 1, 1)), substring(clean, 2))
          }
          # --- PCA embedding (ASW & LISI) ---
          emb <- Embeddings(obj, reduction = red)[, 1:30]
          
          # --- 1. ASW Block ---
          if("ASW" %in% input$eval_metrics) {
            
            dists <- dist(emb)
            if(input$batch_label != "None") {
              sw <- silhouette(as.numeric(factor(obj@meta.data[[input$batch_label]])), dists)[,3]
              detailed_scores <- rbind(detailed_scores, data.frame(Cell=rownames(emb), Score=sw, Integration=display_name, Metric="basw"))
              summary_stats <- rbind(summary_stats, data.frame(Integration=display_name, Method="basw", Median=median(sw)))
            }
            if(input$celltype_label != "None") {
              sw <- silhouette(as.numeric(factor(obj@meta.data[[input$celltype_label]])), dists)[,3]
              detailed_scores <- rbind(detailed_scores, data.frame(Cell=rownames(emb), Score=sw, Integration=display_name, Metric="casw"))
              summary_stats <- rbind(summary_stats, data.frame(Integration=display_name, Method="casw", Median=median(sw)))
            }
          }
          
          # --- 2. LISI Block ---
          if("LISI" %in% input$eval_metrics) {
            
            if(input$batch_label != "None") {
              lisi_res <- compute_lisi(emb, obj@meta.data, input$batch_label)
              detailed_scores <- rbind(detailed_scores, data.frame(Cell=rownames(emb), Score=lisi_res[[input$batch_label]], Integration=display_name, Metric="ilisi"))
              summary_stats <- rbind(summary_stats, data.frame(Integration=display_name, Method="ilisi", Median=median(lisi_res[[input$batch_label]])))
            }
            if(input$celltype_label != "None") {
              lisi_res <- compute_lisi(emb, obj@meta.data, input$celltype_label)
              detailed_scores <- rbind(detailed_scores, data.frame(Cell=rownames(emb), Score=lisi_res[[input$celltype_label]], Integration=display_name, Metric="clisi"))
              summary_stats <- rbind(summary_stats, data.frame(Integration=display_name, Method="clisi", Median=median(lisi_res[[input$celltype_label]])))
            }
          }
          
          # --- 3. GraphLISI Block (High Performance Sparse) ---
          if("GraphLISI" %in% input$eval_metrics) {
            
            target_suffix <- gsub("integrated\\.", "", red)
            # --- NAME MAPPING ---
            # If the reduction is 'pca', the graph is usually just 'snn'
            if(target_suffix == "pca") {
              graph_name <- "snn"
            } else {
              graph_name <- paste0("snn.", target_suffix)
            }
            
            snn_matrix <- obj@graphs[[graph_name]]
            
            # Efficient neighbor extraction
            neighbors <- lapply(1:ncol(snn_matrix), function(i) {
              snn_matrix@i[(snn_matrix@p[i] + 1):snn_matrix@p[i + 1]] + 1
            })
            
            # Helper for local diversity
            get_isi <- function(labels_vec, n_list) {
              sapply(n_list, function(idx) {
                if(length(idx) == 0) return(NA)
                p <- table(labels_vec[idx]) / length(idx)
                1 / sum(p^2)
              })
            }
            
            if(input$batch_label != "None") {
              scores <- get_isi(obj[[input$batch_label, drop=TRUE]], neighbors)
              detailed_scores <- rbind(detailed_scores, data.frame(Cell=rownames(emb), Score=scores, Integration=display_name, Metric="gilisi"))
              summary_stats <- rbind(summary_stats, data.frame(Integration=display_name, Method="gilisi", Median=median(scores, na.rm=TRUE)))
            }
            if(input$celltype_label != "None") {
              scores <- get_isi(obj[[input$celltype_label, drop=TRUE]], neighbors)
              detailed_scores <- rbind(detailed_scores, data.frame(Cell=rownames(emb), Score=scores, Integration=display_name, Metric="gclisi"))
              summary_stats <- rbind(summary_stats, data.frame(Integration=display_name, Method="gclisi", Median=median(scores, na.rm=TRUE)))
            }
          }
        }
      })
      
      # --- RANKING LOGIC ---
      list_metrics <- unique(summary_stats$Method)
      rank_summary <- data.frame()
      
      for(m in list_metrics){
        rank_df <- summary_stats %>% filter(Method == m)
        # ilisi, gilisi, casw: Higher is Better. basw, clisi, gclisi: Lower is Better.
        if(m %in% c("ilisi", "gilisi", "casw")) {
          rank_df <- rank_df %>% arrange(desc(Median))
        } else {
          rank_df <- rank_df %>% arrange(Median)
        }
        rank_df$Rank <- seq(nrow(rank_df), 1) 
        rank_summary <- rbind(rank_summary, rank_df)
      }
      
      return(list(details = detailed_scores, summary = rank_summary))
    })
    
    # --- VISUALIZATIONS ---
    
    output$violin_plots <- renderPlot({
      df <- eval_results()$details
      req(nrow(df) > 0)
      cols_10 <- c("#1F77B4", "#FF7F0E", "#2CA02C","#D62728","#9467BD", "#8C564B", "#E377C2","#7F7F7F", "#BCBD22", "#17BECF")
      
      ggplot(df, aes(x = Integration, y = Score, fill = Integration)) +
        geom_violin(trim = FALSE, show.legend = FALSE) +
        geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, show.legend = FALSE) +
        facet_wrap(~Metric, scales = "free_y", ncol = 3) +
        scale_fill_manual(values = cols_10) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Per-cell Score Distributions", y = "Metric Score")
    })
    
    output$rank_barplot <- renderPlot({
      df <- eval_results()$summary
      req(nrow(df) > 0)
      
      df$Method_group <- ifelse(df$Method %in% c("basw", "ilisi", "gilisi"), 
                                "Batch Mixing", "Cell Type Separation")
      
      cols_6 <- c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F","#EDC948")
      
      ggplot(df, aes(x = Rank, y = Integration, fill = Method)) +
        geom_bar(stat = "identity", color = "white", linewidth = 0.2) +
        facet_grid(. ~ Method_group) +
        scale_fill_manual(values = cols_6) +
        theme_classic() +
        labs(title = "Integration Performance Ranking",
             subtitle = "Higher Rank Score = Better Performance",
             x = "Cumulative Rank Score", y = "Integration Method")
    })
    # --- SCORE TABLE ---
    output$score_table <- renderTable({
      df_details <- eval_results()$details
      req(nrow(df_details) > 0)
      
      df_details %>%
        group_by(Integration, Metric) %>%
        summarize(Median_Score = median(Score, na.rm = TRUE), .groups = 'drop') %>%
        pivot_wider(names_from = Metric, values_from = Median_Score) %>%
        rename(Method = Integration) %>%
        arrange(Method)
    }, digits = 3)
    
   # --- RANK TABLE ###
    output$rank_table <- renderTable({
      # 1. Get the data
      df_rank <- eval_results()$summary
      req(nrow(df_rank) > 0)
      
      # 2. Pivot to Wide format
      # Each integration method gets one row
      rank_pivot <- df_rank %>%
        dplyr::select(Integration, Method, Rank) %>%
        tidyr::pivot_wider(
          names_from = Method, 
          values_from = Rank
        )
      
      # 3. Calculate Total Score safely
      # 'across(where(is.numeric))' automatically ignores the 'Integration' text column
      # and works perfectly whether you have 1 metric or 6 metrics.
      rank_pivot <- rank_pivot %>%
        mutate(Total_Score = rowSums(across(where(is.numeric)), na.rm = TRUE)) %>%
        arrange(desc(Total_Score))
      
      # 4. Final Column Ordering
      # any_of() ensures that if a metric wasn't calculated, the code doesn't break
      rank_pivot <- rank_pivot %>%
        dplyr::select(
          Integration, 
          any_of(c("basw", "ilisi", "gilisi", "casw", "clisi", "gclisi")), 
          Total_Score
        ) %>%
        rename(Method = Integration)
      return(rank_pivot)
    }, digits = 0)
    
  })
}