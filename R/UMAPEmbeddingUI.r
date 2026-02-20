### R/UMAPEmbeddingUI.R

UMAPEmbeddingUI <- function(id) {
  ns <- NS(id)
  
  tabPanel("UMAP Embeddings",
           value = ns("tab3b"),
           sidebarLayout(
             sidebarPanel(
               uiOutput(ns("mz_selector")),
               uiOutput(ns("color_selector")),
               
               # Add color choices UI element
               uiOutput(ns("color_choices")),
               
               sliderInput(
                 ns("point_size"), 
                 "Point Size:", 
                 min = 0.5, 
                 max = 3, 
                 value = 1.5
               ),
               
               sliderInput(
                 ns("color_midpoint"),
                 "Color Scale Midpoint:",
                 min = 0,
                 max = 1,
                 value = 0.5,
                 step = 0.01
               ),
               
               checkboxInput(
                 ns("log_scale"), 
                 "Log scale data", 
                 value = FALSE
               ),
               
               checkboxInput(
                 ns("col_match"),
                 "Match color to original UMAP colors",
                 value = FALSE
               ),
               
               checkboxInput(
                 ns("plot_mz"),
                 "Spatial plot of m/z value?", 
                 value = FALSE
               ),
               tags$hr(),
               tags$div(style = "font-weight:600; margin-bottom:6px;", "Re-run clustering on existing UMAP"),
               checkboxGroupInput(
                 ns("embed_clustering_methods"),
                 "Methods",
                 choices = c(
                   "K-means" = "kmeans",
                   "Hierarchical" = "hierarchical",
                   "DBSCAN" = "dbscan",
                   "HDBSCAN" = "hdbscan",
                   "Spectral" = "spectral",
                   "K-medoids (slow)" = "kmedoids",
                   "Fuzzy C-means" = "fuzzy",
                   "Mclust" = "mclust",
                   "SOM" = "som",
                   "Spherical K-means" = "skmeans"
                 ),
                 selected = c("kmeans", "hierarchical")
               ),
               fluidRow(
                 column(
                   6,
                   numericInput(ns("embed_k_clusters"), "k clusters", value = 12, min = 2)
                 ),
                 column(
                   6,
                   numericInput(ns("embed_minPts"), "minPts", value = 10, min = 2)
                 )
               ),
               numericInput(ns("embed_eps"), "DBSCAN eps", value = 0.5, min = 0.001, step = 0.01),
               checkboxInput(ns("embed_allow_highmem"), "Allow high-memory methods without guard", value = FALSE),
               fluidRow(
                 column(6, numericInput(ns("embed_pairwise_limit"), "Pairwise guard n", value = 12000, min = 2000, step = 500)),
                 column(6, numericInput(ns("embed_spectral_limit"), "Spectral guard n", value = 8000, min = 1000, step = 500))
               ),
               checkboxInput(ns("embed_write_diag"), "Write diagnostics log to working directory", value = FALSE),
               actionButton(ns("rerun_embed_clustering"), "Re-run clustering on selected runs"),
               helpText("Writes clustering outputs to pData columns (e.g., col_kmeans, col_dbscan)."),
               tags$hr(),
               actionButton(ns("embed_store_proc"), "Store embedding clustering in processed data"),
               shinyFiles::shinySaveButton(ns("embed_save_imzml"), "Save imzML (from embeddings)", "Save", filetype = list(""))
             ),
             
             mainPanel(
               DT::dataTableOutput(ns("mytable")),
               imageOutput(ns("umapPlot"), 
                           width = "800px",
                           height = "600px"),
               imageOutput(ns("umapPlot2"),
                           width = "800px", 
                           height = "600px"), 
               plotOutput(ns("phenoplot"),
                          width = "800px",
                          height = "600px")
             )
           ))
}
