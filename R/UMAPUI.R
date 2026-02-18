UMAPUI <- function(id) {
  ns <- NS(id)
  tabPanel(
    "Segmentation - UMAP",
    value = ns("tab3"),
    tags$head(
      tags$style(HTML("
        .umap-sidebar .umap-section {
          margin-bottom: 10px;
          padding: 8px;
          border: 1px solid #ddd;
          border-radius: 4px;
          background: #fafafa;
        }
        .umap-sidebar .umap-section h5 {
          margin-top: 0;
          margin-bottom: 8px;
          font-weight: 600;
        }
        .umap-sidebar .form-group {
          margin-bottom: 6px;
        }
        .umap-sidebar .control-label {
          margin-bottom: 2px;
          font-size: 12px;
        }
        .umap-sidebar .checkbox, .umap-sidebar .radio {
          margin-top: 2px;
          margin-bottom: 2px;
        }
        .umap-sidebar .btn {
          margin-top: 4px;
          margin-right: 6px;
        }
      "))
    ),
    sidebarLayout(
      sidebarPanel(
        class = "umap-sidebar",
        width = 4,
        tags$div(
          class = "umap-section",
          tags$h5("Mode"),
          radioButtons(
            ns("seg_choice"),
            "Segmentation Goal",
            c(
              "Background removal (pixelData)" = "bk_seg",
              "pData/UMAP Editor" = "anat_seg",
              "Fix stray pixels" = "fix_pix"
            )
          ),
          uiOutput(ns("fix_pix_opts"))
        ),
        tags$div(
          class = "umap-section",
          tags$h5("UMAP Parameters"),
          uiOutput(ns("umap_choice"))
        ),
        tags$div(
          class = "umap-section",
          tags$h5("Clustering Parameters"),
          fluidRow(
            column(4, numericInput(ns("k_clustering"), "k", 5, min = 2)),
            column(4, numericInput(ns("eps"), "eps", 0.15, min = 0.0001)),
            column(4, numericInput(ns("minPts"), "minPts", 50L, min = 2))
          ),
          checkboxGroupInput(
            ns("clustering_methods"),
            "Methods",
            choices = list(
              "K-means" = "kmeans",
              "Hierarchical" = "hierarchical",
              "DBSCAN" = "dbscan",
              "HDBSCAN" = "hdbscan",
              "Spectral" = "spectral",
              "K-medoids" = "kmedoids",
              "Fuzzy C-means" = "fuzzy",
              "Mclust" = "mclust",
              "SOM" = "som",
              "Spherical K-means" = "skmeans"
            ),
            selected = c("kmeans")
          ),
          actionButton(ns("action_dbscan"), label = "Re-run color clustering")
        ),
        tags$div(
          class = "umap-section",
          tags$h5("Display and Output"),
          fluidRow(
            column(6, numericInput(ns("cex"), "Tile size (cex)", value = 2, min = 0.1)),
            column(6, numericInput(ns("quant"), "Reduction quantile", value = 0.9, min = 0, max = 0.99))
          ),
          selectInput(
            ns("umap_cols"),
            "Colors for UMAP visualization",
            choices = list(
              "R reduced" = "Reduced2",
              "R colors" = "Reduced",
              "K-means" = "kmeans",
              "Hierarchical" = "hierarchical",
              "DBSCAN" = "dbscan",
              "HDBSCAN" = "hdbscan",
              "Spectral Clustering" = "spectral",
              "K-medoids" = "kmedoids",
              "Fuzzy C-means" = "fuzzy",
              "Model-based Clustering (Mclust)" = "mclust",
              "Self-Organizing Map (SOM)" = "som",
              "Spherical K-means" = "skmeans"
            )
          ),
          uiOutput(ns("Color_choices")),
          actionButton(ns("store_proc"), label = "Store processed data"),
          shinyFiles::shinySaveButton(ns("save_imzml"), "Save imzML File", "Save", filetype = list(""))
        )
      ),
      mainPanel(
        verbatimTextOutput(ns("variables")),
        DT::dataTableOutput(ns("mytable")),
        plot_card_UI(ns("plot_card_umap")),
        textOutput(ns("data_list_table")),
        textOutput(ns("summary")),
        uiOutput(ns("segmentation_plots"))
      )
    )
  )
}
