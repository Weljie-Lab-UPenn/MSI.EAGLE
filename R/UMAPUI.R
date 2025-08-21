UMAPUI <- function(id) {
  ns <- NS(id)
  tabPanel(
    "Segmentation - UMAP",
    value = ns("tab3"),
    sidebarLayout(
      sidebarPanel(
        radioButtons(
          ns('seg_choice'),
          "Segmentation Goal",
          c(
            "Background removal (pixelData)" = 'bk_seg',
            "pData/UMAP Editor" = 'anat_seg',
            "Fix stray pixels" = "fix_pix"
          )
        ),
        uiOutput(ns('fix_pix_opts')),
        uiOutput(ns('umap_choice')),
        
        # **New Checkbox Group Input for Clustering Methods**
        checkboxGroupInput(
          ns("clustering_methods"),
          "Select Clustering Methods:",
          choices = list(
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
          ),
          selected = c("kmeans")  # Default selected methods
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
        uiOutput(ns('Color_choices')),
        
        actionButton(ns("store_proc"), label = "Store processed data"),
        # downloadButton(ns("save_proc"), label = "Save processed data")
        shinyFiles::shinySaveButton(ns("save_imzml"), "Save imzML File", "Save", filetype = list(""))
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
