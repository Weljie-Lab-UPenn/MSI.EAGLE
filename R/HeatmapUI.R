### R/HeatmapUI.R
HeatmapUI <- function(id) {
  ns <- NS(id)
  
  tabPanel(
    "Heatmap and Outlier detection",
    value = ns("tab7"),
    sidebarLayout(
      sidebarPanel(
        uiOutput(ns("hmap_vars"), label = "select_variable_for_heatmap"),
        checkboxGroupInput(
          ns("hmap_params"),
          "Heatmap parameters",
          choices = list(
            "Show only significant peaks" = "hm_sigOnly",
            "Annotate columns" =
              "hm_plotCols",
            "Annotate rows" = "hm_plotRows",
            "Remove legend" = "hmap_legend",
            #"ANOVA "="anova",
            "Scale row" = "rowScale",
            "Scale col" = "colScale",
            "Cluster rows" = "rowclust",
            "Cluster cols" = "colclust",
            "Remove col names" =
              "hm_colNames",
            "Remove row names" =
              "hm_rowNames"
          ),
          selected = c(
            "hm_sigOnly",
            "hm_plotCols",
            "rowclust",
            "colclust",
            "rowScale"
          )
        ),
        uiOutput(ns("hm_labels")),
        numericInput(
          ns("hm_k"),
          label = "kmeans k-value",
          value = NA,
          min = 0
        ),
        numericInput(
          ns("hm_cut_row"),
          label = "cutree row clusters",
          value = NA,
          min = 0
        ),
        numericInput(
          ns("hm_cut_col"),
          label = "cutree col clusters",
          value = NA,
          min = 0
        ),
        selectInput(
          ns("clust_method"),
          label = "clustering method (R hclust)",
          choices = c(
            "ward.D",
            "ward.D2",
            "single",
            "complete",
            "average",
            "mcquitty",
            "median",
            "centroid"
          ),
          selected = "complete",
          multiple = F
        ),
        selectInput(
          ns("clustering_distance_columns"),
          label = "column distance method",
          choices = c(
            "euclidean",
            "maximum",
            "manhattan",
            "canberra",
            "binary",
            "minkowski"
          )
        ),
        selectInput(
          ns("clustering_distance_rows"),
          label = "row distance method",
          choices = c(
            "euclidean",
            "maximum",
            "manhattan",
            "canberra",
            "binary",
            "minkowski"
          )
        ),
        selectInput(
          ns("heatmap_colors"),
          label = "Heatmap colors",
          choices = c(
            "default",
            "viridis","magma", "plasma", "inferno", "cividis" 
          )
          
        ),
        
        uiOutput(ns('hm_sig_select'), label = "Select statistics result for filtering"),
        actionButton(ns("action_hmap"), "Generate heatmap"),
        downloadButton(ns('Export'))
      ),
      
      
      mainPanel(tabsetPanel(
        tabPanel("Heatmap",
                 imageOutput(
                   ns("plot14"), width = "800px", height = "600px"
                 )),
        tabPanel(
          "Outlier detection",
          p('Outlier detection'),
          numericInput(ns('hdo_alpha'), 'alpha value', value =
                         0.01),
          numericInput(ns('hdo_k'), 'k value', value = 10),
          selectInput(
            ns('knnsearchtype'),
            'KNN searchtype',
            choices = c("brute", "kd_tree"),
            selected = "brute"
          ),
          numericInput(
            ns('hdo_p'),
            "Proportion of possible outliers",
            min = 0,
            max = 0.5,
            value = 0.25
          ),
          numericInput(
            ns('tn'),
            "Sample size to calculate emperical threshold",
            value = 50
          ),
          actionButton(ns("action_outliers"), "Find outliers"),
          actionButton(
            ns("action_outliers_heatmap"),
            "Generate heatmap with outliers removed"
          ),
          actionButton(
            ns("action_outliers_save"),
            "Save dataset with outliers removed"
          )
          
        )
        
      ))
      
      
    )
    
  )
  
  
}