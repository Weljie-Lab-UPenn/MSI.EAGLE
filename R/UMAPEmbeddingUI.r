### UMAPEmbeddingUI.R
UMAPEmbeddingUI <- function(id) {
  ns <- NS(id)
  
  tabPanel("UMAP Embeddings",
           value = ns("tab3b"),
           sidebarLayout(
             sidebarPanel(
               # Visualization controls
               uiOutput(ns("mz_selector")),
               uiOutput(ns("color_selector")),
               
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
               )
             ),
             
             mainPanel(
               verbatimTextOutput(ns("variables")),
               # Run selection table
               DT::dataTableOutput(ns("mytable")),
               
               # Visualization outputs
               imageOutput(ns("umapPlot"), 
                           width = "800px",
                           height = "600px"),
               imageOutput(ns("umapPlot2"),
                           width = "800px", 
                           height = "600px"),
               imageOutput(ns("phenoplot"),
                           width = "800px",
                           height = "600px")
             )
           ))
}
