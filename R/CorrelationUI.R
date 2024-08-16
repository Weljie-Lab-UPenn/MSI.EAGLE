### R/CorrelationUI.R
CorrelationUI <- function(id) {
  ns <- NS(id)
  
  tabPanel("Colocalization analysis",
           value="tab8",
           sidebarLayout(
             sidebarPanel(
               uiOutput(ns("corr_vars"), label="select_variable_for_correlation"),
               uiOutput(ns("corr_mz"), label="m/z selection for correlation"),
                 numericInput(ns("corr_n"), "Number of colocalized features?", min=2, value = 10),
               radioButtons(ns("corr_sort"), "Sorted by", choices = list(
                 "correlation"="correlation",
                 "Manders’ colocalization coefficients (M1)"="M1",
                 "Manders’ colocalization coefficients (M2)"="M2"
               ), selected = "correlation"),
               numericInput(ns("corr_plot_n"), "Number of top correlations to plot", min=1, value=3),
               actionButton(ns("action_corr"), "Generate colocalization")
               
               
             ),
             mainPanel(
               imageOutput(ns("plot15"), width="800px", height="600px")
               
             )
             
           )
           
  )
  
  
}