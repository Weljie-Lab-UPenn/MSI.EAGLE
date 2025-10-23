### R/StatsPrepUI.R

# function to create the UI for the statistics tab
StatsPrepUI <- function(id) {
  # create a namespace for the UI
  ns <- NS(id)
  
  
  # create the tabPanel for the statistics tab
  tabPanel("Statistics",
           
           # set the value of the tabPanel to the namespace
           value = ns("tab6"),
           
           # create a sidebarLayout for the tabPanel
           sidebarLayout(
             # create a sidebarPanel for the sidebarLayout
             sidebarPanel(
               # create an uiOutput for the sidebarPanel
               uiOutput(ns("final_data"), label = "Choose target file for analysis"),
               
               # create an actionButton for the sidebarPanel
               actionButton(ns("action_read_stats"), "Read file for stats"),
               
               # create a radioButtons for the sidebarPanel
               radioButtons(
                 ns("stats_test"),
                 "Type of test",
                 choices = list(
                   "Means test (Cardinal)" = "meanstest",
                   "Spatial shrunken centroids (Cardinal)" =
                     "ssctest",
                   "Spatial DGMM (Cardinal)" =
                     "spatialDGMM"
                 )
               ),
               
               #uiOutput('grouping_variables'),
               
               numericInput(
                 ns("FDR_val"),
                 label = "FDR threshold to include in results",
                 value = 0.05,
                 max = 1.0,
                 min = 0.0,
                 step=0.05
               ),
               numericInput(
                 ns("debug_n"),
                 "# Peaks (debugging); set before reading .rds/.imzML file",
                 value = NULL
               ),
               #downloadButton("save_results", label="Save analysis data")
               actionButton(
                 ns("action_demo"), 
                 label="load demo rcc data (CardinalWorkflows)"
               )
             ),
             
             # create a mainPanel for the sidebarLayout
             mainPanel(# create a tabsetPanel for the mainPanel
               # create a verbatimTextOutput for the tabPanel
               verbatimTextOutput(ns("filename_header")),
               
               tabsetPanel(
                 # create a tabPanel for the tabsetPanel
                 tabPanel(
                   "Setup and Run Statistics",
                   
                   # set the value of the tabPanel to the namespace
                   value = ns("tab6a"),
                   
                   
                   
                   # create an uiOutput for the tabPanel
                   uiOutput(ns('phen_cols_stats')),
                   
                   # create an uiOutput for the tabPanel
                   uiOutput(ns("ui"))
                 ),
                 # create a tabPanel for the tabsetPanel
                 tabPanel(
                   "Output Table and Ion Plots",
                   
                   # set the value of the tabPanel to the namespace
                   value = ns("tab6c"),
                   
                   # #ask user for choice of cardinal plot or ggplot
                   # radioButtons(
                   #   ns("plot_choice"),
                   #   "Choose plotting package",
                   #   choices = list(
                   #     "Cardinal / matter" = "cardinal",
                   #     "ggplot2 " = "ggplot"
                   #   )
                   # ),
                   # user check box for which plot to output dynamically
                   uiOutput(ns("plot_choice_ui")),
                   
                   # create an imageOutput for the tabPanel
                   imageOutput(ns("plot11"), width = "800px", height =
                                 "600px"),
                   
                   # create a dataTableOutput for the tabPanel
                   DT::dataTableOutput(ns("stats_table")),
                   
                   wellPanel(
                     # create a textInput for the tabPanel
                     textInput(ns("plot_prefix"), label = "Prefix for plotting / results output", value = "output"),
                     
                     # create an actionButton for the tabPanel
                     actionButton(ns("write_plots"), label = "save significant plots from testing to directory"),
                     
                     # ADD THESE NEW CONTROLS BELOW:
                     tags$hr(),
                     tags$h5("MSI Image Export Settings:"),
                     fluidRow(
                       column(6, selectInput(ns("export_contrast"), "Contrast enhancement", 
                                             c("none", "histogram", "adaptive"), selected = "histogram")),
                       column(6, selectInput(ns("export_colorscale"), "Colorscale", 
                                             c("Spectral", "Cividis", "Viridis", "Inferno", "Plasma", 
                                               "Zissou 1", "Purple-Green", "Berlin", "PiYG", "Grays", 
                                               "Batlow", "turku", "YlOrRd", "Terrain", 
                                               "PrGn", "Green-Brown", "Hawaii", "Cork", "Rocket", "RdYlBu"),
                                             selected = "Spectral"))
                     ),
                     fluidRow(
                       column(6, selectInput(ns("export_smooth"), "Smoothing options", 
                                             c("none", "mean", "gaussian", "bilateral", "adaptive", "diffusion", "guided"),
                                             selected = "none")),
                       column(6, checkboxInput(ns("export_scale"), "Scale multiple images?", value = TRUE))
                     ),
                     checkboxInput(ns("export_dark_bg"), "Dark background?", value = FALSE)
                   ),
                 ),
                 
                 # create a tabPanel for the tabsetPanel
                 tabPanel(
                   "Overview Plots",
                   
                   # set the value of the tabPanel to the namespace
                   value = ns("tab6b"),
                   
                   # create a dataTableOutput for the tabPanel
                   DT::dataTableOutput(ns("mytable_stats_plate")),
                   
                   # create a plot_card_UI for the tabPanel
                   plot_card_UI(ns("plot_card_stats")),
                   
                   # create an imageOutput for the tabPanel
                   #imageOutput(ns("plot12"), width="800px", height="600px"),
                   
                   # create an imageOutput for the tabPanel
                   imageOutput(ns("plot13"), width = "800px", height =
                                 "600px")
                   
                   
                   
                 ),
                 
                 
                 # create a tabPanel for the tabsetPanel
                 tabPanel(
                   "Data Export",
                   
                   # set the value of the tabPanel to the namespace
                   value = ns("tab6d"),
                   
                   # create an uiOutput for the tabPanel
                   #uiOutput(ns("output_factors")),
                   
                   # create an uiOutput for the tabPanel
                   uiOutput(ns("grouping_variables_export")),
                   
                   # create an actionButton for the tabPanel
                   actionButton(
                     ns("data_export"),
                     label = HTML("Export data table with<br/>selected variables")
                   )
                 )
                 
                 
                 
                 
                 
               ))
             
           ))
  
}
