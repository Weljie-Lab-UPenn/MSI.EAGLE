
### R/PeakPickUI.R
PeakPickUI <- function(id) {
  ns <- NS(id)
  
  
  tabPanel(HTML("File Restore and <br/> Overview Analysis"),
           value=ns("tab2"),
           sidebarLayout(
             sidebarPanel(
               
               radioButtons(ns("pp_operation"), "Peak-picking operation",
                            c("Open existing peak-picked dataset"="open_file",
                              "Peak-pick raw files from Data Setup"="pp_raw",
                              "Combine datasets (same coordinates; merge features)"="subset_f",
                              "Combine datasets (same peak list; append runs)"="add_same_pklist"
                              #"Demo data from CardinalWorkflows package"="demo"
                              )
               ),
               uiOutput(ns("operation_help")),
               uiOutput(ns("peak_pick")),
               fluidRow(column(2, verbatimTextOutput(ns("value"))))

               
             ),
             mainPanel(
               verbatimTextOutput(ns("variables")),
               DT::dataTableOutput(ns("peak_pick_selection")),
               textOutput(ns("overview_text")),
               plot_card_UI(ns("card_plot"))
             )
           )
  )
}
