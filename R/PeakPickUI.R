
### R/PeakPickUI.R
PeakPickUI <- function(id) {
  ns <- NS(id)
  
  
  tabPanel(HTML("File Restore and <br/> Overview Analysis"),
           value=ns("tab2"),
           sidebarLayout(
             sidebarPanel(
               
               radioButtons(ns("pp_operation"), "Peakpeaking operation",
                            c("Open previously peakpicked file"="open_file",
                              "Peakpick rawfiles (ie. imzML loaded in Data Setup tab)"="pp_raw",
                              "Add two imagesets (same coordinates)"="subset_f",
                              "Add two imagesets (same peaklist)"="add_same_pklist"
                              #"Demo data from CardinalWorkflows package"="demo"
                              )
               ),
               uiOutput(ns("peak_pick")),
               fluidRow(column(2, verbatimTextOutput(ns("value"))))

               
             ),
             mainPanel(
               DT::dataTableOutput(ns("peak_pick_selection")),
               textOutput(ns("overview_text")),
               plot_card_UI(ns("card_plot"))
             )
           )
  )
}
