
### R/DepthAnalysisUI.R
DepthAnalysisUI <- function(id) {
  ns <- NS(id)
  
  
  tabPanel("Depth Analysis",
           value = ns("tab5"),
           sidebarLayout(
             sidebarPanel(

               
               uiOutput(ns("segmented_data"), label="Choose segmented file as template"),
               actionButton(ns("action_seg"), label = HTML("Read file with coordinates")),
               radioButtons(ns("targeted_pp"), label="Targeted or untargeted analysis", choices = list("Untargeted"="untargeted", "Targeted"="targeted", "Mean spectrum"="mean")),
               uiOutput(ns("pp_options")),
               actionButton(ns("action_seg_run"), label = HTML("Peak pick or bin segmented data")),
               #downloadButton(ns("save_state2"), "Save binned file")
               shinyFiles::shinySaveButton(ns("save_imzml"), "Save imzML and .rds Files", "Save", filetype = list("")),
               #downloadButton(ns("save_state2"), "Save .rds file")
               
             ),
             mainPanel(
               verbatimTextOutput(ns("filename_header")),
               plot_card_UI(ns("card_plot_depth"))
             )
           )
  )
  
}
