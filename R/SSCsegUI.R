### R/SSCsegUI.R
SSCsegUI <- function(id) {
  ns <- NS(id)
  
  tabPanel("Segmentation - SSC",
           value = ns("tab3a"),
           sidebarLayout(
             sidebarPanel(
               checkboxInput(ns('fix_pix_tf'), 
                             "Show fix stray pixel options?",
                             FALSE),
               
               uiOutput(ns('fix_pix_opts')),
               textInput(ns("card_r"), "radius for SSC", value = 2),
               textInput(ns("card_k"), "k (clusters) for SSC", value =
                           5),
               textInput(ns("card_s"), "shrinkage for Cardinal SSC", value =
                           "2^(1:6)"),
               
               selectInput(
                 ns("ssc_method"),
                 "Method for SSC",
                 choices = list("adaptive" = "adaptive", "gaussian" =
                                  "gaussian")
               ),
               
               actionButton(ns("action_ssc"), label = "Start SSC"),
               uiOutput(ns('ssc_model_choices')),
               uiOutput(ns('Color_choices2')),
               selectInput(
                 ns("ssc_cols"),
                 "Colors for ssc visualzation",
                 choices = list("alphabet" = "alphabet", "default" =
                                  NULL)
               ),
               
               actionButton(ns("store_proc2"), label = "Store processed data"),
               shinyFiles::shinySaveButton(ns("save_imzml"), "Save imzML File", "Save", filetype = list(""))
               #downloadButton(ns("save_proc2"), label = "Save processed data")
               
               
               
             ),
             mainPanel(
               DT::dataTableOutput(ns("mytable2")),
               plot_card_UI(ns("card_plot_ssc")),
               imageOutput(ns("plot6a"), width = "800px", height = "600px"),
               imageOutput(ns("plot5a"), width = "800px", height = "600px")
             )
             
           ))
}
