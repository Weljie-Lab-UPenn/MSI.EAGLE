

### R/UMAPUI.R
UMAPUI <- function(id) {
  ns <- NS(id)
  tabPanel("Segmentation - UMAP",
           value = ns("tab3"),
           sidebarLayout(
             sidebarPanel(
               radioButtons(
                 ns('seg_choice'),
                 "Segmentation Goal",
                 c(
                   "Background removal (pixelData)" = 'bk_seg',
                   "Anatomical segmentation (phenotypeData)" =
                     'anat_seg',
                   "Fix stray pixels" = "fix_pix"
                 )
               ),
               uiOutput(ns('fix_pix_opts')),
               uiOutput(ns('umap_choice')),
               
               selectInput(
                 ns("umap_cols"),
                 "Colors for UMAP visualzation",
                 choices = list(
                   "R reduced" = "Reduced2",
                   "R colors" = "Reduced",
                   DBSCAN = "dbscan",
                   HDBSCAN = "hdbscan"
                 )
               ),
               uiOutput(ns('Color_choices')),
               
               
               actionButton(ns("store_proc"), label = "Store processed data"),
               #downloadButton(ns("save_proc"), label = "Save processed data")
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
             
           ))
}
