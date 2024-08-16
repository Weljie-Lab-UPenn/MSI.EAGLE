
### R/PhenoSetupUI.R
PhenoSetupUI <- function(id) {
  ns <- NS(id)
  
  
  tabPanel("Phenotyping",
           value = ns("tab4"),
           sidebarLayout(
             sidebarPanel(
               fileInput(ns('phenotype_file'), "File with Phenotype info (.txt)", accept = "text/plain"),
               radioButtons(ns("img_phen"), "MSI image to phenotype",
                            choices = list(  "Restored from file in 'File restore..' tab"="phen_pp", "Stored data (less common)"="phen_no_bk")),
               actionButton(ns("read_phenotype"), "Read dataset to be phenotyped"),
               radioButtons(ns("phen_method"), "Phenotype method",
                            choices = list("Spectral density"="spec_density", "Periodicity"="period", "Breaks between samples"="breaks", "Manual (x & y limits specified in file)"="manual")),
               uiOutput(ns("phen_method_opts")),
               #radioButtons("phen_type", "Sample to phenotype",
               #             choices = list( "Tissue"="tissue", "Spots"="spots")),
               uiOutput(ns('phen_cols')),
               actionButton(ns("start_phenotype"), "Start phenotyping"),
               checkboxInput(ns("key_tf"), "Hide key / lenged?", value=FALSE),
               
               uiOutput(ns('plot_cols')),
               uiOutput(ns("phen_interaction")),
               actionButton(ns("interaction"), label="Add interaction data"),
               shinyFiles::shinySaveButton(ns("save_imzml"), "Save imzML File", "Save", filetype = list(""))
               
             ),
             mainPanel(
               imageOutput(ns("plot8"), width="800px", height="600px"),
               
               DT::dataTableOutput(ns("phenotype_table_new")),
               DT::dataTableOutput(ns("phenotype_table"))
               #textOutput("data_list_table"),
               #textOutput("summary"),
               #imageOutput("colormaps", width="800px", height="300px")
             )
             
           )
           
  )
  
}
