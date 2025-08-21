

### R/DataSetupUI.R
DataSetupUI <- function(id, ncores, nchunks, rawd, wd) {
  ns <- NS(id)
  
  tabPanel(
    title = "Data Setup",
    value = ns('tab1'),
    sidebarLayout(
      sidebarPanel(
        #p("option to remove pixels from raw data processing"),
        # Copy the line below to make a text input box
        textInput(
          ns("folder"),
          label = p("Raw files directory (ie. imzML file location)"),
          value = rawd
        ),
        
        
        textInput(
          ns("wd"),
          label = p("Working directory (ie. where results are stored)"),
          value = wd
        ),
        
        shinyFiles::shinyDirButton(
          id = ns('wd_new'),
          label = 'Select directory',
          'Please select a new working direcotry',
          FALSE
        ),
        
        textInput(
          ns("regex"),
          label = p("Regex syntax for raw files."),
          value = "(.*)imzML$"
        ),
        uiOutput(ns("par_mode_setup")),
        fluidRow(
          column(6, numericInput(ns("ncores"), "# of cores", ncores)),
          column(6, numericInput(ns("chunks"), "# of chunks", nchunks))
        ),
        # numericInput(ns("ncores"), "# of cores e", ncores),
        # # radioButtons(
        # #   ns("mode"),
        # #   "Ion mode (only for visualization)",
        # #   c("Positive" = "p", "Negative" = "n")
        # # ),
        selectInput(
          ns("pp_params"),
          label = "Peak picking params",
          choices = list("qTof1" = "qtof1", "Hi Res" =
                           "hires"),
          selected = "qtof1"
        ),
        uiOutput(ns("pp_params_display")),
        actionButton(ns("action1"), label = "Import Data"),
        
        
        wellPanel(
          tags$hr(style = "border-top: 2px solid #FF5733;"),  # Change color and thickness
          
          p("Quick Look for Quality Control:"),
          
          
          p("Limits for plotting spectrum to ensure MS looks reasonable"),
        
          numericInput(ns("pix_to_plot"), "% of pixels to randomly sample for plot", 1),
          actionButton(
            ns("action2"),
            label = HTML("Extract sample and <br/> MS Spectrum plot")
          )
        ),
        wellPanel(
          p(),
          tags$hr(style = "border-top: 2px solid #FF5733;"),  # Change color and thickness
          selectInput(
            ns("cardworkdat"),
            label = "Choose sample dataset from CardinalWorkflows",
            choices =
              c(
                "Human Renal Cell Carcinoma (RCC) with background",
                "Human Renal Cell Carcinoma (RCC) (no background)",
                "Whole Pig Fetus Cross-Section"
              )
          ),
          actionButton(ns("action_demo"), label = HTML("Load demo data"))
        )
        
        
        
      ),
      mainPanel(
        p("Setup data directories, some options, and import .imzML files"),
        p(""),
        verbatimTextOutput(ns("variables")),
        DT::dataTableOutput(ns("files")),
        plot_card_UI(ns("card_plot"))
        
      )
    )
  )
  
}
