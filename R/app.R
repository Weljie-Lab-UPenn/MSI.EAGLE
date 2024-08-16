#' @import shiny
#' @import BiocParallel
#' @import Cardinal
#' @import markdown
#' @import shinydashboard
#' @import shinythemes
#' @import uwot

# need to run devtools::document() to ensure imports are added to NAMESPACE



#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[DT]{dataTableOutput}}
#' @rdname MSI.EAGLE
#' @export 
# @importFrom DT dataTableOutput renderDataTable
MSI.EAGLE <- function(...) {
  
  
  #check for directories and set to current wd if not set
  if(!exists("rawd")){
    rawd=getwd()
  }
  
  if(!exists("wd")){
    wd=getwd()
  }
  
  #check for number of cores and set to #detected -2 by default
  if(!exists("ncores")){
    ncores=as.integer(parallel::detectCores())-2
  }
  
  ui <- navbarPage(p(strong("MSI EAGLE")),
                   theme = shinythemes::shinytheme("sandstone"),                   
                   DataSetupUI("tab1", ncores, rawd, wd),
                   PeakPickUI("tab2"),
                   UMAPUI("tab3"),
                   SSCsegUI("tab3a"),
                   PhenoSetupUI("tab4"),
                   DepthAnalysisUI("tab5"),
                   StatsPrepUI("tab6"),
                   HeatmapUI("tab7"),
                   CorrelationUI("tab8"),
                   tabPanel("Help",
                            HTML("<p><b>MSI.EAGLE</b></p>
                                 <p>Please look at the github page for the latest manual,  updates, and contact info.</p>
                                 <a href='http://165.123.67.19:30003/aalim/DESI_Shiny_Processing_script/src/branch/modules', target='_blank'>Click here</a>")
                            ),
                   tabPanel("debug",
                            actionButton("browser", "debug browser()")
                   )
                   
                   
  )
  
  
  server <- function(input, output, session){
    
    setup_values <- DataSetupServer("tab1", rawd = rawd, wd = wd)
    preproc_values <- PeakPickServer("tab2", setup_values)
    UMAPServer("tab3", setup_values, preproc_values)
    SSCsegServer("tab3a", setup_values, preproc_values)
    PhenoServer("tab4", setup_values, preproc_values)
    DepthAnalysisServer("tab5", setup_values)
    proc_values<-StatsPrepServer("tab6", setup_values )
    HeatmapServer("tab7", proc_values)
    CorrelationServer("tab8", proc_values)
    
    observeEvent(input$browser, {
      browser()
    })
    
  }
  
  shinyApp(ui=ui, server=server)
}
