#' @import shiny
#' @import BiocParallel
#' @import Cardinal
#' @import markdown
#' @import shinydashboard
#' @import shinythemes
#' @import uwot
#' @import magrittr

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
  library(magrittr)
  # Allow larger uploads (histology images / polygon files) than Shiny default 5 MB.
  options(shiny.maxRequestSize = 500 * 1024^2)
  
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
  
  #check of number of chunks and set if present
  if(!exists("nchunks")){
    nchunks=20
  }
  
  ui <- navbarPage(p(strong("MSI EAGLE")),
                   theme = shinythemes::shinytheme("sandstone"),                   
                   DataSetupUI("tab1", ncores, nchunks, rawd, wd),
                   PeakPickUI("tab2"),
                   UMAPUI("tab3"),
                   SSCsegUI("tab3a"),
                   UMAPEmbeddingUI("tab3b"),
                   PhenoSetupUI("tab4"),
                   MaskedAnalysisUI("tab5"),
                   StatsPrepUI("tab6"),
                   HeatmapUI("tab7"),
                   CorrelationUI("tab8"),
                   HistologyIntegrationUI("tab9"),
                   tabPanel("Help",
                            HTML("<p><b>MSI.EAGLE</b> is a general MSI analysis platform with optional single-cell and multimodal registration workflows.</p>
                                 <ul>
                                   <li>Routine MSI preprocessing, visualization, clustering, and statistics</li>
                                   <li>Polygon-informed mapping of cell regions into MSI pData</li>
                                   <li>Spatially aware UMAP and clustering workflows</li>
                                   <li>Histology, cluster-image, and polygon registration to MSI space</li>
                                   <li>Integrated statistical testing and reproducible exports</li>
                                 </ul>
                                 <p>For setup and updates, see the project repository:</p>
                                 <a href='https://github.com/Weljie-Lab-UPenn/MSI.EAGLE' target='_blank'>https://github.com/Weljie-Lab-UPenn/MSI.EAGLE</a>")
                            ),
                   tabPanel("debug",
                            actionButton("browser", "debug browser()")
                   )
                   
                   
  )
  
  
  server <- function(input, output, session){
    
    setup_values <- DataSetupServer("tab1", rawd = rawd, wd = wd)
    preproc_values <- PeakPickServer("tab2", setup_values)
    UMAPServer("tab3", setup_values, preproc_values, NULL)
    SSCsegServer("tab3a", setup_values, preproc_values)
    UMAPEmbeddingServer("tab3b", setup_values, preproc_values, NULL)
    PhenoServer("tab4", setup_values, preproc_values)
    MaskedAnalysisServer("tab5", setup_values)
    proc_values<-StatsPrepServer("tab6", setup_values )
    HeatmapServer("tab7", proc_values)
    CorrelationServer("tab8", proc_values, setup_values)
    HistologyIntegrationServer("tab9", setup_values, preproc_values)
    
    observeEvent(input$browser, {
      browser()
    })
    
  }
  
  shinyApp(ui=ui, server=server)
}
