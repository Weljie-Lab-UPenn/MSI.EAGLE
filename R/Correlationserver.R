### R/CorrelationServer.R
CorrelationServer <- function(id,  proc_values) {
  moduleServer(id, function(input, output, session) {
    ns = session$ns #for dyanamic variable namespace
    x5 =  observe({
      proc_values()[["x5"]]
    }) # if more than one value in list; not sure how to do it otherwise
    #x5<- reactive({proc_values()$x5})
    
    
    output$corr_vars <-
      renderUI({
        #https://gist.github.com/wch/4211337
        # If missing input, return to avoid error later in function
        req(proc_values()$x5$stats_results)
        
        
        hmap_choices <-
          colnames(proc_values()$x5$stats_results)[!colnames(proc_values()$x5$stats_results) %in%
                                                     c("mz", "feature")]
        
        list(
          selectInput(
            ns("corr_sig_select"),
            "Choose filtering column for correlation choices",
            choices  = hmap_choices,
            selected = "AdjP"
          ),
          radioButtons(
            ns("corr_sig_direction"),
            label = "Filtering stat direction?",
            choices = list("Ascending" = "ascending",
                           "Descending" = "descending"),
            selected = "descending"
          ),
          numericInput(ns("corr_sig"), "Significance cutoff", value = .05),
          checkboxInput(ns("key_on"), "Colorkey on?", value=TRUE)
        )
      })
    
    output$corr_mz <-
      renderUI({
        #https://gist.github.com/wch/4211337
        # If missing input, return to avoid error later in function
        req(input$corr_sig_select)
        
        
        
        if (input$corr_sig_direction == "descending") {
          sig_mz <-
            proc_values()$x5$stats_results %>% dplyr::filter(.data[[input$corr_sig_select]] <= input$corr_sig) %>% dplyr::select(mz)
        } else {
          sig_mz <-
            proc_values()$x5$stats_results %>% dplyr::filter(.data[[input$corr_sig_select]] >= input$corr_sig) %>% dplyr::select(mz)
        }
        
        
        
        corr_choices <- sig_mz
        
        selectInput(ns("corr_mz"),
                    "Choose m/z value for correlation analysis",
                    choices = corr_choices)
      })
    
    
    
    observeEvent(input$action_corr, {
      output$plot15 <- renderImage({
        req(input$corr_mz)
        
        # A temp file to save the output.
        # This file will be removed later by renderImage
        outfile <- tempfile(fileext = '.png')
        
        png(outfile, width = 800, height = 600)
        
        coloc <- colocalized(
          proc_values()$x5$data_file_selected,
          mz = input$corr_mz,
          n = input$corr_n,
          sort.by = input$corr_sort
        )
        print((coloc))
       # browser()
        
        
        
        p1 <- image(
          proc_values()$x5$data_file_selected,
          mz = coloc$mz[1:input$corr_plot_n],
          layout=c( length(runNames(proc_values()$x5$data_file_selected)), input$corr_plot_n),
          contrast.enhance = "histogram",
          normalize.image = "linear",
          colorscale = pals::parula(255),
          colorkey=input$key_on
        )
        
        
        
        
        print(p1) #, layout = n2mfrow(length(p1$facets)))
        
        
        dev.off()
        
        # Return a list containing the filename
        list(
          src = outfile,
          contentType = 'image/png',
          width = 800,
          height = 600,
          alt = "This is alternate text"
        )
      }, deleteFile = TRUE)
    })
    
    
  })
}