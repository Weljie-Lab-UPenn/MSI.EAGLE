### R/DataSetupServer.R
library(magrittr)
DataSetupServer <- function(id, rawd, wd) {
  moduleServer(id, function(input, output, session) {
    
    ns = session$ns #for dyanamic variable namespace
    
    
    output$pp_params_display <- renderUI({
      
      
      
      # Define the dynamic parameters based on the switch
      # Any NULL in mass range will force auto detection
      params <- switch(
        input$pp_params,
        "qtof1" = list(
          res = NA, 
          tol = 30, 
          mass_range_min = NULL, 
          mass_range_max = NULL
        ),
        "hires" = list(
          res = NA, 
          tol = 10, 
          mass_range_min = NULL, 
          mass_range_max = NULL
        )
      )
      
        
        tagList(
          # Use the parameters to render UI elements dynamically
          numericInput(ns("res"), "resolution (ppm)", params$res),  
          
          fluidRow(
            column(6, 
                   numericInput(ns("tol"), "tolerance", params$tol)),
            column(6, 
                   selectInput(ns("units"), "units", choices = c("ppm", "mz"), selected = "ppm"))
          ),
          fluidRow(
            column(6,
                   numericInput(
                     ns("mass_range_min"),
                     label = p("m/z min for import"),
                     min = 0,
                     value = params$mass_range_min
                   )),
            column(6,
                   numericInput(
                     ns("mass_range_max"),
                     label = p("m/z max for import"),
                     min = 0,
                     value = params$mass_range_max
                   ))
          )
        )
    
      
      
    })
    
    
    # You can access the value of the widget with input$text, e.g.
    output$variables <-
      renderPrint({
        cat(paste("raw folder= \t", input$folder),
            "\n",
            paste("working dir= \t\t", input$wd))
      })
    
    
    root_dir <- wd
    observe({
      #shinyFiles::shinyDirChoose(input, "wd_new", roots = volumes, session = session)
      shinyFiles::shinyDirChoose(input,
                                 "wd_new",
                                 roots = c(wd = '.'),
                                 session = session)
    })
    #create list with both files and outfiles variables
    files <- reactive({
      setwd(input$wd)
      print(paste("wd is", getwd()))
      
      
      #raw files folder
      folder <- input$folder
      
      
      
      
      #read list of files to process; end result should be list of imzML files to process
      files <-
        list.files(path = folder,
                   recursive = F,
                   pattern = input$regex)
      # files <-
      #   sub("\\.imzML", "", files) #needed because Cardinal automatically adds this to the name
      files
      
    })
    
    #working directory
    
    
    observeEvent(input$wd_new, {
      print(getwd())
      print(input$new_wd[[2]])
      message(shinyFiles::parseDirPath(c(wd = root_dir), input$wd_new))
      updateTextInput(session, "wd", value = as.character(shinyFiles::parseDirPath(c(wd =
                                                                                       root_dir), input$wd_new)))
      
    })
    
    #create reactive object for initial data setup
    x1 <- reactiveValues(raw_list = NULL,
                         file_list = NULL,
                         mass.range = NULL)
    
    
    observe({
      if(!input$action_demo){
      
        x1$file_list <- try(as.data.frame(cbind(files = files(), index = 1:length(files()))))
        
       
      } else {
        x1$file_list <- as.data.frame(cbind(files = names(x1$raw_list), index = 1:length(x1$raw_list)))
      }
      
      #browser()
    })
    
    output$files = DT::renderDataTable({
      if(length(dim(x1$file_list))==2){
        #output$files = DT::renderDT({
        DT::datatable(
          #tab(),
          x1$file_list,
          #selection = list(mode = "multiple", selected=c(1:length(files()))),
          selection = 'none',
          caption = "Select files for further processing",
          #editable = TRUE,
          
          extensions = c("Buttons", "Select"),
          options = list(
            dom = 'Bfrtip',
            select = TRUE,
            buttons = list('pageLength', "copy", "selectNone", "selectAll")
          )
        )
        
      }
      
    },
    server = FALSE)
    
    
    files_selected <- reactive({
      ids <- input$files_rows_selected
      files()[ids]
    })
    
    
    output$par_mode_setup <- renderUI({
      # browser()
      #
      if (as.character(Sys.info()['sysname']) == "Windows") {
        pmode = "pc_p"
      } else if (as.character(Sys.info()['sysname']) == "Darwin") {
        pmode = "mac_p"
      } else {
        pmode = "no_p"
      }
      
      radioButtons(
        inputId =  ns("parallel_mode"),
        label =  "Parallelization mode",
        choices = c(
          "None" = "no_p",
          "MulticorParam" = "mac_p",
          "SnowParam" = "pc_p"
        ),
        #selected = "no_p")
        selected = pmode
      ) #, "FutureParam"="fut_p"))
      
      
    })
    
    
    par_mode <- reactive({
      if (!is.null(input$parallel_mode)) {
        par_mode_in <- switch(
          input$parallel_mode,
          "no_p" = SerialParam(),
          "mac_p" = MulticoreParam(workers = input$ncores),
          "pc_p" = SnowParam(workers = input$ncores),
          "fut_p" = FutureParam()
        )
        return(par_mode_in)
      }
    })
    
    
    
     
    observeEvent(input$action1, {
      if (length(files_selected()) < 1) {
        showNotification("No files selected for raw data input! Please select runs from table.", type="error", duration = 10)
        print("no files selected for raw data input!")
        return()
      }
      
      #browser()
      withProgress(message = 'Importing data', value = 1, {
        
        
        #if input$mass_range_min and max are  NULL, set mass.range to NULL
        if (!is.numeric(input$mass_range_min) | !is.numeric(input$mass_range_max)) {
          message("mass range is set to NULL")
          mass.range <- NULL
          x1$mass.range <- NULL
        } else {
          mass.range <- c(input$mass_range_min, input$mass_range_max)
          x1$mass.range <- mass.range
        }

        x1$raw_list <- sapply(files_selected(), function(x)
        {
          y <- Cardinal::readMSIData(
            paste0(input$folder,"//",x),
            #folder=input$folder,
            resolution = input$res,
            units = input$units,
            mass.range = mass.range
          )
          #Cardinal::centroided(y) <- TRUE
          coord(y)$z <- NULL
          return(y)
        })
        
        
        #print(x1$raw_list[[1]])
        #graphics.off()
        names(x1$raw_list) <- files_selected()
        
      })
    })
    
    
    img.dat <- eventReactive(input$action2, {
      
      if (is.null(x1$raw_list))
        return()
       setCardinalBPPARAM(par_mode())
       setCardinalNChunks(input$chunks)
     
      
      withProgress(message = "Extracting sample for plotting", value = 1, {
        tmp.img <- Cardinal::combine(lapply(x1$raw_list, convertMSImagingExperiment2Arrays))
        tmp.img <- convertMSImagingArrays2Experiment(tmp.img, 
                                                     #mass.range = c(input$mass_range_min, input$mass_range_max), 
                                                     mass.range= x1$mass.range,
                                                     resolution = input$res, 
                                                     units = input$units
                                                     )
        #tmp.img <- Cardinal::combine(x1$raw_list[1:length(files())])
        tmp.img<- tmp.img %>% 
           peakProcess(SNR=20, 
                       tolerance = input$tol, 
                       units = input$units, 
                       method="diff",
                       mass.range = x1$mass.range, 
                       sampleSize=input$pix_to_plot/100, filterFreq=0.2) %>% summarizeFeatures()
      })
    })
    
    observe({
      req(img.dat())
      #browser()
      plot_card_server("card_plot", overview_peaks_sel = img.dat())
    })
    
    plot_img.dat <- reactive({
      if (is.null(img.dat()))
        return()
      plot(img.dat())
      
    })
    
    output$plot_ranges <- renderUI({
      aa <- try(img.dat())
      if (class(aa) == "try-error")
        return()
      
      if (is.null(plot_img.dat()))
        return()
      
      sliderInput(
        ns("mass_range_plot"),
        label = p("m/z range for quick MS plot (X)"),
        min = plot_img.dat()$par$xlim[1],
        max = plot_img.dat()$par$xlim[2],
        value = plot_img.dat()$par$xlim,
        step = NULL,
        round = TRUE
      )
      
    })
    
    output$plot_ranges2 <- renderUI({
      if (is.null(plot_img.dat()))
        return()
      
      sliderInput(
        ns("int_range_plot"),
        label = p("intensity range for quick MS plot (Y)"),
        min = plot_img.dat()$channels$y$limits[1],
        max = plot_img.dat()$channels$y$limits[2],
        value = plot_img.dat()$channels$y$limits,
        step = plot_img.dat()$channels$y$limits[2] / 20,
        round = TRUE
      )
    })
    
    output$plot <- renderImage({
      if (is.null(img.dat()))
        return()
      
      # A temp file to save the output.
      # This file will be removed later by renderImage
      outfile <- tempfile(fileext = '.png')
      
      png(outfile, width = 400, height = 300)
      
      
      print(plot(
        img.dat(),
        xlim = input$mass_range_plot,
        ylim = input$int_range_plot
      ))
      
      dev.off()
      
      # Return a list containing the filename
      list(
        src = outfile,
        contentType = 'image/png',
        width = 400,
        height = 300,
        alt = "This is alternate text"
      )
    }, deleteFile = TRUE)
    
    
    
    #plot of
    output$plot2 <- renderImage({
      req(img.dat())
      # A temp file to save the output.
      # This file will be removed later by renderImage
      outfile <- tempfile(fileext = '.png')
      
      png(outfile, width = 800, height = 600)
      
      ion <- switch(input$mode,
                    "p" = 381.2,
                    "n" = 255.2)
      print(
        Cardinal::image(
          img.dat(),
          mz = ion,
          plusminus = 1,
          colorscale = pals::parula(255),
          contrast.enhance = "histogram"
        )
      )
      
      dev.off()
      
      # Return a list containing the filename
      list(
        src = outfile,
        contentType = 'image/png',
        width = 800,
        height = 600,
        alt = "This is alternate text"
      )
      #browser()
    }, deleteFile = TRUE)
    
    observeEvent(input$action_demo, {
      req(input$cardworkdat)
      
      if(!require("CardinalWorkflows", quietly = TRUE)) {
        showNotification("CardinalWorkflows package is required to load demo data", type="error")
        message("CardinalWorkflows package is required to load demo data")
        message("can be installed as follows:\n
                if (!requireNamespace(\"BiocManager\", quietly=TRUE))\n
                install.packages(\"BiocManager\")\n

                BiocManager::install(\"CardinalWorkflows\")")
        return()
      }
      
      
      showNotification( "Loading demo data")
                   
                   {
                     data_in <- switch(
                       input$cardworkdat,
                       "Human Renal Cell Carcinoma (RCC) with background" = "rcc_whole",
                       "Human Renal Cell Carcinoma (RCC) (no background)" = "rcc",
                       "Whole Pig Fetus Cross-Section" = "pig206"
                     )
                     
                     #browser()
                     

                     
                     if (data_in == "pig206") {
                       data(pig206, package = "CardinalWorkflows")
                       #overview_peaks <- as(pig206, "MSImagingExperiment")
                       
                        tmp <- as(pig206, "MSImagingExperiment")
                        x1$raw_list<-lapply((Cardinal::runNames(tmp)), 
                                            function(x) tmp[,Cardinal::run(tmp) %in% x])
                        names(x1$raw_list)<-runNames(tmp)
                        
                     } else if (data_in == "rcc") {
                       data(rcc, package = "CardinalWorkflows")
                       
                       tmp<-as(rcc, "MSImagingExperiment")
                       
                       tmp<-tmp %>%
                         subsetPixels(!is.na(diagnosis))
                       
                       x1$raw_list<-lapply((Cardinal::runNames(tmp)), 
                                           function(x) tmp[,Cardinal::run(tmp) %in% x])
                       names(x1$raw_list)<-runNames(tmp)                       
                       
                       
                       
                       #overview_peaks <- as(rcc, "MSImagingExperiment")
                       
                     } else if (data_in == "rcc_whole") {
                       
                       data(rcc, package = "CardinalWorkflows")
                       
                       tmp<-as(rcc, "MSImagingExperiment")
                       
                       
                       x1$raw_list<-lapply((Cardinal::runNames(tmp)), 
                                           function(x) tmp[,Cardinal::run(tmp) %in% x])
                       names(x1$raw_list)<-runNames(tmp)                       
                       
                       
                     }
                     
                     
                     
                     #x0$overview_peaks <- overview_peaks
      }
      x1$file_list<-names(x1$raw_list)
      
    })
    
    setup_values <- reactive({
      list(
        wd = input$wd,
        rawd = input$rawd,
        ncores = input$ncores,
        chunks = input$chunks,
        x1 = x1,
        par_mode_out = setCardinalBPPARAM(par_mode()),
        par_mode = par_mode(),
        tol = input$tol,
        # used for pick picking. may be better not transferring
        res = input$res, #same as align_tol
        mz_max=input$mass_range_max,
        mz_min=input$mass_range_min,
        units=input$units
      )
    })
    return(setup_values)
    
  })
}
