### R/PeakPickServer.R

PeakPickServer <- function(id, setup_values) {
  moduleServer(id, function(input, output, session) {
    #for dynamic UI
    ns = session$ns
    
    #import parallel mode
    par_mode = reactive({
      setup_values()[["par_mode"]]
    })
    
    #import chunk size
    chunks = reactive({
      setup_values()[["chunks"]]
    })
    
    
    has.new.files <- function() {
      unique(list.files(setup_values()[["wd"]]), recursive=T)
    }
    get.files <- function() {
      list.files(setup_values()[["wd"]], recursive=T)
    }
    
    # store as a reactive instead of output
    my_files <-
      reactivePoll(10, session, checkFunc = has.new.files, valueFunc = get.files)
    
    # any time the reactive changes, update the selectInput
    observeEvent(my_files(),
                 ignoreInit = T,
                 ignoreNULL = T,
                 {
                   print(grep(
                     ".imzML",
                     my_files(),
                     ignore.case = T,
                     value = T
                   ))
                   updateSelectInput(
                     session,
                     ns('peakPickfile'),
                     choices = grep(
                       ".imzML",
                       my_files(),
                       ignore.case = T,
                       value = T
                     )
                   )
                 })
    
    
    #change depending on value of input$peak_pick_status, pp_y is .imzML, pp_old is .rds
    output$pk_file <- renderUI({
      
      # Determine the list of files based on the value of input$peak_pick_status
      file_choices <- switch(
        input$peak_pick_status,
        "pp_old" = grep(".rds", my_files(), ignore.case = TRUE, value = TRUE),
        "pp_y" = grep(".imzML|.rds", my_files(), ignore.case = TRUE, value = TRUE),
        "pp_ref" = grep(".imzML|.rds", my_files(), ignore.case = TRUE, value = TRUE),
        # Add default or other cases if needed
        character(0) # Return an empty character vector if no match
        )
      
      selectInput(
        ns("peakPickfile"),
        "Imageset with peaked peaks",
        file_choices)
      
    })
    
    output$add_file <- renderUI({
      selectInput(
        ns("peakAddfile"),
        "Imageset to add to peak file",
        grep(
          ".rds|.imzML",
          my_files(),
          ignore.case = T,
          value = T
        )
      )
    })
    
    
    
    output$peak_pick <- renderUI({
      switch(
        input$pp_operation,
        "open_file" = list(
          radioButtons(
            ns("peak_pick_status"),
            "Cardinal v3.6+ processed .imzML/.rds file?",
            
            c("Yes" = "pp_y", "No, older .rds" = "pp_old")), 
          uiOutput(ns('pk_file')),
          actionButton(ns("action"), label = HTML("Restore saved file")),
          shinyFiles::shinySaveButton(ns("save_imzml"), "Save imzML File", "Save", filetype = list(""))
        ),
        "pp_raw" = list(
          radioButtons(
            ns("peak_pick_status"),
            "Start from an existing peak_picked file?",
            c(
              "No, fresh pick raw files" = "pp_no",
              "From Reference file" = "pp_ref",
              "From Mean Spectrum" = "pp_mean"
            )
          ),
          uiOutput(ns('pk_file_condition')),
          textInput(
            ns("selected_mz"),
            "Features selected / filtered (by comma or colon)"
          ),
          # downloadButton(
          #   ns("save_selected"),
          #   "Save selected runs from peak picked file"
          # ),
          p(),
          shinyFiles::shinySaveButton(ns("save_imzml"), "Save imzML File", "Save", filetype = list(""))
        ),
        "subset_f" = list(
          uiOutput(ns('pk_file')),
          actionButton(ns("action"), label = HTML("Restore first file")),
          
          uiOutput(ns('add_file')),
          actionButton(ns('action_add_file'), "Add file (same coordinates)"),
          
          
          textInput(
            ns("selected_mz"),
            "Features selected / filtered (by comma or colon)"
          ),
          # downloadButton(ns("save_selected"), "Save selected runs and features"),
          shinyFiles::shinySaveButton(ns("save_imzml"), "Save imzML File", "Save", filetype = list(""))
        ),
        "add_same_pklist" = list(
          uiOutput(ns('pk_file')),
          actionButton(ns("action"), label = HTML("Restore first file")),
          
          uiOutput(ns('add_file')),
          actionButton(
            ns('action_add_file_same'),
            "Add Imageset (same peaklist)"
          ),
          textInput(
            ns("selected_mz"),
            "Features selected / filtered (by comma or colon)"
          ),
          # downloadButton(
          #   ns("save_selected"),
          #   "Save selected runs and features from combined data"
          # ),
          shinyFiles::shinySaveButton(ns("save_imzml"), "Save imzML File", "Save", filetype = list(""))
        )
        # "demo" = list(selectInput(
        #   ns("cardworkdat"),
        #   label = "Choose sample dataset",
        #   choices =
        #     c(
        #       "Human Renal Cell Carcinoma (RCC)",
        #       "Whole Pig Fetus Cross-Section"
        #     )
        #   ),
        #   actionButton(ns("action_demo"), label = HTML("Load demo data"))
        #)
        
        
        
        
      )
    })
    
    
    
    
    output$pk_file_condition <- renderUI({
      req(input$peak_pick_status)
      switch(
        input$peak_pick_status,
        "pp_no" = list(
          # numericInput(
          #   ns("pix_for_peak_picking"),
          #   "% of pixels to randomly sample for peak picking",
          #   value = 10,
          #   min = 1,
          #   max = 100
          # ),
          numericInput(ns("SNR"), "S/N for overview peak picking", 20, min =
                         0),
          numericInput(
            ns("freq_min"),
            "Minimum peak frequency (0-1)",
            0.01,
            min = 0,
            max = 1
          ),
          textInput(
            ns("pp_method"),
            "Peak picking method (“diff”, “sd”, “mad”, “quantile”, “filter”, “cwt”)",
            value = "diff",
            placeholder = "Cardinal pp method, adaptive, mad, or simple"
          ),
          actionButton(ns("action"), label = HTML("Start peak picking"))
        ),
        "pp_ref" = list(uiOutput(ns('pk_file')),
                        actionButton(
                          ns("action"), label = HTML("Start peak binning")
                        )),
        "pp_mean" = list(
          #uiOutput('pk_file'),
          numericInput(ns("SNR"), "S/N for overview peak picking", 10, min =
                         0),
          selectInput(
            ns("pp_method"),
            "Peak picking method",
            choices = c("diff", "sd", "mad", "quantile", "filter", "cwt"),
            selected = "diff"
          ),
          # numericInput(
          #   ns("freq_min"),
          #   "Minimum peak frequency (0-1)",
          #   0.01,
          #   min = 0,
          #   max = 1
          # ),
          actionButton(ns("action"), label = HTML("Start peak binning"))
        )
      )
      
    })
    
    
    
    x0 <- reactiveValues(overview_peaks = NULL)
    
    

    observeEvent(input$action_demo, {
      req(input$cardworkdat)
      
      withProgress(message = "Loading demo data",
                   value = 0.5,
                   detail = "",
                   {
                     data_in <- switch(
                       input$cardworkdat,
                       "Human Renal Cell Carcinoma (RCC)" = "rcc",
                       "Whole Pig Fetus Cross-Section" = "pig206"
                     )
                     
                     #browser()
                     
                     x1 = reactiveValues(
                       raw_list = NULL
                     )
                     
                     if (data_in == "pig206") {
                       data(pig206, package = "CardinalWorkflows")
                       #overview_peaks <- as(pig206, "MSImagingExperiment")
                       
                       x1$raw_list<- as(pig206, "MSImagingExperiment")
                       
                     } else if (data_in == "rcc") {
                       data(rcc, package = "CardinalWorkflows")
                       
                       x1$raw_list<-as(rcc, "MSImagingExperiment")
                       #overview_peaks <- as(rcc, "MSImagingExperiment")
                       
                     }
                     
                     
                     
                     #x0$overview_peaks <- overview_peaks
                   })
      
    })
    
    overview_peaks <- eventReactive(input$action, {
      req(input$peak_pick_status)
      print(input$peak_pick_status)
      
      
      
      x0$overview_peaks<-NULL
      
      #start clock
      ptm<-proc.time()
      
      
      withProgress(message = "Performing peak peaking / binning",
                   value = 0.5,
                   detail = "Can take a while for large datasets...",
                   {
                     if (!is.null(input$peak_pick_status) &
                         input$peak_pick_status == "pp_y") {
                       if (length(input$peakPickfile) < 1) {
                         print("need file to open")
                         return()
                         
                       }
                       
                       
                       #check filetype and open accordingly
                       if (length(grep(".rds", input$peakPickfile, ignore.case = T)) > 0) {
                         overview_peaks <- readRDS(input$peakPickfile)
                         print(overview_peaks)
                       } else if (length(grep(".imzML", input$peakPickfile, ignore.case = T) >0)) {
                         overview_peaks <- readImzML(input$peakPickfile)
                         print(overview_peaks)
                       } else {
                         print("file type not recognized")
                         return()
                       }
                       
                       # overview_peaks <- readImzML(input$peakPickfile)
                       # print(overview_peaks)
                       
                       #browser()
                       
                       #De novo peak picking
                     }else if (input$peak_pick_status == "pp_old") {
                       if (length(input$peakPickfile) < 1) {
                         print("need file to open")
                         return()
                         
                       }
                       
                       old_dat<-readRDS(input$peakPickfile)
                       overview_peaks <- convert_card(old_dat)
                       print(overview_peaks)
                       
                     } else if (input$peak_pick_status == "pp_no") {
                       
                       
                       
                       x1 = setup_values()[["x1"]] # bring in raw_files list
                       
                       if (is.null(x1$raw_list)) {
                         message("Please choose raw data in Data Setup tab for peak picking")
                         return()
                       }
                       
                       setCardinalBPPARAM(par_mode())
                       setCardinalNChunks(chunks())
                       
                       # commented code not relevant to Cardinal 3.6
                       # coord_list_reduced <-
                       #   lapply(x1$raw_list, function(x)
                       #     coord(x)[(sample(
                       #       nrow(coord(x)),
                       #       dim(x)[2] * input$pix_for_peak_picking / 100
                       #     )), ])
                       # 
                       # #function to use raw image list and coordinates to create reduced raw set
                       # select_pix <-
                       #   function(index, raw_list, coord_set) {
                       #     raw_list[[index]][!is.na(prodlim::row.match(
                       #       as.data.frame(coord(raw_list[[index]])),
                       #       as.data.frame(coord_set[[index]])
                       #     ))]
                       #   }
                       # 
                       # #combine reduced images.
                       # test_raw_reduced <-
                       #   try(
                       #   Cardinal::combine(lapply(1:length(x1$raw_list), function(x)
                       #     select_pix(x, x1$raw_list, coord_list_reduced)))
                       #   )
                       # #raw_reduced<-combine(bplapply(1:length(raw_list), function(x) select_pix(x,raw_list, coord_list_reduced)))
                       # 
                       #  if(class(test_raw_reduced) %in% "try-error"){
                       #    message("failed to combine raw files, check files")
                       #    showNotification("failed to combine raw files, check files", type="error")
                       #    return()
                       #  }
                       # 
                       
                       #use Cardinal 3.6 peak processing method
                       
                       msa<-convertMSImagingArrays2Experiment(Cardinal::combine(lapply(x1$raw_list, convertMSImagingExperiment2Arrays)),
                                                                                              #mass.range = c(setup_values()[["mz_max"]], setup_values()[["mz_min"]]), 
                                                                                              x1$mass.range,
                                                                                              resolution = setup_values()[["res"]], 
                                                                                              units = "ppm")
                       
                       mse_queue<- msa |>
                         normalize() |>
                         #smooth() |>
                         #reduceBaseline() |>
                         peakPick(SNR=input$SNR, method=input$method, type="area", tolerance=NA, units="ppm")
                        
                       #plot the middle spectrum?
                       print(plot(mse_queue, i=round(dim(coord(msa))[1]/2,0), linewidth=2))
                       
                       test_mz_reduced<-try(
                          peakAlign(mse_queue, tolerance= setup_values()[["tol"]], units=setup_values()[["units"]]) %>%
                          subsetFeatures( freq > input$freq_min)%>% 
                          summarizeFeatures()
                       )
                          

                         # peakProcess(
                         #   convertMSImagingArrays2Experiment(Cardinal::combine(lapply(x1$raw_list, convertMSImagingExperiment2Arrays)),
                         #                                     mass.range = c(setup_values()[["mz_max"]], setup_values()[["mz_min"]]), 
                         #                                     resolution = setup_values()[["res"]], 
                         #                                     units = "ppm"),
                         #   #combine(x1$raw_list),
                         #   method=input$pp_method,
                         #   SNR=input$SNR,
                         #   #SNR=5,
                         #   filterFreq = input$freq_min,
                         #   tolerance=setup_values()[["tol"]],
                         #   units="ppm",
                         #   type="area",
                         #   sampleSize=input$pix_for_peak_picking/100
                        
                       
                       
                       # 
                       # #create list of peaks for referencing
                       # test_mz_reduced <- try(
                       #   HTS_reproc(
                       #   data1 = test_raw_reduced,
                       #   SN = input$SNR,
                       #   res = setup_values()[["res"]],
                       #   align_tol = setup_values()[["tol"]],
                       #   method = input$pp_method,
                       #   freq.min = input$freq_min
                       #  )
                       # )
                       
                       if(class(test_mz_reduced) %in% "try-error"){
                         message("peak picking failed, check files")
                         showNotification("peakpicking failed, check files")
                         return()
                       }
                       # 
                       # test_mz_reduced
                       # message("Saving mz_list_ref.RData file with peaklist used for peakpicking")
                       # saveRDS(test_mz_reduced, file = paste("mz_list_ref.RData"))
                       # browser()
                       # #create peak list from reference   
                       # proc_list<-lapply(
                       #   x1$raw_list,
                       #   ref_to_peaks_3_6,
                       #   mz_ref = mz(test_mz_reduced),
                       #   tol = setup_values()[["tol"]]
                       # )
                       overview_peaks <- test_mz_reduced
                       #browser()
                       #use mz list from existing file
                     } else if (input$peak_pick_status == "pp_ref") {
                       x1 = setup_values()[["x1"]] # bring in raw_files list
                       if (is.null(x1$raw_list)) {
                         print("Please choose raw data in Data Setup tab for peak picking")
                         return()
                       }
                       
                       if (input$peakPickfile == '') {
                         print("need file to open")
                         return()
                       }
                       #browser()
                       print("binning raw files from reference")
                       #check for rds or imzml and open accordingly
                       if(length(grep(".rds", input$peakPickfile, ignore.case = T, value = T) == 1)) {
                         test_mz_reduced <- readRDS(input$peakPickfile)
                       } else if (length(grep(".imzML", input$peakPickfile, ignore.case = T, value = T) == 1)) {
                         test_mz_reduced <- readImzML(input$peakPickfile)
                       } else {
                         print("file type not recognized")
                         return()
                       }
                       
                       print(test_mz_reduced)
                       setCardinalBPPARAM(par_mode())
                       setCardinalNChunks(chunks())
                       
                       
                       a<-Cardinal::combine(lapply(x1$raw_list, convertMSImagingExperiment2Arrays))
                       
                       overview_peaks <- try(a %>% normalize() %>% peakProcess(
                                                         ref=mz(test_mz_reduced),
                                                         SN=input$SNR,
                                                         type="area",
                                                         tolerance=setup_values()[["tol"]], units=setup_values()[["units"]]) %>% process() %>% summarizeFeatures()
                       )
                       
                       if(class(overview_peaks) %in% "try-error"){
                         message("Peak binning from reference failed, please check files")
                         showNotification("Peak binning from reference failed, please check files", type="error")
                         return()
                       }
                       
                       #calculate mean spectrum for peak list
                     } else if (input$peak_pick_status == "pp_mean") {
                       x1 = setup_values()[["x1"]] # bring in raw_files list
                       if (is.null(x1$raw_list)) {
                         print("Please choose raw data in Data Setup tab for peak picking")
                         return()
                       }
                       
                       
                       
                       
                       if(!is.null(x1$mass.range)) {
                         #check mz_min and mz_max to make sure they are not the same
                         if (setup_values()[["mz_min"]] == setup_values()[["mz_max"]]) {
                           print("mz_min and mz_max are the same, please change")
                           showNotification("mz_min and mz_max are the same, please change", type="error")
                           return()
                         }
                       }
                       
                       
                       #do.call(cbind, x1$raw_list[c(1:3))
                       
                       setCardinalBPPARAM(par_mode())
                       setCardinalNChunks(chunks())
                       print("binning raw files from mean spectrum")
                       test_mz_mean <- try(
                         #for debugging
                         #Cardinal::combine(x1$raw_list[c(1:3,5:8)]) %>%
                         Cardinal::combine(lapply(x1$raw_list, convertMSImagingExperiment2Arrays)) %>%
                         
                           convertMSImagingArrays2Experiment(
                             #mass.range=c(setup_values()[["mz_max"]], setup_values()[["mz_min"]])) %>% 
                             mass.range=x1$mass.range) %>%
                           # normalize() %>% 
                           # process() %>%
                         estimateReferencePeaks( SNR=input$SNR, 
                                                method=input$pp_method)
                       )
                       
                       print(test_mz_mean)
                       #add check to see if 0 peaks picked...
                       
                       #check class of test_mz_mean
                       if(class(test_mz_mean) %in% "try-error") {
                         print("mean spectrum peak picking failed, check input files")
                         showNotification("Mean spectrum peak picking failed, check input files.", type = "error")
                         return()
                       }
                       
                       test_mz_reduced<-try(Cardinal::combine(lapply(x1$raw_list, convertMSImagingExperiment2Arrays)) %>% 
                                              normalize() %>% 
                                              peakProcess(
                                                 ref=mz(test_mz_mean),
                                                 SN=input$SNR,
                                                 type="area",
                                                 tolerance=setup_values()[["tol"]], units=setup_values()[["units"]]) %>% 
                                               summarizeFeatures()
                       )

                       #test_mz_reduced  <- summarizeFeatures(test_mz_reduced)


                       if(class(test_mz_reduced) %in% "try-error") {
                         print("peak picking failed, check input files")
                         showNotification("Mean spectrum peak picking failed, check input files.", type = "error")
                         return()
                       }
                       
                       setCardinalBPPARAM(par_mode())
                       setCardinalNChunks(chunks())
                       
                       overview_peaks <-
                         test_mz_reduced
                       
                       print(test_mz_reduced)
                       #browser() #number of peaks not changing with SN?
                       
                     }
                     
                     #reset UMAP / SSC segmentation parameters
                     x2$mytable_selected <- NULL
                     x2$tf_list <- NULL
                     x2$tf_list_anat <- NULL
                     x2$tf_list_umap <- NULL
                     x2$data_list <- NULL
                     x2$list_proc_img <- NULL
                     x2$pdat_anat <- NULL
                     x2$ssc <- NULL
                     
                     
                     x0$overview_peaks <- overview_peaks
                     #end clock
                     print(proc.time()-ptm)
                     return(overview_peaks)
                     
                   })
      
    
      
      
      
    })
    
    
    
    
    observe(overview_peaks()) #not sure if this is needed or can use x0
    #observe(x0$overview_peaks)
    
    # create subset from peaks within original dataset
    observeEvent(input$action_add_file, {
      #check for .rds or .imzML and open accordingly
      if(length(grep(".rds", input$peakAddfile, ignore.case = T, value = T) == 1)) {
        add_peaks <- readRDS(input$peakAddfile)
        print(add_peaks)
      } else if (length(grep(".imzML", input$peakAddfile, ignore.case = T, value = T) == 1)) {
        add_peaks <- readImzML(input$peakAddfile)
        print(add_peaks)
      } else {
        print("file type not recognized")
        return()
      }
      # add_peaks <- readImzML(input$peakAddfile)
      # print(add_peaks)
      # 
      
      #check to make sure overview_peaks() exists
      
      if (dim(x0$overview_peaks)[2] != dim(add_peaks)[2]) {
        print("inconistent pixel sizes, make sure same coordinate set being used")
        return()
      }
      
      #x<-overview_peaks()[sort(sample(1:3895, 1000)),]
      x <- x0$overview_peaks
      y <- add_peaks
      
      #remove all featureData other than IDs
      fData(x) <- fData(x)[, colnames(fData(x)) %in% "ID"]
      fData(y) <- fData(y)[, colnames(fData(y)) %in% "ID"]
      
      mz_sort <- sort(c(mz(x), mz(y)))
      
      if (sum(mz_sort %in% mz(x)) == length(mz_sort)) {
        message("Datasets do not appear to have unique peaklists, aborting")
        return(NULL)
      }
      #message("debugging file subset features")
      #browser()
      
      print("combining file to already opened file--slow for very large files")
      
      #create first mz, try first file to find mz. This is to avoid errors later with empty dataset
      dat<-(matrix(NA, ncol = ncol(x), nrow = length(mz_sort)))

      # Find indices to insert values from x and y
      idx_x <- match(mz(x), mz_sort)
      idx_y <- match(mz(y), mz_sort)
      
      # Insert values from x and y into the merged dataframe
      dat[idx_x, ] <- (iData(x))
      dat[idx_y, ] <- iData(y)
      
      dat2<-try(Cardinal::MSImagingExperiment( imageData = dat, 
                                     featureData = MassDataFrame(mz_sort), 
                                     pixelData = pData(x0$overview_peaks),
                                     centroided = TRUE)
      )
      
      if(class(dat2)%in% "try-error") {
        print("Combine failed, ensure mass list is not duplicated and coordinates are the same")
        showNotification("Combine failed, ensure mass list is not duplicated and coordinates are the same", type="error")
        return()
      } else {
        dat<-dat2
      }
      
      
      
      
      x0$overview_peaks <- dat
      
    })
    
    observeEvent(input$action_add_file_same, {
      if (input$peakAddfile == '') {
        print("need file to open")
        return()
        
      }
      #check filename and open accordinlgy
      if(length(grep(".rds", input$peakAddfile, ignore.case = T, value = T) == 1)) {
        add_peaks <- readRDS(input$peakAddfile)
        print(add_peaks)
      } else if (length(grep(".imzML", input$peakAddfile, ignore.case = T, value = T) == 1)) {
        add_peaks <- readImzML(input$peakAddfile)
        print(add_peaks)
      } else {
        print("file type not recognized")
        return()
      }
      
      #browser()
      # add_peaks <- readImzML(input$peakAddfile)
      # print(add_peaks)
      # 
      
      #check to make sure overview_peaks() exists and has the same number of peaks
      
      if (!identical(mz(x0$overview_peaks), mz(add_peaks))) {
        print("inconistent mass lists, make sure same peak lists being used")
        return()
      }
      
      #x<-overview_peaks()[sort(sample(1:3895, 1000)),]
      x <- x0$overview_peaks
      y <- add_peaks
      
      
      print("combining files")
      
      
      
      if (!setequal(colnames(pData(x)), colnames(pData(y)))) {
        cols <- union(colnames(pData(x)), colnames(pData(y)))
        cols_overview <- cols[!cols %in% (colnames(pData(x)))]
        cols_selected <- cols[!cols %in% (colnames(pData(y)))]
        
        pData(x)[, cols_overview] <- "NA"
        pData(y)[, cols_selected] <- "NA"
      }
      
      dat<-combine_card(list(x,y))

      x0$overview_peaks <- dat
      
    })
    
    
    #create reactive object which contains info for  (run_table)
    #a color vector for viz / selection, and intermediate imaging dataset (overview_peaks_sel)
    #and the data_list used to store UMAP/SSC results and color coding
    x2 = reactiveValues(
      run_table = NULL,
      #pixel selection table, Overview Analysis
      bkcols = NULL,
      #umap selection colors
      overview_peaks_sel = NULL,
      #imaging dataset after pixels selected in Overview Analysis
      mytable_selected = NULL,
      #imaging dataset used for UMAP (selected by run/rowid)
      pdat_anat = NULL,
      #pdata information corresponding to mytable_selected for annotation
      data_list = NULL,
      #data list containing UMAP data and color clustering schemes
      ssc = NULL,
      #ssc object
      ssc_models = NULL,
      #ssc models
      tf_list = NULL,
      #T/F list to select rows / pixels for segmentation.
      tf_list_anat = NULL,
      #TF list when performing anatomy segmentation
      tf_list_umap = NULL,
      #TF list fixed after UMAP runs for anatomy segmentation
      umap_name = NULL,
      #names of images used for UMAP
      list_proc_img = NULL,
      #list of processed runs after background subtraction to save/continue processing
      rcol_plot = NULL
    )       #variable to record plot status for UMAP colors
    observe({
      req(x0$overview_peaks)
      #browser()
      a <- runNames(isolate(x0$overview_peaks))
      xmax = NULL
      ymax = NULL
      for (i in 1:length(a)) {
        img.dat <-
          x0$overview_peaks %>% subsetPixels(Cardinal::run(x0$overview_peaks) %in% a[i])
        xmax[i] = max(coord(img.dat)$x)
        ymax[i] = max(coord(img.dat)$y)
      }
      
      run_table <-
        data.frame(
          runs = a,
          index = 1:length(a),
          xmin = 0,
          xmax = xmax,
          ymin = 0,
          ymax = ymax
        )
      
      x2$run_table <- run_table
    })
    
    
    
    output$peak_pick_selection = DT::renderDT(
      #output$peak_pick_selection = DT::renderDataTable({
      #datatable(
      req(x2$run_table),
      selection = list(mode = "multiple", selected = c(1:nrow(x2$run_table))),
      caption = HTML(
        "<b/>Select files for further processing and edit pixels to remove noise:"
      ),
      editable = TRUE
    )
    
    
    
    proxy = DT::dataTableProxy('peak_pick_selection')
    
    observeEvent(input$peak_pick_selection_cell_edit, {
      info = input$peak_pick_selection_cell_edit
      str(info)
      i = info$row
      j = info$col
      v = info$value
      x2$run_table[i, j] <-
        isolate(DT::coerceValue(v, x2$run_table[i, j]))
      #replaceData(proxy, run_table(), resetPaging = FALSE)  # important
    })
    
    output$overview_text <- renderText({
      x2$run_table$min
    })
    
    
    observe({
      req(x0$overview_peaks)
      #browser()
      a <- runNames((x0$overview_peaks))
      b <- list()
      for (i in 1:length(a)) {
        tmp <- subsetPixels(x0$overview_peaks,
                            Cardinal::run(x0$overview_peaks) %in% a[i])
        b[[i]] <- subsetPixels(
          tmp,
          x > x2$run_table$xmin[i] ,
          x < (x2$run_table$xmax[i] + 1) ,
          y > x2$run_table$ymin[i],
          y < (x2$run_table$ymax[i] + 1)
        )
      }
      
      #check if freq data is getting lost here
      #browser()
      
      #do.call doesn't always work, so lets use the loop to recombine the data
      
      
      x2$overview_peaks_sel <- combine_card(b)
      

      plot_card_server("card_plot", overview_peaks_sel = x2$overview_peaks_sel)
      
      
    })
    
    
    #save the modified files (subset, added etc) and use the Features selected variable
    # to filter which features are saved.
    
    observeEvent(input$save_imzml, {
      
      req(input$save_imzml)
      
      
      volumes <- c(wd = setup_values()[["wd"]], home = fs::path_home())
      shinyFiles::shinyFileSave(input, "save_imzml", roots = volumes, session = session)
      
      
      save_path <- shinyFiles::parseSavePath(volumes, input$save_imzml)
      if (nrow(save_path) == 0) return(NULL)
      
      #browser()
      filen <- as.character(save_path$datapath)
      # 
      #   #Save state template from here: https://www.r-bloggers.com/2019/06/shiny-application-with-modules-saving-and-restoring-from-rds/
      # output$save_selected <- downloadHandler(
      #   filename = function() {
      #     paste0(getwd(),
      #            "/MSI-proc_hi_SN-peak_picked-subset-",
      #            Sys.Date(),
      #            ".rds")
      #   },
      #   content = function(file) {
          #restore original fData frame since we are only dealing with pixels
          fData(x2$overview_peaks_sel) <- fData(x0$overview_peaks)
          
          
          #dataset feature list
          features <- 1:nrow(x2$overview_peaks_sel)
          #check to see if we're selecting a subset of the features by looking at the mz input field
          if (!is.null(input$selected_mz) &&
              input$selected_mz != "") {
            if (grepl(",|:", input$selected_mz)) {
              features_sel <-
                sort(as.numeric(eval(parse(
                  text = paste0("c(", input$selected_mz, ")")
                ))))
            } else {
              features_sel <-
                sort(as.numeric(unlist(
                  strsplit(input$selected_mz, "\\s+")
                )))
            }
          } else {
            features_sel <- features
          }
          #pk_img <- x0$overview_peaks
          #browser()
          #save subset of data
          pk_img <- x2$overview_peaks_sel[features_sel, ] %>%
            subsetPixels(Cardinal::run(x2$overview_peaks_sel) %in% runNames(x2$overview_peaks_sel)[input$peak_pick_selection_rows_selected])
          #print(pk_img)
          print("saving selected subset of image file")
          print(pk_img)
          print(filen)
          #TODO - why is the ID column replicated??
          #saveRDS(pk_img, file)
          
          
          writeImzML(pk_img, filen)
          saveRDS(pk_img, paste0(filen, ".rds"))
        }
      )
    
    
    
    
    #browser()
    
    preproc_values <- reactive({
      list(x2 = x2,
           x0 = x0)
    })
    return(preproc_values)
    
  })
}
