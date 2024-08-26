### R/DepthAnalysisServer.R
DepthAnalysisServer <- function(id,  setup_values) {
  moduleServer(id, function(input, output, session) {
    # output$table <- DT::renderDataTable(mtcars)
    
    ns = session$ns #for dyanamic variable namespace
    
    #import parallel mode
    par_mode = reactive({
      setup_values()[["par_mode"]]
    })
    
    
    # dir_list<-reactive ({
    #   invalidateLater(20000)
    #   list.files(setup_values()[["wd"]])
    # })
    # replaced dir_list with my_files
    
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
                   #browser()
                   # print(grep(
                   #   ".imzML",
                   #   my_files(),
                   #   ignore.case = T,
                   #   value = T
                   # ))
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
    
    #Setup and values for the Depth processing tab
    x4 <- reactiveValues(
      seg_file = NULL,
      seg_filename = NULL,
      #to keep record of what file is used for coordinates
      depth_peaks = NULL,
      seg_pp_file = NULL
    )
    
    
    output$segmented_data <- renderUI({
      selectInput(
        ns("segmentation_file"),
        "Imageset with segmented coordinates",
        grep(
          "imzML",
          my_files(),
          ignore.case = T,
          value = T
        )
      )
    })
    
    output$pp_options <- renderUI({
      switch(
        input$targeted_pp,
        "untargeted" =
          list(
            numericInput(
              ns("pix_for_peak_picking2"),
              "% of pixels to for peak picking",
              100
            ),
            numericInput(ns("SNR2"), "S/N for comprehensive peak picking", 20),
            numericInput(
              ns("freq_min2"),
              "Minimum peak frequency (0-1",
              0.03,
              min = 0,
              max = 1
            ),
            numericInput(
              ns("tol2"),
              "tolerance for  alignment or peak binning (ppm)",
              30
            ),
            textInput(
              ns("pp_method2"),
              "Peak picking method (“diff”, “sd”, “mad”, “quantile”, “filter”, “cwt”)",
              value = "diff",
              placeholder = "Cardinal pp method, adaptive, mad, or simple"
            )
          ),
        "targeted" =
          list(
            fileInput(
              ns('masses_file'),
              "File with exact mass values info (.txt)",
              placeholder = "Must have 'exact.mass' column",
              accept = "text/plain"
            ),
            #must have column called exact.mass
            numericInput(ns("tol2"), "tolerance for peak binning (ppm)", 15)
          ),
        "mean" = list(
          numericInput(ns("SNR2"), "S/N for comprehensive peak picking", 10),
          textInput(
            ns("pp_method2"),
            "Peak picking method (“diff”, “sd”, “mad”, “quantile”, “filter”, “cwt”)",
            value = "mad",
            placeholder = "Cardinal pp method, adaptive, mad, or simple"
          ),
          numericInput(
            ns("freq_min2"),
            "Minimum peak frequency (0-1",
            0.03,
            min = 0,
            max = 1
          ),
          numericInput(ns("tol2"), "tolerance for peak binning (ppm)", 15)
        )
      )
    })
    
    
    observeEvent(input$action_seg, {
      x1 = setup_values()[["x1"]] # bring in raw_files list
      
      if(is.null(x1$raw_list)){
        message("raw files not selected, please choose from the Data Setup tab first!")
        showNotification("raw files not selected, please choose from the Data Setup tab first!", type="error")
        return()
      }
      
      x4$seg_file <- readImzML(input$segmentation_file)
      print(x4$seg_file)
      x4$seg_filename = input$segmentation_file
      
      
      
      if (is.null(x1$raw_list) | is.null(x4$seg_file)) {
        print(
          "Please choose raw data in Data Setup Tab and input segmentation tempate first for peakpicking"
        )
        return(NULL)
      } else if (sum(paste0(runNames(x4$seg_file), ".imzML") %in% names(x1$raw_list)) < 1) {
        
        print("no raw file compatible with dataset")
        print("runNames(x4$seg_file)")
        print(names(x1$raw_list))
        
        x4$seg_filename = "Choose a file compatible with the raw dataset"
        
      }
      
      
      
    })
    
    output$filename_header <- renderPrint({
      (paste0("working filename= ", req(x4$seg_filename)))
      #if(is.null(input$stats_input_file))
      #  return()
      #print(input$stats_input_file)
      
      #renderPrint(cat(input$stats_input_file))
    })
    
    observeEvent(input$action_seg_run, {
      if (is.null(x4$seg_file)) {
        message("Require template coordinate file first, please select and press")
        message("'Read file with coordinates' button")
        showNotification(
          "Require template coordinate file first, please select and press 'Read file with coordinates' button",
          type = "error"
        )
        
        
        return(NULL)
      }
      
      #start timer
      ptm<-proc.time()
      
      withProgress(message = "setting up coordinates and extracting pixels from raw data",
                   detail = "using existing coordinates",
                   value = 0.2,
                   {
                     req(x4$seg_file)
                     
                     message("setting up coordinates and extracting pixels from raw data")
                     
                     
                     
                     x1 = setup_values()[["x1"]] # bring in raw_files list
                     
                     setCardinalBPPARAM(par_mode())
                     
                     
                     
                     #extract only runs available in raw list
                     
                     names(x1$raw_list) <- gsub(".imzML", "", names(x1$raw_list))
                     
                     seg_file_trimmed<- Cardinal::subsetPixels(x4$seg_file, run %in% names(x1$raw_list))
                     x4$seg_file_trimmed<-seg_file_trimmed
                     
                     #reorder raw list to match file names in segmented file
                     raw_list_ord <-
                       x1$raw_list[runNames(seg_file_trimmed)]
                     
                     # seg_file_ord<-lapply(1:length(runNames(x4$seg_file)),
                     #                      function(x) x4$seg_file[,run(x4$seg_file)%in%names(x1$raw_list)[x]])
                     #
                     # seg_file_ord<-combine(seg_file_ord)
                     #
                     
                     seg_file_ord <- seg_file_trimmed
                     
                     #get coordinates from coordinate dataset
                     coord_list_segmented <-
                       lapply(runNames(seg_file_ord), function(x)
                         coord(seg_file_ord)[Cardinal::run(seg_file_ord) %in% x,])
                     #sum(unlist(lapply(1:16, function(x) dim(coord_list_segmented[[x]])[1])))  # check how many pixels selected
                     names(coord_list_segmented) <-
                       runNames(seg_file_ord)
                     
                     #if only using a subset of the segmented coordinates, reduce to those only in the coordinate set
                     coord_list_segmented <- coord_list_segmented[names(coord_list_segmented) %in% names(raw_list_ord)]
                     
                     
                     
                     if (length(coord_list_segmented) != length(x1$raw_list)) {
                       print(
                         "# of runs in raw data files list and segmented file does not match, please check. Trying to proceed anyway."
                       )
                       
                       if (length(coord_list_segmented) > length(x1$raw_list)) {
                         message(
                           "Mismatch betwen # of runs in coordinate file and # of raw files open. Will try to proceed with open raw file(s). Please check results carefully."
                         )
                         showNotification(
                           "Mismatch betwen # of runs in coordinate file and # of raw files open. Will try to proceed with open raw file(s). Please check results carefully.",
                           type = "warning"
                         )
                       
                        
                         
                         coord_list_segmented <-
                           coord_list_segmented[!is.na(names(raw_list_ord))]
                         raw_list_ord <-
                           x1$raw_list[runNames(seg_file_trimmed)]
                         raw_list_ord <-
                           raw_list_ord[sapply(raw_list_ord, function(i)
                             ! is.null(i))]
                         seg_file_ord <-
                           seg_file_ord %>% subsetPixels(run %in% names(raw_list_ord))
                         
                       } else {
                         message(
                           "Mismatch betwen # of runs in coordinate file and # of raw files open. Will try to proceed by reducing # of open raw file(s). Please check results carefully."
                         )
                         showNotification(
                           "Mismatch betwen # of runs in coordinate file and # of raw files open. Will try to proceed by reducing # of open raw file(s). Please check results carefully.",
                           type = "warning"
                         )
                         
                         #if using a subset of the data, does the NA / NULL entries int he list cause an issue?
                         x1$raw_list <-
                           x1$raw_list[names(coord_list_segmented)]
                         
                         
                         #return()
                         
                       }
                       
                       #return()
                     }
                     
                     #if Untargeted is not selected first, this value is not set, but usually will be 100% anyway
                     if(is.null(input$pix_for_peak_picking2)) {
                       select_ratio=100
                     } else {
                       select_ratio=input$pix_for_peak_picking2
                     }
                     
                     #get coordinates from downsampling
                     coord_list_reduced <-
                       lapply(1:length(coord_list_segmented), function(x)
                         coord_list_segmented[[x]][sample(
                           1:nrow(coord_list_segmented[[x]]),
                           round(
                             nrow(coord_list_segmented[[x]]) *select_ratio / 100
                           )
                         ),])
                     
                     names(coord_list_reduced)<-names(coord_list_segmented)
                     
                     
                     #get coordinates
                     
                     
                     px_expected <-
                       sum(unlist(lapply(1:length(coord_list_segmented), function(x)
                         dim(coord_list_reduced[[x]])[1]))) # check how many pixels selected
                     
                     #function to use raw image list and coordinates to create reduced raw set
                     
                     select_pix <-
                       function(index, raw_list, coord_set) {
                         raw_list[[index]][!is.na(prodlim::row.match(as.data.frame(coord(
                           raw_list[[index]]
                         )), as.data.frame(coord_set[[index]])))]
                       }
                     
                     
                     #most likely point of failure here....
                     
                     if (length(x1$raw_list) == 1) {
                       test_raw_reduced <-
                         select_pix(1, raw_list_ord[runNames(x1$raw_list[[1]])], coord_list_reduced)
                       tmp_seg_coord_list<-list(test_raw_reduced)
                     } else if (length(x1$raw_list) > 1) {
                       tmp_seg_coord_list<-try( bplapply(1:length(raw_list_ord), function(x)
                         select_pix(x, raw_list_ord, coord_list_reduced)))
                       
                       if(class(tmp_seg_coord_list) %in% "try-error") {
                         message("setting coordinate list failed, check names")
                         showNotification("setting coordinate list failed, check names")
                         return()
                       }
                       
                       i=1
                       
                         tmp_names<-unlist(lapply(tmp_seg_coord_list, runNames))
                         names(tmp_seg_coord_list) <- tmp_names
                         test_raw_reduced<-convertMSImagingExperiment2Arrays(tmp_seg_coord_list[[1]])

                         
                         while(i< (length(tmp_seg_coord_list))){
                           k=i+1
                           tmp<-convertMSImagingExperiment2Arrays(tmp_seg_coord_list[[k]])
                           test_raw_reduced<-try(combine(test_raw_reduced, tmp))
                           
                           if(class(test_raw_reduced) %in% "try-error"){
                             message("combining imagesets failed, check data files")
                             showNotification("combining imagesets failed, check data files", type="error")
                             return()
                           }
                           
                           i=i+1
                         }
                       
                      
                     }
                     
                     if (is.null(test_raw_reduced)) {
                       message(
                         "Cannot map coordinates to raw data with multiple files. Perhaps try individual runs."
                       )
                       showNotification(
                         "Cannot map coordinates to raw data with multiple files. Perhaps try individual runs.",
                         type = "error"
                       )
                       return(NULL)
                     }
                     #px_found<-unlist(lapply(1:length(runNames(test_raw_reduced)), function(x) dim(test_raw_reduced[,run(test_raw_reduced)%in%runNames(test_raw_reduced)[x]])[2]))
                     #raw_reduced<-combine(bplapply(1:length(raw_list), function(x) select_pix(x,raw_list, coord_list_reduced)))
                     #browser()
                     if (length(test_raw_reduced) != px_expected) {
                       print("reduced coordinate set does not match input coordinate set")
                       print("Run names from reduced coordinate set:")
                       print(runNames(test_raw_reduced))
                       print("")
                       print("Run names from raw data list:")
                       print(names(x1$raw_list))
                       return()
                     } else {
                       
                       test_raw_reduced<-convertMSImagingArrays2Experiment(test_raw_reduced)
                       pData(test_raw_reduced)<-pData(seg_file_ord)
                       
                      }
                     
                       
                       
                       
                     
                     
                     incProgress(amount = 0.2, message = "Coordinates mapped, now binning / pick picking")
                     #create final peak picked file
                     #browser()
                     if (input$targeted_pp == "untargeted") {
                       #setCardinalBPPARAM(SerialParam())
                       
                       
                       print("performing untargeted analysis")
                       
                       browser()
                       
                       x4$seg_pp_file <-
                         
                         
                         
                         
                         a<-
                         try(HTS_reproc(
                           test_raw_reduced,
                           SN = input$SNR2,
                           res = input$res,
                           align_tol = input$tol2,
                           method = input$pp_method2,
                           freq.min = input$freq_min2
                         ))
                       
                       
                       
                       # mse_queue<- test_raw_reduced |>
                       #   normalize() |>
                       #   #smooth() |>
                       #   #reduceBaseline() |>
                       #   peakPick(SNR=input$SNR2, method=input$method2, type="area", tolerance=NA, units="ppm")
                       # 
                       # #plot the middle spectrum?
                       # print(plot(mse_queue, i=round(dim(coord(msa))[1]/2,0), linewidth=2))
                       # 
                       # b<-try(
                       #   peakAlign(mse_queue, tolerance= input$tol2, units="ppm") %>%
                       #     subsetFeatures( freq > input$freq_min2)%>% 
                       #     summarizeFeatures()
                       # )
                       # 
                       if(class(x4$seg_pp_file) %in% "try-error") {
                         print("peak picking failed, check input files and parameters")
                         showNotification("Peak picking failed, check input files and parameters.", type = "error")
                         return()
                       }
                       
                       #report number of peaks picked
                       print(paste0("Number of peaks picked: ", ncol(x4$seg_pp_file)))
                       
                       
                     } else if (input$targeted_pp == "mean") {
                       #setCardinalBPPARAM(SerialParam())
                       
                       
                       print("using mean spectrum to create peak list and bin raw data")
                       
                       
                       
                       
                       test_mz_mean <- try(
                         #for debugging
                         #Cardinal::combine(x1$raw_list[c(1:3,5:8)]) %>%
                         # convertMSImagingExperiment2Arrays(test_raw_reduced) %>%
                         #   convertMSImagingArrays2Experiment(mass.range=c(setup_values()[["mz_max"]], setup_values()[["mz_min"]])) %>%
                         
                         
                         Cardinal::combine(lapply(tmp_seg_coord_list, convertMSImagingExperiment2Arrays)) %>%
                           convertMSImagingArrays2Experiment(mass.range=c(setup_values()[["mz_max"]], setup_values()[["mz_min"]])) %>%
                           estimateReferencePeaks( SNR=input$SNR2, 
                                                   method=input$pp_method2)
                       )
                       
                       #check class of test_mz_mean
                       if(class(test_mz_mean) %in% "try-error") {
                         print("mean spectrum peak picking failed, check input files")
                         showNotification("Mean spectrum peak picking failed, check input files.", type = "error")
                         return()
                       }
                       
                       incProgress(amount = 0.2, message = "Finished calculating mean spectrum, binning raw file(s)")
                       
                       
                       #test_mz_reduced2<-try(peakProcess(test_raw_reduced, 
                       test_mz_reduced2<-try(peakProcess(Cardinal::combine(lapply(tmp_seg_coord_list, convertMSImagingExperiment2Arrays)), 
                                                        ref=mz(test_mz_mean),
                                                        SN=input$SNR2,
                                                        type="area",
                                                        tolerance=input$tol2, units="ppm") %>% process() %>% summarizeFeatures()
                       )
                       
                       
                       # test_mz_reduced2 <- test_raw_reduced %>%
                       #   summarizeFeatures(FUN = "mean")  %>%
                       #   normalize(method = "tic") %>%
                       #   peakPick(method = input$pp_method2, SNR = input$SNR2) %>%
                       #   peakAlign(ref = "mean",
                       #             tolerance = input$tol2,
                       #             units = "ppm") %>%
                       #   peakFilter(freq.min = input$freq_min2) %>%
                       #   process()
                       # # 
                       # saveRDS(test_mz_reduced2, file = "mean_picked_depth_mz.rds")
                       # 
                       # seg_pp_file <-
                       #   ref_to_peaks(
                       #     test_raw_reduced,
                       #     mz_ref = mz(test_mz_reduced2),
                       #     tol = input$tol2
                       #   )
                       # 
                       x4$seg_pp_file <- test_mz_reduced2
                       
                       print(paste0("Number of peaks picked: ", ncol(x4$seg_pp_file)))
                       
                     } else if (input$targeted_pp == "targeted") {
                       print("performing targeted analysis")
                       print(paste0("Tolerance set to ", input$tol2))
                       file <- input$masses_file
                       ext <- tools::file_ext(file$datapath)
                       
                       req(file)
                       validate(need(ext == "txt", "Please upload a tab delimited .txt file"))
                       
                       neg_masses <-
                         readr::read_delim(file = file$datapath, delim = "\t")
                       
                       
                       
                       #check of exact.mass column
                       if (is.null((neg_masses$exact.mass))) {
                         message("No 'exact.mass' column, exiting. Check your mass list file")
                         showNotification("No 'exact.mass' column, exiting. Check your mass list file",
                                          type = "error")
                         return()
                       }
                       
                       neg_masses$exact.mass <-
                         as.numeric(neg_masses$exact.mass)
                       neg_masses <-
                         neg_masses[order(neg_masses$exact.mass),]
                       
                       neg_ref_mz <-
                         (na.omit(neg_masses$exact.mass))
                       
                       if (length(unique(neg_ref_mz)) != length(neg_ref_mz)) {
                         
                         
                         
                         print("duplicate mz values, please check mass list!!")
                         showNotification("duplicate mz values, please check mass list!! See console for details", type="error")
                         print(names(table(neg_ref_mz))[table(neg_ref_mz) > 1])
                         return()
                       }
                       
                       if (is.null(neg_ref_mz)) {
                         print("No masses from list, exiting")
                         return()
                       }
                       
                       # x4$seg_pp_file <- test_raw_reduced %>%
                       #   normalize(method = "tic") %>%
                       #   peakBin(ref = neg_ref_mz,
                       #           tolerance = input$tol2,
                       #           units = "ppm") %>% process()
                       
                       #browser()
                       
                       tmp.img<-try(peakProcess(test_raw_reduced, 
                                                        ref=neg_ref_mz,
                                                        #SN=input$SNR,
                                                        type="area",
                                                        tolerance=setup_values()[["tol"]], units="ppm") %>% process() %>% summarizeFeatures()
                       )
                       
                       x4$seg_pp_file <-tmp.img
                       
                       #report number of masses binned
                       print(paste0("Number of masses binned: ", ncol(x4$seg_pp_file)))
                       
                     }
                     #add mass annotations?
                     if (input$targeted_pp == "targeted") {
                       featureData(x4$seg_pp_file)$ID <- neg_masses$ID
                     }
                     
                    
                     #if using downsampling, re-bin
                     if (input$pix_for_peak_picking2 < 100) {
                       #combine reduced images.
                       test_raw_reduced <-
                         combine_card(bplapply(1:length(raw_list_ord), function(x)
                           select_pix(x, raw_list_ord, coord_list_segmented)))
                       #raw_reduced<-combine(bplapply(1:length(raw_list), function(x) select_pix(x,raw_list, coord_list_reduced)))
                       
                       
                       test_mz_reduced <- x4$seg_pp_file
                       
                       x4$seg_pp_file <-
                         ref_to_peaks(
                           test_raw_reduced,
                           mz_ref = mz(test_mz_reduced),
                           tol = input$tol2
                         )
                     }
                     #restore pdata
                     #browser()
                     #function to match original pdata to current data coordinates
                     #needed in the event of changes to coordinates externally
                     #pdat is original pData, msddf is the newly created data 
                     pdat_match<-function(pdat, msddf){
                       #msddf<-x4$seg_pp_file
                       #pdat<-pData(seg_file_ord)
                       
                       
                       df1<-coord(pdat)
                       df2<-coord(msddf)
                       
                       df1$run<-run(pdat)
                       df2$run<-run(msddf)
                       
                       
                       # Add a row number to df2 to preserve the original order
                       df2$index <- 1:nrow(df2)
                       
                       # Merge df1 with df2 on 'x' and 'y'
                       merged_df <- merge(df1, df2, by = c("run", "x", "y"))
                       
                       # Sort df1 based on the order of df2
                       ord_vec<-order(match(paste(df1$run, df1$x, df1$y), paste(df2$run, df2$x, df2$y)))
                       
                       sorted_pdat1 <- pdat[ord_vec, ]
                       
                       return(sorted_pdat1)
                       
                     }
                     
                     
                     ordered_pdat<-pdat_match(pdat=pData(seg_file_ord), msddf=x4$seg_pp_file)
                     
                     pData(x4$seg_pp_file) <- ordered_pdat
                     
                     plot_card_server("card_plot_depth", overview_peaks_sel = x4$seg_pp_file)
                     
                   })
      
      print(proc.time()-ptm)
      
    })
    
    
    
    
    
    #read segmented file and plot
    output$plot9 <- renderImage({
      req(x4$seg_pp_file)
      if (is.null(x4$seg_pp_file))
        return()
      
      # A temp file to save the output.
      # This file will be removed later by renderImage
      outfile <- tempfile(fileext = '.png')
      
      png(outfile, width = 800, height = 600)
      
      img.dat <- x4$seg_pp_file
      
      print(
        Cardinal::image(
          x4$seg_pp_file,
          colorscale = pals::parula(255),
          contrast.enhance = "histogram"
        )
      )
      
      
      
      #print(Cardinal::image(mytable_selected(), mz=ion, plusminus=input$plusminus_viz))
      
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
    
    
    
    
    observeEvent(input$save_imzml, {
      
      req(input$save_imzml)
      volumes <- c(wd = setup_values()[["wd"]], home = fs::path_home())
      shinyFiles::shinyFileSave(input, "save_imzml", roots = volumes, session = session)
      
      
      save_path <- shinyFiles::parseSavePath(volumes, input$save_imzml)
      if (nrow(save_path) == 0) return(NULL)
      
      
      filen <- as.character(save_path$datapath)
      
        pk_img <- x4$seg_pp_file
        #pData(pk_img) <- pData(x4$seg_file_trimmed)
        writeImzML(pk_img, filen)
        saveRDS(pk_img, paste0(filen,".rds"))
      }
    )
    
    # output$save_state2 <- downloadHandler(
    #   filename = function() {
    #     paste0(getwd(),
    #            "/MSI-depth_picked_SN-peak_picked-",
    #            Sys.Date(),
    #            ".rds")
    #   },
    #   content = function(filen) {
    #     #
    #     pk_img <- x4$seg_pp_file
    #     #pData(pk_img) <- pData(x4$seg_file_trimmed)
    #     saveRDS(pk_img, filen)
    #   }
    #) 
  })
}
