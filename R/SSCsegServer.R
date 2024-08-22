### R/SSCsegServer.R
SSCsegServer <- function(id, setup_values, preproc_values) {
  moduleServer(id, function(input, output, session) {
    # output$table <- DT::renderDataTable(mtcars)
    
    ns = session$ns #for dyanamic variable namespace
    
    #import parallel mode
    
    par_mode = reactive({
      setup_values()[["par_mode"]]
    })
    
    
    #create table as reactive
    run_table<-reactive({
      
      x2 <- preproc_values()[["x2"]]
      
      req(x2$overview_peaks_sel)
      size_list = unlist(lapply(runNames(x2$overview_peaks_sel), function(x)
        ncol(subsetPixels(x2$overview_peaks_sel, run = x))))
      
      tab<-cbind(
        run = runNames(x2$overview_peaks_sel),
        size = size_list
      )
      
      return(tab)
      
    })
    
    #list runs for selection
    output$mytable2 = DT::renderDataTable({
      
      DT::datatable(
        run_table(),
        #FIX SIZE
        selection = 'none',
        caption = "Choose runs to process for ssc clustering",
        extensions = c("Buttons", "Select"),
        options = list(
          dom = 'Bfrtip',
          select = TRUE,
          buttons = list('pageLength', "copy", "selectNone", "selectAll")
        )
      )
    },
    server = FALSE)
    
    observe({
      x2 <- preproc_values()[["x2"]]
      
      req(x2$ssc)
      #check if stored coordinates match current dataset; not rigorous, could check names in future
      if (!setequal(dim(x2$mytable_selected), c(req(nrow(coord(x2$ssc[[1]]))), req(nrow(x2$ssc[[1]]@featureData))))) {
        print("stored ssc does not match dimensions of currently selected dataset.")
        print("Manually remove stored ssc file and try again.")
        x2$ssc <- NULL
      }
    })
    
    #extract selected datasets from table when row is selected
    observe({
      x2 <- preproc_values()[["x2"]]
      req(x2$overview_peaks_sel)
      req(input$mytable2_rows_selected)
      
     
      
      #note this requires the order of list_proc_img and the table to be the same
      if (!is.null(x2$list_proc_img)) {
        
        #create temp list for reordering
        tmp_list <- x2$list_proc_img[input$mytable2_rows_selected]
        
        x2$mytable_selected <- combine_card(tmp_list)
        x2$tf_list <- rep(TRUE, ncol(x2$mytable_selected)) #reset tf_list based on new data
          
      } else {
        ids <- input$mytable2_rows_selected
        x2$mytable_selected <-
          x2$overview_peaks_sel %>% subsetPixels(
            Cardinal::run(x2$overview_peaks_sel) %in% runNames(x2$overview_peaks_sel)[ids]
          )
      }
    })
    
    observeEvent(input$fix_pix, {
      
      
      x2 <- preproc_values()[["x2"]]
      req(x2$mytable_selected)
      tmp_dat <- fix_pix(
        dat=x2$mytable_selected,
        remove = input$fix_pix_t_f,
        r = input$fix_pix_r,
        n_thresh = input$fix_pix_n_thresh
      )
      if (!is.null(tmp_dat)) {
        x2$mytable_selected <- tmp_dat
        
      }
      print("done fix_pix")
      
    })
    
    output$fix_pix_opts <- renderUI({
      req(input$fix_pix_tf)
      list(
        numericInput(ns("fix_pix_r"), "radius to search for neighbors", 2),
        numericInput(ns("fix_pix_n_thresh"), "min # of neighbors allowed", 3),
        checkboxInput(ns("fix_pix_t_f"), "Actually remove pixels?", value =
                        TRUE),
        actionButton(ns("fix_pix"), HTML("Remove isolated pixels"))
      )
      
    })
    

    
    #Run SSC
    observeEvent(input$action_ssc, {
      x2 <- preproc_values()[["x2"]]
      
      
      
      if (is.null(input$mytable2_rows_selected)) {
        
        showNotification("Select a dataset from table first", type="error", duration=10)
        print("Select a dataset from table first")
        return()
      }
      
      withProgress(message = "Performing Shrunken Spatial Centroids Analsysis", {
        #browser()
        print("starting SSC")
        setCardinalBPPARAM(par_mode())
        
        r = eval(parse(text = input$card_r))
        k = eval(parse(text = input$card_k))
        s_var = eval(parse(text = input$card_s))
        
        ssc_filename = paste(
          "ssc_out_r_",
          paste(r, collapse = ","),
          "_k_",
          paste(k, collapse = ","),
          "_s_",
          paste(s_var, collapse = ","),
          "_runs_",
          paste(input$mytable2_rows_selected, collapse = ",") ,
          ".RData",
          sep = ""
        )
        
        
        if (!file.exists(file = ssc_filename)) {
          #if file doesn't already exist
          
          dat_ssc <-
            try(spatialShrunkenCentroids(
              x2$mytable_selected,
              r = r,
              k = k,
              s = s_var,
              weights = input$ssc_method
            ))
          
          if (class(dat_ssc) == "try-error") {
            message("Cannot run ssc on this dataset. Try removing stray pixels with fix_pix first.")
            
            showNotification(
              "Cannot run ssc on this dataset. Try removing stray pixels with fix_pix first.",
              type = "error",
              duration = 10
            )
            return()
          } else {
            x2$ssc <- dat_ssc
            saveRDS(dat_ssc, file = ssc_filename)
          }
        } else {
          print("Restoring from existing ssc run")
          x2$ssc <- readRDS(file = ssc_filename)
        }
        
      })
      
      
    })
    
    
    output$ssc_model_choices <-
      renderUI({
        #https://gist.github.com/wch/4211337
        x2 <- preproc_values()[["x2"]]
        req(x2$ssc)
        
        
        ssc_models = ResultsList(x2$ssc)

        x2$ssc_models <- names(ssc_models)

        
        # Create the checkboxes for the ssc model to be working with.
        radioButtons(
          ns("ssc_model"),
          "Choose model",
          choices  = x2$ssc_models,
          selected = "none"
        )
      })
    
    observe({
      x2 <- preproc_values()[["x2"]]
      req(x2$ssc)
      req(input$ssc_model)
      
      #create tf_list based on selected model and colors
      
      if (sum(x2$ssc_models %in% input$ssc_model) == 1) {
        model_num <- which(x2$ssc_models %in% input$ssc_model)
        x2$bkcols <- x2$ssc[[model_num]]$class
        x2$tf_list <- x2$bkcols %in% x2$bkcols
      } else {
        print("input ssc model does not match any stored ssc models. Exiting")
        return(NULL)
      }
      
    })
    
    #Plot selected dataset
    observe( {
      x2 <- preproc_values()[["x2"]]
      req(x2$mytable_selected)
      
      if (dim(x2$mytable_selected)[2] > 0 &&
          dim(x2$mytable_selected)[2] == length(x2$tf_list)) {
        img.dat <- x2$mytable_selected %>% subsetPixels(x2$tf_list)
        #browser()
        if (dim(img.dat)[2] > 0) {
          plot_card_server("card_plot_ssc",
                           overview_peaks_sel = img.dat,
                           spatialOnly = TRUE)
        }
      }
      #browser()
      
    })
    
    
    output$Color_choices2 <-
      renderUI({
        #https://gist.github.com/wch/4211337
        # If missing input, return to avoid error later in function
        x2 <- preproc_values()[["x2"]]
        req(x2$ssc)
        # Get the data set with the appropriate name
        dat <- (x2$bkcols)
        cols <- unique(dat)
        
        # Create the checkboxes and select them all by default
        checkboxGroupInput(ns("cols2"),
                           "Choose colors",
                           choices  = cols,
                           selected = cols)
      })
    
    #change T/F list when something unchecked...
    observeEvent(input$cols2, {
      x2 <- preproc_values()[["x2"]]
      req(x2$ssc)
      
      x2$tf_list <- x2$bkcols %in% input$cols2
      
    })
    
    
    
    output$plot6a <- renderImage({
      x2 <- preproc_values()[["x2"]]
      req(x2$ssc)
      req(x2$tf_list)
      req(input$ssc_model)
      req(input$cols2)
      # A temp file to save the output.
      # This file will be removed later by renderImage
      outfile <- tempfile(fileext = '.png')
      
      png(outfile, width = 800, height = 600)
      
      
      cols <- x2$bkcols[x2$tf_list]
      
      
      updateSelectizeInput(session,
                           'Color_choices2',
                           choices = cols,
                           server = TRUE)
      #updateSelectizeInput(session, 'cols2', choices = cols, server = TRUE)
      
      
      img.dat <- x2$mytable_selected %>% subsetPixels(x2$tf_list)
      
      
      if (input$ssc_cols == "alphabet") {
        img.dat$cols<-cols
        p1<-(image(
          img.dat,
          "cols",
          key = T,
          #col = pals::alphabet2(),
        ))
      } else{
        showNotification("no alternate colors")
      }
      
      #browser()
      model_num <- which(x2$ssc_models %in% input$ssc_model)
      
      #browser()
      
      
      #p2<-
      p2<-  plot(x2$ssc, lwd=2, i=model_num, column=input$cols2, main="Weighted Spectrum")
      p3<-plot(x2$ssc, i=model_num, values="statistic", lwd=2, column=input$cols2, main="Statistic")
      
      plot_list<-list(p1,p2,p3)
      
      print(matter::as_facets(plot_list, free="xy"))
      

        #
      
      
      #
      #
      # }
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
    
    observeEvent(input$store_proc2, {
      x2 <- preproc_values()[["x2"]]
      req(x2$tf_list)
      
      #browser()
      
      tmp.img <- x2$mytable_selected %>% subsetPixels(x2$tf_list)
      
      #assign classes to trimmed data
      tmp.img$ssc_cols <- x2$ssc[[input$ssc_model]]$class[x2$tf_list]
      
      all_runs <-
        runNames(x2$overview_peaks_sel) #from the total set being analyzed
      runs <- runNames(tmp.img)
      idx <-
        which(all_runs %in% runs)  #are we assuming that the runs are stored in the same order?
      
      #if x2$list_proc_img is null, create it
      if (is.null(x2$list_proc_img)) {
        #put all runs into a separate list from the originally imported / peakpicked data
        for (i in 1:length(all_runs)) {
          x2$list_proc_img[[i]] <-
            x2$overview_peaks_sel %>% subsetPixels(Cardinal::run(x2$overview_peaks_sel) ==
                                                     all_runs[i])
          #names(x2$list_proc_img[[i]]) <- all_runs[i]
        }
      }
      
      
      
      #replace list elements with updated processed data
      for (i in runs) {
        
        #get names of proc data from proc list
        runnames_all<-unlist(lapply(x2$list_proc_img, runNames))
        
        index<-which(runnames_all %in% i)
        
        x2$list_proc_img[[index]] <-
          tmp.img %>% subsetPixels(Cardinal::run(tmp.img) == i)
        
      }
      
      #ensure the same pData columns
      #find all pData elements across all datasets
      all_pData <-
        unique(unlist(lapply(x2$list_proc_img, function(x)
          colnames(pData(x)))))
      
      #function to check pdata and ensure all good
      chk_pData <- function(dat, pDat_cols = all_pData) {
        col_diff <- setdiff(pDat_cols, colnames(pData(dat)))
        
        if (length(col_diff) == 0)
          return(dat)
        else {
          pData(dat)[col_diff] <- NA
        }
        
        return(dat)
        
      }
      
      x2$list_proc_img <-
        lapply(x2$list_proc_img, function(x)
          chk_pData(x))
      
      #browser()
      #restore the data as a complete imageset
      x2$mytable_selected <-
        combine_card(x2$list_proc_img[input$mytable2_rows_selected])
      
      
      #create tmp data for storage of processed data
      tmp_dat <- x2$list_proc_img  
      
      #reorder data based on ncol size
      #tmp_dat <- tmp_dat[order(unlist(lapply(tmp_dat, ncol)), decreasing = TRUE)]
      
      tmp_dat.img<-combine_card(tmp_dat)
      
      
      #x2$overview_peaks_sel <- tmp_dat.img
      #ssc_tmp_img() <- tmp_dat.img
      
      x2$tf_list <- rep(TRUE, ncol(tmp_dat.img)) #reset tf_list based on new data
      x2$data_list <- NULL
      x2$ssc <- NULL
      #x2$list_proc_img<-NULL
      
      #browser()
      
    })
    
    #create reactive variable for storing processed data
    ssc_tmp_img <- reactive({
      x2<-preproc_values()[["x2"]] 
      return(x2$overview_peaks_sel)
    })
    
    
    
    observeEvent(input$save_imzml, {
      
      req(input$save_imzml)
      if(is.null(x2$list_proc_img)) {
        showNotification("No processed data to save- do you need to store the data?", type="error", duration=10)
        message( "No processed data to save- do you need to store the data?")
          return()
      }
      
      volumes <- c(wd = setup_values()[["wd"]], home = fs::path_home())
      shinyFiles::shinyFileSave(input, "save_imzml", roots = volumes, session = session)
      
      
      save_path <- shinyFiles::parseSavePath(volumes, input$save_imzml)
      if (nrow(save_path) == 0) return(NULL)
      
      browser()
      filen <- as.character(save_path$datapath)
        x2 <- preproc_values()[["x2"]]
        #pk_img <- x0$overview_peaks
        
        #create tmp data for storage of processed data
        tmp_dat <- x2$list_proc_img  
        
        #reorder data based on ncol size
        tmp_dat <- tmp_dat[order(unlist(lapply(tmp_dat, ncol)), decreasing = TRUE)]
        
        tmp_dat.img<-combine_card(tmp_dat)
        
        pk_img <-  tmp_dat.img
        
        a<-fData(pk_img)[unique(colnames(fData(pk_img)))]
        
        fData(pk_img) <- MassDataFrame(mz=a$mz, a %>% as.data.frame() %>% dplyr::select(-mz))
        
        writeImzML(pk_img, filen)
      }
    )
  })
}
