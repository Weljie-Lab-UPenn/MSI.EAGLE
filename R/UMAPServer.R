### R/UMAPServer.R
UMAPServer <- function(id, setup_values, preproc_values) {
  moduleServer(id, function(input, output, session) {
    #for dynamic UI
    ns = session$ns
    
    #import parallel mode
    par_mode = reactive({
      setup_values()[["par_mode"]]
    })
    
    
    
    
    ###START UMAP for background removal
    
    
    
    
    #list runs for selection
    output$mytable = DT::renderDataTable({
      x0 <- preproc_values()[["x0"]]
      req(x0$overview_peaks)
      ovps <- x0$overview_peaks
      #x2<-preproc_values()[["x2"]]
      #req(x2$overview_peaks_sel)
      #ovps<-x2$overview_peaks_sel
      DT::datatable(
        cbind(
          run = runNames(ovps),
          size = lapply(runNames(ovps), function(x)
            ncol(ovps[, Cardinal::run(ovps) %in% x]))
        ),
        #selection = list(mode = "multiple"),
        selection = 'none',
        caption = "Choose runs to process for UMAP clustering",
        extensions = c("Buttons", "Select"),
        options = list(
          dom = 'Bfrtip',
          select = TRUE,
          buttons = list('pageLength', "copy", "selectNone", "selectAll")
        )
      )
    },
    server = FALSE)
    
    #extract selected datasets from table (x2$mytable_selected)
    observeEvent(input$mytable_rows_selected, {
      x2 <- preproc_values()[["x2"]]
     
      if (is.null(x2$list_proc_img)) {
        #only do this if there is no table already created for anatomical segmentation
        req(x2$overview_peaks_sel)
        req(input$mytable_rows_selected)
        
        ids <- input$mytable_rows_selected
        x2$mytable_selected <-
          x2$overview_peaks_sel %>% subsetPixels(
            Cardinal::run(x2$overview_peaks_sel) %in% runNames(x2$overview_peaks_sel)[ids]
          )
      }
      #x2$tf_list<-coord(x2$mytable_selected) %in% coord(x2$mytable_selected)
      #browser()
      # Update tf_list accordingly TEMP, check carefully
      #x2$tf_list <- rep(TRUE, ncol(x2$mytable_selected))
    })
    
    #update the x2$mytable_selected variable if there is a change to the selected rows.
    #use processed data if available... need to think about tf list
    observeEvent(input$mytable_rows_selected, {
      
      x2 <- preproc_values()[["x2"]]
      req(x2$mytable_selected)
      #browser()
      if (input$seg_choice == "anat_seg") {
        if (!is.null(x2$list_proc_img)) {
          print("updating selected table from previously processed phenotype data")
          
          #store runNames from previous mytable_selected to help with tf_list tracking
          
          old_mytable_selected <- x2$mytable_selected
          x2$mytable_selected <-
            combine(x2$list_proc_img[input$mytable_rows_selected])
          
          #update tf list to reflect new dataset. create list of all tf
          
          tf_list_list <-
            lapply(x2$list_proc_img, function(x)
              rep(TRUE, ncol(x)))
          names(tf_list_list) <-
            unlist(lapply(x2$list_proc_img, runNames))
          
          x2$tf_list <-
            unlist(tf_list_list[runNames(x2$mytable_selected)])
          x2$pdat_anat <- (pData(x2$mytable_selected))
          
          
          
        }
      } else if (input$seg_choice == "bk_seg"){
        
          if (!is.null(x2$list_proc_img)) {
            print("updating selected table from previously processed phenotype data")
            
            
            x2$mytable_selected <-
              combine_card(x2$list_proc_img[input$mytable_rows_selected])
            
            x2$tf_list <-
              rep(TRUE, ncol(x2$mytable_selected))
          }
      }
    })
    
    
    #create reactive temporary variable to store dataset being worked on which matches the t/f list
    #start with the inital dataset selected from the table.
    tmp_umap_dat <- reactive({
      x2 <- preproc_values()[["x2"]]
      req(x2$mytable_selected)
      if (!is.null(x2$list_proc_img)) {
        img.dat <-
          combine_card(x2$list_proc_img[input$mytable_rows_selected])
          #do.call(cbind, x2$list_proc_img[input$mytable_rows_selected])
          #Cardinal::combine(x2$list_proc_img[input$mytable_rows_selected])
        #x2$tf_list<-rep(TRUE, ncol(img.dat))
        return(img.dat)
      } else {
        #if there is no processed data, default to the mytable selected original data
        
        x2$tf_list <- rep(TRUE, ncol(x2$mytable_selected))
        return(x2$mytable_selected)
        
      }
    })
    
    
    #create tf_list and tf_list_umap
    observe( {
      #true false list corresponsds to size of selected entries in table.
      x2 <- preproc_values()[["x2"]]
      req(x2$mytable_selected)
      
      #if there is a selection in the table, create a true/false list
      if (dim(x2$mytable_selected)[2] > 0) {
        if (is.null(x2$tf_list) | length(x2$tf_list) == 0) {
          x2$tf_list <-
            !is.na(prodlim::row.match(as.data.frame(coord(
              x2$mytable_selected
            )), as.data.frame(coord(
              x2$mytable_selected
            ))))
        }
        #x2$tf_list<-coord(x2$mytable_selected) %in% coord(x2$mytable_selected)
        x2$tf_list_umap <- NULL
      }
    })
    #
    #
    
    
    
    #create the plot using plot_card_server()
    observe( {
      req(input$mytable_rows_selected)
      x2 <- preproc_values()[["x2"]]
      req(x2$mytable_selected)
      #req(tmp_umap_dat())
      req(x2$tf_list)
      
      if (input$seg_choice == "anat_seg") {
        #look for data in processed image area first
        if (!is.null(x2$list_proc_img)) {
          img.dat <- combine(x2$list_proc_img[input$mytable_rows_selected])
        } else if (!is.null(x2$tf_list_umap) &&
                   !is.null(x2$umap_name)) {
          if (identical(x2$umap_name, runNames(x2$mytable_selected))) {
            img.dat <- x2$mytable_selected %>% subsetPixels(x2$tf_list_umap)
          }
        } else {
          img.dat <- x2$mytable_selected #%>% subsetPixels(x2$tf_list)
        }
      } else {
       # browser()
        
        if(sum(x2$tf_list)==ncol(x2$mytable_selected)) {
          img.dat <- x2$mytable_selected
        } else if (length(x2$tf_list) == ncol(x2$mytable_selected)) { 
          
          img.dat <- x2$mytable_selected %>% subsetPixels(x2$tf_list)
        } else {
          #browser()
          message("tf_list is not the same size as the dataset, exiting")
          showNotification("tf_list is not the same size as the dataset, exiting", type="error")
          return(NULL)
        }

      }
      
      
      
      if (sum(x2$tf_list) == 0) {
        print("nothing selected in tf_list, exiting")
        return(NULL)
        
      }
      
      #if there is already processed data, use that
      if (!is.null(x2$list_proc_img)) {
        img.dat <-
          combine_card(x2$list_proc_img[input$mytable_rows_selected])
      }
      
      plot_card_server("plot_card_umap",
                       overview_peaks_sel = img.dat,
                       spatialOnly = TRUE)
      
    })
    
    #create the plot when dataset updated
    observeEvent(req(tmp_umap_dat()), {
      plot_card_server("plot_card_umap",
                       overview_peaks_sel = tmp_umap_dat(),
                       spatialOnly = TRUE)
      
    })
    
      
    
    
    #run fix_pix to remove isolated pixels
    observeEvent(input$fix_pix, {
      x2 <- preproc_values()[["x2"]]
      req(x2$mytable_selected)
      tmp_dat <- fix_pix(
        x2$mytable_selected,
        remove = input$fix_pix_t_f,
        r = input$fix_pix_r,
        n_thresh = input$fix_pix_n_thresh
      )
      if (!is.null(tmp_dat)) {
        x2$mytable_selected <- tmp_dat
        
      }
      print("done fix_pix")
      
    })
    
    output$fm_params <- renderUI({
      if (input$fastmap == TRUE) {
        list(
          numericInput(ns("fm_r"), "spatial r for fastmap", value = 1),
          numericInput(
            ns("fm_ncomp"),
            "Number of fastmap components",
            value = 20,
            min = 3
          ),
          selectInput(
            ns("fm_method"),
            "Fastmap method",
            choices = c("gaussian", "adaptive")
          ),
          selectInput(
            ns("fm_metric"),
            "Fastmap metric",
            choices = c("average", "neighborhood", "correlation")
          )
        )
      }
      
    })
    
    
    
    
    #Run UMAP
    observeEvent(input$umap1, {
      if (is.null(input$mytable_rows_selected)) {
        showNotification("Select a dataset from table first", type="warning")
        print("Select a dataset from table first")
        return()
      }
      #observe(input$umap1, {
      print("starting UMAP")
      #reset colorscale
      x2 <- preproc_values()[["x2"]]
      
      
      ptm <- proc.time()
      
      withProgress(message = "Performing UMAP analysis and color clustering", value =
                     0.2, {
                       if (input$seg_choice == "anat_seg") {
                         img.dat <- x2$mytable_selected %>% subsetPixels(x2$tf_list)
                         x2$tf_list_anat <- rep(TRUE, ncol(img.dat))
                         x2$tf_list_umap <-
                           x2$tf_list  #store T/F list to store umap pixels
                         x2$umap_name <- runNames(img.dat)
                         x2$pdat_anat <- pData(x2$mytable_selected)
                         
                         
                       } else {
                         #img.dat <- x2$mytable_selected
                         img.dat <- tmp_umap_dat()
                         x2$rcol_plot <- NULL
                       }
                       
                       if (input$fastmap == TRUE) {
                         setCardinalBPPARAM(par_mode())
                         data_list <- try(get_umap(
                           img.dat,
                           outliers = F,
                           nn = input$nn,
                           search_k = input$search_k,
                           min_dist = input$min_dist,
                           n_trees = input$n_trees,
                           set_op_mix_ratio = input$set_op_mix_ratio,
                           pca = eval(parse(text = input$pca_umap)),
                           fastmap = TRUE,
                           fm_r = input$fm_r,
                           fm_method = input$fm_method,
                           fm_metric = input$fm_metric,
                           fm_ncomp = input$fm_ncomp
                         ))
                         if(class(data_list) %in% "try-error"){
                           showNotification("Fastmap umap failed, check parameters", type="error")
                           message("Fastmap umap failed, check parameters")
                           return()
                         }
                         
                       } else {
                         data_list <- get_umap(
                           img.dat = img.dat,
                           outliers = F,
                           #thresh=input$thresh,
                           search_k = input$search_k,
                           min_dist = input$min_dist,
                           set_op_mix_ratio = input$set_op_mix_ratio,
                           pca = eval(parse(text = input$pca_umap))
                         )
                       }
                       
                       
                       
                       if (is.null(data_list$umap_separation$umap_out)) {
                         print("UMAP matrix is null, exiting")
                         return()
                       }
                       
                       incProgress(amount = 0.5,
                                   message = "Clustering",
                                   detail = "Performing dbscan and hdbscan")
                       
                       print("starting dbscan")
                       data_list$dbscan_umap_separation <-
                         dbscan::dbscan(
                           as.matrix(data_list$umap_separation$umap_out),
                           eps = input$eps ,
                           minPts = input$minPts
                         )
                       print("starting hdbscan")
                       
                       
                       if (input$hdbscan_method == "R_hdbscan") {
                         try(data_list$hdbscan_umap_separation <-
                               dbscan::hdbscan(as.matrix(data_list$umap_separation$umap_out) ,
                                               minPts = input$minPts))
                         # } else if (input$hdbscan_method == "py_hdbscan") {
                         #   try(data_list$hdbscan_umap_separation$cluster <-
                         #         hdbscan_python(
                         #           as.matrix(data_list$umap_separation$umap_out) ,
                         #           min_cluster_size = input$minPts,
                         #           min_samples = input$minSamples,
                         #           epsilon = input$epsilon
                         #         ))
                       } else if (input$hdbscan_method == "hdbscan_none") {
                         data_list$hdbscan_umap_separation <- NULL
                         data_list$hdbscan_umap_separation$cluster <-
                           NULL
                         data_list$hdbscan_umap_separation$cluster <-
                           NA
                       }
                       
                       #replace numeric clustering output with text colors to simplify plotting
                       
                       
                       #fine R colors to match polychrome
                       #based on: https://stackoverflow.com/questions/41209395/from-hex-color-code-or-rgb-to-color-name-using-r
                       #define all r colors
                       r_colors <- data.frame(color = colors())
                       r_colors <-
                         cbind(r_colors, t(col2rgb(colors())))
                       
                       pc <-
                         as.data.frame(t(col2rgb(pals::polychrome())))
                       col_data <-
                         cbind(Color = pals::polychrome(), pc)
                       
                       pc_color_name <-
                         sapply(seq_along(col_data$Color),
                                function(i)
                                  r_colors$color[which.min(
                                    (r_colors$red - col_data$red[i]) ^ 2 +
                                      (r_colors$green - col_data$green[i]) ^ 2 +
                                      (r_colors$blue - col_data$blue[i]) ^ 2
                                  )])
                       
                       
                       #ensure that the max number of colors is not exceeded
                       
                       max_col <-
                         function(num_vec_in, color_vec = pc_color_name) {
                           a = num_vec_in
                           
                           if (length(table(a)) / length(color_vec) > 1) {
                             ff <- floor(length(table(a)) / length(color_vec)) + 1
                             col_vector = rep(as.vector(color_vec), ff)
                             print("Warning: more than 36 colors, color replication performed")
                           } else {
                             col_vector = as.vector(color_vec)
                           }
                           
                           names(col_vector) <-
                             c(0:(length(col_vector) - 1))
                           return(col_vector)
                         }
                       
                       dbscan_cols <-
                         max_col(data_list$dbscan_umap_separation$cluster,
                                 pc_color_name)
                       hdbscan_cols <-
                         max_col(data_list$hdbscan_umap_separation$cluster,
                                 pc_color_name)
                       
                       
                       data_list$dbscan_umap_separation$cluster <-
                         dbscan_cols[as.character(data_list$dbscan_umap_separation$cluster)]
                       data_list$hdbscan_umap_separation$cluster <-
                         hdbscan_cols[as.character(data_list$hdbscan_umap_separation$cluster)]
                       
                       
                       print("finished clustering")
                       
                       incProgress(amount = 0.9, detail = "Starting quantile color reduction")
                       
                       print("started color reduction")
                       
                       if (as.character(Sys.info()['sysname']) == "Windows") {
                         #if on windows, use a single core
                         old_nc = setup_values()[["ncores"]]
                         ncores = 1
                       }  else {
                         ncores = setup_values()[["ncores"]]
                       }
                       
                       #reduced<-unlist(future_lapply(data_list$umap_separation$color_scheme, function(x) color.id(x)[1]))
                       reduced <-
                         unlist(
                           parallel::mclapply(data_list$umap_separation$color_scheme, function(x)
                             color.id(x)[1], mc.cores = ncores)
                         )
                       cols_red_list <-
                         unique(names(table(reduced)[table(reduced) > quantile(table(reduced), input$quant)]))
                       #cols_red_90<-unlist(future_lapply(data_list$umap_separation$color_scheme, function(x) color.id(x, palette = cols_red_list)[1]))
                       cols_red_90 <-
                         unlist(
                           parallel::mclapply(data_list$umap_separation$color_scheme, function(x)
                             color.id(x, palette = cols_red_list)[1],  mc.cores = ncores)
                         )
                       print("finished color reduction")
                       
                       if (as.character(Sys.info()['sysname']) == "Windows") {
                         #if on windows, return original ncores
                         ncores = old_nc
                       }  else {
                         ncores = setup_values()[["ncores"]]
                       }
                       
                       
                       data_list$umap_separation$initial_reduced <-
                         reduced
                       data_list$umap_separation$col_reduced <-
                         cols_red_90
                       
                       
                       
                       x2$data_list <- data_list
                       
                     }
                   )
      print(proc.time()-ptm)
    })
    
    #     observe({
    #       x2$data_list <- data_list_shiny()
    #     })
    
    #run clustering only
    observeEvent(input$action_dbscan, {
      x2 <- preproc_values()[["x2"]]
      req(x2$data_list)
      #observe(input$action_dbscan, {
      
      data_list <- x2$data_list
      
      print("Re-running H/DBSCAN")
      
      withProgress(message = "Re-running clustering on existing UMAP", {
        data_list$dbscan_umap_separation <-
          dbscan::dbscan(
            as.matrix(data_list$umap_separation$umap_out),
            eps = input$eps ,
            minPts = input$minPts
          )
        
        if (input$hdbscan_method == "R_hdbscan") {
          try(data_list$hdbscan_umap_separation <-
                dbscan::hdbscan(as.matrix(data_list$umap_separation$umap_out) ,
                                minPts = input$minPts))
        } else if (input$hdbscan_method == "py_hdbscan") {
          try(data_list$hdbscan_umap_separation$cluster <-
                hdbscan_python(
                  as.matrix(data_list$umap_separation$umap_out) ,
                  min_cluster_size = input$minPts,
                  min_samples = input$minSize,
                  epsilon = input$epsilon
                ))
        } else if (input$hdbscan_method == "hdbscan_none") {
          data_list$hdbscan_umap_separation <- NULL
          data_list$hdbscan_umap_separation$cluster <- NULL
          data_list$hdbscan_umap_separation$cluster <- NA
        }
        
        if (as.character(Sys.info()['sysname']) == "Windows") {
          #if on windows, use a single core
          old_nc = setup_values()[["ncores"]]
          ncores = 1
        }  else {
          ncores = setup_values()[["ncores"]]
        }
        
        reduced <-
          unlist(
            parallel::mclapply(data_list$umap_separation$color_scheme, function(x)
              color.id(x)[1], mc.cores = ncores)
          )
        cols_red_list <-
          unique(names(table(reduced)[table(reduced) > quantile(table(reduced), input$quant)]))
        cols_red_90 <-
          unlist(
            parallel::mclapply(data_list$umap_separation$color_scheme, function(x)
              color.id(x, palette = cols_red_list)[1], mc.cores = ncores)
          )
        
        if (as.character(Sys.info()['sysname']) == "Windows") {
          #if on windows,restore value of ncores
          ncores = old_nc
        }
        
        data_list$umap_separation$initial_reduced <- reduced
        data_list$umap_separation$col_reduced <- cols_red_90
        
        
        #replace numeric clustering output with text colors to simplify plotting
        
        
        #find R colors to match polychrome
        #based on: https://stackoverflow.com/questions/41209395/from-hex-color-code-or-rgb-to-color-name-using-r
        #define all r colors
        r_colors <- data.frame(color = colors())
        r_colors <- cbind(r_colors, t(col2rgb(colors())))
        
        pc <- as.data.frame(t(col2rgb(pals::polychrome())))
        col_data <- cbind(Color = pals::polychrome(), pc)
        
        pc_color_name <- sapply(seq_along(col_data$Color),
                                function(i)
                                  r_colors$color[which.min(
                                    (r_colors$red - col_data$red[i]) ^ 2 +
                                      (r_colors$green - col_data$green[i]) ^ 2 +
                                      (r_colors$blue - col_data$blue[i]) ^ 2
                                  )])
        
        
        #ensure that the max number of colors is not exceeded
        
        max_col <- function(num_vec_in, color_vec = pc_color_name) {
          a = num_vec_in
          
          if (length(table(a)) / length(color_vec) > 1) {
            ff <- floor(length(table(a)) / length(color_vec)) + 1
            col_vector = rep(as.vector(color_vec), ff)
            print("Warning: more than 36 colors, color replication performed")
          } else {
            col_vector = as.vector(color_vec)
          }
          
          names(col_vector) <- c(0:(length(col_vector) - 1))
          return(col_vector)
        }
        
        dbscan_cols <-
          max_col(data_list$dbscan_umap_separation$cluster,
                  pc_color_name)
        hdbscan_cols <-
          max_col(data_list$hdbscan_umap_separation$cluster,
                  pc_color_name)
        
        
        data_list$dbscan_umap_separation$cluster <-
          dbscan_cols[as.character(data_list$dbscan_umap_separation$cluster)]
        data_list$hdbscan_umap_separation$cluster <-
          hdbscan_cols[as.character(data_list$hdbscan_umap_separation$cluster)]
        
        
        print("finished clustering")
        
        
        
        
        x2$data_list <- data_list
      })
    })
    
    
    
    #clust<- reactive({
    
    observe({
      x2 <- preproc_values()[["x2"]]
      if (is.null(x2$data_list)) {
        return()
      }
      x2$bkcols <- switch(
        input$umap_cols,
        #get color info
        #"RGB"= data_list_shiny()$umap_separation$color_scheme,
        "Reduced" = x2$data_list$umap_separation$initial_reduced,
        "Reduced2" = x2$data_list$umap_separation$col_reduced,
        "dbscan" = x2$data_list$dbscan_umap_separation$cluster,
        "hdbscan" = x2$data_list$hdbscan_umap_separation$cluster
      )
      
    })
    
    
    output$data_list_table <- renderPrint({
      x2 <- preproc_values()[["x2"]]
      paste(table(x2$bkcols), table(unique(x2$tf_list)))
    })
    
    
    #UMAP embeddings pairwise plot
    output$plot6 <- renderImage({
      x2 <- preproc_values()[["x2"]]
      req(x2$data_list)
      req(x2$tf_list)
      print("plot6 continuing")
      # A temp file to save the output.
      # This file will be removed later by renderImage
      outfile <- tempfile(fileext = '.png')
      
      png(outfile, width = 800, height = 600)
      
      
      if (input$seg_choice == "anat_seg") {
        tf_list = x2$tf_list_anat
      } else {
        tf_list <- x2$tf_list
      }
      
      if (length(tf_list) == 0)
        return(NULL)
      
      cols <- x2$bkcols[tf_list]
      cols[cols == 0] <- 100
      cols[cols %in% NA] <- 100
      
      updateSelectizeInput(session,
                           'Color_choices',
                           choices = cols,
                           server = TRUE)
      
      palette(pals::polychrome())
      print(
        pairs(
          x2$data_list$umap_separation$umap_out[tf_list,],
          col = cols,
          pch = ".",
          cex = input$cex,
          main = runNames(x2$overview_peaks_sel)
        )
      )
      palette("ggplot2")
      
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
    
    
    output$segmentation_plots <- renderUI({
      if (input$seg_choice == "bk_seg") {
        list(
          imageOutput(ns("plot7"), width = "800px", height = "600px"),
          imageOutput(ns("plot6"), width = "800px", height = "600px"),
          imageOutput(ns("colormaps"), width = "800px", height = "300px")
        )
      } else if (input$seg_choice == "anat_seg") {
        list(
          imageOutput(
            ns("plot7_tissue_umap"),
            width = "800px",
            height = "600px"
          ),
          imageOutput(ns("colormaps"), width = "800px", height = "300px"),
          imageOutput(ns("plot7"), width = "800px", height = "600px"),
          imageOutput(ns("plot6"), width = "800px", height = "600px")
          
        )
      }
    })
    
    
    
    output$plot7_tissue_umap <- renderImage({
      
      
      req(input$mytable_rows_selected)
      req(input$seg_pdat_col)
      x2 <- preproc_values()[["x2"]]
      
      x2$rcol_plot <- input$seg_pdat_col
      
      req(x2$rcol_plot)
      req(input$show_dat)
      
      print("plot7_tissue_umap_continuting")
      outfile <- tempfile(fileext = '.png')
      
      png(outfile, width = 800, height = 600)
      
      
      
      
      if (input$seg_choice == "anat_seg") {
        #test for existing proc data w/o umap, then test for umap to visualize
        if (!is.null(x2$list_proc_img) && is.null(x2$umap_name)) {
          img.dat <- combine(x2$list_proc_img[input$mytable_rows_selected])
        } else if (!is.null(x2$tf_list_umap) &&
                   !is.null(x2$umap_name)) {
          if (identical(x2$umap_name, runNames(x2$mytable_selected))) {
            img.dat <- x2$mytable_selected %>% subsetPixels(x2$tf_list_umap)
            if (!is.null(x2$pdat_anat)) {
              pData(img.dat) <- x2$pdat_anat[x2$tf_list_umap,]
            }
          }
        } else {
          #test to see if x2$pdat_anat matches and use that for viz
          img.dat <-
            x2$mytable_selected #%>% subsetPixels(x2$tf_list)
          if (!is.null(x2$pdat_anat) &&
              nrow(x2$pdat_anat) == ncol(img.dat)) {
            pData(img.dat) <- x2$pdat_anat
          }
        }
      } else {
        img.dat <- x2$mytable_selected #%>% subsetPixels(x2$tf_list)
      }
      
      
      
      
      #img.dat<-x2$mytable_selected #x2$overview_peaks_sel  #subset somehow??
      
      
      
      #if(x2$rcol_plot=="Rcol_reduced"){
      cols = as.data.frame(pData(img.dat))[, x2$rcol_plot]
      tf_list <- cols %in% input$show_dat
      
      img.dat <- img.dat %>% subsetPixels(tf_list)
      
      x <-
        coord(img.dat)[, 1]
      y <- coord(img.dat)[, 2]
      runs = Cardinal::run(img.dat)
      aa <- data.frame(x, y, runs)
      cols = (as.data.frame(pData(img.dat))[, x2$rcol_plot])
      
      # } else {
      #   
      #   cols=(as.factor(as.data.frame(pData(img.dat))[,x2$rcol_plot]))
      #   tf_list<-cols%in%input$show_dat
      #
      #   img.dat <- img.dat %>% subsetPixels(tf_list)
      #
      #   x<-coord(img.dat)[,1] ; y<-coord(img.dat)[,2]; runs=Cardinal::run(img.dat)
      #   aa<-data.frame(x,y,runs=req(runs))
      #   cols=as.numeric(as.factor(as.data.frame(pData(img.dat))[,x2$rcol_plot]))+1
      #
      #
      #
      # }
      
      
      
      if (dim(aa)[1] == 0) {
        return(NULL)
      }
      
      #test cols to make sure they will work. from here: https://stackoverflow.com/questions/13289009/check-if-character-string-is-a-valid-color-representation
      areColors <- function(x) {
        sapply(x, function(X) {
          tryCatch(is.matrix(col2rgb(X)), 
                   error = function(e) FALSE)
        })
      }
      
      if(sum(areColors(cols))<length(cols)){
        cat("\npData field for vizualization not recognized as colors, converting to alphabet")
        
        cols_col<-as.numeric(as.factor(cols))
        
        #max_col<-max(as.numeric(as.factor(cols)))
        
        #alph_cols<-palette("alphabet")[1:max_col]
        
        #cols_col<-factor(cols_col, levels=alph_cols)
        
        cols=cols_col
        
        
        p1 <- ggplot2::ggplot(aa, ggplot2::aes(x = x, y = y)) +
          ggplot2::geom_tile(ggplot2::aes(fill = factor(cols))) +
          ggplot2::facet_wrap(. ~ runs) +
          ggplot2::theme_minimal() +
          ggplot2::scale_y_continuous(trans = "reverse")+
          ggplot2::theme(legend.position = "bottom")+
          ggplot2::scale_fill_discrete(name=x2$rcol_plot, labels=levels(factor(as.data.frame(pData(img.dat))[, x2$rcol_plot])))
        
        
      } else {
        
      
        p1 <- ggplot2::ggplot(aa, ggplot2::aes(x = x, y = y)) +
          ggplot2::geom_tile(fill = cols) +
          ggplot2::facet_wrap(. ~ runs) +
          ggplot2::theme_minimal() +
          ggplot2::scale_y_continuous(trans = "reverse")
          
      }
      
      print(p1)
      
      dev.off()
      graphics.off()
      
      # Return a list containing the filename
      list(
        src = outfile,
        contentType = 'image/png',
        width = 800,
        height = 600,
        alt = "This is alternate text"
      )
    }, deleteFile = TRUE)
    
    
    
    #plot to make sure it looks okay
    output$plot7 <- renderImage({
      x2 <- preproc_values()[["x2"]]
      req(x2$mytable_selected)
      req(x2$tf_list)
      
      outfile <- tempfile(fileext = '.png')
      
      png(outfile, width = 800, height = 600)
      
      
      if (input$seg_choice == "anat_seg") {
        tf_list = x2$tf_list_anat
      } else {
        tf_list <- x2$tf_list
      }
      
      
      
      #require(ggplot2)
      
      req(x2$data_list)
      req(x2$bkcols)
      
      # A temp file to save the output.
      # This file will be removed later by renderImage
      #outfile <- tempfile(fileext = '.png')
      
      #png(outfile, width = 800, height = 600)
      
      
      img.dat <-
        x2$data_list$clean_img.dat %>% subsetPixels(tf_list)
      #previous line will revert to original image after saving to x2$overview_peaks_sel
      #img.dat<-x2$mytable_selected %>% subsetPixels(x2$tf_list)
      
      x <-
        coord(img.dat)[, 1]
      y <- coord(img.dat)[, 2]
      runs = Cardinal::run(img.dat)
      req(runs)
      a <- data.frame(x, y, runs)
      cols <- x2$bkcols[tf_list]
      
      cat(names(table(cols)))
      cols[cols == 0] <- 100 #for h/dbscan?
      
      
      
      
      
      #if(input$umap_cols=="Reduced2" | input$umap_cols=="Reduced" ){
      p1 <- ggplot2::ggplot(a, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_tile(fill = cols) +
        ggplot2::facet_wrap(. ~ runs) +
        ggplot2::theme_minimal() +
        ggplot2::scale_y_continuous(trans = "reverse")
      # } else { #for h/dbscan
      #   if(length(table(cols))/length(pals::polychrome()) >1) {
      #     ff<-floor(length(table(cols))/length(pals::polychrome()))+1
      #     col_vector=rep(as.vector(pals::polychrome()), ff)
      #   }else {
      #     col_vector=as.vector(pals::polychrome())
      #   }
      #
      #   #h/dbscan plot
      #   p1<-ggplot2::ggplot(a, ggplot2::aes(x=x, y=y))+
      #     ggplot2::geom_tile(ggplot2::aes( fill=as.factor(cols)))+
      #     ggplot2::facet_wrap(.~runs)+
      #     ggplot2::theme_minimal()+
      #     ggplot2::scale_y_continuous(trans = "reverse")+
      #     ggplot2::scale_fill_manual(values=col_vector) +
      #     ggplot2::guides(colour=FALSE, fill=FALSE)
      # }
      #
      
      print(p1)
      
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
    
    #lets have the user choose which values to keep from clust, then regenerate the dataset based on those values.
    
    output$colormaps <- renderImage({
      if (input$seg_choice == "anat_seg") {
        #only output one colormap for anatomical segmentatin
        
        outfile <- tempfile(fileext = '.png')
        
        png(outfile, width = 800, height = 300)
        
        x2 <- preproc_values()[["x2"]]
        req(x2$mytable_selected)
        req(input$seg_pdat_col)
        req(input$show_dat)
        
        
        #TODO - deal with segmentation from exisitng info, without running UMAP
        img.dat <- x2$mytable_selected
        if (!is.null(x2$pdat_anat) &&
            identical(x2$umap_name, runNames(x2$mytable_selected))) {
          pData(img.dat) <- x2$pdat_anat
        } else if (!is.null(x2$pdat_anat) &&
                   dim(x2$pdat_anat)[1] == ncol(img.dat)) {
          pData(img.dat) <- x2$pdat_anat
        }
        dat <- as.data.frame(pData(img.dat))
        
        cols <-
          dat[, input$seg_pdat_col][dat[, input$seg_pdat_col] %in% input$show_dat]
        a <- as.data.frame(table(cols))
        
        if (input$seg_pdat_col != "Rcol_reduced") {
          a$cols2 <-
            as.numeric(as.factor(as.data.frame(table(cols))$cols)) + 1
        } else {
          a$cols2 <- a$cols
        }
        
        
        
        req(a$cols2)
        
        p1 <-
          ggplot2::ggplot(a, ggplot2::aes(x = cols2, y = log(Freq))) +
          ggplot2::geom_col(fill = a$cols2) +
          ggplot2::geom_text(
            ggplot2::aes(label = cols),
            position = ggplot2::position_dodge(width = 0.9),
            angle = 90,
            size = 10,
            hjust = 1
          ) +
          ggplot2::theme_minimal() +
          ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 1)) +
          ggplot2::ggtitle(input$show_dat)
        #browser() TODO figure out plotting errors during segmentation
        
        print(p1)
        
        
      } else {
        x2 <- preproc_values()[["x2"]]
        req(x2$data_list)
        # A temp file to save the output.
        # This file will be removed later by renderImage
        outfile <- tempfile(fileext = '.png')
        
        png(outfile, width = 800, height = 300)
        
        
        
        
        col_red <- x2$data_list$umap_separation$col_reduced
        col_dbscan <- x2$data_list$dbscan_umap_separation$cluster
        col_hdbscan <- x2$data_list$hdbscan_umap_separation$cluster
        
        #helper function to generate color map plots; a is the color vector
        plot_col_map <- function(a, titles = "plot") {
          #require(ggplot2)
          
          b <- data.frame(table(a))
          colnames(b) <- c("color",  "Freq")
          #if(input$umap_cols=="Reduced2" | input$umap_cols=="Reduced"){
          
          p1 <-
            ggplot2::ggplot(b, ggplot2::aes(x = color, y = log(Freq))) +
            ggplot2::geom_col(fill = b$color) +
            ggplot2::geom_text(
              ggplot2::aes(label = color),
              position = ggplot2::position_dodge(width = 0.9),
              angle = 90,
              size = 10,
              hjust = 1
            ) +
            ggplot2::theme_minimal() +
            ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2)) +
            ggplot2::ggtitle(titles)
          # } else {
          #
          #   browser
          #   if(length(table(a))/length(pals::polychrome()) >1) {
          #     ff<-floor(length(table(a))/length(pals::polychrome()))+1
          #     col_vector=rep(as.vector(pals::polychrome()), ff)
          #     print("Warning: more than 36 colors, color replication performed")
          #   }else {
          #     col_vector=as.vector(pals::polychrome())
          #   }
          #
          #   p1<-ggplot2::ggplot(b, ggplot2::aes(x=color,y=log(Freq)))+
          #     ggplot2::geom_col(ggplot2::aes(fill=color))+
          #     ggplot2::geom_text(ggplot2::aes(label=color), position=ggplot2::position_dodge(width=0.9), angle=90, size=10)+
          #     ggplot2::theme_minimal()+
          #     ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2))+
          #     ggplot2::ggtitle(titles)+
          #     ggplot2::scale_fill_manual(values=col_vector)+
          #     ggplot2::guides(fill="none")
          # }
          
          return(p1)
        }
        
        p1 <- plot_col_map(col_red, titles = "Reduced R colors")
        p2 <- plot_col_map(col_dbscan, titles = "DBSCAN colors")
        
        
        if (sum(is.na(x2$data_list$hdbscan_umap_separation$cluster)) <
            1) {
          p3 <- plot_col_map(col_hdbscan, titles = "HDBSCAN colors")
          print(gridExtra::grid.arrange(p1, p2, p3, nrow = 1))
        } else {
          print(gridExtra::grid.arrange(p1, p2, nrow = 1))
        }
      }
      
      
      
      dev.off()
      
      # Return a list containing the filename
      list(
        src = outfile,
        contentType = 'image/png',
        width = 800,
        height = 300,
        alt = "This is alternate text"
      )
    }, deleteFile = TRUE)
    
    #color choices for umap selection
    output$Color_choices <-
      renderUI({
        #https://gist.github.com/wch/4211337
        # If missing input, return to avoid error later in function
        #if(is.null(x2$bkcols))
        #   return()
        
        # Get the data set with the appropriate name
        x2 <- preproc_values()[["x2"]]
        dat <- (x2$bkcols)
        cols <- unique(dat)
        
        # Create the checkboxes and select them all by default
        checkboxGroupInput(ns("cols"),
                           "Choose colors",
                           choices  = cols,
                           selected = cols)
      })
    
    
    
    #adjust t/f list when colors selected
    observeEvent({
      input$show_dat
      input$cols
      1
    }, {
      #print(input$cols)
      x2 <- preproc_values()[["x2"]]
      req(x2$mytable_selected)
      if (input$seg_choice == "anat_seg") {
        
        
        if (!is.null(x2$pdat_anat)) {
          dat <- as.data.frame(x2$pdat_anat)
        } else {
          dat <- as.data.frame(pData(x2$mytable_selected))
        }
        
        
        # if(is.null(x2$pdat_anat) || !identical(x2$umap_name, runNames(x2$mytable_selected))) {
        #   dat<-as.data.frame(pData(x2$mytable_selected))
        # }else {
        #   dat<-as.data.frame(x2$pdat_anat)
        # }
        
        if (!is.null(input$show_dat)) {
          x2$tf_list <- dat[, input$seg_pdat_col] %in% input$show_dat
          x2$tf_list_anat <- x2$bkcols %in% input$cols
        }
        
      } else {
        x2 <- preproc_values()[["x2"]]
        req(x2$data_list)
        x2$tf_list <- x2$bkcols %in% input$cols
      }
      #print(x2$tf_list)
      #print(x2$bkcols)
      #print(paste(length(x2$tf_list), " is length of tf_list"))
      print(knitr::kable(table(x2$tf_list)))
    })
    
    #save processed data in a list when save button pressed and plot results
    
    observeEvent(input$store_proc, {
      x2 <- preproc_values()[["x2"]]
      
      req(x2$tf_list) #this will tell us  which pixels to use / store??
      print("")
      print("Storing processed data")
      showNotification("Storing processed Data")
      
      
      if (input$seg_choice == "anat_seg") {
        #find colnames only in processed data, and add those to the overview data
        
        if (!is.null(x2$tf_list_anat) 
             &&
             (sum(x2$tf_list_anat) != sum(x2$tf_list))) {
            
          print("Storing UMAP results as processed data")
          
          
          x2$tf_list <- x2$tf_list_anat
          aa <- x2$tf_list_anat
          x2$tf_list_anat <- NULL
          
          req(x2$data_list)
          #x2<-preproc_values()[["x2"]]
          
          if (is.null(x2$overview_peaks_sel$Rcol_reduced)) {
            x2$overview_peaks_sel$Rcol_reduced <- "NA"
          }
          if (is.null(x2$overview_peaks_sel$col_dbscan)) {
            x2$overview_peaks_sel$col_dbscan <- "NA"
          }
          
          if (is.null(x2$overview_peaks_sel$col_hdbscan)) {
            x2$overview_peaks_sel$col_hdbscan <- "NA"
          }
          
          tmp.img <- x2$mytable_selected[x2$tf_list]
          x2$pdat_anat <- x2$pdata_anat[x2$tf_list,]
          
          tmp.img$Rcol_reduced <-
            x2$data_list$umap_separation$col_reduced[x2$tf_list]
          tmp.img$col_dbscan <-
            x2$data_list$dbscan_umap_separation$cluster[x2$tf_list]
          tmp.img$col_hdbscan <-
            x2$data_list$hdbscan_umap_separation$cluster[x2$tf_list]
          
          mytable_selected = tmp.img
          x2$tf_list <- rep(TRUE, ncol(mytable_selected))
        } else {
          mytable_selected <- x2$mytable_selected
          if (dim(pData(mytable_selected))[1] != dim(x2$pdat_anat)[1]) {
            print(paste0(
              "processed data size, ",
              dim(x2$pdat_anat)[1],
              " does not match image size, ",
              dim(pData(mytable_selected))[1]
            ))
            print("choose image matching processed data size before storing")
            print(paste0("possible name: ", x2$umap_name))
            return()
          } else {
            pData(mytable_selected) <- x2$pdat_anat
          }
        }
        
        if (!setequal(colnames(pData(x2$overview_peaks_sel)), colnames(pData(mytable_selected)))) {
          cols <-
            union(colnames(pData(x2$overview_peaks_sel)), colnames(pData(mytable_selected)))
          cols_overview <-
            cols[!cols %in% (colnames(pData(x2$overview_peaks_sel)))]
          cols_selected <-
            cols[!cols %in% (colnames(pData(mytable_selected)))]
          
          pData(x2$overview_peaks_sel)[, cols_overview] <- "NA"
          pData(mytable_selected)[, cols_selected] <- "NA"
        }
        
        
        
        tmp.img <- mytable_selected %>% subsetPixels(x2$tf_list)
        
      } else if (input$seg_choice == "fix_pix") {
        tmp.img <- x2$mytable_selected
        
        all_runs <-
          runNames(x2$overview_peaks_sel) #from the total set being analyzed
        runs <- runNames(tmp.img)
        idx <-
          which(all_runs %in% runs)  #are we assuming that the runs are stored in the same order?
        
        #put all runs into a separate list from the originally imported / peakpicked data
        for (i in 1:length(all_runs)) {
          x2$list_proc_img[[i]] <-
            x2$overview_peaks_sel %>% subsetPixels(Cardinal::run(x2$overview_peaks_sel) ==
                                                     all_runs[i])
          fData(x2$list_proc_img[[i]]) <-
            fData(x2$list_proc_img[[i]])[unique(colnames(fData(x2$list_proc_img[[i]])))] #clean up featuredata
        }
        
        #replace list elements with updated processed data
        for (i in 1:length(runs)) {
          x2$list_proc_img[[idx[i]]] <-
            tmp.img %>% subsetPixels(Cardinal::run(tmp.img) %in% runs[i])
          fData(x2$list_proc_img[[idx[i]]]) <-
            fData(x2$list_proc_img[[idx[i]]])[unique(colnames(fData(x2$list_proc_img[[idx[i]]])))] #clean up featuredata
        }
        
        
        
        
      } else{
        
        
        req(x2$data_list)
        #x2<-preproc_values()[["x2"]]
        
        
        if (is.null(x2$overview_peaks_sel$Rcol_reduced)) {
          x2$overview_peaks_sel$Rcol_reduced <- "NA"
        }
        if (is.null(x2$overview_peaks_sel$col_dbscan)) {
          x2$overview_peaks_sel$col_dbscan <- "NA"
        }
        
        if (is.null(x2$overview_peaks_sel$col_hdbscan)) {
          x2$overview_peaks_sel$col_hdbscan <- "NA"
        }
        
        #input$seg_field add this column to the rest of the data!
        
        #tmp.img<-x2$mytable_selected %>% subsetPixels(x2$tf_list)
        
        
        
        
        tmp.img <- tmp_umap_dat()[x2$tf_list]
        
        tmp.img$Rcol_reduced <-
          x2$data_list$umap_separation$col_reduced[x2$tf_list]
        tmp.img$col_dbscan <-
          x2$data_list$dbscan_umap_separation$cluster[x2$tf_list]
        tmp.img$col_hdbscan <-
          x2$data_list$hdbscan_umap_separation$cluster[x2$tf_list]
        
        #insert plot now?
        
      }
      
      all_runs <-
        runNames(x2$overview_peaks_sel) #from the total set being analyzed
      runs <- runNames(tmp.img)
      idx <-
        which(all_runs %in% runs)  #are we assuming that the runs are stored in the same order?
      
      # if x2$list_proc_img is NULL, then we need to create it
      if (is.null(x2$list_proc_img)) {
        x2$list_proc_img <- vector("list", length = length(all_runs))
               
        #put all runs into a separate list from the originally imported / peakpicked data
        for (i in 1:length(all_runs)) {
          x2$list_proc_img[[i]] <-
            x2$overview_peaks_sel %>% subsetPixels(Cardinal::run(x2$overview_peaks_sel) ==
                                                     all_runs[i])
          
            a<-(fData(x2$list_proc_img[[i]])[unique(colnames(fData(x2$list_proc_img[[i]])))]) #clean up featuredata
            
            #need to use massdatafram explicitly
            fData(x2$list_proc_img[[i]]) <- MassDataFrame(mz=a$mz, a %>% as.data.frame() %>% dplyr::select(-mz))
            
        }
      }
      
      #replace list elements with updated processed data ###TODO
      for (i in 1:length(runs)) {
        x2$list_proc_img[[idx[i]]] <-
          tmp.img %>% subsetPixels(Cardinal::run(tmp.img) %in% runs[i])
         a<- fData(x2$list_proc_img[[idx[i]]])[unique(colnames(fData(x2$list_proc_img[[idx[i]]])))] #clean up featuredata
         
         fData(x2$list_proc_img[[idx[i]]]) <- MassDataFrame(mz=a$mz, a %>% as.data.frame() %>% dplyr::select(-mz))
      }
      
      

      x2$data_list <- NULL

      x2$mytable_selected <-
        combine_card(x2$list_proc_img[input$mytable_rows_selected])
      #browser()
        #do.call(cbind, x2$list_proc_img[input$mytable_rows_selected]) #TODO check carefully!!
      # tmp_umap_dat<-reactive({return(tmp.img)})
      x2$tf_list<-rep(TRUE, ncol(combine_card(x2$list_proc_img[input$mytable_rows_selected])))
      
      
    })
    
    
    
    
    output$umap_choice <- renderUI({
      switch(
        input$seg_choice,
        "bk_seg" = list(
          numericInput(ns("min_dist"), "minimum distance for UMAP", 0.1),
          numericInput(ns("search_k"), "# of node to search", 50),
          
          numericInput(
            ns("nn"),
            "UMAP nearest neighbors",
            min = 2,
            max = 100,
            value = 5
          ),
          numericInput(ns("n_trees"), "Number of trees for UMAP", 50),
          numericInput(
            ns("set_op_mix_ratio"),
            "set_op_mix_ratio",
            min = 0,
            max = 1,
            value = 1
          ),
          textInput(
            ns("pca_umap"),
            "UMAP PCA components (NULL for no PCA)",
            value = NULL
          ),
          #checkboxInput(ns("use_python"), "Use Python?"),
          uiOutput(ns("use_python")),
          numericInput(ns("eps"), "DBSCAN eps", .15),
          numericInput(
            ns("minPts"),
            "H/DBSCAN min points / cluster",
            min = 2,
            value = 50
          ),
          numericInput(
            ns("cex"),
            "cex, tile plot size adjustment",
            min = 0.1,
            value = 2
          ),
          numericInput(
            ns("quant"),
            "R color reduction quantile (0-1)",
            min = 0,
            max = 0.99,
            value = .9
          ),
          checkboxInput(
            ns("fastmap"),
            "Use spatialFastmap? (experimental)",
            value =
              FALSE
          ),
          uiOutput(ns("fm_params")),
          actionButton(ns("umap1"), label = "Start UMAP"),
          actionButton(ns("action_dbscan"), label = "Re-run color clustering")
        ),
        "anat_seg" = list(
          numericInput(ns("min_dist"), "minimum distance for UMAP", 0.01),
          numericInput(
            ns("nn"),
            "UMAP nearest neighbors",
            min = 2,
            max = 100,
            value = 5
          ),
          numericInput(ns("n_trees"), "Number of trees for UMAP", 50),
          numericInput(
            ns("set_op_mix_ratio"),
            "set_op_mix_ratio",
            min = 0,
            max = 1,
            value = 1
          ),
          textInput(
            ns("pca_umap"),
            "UMAP PCA components (NULL for no PCA)",
            value = NULL
          ),
          numericInput(ns("eps"), "DBSCAN eps", .15),
          numericInput(
            ns("minPts"),
            "H/DBSCAN min points / cluster",
            min = 2,
            value = 50L
          ),
          #checkboxInput(ns("use_python"), "Use Python?"),
          uiOutput(ns("use_python")),
          numericInput(
            ns("cex"),
            "cex, tile plot size adjustment",
            min = 0.1,
            value = 2
          ),
          numericInput(
            ns("quant"),
            "R color reduction quantile (0-1)",
            min = 0,
            max = 0.99,
            value = .85
          ),
          checkboxInput(ns("fastmap"), "Use spatialFastmap?", value =
                          FALSE),
          uiOutput(ns("fm_params")),
          actionButton(ns("umap1"), label = "Start UMAP"),
          actionButton(ns("action_dbscan"), label = "Re-run color clustering"),
          uiOutput(ns('anat_pdata')) #pData selection for anatomical features
          
        )
      )
    })
    
    output$use_python <- renderUI({
      # if (input$use_python) {
      #   list(
      #     numericInput(ns("minSize"), "Python HDBSCAN min size", value = 15L),
      #     numericInput(
      #       ns("epsilon"),
      #       "Python HDBSCAN cluster_selection_epsilon",
      #       value = 0.2
      #     ),
      #     selectInput(
      #       ns("hdbscan_method"),
      #       "Method for HDBSCAN calculation",
      #       choices = list(
      #         "R DBSCAN package" = "R_hdbscan",
      #         "python HDBSCAN" = "py_hdbscan",
      #         None = "hdbscan_none"
      #       )
      #     )
      #
      #   )
      # } else {
      list(selectInput(
        ns("hdbscan_method"),
        "Method for HDBSCAN calculation",
        choices = list("R DBSCAN package" = "R_hdbscan", None =
                         "hdbscan_none")
      ))
      #}
    })
    
    output$anat_pdata <- renderUI({
      x2 <- preproc_values()[["x2"]]
      
      #extract pdata
      req(x2$mytable_selected)
      
      if (is.null(x2$pdat_anat)) {
        pdat <- as.data.frame(pData(x2$mytable_selected))
      } else {
        pdat <- as.data.frame(x2$pdat_anat)
      }
      
      v <-
        list(
          selectInput(
            ns("seg_pdat_col"),
            "Select pData field for visualization",
            choices = colnames(pdat),
            selected = "Rcol_reduced"
          ),
          uiOutput(ns("seg_pdat_show")),
          #interface to change names...
          radioButtons(
            ns("seg_rename"),
            "Choose pData field to annotate:",
            choices = c("Replace" = "replace",
                        "Add new field" = "add")
          ),
          uiOutput(ns("seg_field")),
          actionButton(ns("apply_annotation"), "Apply annotation")
          #selectInput("seg_pdat_col", "Select pData field to remove from visualization",
          #             choices = colnames(pdat)),
          #actionButton("remove_seg", "Remove from plot")
          
          
        )
      
      
      return(v)
      
    })
    
    output$seg_field <- renderUI ({
      x2 <- preproc_values()[["x2"]]
      
      
      cnames <-
        unique(c(colnames(as.data.frame(x2$pdat_anat)), colnames(as.data.frame(
          pData(x2$mytable_selected)
        ))))
      choices <-
        cnames[!cnames %in% c("Rcol_reduced",
                              "run",
                              "Plate",
                              "col_dbscan",
                              "col_hdbscan",
                              "x",
                              "y")]
      
      
      
      
      switch(
        req(input$seg_rename),
        "replace" = list(
          selectInput(
            ns("seg_field"),
            "Choose field to annotate",
            choices = choices,
            selected = input$seg_pdat_col
          ),
          uiOutput(ns("current_seg_choices")),
          p(HTML(paste("<b>OR</b>"))),
          textInput(ns("new_seg_value"), "Input new annotation value")
        ),
        "add" = list(
          textInput(ns("seg_field"), "Name of new field to populate"),
          textInput(ns("new_seg_value"), "New annotation value")
        )
      )
      
    })
    
    output$current_seg_choices <- renderUI ({
      x2 <- preproc_values()[["x2"]]
      x0 <- preproc_values()[["x0"]]
      
      
      #generate choices from all possible combinations of stored and existing information
      choices1 <- NULL
      choices2 <- NULL
      choices3 <- NULL
      if (!is.null(x2$pdat_anat) &&
          input$seg_field %in% colnames(x2$pdat_anat)) {
        choices1 = c("", unique(as.data.frame(x2$pdat_anat)[, (input$seg_field)]))
      }
      
      if (input$seg_field %in% colnames(as.data.frame(pData(x0$overview_peaks)))) {
        choices2 = unique(as.data.frame(pData(x0$overview_peaks))[, (input$seg_field)])
      }
      
      if (input$seg_field %in% colnames(as.data.frame(pData(x2$overview_peaks_sel)))) {
        choices3 = unique(as.data.frame(pData(x2$overview_peaks_sel))[, (input$seg_field)])
      }
      
      choices_combined <- unique(c(choices1, choices2, choices3))
      
      
      selectInput(
        ns("current_seg_value"),
        "Choose existing value",
        choices = c("", choices_combined),
        selected = ""
      )
      
      
    })
    
    output$seg_pdat_show <- renderUI ({
      x2 <- preproc_values()[["x2"]]
      req(x2$mytable_selected)
      
      
      #TODO - get values from existing x pdat_anat?
      
      if (!is.null(x2$pdat_anat)) {
        pdat <- as.data.frame(x2$pdat_anat)
      } else {
        pdat <- as.data.frame(pData(x2$mytable_selected))
      }
      
      # if(is.null(x2$pdat_anat) || !identical(runNames(x2$mytable_selected),(x2$umap_name))){
      #   pdat<-as.data.frame(pData(x2$mytable_selected))
      # } else {
      #   pdat<-as.data.frame(x2$pdat_anat)
      # }
      
      selectizeInput(
        ns("show_dat"),
        "Values to show",
        choices = unique(pdat[, input$seg_pdat_col]),
        selected = unique(pdat[, input$seg_pdat_col]),
        multiple = T
      )
    })
    
    output$fix_pix_opts <- renderUI({
      if (input$seg_choice == "fix_pix") {
        list(
          numericInput(ns("fix_pix_r"), "radius to search for neighbors", 2),
          numericInput(ns("fix_pix_n_thresh"), "min # of neighbors allowed", 3),
          checkboxInput(ns("fix_pix_t_f"), "Actually remove pixels?", value =
                          TRUE),
          actionButton(ns("fix_pix"), HTML("Remove isolated pixels"))
        )
      }
    })
    
    #Apply the results of annotations for specific tissue regions
    observeEvent(input$apply_annotation, {
      x2 <- preproc_values()[["x2"]]
      #req(input$seg_field)
      req(x2$mytable_selected)
      
      
      
      
      if (is.null(x2$pdat_anat) ||
          dim(x2$pdat_anat)[1] != dim(pData(x2$mytable_selected))[1]) {
        pdat <- (pData(x2$mytable_selected))
      } else {
        pdat <- x2$pdat_anat
      }
      
      if (input$seg_rename == "replace") {
        
        if (!input$seg_field %in% colnames(pdat)) {
          print(paste0(
            input$seg_field,
            " not found in existing pData, please check"
          ))
          return()
        }
        
        # if(input$seg_field=="Rcol_reduced"){ #check if is valid color to go in color vector??
        #   plotfuncions::
        # }
        
        #find which pixels are common between original dataset selection and umap selection
        #current dataset
        
        
        
        if (!is.null(x2$tf_list_umap) &&
            length(x2$tf_list_anat) != 0) {
          orig_idx <- which(x2$tf_list_umap)
          new_idx <- orig_idx[x2$tf_list_anat]
        } else {
          new_idx <- which(x2$tf_list)
        }
        
        
        #x2$tf_list_anat<-rep(TRUE, ncol(x2$mytable_selected))
        
        aa <- as.data.frame(pdat)[, input$seg_field]
        
        
        #a[1:20]<-"test1_doh"
        
        
        #test for replacement or addition of annotations based on which field is occupied
        
        if (is.null(input$current_seg_value)) {
          if (input$new_seg_value != "") {
            aa[new_idx] <- input$new_seg_value
          } else {
            print("Cannot determine field to apply annotation, please check")
            return(NULL)
          }
          
        } else {
          if (input$current_seg_value != "" && input$new_seg_value == "") {
            aa[new_idx] <- input$current_seg_value
            
          } else if ((input$current_seg_value == "" ||
                      is.null(input$current_seg_value)) &&
                     input$new_seg_value != "") {
            aa[new_idx] <- input$new_seg_value
            
          } else if (input$current_seg_value != "" &&
                     input$new_seg_value != "") {
            print(
              "Cannot have both current and new replacement values to segmentation, choose one!"
            )
            return()
          } else {
            print("cannot deterimine replacement value, check variables.")
            return()
          }
        }
        
        
        
        pdat[, input$seg_field] <- aa
        
        print("Storing annotation information")
        x2$pdat_anat <- pdat
        print(input$seg_field)
        print(table(as.data.frame(x2$pdat_anat)[, input$seg_field]))
        
        #pData(x2$mytable_selected)[,input$seg_field]<-pdat[,input$seg_field]
        
        
        
        
      } else if (input$seg_rename == "add") {
        if (input$seg_field %in% colnames(pdat)) {
          print(paste0(
            input$seg_field,
            " found in existing pData, please use new name"
          ))
          return()
        }
        
        pdat[x2$tf_list, input$seg_field] <- input$new_seg_value
        #pdat[!x2$tf_list,input$seg_field]<-"NA"
        
        print("Applying new annotation information")
        x2$pdat_anat <- pdat
        #pData(x2$mytable_selected)[,input$seg_field]<-pdat[,input$seg_field]
        
      }
      
    })
    
     observeEvent(input$save_imzml, {
      
      req(input$save_imzml)
       #create variables required for saving
       volumes <- c(wd = setup_values()[["wd"]], home = fs::path_home())
       shinyFiles::shinyFileSave(input, "save_imzml", roots = volumes, session = session)
       
       
       
      save_path <- shinyFiles::parseSavePath(volumes, input$save_imzml)
      if (nrow(save_path) == 0) return(NULL)
      
      
      filen <- as.character(save_path$datapath)
    #   
    # 
    # output$save_proc <- downloadHandler(
    #   filename = function() {
    #     if (input$seg_choice == "anat_seg") {
    #       fn = "/MSI-proc-anat_proc-"
    #     } else {
    #       fn = "/MSI-proc-bk_proc-"
    #     }
    #     
    #     paste0(getwd(), fn, Sys.Date(), ".rds")
    #   },
    #   content = function(file) {
        #pk_img <- x0$overview_peaks
        x2 <- preproc_values()[["x2"]]
        
        #pk_img <- Cardinal::combine(x2$list_proc_img)
        
        #issue with combining data, need to reorder with largest first
        pk_img_reorder <- x2$list_proc_img[order(sapply(x2$list_proc_img, ncol), decreasing = TRUE)]
        
        
        #use do.call to cbind the list
        pk_img<-combine_card(pk_img_reorder)
          #do.call(cbind, pk_img_reorder)
        
        
        #may need to restore the order based on the original data
        
        
        
        #cleanup feature data
        
          a<-fData(pk_img)[unique(colnames(fData(pk_img)))]
        
          fData(pk_img) <- MassDataFrame(mz=a$mz, a %>% as.data.frame() %>% dplyr::select(-mz))
        
          writeImzML(pk_img, filen)
        
        
      }
    )
    
    
  })
}
