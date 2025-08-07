### R/StatsPrepServer.R

StatsPrepServer <- function(id,  setup_values) {
  moduleServer(id, function(input, output, session) {
    # output$table <- DT::renderDataTable(mtcars)
    
    ns = session$ns #for dyanamic variable namespace
    
    output$variables <- renderPrint({
      setup_vals <- setup_values()
      paste("working directory=", setup_vals$wd)
    })
    
    #import parallel mode
    par_mode = reactive({
      setup_values()[["par_mode"]]
    })
    
    
    
    has.new.files <- function() {
      unique(list.files(setup_values()[["wd"]], recursive = T))
    }
    get.files <- function() {
      list.files(setup_values()[["wd"]], recursive = T)
    }
    
    # store as a reactive instead of output
    my_files <-
      reactivePoll(10, session, checkFunc = has.new.files, valueFunc = get.files)
    
    # any time the reactive changes, update the selectInput
    observeEvent(my_files(),
                 ignoreInit = T,
                 ignoreNULL = T,
                 {
                   #print(grep("rds|RData", my_files(), ignore.case =T, value=T  ))
                   updateSelectInput(
                     session,
                     ns('peakPickfile'),
                     choices = grep(
                       "rds|RData|imzML$",
                       my_files(),
                       ignore.case = T,
                       value = T
                     )
                   )
                 })
    
    
    
    ####STATISTICS####
    
    output$ui <- renderUI({
      #req(input$phen_cols_stats)
      switch(
        input$stats_test,
        "ssctest" = list(
          uiOutput(ns('test_membership')),
          uiOutput(ns('ssc_params')),
          uiOutput(ns("phen_interaction_stats")),
          uiOutput(ns('ssc_factors')),
          actionButton(ns("run_test"), label = "Run test")
        ),
        "meanstest" = list(
          uiOutput(ns('test_membership')),
          uiOutput(ns("phen_interaction_stats")),
          actionButton(ns("run_test"), label = "Run test")
        ),
        "spatialDGMM" = list(
          uiOutput(ns('test_membership')),
          uiOutput(ns('ssc_params')),
          uiOutput(ns("phen_interaction_stats")),
          checkboxInput(ns("var_filt"), "Filter by variance?", value = TRUE),
          numericInput(ns("var_thresh"), "Quantile filter for variance-based threshold", 0.8),
          actionButton(ns("run_test"), label = "Run test"),
          textInput(ns("plot_prefix"), label="Prefix for plotting / results output", value="output"),
          actionButton(ns("write_plots"), label = "save significant plots directory 'plots'")
        )
      )
    })
    
    
    #Setup and values for the statistics tab
    x5 <- reactiveValues(
      data_file = NULL,
      data_file_selected = NULL,
      mytable_stats_plate_selected = NULL,
      phen_options = NULL,
      phen_cols_stats = NULL,
      #from input$phen_col_stats, the test variable
      test_result = NULL,
      stats_results = NULL,
      plot_list = NULL,
      tf_list = NULL,
      #True/False list for selecting df to visualize
      test_result_feature_test = NULL,
      source_file = NULL,
      model_restore_path = NULL,
      dat_long_tech_avg = NULL,
      #stored initial data to create plots for visualization
      size_ok = NULL,
      groupsx = NULL,
      #full set of grouping variables for plot / output / grouping
      group_var = NULL,
      # names of groupign variables
      stats_table_filtered = NULL, #filtered stats table final form
      mytable_stats_plate_rows_selected = NULL #selected rows from stats table
    )
    
    output$final_data <- renderUI({
      selectInput(
        ns("stats_input_file"),
        "Imageset for analysis",
        grep(
          "imzML$|rds|RData",
          my_files(),
          ignore.case = T,
          value = T
        )
      )
    })
    
    observeEvent(input$action_read_stats, {
      req(input$stats_input_file)
      
      
      
      
      
      if(length(grep("rds", input$stats_input_file)>0)) {
        
        x5$data_file <- readRDS(input$stats_input_file)
      } else {
        x5$data_file <- readMSIData(input$stats_input_file)
        #x5$data_file <- readImzML(input$stats_input_file)
      }
      print(x5$data_file)
      
      if (!is.na(input$debug_n)) {
        message("selecting peaks for debugging")
        
        if (dim(x5$data_file)[1] < input$debug_n) {
          message("more debug peaks than features, keeping all peaks")
        } else {
          x5$data_file <- x5$data_file[1:input$debug_n, ]
          print(x5$data_file)
        }
      }
      
      
      x5$source_file <- input$stats_input_file
    })
    
    
    observeEvent(input$action_demo, {
      
      withProgress(message = "Loading rcc demo data, including peak picking. See CardinalWorkflows classification vignette for details.",
                   value = 0.5,
                   detail = "",
                   {
                     data(rcc, package = "CardinalWorkflows")
                     
                     rcc<-CardinalWorkflows::exampleMSIData("rcc")
                     
                     rcc_peaks <- rcc |>
                       normalize(method="tic") |>
                       peakProcess(SNR=3, filterFreq=FALSE,
                                   tolerance=0.5, units="mz")
                     
                     rcc_nobg <- subsetPixels(rcc_peaks, !is.na(diagnosis))
                     
                     overview_peaks <- rcc_nobg
                     
                     
                     x5$data_file <- overview_peaks
                     
                     
                   })
      
    })
    
    output$filename_header <- renderPrint({
      (paste0("working filename= ", req(x5$source_file)))
    })
    
    
    output$ssc_params <- renderUI({
      if (input$stats_test %in% c("ssctest")) {
        list(
          textInput(ns('sscr'), 'SSC r value', "1"),
          textInput(ns("sscs"), 'SSC s value ', "c(0, (4)^(1:2))"),
          textInput(ns("ssck"), 'SSC k value ', "3")
        )
      } else if (input$stats_test %in% c("spatialDGMM")) {
        list(textInput(ns('sscr'), 'DGMM r value', "1"),
             textInput(ns("sscs"), 'DGMM k value ', "5"))
      } else {
        return()
      }
    })
    
    output$phen_cols_stats <-
      renderUI({
        #https://gist.github.com/wch/4211337
        # If missing input, return to avoid error later in function
        if (is.null(x5$data_file))
          return()
        
        # Get the data set with the appropriate name
        
        dat <- (as.data.frame(pData(x5$data_file)))
        x5$phen_options <- colnames(dat)
        
        
        # Create the checkboxes and select them all by default
        selectInput(ns("phen_cols_stats"),
                    "Choose variable to test",
                    choices  = x5$phen_options)
      })
    
    output$grouping_variables <- renderUI({
      selectInput(
        ns("grouping_variables"),
        "Choose grouping variables",
        choices  = unique(c("none", x5$phen_options)),
        multiple = T
      )
      
    })
    
    output$grouping_variables_export <- renderUI({
      selectInput(
        ns("grouping_variables_export"),
        "Choose grouping variables",
        choices  = unique(c("none", x5$phen_options)),
        multiple = T
      )
      
    })
    
    
    output$test_membership <-
      renderUI({
        
        #https://gist.github.com/wch/4211337
        # If missing input, return to avoid error later in function
        if (is.null(x5$data_file))
          return()
        if (is.null(input$phen_cols_stats))
          return()
        
        #store test variable in x5 for other modules to use
        x5$phen_cols_stats <- input$phen_cols_stats
        
        
        
        # Get the members from variable being tested
        
        members <-
          unique(as.data.frame(pData(x5$data_file))[, input$phen_cols_stats])
        
        
        
        # Create the checkboxes and select them all by default
        selectizeInput(
          ns("test_membership"),
          "Choose members to include",
          choices  = members,
          selected = members,
          multiple=TRUE
        )
      })
    
    
    
    output$phen_interaction_stats <-
      renderUI({
        #https://gist.github.com/wch/4211337
        # If missing input, return to avoid error later in function
        if (is.null(x5$data_file))
          return()
        
        # Get the data set with the appropriate name
        # Create the checkboxes and select them all by default
        
        selectInput(
          ns("phen_interaction_stats"),
          "Choose sample definition variable(s)",
          choices  = c("none", 1, x5$phen_options),
          multiple = T
        )
      })
    
    output$output_factors <-
      renderUI({
        #https://gist.github.com/wch/4211337
        # If missing input, return to avoid error later in function
        if (is.null(x5$data_file))
          return()
        
        # Get the data set with the appropriate name
        # Create the checkboxes and select them all by default
        
        selectInput(
          ns("output_factors"),
          "Choose variables for modeling",
          choices  = c(x5$phen_options),
          multiple = TRUE,
          selected = input$phen_cols_stats
        )
      })
    
    output$ssc_factors <- renderUI({
      req(x5$data_file)
      req(input$phen_cols_stats)
      list(
        #checkbox for MIL
        checkboxInput(ns("ssc_MIL"), "Cross-validation with MIL", value = FALSE),
        numericInput(ns("ssc_MIL_folds"), "Number of folds/bags for MIL", 3),
        selectInput(
          ns("ssc_fold_vars"),
          "Choose variables for creating folds/bags",
          choices  =  x5$phen_options,
          multiple = TRUE,
          selected = ""
        )
      )
    })
    
    
    observeEvent(input$run_test, {
      message("Running stats test")
      gc()
      
      # Clear previous test results
      x5$test_result <- NULL
      x5$test_result_feature_test <- NULL 
      x5$stats_results <- NULL
      x5$data_file_selected <- NULL
      x5$dat_long_tech_avg <- NULL
      x5$stats_table_filtered <- NULL
      x5$size_ok <- NULL
      x5$groupsx <- NULL
      x5$plot_list <- NULL
      
      
      setCardinalBPPARAM(par_mode())
      setCardinalNChunks(setup_values()[["chunks"]])
      
      if (is.null(x5$data_file)) {
        print(
          "Select dataset to analyze and then press 'Read file for stats' button in the left panel"
        )
        return()
      }
      req(input$phen_cols_stats)
      
      withProgress(message = "Running stats test", value = 0.2, {
        print("running test")
        
        
        #create selected data file
        #For all cases!
        #create selected data files based on chosen user input.
        
        select_vec <-
          as.data.frame(pData(x5$data_file))[, input$phen_cols_stats] %in% input$test_membership
        
        x5$data_file_selected = x5$data_file[, select_vec]
        
        print("checking groups")
        
        group_var = input$phen_interaction_stats
        
        x5$group_var <- group_var
        
        if (sum(group_var %in% "1") == 1) {
          groupsx = 1
        } else if (sum(group_var %in% "none") > 0) {
          message("must have grouping for this test. Exiting.")
          return(NULL)
        } else{
          grouping = as.data.frame(pData(x5$data_file_selected))[, group_var]
          if (is.factor(grouping)) {
            groupsx = droplevels(as.data.frame(pData(x5$data_file_selected))[, group_var])
          } else {
            message(
              "grouping variable is not a factor-- do you need to create an interaction term?"
            )
            
            if (length(grouping) == 0) {
              message("must have value for grouping, or 'none' (for ANOVA)")
              showNotification(
                "must have value for grouping",
                type = "error",
                duration = 10
              )
              return(NULL)
            }
            message("continuing by creating a factor. Check the output carefully!")
            groupsx = droplevels(as.factor(interaction(grouping)))
          }
        }
        x5$groupsx <- groupsx
        print("groups done")
        
        if (input$stats_test == "meanstest") {
          
          
          select_vec <-
            as.data.frame(pData(x5$data_file))[, input$phen_cols_stats] %in% input$test_membership
          
          x5$data_file_selected = x5$data_file[, select_vec]
          
          
          
          message("performing means test")
          
          
          if (input$phen_cols_stats == "run") {
            pData(x5$data_file_selected)$run <- run(x5$data_file_selected)
          }
          
          
          #convert input$phen_cols_stats to factor in dataset
          #check if numbers
          tmp<-as.data.frame(pData(x5$data_file_selected)[,input$phen_cols_stats])
          if(all(grepl("^[0-9]+$", tmp[,1]))) {
            pData(x5$data_file_selected)[,input$phen_cols_stats]<-droplevels(factor(paste0("X.", tmp[,1])))
          } else {
            pData(x5$data_file_selected)[,input$phen_cols_stats]<-droplevels(factor(as.data.frame(pData(x5$data_file_selected))[,input$phen_cols_stats]))
          }
          
          #check how many groupsx are NA; if any, remove them
          na_vec<-is.na(x5$groupsx)
          
          message("removing pixels with NA sample values")
          
          
          dat<-x5$data_file_selected[,!na_vec]
          groupsx<-(as.character(x5$groupsx[!na_vec]))
          
          x5$data_file_selected<-dat
          x5$groupsx<-groupsx
          mt <- 
            try(meansTest(x5$data_file_selected[, ],
                          as.formula(paste0("~", input$phen_cols_stats)),
                          samples =
                            droplevels(as.factor(x5$groupsx))))
          
          if(class(mt) %in% "try-error") {
            print("meanstest failed. check data and data size")
            showNotification("meanstest failed. check data and data size", type="error")
            print(table(interaction(groupsx, as.data.frame(pData(x5$data_file_selected))[,input$phen_cols_stats])))
            message("means test failed")
            return()
          }
          
          
          if (class(mt@listData[[1]]$model) %in% "try-error") {
            print("cannot create meanstest summary, check grouping variable(s)")
            stop("means test failed")
          }
          
          incProgress(amount = 0.5, message = "extracting top features... can take a while for large datasets")
          print("extracting top features... can take a while for large datasets")
          
          
          stats_results <-
            as.data.frame(topFeatures(mt, n=length(mt), p.adjust="BH"))
          x5$stats_results <- stats_results
          
          x5$test_result <- mt
          
          x5$test_result_feature_test <- NULL
          bpstop(par_mode())
          
          print(
            "Means test complete, check Output table tab for table and check FDR cutoff if results not visible."
          )
          showNotification(
            "Means test complete, check Output table tab for table and check FDR cutoff if results not visible."
          )
          
          
          
        } else if (input$stats_test == "ssctest") {
          
          
          
          if (input$ssc_MIL) {
            if(identical(input$ssc_fold_vars, input$phen_cols_stats)) {
              showNotification("MIL fold variable cannot be the same as the test variable, exiting", type="error")
              message("MIL fold variable cannot be the same as the test variable, exiting")
              return(NULL)
            }
          }
          
          select_vec <-
            as.data.frame(pData(x5$data_file))[, input$phen_cols_stats] %in% input$test_membership
          
          x5$data_file_selected = x5$data_file[, select_vec]
          
          
          sscr = as.numeric(unlist(strsplit(input$sscr, split = ",")))
          sscs = as.numeric(eval(parse(text=input$sscs)))
          ssck = as.numeric(unlist(strsplit(input$ssck, split = ",")))
          
          if (sum(is.na(c(sscr, sscs, ssck))) > 0) {
            showNotification(paste0(
              "NA value detected in r, s or k. Values are:",
              paste(c(sscr, sscs, ssck), collapse = ",")
            ))
            return(NULL)
          }
          
          
          a <-
            droplevels(as.data.frame(pData(x5$data_file_selected)))
          dat <- x5$data_file_selected
          
          pData(dat) <- PositionDataFrame(coord(dat), run = a$run,  a[,!colnames(a) %in% c("x", "y", "run")])
          y = a[, input$phen_cols_stats]
          
          
          if (!input$ssc_MIL) {
            myfold <- groupsx
            
            res <-
              spatialShrunkenCentroids(dat,
                                       y,
                                       r = sscr,
                                       s = sscs,
                                       k = ssck)
            
            print((res))
            
            #create top features which will be a list
            tf_ssc <-
              topFeatures(res, n = dim(dat)[1])
            
            #create one dataframe from list. Use names of res as new column
            tf_ssc <-
              do.call(
                rbind,
                lapply(1:length(tf_ssc), function(i)
                  cbind(tf_ssc[[i]], model = names(res)[i])
                )
              ) 
            
            tf_ssc <- as.data.frame(tf_ssc[order(tf_ssc$statistic, decreasing=T), ])
            
            tf_ssc<- tf_ssc %>% dplyr::mutate(
              statistic = round(statistic, 2),
              centers = round(centers, 2),
              sd = round(sd, 2)
            )
            
          } else if (input$ssc_MIL) {
            
            
            if (length(input$ssc_fold_vars) == 0) {
              showNotification("No fold variables selected for MIL, exiting")
              return(NULL)
            }
            
            nfold=as.numeric(input$ssc_MIL_folds)
            if (nfold < 2) {
              showNotification("Number of folds must be greater than 1, exiting")
              return(NULL)
            }
            
            incProgress(amount=0.3, message="Initial SSC setup done, now performing cross-validation. Check console for details.")
            
            aa<-as.data.frame(pData(x5$data_file_selected))[, c(input$ssc_fold_vars, input$ssc_fold_vars, x5$group_var)]
            
            create_folds <- function(df2, group_col, run_col, f) {
              
              df2<-as.data.frame(df2)
              
              # Set the number of folds
              f=input$ssc_MIL_folds
              
              # Create a new column to store fold assignments
              df2$fold <- NA
              
              # Combine group_col and run_col, and remove duplicates
              all_group_cols <- unique(c(run_col))
              
              
              # Group data by 'Group.time.point' and 'run'
              grouped_data <- df2 %>% dplyr::group_by(!!!dplyr::syms(all_group_cols))
              
              # Assign folds using a round-robin assignment
              df2 <- grouped_data %>%
                dplyr::mutate(fold = (dplyr::cur_group_id() - 1) %% f + 1) %>%
                dplyr::ungroup()
              
              
              return(df2)
            }
            
            myfold <- try(create_folds(aa, input$phen_cols_stats, input$ssc_fold_vars, input$ssc_MIL_folds))
            
            if(class(myfold)[1]=="try-error") {
              showNotification("MIL fold creation failed, check variables and try again")
              return(NULL)
            }
            
            cv_ssc <-
              try(crossValidate(spatialShrunkenCentroids,
                                x = dat,
                                y = y,
                                folds = myfold$fold,
                                bags= groupsx,
                                r = sscr,
                                s = sscs
              ))
            
            
            
            if (class(cv_ssc) == "try-error") {
              showNotification(
                "SSC not able to complete-- try changing Grouping variable to create folds or fixing stray pixels in the Segmentation - UMAP tab, save new file, and try again",
                duration = 5
              )
              on.exit(bpstop(par_mode()), add = TRUE)
              return()
            }
            on.exit(bpstop(par_mode()), add = TRUE)
            
            a<-cv_ssc$average
            #choose best model
            max_macro<-max(a[,"MacroRecall"], na.rm=T)
            best_model<-which(a[,"MacroRecall"]==max_macro)
            if(length(best_model)==1) {
              print(a)
              mod_name<-rownames(a)[best_model]
              
              # Extract r and s values from the string
              params <- strsplit(mod_name, ",")[[1]]
              r_value <- as.numeric(gsub("r=", "", params[1]))
              s_value <- as.numeric(gsub("s=", "", params[2]))
              
              
              res <- spatialShrunkenCentroids(dat,
                                              y=y, r=r_value, s=s_value)
              
              
              
              topFeatures(res, n = 100)
              
            } else {
              showNotification("No best model found from CV check folds, exiting")
              print("interaction of groups and fold selection here")
              print(table(interaction(groupsx, dat$Group.time.point)))
              print(a)
              return(NULL)
            }
            
            tf_ssc <-
              topFeatures(res, n = dim(featureData(res))[1])
            
            tf_ssc<- as.data.frame(tf_ssc) %>% dplyr::mutate(
              statistic = round(statistic, 2),
              centers = round(centers, 2),
              sd = round(sd, 2), model = mod_name
            )
            
          }
          
          x5$stats_results <- as.data.frame(tf_ssc)
          
          if(class(res)=="SpatialShrunkenCentroids") {
            x5$test_result <- list(res)
            names(x5$test_result) <- mod_name
          }
          x5$test_result <- res
          bpstop(par_mode())
          print(
            "SSC test complete, check Output table tab for table and check FDR cutoff if results not visible."
          )
          
        } else if (input$stats_test == "spatialDGMM") {
          
          
          sscr = as.numeric(unlist(strsplit(input$sscr, split = ",")))
          sscs = as.numeric(unlist(strsplit(input$sscs, split = ",")))
          
          #add filter for variance
          if(input$var_filt) {
            var_red_peaks <- summarizeFeatures(x5$data_file_selected, stat=c(Variance="var"))
            
            var_red_peaks <- subsetFeatures(var_red_peaks, Variance >= quantile(Variance, input$var_thresh))
          } else {
            var_red_peaks <- x5$data_file_selected
          }
          
          
          print(var_red_peaks)
          
          
          dgmm <-
            try(spatialDGMM(
              var_red_peaks,
              r = sscr,
              k = sscs,
              groups = droplevels(groupsx),
              weights="gaussian"
            ))
          
          on.exit(bpstop(par_mode()), add = TRUE)
          
          
          if (class(dgmm) == "try-error") {
            print("check variables")
            showNotification(
              "DGMM not able to complete-- check variables, add or change variance filter and try again",
              duration = 5
            )
            
            print("DGMM failed")
            return()
            
          }
          
          x5$test_result <- dgmm
          
          print("spatialDGMM finished. Starting means test.")
          showNotification("spatialDGMM finished. Starting means test.")
          
          mtest <-
            meansTest(dgmm, as.formula(paste0("~ ", input$phen_cols_stats)))
          x5$test_result_feature_test <-
            mtest #in order to save later
          
          tf<-  as.data.frame(topFeatures(mtest, n = length(mtest), p.adjust = "BH"))
          
          #add ID from fData if it exists
          fdat<-as.data.frame(fData(var_red_peaks))
          if("ID" %in% colnames(fdat)) {
            tf <- tf %>%
              dplyr::mutate(ID = fdat$ID[match(tf$mz, fdat$mz)]) %>%  # Add the new column
              dplyr::relocate(ID, .after = mz)  # Relocate it after the "mz" column
          }
          
          x5$stats_results <- tf
          
          print((x5$test_result))
          bpstop(par_mode())
          print(
            "Spatial DGMM test complete, check Output table tab for table and check FDR cutoff if results not visible."
          )
        }
      })
    })
    
    observeEvent(input$data_export, {
      req(input$grouping_variables_export)
      
      showNotification("Starting data export.")
      
      numdat <- as.matrix(spectra(x5$data_file))
      rownames(numdat) <- mz(x5$data_file)
      pdat <- as.data.frame(pData(x5$data_file))
      
      dat <- cbind(pdat, t(numdat))
      
      dat_long_tech_avg <-
        (dat) %>%  dplyr::group_by_at(c(input$grouping_variables_export)) %>%
        dplyr::summarize(across(where(is.numeric), mean), .groups = 'keep')
      
      file_name=paste0("MSI.EAGLE_data_table_export_", Sys.Date(), ".tsv")
      
      print("Exporting datatable summarized by selected variables.")
      showNotification(paste("Exporting datatable in the working directory summarized by selected variables. Filename is: ", file_name, "\n"), duration = 10)
      
      write.table(
        dat_long_tech_avg,
        file = paste0("data_table_export_", Sys.Date(), ".txt"),
        sep = "\t",
        row.names = F
      )
      
    })
    
    
    observeEvent(input$mummichog, {
      req(x5$stats_results)
      
      
      if (input$stats_test %in% c("meanstest", "spatialDGMM")) {
        dat <- x5$stats_results %>%
          dplyr::mutate_at(dplyr::vars(pvalue, fdr), ~ (round(., 5))) %>% 
          dplyr::filter(fdr <= input$FDR_val)
      } else {
        print("not yet")
        showNotification("Not yet implemented for this test")
        return()
      }
      print("Writing Mummichog import file.")
      mchog_all <-
        cbind(
          mz = dat$mz,
          rt = dat[, "i"],
          p.value = dat[, "pvalue"],
          LR = dat[, "statistic"],
          fdr = dat[, "fdr"]
        )
      write.table(
        mchog_all,
        file = paste0(
          "mummichog_input_all_FDR_",
          input$FDR_val,
          "_",
          Sys.Date(),
          ".txt"
        ),
        sep = "\t",
        row.names = F
      )
      
    })
    
    observeEvent(input$metaboanalyst, {
      req(x5$stats_results)
      
      if (input$stats_test %in% c("meanstest", "spatialDGMM")) {
        dat <- x5$stats_results %>%
          dplyr::mutate_at(dplyr::vars(pvalue, fdr), ~ (round(., 5))) %>%
          dplyr::filter(fdr <= input$FDR_val)
      } else {
        print("not yet")
        showNotification("Not yet implemented for this test")
      }
      
      print("Writing Metaboanalyst input file.")
      hmdb_all_mc2 <-
        as.data.frame(cbind(
          m.z = dat$mz ,
          p.value = dat[, "pvalue"],
          t.score = dat[, "statistic"]
        ))
      hmdb_all_mc1 <- cbind(m.z = dat$mz[order(dat$pvalue)])
      write.table(
        hmdb_all_mc2,
        file = paste0(
          "hmdb_input_mc2_all_FDR_",
          input$FDR_val,
          "_",
          Sys.Date(),
          ".txt"
        ),
        sep = "\t",
        row.names = F,
        quote = F
      )
      write.table(
        hmdb_all_mc1,
        file = paste0("hmdb_input_mc1_", "_", Sys.Date(), ".txt"),
        sep = "\t",
        row.names = F,
        quote = F
      )
    })
    
    
    observeEvent(input$save_stats_models, {
      req(x5$test_result)
      req(x5$data_file_selected)
      
      input_image_dataset <- x5$data_file_selected
      parent_image_dataset <- x5$source_file
      
      if (input$stats_test == "spatialDGMM") {
        print("Saving model")
        dgmm_model <- x5$test_result
        dgmm_seg_test <- x5$test_result_feature_test
        
        dgmm_stats_result <- x5$stats_results
        dgmm_dat_long_tech_avg <- x5$dat_long_tech_avg
        
        save(
          input_image_dataset,
          dgmm_dat_long_tech_avg,
          dgmm_model,
          dgmm_seg_test,
          dgmm_stats_result,
          parent_image_dataset,
          file = paste0(
            "DGMM_",
            parent_image_dataset,
            "_",
            Sys.Date(),
            ".RData"
          )
        )
        print("DGMM model saved in working directory")
        
        
      } else if (input$stats_test == "meanstest") {
        
        print("Saving model")
        means_model <- x5$test_result
        means_test <- x5$test_result_feature_test
        means_test_stats_result <- x5$stats_results
        means_test_dat_long_tech_avg <- x5$dat_long_tech_avg
        phen_options <- x5$phen_options
        phen_cols_stats <- input$phen_cols_stats
        phen_interaction_stats <- input$phen_interaction_stats
        
        
        
        save(
          input_image_dataset,
          means_test_dat_long_tech_avg,
          means_model,
          means_test,
          parent_image_dataset,
          means_test_stats_result,
          phen_options,
          phen_cols_stats,
          phen_interaction_stats,
          file = paste0(
            "Means_",
            parent_image_dataset,
            "_",
            Sys.Date(),
            ".RData"
          )
        )
        
        print("meanstest model saved in working directory")
        
      } else {
        print("not yet")
      }
      
    })
    
    
    
    observeEvent(input$restore_stats_models, {
      
      tryPath <- tryCatch(
        file.choose()
        ,
        error = function(e) {
          e
        }
      )
      
      if (inherits(tryPath, "error")) {
        x5$model_restore_path <- NULL
      } else {
        x5$model_restore_path <- tryPath
      }
      
      print(paste("restoring", x5$model_restore_path))
      
      if (!is.null(x5$model_restore_path)) {
        load(x5$model_restore_path)
      } else {
        print("Path not valid, choose a valid model file")
        return()
      }
      
      if (grepl("means", x5$model_restore_path, ignore.case = T)) {
        selected_test <- "meanstest"
        x5$data_file_selected <- input_image_dataset
        x5$test_result <- means_model
        x5$test_result_feature_test <- means_test
        x5$stats_results <- means_test_stats_result
        x5$phen_options <- phen_options
        
        updateSelectInput(
          session = session,
          ns("phen_cols_stats"),
          "Choose variable to test",
          choices  = x5$phen_options,
          selected = phen_cols_stats
        )
        
        updateSelectInput(
          session = session,
          ns("phen_interaction_stats"),
          "Choose sample definition variable(s)",
          choices  = c("none", 1, x5$phen_options),
          selected = input$phen_cols_stats
        )
        
        message(paste("Test variable is:", phen_cols_stats))
        message(paste(
          "Grouping variable(s):",
          paste(phen_interaction_stats, collapse = " ")
        ))
        
        
      } else if (grepl("dgmm", x5$model_restore_path, ignore.case = T)) {
        selected_test <- "spatialDGMM"
        x5$data_file_selected <- input_image_dataset
        x5$test_result <- dgmm_model
        x5$test_result_feature_test <- dgmm_seg_test
        x5$stats_results <- dgmm_stats_result
        
      } else {
        print("not yet")
      }
      
      updateRadioButtons(
        session,
        "stats_test",
        "Type of test",
        choices = list(
          "Means test" = "meanstest",
          "Spatial shrunken centroids" =
            "ssctest",
          "Spatial DGMM" = "spatialDGMM"
        ),
        selected = selected_test
      )
      
      
    })
    
    observeEvent(input$write_plots, {
      req(x5$stats_results)
      req(x5$stats_table_filtered)
      
      wd <- getwd()
      dir.create(file.path(wd, "plots"), showWarnings = FALSE)
      
      if (input$stats_test %in% c("meanstest", "spatialDGMM", "ssctest")) {
        nplots<-dim(x5$stats_table_filtered)[1]
        
        
        withProgress(
          message = paste0( "Saving ", 
                            nplots,
                            " plots"),
          detail = paste0(
            "folder in ",
            wd,
            ". Can be slow. Only way to cancel is restarting!"
          ),
          {
            
            
            for(i in (1:dim(x5$stats_table_filtered)[1])){
              
              
              
              p1<- try(plot_stats_results(x5 = x5,
                                          stats_table_rows_selected = i,
                                          stats_test = input$stats_test,
                                          phen_cols_stats = input$phen_cols_stats,
                                          group_var = x5$group_var,
                                          plot_choice = input$plot_choice,
              ))
              
              
              if (class(p1)[1] == "try-error") {
                print("Error in plotting")
                showNotification("Error in plotting, check variables and try again")
                return()
              }
              
              
              incProgress(amount = 1 / dim(x5$stats_table_filtered)[1])
              pdf(
                file = paste(
                  "plots/",
                  input$plot_prefix,
                  "_",
                  input$plot_choice,
                  "_",
                  x5$stats_table_filtered$i[i],
                  "_",
                  round(x5$stats_table_filtered$mz[i], digits = 4),
                  ".pdf",
                  sep = ""
                ),
                width = 6,
                height = 4
              )
              
              print(p1)
              
              dev.off()
            }
            
          }
        )
      } else{
        print("not yet")
        showNotification("Plots not available for this statistical test")
      }
    })
    
    output$mytable_stats_plate = DT::renderDataTable({
      req(x5$data_file_selected)
      DT::datatable(
        cbind(run = runNames(x5$data_file_selected), "run"),
        selection = list(mode = "multiple", selected = c(1:length(
          runNames(x5$data_file_selected)
        ))),
        caption = "Choose runs to visualize"
      )
    })
    
    #extract selected datasets from table
    
    observe({
      if (is.null(x5$data_file_selected))
        return()
      
      req(input$mytable_stats_plate_rows_selected)
      ids <- input$mytable_stats_plate_rows_selected
      x5$mytable_stats_plate_selected <-
        x5$data_file_selected %>% subsetPixels(Cardinal::run(x5$data_file_selected) %in% runNames(x5$data_file_selected)[ids])
      x5$mytable_stats_plate_rows_selected <- ids
    })
    
    observe({
      req(x5$data_file_selected)
      req(x5$mytable_stats_plate_selected)
      
      x5$tf_list <-
        !is.na(prodlim::row.match(as.data.frame(
          coord(x5$mytable_stats_plate_selected)
        ), as.data.frame(
          coord(x5$mytable_stats_plate_selected)
        )))
    })
    
    observe({
      req(x5$mytable_stats_plate_selected)
      req(x5$tf_list)
      
      img.dat <-
        x5$mytable_stats_plate_selected %>% subsetPixels(x5$tf_list)
      
      plot_card_server("plot_card_stats",
                       overview_peaks_sel = img.dat,
                       spatialOnly = TRUE)
      
    })
    
    
    output$stats_table = DT::renderDataTable({
      req(x5$stats_results)
      dat <- x5$stats_results
      
      #apply FDR filter and round
      if (input$stats_test %in% c("meanstest", "spatialDGMM")) {
        dat <- x5$stats_results %>%
          dplyr::filter(fdr <= input$FDR_val) %>%
          dplyr::mutate(
            pvalue = signif(pvalue, digits = 3),
            fdr = signif(fdr,  digits = 3),
            statistic = signif(statistic,  digits = 4)
          )
        
      } else {
        dat <- dat
        
      }
      
      
      if (input$stats_test %in% c("meanstest", "spatialDGMM", "ssctest")) {
        
        labels = (unique(as.data.frame(pData(
          x5$data_file_selected
        ))[, input$phen_cols_stats]))
        
        if (dim(dat)[1] == 0) {
          showNotification(
            "No features selected, check FDR threshold (set to 1 for troubleshooting",
            duration = 10,
            type = "error"
          )
          return(NULL)
        }
        
        # compute FC
        if (length(labels) == 2) {
          mean_spectra <-
            as.data.frame(
              summarizeFeatures(
                x5$data_file_selected %>% subsetFeatures(mz %in% dat$mz),
                groups = droplevels(as.factor(as.data.frame(
                  pData(x5$data_file_selected), na.rm=TRUE
                )[, input$phen_cols_stats]))) %>% fData()
              
            )
          
          
          #reorder to match x5$stats_results
          new_mean <-summarizeFeatures(
            x5$data_file_selected %>% 
              subsetFeatures(mz %in% dat$mz)
          ) %>% 
            fData %>% 
            as.data.frame %>% 
            dplyr::select(mean, mz)
          
          mean_spectra$mean<-new_mean$mean
          
          #reorder to match x5$stats_results  
          mean_spectra <-
            mean_spectra[match(dat$mz, mean_spectra$mz), ]
          
          vars<- colnames(mean_spectra)[colnames(mean_spectra) %in% paste0(as.character(labels), ".mean")]
          
          if(length(vars)!=2){
            print("Two group labels not found, checking for numeric Phenotype test variables")
            
            vars<- colnames(mean_spectra)[colnames(mean_spectra) %in% paste0("X", as.character(labels), ".mean")]
            
            if(length(vars)!=2) {
              showNotification("Two group labels not found, check Phenotype test variables", type="error")
              return()
            }
          }
          
          log2FC = try(formatC(log2(mean_spectra[, vars[2]] / mean_spectra[, vars[1]]), format="fg", digits=3))
          
          if(class(log2FC)=="try-error") {
            print('log2FC not computed, check data columns not named "mz", "count", "freq", "ID"')
            return()
          }
          
          
          if (dim(dat)[1] == length(log2FC)) {
            lab <- paste0("log2FC(", labels[2], "/", labels[1], ")")
            dat <- cbind(dat, log2FC)
            names(dat)[names(dat) == "log2FC"] <- lab
            dat <- cbind(dat, round(mean_spectra[, vars], 4))
            dat <-
              cbind(dat, max_mean_group = colnames(mean_spectra[, vars])[max.col(mean_spectra[, vars])])
            dat$mz <- round(dat$mz, 4)
            
            
            #if ID exists in mean_spectra, add it to dat, and reorder to have i, mz and ID first
            if(!is.null(mean_spectra$ID)){
              #if ID doesn't exist in DAT already, add it
              if(!"ID" %in% colnames(dat)){
                dat <- cbind(ID=mean_spectra$ID, dat)
              }
            }
            
            #check for mean column in mean_spectr and append it to dat
            if(!"mean" %in% colnames(dat)){
              dat <- cbind(dat, mean=round(mean_spectra$mean, 4))
            }
            
            
            #check for ID column, and if it exists, move it to first column after mz
            if("ID" %in% colnames(dat)){
              dat <- dat %>% dplyr::select(i, mz, ID, everything())
            }
            
          } else {
            print("stats table and fold change results have different lengths, not adding FC")
          }
          
        } else {
          print("more than two labels in test data, FC not computed")
          
          #calculate group statistics
          groups_clean<-as.factor(as.data.frame(
            pData(
              x5$data_file_selected %>% subsetFeatures(mz %in% dat$mz)
            )
          )[, input$phen_cols_stats])
          
          groups_clean<-droplevels(groups_clean)
          
          mean_spectra <-
            as.data.frame(fData(
              summarizeFeatures(
                x5$data_file_selected %>% subsetFeatures(mz %in% dat$mz),
                groups = groups_clean
              )))
          
          pattern <- "^mz$|^mean$|^mean\\..*$|^ID\\..|^ID$"
          
          # Use grep to identify columns that match the pattern
          cols_to_exclude <- grep(pattern, names(mean_spectra))
          
          # Exclude those columns using negative indexing
          df_filtered <- mean_spectra[, -cols_to_exclude]
          
          temp_means <-
            cbind( 
              round(df_filtered[], 4),
              max_mean_group = colnames(df_filtered)[max.col(df_filtered)]
            )
          if(!is.null(mean_spectra$ID)){
            temp_means <- cbind(mz=round(mean_spectra$mz,4), ID=mean_spectra$ID, temp_means)
          } else {
            temp_means <- cbind(mz=round(mean_spectra$mz,4), temp_means)
          }
          
          
          dat$mz <- round(dat$mz, 4)
          
          dat <- merge(dat, temp_means, by=c("mz"))
          
          dat2<-dat
          #if ID.x and ID.y are the same, remove ID.y and rename ID.x to ID
          if("ID.x" %in% colnames(dat2) & "ID.y" %in% colnames(dat2)){
            dat2$ID <- ifelse(dat2$ID.x == dat2$ID.y, dat2$ID.x, NA)
            dat2 <- dat2[, !grepl("ID.y", colnames(dat2))]
            dat2 <- dat2[, !grepl("ID.x", colnames(dat2))]
            #relocate ID to first column after mz using dplyr
            dat2 <- dat2 %>% dplyr::select(i, mz, ID, everything())
            dat<-dat2
          }
          
          x5$stats_table_filtered <- dat
          
        }
        
        
      }
      
      #store filtered table for later use
      x5$stats_table_filtered <- dat
      
      DT::datatable(
        dat,
        selection = list(mode = "multiple", selected = "none"),
        caption = "Results",
        editable = FALSE,
        extensions = c("Buttons"),
        options = list(
          dom = 'Bfrtipl',
          buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
        )
      )
    }, server = FALSE)
    
    
    output$plot13 <- renderImage({
      req(x5$data_file_selected)
      req(x5$tf_list)
      req(input$stats_table_rows_selected)
      
      # A temp file to save the output.
      # This file will be removed later by renderImage
      outfile <- tempfile(fileext = '.png')
      
      png(outfile, width = 800, height = 600)
      
      
      plate1 <-
        x5$data_file_selected[, Cardinal::run(x5$data_file_selected) %in% runNames(x5$data_file_selected)[1]]
      
      
      
      plate_sel <-
        x5$mytable_stats_plate_selected[,x5$tf_list]
      
      nplots = length(runNames(plate_sel)) + 1
      
      p2 <-
        image(
          plate_sel,
          input$phen_cols_stats,
          key = T,
          col = ggsci::pal_npg()(10)
        )
      
      print(p2)
      
      
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
    
    #create user interface choices for plotting based on model type
    output$plot_choice_ui <- renderUI({
      req(x5$stats_results)
      
      if (input$stats_test %in% c("meanstest")) {
        plot_choices <- list(
          "Cardinal / matter" = "cardinal",
          "ggplot2 " = "ggplot",
          "MSI image" = "msi_image",
          "Test Phenotype" = "groupings"
        )
      } else if (input$stats_test %in% c("ssctest")) {
        plot_choices <- list(
          "Means Plot" = "means_plot",
          "SSC ion image" = "ion_image",
          "T-statistic" = "t_statistic",
          "Test Group" = "groupings",
          "MSI image" = "msi_image"
        )
      } else if (input$stats_test %in% c("spatialDGMM")) {
        plot_choices <- list(
          "DGMM means test results" = "dgmm_means_test",
          "DGMM segment ranks image" = "dgmm_ranks",
          "DGMM estimated segment parameters" = "dgmm_params",
          "MSI image" = "msi_image"
        )
      } else {
        plot_choices <- list()
      }
      
      selectInput(
        ns("plot_choice"),
        "Choose plot type:",
        choices = plot_choices
      )
    })
    
    
    
    
    #Plots selected statistical results
    output$plot11 <- renderImage({
      req(input$stats_table_rows_selected)
      
      
      # Setup to ensure proper cleanup
      on.exit({
        if (!is.null(par_mode())) {
          bpstop(par_mode())  # Ensure parallel workers are stopped
        }
      }, add = TRUE)
      
      # A temp file to save the output.
      # This file will be removed later by renderImage
      outfile <- tempfile(fileext = '.png')
      
      png(outfile, width = 800, height = 600)
      
      if(input$stats_test %in% c("meanstest", "spatialDGMM", "ssctest")) {
        #plot selected statistical results
        setCardinalBPPARAM(par_mode())
        setCardinalNChunks(setup_values()[["chunks"]])
        
        
        
        #plot selected statistical results
        p1 <- plot_stats_results(x5 = x5,
                                 stats_table_rows_selected = input$stats_table_rows_selected,
                                 stats_test = input$stats_test,
                                 phen_cols_stats = input$phen_cols_stats,
                                 group_var = x5$group_var,
                                 plot_choice = input$plot_choice,
                                 chunks = setup_values()[["chunks"]]
        )
        print(p1)
        
      } else{
        print("Plot not supported yet")
      }
      
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
    
    proc_values <- reactive({
      list(x5 = x5,
           par_mode = par_mode)
    })
    return(proc_values)
    
    
  })
}