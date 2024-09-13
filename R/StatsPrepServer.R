### R/StatsPrepServer.R
StatsPrepServer <- function(id,  setup_values) {
  moduleServer(id, function(input, output, session) {
    # output$table <- DT::renderDataTable(mtcars)
    
    ns = session$ns #for dyanamic variable namespace
    
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
          #uiOutput(ns("output_factors")),
          uiOutput(ns('ssc_factors')),
          #check if we need this??
          actionButton(ns("run_test"), label = "Run test")
          ),
        "meanstest" = list(
          uiOutput(ns('test_membership')),
          uiOutput(ns("phen_interaction_stats")),
          actionButton(ns("run_test"), label = "Run test"),
          uiOutput(ns("output_factors"))
          #uiOutput(ns('anova_factors')),
          #actionButton(ns("save_stats_models"), label =
          #               "Save means test model"),
          #actionButton(ns("restore_stats_models"), "Restore means test model"))
          ),
        "spatialDGMM" = list(
          uiOutput(ns('test_membership')),
          uiOutput(ns('ssc_params')),
          uiOutput(ns("phen_interaction_stats")),
          uiOutput(ns("output_factors")),
          checkboxInput(ns("var_filt"), "Filter by variance?", value = TRUE),
          numericInput(ns("var_thresh"), "Quantile filter for variance-based threshold", 0.8),
          #checkboxInput(ns("dgmm_means_test"), "Perform means test on DGMM results?", value = TRUE),
          actionButton(ns("run_test"), label = "Run test"),
          #textInput(ns("plot_prefix"), label="Prefix for plotting / results output", value="output"),
          actionButton(ns("write_plots"), label =
                         "save significant plots directory 'plots'"),
          actionButton(ns("save_stats_models"), "Save DGMM model"),
          actionButton(ns("restore_stats_models"), "Restore DGMM model")
          ),
        #https://www.datanovia.com/en/lessons/mixed-anova-in-r/ some ANOVA information
        #https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/
        #https://www.datanovia.com/en/lessons/anova-in-r/
        "anova" = list(
          radioButtons(
            ns("anova_type"),
            "Type of ANOVA",
            choices = list(
              "One way (Y ~ A)" = "anova_oneway",
              "Two way additive, Independent (Y ~ A + B)" =
                "anova_add",
              "Two way Interaction (Y ~ A * B)" =
                "anova_interaction",
              "Repeated measures, Between subjects interaction (Y ~ A*B +Error(C)" =
                "between_interacting" #,
              #"Repeated measures, independent (Y ~ A+B + Error (C)"="between_independent")
            )
          ),
          #uiOutput("phen_interaction_stats"),
          
          uiOutput(ns('grouping_variables')),
          uiOutput(ns("output_factors")),
          uiOutput(ns('anova_factors')),
          uiOutput(ns('test_membership')),
          actionButton(ns("run_test"), label = "Run test")
          #actionButton(ns("save_stats_models"), label =
           #              "Save ANOVA test model"),
          #actionButton(ns("restore_stats_models"), "Restore ANOVA test model")
        )
      )
      #
      #
      #
      # ,
     
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
      anova_vars = NULL,
      tf_list = NULL,
      #True/False list for selecting df to visualize
      test_result_feature_test = NULL,
      source_file = NULL,
      model_restore_path = NULL,
      dat_long_tech_avg = NULL,
      #stored initial data to create plots for visualization
      size_ok = NULL,
      #check membership size for ANOVA
      groupsx = NULL,
      #full set of grouping variables for plot / output / grouping
      group_var = NULL,
      # names of groupign variables
      stats_table_filtered = NULL #filtered stats table final form
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
      #if(is.null(input$stats_input_file))
      #  return()
      #print(input$stats_input_file)
      
      #renderPrint(cat(input$stats_input_file))
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
        
        if (input$stats_test == "anova" &&
            #!is.null(input$aov_vars2) &&
            input$anova_type == "anova_oneway") {
          #TODO .. fix this funkiness!!!
          
          members <-
            unique(as.data.frame(pData(x5$data_file))[, input$aov_vars1])
          
        } else if (input$stats_test == "anova" &&
                   !is.null(input$aov_vars2))  {
          #members<-paste(aov_members[,1], aov_members[,2], sep="____")
          members <-
            c(
              paste(input$aov_vars1, unique(as.data.frame(
                pData(x5$data_file)
              )[, input$aov_vars1]), sep = "____"),
              paste(input$aov_vars2, unique(as.data.frame(
                pData(x5$data_file)
              )[, input$aov_vars2]), sep = "____")
            )
          
          
        } else {
          members <-
            unique(as.data.frame(pData(x5$data_file))[, input$phen_cols_stats])
          
        }
        
        
        
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
    
    
    output$anova_factors <- renderUI({
      x5$anova_vars <- input$output_factors
      if (is.null(x5$anova_vars))
        return()
      
      v <- list()
      
      if (is.null(input$anova_type)) {
        anova_choices = c("Fixed", "Ind. Variable 1", "Ind. Variable 2")
        
        for (i in 1:length(anova_choices)) {
          v[[i]] <- #box(width = NULL, background = "blue",
            #    title = h6(anova_choices[i], style = "display:inline; font-weight:bold"),
            selectInput(ns(eval(paste0(
              "aov_vars", i
            ))),
            label = NULL,
            choices = c(x5$anova_vars))
          #)
          
        }
      } else if (input$anova_type == "anova_oneway") {
        anova_choices = c("Y Variable for one-way ANOVA")
        for (i in 1:length(anova_choices)) {
          v[[1]] = box(
            width = NULL,
            background = "blue",
            title = h6(anova_choices[i], style = "display:inline; font-weight:bold"),
            selectInput(
              ns(paste0("aov_vars", i)),
              label = NULL,
              choices = c(input$phen_cols_stats)
            )
          )
        }
        
        
        
        
        # for (i in 1:length(anova_choices)){
        #
        #   v[[i]] <- box(width = NULL, background = "blue",
        #                 title = h6(anova_choices[i], style = "display:inline; font-weight:bold"),
        #                 selectInput(ns(paste0("aov_vars",i)), label = NULL,choices = x5$anova_vars[!c(x5$anova_vars)%in%input$phen_cols_stats])
        
        
      } else if (input$anova_type == "anova_add") {
        anova_choices = c("Ind V1", "Ind V2")
        
        for (i in 1:length(anova_choices)) {
          v[[i]] <- box(
            width = NULL,
            background = "blue",
            title = h6(anova_choices[i], style = "display:inline; font-weight:bold"),
            selectInput(
              ns(paste0("aov_vars", i)),
              label = NULL,
              choices = c(x5$anova_vars)
            )
          )
        }
      } else if (input$anova_type == "anova_interaction") {
        anova_choices = c("Int V1", "Int V2")
        
        for (i in 1:length(anova_choices)) {
          v[[i]] <- box(
            width = NULL,
            background = "blue",
            title = h6(anova_choices[i], style = "display:inline; font-weight:bold"),
            selectInput(
              ns(paste0("aov_vars", i)),
              label = NULL,
              choices = c(x5$anova_vars)
            )
          )
        }
      } else if (input$anova_type == "between_interacting") {
        anova_choices = c("Variable 1", "Variable 2", "Subject")
        
        for (i in 1:length(anova_choices)) {
          v[[i]] <- box(
            width = NULL,
            background = "blue",
            title = h6(anova_choices[i], style = "display:inline; font-weight:bold"),
            selectInput(
              ns(paste0("aov_vars", i)),
              label = NULL,
              choices = c(x5$anova_vars)
            )
          )
        }
      } else if (input$anova_type == "between_independent") {
        anova_choices = c("Subject", "Ind. Variable 1", "Ind. Variable 2")
        
        for (i in 1:length(anova_choices)) {
          v[[i]] <- box(
            width = NULL,
            background = "blue",
            title = h6(anova_choices[i], style = "display:inline; font-weight:bold"),
            selectInput(
              ns(paste0("aov_vars", i)),
              label = NULL,
              choices = c(x5$anova_vars)
            )
          )
        }
      }
      return(v)
      
    })
    
    
    
    
    
    observeEvent(input$run_test, {
      message("Running stats test")
      gc()
      if (is.null(x5$data_file)) {
        print(
          "Select dataset to analyze and then press 'Read file for stats' button in the left panel"
        )
        return()
      }
      req(input$phen_cols_stats)
      
      # Setup to ensure proper cleanup
      # on.exit({
      #   if (!is.null(par_mode())) {
      #     bpstop(par_mode())  # Ensure parallel workers are stopped
      #   }
      # }, add = TRUE)
      
      withProgress(message = "Running stats test", value = 0.2, {
        print("running test")
        
        
        #create selected data file. The ANOVA variables necessiate a slightly more complex selection
        if (!is.null(input$anova_type) &&
            input$anova_type %in% c("anova_add",
                                    "anova_interaction",
                                    "between_interacting")) {
          if (!input$phen_cols_stats %in% c(input$aov_vars1, input$aov_vars2)) {
            print("test variable must be one of the ANOVA variables. Quitting.")
            return()
          }
          
          
          members <-
            data.frame(t = input$test_membership) %>% tidyr::separate(
              col = t,
              into = c("coln", "membs"),
              sep = "____"
            )
          
          #use do.call to extract?
          
          select_list <-
            lapply(1:dim(members)[1], function(i)
              ((
                as.data.frame(pData(x5$data_file))[, members[i, 1]] %in% members[i, 2]
              )))
          
          select_vec <-
            (Reduce("+", select_list)) == length(unique(members[, 1]))
          
          x5$data_file_selected = x5$data_file[, select_vec]
        } else {
          #For all other cases!
          #create selected data files based on chosen user input.
          
          select_vec <-
            as.data.frame(pData(x5$data_file))[, input$phen_cols_stats] %in% input$test_membership
          
          x5$data_file_selected = x5$data_file[, select_vec]
          #x5$data_file_selected = x5$data_file[1:100,as.data.frame(pData(x5$data_file))[,input$phen_cols_stats]%in%input$test_membership]
          
        }
        
        print("checking groups")
        
        
        
        #TODO check this..!
        if (input$stats_test == "anova") {
          group_var = input$grouping_variables
        } else {
          group_var = input$phen_interaction_stats
        }
        
        # if(!is.null(group_var) && group_var=='none'){
        #   group_var=input$phen_cols_stats
        # }
        
        x5$group_var <- group_var
        
        # if(sum(group_var%in%input$phen_cols_stats)>0) {
        #   print("Test variable in Grouping variable(s), exiting")
        #   return()
        # }else
        if (sum(group_var %in% "1") == 1) {
          groupsx = 1
        } else if (input$stats_test == "anova" &&
                   sum(group_var %in% "none") > 0) {
          groupsx = "none"
        } else if (input$stats_test != "anova" &&
                   sum(group_var %in% "none") > 0) {
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
                "must have value for grouping, or can have 'none' (for ANOVA)",
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
          pData(x5$data_file_selected)[,input$phen_cols_stats]<-droplevels(factor(as.data.frame(pData(x5$data_file_selected))[,input$phen_cols_stats]))
          
          mt <- 
            try(meansTest(x5$data_file_selected[, ],
                      as.formula(paste0("~", input$phen_cols_stats)),
                      samples =
                        droplevels(as.factor(groupsx))))
          #mt2<-meansTest(x5$data_file_selected, as.formula(paste("~", input$phen_cols_stats)), groups=droplevels(x5$data_file_selected$Plate.Group))
          
          on.exit(bpstop(par_mode()), add = TRUE)
          
          if(class(mt) %in% "try-error") {
            print("meanstest failed. check data and data size")
            showNotification("meanstest failed. check data and data size", type="error")
            print(table(interaction(groupsx, as.data.frame(pData(x5$data_file_selected))[,input$phen_cols_stats])))
            stop("means test failed")
          }
          
          
          if (class(mt@listData[[1]]$model) %in% "try-error") {
            print("cannot create meanstest summary, check grouping variable(s)")
            stop("means test failed")
          }
          
          # check summary results are not NaN
          
          
          ntests = max(round(length(mt@listData) / 100), 1)
          
          # if (sum(is.nan(summary(mt[1:ntests])$PValue)) == ntests) {
          #   message("first 1% of PValues are NaN, checking median PValue")
          #   
          #   if (is.nan(summary(mt[median(1:length(resultData(mt)))])$PValue)) {
          #     message("Median value Nan, check grouping variable!  ")
          #     return(NULL)
          #     
          #   } else {
          #     message("median value non NaN, continuing")
          #   }
          # }
          incProgress(amount = 0.5, message = "extracting top features... can take a while for large datasets")
          print("extracting top features... can take a while for large datasets")
          
          
          stats_results <-
            as.data.frame(topFeatures(mt, n=length(mt), p.adjust="BH"))
          x5$stats_results <- stats_results
          
          #x5$stats_results<-Cardinal::subset(x5$stats_results, AdjP<as.numeric(input$FDR_val))
          #check for nan issue
          
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
                                       #x5$data_file_selected,
                                       y,
                                       # droplevels(as.data.frame(pData(x5$data_file_selected))[, input$phen_cols_stats]),
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
            #res <- spatialShrunkenCentroids(x=x5$data_file_selected, groups=groupsx, r=sscr, s=sscs, k=ssck)
                    
          
          #myfold=droplevels(run(x5$data_file_selected))
          
          # tf_list<-te %in% c("LF.PARAFILM", "UF.PARAFILM")
          # dat_trim<-dat[, tf_list]
          # aa <- droplevels(as.data.frame(pData(dat_trim)))
          # pData(dat_trim)<-PositionDataFrame(coord(dat_trim), run = aa$run,  aa)
          # myfold=run(dat_trim)
          #
          
          #work on this....
          #check for MIL
          
          
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
              
              # df2=aa
              # group_col=input$phen_cols_stats
              # run_col=input$ssc_fold_vars
          
            
            
              
              df2<-as.data.frame(df2)
              
              # Assuming your data is in a data frame named 'df'
              
              # Set the number of folds
              f=input$ssc_MIL_folds
              
              # Create a new column to store fold assignments
              df2$fold <- NA
              
              # Combine group_col and run_col, and remove duplicates
              #all_group_cols <- unique(c(group_col, run_col, x5$group_var))
              all_group_cols <- unique(c(run_col))
  
              
              # Group data by 'Group.time.point' and 'run'
              grouped_data <- df2 %>% dplyr::group_by(!!!dplyr::syms(all_group_cols))
              
              # Assign folds in a round-robin fashion
              # for (i in 1:nrow(grouped_data)) {
              #   group <- grouped_data[i, ]
              #   fold_assignment <- (i - 1) %% f + 1  # Round-robin assignment
              #   df2$fold[df2[[group_col]] == group[[group_col]] & df2[[run_col]] == group[[run_col]]] <- fold_assignment
              # }
              # 
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
                                #.fun = "spatialShrunkenCentroids",
                                folds = myfold$fold,
                                bags= groupsx, #pData(dat)[, input$phen_interaction_stats],
                                #r=sscr, s=sscs, .fold=run(x5$data_file_selected))
                                r = sscr,
                                s = sscs
              )) #,
            #k = ssck,
            # .process = FALSE,
            # .processControl = list(
            #   SNR = 3,
            #   tolerance = 15,
            #   units = "ppm"
            # )
            # #BPPARAM=SerialParam()
            #))
            
            
            
            if (class(cv_ssc) == "try-error") {
              #browser()
              showNotification(
                "SSC not able to complete-- try changing Grouping variable to create folds or fixing stray pixels in the Segmentation - UMAP tab, save new file, and try again",
                duration = 5
              )
              on.exit(bpstop(par_mode()), add = TRUE)
              #stop("SSC failed")
              return()
            }
            on.exit(bpstop(par_mode()), add = TRUE)
            #print((res))
            #print((cv_ssc))
            #graphics.off()
            #print(plot(summary(cv_ssc), Accuracy ~ s, type='b')) #move to main area eventually
            
            a<-cv_ssc$average
            #choose best model
            max_macro<-max(a[,"MacroRecall"], na.rm=T)
            best_model<-which(a[,"MacroRecall"]==max_macro)
            if(length(best_model)==1) {
              #browser()
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
              #message(gridExtra::tableGrob(table(interaction(groupsx, dat$Group.time.point))))
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
          
          #test_dgmm<-try(spatialDGMM(x5$data_file_selected[1,], r=sscr, k=sscs, groups=droplevels(as.data.frame(pData(x5$data_file_selected))[,input$phen_interaction_stats])))
          
          
          #add filter for variance
          if(input$var_filt) {
            var_red_peaks <- summarizeFeatures(x5$data_file_selected, stat=c(Variance="var"))
            
            #plot(var_red_peaks, "Variance", xlab="m/z", ylab="Intensity")
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
              weights="gaussian" #add option for weights at somepoint?
            ))
          
          on.exit(bpstop(par_mode()), add = TRUE)
          
          
          if (class(dgmm) == "try-error") {
            #table(as.data.frame(na.omit(dat_long_tech_avg)) %>% subset(mz == mzs[1]) %>% dplyr::select(.data[[input$aov_vars1]]))
            #table(as.data.frame(na.omit(dat_long_tech_avg)) %>% subset(mz ==
            #                                                             mzs[1]) %>% dplyr::select(.data[[input$aov_vars2]]))
            
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
          
          
          #topFeatures(stest, p.adjust="fdr", AdjP < .1)
          
          tf<-  as.data.frame(topFeatures(mtest, n = length(mtest), p.adjust = "BH"))
          #browser()
          
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
        }  else if (input$stats_test == "MIL") {
          browser()
          sscr = as.numeric(unlist(strsplit(input$sscr, split = ",")))
          sscs = as.numeric(unlist(strsplit(input$sscs, split = ",")))
          
          #test_dgmm<-try(spatialDGMM(x5$data_file_selected[1,], r=sscr, k=sscs, groups=droplevels(as.data.frame(pData(x5$data_file_selected))[,input$phen_interaction_stats])))
          
          
          #add filter for variance
          if(input$var_filt) {
            var_red_peaks <- summarizeFeatures(x5$data_file_selected, stat=c(Variance="var"))
            
            #plot(var_red_peaks, "Variance", xlab="m/z", ylab="Intensity")
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
              weights="adaptive" #add option for weights at somepoint?
            ))
          
          on.exit(bpstop(par_mode()), add = TRUE)
          
          
          if (class(dgmm) == "try-error") {
            #table(as.data.frame(na.omit(dat_long_tech_avg)) %>% subset(mz == mzs[1]) %>% dplyr::select(.data[[input$aov_vars1]]))
            #table(as.data.frame(na.omit(dat_long_tech_avg)) %>% subset(mz ==
            #                                                             mzs[1]) %>% dplyr::select(.data[[input$aov_vars2]]))
            
            print("check variables")
            showNotification(
              "DGMM not able to complete-- check variables, add or change variance filter and try again",
              duration = 5
            )
            bpstop(par_mode())
            stop("DGMM failed")
            
          }
          
          x5$test_result <- dgmm
          
          print("spatialDGMM finished. Starting means test.")
          showNotification("spatialDGMM finished. Starting means test.")
          
          mtest <-
            meansTest(dgmm, as.formula(paste0("~ ", input$phen_cols_stats)))
          x5$test_result_feature_test <-
            mtest #in order to save later
          
          
          #topFeatures(stest, p.adjust="fdr", AdjP < .1)
          
          tf<-  as.data.frame(topFeatures(mtest, n = length(mtest), p.adjust = "BH"))
          #browser()
          
          #add ID from fData if it exists
          fdat<-as.data.frame(fData(var_red_peaks))
          if("ID" %in% colnames(fdat)) {
            tf <- tf %>%
              dplyr::mutate(ID = fdat$ID) %>%  # Add the new column
              dplyr::relocate(ID, .after = mz)  # Relocate it after the "mz" column
          }
          
          
          
          
          x5$stats_results <- tf
          
          print(summary(x5$test_result))
          bpstop(par_mode())
          print(
            "Spatial DGMM test complete, check Output table tab for table and check FDR cutoff if results not visible."
          )
        } else if (input$stats_test == "anova") {
          
          req(input$aov_vars1)
          
          
          numdat <- as.matrix(spectra(x5$data_file_selected))
          rownames(numdat) <- mz(x5$data_file_selected)
          pdat <- as.data.frame(pData(x5$data_file_selected))
          
          dat <- cbind(pdat, t(numdat))
          
          dat_long <-
            (
              tidyr::pivot_longer(
                dat,
                cols = (ncol(pdat) + 1):ncol(dat),
                names_to = "mz",
                values_to = "response"
              )
            )
          
          
          #dat_long_tech_avg<- dat_long %>% group_by(eval(input$phen_interaction_stats), mz) %>%
          if (sum(input$grouping_variables %in% "none") == 1) {
            print("'none' selected for grouping, so no groupings will be used. Check carefully")
            dat_long_tech_avg <-
              dat_long %>% dplyr::group_by_at(c(unique(c(
                input$output_factors
              )), "mz")) %>%
              dplyr::summarize(tech_avg = mean(response),
                               .groups = "keep")
          } else {
            dat_long_tech_avg <-
              dat_long %>% dplyr::group_by_at(c(unique(
                c(input$grouping_variables, input$output_factors)
              ), "mz")) %>%
              dplyr::summarize(tech_avg = mean(response),
                               .groups = "keep")
          }
          
          
          
          x5$dat_long_tech_avg <- dat_long_tech_avg #all mzs?
          mzs <- fData(x5$data_file_selected)
          feats <- features(x5$data_file_selected)
          
          if (is.null(dat_long_tech_avg[, input$phen_cols_stats])) {
            print("Test variable not found!")
            return()
          }
          
          
          
          
          #x5$test_result<-dat_long_tech_avg  #TODO check
          
          #calculate TIC /pixel; requires specifying color vector in output results
          #usually used for tissue proportions TODO: check
          #test for sample_id
          if (!is.null(dat_long$sample_id)) {
            nmz = length(unique(dat_long$mz))
            dat_long_col_proportion <- dat_long_tech_avg %>%
              dplyr::group_by(sample_id) %>% summarize(
                group = group,
                run = run,
                sample_id = sample_id,
                region_col = region_col,
                n_color_by_sample = n_color_by_sample,
                sum_pix_sample = sum(n_color_by_sample),
                .groups = "keep"
              ) %>%
              dplyr::mutate(col_proportion = n_color_by_sample / sum_pix_sample)
          } else {
            print("no 'sample_id' thus proportions not calculated")
          }
          
          
          #eventually use lapply or bplapply???
          #create vector of mz values
          
          
          if (input$anova_type == "anova_oneway") {
            #this order for input$aov_vars: anova_choices=c("Subject", "Within condition", "Within time")
            
            
            #test for more than one variable / group / mz
            x = dat_long_tech_avg$mz[1]
            a <-
              as.data.frame(na.omit(dat_long_tech_avg))[as.data.frame(na.omit(dat_long_tech_avg))$mz %in% x, ]
            #a<-as.data.frame(na.omit(dat_long_tech_avg)) %>% subset(mz==x)
            if (sum(table(a[eval(input$aov_vars1)])) == dim(unique(a[eval(input$aov_vars1)]))[1]) {
              print("only one value per group, exiting. Check groups")
              return()
            }
            
            
            #create models for each mz value
            fm = as.formula(paste0("tech_avg ~", input$aov_vars1))
            anova_list <-
              lapply(mz(mzs), function(x)
                as.data.frame(na.omit(dat_long_tech_avg))[dat_long_tech_avg$mz %in% x, ] %>% rstatix::anova_test(formula =
                                                                                                                   fm))
            
            names(anova_list) <- feats
            x5$test_result <- anova_list
            
            pvals <- unlist(lapply(anova_list, '[[', 5))
            adjustp <- p.adjust(pvals)
            
            stats_results <-
              cbind(mz = mz(mzs),
                    mzs,
                    feature = feats,
                    pvals,
                    adjustp)
            names(stats_results) <-
              c("mz", "feature", "pvals", "pvals.BH")
            
            x5$stats_results <-
              as.data.frame(cbind(as.data.frame(mzs), feature = feats, pvals, adjustp))
            
            
            
            
            
          } else if (input$anova_type == "anova_add") {
            if (input$aov_vars1 == input$aov_vars2) {
              print("Both variables cannot be the same for two-way ANOVA")
              return()
            }
            
            
            #plots modelled from here
            #https://csdaw.github.io/ggprism/articles/pvalues.html
            
            
            fm = as.formula(paste0(
              "tech_avg ~",
              input$aov_vars1,
              "+",
              input$aov_vars2
            ))
            
            anova_list <- try(bplapply(mz(mzs), function(x)
              as.data.frame(na.omit(dat_long_tech_avg))[dat_long_tech_avg$mz %in% x, ] %>% rstatix::anova_test(fm)))
            
            if (class(anova_list) == "try-error") {
              
              table(as.data.frame(na.omit(dat_long_tech_avg)) %>% subset(mz ==
                                                                           mzs[1]) %>% select(.data[[input$aov_vars1]]))
              table(as.data.frame(na.omit(dat_long_tech_avg)) %>% subset(mz ==
                                                                           mzs[1]) %>% select(.data[[input$aov_vars2]]))
              
              print("check group sizes")
              stop("ANOVA failed")
            }
            
            
            names(anova_list) <- feats
            x5$test_result <- anova_list
            
            pvals <-
              do.call(rbind.data.frame, (lapply(anova_list, '[[', 5)))
            colnames(pvals) <- paste0(anova_list[[1]]$Effect, ".p")
            adjustp <- apply(pvals, 2, p.adjust)
            colnames(adjustp) <-
              sapply(colnames(pvals), function(x)
                paste0(x, ".BH"))
            
            stats_results <-
              cbind(mz = mz(mzs),
                    feature = feats,
                    pvals,
                    adjustp)
            #names(stats_results)<-c("mz", "feature", "pvals", "pvals.BH")
            
            x5$stats_results <- as.data.frame(stats_results)
            
            
            
            
            
          } else if (input$anova_type == "anova_interaction") {
            if (input$aov_vars1 == input$aov_vars2) {
              print("Both variables cannot be the same for two-way ANOVA")
              return()
            }
            
            #formula for simple two way interaction ANOVA
            fm = as.formula(paste0(
              "tech_avg ~",
              input$aov_vars1,
              "*",
              input$aov_vars2
            ))
            
            #run a single test to check for issues
            
            
            test_aov <-
              try(as.data.frame(na.omit(dat_long_tech_avg))[dat_long_tech_avg$mz %in%
                                                              mz(mzs)[1], ] %>% rstatix::anova_test(fm))
            
            if (sum(class(test_aov) %in% "try-error") > 0) {
              
              table(as.data.frame(na.omit(dat_long_tech_avg)) %>% subset(mz ==
                                                                           mz(mzs)[1]) %>% select(.data[[input$aov_vars1]]))
              table(as.data.frame(na.omit(dat_long_tech_avg)) %>% subset(mz ==
                                                                           mz(mzs)[1]) %>% select(.data[[input$aov_vars2]]))
              
              print("check group sizes")
              stop("ANOVA failed")
            }
            
            nm = table(attributes(test_aov)$args$data[, input$phen_cols_stats])
            
            if (min(nm) < 3) {
              
              showModal(modalDialog(
                title = "Sample size warning!",
                HTML(paste(
                  names(nm), "n=  ", print(nm), "<p>"
                )),
                footer = tagList(
                  modalButton("Cancel"),
                  actionButton("ok", "Continue")
                  # actionButton("tryagain", "try_again")
                )
              ))
              
              observeEvent(input$ok, {
                x5$size_ok <- T
                removeModal()
              })
              
            } else {
              x5$size_ok <- T
            }
            
            
            observeEvent (x5$size_ok, {
              anova_list <-
                lapply(mz(mzs), function(x)
                  as.data.frame(na.omit(dat_long_tech_avg))[dat_long_tech_avg$mz %in% x, ] %>% rstatix::anova_test(fm))
              names(anova_list) <- feats[]
              x5$test_result <- anova_list
              
              pvals <-
                do.call(rbind.data.frame, (lapply(anova_list, '[[', 5)))
              colnames(pvals) <-
                paste0(anova_list[[1]]$Effect, ".p")
              adjustp <- apply(pvals, 2, p.adjust)
              colnames(adjustp) <-
                sapply(colnames(pvals), function(x)
                  paste0(x, ".BH"))
              
              stats_results <-
                cbind(mz = mz(mzs),
                      mzs,
                      feature = feats,
                      pvals,
                      adjustp)
              #names(stats_results)<-c("mz", "feature", "pvals", "pvals.BH")
              
              x5$stats_results <- as.data.frame(stats_results)
            })
            
            
            
            
          } else if (input$anova_type == "between_interacting") {
            #order of variables is: anova_choices=c( "Variable 1", "Variable 2", "Subject")
            
            
            if (max(table(
              c(
                input$aov_vars1,
                input$aov_vars2,
                input$aov_vars3
              )
            )) > 1) {
              print("Require unique subject and predictor variables")
              return()
            }
            
            #formula for simple two way interaction ANOVA
            fm = as.formula(
              paste0(
                "tech_avg ~",
                input$aov_vars1,
                "*",
                input$aov_vars2,
                "+ Error(",
                input$aov_vars3,
                ")"
              )
            )
            
            
            #run a single test to check for issues
            library(rstatix)
            
            
            
            
            two_way_rep <-
              function(dat,
                       mz_sel,
                       subject,
                       between,
                       within) {
                print(mz_sel)
                #   tryCatch( expr = {
                (na.omit(dat)) %>%
                  subset(mz == mz_sel) %>%
                  dplyr::group_by_(between, within, subject) %>%
                  dplyr::summarize(tech_avg = mean(tech_avg),
                                   .groups = "keep") %>%
                  dplyr::ungroup() %>%
                  rstatix::anova_test(
                    dv = tech_avg,
                    wid = subject,
                    between = between,
                    within = within
                  )
              } #,
            #   error=function(e){
            #     message("Encountered an error, check variables")
            #     print(e)
            #     return()
            #   }
            #   )
            #
            # }
            
            test_aov <- try(two_way_rep(
              dat = dat_long_tech_avg,
              mz_sel = mzs@mz[1],
              subject = input$aov_vars3,
              between = input$aov_vars1,
              within = input$aov_vars2
            ))
            
            if (class(test_aov) == "try-error") {
              
              table(as.data.frame(na.omit(dat_long_tech_avg)) %>% subset(mz ==
                                                                           mzs@mz[1]) %>% select(.data[[input$aov_vars1]]))
              table(as.data.frame(na.omit(dat_long_tech_avg)) %>% subset(mz ==
                                                                           mzs@mz[1]) %>% select(.data[[input$aov_vars2]]))
              table(as.data.frame(na.omit(dat_long_tech_avg)) %>% subset(mz ==
                                                                           mzs@mz[1]) %>% select(.data[[input$aov_vars3]]))
              print("check group sizes")
              stop("ANOVA failed")
            }
            
            
            
            
            nm = table(attributes(test_aov)$args$data[, input$phen_cols_stats])
            
            if (min(nm) < 3) {
              
              showModal(modalDialog(
                title = "Sample size warning!",
                HTML(paste(
                  names(nm), "n=  ", print(nm), "<p>"
                )),
                footer = tagList(
                  modalButton("Cancel"),
                  actionButton("ok", "Continue")
                  # actionButton("tryagain", "try_again")
                )
              ))
              
              observeEvent(input$ok, {
                x5$size_ok <- T
                removeModal()
              })
              
            } else {
              x5$size_ok <- T
            }
            
            
            observeEvent (x5$size_ok, {
              anova_list <-
                future_lapply(mzs[], function(x)
                  two_way_rep(
                    dat = dat_long_tech_avg,
                    mz = x,
                    subject = input$aov_vars3,
                    between = input$aov_vars1,
                    within = input$aov_vars2
                  ))
              names(anova_list) <- feats[]
              x5$test_result <- anova_list
              
              pvals <-
                do.call(rbind.data.frame, (lapply(anova_list, '[[', 5)))
              colnames(pvals) <-
                paste0(anova_list[[1]]$Effect, ".p")
              adjustp <- apply(pvals, 2, p.adjust)
              colnames(adjustp) <-
                sapply(colnames(pvals), function(x)
                  paste0(x, ".BH"))
              
              stats_results <-
                cbind(mz = mzs,
                      feature = feats,
                      pvals,
                      adjustp)
              #names(stats_results)<-c("mz", "feature", "pvals", "pvals.BH")
              
              x5$stats_results <- as.data.frame(stats_results)
            })
            
            
          }
          
          print(
            "ANOVA test complete, Output table tab for table and check FDR cutoff if results not visible."
          )
          
          
          
          
          #p_val(dat_long_tech_avg, mzs[2])
          
          #
          #
          #adjustp<-apply(pvals, 2, p.adjust)
          # colnames(adjustp)<-paste0(colnames(pvals),".BH")
          #
          # 
          # x5$stats_results<-as.data.frame(cbind(mz=mzs, feature=feats, pvals, adjustp))
          
        } #end of total anova testing
      })
    })
    
    observeEvent(input$data_export, {
      #req(input$output_factors)
      req(input$grouping_variables_export)
      
      
      
      numdat <- as.matrix(spectra(x5$data_file))
      rownames(numdat) <- mz(x5$data_file)
      pdat <- as.data.frame(pData(x5$data_file))
      
      dat <- cbind(pdat, t(numdat))
      
      dat_long_tech_avg <-
        (dat) %>%  dplyr::group_by_at(c(input$grouping_variables_export)) %>%
        dplyr::summarize(across(where(is.numeric), mean), .groups = 'keep')
      
      file_name=paste0("data_table_export_", Sys.Date(), ".txt")
      
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
      
      browser()
      if (input$stats_test %in% c("meanstest", "spatialDGMM")) {
        dat <- x5$stats_results %>%
          dplyr::mutate_at(dplyr::vars(PValue, AdjP), ~ (round(., 5))) %>% 
          dplyr::filter(AdjP <= input$FDR_val)
      } else if (input$stats_test %in% c("anova")) {
        dat <- x5$stats_results %>%
          print("not yet- anova")
        return()
      } else {
        print("not yet")
        return()
        #dat<-x5$stats_results
      }
      print("Writing Mummichog import file.")
      mchog_all <-
        cbind(
          mz = dat$mz,
          rt = dat[, "feature"],
          p.value = dat[, "PValue"],
          LR = dat[, "LR"],
          fdr = dat[, "AdjP"]
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
      
      browser()
      
      if (input$stats_test %in% c("meanstest", "spatialDGMM")) {
        dat <- x5$stats_results %>%
          dplyr::mutate_at(dplyr::vars(PValue, AdjP), ~ (round(., 5))) %>%
          dplyr::filter(AdjP <= input$FDR_val)
      } else {
        print("not yet")
        #dat<-x5$stats_results
      }
      
      print("Writing Metaboanalyst input file.")
      hmdb_all_mc2 <-
        as.data.frame(cbind(
          m.z = dat$mz ,
          p.value = dat[, "PValue"],
          t.score = dat[, "LR"]
        ))
      hmdb_all_mc1 <- cbind(m.z = dat$mz[order(dat$PValue)])
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
      
      # if(is.null(x5$dat_long_tech_avg)){
      #   print("need to generate meansplots first")
      # }
      
      
      
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
        
      } else if (input$stats_test == "anova") {
        print("Saving model")
        anova_list <- x5$test_result
        anova_test_stats_result <- x5$stats_results
        aov_vars = list(input$aov_vars1, input$aov_vars2, input$aov_vars3)
        aov_test_var <- input$phen_cols_stats
        aov_grouping <- input$grouping_variables
        aov_type <- input$anova_type
        aov_dat_long_tech_avg <- x5$dat_long_tech_avg
        
        
        
        save(
          input_image_dataset,
          anova_list,
          parent_image_dataset,
          anova_test_stats_result,
          aov_vars,
          aov_test_var,
          aov_grouping,
          aov_type,
          aov_dat_long_tech_avg,
          file = paste0(
            "ANOVA_",
            parent_image_dataset,
            "_",
            Sys.Date(),
            ".RData"
          )
        )
        
        print("Anova models saved in working directory")
        
        
        
      } else {
        print("not yet")
      }
      
    })
    
    
    
    observeEvent(input$restore_stats_models, {
      
      
      #req(x5$test_result)
      #req(x5$data_file_selected)
      #modeled from here: https://stackoverflow.com/questions/51191701/r-shiny-fileinput-large-files
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
      
      #restore anova models?? test for test type in filename etc
      
      if (grepl("means", x5$model_restore_path, ignore.case = T)) {
        selected_test <- "meanstest"
        x5$data_file_selected <- input_image_dataset
        #input$stats_input_file<-parent_image_dataset
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
        #input$stats_input_file<-parent_image_dataset
        x5$test_result <- dgmm_model
        x5$test_result_feature_test <- dgmm_seg_test
        x5$stats_results <- dgmm_stats_result
        
      } else if (grepl("ANOVA", x5$model_restore_path, ignore.case = T)) {
        selected_test <- "anova"
        x5$data_file_selected <- input_image_dataset
        #input$stats_input_file<-parent_image_dataset
        
        
        x5$test_result <- anova_list
        x5$stats_results <- anova_test_stats_result
        x5$dat_long_tech_avg <- aov_dat_long_tech_avg
        updateRadioButtons(inputId = ns("anova_type"),
                           selected = aov_type) #need full list?
        updateSelectInput(
          session = session,
          inputId = ns("grouping_variables"),
          "Choose grouping variables",
          choices  = unique(c("none", x5$phen_options)),
          selected = aov_grouping
        )
        updateSelectInput(
          session = session,
          inputId = ns("grouping_variables_export"),
          "Choose grouping variables",
          choices  = unique(c("none", x5$phen_options)),
          selected = aov_grouping
        )
        updateSelectInput(
          session,
          ns("phen_cols_stats"),
          "Choose variable to test",
          choices  = x5$phen_options,
          selected = aov_test_var
        )
        updateSelectInput(session, ns("aov_vars1"), selected = aov_vars[[1]])
        updateSelectInput(session, ns("aov_vars2"), selected = aov_vars[[2]])
        updateSelectInput(session, ns("aov_vars3"), selected = aov_vars[[3]])
        
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
          "Spatial DGMM" = "spatialDGMM",
          "ANOVA / t-test lm" = "anova"
        ),
        selected = selected_test
      )
      
      
    })
    
    observeEvent(input$write_plots, {
      req(x5$stats_results)
      
      
      
      wd <- getwd()
      dir.create(file.path(wd, "plots"), showWarnings = FALSE)
      #setwd(file.path(wd, "plots"))
      
      
      
      if (input$stats_test %in% c("meanstest", "spatialDGMM")) {
        withProgress(
          message = "Saving plots",
          detail = paste0(
            "folder in ",
            wd,
            ". Can be slow. Only way to cancel is restarting!"
          ),
          {
            req(x5$stats_results)
            
            stats_results_trimmed <-
              subset(x5$stats_results, AdjP <= input$FDR_val)
            
            
            for (m in 1:length(stats_results_trimmed$feature)) {
              incProgress(amount = 1 / length(stats_results_trimmed$feature))
              pdf(
                file = paste(
                  "plots/",
                  input$plot_prefix,
                  "_",
                  m,
                  "_",
                  round(stats_results_trimmed$mz[m], digits = 4),
                  ".pdf",
                  sep = ""
                ),
                # The directory you want to save the file in
                width = 4,
                # The width of the plot in inches
                height = 4
              ) # The height of the plot in inches
              
              
              #for means test
              if (input$stats_test == "meanstest") {
                p1 <-
                  image(
                    x5$test_result,
                    model = list(x5$stats_results$feature[m]),
                    normalize.image = 'linear',
                    contrast.enhance = "histogram",
                    key = F
                  )
              } else {
                #for DGMM
                p1 <-
                  image(
                    x5$test_result_feature_test,
                    model = list(x5$stats_results$feature[m]),
                    values = "mapping",
                    contrast.enhance = "histogram",
                    key = F
                  )
              }
              nplots <- length(p1$dpages) + 2
              
              
              print(p1, layout = n2mfrow(nplots))
              print(plot(
                x5$test_result,
                model = list(x5$stats_results$feature[m]),
                key = F
              ), layout = FALSE)
              if (input$stats_test == "spatialDGMM") {
                print(plot(
                  x5$test_result_feature_test,
                  model = list(x5$stats_results$feature[m])
                ),
                layout = FALSE)
              }
              title(paste(
                "mz= ",
                round(x5$stats_results$mz[m], 4),
                " FDR= ",
                round(x5$stats_results$AdjP[m], 2)
              ))
              
              
              
              
              #print(plot(x5$test_result, model=list(stats_results_trimmed$feature[m])))
              #title(paste("mz= ",round(x5$stats_results$mz[m],4)," FDR= ", round(x5$stats_results$AdjP[m],2)))
              # Step 3: Run dev.off() to create the file!
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
        # , size=ncol(x0$overview_peaks)), #FIX SIZE
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
    
    
    # #plot to make sure it looks okay
    # output$plot12 <- renderImage({
    #   req(x5$mytable_stats_plate_selected)
    #   req(x5$tf_list)
    #   # A temp file to save the output.
    #   # This file will be removed later by renderImage
    #   outfile <- tempfile(fileext = '.png')
    #
    #   png(outfile, width = 800, height = 600)
    #
    #
    #   img.dat<-x5$mytable_stats_plate_selected %>% subsetPixels(x5$tf_list)
    #
    #
    #
    #   #TODO - use an external module for visualization here
    #   # if (input$ion_viz=="viz_all") {
    #   #   
    #   #
    #   #  plusminus=input$plusminus_viz
    #     plusminus=1
    #     print(Cardinal::image(img.dat, plusminus=input$plusminus_viz, colorscale=pals::parula(255), contrast.enhance="histogram"))
    #
    #   #
    #   # } else if (input$ion_viz=="custom"){
    #   #   ion=input$mz_viz
    #   #   plusminus=input$plusminus_viz
    #   #   print(Cardinal::image(img.dat, mz=ion, plusminus=input$plusminus_viz, colorscale=pals::parula(255), contrast.enhance="histogram"))
    #   #
    #   # }
    #
    #   #print(Cardinal::image(mytable_selected(), mz=ion, plusminus=input$plusminus_viz))
    #
    #   dev.off()
    #
    #   # Return a list containing the filename
    #   list(src = outfile,
    #        contentType = 'image/png',
    #        width = 800,
    #        height = 600,
    #        alt = "This is alternate text")
    # }, deleteFile = TRUE)
    #
    
    
    output$stats_table = DT::renderDataTable({
      req(x5$stats_results)
      dat <- x5$stats_results
      

      # if (input$stats_test %in% c("meanstest", "spatialDGMM")) {
      #   dat$AdjP <- round(dat$AdjP, 4)
      #   dat$PValue <- round(dat$PValue, 4)
      # }
      #
      
      
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
        ))[, input$phen_cols_stats])) #may mess up consistent coloring with droplevels?
        
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
            
          #browser()
          #reorder to match x5$stats_results
          mean_spectra <-
            mean_spectra[match(dat$mz, mean_spectra$mz), ]
          
          vars=colnames(mean_spectra)[!colnames(mean_spectra) %in% c("mz", "count", "freq", "ID")]
          #remove vars that end in .1
          vars <- vars[!grepl("\\.1$", vars)]
          
          
          log2FC = try(formatC(log2(mean_spectra[, vars[2]] / mean_spectra[, vars[1]]), format="fg", digits=3))
          
          if(class(log2FC)=="try-error") {
            print('log2FC not computed, check data columns not named "mz", "count", "freq", "ID"')
            return()
          }
          
          
          if (dim(dat)[1] == length(log2FC)) {
            lab <- paste0("log2FC(", labels[2], "/", labels[1], ")")
            dat <- cbind(dat, log2FC)
            names(dat)[names(dat) == "log2FC"] <- lab
            dat <- cbind(dat, round(mean_spectra[, vars], 3))
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
          
          #find max from each row, only from columns not named mz or mean in temp_means
          

          #create temp dataset with only columns of interest
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
      
      
      
      
      
      #add fold change data if available for meanstest, ssctest, and spatialDGMM
      #with condition that only two groups are in the test group for now
      
      
      
      
      x5$stats_table_filtered <- dat
      
      DT::datatable(
        dat,
        selection = list(mode = "multiple", selected = "none"),
        caption = "Results",
        editable = FALSE,
        #server = FALSE,
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

      
      
      #print(p1, layout=n2mfrow(nplots))
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
               "ggplot2 " = "ggplot"
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
        #selected = "ggplot"
      )
    })
    
    
    
    
    #Plots selected statistical results
    output$plot11 <- renderImage({
      req(input$stats_table_rows_selected)
      if (is.null(input$stats_table_rows_selected))
        return()
      
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
      
      
      dat <- x5$stats_table_filtered
      
      #get row from orignial dataset corresponding to selected data
      m =  which(x5$stats_results$i %in% dat[input$stats_table_rows_selected, ]$i)
      if (input$stats_test == "ssctest") {
        m_model<-dat[input$stats_table_rows_selected, ]$model
        
        #find which rows of x5$stats_result match both ion and model for each value of input$stats_table_rows_selected
        idx<-NULL
        for(i in 1:length(m)){
          if(sum(x5$stats_results$i %in% dat[input$stats_table_rows_selected[i], ]$i & x5$stats_results$model %in% m_model[i])>0) {
            idx[i] <- which(x5$stats_results$i %in% dat[input$stats_table_rows_selected[i], ]$i & x5$stats_results$model %in% m_model[i])
          }
        }
        
        if(length(idx)>1){
         showNotification("Multiple ions selected, only first ion will be plotted", duration = 10)
         m = idx[1]
        }
      }
      
      if (input$stats_test %in% c("meanstest")) {
        req(x5$stats_results)
        
        mycols = ggsci::pal_npg()(as.data.frame(pData(x5$data_file_selected))[, input$phen_cols_stats] %>%
                                    factor() %>% droplevels() %>%
                                    levels() %>% length())
        
        
        #print base plot
        nplots=length(input$stats_table_rows_selected)
        mplot <- F
        
        
        #going to ignore mean plots for now, will come back if needed later
        # if (!is.null(input$output_factors)) {
        #   print("adding means plot across multiple conditions")
        #   mplot <- T
        #   nplots <- nplots + 1
        # }
        
       
        
        # p1 <-
        #   plot(
        #     x5$test_result,
        #     i = c(dat[input$stats_table_rows_selected, ]$i),
        #     col = mycols,
        #     las = 0,
        #     fill=T,
        #     panel.first = NULL,
        #     panel.last=NULL
        #   )
        # print(p1)

        # 
        
        #adjust a and b to work with multiple ions if m is longer than 1
        
        
        if(input$plot_choice == "ggplot") {
          #browser()
          #plots <- list()
          
        #extract all spectra and metadata for selected ions
          a<-lapply(m, function(x) {
            spectra(
            x5$data_file_selected[
              which(
                mz(x5$data_file_selected
                   ) %in% x5$stats_results$mz[x]
              ), 
              ]
            )
          })
          
          i_vals<-x5$stats_results$i[m]
          mz_vals<-x5$stats_results$mz[m]
          fdr_vals<-x5$stats_results$fdr[m]
          
          names(a) <- paste0("mz=", round(mz_vals, 4), " i=", i_vals, " FDR=", round(fdr_vals, 3))
          
          
          b<-pData(
            x5$data_file_selected
            )[, c(input$phen_cols_stats, x5$group_var)]
          
          dat_comb<-t(do.call(rbind,a))
          colnames(dat_comb)<-names(a)
          
          dat_comb<-cbind(dat_comb, b)
          
          summarized_df <- dat_comb %>% as.data.frame() %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(c(input$phen_cols_stats, x5$group_var)))) %>%  # Group by the categorical columns
            dplyr::summarize(dplyr::across(dplyr::everything(), mean, na.rm = TRUE))  # Summarize by calculating the mean for each ion
 
          long_format_df <- summarized_df %>%
            tidyr::pivot_longer(
              cols = -c(input$phen_cols_stats, x5$group_var),  # Exclude the grouping columns
              names_to = "ion",  # New column to hold the names of the numeric columns
              values_to = "value"  # New column to hold the values of the numeric columns
            ) %>% na.omit()
          
          n_samples <- long_format_df %>%
            dplyr::group_by(ion, !!rlang::sym(input$phen_cols_stats)) %>%
            dplyr::summarize(n = dplyr::n(), 
                             min_value = min(value, na.rm = TRUE),  # Calculate max value per ion and group
                             .groups = 'drop')
          
          
            p1<- ggplot2::ggplot(long_format_df, 
                          ggplot2::aes(x = !!ggplot2::sym(input$phen_cols_stats), 
                                       y = value,
                                       fill=!!ggplot2::sym(input$phen_cols_stats)))+
            
            ggplot2::geom_jitter(alpha = 0.5, width=0.2) +
            ggplot2::geom_boxplot(size = 1, alpha=0.8)+
            #ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray") + 
            ggplot2::theme_minimal() +
            ggplot2::labs(
              x = "",
              y = "Mean (normalized intensity, a.u.)",
              #title = title_t
            ) +
            ggprism::theme_prism()+
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))+
            #remove legend
            ggplot2::theme(legend.position = "none")+
            ggplot2::facet_wrap(~ion, scales = "free_y")+
            ggplot2::scale_fill_manual(values = mycols)+
            #add number of samples as annotation to plot
            ggplot2::geom_text(
              data = n_samples,
              ggplot2::aes(x = !!ggplot2::sym(input$phen_cols_stats),
                           y = 0.95*(min_value), 
                           label = paste0("n=", n)),
              vjust = -0.5
            )
          
          
          # 
          # 
          # plots<-lapply(m, function(x) {
          # #for (i in 1:length(m)) {
          #   
          #   #get spectra for currently selected table ions (m)
          #   a <-
          #     as.matrix(
          #       spectra(
          #         x5$data_file_selected[
          #           which(
          #             mz(x5$data_file_selected
          #                ) %in%x5$stats_results$mz[x]
          #           ), 
          #           ]
          #         )
          #       )
          #   
          #   #get metadata (pdata) for currently selected table ions (m)
          #   b <-
          #     as.data.frame(
          #       pData(
          #         x5$data_file_selected[
          #           which(
          #             mz(
          #               x5$data_file_selected
          #               ) %in% x5$stats_results$mz[x]
          #                   
          #             ), 
          #           ]
          #         )
          #       )
          #   
          #    
          #   
          #   dat2 <-  cbind(a=(a), b)
          #   
          #   #dat_long_tech_avg<- as.data.frame(na.omit(dat)) %>% dplyr::group_by_at(c(input$aov_vars1, input$aov_vars2, input$aov_vars3)) %>%
          #   dat_long_tech_avg <-
          #     dat2 %>% dplyr::group_by_at(c(input$phen_cols_stats, x5$group_var)) %>%
          #     dplyr::summarize(tech_avg = mean(a),
          #                      .groups = "keep")  %>% na.omit()
          #   
          #   #x5$dat_long_tech_avg<-dat_long_tech_avg
          #   
          #   fm <-
          #     as.formula(paste0("tech_avg~", input$phen_cols_stats))
          #   
          #   df_summary <- dat_long_tech_avg %>%
          #     dplyr::group_by(dplyr::across(dplyr::all_of(input$phen_cols_stats))) %>%
          #     dplyr::summarize(
          #       mean_tech_avg = mean(tech_avg),
          #       se_tech_avg = sd(tech_avg) / sqrt(dplyr::n())
          #     )
          #   
          #   # View the summary data
          #   df_summary
          #   
          #   title_t=paste(
          #     "mz= ",
          #     round(x5$stats_results$mz[x], 4),
          #     " FDR= ",
          #     formatC(x5$stats_results$fdr[x], format = "g", digits=2)
          #   )
          # 
          #   
          #   dat_long_tech_avg<-as.data.frame(dat_long_tech_avg)
          #   
          #   #check if factor
          #   dat_long_tech_avg[,input$phen_cols_stats] <- as.factor(dat_long_tech_avg[,input$phen_cols_stats])
          #   names(mycols) <- levels(as.data.frame(dat_long_tech_avg)[,input$phen_cols_stats])
          #   
          #   #calculate the number of samples in each group
          #   n_samples <- dat_long_tech_avg %>% 
          #     dplyr::group_by_at(c(input$phen_cols_stats)) %>% 
          #     dplyr::summarize(n=dplyr::n(), .groups = 'drop')
          #   
          #   
          # ggplot2::ggplot(dat_long_tech_avg, 
          #                               ggplot2::aes(x = !!ggplot2::sym(input$phen_cols_stats), 
          #                                            y = tech_avg,
          #                                            fill=!!ggplot2::sym(input$phen_cols_stats)))+
          #      
          #     ggplot2::geom_jitter(alpha = 0.5, width=0.2) +
          #     ggplot2::geom_boxplot(size = 1, alpha=0.8)+
          #     #ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray") + 
          #     ggplot2::theme_minimal() +
          #     ggplot2::labs(
          #       x = "",
          #       y = "Mean (normalized intensity, a.u.)",
          #       title = title_t
          #     ) +
          #     ggprism::theme_prism()+
          #     ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))+
          #     #remove legend
          #     ggplot2::theme(legend.position = "none") +
          #     ggplot2::scale_fill_manual(values = mycols)+
          #     #add number of samples as annotation to plot
          #     ggplot2::geom_text(
          #       data = n_samples,
          #       ggplot2::aes(x = !!ggplot2::sym(input$phen_cols_stats),
          #                    y = min(dat_long_tech_avg$tech_avg) - 0.05/mean(dat_long_tech_avg$tech_avg, na.rm=T), 
          #                    label = paste0("n=", n)),
          #       vjust = -0.5
          #     )
          # 
          #   
          #   
          # })
          # 
          # #browser()
          # # Number of plots
          # num <- length(plots)
          # 
          # # Calculate number of columns and rows dynamically
          # ncol_var <- ceiling(sqrt(num))
          # nrow_var <- ceiling(num / ncol_var)
          # 
          #print(gridExtra::grid.arrange(grobs = plots, ncol=ncol_var, nrow=nrow_var))
          print(p1)
            bpstop(par_mode())
        } else {
          #browser()
          m=input$stats_table_rows_selected
          #define mz and features of interest
          mz_vals=x5$stats_results$mz[m]
          i_vals<-x5$stats_results$i[m]
          fdr_vals<-x5$stats_results$fdr[m]
          names(i_vals)<-(paste0("mz= ",mz_vals, " FDR= ", round(fdr_vals,3)))
          
          
          p1 <- plot(x5$test_result, i = i_vals, col = mycols, las = 0, fill=T, free="y")
          print(p1)
          bpstop(par_mode())
          
        }
          # x=df_summary[[input$phen_cols_stats]]
          # y=df_summary$mean_tech_avg
          # ymin=df_summary$mean_tech_avg-df_summary$se_tech_avg
          # ymax=df_summary$mean_tech_avg+df_summary$se_tech_avg
          # 
          # 
          # 
          #  p2<-vizi( x=x, y = y,
          #   ymin = ymin, ymax=ymax) 
          #  
          # 
          #   p2<- add_mark(p2,
          #       "intervals",
          #    ) 
          #   
          # 
          # 
          # #
          # browser()
          # matter::as_facets(p1, p2, ncol = 1)
          # 
          # 
          
        #   print(
        #     gplots::plotmeans(
        #       fm,
        #       dat_long_tech_avg,
        #       mean.labels = T,
        #       digits = 1,
        #       las = 2
        #     ),
        #     layout = FALSE
        #   )
         # } else {
         #  print(p1)
         #}
      } else if (input$stats_test %in% c("spatialDGMM")) {
        req(x5$stats_results)
        
        
        mycols = ggsci::pal_npg()(as.data.frame(pData(x5$data_file_selected))[, input$phen_cols_stats] %>%
                                    factor() %>% droplevels() %>%
                                    levels() %>% length())
        
        #define mz and features of interest
        mz_vals=x5$stats_results$mz[m]
        i_vals<-x5$stats_results$i[m]
        fdr_vals<-x5$stats_results$fdr[m]
        names(i_vals)<-(paste0("mz= ",round(mz_vals, 4), " FDR= ", round(fdr_vals,3)))
        
        #if multiple ions are selected, create idx for the first one and show a message
        if(length(m)>1){
          showNotification("Multiple ions selected, only first ion will be plotted", duration = 10)
          idx = m[1]
        } else {
          idx = m
        }
        
        
        
        #dgmm plot
        p1<-plot(x5$test_result, i=i_vals, fill=T, free="xy")
        #means test plot
        p2<-plot(x5$test_result_feature_test, i=i_vals, col=mycols, fill=T, free="xy")
        p3<-image(x5$test_result, i=i_vals, smooth="bilateral", enhance="adaptive", scale=TRUE)
        p4 <-
          image(
            x5$data_file_selected,
            mz = (x5$stats_results$mz[idx]),
            
            enhance = "histogram",
            scale = T, free="xy"
            #col=mycols
          )
        
        
        #choose image to plot using switch from input$plot_choice
        plot_choice <- switch(input$plot_choice,
                              "dgmm_means_test" = p2,
                              "dgmm_ranks" = p3,
                              "dgmm_params" = p1,
                              "msi_image" = p4
        )
        
        
        
        
        print(plot_choice)
        
        #one day we could combine these...
        # matter::as_facets( p1, p3, free="xy"
        #                    )
        # 
        # 
        # 
        # if (mplot == T) {
        #   
        #   #require(gplots)
        #   #mzdat=subset(x5$test_result, mz==x5$stats_results$mz[m])
        #   a <-
        #     as.matrix(spectra(x5$data_file_selected[which(mz(x5$data_file_selected) %in%
        #                                                   x5$stats_results$mz[m]), ]))
        #   b <-
        #     pData(x5$data_file_selected[which(mz(x5$data_file_selected) %in% x5$stats_results$mz[m]), ])
        #   dat = cbind(run = Cardinal::run(x5$data_file_selected),
        #               b,
        #               value = a[1, ])
        #   
        #   dat_long_tech_avg <-
        #     as.data.frame(dat) %>% dplyr::group_by_at(c(
        #       "run",
        #       input$aov_vars1,
        #       input$aov_vars2,
        #       input$aov_vars3
        #     )) %>%
        #     dplyr::summarize(tech_avg = mean(value),
        #                      .groups = "keep")
        #   
        #   #x5$dat_long_tech_avg<-dat_long_tech_avg #single m/z?
        #   
        #   fm <-
        #     as.formula(
        #       paste0(
        #         "tech_avg~interaction(",
        #         input$aov_vars2,
        #         ",",
        #         input$aov_vars3,
        #         ")"
        #       )
        #     )
        #   
        #   print(
        #     gplots::plotmeans(
        #       fm,
        #       dat_long_tech_avg,
        #       mean.labels = T,
        #       digits = 1
        #     ),
        #     layout = FALSE
        #   )
        # }
        # 
        
        
      } else if (input$stats_test == "anova") {
        
        #dev.new()
        
        #library(gplots)
        #library(rstatix)
        
        #m=input$stats_table_rows_selected
        #print(p1, layout=n2mfrow(nplots))
        res.test <- x5$test_result[[m]]
        print(rstatix::get_anova_table(res.test))
        dat <-
          as.data.frame(subset(
            na.omit(x5$dat_long_tech_avg),
            mz == x5$stats_results[m, "mz"]
          ))
        
        
        
        
        if (input$anova_type == "anova_oneway") {
          if (!input$aov_vars1 %in% colnames(dat)) {
            print("Check the Group variable, different from comparison variable")
            return()
          }
          
          
          dat[, input$aov_vars1] <-
            as.factor(dat[, input$aov_vars1])
          
          
          library(ggpubr)
          bxp <- ggboxplot(
            dat,
            x = input$aov_vars1,
            y = 'tech_avg',
            #color = input$aov_vars1, add="jitter", palette = "jco"
            color = input$aov_vars1,
            add = "jitter",
            palette = "npg",
            size = 1.3,
            font.label = list(size = 15, color = "black")
          )
          #pair wise comparision
          pwc <-
            (dat) %>% rstatix::tukey_hsd(attributes(res.test)$args$formula)
          pwc <-
            (pwc) %>% rstatix::add_xy_position(x = eval(input$aov_vars1))
          
          print(as.data.frame(pwc))
          
          
          
          p1 <- bxp +
            stat_pvalue_manual(pwc,
                               tip.length = 0,
                               hide.ns = TRUE) +
            labs(
              title = paste0(
                "m/z= ",
                round(x5$stats_results[m, "mz"], 4),
                ",  FDR= ",
                x5$stats_results[m, "adjustp"]
              ),
              subtitle = rstatix::get_test_label(res.test, detailed = TRUE) ,
              caption = rstatix::get_pwc_label(pwc),
              size = 15
            )
          
          print(p1)
          
        } else if (input$anova_type %in% c("anova_add",
                                           "anova_interaction",
                                           "between_interacting")) {
          
          
          #library(ggpubr)
          #library(ggprism)
          
          dat[, input$aov_vars1] <-
            as.factor(dat[, input$aov_vars1])
          dat[, input$aov_vars2] <-
            as.factor(dat[, input$aov_vars2])
          
          
          if (input$anova_type %in% c("between_interacting")) {
            dat[, input$aov_vars3] <- as.factor(dat[, input$aov_vars3])
          }
          
          
          
          
          
          
          plot_means_p <- function(dat, v1, v2, cols) {
            fm1 = as.formula(paste0("tech_avg ~ ", v2))
            fm2 = as.formula(paste0("tech_avg ~ ", v1))
            
            
            
            
            df_p_val1 <- (dat) %>%
              rstatix::group_by(.data[[v1]]) %>%
              as.data.frame() %>%
              rstatix::t_test(fm1) %>%
              rstatix::adjust_pvalue(p.col = "p", method = "fdr") %>%
              rstatix::add_significance(p.col = "p.adj") %>%
              rstatix::add_xy_position(x = v1, dodge = 0.8) # important for positioning!
            
            df_p_val2 <- rstatix::t_test(dat, fm2,
                                         p.adjust.method = "fdr") %>%
              rstatix::add_xy_position()
            
            if (is.null(df_p_val2$p.adj)) {
              df_p_val2 <- df_p_val2 %>% dplyr::mutate(p.adj = p)
            }
            
            
            # bxp <- ggboxplot(
            #   dat, x = v1, y = 'tech_avg',
            #   color = v2, add="jitter", palette = "jco"
            #                 )+theme_prism()
            
            p1 <-
              ggplot2::ggplot(dat, ggplot2::aes(x = .data[[v1]], y = tech_avg)) +
              ggplot2::geom_boxplot(ggplot2::aes(fill = .data[[v2]])) +
              ggplot2::scale_fill_manual(values = as.vector(cols)) +
              ggprism::theme_prism(base_family = "Arial") +
              
              ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25)) + #, hjust = 1))+
              ggplot2::ylab("response") +
              ggplot2::coord_cartesian(ylim = c(0, max(df_p_val2$y.position *
                                                         1.2))) +
              ggprism::add_pvalue(
                df_p_val1,
                xmin = "xmin",
                xmax = "xmax",
                label = "p = {p.adj}",
                tip.length = 0
              ) +
              ggprism::add_pvalue(
                df_p_val2,
                label = "p = {p.adj}",
                tip.length = 0.01,
                bracket.nudge.y = df_p_val2$y.position * .1,
                step.increase = 0.015
              )
            
            
            
            return(p1)
          }
          
          #dat<-dat %>% filter(!Plate.Group=="N9.healthy control")
          
          dat <- droplevels(dat)
          
          match_cols <-
            ggsci::pal_npg()(length(levels(droplevels(dat[, input$aov_vars1]))))
          contrast_cols <-
            pals::alphabet2((length(levels(
              droplevels(dat[, input$aov_vars2])
            ))))
          
          p1 <-
            plot_means_p(
              dat,
              v1 = input$aov_vars1,
              v2 = input$aov_vars2,
              cols = contrast_cols
            )
          p2 <-
            plot_means_p(
              dat,
              v1 = input$aov_vars2,
              v2 = input$aov_vars1,
              cols = match_cols
            )
          
          
          
          
          
          #add anova plots? qq plots, anything else?
          sig_cols <-
            grep("BH|mz|ID",
                 colnames(x5$stats_results[m, ]),
                 value = T)
          sig_vals <-
            x5$stats_results[m, c(colnames(x5$stats_results) %in% sig_cols)]
          sig_vals$mz <- round(sig_vals$mz, 4)
          
          table_vals <-
            table(interaction(dat[, input$aov_vars1], dat[, input$aov_vars2]))
          
          print(
            gridExtra::grid.arrange(
              p1,
              p2,
              gridExtra::tableGrob(sig_vals),
              gridExtra::tableGrob(t(as.data.frame(
                table_vals
              ))),
              nrow = 3,
              heights = grid::unit(c(5, 1, 1), c("in", "in")),
              layout_matrix = rbind(c(1, 1, 2, 2), c(NA, 3, 3, NA), c(NA, 4, 4, NA))
            ),
          )
          
          
          
          #fm3<-as.formula(paste0("tech_avg~interaction(",input$aov_vars1,",",input$aov_vars2,")"))
          #plotmeans(formula=fm3, data=dat)
          
        } else if (input$anova_type %in% c("anova_in")) {
          #library(ggpubr)
          
          
          
          bxp <- ggpubr::ggboxplot(
            dat,
            x = input$aov_vars2,
            y = 'tech_avg',
            color = input$aov_vars1,
            add = "jitter",
            palette = "jco"
          )
          qq <-
            ggpubr::ggqqplot(dat, "tech_avg", ggtheme = theme_bw()) +
            ggplot2::facet_grid(paste0(input$aov_vars2, " ~ ", input$aov_vars1))
          
          #library(emmeans)
          fm2 <-
            as.formula(paste0("tech_avg ~ ", input$aov_vars1))
          pwc <-
            dat %>% dplyr::group_by(.data[[input$aov_vars2]]) %>%
            #group_by( dat$time.point) %>%
            rstatix::emmeans_test(fm2, p.adjust.method = "fdr")
          pwc
          
          pwc <-
            pwc %>% rstatix::add_xy_position(x = eval(input$aov_vars1)) #TODO check
          
          
          p1 <- bxp +
            ggpubr::stat_pvalue_manual(pwc) +
            ggplot2::labs(
              subtitle = rstatix::get_test_label(res.test, detailed = TRUE),
              caption = rstatix::get_pwc_label(pwc)
            )
        }
        
        
        
        
      } else if (input$stats_test == "ssctest") {
        
        
        
        
        #require(gplots)
        #mzdat=subset(x5$test_result, mz==x5$stats_results$mz[m])
        a <-
          as.matrix(spectra(
            subsetFeatures(x5$data_file_selected,
                           mz == x5$stats_results$mz[m])
          ))
        b <-
          pData(subsetFeatures(x5$data_file_selected,
                               mz == x5$stats_results$mz[m]))
        dat = cbind(run = Cardinal::run(x5$data_file_selected),
                    b,
                    value = a)
        
        dat_long_tech_avg <-
          as.data.frame(dat) %>% dplyr::group_by_at(unique(
            c(
              "run",
              input$phen_cols_stats
            )
          )) %>%
          dplyr::summarize(tech_avg = mean(value),
                           .groups = "keep")
        
        #x5$dat_long_tech_avg<-dat_long_tech_avg #single m/z
        
        
        fm <-
          as.formula(
            paste0(
              "tech_avg~(",
              input$phen_cols_stats,")"))
        
        #choose one model to plot, hightest accuracy
        idx<-input$stats_table_rows_selected
        ssc_model <- x5$stats_table_filtered[idx, "model"]
        
        ncolors <- length(levels(as.factor(as.data.frame(
          pData(x5$data_file_selected)[, input$phen_cols_stats]
        )[,1])))
        
        
        mycols = ggsci::pal_npg()(ncolors)
        
        
        
        if(is.null(names(x5$test_result))){
          p1 <-
            image(x5$test_result[[1]],
                  #model = ssc_model,
                  key = F,
                  col=mycols)
        } else {
          p1 <-
            image(x5$test_result,
                  model = ssc_model,
                  key = F,
                  col=mycols)
          
        }
        
        if(length(mz)>1){
          showNotification("Multiple ions selected, only first ion will be plotted", duration = 10)
          idx = idx[1]
        }
        
        p2 <-
          image(
            x5$data_file_selected,
            mz = (x5$stats_results$mz[idx]),
            model=ssc_model,
            
            enhance = "histogram",
            scale = T, free="xy"
            #col=mycols
          )
        
        fm2 <-
          (paste0(input$phen_cols_stats )) #, "~x*y"))
        
        #if more than one model present, create a matter vizi facet with all models
        
       if(class(x5$test_result)=="SpatialShrunkenCentroids") {
         x5$test_result<-list(x5$test_result)
       }
        
        if(length(x5$test_result)>1){
          plot_list <- list()
          nmodels=length(x5$test_result)
          for(i in 1:nmodels){
            plot_list[[i]] <-
              plot(x5$test_result[[i]], type="statistic", linewidth=2, col=mycols)
            
          }
          names(plot_list) <- names(x5$test_result)
        } else {
          plot_list <- plot(x5$test_result[[1]], type="statistic", linewidth=2, col=mycols)
        }
        p3 <-
          matter::as_facets(plot_list, ncol = 1)
        
        p4 <-
          image(x5$data_file_selected[, Cardinal::run(x5$data_file_selected) %in% runNames(x5$data_file_selected)[]],
                fm2,
                key = T,
                col = mycols)
        
        plot_out<-switch(
          input$plot_choice,
          "ion_image" = p1,
          "means_plot" = gplots::plotmeans(fm, dat_long_tech_avg),
          "t_statistic" = p3,
          "groupings" = p4,
          "msi_image" = p2,
        )
        
        print(plot_out)
        
      } else{
        print("Plot not supported yet")
      }
      
      
      
      
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
    
    proc_values <- reactive({
      list(x5 = x5,
           par_mode = par_mode)
    })
    return(proc_values)
    
    
  })
}
