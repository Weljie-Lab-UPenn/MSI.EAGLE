### R/PeakPickServer.R

PeakPickServer <- function(id, setup_values) {
  moduleServer(id, function(input, output, session) {
    #for dynamic UI
    ns = session$ns
    
    output$variables <- renderPrint({
      setup_vals <- setup_values()
      paste("working directory=", setup_vals$wd)
    })
    
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

    listed_files <- reactive({
      files <- my_files()
      files[!is.na(files) & nzchar(files)]
    })
    is_rds_file <- function(path) {
      grepl("\\.(rds|RData)$", basename(path), ignore.case = TRUE)
    }
    is_imzml_file <- function(path) {
      grepl("\\.imzML$", basename(path), ignore.case = TRUE)
    }
    rds_files <- reactive({
      files <- listed_files()
      files[is_rds_file(files)]
    })
    imzml_or_rds_files <- reactive({
      files <- listed_files()
      files[is_rds_file(files) | is_imzml_file(files)]
    })
    resolve_wd_path <- function(path) {
      if (is.null(path) || !nzchar(path)) return(path)
      if (grepl("^(/|[A-Za-z]:[\\\\/])", path)) return(path)
      file.path(setup_values()[["wd"]], path)
    }
    log_file_load <- function(context, selected_name, resolved_path, obj = NULL) {
      message(sprintf(
        "[PeakPick] Loaded %s | selected='%s' | path='%s'",
        context,
        as.character(selected_name),
        as.character(resolved_path)
      ))
      if (!is.null(obj) && !inherits(obj, "try-error")) {
        n_feat <- try(nrow(obj), silent = TRUE)
        n_pix <- try(ncol(obj), silent = TRUE)
        if (!inherits(n_feat, "try-error") && !inherits(n_pix, "try-error")) {
          message(sprintf("[PeakPick] Summary %s: features=%s, pixels=%s", context, as.character(n_feat), as.character(n_pix)))
        }
      }
    }
    
    # Detect common Cardinal/matter serialization mismatches early so we can fail gracefully.
    validate_msi_object <- function(obj, context = "dataset", notify = TRUE) {
      if (is.null(obj)) {
        if (notify) {
          showNotification(sprintf("Could not load %s: object is NULL.", context), type = "error", duration = 10)
        }
        return(FALSE)
      }
      npx <- try(ncol(obj), silent = TRUE)
      if (inherits(npx, "try-error") || is.null(npx) || length(npx) != 1 || !is.finite(npx) || npx < 1) {
        if (notify) {
          showNotification(
            sprintf("Could not use %s: object does not look like a valid MSI imaging dataset.", context),
            type = "error",
            duration = 10
          )
        }
        return(FALSE)
      }
      test_idx <- rep(FALSE, npx)
      test_idx[1] <- TRUE
      test_subset <- try(Cardinal::subsetPixels(obj, test_idx), silent = TRUE)
      if (inherits(test_subset, "try-error")) {
        err_txt <- as.character(test_subset)
        hint <- if (grepl("slot of name \"refs\"|class \"atoms\"", err_txt, ignore.case = TRUE)) {
          "This file appears incompatible with the installed Cardinal/matter versions. Try loading the source .imzML or re-saving the .rds in this environment."
        } else {
          "Please try the source .imzML file or regenerate the .rds."
        }
        if (notify) {
          showNotification(
            sprintf("Failed to load %s for pixel operations. %s", context, hint),
            type = "error",
            duration = 12
          )
        }
        message("Validation failure in PeakPickServer for ", context, ": ", err_txt)
        return(FALSE)
      }
      TRUE
    }

    parse_feature_index_selection <- function(selection_txt, n_features) {
      if (is.null(selection_txt) || !nzchar(trimws(selection_txt))) {
        return(seq_len(n_features))
      }
      tokens <- unlist(strsplit(trimws(selection_txt), "[,\\s]+"))
      tokens <- tokens[nzchar(tokens)]
      if (length(tokens) == 0) {
        return(seq_len(n_features))
      }

      idx <- integer(0)
      invalid <- character(0)
      for (tok in tokens) {
        if (grepl("^\\d+$", tok)) {
          idx <- c(idx, as.integer(tok))
        } else if (grepl("^\\d+:\\d+$", tok)) {
          bounds <- as.integer(strsplit(tok, ":", fixed = TRUE)[[1]])
          idx <- c(idx, seq.int(bounds[1], bounds[2]))
        } else {
          invalid <- c(invalid, tok)
        }
      }

      if (length(invalid) > 0) {
        stop(sprintf(
          "Invalid feature selector token(s): %s. Use i indices like 1, 5, 20:40.",
          paste(unique(invalid), collapse = ", ")
        ))
      }

      idx <- sort(unique(idx))
      idx <- idx[is.finite(idx) & idx >= 1 & idx <= n_features]
      if (length(idx) == 0) {
        stop(sprintf("No valid feature indices remain after filtering to 1..%d.", n_features))
      }
      idx
    }

    combine_with_cardinal_fallback <- function(img_list, context = "combine") {
      combined <- try(Cardinal::combine(img_list), silent = TRUE)
      if (!inherits(combined, "try-error")) {
        message(sprintf("[PeakPick] %s: using Cardinal::combine", context))
        return(combined)
      }
      message(sprintf(
        "[PeakPick] %s: Cardinal::combine failed; using combine_card fallback. Error: %s",
        context, as.character(combined)
      ))
      combine_card(img_list)
    }
    
    # any time the reactive changes, update the selectInput
    observeEvent(my_files(),
                 ignoreInit = T,
                 ignoreNULL = T,
                 {
                   updateSelectInput(
                     session,
                     ns('peakPickfile'),
                     choices = imzml_or_rds_files()
                   )
                 })
    
    
    #change depending on value of input$peak_pick_status, pp_y is .imzML, pp_old is .rds
    output$pk_file <- renderUI({
      
      # Determine the list of files based on the value of input$peak_pick_status
      file_choices <- switch(
        input$peak_pick_status,
        "pp_old" = rds_files(),
        "pp_y" = imzml_or_rds_files(),
        "pp_ref" = imzml_or_rds_files(),
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
        imzml_or_rds_files()
      )
    })
    
    
    
    output$operation_help <- renderUI({
      req(input$pp_operation)
      switch(
        input$pp_operation,
        "open_file" = tags$div(
          tags$strong("Open existing dataset"),
          tags$ul(
            tags$li("Restore a peak-picked .imzML/.rds dataset."),
            tags$li("Edit xmin/xmax/ymin/ymax in the run table to crop per run."),
            tags$li("Save exports the cropped runs plus optional feature index filter (i values).")
          )
        ),
        "pp_raw" = tags$div(
          tags$strong("Peak-pick raw files"),
          tags$ul(
            tags$li("Start from raw files in Data Setup, from a reference file, or from mean-spectrum peaks."),
            tags$li("After processing, use the run table limits to crop."),
            tags$li("Optional feature index filter applies on save (i values).")
          )
        ),
        "subset_f" = tags$div(
          tags$strong("Combine by coordinates"),
          tags$ul(
            tags$li("Use when two datasets share the same coordinate grid but have different peak lists."),
            tags$li("The merged feature list can be filtered by i indices at save."),
            tags$li("Run-table limits still define cropped export boundaries.")
          )
        ),
        "add_same_pklist" = tags$div(
          tags$strong("Combine by peak list"),
          tags$ul(
            tags$li("Use when datasets have the same mz/feature list and you want to append runs/pixels."),
            tags$li("Optional feature index filter (i) is available at save."),
            tags$li("Run-table x/y limits define cropped output.")
          )
        )
      )
    })

    output$peak_pick <- renderUI({
      feature_select_ui <- tagList(
        textInput(
          ns("selected_mz"),
          "Optional feature index filter for save (i values)",
          placeholder = "Examples: 1,5,20:40 or 1 5 20:40. Leave blank to save all features."
        ),
        helpText("Applies to save in all modes, including single-file open/restore.")
      )
      save_button_ui <- shinyFiles::shinySaveButton(
        ns("save_imzml"),
        "Save cropped / filtered imzML",
        "Save",
        filetype = list("")
      )

      switch(
        input$pp_operation,
        "open_file" = list(
          radioButtons(
            ns("peak_pick_status"),
            "Cardinal v3.6+ processed .imzML/.rds file?",
            c("Yes" = "pp_y", "No, older .rds" = "pp_old")
          ),
          uiOutput(ns('pk_file')),
          actionButton(ns("action"), label = HTML("Restore saved file")),
          feature_select_ui,
          save_button_ui
        ),
        "pp_raw" = list(
          radioButtons(
            ns("peak_pick_status"),
            "Start from an existing peak-picked file?",
            c(
              "No, fresh pick raw files" = "pp_no",
              "From reference file" = "pp_ref",
              "From mean spectrum" = "pp_mean"
            )
          ),
          uiOutput(ns('pk_file_condition')),
          feature_select_ui,
          p(),
          save_button_ui
        ),
        "subset_f" = list(
          uiOutput(ns('pk_file')),
          actionButton(ns("action"), label = HTML("Restore first file")),
          uiOutput(ns('add_file')),
          actionButton(ns('action_add_file'), "Add file (same coordinates)"),
          feature_select_ui,
          save_button_ui
        ),
        "add_same_pklist" = list(
          uiOutput(ns('pk_file')),
          actionButton(ns("action"), label = HTML("Restore first file")),
          uiOutput(ns('add_file')),
          actionButton(
            ns('action_add_file_same'),
            "Add imageset (same peak list)"
          ),
          feature_select_ui,
          save_button_ui
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
      if (is.null(input$cardworkdat) || !nzchar(input$cardworkdat)) {
        showNotification("Choose a demo dataset before loading demo data.", type = "warning", duration = 7)
        return()
      }
      if (!requireNamespace("CardinalWorkflows", quietly = TRUE)) {
        showNotification("CardinalWorkflows package is required to load demo data.", type = "error", duration = 8)
        return()
      }
      
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
                         showNotification("Choose an input file first.", type = "warning", duration = 6)
                         return()
                         
                       }
                       
                       
                       pick_path <- resolve_wd_path(input$peakPickfile)
                       #check filetype and open accordingly
                       if (is_rds_file(input$peakPickfile)) {
                         message(sprintf("[PeakPick] Loading main dataset from .rds | selected='%s' | path='%s'", input$peakPickfile, pick_path))
                         overview_peaks <- try(readRDS(pick_path), silent = TRUE)
                         print(overview_peaks)
                       } else if (is_imzml_file(input$peakPickfile)) {
                         message(sprintf("[PeakPick] Loading main dataset from .imzML | selected='%s' | path='%s'", input$peakPickfile, pick_path))
                         overview_peaks <- try(readImzML(pick_path), silent = TRUE)
                         print(overview_peaks)
                       } else {
                         print("file type not recognized")
                         showNotification("Input file type not recognized. Use .imzML or .rds/.RData.", type = "error", duration = 8)
                         return()
                       }
                       if (inherits(overview_peaks, "try-error")) {
                         showNotification("Failed to read selected file. Try the source .imzML file if this is an .rds compatibility issue.", type = "error", duration = 10)
                         return()
                       }
                       if (!validate_msi_object(overview_peaks, context = basename(input$peakPickfile), notify = TRUE)) {
                         return()
                       }
                       log_file_load("main dataset", input$peakPickfile, pick_path, overview_peaks)
                       
                       # overview_peaks <- readImzML(input$peakPickfile)
                       # print(overview_peaks)
                       
                       #browser()
                       
                       #De novo peak picking
                     }else if (input$peak_pick_status == "pp_old") {
                       if (length(input$peakPickfile) < 1) {
                         print("need file to open")
                         showNotification("Choose an .rds file first.", type = "warning", duration = 6)
                         return()
                         
                       }
                       
                       pick_path <- resolve_wd_path(input$peakPickfile)
                       message(sprintf("[PeakPick] Loading legacy dataset from .rds | selected='%s' | path='%s'", input$peakPickfile, pick_path))
                       old_dat <- try(readRDS(pick_path), silent = TRUE)
                       if (inherits(old_dat, "try-error")) {
                         showNotification("Failed to read legacy .rds file. Please verify the file or use .imzML input.", type = "error", duration = 10)
                         return()
                       }
                       overview_peaks <- try(convert_card(old_dat), silent = TRUE)
                       if (inherits(overview_peaks, "try-error")) {
                         showNotification("Failed to convert legacy .rds file. Try loading from .imzML instead.", type = "error", duration = 10)
                         return()
                       }
                       if (!validate_msi_object(overview_peaks, context = basename(input$peakPickfile), notify = TRUE)) {
                         return()
                       }
                       log_file_load("legacy converted dataset", input$peakPickfile, pick_path, overview_peaks)
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
                                                                                              mass.range = c(setup_values()[["mz_min"]], setup_values()[["mz_max"]]), 
                                                                                              #x1$mass.range,
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
                         showNotification("Choose a reference file before starting peak binning.", type = "warning", duration = 6)
                         return()
                       }
                       #browser()
                       print("binning raw files from reference")
                       #check for rds or imzml and open accordingly
                       pick_path <- resolve_wd_path(input$peakPickfile)
                       if (is_rds_file(input$peakPickfile)) {
                         message(sprintf("[PeakPick] Loading reference peak list from .rds | selected='%s' | path='%s'", input$peakPickfile, pick_path))
                         test_mz_reduced <- readRDS(pick_path)
                       } else if (is_imzml_file(input$peakPickfile)) {
                         message(sprintf("[PeakPick] Loading reference peak list from .imzML | selected='%s' | path='%s'", input$peakPickfile, pick_path))
                         test_mz_reduced <- readImzML(pick_path)
                       } else {
                         print("file type not recognized")
                         showNotification("Reference file type not recognized. Use .imzML or .rds/.RData.", type = "error", duration = 8)
                         return()
                       }
                       log_file_load("reference peak list", input$peakPickfile, pick_path, test_mz_reduced)
                       
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
                             mass.range=c(setup_values()[["mz_min"]], setup_values()[["mz_max"]])) %>% 
                             #mass.range=x1$mass.range) %>%
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
      if (is.null(input$peakAddfile) || !nzchar(input$peakAddfile)) {
        showNotification("Choose a file to add first.", type = "warning", duration = 6)
        return()
      }
      if (is.null(x0$overview_peaks)) {
        showNotification("No base dataset loaded. Restore or create a peak-picked dataset first.", type = "warning", duration = 7)
        return()
      }
      #check for .rds or .imzML and open accordingly
      add_path <- resolve_wd_path(input$peakAddfile)
      if (is_rds_file(input$peakAddfile)) {
        message(sprintf("[PeakPick] Loading file-to-add from .rds | selected='%s' | path='%s'", input$peakAddfile, add_path))
        add_peaks <- try(readRDS(add_path), silent = TRUE)
        print(add_peaks)
      } else if (is_imzml_file(input$peakAddfile)) {
        message(sprintf("[PeakPick] Loading file-to-add from .imzML | selected='%s' | path='%s'", input$peakAddfile, add_path))
        add_peaks <- try(readImzML(add_path), silent = TRUE)
        print(add_peaks)
      } else {
        print("file type not recognized")
        showNotification("Input file type not recognized. Use .imzML or .rds/.RData.", type = "error", duration = 8)
        return()
      }
      if (inherits(add_peaks, "try-error")) {
        showNotification("Could not read file to add. Try .imzML if this .rds is not compatible.", type = "error", duration = 10)
        return()
      }
      if (!validate_msi_object(add_peaks, context = basename(input$peakAddfile), notify = TRUE)) {
        return()
      }
      log_file_load("file-to-add", input$peakAddfile, add_path, add_peaks)
      # add_peaks <- readImzML(input$peakAddfile)
      # print(add_peaks)
      # 
      
      #check to make sure overview_peaks() exists
      
      if (dim(x0$overview_peaks)[2] != dim(add_peaks)[2]) {
        print("inconistent pixel sizes, make sure same coordinate set being used")
        showNotification("Pixel dimensions do not match. Use files with the same coordinate set.", type = "error", duration = 8)
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
      if (is.null(x0$overview_peaks)) {
        showNotification("No base dataset loaded. Restore or create a peak-picked dataset first.", type = "warning", duration = 7)
        return()
      }
      if (input$peakAddfile == '') {
        print("need file to open")
        showNotification("Choose a file to add first.", type = "warning", duration = 6)
        return()
        
      }
      #check filename and open accordinlgy
      add_path <- resolve_wd_path(input$peakAddfile)
      if (is_rds_file(input$peakAddfile)) {
        message(sprintf("[PeakPick] Loading file-to-add (same peaklist) from .rds | selected='%s' | path='%s'", input$peakAddfile, add_path))
        add_peaks <- try(readRDS(add_path), silent = TRUE)
        print(add_peaks)
      } else if (is_imzml_file(input$peakAddfile)) {
        message(sprintf("[PeakPick] Loading file-to-add (same peaklist) from .imzML | selected='%s' | path='%s'", input$peakAddfile, add_path))
        add_peaks <- try(readImzML(add_path), silent = TRUE)
        print(add_peaks)
      } else {
        print("file type not recognized")
        showNotification("Input file type not recognized. Use .imzML or .rds/.RData.", type = "error", duration = 8)
        return()
      }
      if (inherits(add_peaks, "try-error")) {
        showNotification("Could not read file to add. Try .imzML if this .rds is not compatible.", type = "error", duration = 10)
        return()
      }
      if (!validate_msi_object(add_peaks, context = basename(input$peakAddfile), notify = TRUE)) {
        return()
      }
      log_file_load("file-to-add same-peaklist", input$peakAddfile, add_path, add_peaks)
      
      #browser()
      # add_peaks <- readImzML(input$peakAddfile)
      # print(add_peaks)
      # 
      
      #check to make sure overview_peaks() exists and has the same number of peaks
      
      if (!identical(mz(x0$overview_peaks), mz(add_peaks))) {
        print("inconistent mass lists, make sure same peak lists being used")
        showNotification("Mass lists differ between files. Use datasets with the same peak list for this option.", type = "error", duration = 8)
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
      
      dat <- try(combine_with_cardinal_fallback(list(x, y), context = "same-peaklist combine"), silent = TRUE)
      if (inherits(dat, "try-error")) {
        showNotification("Combining same-peaklist datasets failed. Check file compatibility and run metadata.", type = "error", duration = 10)
        message("PeakPick same-peaklist combine failure: ", as.character(dat))
        return()
      }

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
      a <- runNames(isolate(x0$overview_peaks))
      cd <- try(as.data.frame(coord(x0$overview_peaks)), silent = TRUE)
      runs_px <- try(as.character(Cardinal::run(x0$overview_peaks)), silent = TRUE)
      if (inherits(cd, "try-error") || inherits(runs_px, "try-error") || is.null(cd$x) || is.null(cd$y)) {
        showNotification("Could not parse coordinates from loaded dataset. Try loading the source .imzML file.", type = "error", duration = 10)
        x2$run_table <- NULL
        return()
      }
      if (length(runs_px) != nrow(cd)) {
        showNotification("Loaded dataset has mismatched run/coordinate lengths. Try loading the source .imzML file.", type = "error", duration = 10)
        x2$run_table <- NULL
        return()
      }
      xmax <- numeric(length(a))
      ymax <- numeric(length(a))
      for (i in seq_along(a)) {
        idx <- runs_px %in% a[i]
        if (!any(idx)) {
          xmax[i] <- 0
          ymax[i] <- 0
        } else {
          xmax[i] <- suppressWarnings(max(cd$x[idx], na.rm = TRUE))
          ymax[i] <- suppressWarnings(max(cd$y[idx], na.rm = TRUE))
          if (!is.finite(xmax[i])) xmax[i] <- 0
          if (!is.finite(ymax[i])) ymax[i] <- 0
        }
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
        "<b>Run table:</b> select runs to export. Edit xmin/xmax/ymin/ymax to crop each run before saving."
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
      "Tip: saved output uses current run selections, x/y crop bounds from this table, and optional feature index filter (i values)."
    })
    
    
    observe({
      req(x0$overview_peaks)
      req(x2$run_table)
      a <- runNames((x0$overview_peaks))
      cd <- try(as.data.frame(coord(x0$overview_peaks)), silent = TRUE)
      runs_px <- try(as.character(Cardinal::run(x0$overview_peaks)), silent = TRUE)
      if (inherits(cd, "try-error") || inherits(runs_px, "try-error") || is.null(cd$x) || is.null(cd$y)) {
        showNotification("Could not subset loaded dataset for plotting. Try loading the source .imzML file.", type = "error", duration = 10)
        x2$overview_peaks_sel <- NULL
        return()
      }
      keep <- rep(FALSE, nrow(cd))
      for (i in seq_along(a)) {
        rt_i <- match(a[i], x2$run_table$runs)
        if (is.na(rt_i)) next
        x_min <- suppressWarnings(as.numeric(x2$run_table$xmin[rt_i]))
        x_max <- suppressWarnings(as.numeric(x2$run_table$xmax[rt_i]))
        y_min <- suppressWarnings(as.numeric(x2$run_table$ymin[rt_i]))
        y_max <- suppressWarnings(as.numeric(x2$run_table$ymax[rt_i]))
        if (!is.finite(x_min)) x_min <- 0
        if (!is.finite(x_max)) x_max <- max(cd$x[runs_px %in% a[i]], na.rm = TRUE)
        if (!is.finite(y_min)) y_min <- 0
        if (!is.finite(y_max)) y_max <- max(cd$y[runs_px %in% a[i]], na.rm = TRUE)
        keep <- keep | (
          runs_px %in% a[i] &
            cd$x > x_min &
            cd$x < (x_max + 1) &
            cd$y > y_min &
            cd$y < (y_max + 1)
        )
      }
      if (!any(keep)) {
        showNotification("No pixels remain after current run-boundary filters.", type = "warning", duration = 8)
        x2$overview_peaks_sel <- NULL
        return()
      }
      selected_img <- try(Cardinal::subsetPixels(x0$overview_peaks, keep), silent = TRUE)
      if (inherits(selected_img, "try-error")) {
        err_txt <- as.character(selected_img)
        showNotification("Failed to subset loaded dataset for plotting. Try loading the source .imzML file.", type = "error", duration = 10)
        message("PeakPick subset failure while building overview selection: ", err_txt)
        x2$overview_peaks_sel <- NULL
        return()
      }
      
      #check if freq data is getting lost here
      #browser()
      
      #do.call doesn't always work, so lets use the loop to recombine the data
      
      
      x2$overview_peaks_sel <- selected_img
      

      plot_card_server("card_plot", overview_peaks_sel = x2$overview_peaks_sel)
      
      
    })
    
    
    # Save current cropped dataset (run-table x/y limits), selected runs, and optional feature indices.
    observeEvent(input$save_imzml, {
      req(input$save_imzml)
      if (is.null(x2$overview_peaks_sel) || is.null(x0$overview_peaks)) {
        showNotification("No peak-picked dataset available to save. Load or create a dataset first.", type = "warning", duration = 8)
        return(NULL)
      }
      if (is.null(input$peak_pick_selection_rows_selected) || length(input$peak_pick_selection_rows_selected) == 0) {
        showNotification("Select at least one run in the peak-pick table before saving.", type = "warning", duration = 8)
        return(NULL)
      }

      volumes <- c(wd = setup_values()[["wd"]], home = fs::path_home())
      shinyFiles::shinyFileSave(input, "save_imzml", roots = volumes, session = session)
      save_path <- shinyFiles::parseSavePath(volumes, input$save_imzml)
      if (nrow(save_path) == 0) return(NULL)

      filen <- as.character(save_path$datapath)

      # restore full feature metadata before optional feature filtering
      fData(x2$overview_peaks_sel) <- fData(x0$overview_peaks)

      features_sel <- tryCatch(
        parse_feature_index_selection(input$selected_mz, nrow(x2$overview_peaks_sel)),
        error = function(e) {
          showNotification(conditionMessage(e), type = "error", duration = 10)
          return(NULL)
        }
      )
      if (is.null(features_sel)) return(NULL)

      run_keep <- unique(as.character(x2$run_table$runs[input$peak_pick_selection_rows_selected]))
      run_keep <- run_keep[nzchar(run_keep)]
      if (length(run_keep) == 0) {
        showNotification("No valid runs selected to save.", type = "warning", duration = 8)
        return(NULL)
      }

      pk_img <- try({
        x2$overview_peaks_sel[features_sel, ] %>%
          subsetPixels(Cardinal::run(x2$overview_peaks_sel) %in% run_keep)
      }, silent = TRUE)
      if (inherits(pk_img, "try-error")) {
        showNotification("Could not subset dataset for save. This file may be incompatible; try exporting from .imzML source.", type = "error", duration = 10)
        message("PeakPick save subset failure: ", as.character(pk_img))
        return(NULL)
      }

      writeImzML(pk_img, filen)
      saveRDS(pk_img, paste0(filen, ".rds"))

      n_runs_saved <- length(unique(as.character(Cardinal::run(pk_img))))
      showNotification(
        sprintf("Saved cropped dataset: %d run(s), %d feature(s), %d pixel(s).", n_runs_saved, nrow(pk_img), ncol(pk_img)),
        type = "message",
        duration = 8
      )
    })
    
    
    
    
    #browser()
    
    preproc_values <- reactive({
      list(x2 = x2,
           x0 = x0)
    })
    return(preproc_values)
    
  })
}
