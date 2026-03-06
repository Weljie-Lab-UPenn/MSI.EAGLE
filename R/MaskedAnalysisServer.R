### R/MaskedAnalysisServer.R
MaskedAnalysisServer <- function(id,  setup_values) {
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

    segmentation_file_choices <- reactive({
      files <- my_files()
      files <- files[!is.na(files) & nzchar(files)]
      files[grepl("\\.(imzML|rds)$", basename(files), ignore.case = TRUE)]
    })
    resolve_wd_path <- function(path) {
      if (is.null(path) || !nzchar(path)) return(path)
      if (grepl("^(/|[A-Za-z]:[\\\\/])", path)) return(path)
      file.path(setup_values()[["wd"]], path)
    }

    extract_imzml_template_pdata <- function(seg_path) {
      parsed <- try(CardinalIO::parseImzML(seg_path), silent = TRUE)
      if (inherits(parsed, "try-error")) {
        stop(sprintf("parseImzML failed: %s", as.character(parsed)))
      }

      pos <- try(as.data.frame(parsed$run$spectrumList$positions, stringsAsFactors = FALSE), silent = TRUE)
      if (inherits(pos, "try-error") || nrow(pos) == 0L) {
        stop("No position metadata found in imzML run/spectrumList.")
      }

      pos_names <- names(pos)
      x_col <- pos_names[grep("position\\s*x|^x$", pos_names, ignore.case = TRUE)][1]
      y_col <- pos_names[grep("position\\s*y|^y$", pos_names, ignore.case = TRUE)][1]
      z_col <- pos_names[grep("position\\s*z|^z$", pos_names, ignore.case = TRUE)][1]
      if (!is.character(x_col) || !nzchar(x_col) || !is.character(y_col) || !nzchar(y_col)) {
        stop("Could not identify X/Y coordinate columns in imzML position metadata.")
      }

      coord_df <- data.frame(
        x = suppressWarnings(as.integer(pos[[x_col]])),
        y = suppressWarnings(as.integer(pos[[y_col]])),
        stringsAsFactors = FALSE
      )
      if (is.character(z_col) && nzchar(z_col)) {
        z_vals <- suppressWarnings(as.integer(pos[[z_col]]))
        if (any(is.finite(z_vals))) {
          coord_df$z <- z_vals
        }
      }

      n_pix <- nrow(coord_df)
      base_run <- tools::file_path_sans_ext(basename(seg_path))
      run_vec <- rep(base_run, n_pix)

      pdata_sidecar_path <- sub("\\.imzML$", ".pdata", seg_path, ignore.case = TRUE)
      extra_df <- NULL
      if (file.exists(pdata_sidecar_path)) {
        sidecar <- try(read.delim(pdata_sidecar_path, sep = "", stringsAsFactors = FALSE, check.names = FALSE), silent = TRUE)
        if (!inherits(sidecar, "try-error") && nrow(sidecar) == n_pix) {
          sidecar <- as.data.frame(sidecar, stringsAsFactors = FALSE, check.names = FALSE)
          if (!is.null(rownames(sidecar)) && !is.null(rownames(pos)) &&
              all(rownames(pos) %in% rownames(sidecar))) {
            sidecar <- sidecar[rownames(pos), , drop = FALSE]
          }
          if ("run" %in% names(sidecar)) {
            rv <- as.character(sidecar$run)
            rv[is.na(rv) | !nzchar(rv)] <- base_run
            run_vec <- rv
          }
          drop_cols <- names(sidecar) %in% c("x", "y", "z", "run", x_col, y_col, z_col)
          extra_df <- sidecar[, !drop_cols, drop = FALSE]
        }
      }

      pdat_args <- list(coord = coord_df, run = run_vec, row.names = FALSE)
      if (!is.null(extra_df) && ncol(extra_df) > 0) {
        pdat_args <- c(pdat_args, as.list(extra_df))
      }
      pdat <- do.call(Cardinal::PositionDataFrame, pdat_args)
      pdat
    }

    extract_template_pdata <- function(seg_path) {
      if (grepl("\\.imzML$", basename(seg_path), ignore.case = TRUE)) {
        # Fast path: parse XML + .pdata metadata only; do not load spectra intensities.
        return(extract_imzml_template_pdata(seg_path))
      }
      if (grepl("\\.rds$", basename(seg_path), ignore.case = TRUE)) {
        obj <- readRDS(seg_path)
        if (inherits(obj, "MSImagingArrays")) {
          obj <- convertMSImagingArrays2Experiment(obj)
        }
        if (inherits(obj, "MSImagingExperiment")) {
          return(pData(obj))
        }
        if (inherits(obj, "PositionDataFrame")) {
          return(obj)
        }
        if (is.data.frame(obj) && all(c("x", "y") %in% names(obj))) {
          run_vec <- if ("run" %in% names(obj)) as.character(obj$run) else rep("run0", nrow(obj))
          extra_df <- obj[, !(names(obj) %in% c("x", "y", "z", "run")), drop = FALSE]
          coord_df <- obj[, c("x", "y", intersect("z", names(obj))), drop = FALSE]
          pdat_args <- list(coord = coord_df, run = run_vec, row.names = FALSE)
          if (ncol(extra_df) > 0) pdat_args <- c(pdat_args, as.list(extra_df))
          return(do.call(Cardinal::PositionDataFrame, pdat_args))
        }
        stop("Selected .rds file does not contain an MSI object or coordinate table.")
      }
      stop("Please select a valid .imzML or .rds file.")
    }
    
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
                     ns("segmentation_file"),
                     choices = segmentation_file_choices()
                   )
                 })
    
    #Setup and values for the Masked processing tab
    x4 <- reactiveValues(
      seg_file = NULL,
      seg_pdata = NULL,
      seg_filename = NULL,
      #to keep record of what file is used for coordinates
      masked_peaks = NULL,
      seg_pp_file = NULL
    )
    
    
    output$segmented_data <- renderUI({
      selectInput(
        ns("segmentation_file"),
        "Imageset with masked coordinates",
        segmentation_file_choices()
      )
    })
    
    output$pp_options <- renderUI({
      switch(
        input$targeted_pp,
        "untargeted" =
          #list(
            # numericInput(
            #   ns("pix_for_peak_picking2"),
            #   "% of pixels to for peak picking",
            #   100
            # ),
            list(
              fluidRow(
                column(6, 
                       numericInput(
                         ns("SNR2"),
                         "S/N for comprehensive peak picking",
                         20
                       )
                ),
                column(6, 
                       numericInput(
                         ns("freq_min2"),
                         "Minimum peak frequency (0-1",
                         0.03,
                         min = 0,
                         max = 1
                       )
                )
              ),
              fluidRow( column(6, numericInput(
                ns("tol2"),
                "tolerance",
                30
                )
              ),
              column(6, selectInput(
                ns("units2"),
                "units",
                choices=c("ppm", "mz"),
                selected="ppm"
                )
              )
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
            numericInput(ns("tol2"), "tolerance for peak binning (ppm)", 25)
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
      
      if (is.null(input$segmentation_file) || !nzchar(input$segmentation_file)) {
        showNotification("Choose a coordinate template file (.imzML or .rds) first.", type = "warning", duration = 7)
        return()
      }
      
      
      if(is.null(x1$raw_list)){
        message("raw files not selected, please choose from the Data Setup tab first!")
        showNotification("raw files not selected, please choose from the Data Setup tab first!", type="error")
        return()
      }
      
      seg_path <- resolve_wd_path(input$segmentation_file)
      x4$seg_pdata <- try(extract_template_pdata(seg_path), silent = TRUE)
      if (inherits(x4$seg_pdata, "try-error") || is.null(x4$seg_pdata) || nrow(x4$seg_pdata) == 0L) {
        showNotification("Could not extract coordinate template from selected file.", type = "error", duration = 8)
        message("MaskedAnalysis template extraction failed: ", as.character(x4$seg_pdata))
        x4$seg_pdata <- NULL
        return()
      }
      print(sprintf(
        "[MaskedAnalysis] Loaded coordinate template: pixels=%d runs=%d",
        nrow(x4$seg_pdata), length(unique(as.character(Cardinal::run(x4$seg_pdata))))
      ))
      x4$seg_filename = input$segmentation_file
      
      
      
      if (is.null(x1$raw_list) | is.null(x4$seg_pdata)) {
        print(
          "Please choose raw data in Data Setup Tab and input segmentation tempate first for peakpicking"
        )
        return(NULL)
        
      } else if (length(unique(as.character(Cardinal::run(x4$seg_pdata)))) == 1 &&  length(names(x1$raw_list)) == 1) {
        message("only one run in segmented file and one in raw file list, will proceed")
        message("Ensure correct datasets being used.")
        
      } else if (sum(unique(as.character(Cardinal::run(x4$seg_pdata))) %in% names(x1$raw_list)) < 1) {
        
        print("no raw file compatible with dataset")
        print("run names in template coordinate file")
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
      if (is.null(x4$seg_pdata)) {
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
      #browser()
      withProgress(message = "setting up coordinates and extracting pixels from raw data",
                   detail = "using existing coordinates",
                   value = 0.2,
                   {
                     req(x4$seg_pdata)
                     
                     message("setting up coordinates and extracting pixels from raw data")
                     
                     
                    
                     x1 = setup_values()[["x1"]] # bring in raw_files list
                     
                     setCardinalBPPARAM(par_mode())
                     setCardinalNChunks(setup_values()[["chunks"]])
                     
                     #browser()
                     
                     #extract only runs available in raw list
                     
                     names(x1$raw_list) <- gsub("\\.imzML$", "", names(x1$raw_list), ignore.case = TRUE)
                     
                     #check to make sure names in x4$seg_file are in names of x1$raw_list
                     if (sum(unique(as.character(Cardinal::run(x4$seg_pdata))) %in% names(x1$raw_list)) < 1) {
                       
                         #if there is only one run in coordinate template and one in raw_list, continue with warning
                         if (length(unique(as.character(Cardinal::run(x4$seg_pdata)))) == 1 & length(names(x1$raw_list)) == 1) {
                           message("only one run in segmented file and one in raw file list, will proceed")
                           showNotification(
                             "only one run in segmented file and one in raw file list, will proceed. Ensure correct datasets being used.",
                             type = "warning"
                           )
                           
                           Cardinal::run(x4$seg_pdata) <- factor(
                             rep(names(x1$raw_list), nrow(x4$seg_pdata)),
                             levels = names(x1$raw_list)
                           )
                           
                           
  
                         } else {
                         
                         
                         message("no raw file compatible with dataset, names shown below:")
                         showNotification("no raw file compatible with dataset", type = "error")
                         print(unique(as.character(Cardinal::run(x4$seg_pdata))))
                         print(names(x1$raw_list))
                         return()
                       }
                     } 
                     
                     seg_file_trimmed <- x4$seg_pdata[as.character(Cardinal::run(x4$seg_pdata)) %in% names(x1$raw_list), ]
                     
                     
                     x4$seg_file_trimmed<-seg_file_trimmed
                     
                     #reorder raw list to match file names in segmented file
                     seg_runs <- unique(as.character(Cardinal::run(seg_file_trimmed)))
                     raw_list_ord <-
                       x1$raw_list[seg_runs]

                     # Keep only runs that were actually found in the raw list.
                     is_valid_raw <- sapply(raw_list_ord, function(z) !is.null(z))
                     valid_raw_idx <- which(!is.na(names(raw_list_ord)) & is_valid_raw)
                     if (length(valid_raw_idx) < 1) {
                       showNotification("No valid raw runs matched the coordinate template after filtering.", type = "error")
                       return()
                     }
                     raw_list_ord <- raw_list_ord[valid_raw_idx]
                     seg_runs <- names(raw_list_ord)
                     
                     # seg_file_ord<-lapply(1:length(runNames(x4$seg_file)),
                     #                      function(x) x4$seg_file[,run(x4$seg_file)%in%names(x1$raw_list)[x]])
                     #
                     # seg_file_ord<-combine(seg_file_ord)
                     #
                     
                     seg_file_ord <- seg_file_trimmed
                     
                     #get coordinates from coordinate dataset
                     coord_list_segmented <-
                       lapply(seg_runs, function(x)
                         coord(seg_file_ord)[as.character(Cardinal::run(seg_file_ord)) %in% x,])
                     #sum(unlist(lapply(1:16, function(x) dim(coord_list_segmented[[x]])[1])))  # check how many pixels selected
                     names(coord_list_segmented) <-
                       seg_runs
                     
                     #if only using a subset of the masked coordinates, reduce to those only in the coordinate set
                     coord_list_segmented <- coord_list_segmented[names(coord_list_segmented) %in% names(raw_list_ord)]
                     if (length(coord_list_segmented) < 1) {
                       showNotification("Coordinate template has no runs overlapping with the selected raw data.", type = "error")
                       return()
                     }
                     
                     
                     
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
                           x1$raw_list[unique(as.character(Cardinal::run(seg_file_trimmed)))]
                         raw_list_ord <-
                           raw_list_ord[sapply(raw_list_ord, function(i)
                             ! is.null(i))]
                         seg_file_ord <-
                           seg_file_ord[as.character(Cardinal::run(seg_file_ord)) %in% names(raw_list_ord), ]
                         
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
                     
                     coord_map_ptm <- proc.time()

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
                         raw_img <- raw_list[[index]]
                         raw_coord <- as.data.frame(coord(raw_img))
                         target_coord <- as.data.frame(coord_set[[index]])

                         if (nrow(raw_coord) == nrow(target_coord)) {
                           same_xy <- suppressWarnings(
                             all(raw_coord$x == target_coord$x & raw_coord$y == target_coord$y, na.rm = FALSE)
                           )
                           if (isTRUE(same_xy)) {
                             return(raw_img)
                           }
                         }

                         idx <- prodlim::row.match(raw_coord, target_coord)
                         raw_img[!is.na(idx)]
                       }
                     
                     
                     #most likely point of failure here....
                     get_pixel_count <- function(obj) {
                       if (is.null(obj)) return(NA_integer_)

                       n_try <- suppressWarnings(try(as.integer(ncol(obj)), silent = TRUE))
                       if (!inherits(n_try, "try-error") && length(n_try) == 1 && is.finite(n_try)) {
                         return(n_try)
                       }

                       n_try <- suppressWarnings(try(as.integer(length(obj)), silent = TRUE))
                       if (!inherits(n_try, "try-error") && length(n_try) == 1 && is.finite(n_try)) {
                         return(n_try)
                       }

                       c_try <- suppressWarnings(try(as.data.frame(coord(obj)), silent = TRUE))
                       if (!inherits(c_try, "try-error")) {
                         return(as.integer(nrow(c_try)))
                       }

                       NA_integer_
                     }

                     raw_npix <- vapply(raw_list_ord, get_pixel_count, integer(1))
                     coord_npix <- vapply(coord_list_reduced, function(obj) {
                       suppressWarnings(as.integer(nrow(obj)))
                     }, integer(1))

                     message(sprintf(
                       "[MaskedAnalysis] Raw object classes: %s | raw_npix=%s | coord_npix=%s",
                       paste(vapply(raw_list_ord, function(z) class(z)[1], character(1)), collapse = ","),
                       paste(raw_npix, collapse = ","),
                       paste(coord_npix, collapse = ",")
                     ))

                     full_cover <- isTRUE(select_ratio >= 100) &&
                       length(raw_list_ord) == length(coord_list_reduced) &&
                       length(raw_npix) > 0 &&
                       all(is.finite(raw_npix)) &&
                       all(is.finite(coord_npix)) &&
                       all(raw_npix == coord_npix)

                     if (full_cover) {
                       message("[MaskedAnalysis] Coordinate template matches full raw pixel coverage; skipping row-match remap.")
                       tmp_seg_coord_list <- raw_list_ord
                     } else if (length(raw_list_ord) == 1) {
                       tmp_seg_coord_list <- list(select_pix(1, raw_list_ord, coord_list_reduced))
                     } else if (length(raw_list_ord) > 1) {
                       tmp_seg_coord_list <- try(bplapply(seq_along(raw_list_ord), function(x)
                         select_pix(x, raw_list_ord, coord_list_reduced)))
                       if (inherits(tmp_seg_coord_list, "try-error")) {
                         message("setting coordinate list failed, check names")
                         showNotification("setting coordinate list failed, check names")
                         return()
                       }
                     } else {
                       tmp_seg_coord_list <- NULL
                     }

                     if (!is.null(tmp_seg_coord_list) && !inherits(tmp_seg_coord_list, "try-error")) {
                       tmp_names <- unlist(lapply(tmp_seg_coord_list, runNames))
                       if (length(tmp_names) == length(tmp_seg_coord_list)) {
                         names(tmp_seg_coord_list) <- tmp_names
                       }
                     }

                     if (!is.null(tmp_seg_coord_list) && !inherits(tmp_seg_coord_list, "try-error") && length(tmp_seg_coord_list) >= 1) {
                       test_raw_reduced <- convertMSImagingExperiment2Arrays(tmp_seg_coord_list[[1]])
                       if (length(tmp_seg_coord_list) > 1) {
                         i <- 1
                         while (i < length(tmp_seg_coord_list)) {
                           k <- i + 1
                           tmp <- convertMSImagingExperiment2Arrays(tmp_seg_coord_list[[k]])
                           test_raw_reduced <- try(combine(test_raw_reduced, tmp))

                           if (inherits(test_raw_reduced, "try-error")) {
                             message("combining imagesets failed, check data files")
                             showNotification("combining imagesets failed, check data files", type = "error")
                             return()
                           }
                           i <- i + 1
                         }
                       }
                     } else {
                       test_raw_reduced <- NULL
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

                     coord_map_elapsed <- proc.time() - coord_map_ptm
                     message(sprintf("[MaskedAnalysis] Coordinate mapping stage elapsed: %.2fs", unname(coord_map_elapsed["elapsed"])))
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
                       #browser()
                       test_raw_reduced<-convertMSImagingArrays2Experiment(test_raw_reduced)
                       #pData(test_raw_reduced)<-pData(seg_file_ord)
                       
                      }
                     
                       
                       
                       
                     
                     
                     incProgress(amount = 0.2, message = "Coordinates mapped, now binning / pick picking")
                     #create final peak picked file
                     #browser()
                     if (input$targeted_pp == "untargeted") {
                       #setCardinalBPPARAM(SerialParam())
                       
                       
                       print("performing untargeted analysis")
                       
                       #browser()
                       
                       x4$seg_pp_file <-
                         
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
                       print(paste0("Number of peaks picked: ", nrow(x4$seg_pp_file)))
                       
                       
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
                       
                       
                       print(paste0("Number of peaks picked: ", nrow(x4$seg_pp_file)))
                       
                     } else if (input$targeted_pp == "targeted") {
                       print("performing targeted analysis")
                       print(paste0("Tolerance set to ", input$tol2))
                       file <- input$masses_file
                       ext <- tools::file_ext(file$datapath)
                       
                       req(file)
                       shiny::validate(shiny::need(ext == "txt", "Please upload a tab delimited .txt file"))
                       
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
                       gc()
                       targeted_ptm <- proc.time()
                       tmp.img<-try(test_raw_reduced |> normalize() |> peakProcess( 
                                                ref=neg_ref_mz,
                                                #SN=input$SNR,
                                                type="area",
                                                tolerance=input$tol2, units="ppm") %>% process() %>% summarizeFeatures()
                       )
                       targeted_elapsed <- proc.time() - targeted_ptm
                       message(sprintf("[MaskedAnalysis] Targeted binning stage elapsed: %.2fs", unname(targeted_elapsed["elapsed"])))
                       
                       
                       
                       
                       
                         x4$seg_pp_file <-tmp.img
                       
                       #report number of masses binned
                       print(paste0("Number of masses binned: ", nrow(x4$seg_pp_file)))
                       
                     }
                     #add mass annotations?
                     if (input$targeted_pp == "targeted") {
                       featureData(x4$seg_pp_file)$ID <- neg_masses$ID
                     }
                     
                    ##removed for now...
                     # #if using downsampling, re-bin
                     # if (input$pix_for_peak_picking2 < 100) {
                     #   #combine reduced images.
                     #   test_raw_reduced <-
                     #     combine_card(bplapply(1:length(raw_list_ord), function(x)
                     #       select_pix(x, raw_list_ord, coord_list_segmented)))
                     #   #raw_reduced<-combine(bplapply(1:length(raw_list), function(x) select_pix(x,raw_list, coord_list_reduced)))
                     #   
                     #   
                     #   test_mz_reduced <- x4$seg_pp_file
                     #   
                     #   x4$seg_pp_file <-
                     #     ref_to_peaks(
                     #       test_raw_reduced,
                     #       mz_ref = mz(test_mz_reduced),
                     #       tol = input$tol2
                     #     )
                     # }
                     
                     
                     #restore pdata
                     #function to match original pdata to current data coordinates
                     #needed in the event of changes to coordinates externally
                     #pdat is original pData, msddf is the newly created data 
                     pdat_match<-function(pdat, msddf){
                       #msddf<-x4$seg_pp_file
                       #pdat<-pData(seg_file_ord)
                       
                       
                       df1 <- as.data.frame(coord(pdat))
                       df2 <- as.data.frame(coord(msddf))
                       
                       df1$run <- as.character(run(pdat))
                       df2$run <- as.character(run(msddf))
                       
                       key1 <- paste(df1$run, df1$x, df1$y, sep = "\r")
                       key2 <- paste(df2$run, df2$x, df2$y, sep = "\r")
                       ord_vec <- match(key2, key1)
                       if (anyNA(ord_vec)) {
                         message(sprintf(
                           "[MaskedAnalysis] pData remap: %d/%d coordinates did not match template; preserving original row order for unmatched entries.",
                           sum(is.na(ord_vec)), length(ord_vec)
                         ))
                         fallback <- seq_len(nrow(df2))
                         fallback[!is.na(ord_vec)] <- ord_vec[!is.na(ord_vec)]
                         ord_vec <- fallback
                       }
                       sorted_pdat1 <- pdat[ord_vec, ]
                       
                       return(sorted_pdat1)
                       
                     }
                     
                     
                     pdat_ptm <- proc.time()
                     template_pdat <- NULL
                     if (is(seg_file_ord, "PositionDataFrame")) {
                       template_pdat <- seg_file_ord
                     } else {
                       template_pdat <- try(pData(seg_file_ord), silent = TRUE)
                     }
                     if (inherits(template_pdat, "try-error") || is.null(template_pdat)) {
                       showNotification("Could not recover template pData for remapping after masked processing.", type = "error")
                       message("[MaskedAnalysis] template_pdata recovery failed")
                       return()
                     }

                     ordered_pdat <- pdat_match(pdat = template_pdat, msddf = x4$seg_pp_file)
                     pdat_elapsed <- proc.time() - pdat_ptm
                     message(sprintf("[MaskedAnalysis] pData remap stage elapsed: %.2fs", unname(pdat_elapsed["elapsed"])))
                     
                     pData(x4$seg_pp_file) <- ordered_pdat
                     
                     plot_card_server("card_plot_masked", overview_peaks_sel = x4$seg_pp_file)
                     
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
      if (is.null(x4$seg_pp_file)) {
        showNotification("No segmented/peak-picked dataset to save. Run 'Peak pick or bin segmented data' first.", type = "error", duration = 8)
        return()
      }
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
