### R/CorrelationServer.R
CorrelationServer <- function(id, proc_values, setup_values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns # For dynamic variable namespace
    
    # Reactive expression to access x5 from proc_values
    x5 <- reactive({
      req(proc_values()[["x5"]])
      proc_values()[["x5"]]
    })
    
    # Functions to check and get files from the working directory
    has_new_files <- function() {
      unique(list.files(req(setup_values())[["wd"]], recursive = TRUE))
    }
    
    get_files <- function() {
      list.files(req(setup_values())[["wd"]], recursive = TRUE)
    }
    
    # Reactive polling to detect new files every 10 seconds (10000 ms)
    my_files <- reactivePoll(10000, session, checkFunc = has_new_files, valueFunc = get_files)
    valid_file_choices <- reactive({
      files <- my_files()
      files <- files[!is.na(files) & nzchar(files)]
      files[grepl("\\.(imzML|rds)$", basename(files), ignore.case = TRUE)]
    })
    
    # Update the 'correlation_file' selectInput when new files are detected
    observeEvent(my_files(), ignoreInit = TRUE, ignoreNULL = TRUE, {
      updateSelectInput(session, ns('correlation_file'), choices = valid_file_choices())
    })
    
    # ReactiveValues to store data_file and hmap_choices
    x6 <- reactiveValues(
      data_file = NULL,
      hmap_choices = NULL,
      coloc_raw = NULL,
      coloc_sorted = NULL,
      coloc_data_file = NULL,
      coloc_seed_mz = NULL,
      coloc_sort_metric = NULL,
      coloc_full_table = FALSE
    )

    sort_coloc_table <- function(df, metric = NULL, direction = "desc") {
      if (is.null(df)) return(NULL)
      out <- as.data.frame(df, stringsAsFactors = FALSE)
      if (!is.null(metric) && nzchar(metric) && metric %in% colnames(out) && nrow(out) > 1) {
        ord_key <- tryCatch(xtfrm(out[[metric]]), error = function(e) seq_len(nrow(out)))
        ord <- order(ord_key, decreasing = identical(direction, "desc"), na.last = TRUE)
        out <- out[ord, , drop = FALSE]
      }
      rownames(out) <- NULL
      out
    }

    format_sig4_label <- function(x) {
      out <- suppressWarnings(as.numeric(x))
      lbl <- as.character(x)
      if (length(lbl) != length(out)) {
        lbl <- vapply(out, function(v) {
          if (is.finite(v)) format(v, digits = 15, trim = TRUE) else NA_character_
        }, character(1))
      }
      ok <- is.finite(out)
      if (any(ok)) {
        lbl[ok] <- vapply(out[ok], function(v) format(signif(v, 4), digits = 4, trim = TRUE), character(1))
      }
      lbl
    }

    format_mz_fixed4_label <- function(x) {
      out <- suppressWarnings(as.numeric(x))
      lbl <- as.character(x)
      ok <- is.finite(out)
      if (any(ok)) {
        lbl[ok] <- sprintf("%.4f", out[ok])
      }
      lbl
    }

    make_mz_select_choices <- function(mz_values) {
      mz_num <- suppressWarnings(as.numeric(mz_values))
      mz_num <- mz_num[is.finite(mz_num)]
      if (length(mz_num) == 0) return(character(0))
      vals <- vapply(mz_num, function(v) format(v, digits = 15, trim = TRUE), character(1))
      lbls <- format_mz_fixed4_label(mz_num)
      stats::setNames(vals, lbls)
    }

    refresh_coloc_sorted <- function() {
      x6$coloc_sorted <- sort_coloc_table(
        x6$coloc_raw,
        metric = x6$coloc_sort_metric,
        direction = if (!is.null(input$corr_table_order)) input$corr_table_order else "desc"
      )
      invisible(NULL)
    }

    get_selected_coloc_row <- function() {
      idx <- suppressWarnings(as.integer(input$corr_table_rows_selected))
      idx <- idx[is.finite(idx)]
      if (length(idx) == 0 || is.null(x6$coloc_sorted) || nrow(x6$coloc_sorted) == 0) return(NULL)
      idx <- idx[[1]]
      if (idx < 1 || idx > nrow(x6$coloc_sorted)) return(NULL)
      x6$coloc_sorted[idx, , drop = FALSE]
    }

    run_colocalization <- function(seed_mz_override = NULL) {
      if (identical(input$corr_source, "from_stats")) {
        if (is.null(x6$data_file)) {
          showNotification("No stats-backed dataset available. Run/read stats data first.", type = "warning", duration = 7)
          return(invisible(NULL))
        }
      } else if (identical(input$corr_source, "from_file")) {
        if (is.null(x6$data_file)) {
          showNotification("No file loaded. Click 'Read Selected File' first.", type = "warning", duration = 7)
          return(invisible(NULL))
        }
      } else {
        showNotification("Invalid data source selected.", type = "error")
        return(invisible(NULL))
      }

      selected_mz <- if (!is.null(seed_mz_override)) as.character(seed_mz_override) else input$corr_mz
      if (is.null(selected_mz) || !nzchar(selected_mz)) {
        showNotification("Choose a valid m/z value before generating colocalization.", type = "warning", duration = 7)
        return(invisible(NULL))
      }

      mz_num <- suppressWarnings(as.numeric(selected_mz))
      if (!is.finite(mz_num)) {
        showNotification("Selected m/z value is not numeric.", type = "error")
        return(invisible(NULL))
      }

      plot_n_req <- suppressWarnings(as.integer(input$corr_plot_n))
      if (!is.finite(plot_n_req) || plot_n_req < 1) plot_n_req <- 3L
      # Without a separate "number of features" control, compute a modest table by
      # default and use the full-table toggle when the user wants all rows.
      n_req <- if (isTRUE(input$corr_full_table)) Inf else max(25L, plot_n_req)

      coloc <- tryCatch({
        colocalized(
          x6$data_file,
          mz = mz_num,
          n = n_req,
          sort.by = input$corr_sort
        )
      }, error = function(e) {
        showNotification(paste("Error in correlation computation:", e$message), type = "error", duration = 8)
        NULL
      })
      if (is.null(coloc)) return(invisible(NULL))

      coloc_df <- tryCatch(as.data.frame(coloc), error = function(e) NULL)
      if (is.null(coloc_df) || nrow(coloc_df) == 0) {
        showNotification("Correlation returned no rows.", type = "warning", duration = 6)
        return(invisible(NULL))
      }
      if (!"mz" %in% colnames(coloc_df)) {
        showNotification("Correlation result does not contain 'mz' values.", type = "error", duration = 8)
        return(invisible(NULL))
      }

      x6$coloc_raw <- coloc_df
      x6$coloc_data_file <- x6$data_file
      x6$coloc_seed_mz <- mz_num
      x6$coloc_sort_metric <- input$corr_sort
      x6$coloc_full_table <- isTRUE(input$corr_full_table)
      refresh_coloc_sorted()

      showNotification(
        sprintf("Colocalization computed for m/z %.6f (%d rows).", mz_num, nrow(coloc_df)),
        type = "message",
        duration = 5
      )
      invisible(NULL)
    }
    
    ### Render UI for Source Selection ###
    output$source_ui <- renderUI({
      if (input$corr_source == "from_file") {
        tagList(
          selectInput(
            ns("correlation_file"),
            "Select File for Correlation Analysis",
            choices = valid_file_choices(),
            selected = NULL
          ),
          actionButton(ns("action_seg"), label = HTML("Read Selected File"))
        )
      } else {
        # If 'from_stats' is selected, no additional UI elements are needed here
        NULL
      }
    })
    
    ### Observe Event for Reading File ###
    observeEvent(input$action_seg, {
      if (is.null(input$correlation_file) || !nzchar(input$correlation_file)) {
        showNotification("Choose a .imzML or .rds file before clicking 'Read Selected File'.", type = "warning", duration = 7)
        return(NULL)
      }
      
      # Construct full file path
      file_path <- file.path(setup_values()[["wd"]], input$correlation_file)
      
      # Read the selected file based on its extension
      data_file <- NULL
      if (grepl("\\.rds$", basename(input$correlation_file), ignore.case = TRUE)) {
        data_file <- tryCatch({
          readRDS(file_path)
        }, error = function(e) {
          showNotification(paste("Error reading .rds file:", e$message), type = "error")
          NULL
        })
      } else if (grepl("\\.imzML$", basename(input$correlation_file), ignore.case = TRUE)) {
        data_file <- tryCatch({
          readMSIData(file_path)
        }, error = function(e) {
          showNotification(paste("Error reading .imzML file:", e$message), type = "error")
          NULL
        })
      } else {
        showNotification("Unsupported file type selected.", type = "error")
      }
      
      if (!is.null(data_file)) {
        x6$data_file <- data_file
        showNotification("Data file loaded successfully.", type = "message")
      }
    })
    
    ### Observe Event for Source Selection ###
    observeEvent(input$corr_source, {
      if (input$corr_source == "from_stats") {
        if (is.null(x5()$data_file_selected)) {
          x6$data_file <- NULL
          x6$hmap_choices <- NULL
          showNotification("No stats dataset available. Run a stats test first, then switch source to 'From stats'.", type = "warning", duration = 8)
          return(NULL)
        }
        if (is.null(x5()$stats_results)) {
          x6$data_file <- NULL
          x6$hmap_choices <- NULL
          showNotification("No stats results found. Run a stats test first.", type = "warning", duration = 8)
          return(NULL)
        }
        x6$data_file <- x5()$data_file_selected
        
        # Extract hmap_choices from stats_results, excluding 'mz' and 'feature'
        x6$hmap_choices <- setdiff(colnames(x5()$stats_results), c("mz", "feature"))
      } else {
        # If 'from_file' is selected, clear hmap_choices
        x6$hmap_choices <- NULL
      }
    }, ignoreNULL = FALSE) # Trigger on initial load as well
    
    ### Render UI for Correlation Variables ###
    output$corr_vars <- renderUI({
      if (input$corr_source == "from_stats") {
        req(x6$hmap_choices)
        tagList(
          selectInput(
            ns("corr_sig_select"),
            "Choose Filtering Column for Correlation Choices",
            choices = x6$hmap_choices,
            selected = "fdr"
          ),
          radioButtons(
            ns("corr_sig_direction"),
            label = "Filtering Statistic Direction:",
            choices = list("Ascending (≤)" = "ascending", "Descending (≥)" = "descending"),
            selected = "descending"
          ),
          numericInput(ns("corr_sig"), "Significance Cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
          checkboxInput(ns("key_on"), "Display Color Key", value = TRUE)
        )
      } else {
        # If 'from_file', no additional correlation variables are needed here
        NULL
      }
    })
    
    ### Render UI for m/z Selection ###
    output$corr_mz <- renderUI({
      if (input$corr_source == "from_stats") {
        req(x5()$stats_results)
        req(input$corr_sig_select)
        req(input$corr_sig_direction)
        req(input$corr_sig)
        
        # Filter m/z based on statistical test results
        if (input$corr_sig_direction == "descending") {
          sig_mz <- x5()$stats_results %>%
            dplyr::filter(.data[[input$corr_sig_select]] >= input$corr_sig) %>%
            dplyr::pull(mz)
        } else {
          sig_mz <- x5()$stats_results %>%
            dplyr::filter(.data[[input$corr_sig_select]] <= input$corr_sig) %>%
            dplyr::pull(mz)
        }
        
        # Remove duplicates and sort
        sig_mz <- unique(sort(sig_mz))
        
        mz_choices <- make_mz_select_choices(sig_mz)

        selectInput(
          ns("corr_mz"),
          "Choose m/z Value for Correlation Analysis",
          choices = mz_choices,
          selected = if (!is.null(input$corr_mz) && nzchar(input$corr_mz)) input$corr_mz else NULL
        )
      } else if (input$corr_source == "from_file") {
        req(x6$data_file)
        
          mz_values <- mz(x6$data_file)
          mz_choices <- make_mz_select_choices(mz_values)
        
        selectInput(
          ns("corr_mz"),
          "Choose m/z Value for Correlation Analysis",
          choices = mz_choices,
          selected = if (!is.null(input$corr_mz) && nzchar(input$corr_mz)) input$corr_mz else NULL
        )
      } else {
        NULL
      }
    })
    
    ### Correlation Execution / Result Controls ###
    observeEvent(input$action_corr, {
      run_colocalization()
    })

    observeEvent(input$corr_table_order, {
      if (!is.null(x6$coloc_raw)) refresh_coloc_sorted()
    }, ignoreInit = TRUE)

    observeEvent(input$corr_use_selected, {
      row <- get_selected_coloc_row()
      if (is.null(row) || !"mz" %in% colnames(row)) {
        showNotification("Select a row in the correlation table first.", type = "warning", duration = 6)
        return(NULL)
      }
      mz_val <- as.character(row$mz[[1]])
      updateSelectInput(session, "corr_mz", selected = mz_val)
      showNotification(sprintf("Set correlation seed m/z to %s (results not recomputed yet).", mz_val), type = "message", duration = 5)
    })

    observeEvent(input$corr_rerun_selected, {
      row <- get_selected_coloc_row()
      if (is.null(row) || !"mz" %in% colnames(row)) {
        showNotification("Select a row in the correlation table first.", type = "warning", duration = 6)
        return(NULL)
      }
      mz_val <- as.character(row$mz[[1]])
      updateSelectInput(session, "corr_mz", selected = mz_val)
      run_colocalization(seed_mz_override = mz_val)
    })

    output$corr_table <- DT::renderDT({
      req(x6$coloc_sorted)
      df <- x6$coloc_sorted
      dt <- DT::datatable(
        df,
        rownames = FALSE,
        selection = "multiple",
        options = list(
          pageLength = 15,
          scrollX = TRUE
        )
      )
      num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
      signif_cols <- setdiff(num_cols, c("mz", "i"))
      if (length(signif_cols) > 0) {
        dt <- DT::formatSignif(dt, columns = signif_cols, digits = 4)
      }
      if ("mz" %in% colnames(df)) {
        dt <- DT::formatRound(dt, columns = "mz", digits = 4)
      }
      if ("i" %in% colnames(df)) {
        dt <- DT::formatRound(dt, columns = "i", digits = 0)
      }
      dt
    }, server = TRUE)

    ### Render Correlation Plot (from cached results only) ###
    output$plot15 <- renderImage({
      req(x6$coloc_sorted)
      req(x6$coloc_data_file)
      req(input$corr_plot_n)

      coloc_df <- x6$coloc_sorted
      if (!"mz" %in% colnames(coloc_df) || nrow(coloc_df) == 0) {
        return(NULL)
      }

      sel_idx <- suppressWarnings(as.integer(input$corr_table_rows_selected))
      sel_idx <- sel_idx[is.finite(sel_idx) & sel_idx >= 1 & sel_idx <= nrow(coloc_df)]

      if (length(sel_idx) > 0) {
        sel_idx <- unique(sel_idx)
        plot_mz <- suppressWarnings(as.numeric(coloc_df$mz[sel_idx]))
      } else {
        plot_n <- suppressWarnings(as.integer(input$corr_plot_n))
        if (!is.finite(plot_n) || plot_n < 1) plot_n <- 1L
        plot_n <- min(plot_n, nrow(coloc_df))
        plot_mz <- suppressWarnings(as.numeric(coloc_df$mz[seq_len(plot_n)]))
      }
      plot_mz <- plot_mz[is.finite(plot_mz)]
      if (isTRUE(input$corr_plot_include_source) && is.finite(x6$coloc_seed_mz)) {
        plot_mz <- unique(c(as.numeric(x6$coloc_seed_mz), plot_mz))
      }
      if (length(plot_mz) == 0) {
        showNotification("No valid m/z values available in cached correlation results for plotting.", type = "warning", duration = 6)
        return(NULL)
      }

      n_runs <- suppressWarnings(length(runNames(x6$coloc_data_file)))
      if (!is.finite(n_runs) || n_runs < 1) n_runs <- 1L

      # Keep run-by-feature layout for multi-run data; use a denser grid for single-run
      # datasets so labels stay readable when plotting many correlated ions.
      if (n_runs == 1L) {
        ncol_layout <- max(1L, min(4L, ceiling(sqrt(length(plot_mz)))))
        nrow_layout <- max(1L, ceiling(length(plot_mz) / ncol_layout))
      } else {
        nrow_layout <- n_runs
        ncol_layout <- length(plot_mz)
      }

      panel_w <- if (n_runs == 1L) 360L else 320L
      panel_h <- if (n_runs == 1L) 300L else 280L
      img_width <- max(900L, min(2800L, as.integer(panel_w * ncol_layout + if (isTRUE(input$key_on)) 180L else 80L)))
      img_height <- max(650L, min(2400L, as.integer(panel_h * nrow_layout + 110L)))

      outfile <- tempfile(fileext = ".png")
      grDevices::png(outfile, width = img_width, height = img_height, res = 120)
      on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

      img_args <- list(
        x = x6$coloc_data_file,
        mz = plot_mz,
        layout = c(nrow_layout, ncol_layout),
        scale = TRUE,
        col = function(n) rev(cpal("Spectral")(n)),
        colorkey = isTRUE(input$key_on)
      )
      if (isTRUE(input$corr_enhance_hist)) {
        img_args$enhance <- "histogram"
      }
      if (isTRUE(input$corr_smooth_gaussian)) {
        img_args$smooth <- "gaussian"
      }

      p1 <- tryCatch({
        do.call(image, img_args)
      }, error = function(e) {
        showNotification(paste("Error in plotting:", e$message), type = "error", duration = 8)
        NULL
      })

      if (!is.null(p1)) {
        font_scale <- min(1.35, 1 + 0.06 * max(0, length(plot_mz) - 1))
        if (inherits(p1, "trellis")) {
          print(
            p1,
            par.settings = list(
              fontsize = list(
                text = 12 * font_scale,
                points = 11 * font_scale
              )
            )
          )
        } else {
          print(p1)
        }
      }

      list(
        src = outfile,
        contentType = "image/png",
        width = img_width,
        height = img_height,
        alt = "Correlation Heatmap"
      )
    }, deleteFile = TRUE)
  })
}
