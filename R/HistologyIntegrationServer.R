HistologyIntegrationServer <- function(id, setup_values, preproc_values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    xh <- reactiveValues(
      mapped_obj = NULL,
      cluster_lookup = NULL,
      mapped_column = NULL,
      mapping_source = NULL,
      opt_xy_candidates = NULL,
      restore_msi_mode = NULL,
      restore_mz_select = NULL,
      restore_mz_value = NULL,
      restore_rgb_select = NULL,
      restore_rgb_values = NULL,
      restore_rgb_auto_n = NULL,
      rgb_mz_applied = NULL,
      restore_overlay_pdata_field = NULL,
      stat_fit_grid = NULL,
      stat_fit_candidates = NULL,
      stat_fit_summary = NULL,
      polygon_cluster_result = NULL,
      polygon_label_field_before_cluster = NULL
    )
    overlay_last_recorded <- reactiveVal(NULL)

    output$msi_upload_ui <- renderUI({
      if (identical(input$msi_source, "upload")) {
        fileInput(
          ns("msi_upload"),
          "Upload MSI data (.rds)",
          accept = c(".rds")
        )
      } else {
        NULL
      }
    })

    output$overlay_info_ui <- renderUI({
      if (isTRUE(input$show_overlay_info)) {
        verbatimTextOutput(ns("overlay_info"))
      } else {
        NULL
      }
    })

    output$msi_mode_controls <- renderUI({
      req(input$msi_plot_mode)
      mode <- tolower(trimws(as.character(input$msi_plot_mode)[1]))

      if (identical(mode, "rgb")) {
        return(
          tagList(
            selectizeInput(
              ns("rgb_mz_select"),
              "RGB m/z channels (pick 2 or 3)",
              choices = character(0),
              selected = NULL,
              multiple = TRUE,
              options = list(
                maxItems = 3,
                placeholder = "Select 2 or 3 m/z values (R,G,B)"
              )
            ),
            fluidRow(
              column(
                4,
                actionButton(ns("rgb_add_from_mz"), "Add selected m/z")
              ),
              column(
                4,
                actionButton(ns("rgb_clear_mz"), "Clear RGB channels")
              ),
              column(
                4,
                actionButton(ns("rgb_apply_mz"), "Apply RGB channels")
              )
            ),
            textOutput(ns("rgb_channels_summary")),
            fluidRow(
              column(
                6,
                radioButtons(
                  ns("rgb_auto_n"),
                  "Auto-pick channels",
                  choices = c("2" = 2, "3" = 3),
                  selected = 3,
                  inline = TRUE
                )
              ),
              column(
                6,
                actionButton(ns("suggest_rgb_mz"), "Suggest RGB ions")
              )
            ),
            fluidRow(
              column(
                6,
                selectInput(
                  ns("rgb_render_mode"),
                  "RGB rendering",
                  choices = c("Dominant channel (pure RGB)" = "dominant", "Additive blend" = "additive"),
                  selected = "dominant"
                )
              ),
              column(
                6,
                sliderInput(
                  ns("rgb_bg_cutoff"),
                  "RGB background cutoff",
                  min = 0,
                  max = 1,
                  value = 0.05,
                  step = 0.01
                )
              )
            ),
            tags$small("Smart picker prioritizes high-contrast ions with low redundancy between channels.")
          )
        )
      }

      if (identical(mode, "pdata")) {
        obj <- msi_for_pdata()
        cols <- colnames(as.data.frame(Cardinal::pData(obj)))
        if (is.null(cols) || length(cols) == 0) {
          return(tags$small("No pData fields available in current MSI object."))
        }
        restore_field <- isolate(xh$restore_overlay_pdata_field)
        default_col <- if (!is.null(restore_field) && restore_field %in% cols) {
          restore_field
        } else if ("histo_cluster" %in% cols) {
          "histo_cluster"
        } else {
          cols[1]
        }
        return(
          selectInput(
            ns("overlay_pdata_field"),
            "pData field for MSI base plot",
            choices = cols,
            selected = default_col
          )
        )
      }

      tags$small("Single-ion mode uses the selected m/z value above.")
    })

    output$mapping_source_ui <- renderUI({
      req(input$mapping_source)
      if (identical(input$mapping_source, "cluster")) {
        tagList(
          textInput(ns("cluster_pdata_col"), "pData column name", value = "histo_cluster_mapped"),
          numericInput(ns("cluster_alpha_threshold"), "Alpha threshold for mapped clusters", value = 0.05, min = 0, max = 1, step = 0.01)
        )
      } else {
        tagList(
          uiOutput(ns("polygon_label_field_ui")),
          textInput(ns("polygon_pdata_col"), "pData column name", value = "polygon_region"),
          selectInput(ns("polygon_overlap_rule"), "If multiple polygons hit one pixel", choices = c("First match" = "first", "All matches (semicolon-separated)" = "all"), selected = "first")
        )
      }
    })

    observeEvent(input$mapping_source, {
      if (identical(input$mapping_source, "polygon")) {
        updateRadioButtons(session, "overlay_layer", selected = "polygon")
      } else {
        updateRadioButtons(session, "overlay_layer", selected = "cluster")
      }
    }, ignoreInit = TRUE)

    msi_data <- reactive({
      if (identical(input$msi_source, "upload")) {
        req(input$msi_upload)
        file <- input$msi_upload
        ext <- tolower(tools::file_ext(file$name))
        validate(need(ext == "rds", "Please upload an .rds MSI file."))

        out <- try(readRDS(file$datapath), silent = TRUE)
        validate(need(!inherits(out, "try-error"), "Could not read uploaded .rds MSI file."))
        return(out)
      }

      x2 <- preproc_values()[["x2"]]
      validate(need(!is.null(x2$overview_peaks_sel), "No current MSI data found. Load or peak-pick data first in Overview Analysis."))
      x2$overview_peaks_sel
    })

    msi_for_pdata <- reactive({
      if (!is.null(xh$mapped_obj)) {
        xh$mapped_obj
      } else {
        msi_data()
      }
    })

    refresh_mz_ion_inputs <- function(msi_obj) {
      mzv <- try(Cardinal::mz(msi_obj), silent = TRUE)
      if (inherits(mzv, "try-error") || length(mzv) == 0) {
        return(invisible(FALSE))
      }
      mz_num <- suppressWarnings(as.numeric(mzv))

      idx_vals <- as.character(seq_along(mzv))
      labels <- format(mzv, digits = 10, scientific = FALSE, trim = TRUE)
      choices <- idx_vals
      names(choices) <- labels

      cur <- isolate(input$mz_select)
      sel <- if (!is.null(cur) && cur %in% idx_vals) cur else idx_vals[1]
      restore_sel <- isolate(xh$restore_mz_select)
      if (!is.null(restore_sel)) {
        restore_sel <- as.character(restore_sel)[1]
        if (!is.na(restore_sel) && nzchar(restore_sel) && restore_sel %in% idx_vals) {
          sel <- restore_sel
        }
        xh$restore_mz_select <- NULL
      }
      restore_mz_val <- suppressWarnings(as.numeric(isolate(xh$restore_mz_value)))
      if (length(restore_mz_val) == 0L) restore_mz_val <- NA_real_
      if (isTRUE(is.finite(restore_mz_val)) && any(is.finite(mz_num))) {
        fin <- which(is.finite(mz_num))
        near <- fin[which.min(abs(mz_num[fin] - restore_mz_val))]
        if (isTRUE(length(near) == 1L && is.finite(near))) {
          sel <- idx_vals[near]
        }
        xh$restore_mz_value <- NULL
      }
      updateSelectizeInput(session, "mz_select", choices = choices, selected = sel, server = TRUE)

      mode_now <- tolower(trimws(as.character(input$msi_plot_mode)[1]))
      if (identical(mode_now, "rgb")) {
        cur_rgb <- isolate(input$rgb_mz_select)
        rgb_sel <- intersect(as.character(cur_rgb), idx_vals)
        restore_rgb_idx <- isolate(xh$restore_rgb_select)
        if (!is.null(restore_rgb_idx)) {
          rr <- intersect(as.character(restore_rgb_idx), idx_vals)
          if (length(rr) > 0) {
            rgb_sel <- rr[seq_len(min(3L, length(rr)))]
          }
          xh$restore_rgb_select <- NULL
        }
        restore_rgb_vals <- isolate(xh$restore_rgb_values)
        if (!is.null(restore_rgb_vals) && length(restore_rgb_vals) > 0 && any(is.finite(mz_num))) {
          fin <- which(is.finite(mz_num))
          rr <- unique(vapply(
            as.numeric(restore_rgb_vals),
            function(target) {
              fin[which.min(abs(mz_num[fin] - target))]
            },
            integer(1)
          ))
          rr <- rr[is.finite(rr) & rr >= 1L & rr <= length(idx_vals)]
          if (length(rr) > 0) {
            rgb_sel <- idx_vals[rr][seq_len(min(3L, length(rr)))]
          }
          xh$restore_rgb_values <- NULL
        }
        if (length(rgb_sel) < 2) {
          rgb_sel <- idx_vals[seq_len(min(2, length(idx_vals)))]
        }
        updateSelectizeInput(session, "rgb_mz_select", choices = choices, selected = rgb_sel, server = TRUE)
        # Keep an applied RGB selection so UI updates don't trigger repeated heavy re-renders.
        applied_now <- intersect(as.character(isolate(xh$rgb_mz_applied)), idx_vals)
        if (length(applied_now) < 2) {
          xh$rgb_mz_applied <- rgb_sel
        } else {
          xh$rgb_mz_applied <- applied_now[seq_len(min(3L, length(applied_now)))]
        }
        restore_n <- suppressWarnings(as.integer(isolate(xh$restore_rgb_auto_n)))
        if (length(restore_n) == 0L) restore_n <- NA_integer_
        if (isTRUE(is.finite(restore_n) && restore_n %in% c(2L, 3L))) {
          updateSelectInput(session, "rgb_auto_n", selected = restore_n)
          xh$restore_rgb_auto_n <- NULL
        } else {
          updateSelectInput(session, "rgb_auto_n", selected = if (length(rgb_sel) >= 3) 3 else 2)
        }
      }
      invisible(TRUE)
    }

    observeEvent(msi_data(), {
      obj <- try(msi_data(), silent = TRUE)
      if (inherits(obj, "try-error") || is.null(obj)) return()
      refresh_mz_ion_inputs(obj)
    }, ignoreInit = FALSE)

    observeEvent(input$msi_plot_mode, {
      mode <- tolower(trimws(as.character(input$msi_plot_mode)[1]))
      if (!identical(mode, "rgb")) return()
      session$onFlushed(function() {
        obj <- try(msi_data(), silent = TRUE)
        if (!inherits(obj, "try-error") && !is.null(obj)) {
          refresh_mz_ion_inputs(obj)
        }
      }, once = TRUE)
    }, ignoreInit = TRUE)

    observeEvent(input$rgb_add_from_mz, {
      obj <- try(msi_for_pdata(), silent = TRUE)
      if (inherits(obj, "try-error")) {
        showNotification("Load MSI data first.", type = "warning", duration = 5)
        return()
      }
      mzv <- suppressWarnings(as.numeric(Cardinal::mz(obj)))
      if (length(mzv) == 0) {
        showNotification("No m/z values available to add.", type = "warning", duration = 5)
        return()
      }

      idx_vals <- as.character(seq_along(mzv))
      cur <- intersect(as.character(isolate(input$rgb_mz_select)), idx_vals)
      add_one <- as.character(isolate(input$mz_select))[1]
      if (length(add_one) == 0 || is.na(add_one) || !nzchar(add_one) || !add_one %in% idx_vals) {
        showNotification("Pick an m/z in 'Select m/z' first, then click add.", type = "message", duration = 5)
        return()
      }

      next_sel <- unique(c(cur, add_one))
      if (length(next_sel) > 3) {
        next_sel <- utils::tail(next_sel, 3)
        showNotification("RGB supports up to 3 channels. Keeping the most recent 3 selections.", type = "message", duration = 6)
      }

      labels <- format(mzv, digits = 10, scientific = FALSE, trim = TRUE)
      choices <- idx_vals
      names(choices) <- labels
      updateSelectizeInput(session, "rgb_mz_select", choices = choices, selected = next_sel, server = TRUE)
    }, ignoreInit = TRUE)

    observeEvent(input$rgb_clear_mz, {
      obj <- try(msi_for_pdata(), silent = TRUE)
      if (inherits(obj, "try-error")) return()
      mzv <- suppressWarnings(as.numeric(Cardinal::mz(obj)))
      if (length(mzv) == 0) return()
      idx_vals <- as.character(seq_along(mzv))
      labels <- format(mzv, digits = 10, scientific = FALSE, trim = TRUE)
      choices <- idx_vals
      names(choices) <- labels
      updateSelectizeInput(session, "rgb_mz_select", choices = choices, selected = character(0), server = TRUE)
      showNotification("RGB channel selection cleared. Click 'Apply RGB channels' to refresh the plot.", type = "message", duration = 5)
    }, ignoreInit = TRUE)

    observeEvent(input$rgb_apply_mz, {
      obj <- try(msi_for_pdata(), silent = TRUE)
      if (inherits(obj, "try-error")) {
        showNotification("Load MSI data first.", type = "warning", duration = 5)
        return()
      }
      mzv <- suppressWarnings(as.numeric(Cardinal::mz(obj)))
      if (length(mzv) == 0) {
        showNotification("No m/z values available to apply.", type = "warning", duration = 5)
        return()
      }
      idx_vals <- as.character(seq_along(mzv))
      sel <- intersect(as.character(isolate(input$rgb_mz_select)), idx_vals)
      sel <- unique(sel)
      if (length(sel) > 3) sel <- sel[seq_len(3L)]
      xh$rgb_mz_applied <- sel
      if (length(sel) %in% c(2L, 3L)) {
        showNotification("Applied RGB channel selection.", type = "message", duration = 4)
      } else {
        showNotification("Select 2 or 3 RGB channels, then click Apply.", type = "warning", duration = 6)
      }
    }, ignoreInit = TRUE)

    output$rgb_channels_summary <- renderText({
      req(input$msi_plot_mode)
      if (!identical(tolower(trimws(as.character(input$msi_plot_mode)[1])), "rgb")) {
        return("")
      }

      obj <- try(msi_for_pdata(), silent = TRUE)
      if (inherits(obj, "try-error")) return("RGB channels: no MSI loaded")

      mzv <- suppressWarnings(as.numeric(Cardinal::mz(obj)))
      if (length(mzv) == 0) return("RGB channels: no m/z values available")

      idx_pending <- suppressWarnings(as.integer(input$rgb_mz_select))
      idx_pending <- idx_pending[is.finite(idx_pending) & idx_pending >= 1L & idx_pending <= length(mzv)]
      idx_pending <- unique(idx_pending)
      idx_applied <- suppressWarnings(as.integer(isolate(xh$rgb_mz_applied)))
      idx_applied <- idx_applied[is.finite(idx_applied) & idx_applied >= 1L & idx_applied <= length(mzv)]
      idx_applied <- unique(idx_applied)
      if (length(idx_pending) == 0 && length(idx_applied) == 0) return("RGB channels: none selected")

      fmt_ch <- function(idx) {
        if (length(idx) == 0) return("none")
        ch <- c("R", "G", "B")[seq_len(min(3, length(idx)))]
        vals <- format(mzv[idx], digits = 10, scientific = FALSE, trim = TRUE)
        paste(sprintf("%s=%s", ch, vals), collapse = ", ")
      }
      pending_txt <- fmt_ch(idx_pending)
      applied_txt <- fmt_ch(idx_applied)
      if (!identical(pending_txt, applied_txt)) {
        paste0("RGB channels (selected): ", pending_txt, " | applied: ", applied_txt)
      } else {
        paste0("RGB channels: ", applied_txt)
      }
    })

    observeEvent(input$suggest_rgb_mz, {
      req(msi_for_pdata())
      validate(need(identical(tolower(trimws(as.character(input$msi_plot_mode)[1])), "rgb"), "Switch MSI display mode to RGB before auto-picking ions."))

      obj <- msi_for_pdata()
      mz_axis <- suppressWarnings(as.numeric(Cardinal::mz(obj)))
      validate(need(length(mz_axis) >= 2, "MSI data needs at least 2 m/z values to build RGB channels."))

      n_pick <- suppressWarnings(as.integer(input$rgb_auto_n))
      if (!n_pick %in% c(2L, 3L)) n_pick <- 3L

      withProgress(message = "Selecting RGB ions", value = 0.2, {
        incProgress(0.4, detail = "Scoring channel contrast and redundancy")
        idx <- try(suggest_rgb_indices(obj, n_channels = n_pick), silent = TRUE)
        if (inherits(idx, "try-error")) idx <- integer(0)
        idx <- as.integer(idx[is.finite(idx)])
        idx <- idx[idx >= 1L & idx <= length(mz_axis)]

        if (length(idx) < 2L) {
          showNotification("Could not generate a robust RGB suggestion from this dataset. Please pick channels manually.", type = "warning", duration = 7)
          return()
        }

        idx <- idx[seq_len(min(length(idx), n_pick))]
        incProgress(0.4, detail = "Applying selection")
        idx_vals <- as.character(seq_along(mz_axis))
        labels <- format(mz_axis, digits = 10, scientific = FALSE, trim = TRUE)
        choices <- idx_vals
        names(choices) <- labels
        updateSelectizeInput(session, "rgb_mz_select", choices = choices, selected = as.character(idx), server = TRUE)

        mz_lab <- format(mz_axis[idx], digits = 10, scientific = FALSE, trim = TRUE)
        showNotification(
          sprintf("Suggested RGB ions: %s", paste(mz_lab, collapse = ", ")),
          type = "message",
          duration = 7
        )
      })
    }, ignoreInit = TRUE)

    histology_image <- reactive({
      req(input$histology_upload)
      file <- input$histology_upload
      out <- try(magick::image_read(file$datapath), silent = TRUE)
      validate(need(!inherits(out, "try-error"), "Could not read histology image."))
      out
    })

    cluster_overlay_image <- reactive({
      req(input$cluster_overlay_upload)
      file <- input$cluster_overlay_upload
      out <- try(magick::image_read(file$datapath), silent = TRUE)
      validate(need(!inherits(out, "try-error"), "Could not read cluster-color overlay image."))
      out
    })

    polygon_data <- reactive({
      req(input$polygon_file)
      validate(need(requireNamespace("sf", quietly = TRUE), "Package 'sf' is required for polygon mapping."))
      poly <- try(sf::st_read(input$polygon_file$datapath, quiet = TRUE), silent = TRUE)
      validate(need(!inherits(poly, "try-error"), "Could not read polygon file. Please provide a valid GeoJSON."))
      poly
    })

    output$polygon_label_field_ui <- renderUI({
      req(input$mapping_source)
      if (!identical(input$mapping_source, "polygon")) return(NULL)
      req(polygon_data())

      poly <- polygon_data()
      fields <- colnames(as.data.frame(poly))
      fields <- fields[!fields %in% attr(poly, "sf_column")]
      choices <- c("row_index", fields)
      if (!is.null(xh$polygon_cluster_result) && !is.null(xh$polygon_cluster_result$cluster_label)) {
        choices <- c("clustered_polygon_class", choices)
      }
      cur <- isolate(input$polygon_label_field)
      default <- if (!is.null(cur) && cur %in% choices) cur else if ("classification" %in% choices) "classification" else choices[1]
      selectInput(ns("polygon_label_field"), "Polygon label field", choices = choices, selected = default)
    })

    output$edge_fit_pdata_field_ui <- renderUI({
      req(input$edge_fit_signal_source)
      if (!identical(input$edge_fit_signal_source, "pdata")) return(NULL)
      obj <- try(msi_for_pdata(), silent = TRUE)
      if (inherits(obj, "try-error") || is.null(obj)) {
        return(tags$small("Load MSI data to choose a pData field for Edge Fit optimization."))
      }
      pd <- try(as.data.frame(Cardinal::pData(obj)), silent = TRUE)
      if (inherits(pd, "try-error") || is.null(pd) || ncol(pd) == 0L) {
        return(tags$small("No pData fields available in current MSI object."))
      }
      cols <- colnames(pd)
      cur <- isolate(input$edge_fit_pdata_field)
      default <- if (!is.null(cur) && cur %in% cols) {
        cur
      } else if ("polygon_cluster_class" %in% cols) {
        "polygon_cluster_class"
      } else if ("Rcol_reduced" %in% cols) {
        "Rcol_reduced"
      } else {
        cols[1]
      }
      selectInput(ns("edge_fit_pdata_field"), "pData field for Edge Fit optimization", choices = cols, selected = default)
    })

    parse_stat_fit_group_field <- function(sel = NULL) {
      s <- sel
      if (is.null(s)) s <- input$stat_fit_group_field
      s <- as.character(s)[1]
      if (is.na(s) || !nzchar(s)) {
        return(list(source = "pdata", field = "polygon_cluster_class", raw = ""))
      }
      if (startsWith(s, "pdata::")) {
        return(list(source = "pdata", field = sub("^pdata::", "", s), raw = s))
      }
      if (startsWith(s, "polygon::")) {
        return(list(source = "polygon", field = sub("^polygon::", "", s), raw = s))
      }
      list(source = "auto", field = s, raw = s)
    }

    stat_fit_group_field_uses_pdata <- function(sel = NULL) {
      spec <- parse_stat_fit_group_field(sel)
      if (identical(spec$source, "pdata")) return(TRUE)
      if (identical(spec$source, "polygon")) return(FALSE)
      obj_try <- try(msi_for_pdata(), silent = TRUE)
      if (inherits(obj_try, "try-error") || is.null(obj_try)) return(FALSE)
      pd_cols <- try(colnames(as.data.frame(Cardinal::pData(obj_try))), silent = TRUE)
      if (inherits(pd_cols, "try-error") || length(pd_cols) == 0L) return(FALSE)
      spec$field %in% pd_cols
    }

    output$stat_fit_group_field_ui <- renderUI({
      req(input$stat_fit_metric_mode)
      mode <- tolower(trimws(as.character(input$stat_fit_metric_mode)[1]))
      if (!identical(mode, "polygon_cluster_groups")) return(NULL)

      choices <- character(0)
      labels <- character(0)

      obj_try <- try(msi_for_pdata(), silent = TRUE)
      if (!inherits(obj_try, "try-error") && !is.null(obj_try)) {
        pd_cols <- try(colnames(as.data.frame(Cardinal::pData(obj_try))), silent = TRUE)
        if (!inherits(pd_cols, "try-error") && length(pd_cols) > 0) {
          pd_vals <- paste0("pdata::", pd_cols)
          choices <- c(choices, pd_vals)
          labels <- c(labels, paste0(pd_cols, " (pData)"))
        }
      }

      if (!is.null(input$polygon_file) && nzchar(input$polygon_file$name)) {
        poly_try <- try(polygon_data(), silent = TRUE)
        if (!inherits(poly_try, "try-error") && !is.null(poly_try)) {
          poly <- poly_try
          fields <- colnames(as.data.frame(poly))
          fields <- fields[!fields %in% attr(poly, "sf_column")]
          poly_choices <- c("row_index", fields)
          if (!is.null(xh$polygon_cluster_result) && !is.null(xh$polygon_cluster_result$cluster_label)) {
            poly_choices <- c("clustered_polygon_class", poly_choices)
          }
          poly_choices <- unique(poly_choices)
          poly_vals <- paste0("polygon::", poly_choices)
          choices <- c(choices, poly_vals)
          labels <- c(labels, paste0(poly_choices, " (polygon)"))
        }
      }

      if (length(choices) == 0L) {
        return(tags$small("Load an MSI dataset with mapped polygon groups or a polygon GeoJSON to select a Stat Fit group field."))
      }

      sel_cur <- isolate(input$stat_fit_group_field)
      defaults <- c("pdata::polygon_cluster_class", "polygon::polygon_cluster_class", "polygon::clustered_polygon_class")
      sel_default <- if (!is.null(sel_cur) && sel_cur %in% choices) {
        sel_cur
      } else {
        hit <- defaults[defaults %in% choices]
        if (length(hit) > 0L) hit[1] else choices[1]
      }

      vals <- choices
      names(vals) <- labels
      selectInput(
        ns("stat_fit_group_field"),
        "Group field for cluster-group metric",
        choices = vals,
        selected = sel_default
      )
    })

    normalize_join_key <- function(x) {
      x <- trimws(as.character(x))
      x[is.na(x)] <- ""
      make.names(x, unique = FALSE)
    }

    read_polygon_measurement_table <- function(file_path, file_name = NULL) {
      header_line <- ""
      header_line <- tryCatch(readLines(file_path, n = 1L, warn = FALSE), error = function(e) "")
      header_line <- if (length(header_line) > 0) header_line[1] else ""
      ext <- tolower(tools::file_ext(if (is.null(file_name)) file_path else file_name))
      if (identical(ext, "tsv") || grepl("\t", header_line, fixed = TRUE)) {
        out <- try(utils::read.delim(file_path, check.names = FALSE, stringsAsFactors = FALSE), silent = TRUE)
      } else {
        out <- try(utils::read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE), silent = TRUE)
        if (inherits(out, "try-error")) {
          out <- try(utils::read.delim(file_path, check.names = FALSE, stringsAsFactors = FALSE), silent = TRUE)
        }
      }
      if (inherits(out, "try-error")) return(NULL)
      as.data.frame(out, stringsAsFactors = FALSE)
    }

    add_ratio_if_present <- function(df, new_name, num_col, den_col) {
      if (!all(c(num_col, den_col) %in% names(df))) return(df)
      num <- suppressWarnings(as.numeric(df[[num_col]]))
      den <- suppressWarnings(as.numeric(df[[den_col]]))
      out <- rep(NA_real_, length(num))
      ok <- is.finite(num) & is.finite(den) & den != 0
      out[ok] <- num[ok] / den[ok]
      df[[new_name]] <- out
      df
    }

    augment_polygon_measurements <- function(df) {
      if (is.null(df) || nrow(df) == 0) return(df)
      df <- add_ratio_if_present(df, "Nucleus_Hematoxylin_Eosin_Ratio", "Nucleus: Hematoxylin OD mean", "Nucleus: Eosin OD mean")
      df <- add_ratio_if_present(df, "Cell_Hematoxylin_Eosin_Ratio", "Cell: Hematoxylin OD mean", "Cell: Eosin OD mean")
      df <- add_ratio_if_present(df, "Cytoplasm_Hematoxylin_Eosin_Ratio", "Cytoplasm: Hematoxylin OD mean", "Cytoplasm: Eosin OD mean")
      df <- add_ratio_if_present(df, "Nucleus_Aspect_Ratio", "Nucleus: Max caliper", "Nucleus: Min caliper")
      df <- add_ratio_if_present(df, "Nucleus_Area_Perimeter_Ratio", "Nucleus: Area", "Nucleus: Perimeter")
      df <- add_ratio_if_present(df, "Cytoplasm_Nucleus_Hematoxylin_Ratio", "Cytoplasm: Hematoxylin OD mean", "Nucleus: Hematoxylin OD mean")

      if (all(c("Cell: Area", "Nucleus: Area") %in% names(df))) {
        cell_area <- suppressWarnings(as.numeric(df[["Cell: Area"]]))
        nuc_area <- suppressWarnings(as.numeric(df[["Nucleus: Area"]]))
        out <- rep(NA_real_, length(cell_area))
        ok <- is.finite(cell_area) & is.finite(nuc_area) & nuc_area != 0
        out[ok] <- (cell_area[ok] - nuc_area[ok]) / nuc_area[ok]
        df[["Cytoplasm_Nucleus_Area_Ratio"]] <- out
      }
      df
    }

    extract_polygon_measurements_from_sf <- function(poly) {
      if (is.null(poly) || nrow(poly) == 0) return(NULL)
      df0 <- try(as.data.frame(sf::st_drop_geometry(poly)), silent = TRUE)
      if (inherits(df0, "try-error") || is.null(df0)) return(NULL)
      df0 <- as.data.frame(df0, stringsAsFactors = FALSE, check.names = FALSE)

      if (!"measurements" %in% names(df0)) {
        return(df0)
      }

      base_df <- df0[, setdiff(names(df0), "measurements"), drop = FALSE]
      meas_col <- df0[["measurements"]]
      n <- nrow(df0)

      parsed_rows <- vector("list", n)
      all_keys <- character(0)

      if (requireNamespace("jsonlite", quietly = TRUE)) {
        for (i in seq_len(n)) {
          mi <- meas_col[[i]]
          if (is.null(mi) || (length(mi) == 1L && (is.na(mi) || !nzchar(trimws(as.character(mi)))))) {
            parsed_rows[[i]] <- list()
            next
          }
          one <- try({
            if (is.character(mi) && length(mi) == 1L) {
              jsonlite::fromJSON(mi, simplifyVector = TRUE)
            } else {
              mi
            }
          }, silent = TRUE)
          if (inherits(one, "try-error") || is.null(one)) {
            parsed_rows[[i]] <- list()
            next
          }
          if (is.data.frame(one)) one <- as.list(one[1, , drop = FALSE])
          if (!is.list(one)) one <- list(value = one)
          one <- lapply(one, function(v) {
            if (length(v) == 0L || is.null(v)) return(NA_real_)
            v <- v[[1]]
            suppressWarnings(as.numeric(v))
          })
          nm <- names(one)
          if (is.null(nm)) {
            parsed_rows[[i]] <- list()
          } else {
            keep <- nzchar(nm)
            one <- one[keep]
            parsed_rows[[i]] <- one
            all_keys <- unique(c(all_keys, names(one)))
          }
        }
      }

      if (length(all_keys) == 0L) {
        return(base_df)
      }

      meas_mat <- matrix(NA_real_, nrow = n, ncol = length(all_keys))
      colnames(meas_mat) <- all_keys
      for (i in seq_len(n)) {
        rowi <- parsed_rows[[i]]
        if (length(rowi) == 0L) next
        nm <- intersect(names(rowi), all_keys)
        if (length(nm) == 0L) next
        meas_mat[i, nm] <- suppressWarnings(as.numeric(unlist(rowi[nm], use.names = FALSE)))
      }
      cbind(base_df, as.data.frame(meas_mat, check.names = FALSE, stringsAsFactors = FALSE))
    }

    sanitize_polygon_sf_planar <- function(poly_sf) {
      if (is.null(poly_sf) || nrow(poly_sf) == 0) return(poly_sf)
      poly <- poly_sf
      # Image-space polygons should be treated as planar geometry.
      poly <- tryCatch(sf::st_set_crs(poly, NA), error = function(e) poly)
      poly <- tryCatch(sf::st_make_valid(poly), error = function(e) poly)
      gtypes <- try(as.character(sf::st_geometry_type(poly, by_geometry = TRUE)), silent = TRUE)
      if (!inherits(gtypes, "try-error") && length(gtypes) > 0 && any(!gtypes %in% c("POLYGON", "MULTIPOLYGON"))) {
        poly <- tryCatch(sf::st_collection_extract(poly, "POLYGON", warn = FALSE), error = function(e) poly)
      }
      # Drop rows with empty geometries after repair/extract.
      empt <- try(sf::st_is_empty(poly), silent = TRUE)
      if (!inherits(empt, "try-error") && length(empt) == nrow(poly)) {
        empt[is.na(empt)] <- TRUE
        poly <- poly[!empt, , drop = FALSE]
      }
      poly
    }

    polygon_geometry_features <- function(poly) {
      n <- nrow(poly)
      if (n == 0) return(data.frame(.poly_row = integer(0)))
      poly_work <- poly
      poly_work$.orig_poly_row <- seq_len(n)
      poly_work <- sanitize_polygon_sf_planar(poly_work)
      validate(need(nrow(poly_work) > 0, "No valid polygon geometries remain after repair."))
      if (!".orig_poly_row" %in% names(poly_work)) {
        poly_work$.orig_poly_row <- seq_len(nrow(poly_work))
      }
      idx_map <- suppressWarnings(as.integer(poly_work$.orig_poly_row))
      idx_map <- idx_map[is.finite(idx_map) & idx_map >= 1L & idx_map <= n]
      validate(need(length(idx_map) > 0, "No valid polygon rows remain after repair."))

      geom <- sf::st_geometry(poly_work)
      n_work <- nrow(poly_work)

      area <- rep(NA_real_, n)
      perim <- rep(NA_real_, n)
      old_s2 <- try(sf::sf_use_s2(), silent = TRUE)
      if (!inherits(old_s2, "try-error")) {
        suppressWarnings(try(sf::sf_use_s2(FALSE), silent = TRUE))
        on.exit(suppressWarnings(try(sf::sf_use_s2(old_s2), silent = TRUE)), add = TRUE)
      }
      area_w <- rep(NA_real_, n_work)
      perim_w <- rep(NA_real_, n_work)
      for (i in seq_len(n_work)) {
        gi <- geom[i]
        area_w[i] <- tryCatch(suppressWarnings(as.numeric(sf::st_area(gi))), error = function(e) NA_real_)
        perim_w[i] <- tryCatch(suppressWarnings(as.numeric(sf::st_length(sf::st_boundary(gi)))), error = function(e) NA_real_)
      }
      area_w[!is.finite(area_w)] <- NA_real_
      perim_w[!is.finite(perim_w)] <- NA_real_
      area[idx_map] <- area_w
      perim[idx_map] <- perim_w

      bbox_vals <- lapply(seq_len(n_work), function(i) {
        bb <- sf::st_bbox(geom[i])
        c(
          xmin = as.numeric(bb[["xmin"]]),
          ymin = as.numeric(bb[["ymin"]]),
          xmax = as.numeric(bb[["xmax"]]),
          ymax = as.numeric(bb[["ymax"]])
        )
      })
      bbm <- do.call(rbind, bbox_vals)
      bb_w_w <- bbm[, "xmax"] - bbm[, "xmin"]
      bb_h_w <- bbm[, "ymax"] - bbm[, "ymin"]
      bb_area_w <- bb_w_w * bb_h_w
      bb_w <- rep(NA_real_, n)
      bb_h <- rep(NA_real_, n)
      bb_area <- rep(NA_real_, n)
      bb_w[idx_map] <- bb_w_w
      bb_h[idx_map] <- bb_h_w
      bb_area[idx_map] <- bb_area_w

      circularity <- rep(NA_real_, n)
      ok_cp <- is.finite(area) & is.finite(perim) & perim > 0
      circularity[ok_cp] <- 4 * pi * area[ok_cp] / (perim[ok_cp]^2)

      bbox_aspect <- rep(NA_real_, n)
      ok_ar <- is.finite(bb_w) & is.finite(bb_h) & pmin(bb_w, bb_h) > 0
      bbox_aspect[ok_ar] <- pmax(bb_w[ok_ar], bb_h[ok_ar]) / pmin(bb_w[ok_ar], bb_h[ok_ar])

      fill_ratio <- rep(NA_real_, n)
      ok_fr <- is.finite(area) & is.finite(bb_area) & bb_area > 0
      fill_ratio[ok_fr] <- area[ok_fr] / bb_area[ok_fr]

      n_vertices_w <- vapply(seq_len(n_work), function(i) {
        cc <- try(sf::st_coordinates(geom[i]), silent = TRUE)
        if (inherits(cc, "try-error") || is.null(dim(cc))) return(NA_real_)
        if (!("L1" %in% colnames(cc))) return(nrow(cc))
        sum(!duplicated(cc[, c("L1", "X", "Y"), drop = FALSE]))
      }, numeric(1))
      n_vertices <- rep(NA_real_, n)
      n_vertices[idx_map] <- n_vertices_w

      ctr <- try(suppressWarnings(sf::st_centroid(geom)), silent = TRUE)
      ctr_xy <- try(sf::st_coordinates(ctr), silent = TRUE)
      if (inherits(ctr_xy, "try-error") || is.null(dim(ctr_xy)) || nrow(ctr_xy) != n_work) {
        cx <- rep(NA_real_, n)
        cy <- rep(NA_real_, n)
      } else {
        cx <- rep(NA_real_, n)
        cy <- rep(NA_real_, n)
        cx[idx_map] <- suppressWarnings(as.numeric(ctr_xy[, "X"]))
        cy[idx_map] <- suppressWarnings(as.numeric(ctr_xy[, "Y"]))
      }

      data.frame(
        .poly_row = seq_len(n),
        geom_area = area,
        geom_perimeter = perim,
        geom_circularity = circularity,
        geom_bbox_width = bb_w,
        geom_bbox_height = bb_h,
        geom_bbox_aspect = bbox_aspect,
        geom_fill_ratio = fill_ratio,
        geom_n_vertices = n_vertices,
        geom_centroid_x = cx,
        geom_centroid_y = cy,
        stringsAsFactors = FALSE
      )
    }

    polygon_cluster_measurements_raw <- reactive({
      poly_df <- NULL
      if (!is.null(input$polygon_file) && !is.null(input$polygon_file$datapath)) {
        poly_obj <- try(polygon_data(), silent = TRUE)
        if (!inherits(poly_obj, "try-error") && !is.null(poly_obj)) {
          poly_df <- extract_polygon_measurements_from_sf(poly_obj)
        }
      }

      if (!is.null(input$polygon_cluster_table) && !is.null(input$polygon_cluster_table$datapath)) {
        ext_df <- read_polygon_measurement_table(input$polygon_cluster_table$datapath, input$polygon_cluster_table$name)
        validate(need(!is.null(ext_df) && nrow(ext_df) > 0, "Could not read external polygon measurements table."))
        return(ext_df)
      }

      poly_df
    })

    polygon_cluster_measurements_aug <- reactive({
      df <- polygon_cluster_measurements_raw()
      if (is.null(df)) return(NULL)
      augment_polygon_measurements(df)
    })

    guess_polygon_cluster_id_col <- function(df) {
      if (is.null(df) || ncol(df) == 0) return(NULL)
      nms <- names(df)
      low <- tolower(nms)
      preferred <- c("poly_name", "polygon_id", "polygonid", "id", "name", "object id", "objectid", "classification")
      hit <- match(preferred, low)
      hit <- hit[!is.na(hit)]
      if (length(hit) > 0) return(nms[hit[1]])
      chr_like <- which(vapply(df, function(z) is.character(z) || is.factor(z), logical(1)))
      if (length(chr_like) > 0) return(nms[chr_like[1]])
      NULL
    }

    default_polygon_cluster_feature_cols <- function(df) {
      if (is.null(df) || nrow(df) == 0) return(character(0))
      is_num <- vapply(df, function(z) is.numeric(z) || is.integer(z), logical(1))
      nms <- names(df)[is_num]
      if (length(nms) == 0) return(character(0))
      low <- tolower(nms)
      drop_pat <- "(^x$|^y$|centroid|coord|tile|row|column|^i$|index|label|class|cluster|fold|parent|child)"
      keep <- !grepl(drop_pat, low)
      out <- nms[keep]
      if (length(out) == 0) out <- nms
      out
    }

    output$polygon_cluster_id_field_ui <- renderUI({
      df <- try(polygon_cluster_measurements_aug(), silent = TRUE)
      if (inherits(df, "try-error") || is.null(df)) {
        return(tags$small("No polygon attributes/measurements detected yet. Clustering can still use geometry-only features."))
      }
      default_id <- guess_polygon_cluster_id_col(df)
      choices <- names(df)
      selectInput(
        ns("polygon_cluster_id_field"),
        "Attribute/CSV polygon ID column",
        choices = choices,
        selected = if (!is.null(default_id) && default_id %in% choices) default_id else choices[1]
      )
    })

    output$polygon_cluster_features_ui <- renderUI({
      df <- try(polygon_cluster_measurements_aug(), silent = TRUE)
      if (inherits(df, "try-error") || is.null(df)) return(NULL)
      num_cols <- names(df)[vapply(df, function(z) is.numeric(z) || is.integer(z), logical(1))]
      if (length(num_cols) == 0) {
        return(tags$small("No numeric measurement columns detected in polygon attributes / external table."))
      }
      default_cols <- default_polygon_cluster_feature_cols(df)
      default_cols <- intersect(default_cols, num_cols)
      selectizeInput(
        ns("polygon_cluster_feature_cols"),
        "Measurement features (numeric)",
        choices = num_cols,
        selected = default_cols,
        multiple = TRUE,
        options = list(placeholder = "Select measurement features (optional in auto mode)")
      )
    })

    output$polygon_cluster_join_field_ui <- renderUI({
      req(polygon_data())
      poly <- polygon_data()
      fields <- colnames(as.data.frame(poly))
      fields <- fields[!fields %in% attr(poly, "sf_column")]
      choices <- c("row_index", fields)
      cur_join <- isolate(input$polygon_cluster_join_field)
      default <- if (!is.null(cur_join) && cur_join %in% choices) {
        cur_join
      } else if ("id" %in% choices) {
        "id"
      } else if ("classification" %in% choices) {
        "classification"
      } else {
        choices[1]
      }
      selectInput(ns("polygon_cluster_join_field"), "Polygon ID field (GeoJSON)", choices = choices, selected = default)
    })

    observeEvent(input$run_polygon_clustering, {
      req(polygon_data())
      poly <- polygon_data()
      validate(need(nrow(poly) >= 2, "Need at least 2 polygons for clustering."))

      geom_tbl <- polygon_geometry_features(poly)
      n_poly <- nrow(geom_tbl)

      join_field <- as.character(input$polygon_cluster_join_field)[1]
      if (length(join_field) == 0 || is.na(join_field) || !nzchar(join_field)) join_field <- "row_index"
      if (identical(join_field, "row_index") || !join_field %in% colnames(poly)) {
        poly_key_raw <- as.character(seq_len(n_poly))
      } else {
        poly_key_raw <- as.character(poly[[join_field]])
      }
      poly_key_norm <- normalize_join_key(poly_key_raw)

      meas_df <- NULL
      meas_feat_aligned <- NULL
      matched_meas <- rep(FALSE, n_poly)
      csv_feature_cols_used <- character(0)
      csv_id_col <- NULL

      mode <- tolower(trimws(as.character(input$polygon_cluster_feature_mode)[1]))
      if (!mode %in% c("auto", "measurements", "measurements_plus_geometry", "geometry")) mode <- "auto"

      if (!identical(mode, "geometry")) {
        meas_df <- try(polygon_cluster_measurements_aug(), silent = TRUE)
        if (inherits(meas_df, "try-error")) meas_df <- NULL
      }

      if (!is.null(meas_df) && nrow(meas_df) > 0) {
        csv_id_col <- as.character(input$polygon_cluster_id_field)[1]
        if (length(csv_id_col) == 0 || is.na(csv_id_col) || !nzchar(csv_id_col) || !csv_id_col %in% names(meas_df)) {
          csv_id_col <- guess_polygon_cluster_id_col(meas_df)
        }
        join_mode <- as.character(input$polygon_cluster_join_mode)[1]
        if (!join_mode %in% c("exact", "row_order")) join_mode <- "exact"

        feat_cols <- intersect(as.character(input$polygon_cluster_feature_cols), names(meas_df))
        if (length(feat_cols) == 0) feat_cols <- default_polygon_cluster_feature_cols(meas_df)
        feat_cols <- feat_cols[vapply(meas_df[feat_cols], function(z) is.numeric(z) || is.integer(z), logical(1))]
        csv_feature_cols_used <- feat_cols

        if (length(feat_cols) > 0) {
          feat_template <- as.data.frame(matrix(NA_real_, nrow = n_poly, ncol = length(feat_cols)))
          names(feat_template) <- feat_cols

          if (identical(join_mode, "row_order")) {
            nn <- min(n_poly, nrow(meas_df))
            if (nn > 0) {
              feat_template[seq_len(nn), ] <- lapply(meas_df[seq_len(nn), feat_cols, drop = FALSE], function(z) suppressWarnings(as.numeric(z)))
              matched_meas[seq_len(nn)] <- TRUE
            }
          } else {
            validate(need(!is.null(csv_id_col) && csv_id_col %in% names(meas_df), "Select a valid CSV polygon ID column for exact matching."))
            meas_key_norm <- normalize_join_key(meas_df[[csv_id_col]])
            m_idx <- match(poly_key_norm, meas_key_norm)
            hit <- !is.na(m_idx)
            if (any(hit)) {
              feat_template[hit, ] <- lapply(meas_df[m_idx[hit], feat_cols, drop = FALSE], function(z) suppressWarnings(as.numeric(z)))
              matched_meas[hit] <- TRUE
            }
          }
          meas_feat_aligned <- feat_template
        }
      }

      geom_features <- geom_tbl[, setdiff(names(geom_tbl), ".poly_row"), drop = FALSE]
      if ("geom_centroid_x" %in% names(geom_features)) geom_features$geom_centroid_x <- NULL
      if ("geom_centroid_y" %in% names(geom_features)) geom_features$geom_centroid_y <- NULL

      feature_df <- NULL
      if (identical(mode, "geometry")) {
        feature_df <- geom_features
      } else if (identical(mode, "measurements")) {
        validate(need(!is.null(meas_feat_aligned) && ncol(meas_feat_aligned) > 0, "No numeric measurement features available. Upload a QuPath table or switch to geometry-only mode."))
        feature_df <- meas_feat_aligned
      } else if (identical(mode, "measurements_plus_geometry")) {
        validate(need(!is.null(meas_feat_aligned) && ncol(meas_feat_aligned) > 0, "No numeric measurement features available for measurements+geometry mode."))
        feature_df <- cbind(meas_feat_aligned, geom_features)
      } else {
        if (!is.null(meas_feat_aligned) && ncol(meas_feat_aligned) > 0 && any(matched_meas)) {
          feature_df <- meas_feat_aligned
        } else {
          feature_df <- geom_features
        }
      }

      validate(need(!is.null(feature_df) && ncol(feature_df) > 0, "No usable polygon features found for clustering."))
      feature_df <- as.data.frame(feature_df, stringsAsFactors = FALSE)
      feature_df[] <- lapply(feature_df, function(z) suppressWarnings(as.numeric(z)))

      keep_col <- vapply(feature_df, function(z) {
        zf <- z[is.finite(z)]
        length(zf) >= 2 && stats::sd(zf) > 0
      }, logical(1))
      feature_df <- feature_df[, keep_col, drop = FALSE]
      validate(need(ncol(feature_df) > 0, "All candidate features are constant or missing after filtering."))

      eligible <- rowSums(is.finite(as.matrix(feature_df))) > 0
      if (isTRUE(input$polygon_cluster_drop_unmatched) && any(!matched_meas) && !identical(mode, "geometry")) {
        eligible <- eligible & matched_meas
      }
      validate(need(sum(eligible) >= 2, "Need at least 2 polygons with usable features after filtering/matching."))

      X <- feature_df[eligible, , drop = FALSE]
      for (j in seq_len(ncol(X))) {
        colj <- X[[j]]
        med <- stats::median(colj[is.finite(colj)], na.rm = TRUE)
        if (!is.finite(med)) med <- 0
        colj[!is.finite(colj)] <- med
        X[[j]] <- colj
      }

      X_mat <- as.matrix(X)
      if (isTRUE(input$polygon_cluster_scale)) {
        X_mat <- scale(X_mat)
        X_mat[!is.finite(X_mat)] <- 0
      }

      pca_max <- min(ncol(X_mat), nrow(X_mat) - 1L)
      pca_req <- suppressWarnings(as.integer(input$polygon_cluster_pca_dims))
      if (!is.finite(pca_req) || pca_req < 1) pca_req <- 4L
      pca_use <- min(max(1L, pca_req), max(1L, pca_max))

      pca_obj <- NULL
      embed_mat <- X_mat
      if (ncol(X_mat) > 1 && nrow(X_mat) > 2) {
        pca_obj <- try(stats::prcomp(X_mat, center = FALSE, scale. = FALSE), silent = TRUE)
        if (!inherits(pca_obj, "try-error") && !is.null(pca_obj$x)) {
          pcs_avail <- min(ncol(pca_obj$x), pca_use)
          embed_mat <- pca_obj$x[, seq_len(pcs_avail), drop = FALSE]
        }
      }

      k_req <- suppressWarnings(as.integer(input$polygon_cluster_k))
      if (!is.finite(k_req) || k_req < 2L) k_req <- 4L
      n_unique_rows <- try(nrow(unique(as.data.frame(embed_mat))), silent = TRUE)
      if (inherits(n_unique_rows, "try-error") || !is.finite(n_unique_rows)) n_unique_rows <- nrow(embed_mat)
      k_use <- min(k_req, nrow(embed_mat), as.integer(n_unique_rows))
      validate(need(k_use >= 2L, "Need at least 2 eligible polygons for clustering."))
      if (k_use < k_req) {
        showNotification(sprintf("Reduced k from %d to %d (number of eligible polygons).", k_req, k_use), type = "message", duration = 6)
      }

      seed_val <- suppressWarnings(as.integer(input$polygon_cluster_seed))
      if (!is.finite(seed_val)) seed_val <- 123L
      nstart_val <- suppressWarnings(as.integer(input$polygon_cluster_nstart))
      if (!is.finite(nstart_val) || nstart_val < 1L) nstart_val <- 25L
      keep_prop <- suppressWarnings(as.numeric(input$polygon_cluster_keep_prop))
      if (!is.finite(keep_prop)) keep_prop <- 1
      keep_prop <- clamp(keep_prop, 0.01, 1)

      set.seed(seed_val)
      km <- try(stats::kmeans(embed_mat, centers = k_use, nstart = nstart_val), silent = TRUE)
      validate(need(!inherits(km, "try-error"), "k-means clustering failed for the selected polygon features."))

      dig <- max(2L, nchar(as.character(k_use)))
      label_levels <- paste0("cluster_", formatC(seq_len(k_use), width = dig, flag = "0"))
      labels_eligible_raw <- label_levels[as.integer(km$cluster)]
      cluster_num_eligible <- as.integer(km$cluster)

      centers_mat <- as.matrix(km$centers)
      if (is.null(dim(centers_mat))) {
        centers_mat <- matrix(centers_mat, nrow = k_use, ncol = ncol(embed_mat), byrow = TRUE)
      }
      dist_sq <- matrix(Inf, nrow = nrow(embed_mat), ncol = k_use)
      for (j in seq_len(k_use)) {
        ctr_j <- matrix(centers_mat[j, , drop = TRUE], nrow = nrow(embed_mat), ncol = ncol(embed_mat), byrow = TRUE)
        dj <- rowSums((embed_mat - ctr_j)^2)
        dj[!is.finite(dj)] <- Inf
        dist_sq[, j] <- dj
      }
      d_assigned <- dist_sq[cbind(seq_len(nrow(dist_sq)), cluster_num_eligible)]
      d_other <- vapply(seq_len(nrow(dist_sq)), function(i) {
        jj <- seq_len(k_use)
        jj <- jj[jj != cluster_num_eligible[i]]
        if (length(jj) == 0L) return(Inf)
        suppressWarnings(min(dist_sq[i, jj], na.rm = TRUE))
      }, numeric(1))
      distinctiveness <- d_other - d_assigned
      distinctiveness[!is.finite(distinctiveness)] <- -Inf

      keep_mask_eligible <- rep(TRUE, length(cluster_num_eligible))
      if (keep_prop < 0.9999) {
        keep_mask_eligible <- rep(FALSE, length(cluster_num_eligible))
        for (cl_id in sort(unique(cluster_num_eligible))) {
          ii <- which(cluster_num_eligible == cl_id)
          if (length(ii) == 0L) next
          n_keep <- max(1L, min(length(ii), ceiling(length(ii) * keep_prop)))
          ord <- order(distinctiveness[ii], decreasing = TRUE, na.last = TRUE)
          keep_mask_eligible[ii[ord[seq_len(n_keep)]]] <- TRUE
        }
      }
      labels_eligible <- labels_eligible_raw
      labels_eligible[!keep_mask_eligible] <- "unclustered_polygon"

      cluster_label_all <- rep(NA_character_, n_poly)
      cluster_num_all <- rep(NA_integer_, n_poly)
      cluster_label_all[eligible] <- labels_eligible
      cluster_num_all[eligible] <- cluster_num_eligible
      cluster_label_all[is.na(cluster_label_all)] <- "unclustered_polygon"

      pca_df <- NULL
      if (!is.null(embed_mat) && nrow(embed_mat) == sum(eligible)) {
        embed_df <- as.data.frame(embed_mat, stringsAsFactors = FALSE)
        if (ncol(embed_df) == 1) names(embed_df) <- "PC1"
        if (ncol(embed_df) >= 2) names(embed_df)[1:2] <- c("PC1", "PC2")
        pca_df <- data.frame(
          polygon_row = which(eligible),
          cluster_label = labels_eligible,
          cluster_label_raw = labels_eligible_raw,
          cluster_id = cluster_num_eligible,
          distinctiveness_score = distinctiveness,
          retained = keep_mask_eligible,
          embed_df,
          stringsAsFactors = FALSE
        )
      }

      total_counts <- tabulate(cluster_num_eligible, nbins = k_use)
      kept_counts <- tabulate(cluster_num_eligible[keep_mask_eligible], nbins = k_use)
      counts_df <- data.frame(
        cluster_id = seq_len(k_use),
        cluster_label = label_levels,
        n_polygons = as.integer(kept_counts),
        n_total = as.integer(total_counts),
        n_filtered_out = as.integer(total_counts - kept_counts),
        keep_prop_used = rep(keep_prop, k_use),
        stringsAsFactors = FALSE
      )

      preview_df <- data.frame(
        polygon_row = seq_len(n_poly),
        polygon_key = poly_key_raw,
        polygon_key_norm = poly_key_norm,
        matched_measurements = matched_meas,
        retained = FALSE,
        distinctiveness_score = NA_real_,
        cluster_label_raw = "unclustered_polygon",
        cluster_label = cluster_label_all,
        stringsAsFactors = FALSE
      )
      preview_df$retained[eligible] <- keep_mask_eligible
      preview_df$distinctiveness_score[eligible] <- distinctiveness
      preview_df$cluster_label_raw[eligible] <- labels_eligible_raw
      preview_df <- cbind(preview_df, geom_features[, intersect(c("geom_area", "geom_circularity", "geom_bbox_aspect"), names(geom_features)), drop = FALSE])
      preview_df <- preview_df[seq_len(min(500L, nrow(preview_df))), , drop = FALSE]

      feature_cols_final <- colnames(feature_df)
      xh$polygon_cluster_result <- list(
        cluster_label = cluster_label_all,
        cluster_id = cluster_num_all,
        eligible = eligible,
        matched_measurements = matched_meas,
        polygon_key = poly_key_raw,
        polygon_key_norm = poly_key_norm,
        polygon_key_field = join_field,
        csv_id_field = csv_id_col,
        csv_feature_cols = csv_feature_cols_used,
        feature_cols_used = feature_cols_final,
        feature_mode = mode,
        k = k_use,
        keep_prop = keep_prop,
        n_polygons = n_poly,
        n_eligible = sum(eligible),
        n_retained = sum(eligible & cluster_label_all != "unclustered_polygon", na.rm = TRUE),
        n_matched_measurements = sum(matched_meas),
        counts = counts_df,
        preview = preview_df,
        pca = pca_df
      )

      showNotification(
        sprintf(
          "Polygon clustering complete: %d clusters across %d/%d polygons; kept %d (prop=%.2f), matched=%d.",
          k_use, sum(eligible), n_poly, sum(keep_mask_eligible), keep_prop, sum(matched_meas)
        ),
        type = "message",
        duration = 7
      )
    }, ignoreInit = TRUE)

    observeEvent(input$use_polygon_clusters_for_labels, {
      if (is.null(xh$polygon_cluster_result) || is.null(xh$polygon_cluster_result$cluster_label)) {
        showNotification("Run polygon clustering first.", type = "warning", duration = 6)
        return()
      }
      cur_label_field <- as.character(input$polygon_label_field)[1]
      if (!is.na(cur_label_field) && nzchar(cur_label_field) && !identical(cur_label_field, "clustered_polygon_class")) {
        xh$polygon_label_field_before_cluster <- cur_label_field
      }
      out_col <- trimws(as.character(input$polygon_cluster_pdata_col)[1])
      if (is.na(out_col) || !nzchar(out_col)) out_col <- "polygon_cluster_class"
      out_col <- make.names(out_col)
      updateRadioButtons(session, "mapping_source", selected = "polygon")
      updateRadioButtons(session, "overlay_layer", selected = "polygon")
      updateCheckboxInput(session, "polygon_color_by_label", value = TRUE)
      session$onFlushed(function() {
        updateTextInput(session, "polygon_pdata_col", value = out_col)
        updateSelectInput(session, "polygon_label_field", selected = "clustered_polygon_class")
      }, once = TRUE)
      showNotification(
        sprintf("Polygon label field set to clustered classes. Mapping target pData column set to '%s' (Mapping & Save panel). Click 'Map selected source to pData' next.", out_col),
        type = "message",
        duration = 10
      )
    }, ignoreInit = TRUE)

    output$polygon_cluster_summary <- renderPrint({
      res <- xh$polygon_cluster_result
      if (is.null(res)) {
        cat("No polygon clustering results yet.\n")
        cat("Run clustering in this tab using QuPath measurements or geometry-only features.\n")
        return(invisible(NULL))
      }
      info <- list(
        feature_mode = res$feature_mode,
        polygon_id_field = res$polygon_key_field,
        csv_id_field = if (!is.null(res$csv_id_field)) res$csv_id_field else NA_character_,
        n_polygons = res$n_polygons,
        n_eligible = res$n_eligible,
        n_retained = res$n_retained,
        keep_prop = res$keep_prop,
        n_matched_measurements = res$n_matched_measurements,
        k = res$k,
        n_features_used = length(res$feature_cols_used),
        features_used_head = utils::head(res$feature_cols_used, 12)
      )
      print(info)
    })

    output$polygon_cluster_counts_table <- DT::renderDataTable({
      req(xh$polygon_cluster_result)
      counts_df <- xh$polygon_cluster_result$counts
      DT::datatable(counts_df, options = list(dom = "tip", pageLength = 8), rownames = FALSE)
    })

    output$polygon_cluster_preview_table <- DT::renderDataTable({
      req(xh$polygon_cluster_result)
      DT::datatable(xh$polygon_cluster_result$preview, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
    })

    output$polygon_cluster_plot <- renderPlot({
      req(xh$polygon_cluster_result)
      res <- xh$polygon_cluster_result
      pca_df <- res$pca
      if (is.null(pca_df) || !is.data.frame(pca_df) || nrow(pca_df) == 0) {
        cnt <- res$counts
        graphics::barplot(cnt$n_polygons, names.arg = cnt$cluster_label, las = 2, cex.names = 0.8, ylab = "n polygons", main = "Polygon clusters")
        return(invisible(NULL))
      }

      if (!all(c("PC1", "PC2") %in% names(pca_df))) {
        cnt <- res$counts
        graphics::barplot(cnt$n_polygons, names.arg = cnt$cluster_label, las = 2, cex.names = 0.8, ylab = "n polygons", main = "Polygon clusters")
        return(invisible(NULL))
      }

      pal <- get_discrete_palette(length(unique(pca_df$cluster_label)), input$cluster_palette)
      labs <- sort(unique(pca_df$cluster_label[pca_df$cluster_label != "unclustered_polygon"]))
      if (length(labs) == 0) labs <- sort(unique(pca_df$cluster_label))
      pal <- pal[seq_len(length(labs))]
      names(pal) <- labs

      p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2))
      if ("retained" %in% names(pca_df) && any(!pca_df$retained, na.rm = TRUE)) {
        p <- p + ggplot2::geom_point(
          data = pca_df[!pca_df$retained, , drop = FALSE],
          color = "grey75",
          alpha = 0.35,
          size = 0.8
        )
      }
      p +
        ggplot2::geom_point(
          data = if ("retained" %in% names(pca_df)) pca_df[pca_df$retained, , drop = FALSE] else pca_df,
          ggplot2::aes(color = cluster_label_raw),
          alpha = 0.8,
          size = 1.1
        ) +
        ggplot2::scale_color_manual(values = pal) +
        ggplot2::labs(title = "Polygon clustering (PCA space)", subtitle = sprintf("Retained per cluster: %.2f", if (!is.null(res$keep_prop)) res$keep_prop else 1), color = "Cluster") +
        ggplot2::theme_minimal(base_size = 11)
    })

    suggest_rgb_indices <- function(obj, n_channels = 3L, max_cells = 1.5e7) {
      n_channels <- suppressWarnings(as.integer(n_channels))
      if (!n_channels %in% c(2L, 3L)) n_channels <- 3L

      mz_axis <- suppressWarnings(as.numeric(Cardinal::mz(obj)))
      n_feat <- length(mz_axis)
      n_pix <- ncol(obj)
      if (n_feat < 2L || n_pix < 2L) return(integer(0))

      max_cells <- suppressWarnings(as.numeric(max_cells))
      if (!is.finite(max_cells) || max_cells <= 1e6) max_cells <- 1.5e7

      sample_pix <- floor(max_cells / max(1L, n_feat))
      sample_pix <- max(200L, sample_pix)
      sample_pix <- min(n_pix, sample_pix)
      if (sample_pix < 2L) sample_pix <- min(n_pix, 2L)

      pix_idx <- if (sample_pix < n_pix) sort(sample.int(n_pix, sample_pix)) else seq_len(n_pix)
      sp <- try(as.matrix(Cardinal::spectra(obj)[, pix_idx, drop = FALSE]), silent = TRUE)
      if (inherits(sp, "try-error") || nrow(sp) < 2L || ncol(sp) < 2L) return(integer(0))

      storage.mode(sp) <- "double"
      sp[!is.finite(sp)] <- NA_real_
      sp <- log1p(pmax(sp, 0))

      score <- apply(sp, 1, stats::sd, na.rm = TRUE)
      nz <- rowMeans(sp > 0, na.rm = TRUE)
      score[!is.finite(score)] <- 0
      nz[!is.finite(nz)] <- 0
      score <- score * (0.25 + 0.75 * pmin(1, nz * 2))

      if (length(mz_axis) >= length(score)) {
        valid_mz <- is.finite(mz_axis[seq_len(length(score))])
        score[!valid_mz] <- -Inf
      }

      ord <- order(score, decreasing = TRUE, na.last = NA)
      ord <- ord[is.finite(score[ord]) & score[ord] > 0]
      if (length(ord) < 2L) {
        ord <- which(is.finite(score))
        ord <- ord[order(score[ord], decreasing = TRUE)]
      }
      if (length(ord) < 2L) return(integer(0))

      k <- min(
        length(ord),
        max(20L, min(200L, as.integer(ceiling(sqrt(length(ord)) * 8))))
      )
      cand <- ord[seq_len(k)]
      cand_mat <- sp[cand, , drop = FALSE]
      keep <- rowSums(is.finite(cand_mat)) >= 5L
      cand <- cand[keep]
      cand_mat <- cand_mat[keep, , drop = FALSE]
      if (length(cand) < 2L) {
        return(ord[seq_len(min(length(ord), n_channels))])
      }

      corr <- suppressWarnings(stats::cor(t(cand_mat), use = "pairwise.complete.obs"))
      if (is.null(dim(corr))) {
        return(cand[seq_len(min(length(cand), n_channels))])
      }
      corr[!is.finite(corr)] <- 1
      diag(corr) <- 1

      cand_score <- score[cand]
      selected <- cand[which.max(cand_score)]
      while (length(selected) < n_channels && length(selected) < length(cand)) {
        sel_pos <- match(selected, cand)
        diversity <- 1 - rowMeans(abs(corr[, sel_pos, drop = FALSE]), na.rm = TRUE)
        diversity[!is.finite(diversity)] <- 0
        rank_score <- cand_score * pmax(0, diversity)
        rank_score[sel_pos] <- -Inf
        next_pos <- which.max(rank_score)
        if (!is.finite(rank_score[next_pos]) || rank_score[next_pos] == -Inf) break
        selected <- c(selected, cand[next_pos])
      }

      selected <- unique(selected)
      if (length(selected) < 2L) {
        selected <- unique(c(selected, ord[seq_len(min(2L, length(ord)))]))
      }
      selected[seq_len(min(length(selected), n_channels))]
    }

    transform_intensity <- function(x, trans) {
      if (identical(trans, "sqrt")) return(sqrt(pmax(x, 0)))
      if (identical(trans, "log1p")) return(log1p(pmax(x, 0)))
      if (identical(trans, "asinh")) return(asinh(x))
      x
    }

    clamp <- function(x, lo, hi) {
      pmax(lo, pmin(hi, x))
    }

    parse_numeric_tokens <- function(x) {
      if (is.null(x) || length(x) == 0) return(numeric(0))
      txt <- as.character(x)[1]
      if (is.null(txt) || is.na(txt) || !nzchar(trimws(txt))) return(numeric(0))
      toks <- trimws(unlist(strsplit(txt, "[,;[:space:]]+", perl = TRUE)))
      toks <- toks[nzchar(toks)]
      if (length(toks) == 0) return(numeric(0))
      vals <- suppressWarnings(as.numeric(toks))
      vals[is.finite(vals)]
    }

    safe_color <- function(col, fallback = "#73FFFF") {
      col <- as.character(col)[1]
      if (is.null(col) || is.na(col) || !nzchar(trimws(col))) return(fallback)
      ok <- try(grDevices::col2rgb(col), silent = TRUE)
      if (inherits(ok, "try-error")) fallback else col
    }

    normalize_labels <- function(x) {
      x <- as.character(x)
      x <- iconv(x, from = "", to = "UTF-8", sub = "")
      x <- gsub("[[:cntrl:]\u200B\u200C\u200D\uFEFF]", "", x, perl = TRUE)
      x <- gsub("\\s+", " ", x, perl = TRUE)
      x <- trimws(x)
      x[is.na(x) | !nzchar(x)] <- "NA"
      x
    }

    are_valid_colors <- function(x) {
      x <- as.character(x)
      ok <- rep(FALSE, length(x))
      keep <- !is.na(x) & nzchar(x)
      if (any(keep)) {
        ok[keep] <- vapply(
          x[keep],
          function(val) {
            tryCatch({
              grDevices::col2rgb(val)
              TRUE
            }, error = function(e) FALSE)
          },
          logical(1)
        )
      }
      ok
    }

    rescale01 <- function(x, enhance = TRUE) {
      out <- rep(NA_real_, length(x))
      ok <- is.finite(x)
      if (!any(ok)) {
        return(rep(0, length(x)))
      }

      x_work <- x
      if (isTRUE(enhance)) {
        q <- stats::quantile(x_work[ok], probs = c(0.01, 0.99), na.rm = TRUE, names = FALSE, type = 8)
        if (all(is.finite(q)) && q[2] > q[1]) {
          x_work <- pmin(pmax(x_work, q[1]), q[2])
          lo <- q[1]
          hi <- q[2]
        } else {
          rr <- range(x_work[ok], na.rm = TRUE)
          lo <- rr[1]
          hi <- rr[2]
        }
      } else {
        rr <- range(x_work[ok], na.rm = TRUE)
        lo <- rr[1]
        hi <- rr[2]
      }

      if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
        out[ok] <- 0
        return(out)
      }

      out[ok] <- (x_work[ok] - lo) / (hi - lo)
      pmin(1, pmax(0, out))
    }

    gaussian_smooth_matrix <- function(mat, sigma) {
      sigma <- as.numeric(sigma)
      if (!is.finite(sigma) || sigma <= 0) {
        return(mat)
      }

      radius <- max(1L, as.integer(ceiling(3 * sigma)))
      x <- seq.int(-radius, radius)
      k <- exp(-(x * x) / (2 * sigma * sigma))
      k <- k / sum(k)

      smooth_vec <- function(v) {
        w <- as.numeric(!is.na(v))
        v0 <- v
        v0[is.na(v0)] <- 0
        num <- as.numeric(stats::filter(v0, k, sides = 2))
        den <- as.numeric(stats::filter(w, k, sides = 2))
        out <- num / den
        out[!is.finite(out) | den <= 0] <- NA_real_
        out
      }

      tmp <- t(apply(mat, 1, smooth_vec))
      out <- apply(tmp, 2, smooth_vec)
      if (!is.matrix(out)) {
        out <- matrix(out, nrow = nrow(mat), ncol = ncol(mat))
      }
      out
    }

    get_msi_palette <- function(name) {
      pal <- try(grDevices::hcl.colors(256, name), silent = TRUE)
      if (inherits(pal, "try-error") || length(pal) < 2) {
        grDevices::hcl.colors(256, "Inferno")
      } else {
        pal
      }
    }

    get_discrete_palette <- function(n, name) {
      if (n <= 0) return(character(0))
      pal <- try(grDevices::hcl.colors(max(n, 3), name), silent = TRUE)
      if (inherits(pal, "try-error") || length(pal) < n) {
        pal <- grDevices::hcl.colors(max(n, 3), "Set 2")
      }
      pal[seq_len(n)]
    }

    resolve_polygon_axis_mode <- function(poly_sf, mode = "auto", hist_img = NULL) {
      mode_chr <- tolower(trimws(if (is.null(mode) || length(mode) == 0) "auto" else as.character(mode[[1]])))
      if (mode_chr %in% c("xy", "yx")) return(mode_chr)
      if (is.null(hist_img)) return("xy")

      info <- try(magick::image_info(hist_img)[1, ], silent = TRUE)
      if (inherits(info, "try-error")) return("xy")
      img_w <- as.numeric(info$width)
      img_h <- as.numeric(info$height)
      if (!is.finite(img_w) || !is.finite(img_h) || img_w <= 0 || img_h <= 0) return("xy")

      bb <- try(sf::st_bbox(poly_sf), silent = TRUE)
      if (inherits(bb, "try-error")) return("xy")
      w <- as.numeric(bb$xmax - bb$xmin)
      h <- as.numeric(bb$ymax - bb$ymin)
      if (!is.finite(w) || !is.finite(h) || w <= 0 || h <= 0) return("xy")

      ar_img <- img_w / img_h
      d_xy <- abs(log((w / h) / ar_img))
      d_yx <- abs(log((h / w) / ar_img))
      if (is.finite(d_yx) && is.finite(d_xy) && d_yx < d_xy) "yx" else "xy"
    }

    get_histology_image_optional <- function() {
      if (is.null(input$histology_upload)) return(NULL)
      out <- try(histology_image(), silent = TRUE)
      if (inherits(out, "try-error")) NULL else out
    }

    get_overlay_source_dims <- function() {
      img <- get_histology_image_optional()
      if (!is.null(img)) {
        info <- try(magick::image_info(img)[1, ], silent = TRUE)
        if (!inherits(info, "try-error")) {
          w <- as.numeric(info$width)
          h <- as.numeric(info$height)
          if (is.finite(w) && is.finite(h) && w > 0 && h > 0) {
            return(list(width = w, height = h, source = "histology"))
          }
        }
      }

      if (!is.null(input$cluster_overlay_upload)) {
        img2 <- try(cluster_overlay_image(), silent = TRUE)
        if (!inherits(img2, "try-error") && !is.null(img2)) {
          info2 <- try(magick::image_info(img2)[1, ], silent = TRUE)
          if (!inherits(info2, "try-error")) {
            w2 <- as.numeric(info2$width)
            h2 <- as.numeric(info2$height)
            if (is.finite(w2) && is.finite(h2) && w2 > 0 && h2 > 0) {
              return(list(width = w2, height = h2, source = "cluster"))
            }
          }
        }
      }
      NULL
    }

    to_bool <- function(x, default = FALSE) {
      if (is.logical(x) && length(x) > 0 && !is.na(x[1])) return(isTRUE(x[1]))
      sx <- tolower(trimws(as.character(x)[1]))
      if (sx %in% c("true", "t", "1", "yes", "y", "on")) return(TRUE)
      if (sx %in% c("false", "f", "0", "no", "n", "off")) return(FALSE)
      default
    }

    current_registration_params <- reactive({
      obj <- try(msi_for_pdata(), silent = TRUE)
      mz_axis <- numeric(0)
      if (!inherits(obj, "try-error")) {
        mz_axis <- suppressWarnings(as.numeric(Cardinal::mz(obj)))
      }

      mz_sel_idx <- suppressWarnings(as.integer(input$mz_select))
      mz_sel_val <- NA_real_
      if (length(mz_sel_idx) == 1L && is.finite(mz_sel_idx) && mz_sel_idx >= 1L && mz_sel_idx <= length(mz_axis)) {
        mz_sel_val <- mz_axis[mz_sel_idx]
      }

      rgb_sel_source <- isolate(xh$rgb_mz_applied)
      if (is.null(rgb_sel_source) || length(rgb_sel_source) == 0) rgb_sel_source <- input$rgb_mz_select
      rgb_sel_idx <- suppressWarnings(as.integer(rgb_sel_source))
      rgb_sel_idx <- rgb_sel_idx[is.finite(rgb_sel_idx) & rgb_sel_idx >= 1L & rgb_sel_idx <= length(mz_axis)]
      rgb_sel_idx <- unique(rgb_sel_idx)
      rgb_sel_vals <- if (length(rgb_sel_idx) > 0) mz_axis[rgb_sel_idx] else numeric(0)

      list(
        msi_plot_mode = input$msi_plot_mode,
        mz_select = input$mz_select,
        mz_value = mz_sel_val,
        rgb_mz_select = paste(as.character(rgb_sel_idx), collapse = ","),
        rgb_mz_values = paste(format(rgb_sel_vals, digits = 10, scientific = FALSE, trim = TRUE), collapse = ","),
        rgb_auto_n = input$rgb_auto_n,
        overlay_pdata_field = input$overlay_pdata_field,
        overlay_scale_mode = input$overlay_scale_mode,
        histology_um_per_px = input$histology_um_per_px,
        msi_um_per_px = input$msi_um_per_px,
        histology_resample_factor = input$histology_resample_factor,
        scale_x = input$scale_x,
        scale_y = input$scale_y,
        rotate_deg = input$rotate_deg,
        translate_x = input$translate_x,
        translate_y = input$translate_y,
        flip_histology_y = isTRUE(input$flip_histology_y),
        histology_alpha = input$histology_alpha,
        cluster_alpha = input$cluster_alpha,
        polygon_axis_mode = input$polygon_axis_mode,
        polygon_color_by_label = isTRUE(input$polygon_color_by_label),
        polygon_outline_color = input$polygon_outline_color,
        polygon_linewidth = input$polygon_linewidth,
        intensity_transform = input$intensity_transform,
        msi_palette = input$msi_palette,
        cluster_palette = input$cluster_palette,
        rgb_render_mode = input$rgb_render_mode,
        rgb_bg_cutoff = input$rgb_bg_cutoff,
        optimize_edge_band = input$optimize_edge_band,
        stat_fit_metric_mode = input$stat_fit_metric_mode,
        stat_fit_group_field = input$stat_fit_group_field,
        stat_fit_outside_mode = input$stat_fit_outside_mode,
        stat_fit_objective = input$stat_fit_objective,
        stat_fit_bbox_pad = input$stat_fit_bbox_pad,
        stat_fit_top_n = input$stat_fit_top_n,
        show_fit_info = isTRUE(input$show_fit_info),
        enhance_contrast = isTRUE(input$enhance_contrast),
        gaussian_smooth = isTRUE(input$gaussian_smooth),
        gaussian_sigma = input$gaussian_sigma
      )
    })

    apply_resolution_scale <- function(show_message = FALSE) {
      h <- suppressWarnings(as.numeric(input$histology_um_per_px))
      m <- suppressWarnings(as.numeric(input$msi_um_per_px))
      f <- suppressWarnings(as.numeric(input$histology_resample_factor))
      if (!is.finite(f) || f <= 0) f <- 1
      if (!is.finite(h) || !is.finite(m) || h <= 0 || m <= 0) {
        if (isTRUE(show_message)) {
          showNotification("Resolution values must be positive numbers.", type = "error", duration = 6)
        }
        return(FALSE)
      }
      ratio <- (h * f) / m
      ratio_use <- clamp(ratio, 0.001, 50)
      updateSliderInput(session, "scale_x", value = ratio_use)
      updateSliderInput(session, "scale_y", value = ratio_use)
      updateNumericInput(session, "scale_x_num", value = ratio_use)
      updateNumericInput(session, "scale_y_num", value = ratio_use)
      if (isTRUE(show_message)) {
        showNotification(sprintf("Scale set to %.6f ((histology_um_per_px * export_factor) / msi_um_per_px).", ratio), type = "message", duration = 6)
      }
      TRUE
    }

    make_msi_raster <- reactive({
      req(msi_for_pdata())

      obj <- msi_for_pdata()
      cd <- as.data.frame(Cardinal::coord(obj))
      validate(need(all(c("x", "y") %in% names(cd)), "MSI coordinates must contain x and y."))
      validate(need(nrow(cd) > 0, "MSI data has no coordinates."))

      x_norm <- as.integer(cd$x - min(cd$x, na.rm = TRUE) + 1L)
      y_norm <- as.integer(cd$y - min(cd$y, na.rm = TRUE) + 1L)
      nx <- max(x_norm, na.rm = TRUE)
      ny <- max(y_norm, na.rm = TRUE)
      row_idx <- ny - y_norm + 1L

      smooth_sigma <- suppressWarnings(as.numeric(input$gaussian_sigma))
      if (!is.finite(smooth_sigma) || smooth_sigma <= 0) smooth_sigma <- 1

      build_numeric_signal <- function(vals, apply_transform = TRUE, apply_enhance = TRUE) {
        vals <- as.numeric(vals)
        if (isTRUE(apply_transform)) {
          vals <- transform_intensity(vals, input$intensity_transform)
        }

        mat <- matrix(NA_real_, nrow = ny, ncol = nx)
        mat[cbind(row_idx, x_norm)] <- vals

        if (isTRUE(input$gaussian_smooth)) {
          mat <- gaussian_smooth_matrix(mat, smooth_sigma)
        }

        sig <- matrix(
          rescale01(as.vector(mat), enhance = isTRUE(input$enhance_contrast)),
          nrow = ny,
          ncol = nx
        )
        if (!isTRUE(apply_enhance)) {
          sig <- matrix(
            rescale01(as.vector(mat), enhance = FALSE),
            nrow = ny,
            ncol = nx
          )
        }
        sig
      }

      build_numeric_raster <- function(vals, apply_transform = TRUE) {
        mat_scaled <- build_numeric_signal(vals, apply_transform = apply_transform, apply_enhance = TRUE)

        pal <- get_msi_palette(input$msi_palette)
        bin <- pmax(1L, pmin(256L, as.integer(mat_scaled * 255) + 1L))
        col_mat <- matrix(pal[bin], nrow = ny, ncol = nx)
        col_mat[is.na(mat_scaled)] <- "#00000000"
        as.raster(col_mat)
      }

      build_categorical_signal <- function(labels) {
        labs <- normalize_labels(labels)
        lab_mat <- matrix(NA_character_, nrow = ny, ncol = nx)
        lab_mat[cbind(row_idx, x_norm)] <- labs

        valid <- !is.na(lab_mat)
        right <- cbind(lab_mat[, -1, drop = FALSE], NA_character_)
        left <- cbind(NA_character_, lab_mat[, -ncol(lab_mat), drop = FALSE])
        down <- rbind(lab_mat[-1, , drop = FALSE], rep(NA_character_, ncol(lab_mat)))
        up <- rbind(rep(NA_character_, ncol(lab_mat)), lab_mat[-nrow(lab_mat), , drop = FALSE])

        b <- matrix(0, nrow = ny, ncol = nx)
        b <- b + ((lab_mat != right) & valid & !is.na(right))
        b <- b + ((lab_mat != left) & valid & !is.na(left))
        b <- b + ((lab_mat != down) & valid & !is.na(down))
        b <- b + ((lab_mat != up) & valid & !is.na(up))
        b[!is.finite(b)] <- 0
        b <- b / 4
        b[!valid] <- 0

        if (isTRUE(input$gaussian_smooth)) {
          b <- gaussian_smooth_matrix(b, smooth_sigma)
        }

        b <- matrix(rescale01(as.vector(b), enhance = FALSE), nrow = ny, ncol = nx)
        b[!is.finite(b)] <- 0
        b
      }

      build_rgb_raster <- function(vals_r, vals_g, vals_b) {
        to_channel <- function(vals) {
          vals <- transform_intensity(as.numeric(vals), input$intensity_transform)
          mat <- matrix(NA_real_, nrow = ny, ncol = nx)
          mat[cbind(row_idx, x_norm)] <- vals
          if (isTRUE(input$gaussian_smooth)) {
            mat <- gaussian_smooth_matrix(mat, smooth_sigma)
          }
          matrix(
            rescale01(as.vector(mat), enhance = isTRUE(input$enhance_contrast)),
            nrow = ny,
            ncol = nx
          )
        }

        r01 <- to_channel(vals_r)
        g01 <- to_channel(vals_g)
        b01 <- to_channel(vals_b)
        rgb_mode <- tolower(trimws(as.character(input$rgb_render_mode)[1]))
        if (!rgb_mode %in% c("additive", "dominant")) rgb_mode <- "dominant"
        bg_cut <- suppressWarnings(as.numeric(input$rgb_bg_cutoff))
        if (!is.finite(bg_cut)) bg_cut <- 0.05
        bg_cut <- clamp(bg_cut, 0, 1)

        r <- ifelse(is.finite(r01), r01, 0)
        g <- ifelse(is.finite(g01), g01, 0)
        b <- ifelse(is.finite(b01), b01, 0)
        maxv <- pmax(r, g, b)
        mask <- (is.finite(r01) | is.finite(g01) | is.finite(b01)) & maxv >= bg_cut

        if (identical(rgb_mode, "dominant")) {
          rgb_stack <- cbind(as.vector(r), as.vector(g), as.vector(b))
          dom <- max.col(rgb_stack, ties.method = "first")
          inten <- as.vector(maxv)
          rr_v <- gg_v <- bb_v <- rep(0, length(dom))
          rr_v[dom == 1] <- inten[dom == 1]
          gg_v[dom == 2] <- inten[dom == 2]
          bb_v[dom == 3] <- inten[dom == 3]
          rr <- matrix(rr_v, nrow = ny, ncol = nx)
          gg <- matrix(gg_v, nrow = ny, ncol = nx)
          bb <- matrix(bb_v, nrow = ny, ncol = nx)
        } else {
          rr <- r
          gg <- g
          bb <- b
        }

        rr <- as.integer(pmax(0, pmin(255, round(rr * 255))))
        gg <- as.integer(pmax(0, pmin(255, round(gg * 255))))
        bb <- as.integer(pmax(0, pmin(255, round(bb * 255))))

        col_hex <- matrix(sprintf("#%02X%02X%02X", rr, gg, bb), nrow = ny, ncol = nx)
        col_hex[!mask] <- "#00000000"
        as.raster(col_hex)
      }

      mode_req <- tolower(trimws(as.character(input$msi_plot_mode)[1]))
      if (!mode_req %in% c("mz", "rgb", "pdata")) mode_req <- "mz"

      # 1) pData mode
      if (identical(mode_req, "pdata")) {
        pd <- as.data.frame(Cardinal::pData(obj))
        field <- as.character(input$overlay_pdata_field)[1]
        if (!is.null(field) && length(field) > 0 && !is.na(field) && nzchar(field) && field %in% colnames(pd)) {
          vals <- pd[[field]]

          if (is.numeric(vals) || is.integer(vals)) {
            ras <- build_numeric_raster(vals, apply_transform = FALSE)
            sig_opt <- build_numeric_signal(vals, apply_transform = FALSE, apply_enhance = FALSE)
            return(list(
              raster = ras,
              nx = nx,
              ny = ny,
              mz_selected = NA_real_,
              mz_index = NA_integer_,
              rgb_mz = NULL,
              pdata_field = field,
              mode = "pdata",
              mode_requested = mode_req,
              display_label = paste0("pData: ", field),
              opt_signal = sig_opt,
              x_norm = x_norm,
              y_norm = y_norm,
              row_idx = row_idx
            ))
          }

          labs <- normalize_labels(vals)
          sig_opt <- build_categorical_signal(vals)
          lev <- unique(labs)
          col_map <- stats::setNames(rep("grey70", length(lev)), lev)
          valid_col <- are_valid_colors(lev)
          if (any(valid_col)) col_map[valid_col] <- lev[valid_col]
          if (any(!valid_col)) col_map[!valid_col] <- get_discrete_palette(sum(!valid_col), input$cluster_palette)

          lab_mat <- matrix(NA_character_, nrow = ny, ncol = nx)
          lab_mat[cbind(row_idx, x_norm)] <- labs
          col_mat <- matrix("#00000000", nrow = ny, ncol = nx)
          keep <- !is.na(lab_mat)
          col_mat[keep] <- col_map[lab_mat[keep]]

          return(list(
            raster = as.raster(col_mat),
            nx = nx,
            ny = ny,
            mz_selected = NA_real_,
            mz_index = NA_integer_,
            rgb_mz = NULL,
            pdata_field = field,
            mode = "pdata",
            mode_requested = mode_req,
            display_label = paste0("pData: ", field),
            opt_signal = sig_opt,
            x_norm = x_norm,
            y_norm = y_norm,
            row_idx = row_idx
          ))
        }
      }

      mzv <- try(Cardinal::mz(obj), silent = TRUE)
      mzv <- suppressWarnings(as.numeric(mzv))
      mzv <- mzv[is.finite(mzv)]
      validate(need(length(mzv) > 0, "MSI data has no m/z axis for selected display mode."))

      # 2) RGB mode
      if (identical(mode_req, "rgb")) {
        idx_rgb <- suppressWarnings(as.integer(xh$rgb_mz_applied))
        idx_rgb <- unique(idx_rgb[is.finite(idx_rgb) & idx_rgb >= 1L & idx_rgb <= length(mzv)])
        if (length(idx_rgb) %in% c(2, 3)) {
          sp <- try(as.matrix(Cardinal::spectra(obj)[idx_rgb, , drop = FALSE]), silent = TRUE)
          if (!inherits(sp, "try-error") && nrow(sp) == length(idx_rgb) && ncol(sp) == ncol(obj)) {
            vals_r <- sp[1, ]
            vals_g <- sp[2, ]
            vals_b <- if (length(idx_rgb) == 3) sp[3, ] else rep(0, ncol(obj))
            ras <- build_rgb_raster(vals_r, vals_g, vals_b)
            rgb_mz <- mzv[idx_rgb]
            rgb_label <- if (length(rgb_mz) == 2) {
              sprintf("RGB: R=%.5f, G=%.5f, B=0", rgb_mz[1], rgb_mz[2])
            } else {
              sprintf("RGB: R=%.5f, G=%.5f, B=%.5f", rgb_mz[1], rgb_mz[2], rgb_mz[3])
            }
            return(list(
              raster = ras,
              nx = nx,
              ny = ny,
              mz_selected = NA_real_,
              mz_index = NA_integer_,
              rgb_mz = rgb_mz,
              pdata_field = NULL,
              mode = "rgb",
              mode_requested = mode_req,
              display_label = rgb_label,
              opt_signal = NULL,
              x_norm = x_norm,
              y_norm = y_norm,
              row_idx = row_idx
            ))
          }
        }
        # In RGB mode, do not fall back to single-ion rendering while the user is
        # editing channel selections. This avoids repeated expensive redraws from
        # changes to the standalone m/z selector used only for channel picking.
        empty_mat <- matrix("#00000000", nrow = ny, ncol = nx)
        return(list(
          raster = as.raster(empty_mat),
          nx = nx,
          ny = ny,
          mz_selected = NA_real_,
          mz_index = NA_integer_,
          rgb_mz = numeric(0),
          pdata_field = NULL,
          mode = "rgb",
          mode_requested = mode_req,
          display_label = "RGB: select 2-3 channels and click Apply RGB channels",
          opt_signal = NULL,
          x_norm = x_norm,
          y_norm = y_norm,
          row_idx = row_idx
        ))
      }

      # 3) Single m/z mode (default fallback)
      idx <- suppressWarnings(as.integer(input$mz_select))
      if (!is.finite(idx) || idx < 1L || idx > length(mzv)) idx <- 1L

      vals <- try(as.numeric(Cardinal::spectra(obj)[idx, ]), silent = TRUE)
      validate(need(!inherits(vals, "try-error"), "Could not extract MSI intensities for selected m/z."))
      ras <- build_numeric_raster(vals, apply_transform = TRUE)

      list(
        raster = ras,
        nx = nx,
        ny = ny,
        mz_selected = mzv[idx],
        mz_index = idx,
        rgb_mz = NULL,
        pdata_field = NULL,
        mode = "mz",
        mode_requested = mode_req,
        display_label = sprintf("m/z %.5f", mzv[idx]),
        opt_signal = NULL,
        x_norm = x_norm,
        y_norm = y_norm,
        row_idx = row_idx
      )
    })

    image_to_raster_rgba <- function(img, alpha_scale = 1) {
      dat <- magick::image_data(img, channels = "rgba")
      h <- dim(dat)[2]
      w <- dim(dat)[3]

      r <- matrix(as.character(dat[1, , ]), nrow = h, ncol = w)
      g <- matrix(as.character(dat[2, , ]), nrow = h, ncol = w)
      b <- matrix(as.character(dat[3, , ]), nrow = h, ncol = w)
      a_hex <- matrix(as.character(dat[4, , ]), nrow = h, ncol = w)

      a_num <- suppressWarnings(strtoi(a_hex, base = 16L))
      if (anyNA(a_num)) a_num[is.na(a_num)] <- 255L
      a_num <- as.integer(pmax(0, pmin(255, round(a_num * alpha_scale))))
      a <- sprintf("%02X", a_num)

      as.raster(matrix(paste0("#", r, g, b, a), nrow = h, ncol = w))
    }

    transform_overlay_image <- function(img, use_point_filter = FALSE, alpha_scale = 1, scale_correction = NULL) {
      req(make_msi_raster())
      img_work <- img

      if (isTRUE(input$flip_histology_y)) img_work <- magick::image_flip(img_work)

      info0 <- magick::image_info(img_work)[1, ]
      nx <- make_msi_raster()$nx
      ny <- make_msi_raster()$ny

      scale_mode <- tolower(trimws(as.character(input$overlay_scale_mode)[1]))
      if (!scale_mode %in% c("absolute", "fit")) scale_mode <- "absolute"
      if (is.null(scale_correction)) {
        scale_correction <- list(fx = 1, fy = 1)
      }
      corr_x <- suppressWarnings(as.numeric(scale_correction$fx))
      corr_y <- suppressWarnings(as.numeric(scale_correction$fy))
      if (!is.finite(corr_x) || corr_x <= 0) corr_x <- 1
      if (!is.finite(corr_y) || corr_y <= 0) corr_y <- 1
      if (identical(scale_mode, "fit")) {
        fit_scale <- min(nx / info0$width, ny / info0$height)
        target_w <- max(1L, as.integer(round(info0$width * fit_scale * input$scale_x)))
        target_h <- max(1L, as.integer(round(info0$height * fit_scale * input$scale_y)))
      } else {
        target_w <- max(1L, as.integer(round(info0$width * input$scale_x * corr_x)))
        target_h <- max(1L, as.integer(round(info0$height * input$scale_y * corr_y)))
      }

      if (isTRUE(use_point_filter)) {
        img_work <- magick::image_resize(img_work, paste0(target_w, "x", target_h, "!"), filter = "point")
      } else {
        img_work <- magick::image_resize(img_work, paste0(target_w, "x", target_h, "!"))
      }
      img_work <- magick::image_background(img_work, "none")
      if (!is.null(input$rotate_deg) && is.finite(input$rotate_deg) && input$rotate_deg != 0) {
        img_work <- magick::image_rotate(img_work, input$rotate_deg)
      }

      ras <- image_to_raster_rgba(img_work, alpha_scale = alpha_scale)
      list(
        raster = ras,
        width = ncol(ras),
        height = nrow(ras),
        source_width = as.numeric(info0$width),
        source_height = as.numeric(info0$height),
        scale_correction_x = corr_x,
        scale_correction_y = corr_y
      )
    }

    transformed_overlay <- reactive({
      req(input$overlay_layer)
      req(make_msi_raster())
      scale_corr <- list(fx = 1, fy = 1)

      if (identical(input$overlay_layer, "polygon")) {
        req(input$polygon_file)
        req(polygon_data())

        msi <- make_msi_raster()
        poly <- polygon_data()
        poly$map_label <- get_polygon_labels(poly, input$polygon_label_field)
        axis_mode <- resolve_polygon_axis_mode(poly, input$polygon_axis_mode, get_histology_image_optional())
        src_dims <- get_overlay_source_dims()
        poly_t <- transform_polygon_sf(
          poly_sf = poly,
          nx = msi$nx,
          ny = msi$ny,
          scale_x = input$scale_x,
          scale_y = input$scale_y,
          translate_x = input$translate_x,
          translate_y = input$translate_y,
          rotate_deg = input$rotate_deg,
          flip_y = isTRUE(input$flip_histology_y),
          swap_xy = identical(axis_mode, "yx"),
          scale_mode = input$overlay_scale_mode,
          source_width = if (!is.null(src_dims)) src_dims$width else NA_real_,
          source_height = if (!is.null(src_dims)) src_dims$height else NA_real_
        )

        return(list(
          layer = "polygon",
          polygons = poly_t,
          alpha_used = NA_real_,
          axis_mode = axis_mode,
          polygon_source_dim = if (!is.null(src_dims)) paste0(src_dims$width, "x", src_dims$height, " (", src_dims$source, ")") else NA_character_
        ))
      }

      if (identical(input$overlay_layer, "cluster")) {
        validate(need(!is.null(input$cluster_overlay_upload), "Upload a cluster-color image or switch overlay to Histology image."))
        alpha_use <- if (is.finite(input$cluster_alpha)) input$cluster_alpha else 0.7
        out <- transform_overlay_image(cluster_overlay_image(), use_point_filter = TRUE, alpha_scale = alpha_use, scale_correction = scale_corr)
        out$alpha_used <- alpha_use
        out$layer <- "cluster"
        return(out)
      }

      alpha_use <- if (is.finite(input$histology_alpha)) input$histology_alpha else 0.5
      out <- transform_overlay_image(histology_image(), use_point_filter = FALSE, alpha_scale = alpha_use, scale_correction = scale_corr)
      out$alpha_used <- alpha_use
      out$layer <- "histology"
      out
    })

    cluster_mapping_payload <- reactive({
      req(make_msi_raster())
      req(input$cluster_overlay_upload)

      scale_corr <- list(fx = 1, fy = 1)
      tr <- transform_overlay_image(cluster_overlay_image(), use_point_filter = TRUE, alpha_scale = 1, scale_correction = scale_corr)
      msi <- make_msi_raster()

      nx <- msi$nx
      ny <- msi$ny
      w <- tr$width
      h <- tr$height

      cx <- (nx / 2) + input$translate_x
      cy <- (ny / 2) - input$translate_y
      xleft <- cx - w / 2
      xright <- cx + w / 2
      ybottom <- cy - h / 2
      ytop <- cy + h / 2

      xpix <- msi$x_norm
      ypix <- msi$y_norm

      u <- (xpix - xleft) / (xright - xleft)
      v <- (ytop - ypix) / (ytop - ybottom)
      inside <- is.finite(u) & is.finite(v) & u >= 0 & u <= 1 & v >= 0 & v <= 1

      col_idx <- pmin(w, pmax(1L, as.integer(floor(u * (w - 1L)) + 1L)))
      row_idx <- pmin(h, pmax(1L, as.integer(floor(v * (h - 1L)) + 1L)))

      col_mat <- as.matrix(tr$raster)
      sampled <- rep(NA_character_, length(xpix))
      sampled[inside] <- col_mat[cbind(row_idx[inside], col_idx[inside])]

      list(sampled_rgba = sampled)
    })

    extract_reference_cluster_colors <- function(img, max_colors = 128L) {
      dat <- magick::image_data(img, channels = "rgba")
      r <- as.character(dat[1, , ])
      g <- as.character(dat[2, , ])
      b <- as.character(dat[3, , ])
      a <- as.character(dat[4, , ])
      a_num <- suppressWarnings(strtoi(a, base = 16L))
      keep <- !is.na(a_num) & a_num > 0L
      if (!any(keep)) return(character(0))

      rgb_hex <- paste0("#", r, g, b)
      tab <- sort(table(rgb_hex[keep]), decreasing = TRUE)
      cols <- names(tab)
      if (length(cols) > max_colors) cols <- cols[seq_len(max_colors)]
      cols
    }

    nearest_palette_map <- function(query_hex, ref_hex) {
      query_hex <- as.character(query_hex)
      out <- rep(NA_character_, length(query_hex))
      ok <- !is.na(query_hex)
      if (!any(ok) || length(ref_hex) == 0) return(out)

      uq <- unique(query_hex[ok])
      ref_rgb <- t(grDevices::col2rgb(ref_hex))
      lut <- character(length(uq))
      names(lut) <- uq

      for (u in uq) {
        if (u %in% ref_hex) {
          lut[u] <- u
        } else {
          uu <- as.numeric(grDevices::col2rgb(u))
          d <- rowSums((ref_rgb - matrix(uu, nrow = nrow(ref_rgb), ncol = 3, byrow = TRUE))^2)
          lut[u] <- ref_hex[which.min(d)]
        }
      }

      out[ok] <- lut[query_hex[ok]]
      out
    }

    build_position_pdata <- function(obj, pd_df) {
      run_vec <- if ("run" %in% colnames(pd_df)) pd_df$run else Cardinal::run(obj)
      extra_pd <- pd_df[, !colnames(pd_df) %in% c("x", "y", "run"), drop = FALSE]
      new_pd <- try(
        Cardinal::PositionDataFrame(
          run = run_vec,
          coord = Cardinal::coord(obj),
          extra_pd
        ),
        silent = TRUE
      )
      validate(need(!inherits(new_pd, "try-error"), "Could not construct PositionDataFrame for mapped pData."))
      new_pd
    }

    normalize_crs <- function(crs_obj) {
      if (inherits(crs_obj, "try-error") || is.null(crs_obj)) return(NULL)
      if (inherits(crs_obj, "crs")) return(crs_obj)
      if (is.numeric(crs_obj) || is.character(crs_obj)) return(crs_obj)
      NULL
    }

    get_polygon_labels <- function(poly, label_field) {
      n <- nrow(poly)
      if (n == 0) return(character(0))
      if (identical(label_field, "clustered_polygon_class")) {
        cl <- NULL
        res <- xh$polygon_cluster_result
        if (!is.null(res)) {
          # Prefer stable-ID mapping when the clustering run stored a polygon key field.
          key_field <- as.character(res$polygon_key_field)[1]
          key_ref <- res$polygon_key_norm
          cl_ref <- as.character(res$cluster_label)
          if (is.null(key_ref) && !is.null(res$polygon_key)) {
            key_ref <- normalize_join_key(res$polygon_key)
          }
          if (!is.null(key_field) && nzchar(key_field) && !identical(key_field, "row_index") &&
              key_field %in% colnames(poly) && !is.null(key_ref) && length(key_ref) == length(cl_ref) && length(key_ref) > 0) {
            key_now <- normalize_join_key(poly[[key_field]])
            m <- match(key_now, key_ref)
            cl <- rep(NA_character_, n)
            hit <- !is.na(m)
            if (any(hit)) cl[hit] <- cl_ref[m[hit]]
          } else {
            cl <- cl_ref
          }
        }
        cl <- as.character(cl)
        if (length(cl) != n) {
          cl <- c(cl, rep(NA_character_, max(0L, n - length(cl))))
          cl <- cl[seq_len(n)]
        }
        cl[is.na(cl) | trimws(cl) == ""] <- "unclustered_polygon"
        return(make.names(cl, unique = FALSE))
      }
      if (is.null(label_field) || identical(label_field, "row_index") || !label_field %in% colnames(poly)) {
        return(paste0("polygon_", sprintf("%03d", seq_len(n))))
      }
      vals <- as.character(poly[[label_field]])
      vals[is.na(vals) | trimws(vals) == ""] <- paste0("polygon_", sprintf("%03d", which(is.na(vals) | trimws(vals) == "")))
      make.names(vals, unique = TRUE)
    }

    is_cell_like_polygon_label <- function(lbl) {
      s <- tolower(trimws(as.character(lbl)))
      s[is.na(s)] <- ""
      # Common non-cell/background annotation tokens from exported polygon workflows.
      non_cell <- grepl(
        "^(line[0-9_.-]*|background|bg|outside|outer|border|boundary|frame|artifact|noncell|notcell|mask|roi|region)$",
        s
      )
      nzchar(s) & !non_cell
    }

    transform_polygon_sf <- function(poly_sf, nx, ny, scale_x, scale_y, translate_x, translate_y, rotate_deg, flip_y = FALSE, swap_xy = FALSE, scale_mode = "absolute", source_width = NA_real_, source_height = NA_real_) {
      geom_type <- as.character(sf::st_geometry_type(poly_sf))
      keep <- geom_type %in% c("POLYGON", "MULTIPOLYGON")
      poly <- poly_sf[keep, , drop = FALSE]
      validate(need(nrow(poly) > 0, "No polygon geometries found in input file."))

      bb <- sf::st_bbox(poly)
      if (isTRUE(swap_xy)) {
        bb_xmin <- as.numeric(bb$ymin)
        bb_xmax <- as.numeric(bb$ymax)
        bb_ymin <- as.numeric(bb$xmin)
        bb_ymax <- as.numeric(bb$xmax)
      } else {
        bb_xmin <- as.numeric(bb$xmin)
        bb_xmax <- as.numeric(bb$xmax)
        bb_ymin <- as.numeric(bb$ymin)
        bb_ymax <- as.numeric(bb$ymax)
      }
      bb_w <- as.numeric(bb_xmax - bb_xmin)
      bb_h <- as.numeric(bb_ymax - bb_ymin)
      validate(need(is.finite(bb_w) && is.finite(bb_h) && bb_w > 0 && bb_h > 0, "Polygon bounding box is invalid."))

      has_source_dims <- is.finite(source_width) && is.finite(source_height) && source_width > 0 && source_height > 0
      if (has_source_dims) {
        if (isTRUE(swap_xy)) {
          frame_w <- as.numeric(source_height)
          frame_h <- as.numeric(source_width)
        } else {
          frame_w <- as.numeric(source_width)
          frame_h <- as.numeric(source_height)
        }
      } else {
        frame_w <- bb_w
        frame_h <- bb_h
      }
      validate(need(is.finite(frame_w) && is.finite(frame_h) && frame_w > 0 && frame_h > 0, "Polygon reference frame is invalid."))

      mode <- tolower(trimws(as.character(scale_mode)[1]))
      if (!mode %in% c("absolute", "fit")) mode <- "absolute"
      if (identical(mode, "fit")) {
        fit_scale <- min(nx / frame_w, ny / frame_h)
        sx <- fit_scale * scale_x
        sy <- fit_scale * scale_y
      } else {
        sx <- as.numeric(scale_x)
        sy <- as.numeric(scale_y)
      }
      validate(need(is.finite(sx) && is.finite(sy) && sx > 0 && sy > 0, "Polygon scale must be positive."))

      cx <- (nx / 2) + translate_x
      cy <- (ny / 2) - translate_y
      theta <- rotate_deg * pi / 180
      ct <- cos(theta)
      st <- sin(theta)

      in_frame <- if (has_source_dims) {
        bb_xmin >= -1 && bb_ymin >= -1 && bb_xmax <= (frame_w + 1) && bb_ymax <= (frame_h + 1)
      } else {
        FALSE
      }
      # Keep polygon and image on the same scale path:
      # never auto-normalize/stretch polygon bbox to frame.
      use_normalized_frame <- FALSE
      # If coordinates are offset but otherwise image-space sized, shift into frame.
      use_shifted_origin <- has_source_dims && !in_frame
      # Keep polygon and image in the same underlying source-image frame.
      # Do not implicitly switch polygon frame size based on bbox dimensions.
      use_bbox_frame <- FALSE
      geom_frame_w <- frame_w
      geom_frame_h <- frame_h
      w <- geom_frame_w * sx
      h <- geom_frame_h * sy

      rebuild_one <- function(one_row) {
        cc <- sf::st_coordinates(one_row)
        validate(need(nrow(cc) > 0, "Polygon has empty coordinates after read."))
        cc_names <- colnames(cc)
        x_col <- if ("X" %in% cc_names) {
          "X"
        } else if ("x" %in% cc_names) {
          "x"
        } else {
          NA_character_
        }
        y_col <- if ("Y" %in% cc_names) {
          "Y"
        } else if ("y" %in% cc_names) {
          "y"
        } else {
          NA_character_
        }
        if (!is.finite(match(x_col, cc_names)) || !is.finite(match(y_col, cc_names))) {
          xy_cols <- cc_names[cc_names %in% c("X", "Y", "x", "y")]
          if (length(xy_cols) >= 2) {
            x_col <- xy_cols[1]
            y_col <- xy_cols[2]
          }
        }
        validate(need(is.character(x_col) && nzchar(x_col) && x_col %in% cc_names, "Polygon coordinates missing X column after transformation."))
        validate(need(is.character(y_col) && nzchar(y_col) && y_col %in% cc_names, "Polygon coordinates missing Y column after transformation."))

        x_src <- cc[, x_col]
        y_src <- cc[, y_col]
        if (isTRUE(swap_xy)) {
          tmp <- x_src
          x_src <- y_src
          y_src <- tmp
        }

        if (use_normalized_frame) {
          x_base <- (x_src - bb_xmin) / bb_w * frame_w
          y_base <- (y_src - bb_ymin) / bb_h * frame_h
        } else if (use_shifted_origin && !use_bbox_frame) {
          x_base <- (x_src - bb_xmin) / bb_w * frame_w
          y_base <- (y_src - bb_ymin) / bb_h * frame_h
        } else if (use_shifted_origin) {
          x_base <- x_src - bb_xmin
          y_base <- y_src - bb_ymin
        } else if (has_source_dims) {
          x_base <- x_src
          y_base <- y_src
        } else {
          x_base <- x_src - bb_xmin
          y_base <- y_src - bb_ymin
        }

        x0 <- x_base * sx
        if (isTRUE(flip_y)) {
          y0 <- y_base * sy
        } else {
          # Convert top-left image-style Y coordinates into plot Y coordinates.
          y0 <- (geom_frame_h - y_base) * sy
        }

        xr <- x0 - w / 2
        yr <- y0 - h / 2
        x1 <- cx + (ct * xr - st * yr)
        y1 <- cy + (st * xr + ct * yr)

        cc_new <- cbind(x1, y1, cc[, setdiff(colnames(cc), c("X", "Y")), drop = FALSE])

        if ("L3" %in% colnames(cc_new)) {
          polys <- lapply(unique(cc_new[, "L3"]), function(l3) {
            sub <- cc_new[cc_new[, "L3"] == l3, , drop = FALSE]
            rings <- lapply(unique(sub[, "L2"]), function(l2) {
              ring <- as.matrix(sub[sub[, "L2"] == l2, c("x1", "y1"), drop = FALSE])
              if (nrow(ring) > 2 && !all(ring[1, ] == ring[nrow(ring), ])) ring <- rbind(ring, ring[1, ])
              ring
            })
            sf::st_polygon(rings)
          })
          if (length(polys) == 1) polys[[1]] else sf::st_multipolygon(polys)
        } else if ("L2" %in% colnames(cc_new)) {
          rings <- lapply(unique(cc_new[, "L2"]), function(l2) {
            ring <- as.matrix(cc_new[cc_new[, "L2"] == l2, c("x1", "y1"), drop = FALSE])
            if (nrow(ring) > 2 && !all(ring[1, ] == ring[nrow(ring), ])) ring <- rbind(ring, ring[1, ])
            ring
          })
          sf::st_polygon(rings)
        } else {
          ring <- as.matrix(cc_new[, c("x1", "y1"), drop = FALSE])
          if (nrow(ring) > 2 && !all(ring[1, ] == ring[nrow(ring), ])) ring <- rbind(ring, ring[1, ])
          sf::st_polygon(list(ring))
        }
      }

      geoms <- lapply(seq_len(nrow(poly)), function(i) {
        try(rebuild_one(poly[i, ]), silent = TRUE)
      })
      good <- !vapply(geoms, inherits, logical(1), "try-error")
      poly <- poly[good, , drop = FALSE]
      geoms <- geoms[good]
      validate(need(length(geoms) > 0, "Could not transform polygon geometries."))

      # Geometry is transformed into MSI pixel space; keep CRS unset.
      sfc <- sf::st_sfc(geoms)
      poly <- sf::st_set_geometry(poly, sfc)
      poly <- sf::st_make_valid(poly)
      gtypes <- unique(as.character(sf::st_geometry_type(poly, by_geometry = TRUE)))
      if (any(!gtypes %in% c("POLYGON", "MULTIPOLYGON"))) {
        # Keep only polygonal geometries; st_make_valid can emit mixed collections.
        poly <- sf::st_collection_extract(poly, "POLYGON", warn = FALSE)
      }
      validate(need(nrow(poly) > 0, "No valid polygon geometries remain after transformation."))
      poly
    }

    msi_signal_matrix <- function(msi_obj) {
      validate(need(!is.null(msi_obj$raster), "MSI raster is unavailable for optimization."))
      ny <- as.integer(msi_obj$ny)
      nx <- as.integer(msi_obj$nx)
      validate(need(is.finite(nx) && nx > 0 && is.finite(ny) && ny > 0, "Invalid MSI canvas for optimization."))

      # Prefer palette-independent optimization signal when provided
      # (used for pData mode so fit is not color-map dependent).
      if (!is.null(msi_obj$opt_signal) && is.matrix(msi_obj$opt_signal)) {
        sig <- msi_obj$opt_signal
        if (nrow(sig) == ny && ncol(sig) == nx) {
          sig <- suppressWarnings(matrix(as.numeric(sig), nrow = ny, ncol = nx))
          sig[!is.finite(sig)] <- 0
          q <- suppressWarnings(stats::quantile(sig[sig > 0], probs = 0.75, na.rm = TRUE, names = FALSE, type = 8))
          if (is.finite(q) && q > 0) sig <- pmax(sig - q, 0)
          mx <- suppressWarnings(max(sig, na.rm = TRUE))
          if (is.finite(mx) && mx > 0) sig <- sig / mx
          return(sig)
        }
      }

      hex <- matrix(as.character(msi_obj$raster), nrow = ny, ncol = nx)
      rr <- suppressWarnings(strtoi(substr(hex, 2, 3), base = 16L))
      gg <- suppressWarnings(strtoi(substr(hex, 4, 5), base = 16L))
      bb <- suppressWarnings(strtoi(substr(hex, 6, 7), base = 16L))
      aa <- rep(255L, length(rr))
      has_alpha <- !is.na(hex) & nchar(hex) >= 9
      aa[has_alpha] <- suppressWarnings(strtoi(substr(hex[has_alpha], 8, 9), base = 16L))

      rr[!is.finite(rr)] <- 0
      gg[!is.finite(gg)] <- 0
      bb[!is.finite(bb)] <- 0
      aa[!is.finite(aa)] <- 255

      # For RGB overlays, use channel-order agnostic signal so optimization
      # does not depend on whether a feature is in R, G, or B.
      if (identical(msi_obj$mode, "rgb")) {
        chroma <- (abs(rr - gg) + abs(rr - bb) + abs(gg - bb)) / 3
        sig <- pmax(rr, gg, bb) + 0.25 * chroma
      } else {
        sig <- 0.2126 * rr + 0.7152 * gg + 0.0722 * bb
      }
      sig[aa <= 0] <- 0

      ok <- is.finite(sig) & sig > 0
      if (any(ok)) {
        q <- suppressWarnings(stats::quantile(sig[ok], probs = 0.75, na.rm = TRUE, names = FALSE, type = 8))
        if (is.finite(q) && q > 0) sig <- pmax(sig - q, 0)
      }
      mx <- suppressWarnings(max(sig, na.rm = TRUE))
      if (is.finite(mx) && mx > 0) sig <- sig / mx
      matrix(sig, nrow = ny, ncol = nx)
    }

    polygon_mask_matrix <- function(msi_obj, tx, ty) {
      req(polygon_data())
      poly <- polygon_data()
      poly$map_label <- get_polygon_labels(poly, input$polygon_label_field)
      axis_mode <- resolve_polygon_axis_mode(poly, input$polygon_axis_mode, get_histology_image_optional())
      src_dims <- get_overlay_source_dims()
      poly_t <- transform_polygon_sf(
        poly_sf = poly,
        nx = msi_obj$nx,
        ny = msi_obj$ny,
        scale_x = input$scale_x,
        scale_y = input$scale_y,
        translate_x = tx,
        translate_y = ty,
        rotate_deg = input$rotate_deg,
        flip_y = isTRUE(input$flip_histology_y),
        swap_xy = identical(axis_mode, "yx"),
        scale_mode = input$overlay_scale_mode,
        source_width = if (!is.null(src_dims)) src_dims$width else NA_real_,
        source_height = if (!is.null(src_dims)) src_dims$height else NA_real_
      )

      grid <- expand.grid(
        x = seq_len(as.integer(msi_obj$nx)),
        y = seq_len(as.integer(msi_obj$ny)),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
      poly_crs <- normalize_crs(try(sf::st_crs(poly_t), silent = TRUE))
      if (is.null(poly_crs)) {
        pts_sf <- sf::st_as_sf(grid, coords = c("x", "y"), remove = FALSE)
      } else {
        pts_sf <- sf::st_as_sf(grid, coords = c("x", "y"), crs = poly_crs, remove = FALSE)
      }
      hit <- sf::st_intersects(pts_sf, poly_t)
      hit_mask <- lengths(hit) > 0
      row_idx <- as.integer(msi_obj$ny - grid$y + 1L)

      mask <- matrix(FALSE, nrow = as.integer(msi_obj$ny), ncol = as.integer(msi_obj$nx))
      mask[cbind(row_idx, as.integer(grid$x))] <- hit_mask
      mask
    }

    polygon_group_label_code_matrix <- function(
      msi_obj,
      tx,
      ty,
      label_field = "clustered_polygon_class",
      exclude_labels = c("outside_polygon", "outside", "unassigned", "unclustered_polygon")
    ) {
      req(polygon_data())
      poly <- polygon_data()
      poly$map_label <- get_polygon_labels(poly, label_field)
      axis_mode <- resolve_polygon_axis_mode(poly, input$polygon_axis_mode, get_histology_image_optional())
      src_dims <- get_overlay_source_dims()
      poly_t <- transform_polygon_sf(
        poly_sf = poly,
        nx = msi_obj$nx,
        ny = msi_obj$ny,
        scale_x = input$scale_x,
        scale_y = input$scale_y,
        translate_x = tx,
        translate_y = ty,
        rotate_deg = input$rotate_deg,
        flip_y = isTRUE(input$flip_histology_y),
        swap_xy = identical(axis_mode, "yx"),
        scale_mode = input$overlay_scale_mode,
        source_width = if (!is.null(src_dims)) src_dims$width else NA_real_,
        source_height = if (!is.null(src_dims)) src_dims$height else NA_real_
      )

      grid <- expand.grid(
        x = seq_len(as.integer(msi_obj$nx)),
        y = seq_len(as.integer(msi_obj$ny)),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
      poly_crs <- normalize_crs(try(sf::st_crs(poly_t), silent = TRUE))
      if (is.null(poly_crs)) {
        pts_sf <- sf::st_as_sf(grid, coords = c("x", "y"), remove = FALSE)
      } else {
        pts_sf <- sf::st_as_sf(grid, coords = c("x", "y"), crs = poly_crs, remove = FALSE)
      }
      hit <- sf::st_intersects(pts_sf, poly_t)

      raw_labels <- as.character(poly_t$map_label)
      is_cell_poly <- is_cell_like_polygon_label(raw_labels)
      poly_area <- suppressWarnings(as.numeric(sf::st_area(poly_t)))
      poly_area[!is.finite(poly_area)] <- 0
      canvas_area <- as.numeric(msi_obj$nx) * as.numeric(msi_obj$ny)
      huge_poly <- is.finite(poly_area) & poly_area > (0.5 * canvas_area)
      if (any(huge_poly) && any(is_cell_poly & !huge_poly)) {
        is_cell_poly <- is_cell_poly & !huge_poly
      }
      if (!any(is_cell_poly)) is_cell_poly <- rep(TRUE, length(raw_labels))

      label_grid <- rep(NA_character_, nrow(grid))
      has_hit <- lengths(hit) > 0L
      if (any(has_hit)) {
        label_grid[has_hit] <- vapply(hit[has_hit], function(ix) {
          ix_use <- ix[is_cell_poly[ix]]
          if (length(ix_use) == 0L) return(NA_character_)
          labs <- as.character(raw_labels[ix_use])
          labs <- labs[!is.na(labs) & nzchar(trimws(labs))]
          if (length(labs) == 0L) return(NA_character_)
          labs[1]
        }, character(1))
      }

      ex_norm <- tolower(trimws(as.character(exclude_labels)))
      ex_norm <- ex_norm[is.finite(nchar(ex_norm))]
      lab_norm <- tolower(trimws(label_grid))
      lab_norm[is.na(lab_norm)] <- ""
      valid_lab <- nzchar(lab_norm) & !(lab_norm %in% ex_norm)
      group_levels <- sort(unique(label_grid[valid_lab]))

      code_vec <- integer(nrow(grid))
      if (length(group_levels) > 0L) {
        m <- match(label_grid, group_levels)
        keep <- valid_lab & !is.na(m)
        code_vec[keep] <- as.integer(m[keep])
      }

      row_idx <- as.integer(msi_obj$ny - grid$y + 1L)
      code_mat <- matrix(0L, nrow = as.integer(msi_obj$ny), ncol = as.integer(msi_obj$nx))
      code_mat[cbind(row_idx, as.integer(grid$x))] <- code_vec

      list(
        code_mat = code_mat,
        group_levels = group_levels,
        n_groups = length(group_levels),
        n_assigned = as.integer(sum(code_vec > 0L, na.rm = TRUE))
      )
    }

    pdata_group_label_code_matrix <- function(
      msi_obj,
      field,
      exclude_labels = c("outside_polygon", "outside", "unassigned", "unclustered_polygon")
    ) {
      obj <- msi_for_pdata()
      pd <- as.data.frame(Cardinal::pData(obj))
      if (!field %in% names(pd)) {
        stop(sprintf("pData field '%s' not found.", field))
      }
      if (nrow(pd) != length(msi_obj$x_norm)) {
        stop(sprintf("pData length mismatch for Stat Fit group labels (%d rows vs %d pixels).", nrow(pd), length(msi_obj$x_norm)))
      }

      labels <- as.character(pd[[field]])
      labels[is.na(labels)] <- ""
      labels <- trimws(labels)
      if (any(grepl(";", labels, fixed = TRUE), na.rm = TRUE)) {
        labels <- trimws(sub(";.*$", "", labels))
      }

      ex_norm <- tolower(trimws(as.character(exclude_labels)))
      ex_norm <- ex_norm[!is.na(ex_norm) & nzchar(ex_norm)]
      lab_norm <- tolower(labels)
      lab_norm[is.na(lab_norm)] <- ""
      valid_lab <- nzchar(lab_norm) & !(lab_norm %in% ex_norm)
      group_levels <- sort(unique(labels[valid_lab]))

      code_vec <- integer(length(labels))
      if (length(group_levels) > 0L) {
        m <- match(labels, group_levels)
        keep <- valid_lab & !is.na(m)
        code_vec[keep] <- as.integer(m[keep])
      }

      code_mat <- matrix(0L, nrow = as.integer(msi_obj$ny), ncol = as.integer(msi_obj$nx))
      code_mat[cbind(as.integer(msi_obj$row_idx), as.integer(msi_obj$x_norm))] <- code_vec

      list(
        code_mat = code_mat,
        group_levels = group_levels,
        n_groups = length(group_levels),
        n_assigned = as.integer(sum(code_vec > 0L, na.rm = TRUE))
      )
    }

    stat_fit_group_label_code_matrix <- function(
      msi_obj,
      tx,
      ty,
      group_field_sel = NULL,
      exclude_labels = c("outside_polygon", "outside", "unassigned", "unclustered_polygon")
    ) {
      spec <- parse_stat_fit_group_field(group_field_sel)
      if (identical(spec$source, "pdata")) {
        out <- pdata_group_label_code_matrix(msi_obj = msi_obj, field = spec$field, exclude_labels = exclude_labels)
        out$group_source <- "pdata"
        out$group_field <- spec$field
        return(out)
      }
      if (identical(spec$source, "polygon")) {
        out <- polygon_group_label_code_matrix(msi_obj = msi_obj, tx = tx, ty = ty, label_field = spec$field, exclude_labels = exclude_labels)
        out$group_source <- "polygon"
        out$group_field <- spec$field
        return(out)
      }

      # Auto fallback: prefer pData if present, otherwise polygon field.
      obj_try <- try(msi_for_pdata(), silent = TRUE)
      if (!inherits(obj_try, "try-error") && !is.null(obj_try)) {
        pd_cols <- try(colnames(as.data.frame(Cardinal::pData(obj_try))), silent = TRUE)
        if (!inherits(pd_cols, "try-error") && spec$field %in% pd_cols) {
          out <- pdata_group_label_code_matrix(msi_obj = msi_obj, field = spec$field, exclude_labels = exclude_labels)
          out$group_source <- "pdata"
          out$group_field <- spec$field
          return(out)
        }
      }
      out <- polygon_group_label_code_matrix(msi_obj = msi_obj, tx = tx, ty = ty, label_field = spec$field, exclude_labels = exclude_labels)
      out$group_source <- "polygon"
      out$group_field <- spec$field
      out
    }

    edge_fit_signal_from_pdata_field <- function(msi_obj, field) {
      obj <- msi_for_pdata()
      pd <- as.data.frame(Cardinal::pData(obj))
      if (!field %in% names(pd)) {
        stop(sprintf("pData field '%s' not found for Edge Fit optimization.", field))
      }

      ny <- as.integer(msi_obj$ny)
      nx <- as.integer(msi_obj$nx)
      row_idx <- as.integer(msi_obj$row_idx)
      x_norm <- as.integer(msi_obj$x_norm)
      vals <- pd[[field]]

      # Treat low-cardinality integer-like numerics as categorical labels.
      vals_num <- suppressWarnings(as.numeric(vals))
      is_num <- is.numeric(vals) || is.integer(vals)
      unique_n <- length(unique(vals[!is.na(vals)]))
      integer_like <- is_num && all(is.na(vals_num) | abs(vals_num - round(vals_num)) < 1e-8)
      as_categorical <- (!is_num) || (integer_like && unique_n <= 128L)

      if (!as_categorical) {
        mat <- matrix(NA_real_, nrow = ny, ncol = nx)
        mat[cbind(row_idx, x_norm)] <- vals_num
        if (isTRUE(input$gaussian_smooth)) {
          smooth_sigma <- suppressWarnings(as.numeric(input$gaussian_sigma))
          if (!is.finite(smooth_sigma) || smooth_sigma <= 0) smooth_sigma <- 1
          mat <- gaussian_smooth_matrix(mat, smooth_sigma)
        }
        mat[!is.finite(mat)] <- 0
        return(list(signal = mat, field = field, type = "numeric"))
      }

      labs <- normalize_labels(vals)
      lab_mat <- matrix(NA_character_, nrow = ny, ncol = nx)
      lab_mat[cbind(row_idx, x_norm)] <- labs
      valid <- !is.na(lab_mat)
      right <- cbind(lab_mat[, -1, drop = FALSE], NA_character_)
      left <- cbind(NA_character_, lab_mat[, -ncol(lab_mat), drop = FALSE])
      down <- rbind(lab_mat[-1, , drop = FALSE], rep(NA_character_, ncol(lab_mat)))
      up <- rbind(rep(NA_character_, ncol(lab_mat)), lab_mat[-nrow(lab_mat), , drop = FALSE])

      b <- matrix(0, nrow = ny, ncol = nx)
      b <- b + ((lab_mat != right) & valid & !is.na(right))
      b <- b + ((lab_mat != left) & valid & !is.na(left))
      b <- b + ((lab_mat != down) & valid & !is.na(down))
      b <- b + ((lab_mat != up) & valid & !is.na(up))
      b[!is.finite(b)] <- 0
      b[!valid] <- 0
      if (isTRUE(input$gaussian_smooth)) {
        smooth_sigma <- suppressWarnings(as.numeric(input$gaussian_sigma))
        if (!is.finite(smooth_sigma) || smooth_sigma <= 0) smooth_sigma <- 1
        b <- gaussian_smooth_matrix(b, smooth_sigma)
      }
      b[!is.finite(b)] <- 0
      list(signal = b, field = field, type = "categorical")
    }

    edge_fit_signal_matrix <- function(msi_obj) {
      source_mode <- tolower(trimws(as.character(input$edge_fit_signal_source)[1]))
      if (!source_mode %in% c("current", "pdata")) source_mode <- "current"
      if (!identical(source_mode, "pdata")) {
        return(list(
          signal = msi_signal_matrix(msi_obj),
          source = "current",
          field = if (!is.null(msi_obj$pdata_field)) as.character(msi_obj$pdata_field) else NA_character_,
          type = if (!is.null(msi_obj$mode)) as.character(msi_obj$mode) else NA_character_
        ))
      }
      field <- as.character(input$edge_fit_pdata_field)[1]
      if (is.na(field) || !nzchar(field)) {
        stop("Select a pData field for Edge Fit optimization.")
      }
      out <- edge_fit_signal_from_pdata_field(msi_obj = msi_obj, field = field)
      list(
        signal = out$signal,
        source = "pdata",
        field = out$field,
        type = out$type
      )
    }

    shift_bool <- function(mat, dr = 0L, dc = 0L) {
      ny <- nrow(mat)
      nx <- ncol(mat)
      out <- matrix(FALSE, nrow = ny, ncol = nx)
      src_r <- seq_len(ny)
      src_c <- seq_len(nx)
      tgt_r <- src_r + as.integer(dr)
      tgt_c <- src_c + as.integer(dc)
      keep_r <- which(tgt_r >= 1L & tgt_r <= ny)
      keep_c <- which(tgt_c >= 1L & tgt_c <= nx)
      if (length(keep_r) == 0L || length(keep_c) == 0L) return(out)
      out[tgt_r[keep_r], tgt_c[keep_c]] <- mat[src_r[keep_r], src_c[keep_c], drop = FALSE]
      out
    }

    shift_int <- function(mat, dr = 0L, dc = 0L, fill = 0L) {
      ny <- nrow(mat)
      nx <- ncol(mat)
      out <- matrix(as.integer(fill), nrow = ny, ncol = nx)
      src_r <- seq_len(ny)
      src_c <- seq_len(nx)
      tgt_r <- src_r + as.integer(dr)
      tgt_c <- src_c + as.integer(dc)
      keep_r <- which(tgt_r >= 1L & tgt_r <= ny)
      keep_c <- which(tgt_c >= 1L & tgt_c <= nx)
      if (length(keep_r) == 0L || length(keep_c) == 0L) return(out)
      out[tgt_r[keep_r], tgt_c[keep_c]] <- mat[src_r[keep_r], src_c[keep_c], drop = FALSE]
      out
    }

    dilate8 <- function(mat) {
      mat <- mat & TRUE
      out <- mat
      for (dr in -1L:1L) {
        for (dc in -1L:1L) {
          if (dr == 0L && dc == 0L) next
          out <- out | shift_bool(mat, dr = dr, dc = dc)
        }
      }
      out
    }

    erode8 <- function(mat) {
      mat <- mat & TRUE
      out <- mat
      for (dr in -1L:1L) {
        for (dc in -1L:1L) {
          if (dr == 0L && dc == 0L) next
          out <- out & shift_bool(mat, dr = dr, dc = dc)
        }
      }
      out
    }

    idx_to_coords <- function(idx, nr) {
      idx <- as.integer(idx)
      if (length(idx) == 0L) {
        return(list(r = integer(0), c = integer(0), n_total = 0L))
      }
      list(
        r = ((idx - 1L) %% nr) + 1L,
        c = ((idx - 1L) %/% nr) + 1L,
        n_total = length(idx)
      )
    }

    build_mask_bands <- function(mask, band_px = 4L) {
      mask <- (mask & TRUE)
      ny <- nrow(mask)
      nx <- ncol(mask)
      band_px <- suppressWarnings(as.integer(band_px))
      if (!is.finite(band_px) || band_px < 1L) band_px <- 4L
      band_px <- as.integer(min(50L, max(1L, band_px)))

      inside_rings <- vector("list", band_px)
      outside_rings <- vector("list", band_px)

      inside_work <- mask
      grown <- mask
      for (k in seq_len(band_px)) {
        eroded <- erode8(inside_work)
        in_ring <- inside_work & !eroded
        inside_work <- eroded

        dilated <- dilate8(grown)
        out_ring <- dilated & !grown
        grown <- dilated

        inside_rings[[k]] <- idx_to_coords(which(in_ring), ny)
        outside_rings[[k]] <- idx_to_coords(which(out_ring), ny)
      }

      list(
        ny = ny,
        nx = nx,
        band_px = band_px,
        inside_all = idx_to_coords(which(mask), ny),
        inside_core = idx_to_coords(which(inside_work), ny),
        boundary = idx_to_coords(which(mask & !erode8(mask)), ny),
        inside_rings = inside_rings,
        outside_rings = outside_rings
      )
    }

    sample_shifted_stats <- function(coord_obj, signal, dx = 0L, dy = 0L, min_n = 5L) {
      if (is.null(coord_obj) || length(coord_obj$r) == 0L) {
        return(list(ok = FALSE, mean = NA_real_, q90 = NA_real_, n = 0L, cov = 0))
      }
      ny <- nrow(signal)
      nx <- ncol(signal)
      rr <- coord_obj$r + as.integer(dy)
      cc <- coord_obj$c + as.integer(dx)
      keep <- rr >= 1L & rr <= ny & cc >= 1L & cc <= nx
      n_keep <- sum(keep)
      cov <- n_keep / max(1L, as.integer(coord_obj$n_total))
      if (!is.finite(n_keep) || n_keep < min_n) {
        return(list(ok = FALSE, mean = NA_real_, q90 = NA_real_, n = as.integer(n_keep), cov = cov))
      }
      vals <- as.numeric(signal[cbind(rr[keep], cc[keep])])
      vals <- vals[is.finite(vals)]
      n_vals <- length(vals)
      if (n_vals < min_n) {
        return(list(ok = FALSE, mean = NA_real_, q90 = NA_real_, n = as.integer(n_vals), cov = cov))
      }
      list(
        ok = TRUE,
        mean = suppressWarnings(mean(vals, na.rm = TRUE)),
        q90 = suppressWarnings(stats::quantile(vals, probs = 0.9, na.rm = TRUE, names = FALSE, type = 8)),
        n = as.integer(n_vals),
        cov = cov
      )
    }

    score_mask_shift <- function(mask_cache, signal, dx = 0L, dy = 0L) {
      if (is.matrix(mask_cache)) {
        mask_cache <- build_mask_bands(mask_cache, band_px = 4L)
      }
      if (is.null(mask_cache) || is.null(mask_cache$inside_all)) return(-Inf)
      dx <- as.integer(round(dx))
      dy <- as.integer(round(dy))

      in_all <- sample_shifted_stats(mask_cache$inside_all, signal, dx = dx, dy = dy, min_n = 10L)
      in_relaxed <- if (isTRUE(in_all$ok)) {
        in_all
      } else {
        sample_shifted_stats(mask_cache$inside_all, signal, dx = dx, dy = dy, min_n = 1L)
      }
      if (!isTRUE(in_relaxed$ok)) return(-Inf)
      if (!is.finite(in_relaxed$cov) || in_relaxed$cov <= 0) return(-Inf)
      coverage_all <- pmax(0, pmin(1, in_relaxed$cov))
      if (coverage_all < 0.01) return(-Inf)

      band_n <- min(
        length(mask_cache$inside_rings),
        length(mask_cache$outside_rings)
      )
      if (!is.finite(band_n) || band_n < 1) {
        band_n <- 0L
      }

      # Boundary-gradient score: compare inside and outside signal across ring depth.
      w <- if (band_n > 0L) 1 / seq_len(band_n) else numeric(0)
      grad_num <- 0
      grad_den <- 0
      cov_num <- 0
      cov_den <- 0
      out_sum <- 0
      out_n <- 0
      if (band_n > 0L) {
        for (k in seq_len(band_n)) {
          in_k <- sample_shifted_stats(mask_cache$inside_rings[[k]], signal, dx = dx, dy = dy, min_n = 3L)
          out_k <- sample_shifted_stats(mask_cache$outside_rings[[k]], signal, dx = dx, dy = dy, min_n = 3L)
          if (isTRUE(in_k$ok) && isTRUE(out_k$ok)) {
            wk <- w[k]
            grad_num <- grad_num + wk * (in_k$mean - out_k$mean)
            grad_den <- grad_den + wk
            cov_num <- cov_num + wk * min(in_k$cov, out_k$cov)
            cov_den <- cov_den + wk
            out_sum <- out_sum + out_k$mean * out_k$n
            out_n <- out_n + out_k$n
          }
        }
      }

      out_mean <- if (out_n > 0) out_sum / out_n else suppressWarnings(mean(as.numeric(signal), na.rm = TRUE))
      if (!is.finite(out_mean)) out_mean <- 0

      has_grad <- (grad_den > 0 && cov_den > 0)
      grad_contrast <- if (has_grad) grad_num / grad_den else 0
      ring_cov <- if (has_grad) cov_num / cov_den else coverage_all * 0.5
      if (!is.finite(ring_cov)) ring_cov <- 0

      global_contrast <- in_relaxed$mean - out_mean
      q_boost <- in_relaxed$q90 - out_mean
      if (!is.finite(q_boost)) q_boost <- 0
      core <- sample_shifted_stats(mask_cache$inside_core, signal, dx = dx, dy = dy, min_n = 5L)
      core_term <- if (isTRUE(core$ok)) (core$mean - in_relaxed$mean) else 0

      coverage <- pmax(0.01, pmin(1, if (has_grad) pmin(coverage_all, pmax(0, ring_cov)) else coverage_all))

      base <- 0.60 * global_contrast + 0.25 * grad_contrast + 0.10 * q_boost + 0.05 * core_term
      score <- base - 1.20 * (1 - coverage)^2
      if (!is.finite(score)) {
        score <- global_contrast - 1.50 * (1 - coverage_all)^2
      }
      if (!is.finite(score)) return(-Inf)
      score
    }

    build_edge_distance_cache <- function(signal, q = 0.85) {
      ny <- nrow(signal)
      nx <- ncol(signal)
      if (!is.finite(ny) || !is.finite(nx) || ny < 2 || nx < 2) {
        return(NULL)
      }

      right <- cbind(signal[, -1, drop = FALSE], signal[, nx, drop = FALSE])
      down <- rbind(signal[-1, , drop = FALSE], signal[ny, , drop = FALSE])
      grad <- abs(right - signal) + abs(down - signal)
      grad[!is.finite(grad)] <- 0

      g_ok <- grad[is.finite(grad) & grad > 0]
      if (length(g_ok) < 10L) return(NULL)
      thr <- suppressWarnings(stats::quantile(g_ok, probs = q, na.rm = TRUE, names = FALSE, type = 8))
      if (!is.finite(thr)) {
        thr <- suppressWarnings(mean(g_ok, na.rm = TRUE))
      }
      if (!is.finite(thr)) return(NULL)

      edge <- grad >= thr
      n_edge <- suppressWarnings(sum(edge, na.rm = TRUE))
      if (!is.finite(n_edge) || n_edge < 5L) return(NULL)

      # Two-pass chamfer distance transform from MSI edges.
      inf <- 1e6
      d <- matrix(inf, nrow = ny, ncol = nx)
      d[edge] <- 0
      wdiag <- sqrt(2)

      for (r in seq_len(ny)) {
        for (c in seq_len(nx)) {
          v <- d[r, c]
          if (r > 1L) v <- min(v, d[r - 1L, c] + 1)
          if (c > 1L) v <- min(v, d[r, c - 1L] + 1)
          if (r > 1L && c > 1L) v <- min(v, d[r - 1L, c - 1L] + wdiag)
          if (r > 1L && c < nx) v <- min(v, d[r - 1L, c + 1L] + wdiag)
          d[r, c] <- v
        }
      }
      for (r in seq.int(ny, 1L, by = -1L)) {
        for (c in seq.int(nx, 1L, by = -1L)) {
          v <- d[r, c]
          if (r < ny) v <- min(v, d[r + 1L, c] + 1)
          if (c < nx) v <- min(v, d[r, c + 1L] + 1)
          if (r < ny && c < nx) v <- min(v, d[r + 1L, c + 1L] + wdiag)
          if (r < ny && c > 1L) v <- min(v, d[r + 1L, c - 1L] + wdiag)
          d[r, c] <- v
        }
      }
      d <- pmin(d, 50)
      list(dist = d, n_edge = as.integer(n_edge), threshold = thr)
    }

    score_edge_shift <- function(mask_cache, edge_cache, dx = 0L, dy = 0L) {
      if (is.null(mask_cache) || is.null(mask_cache$boundary) || is.null(edge_cache) || is.null(edge_cache$dist)) {
        return(-Inf)
      }
      dx <- as.integer(round(dx))
      dy <- as.integer(round(dy))

      b <- sample_shifted_stats(mask_cache$boundary, edge_cache$dist, dx = dx, dy = dy, min_n = 8L)
      if (!isTRUE(b$ok)) {
        b <- sample_shifted_stats(mask_cache$boundary, edge_cache$dist, dx = dx, dy = dy, min_n = 1L)
      }
      if (!isTRUE(b$ok)) return(-Inf)
      if (!is.finite(b$mean)) return(-Inf)

      cov <- pmax(0.01, pmin(1, if (is.finite(b$cov)) b$cov else 0))
      # Lower distance to MSI edges is better; include coverage penalty.
      sc <- (-b$mean) - 1.10 * (1 - cov)^2
      if (!is.finite(sc)) return(-Inf)
      sc
    }

    get_intensity_matrix_fast <- function(msi_obj) {
      X <- try(as.matrix(Cardinal::spectra(msi_obj)), silent = TRUE)
      if (inherits(X, "try-error") || is.null(X)) {
        X <- NULL
      }
      if (is.null(X) && requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        X <- try(as.matrix(SummarizedExperiment::assay(msi_obj)), silent = TRUE)
        if (inherits(X, "try-error")) X <- NULL
      }
      if (is.null(X)) {
        stop("Could not extract MSI intensity matrix for statistical fit.")
      }
      if (nrow(X) == ncol(msi_obj) && ncol(X) == nrow(msi_obj)) {
        X <- t(X)
      }
      if (nrow(X) != nrow(msi_obj) || ncol(X) != ncol(msi_obj)) {
        stop(sprintf(
          "Intensity matrix dimensions mismatch: got %d x %d, expected %d x %d (features x pixels).",
          nrow(X), ncol(X), nrow(msi_obj), ncol(msi_obj)
        ))
      }
      storage.mode(X) <- "double"
      X
    }

    select_top_variable_features <- function(X, max_features = 200L) {
      n_feat <- nrow(X)
      if (!is.finite(n_feat) || n_feat < 1) return(integer(0))
      k <- suppressWarnings(as.integer(max_features))
      if (!is.finite(k) || k < 1L || k >= n_feat) return(seq_len(n_feat))

      Xz <- X
      Xz[!is.finite(Xz)] <- 0
      nobs <- rowSums(is.finite(X))
      s <- rowSums(Xz)
      ss <- rowSums(Xz * Xz)
      vv <- (ss - (s * s) / pmax(nobs, 1L)) / pmax(nobs - 1L, 1L)
      vv[!is.finite(vv)] <- 0
      ord <- order(vv, decreasing = TRUE)
      ord[seq_len(min(k, length(ord)))]
    }

    local_outside_mask <- function(mask_inside, buffer_px = 3L) {
      b <- suppressWarnings(as.integer(buffer_px))
      if (!is.finite(b) || b < 1L) b <- 3L
      b <- as.integer(min(25L, max(1L, b)))
      grown <- mask_inside
      for (i in seq_len(b)) {
        grown <- dilate8(grown)
      }
      grown & !mask_inside
    }

    bbox_outside_mask <- function(mask_inside, pad_px = 25L) {
      b <- suppressWarnings(as.integer(pad_px))
      if (!is.finite(b) || b < 0L) b <- 25L
      b <- as.integer(min(500L, max(0L, b)))

      idx <- which(mask_inside)
      if (length(idx) == 0L) {
        return(matrix(FALSE, nrow = nrow(mask_inside), ncol = ncol(mask_inside)))
      }
      rc <- arrayInd(idx, .dim = dim(mask_inside))
      rmin <- max(1L, min(rc[, 1]) - b)
      rmax <- min(nrow(mask_inside), max(rc[, 1]) + b)
      cmin <- max(1L, min(rc[, 2]) - b)
      cmax <- min(ncol(mask_inside), max(rc[, 2]) + b)

      box <- matrix(FALSE, nrow = nrow(mask_inside), ncol = ncol(mask_inside))
      box[rmin:rmax, cmin:cmax] <- TRUE
      box & !mask_inside
    }

    stat_ttest_count_sig <- function(
      X,
      idx_inside,
      idx_outside,
      alpha = 0.1,
      use_adjusted = TRUE,
      use_abs_lfc = TRUE
    ) {
      if (length(idx_inside) < 2L || length(idx_outside) < 2L) {
        return(list(n_sig = 0L, score = -Inf, med_lfc = NA_real_, mean_logp = NA_real_))
      }

      X1 <- X[, idx_inside, drop = FALSE]
      X2 <- X[, idx_outside, drop = FALSE]

      n1 <- rowSums(is.finite(X1))
      n2 <- rowSums(is.finite(X2))

      X1z <- X1
      X2z <- X2
      X1z[!is.finite(X1z)] <- 0
      X2z[!is.finite(X2z)] <- 0

      s1 <- rowSums(X1z)
      s2 <- rowSums(X2z)
      ss1 <- rowSums(X1z * X1z)
      ss2 <- rowSums(X2z * X2z)

      mean1 <- s1 / pmax(n1, 1L)
      mean2 <- s2 / pmax(n2, 1L)
      var1 <- (ss1 - (s1 * s1) / pmax(n1, 1L)) / pmax(n1 - 1L, 1L)
      var2 <- (ss2 - (s2 * s2) / pmax(n2, 1L)) / pmax(n2 - 1L, 1L)

      se2 <- (var1 / pmax(n1, 1L)) + (var2 / pmax(n2, 1L))
      tstat <- (mean1 - mean2) / sqrt(se2)
      df <- (se2 * se2) / (
        (var1 * var1) / (pmax(n1, 1L) * pmax(n1, 1L) * pmax(n1 - 1L, 1L)) +
          (var2 * var2) / (pmax(n2, 1L) * pmax(n2, 1L) * pmax(n2 - 1L, 1L))
      )

      bad <- !is.finite(tstat) | !is.finite(df) | df <= 0 | !is.finite(se2) | se2 <= 0
      pval <- rep(NA_real_, length(tstat))
      pval[!bad] <- 2 * stats::pt(abs(tstat[!bad]), df = df[!bad], lower.tail = FALSE)

      score_p <- if (isTRUE(use_adjusted)) stats::p.adjust(pval, method = "BH") else pval
      score_p[!is.finite(score_p)] <- NA_real_

      lfc <- log2((mean1 + 1e-12) / (mean2 + 1e-12))
      eff <- if (isTRUE(use_abs_lfc)) abs(lfc) else lfc
      eff[!is.finite(eff)] <- NA_real_

      sig <- is.finite(score_p) & (score_p < alpha)
      n_sig <- as.integer(sum(sig, na.rm = TRUE))
      med_lfc <- if (any(sig, na.rm = TRUE)) suppressWarnings(stats::median(eff[sig], na.rm = TRUE)) else 0
      mean_logp <- if (any(sig, na.rm = TRUE)) suppressWarnings(mean(-log10(score_p[sig] + 1e-300), na.rm = TRUE)) else 0
      if (!is.finite(med_lfc)) med_lfc <- 0
      if (!is.finite(mean_logp)) mean_logp <- 0

      abs_t <- abs(tstat)
      abs_t <- abs_t[is.finite(abs_t)]
      if (length(abs_t) > 0L) {
        abs_t <- sort(abs_t, decreasing = TRUE)
        top_k <- max(1L, floor(length(abs_t) * 0.1))
        mean_abs_t_top <- suppressWarnings(mean(abs_t[seq_len(top_k)], na.rm = TRUE))
      } else {
        mean_abs_t_top <- 0
      }
      if (!is.finite(mean_abs_t_top)) mean_abs_t_top <- 0

      # Primary objective tracks "count significant features".
      # Add signal-strength tie-breakers so scores still vary when n_sig is low.
      score <- as.numeric(n_sig) + 0.25 * med_lfc + 0.1 * mean_logp + 0.05 * mean_abs_t_top
      list(
        n_sig = n_sig,
        score = score,
        med_lfc = med_lfc,
        mean_logp = mean_logp,
        mean_abs_t_top = mean_abs_t_top
      )
    }

    stat_group_anova_core <- function(X, idx_pixels, group_index, use_adjusted = TRUE) {
      if (length(idx_pixels) < 3L || length(group_index) != length(idx_pixels)) return(NULL)

      idx_pixels <- suppressWarnings(as.integer(idx_pixels))
      group_index <- suppressWarnings(as.integer(group_index))
      ok <- is.finite(idx_pixels) & idx_pixels >= 1L & idx_pixels <= ncol(X) &
        is.finite(group_index) & group_index >= 1L
      idx_pixels <- idx_pixels[ok]
      group_index <- group_index[ok]
      if (length(idx_pixels) < 3L) return(NULL)

      lev <- sort(unique(group_index))
      if (length(lev) < 2L) return(NULL)
      group_index <- match(group_index, lev)
      K <- length(lev)
      if (length(group_index) <= K) return(NULL)

      Xg <- X[, idx_pixels, drop = FALSE]
      nf <- nrow(Xg)
      n_by <- matrix(0, nrow = nf, ncol = K)
      s_by <- matrix(0, nrow = nf, ncol = K)
      ss_by <- matrix(0, nrow = nf, ncol = K)

      for (j in seq_len(K)) {
        jj <- which(group_index == j)
        if (length(jj) == 0L) next
        Xj <- Xg[, jj, drop = FALSE]
        fin <- is.finite(Xj)
        Xjz <- Xj
        Xjz[!fin] <- 0
        n_by[, j] <- rowSums(fin)
        s_by[, j] <- rowSums(Xjz)
        ss_by[, j] <- rowSums(Xjz * Xjz)
      }

      n_tot <- rowSums(n_by)
      s_tot <- rowSums(s_by)
      ss_tot <- rowSums(ss_by)
      k_feat <- rowSums(n_by > 0L)

      term_between <- rowSums(ifelse(n_by > 0L, (s_by * s_by) / pmax(n_by, 1L), 0))
      ss_between <- term_between - ifelse(n_tot > 0L, (s_tot * s_tot) / pmax(n_tot, 1L), 0)
      ss_within <- ss_tot - term_between
      ss_between[!is.finite(ss_between)] <- NA_real_
      ss_within[!is.finite(ss_within)] <- NA_real_
      ss_between <- pmax(ss_between, 0)
      ss_within <- pmax(ss_within, 0)

      df1 <- k_feat - 1L
      df2 <- n_tot - k_feat
      ms_between <- ss_between / pmax(df1, 1L)
      ms_within <- ss_within / pmax(df2, 1L)
      f_stat <- ms_between / ms_within

      bad <- !is.finite(f_stat) | !is.finite(df1) | !is.finite(df2) |
        df1 <= 0L | df2 <= 0L | !is.finite(ms_within) | ms_within <= 0
      pval <- rep(NA_real_, length(f_stat))
      pval[!bad] <- stats::pf(f_stat[!bad], df1 = df1[!bad], df2 = df2[!bad], lower.tail = FALSE)
      padj <- stats::p.adjust(pval, method = "BH")
      score_p <- if (isTRUE(use_adjusted)) padj else pval
      score_p[!is.finite(score_p)] <- NA_real_

      eta2 <- ss_between / pmax(ss_between + ss_within, 1e-300)
      eta2[!is.finite(eta2)] <- NA_real_

      mean_by <- s_by / pmax(n_by, 1L)
      mean_by[n_by <= 0L] <- NA_real_
      mean_max <- suppressWarnings(apply(mean_by, 1L, function(v) if (all(!is.finite(v))) NA_real_ else max(v, na.rm = TRUE)))
      mean_min <- suppressWarnings(apply(mean_by, 1L, function(v) if (all(!is.finite(v))) NA_real_ else min(v, na.rm = TRUE)))
      mean_range <- mean_max - mean_min
      mean_range[!is.finite(mean_range)] <- NA_real_

      list(
        p_value = as.numeric(pval),
        p_adj = as.numeric(padj),
        score_p = as.numeric(score_p),
        f_stat = as.numeric(f_stat),
        eta2 = as.numeric(eta2),
        mean_range = as.numeric(mean_range),
        n_groups = as.integer(k_feat),
        n_pixels = as.integer(n_tot)
      )
    }

    stat_group_anova_count_sig <- function(
      X,
      idx_pixels,
      group_index,
      alpha = 0.1,
      use_adjusted = TRUE
    ) {
      core <- stat_group_anova_core(
        X = X,
        idx_pixels = idx_pixels,
        group_index = group_index,
        use_adjusted = use_adjusted
      )
      if (is.null(core)) {
        return(list(n_sig = 0L, score = -Inf, med_lfc = NA_real_, mean_logp = NA_real_, mean_abs_t_top = NA_real_))
      }

      score_p <- core$score_p
      sig <- is.finite(score_p) & (score_p < alpha)
      n_sig <- as.integer(sum(sig, na.rm = TRUE))
      med_eta2 <- if (any(sig, na.rm = TRUE)) suppressWarnings(stats::median(core$eta2[sig], na.rm = TRUE)) else 0
      mean_logp <- if (any(sig, na.rm = TRUE)) suppressWarnings(mean(-log10(score_p[sig] + 1e-300), na.rm = TRUE)) else 0
      if (!is.finite(med_eta2)) med_eta2 <- 0
      if (!is.finite(mean_logp)) mean_logp <- 0

      f_ok <- core$f_stat[is.finite(core$f_stat) & core$f_stat > 0]
      if (length(f_ok) > 0L) {
        f_ok <- sort(f_ok, decreasing = TRUE)
        top_k <- max(1L, floor(length(f_ok) * 0.1))
        mean_top_strength <- suppressWarnings(mean(log10(f_ok[seq_len(top_k)] + 1), na.rm = TRUE))
      } else {
        mean_top_strength <- 0
      }
      if (!is.finite(mean_top_strength)) mean_top_strength <- 0

      score <- as.numeric(n_sig) + 0.25 * med_eta2 + 0.1 * mean_logp + 0.05 * mean_top_strength
      list(
        n_sig = n_sig,
        score = score,
        med_lfc = med_eta2,
        mean_logp = mean_logp,
        mean_abs_t_top = mean_top_strength
      )
    }

    stat_group_anova_feature_table <- function(
      X,
      idx_pixels,
      group_index,
      mz_values = NULL,
      use_adjusted = TRUE
    ) {
      core <- stat_group_anova_core(
        X = X,
        idx_pixels = idx_pixels,
        group_index = group_index,
        use_adjusted = use_adjusted
      )
      if (is.null(core)) return(data.frame())

      if (is.null(mz_values) || length(mz_values) != nrow(X)) {
        mz_values <- rep(NA_real_, nrow(X))
      }
      neglogp <- -log10(core$score_p + 1e-300)
      neglogp[!is.finite(neglogp)] <- 0

      out <- data.frame(
        feature_index = seq_len(nrow(X)),
        mz = suppressWarnings(as.numeric(mz_values)),
        n_groups = as.integer(core$n_groups),
        n_pixels = as.integer(core$n_pixels),
        F_stat = as.numeric(core$f_stat),
        eta2 = as.numeric(core$eta2),
        mean_range = as.numeric(core$mean_range),
        p_value = as.numeric(core$p_value),
        p_adj = as.numeric(core$p_adj),
        score_p = as.numeric(core$score_p),
        neglog10_score_p = as.numeric(neglogp),
        stringsAsFactors = FALSE
      )
      out <- out[is.finite(out$score_p) & is.finite(out$eta2), , drop = FALSE]
      out
    }

    stat_ttest_feature_table <- function(
      X,
      idx_inside,
      idx_outside,
      mz_values = NULL,
      use_adjusted = TRUE
    ) {
      if (length(idx_inside) < 2L || length(idx_outside) < 2L) {
        return(data.frame())
      }

      X1 <- X[, idx_inside, drop = FALSE]
      X2 <- X[, idx_outside, drop = FALSE]

      n1 <- rowSums(is.finite(X1))
      n2 <- rowSums(is.finite(X2))

      X1z <- X1
      X2z <- X2
      X1z[!is.finite(X1z)] <- 0
      X2z[!is.finite(X2z)] <- 0

      s1 <- rowSums(X1z)
      s2 <- rowSums(X2z)
      ss1 <- rowSums(X1z * X1z)
      ss2 <- rowSums(X2z * X2z)

      mean1 <- s1 / pmax(n1, 1L)
      mean2 <- s2 / pmax(n2, 1L)
      var1 <- (ss1 - (s1 * s1) / pmax(n1, 1L)) / pmax(n1 - 1L, 1L)
      var2 <- (ss2 - (s2 * s2) / pmax(n2, 1L)) / pmax(n2 - 1L, 1L)

      se2 <- (var1 / pmax(n1, 1L)) + (var2 / pmax(n2, 1L))
      tstat <- (mean1 - mean2) / sqrt(se2)
      df <- (se2 * se2) / (
        (var1 * var1) / (pmax(n1, 1L) * pmax(n1, 1L) * pmax(n1 - 1L, 1L)) +
          (var2 * var2) / (pmax(n2, 1L) * pmax(n2, 1L) * pmax(n2 - 1L, 1L))
      )

      bad <- !is.finite(tstat) | !is.finite(df) | df <= 0 | !is.finite(se2) | se2 <= 0
      pval <- rep(NA_real_, length(tstat))
      pval[!bad] <- 2 * stats::pt(abs(tstat[!bad]), df = df[!bad], lower.tail = FALSE)
      padj <- stats::p.adjust(pval, method = "BH")
      score_p <- if (isTRUE(use_adjusted)) padj else pval

      lfc <- log2((mean1 + 1e-12) / (mean2 + 1e-12))
      abs_lfc <- abs(lfc)
      neglogp <- -log10(score_p + 1e-300)
      neglogp[!is.finite(neglogp)] <- 0

      if (is.null(mz_values) || length(mz_values) != nrow(X)) {
        mz_values <- rep(NA_real_, nrow(X))
      }

      out <- data.frame(
        feature_index = seq_len(nrow(X)),
        mz = suppressWarnings(as.numeric(mz_values)),
        n_inside = as.integer(n1),
        n_outside = as.integer(n2),
        mean_inside = as.numeric(mean1),
        mean_outside = as.numeric(mean2),
        log2FC = as.numeric(lfc),
        abs_log2FC = as.numeric(abs_lfc),
        t_stat = as.numeric(tstat),
        p_value = as.numeric(pval),
        p_adj = as.numeric(padj),
        score_p = as.numeric(score_p),
        neglog10_score_p = as.numeric(neglogp),
        stringsAsFactors = FALSE
      )
      out <- out[is.finite(out$score_p) & is.finite(out$abs_log2FC), , drop = FALSE]
      out
    }

    output$optimize_xy_preview_ui <- renderUI({
      cand <- xh$opt_xy_candidates
      if (is.null(cand) || nrow(cand) == 0) {
        return(tags$small("Auto-fit preview: run Auto-fit XY to generate up to 5 candidates."))
      }
      labels <- sprintf(
        "#%d dX=%.2f dY=%.2f score=%.4f",
        cand$rank, cand$dX, cand$dY, cand$score
      )
      vals <- as.character(cand$rank)
      names(vals) <- labels
      selectInput(
        ns("optimize_xy_choice"),
        "Auto-fit candidates (top 5)",
        choices = vals,
        selected = vals[1]
      )
    })

    observeEvent(input$apply_optimize_xy_choice, {
      cand <- xh$opt_xy_candidates
      if (is.null(cand) || nrow(cand) == 0) {
        showNotification("No auto-fit candidates available yet. Run Auto-fit XY first.", type = "warning", duration = 6)
        return()
      }
      pick_rank <- suppressWarnings(as.integer(input$optimize_xy_choice))
      if (!is.finite(pick_rank)) pick_rank <- cand$rank[1]
      row <- cand[cand$rank == pick_rank, , drop = FALSE]
      if (nrow(row) == 0) row <- cand[1, , drop = FALSE]

      updateSliderInput(session, "translate_x", value = row$translate_x[1])
      updateSliderInput(session, "translate_y", value = row$translate_y[1])
      updateNumericInput(session, "translate_x_num", value = row$translate_x[1])
      updateNumericInput(session, "translate_y_num", value = row$translate_y[1])

      showNotification(
        sprintf(
          "Applied candidate #%d: dX=%.2f, dY=%.2f, score=%.4f",
          row$rank[1], row$dX[1], row$dY[1], row$score[1]
        ),
        type = "message",
        duration = 6
      )
    }, ignoreInit = TRUE)

    output$stat_fit_preview_ui <- renderUI({
      cand <- xh$stat_fit_candidates
      if (is.null(cand) || nrow(cand) == 0) {
        return(tags$small("Stat-fit preview: run Stat Fit to generate up to 5 candidates."))
      }
      labels <- sprintf(
        "#%d dX=%d dY=%d n_sig=%d score=%.2f",
        cand$rank, cand$dX, cand$dY, cand$n_sig, cand$score
      )
      vals <- as.character(cand$rank)
      names(vals) <- labels
      selectInput(
        ns("stat_fit_choice"),
        "Stat-fit candidates (top 5)",
        choices = vals,
        selected = vals[1]
      )
    })

    observeEvent(input$apply_stat_fit_choice, {
      cand <- xh$stat_fit_candidates
      if (is.null(cand) || nrow(cand) == 0) {
        showNotification("No stat-fit candidates available yet. Run Stat Fit first.", type = "warning", duration = 6)
        return()
      }
      pick_rank <- suppressWarnings(as.integer(input$stat_fit_choice))
      if (!is.finite(pick_rank)) pick_rank <- cand$rank[1]
      row <- cand[cand$rank == pick_rank, , drop = FALSE]
      if (nrow(row) == 0) row <- cand[1, , drop = FALSE]

      updateSliderInput(session, "translate_x", value = row$translate_x[1])
      updateSliderInput(session, "translate_y", value = row$translate_y[1])
      updateNumericInput(session, "translate_x_num", value = row$translate_x[1])
      updateNumericInput(session, "translate_y_num", value = row$translate_y[1])

      showNotification(
        sprintf(
          "Applied stat-fit candidate #%d: dX=%d, dY=%d, n_sig=%d, score=%.2f",
          row$rank[1], row$dX[1], row$dY[1], row$n_sig[1], row$score[1]
        ),
        type = "message",
        duration = 6
      )
    }, ignoreInit = TRUE)

    observeEvent(input$stat_fit_metric_mode, {
      mode <- tolower(trimws(as.character(input$stat_fit_metric_mode)[1]))
      if (identical(mode, "polygon_cluster_groups")) {
        updateSelectInput(session, "stat_fit_objective", selected = "max")
      } else if (identical(mode, "inside_outside")) {
        updateSelectInput(session, "stat_fit_objective", selected = "min")
      }
    }, ignoreInit = TRUE)

    observeEvent(input$stat_fit_max_info, {
      req(msi_for_pdata())
      stat_metric_mode_pre <- tolower(trimws(as.character(input$stat_fit_metric_mode)[1]))
      if (!stat_metric_mode_pre %in% c("inside_outside", "polygon_cluster_groups")) stat_metric_mode_pre <- "inside_outside"
      stat_group_spec_pre <- parse_stat_fit_group_field()
      need_polygon_file <- !(
        identical(stat_metric_mode_pre, "polygon_cluster_groups") &&
          isTRUE(stat_fit_group_field_uses_pdata(stat_group_spec_pre$raw))
      )
      if (isTRUE(need_polygon_file) && (is.null(input$polygon_file) || !nzchar(input$polygon_file$name))) {
        showNotification("Load a polygon file before running max info.", type = "warning", duration = 7)
        return()
      }
      nid_work <- showNotification(
        "Generating Stat Fit max-info output (feature ranking) ...",
        type = "message",
        duration = NULL,
        closeButton = FALSE
      )
      on.exit({
        if (!is.null(nid_work)) {
          try(removeNotification(nid_work), silent = TRUE)
        }
      }, add = TRUE)

      tx <- suppressWarnings(as.numeric(input$translate_x))
      ty <- suppressWarnings(as.numeric(input$translate_y))
      if (!is.finite(tx)) tx <- 0
      if (!is.finite(ty)) ty <- 0

      outside_mode <- tolower(trimws(as.character(input$stat_fit_outside_mode)[1]))
      if (!outside_mode %in% c("local", "bbox", "global")) outside_mode <- "bbox"
      stat_metric_mode <- tolower(trimws(as.character(input$stat_fit_metric_mode)[1]))
      if (!stat_metric_mode %in% c("inside_outside", "polygon_cluster_groups")) stat_metric_mode <- "inside_outside"
      stat_group_spec <- parse_stat_fit_group_field()
      stat_objective <- tolower(trimws(as.character(input$stat_fit_objective)[1]))
      if (!stat_objective %in% c("min", "max")) stat_objective <- if (identical(stat_metric_mode, "polygon_cluster_groups")) "max" else "min"

      buf <- suppressWarnings(as.integer(input$stat_fit_buffer_px))
      bbox_pad <- suppressWarnings(as.integer(input$stat_fit_bbox_pad))
      min_pix <- suppressWarnings(as.integer(input$stat_fit_min_pixels))
      max_ions <- suppressWarnings(as.integer(input$stat_fit_max_ions))
      top_n <- suppressWarnings(as.integer(input$stat_fit_top_n))
      use_adjusted <- isTRUE(input$stat_fit_use_adjusted)

      if (!is.finite(buf) || buf < 1L) buf <- 3L
      if (!is.finite(bbox_pad) || bbox_pad < 0L) bbox_pad <- 25L
      if (!is.finite(min_pix) || min_pix < 2L) min_pix <- 10L
      if (!is.finite(max_ions) || max_ions < 20L) max_ions <- 200L
      if (!is.finite(top_n) || top_n < 1L) top_n <- 10L
      buf <- as.integer(min(25L, max(1L, buf)))
      bbox_pad <- as.integer(min(500L, max(0L, bbox_pad)))
      min_pix <- as.integer(min(1000L, max(2L, min_pix)))
      max_ions <- as.integer(min(2000L, max(20L, max_ions)))
      top_n <- as.integer(min(100L, max(1L, top_n)))

      msi <- make_msi_raster()
      pix_rows <- as.integer(msi$row_idx)
      pix_cols <- as.integer(msi$x_norm)
      valid_pix <- matrix(FALSE, nrow = as.integer(msi$ny), ncol = as.integer(msi$nx))
      valid_pix[cbind(pix_rows, pix_cols)] <- TRUE
      idx_in <- integer(0)
      idx_out <- integer(0)
      grp_idx <- integer(0)
      grp_levels_present <- character(0)
      n_in <- 0L
      n_out <- 0L
      n_grp <- 0L

      if (identical(stat_metric_mode, "polygon_cluster_groups")) {
        lab0 <- try(
          stat_fit_group_label_code_matrix(
            msi_obj = msi,
            tx = tx,
            ty = ty,
            group_field_sel = stat_group_spec$raw,
            exclude_labels = c("outside_polygon", "outside", "unassigned", "unclustered_polygon")
          ),
          silent = TRUE
        )
        if (inherits(lab0, "try-error") || is.null(lab0)) {
          showNotification("Could not build cluster-group labels for current transform. Check selected group field / polygon file.", type = "error", duration = 8)
          return()
        }
        grp_vec_full <- as.integer(lab0$code_mat[cbind(pix_rows, pix_cols)])
        valid_grp <- is.finite(grp_vec_full) & grp_vec_full > 0L
        if (!any(valid_grp)) {
          showNotification(
            sprintf("No assigned labels in '%s' (source: %s) at current transform after exclusions.", lab0$group_field, lab0$group_source),
            type = "warning",
            duration = 8
          )
          return()
        }
        tab_grp <- table(grp_vec_full[valid_grp])
        keep_codes <- suppressWarnings(as.integer(names(tab_grp)[tab_grp >= min_pix]))
        keep_codes <- keep_codes[is.finite(keep_codes) & keep_codes > 0L]
        if (length(keep_codes) < 2L) {
          showNotification(
            sprintf("Need at least two groups in '%s' (source: %s) with >= %d pixels at current transform.", lab0$group_field, lab0$group_source, min_pix),
            type = "warning",
            duration = 8
          )
          return()
        }
        idx_in <- which(valid_grp & (grp_vec_full %in% keep_codes))
        grp_idx <- match(grp_vec_full[idx_in], keep_codes)
        grp_levels_present <- as.character(lab0$group_levels[keep_codes])
        n_in <- as.integer(length(idx_in))
        n_out <- as.integer(sum(!valid_grp | !(grp_vec_full %in% keep_codes), na.rm = TRUE))
        n_grp <- as.integer(length(unique(grp_idx)))
      } else {
        mask_i <- try(polygon_mask_matrix(msi, tx = tx, ty = ty), silent = TRUE)
        if (inherits(mask_i, "try-error") || is.null(mask_i)) {
          showNotification("Could not build polygon mask for current transform.", type = "error", duration = 8)
          return()
        }

        outside_i <- switch(
          outside_mode,
          local = local_outside_mask(mask_i, buffer_px = buf),
          bbox = bbox_outside_mask(mask_i, pad_px = bbox_pad),
          global = !mask_i,
          bbox_outside_mask(mask_i, pad_px = bbox_pad)
        )
        outside_i <- outside_i & !mask_i
        outside_i <- outside_i & valid_pix

        inside_vec <- mask_i[cbind(pix_rows, pix_cols)]
        outside_vec <- outside_i[cbind(pix_rows, pix_cols)]
        n_in <- as.integer(sum(inside_vec, na.rm = TRUE))
        n_out <- as.integer(sum(outside_vec, na.rm = TRUE))
        if (n_in < min_pix || n_out < min_pix) {
          showNotification(
            sprintf("Not enough pixels for feature ranking at current transform (inside=%d, outside=%d; min=%d).", n_in, n_out, min_pix),
            type = "warning",
            duration = 8
          )
          return()
        }
        idx_in <- which(inside_vec)
        idx_out <- which(outside_vec)
      }

      obj <- msi_for_pdata()
      Xfull <- try(get_intensity_matrix_fast(obj), silent = TRUE)
      if (inherits(Xfull, "try-error") || is.null(Xfull)) {
        showNotification("Could not extract MSI intensities.", type = "error", duration = 8)
        return()
      }
      feat_idx <- select_top_variable_features(Xfull, max_features = max_ions)
      X <- Xfull[feat_idx, , drop = FALSE]
      rm(Xfull)

      mz_axis <- suppressWarnings(as.numeric(Cardinal::mz(obj)))
      mz_vals <- rep(NA_real_, length(feat_idx))
      if (length(mz_axis) > 0) {
        good_idx <- feat_idx[feat_idx >= 1L & feat_idx <= length(mz_axis)]
        mz_vals[feat_idx >= 1L & feat_idx <= length(mz_axis)] <- mz_axis[good_idx]
      }

      det <- rowMeans(is.finite(X) & X > 0, na.rm = TRUE)
      keep <- is.finite(det) & det >= 0.05
      if (sum(keep) >= 20L) {
        X <- X[keep, , drop = FALSE]
        mz_vals <- mz_vals[keep]
      }
      if (nrow(X) < 10L) {
        showNotification("Too few informative ions after filtering.", type = "warning", duration = 8)
        return()
      }

      if (!identical(input$intensity_transform, "none")) {
        X <- transform_intensity(X, input$intensity_transform)
      }

      if (identical(stat_metric_mode, "polygon_cluster_groups")) {
        tbl <- stat_group_anova_feature_table(
          X,
          idx_pixels = idx_in,
          group_index = grp_idx,
          mz_values = mz_vals,
          use_adjusted = use_adjusted
        )
      } else {
        tbl <- stat_ttest_feature_table(
          X,
          idx_inside = idx_in,
          idx_outside = idx_out,
          mz_values = mz_vals,
          use_adjusted = use_adjusted
        )
      }
      if (nrow(tbl) == 0L) {
        showNotification("No valid feature statistics were computed.", type = "warning", duration = 7)
        return()
      }

      if (identical(stat_metric_mode, "polygon_cluster_groups")) {
        tbl$fit_score <- tbl$eta2 * tbl$neglog10_score_p
        tbl$inform_abs <- tbl$mean_range * tbl$neglog10_score_p
      } else {
        tbl$fit_score <- if (identical(stat_objective, "min")) {
          (-tbl$log2FC) * tbl$neglog10_score_p
        } else {
          tbl$log2FC * tbl$neglog10_score_p
        }
        tbl$inform_abs <- tbl$abs_log2FC * tbl$neglog10_score_p
      }
      tbl$fit_score[!is.finite(tbl$fit_score)] <- -Inf
      tbl$inform_abs[!is.finite(tbl$inform_abs)] <- 0
      tbl <- tbl[order(tbl$fit_score, tbl$inform_abs, decreasing = TRUE), , drop = FALSE]
      top_tbl <- utils::head(tbl, top_n)
      # Display order in console: ascending p-value for easier interpretation,
      # while preserving the existing fit-score ranking used to choose the top set.
      if (nrow(top_tbl) > 1L) {
        p_ord <- order(
          ifelse(is.finite(top_tbl$p_value), top_tbl$p_value, Inf),
          ifelse(is.finite(top_tbl$p_adj), top_tbl$p_adj, Inf),
          -ifelse(is.finite(top_tbl$inform_abs), top_tbl$inform_abs, 0),
          na.last = TRUE
        )
        top_tbl <- top_tbl[p_ord, , drop = FALSE]
      }

      cat("\n[Stat Fit Max Info] Top features for current transform\n")
      if (identical(stat_metric_mode, "polygon_cluster_groups")) {
        cat(sprintf(
          "transform: translate_x=%.3f translate_y=%.3f | metric=polygon_cluster_groups | field=%s (%s) | objective=%s | groups=%d | grouped_pixels=%d excluded=%d\n",
          tx, ty, stat_group_spec$field, ifelse(is.null(lab0$group_source), "unknown", lab0$group_source), stat_objective, n_grp, n_in, n_out
        ))
        if (length(grp_levels_present) > 0L) {
          cat("groups used: ", paste(grp_levels_present, collapse = ", "), "\n", sep = "")
        }
        print(top_tbl[, c(
          "feature_index", "mz", "fit_score", "neglog10_score_p", "eta2", "mean_range",
          "F_stat", "n_groups", "n_pixels", "p_value", "p_adj"
        ), drop = FALSE], row.names = FALSE)
      } else {
        cat(sprintf(
          "transform: translate_x=%.3f translate_y=%.3f | metric=inside_outside | objective=%s | outside_mode=%s | inside=%d outside=%d\n",
          tx, ty, stat_objective, outside_mode, n_in, n_out
        ))
        print(top_tbl[, c(
          "feature_index", "mz", "fit_score", "neglog10_score_p", "log2FC", "abs_log2FC",
          "mean_inside", "mean_outside", "t_stat", "p_value", "p_adj"
        ), drop = FALSE], row.names = FALSE)
      }
      flush.console()

      showNotification(sprintf("Printed top %d informative features to console.", nrow(top_tbl)), type = "message", duration = 6)
    }, ignoreInit = TRUE)

    observeEvent(input$run_stat_fit, {
      req(msi_for_pdata())
      stat_metric_mode_pre <- tolower(trimws(as.character(input$stat_fit_metric_mode)[1]))
      if (!stat_metric_mode_pre %in% c("inside_outside", "polygon_cluster_groups")) stat_metric_mode_pre <- "inside_outside"
      stat_group_spec_pre <- parse_stat_fit_group_field()
      need_polygon_file <- !(
        identical(stat_metric_mode_pre, "polygon_cluster_groups") &&
          isTRUE(stat_fit_group_field_uses_pdata(stat_group_spec_pre$raw))
      )
      if (isTRUE(need_polygon_file) && (is.null(input$polygon_file) || !nzchar(input$polygon_file$name))) {
        showNotification("Load a polygon file before running Stat Fit.", type = "warning", duration = 7)
        return()
      }

      tx0 <- suppressWarnings(as.numeric(input$translate_x))
      ty0 <- suppressWarnings(as.numeric(input$translate_y))
      if (!is.finite(tx0)) tx0 <- 0
      if (!is.finite(ty0)) ty0 <- 0

      rng <- suppressWarnings(as.integer(input$stat_fit_range))
      stp <- suppressWarnings(as.integer(input$stat_fit_step))
      buf <- suppressWarnings(as.integer(input$stat_fit_buffer_px))
      bbox_pad <- suppressWarnings(as.integer(input$stat_fit_bbox_pad))
      min_pix <- suppressWarnings(as.integer(input$stat_fit_min_pixels))
      max_ions <- suppressWarnings(as.integer(input$stat_fit_max_ions))
      alpha <- suppressWarnings(as.numeric(input$stat_fit_alpha))
      use_adjusted <- isTRUE(input$stat_fit_use_adjusted)
      use_abs_lfc <- isTRUE(input$stat_fit_use_abs_lfc)
      stat_metric_mode <- tolower(trimws(as.character(input$stat_fit_metric_mode)[1]))
      if (!stat_metric_mode %in% c("inside_outside", "polygon_cluster_groups")) stat_metric_mode <- "inside_outside"
      stat_group_spec <- parse_stat_fit_group_field()
      outside_mode <- tolower(trimws(as.character(input$stat_fit_outside_mode)[1]))
      if (!outside_mode %in% c("local", "bbox", "global")) outside_mode <- "bbox"
      stat_objective <- tolower(trimws(as.character(input$stat_fit_objective)[1]))
      if (!stat_objective %in% c("min", "max")) stat_objective <- if (identical(stat_metric_mode, "polygon_cluster_groups")) "max" else "min"

      if (!is.finite(rng) || rng < 1L) rng <- 20L
      if (!is.finite(stp) || stp < 1L) stp <- 2L
      if (!is.finite(buf) || buf < 1L) buf <- 3L
      if (!is.finite(bbox_pad) || bbox_pad < 0L) bbox_pad <- 25L
      if (!is.finite(min_pix) || min_pix < 2L) min_pix <- 10L
      if (!is.finite(max_ions) || max_ions < 20L) max_ions <- 200L
      if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) alpha <- 0.1
      rng <- as.integer(min(300L, max(1L, rng)))
      stp <- as.integer(min(50L, max(1L, stp)))
      buf <- as.integer(min(25L, max(1L, buf)))
      bbox_pad <- as.integer(min(500L, max(0L, bbox_pad)))
      min_pix <- as.integer(min(1000L, max(2L, min_pix)))
      max_ions <- as.integer(min(2000L, max(20L, max_ions)))

      msi <- make_msi_raster()
      withProgress(message = "Stat-fit grid search", value = 0.05, {
        mask0 <- NULL
        group_code0 <- NULL
        group_levels0 <- character(0)
        if (identical(stat_metric_mode, "polygon_cluster_groups")) {
          incProgress(0.10, detail = "Rasterizing polygon cluster-group labels")
          lab0 <- try(
            stat_fit_group_label_code_matrix(
              msi_obj = msi,
              tx = tx0,
              ty = ty0,
              group_field_sel = stat_group_spec$raw,
              exclude_labels = c("outside_polygon", "outside", "unassigned", "unclustered_polygon")
            ),
            silent = TRUE
          )
          if (inherits(lab0, "try-error") || is.null(lab0)) {
            xh$stat_fit_grid <- NULL
            xh$stat_fit_candidates <- NULL
            xh$stat_fit_summary <- NULL
            showNotification("Could not build cluster-group labels for Stat Fit. Check selected group field / polygon file (or use pData field).", type = "error", duration = 8)
            return()
          }
          group_code0 <- lab0$code_mat
          group_levels0 <- as.character(lab0$group_levels)
          if (length(group_levels0) < 2L) {
            xh$stat_fit_grid <- NULL
            xh$stat_fit_candidates <- NULL
            xh$stat_fit_summary <- NULL
            showNotification(
              sprintf(
                "Stat Fit group metric needs at least two groups in '%s' (source: %s), excluding outside/unassigned/unclustered.",
                lab0$group_field,
                lab0$group_source
              ),
              type = "warning",
              duration = 8
            )
            return()
          }
        } else {
          incProgress(0.10, detail = "Rasterizing base polygon mask")
          mask0 <- try(polygon_mask_matrix(msi, tx = tx0, ty = ty0), silent = TRUE)
          if (inherits(mask0, "try-error") || is.null(mask0)) {
            xh$stat_fit_grid <- NULL
            xh$stat_fit_candidates <- NULL
            xh$stat_fit_summary <- NULL
            showNotification("Could not build polygon mask for Stat Fit. Check polygon scaling/position.", type = "error", duration = 8)
            return()
          }

          n_mask <- suppressWarnings(sum(mask0, na.rm = TRUE))
          if (!is.finite(n_mask) || n_mask < 10L) {
            xh$stat_fit_grid <- NULL
            xh$stat_fit_candidates <- NULL
            xh$stat_fit_summary <- NULL
            showNotification("Polygon mask is too small for Stat Fit.", type = "warning", duration = 8)
            return()
          }
        }

        incProgress(0.15, detail = "Extracting MSI intensity matrix")
        obj <- msi_for_pdata()
        Xfull <- try(get_intensity_matrix_fast(obj), silent = TRUE)
        if (inherits(Xfull, "try-error") || is.null(Xfull)) {
          xh$stat_fit_grid <- NULL
          xh$stat_fit_candidates <- NULL
          xh$stat_fit_summary <- NULL
          showNotification("Could not extract MSI intensities for Stat Fit.", type = "error", duration = 8)
          return()
        }

        feat_idx <- select_top_variable_features(Xfull, max_features = max_ions)
        X <- Xfull[feat_idx, , drop = FALSE]
        rm(Xfull)
        det <- rowMeans(is.finite(X) & X > 0, na.rm = TRUE)
        keep <- is.finite(det) & det >= 0.05
        if (sum(keep) >= 20L) {
          X <- X[keep, , drop = FALSE]
        }
        if (nrow(X) < 10L) {
          xh$stat_fit_grid <- NULL
          xh$stat_fit_candidates <- NULL
          xh$stat_fit_summary <- NULL
          showNotification("Too few informative ions after filtering for Stat Fit.", type = "warning", duration = 8)
          return()
        }

        if (!identical(input$intensity_transform, "none")) {
          X <- transform_intensity(X, input$intensity_transform)
        }

        pix_rows <- as.integer(msi$row_idx)
        pix_cols <- as.integer(msi$x_norm)
        valid_pix <- matrix(FALSE, nrow = as.integer(msi$ny), ncol = as.integer(msi$nx))
        valid_pix[cbind(pix_rows, pix_cols)] <- TRUE
        dx_vals <- seq.int(-rng, rng, by = stp)
        dy_vals <- seq.int(-rng, rng, by = stp)
        grid <- expand.grid(dx = dx_vals, dy = dy_vals, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
        grid$n_inside <- NA_integer_
        grid$n_outside <- NA_integer_
        grid$n_groups <- NA_integer_
        grid$n_sig <- NA_integer_
        grid$score <- NA_real_
        grid$med_lfc <- NA_real_
        grid$mean_logp <- NA_real_
        grid$mean_abs_t_top <- NA_real_

        n_total <- nrow(grid)
        step_update <- max(1L, as.integer(n_total / 80L))

        eval_stat_fit_candidate <- function(i) {
          tryCatch({
            dx <- as.integer(grid$dx[i])
            dy <- as.integer(grid$dy[i])
            if (identical(stat_metric_mode, "polygon_cluster_groups")) {
              lab_i <- shift_int(group_code0, dr = dy, dc = dx, fill = 0L)
              grp_vec <- as.integer(lab_i[cbind(pix_rows, pix_cols)])
              valid_grp <- is.finite(grp_vec) & grp_vec > 0L
              n_assigned <- as.integer(sum(valid_grp, na.rm = TRUE))
              if (n_assigned < (2L * min_pix)) {
                return(list(
                  n_inside = n_assigned,
                  n_outside = as.integer(length(grp_vec) - n_assigned),
                  n_groups = as.integer(sum(table(grp_vec[valid_grp]) >= min_pix)),
                  n_sig = NA_integer_,
                  score = NA_real_,
                  med_lfc = NA_real_,
                  mean_logp = NA_real_,
                  mean_abs_t_top = NA_real_,
                  error = NULL
                ))
              }

              tab_grp <- table(grp_vec[valid_grp])
              keep_codes <- suppressWarnings(as.integer(names(tab_grp)[tab_grp >= min_pix]))
              keep_codes <- keep_codes[is.finite(keep_codes) & keep_codes > 0L]
              n_groups_kept <- as.integer(length(keep_codes))
              if (n_groups_kept < 2L) {
                return(list(
                  n_inside = n_assigned,
                  n_outside = as.integer(length(grp_vec) - n_assigned),
                  n_groups = n_groups_kept,
                  n_sig = NA_integer_,
                  score = NA_real_,
                  med_lfc = NA_real_,
                  mean_logp = NA_real_,
                  mean_abs_t_top = NA_real_,
                  error = NULL
                ))
              }

              idx_grp <- which(valid_grp & (grp_vec %in% keep_codes))
              grp_idx <- match(grp_vec[idx_grp], keep_codes)
              mt <- stat_group_anova_count_sig(
                X,
                idx_pixels = idx_grp,
                group_index = grp_idx,
                alpha = alpha,
                use_adjusted = use_adjusted
              )
              list(
                n_inside = as.integer(length(idx_grp)),
                n_outside = as.integer(length(grp_vec) - length(idx_grp)),
                n_groups = n_groups_kept,
                n_sig = as.integer(mt$n_sig),
                score = as.numeric(mt$score),
                med_lfc = as.numeric(mt$med_lfc),
                mean_logp = as.numeric(mt$mean_logp),
                mean_abs_t_top = as.numeric(mt$mean_abs_t_top),
                error = NULL
              )
            } else {
              mask_i <- shift_bool(mask0, dr = dy, dc = dx)
              outside_i <- switch(
                outside_mode,
                local = local_outside_mask(mask_i, buffer_px = buf),
                bbox = bbox_outside_mask(mask_i, pad_px = bbox_pad),
                global = !mask_i,
                bbox_outside_mask(mask_i, pad_px = bbox_pad)
              )
              outside_i <- outside_i & !mask_i
              outside_i <- outside_i & valid_pix

              inside_vec <- mask_i[cbind(pix_rows, pix_cols)]
              outside_vec <- outside_i[cbind(pix_rows, pix_cols)]
              n_in <- as.integer(sum(inside_vec, na.rm = TRUE))
              n_out <- as.integer(sum(outside_vec, na.rm = TRUE))

              if (n_in < min_pix || n_out < min_pix) {
                return(list(
                  n_inside = n_in,
                  n_outside = n_out,
                  n_groups = NA_integer_,
                  n_sig = NA_integer_,
                  score = NA_real_,
                  med_lfc = NA_real_,
                  mean_logp = NA_real_,
                  mean_abs_t_top = NA_real_,
                  error = NULL
                ))
              }

              idx_in <- which(inside_vec)
              idx_out <- which(outside_vec)
              mt <- stat_ttest_count_sig(
                X,
                idx_inside = idx_in,
                idx_outside = idx_out,
                alpha = alpha,
                use_adjusted = use_adjusted,
                use_abs_lfc = use_abs_lfc
              )
              list(
                n_inside = n_in,
                n_outside = n_out,
                n_groups = NA_integer_,
                n_sig = as.integer(mt$n_sig),
                score = as.numeric(mt$score),
                med_lfc = as.numeric(mt$med_lfc),
                mean_logp = as.numeric(mt$mean_logp),
                mean_abs_t_top = as.numeric(mt$mean_abs_t_top),
                error = NULL
              )
            }
          }, error = function(e) {
            list(
              n_inside = NA_integer_,
              n_outside = NA_integer_,
              n_groups = NA_integer_,
              n_sig = NA_integer_,
              score = NA_real_,
              med_lfc = NA_real_,
              mean_logp = NA_real_,
              mean_abs_t_top = NA_real_,
              error = conditionMessage(e)
            )
          })
        }

        apply_stat_fit_result <- function(i, res) {
          grid$n_inside[i] <<- as.integer(res$n_inside)
          grid$n_outside[i] <<- as.integer(res$n_outside)
          grid$n_groups[i] <<- as.integer(res$n_groups)
          grid$n_sig[i] <<- as.integer(res$n_sig)
          grid$score[i] <<- as.numeric(res$score)
          grid$med_lfc[i] <<- as.numeric(res$med_lfc)
          grid$mean_logp[i] <<- as.numeric(res$mean_logp)
          grid$mean_abs_t_top[i] <<- as.numeric(res$mean_abs_t_top)
        }

        stat_fit_errors <- character(0)
        note_stat_fit_error <- function(i, res) {
          if (!is.null(res$error) && nzchar(res$error)) {
            stat_fit_errors <<- c(stat_fit_errors, sprintf("idx=%d dX=%d dY=%d: %s", i, grid$dx[i], grid$dy[i], res$error))
          }
        }

        ncores_req <- suppressWarnings(as.integer(tryCatch(setup_values()[["ncores"]], error = function(e) NA_integer_)))
        if (!is.finite(ncores_req) || ncores_req < 1L) ncores_req <- 1L
        ncores_detect <- suppressWarnings(as.integer(tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA_integer_)))
        if (!is.finite(ncores_detect) || ncores_detect < 1L) ncores_detect <- ncores_req
        ncores_use <- as.integer(max(1L, min(ncores_req, ncores_detect)))
        can_parallel <- (.Platform$OS.type != "windows") &&
          requireNamespace("parallel", quietly = TRUE) &&
          ncores_use > 1L &&
          n_total >= (2L * ncores_use)

        message(sprintf(
          "[Histology Stat Fit] Grid candidates=%d | mode=%s | cores=%d (%s)",
          n_total,
          if (can_parallel) "parallel" else "serial",
          ncores_use,
          .Platform$OS.type
        ))

        if (can_parallel) {
          incProgress(0, detail = sprintf("Scoring translations in parallel (%d cores)", ncores_use))
          chunk_size <- max(16L, 8L * ncores_use)
          chunks <- split(seq_len(n_total), ceiling(seq_len(n_total) / chunk_size))
          parallel_failed <- FALSE

          for (chunk_idx in chunks) {
            res_chunk <- try(
              parallel::mclapply(
                chunk_idx,
                eval_stat_fit_candidate,
                mc.cores = ncores_use,
                mc.preschedule = TRUE
              ),
              silent = TRUE
            )

            if (inherits(res_chunk, "try-error") || length(res_chunk) != length(chunk_idx)) {
              parallel_failed <- TRUE
              msg_pf <- if (inherits(res_chunk, "try-error")) as.character(res_chunk) else "unexpected result length"
              message("[Histology Stat Fit] Parallel evaluation failed; falling back to serial for remaining grid: ", msg_pf)
              start_i <- min(chunk_idx)
              break
            }

            for (k in seq_along(chunk_idx)) {
              i <- chunk_idx[k]
              res <- res_chunk[[k]]
              apply_stat_fit_result(i, res)
              note_stat_fit_error(i, res)
            }

            incProgress(0.70 * length(chunk_idx) / n_total)
          }

          if (isTRUE(parallel_failed)) {
            rem_idx <- which(!is.finite(grid$n_inside) & !is.finite(grid$n_outside))
            if (length(rem_idx) == 0L) {
              rem_idx <- seq.int(start_i, n_total)
            }
            for (i in rem_idx) {
              res <- eval_stat_fit_candidate(i)
              apply_stat_fit_result(i, res)
              note_stat_fit_error(i, res)
              if (i %% step_update == 0L) {
                incProgress(0.70 / ceiling(n_total / step_update))
              }
            }
          }
        } else {
          incProgress(0, detail = "Scoring translations")
          for (i in seq_len(n_total)) {
            res <- eval_stat_fit_candidate(i)
            apply_stat_fit_result(i, res)
            note_stat_fit_error(i, res)
            if (i %% step_update == 0L) {
              incProgress(0.70 / ceiling(n_total / step_update))
            }
          }
        }

        if (length(stat_fit_errors) > 0L) {
          message(sprintf("[Histology Stat Fit] %d grid evaluations returned errors (showing up to 5):", length(stat_fit_errors)))
          for (em in utils::head(stat_fit_errors, 5L)) message("  ", em)
        }

        cand <- grid[is.finite(grid$score), , drop = FALSE]
        if (nrow(cand) == 0L) {
          xh$stat_fit_grid <- grid
          xh$stat_fit_candidates <- NULL
          xh$stat_fit_summary <- NULL
          showNotification("Stat Fit found no valid candidates. Try smaller range, larger step, or lower min pixels.", type = "warning", duration = 8)
          return()
        }

        if (identical(stat_objective, "min")) {
          cand <- cand[order(cand$score, cand$n_sig, decreasing = c(FALSE, FALSE)), , drop = FALSE]
        } else {
          cand <- cand[order(cand$score, cand$n_sig, decreasing = c(TRUE, TRUE)), , drop = FALSE]
        }
        cand$translate_x <- clamp(tx0 + cand$dx, -1000, 1000)
        cand$translate_y <- clamp(ty0 + cand$dy, -1000, 1000)
        cand$dX <- as.integer(cand$dx)
        cand$dY <- as.integer(cand$dy)
        cand$rank <- seq_len(nrow(cand))

        top_n <- min(5L, nrow(cand))
        xh$stat_fit_grid <- grid
        xh$stat_fit_candidates <- cand[seq_len(top_n), c(
          "rank", "dX", "dY", "n_sig", "score", "med_lfc", "mean_logp", "mean_abs_t_top",
          "n_inside", "n_outside", "n_groups", "translate_x", "translate_y"
        ), drop = FALSE]

        best <- cand[1, , drop = FALSE]
        xh$stat_fit_summary <- list(
          best_dX = best$dX[1],
          best_dY = best$dY[1],
          best_translate_x = best$translate_x[1],
          best_translate_y = best$translate_y[1],
          best_n_sig = best$n_sig[1],
          best_n_groups = if ("n_groups" %in% names(best)) best$n_groups[1] else NA_integer_,
          best_score = best$score[1],
          best_med_lfc = best$med_lfc[1],
          best_mean_logp = best$mean_logp[1],
          best_mean_abs_t_top = best$mean_abs_t_top[1],
          n_candidates = nrow(cand),
          n_ions_tested = nrow(X),
          range = rng,
          step = stp,
          metric_mode = stat_metric_mode,
          metric_group_var = if (identical(stat_metric_mode, "polygon_cluster_groups")) if (!is.null(lab0$group_field)) lab0$group_field else stat_group_spec$field else NA_character_,
          metric_group_source = if (identical(stat_metric_mode, "polygon_cluster_groups")) if (!is.null(lab0$group_source)) lab0$group_source else stat_group_spec$source else NA_character_,
          n_group_levels_base = if (identical(stat_metric_mode, "polygon_cluster_groups")) length(group_levels0) else NA_integer_,
          objective = stat_objective,
          outside_mode = outside_mode,
          outside_buffer_px = buf,
          bbox_pad_px = bbox_pad,
          min_pixels_per_group = min_pix,
          alpha = alpha,
          adjusted_p = use_adjusted
        )

        showNotification(
          sprintf(
            "Stat Fit complete (%s; %s). Best dX=%d, dY=%d (n_sig=%d, score=%.2f).",
            if (identical(stat_objective, "min")) "min score" else "max score",
            if (identical(stat_metric_mode, "polygon_cluster_groups")) sprintf("%s groups (%s)", if (!is.null(lab0$group_field)) lab0$group_field else stat_group_spec$field, if (!is.null(lab0$group_source)) lab0$group_source else stat_group_spec$source) else "inside/outside",
            best$dX[1], best$dY[1], best$n_sig[1], best$score[1]
          ),
          type = "message",
          duration = 8
        )
      })
    }, ignoreInit = TRUE)

    output$stat_fit_heatmap_ui <- renderUI({
      if (requireNamespace("plotly", quietly = TRUE)) {
        plotly::plotlyOutput(ns("stat_fit_heatmap"), height = "240px")
      } else {
        plotOutput(ns("stat_fit_heatmap"), height = "240px")
      }
    })

    get_stat_fit_objective <- function() {
      obj <- NULL
      if (!is.null(xh$stat_fit_summary) && !is.null(xh$stat_fit_summary$objective)) {
        obj <- tolower(trimws(as.character(xh$stat_fit_summary$objective)[1]))
      }
      if (!obj %in% c("min", "max")) {
        obj <- tolower(trimws(as.character(input$stat_fit_objective)[1]))
      }
      if (!obj %in% c("min", "max")) obj <- "min"
      obj
    }

    get_stat_fit_best_row <- function(g2) {
      obj <- get_stat_fit_objective()
      idx <- if (identical(obj, "min")) which.min(g2$score) else which.max(g2$score)
      g2[idx, , drop = FALSE]
    }

    build_stat_fit_heatmap <- function(g2, best, with_tooltip = FALSE) {
      if (isTRUE(with_tooltip)) {
        g2$.tip <- sprintf(
          "dX: %d<br>dY: %d<br>n_sig: %d<br>score: %.2f<br>n_inside: %d<br>n_outside: %d",
          as.integer(g2$dx), as.integer(g2$dy), as.integer(g2$n_sig), as.numeric(g2$score),
          as.integer(g2$n_inside), as.integer(g2$n_outside)
        )
        p <- ggplot2::ggplot(g2, ggplot2::aes(x = dx, y = dy, fill = n_sig, text = .tip))
      } else {
        p <- ggplot2::ggplot(g2, ggplot2::aes(x = dx, y = dy, fill = n_sig))
      }

      p +
        ggplot2::geom_tile() +
        ggplot2::geom_point(data = best, ggplot2::aes(x = dx, y = dy), inherit.aes = FALSE, shape = 4, size = 4, stroke = 1.2, color = "red") +
        ggplot2::coord_equal() +
        ggplot2::scale_fill_viridis_c(option = "plasma", na.value = "grey85") +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "Stat-fit heatmap", subtitle = "Fill = number of significant ions", x = "dX", y = "dY", fill = "n_sig")
    }

    if (requireNamespace("plotly", quietly = TRUE)) {
      output$stat_fit_heatmap <- plotly::renderPlotly({
        if (!isTRUE(input$show_fit_info)) return(NULL)
        grid <- xh$stat_fit_grid
        validate(need(!is.null(grid) && nrow(grid) > 0, "Run Stat Fit to view diagnostic heatmap."))
        g2 <- grid[is.finite(grid$n_sig), , drop = FALSE]
        validate(need(nrow(g2) > 0, "No valid Stat Fit points to plot."))
        best <- get_stat_fit_best_row(g2)
        p <- build_stat_fit_heatmap(g2, best, with_tooltip = TRUE)
        plotly::ggplotly(p, tooltip = "text")
      })
    } else {
      output$stat_fit_heatmap <- renderPlot({
        if (!isTRUE(input$show_fit_info)) return(invisible(NULL))
        grid <- xh$stat_fit_grid
        validate(need(!is.null(grid) && nrow(grid) > 0, "Run Stat Fit to view diagnostic heatmap."))
        g2 <- grid[is.finite(grid$n_sig), , drop = FALSE]
        validate(need(nrow(g2) > 0, "No valid Stat Fit points to plot."))
        best <- get_stat_fit_best_row(g2)
        build_stat_fit_heatmap(g2, best, with_tooltip = FALSE)
      })
    }

    output$stat_fit_contour <- renderPlot({
      if (!isTRUE(input$show_fit_info)) return(invisible(NULL))
      grid <- xh$stat_fit_grid
      validate(need(!is.null(grid) && nrow(grid) > 0, "Run Stat Fit to view contour diagnostics."))
      g2 <- grid[is.finite(grid$score), , drop = FALSE]
      validate(need(nrow(g2) > 0, "No valid Stat Fit score values to plot."))
      best <- get_stat_fit_best_row(g2)
      ggplot2::ggplot(g2, ggplot2::aes(x = dx, y = dy, z = score)) +
        ggplot2::geom_contour(bins = 12, linewidth = 0.5) +
        ggplot2::geom_point(data = best, ggplot2::aes(x = dx, y = dy), shape = 4, size = 4, stroke = 1.2, color = "red") +
        ggplot2::coord_equal() +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "Stat-fit contours", subtitle = "Contours = composite statistical score", x = "dX", y = "dY")
    })

    output$stat_fit_summary <- renderPrint({
      if (!isTRUE(input$show_fit_info)) return(invisible(NULL))
      s <- xh$stat_fit_summary
      if (is.null(s)) {
        cat("Run Stat Fit to populate summary.\n")
      } else {
        print(s)
      }
    })

    observeEvent(input$optimize_xy, {
      req(msi_for_pdata())
      if (is.null(input$polygon_file) || !nzchar(input$polygon_file$name)) {
        showNotification("Load a polygon file before running Auto-fit XY.", type = "warning", duration = 7)
        return()
      }

      tx0 <- suppressWarnings(as.numeric(input$translate_x))
      ty0 <- suppressWarnings(as.numeric(input$translate_y))
      if (!is.finite(tx0)) tx0 <- 0
      if (!is.finite(ty0)) ty0 <- 0

      rng <- suppressWarnings(as.integer(input$optimize_xy_range))
      stp <- suppressWarnings(as.integer(input$optimize_xy_step))
      band_px <- suppressWarnings(as.integer(input$optimize_edge_band))
      if (!is.finite(rng) || rng < 1L) rng <- 40L
      if (!is.finite(stp) || stp < 1L) stp <- 2L
      if (!is.finite(band_px) || band_px < 1L) band_px <- 4L
      rng <- as.integer(min(500L, max(1L, rng)))
      stp <- as.integer(min(50L, max(1L, stp)))
      band_px <- as.integer(min(20L, max(1L, band_px)))

      msi <- make_msi_raster()
      withProgress(message = "Auto-fitting XY translation", value = 0.1, {
        incProgress(0.2, detail = "Preparing MSI target map")
        sig_info <- try(edge_fit_signal_matrix(msi), silent = TRUE)
        if (inherits(sig_info, "try-error") || is.null(sig_info) || !is.matrix(sig_info$signal)) {
          showNotification("Could not build Edge Fit optimization signal. Check selected pData field or display mode.", type = "error", duration = 8)
          return()
        }
        signal <- sig_info$signal
        message(sprintf(
          "[Histology Edge Fit] Optimization signal source=%s field=%s type=%s",
          if (!is.null(sig_info$source)) as.character(sig_info$source) else "unknown",
          if (!is.null(sig_info$field) && nzchar(as.character(sig_info$field))) as.character(sig_info$field) else "n/a",
          if (!is.null(sig_info$type) && nzchar(as.character(sig_info$type))) as.character(sig_info$type) else "n/a"
        ))
        edge_cache <- build_edge_distance_cache(signal)

        incProgress(0.2, detail = "Rasterizing polygon mask")

        # Build mask from multiple translation anchors so optimization can recover
        # from poor initial placement (where a single anchor may clip most polygons).
        anchor_vals <- unique(as.integer(round(c(-rng, -rng / 2, 0, rng / 2, rng))))
        anchor_grid <- expand.grid(
          dx = anchor_vals,
          dy = anchor_vals,
          KEEP.OUT.ATTRS = FALSE,
          stringsAsFactors = FALSE
        )

        mask0 <- NULL
        anchor_dx <- 0L
        anchor_dy <- 0L
        best_anchor_area <- -Inf
        for (i in seq_len(nrow(anchor_grid))) {
          adx <- as.integer(anchor_grid$dx[i])
          ady <- as.integer(anchor_grid$dy[i])
          m_try <- try(
            polygon_mask_matrix(msi, tx = tx0 + adx, ty = ty0 + ady),
            silent = TRUE
          )
          if (inherits(m_try, "try-error") || is.null(m_try)) next
          n_try <- suppressWarnings(sum(m_try, na.rm = TRUE))
          if (is.finite(n_try) && n_try > best_anchor_area) {
            best_anchor_area <- n_try
            mask0 <- m_try
            anchor_dx <- adx
            anchor_dy <- ady
          }
        }

        if (is.null(mask0)) {
          showNotification("Could not build polygon mask for optimization. Check polygon geometry.", type = "error", duration = 8)
          return()
        }
        n_mask <- sum(mask0, na.rm = TRUE)
        if (!is.finite(n_mask) || n_mask < 10) {
          showNotification("Polygon mask is too small for stable XY optimization.", type = "warning", duration = 8)
          return()
        }

        mask_cache <- build_mask_bands(mask0, band_px = band_px)

        dx_vals <- seq.int(-rng, rng, by = stp)
        dy_vals <- seq.int(-rng, rng, by = stp)

        ncores_req <- suppressWarnings(as.integer(tryCatch(setup_values()[["ncores"]], error = function(e) NA_integer_)))
        if (!is.finite(ncores_req) || ncores_req < 1L) ncores_req <- 1L
        ncores_detect <- suppressWarnings(as.integer(tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA_integer_)))
        if (!is.finite(ncores_detect) || ncores_detect < 1L) ncores_detect <- ncores_req
        ncores_use <- as.integer(max(1L, min(ncores_req, ncores_detect)))
        can_parallel_edge <- (.Platform$OS.type != "windows") &&
          requireNamespace("parallel", quietly = TRUE) &&
          ncores_use > 1L
        message(sprintf(
          "[Histology Edge Fit] Grid scoring mode=%s cores=%d",
          if (can_parallel_edge) "parallel" else "serial",
          ncores_use
        ))

        score_edge_pair <- function(dx, dy) {
          sc <- score_mask_shift(mask_cache, signal, dx = dx - anchor_dx, dy = dy - anchor_dy)
          sc_edge <- score_edge_shift(mask_cache, edge_cache, dx = dx - anchor_dx, dy = dy - anchor_dy)
          list(
            dX = as.integer(dx),
            dY = as.integer(dy),
            raw_score = as.numeric(sc),
            edge_score = as.numeric(sc_edge)
          )
        }

        score_edge_grid <- function(dx_set, dy_set, progress_weight = 0, detail = NULL) {
          if (!is.null(detail)) incProgress(0, detail = detail)
          grid_pairs <- expand.grid(
            dX = as.integer(dx_set),
            dY = as.integer(dy_set),
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE
          )
          if (nrow(grid_pairs) == 0L) return(data.frame())

          res_list <- vector("list", nrow(grid_pairs))
          if (isTRUE(can_parallel_edge) && nrow(grid_pairs) >= (2L * ncores_use)) {
            chunk_size <- max(64L, 16L * ncores_use)
            chunks <- split(seq_len(nrow(grid_pairs)), ceiling(seq_len(nrow(grid_pairs)) / chunk_size))
            parallel_failed <- FALSE
            for (chunk_idx in chunks) {
              chunk_df <- grid_pairs[chunk_idx, , drop = FALSE]
              out_chunk <- try(
                parallel::mclapply(
                  seq_len(nrow(chunk_df)),
                  function(k) score_edge_pair(dx = chunk_df$dX[k], dy = chunk_df$dY[k]),
                  mc.cores = ncores_use,
                  mc.preschedule = TRUE
                ),
                silent = TRUE
              )
              if (inherits(out_chunk, "try-error") || length(out_chunk) != length(chunk_idx)) {
                parallel_failed <- TRUE
                message("[Histology Edge Fit] Parallel grid scoring failed; falling back to serial: ", if (inherits(out_chunk, "try-error")) as.character(out_chunk) else "unexpected result length")
                break
              }
              res_list[chunk_idx] <- out_chunk
              if (progress_weight > 0) incProgress(progress_weight * length(chunk_idx) / nrow(grid_pairs))
            }
            if (isTRUE(parallel_failed)) {
              rem <- which(vapply(res_list, is.null, logical(1)))
              for (ii in rem) {
                res_list[[ii]] <- score_edge_pair(dx = grid_pairs$dX[ii], dy = grid_pairs$dY[ii])
                if (progress_weight > 0 && (ii %% max(1L, floor(nrow(grid_pairs) / 40L)) == 0L)) {
                  incProgress(progress_weight / max(1L, ceiling(nrow(grid_pairs) / max(1L, floor(nrow(grid_pairs) / 40L)))))
                }
              }
            }
          } else {
            step_update_local <- max(1L, as.integer(nrow(grid_pairs) / 40L))
            for (ii in seq_len(nrow(grid_pairs))) {
              res_list[[ii]] <- score_edge_pair(dx = grid_pairs$dX[ii], dy = grid_pairs$dY[ii])
              if (progress_weight > 0 && (ii %% step_update_local == 0L)) {
                incProgress(progress_weight / ceiling(nrow(grid_pairs) / step_update_local))
              }
            }
          }

          out <- do.call(rbind, lapply(res_list, function(z) {
            if (is.null(z)) return(NULL)
            data.frame(
              dX = as.integer(z$dX),
              dY = as.integer(z$dY),
              raw_score = as.numeric(z$raw_score),
              edge_score = as.numeric(z$edge_score),
              stringsAsFactors = FALSE
            )
          }))
          if (is.null(out) || nrow(out) == 0L) out <- data.frame()
          out
        }

        cand_dx <- integer(0)
        cand_dy <- integer(0)
        cand_sc <- numeric(0)
        cand_edge_sc <- numeric(0)

        collapse_candidates <- function(dx_vec, dy_vec, sc_vec, edge_sc_vec, rng_local, stp_local) {
          if (length(sc_vec) == 0L || length(edge_sc_vec) == 0L) return(data.frame())
          key <- paste(dx_vec, dy_vec, sep = ",")
          key_levels <- unique(key)
          max_sc <- vapply(key_levels, function(k) {
            vv <- sc_vec[key == k]
            if (all(!is.finite(vv))) return(NA_real_)
            max(vv, na.rm = TRUE)
          }, numeric(1))
          max_edge_sc <- vapply(key_levels, function(k) {
            vv <- edge_sc_vec[key == k]
            if (all(!is.finite(vv))) return(NA_real_)
            max(vv, na.rm = TRUE)
          }, numeric(1))
          parts <- do.call(rbind, strsplit(key_levels, ",", fixed = TRUE))
          out <- data.frame(
            dX = as.integer(parts[, 1]),
            dY = as.integer(parts[, 2]),
            raw_score = as.numeric(max_sc),
            edge_score = as.numeric(max_edge_sc),
            stringsAsFactors = FALSE
          )
          out <- out[is.finite(out$raw_score) | is.finite(out$edge_score), , drop = FALSE]
          if (nrow(out) == 0L) return(out)

          robust_z <- function(x) {
            x <- as.numeric(x)
            x[!is.finite(x)] <- NA_real_
            keep <- is.finite(x)
            if (!any(keep)) return(rep(0, length(x)))
            med <- suppressWarnings(stats::median(x[keep], na.rm = TRUE))
            madv <- suppressWarnings(stats::mad(x[keep], center = med, constant = 1, na.rm = TRUE))
            if (!is.finite(madv) || madv <= 1e-6) return(rep(0, length(x)))
            z <- (x - med) / madv
            z[!is.finite(z)] <- 0
            z
          }

          rng_used <- max(1, as.integer(rng_local))
          stp_used <- max(1, as.integer(stp_local))
          edge_band <- max(1, as.integer(2L * stp_used))

          out$dist_norm <- sqrt((out$dX / rng_used)^2 + (out$dY / rng_used)^2)
          out$dist_norm <- pmax(0, pmin(out$dist_norm, 2))
          out$edge_margin <- pmin(rng_used - abs(out$dX), rng_used - abs(out$dY))
          out$edge_norm <- pmax(0, pmin(1, out$edge_margin / edge_band))

          # Soft priors to reduce boundary drift:
          # - small displacement prior around current translation
          # - explicit penalty for hugging search bounds
          pen_shift <- 0.10 * (out$dist_norm^2)
          pen_edge <- 0.25 * ((1 - out$edge_norm)^2)
          z_contrast <- robust_z(out$raw_score)
          z_edge <- robust_z(out$edge_score)
          out$score <- 0.55 * z_contrast + 0.45 * z_edge - pen_shift - pen_edge
          out
        }

        pick_best_candidate <- function(df, rng_local, stp_local) {
          if (is.null(df) || nrow(df) == 0L) return(NULL)
          ord <- order(df$score, decreasing = TRUE)
          df_ord <- df[ord, , drop = FALSE]
          best <- df_ord[1, , drop = FALSE]

          edge_thr <- max(1L, as.integer(stp_local))
          near_edge <- (abs(best$dX[1]) >= (rng_local - edge_thr)) || (abs(best$dY[1]) >= (rng_local - edge_thr))
          if (isTRUE(near_edge)) {
            interior <- df_ord[
              abs(df_ord$dX) <= (rng_local - edge_thr) &
                abs(df_ord$dY) <= (rng_local - edge_thr),
              ,
              drop = FALSE
            ]
            if (nrow(interior) > 0L) {
              alt <- interior[which.max(interior$score), , drop = FALSE]
              tol <- max(0.02, 0.05 * abs(best$score[1]))
              if (is.finite(alt$score[1]) && alt$score[1] >= (best$score[1] - tol)) {
                best <- alt
              }
            }
          }
          best
        }

        coarse_df <- score_edge_grid(dx_vals, dy_vals, progress_weight = 0.45, detail = "Scoring coarse XY grid")
        if (nrow(coarse_df) > 0L) {
          keep_coarse <- is.finite(coarse_df$raw_score) | is.finite(coarse_df$edge_score)
          coarse_df <- coarse_df[keep_coarse, , drop = FALSE]
          cand_dx <- c(cand_dx, as.integer(coarse_df$dX))
          cand_dy <- c(cand_dy, as.integer(coarse_df$dY))
          cand_sc <- c(cand_sc, as.numeric(coarse_df$raw_score))
          cand_edge_sc <- c(cand_edge_sc, as.numeric(coarse_df$edge_score))
        }

        cand_df <- collapse_candidates(cand_dx, cand_dy, cand_sc, cand_edge_sc, rng_local = rng, stp_local = stp)
        best_row <- pick_best_candidate(cand_df, rng_local = rng, stp_local = stp)
        if (is.null(best_row) || nrow(best_row) == 0L || !is.finite(best_row$score[1])) {
          xh$opt_xy_candidates <- NULL
          showNotification(
            "Auto-fit XY could not find a valid translation score in this range. Try a smaller step, different MSI mode/ion, or check polygon scaling.",
            type = "warning",
            duration = 8
          )
          return()
        }

        # Local refinement at 1-pixel resolution near best stabilized coarse candidate.
        best_dx <- as.integer(best_row$dX[1])
        best_dy <- as.integer(best_row$dY[1])
        fine_dx <- seq.int(best_dx - stp, best_dx + stp, by = 1L)
        fine_dy <- seq.int(best_dy - stp, best_dy + stp, by = 1L)
        fine_dx <- fine_dx[fine_dx >= -rng & fine_dx <= rng]
        fine_dy <- fine_dy[fine_dy >= -rng & fine_dy <= rng]
        fine_df <- score_edge_grid(fine_dx, fine_dy, progress_weight = 0.08, detail = "Scoring local XY refinement")
        if (nrow(fine_df) > 0L) {
          keep_fine <- is.finite(fine_df$raw_score) | is.finite(fine_df$edge_score)
          fine_df <- fine_df[keep_fine, , drop = FALSE]
          cand_dx <- c(cand_dx, as.integer(fine_df$dX))
          cand_dy <- c(cand_dy, as.integer(fine_df$dY))
          cand_sc <- c(cand_sc, as.numeric(fine_df$raw_score))
          cand_edge_sc <- c(cand_edge_sc, as.numeric(fine_df$edge_score))
        }

        cand_df <- collapse_candidates(cand_dx, cand_dy, cand_sc, cand_edge_sc, rng_local = rng, stp_local = stp)
        best_row <- pick_best_candidate(cand_df, rng_local = rng, stp_local = stp)
        if (is.null(best_row) || nrow(best_row) == 0L || !is.finite(best_row$score[1])) {
          xh$opt_xy_candidates <- NULL
          showNotification(
            "Auto-fit XY could not find a valid translation score in this range. Try a smaller step, different MSI mode/ion, or check polygon scaling.",
            type = "warning",
            duration = 8
          )
          return()
        }
        incProgress(0.12, detail = "Exact boundary snap refinement")

        robust_z2 <- function(x) {
          x <- as.numeric(x)
          x[!is.finite(x)] <- NA_real_
          keep <- is.finite(x)
          if (!any(keep)) return(rep(0, length(x)))
          med <- suppressWarnings(stats::median(x[keep], na.rm = TRUE))
          madv <- suppressWarnings(stats::mad(x[keep], center = med, constant = 1, na.rm = TRUE))
          if (!is.finite(madv) || madv <= 1e-6) return(rep(0, length(x)))
          z <- (x - med) / madv
          z[!is.finite(z)] <- 0
          z
        }

        score_translation_exact <- function(tx_abs, ty_abs) {
          m_try <- try(polygon_mask_matrix(msi, tx = tx_abs, ty = ty_abs), silent = TRUE)
          if (inherits(m_try, "try-error") || is.null(m_try)) {
            return(list(ok = FALSE, raw_score = NA_real_, edge_score = NA_real_, local_score = -Inf, n = 0L))
          }
          n_try <- suppressWarnings(sum(m_try, na.rm = TRUE))
          if (!is.finite(n_try) || n_try < 10) {
            return(list(ok = FALSE, raw_score = NA_real_, edge_score = NA_real_, local_score = -Inf, n = 0L))
          }
          cache_try <- build_mask_bands(m_try, band_px = band_px)
          sc_raw <- score_mask_shift(cache_try, signal, dx = 0, dy = 0)
          sc_edge <- score_edge_shift(cache_try, edge_cache, dx = 0, dy = 0)
          if (!is.finite(sc_raw) && !is.finite(sc_edge)) {
            return(list(ok = FALSE, raw_score = sc_raw, edge_score = sc_edge, local_score = -Inf, n = as.integer(n_try)))
          }
          sc_local <- if (is.finite(sc_edge)) sc_edge else -2
          if (is.finite(sc_raw)) sc_local <- sc_local + 0.10 * sc_raw
          list(ok = TRUE, raw_score = sc_raw, edge_score = sc_edge, local_score = sc_local, n = as.integer(n_try))
        }

        seed_df <- cand_df[order(cand_df$score, decreasing = TRUE), , drop = FALSE]
        if (nrow(seed_df) > 12L) seed_df <- seed_df[seq_len(12L), , drop = FALSE]
        if (nrow(seed_df) == 0L) {
          seed_df <- data.frame(dX = as.numeric(best_row$dX[1]), dY = as.numeric(best_row$dY[1]), stringsAsFactors = FALSE)
        }
        seed_df$translate_x <- clamp(tx0 + as.numeric(seed_df$dX), -1000, 1000)
        seed_df$translate_y <- clamp(ty0 + as.numeric(seed_df$dY), -1000, 1000)

        exact_rows <- vector("list", nrow(seed_df))
        for (i in seq_len(nrow(seed_df))) {
          ev <- score_translation_exact(seed_df$translate_x[i], seed_df$translate_y[i])
          if (!isTRUE(ev$ok)) next
          exact_rows[[i]] <- data.frame(
            translate_x = as.numeric(seed_df$translate_x[i]),
            translate_y = as.numeric(seed_df$translate_y[i]),
            dX = as.numeric(seed_df$translate_x[i] - tx0),
            dY = as.numeric(seed_df$translate_y[i] - ty0),
            raw_score = as.numeric(ev$raw_score),
            edge_score = as.numeric(ev$edge_score),
            local_score = as.numeric(ev$local_score),
            stringsAsFactors = FALSE
          )
        }
        exact_rows <- Filter(Negate(is.null), exact_rows)
        exact_df <- if (length(exact_rows) > 0) do.call(rbind, exact_rows) else data.frame()

        if (nrow(exact_df) > 0) {
          edge_band2 <- max(1, as.integer(2L * stp))
          exact_df$dist_norm <- sqrt((exact_df$dX / max(1, rng))^2 + (exact_df$dY / max(1, rng))^2)
          exact_df$dist_norm <- pmax(0, pmin(exact_df$dist_norm, 2))
          exact_df$edge_margin <- pmin(rng - abs(exact_df$dX), rng - abs(exact_df$dY))
          exact_df$edge_norm <- pmax(0, pmin(1, exact_df$edge_margin / edge_band2))
          pen_shift2 <- 0.06 * (exact_df$dist_norm^2)
          pen_edge2 <- 0.18 * ((1 - exact_df$edge_norm)^2)
          z_raw2 <- robust_z2(exact_df$raw_score)
          z_edge2 <- robust_z2(exact_df$edge_score)
          z_local2 <- robust_z2(exact_df$local_score)
          exact_df$score_exact <- 0.50 * z_edge2 + 0.35 * z_raw2 + 0.15 * z_local2 - pen_shift2 - pen_edge2
          exact_df <- exact_df[order(exact_df$score_exact, decreasing = TRUE), , drop = FALSE]
        }

        if (nrow(exact_df) == 0L) {
          tx_seed <- clamp(tx0 + as.numeric(best_row$dX[1]), -1000, 1000)
          ty_seed <- clamp(ty0 + as.numeric(best_row$dY[1]), -1000, 1000)
          ev_seed <- score_translation_exact(tx_seed, ty_seed)
          exact_df <- data.frame(
            translate_x = tx_seed,
            translate_y = ty_seed,
            dX = tx_seed - tx0,
            dY = ty_seed - ty0,
            raw_score = as.numeric(ev_seed$raw_score),
            edge_score = as.numeric(ev_seed$edge_score),
            local_score = as.numeric(ev_seed$local_score),
            score_exact = as.numeric(ev_seed$local_score),
            stringsAsFactors = FALSE
          )
        }

        # Strong local snap around best exact seed.
        tx_cur <- as.numeric(exact_df$translate_x[1])
        ty_cur <- as.numeric(exact_df$translate_y[1])
        ev_cur <- score_translation_exact(tx_cur, ty_cur)
        sc_cur <- if (isTRUE(ev_cur$ok) && is.finite(ev_cur$local_score)) ev_cur$local_score else as.numeric(exact_df$local_score[1])
        dirs <- matrix(
          c(
            1, 0,
            -1, 0,
            0, 1,
            0, -1,
            1, 1,
            1, -1,
            -1, 1,
            -1, -1
          ),
          ncol = 2,
          byrow = TRUE
        )
        for (iter in seq_len(25L)) {
          tx_c <- clamp(tx_cur + dirs[, 1], -1000, 1000)
          ty_c <- clamp(ty_cur + dirs[, 2], -1000, 1000)
          key <- paste(tx_c, ty_c, sep = ",")
          keep <- !duplicated(key)
          tx_c <- tx_c[keep]
          ty_c <- ty_c[keep]
          best_tx <- tx_cur
          best_ty <- ty_cur
          best_sc <- sc_cur
          best_ev <- ev_cur
          for (j in seq_along(tx_c)) {
            ev_j <- score_translation_exact(tx_c[j], ty_c[j])
            if (!isTRUE(ev_j$ok) || !is.finite(ev_j$local_score)) next
            if (!is.finite(best_sc) || ev_j$local_score > (best_sc + 1e-4)) {
              best_sc <- ev_j$local_score
              best_tx <- tx_c[j]
              best_ty <- ty_c[j]
              best_ev <- ev_j
            }
          }
          if (isTRUE(all.equal(best_tx, tx_cur)) && isTRUE(all.equal(best_ty, ty_cur))) break
          tx_cur <- best_tx
          ty_cur <- best_ty
          sc_cur <- best_sc
          ev_cur <- best_ev
        }

        # Basin-jump scan (1 px) to recover from nearby local minima
        # that can remain offset by ~10-30 pixels in one axis.
        jump_span <- as.integer(min(80L, max(20L, round(rng * 0.8))))
        if (jump_span >= 2L) {
          best_tx <- tx_cur
          best_ty <- ty_cur
          best_sc <- sc_cur
          best_ev <- ev_cur

          # Scan Y around current X.
          for (ddy in seq.int(-jump_span, jump_span, by = 1L)) {
            ty_try <- clamp(ty_cur + ddy, -1000, 1000)
            ev_try <- score_translation_exact(tx_cur, ty_try)
            if (!isTRUE(ev_try$ok) || !is.finite(ev_try$local_score)) next
            if (!is.finite(best_sc) || ev_try$local_score > (best_sc + 1e-4)) {
              best_tx <- tx_cur
              best_ty <- ty_try
              best_sc <- ev_try$local_score
              best_ev <- ev_try
            }
          }

          # Scan X around improved Y.
          for (ddx in seq.int(-jump_span, jump_span, by = 1L)) {
            tx_try <- clamp(best_tx + ddx, -1000, 1000)
            ev_try <- score_translation_exact(tx_try, best_ty)
            if (!isTRUE(ev_try$ok) || !is.finite(ev_try$local_score)) next
            if (!is.finite(best_sc) || ev_try$local_score > (best_sc + 1e-4)) {
              best_tx <- tx_try
              best_sc <- ev_try$local_score
              best_ev <- ev_try
            }
          }

          tx_cur <- best_tx
          ty_cur <- best_ty
          sc_cur <- best_sc
          ev_cur <- best_ev

          # One final short local snap around jumped point.
          for (iter in seq_len(10L)) {
            tx_c <- clamp(tx_cur + dirs[, 1], -1000, 1000)
            ty_c <- clamp(ty_cur + dirs[, 2], -1000, 1000)
            key <- paste(tx_c, ty_c, sep = ",")
            keep <- !duplicated(key)
            tx_c <- tx_c[keep]
            ty_c <- ty_c[keep]
            best_tx2 <- tx_cur
            best_ty2 <- ty_cur
            best_sc2 <- sc_cur
            best_ev2 <- ev_cur
            for (j in seq_along(tx_c)) {
              ev_j <- score_translation_exact(tx_c[j], ty_c[j])
              if (!isTRUE(ev_j$ok) || !is.finite(ev_j$local_score)) next
              if (!is.finite(best_sc2) || ev_j$local_score > (best_sc2 + 1e-4)) {
                best_sc2 <- ev_j$local_score
                best_tx2 <- tx_c[j]
                best_ty2 <- ty_c[j]
                best_ev2 <- ev_j
              }
            }
            if (isTRUE(all.equal(best_tx2, tx_cur)) && isTRUE(all.equal(best_ty2, ty_cur))) break
            tx_cur <- best_tx2
            ty_cur <- best_ty2
            sc_cur <- best_sc2
            ev_cur <- best_ev2
          }
        }

        incProgress(0.08, detail = "Applying best XY translation")
        tx_new <- clamp(tx_cur, -1000, 1000)
        ty_new <- clamp(ty_cur, -1000, 1000)
        best_dx <- tx_new - tx0
        best_dy <- ty_new - ty0
        best_score <- if (isTRUE(ev_cur$ok)) as.numeric(ev_cur$local_score) else NA_real_
        best_raw_score <- if (isTRUE(ev_cur$ok)) as.numeric(ev_cur$raw_score) else NA_real_
        best_edge_score <- if (isTRUE(ev_cur$ok)) as.numeric(ev_cur$edge_score) else NA_real_
        updateSliderInput(session, "translate_x", value = tx_new)
        updateSliderInput(session, "translate_y", value = ty_new)
        updateNumericInput(session, "translate_x_num", value = tx_new)
        updateNumericInput(session, "translate_y_num", value = ty_new)

        # Build and store top-5 candidate preview from exact refinement if available.
        if (!is.null(exact_df) && nrow(exact_df) > 0) {
          prev <- exact_df[order(exact_df$score_exact, decreasing = TRUE), , drop = FALSE]
          if (nrow(prev) > 5L) prev <- prev[seq_len(5L), , drop = FALSE]
          prev$score <- prev$score_exact
          prev$rank <- seq_len(nrow(prev))
          xh$opt_xy_candidates <- prev[, c("rank", "dX", "dY", "score", "raw_score", "edge_score", "translate_x", "translate_y"), drop = FALSE]
        } else if (!is.null(cand_df) && nrow(cand_df) > 0) {
          cand_df <- cand_df[order(cand_df$score, decreasing = TRUE), , drop = FALSE]
          if (nrow(cand_df) > 5L) cand_df <- cand_df[seq_len(5L), , drop = FALSE]
          cand_df$translate_x <- clamp(tx0 + cand_df$dX, -1000, 1000)
          cand_df$translate_y <- clamp(ty0 + cand_df$dY, -1000, 1000)
          cand_df$rank <- seq_len(nrow(cand_df))
          xh$opt_xy_candidates <- cand_df[, c("rank", "dX", "dY", "score", "raw_score", "edge_score", "translate_x", "translate_y"), drop = FALSE]
        } else {
          xh$opt_xy_candidates <- NULL
        }

        showNotification(
          sprintf(
            "Auto-fit XY applied: dX=%.2f, dY=%.2f (new X=%.2f, Y=%.2f; local=%.4f, contrast=%.4f, edge=%.4f; anchor dX=%d,dY=%d; band=%d px). Top-5 candidates are available below.",
            best_dx, best_dy, tx_new, ty_new, best_score, best_raw_score, best_edge_score, anchor_dx, anchor_dy, band_px
          ),
          type = "message",
          duration = 8
        )
      })
    }, ignoreInit = TRUE)

    observeEvent(input$map_to_pdata, {
      if (is.null(msi_data())) {
        showNotification("Load an MSI dataset first before mapping overlays to pData.", type = "warning", duration = 7)
        return()
      }
      if (is.null(input$mapping_source) || !nzchar(input$mapping_source)) {
        showNotification("Choose a mapping source (polygon or cluster image) first.", type = "warning", duration = 7)
        return()
      }

      if (identical(input$mapping_source, "cluster")) {
        if (is.null(input$cluster_overlay_upload)) {
          showNotification("Upload a cluster-color image first or switch mapping source to polygon.", type = "warning", duration = 8)
          return()
        }

        col_name <- trimws(input$cluster_pdata_col)
        validate(need(nchar(col_name) > 0, "Please enter a pData column name."))
        col_name <- make.names(col_name)

        payload <- cluster_mapping_payload()
        sampled <- payload$sampled_rgba
        validate(need(length(sampled) > 0, "No cluster pixels were sampled from overlay."))

        alpha_thresh <- clamp(as.numeric(input$cluster_alpha_threshold), 0, 1)
        a_byte <- rep(0L, length(sampled))
        has_alpha <- !is.na(sampled) & nchar(sampled) == 9
        a_byte[has_alpha] <- suppressWarnings(strtoi(substr(sampled[has_alpha], 8, 9), base = 16L))
        alpha_val <- a_byte / 255

        rgb_hex <- rep(NA_character_, length(sampled))
        rgb_hex[!is.na(sampled)] <- substr(sampled[!is.na(sampled)], 1, 7)
        rgb_hex[alpha_val <= alpha_thresh] <- NA_character_

        ref_cols <- extract_reference_cluster_colors(cluster_overlay_image(), max_colors = 128L)
        validate(need(length(ref_cols) > 0, "Could not detect cluster colors from cluster overlay image."))

        mapped_hex <- nearest_palette_map(rgb_hex, ref_cols)
        uniq_ref <- ref_cols
        cluster_labels <- paste0("cluster_", sprintf("%03d", seq_along(uniq_ref)))
        names(cluster_labels) <- uniq_ref
        mapped_lab <- cluster_labels[mapped_hex]

        obj <- msi_for_pdata()
        pd <- as.data.frame(Cardinal::pData(obj))
        validate(need(nrow(pd) == length(mapped_lab), "Length mismatch between pData and mapped clusters."))
        pd[[col_name]] <- factor(mapped_lab, levels = cluster_labels)

        Cardinal::pData(obj) <- build_position_pdata(obj, pd)
        xh$mapped_obj <- obj
        xh$mapped_column <- col_name
        xh$mapping_source <- "cluster"

        lookup <- data.frame(
          label = cluster_labels,
          color_hex = uniq_ref,
          stringsAsFactors = FALSE
        )
        counts <- table(pd[[col_name]])
        lookup$n_pixels <- as.integer(counts[lookup$label])
        lookup$n_pixels[is.na(lookup$n_pixels)] <- 0L
        xh$cluster_lookup <- lookup

        showNotification(
          sprintf("Mapped cluster overlay to pData column '%s' (%d non-NA assignments).", col_name, sum(!is.na(pd[[col_name]]))),
          type = "message",
          duration = 6
        )
        return()
      }

      if (is.null(input$polygon_file) || !nzchar(input$polygon_file$name)) {
        showNotification("Choose a polygon file first or switch mapping source to cluster image.", type = "warning", duration = 8)
        return()
      }
      if (is.null(polygon_data())) {
        showNotification("Polygon file could not be read. Check file format and retry.", type = "error", duration = 8)
        return()
      }

      col_name <- trimws(input$polygon_pdata_col)
      validate(need(nchar(col_name) > 0, "Please enter a pData column name."))
      col_name <- make.names(col_name)

      msi <- make_msi_raster()
      poly <- polygon_data()
      current_label_field <- as.character(input$polygon_label_field)[1]
      if (is.na(current_label_field) || !nzchar(current_label_field)) current_label_field <- "row_index"
      companion_base_field <- NULL
      if (identical(current_label_field, "clustered_polygon_class") && !is.null(xh$polygon_cluster_result)) {
        fields_poly <- colnames(as.data.frame(poly))
        fields_poly <- fields_poly[!fields_poly %in% attr(poly, "sf_column")]
        companion_base_field <- as.character(xh$polygon_label_field_before_cluster)[1]
        if (is.na(companion_base_field) || !nzchar(companion_base_field) || identical(companion_base_field, "clustered_polygon_class")) {
          companion_base_field <- if ("classification" %in% fields_poly) {
            "classification"
          } else if ("id" %in% fields_poly) {
            "id"
          } else {
            "row_index"
          }
        }
        poly$map_label_original <- get_polygon_labels(poly, companion_base_field)
      }
      poly$map_label <- get_polygon_labels(poly, current_label_field)
      axis_mode <- resolve_polygon_axis_mode(poly, input$polygon_axis_mode, get_histology_image_optional())
      src_dims <- get_overlay_source_dims()
      poly_t <- transform_polygon_sf(
        poly_sf = poly,
        nx = msi$nx,
        ny = msi$ny,
        scale_x = input$scale_x,
        scale_y = input$scale_y,
        translate_x = input$translate_x,
        translate_y = input$translate_y,
        rotate_deg = input$rotate_deg,
        flip_y = isTRUE(input$flip_histology_y),
        swap_xy = identical(axis_mode, "yx"),
        scale_mode = input$overlay_scale_mode,
        source_width = if (!is.null(src_dims)) src_dims$width else NA_real_,
        source_height = if (!is.null(src_dims)) src_dims$height else NA_real_
      )

      pts_df <- data.frame(
        x = msi$x_norm,
        y = msi$y_norm,
        idx = seq_along(msi$x_norm)
      )
      poly_crs <- normalize_crs(try(sf::st_crs(poly_t), silent = TRUE))
      if (is.null(poly_crs)) {
        pts_sf <- sf::st_as_sf(pts_df, coords = c("x", "y"), remove = FALSE)
      } else {
        pts_sf <- sf::st_as_sf(pts_df, coords = c("x", "y"), crs = poly_crs, remove = FALSE)
      }
      hit <- sf::st_intersects(pts_sf, poly_t)

      outside_label <- "outside_polygon"
      mapped_lab <- rep(outside_label, nrow(pts_df))
      poly_uid <- paste0("polygon_", sprintf("%05d", seq_len(nrow(poly_t))))
      mapped_uid <- rep(outside_label, nrow(pts_df))
      has_hit <- lengths(hit) > 0
      raw_labels <- as.character(poly_t$map_label)
      is_cell_poly <- is_cell_like_polygon_label(raw_labels)
      # Guard against full-canvas/non-cell polygons dominating intersections.
      poly_area <- suppressWarnings(as.numeric(sf::st_area(poly_t)))
      poly_area[!is.finite(poly_area)] <- 0
      canvas_area <- as.numeric(msi$nx) * as.numeric(msi$ny)
      huge_poly <- is.finite(poly_area) & poly_area > (0.5 * canvas_area)
      if (any(huge_poly) && any(is_cell_poly & !huge_poly)) {
        is_cell_poly <- is_cell_poly & !huge_poly
      }
      if (!any(is_cell_poly)) {
        # Fallback: if heuristics remove everything, use all polygons.
        is_cell_poly <- rep(TRUE, length(raw_labels))
      }

      hit_cell <- logical(nrow(pts_df))
      if (nrow(pts_df) > 0) {
        hit_cell <- vapply(hit, function(ix) {
          any(is_cell_poly[ix])
        }, logical(1))
      }
      if (any(has_hit)) {
        if (identical(input$polygon_overlap_rule, "all")) {
          mapped_lab[has_hit] <- vapply(hit[has_hit], function(ix) {
            ix_use <- ix[is_cell_poly[ix]]
            if (length(ix_use) == 0L) return(outside_label)
            paste(unique(poly_t$map_label[ix_use]), collapse = ";")
          }, character(1))
          mapped_uid[has_hit] <- vapply(hit[has_hit], function(ix) {
            ix_use <- ix[is_cell_poly[ix]]
            if (length(ix_use) == 0L) return(outside_label)
            paste(unique(poly_uid[ix_use]), collapse = ";")
          }, character(1))
        } else {
          mapped_lab[has_hit] <- vapply(hit[has_hit], function(ix) {
            ix_use <- ix[is_cell_poly[ix]]
            if (length(ix_use) == 0L) return(outside_label)
            as.character(poly_t$map_label[ix_use[1]])
          }, character(1))
          mapped_uid[has_hit] <- vapply(hit[has_hit], function(ix) {
            ix_use <- ix[is_cell_poly[ix]]
            if (length(ix_use) == 0L) return(outside_label)
            as.character(poly_uid[ix_use[1]])
          }, character(1))
        }
      }

      map_polygon_labels_from_vector <- function(label_vec_full, outside = "outside_polygon") {
        out <- rep(outside, nrow(pts_df))
        if (!any(has_hit)) return(out)
        if (identical(input$polygon_overlap_rule, "all")) {
          out[has_hit] <- vapply(hit[has_hit], function(ix) {
            ix_use <- ix[is_cell_poly[ix]]
            if (length(ix_use) == 0L) return(outside)
            paste(unique(as.character(label_vec_full[ix_use])), collapse = ";")
          }, character(1))
        } else {
          out[has_hit] <- vapply(hit[has_hit], function(ix) {
            ix_use <- ix[is_cell_poly[ix]]
            if (length(ix_use) == 0L) return(outside)
            as.character(label_vec_full[ix_use[1]])
          }, character(1))
        }
        out
      }

      mapped_lab_original <- NULL
      orig_col_name <- NULL
      if (!is.null(companion_base_field) && "map_label_original" %in% colnames(poly_t)) {
        mapped_lab_original <- map_polygon_labels_from_vector(poly_t$map_label_original, outside = outside_label)
        orig_col_name <- "polygon_region"
        if (identical(orig_col_name, col_name)) {
          orig_col_name <- make.names(paste0(col_name, "_original"))
        }
      }

      obj <- msi_for_pdata()
      pd <- as.data.frame(Cardinal::pData(obj))
      validate(need(nrow(pd) == length(mapped_lab), "Length mismatch between pData and polygon mapping."))
      pd[[col_name]] <- mapped_lab
      if (!is.null(mapped_lab_original) && length(mapped_lab_original) == nrow(pd) && !is.null(orig_col_name)) {
        pd[[orig_col_name]] <- mapped_lab_original
      }
      uid_col <- "polygon_cell_id"
      pd[[uid_col]] <- mapped_uid
      binary_col <- "polygon_is_cell"
      pd[[binary_col]] <- as.integer(hit_cell)

      Cardinal::pData(obj) <- build_position_pdata(obj, pd)
      xh$mapped_obj <- obj
      xh$mapped_column <- col_name
      xh$mapping_source <- "polygon"

      tt <- sort(table(mapped_lab), decreasing = TRUE)
      if (length(tt) > 0 && "" %in% names(tt)) tt <- tt[names(tt) != ""]
      lookup <- data.frame(
        label = names(tt),
        n_pixels = as.integer(tt),
        stringsAsFactors = FALSE
      )
      xh$cluster_lookup <- lookup

      if (!is.null(orig_col_name) && orig_col_name %in% names(pd)) {
        showNotification(
          sprintf(
            "Mapped polygons to cluster column '%s', original-label column '%s', '%s', and binary '%s' (inside=%d, outside=%d).",
            col_name,
            orig_col_name,
            uid_col,
            binary_col,
            sum(pd[[binary_col]] == 1, na.rm = TRUE),
            sum(pd[[binary_col]] == 0, na.rm = TRUE)
          ),
          type = "message",
          duration = 8
        )
      } else {
        showNotification(
          sprintf(
            "Mapped polygons to '%s', '%s', and binary '%s' (inside=%d, outside=%d).",
            col_name,
            uid_col,
            binary_col,
            sum(pd[[binary_col]] == 1, na.rm = TRUE),
            sum(pd[[binary_col]] == 0, na.rm = TRUE)
          ),
          type = "message",
          duration = 6
        )
      }
    })

    output$overlay_plot <- renderPlot({
      req(make_msi_raster())
      req(transformed_overlay())

      msi <- make_msi_raster()
      ov <- transformed_overlay()
      msi_label <- if (!is.null(msi$display_label) && nzchar(msi$display_label)) {
        msi$display_label
      } else if (is.finite(msi$mz_selected)) {
        sprintf("m/z %.5f", msi$mz_selected)
      } else {
        "MSI"
      }

      graphics::par(mar = c(0.5, 0.5, 1.5, 0.5))
      graphics::plot.new()
      graphics::plot.window(xlim = c(1, msi$nx), ylim = c(1, msi$ny), asp = 1, xaxs = "i", yaxs = "i")
      graphics::rasterImage(msi$raster, 1, 1, msi$nx, msi$ny, interpolate = FALSE)

      if (identical(ov$layer, "polygon")) {
        poly <- ov$polygons
        label_vec <- as.character(poly$map_label)
        label_vec[is.na(label_vec) | trimws(label_vec) == ""] <- "polygon"
        uniq_labels <- unique(label_vec)
        use_palette <- isTRUE(input$polygon_color_by_label)
        base_poly_col <- safe_color(input$polygon_outline_color, "#73FFFF")
        poly_lwd <- suppressWarnings(as.numeric(input$polygon_linewidth))
        if (!is.finite(poly_lwd) || poly_lwd <= 0) poly_lwd <- 1
        if (use_palette) {
          pal <- get_discrete_palette(length(uniq_labels), input$cluster_palette)
          names(pal) <- uniq_labels
        }

        for (lab in uniq_labels) {
          ii <- which(label_vec == lab)
          if (length(ii) == 0) next
          border_col <- if (use_palette) pal[lab] else base_poly_col
          graphics::plot(
            sf::st_geometry(poly[ii, , drop = FALSE]),
            add = TRUE,
            col = NA,
            border = border_col,
            lwd = poly_lwd,
            axes = FALSE,
            reset = FALSE
          )
        }
        graphics::title(main = sprintf("Polygon Overlay on MSI (%s)", msi_label))
      } else {
        cx <- (msi$nx / 2) + input$translate_x
        cy <- (msi$ny / 2) - input$translate_y
        xleft <- cx - ov$width / 2
        xright <- cx + ov$width / 2
        ybottom <- cy - ov$height / 2
        ytop <- cy + ov$height / 2
        graphics::rasterImage(ov$raster, xleft, ybottom, xright, ytop, interpolate = TRUE)
        if (identical(ov$layer, "cluster")) {
          graphics::title(main = sprintf("Cluster Overlay on MSI (%s)", msi_label))
        } else {
          graphics::title(main = sprintf("Histology Overlay on MSI (%s)", msi_label))
        }
      }
      graphics::box()
      rec <- try(grDevices::recordPlot(), silent = TRUE)
      if (!inherits(rec, "try-error")) {
        overlay_last_recorded(rec)
      }
    })

    output$download_overlay_pdf <- downloadHandler(
      filename = function() {
        paste0("histology_overlay_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
      },
      content = function(file) {
        rec <- overlay_last_recorded()
        if (is.null(rec)) {
          showNotification("No overlay plot available yet. Render the overlay first, then download PDF.", type = "warning", duration = 5)
          stop("No rendered overlay plot available for PDF export.")
        }

        width_px <- 1400
        height_px <- 900
        grDevices::pdf(file = file, width = width_px / 72, height = height_px / 72, onefile = TRUE)
        on.exit(grDevices::dev.off(), add = TRUE)
        grDevices::replayPlot(rec)
      }
    )

    output$overlay_info <- renderPrint({
      req(make_msi_raster())
      req(transformed_overlay())
      msi <- make_msi_raster()
      ov <- transformed_overlay()
      poly_corr <- list(fx = 1, fy = 1, source = "disabled")

      list(
        working_directory = setup_values()[["wd"]],
        overlay_shown = input$overlay_layer,
        msi_canvas = sprintf("%d x %d", msi$nx, msi$ny),
        msi_display_mode = msi$mode,
        msi_display_label = msi$display_label,
        mz_selected = if (is.finite(msi$mz_selected)) round(msi$mz_selected, 6) else NA_real_,
        mz_index = if (is.finite(msi$mz_index)) as.integer(msi$mz_index) else NA_integer_,
        rgb_mz = if (!is.null(msi$rgb_mz) && length(msi$rgb_mz) > 0) paste(round(msi$rgb_mz, 6), collapse = ", ") else NA_character_,
        overlay_pdata_field = if (!is.null(msi$pdata_field)) msi$pdata_field else NA_character_,
        alpha_used = if (is.finite(ov$alpha_used)) round(ov$alpha_used, 3) else NA_real_,
        alpha_histology = input$histology_alpha,
        alpha_cluster = input$cluster_alpha,
        overlay_scale_mode = input$overlay_scale_mode,
        histology_um_per_px = input$histology_um_per_px,
        msi_um_per_px = input$msi_um_per_px,
        histology_resample_factor = input$histology_resample_factor,
        effective_histology_um_per_px = if (is.finite(input$histology_um_per_px) && is.finite(input$histology_resample_factor) && input$histology_resample_factor > 0) input$histology_um_per_px * input$histology_resample_factor else NA_real_,
        suggested_scale_ratio = if (is.finite(input$histology_um_per_px) && is.finite(input$histology_resample_factor) && is.finite(input$msi_um_per_px) && input$msi_um_per_px > 0) (input$histology_um_per_px * input$histology_resample_factor) / input$msi_um_per_px else NA_real_,
        scale_correction_x = if (is.finite(poly_corr$fx)) poly_corr$fx else 1,
        scale_correction_y = if (is.finite(poly_corr$fy)) poly_corr$fy else 1,
        scale_correction_mode = "disabled (physical scaling only)",
        source_image_dim = if (!is.null(ov$source_width) && !is.null(ov$source_height)) paste0(ov$source_width, " x ", ov$source_height) else NA_character_,
        transformed_overlay_dim = if (!is.null(ov$width) && !is.null(ov$height)) paste0(ov$width, " x ", ov$height) else NA_character_,
        suggested_export_factor_fit_x = if (!is.null(ov$source_width) && is.finite(ov$source_width) && ov$source_width > 0 && is.finite(input$histology_um_per_px) && input$histology_um_per_px > 0 && is.finite(input$msi_um_per_px) && input$msi_um_per_px > 0) {
          msi$nx / (ov$source_width * (input$histology_um_per_px / input$msi_um_per_px))
        } else NA_real_,
        suggested_export_factor_fit_y = if (!is.null(ov$source_height) && is.finite(ov$source_height) && ov$source_height > 0 && is.finite(input$histology_um_per_px) && input$histology_um_per_px > 0 && is.finite(input$msi_um_per_px) && input$msi_um_per_px > 0) {
          msi$ny / (ov$source_height * (input$histology_um_per_px / input$msi_um_per_px))
        } else NA_real_,
        msi_palette = input$msi_palette,
        rgb_render_mode = input$rgb_render_mode,
        rgb_bg_cutoff = input$rgb_bg_cutoff,
        enhance_contrast = isTRUE(input$enhance_contrast),
        gaussian_smooth = isTRUE(input$gaussian_smooth),
        gaussian_sigma = input$gaussian_sigma,
        polygon_n = if (!is.null(ov$polygons)) nrow(ov$polygons) else NULL,
        polygon_axis_mode = if (!is.null(ov$axis_mode)) ov$axis_mode else NA_character_,
        polygon_source_dim = if (!is.null(ov$polygon_source_dim)) ov$polygon_source_dim else NA_character_,
        polygon_color = input$polygon_outline_color,
        polygon_linewidth = input$polygon_linewidth,
        polygon_color_by_label = isTRUE(input$polygon_color_by_label),
        mapping_source = xh$mapping_source,
        mapped_pdata_column = xh$mapped_column
      )
    })

    output$pdata_field_ui <- renderUI({
      obj <- msi_for_pdata()
      cols <- colnames(as.data.frame(Cardinal::pData(obj)))
      if (is.null(cols) || length(cols) == 0) return(NULL)

      default_col <- if ("histo_cluster" %in% cols) "histo_cluster" else cols[1]
      selectInput(ns("pdata_field"), "Visualize pData field", choices = cols, selected = default_col)
    })

    output$pdata_plot <- renderPlot({
      req(msi_for_pdata())
      req(input$pdata_field)

      obj <- msi_for_pdata()
      pd <- as.data.frame(Cardinal::pData(obj))
      validate(need(input$pdata_field %in% names(pd), "Selected pData field is not available."))

      cd <- as.data.frame(Cardinal::coord(obj))
      x <- as.integer(cd$x - min(cd$x, na.rm = TRUE) + 1L)
      y <- as.integer(cd$y - min(cd$y, na.rm = TRUE) + 1L)
      y_plot <- max(y, na.rm = TRUE) - y + 1L

      v <- pd[[input$pdata_field]]
      df <- data.frame(x = x, y = y_plot, value = v)
      hide_leg <- isTRUE(input$hide_pdata_legend)

      if (is.numeric(v)) {
        p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = value)) +
          ggplot2::geom_raster() +
          ggplot2::scale_fill_gradientn(colors = get_msi_palette(input$msi_palette), na.value = "transparent") +
          ggplot2::coord_equal() +
          ggplot2::labs(fill = input$pdata_field, x = "x", y = "y", title = paste("pData:", input$pdata_field)) +
          ggplot2::theme_minimal()
        if (hide_leg) p <- p + ggplot2::theme(legend.position = "none")
        p
      } else {
        vf <- as.factor(v)
        lev <- levels(vf)
        pal <- get_discrete_palette(length(lev), input$cluster_palette)
        names(pal) <- lev
        p <- ggplot2::ggplot(data.frame(x = x, y = y_plot, value = vf), ggplot2::aes(x = x, y = y, fill = value)) +
          ggplot2::geom_raster() +
          ggplot2::scale_fill_manual(values = pal, na.value = "transparent") +
          ggplot2::coord_equal() +
          ggplot2::labs(fill = input$pdata_field, x = "x", y = "y", title = paste("pData:", input$pdata_field)) +
          ggplot2::theme_minimal()
        if (hide_leg) p <- p + ggplot2::theme(legend.position = "none")
        p
      }
    })

    output$pdata_preview <- DT::renderDataTable({
      obj <- msi_for_pdata()
      pd <- as.data.frame(Cardinal::pData(obj))
      DT::datatable(utils::head(pd, 200), options = list(scrollX = TRUE, pageLength = 8))
    })

    output$cluster_lookup_table <- DT::renderDataTable({
      req(xh$cluster_lookup)
      DT::datatable(xh$cluster_lookup, options = list(pageLength = 8), rownames = FALSE)
    })

    output$download_registration_params <- downloadHandler(
      filename = function() {
        paste0("registration_params_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
      },
      content = function(file) {
        p <- current_registration_params()
        df <- data.frame(
          parameter = names(p),
          value = vapply(p, function(x) as.character(x)[1], character(1)),
          stringsAsFactors = FALSE
        )
        utils::write.table(df, file = file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )

    read_registration_params <- function(path) {
      dat <- try(utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE), silent = TRUE)
      if (!inherits(dat, "try-error") && all(c("parameter", "value") %in% colnames(dat))) {
        return(as.list(setNames(as.character(dat$value), as.character(dat$parameter))))
      }

      dat <- try(utils::read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE), silent = TRUE)
      if (!inherits(dat, "try-error") && all(c("parameter", "value") %in% colnames(dat))) {
        return(as.list(setNames(as.character(dat$value), as.character(dat$parameter))))
      }

      lines <- try(readLines(path, warn = FALSE), silent = TRUE)
      if (inherits(lines, "try-error")) return(NULL)
      kv <- strsplit(lines, "=", fixed = TRUE)
      ok <- lengths(kv) >= 2
      if (!any(ok)) return(NULL)
      keys <- trimws(vapply(kv[ok], `[`, character(1), 1))
      vals <- trimws(vapply(kv[ok], function(x) paste(x[-1], collapse = "="), character(1)))
      if (length(keys) == 0) return(NULL)
      as.list(setNames(vals, keys))
    }

    observeEvent(input$registration_params_upload, {
      req(input$registration_params_upload)
      kv <- read_registration_params(input$registration_params_upload$datapath)
      validate(need(!is.null(kv) && length(kv) > 0, "Could not parse registration parameter file."))
      kv <- as.list(kv)

      if ("msi_plot_mode" %in% names(kv) && !is.null(kv[["msi_plot_mode"]])) {
        mode <- tolower(trimws(as.character(kv[["msi_plot_mode"]])[1]))
        if (mode %in% c("mz", "rgb", "pdata")) {
          updateSelectInput(session, "msi_plot_mode", selected = mode)
          xh$restore_msi_mode <- mode
        }
      }
      if ("mz_value" %in% names(kv) && !is.null(kv[["mz_value"]])) {
        v <- suppressWarnings(as.numeric(kv[["mz_value"]]))
        if (is.finite(v)) {
          # Value-based restore is authoritative; index is fallback-only.
          xh$restore_mz_value <- v
          xh$restore_mz_select <- NULL
        }
      }
      if ("mz_select" %in% names(kv) && !is.null(kv[["mz_select"]]) && is.null(xh$restore_mz_value)) {
        v <- trimws(as.character(kv[["mz_select"]])[1])
        if (!is.na(v) && nzchar(v)) xh$restore_mz_select <- v
      }
      if ("rgb_mz_values" %in% names(kv) && !is.null(kv[["rgb_mz_values"]])) {
        vals <- parse_numeric_tokens(kv[["rgb_mz_values"]])
        if (length(vals) > 0) {
          # Value-based restore is authoritative; index is fallback-only.
          xh$restore_rgb_values <- vals
          xh$restore_rgb_select <- NULL
        }
      }
      if ("rgb_mz_select" %in% names(kv) && !is.null(kv[["rgb_mz_select"]]) && is.null(xh$restore_rgb_values)) {
        toks <- trimws(unlist(strsplit(as.character(kv[["rgb_mz_select"]])[1], "[,;[:space:]]+", perl = TRUE)))
        toks <- toks[nzchar(toks)]
        if (length(toks) > 0) xh$restore_rgb_select <- toks
      }
      if ("rgb_auto_n" %in% names(kv) && !is.null(kv[["rgb_auto_n"]])) {
        v <- suppressWarnings(as.integer(kv[["rgb_auto_n"]]))
        if (is.finite(v) && v %in% c(2L, 3L)) xh$restore_rgb_auto_n <- v
      }
      if ("overlay_pdata_field" %in% names(kv) && !is.null(kv[["overlay_pdata_field"]])) {
        v <- trimws(as.character(kv[["overlay_pdata_field"]])[1])
        if (!is.na(v) && nzchar(v)) {
          xh$restore_overlay_pdata_field <- v
          updateSelectInput(session, "overlay_pdata_field", selected = v)
        }
      }

      if (!is.null(kv[["overlay_scale_mode"]])) {
        mode <- tolower(trimws(kv[["overlay_scale_mode"]]))
        if (mode %in% c("absolute", "fit")) updateSelectInput(session, "overlay_scale_mode", selected = mode)
      }

      if (!is.null(kv[["histology_um_per_px"]])) {
        v <- suppressWarnings(as.numeric(kv[["histology_um_per_px"]]))
        if (is.finite(v) && v > 0) updateNumericInput(session, "histology_um_per_px", value = v)
      }
      if (!is.null(kv[["msi_um_per_px"]])) {
        v <- suppressWarnings(as.numeric(kv[["msi_um_per_px"]]))
        if (is.finite(v) && v > 0) updateNumericInput(session, "msi_um_per_px", value = v)
      }
      if (!is.null(kv[["histology_resample_factor"]])) {
        v <- suppressWarnings(as.numeric(kv[["histology_resample_factor"]]))
        if (is.finite(v) && v > 0) updateNumericInput(session, "histology_resample_factor", value = v)
      }

      if (!is.null(kv[["scale_x"]])) {
        v <- suppressWarnings(as.numeric(kv[["scale_x"]]))
        if (is.finite(v)) {
          updateSliderInput(session, "scale_x", value = clamp(v, 0.001, 50))
          updateNumericInput(session, "scale_x_num", value = clamp(v, 0.001, 50))
        }
      }
      if (!is.null(kv[["scale_y"]])) {
        v <- suppressWarnings(as.numeric(kv[["scale_y"]]))
        if (is.finite(v)) {
          updateSliderInput(session, "scale_y", value = clamp(v, 0.001, 50))
          updateNumericInput(session, "scale_y_num", value = clamp(v, 0.001, 50))
        }
      }
      if (!is.null(kv[["rotate_deg"]])) {
        v <- suppressWarnings(as.numeric(kv[["rotate_deg"]]))
        if (is.finite(v)) {
          updateSliderInput(session, "rotate_deg", value = clamp(v, -180, 180))
          updateNumericInput(session, "rotate_deg_num", value = clamp(v, -180, 180))
        }
      }
      if (!is.null(kv[["translate_x"]])) {
        v <- suppressWarnings(as.numeric(kv[["translate_x"]]))
        if (is.finite(v)) {
          updateSliderInput(session, "translate_x", value = clamp(v, -1000, 1000))
          updateNumericInput(session, "translate_x_num", value = clamp(v, -1000, 1000))
        }
      }
      if (!is.null(kv[["translate_y"]])) {
        v <- suppressWarnings(as.numeric(kv[["translate_y"]]))
        if (is.finite(v)) {
          updateSliderInput(session, "translate_y", value = clamp(v, -1000, 1000))
          updateNumericInput(session, "translate_y_num", value = clamp(v, -1000, 1000))
        }
      }
      if (!is.null(kv[["histology_alpha"]])) {
        v <- suppressWarnings(as.numeric(kv[["histology_alpha"]]))
        if (is.finite(v)) {
          updateSliderInput(session, "histology_alpha", value = clamp(v, 0, 1))
          updateNumericInput(session, "histology_alpha_num", value = clamp(v, 0, 1))
        }
      }
      if (!is.null(kv[["cluster_alpha"]])) {
        v <- suppressWarnings(as.numeric(kv[["cluster_alpha"]]))
        if (is.finite(v)) {
          updateSliderInput(session, "cluster_alpha", value = clamp(v, 0, 1))
          updateNumericInput(session, "cluster_alpha_num", value = clamp(v, 0, 1))
        }
      }
      if (!is.null(kv[["gaussian_sigma"]])) {
        v <- suppressWarnings(as.numeric(kv[["gaussian_sigma"]]))
        if (is.finite(v) && v > 0) updateNumericInput(session, "gaussian_sigma", value = v)
      }

      if (!is.null(kv[["flip_histology_y"]])) updateCheckboxInput(session, "flip_histology_y", value = to_bool(kv[["flip_histology_y"]]))
      if (!is.null(kv[["polygon_color_by_label"]])) updateCheckboxInput(session, "polygon_color_by_label", value = to_bool(kv[["polygon_color_by_label"]]))
      if (!is.null(kv[["enhance_contrast"]])) updateCheckboxInput(session, "enhance_contrast", value = to_bool(kv[["enhance_contrast"]], TRUE))
      if (!is.null(kv[["gaussian_smooth"]])) updateCheckboxInput(session, "gaussian_smooth", value = to_bool(kv[["gaussian_smooth"]], TRUE))

      if (!is.null(kv[["polygon_axis_mode"]])) {
        v <- tolower(trimws(kv[["polygon_axis_mode"]]))
        if (v %in% c("auto", "xy", "yx")) updateSelectInput(session, "polygon_axis_mode", selected = v)
      }
      if (!is.null(kv[["polygon_outline_color"]])) {
        col_val <- safe_color(kv[["polygon_outline_color"]], "#73FFFF")
        if (requireNamespace("colourpicker", quietly = TRUE)) {
          colourpicker::updateColourInput(session, "polygon_outline_color", value = col_val)
        } else {
          updateTextInput(session, "polygon_outline_color", value = col_val)
        }
      }
      if (!is.null(kv[["polygon_linewidth"]])) {
        v <- suppressWarnings(as.numeric(kv[["polygon_linewidth"]]))
        if (is.finite(v)) updateSliderInput(session, "polygon_linewidth", value = clamp(v, 0.2, 6))
      }
      if (!is.null(kv[["intensity_transform"]])) {
        v <- as.character(kv[["intensity_transform"]])
        if (v %in% c("none", "sqrt", "log1p", "asinh")) updateSelectInput(session, "intensity_transform", selected = v)
      }
      if ("rgb_render_mode" %in% names(kv) && !is.null(kv[["rgb_render_mode"]])) {
        v <- tolower(trimws(as.character(kv[["rgb_render_mode"]])))
        if (v %in% c("additive", "dominant")) updateSelectInput(session, "rgb_render_mode", selected = v)
      }
            if ("rgb_bg_cutoff" %in% names(kv) && !is.null(kv[["rgb_bg_cutoff"]])) {
              v <- suppressWarnings(as.numeric(kv[["rgb_bg_cutoff"]]))
              if (is.finite(v)) updateSliderInput(session, "rgb_bg_cutoff", value = clamp(v, 0, 1))
            }
            if ("optimize_edge_band" %in% names(kv) && !is.null(kv[["optimize_edge_band"]])) {
              v <- suppressWarnings(as.numeric(kv[["optimize_edge_band"]]))
              if (is.finite(v)) updateNumericInput(session, "optimize_edge_band", value = clamp(v, 1, 20))
            }
            if ("stat_fit_outside_mode" %in% names(kv) && !is.null(kv[["stat_fit_outside_mode"]])) {
              v <- tolower(trimws(as.character(kv[["stat_fit_outside_mode"]])[1]))
              if (v %in% c("bbox", "local", "global")) updateSelectInput(session, "stat_fit_outside_mode", selected = v)
            }
            if ("stat_fit_metric_mode" %in% names(kv) && !is.null(kv[["stat_fit_metric_mode"]])) {
              v <- tolower(trimws(as.character(kv[["stat_fit_metric_mode"]])[1]))
              if (v %in% c("inside_outside", "polygon_cluster_groups")) updateSelectInput(session, "stat_fit_metric_mode", selected = v)
            }
            if ("stat_fit_group_field" %in% names(kv) && !is.null(kv[["stat_fit_group_field"]])) {
              v <- as.character(kv[["stat_fit_group_field"]])[1]
              if (!is.na(v) && nzchar(v)) updateSelectInput(session, "stat_fit_group_field", selected = v)
            }
            if ("stat_fit_objective" %in% names(kv) && !is.null(kv[["stat_fit_objective"]])) {
              v <- tolower(trimws(as.character(kv[["stat_fit_objective"]])[1]))
              if (v %in% c("min", "max")) updateSelectInput(session, "stat_fit_objective", selected = v)
            }
            if ("stat_fit_bbox_pad" %in% names(kv) && !is.null(kv[["stat_fit_bbox_pad"]])) {
              v <- suppressWarnings(as.numeric(kv[["stat_fit_bbox_pad"]]))
              if (is.finite(v)) updateNumericInput(session, "stat_fit_bbox_pad", value = clamp(v, 0, 500))
            }
            if ("stat_fit_top_n" %in% names(kv) && !is.null(kv[["stat_fit_top_n"]])) {
              v <- suppressWarnings(as.numeric(kv[["stat_fit_top_n"]]))
              if (is.finite(v)) updateNumericInput(session, "stat_fit_top_n", value = clamp(v, 1, 100))
            }
            if ("show_fit_info" %in% names(kv) && !is.null(kv[["show_fit_info"]])) {
              updateCheckboxInput(session, "show_fit_info", value = to_bool(kv[["show_fit_info"]], TRUE))
            }
            if (!is.null(kv[["msi_palette"]])) updateSelectInput(session, "msi_palette", selected = as.character(kv[["msi_palette"]]))
            if (!is.null(kv[["cluster_palette"]])) updateSelectInput(session, "cluster_palette", selected = as.character(kv[["cluster_palette"]]))

      # Re-apply ion selectors after the UI has re-rendered (e.g., when mode switches to RGB).
      session$onFlushed(function() {
        obj_now <- try(msi_data(), silent = TRUE)
        if (!inherits(obj_now, "try-error") && !is.null(obj_now)) {
          refresh_mz_ion_inputs(obj_now)
        }
      }, once = TRUE)

      showNotification("Loaded registration parameters.", type = "message", duration = 5)
    })

    observeEvent(input$set_scale_from_resolution, {
      apply_resolution_scale(show_message = TRUE)
    })

    observeEvent(input$overlay_scale_mode, {
      req(input$overlay_scale_mode)
      if (identical(input$overlay_scale_mode, "absolute")) {
        apply_resolution_scale(show_message = FALSE)
      }
    }, ignoreInit = FALSE)

    observeEvent(input$histology_resample_factor, {
      req(input$histology_resample_factor)
      if (identical(input$overlay_scale_mode, "absolute")) {
        apply_resolution_scale(show_message = FALSE)
      }
    }, ignoreInit = TRUE)

    nudge_translation <- function(dx = 0, dy = 0) {
      tx <- suppressWarnings(as.numeric(input$translate_x))
      ty <- suppressWarnings(as.numeric(input$translate_y))
      if (!is.finite(tx)) tx <- 0
      if (!is.finite(ty)) ty <- 0
      nx <- clamp(tx + dx, -1000, 1000)
      ny <- clamp(ty + dy, -1000, 1000)
      updateSliderInput(session, "translate_x", value = nx)
      updateSliderInput(session, "translate_y", value = ny)
      updateNumericInput(session, "translate_x_num", value = nx)
      updateNumericInput(session, "translate_y_num", value = ny)
    }

    move_step <- reactive({
      if (isTRUE(input$fine_move)) 1 else 5
    })

    observeEvent(input$move_left, {
      nudge_translation(dx = -move_step(), dy = 0)
    }, ignoreInit = TRUE)

    observeEvent(input$move_right, {
      nudge_translation(dx = move_step(), dy = 0)
    }, ignoreInit = TRUE)

    observeEvent(input$move_up, {
      # Positive translate_y moves overlay down in current coordinate convention.
      nudge_translation(dx = 0, dy = -move_step())
    }, ignoreInit = TRUE)

    observeEvent(input$move_down, {
      nudge_translation(dx = 0, dy = move_step())
    }, ignoreInit = TRUE)

    observeEvent(input$save_mapped_imzml, {
      req(input$save_mapped_imzml)
      if (is.null(xh$mapped_obj)) {
        showNotification("No mapped MSI object available. Run 'Map selected source to pData' first.", type = "warning", duration = 7)
        return()
      }

      volumes <- c(wd = setup_values()[["wd"]], home = fs::path_home())
      shinyFiles::shinyFileSave(input, "save_mapped_imzml", roots = volumes, session = session)
      save_path <- shinyFiles::parseSavePath(volumes, input$save_mapped_imzml)
      if (nrow(save_path) == 0) return(NULL)

      filen <- as.character(save_path$datapath)
      # Keep legacy Cardinal folder naming format (no .imzML suffix in folder name).
      filen <- sub("\\.imzML$", "", filen, ignore.case = TRUE)

      obj <- xh$mapped_obj
      n_feat <- suppressWarnings(as.integer(try(nrow(obj), silent = TRUE)))
      n_pix <- suppressWarnings(as.integer(try(ncol(obj), silent = TRUE)))
      est_bytes <- NA_real_
      if (is.finite(n_feat) && is.finite(n_pix) && n_feat > 0 && n_pix > 0) {
        # Rough lower-bound estimate for one double-precision intensity array.
        est_bytes <- as.numeric(n_feat) * as.numeric(n_pix) * 8
      }
      if (is.finite(est_bytes)) {
        message(sprintf(
          "[HistologyIntegration] Saving mapped imzML: features=%d pixels=%d rough intensity bytes=%.2f GB",
          n_feat, n_pix, est_bytes / (1024^3)
        ))
      } else {
        message("[HistologyIntegration] Saving mapped imzML: could not estimate output size.")
      }

      app_chunks <- suppressWarnings(as.integer(try(setup_values()[["chunks"]], silent = TRUE)))
      if (!is.finite(app_chunks) || app_chunks < 1L) app_chunks <- 20L
      write_chunks <- if (is.finite(n_pix) && n_pix > 0) {
        as.integer(max(app_chunks, min(500L, ceiling(n_pix / 2000))))
      } else {
        as.integer(max(app_chunks, 100L))
      }
      if (!is.finite(write_chunks) || write_chunks < 1L) write_chunks <- 100L

      showNotification(
        sprintf("Saving mapped imzML (serial write, %d chunks). This may take time.", write_chunks),
        type = "message",
        duration = 8
      )

      old_chunks <- try(Cardinal::getCardinalNChunks(), silent = TRUE)
      old_bp <- try(Cardinal::getCardinalBPPARAM(), silent = TRUE)
      on.exit({
        # Restore app defaults if available; otherwise restore previous values if we captured them.
        bp_restore <- try(setup_values()[["par_mode"]], silent = TRUE)
        if (!inherits(bp_restore, "try-error") && !is.null(bp_restore)) {
          try(Cardinal::setCardinalBPPARAM(bp_restore), silent = TRUE)
        } else if (!inherits(old_bp, "try-error") && !is.null(old_bp)) {
          try(Cardinal::setCardinalBPPARAM(old_bp), silent = TRUE)
        }
        chunks_restore <- suppressWarnings(as.integer(try(setup_values()[["chunks"]], silent = TRUE)))
        if (is.finite(chunks_restore) && chunks_restore > 0L) {
          try(Cardinal::setCardinalNChunks(chunks_restore), silent = TRUE)
        } else if (!inherits(old_chunks, "try-error")) {
          try(Cardinal::setCardinalNChunks(old_chunks), silent = TRUE)
        }
      }, add = TRUE)

      if (requireNamespace("BiocParallel", quietly = TRUE)) {
        try(Cardinal::setCardinalBPPARAM(BiocParallel::SerialParam()), silent = TRUE)
      }
      try(Cardinal::setCardinalNChunks(write_chunks), silent = TRUE)

      write_ok <- try(writeImzML(obj, filen), silent = TRUE)
      if (inherits(write_ok, "try-error")) {
        err_txt <- gsub("\\s+", " ", as.character(write_ok))
        message("[HistologyIntegration] writeImzML failed for mapped object: ", err_txt)
        showNotification(
          paste0(
            "Mapped imzML save failed. Serial write was used, so likely causes are disk space or a partial/locked output folder. ",
            "Check free space on the target drive (often several GB required) and delete any partial output folder before retrying."
          ),
          type = "error",
          duration = 12
        )
        return()
      }

      rds_ok <- try(saveRDS(obj, paste0(filen, ".rds")), silent = TRUE)
      if (inherits(rds_ok, "try-error")) {
        message("[HistologyIntegration] saveRDS companion file failed: ", as.character(rds_ok))
        showNotification("imzML saved, but companion .rds save failed.", type = "warning", duration = 8)
        return()
      }

      showNotification("Saved mapped imzML and companion .rds file.", type = "message", duration = 6)
    })

    observeEvent(input$histology_alpha, {
      if (is.null(input$histology_alpha_num) || !isTRUE(all.equal(input$histology_alpha, input$histology_alpha_num))) {
        updateNumericInput(session, "histology_alpha_num", value = input$histology_alpha)
      }
    }, ignoreInit = TRUE)

    observeEvent(input$cluster_alpha, {
      if (is.null(input$cluster_alpha_num) || !isTRUE(all.equal(input$cluster_alpha, input$cluster_alpha_num))) {
        updateNumericInput(session, "cluster_alpha_num", value = input$cluster_alpha)
      }
    }, ignoreInit = TRUE)

    observeEvent(input$scale_x, {
      if (is.null(input$scale_x_num) || !isTRUE(all.equal(input$scale_x, input$scale_x_num))) {
        updateNumericInput(session, "scale_x_num", value = input$scale_x)
      }
    }, ignoreInit = TRUE)

    observeEvent(input$scale_y, {
      if (is.null(input$scale_y_num) || !isTRUE(all.equal(input$scale_y, input$scale_y_num))) {
        updateNumericInput(session, "scale_y_num", value = input$scale_y)
      }
    }, ignoreInit = TRUE)

    observeEvent(input$rotate_deg, {
      if (is.null(input$rotate_deg_num) || !isTRUE(all.equal(input$rotate_deg, input$rotate_deg_num))) {
        updateNumericInput(session, "rotate_deg_num", value = input$rotate_deg)
      }
    }, ignoreInit = TRUE)

    observeEvent(input$translate_x, {
      if (is.null(input$translate_x_num) || !isTRUE(all.equal(input$translate_x, input$translate_x_num))) {
        updateNumericInput(session, "translate_x_num", value = input$translate_x)
      }
    }, ignoreInit = TRUE)

    observeEvent(input$translate_y, {
      if (is.null(input$translate_y_num) || !isTRUE(all.equal(input$translate_y, input$translate_y_num))) {
        updateNumericInput(session, "translate_y_num", value = input$translate_y)
      }
    }, ignoreInit = TRUE)

    observeEvent(input$histology_alpha_num, {
      req(input$histology_alpha_num)
      val <- clamp(input$histology_alpha_num, 0, 1)
      if (!isTRUE(all.equal(val, input$histology_alpha))) {
        updateSliderInput(session, "histology_alpha", value = val)
      }
    }, ignoreInit = TRUE)

    observeEvent(input$cluster_alpha_num, {
      req(input$cluster_alpha_num)
      val <- clamp(input$cluster_alpha_num, 0, 1)
      if (!isTRUE(all.equal(val, input$cluster_alpha))) {
        updateSliderInput(session, "cluster_alpha", value = val)
      }
    }, ignoreInit = TRUE)

    observeEvent(input$scale_x_num, {
      req(input$scale_x_num)
      val <- clamp(input$scale_x_num, 0.001, 50)
      if (!isTRUE(all.equal(val, input$scale_x))) {
        updateSliderInput(session, "scale_x", value = val)
      }
    }, ignoreInit = TRUE)

    observeEvent(input$scale_y_num, {
      req(input$scale_y_num)
      val <- clamp(input$scale_y_num, 0.001, 50)
      if (!isTRUE(all.equal(val, input$scale_y))) {
        updateSliderInput(session, "scale_y", value = val)
      }
    }, ignoreInit = TRUE)

    observeEvent(input$rotate_deg_num, {
      req(input$rotate_deg_num)
      val <- clamp(input$rotate_deg_num, -180, 180)
      if (!isTRUE(all.equal(val, input$rotate_deg))) {
        updateSliderInput(session, "rotate_deg", value = val)
      }
    }, ignoreInit = TRUE)

    observeEvent(input$translate_x_num, {
      req(input$translate_x_num)
      val <- clamp(input$translate_x_num, -1000, 1000)
      if (!isTRUE(all.equal(val, input$translate_x))) {
        updateSliderInput(session, "translate_x", value = val)
      }
    }, ignoreInit = TRUE)

    observeEvent(input$translate_y_num, {
      req(input$translate_y_num)
      val <- clamp(input$translate_y_num, -1000, 1000)
      if (!isTRUE(all.equal(val, input$translate_y))) {
        updateSliderInput(session, "translate_y", value = val)
      }
    }, ignoreInit = TRUE)

    observeEvent(input$reset_transform, {
      updateSliderInput(session, "histology_alpha", value = 0.5)
      updateSliderInput(session, "cluster_alpha", value = 0.7)
      updateSliderInput(session, "scale_x", value = 1)
      updateSliderInput(session, "scale_y", value = 1)
      updateSliderInput(session, "rotate_deg", value = 0)
      updateSliderInput(session, "translate_x", value = 0)
      updateSliderInput(session, "translate_y", value = 0)

      updateNumericInput(session, "histology_alpha_num", value = 0.5)
      updateNumericInput(session, "cluster_alpha_num", value = 0.7)
      updateNumericInput(session, "scale_x_num", value = 1)
      updateNumericInput(session, "scale_y_num", value = 1)
      updateNumericInput(session, "rotate_deg_num", value = 0)
      updateNumericInput(session, "translate_x_num", value = 0)
      updateNumericInput(session, "translate_y_num", value = 0)

      updateCheckboxInput(session, "flip_histology_y", value = FALSE)
      updateSelectInput(session, "polygon_axis_mode", selected = "yx")
      updateSelectInput(session, "overlay_scale_mode", selected = "absolute")
      updateCheckboxInput(session, "fine_move", value = FALSE)
      updateCheckboxInput(session, "polygon_color_by_label", value = FALSE)
      updateCheckboxInput(session, "enhance_contrast", value = TRUE)
      updateCheckboxInput(session, "gaussian_smooth", value = TRUE)
      updateCheckboxInput(session, "show_fit_info", value = TRUE)
      updateSliderInput(session, "polygon_linewidth", value = 1)
            updateSelectInput(session, "rgb_render_mode", selected = "dominant")
            updateSliderInput(session, "rgb_bg_cutoff", value = 0.05)
            updateNumericInput(session, "optimize_edge_band", value = 4)
            updateSelectInput(session, "stat_fit_metric_mode", selected = "inside_outside")
            updateSelectInput(session, "stat_fit_outside_mode", selected = "bbox")
            updateSelectInput(session, "stat_fit_objective", selected = "min")
            updateNumericInput(session, "stat_fit_bbox_pad", value = 25)
            updateNumericInput(session, "stat_fit_top_n", value = 10)
            updateNumericInput(session, "gaussian_sigma", value = 1)
            updateNumericInput(session, "histology_resample_factor", value = 1)
            apply_resolution_scale(show_message = FALSE)
        })
  })
}
