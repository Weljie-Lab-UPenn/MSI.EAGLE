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
      histology_fit_grid = NULL,
      histology_fit_candidates = NULL,
      histology_fit_summary = NULL,
      histology_fit_signature = NULL,
      histology_fit_target_info = NULL,
      histology_fit_relation_auto = TRUE,
      histology_fit_relation_last_default = NULL,
      histology_fit_relation_updating = FALSE,
      polygon_cluster_result = NULL,
      polygon_cluster_profile_table = NULL,
      polygon_label_field_before_cluster = NULL,
      loaded_registration_reference_frame = NULL
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
        cur_field <- isolate(input$overlay_pdata_field)
        restore_field <- isolate(xh$restore_overlay_pdata_field)
        default_col <- if (!is.null(cur_field) && length(cur_field) > 0L && !is.na(cur_field[1]) && nzchar(cur_field[1]) && cur_field[1] %in% cols) {
          cur_field[1]
        } else if (!is.null(restore_field) && restore_field %in% cols) {
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
          selectInput(ns("polygon_overlap_rule"), "If multiple polygons hit one pixel", choices = c("First match" = "first", "All matches (semicolon-separated)" = "all"), selected = "first"),
          if (!is.null(input$nucleus_polygon_file) && nzchar(input$nucleus_polygon_file$name)) {
            tags$small("Nucleus GeoJSON detected: mapping will also write nucleus/cytoplasm companion columns.")
          } else {
            tags$small("Optional: upload a nucleus GeoJSON to annotate nucleus-versus-cytoplasm pixels during mapping.")
          }
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

    nucleus_polygon_data <- reactive({
      req(input$nucleus_polygon_file)
      validate(need(requireNamespace("sf", quietly = TRUE), "Package 'sf' is required for polygon mapping."))
      poly <- try(sf::st_read(input$nucleus_polygon_file$datapath, quiet = TRUE), silent = TRUE)
      validate(need(!inherits(poly, "try-error"), "Could not read nucleus polygon file. Please provide a valid GeoJSON."))
      poly
    })

    roi_polygon_data <- reactive({
      req(input$roi_polygon_file)
      validate(need(requireNamespace("sf", quietly = TRUE), "Package 'sf' is required for polygon mapping."))
      poly <- try(sf::st_read(input$roi_polygon_file$datapath, quiet = TRUE), silent = TRUE)
      validate(need(!inherits(poly, "try-error"), "Could not read ROI polygon file. Please provide a valid GeoJSON."))
      poly
    })

    histology_metadata <- reactive({
      if (is.null(input$histology_metadata_file) || is.null(input$histology_metadata_file$datapath)) return(NULL)
      read_histology_metadata_sidecar(input$histology_metadata_file$datapath)
    })

    roi_anchor_info <- reactive({
      explicit_roi <- NULL
      roi_source <- "none"
      if (!is.null(input$roi_polygon_file) && nzchar(input$roi_polygon_file$name)) {
        explicit_roi <- roi_polygon_data()
        roi_source <- "separate_roi_file"
      } else if (!is.null(input$polygon_file) && nzchar(input$polygon_file$name)) {
        poly <- polygon_data()
        roi_guess <- detect_qupath_roi_anchor(poly)
        if (!is.null(roi_guess) && !is.null(roi_guess$polygon) && nrow(roi_guess$polygon) > 0L) {
          explicit_roi <- roi_guess$polygon
          roi_source <- as.character(roi_guess$source %||% "embedded_annotation")
        }
      }
      bbox <- polygon_geometry_bbox(explicit_roi)
      list(
        polygon = explicit_roi,
        source = roi_source,
        bbox = bbox,
        present = !is.null(explicit_roi) && nrow(explicit_roi) > 0L
      )
    })

    registration_transform <- reactive({
      hist_img <- if (!is.null(input$histology_upload) && nzchar(input$histology_upload$name)) {
        try(histology_image(), silent = TRUE)
      } else {
        NULL
      }
      if (inherits(hist_img, "try-error")) hist_img <- NULL
      roi_info <- roi_anchor_info()
      overlay_frame <- overlay_frame_from_inputs(
        hist_img = hist_img,
        roi_poly = if (isTRUE(roi_info$present)) roi_info$polygon else NULL,
        metadata = histology_metadata(),
        histology_resample_factor = input$histology_resample_factor,
        polygon_file_name = if (!is.null(input$polygon_file)) input$polygon_file$name else NULL,
        histology_file_name = if (!is.null(input$histology_upload)) input$histology_upload$name else NULL
      )
      txy <- current_translate_xy()
      list(
        slide_key = overlay_frame$slide_key,
        overlay_source_origin_x = overlay_frame$origin_x,
        overlay_source_origin_y = overlay_frame$origin_y,
        overlay_source_width = overlay_frame$source_width,
        overlay_source_height = overlay_frame$source_height,
        overlay_image_width = overlay_frame$image_width,
        overlay_image_height = overlay_frame$image_height,
        overlay_downsample_factor = overlay_frame$overlay_downsample_factor,
        scale_x = suppressWarnings(as.numeric(input$scale_x)),
        scale_y = suppressWarnings(as.numeric(input$scale_y)),
        rotate_deg = suppressWarnings(as.numeric(input$rotate_deg)),
        translate_x = suppressWarnings(as.numeric(txy$tx)),
        translate_y = suppressWarnings(as.numeric(txy$ty)),
        orientation_override = list(
          flip_histology_y = isTRUE(input$flip_histology_y),
          polygon_axis_mode = tolower(trimws(as.character(input$polygon_axis_mode)[1]))
        ),
        roi_anchor_id = overlay_frame$roi_source,
        roi_anchor_bbox = overlay_frame$roi_anchor_bbox,
        roi_anchor_present = isTRUE(overlay_frame$roi_anchor_present),
        overlay_source_frame_type = overlay_frame$source_frame_type,
        frame_status = overlay_frame$frame_status
      )
    })

    output$registration_frame_status_ui <- renderUI({
      tr <- registration_transform()
      bbox <- tr$roi_anchor_bbox
      bbox_txt <- if (!is.null(bbox) && all(is.finite(bbox))) {
        sprintf("ROI bbox: [%.0f, %.0f] to [%.0f, %.0f]", bbox[["xmin"]], bbox[["ymin"]], bbox[["xmax"]], bbox[["ymax"]])
      } else {
        "ROI bbox: not available"
      }
      tags$div(
        tags$small(
          sprintf(
            "Frame status: %s | slide key: %s | ROI source: %s | downsample: %s",
            tr$frame_status %||% "legacy_best_effort",
            tr$slide_key %||% "n/a",
            tr$roi_anchor_id %||% "none",
            if (is.finite(tr$overlay_downsample_factor)) format(tr$overlay_downsample_factor, trim = TRUE) else "n/a"
          )
        ),
        tags$br(),
        tags$small(bbox_txt)
      )
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

    output$histology_fit_pdata_field_ui <- renderUI({
      req(input$histology_fit_signal_source)
      source_mode <- tolower(trimws(as.character(input$histology_fit_signal_source)[1]))
      if (identical(source_mode, "multi")) {
        cur_n <- suppressWarnings(as.integer(isolate(input$histology_fit_multi_n)))
        if (length(cur_n) != 1L || !isTRUE(is.finite(cur_n)) || cur_n < 2L) cur_n <- 12L
        return(tagList(
          numericInput(ns("histology_fit_multi_n"), "Top ions for fused target", value = clamp(cur_n, 2L, 64L), min = 2, max = 64, step = 1),
          tags$small("Build a fused MSI target from the top spatially informative ions rather than the current single-ion display.")
        ))
      }
      if (identical(source_mode, "pca")) {
        cur_component <- suppressWarnings(as.integer(isolate(input$histology_fit_pca_component)))
        cur_n_components <- suppressWarnings(as.integer(isolate(input$histology_fit_pca_n_components)))
        cur_n_features <- suppressWarnings(as.integer(isolate(input$histology_fit_pca_n_features)))
        if (length(cur_component) != 1L || !isTRUE(is.finite(cur_component)) || cur_component < 1L) cur_component <- 1L
        if (length(cur_n_components) != 1L || !isTRUE(is.finite(cur_n_components)) || cur_n_components < 1L) cur_n_components <- 1L
        if (length(cur_n_features) != 1L || !isTRUE(is.finite(cur_n_features)) || cur_n_features < 4L) cur_n_features <- 32L
        return(tagList(
          numericInput(ns("histology_fit_pca_component"), "PCA component", value = clamp(cur_component, 1L, 64L), min = 1, max = 64, step = 1),
          numericInput(ns("histology_fit_pca_n_components"), "PCs to combine", value = clamp(cur_n_components, 1L, 8L), min = 1, max = 8, step = 1),
          numericInput(ns("histology_fit_pca_n_features"), "Informative ions for PCA", value = clamp(cur_n_features, 4L, 128L), min = 4, max = 128, step = 4),
          tags$small("Build a self-contained MSI target from one or more low-order PCA components computed on a bounded set of informative ions.")
        ))
      }
      if (!identical(source_mode, "pdata")) return(NULL)
      obj <- try(msi_for_pdata(), silent = TRUE)
      if (inherits(obj, "try-error") || is.null(obj)) {
        return(tags$small("Load MSI data to choose a pData field for Histology Fit."))
      }
      pd <- try(as.data.frame(Cardinal::pData(obj)), silent = TRUE)
      if (inherits(pd, "try-error") || is.null(pd) || ncol(pd) == 0L) {
        return(tags$small("No pData fields available in current MSI object."))
      }
      cols <- colnames(pd)
      cur <- isolate(input$histology_fit_pdata_field)
      default <- if (!is.null(cur) && cur %in% cols) {
        cur
      } else if ("polygon_is_cell" %in% cols) {
        "polygon_is_cell"
      } else if ("polygon_region" %in% cols) {
        "polygon_region"
      } else {
        cols[1]
      }
      selectInput(ns("histology_fit_pdata_field"), "pData field for Histology Fit", choices = cols, selected = default)
    })

    observeEvent(list(input$histology_fit_signal_source, input$msi_plot_mode), {
      rel_default <- default_histology_fit_relation(
        source_mode = input$histology_fit_signal_source,
        msi_mode = input$msi_plot_mode
      )
      prev_default <- isolate(xh$histology_fit_relation_last_default)
      cur_rel <- resolve_histology_intensity_relation(input$histology_fit_intensity_relation)
      should_apply <- isTRUE(xh$histology_fit_relation_auto) ||
        is.null(prev_default) ||
        identical(cur_rel, prev_default)

      if (isTRUE(should_apply) && !identical(cur_rel, rel_default)) {
        xh$histology_fit_relation_updating <- TRUE
        updateSelectInput(session, "histology_fit_intensity_relation", selected = rel_default)
      }
      xh$histology_fit_relation_auto <- isTRUE(should_apply)
      xh$histology_fit_relation_last_default <- rel_default
    }, ignoreInit = FALSE)

    observeEvent(input$histology_fit_intensity_relation, {
      if (isTRUE(xh$histology_fit_relation_updating)) {
        xh$histology_fit_relation_updating <- FALSE
        return()
      }
      cur_rel <- resolve_histology_intensity_relation(input$histology_fit_intensity_relation)
      xh$histology_fit_relation_auto <- identical(cur_rel, isolate(xh$histology_fit_relation_last_default))
    }, ignoreInit = TRUE)

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
      defaults <- c("polygon::clustered_polygon_class", "polygon::polygon_cluster_class", "pdata::polygon_cluster_class")
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

    summarize_polygon_cluster_features <- function(res, top_n = 10L) {
      if (is.null(res) || is.null(res$cluster_label) || is.null(res$feature_matrix_all)) {
        return(data.frame())
      }
      feature_df <- as.data.frame(res$feature_matrix_all, stringsAsFactors = FALSE)
      if (nrow(feature_df) == 0 || ncol(feature_df) == 0) return(data.frame())
      labels_all <- as.character(res$cluster_label)
      if (length(labels_all) != nrow(feature_df)) return(data.frame())

      retained <- res$retained_all
      if (is.null(retained) || length(retained) != length(labels_all)) {
        retained <- labels_all != "unclustered_polygon"
      }
      drop_labels <- c("", "outside_polygon", "outside", "unassigned", "unclustered_polygon")
      keep <- as.logical(retained) &
        !is.na(labels_all) & nzchar(labels_all) & !(labels_all %in% drop_labels)
      keep <- as.logical(keep)
      keep[is.na(keep)] <- FALSE
      if (!any(keep)) return(data.frame())

      labels <- labels_all[keep]
      X <- feature_df[keep, , drop = FALSE]

      num_cols <- names(X)[vapply(X, function(z) is.numeric(z) || is.integer(z), logical(1))]
      if (length(num_cols) == 0) return(data.frame())
      meas_cols <- intersect(as.character(res$measurement_feature_cols), num_cols)
      use_measurements <- length(meas_cols) > 0
      feat_cols <- if (use_measurements) meas_cols else num_cols
      if (length(feat_cols) == 0) return(data.frame())

      grp_levels <- sort(unique(labels))
      if (length(grp_levels) < 2L) return(data.frame())

      top_n <- suppressWarnings(as.integer(top_n))
      if (!is.finite(top_n) || top_n < 1L) top_n <- 10L
      top_n <- as.integer(min(100L, max(1L, top_n)))

      rows <- vector("list", 0)
      for (cl in grp_levels) {
        idx_cl <- which(labels == cl)
        idx_ot <- which(labels != cl)
        if (length(idx_cl) < 2L || length(idx_ot) < 2L) next
        for (fn in feat_cols) {
          v <- suppressWarnings(as.numeric(X[[fn]]))
          x1 <- v[idx_cl]
          x0 <- v[idx_ot]
          x1 <- x1[is.finite(x1)]
          x0 <- x0[is.finite(x0)]
          n1 <- length(x1)
          n0 <- length(x0)
          if (n1 < 2L || n0 < 2L) next

          m1 <- suppressWarnings(mean(x1, na.rm = TRUE))
          m0 <- suppressWarnings(mean(x0, na.rm = TRUE))
          med1 <- suppressWarnings(stats::median(x1, na.rm = TRUE))
          med0 <- suppressWarnings(stats::median(x0, na.rm = TRUE))
          delta <- m1 - m0
          sd_all <- suppressWarnings(stats::sd(c(x1, x0), na.rm = TRUE))
          effect <- if (is.finite(sd_all) && sd_all > 0) delta / sd_all else NA_real_

          pv <- NA_real_
          tt <- try(stats::t.test(x1, x0), silent = TRUE)
          if (!inherits(tt, "try-error") && !is.null(tt$p.value)) pv <- as.numeric(tt$p.value)

          rows[[length(rows) + 1L]] <- data.frame(
            cluster_label = cl,
            feature = fn,
            n_cluster = as.integer(n1),
            n_other = as.integer(n0),
            mean_cluster = as.numeric(m1),
            mean_other = as.numeric(m0),
            median_cluster = as.numeric(med1),
            median_other = as.numeric(med0),
            delta_mean = as.numeric(delta),
            effect_size = as.numeric(effect),
            p_value = as.numeric(pv),
            direction = ifelse(is.finite(delta) && delta >= 0, "higher_in_cluster", "lower_in_cluster"),
            feature_source = if (isTRUE(use_measurements)) "measurement" else "all_numeric_features",
            stringsAsFactors = FALSE
          )
        }
      }
      if (length(rows) == 0) return(data.frame())

      out <- do.call(rbind, rows)
      out$p_adj_cluster <- NA_real_
      for (cl in unique(out$cluster_label)) {
        ii <- which(out$cluster_label == cl & is.finite(out$p_value))
        if (length(ii) > 0) {
          out$p_adj_cluster[ii] <- stats::p.adjust(out$p_value[ii], method = "BH")
        }
      }
      out$score <- abs(ifelse(is.finite(out$effect_size), out$effect_size, 0)) *
        ifelse(is.finite(out$p_adj_cluster), -log10(out$p_adj_cluster + 1e-300), 1)
      out$score[!is.finite(out$score)] <- 0

      rank_out <- vector("list", 0)
      for (cl in unique(out$cluster_label)) {
        dd <- out[out$cluster_label == cl, , drop = FALSE]
        dd <- dd[order(dd$score, abs(dd$delta_mean), decreasing = TRUE), , drop = FALSE]
        dd$rank_in_cluster <- seq_len(nrow(dd))
        dd <- dd[dd$rank_in_cluster <= top_n, , drop = FALSE]
        rank_out[[length(rank_out) + 1L]] <- dd
      }
      out_top <- do.call(rbind, rank_out)
      out_top <- out_top[order(out_top$cluster_label, out_top$rank_in_cluster), , drop = FALSE]
      rownames(out_top) <- NULL
      out_top
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
      xh$polygon_cluster_profile_table <- NULL
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

      preview_df_all <- data.frame(
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
      preview_df_all$retained[eligible] <- keep_mask_eligible
      preview_df_all$distinctiveness_score[eligible] <- distinctiveness
      preview_df_all$cluster_label_raw[eligible] <- labels_eligible_raw
      preview_df_all <- cbind(preview_df_all, geom_features[, intersect(c("geom_area", "geom_circularity", "geom_bbox_aspect"), names(geom_features)), drop = FALSE])
      preview_df <- preview_df_all[seq_len(min(500L, nrow(preview_df_all))), , drop = FALSE]

      retained_all <- rep(FALSE, n_poly)
      retained_all[eligible] <- keep_mask_eligible
      cluster_label_raw_all <- rep("unclustered_polygon", n_poly)
      cluster_label_raw_all[eligible] <- labels_eligible_raw

      feature_cols_final <- colnames(feature_df)
      xh$polygon_cluster_result <- list(
        cluster_label = cluster_label_all,
        cluster_label_raw_all = cluster_label_raw_all,
        cluster_id = cluster_num_all,
        eligible = eligible,
        retained_all = retained_all,
        matched_measurements = matched_meas,
        polygon_key = poly_key_raw,
        polygon_key_norm = poly_key_norm,
        polygon_key_field = join_field,
        csv_id_field = csv_id_col,
        csv_feature_cols = csv_feature_cols_used,
        measurement_feature_cols = csv_feature_cols_used,
        feature_cols_used = feature_cols_final,
        feature_matrix_all = feature_df,
        feature_mode = mode,
        k = k_use,
        keep_prop = keep_prop,
        n_polygons = n_poly,
        n_eligible = sum(eligible),
        n_retained = sum(retained_all & cluster_label_all != "unclustered_polygon", na.rm = TRUE),
        n_matched_measurements = sum(matched_meas),
        n_outliers_removed = 0L,
        outlier_sd_threshold = NA_real_,
        counts = counts_df,
        preview_all = preview_df_all,
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

    observeEvent(input$remove_polygon_cluster_outliers, {
      res <- xh$polygon_cluster_result
      if (is.null(res) || is.null(res$pca) || !is.data.frame(res$pca) || nrow(res$pca) < 3) {
        showNotification("Run polygon clustering first (with PCA output).", type = "warning", duration = 6)
        return()
      }
      pca_df <- res$pca
      pc_cols <- grep("^PC[0-9]+$", names(pca_df), value = TRUE)
      if (length(pc_cols) == 0) {
        showNotification("No PCA components found for outlier filtering.", type = "warning", duration = 6)
        return()
      }
      sd_thr <- suppressWarnings(as.numeric(input$polygon_cluster_outlier_sd))
      if (!is.finite(sd_thr) || sd_thr <= 0) sd_thr <- 2

      use_idx <- which((pca_df$retained %in% TRUE) & is.finite(pca_df$polygon_row))
      if (length(use_idx) < 3) {
        showNotification("Need at least 3 retained polygons to remove PCA outliers.", type = "warning", duration = 6)
        return()
      }
      pc_mat <- as.matrix(pca_df[use_idx, pc_cols, drop = FALSE])
      storage.mode(pc_mat) <- "double"
      mu <- colMeans(pc_mat, na.rm = TRUE)
      sig <- apply(pc_mat, 2, stats::sd, na.rm = TRUE)
      good <- is.finite(sig) & sig > 0
      if (!any(good)) {
        showNotification("All PCA components have zero variance; no outliers to remove.", type = "message", duration = 6)
        return()
      }
      z <- matrix(0, nrow = nrow(pc_mat), ncol = ncol(pc_mat))
      colnames(z) <- colnames(pc_mat)
      z[, good] <- sweep(sweep(pc_mat[, good, drop = FALSE], 2, mu[good], "-"), 2, sig[good], "/")
      out_local <- rowSums(abs(z[, good, drop = FALSE]) > sd_thr, na.rm = TRUE) > 0
      if (!any(out_local)) {
        showNotification(sprintf("No PCA outliers detected at %.2f SD.", sd_thr), type = "message", duration = 6)
        return()
      }

      drop_rows <- unique(as.integer(pca_df$polygon_row[use_idx][out_local]))
      drop_rows <- drop_rows[is.finite(drop_rows)]
      drop_rows <- drop_rows[drop_rows >= 1 & drop_rows <= length(res$cluster_label)]
      if (length(drop_rows) == 0) {
        showNotification("No valid polygon rows to remove.", type = "warning", duration = 6)
        return()
      }

      retained_all <- res$retained_all
      if (is.null(retained_all) || length(retained_all) != length(res$cluster_label)) {
        retained_all <- res$cluster_label != "unclustered_polygon"
      }
      retained_all[drop_rows] <- FALSE
      res$cluster_label[drop_rows] <- "unclustered_polygon"
      res$retained_all <- retained_all

      if ("polygon_row" %in% names(pca_df)) {
        map_idx <- as.integer(pca_df$polygon_row)
        keep_map <- is.finite(map_idx) & map_idx >= 1 & map_idx <= length(retained_all)
        pca_df$retained <- FALSE
        pca_df$cluster_label <- "unclustered_polygon"
        pca_df$retained[keep_map] <- retained_all[map_idx[keep_map]]
        pca_df$cluster_label[keep_map] <- res$cluster_label[map_idx[keep_map]]
        if (!is.null(res$cluster_label_raw_all)) {
          pca_df$cluster_label_raw[keep_map] <- res$cluster_label_raw_all[map_idx[keep_map]]
        }
        res$pca <- pca_df
      }

      k_use <- suppressWarnings(as.integer(res$k))
      if (!is.finite(k_use) || k_use < 1) {
        k_use <- suppressWarnings(max(res$cluster_id, na.rm = TRUE))
      }
      if (is.finite(k_use) && k_use >= 1) {
        dig <- max(2L, nchar(as.character(k_use)))
        label_levels <- if (!is.null(res$counts) && is.data.frame(res$counts) && nrow(res$counts) >= k_use) {
          as.character(res$counts$cluster_label[seq_len(k_use)])
        } else {
          paste0("cluster_", formatC(seq_len(k_use), width = dig, flag = "0"))
        }
        cluster_id <- suppressWarnings(as.integer(res$cluster_id))
        total_counts <- tabulate(cluster_id[is.finite(cluster_id)], nbins = k_use)
        kept_counts <- tabulate(cluster_id[is.finite(cluster_id) & retained_all & res$cluster_label != "unclustered_polygon"], nbins = k_use)
        res$counts <- data.frame(
          cluster_id = seq_len(k_use),
          cluster_label = label_levels,
          n_polygons = as.integer(kept_counts),
          n_total = as.integer(total_counts),
          n_filtered_out = as.integer(total_counts - kept_counts),
          keep_prop_used = rep(if (is.null(res$keep_prop)) 1 else res$keep_prop, k_use),
          stringsAsFactors = FALSE
        )
      }

      if (!is.null(res$preview_all) && is.data.frame(res$preview_all) && nrow(res$preview_all) == length(res$cluster_label)) {
        res$preview_all$retained <- retained_all
        res$preview_all$cluster_label <- res$cluster_label
        if (!is.null(res$cluster_label_raw_all) && length(res$cluster_label_raw_all) == nrow(res$preview_all)) {
          res$preview_all$cluster_label_raw <- res$cluster_label_raw_all
        }
        res$preview <- res$preview_all[seq_len(min(500L, nrow(res$preview_all))), , drop = FALSE]
      } else if (!is.null(res$preview) && is.data.frame(res$preview) && "polygon_row" %in% names(res$preview)) {
        pidx <- suppressWarnings(as.integer(res$preview$polygon_row))
        keep <- is.finite(pidx) & pidx >= 1 & pidx <= length(res$cluster_label)
        res$preview$retained[keep] <- retained_all[pidx[keep]]
        res$preview$cluster_label[keep] <- res$cluster_label[pidx[keep]]
      }

      res$n_outliers_removed <- as.integer(if (is.null(res$n_outliers_removed)) 0 else res$n_outliers_removed) + length(drop_rows)
      res$outlier_sd_threshold <- sd_thr
      res$n_retained <- sum(retained_all & res$cluster_label != "unclustered_polygon", na.rm = TRUE)
      xh$polygon_cluster_result <- res
      xh$polygon_cluster_profile_table <- NULL

      showNotification(
        sprintf(
          "Removed %d PCA outlier polygons at %.2f SD. Retained polygons: %d.",
          length(drop_rows), sd_thr, res$n_retained
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
        outlier_sd_threshold = if (!is.null(res$outlier_sd_threshold)) res$outlier_sd_threshold else NA_real_,
        n_outliers_removed = if (!is.null(res$n_outliers_removed)) res$n_outliers_removed else 0L,
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

    observeEvent(input$summarize_polygon_clusters, {
      res <- xh$polygon_cluster_result
      if (is.null(res) || is.null(res$cluster_label)) {
        showNotification("Run polygon clustering first.", type = "warning", duration = 6)
        return()
      }
      top_n <- suppressWarnings(as.integer(input$polygon_cluster_profile_top_n))
      if (!is.finite(top_n) || top_n < 1L) top_n <- 10L

      prof <- try(summarize_polygon_cluster_features(res, top_n = top_n), silent = TRUE)
      if (inherits(prof, "try-error")) {
        showNotification("Could not compute cluster feature summary.", type = "error", duration = 7)
        return()
      }
      if (is.null(prof) || nrow(prof) == 0) {
        xh$polygon_cluster_profile_table <- NULL
        showNotification("No cluster-defining feature summary available for the current clustering result.", type = "warning", duration = 7)
        return()
      }
      xh$polygon_cluster_profile_table <- prof
      src_used <- if ("feature_source" %in% names(prof)) unique(as.character(prof$feature_source)) else "unknown"
      message(sprintf(
        "[Polygon Clustering] Cluster feature summary generated: rows=%d | clusters=%d | top_n=%d | source=%s",
        nrow(prof),
        length(unique(prof$cluster_label)),
        top_n,
        src_used[1]
      ))
      if (identical(src_used[1], "all_numeric_features")) {
        showNotification(
          "No measurement feature columns were available; summary used all numeric polygon features (including geometry).",
          type = "warning",
          duration = 8
        )
      }
      showNotification(
        sprintf(
          "Cluster feature summary generated (%d rows; %d clusters).",
          nrow(prof), length(unique(prof$cluster_label))
        ),
        type = "message",
        duration = 6
      )
    }, ignoreInit = TRUE)

    output$polygon_cluster_profile_table <- DT::renderDataTable({
      req(xh$polygon_cluster_profile_table)
      tbl <- xh$polygon_cluster_profile_table
      num_cols <- intersect(
        c("mean_cluster", "mean_other", "median_cluster", "median_other", "delta_mean", "effect_size", "p_value", "p_adj_cluster", "score"),
        names(tbl)
      )
      for (nm in num_cols) {
        tbl[[nm]] <- suppressWarnings(round(as.numeric(tbl[[nm]]), 5))
      }
      DT::datatable(
        tbl,
        options = list(pageLength = 12, scrollX = TRUE),
        rownames = FALSE
      )
    })

    output$download_polygon_cluster_profile <- downloadHandler(
      filename = function() {
        paste0("polygon_cluster_feature_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
      },
      content = function(file) {
        tbl <- xh$polygon_cluster_profile_table
        if (is.null(tbl) || nrow(tbl) == 0) {
          stop("No cluster feature summary available. Click 'Summarize cluster features' first.")
        }
        utils::write.csv(tbl, file, row.names = FALSE)
      }
    )

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

    stable_subsample_indices <- function(n_total, n_keep) {
      n_total <- suppressWarnings(as.integer(n_total))
      n_keep <- suppressWarnings(as.integer(n_keep))
      if (!is.finite(n_total) || n_total <= 0L) return(integer(0))
      if (!is.finite(n_keep) || n_keep >= n_total) return(seq_len(n_total))
      n_keep <- max(1L, min(n_total, n_keep))
      idx <- unique(as.integer(ceiling(((seq_len(n_keep) - 0.5) * n_total) / n_keep)))
      idx <- idx[is.finite(idx) & idx >= 1L & idx <= n_total]
      if (length(idx) < n_keep) {
        miss <- setdiff(seq_len(n_total), idx)
        if (length(miss) > 0L) {
          fill_n <- min(length(miss), n_keep - length(idx))
          fill_idx <- unique(as.integer(ceiling(((seq_len(fill_n) - 0.5) * length(miss)) / fill_n)))
          fill_idx <- fill_idx[is.finite(fill_idx) & fill_idx >= 1L & fill_idx <= length(miss)]
          idx <- sort(unique(c(idx, miss[fill_idx])))
        }
      }
      if (length(idx) > n_keep) idx <- idx[seq_len(n_keep)]
      sort(unique(idx))
    }

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

      pix_idx <- if (sample_pix < n_pix) stable_subsample_indices(n_pix, sample_pix) else seq_len(n_pix)
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

    build_numeric_signal_matrix <- function(
      vals,
      msi_obj,
      intensity_transform = "none",
      apply_transform = TRUE,
      gaussian_smooth = FALSE,
      smooth_sigma = 1,
      enhance = FALSE
    ) {
      ny <- as.integer(msi_obj$ny)
      nx <- as.integer(msi_obj$nx)
      row_idx <- as.integer(msi_obj$row_idx)
      x_norm <- as.integer(msi_obj$x_norm)

      vals <- as.numeric(vals)
      if (isTRUE(apply_transform)) {
        vals <- transform_intensity(vals, intensity_transform)
      }

      mat <- matrix(NA_real_, nrow = ny, ncol = nx)
      mat[cbind(row_idx, x_norm)] <- vals
      valid_mat <- is.finite(mat)

      if (!is.finite(smooth_sigma) || smooth_sigma <= 0) smooth_sigma <- 1
      if (isTRUE(gaussian_smooth)) {
        mat <- gaussian_smooth_matrix(mat, smooth_sigma)
      }

      sig <- matrix(
        rescale01(as.vector(mat), enhance = isTRUE(enhance)),
        nrow = ny,
        ncol = nx
      )
      sig[!is.finite(sig)] <- 0
      list(signal = sig, valid = valid_mat)
    }

    score_spatial_signal_map <- function(signal, valid_mask = NULL) {
      if (is.null(signal) || !is.matrix(signal) || nrow(signal) < 2L || ncol(signal) < 2L) {
        return(-Inf)
      }
      sig <- suppressWarnings(matrix(as.numeric(signal), nrow = nrow(signal), ncol = ncol(signal)))
      sig[!is.finite(sig)] <- 0
      if (is.null(valid_mask)) {
        valid <- is.finite(signal)
      } else {
        valid <- matrix(as.logical(valid_mask), nrow = nrow(signal), ncol = ncol(signal))
      }
      valid[!is.finite(valid)] <- FALSE

      vals <- sig[valid & is.finite(sig)]
      if (length(vals) < 25L) return(-Inf)

      right_grad <- abs(sig[, -1, drop = FALSE] - sig[, -ncol(sig), drop = FALSE])
      down_grad <- abs(sig[-1, , drop = FALSE] - sig[-nrow(sig), , drop = FALSE])
      right_valid <- valid[, -1, drop = FALSE] & valid[, -ncol(valid), drop = FALSE]
      down_valid <- valid[-1, , drop = FALSE] & valid[-nrow(valid), , drop = FALSE]
      grad_vals <- c(right_grad[right_valid], down_grad[down_valid])
      grad_vals <- grad_vals[is.finite(grad_vals)]
      if (length(grad_vals) < 10L) return(-Inf)

      qv <- suppressWarnings(stats::quantile(vals, probs = c(0.50, 0.90, 0.99), na.rm = TRUE, names = FALSE, type = 8))
      qg <- suppressWarnings(stats::quantile(grad_vals, probs = c(0.75, 0.95), na.rm = TRUE, names = FALSE, type = 8))
      if (!all(is.finite(qv)) || !all(is.finite(qg))) return(-Inf)

      contrast_score <- (qv[2] - qv[1]) + 0.5 * (qv[3] - qv[2])
      edge_score <- qg[2] + 0.5 * qg[1]
      score <- 0.60 * edge_score + 0.40 * contrast_score
      if (!is.finite(score)) return(-Inf)
      score
    }

    score_spatial_coherence_map <- function(signal, valid_mask = NULL, sigma = 2) {
      if (is.null(signal) || !is.matrix(signal) || nrow(signal) < 2L || ncol(signal) < 2L) {
        return(-Inf)
      }
      sig <- suppressWarnings(matrix(as.numeric(signal), nrow = nrow(signal), ncol = ncol(signal)))
      sig[!is.finite(sig)] <- 0
      if (is.null(valid_mask)) {
        valid <- is.finite(signal)
      } else {
        valid <- matrix(as.logical(valid_mask), nrow = nrow(signal), ncol = ncol(signal))
      }
      valid[!is.finite(valid)] <- FALSE
      if (sum(valid, na.rm = TRUE) < 25L) return(-Inf)

      sigma <- suppressWarnings(as.numeric(sigma))
      if (!is.finite(sigma) || sigma <= 0) sigma <- 2

      sig_masked <- sig
      sig_masked[!valid] <- NA_real_
      smooth <- gaussian_smooth_matrix(sig_masked, sigma = sigma)
      keep <- valid & is.finite(smooth)
      if (sum(keep, na.rm = TRUE) < 25L) return(-Inf)

      raw_vals <- sig[keep]
      smooth_vals <- smooth[keep]
      detail_vals <- raw_vals - smooth_vals

      corr <- suppressWarnings(stats::cor(raw_vals, smooth_vals, use = "pairwise.complete.obs", method = "spearman"))
      if (!is.finite(corr)) corr <- 0
      corr <- max(0, min(1, corr))

      var_smooth <- suppressWarnings(stats::var(smooth_vals, na.rm = TRUE))
      var_detail <- suppressWarnings(stats::var(detail_vals, na.rm = TRUE))
      if (!is.finite(var_smooth) || var_smooth < 0) var_smooth <- 0
      if (!is.finite(var_detail) || var_detail < 0) var_detail <- 0
      low_freq_frac <- var_smooth / max(1e-8, var_smooth + var_detail)
      low_freq_frac <- max(0, min(1, low_freq_frac))

      qv <- suppressWarnings(stats::quantile(raw_vals, probs = c(0.50, 0.90, 0.99), na.rm = TRUE, names = FALSE, type = 8))
      contrast <- if (all(is.finite(qv))) (qv[2] - qv[1]) + 0.35 * (qv[3] - qv[2]) else 0
      if (!is.finite(contrast)) contrast <- 0

      structure_score <- score_spatial_signal_map(signal, valid_mask = valid)
      if (!is.finite(structure_score)) structure_score <- 0

      0.45 * corr + 0.30 * low_freq_frac + 0.15 * contrast + 0.10 * structure_score
    }

    suggest_histology_fusion_indices <- function(
      obj,
      msi_obj = NULL,
      n_features = 12L,
      candidate_pool = NULL,
      max_cells = 1.5e7,
      intensity_transform = "none"
    ) {
      n_features <- suppressWarnings(as.integer(n_features))
      if (!is.finite(n_features) || n_features < 2L) n_features <- 12L
      n_features <- as.integer(min(64L, max(2L, n_features)))

      mz_axis <- suppressWarnings(as.numeric(Cardinal::mz(obj)))
      n_feat <- length(mz_axis)
      n_pix <- ncol(obj)
      if (n_feat < 2L || n_pix < 2L) return(integer(0))

      candidate_pool <- suppressWarnings(as.integer(candidate_pool))
      if (!is.finite(candidate_pool) || candidate_pool < n_features) {
        candidate_pool <- max(24L, 4L * n_features)
      }
      candidate_pool <- as.integer(min(n_feat, max(n_features, candidate_pool)))

      max_cells <- suppressWarnings(as.numeric(max_cells))
      if (!is.finite(max_cells) || max_cells <= 1e6) max_cells <- 1.5e7

      sample_pix <- floor(max_cells / max(1L, n_feat))
      sample_pix <- max(200L, sample_pix)
      sample_pix <- min(n_pix, sample_pix)
      if (sample_pix < 2L) sample_pix <- min(n_pix, 2L)

      pix_idx <- if (sample_pix < n_pix) stable_subsample_indices(n_pix, sample_pix) else seq_len(n_pix)
      sp <- try(as.matrix(Cardinal::spectra(obj)[, pix_idx, drop = FALSE]), silent = TRUE)
      if (inherits(sp, "try-error") || nrow(sp) < 2L || ncol(sp) < 2L) return(integer(0))

      storage.mode(sp) <- "double"
      sp[!is.finite(sp)] <- NA_real_
      if (!identical(intensity_transform, "none")) {
        sp <- transform_intensity(sp, intensity_transform)
      }

      score <- apply(sp, 1, stats::sd, na.rm = TRUE)
      nz <- rowMeans(is.finite(sp) & (abs(sp) > 0), na.rm = TRUE)
      score[!is.finite(score)] <- 0
      nz[!is.finite(nz)] <- 0
      score <- score * (0.25 + 0.75 * pmin(1, nz * 2))

      if (length(mz_axis) >= length(score)) {
        valid_mz <- is.finite(mz_axis[seq_len(length(score))])
        score[!valid_mz] <- -Inf
      }

      ord <- order(score, decreasing = TRUE, na.last = NA)
      ord <- ord[is.finite(score[ord]) & score[ord] > 0]
      if (length(ord) == 0L) {
        ord <- which(is.finite(score))
        ord <- ord[order(score[ord], decreasing = TRUE)]
      }
      if (length(ord) == 0L) return(integer(0))
      if (length(ord) <= candidate_pool) return(ord)

      cand_k <- min(length(ord), max(candidate_pool, as.integer(ceiling(sqrt(length(ord)) * 6))))
      cand <- ord[seq_len(cand_k)]

      robust_z_local <- function(x) {
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

      if (!is.null(msi_obj) && length(cand) > 1L) {
        eval_n <- min(length(cand), max(candidate_pool, min(96L, 4L * n_features)))
        eval_idx <- cand[seq_len(eval_n)]
        sp_eval <- try(as.matrix(Cardinal::spectra(obj)[eval_idx, , drop = FALSE]), silent = TRUE)
        if (!inherits(sp_eval, "try-error") && !is.null(sp_eval)) {
          if (nrow(sp_eval) != length(eval_idx) && ncol(sp_eval) == length(eval_idx)) {
            sp_eval <- t(sp_eval)
          }
          if (nrow(sp_eval) == length(eval_idx)) {
            storage.mode(sp_eval) <- "double"
            if (!identical(intensity_transform, "none")) {
              sp_eval <- transform_intensity(sp_eval, intensity_transform)
            }
            sp_eval[!is.finite(sp_eval)] <- 0

            spatial_score <- rep(NA_real_, length(eval_idx))
            for (ii in seq_along(eval_idx)) {
              sig_payload <- build_numeric_signal_matrix(
                vals = sp_eval[ii, ],
                msi_obj = msi_obj,
                intensity_transform = "none",
                apply_transform = FALSE,
                gaussian_smooth = FALSE,
                enhance = FALSE
              )
              spatial_score[ii] <- score_spatial_coherence_map(sig_payload$signal, valid_mask = sig_payload$valid, sigma = 2)
            }

            var_z <- robust_z_local(score[eval_idx])
            spatial_z <- robust_z_local(spatial_score)
            combined <- 0.25 * var_z + 0.75 * spatial_z
            bad_spatial <- is.finite(spatial_z) & spatial_z < -0.75
            combined[bad_spatial] <- combined[bad_spatial] - 2.0
            score[eval_idx] <- combined

            ord <- order(score, decreasing = TRUE, na.last = NA)
            ord <- ord[is.finite(score[ord])]
            if (length(ord) > 0L) {
              cand <- ord[seq_len(min(length(ord), cand_k))]
            }
          }
        }
      }

      cand_mat <- sp[cand, , drop = FALSE]
      keep <- rowSums(is.finite(cand_mat)) >= 5L
      cand <- cand[keep]
      cand_mat <- cand_mat[keep, , drop = FALSE]
      if (length(cand) <= candidate_pool || nrow(cand_mat) <= 1L) {
        return(cand[seq_len(min(length(cand), candidate_pool))])
      }

      corr <- suppressWarnings(stats::cor(t(cand_mat), use = "pairwise.complete.obs"))
      if (is.null(dim(corr))) {
        return(cand[seq_len(min(length(cand), candidate_pool))])
      }
      corr[!is.finite(corr)] <- 1
      diag(corr) <- 1

      cand_score <- score[cand]
      selected <- cand[which.max(cand_score)]
      while (length(selected) < candidate_pool && length(selected) < length(cand)) {
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
      selected[seq_len(min(length(selected), candidate_pool))]
    }

    clamp <- function(x, lo, hi) {
      pmax(lo, pmin(hi, x))
    }

    current_translate_xy <- function() {
      tx_num <- suppressWarnings(as.numeric(input$translate_x_num))
      ty_num <- suppressWarnings(as.numeric(input$translate_y_num))
      tx_sl <- suppressWarnings(as.numeric(input$translate_x))
      ty_sl <- suppressWarnings(as.numeric(input$translate_y))
      tx <- if (is.finite(tx_num)) tx_num else tx_sl
      ty <- if (is.finite(ty_num)) ty_num else ty_sl
      if (!is.finite(tx)) tx <- 0
      if (!is.finite(ty)) ty <- 0
      list(tx = tx, ty = ty)
    }

    sig_chr <- function(x) {
      if (is.null(x) || length(x) == 0L) return("")
      vals <- trimws(as.character(x))
      vals <- vals[!is.na(vals) & nzchar(vals)]
      if (length(vals) == 0L) "" else paste(vals, collapse = "|")
    }

    sig_num <- function(x, digits = 6L) {
      val <- suppressWarnings(as.numeric(x)[1])
      if (!is.finite(val)) return(NA_real_)
      round(val, digits = digits)
    }

    current_histology_fit_signature <- function() {
      txy <- current_translate_xy()
      hist_name <- ""
      if (!is.null(input$histology_upload) && !is.null(input$histology_upload$name)) {
        hist_name <- sig_chr(input$histology_upload$name)
      }
      source_mode <- sig_chr(input$histology_fit_signal_source)
      sig <- list(
        translate_x = sig_num(txy$tx),
        translate_y = sig_num(txy$ty),
        rotate_deg = sig_num(input$rotate_deg),
        scale_x = sig_num(input$scale_x),
        scale_y = sig_num(input$scale_y),
        overlay_scale_mode = sig_chr(input$overlay_scale_mode),
        flip_histology_y = isTRUE(input$flip_histology_y),
        histology_feature_mode = sig_chr(input$histology_feature_mode),
        histology_fit_signal_source = source_mode,
        histology_fit_intensity_relation = sig_chr(input$histology_fit_intensity_relation),
        histology_fit_range = sig_num(input$histology_fit_range, digits = 0L),
        histology_fit_step = sig_num(input$histology_fit_step, digits = 0L),
        histology_upload = hist_name,
        intensity_transform = sig_chr(input$intensity_transform),
        gaussian_smooth = isTRUE(input$gaussian_smooth),
        gaussian_sigma = sig_num(input$gaussian_sigma, digits = 4L)
      )
      if (identical(source_mode, "pdata")) {
        sig$histology_fit_pdata_field <- sig_chr(input$histology_fit_pdata_field)
      }
      if (identical(source_mode, "multi")) {
        sig$histology_fit_multi_n <- sig_num(input$histology_fit_multi_n, digits = 0L)
      }
      if (identical(source_mode, "pca")) {
        sig$histology_fit_pca_component <- sig_num(input$histology_fit_pca_component, digits = 0L)
        sig$histology_fit_pca_n_components <- sig_num(input$histology_fit_pca_n_components, digits = 0L)
        sig$histology_fit_pca_n_features <- sig_num(input$histology_fit_pca_n_features, digits = 0L)
      }
      if (identical(source_mode, "current")) {
        mode_now <- sig_chr(input$msi_plot_mode)
        sig$msi_plot_mode <- mode_now
        if (identical(mode_now, "pdata")) {
          sig$overlay_pdata_field <- sig_chr(input$overlay_pdata_field)
        } else if (identical(mode_now, "rgb")) {
          rgb_applied <- isolate(xh$rgb_mz_applied)
          sig$rgb_mz_select <- if (!is.null(rgb_applied) && length(rgb_applied) > 0L) {
            sig_chr(rgb_applied)
          } else {
            sig_chr(input$rgb_mz_select)
          }
          sig$rgb_render_mode <- sig_chr(input$rgb_render_mode)
          sig$rgb_bg_cutoff <- sig_num(input$rgb_bg_cutoff, digits = 4L)
        } else {
          sig$mz_select <- sig_chr(input$mz_select)
        }
      }
      sig
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

    signal_matrix_to_raster <- function(signal, palette_name = NULL) {
      if (is.null(signal) || !is.matrix(signal) || nrow(signal) == 0L || ncol(signal) == 0L) {
        return(NULL)
      }
      mat <- suppressWarnings(matrix(as.numeric(signal), nrow = nrow(signal), ncol = ncol(signal)))
      valid <- is.finite(mat)
      if (!any(valid)) {
        return(as.raster(matrix("#00000000", nrow = nrow(mat), ncol = ncol(mat))))
      }
      mat01 <- matrix(rescale01(as.vector(mat), enhance = FALSE), nrow = nrow(mat), ncol = ncol(mat))
      pal <- get_msi_palette(if (is.null(palette_name)) input$msi_palette else palette_name)
      bin <- pmax(1L, pmin(256L, as.integer(mat01 * 255) + 1L))
      col_mat <- matrix(pal[bin], nrow = nrow(mat01), ncol = ncol(mat01))
      col_mat[!valid] <- "#00000000"
      as.raster(col_mat)
    }

    get_discrete_palette <- function(n, name) {
      if (n <= 0) return(character(0))
      name_chr <- trimws(if (is.null(name) || length(name) == 0) "Alphabet" else as.character(name[[1]]))
      if (tolower(name_chr) == "alphabet") {
        base_n <- max(26L, n)
        hues <- seq(15, 375, length.out = base_n + 1L)[seq_len(base_n)]
        # Use alternating luminance so adjacent categories stay distinguishable.
        l_vals <- rep(c(62, 74), length.out = base_n)
        c_vals <- rep(c(85, 70), length.out = base_n)
        pal <- grDevices::hcl(h = hues, c = c_vals, l = l_vals)
      } else {
        pal <- try(grDevices::hcl.colors(max(n, 3), name_chr), silent = TRUE)
      }
      if (inherits(pal, "try-error") || length(pal) < n) {
        pal <- grDevices::hcl.colors(max(n, 3), "Set 2")
      }
      pal[seq_len(n)]
    }

    pdata_values_are_categorical <- function(vals, max_levels = 128L) {
      vals_num <- suppressWarnings(as.numeric(vals))
      is_num <- is.numeric(vals) || is.integer(vals)
      if (is.logical(vals)) return(TRUE)
      if (!is_num) return(TRUE)
      max_levels <- suppressWarnings(as.integer(max_levels))
      if (!is.finite(max_levels) || max_levels < 2L) max_levels <- 128L
      unique_n <- length(unique(vals[!is.na(vals)]))
      integer_like <- all(is.na(vals_num) | abs(vals_num - round(vals_num)) < 1e-8)
      isTRUE(integer_like) && isTRUE(unique_n <= max_levels)
    }

    resolve_polygon_axis_mode <- function(poly_sf, mode = "auto", hist_img = NULL) {
      mode_chr <- tolower(trimws(if (is.null(mode) || length(mode) == 0) "auto" else as.character(mode[[1]])))
      if (mode_chr %in% c("xy", "yx")) return(mode_chr)
      reg_tr <- try(registration_transform(), silent = TRUE)
      if (!inherits(reg_tr, "try-error") && !is.null(reg_tr)) {
        # In the ROI/metadata-backed workflow, QuPath image exports and polygon
        # coordinates are assumed to share the same source orientation.
        if (isTRUE(reg_tr$roi_anchor_present) || identical(reg_tr$overlay_source_frame_type, "whole_slide")) {
          return("xy")
        }
      }
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

    wrap_rotation_deg <- function(deg) {
      deg <- suppressWarnings(as.numeric(deg))
      if (!is.finite(deg)) return(NA_real_)
      out <- ((deg + 180) %% 360) - 180
      if (isTRUE(all.equal(out, -180))) out <- 180
      out
    }

    current_msi_coord_frame <- function(obj = NULL) {
      if (is.null(obj)) {
        obj <- try(msi_for_pdata(), silent = TRUE)
      }
      if (inherits(obj, "try-error") || is.null(obj)) return(NULL)
      cd <- try(as.data.frame(Cardinal::coord(obj)), silent = TRUE)
      if (inherits(cd, "try-error") || !is.data.frame(cd) || nrow(cd) == 0L || !all(c("x", "y") %in% names(cd))) {
        return(NULL)
      }
      x <- suppressWarnings(as.numeric(cd$x))
      y <- suppressWarnings(as.numeric(cd$y))
      if (!any(is.finite(x)) || !any(is.finite(y))) return(NULL)
      xmin <- suppressWarnings(min(x, na.rm = TRUE))
      xmax <- suppressWarnings(max(x, na.rm = TRUE))
      ymin <- suppressWarnings(min(y, na.rm = TRUE))
      ymax <- suppressWarnings(max(y, na.rm = TRUE))
      if (!all(is.finite(c(xmin, xmax, ymin, ymax)))) return(NULL)
      nx <- suppressWarnings(as.integer(round(xmax - xmin + 1)))
      ny <- suppressWarnings(as.integer(round(ymax - ymin + 1)))
      if (!is.finite(nx) || nx < 1L || !is.finite(ny) || ny < 1L) return(NULL)
      list(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        nx = nx,
        ny = ny
      )
    }

    registration_reference_coord_frame <- function(obj = NULL) {
      current_frame <- current_msi_coord_frame(obj)
      saved_frame <- isolate(xh$loaded_registration_reference_frame)
      if (current_frame_is_saved_subset(current_frame, saved_frame)) {
        saved_frame
      } else {
        current_frame
      }
    }

    current_registration_view_coord_frame <- function(obj = NULL) {
      current_frame <- current_msi_coord_frame(obj)
      saved_frame <- isolate(xh$loaded_registration_reference_frame)
      if (!current_frame_is_saved_subset(current_frame, saved_frame)) {
        return(current_frame)
      }
      list(
        xmin = current_frame$xmin,
        xmax = current_frame$xmax,
        ymin = saved_frame$ymax - current_frame$ymax + saved_frame$ymin,
        ymax = saved_frame$ymax - current_frame$ymin + saved_frame$ymin,
        nx = current_frame$nx,
        ny = current_frame$ny
      )
    }

    registration_pixel_source_xy <- function(msi_obj) {
      x_src <- suppressWarnings(as.numeric(msi_obj$x_source))
      y_src <- suppressWarnings(as.numeric(msi_obj$y_source))
      current_frame <- current_msi_coord_frame()
      saved_frame <- isolate(xh$loaded_registration_reference_frame)
      if (current_frame_is_saved_subset(current_frame, saved_frame)) {
        y_src <- saved_frame$ymax - y_src + saved_frame$ymin
      }
      list(x = x_src, y = y_src)
    }

    current_registration_params <- reactive({
      obj <- try(msi_for_pdata(), silent = TRUE)
      mz_axis <- numeric(0)
      coord_frame <- NULL
      reg_tr <- try(registration_transform(), silent = TRUE)
      if (inherits(reg_tr, "try-error")) reg_tr <- NULL
      if (!inherits(obj, "try-error")) {
        mz_axis <- suppressWarnings(as.numeric(Cardinal::mz(obj)))
        coord_frame <- registration_reference_coord_frame(obj)
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
        registration_canonical_frame = "msi_source_coords_v1",
        registration_rotation_convention = "shared_overlay_rotation_v2",
        registration_translation_frame = "overlay_center_in_msi_source_coords_v3",
        registration_saved_overlay_layer = input$overlay_layer,
        slide_key = if (!is.null(reg_tr)) reg_tr$slide_key else NA_character_,
        overlay_source_frame_type = if (!is.null(reg_tr)) reg_tr$overlay_source_frame_type else NA_character_,
        histology_export_origin_x = if (!is.null(reg_tr)) reg_tr$overlay_source_origin_x else NA_real_,
        histology_export_origin_y = if (!is.null(reg_tr)) reg_tr$overlay_source_origin_y else NA_real_,
        histology_export_width = if (!is.null(reg_tr)) reg_tr$overlay_source_width else NA_real_,
        histology_export_height = if (!is.null(reg_tr)) reg_tr$overlay_source_height else NA_real_,
        histology_downsample_factor = if (!is.null(reg_tr)) reg_tr$overlay_downsample_factor else NA_real_,
        roi_anchor_present = if (!is.null(reg_tr)) isTRUE(reg_tr$roi_anchor_present) else FALSE,
        roi_anchor_bbox = if (!is.null(reg_tr) && !is.null(reg_tr$roi_anchor_bbox)) paste(as.numeric(reg_tr$roi_anchor_bbox), collapse = ",") else NA_character_,
        orientation_override_active = if (!is.null(reg_tr)) {
          isTRUE(reg_tr$orientation_override$flip_histology_y) ||
            (!is.null(reg_tr$orientation_override$polygon_axis_mode) && !identical(reg_tr$orientation_override$polygon_axis_mode, "auto"))
        } else {
          FALSE
        },
        frame_status = if (!is.null(reg_tr)) reg_tr$frame_status else NA_character_,
        msi_coord_frame = if (!is.null(coord_frame)) "absolute_cardinal_coords_v1" else NA_character_,
        msi_coord_xmin = if (!is.null(coord_frame)) coord_frame$xmin else NA_real_,
        msi_coord_xmax = if (!is.null(coord_frame)) coord_frame$xmax else NA_real_,
        msi_coord_ymin = if (!is.null(coord_frame)) coord_frame$ymin else NA_real_,
        msi_coord_ymax = if (!is.null(coord_frame)) coord_frame$ymax else NA_real_,
        overlay_center_source_x = if (!is.null(coord_frame)) ((coord_frame$xmin + coord_frame$xmax) / 2) + reg_tr$translate_x else NA_real_,
        overlay_center_source_y = if (!is.null(coord_frame)) ((coord_frame$ymin + coord_frame$ymax) / 2) + reg_tr$translate_y else NA_real_,
        msi_canvas_nx = if (!is.null(coord_frame)) coord_frame$nx else NA_integer_,
        msi_canvas_ny = if (!is.null(coord_frame)) coord_frame$ny else NA_integer_,
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
        histology_fit_range = input$histology_fit_range,
        histology_fit_step = input$histology_fit_step,
        histology_feature_mode = input$histology_feature_mode,
        histology_fit_signal_source = input$histology_fit_signal_source,
        histology_fit_pdata_field = input$histology_fit_pdata_field,
        histology_fit_multi_n = input$histology_fit_multi_n,
        histology_fit_pca_component = input$histology_fit_pca_component,
        histology_fit_pca_n_components = input$histology_fit_pca_n_components,
        histology_fit_pca_n_features = input$histology_fit_pca_n_features,
        histology_fit_intensity_relation = input$histology_fit_intensity_relation,
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

      x_source <- suppressWarnings(as.numeric(cd$x))
      y_source <- suppressWarnings(as.numeric(cd$y))
      xmin <- min(x_source, na.rm = TRUE)
      xmax <- max(x_source, na.rm = TRUE)
      ymin <- min(y_source, na.rm = TRUE)
      ymax <- max(y_source, na.rm = TRUE)
      x_norm <- as.integer(x_source - xmin + 1L)
      y_norm <- as.integer(y_source - ymin + 1L)
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
          use_categorical <- pdata_values_are_categorical(vals)

          if (!isTRUE(use_categorical)) {
            ras <- build_numeric_raster(vals, apply_transform = FALSE)
            sig_opt <- build_numeric_signal(vals, apply_transform = FALSE, apply_enhance = FALSE)
            return(list(
              raster = ras,
              nx = nx,
              ny = ny,
              xmin = xmin,
              xmax = xmax,
              ymin = ymin,
              ymax = ymax,
              mz_selected = NA_real_,
              mz_index = NA_integer_,
              rgb_mz = NULL,
              pdata_field = field,
              mode = "pdata",
              mode_requested = mode_req,
              display_label = paste0("pData: ", field),
              opt_signal = sig_opt,
              x_source = x_source,
              y_source = y_source,
              x_display = x_norm,
              y_display = row_idx,
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
            xmin = xmin,
            xmax = xmax,
            ymin = ymin,
            ymax = ymax,
            mz_selected = NA_real_,
            mz_index = NA_integer_,
            rgb_mz = NULL,
            pdata_field = field,
            mode = "pdata",
            mode_requested = mode_req,
            display_label = paste0("pData: ", field),
            opt_signal = sig_opt,
            x_source = x_source,
            y_source = y_source,
            x_display = x_norm,
            y_display = row_idx,
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
              xmin = xmin,
              xmax = xmax,
              ymin = ymin,
              ymax = ymax,
              mz_selected = NA_real_,
              mz_index = NA_integer_,
              rgb_mz = rgb_mz,
              pdata_field = NULL,
              mode = "rgb",
              mode_requested = mode_req,
              display_label = rgb_label,
              opt_signal = NULL,
              x_source = x_source,
              y_source = y_source,
              x_display = x_norm,
              y_display = row_idx,
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
          xmin = xmin,
          xmax = xmax,
          ymin = ymin,
          ymax = ymax,
          mz_selected = NA_real_,
          mz_index = NA_integer_,
          rgb_mz = numeric(0),
          pdata_field = NULL,
          mode = "rgb",
          mode_requested = mode_req,
          display_label = "RGB: select 2-3 channels and click Apply RGB channels",
          opt_signal = NULL,
          x_source = x_source,
          y_source = y_source,
          x_display = x_norm,
          y_display = row_idx,
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
      sig_opt <- build_numeric_signal(vals, apply_transform = TRUE, apply_enhance = FALSE)

      list(
        raster = ras,
        nx = nx,
        ny = ny,
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        mz_selected = mzv[idx],
        mz_index = idx,
        rgb_mz = NULL,
        pdata_field = NULL,
        mode = "mz",
        mode_requested = mode_req,
        display_label = sprintf("m/z %.5f", mzv[idx]),
        opt_signal = sig_opt,
        x_source = x_source,
        y_source = y_source,
        x_display = x_norm,
        y_display = row_idx,
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
      txy_base <- current_translate_xy()
      tr_spec <- registration_transform()

      if (isTRUE(effective_overlay_flip_y())) img_work <- magick::image_flip(img_work)

      info0 <- magick::image_info(img_work)[1, ]
      msi <- make_msi_raster()
      nx <- msi$nx
      ny <- msi$ny
      coord_frame <- registration_reference_coord_frame()

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
      rot_deg <- effective_overlay_rotate_deg()
      if (is.finite(rot_deg) && rot_deg != 0) {
        img_work <- magick::image_rotate(img_work, rot_deg)
      }

      ras <- image_to_raster_rgba(img_work, alpha_scale = alpha_scale)
      src_center_x <- if (!is.null(coord_frame)) ((coord_frame$xmin + coord_frame$xmax) / 2) + tr_spec$translate_x else (nx / 2) + tr_spec$translate_x
      src_center_y <- if (!is.null(coord_frame)) ((coord_frame$ymin + coord_frame$ymax) / 2) + tr_spec$translate_y else (ny / 2) + tr_spec$translate_y
      list(
        raster = ras,
        width = ncol(ras),
        height = nrow(ras),
        translate_x_base = as.numeric(txy_base$tx),
        translate_y_base = as.numeric(txy_base$ty),
        xleft_source = src_center_x - (ncol(ras) / 2),
        xright_source = src_center_x + (ncol(ras) / 2),
        ybottom_source = src_center_y - (nrow(ras) / 2),
        ytop_source = src_center_y + (nrow(ras) / 2),
        source_width = as.numeric(info0$width),
        source_height = as.numeric(info0$height),
        scale_correction_x = corr_x,
        scale_correction_y = corr_y
      )
    }

    effective_overlay_rotate_deg <- function() {
      deg <- suppressWarnings(as.numeric(input$rotate_deg))
      if (!is.finite(deg)) deg <- 0
      # Histology-family overlays and polygon exports share the same 90-degree
      # source-frame offset, so keep the user-facing rotate_deg in a shared
      # registration space and apply the source offset internally.
      deg + 90
    }

    overlay_uses_qupath_source_frame <- function(tr_spec = NULL) {
      if (is.null(tr_spec)) {
        tr_spec <- try(registration_transform(), silent = TRUE)
      }
      if (inherits(tr_spec, "try-error") || is.null(tr_spec)) return(FALSE)
      isTRUE(tr_spec$roi_anchor_present) ||
        identical(tr_spec$overlay_source_frame_type, "whole_slide") ||
        identical(tr_spec$overlay_source_frame_type, "roi_relative") ||
        identical(tr_spec$overlay_source_frame_type, "image_bbox")
    }

    effective_polygon_rotate_deg <- function() {
      deg <- suppressWarnings(as.numeric(input$rotate_deg))
      if (!is.finite(deg)) deg <- 0
      deg
    }

    effective_overlay_flip_y <- function() {
      isTRUE(input$flip_histology_y)
    }

    effective_polygon_flip_y <- function(tr_spec = NULL) {
      if (is.null(tr_spec)) {
        tr_spec <- try(registration_transform(), silent = TRUE)
      }
      base_flip <- isTRUE(overlay_uses_qupath_source_frame(tr_spec))
      user_flip <- isTRUE(input$flip_histology_y)
      xor(base_flip, user_flip)
    }

    transformed_overlay <- reactive({
      req(input$overlay_layer)
      req(make_msi_raster())
      scale_corr <- list(fx = 1, fy = 1)

      if (identical(input$overlay_layer, "polygon")) {
        req(input$polygon_file)
        req(polygon_data())

        poly <- polygon_data()
        poly$map_label <- get_polygon_labels(poly, input$polygon_label_field)
        poly_t <- transform_polygons_for_current_registration(poly, tx = input$translate_x, ty = input$translate_y)

        return(list(
          layer = "polygon",
          polygons = poly_t$display,
          polygons_source = poly_t$source,
          alpha_used = NA_real_,
          axis_mode = poly_t$axis_mode,
          polygon_source_dim = {
            tr_spec <- registration_transform()
            if (is.finite(tr_spec$overlay_source_width) && is.finite(tr_spec$overlay_source_height)) {
              paste0(tr_spec$overlay_source_width, "x", tr_spec$overlay_source_height, " (registration_frame)")
            } else {
              NA_character_
            }
          }
        ))
      }

      if (identical(input$overlay_layer, "combined")) {
        hist_out <- NULL
        poly_out <- NULL
        if (!is.null(input$histology_upload) && nzchar(input$histology_upload$name)) {
          alpha_use <- if (is.finite(input$histology_alpha)) input$histology_alpha else 0.5
          hist_out <- transform_overlay_image(histology_image(), use_point_filter = FALSE, alpha_scale = alpha_use, scale_correction = scale_corr)
          hist_out$alpha_used <- alpha_use
        }
        if (!is.null(input$polygon_file) && nzchar(input$polygon_file$name)) {
          poly <- polygon_data()
          poly$map_label <- get_polygon_labels(poly, input$polygon_label_field)
          poly_out <- transform_polygons_for_current_registration(poly, tx = input$translate_x, ty = input$translate_y)
        }
        validate(need(!is.null(hist_out) || !is.null(poly_out), "Load a histology image or polygon file for combined validation."))
        return(list(
          layer = "combined",
          histology = hist_out,
          polygons = if (!is.null(poly_out)) poly_out$display else NULL,
          polygons_source = if (!is.null(poly_out)) poly_out$source else NULL,
          axis_mode = if (!is.null(poly_out)) poly_out$axis_mode else NA_character_,
          alpha_used = if (!is.null(hist_out)) hist_out$alpha_used else NA_real_
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

    sample_transformed_overlay_to_msi <- function(msi_obj, tr, tx, ty) {
      w <- as.integer(tr$width)
      h <- as.integer(tr$height)
      validate(need(is.finite(w) && w > 0 && is.finite(h) && h > 0, "Invalid transformed overlay dimensions."))

      if (!all(is.finite(c(tr$xleft_source, tr$xright_source, tr$ybottom_source, tr$ytop_source)))) {
        stop("Transformed overlay is missing canonical source bounds.")
      }

      tx_base <- suppressWarnings(as.numeric(tr$translate_x_base))
      ty_base <- suppressWarnings(as.numeric(tr$translate_y_base))
      tx_now <- suppressWarnings(as.numeric(tx))
      ty_now <- suppressWarnings(as.numeric(ty))
      if (!is.finite(tx_base)) tx_base <- 0
      if (!is.finite(ty_base)) ty_base <- 0
      if (!is.finite(tx_now)) tx_now <- tx_base
      if (!is.finite(ty_now)) ty_now <- ty_base
      dx_shift <- tx_now - tx_base
      dy_shift <- ty_now - ty_base

      xleft <- tr$xleft_source + dx_shift
      xright <- tr$xright_source + dx_shift
      ybottom <- tr$ybottom_source + dy_shift
      ytop <- tr$ytop_source + dy_shift

      pix_src <- registration_pixel_source_xy(msi_obj)
      xpix <- suppressWarnings(as.numeric(pix_src$x))
      ypix <- suppressWarnings(as.numeric(pix_src$y))

      u <- (xpix - xleft) / (xright - xleft)
      v <- (ytop - ypix) / (ytop - ybottom)
      inside <- is.finite(u) & is.finite(v) & u >= 0 & u <= 1 & v >= 0 & v <= 1

      col_idx <- pmin(w, pmax(1L, as.integer(floor(u * (w - 1L)) + 1L)))
      row_idx <- pmin(h, pmax(1L, as.integer(floor(v * (h - 1L)) + 1L)))

      col_mat <- as.matrix(tr$raster)
      sampled <- rep(NA_character_, length(xpix))
      sampled[inside] <- col_mat[cbind(row_idx[inside], col_idx[inside])]

      list(
        sampled_rgba = sampled,
        inside = inside,
        width = w,
        height = h
      )
    }

    source_canvas_grid_df <- function(msi_obj) {
      coord_frame <- current_registration_view_coord_frame()
      validate(need(!is.null(coord_frame), "MSI source coordinate frame is unavailable."))
      x_seq <- seq(coord_frame$xmin, coord_frame$xmax, by = 1)
      y_seq <- seq(coord_frame$ymax, coord_frame$ymin, by = -1)
      expand.grid(
        x = x_seq,
        y = y_seq,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
    }

    overlay_source_bbox_to_display <- function(ov) {
      coord_frame <- current_registration_view_coord_frame()
      if (is.null(coord_frame) || is.null(ov) ||
          !all(is.finite(c(ov$xleft_source, ov$xright_source, ov$ybottom_source, ov$ytop_source)))) {
        return(NULL)
      }
      source_bbox_to_display_bbox(
        c(
          xmin = ov$xleft_source,
          xmax = ov$xright_source,
          ymin = ov$ybottom_source,
          ymax = ov$ytop_source
        ),
        coord_frame = coord_frame
      )
    }

    transform_polygons_for_current_registration <- function(poly_sf, tx, ty, axis_mode = NULL) {
      msi <- make_msi_raster()
      tr_spec <- registration_transform()
      if (is.null(axis_mode)) {
        axis_mode <- resolve_polygon_axis_mode(poly_sf, input$polygon_axis_mode, get_histology_image_optional())
      }
      polygon_rotate <- effective_polygon_rotate_deg()
      polygon_flip_y <- effective_polygon_flip_y(tr_spec)
      poly_src <- transform_polygon_sf(
        poly_sf = poly_sf,
        nx = msi$nx,
        ny = msi$ny,
        scale_x = input$scale_x,
        scale_y = input$scale_y,
        translate_x = tx,
        translate_y = ty,
        rotate_deg = polygon_rotate,
        flip_y = polygon_flip_y,
        swap_xy = identical(axis_mode, "yx"),
        scale_mode = input$overlay_scale_mode,
        source_width = tr_spec$overlay_source_width,
        source_height = tr_spec$overlay_source_height,
        source_origin_x = tr_spec$overlay_source_origin_x,
        source_origin_y = tr_spec$overlay_source_origin_y,
        downsample_factor = tr_spec$overlay_downsample_factor,
        coord_frame = registration_reference_coord_frame()
      )
      poly_display <- source_sf_to_display_sf(poly_src, current_registration_view_coord_frame())
      list(source = poly_src, display = poly_display, axis_mode = axis_mode)
    }

    decode_rgba_hex <- function(rgba_hex) {
      rgba_hex <- as.character(rgba_hex)
      n <- length(rgba_hex)
      out <- list(
        r = rep(0, n),
        g = rep(0, n),
        b = rep(0, n),
        a = rep(0, n)
      )
      keep <- !is.na(rgba_hex) & nchar(rgba_hex) >= 7
      if (!any(keep)) return(out)

      rr <- suppressWarnings(strtoi(substr(rgba_hex[keep], 2, 3), base = 16L))
      gg <- suppressWarnings(strtoi(substr(rgba_hex[keep], 4, 5), base = 16L))
      bb <- suppressWarnings(strtoi(substr(rgba_hex[keep], 6, 7), base = 16L))
      aa <- rep(255L, sum(keep))
      has_alpha <- nchar(rgba_hex[keep]) >= 9
      aa[has_alpha] <- suppressWarnings(strtoi(substr(rgba_hex[keep][has_alpha], 8, 9), base = 16L))

      rr[!is.finite(rr)] <- 0L
      gg[!is.finite(gg)] <- 0L
      bb[!is.finite(bb)] <- 0L
      aa[!is.finite(aa)] <- 255L

      out$r[keep] <- rr / 255
      out$g[keep] <- gg / 255
      out$b[keep] <- bb / 255
      out$a[keep] <- aa / 255
      out
    }

    build_histology_feature_from_sampled <- function(sampled_rgba, msi_obj, feature_mode = "hematoxylin", gaussian_smooth = FALSE, smooth_sigma = 1) {
      ny <- as.integer(msi_obj$ny)
      nx <- as.integer(msi_obj$nx)
      feat_mode <- tolower(trimws(as.character(feature_mode)[1]))
      if (!feat_mode %in% c("hematoxylin", "darkness", "purple")) feat_mode <- "hematoxylin"

      ch <- decode_rgba_hex(sampled_rgba)
      valid_vec <- !is.na(sampled_rgba) & is.finite(ch$a) & ch$a > 0
      if (!any(valid_vec)) {
        return(list(signal = matrix(0, nrow = ny, ncol = nx), valid = matrix(FALSE, nrow = ny, ncol = nx)))
      }

      lum <- 0.2126 * ch$r + 0.7152 * ch$g + 0.0722 * ch$b
      dark <- pmax(0, 1 - lum)
      purple <- pmax(0, ((ch$r + ch$b) / 2) - ch$g)

      feat <- switch(
        feat_mode,
        darkness = dark,
        purple = 0.80 * purple + 0.20 * dark,
        0.65 * dark + 0.35 * purple
      )
      feat[!valid_vec] <- 0

      mat <- matrix(0, nrow = ny, ncol = nx)
      mat[cbind(as.integer(msi_obj$row_idx), as.integer(msi_obj$x_norm))] <- feat

      valid_mat <- matrix(FALSE, nrow = ny, ncol = nx)
      valid_mat[cbind(as.integer(msi_obj$row_idx), as.integer(msi_obj$x_norm))] <- valid_vec

      if (!is.finite(smooth_sigma) || smooth_sigma <= 0) smooth_sigma <- 1
      if (isTRUE(gaussian_smooth)) {
        mat <- gaussian_smooth_matrix(mat, smooth_sigma)
      }

      sig <- matrix(
        rescale01(as.vector(mat), enhance = FALSE),
        nrow = ny,
        ncol = nx
      )
      sig[!is.finite(sig)] <- 0

      list(signal = sig, valid = valid_mat)
    }

    build_edge_mask <- function(signal, q = 0.85, valid_mask = NULL) {
      ny <- nrow(signal)
      nx <- ncol(signal)
      if (!is.finite(ny) || !is.finite(nx) || ny < 2 || nx < 2) {
        return(NULL)
      }

      sig <- suppressWarnings(matrix(as.numeric(signal), nrow = ny, ncol = nx))
      sig[!is.finite(sig)] <- 0
      if (is.null(valid_mask)) {
        valid <- is.finite(signal)
      } else {
        valid <- matrix(as.logical(valid_mask), nrow = ny, ncol = nx)
      }
      valid[!is.finite(valid)] <- FALSE

      right <- cbind(sig[, -1, drop = FALSE], sig[, nx, drop = FALSE])
      down <- rbind(sig[-1, , drop = FALSE], sig[ny, , drop = FALSE])
      right_valid <- cbind(valid[, -1, drop = FALSE], FALSE)
      down_valid <- rbind(valid[-1, , drop = FALSE], rep(FALSE, nx))

      grad <- abs(right - sig) + abs(down - sig)
      grad[!(valid & right_valid & down_valid)] <- 0
      grad[!is.finite(grad)] <- 0

      g_ok <- grad[valid & is.finite(grad) & grad > 0]
      if (length(g_ok) < 10L) return(NULL)
      thr <- suppressWarnings(stats::quantile(g_ok, probs = q, na.rm = TRUE, names = FALSE, type = 8))
      if (!is.finite(thr)) {
        thr <- suppressWarnings(mean(g_ok, na.rm = TRUE))
      }
      if (!is.finite(thr)) return(NULL)

      edge <- grad >= thr
      edge[!valid] <- FALSE
      n_edge <- suppressWarnings(sum(edge, na.rm = TRUE))
      if (!is.finite(n_edge) || n_edge < 5L) return(NULL)
      list(edge = edge, n_edge = as.integer(n_edge), threshold = thr)
    }

    normalized_mutual_information <- function(x, y, nbins = 24L) {
      x <- as.numeric(x)
      y <- as.numeric(y)
      keep <- is.finite(x) & is.finite(y)
      if (sum(keep) < 20L) return(NA_real_)

      xk <- x[keep]
      yk <- y[keep]
      xr <- range(xk, na.rm = TRUE)
      yr <- range(yk, na.rm = TRUE)
      if (!all(is.finite(xr)) || !all(is.finite(yr)) || diff(xr) <= 1e-8 || diff(yr) <= 1e-8) {
        return(NA_real_)
      }

      nbins <- suppressWarnings(as.integer(nbins))
      if (!is.finite(nbins) || nbins < 4L) nbins <- 24L
      nbins <- as.integer(min(64L, max(4L, nbins)))

      bx <- pmin(nbins, pmax(1L, as.integer(floor((xk - xr[1]) / max(1e-8, diff(xr)) * nbins) + 1L)))
      by <- pmin(nbins, pmax(1L, as.integer(floor((yk - yr[1]) / max(1e-8, diff(yr)) * nbins) + 1L)))

      joint <- table(bx, by)
      pxy <- joint / sum(joint)
      px <- rowSums(pxy)
      py <- colSums(pxy)

      nz <- which(pxy > 0, arr.ind = TRUE)
      if (nrow(nz) == 0L) return(NA_real_)
      mi <- sum(pxy[nz] * log(pxy[nz] / (px[nz[, 1]] * py[nz[, 2]])))
      hx <- -sum(px[px > 0] * log(px[px > 0]))
      hy <- -sum(py[py > 0] * log(py[py > 0]))
      den <- hx + hy
      if (!is.finite(mi) || !is.finite(den) || den <= 0) return(NA_real_)
      2 * mi / den
    }

    signed_rank_correlation <- function(x, y) {
      x <- as.numeric(x)
      y <- as.numeric(y)
      keep <- is.finite(x) & is.finite(y)
      if (sum(keep) < 20L) return(NA_real_)

      xk <- x[keep]
      yk <- y[keep]
      xr <- range(xk, na.rm = TRUE)
      yr <- range(yk, na.rm = TRUE)
      if (!all(is.finite(xr)) || !all(is.finite(yr)) || diff(xr) <= 1e-8 || diff(yr) <= 1e-8) {
        return(NA_real_)
      }

      sc <- suppressWarnings(stats::cor(xk, yk, method = "spearman"))
      if (!is.finite(sc)) return(NA_real_)
      max(-1, min(1, as.numeric(sc)))
    }

    categorical_signal_association <- function(x, labels, min_group_n = 5L) {
      x <- as.numeric(x)
      labs <- trimws(as.character(labels))
      keep <- is.finite(x) & !is.na(labs) & nzchar(labs)
      if (sum(keep) < 25L) return(NA_real_)

      xk <- x[keep]
      lk <- labs[keep]
      group_n <- table(lk)
      min_group_n <- suppressWarnings(as.integer(min_group_n))
      if (!is.finite(min_group_n) || min_group_n < 2L) min_group_n <- 5L
      keep_groups <- names(group_n[group_n >= min_group_n])
      if (length(keep_groups) < 2L) return(NA_real_)

      group_keep <- lk %in% keep_groups
      xk <- xk[group_keep]
      lk <- lk[group_keep]
      if (length(xk) < 25L) return(NA_real_)

      grand_mean <- suppressWarnings(mean(xk, na.rm = TRUE))
      ss_tot <- suppressWarnings(sum((xk - grand_mean)^2, na.rm = TRUE))
      if (!is.finite(ss_tot) || ss_tot <= 1e-10) return(NA_real_)

      split_x <- split(xk, lk, drop = TRUE)
      ss_between <- sum(vapply(split_x, function(v) {
        v <- v[is.finite(v)]
        if (length(v) == 0L) return(0)
        length(v) * (mean(v, na.rm = TRUE) - grand_mean)^2
      }, numeric(1)), na.rm = TRUE)

      eta2 <- ss_between / ss_tot
      if (!is.finite(eta2)) return(NA_real_)
      max(0, min(1, as.numeric(eta2)))
    }

    resolve_histology_intensity_relation <- function(x = NULL) {
      rel <- x
      if (is.null(rel)) rel <- input$histology_fit_intensity_relation
      rel <- tolower(trimws(as.character(rel)[1]))
      if (!rel %in% c("direct", "inverse", "either")) rel <- "either"
      rel
    }

    histology_intensity_relation_label <- function(x = NULL) {
      rel <- resolve_histology_intensity_relation(x)
      switch(
        rel,
        direct = "direct",
        inverse = "inverse",
        "either / unsigned"
      )
    }

    default_histology_fit_relation <- function(source_mode = NULL, msi_mode = NULL) {
      source_mode <- tolower(trimws(if (is.null(source_mode) || length(source_mode) == 0L) "current" else as.character(source_mode)[1]))
      if (!source_mode %in% c("current", "pdata", "multi", "pca")) source_mode <- "current"
      msi_mode <- tolower(trimws(if (is.null(msi_mode) || length(msi_mode) == 0L) "mz" else as.character(msi_mode)[1]))
      if (!msi_mode %in% c("mz", "rgb", "pdata")) msi_mode <- "mz"

      if (identical(source_mode, "pdata") || (identical(source_mode, "current") && identical(msi_mode, "pdata"))) {
        "either"
      } else {
        "inverse"
      }
    }

    score_edge_mask_against_distance <- function(edge_mask, edge_cache) {
      if (is.null(edge_mask) || is.null(edge_cache) || is.null(edge_cache$dist)) return(NA_real_)
      keep <- edge_mask & is.finite(edge_cache$dist)
      n_keep <- suppressWarnings(sum(keep, na.rm = TRUE))
      if (!is.finite(n_keep) || n_keep < 10L) return(NA_real_)
      d <- as.numeric(edge_cache$dist[keep])
      d <- d[is.finite(d)]
      if (length(d) < 10L) return(NA_real_)
      cov <- min(1, length(d) / max(25, suppressWarnings(as.integer(edge_cache$n_edge))))
      (-mean(d)) - 0.60 * (1 - cov)^2
    }

    histology_feature_payload <- function(msi_obj, tr_hist, tx, ty, feature_mode = "hematoxylin", gaussian_smooth = FALSE, smooth_sigma = 1) {
      sampled <- sample_transformed_overlay_to_msi(msi_obj = msi_obj, tr = tr_hist, tx = tx, ty = ty)
      feat <- build_histology_feature_from_sampled(
        sampled_rgba = sampled$sampled_rgba,
        msi_obj = msi_obj,
        feature_mode = feature_mode,
        gaussian_smooth = gaussian_smooth,
        smooth_sigma = smooth_sigma
      )
      feat$inside <- sampled$inside
      feat
    }

    cluster_mapping_payload <- reactive({
      req(make_msi_raster())
      req(input$cluster_overlay_upload)

      scale_corr <- list(fx = 1, fy = 1)
      tr <- transform_overlay_image(cluster_overlay_image(), use_point_filter = TRUE, alpha_scale = 1, scale_correction = scale_corr)
      msi <- make_msi_raster()
      sampled <- sample_transformed_overlay_to_msi(msi_obj = msi, tr = tr, tx = input$translate_x, ty = input$translate_y)
      list(sampled_rgba = sampled$sampled_rgba)
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

    flip_values_vertically_to_match_display <- function(vals, msi_obj) {
      vals <- as.vector(vals)
      x_src <- suppressWarnings(as.integer(msi_obj$x_source))
      y_src <- suppressWarnings(as.integer(msi_obj$y_source))
      if (length(vals) != length(x_src) || length(vals) != length(y_src)) return(vals)
      if (!all(is.finite(x_src)) || !all(is.finite(y_src))) return(vals)

      ymin <- min(y_src, na.rm = TRUE)
      ymax <- max(y_src, na.rm = TRUE)
      key_now <- paste(x_src, y_src, sep = ":")
      key_flip <- paste(x_src, ymax - y_src + ymin, sep = ":")
      m <- match(key_flip, key_now)
      ok <- is.finite(m) & m >= 1L & m <= length(vals)
      out <- vals
      out[ok] <- vals[m[ok]]
      out
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

    match_nucleus_polygons_to_cells <- function(cell_poly_t, nucleus_poly_t, cell_keep = NULL) {
      if (is.null(cell_poly_t) || nrow(cell_poly_t) == 0 || is.null(nucleus_poly_t) || nrow(nucleus_poly_t) == 0) {
        return(rep(NA_integer_, if (is.null(nucleus_poly_t)) 0L else nrow(nucleus_poly_t)))
      }

      cell_keep <- rep_len(if (is.null(cell_keep)) TRUE else as.logical(cell_keep), nrow(cell_poly_t))
      cell_keep[is.na(cell_keep)] <- FALSE
      if (!any(cell_keep)) cell_keep[] <- TRUE

      cell_idx_all <- seq_len(nrow(cell_poly_t))
      cell_idx_use <- cell_idx_all[cell_keep]
      cell_use <- cell_poly_t[cell_keep, , drop = FALSE]

      cell_area <- suppressWarnings(as.numeric(sf::st_area(cell_use)))
      cell_area[!is.finite(cell_area)] <- Inf

      cell_pts <- try(suppressWarnings(sf::st_point_on_surface(cell_use)), silent = TRUE)
      cell_xy <- try(sf::st_coordinates(cell_pts), silent = TRUE)
      nuc_pts <- try(suppressWarnings(sf::st_point_on_surface(nucleus_poly_t)), silent = TRUE)
      if (inherits(nuc_pts, "try-error")) {
        nuc_pts <- try(suppressWarnings(sf::st_centroid(nucleus_poly_t)), silent = TRUE)
      }
      nuc_xy <- try(sf::st_coordinates(nuc_pts), silent = TRUE)

      choose_best_local <- function(ix, i_nuc = NA_integer_) {
        ix <- suppressWarnings(as.integer(ix))
        ix <- ix[is.finite(ix) & ix >= 1L & ix <= nrow(cell_use)]
        if (length(ix) == 0L) return(NA_integer_)
        if (length(ix) == 1L) return(cell_idx_use[ix])

        area_ix <- cell_area[ix]
        if (any(is.finite(area_ix))) {
          best_area <- min(area_ix, na.rm = TRUE)
          ix_area <- ix[is.finite(area_ix) & abs(area_ix - best_area) <= 1e-8]
          if (length(ix_area) == 1L) return(cell_idx_use[ix_area])
          if (length(ix_area) > 0L) ix <- ix_area
        }

        if (!inherits(cell_xy, "try-error") && !inherits(nuc_xy, "try-error") &&
            nrow(cell_xy) >= nrow(cell_use) && is.finite(i_nuc) && i_nuc >= 1L && i_nuc <= nrow(nucleus_poly_t) &&
            nrow(nuc_xy) >= i_nuc) {
          dx <- cell_xy[ix, 1] - nuc_xy[i_nuc, 1]
          dy <- cell_xy[ix, 2] - nuc_xy[i_nuc, 2]
          d2 <- dx * dx + dy * dy
          d2[!is.finite(d2)] <- Inf
          return(cell_idx_use[ix[which.min(d2)]])
        }

        cell_idx_use[ix[1]]
      }

      matched <- rep(NA_integer_, nrow(nucleus_poly_t))

      hit_within <- try(sf::st_within(nuc_pts, cell_use), silent = TRUE)
      if (!inherits(hit_within, "try-error")) {
        for (i in seq_along(hit_within)) {
          matched[i] <- choose_best_local(hit_within[[i]], i_nuc = i)
        }
      }

      need_intersection <- which(!is.finite(matched))
      if (length(need_intersection) > 0L) {
        hit_intersects <- try(sf::st_intersects(nucleus_poly_t[need_intersection, , drop = FALSE], cell_use), silent = TRUE)
        if (!inherits(hit_intersects, "try-error")) {
          for (k in seq_along(need_intersection)) {
            i <- need_intersection[k]
            matched[i] <- choose_best_local(hit_intersects[[k]], i_nuc = i)
          }
        }
      }

      need_nearest <- which(!is.finite(matched))
      if (length(need_nearest) > 0L && !inherits(cell_xy, "try-error") && !inherits(nuc_xy, "try-error") &&
          nrow(cell_xy) >= nrow(cell_use) && nrow(nuc_xy) >= max(need_nearest)) {
        for (i in need_nearest) {
          dx <- cell_xy[, 1] - nuc_xy[i, 1]
          dy <- cell_xy[, 2] - nuc_xy[i, 2]
          d2 <- dx * dx + dy * dy
          d2[!is.finite(d2)] <- Inf
          matched[i] <- cell_idx_use[which.min(d2)]
        }
      }

      matched
    }

    transform_polygon_sf <- function(poly_sf, nx, ny, scale_x, scale_y, translate_x, translate_y, rotate_deg, flip_y = FALSE, swap_xy = FALSE, scale_mode = "absolute", source_width = NA_real_, source_height = NA_real_, source_origin_x = NA_real_, source_origin_y = NA_real_, downsample_factor = 1, coord_frame = NULL) {
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
      downsample_factor <- suppressWarnings(as.numeric(downsample_factor)[1])
      if (!is.finite(downsample_factor) || downsample_factor <= 0) downsample_factor <- 1
      frame_w_img <- frame_w / downsample_factor
      frame_h_img <- frame_h / downsample_factor

      mode <- tolower(trimws(as.character(scale_mode)[1]))
      if (!mode %in% c("absolute", "fit")) mode <- "absolute"
      if (identical(mode, "fit")) {
        fit_scale <- min(nx / frame_w_img, ny / frame_h_img)
        sx <- fit_scale * scale_x
        sy <- fit_scale * scale_y
      } else {
        sx <- as.numeric(scale_x)
        sy <- as.numeric(scale_y)
      }
      validate(need(is.finite(sx) && is.finite(sy) && sx > 0 && sy > 0, "Polygon scale must be positive."))

      center_frame <- coord_frame
      if (is.null(center_frame)) center_frame <- current_msi_coord_frame()
      cx <- if (!is.null(center_frame)) ((center_frame$xmin + center_frame$xmax) / 2) + translate_x else (nx / 2) + translate_x
      cy <- if (!is.null(center_frame)) ((center_frame$ymin + center_frame$ymax) / 2) + translate_y else (ny / 2) + translate_y
      theta <- rotate_deg * pi / 180
      ct <- cos(theta)
      st <- sin(theta)

      use_source_origin <- is.finite(source_origin_x) && is.finite(source_origin_y)
      in_frame <- if (has_source_dims) {
        bb_xmin >= -1 && bb_ymin >= -1 && bb_xmax <= (frame_w + 1) && bb_ymax <= (frame_h + 1)
      } else {
        FALSE
      }
      use_source_frame <- use_source_origin || (has_source_dims && in_frame)
      geom_frame_w <- if (use_source_frame) frame_w_img else (bb_w / downsample_factor)
      geom_frame_h <- if (use_source_frame) frame_h_img else (bb_h / downsample_factor)
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

        if (use_source_origin) {
          x_base <- x_src - source_origin_x
          y_base <- y_src - source_origin_y
        } else if (has_source_dims && in_frame) {
          x_base <- x_src
          y_base <- y_src
        } else {
          x_base <- x_src - bb_xmin
          y_base <- y_src - bb_ymin
        }

        x_frame <- x_base / downsample_factor
        y_frame <- y_base / downsample_factor
        x0 <- x_frame * sx
        if (isTRUE(flip_y)) {
          y0 <- (geom_frame_h - y_frame) * sy
        } else {
          y0 <- y_frame * sy
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
      poly_t <- transform_polygons_for_current_registration(poly, tx = tx, ty = ty)$source

      # For inside/outside optimization, suppress obvious non-cell and
      # overly large polygons that can make the mask nearly invariant.
      raw_labels <- as.character(poly_t$map_label)
      keep_mask <- rep(TRUE, nrow(poly_t))
      cell_like <- is_cell_like_polygon_label(raw_labels)
      if (any(cell_like, na.rm = TRUE)) {
        keep_mask <- keep_mask & cell_like
      }
      poly_area <- suppressWarnings(as.numeric(sf::st_area(poly_t)))
      poly_area[!is.finite(poly_area)] <- 0
      canvas_area <- as.numeric(msi_obj$nx) * as.numeric(msi_obj$ny)
      huge <- is.finite(poly_area) & (poly_area > (0.5 * canvas_area))
      if (any(huge) && any(keep_mask & !huge)) {
        keep_mask <- keep_mask & !huge
      }
      if (any(keep_mask)) {
        poly_t <- poly_t[keep_mask, , drop = FALSE]
      }
      validate(need(nrow(poly_t) > 0, "No valid polygons remain for inside/outside mask (after filtering non-cell/oversized polygons)."))

      grid <- source_canvas_grid_df(msi_obj)
      poly_crs <- normalize_crs(try(sf::st_crs(poly_t), silent = TRUE))
      if (is.null(poly_crs)) {
        pts_sf <- sf::st_as_sf(grid, coords = c("x", "y"), remove = FALSE)
      } else {
        pts_sf <- sf::st_as_sf(grid, coords = c("x", "y"), crs = poly_crs, remove = FALSE)
      }
      hit <- sf::st_intersects(pts_sf, poly_t)
      hit_mask <- lengths(hit) > 0
      coord_frame <- current_registration_view_coord_frame()
      row_idx <- as.integer(coord_frame$ymax - grid$y + 1L)
      col_idx <- as.integer(grid$x - coord_frame$xmin + 1L)

      mask <- matrix(FALSE, nrow = as.integer(msi_obj$ny), ncol = as.integer(msi_obj$nx))
      mask[cbind(row_idx, col_idx)] <- hit_mask
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
      poly_t <- transform_polygons_for_current_registration(poly, tx = tx, ty = ty)$source

      grid <- source_canvas_grid_df(msi_obj)
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

      coord_frame <- current_registration_view_coord_frame()
      row_idx <- as.integer(coord_frame$ymax - grid$y + 1L)
      col_idx <- as.integer(grid$x - coord_frame$xmin + 1L)
      code_mat <- matrix(0L, nrow = as.integer(msi_obj$ny), ncol = as.integer(msi_obj$nx))
      code_mat[cbind(row_idx, col_idx)] <- code_vec

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
      shift_code_to_anchor <- function(code_mat, tx_use, ty_use) {
        tx_i <- suppressWarnings(as.integer(round(tx_use)))
        ty_i <- suppressWarnings(as.integer(round(ty_use)))
        if (!is.finite(tx_i)) tx_i <- 0L
        if (!is.finite(ty_i)) ty_i <- 0L
        if (tx_i != 0L || ty_i != 0L) {
          code_mat <- shift_int(code_mat, dr = ty_i, dc = tx_i, fill = 0L)
        }
        list(code_mat = code_mat, tx = tx_i, ty = ty_i)
      }

      spec <- parse_stat_fit_group_field(group_field_sel)
      if (identical(spec$source, "pdata")) {
        out <- pdata_group_label_code_matrix(msi_obj = msi_obj, field = spec$field, exclude_labels = exclude_labels)
        anchored <- shift_code_to_anchor(out$code_mat, tx, ty)
        out$code_mat <- anchored$code_mat
        out$anchor_shift_x <- anchored$tx
        out$anchor_shift_y <- anchored$ty
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
          anchored <- shift_code_to_anchor(out$code_mat, tx, ty)
          out$code_mat <- anchored$code_mat
          out$anchor_shift_x <- anchored$tx
          out$anchor_shift_y <- anchored$ty
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

    edge_fit_signal_from_pdata_field <- function(msi_obj, field, binary_mode = c("categorical", "mask", "label_map")) {
      obj <- msi_for_pdata()
      pd <- as.data.frame(Cardinal::pData(obj))
      if (!field %in% names(pd)) {
        stop(sprintf("pData field '%s' not found for Edge Fit optimization.", field))
      }
      binary_mode <- match.arg(binary_mode)

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
      uniq_num <- sort(unique(vals_num[is.finite(vals_num)]))
      binary_like <- (is.logical(vals) || integer_like) &&
        length(uniq_num) <= 2L &&
        all(uniq_num %in% c(0, 1))
      as_categorical <- (!is_num) || (integer_like && unique_n <= 128L)

      if (binary_mode %in% c("mask", "label_map") && isTRUE(binary_like)) {
        sig_payload <- build_numeric_signal_matrix(
          vals = vals_num,
          msi_obj = msi_obj,
          intensity_transform = "none",
          apply_transform = FALSE,
          gaussian_smooth = isTRUE(input$gaussian_smooth),
          smooth_sigma = suppressWarnings(as.numeric(input$gaussian_sigma)),
          enhance = FALSE
        )
        return(list(signal = sig_payload$signal, field = field, type = "binary_mask"))
      }

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
      lev <- unique(labs)
      lev <- lev[!is.na(lev)]
      col_map <- stats::setNames(rep("grey70", length(lev)), lev)
      valid_col <- are_valid_colors(lev)
      if (any(valid_col)) col_map[valid_col] <- lev[valid_col]
      if (any(!valid_col)) col_map[!valid_col] <- get_discrete_palette(sum(!valid_col), input$cluster_palette)

      if (isTRUE(binary_mode == "label_map")) {
        lev_ord <- sort(unique(labs[!is.na(labs)]))
        code_map <- stats::setNames(seq_along(lev_ord), lev_ord)
        code_mat <- matrix(NA_real_, nrow = ny, ncol = nx)
        keep_lab <- !is.na(lab_mat)
        code_mat[keep_lab] <- unname(code_map[lab_mat[keep_lab]])
        return(list(
          signal = code_mat,
          field = field,
          type = "categorical_label_map",
          label_matrix = lab_mat,
          color_map = col_map
        ))
      }

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

    histology_fit_signal_from_fused_ions <- function(msi_obj, n_ions = 12L) {
      obj <- msi_for_pdata()
      n_ions <- suppressWarnings(as.integer(n_ions))
      if (length(n_ions) != 1L || !isTRUE(is.finite(n_ions)) || n_ions < 2L) n_ions <- 12L
      n_ions <- as.integer(min(64L, max(2L, n_ions)))

      intensity_transform <- tolower(trimws(as.character(input$intensity_transform)[1]))
      if (!intensity_transform %in% c("none", "sqrt", "log1p", "asinh")) intensity_transform <- "none"
      smooth_sigma <- suppressWarnings(as.numeric(input$gaussian_sigma))
      if (!is.finite(smooth_sigma) || smooth_sigma <= 0) smooth_sigma <- 1
      use_smooth <- isTRUE(input$gaussian_smooth)

      cand_n <- as.integer(min(max(24L, 4L * n_ions), max(24L, nrow(obj))))
      cand_idx <- suggest_histology_fusion_indices(
        obj = obj,
        msi_obj = msi_obj,
        n_features = n_ions,
        candidate_pool = cand_n,
        intensity_transform = intensity_transform
      )
      if (length(cand_idx) == 0L) {
        mz_axis_fallback <- suppressWarnings(as.numeric(Cardinal::mz(obj)))
        cand_idx <- seq_len(min(length(mz_axis_fallback), n_ions))
        cand_idx <- cand_idx[is.finite(cand_idx) & cand_idx >= 1L]
        if (length(cand_idx) == 0L) {
          stop("Could not identify informative ions for fused Histology Fit target.")
        }
        message(sprintf("[Histology Image Fit] Fused target fallback: using first %d m/z features.", length(cand_idx)))
      }

      sp <- try(as.matrix(Cardinal::spectra(obj)[cand_idx, , drop = FALSE]), silent = TRUE)
      if (inherits(sp, "try-error") || is.null(sp) || nrow(sp) == 0L || ncol(sp) == 0L) {
        stop("Could not extract MSI intensities for fused Histology Fit target.")
      }
      if (nrow(sp) != length(cand_idx) && ncol(sp) == length(cand_idx)) {
        sp <- t(sp)
      }
      if (nrow(sp) != length(cand_idx)) {
        stop("Unexpected MSI matrix dimensions while building fused Histology Fit target.")
      }
      storage.mode(sp) <- "double"

      sig_list <- vector("list", length(cand_idx))
      sig_score <- rep(-Inf, length(cand_idx))
      for (ii in seq_along(cand_idx)) {
        sig_payload <- build_numeric_signal_matrix(
          vals = sp[ii, ],
          msi_obj = msi_obj,
          intensity_transform = intensity_transform,
          apply_transform = TRUE,
          gaussian_smooth = use_smooth,
          smooth_sigma = smooth_sigma,
          enhance = FALSE
        )
        sig_list[[ii]] <- sig_payload$signal
        sig_score[ii] <- combine_weighted_scores(
          c(
            score_spatial_coherence_map(sig_payload$signal, valid_mask = sig_payload$valid, sigma = 2),
            score_spatial_signal_map(sig_payload$signal, valid_mask = sig_payload$valid)
          ),
          c(0.75, 0.25)
        )
      }

      keep <- is.finite(sig_score)
      if (!any(keep)) {
        keep <- rep(TRUE, length(sig_list))
        sig_score <- rep(1, length(sig_list))
        message(sprintf("[Histology Image Fit] Fused target fallback: spatial ranking was not informative; using equal weights across %d ion maps.", length(sig_list)))
      }
      cand_idx <- cand_idx[keep]
      sig_list <- sig_list[keep]
      sig_score <- sig_score[keep]

      ord <- order(sig_score, decreasing = TRUE, na.last = NA)
      top_k <- min(length(ord), n_ions)
      ord <- ord[seq_len(top_k)]
      sel_idx <- cand_idx[ord]
      sel_scores <- sig_score[ord]
      sel_signals <- sig_list[ord]

      w <- as.numeric(sel_scores)
      w[!is.finite(w)] <- 0
      w <- w - min(w, na.rm = TRUE)
      if (!all(is.finite(w)) || sum(w) <= 1e-12) {
        w <- rep(1, length(sel_signals))
      } else {
        w <- w + 1e-6
      }
      w <- w / sum(w)

      fused <- matrix(0, nrow = as.integer(msi_obj$ny), ncol = as.integer(msi_obj$nx))
      for (ii in seq_along(sel_signals)) {
        fused <- fused + (w[ii] * sel_signals[[ii]])
      }
      fused[!is.finite(fused)] <- 0
      ok <- is.finite(fused) & fused > 0
      if (any(ok)) {
        q <- suppressWarnings(stats::quantile(fused[ok], probs = 0.75, na.rm = TRUE, names = FALSE, type = 8))
        if (is.finite(q) && q > 0) fused <- pmax(fused - q, 0)
      }
      mx <- suppressWarnings(max(fused, na.rm = TRUE))
      if (is.finite(mx) && mx > 0) fused <- fused / mx

      mz_axis <- suppressWarnings(as.numeric(Cardinal::mz(obj)))
      sel_mz <- rep(NA_real_, length(sel_idx))
      if (length(mz_axis) >= max(sel_idx)) {
        sel_mz <- mz_axis[sel_idx]
      }
      mz_txt <- format(sel_mz[is.finite(sel_mz)], digits = 10, scientific = FALSE, trim = TRUE)
      mz_preview <- if (length(mz_txt) > 0L) {
        paste(utils::head(mz_txt, 6L), collapse = ", ")
      } else {
        "n/a"
      }
      if (length(mz_txt) > 6L) {
        mz_preview <- paste0(mz_preview, ", ...")
      }

      list(
        signal = fused,
        field = NA_character_,
        type = "mz_fusion",
        label = sprintf("Fused MSI (%d ions)", length(sel_idx)),
        detail = mz_preview,
        mz_values = sel_mz,
        feature_index = sel_idx,
        feature_count = length(sel_idx)
      )
    }

    histology_fit_signal_from_pca <- function(msi_obj, component = 1L, n_components = 1L, n_features = 32L) {
      obj <- msi_for_pdata()
      component <- suppressWarnings(as.integer(component))
      n_components <- suppressWarnings(as.integer(n_components))
      n_features <- suppressWarnings(as.integer(n_features))
      if (length(component) != 1L || !isTRUE(is.finite(component)) || component < 1L) component <- 1L
      if (length(n_components) != 1L || !isTRUE(is.finite(n_components)) || n_components < 1L) n_components <- 1L
      if (length(n_features) != 1L || !isTRUE(is.finite(n_features)) || n_features < 4L) n_features <- 32L
      component <- as.integer(min(64L, max(1L, component)))
      n_components <- as.integer(min(8L, max(1L, n_components)))
      n_features <- as.integer(min(128L, max(4L, n_features)))

      intensity_transform <- tolower(trimws(as.character(input$intensity_transform)[1]))
      if (!intensity_transform %in% c("none", "sqrt", "log1p", "asinh")) intensity_transform <- "none"
      smooth_sigma <- suppressWarnings(as.numeric(input$gaussian_sigma))
      if (!is.finite(smooth_sigma) || smooth_sigma <= 0) smooth_sigma <- 1
      use_smooth <- isTRUE(input$gaussian_smooth)

      cand_n <- as.integer(min(max(48L, 4L * n_features), max(48L, nrow(obj))))
      sel_idx <- suggest_histology_fusion_indices(
        obj = obj,
        msi_obj = msi_obj,
        n_features = n_features,
        candidate_pool = cand_n,
        intensity_transform = intensity_transform
      )
      if (length(sel_idx) < 2L) {
        mz_axis_fallback <- suppressWarnings(as.numeric(Cardinal::mz(obj)))
        sel_idx <- seq_len(min(length(mz_axis_fallback), max(2L, n_features)))
        sel_idx <- sel_idx[is.finite(sel_idx) & sel_idx >= 1L]
        if (length(sel_idx) < 2L) {
          stop("Could not identify enough informative ions for PCA Histology Fit target.")
        }
        message(sprintf("[Histology Image Fit] PCA target fallback: using first %d m/z features.", length(sel_idx)))
      }

      sp <- try(as.matrix(Cardinal::spectra(obj)[sel_idx, , drop = FALSE]), silent = TRUE)
      if (inherits(sp, "try-error") || is.null(sp) || nrow(sp) < 2L || ncol(sp) < 2L) {
        stop("Could not extract MSI intensities for PCA Histology Fit target.")
      }
      if (nrow(sp) != length(sel_idx) && ncol(sp) == length(sel_idx)) {
        sp <- t(sp)
      }
      if (nrow(sp) != length(sel_idx)) {
        stop("Unexpected MSI matrix dimensions while building PCA Histology Fit target.")
      }

      storage.mode(sp) <- "double"
      if (!identical(intensity_transform, "none")) {
        sp <- transform_intensity(sp, intensity_transform)
      }
      sp[!is.finite(sp)] <- 0

      feat_mean <- rowMeans(sp)
      sp_centered <- sp - feat_mean
      feat_scale <- sqrt(rowMeans(sp_centered^2))
      feat_scale[!is.finite(feat_scale) | feat_scale <= 1e-8] <- 1
      sp_scaled <- sp_centered / feat_scale

      cov_mat <- tcrossprod(sp_scaled) / max(1, ncol(sp_scaled) - 1L)
      eig <- try(eigen(cov_mat, symmetric = TRUE), silent = TRUE)
      if (inherits(eig, "try-error") || is.null(eig$values) || is.null(eig$vectors)) {
        stop("Could not compute PCA for Histology Fit target.")
      }

      keep_comp <- is.finite(eig$values) & (eig$values > 1e-10)
      if (!any(keep_comp)) {
        stop("PCA target did not yield any finite components.")
      }
      values <- as.numeric(eig$values[keep_comp])
      vectors <- eig$vectors[, keep_comp, drop = FALSE]
      avg_signal <- suppressWarnings(colMeans(sp))
      mz_axis <- suppressWarnings(as.numeric(Cardinal::mz(obj)))
      sel_mz <- rep(NA_real_, length(sel_idx))
      if (length(mz_axis) >= max(sel_idx)) {
        sel_mz <- mz_axis[sel_idx]
      }

      comp_ids <- seq.int(component, min(ncol(vectors), component + n_components - 1L))
      if (length(comp_ids) == 0L) comp_ids <- 1L

      signal_list <- vector("list", length(comp_ids))
      weight_vec <- rep(0, length(comp_ids))
      var_vec <- rep(NA_real_, length(comp_ids))
      comp_structure_vec <- rep(NA_real_, length(comp_ids))
      comp_coherence_vec <- rep(NA_real_, length(comp_ids))
      load_mat <- matrix(0, nrow = length(sel_idx), ncol = length(comp_ids))
      for (ii in seq_along(comp_ids)) {
        comp_use <- comp_ids[ii]
        load_vec <- as.numeric(vectors[, comp_use])
        score_vec <- as.numeric(crossprod(load_vec, sp_scaled))
        orient <- suppressWarnings(stats::cor(score_vec, avg_signal, use = "pairwise.complete.obs", method = "spearman"))
        if ((is.finite(orient) && orient < 0) || (!is.finite(orient) && sum(load_vec, na.rm = TRUE) < 0)) {
          load_vec <- -load_vec
          score_vec <- -score_vec
        }
        load_mat[, ii] <- load_vec

        sig_payload <- build_numeric_signal_matrix(
          vals = score_vec,
          msi_obj = msi_obj,
          intensity_transform = "none",
          apply_transform = FALSE,
          gaussian_smooth = use_smooth,
          smooth_sigma = smooth_sigma,
          enhance = FALSE
        )
        signal_list[[ii]] <- sig_payload$signal
        var_frac <- values[comp_use] / sum(values)
        var_vec[ii] <- var_frac
        comp_structure_vec[ii] <- score_spatial_signal_map(sig_payload$signal, valid_mask = sig_payload$valid)
        comp_coherence_vec[ii] <- score_spatial_coherence_map(sig_payload$signal, valid_mask = sig_payload$valid, sigma = 2)
        spatial_weight <- combine_weighted_scores(
          c(comp_coherence_vec[ii], comp_structure_vec[ii]),
          c(0.75, 0.25)
        )
        spatial_gain <- if (is.finite(spatial_weight)) max(spatial_weight, 0) else 0
        spatial_penalty <- if (is.finite(spatial_weight) && spatial_weight < -0.25) 0.20 else 1.0
        weight_vec[ii] <- max(1e-6, if (is.finite(var_frac)) var_frac else 0) * spatial_penalty * (1 + spatial_gain)
      }

      if (!any(is.finite(weight_vec) & weight_vec > 0)) {
        weight_vec <- rep(1, length(signal_list))
      }
      weight_vec[!is.finite(weight_vec) | weight_vec <= 0] <- 0
      if (sum(weight_vec) <= 0) weight_vec <- rep(1, length(signal_list))
      weight_vec <- weight_vec / sum(weight_vec)

      pca_signal <- matrix(0, nrow = as.integer(msi_obj$ny), ncol = as.integer(msi_obj$nx))
      for (ii in seq_along(signal_list)) {
        pca_signal <- pca_signal + (weight_vec[ii] * signal_list[[ii]])
      }
      pca_signal <- matrix(rescale01(as.vector(pca_signal), enhance = FALSE), nrow = nrow(pca_signal), ncol = ncol(pca_signal))
      pca_signal[!is.finite(pca_signal)] <- 0

      top_loading_text <- function(loadings, mz_values, n_show = 3L, decreasing = TRUE) {
        ord <- order(loadings, decreasing = decreasing, na.last = NA)
        ord <- ord[seq_len(min(length(ord), n_show))]
        vals <- mz_values[ord]
        vals <- vals[is.finite(vals)]
        if (length(vals) == 0L) return("n/a")
        paste(format(vals, digits = 10, scientific = FALSE, trim = TRUE), collapse = ", ")
      }

      load_summary <- if (ncol(load_mat) == length(weight_vec) && length(weight_vec) > 0L) {
        as.numeric(drop(load_mat %*% weight_vec))
      } else if (ncol(load_mat) > 0L) {
        as.numeric(load_mat[, 1])
      } else {
        rep(NA_real_, length(sel_idx))
      }

      comp_label <- if (length(comp_ids) == 1L) {
        sprintf("PC%d", comp_ids[1])
      } else {
        sprintf("PC%d-%d", min(comp_ids), max(comp_ids))
      }
      var_txt <- paste(sprintf("PC%d %.1f%%", comp_ids, 100 * ifelse(is.finite(var_vec), var_vec, 0)), collapse = ", ")
      detail_txt <- sprintf(
        "%s | %s | loadings + %s | loadings - %s",
        comp_label,
        var_txt,
        top_loading_text(load_summary, sel_mz, decreasing = TRUE),
        top_loading_text(load_summary, sel_mz, decreasing = FALSE)
      )

      list(
        signal = pca_signal,
        field = NA_character_,
        type = "pca_component",
        label = sprintf("MSI PCA (%s)", comp_label),
        detail = detail_txt,
        mz_values = sel_mz,
        feature_index = sel_idx,
        feature_count = length(sel_idx),
        component_ids = comp_ids,
        component_variance = var_vec,
        component_weights = weight_vec,
        top_positive_mz = top_loading_text(load_summary, sel_mz, n_show = 8L, decreasing = TRUE),
        top_negative_mz = top_loading_text(load_summary, sel_mz, n_show = 8L, decreasing = FALSE),
        selected_mz_preview = {
          ord_abs <- order(abs(load_summary), decreasing = TRUE, na.last = NA)
          mz_preview_vals <- sel_mz[ord_abs]
          mz_preview_vals <- mz_preview_vals[is.finite(mz_preview_vals)]
          if (length(mz_preview_vals) == 0L) {
            "n/a"
          } else {
            paste(format(utils::head(mz_preview_vals, 12L), digits = 10, scientific = FALSE, trim = TRUE), collapse = ", ")
          }
        }
      )
    }

    resolve_optimization_signal_matrix <- function(msi_obj, source_mode = "current", field = NULL, context_label = "optimization") {
      source_mode <- tolower(trimws(as.character(source_mode)[1]))
      if (!source_mode %in% c("current", "pdata", "multi", "pca")) source_mode <- "current"

      histology_pdata_target_from_field <- function(msi_obj, field, label_txt = NULL) {
        out_preview <- edge_fit_signal_from_pdata_field(
          msi_obj = msi_obj,
          field = field,
          binary_mode = "label_map"
        )
        out_type <- if (!is.null(out_preview$type) && nzchar(as.character(out_preview$type)[1])) {
          as.character(out_preview$type)[1]
        } else {
          NA_character_
        }
        if (!identical(out_type, "categorical_label_map")) {
          return(list(
            signal = out_preview$signal,
            field = out_preview$field,
            type = out_type,
            label = if (!is.null(label_txt) && nzchar(as.character(label_txt)[1])) as.character(label_txt)[1] else paste0("pData: ", field),
            detail = NA_character_,
            label_matrix = if (!is.null(out_preview$label_matrix)) out_preview$label_matrix else NULL,
            color_map = if (!is.null(out_preview$color_map)) out_preview$color_map else NULL
          ))
        }

        out_score <- edge_fit_signal_from_pdata_field(
          msi_obj = msi_obj,
          field = field,
          binary_mode = "categorical"
        )
        list(
          signal = out_score$signal,
          field = out_score$field,
          type = "categorical_spatial_map",
          label = if (!is.null(label_txt) && nzchar(as.character(label_txt)[1])) as.character(label_txt)[1] else paste0("pData: ", field),
          detail = "Scored as all-label boundary-density map across the full categorical field.",
          label_matrix = if (!is.null(out_preview$label_matrix)) out_preview$label_matrix else NULL,
          color_map = if (!is.null(out_preview$color_map)) out_preview$color_map else NULL
        )
      }

      if (identical(source_mode, "current")) {
        use_histology_label_mode <- grepl("histology", tolower(context_label), fixed = TRUE) &&
          !is.null(msi_obj$mode) &&
          identical(as.character(msi_obj$mode)[1], "pdata") &&
          !is.null(msi_obj$pdata_field) &&
          nzchar(as.character(msi_obj$pdata_field)[1])
        if (isTRUE(use_histology_label_mode)) {
          out <- histology_pdata_target_from_field(
            msi_obj = msi_obj,
            field = as.character(msi_obj$pdata_field)[1],
            label_txt = if (!is.null(msi_obj$display_label)) as.character(msi_obj$display_label) else "Current MSI display"
          )
          return(list(
            signal = out$signal,
            source = "current",
            field = out$field,
            type = out$type,
            label = out$label,
            detail = out$detail,
            label_matrix = if (!is.null(out$label_matrix)) out$label_matrix else NULL,
            color_map = if (!is.null(out$color_map)) out$color_map else NULL
          ))
        }
        return(list(
          signal = msi_signal_matrix(msi_obj),
          source = "current",
          field = if (!is.null(msi_obj$pdata_field)) as.character(msi_obj$pdata_field) else NA_character_,
          type = if (!is.null(msi_obj$mode)) as.character(msi_obj$mode) else NA_character_,
          label = if (!is.null(msi_obj$display_label)) as.character(msi_obj$display_label) else "Current MSI display",
          detail = NA_character_
        ))
      }
      if (identical(source_mode, "multi")) {
        n_ions <- suppressWarnings(as.integer(input$histology_fit_multi_n))
        if (length(n_ions) != 1L || !isTRUE(is.finite(n_ions)) || n_ions < 2L) n_ions <- 12L
        out <- histology_fit_signal_from_fused_ions(msi_obj = msi_obj, n_ions = n_ions)
        return(list(
          signal = out$signal,
          source = "multi",
          field = out$field,
          type = out$type,
          label = out$label,
          detail = out$detail,
          mz_values = out$mz_values,
          feature_index = out$feature_index,
          feature_count = out$feature_count
        ))
      }
      if (identical(source_mode, "pca")) {
        comp_id <- suppressWarnings(as.integer(input$histology_fit_pca_component))
        n_comp <- suppressWarnings(as.integer(input$histology_fit_pca_n_components))
        n_features <- suppressWarnings(as.integer(input$histology_fit_pca_n_features))
        if (length(comp_id) != 1L || !isTRUE(is.finite(comp_id)) || comp_id < 1L) comp_id <- 1L
        if (length(n_comp) != 1L || !isTRUE(is.finite(n_comp)) || n_comp < 1L) n_comp <- 1L
        if (length(n_features) != 1L || !isTRUE(is.finite(n_features)) || n_features < 4L) n_features <- 32L
        out <- histology_fit_signal_from_pca(msi_obj = msi_obj, component = comp_id, n_components = n_comp, n_features = n_features)
        return(list(
          signal = out$signal,
          source = "pca",
          field = out$field,
          type = out$type,
          label = out$label,
          detail = out$detail,
          mz_values = out$mz_values,
          feature_index = out$feature_index,
          feature_count = out$feature_count
        ))
      }
      field <- as.character(field)[1]
      if (is.na(field) || !nzchar(field)) {
        stop(sprintf("Select a pData field for %s.", context_label))
      }
      is_histology_context <- grepl("histology", tolower(context_label), fixed = TRUE)
      if (isTRUE(is_histology_context)) {
        out <- histology_pdata_target_from_field(msi_obj = msi_obj, field = field, label_txt = paste0("pData: ", field))
      } else {
        out <- edge_fit_signal_from_pdata_field(msi_obj = msi_obj, field = field, binary_mode = "categorical")
      }
      list(
        signal = out$signal,
        source = "pdata",
        field = out$field,
        type = out$type,
        label = if (!is.null(out$label)) out$label else paste0("pData: ", out$field),
        detail = if (!is.null(out$detail)) out$detail else NA_character_,
        label_matrix = if (!is.null(out$label_matrix)) out$label_matrix else NULL,
        color_map = if (!is.null(out$color_map)) out$color_map else NULL
      )
    }

    edge_fit_signal_matrix <- function(msi_obj) {
      resolve_optimization_signal_matrix(
        msi_obj = msi_obj,
        source_mode = input$edge_fit_signal_source,
        field = input$edge_fit_pdata_field,
        context_label = "Edge Fit optimization"
      )
    }

    histology_fit_signal_matrix <- function(msi_obj) {
      resolve_optimization_signal_matrix(
        msi_obj = msi_obj,
        source_mode = input$histology_fit_signal_source,
        field = input$histology_fit_pdata_field,
        context_label = "Histology Fit optimization"
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

    build_histology_mask_fit_cache <- function(mask, band_px = 6L) {
      mask <- matrix(as.logical(mask), nrow = nrow(mask), ncol = ncol(mask))
      mask[!is.finite(mask)] <- FALSE
      bands <- build_mask_bands(mask, band_px = band_px)
      center_weight <- matrix(0, nrow = nrow(mask), ncol = ncol(mask))
      boundary_cache <- build_mask_boundary_distance_cache(mask)
      if (!is.null(boundary_cache) && !is.null(boundary_cache$dist)) {
        dist_in <- suppressWarnings(matrix(as.numeric(boundary_cache$dist), nrow = nrow(mask), ncol = ncol(mask)))
        dist_in[!mask | !is.finite(dist_in)] <- 0
        d_keep <- dist_in[mask & is.finite(dist_in) & dist_in > 0]
        if (length(d_keep) > 0L) {
          d_norm <- suppressWarnings(stats::quantile(d_keep, probs = 0.95, na.rm = TRUE, names = FALSE, type = 8))
          if (!is.finite(d_norm) || d_norm <= 0) d_norm <- suppressWarnings(max(d_keep, na.rm = TRUE))
          if (is.finite(d_norm) && d_norm > 0) {
            center_weight[mask] <- pmin(1, dist_in[mask] / d_norm)
          }
        }
      }
      bands$mask <- mask
      bands$center_weight <- center_weight
      bands
    }

    sample_shifted_weighted_signal <- function(coord_obj, signal, weight, dx = 0L, dy = 0L, min_n = 5L) {
      if (is.null(coord_obj) || length(coord_obj$r) == 0L || is.null(weight) || !is.matrix(weight)) {
        return(list(ok = FALSE, mean = NA_real_, n = 0L, cov = 0))
      }
      ny <- nrow(signal)
      nx <- ncol(signal)
      rr <- coord_obj$r + as.integer(dy)
      cc <- coord_obj$c + as.integer(dx)
      keep <- rr >= 1L & rr <= ny & cc >= 1L & cc <= nx
      n_keep <- sum(keep)
      cov <- n_keep / max(1L, as.integer(coord_obj$n_total))
      if (!is.finite(n_keep) || n_keep < min_n) {
        return(list(ok = FALSE, mean = NA_real_, n = as.integer(n_keep), cov = cov))
      }

      vals <- as.numeric(signal[cbind(rr[keep], cc[keep])])
      ww <- as.numeric(weight[cbind(rr[keep], cc[keep])])
      good <- is.finite(vals) & is.finite(ww) & ww > 0
      if (sum(good) < min_n) {
        return(list(ok = FALSE, mean = NA_real_, n = as.integer(sum(good)), cov = cov))
      }
      ww <- ww[good]
      vals <- vals[good]
      wsum <- sum(ww)
      if (!is.finite(wsum) || wsum <= 0) {
        return(list(ok = FALSE, mean = NA_real_, n = as.integer(length(vals)), cov = cov))
      }
      list(
        ok = TRUE,
        mean = sum(vals * ww) / wsum,
        n = as.integer(length(vals)),
        cov = cov
      )
    }

    score_mask_centrality_shift <- function(mask_cache, signal, dx = 0L, dy = 0L) {
      if (is.null(mask_cache) || is.null(mask_cache$inside_all)) return(-Inf)
      dx <- as.integer(round(dx))
      dy <- as.integer(round(dy))

      inside_all <- sample_shifted_stats(mask_cache$inside_all, signal, dx = dx, dy = dy, min_n = 10L)
      if (!isTRUE(inside_all$ok)) {
        inside_all <- sample_shifted_stats(mask_cache$inside_all, signal, dx = dx, dy = dy, min_n = 1L)
      }
      if (!isTRUE(inside_all$ok)) return(-Inf)

      center_term <- if (!is.null(mask_cache$center_weight)) {
        sample_shifted_weighted_signal(mask_cache$inside_all, signal, mask_cache$center_weight, dx = dx, dy = dy, min_n = 5L)
      } else {
        list(ok = FALSE, mean = NA_real_, cov = inside_all$cov)
      }
      core_term <- sample_shifted_stats(mask_cache$inside_core, signal, dx = dx, dy = dy, min_n = 5L)
      boundary_in <- if (length(mask_cache$inside_rings) >= 1L) {
        sample_shifted_stats(mask_cache$inside_rings[[1]], signal, dx = dx, dy = dy, min_n = 3L)
      } else {
        list(ok = FALSE, mean = NA_real_, cov = 0)
      }
      boundary_out <- if (length(mask_cache$outside_rings) >= 1L) {
        sample_shifted_stats(mask_cache$outside_rings[[1]], signal, dx = dx, dy = dy, min_n = 3L)
      } else {
        list(ok = FALSE, mean = NA_real_, cov = 0)
      }

      center_mean <- if (isTRUE(center_term$ok)) center_term$mean else if (isTRUE(core_term$ok)) core_term$mean else inside_all$mean
      boundary_in_mean <- if (isTRUE(boundary_in$ok)) boundary_in$mean else inside_all$mean
      boundary_out_mean <- if (isTRUE(boundary_out$ok)) boundary_out$mean else suppressWarnings(mean(as.numeric(signal), na.rm = TRUE))
      if (!is.finite(boundary_out_mean)) boundary_out_mean <- 0

      centeredness <- center_mean - boundary_in_mean
      inside_vs_outside <- center_mean - boundary_out_mean
      core_boost <- if (isTRUE(core_term$ok)) core_term$mean - inside_all$mean else 0
      global_contrast <- inside_all$mean - boundary_out_mean
      coverage <- pmax(
        0.01,
        pmin(
          1,
          min(
            if (is.finite(inside_all$cov)) inside_all$cov else 0,
            if (isTRUE(center_term$ok) && is.finite(center_term$cov)) center_term$cov else if (is.finite(inside_all$cov)) inside_all$cov else 0
          )
        )
      )

      base <- 0.45 * centeredness + 0.30 * inside_vs_outside + 0.15 * core_boost + 0.10 * global_contrast
      score <- base - 1.10 * (1 - coverage)^2
      if (!is.finite(score)) return(-Inf)
      score
    }

    build_edge_distance_cache <- function(signal, q = 0.85) {
      ny <- nrow(signal)
      nx <- ncol(signal)
      if (!is.finite(ny) || !is.finite(nx) || ny < 2 || nx < 2) {
        return(NULL)
      }
      edge_info <- build_edge_mask(signal, q = q)
      if (is.null(edge_info) || is.null(edge_info$edge)) return(NULL)
      out <- build_distance_cache_from_edge_mask(edge_info$edge, max_dist = 50L)
      if (is.null(out)) return(NULL)
      out$threshold <- edge_info$threshold
      out
    }

    histology_fit_crop_window <- function(msi_obj, tr_hist, tx, ty, search_range = 0L, pad_extra = 2L) {
      nx <- as.integer(msi_obj$nx)
      ny <- as.integer(msi_obj$ny)
      full_rows <- seq_len(ny)
      full_cols <- seq_len(nx)

      w <- suppressWarnings(as.numeric(tr_hist$width))
      h <- suppressWarnings(as.numeric(tr_hist$height))
      if (!is.finite(w) || w <= 0 || !is.finite(h) || h <= 0) {
        return(list(rows = full_rows, cols = full_cols, label = sprintf("%d x %d", length(full_rows), length(full_cols))))
      }

      pad <- suppressWarnings(as.integer(ceiling(search_range)))
      if (!is.finite(pad) || pad < 0L) pad <- 0L
      pad_extra <- suppressWarnings(as.integer(ceiling(pad_extra)))
      if (!is.finite(pad_extra) || pad_extra < 0L) pad_extra <- 2L
      pad <- pad + pad_extra

      cx <- (nx / 2) + as.numeric(tx)
      cy <- (ny / 2) + as.numeric(ty)
      x_lo <- floor(cx - (w / 2) - pad)
      x_hi <- ceiling(cx + (w / 2) + pad)
      y_lo <- floor(cy - (h / 2) - pad)
      y_hi <- ceiling(cy + (h / 2) + pad)

      x_norm <- as.integer(msi_obj$x_norm)
      y_norm <- as.integer(msi_obj$y_norm)
      pix_keep <- is.finite(x_norm) & is.finite(y_norm) &
        x_norm >= x_lo & x_norm <= x_hi &
        y_norm >= y_lo & y_norm <= y_hi

      if (!any(pix_keep)) {
        return(list(rows = full_rows, cols = full_cols, label = sprintf("%d x %d", length(full_rows), length(full_cols))))
      }

      row_idx <- as.integer(msi_obj$row_idx[pix_keep])
      col_idx <- as.integer(msi_obj$x_norm[pix_keep])
      row_min <- max(1L, min(row_idx, na.rm = TRUE))
      row_max <- min(ny, max(row_idx, na.rm = TRUE))
      col_min <- max(1L, min(col_idx, na.rm = TRUE))
      col_max <- min(nx, max(col_idx, na.rm = TRUE))

      rows <- seq.int(row_min, row_max)
      cols <- seq.int(col_min, col_max)
      list(rows = rows, cols = cols, label = sprintf("%d x %d", length(rows), length(cols)))
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

    binary_signal_to_mask <- function(signal, threshold = NULL) {
      if (is.null(signal) || !is.matrix(signal)) {
        return(matrix(FALSE, nrow = 0L, ncol = 0L))
      }
      sig <- suppressWarnings(matrix(as.numeric(signal), nrow = nrow(signal), ncol = ncol(signal)))
      keep <- is.finite(sig)
      out <- matrix(FALSE, nrow = nrow(sig), ncol = ncol(sig))
      if (!any(keep)) return(out)

      sig_keep <- sig[keep]
      sig_max <- suppressWarnings(max(sig_keep, na.rm = TRUE))
      if (!is.finite(sig_max) || sig_max <= 0) return(out)

      if (length(threshold) != 1L || !is.finite(threshold)) {
        sig_rng <- suppressWarnings(range(sig_keep, na.rm = TRUE))
        if (!all(is.finite(sig_rng))) return(out)
        threshold <- if (diff(sig_rng) <= 1e-8) {
          sig_rng[2] / 2
        } else {
          mean(sig_rng)
        }
      }
      out[keep] <- sig[keep] > as.numeric(threshold)
      out
    }

    combine_weighted_scores <- function(values, weights = NULL) {
      vals <- as.numeric(values)
      if (length(vals) == 0L) return(NA_real_)
      if (is.null(weights) || length(weights) != length(vals)) {
        weights <- rep(1, length(vals))
      }
      weights <- as.numeric(weights)
      keep <- is.finite(vals) & is.finite(weights) & weights > 0
      if (!any(keep)) return(NA_real_)
      sum(vals[keep] * weights[keep]) / sum(weights[keep])
    }

    build_structure_signal_map <- function(signal, valid_mask = NULL, sigma = 4, positive_only = FALSE) {
      if (is.null(signal) || !is.matrix(signal) || nrow(signal) < 2L || ncol(signal) < 2L) {
        return(NULL)
      }

      ny <- nrow(signal)
      nx <- ncol(signal)
      sig <- suppressWarnings(matrix(as.numeric(signal), nrow = ny, ncol = nx))
      sig[!is.finite(sig)] <- 0

      if (is.null(valid_mask)) {
        valid <- is.finite(signal)
      } else {
        valid <- matrix(as.logical(valid_mask), nrow = ny, ncol = nx)
      }
      valid[!is.finite(valid)] <- FALSE
      if (sum(valid, na.rm = TRUE) < 25L) return(NULL)

      sigma <- suppressWarnings(as.numeric(sigma))
      if (!is.finite(sigma) || sigma <= 0) sigma <- 4

      sig_masked <- sig
      sig_masked[!valid] <- NA_real_
      smooth <- gaussian_smooth_matrix(sig_masked, sigma = sigma)
      detail <- sig - smooth
      detail[!valid] <- NA_real_
      detail_work <- if (isTRUE(positive_only)) pmax(detail, 0) else abs(detail)

      right <- cbind(sig[, -1, drop = FALSE], sig[, nx, drop = FALSE])
      down <- rbind(sig[-1, , drop = FALSE], sig[ny, , drop = FALSE])
      right_valid <- cbind(valid[, -1, drop = FALSE], FALSE)
      down_valid <- rbind(valid[-1, , drop = FALSE], rep(FALSE, nx))
      grad <- abs(right - sig) + abs(down - sig)
      grad[!(valid & (right_valid | down_valid))] <- NA_real_

      detail_sig <- matrix(rescale01(as.vector(detail_work), enhance = FALSE), nrow = ny, ncol = nx)
      grad_sig <- matrix(rescale01(as.vector(grad), enhance = FALSE), nrow = ny, ncol = nx)
      out <- (0.70 * detail_sig) + (0.30 * grad_sig)
      out[!is.finite(out)] <- 0
      out[!valid] <- 0

      list(
        signal = out,
        valid = valid,
        detail = detail,
        gradient = grad
      )
    }

    build_distance_cache_from_edge_mask <- function(edge_mask, max_dist = 50L) {
      if (is.null(edge_mask) || !is.matrix(edge_mask) || nrow(edge_mask) < 2L || ncol(edge_mask) < 2L) {
        return(NULL)
      }
      edge <- matrix(as.logical(edge_mask), nrow = nrow(edge_mask), ncol = ncol(edge_mask))
      edge[!is.finite(edge)] <- FALSE
      n_edge <- suppressWarnings(sum(edge, na.rm = TRUE))
      if (!is.finite(n_edge) || n_edge < 5L) return(NULL)

      ny <- nrow(edge)
      nx <- ncol(edge)
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

      max_dist <- suppressWarnings(as.numeric(max_dist))
      if (!is.finite(max_dist) || max_dist < 1) max_dist <- 50
      d <- pmin(d, max_dist)
      list(dist = d, edge = edge, n_edge = as.integer(n_edge))
    }

    build_mask_boundary_distance_cache <- function(mask) {
      if (is.null(mask) || !is.matrix(mask)) return(NULL)
      mask <- matrix(as.logical(mask), nrow = nrow(mask), ncol = ncol(mask))
      mask[!is.finite(mask)] <- FALSE
      if (sum(mask, na.rm = TRUE) < 10L) return(NULL)
      boundary <- mask & !erode8(mask)
      build_distance_cache_from_edge_mask(boundary)
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
      if (length(pick_rank) != 1L || !is.finite(pick_rank)) pick_rank <- cand$rank[1]
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

    output$histology_fit_preview_ui <- renderUI({
      cand <- xh$histology_fit_candidates
      if (is.null(cand) || nrow(cand) == 0) {
        return(tags$small("Histology-fit preview: run Histology Fit to generate up to 5 candidates."))
      }
      labels <- sprintf(
        "%s%s#%d dX=%.2f dY=%.2f score=%.4f | int=%.4f | corr=%.4f | edge=%.4f%s",
        ifelse(!is.null(cand$selected) & cand$selected, "[selected] ", ""),
        ifelse(!is.null(cand$start_position) & cand$start_position, "[start] ", ""),
        cand$rank, cand$dX, cand$dY, cand$score, cand$intensity_score, cand$signed_corr, cand$edge_score,
        ifelse(!is.null(cand$peak_score) & cand$peak_score, " | raw peak", "")
      )
      vals <- as.character(cand$rank)
      names(vals) <- labels
      selectInput(
        ns("histology_fit_choice"),
        "Histology-fit candidates (top 5)",
        choices = vals,
        selected = vals[1]
      )
    })

    observeEvent(input$apply_histology_fit_choice, {
      cand <- xh$histology_fit_candidates
      if (is.null(cand) || nrow(cand) == 0) {
        showNotification("No histology-fit candidates available yet. Run Histology Fit first.", type = "warning", duration = 6)
        return()
      }
      pick_rank <- suppressWarnings(as.integer(input$histology_fit_choice))
      if (length(pick_rank) != 1L || !is.finite(pick_rank)) pick_rank <- cand$rank[1]
      row <- cand[cand$rank == pick_rank, , drop = FALSE]
      if (nrow(row) == 0) row <- cand[1, , drop = FALSE]

      updateSliderInput(session, "translate_x", value = row$translate_x[1])
      updateSliderInput(session, "translate_y", value = row$translate_y[1])
      updateNumericInput(session, "translate_x_num", value = row$translate_x[1])
      updateNumericInput(session, "translate_y_num", value = row$translate_y[1])

      showNotification(
        sprintf(
          "Applied histology-fit candidate #%d: dX=%.2f, dY=%.2f, score=%.4f",
          row$rank[1], row$dX[1], row$dY[1], row$score[1]
        ),
        type = "message",
        duration = 6
      )
    }, ignoreInit = TRUE)

    output$histology_fit_summary <- renderPrint({
      sm <- xh$histology_fit_summary
      if (is.null(sm)) {
        cat("No histology-fit results yet.\n")
        cat("Run Histology Fit to score local XY candidates around the current transform.\n")
        return(invisible(NULL))
      }
      print(sm)
      cand <- xh$histology_fit_candidates
      if (!is.null(cand) && nrow(cand) > 0L) {
        cat("\nTop histology-fit candidates:\n")
        print(cand)
      }
    })

    output$histology_fit_target_info_ui <- renderUI({
      info <- xh$histology_fit_target_info
      if (is.null(info) || is.null(info$signal) || !is.matrix(info$signal)) {
        return(tags$small("Run Histology Fit to preview the exact target signal image used for scoring."))
      }
      detail_txt <- if (!is.null(info$detail) && nzchar(as.character(info$detail)[1])) as.character(info$detail)[1] else NULL
      feature_txt <- if (!is.null(info$feature_count) && is.finite(as.numeric(info$feature_count)[1])) {
        sprintf("Features: %d", as.integer(info$feature_count)[1])
      } else {
        NULL
      }
      selected_mz_txt <- if (!is.null(info$selected_mz_preview) && nzchar(as.character(info$selected_mz_preview)[1])) {
        sprintf("Selected ions preview: %s", as.character(info$selected_mz_preview)[1])
      } else {
        NULL
      }
      pca_weight_txt <- if (!is.null(info$component_ids) && !is.null(info$component_weights)) {
        comp_ids <- suppressWarnings(as.integer(info$component_ids))
        comp_w <- suppressWarnings(as.numeric(info$component_weights))
        keep <- is.finite(comp_ids) & is.finite(comp_w)
        if (any(keep)) {
          comp_ids <- comp_ids[keep]
          comp_w <- comp_w[keep]
          sprintf(
            "Component weights: %s",
            paste(sprintf("PC%d %.2f", comp_ids, comp_w), collapse = ", ")
          )
        } else {
          NULL
        }
      } else {
        NULL
      }
      loading_txt <- if (!is.null(info$top_positive_mz) && nzchar(as.character(info$top_positive_mz)[1])) {
        sprintf(
          "Blended PCA loadings: + %s | - %s",
          as.character(info$top_positive_mz)[1],
          if (!is.null(info$top_negative_mz) && nzchar(as.character(info$top_negative_mz)[1])) as.character(info$top_negative_mz)[1] else "n/a"
        )
      } else {
        NULL
      }
      tags$div(
        tags$strong(if (!is.null(info$label) && nzchar(as.character(info$label)[1])) as.character(info$label)[1] else "Histology-fit target"),
        if (!is.null(detail_txt)) tags$p(style = "margin: 6px 0 0 0;", detail_txt),
        tags$p(style = "margin: 6px 0 0 0;", sprintf("Signal type: %s", if (!is.null(info$type) && nzchar(as.character(info$type)[1])) as.character(info$type)[1] else "n/a")),
        if (!is.null(feature_txt)) tags$p(style = "margin: 6px 0 0 0;", feature_txt),
        if (!is.null(pca_weight_txt)) tags$p(style = "margin: 6px 0 0 0;", pca_weight_txt),
        if (!is.null(loading_txt)) tags$p(style = "margin: 6px 0 0 0;", loading_txt),
        if (!is.null(selected_mz_txt)) tags$p(style = "margin: 6px 0 0 0;", selected_mz_txt)
      )
    })

    draw_histology_fit_target_plot <- function(info, main = NULL) {
      if (is.null(info) || is.null(info$signal) || !is.matrix(info$signal)) return(invisible(NULL))
      sig <- info$signal
      if (!is.null(info$type) &&
          identical(as.character(info$type)[1], "categorical_label_map") &&
          !is.null(info$label_matrix) &&
          is.matrix(info$label_matrix) &&
          !is.null(info$color_map)) {
        lab_mat <- info$label_matrix
        col_map <- as.character(info$color_map)
        if (!is.null(names(info$color_map))) names(col_map) <- names(info$color_map)
        col_mat <- matrix("#00000000", nrow = nrow(lab_mat), ncol = ncol(lab_mat))
        keep <- !is.na(lab_mat)
        col_mat[keep] <- unname(col_map[lab_mat[keep]])
        rast <- as.raster(col_mat)
      } else {
        palette_use <- if (!is.null(info$type) &&
          identical(as.character(info$type)[1], "categorical_spatial_map")) {
          "Inferno"
        } else {
          input$msi_palette
        }
        rast <- signal_matrix_to_raster(sig, palette_name = palette_use)
      }
      if (is.null(rast)) return(invisible(NULL))
      if (is.null(main) || !nzchar(trimws(as.character(main)[1]))) {
        main <- if (!is.null(info$label) && nzchar(as.character(info$label)[1])) {
          sprintf("Histology-fit target (%s)", as.character(info$label)[1])
        } else {
          "Histology-fit target"
        }
      }
      graphics::par(mar = c(0.5, 0.5, 2.2, 0.5))
      graphics::plot.new()
      graphics::plot.window(xlim = c(1, ncol(sig)), ylim = c(nrow(sig), 1), asp = 1, xaxs = "i", yaxs = "i")
      graphics::rasterImage(rast, 1, 1, ncol(sig), nrow(sig), interpolate = FALSE)
      graphics::box(col = "grey70")
      graphics::title(main = main)
    }

    output$histology_fit_target_plot <- renderPlot({
      if (!isTRUE(input$show_fit_info)) return(invisible(NULL))
      info <- xh$histology_fit_target_info
      validate(need(!is.null(info) && !is.null(info$signal) && is.matrix(info$signal), "Run Histology Fit to preview the target signal image."))
      draw_histology_fit_target_plot(info)
    })

    output$download_histology_fit_target_png <- downloadHandler(
      filename = function() {
        info <- xh$histology_fit_target_info
        stem <- if (!is.null(info$label) && nzchar(as.character(info$label)[1])) {
          gsub("[^A-Za-z0-9._-]+", "_", as.character(info$label)[1])
        } else {
          "histology_fit_target"
        }
        sprintf("%s_%s.png", stem, format(Sys.time(), "%Y%m%d_%H%M%S"))
      },
      content = function(file) {
        info <- xh$histology_fit_target_info
        if (is.null(info) || is.null(info$signal) || !is.matrix(info$signal)) {
          stop("No Histology Fit target image is available yet. Run Histology Fit first.")
        }
        grDevices::png(file, width = 1400, height = 1100, res = 170)
        on.exit(grDevices::dev.off(), add = TRUE)
        draw_histology_fit_target_plot(info)
      }
    )

    observeEvent(input$run_histology_fit, {
      req(msi_for_pdata())
      histology_fit_notify <- function(msg, type = "message", duration = 7) {
        message(sprintf("[Histology Image Fit] %s", msg))
        showNotification(msg, type = type, duration = duration)
      }
      if (is.null(input$histology_upload) || !nzchar(input$histology_upload$name)) {
        histology_fit_notify("Load a histology image before running Histology Fit.", type = "warning", duration = 7)
        return()
      }
      xh$histology_fit_grid <- NULL
      xh$histology_fit_candidates <- NULL
      xh$histology_fit_summary <- NULL
      xh$histology_fit_signature <- NULL
      xh$histology_fit_target_info <- NULL

      txy <- current_translate_xy()
      tx0 <- as.numeric(txy$tx)
      ty0 <- as.numeric(txy$ty)

      rng <- suppressWarnings(as.integer(input$histology_fit_range))
      stp <- suppressWarnings(as.integer(input$histology_fit_step))
      if (!is.finite(rng) || rng < 1L) rng <- 30L
      if (!is.finite(stp) || stp < 1L) stp <- 2L
      rng <- as.integer(min(500L, max(1L, rng)))
      stp <- as.integer(min(50L, max(1L, stp)))
      message(sprintf(
        "[Histology Image Fit] Starting local XY refinement at X=%.2f Y=%.2f range=%d step=%d",
        tx0,
        ty0,
        rng,
        stp
      ))
      histology_fit_relation <- resolve_histology_intensity_relation(input$histology_fit_intensity_relation)

      msi <- make_msi_raster()
      feat_mode <- tolower(trimws(as.character(input$histology_feature_mode)[1]))
      if (!feat_mode %in% c("hematoxylin", "darkness", "purple")) feat_mode <- "hematoxylin"
      histology_use_smooth <- isTRUE(input$gaussian_smooth)
      histology_smooth_sigma <- suppressWarnings(as.numeric(input$gaussian_sigma))
      if (!is.finite(histology_smooth_sigma) || histology_smooth_sigma <= 0) histology_smooth_sigma <- 1

      ncores_req <- suppressWarnings(as.integer(tryCatch(setup_values()[["ncores"]], error = function(e) NA_integer_)))
      if (!is.finite(ncores_req) || ncores_req < 1L) ncores_req <- 1L
      ncores_detect <- suppressWarnings(as.integer(tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA_integer_)))
      if (!is.finite(ncores_detect) || ncores_detect < 1L) ncores_detect <- ncores_req
      ncores_use <- as.integer(max(1L, min(ncores_req, ncores_detect)))
      can_parallel_histology <- (.Platform$OS.type != "windows") &&
        requireNamespace("parallel", quietly = TRUE) &&
        ncores_use > 1L

      withProgress(message = "Running Histology Fit", value = 0.08, {
        incProgress(0.18, detail = "Preparing MSI target map")
        sig_info <- try(histology_fit_signal_matrix(msi), silent = TRUE)
        if (inherits(sig_info, "try-error") || is.null(sig_info) || !is.matrix(sig_info$signal)) {
          xh$histology_fit_target_info <- NULL
          histology_fit_notify("Could not build Histology Fit target signal. Check current display mode or selected pData field.", type = "error", duration = 8)
          return()
        }
        xh$histology_fit_target_info <- list(
          signal = sig_info$signal,
          source = if (!is.null(sig_info$source)) sig_info$source else NA_character_,
          field = if (!is.null(sig_info$field)) sig_info$field else NA_character_,
          type = if (!is.null(sig_info$type)) sig_info$type else NA_character_,
          label = if (!is.null(sig_info$label)) sig_info$label else NA_character_,
          detail = if (!is.null(sig_info$detail)) sig_info$detail else NA_character_,
          feature_count = if (!is.null(sig_info$feature_count)) sig_info$feature_count else NA_integer_,
          mz_values = if (!is.null(sig_info$mz_values)) sig_info$mz_values else NULL,
          component_ids = if (!is.null(sig_info$component_ids)) sig_info$component_ids else NULL,
          component_variance = if (!is.null(sig_info$component_variance)) sig_info$component_variance else NULL,
          component_weights = if (!is.null(sig_info$component_weights)) sig_info$component_weights else NULL,
          top_positive_mz = if (!is.null(sig_info$top_positive_mz)) sig_info$top_positive_mz else NULL,
          top_negative_mz = if (!is.null(sig_info$top_negative_mz)) sig_info$top_negative_mz else NULL,
          selected_mz_preview = if (!is.null(sig_info$selected_mz_preview)) sig_info$selected_mz_preview else NULL,
          label_matrix = if (!is.null(sig_info$label_matrix)) sig_info$label_matrix else NULL,
          color_map = if (!is.null(sig_info$color_map)) sig_info$color_map else NULL
        )
        signal <- suppressWarnings(matrix(as.numeric(sig_info$signal), nrow = nrow(sig_info$signal), ncol = ncol(sig_info$signal)))
        signal[!is.finite(signal)] <- 0
        message(sprintf(
          "[Histology Image Fit] Target signal source=%s field=%s type=%s label=%s detail=%s feature=%s",
          if (!is.null(sig_info$source)) as.character(sig_info$source) else "unknown",
          if (!is.null(sig_info$field) && nzchar(as.character(sig_info$field))) as.character(sig_info$field) else "n/a",
          if (!is.null(sig_info$type) && nzchar(as.character(sig_info$type))) as.character(sig_info$type) else "n/a",
          if (!is.null(sig_info$label) && nzchar(as.character(sig_info$label))) as.character(sig_info$label) else "n/a",
          if (!is.null(sig_info$detail) && nzchar(as.character(sig_info$detail))) as.character(sig_info$detail) else "n/a",
          feat_mode
        ))

        incProgress(0.14, detail = "Caching transformed histology image")
        tr_hist <- try(
          transform_overlay_image(
            histology_image(),
            use_point_filter = FALSE,
            alpha_scale = 1,
            scale_correction = list(fx = 1, fy = 1)
          ),
          silent = TRUE
        )
        if (inherits(tr_hist, "try-error") || is.null(tr_hist) || is.null(tr_hist$raster)) {
          xh$histology_fit_target_info <- NULL
          histology_fit_notify("Could not transform histology image for Histology Fit.", type = "error", duration = 8)
          return()
        }
        fit_window <- histology_fit_crop_window(msi, tr_hist, tx = tx0, ty = ty0, search_range = rng, pad_extra = stp)
        signal_use <- signal[fit_window$rows, fit_window$cols, drop = FALSE]
        target_type <- if (!is.null(sig_info$type) && nzchar(as.character(sig_info$type)[1])) {
          as.character(sig_info$type)[1]
        } else {
          NA_character_
        }
        label_use <- if (!is.null(sig_info$type) &&
          identical(as.character(sig_info$type)[1], "categorical_label_map") &&
          !is.null(sig_info$label_matrix) &&
          is.matrix(sig_info$label_matrix)) {
          sig_info$label_matrix[fit_window$rows, fit_window$cols, drop = FALSE]
        } else {
          NULL
        }
        signal_mask <- NULL
        binary_mask_cache <- NULL
        binary_mask_reason <- NA_character_
        binary_mask_n <- NA_integer_
        binary_mask_total <- as.integer(length(signal_use))
        binary_mask_frac <- NA_real_
        if (identical(target_type, "binary_mask")) {
          signal_mask <- binary_signal_to_mask(signal_use)
          n_mask <- suppressWarnings(sum(signal_mask, na.rm = TRUE))
          binary_mask_n <- suppressWarnings(as.integer(n_mask))
          if (is.finite(n_mask) && is.finite(binary_mask_total) && binary_mask_total > 0L) {
            binary_mask_frac <- n_mask / binary_mask_total
          }
          if (!is.finite(n_mask) || n_mask <= 0L) {
            binary_mask_reason <- "mask has no positive pixels in the scoring window"
          } else if (n_mask < 25L) {
            binary_mask_reason <- sprintf("mask is too small in the scoring window (%d positive pixels)", as.integer(n_mask))
          } else if (n_mask > (length(signal_mask) - 25L)) {
            binary_mask_reason <- sprintf(
              "mask nearly fills the scoring window (%d / %d positive pixels)",
              as.integer(n_mask),
              as.integer(length(signal_mask))
            )
          } else {
            binary_mask_cache <- build_histology_mask_fit_cache(signal_mask, band_px = 6L)
            if (is.null(binary_mask_cache)) {
              binary_mask_reason <- "mask cache could not be built"
            }
          }
          if (!is.null(binary_mask_cache)) {
            message(sprintf(
              "[Histology Image Fit] Binary mask stats window=%s positive_pixels=%d total_pixels=%d fraction=%.4f mask_scoring=enabled",
              fit_window$label,
              if (is.finite(binary_mask_n)) binary_mask_n else -1L,
              if (is.finite(binary_mask_total)) binary_mask_total else -1L,
              if (is.finite(binary_mask_frac)) binary_mask_frac else NA_real_
            ))
          } else {
            histology_fit_notify(
              sprintf(
                "Binary-mask Histology Fit could not initialize mask scoring for window %s. %s.",
                fit_window$label,
                if (!is.na(binary_mask_reason) && nzchar(binary_mask_reason)) binary_mask_reason else "Mask statistics were invalid"
              ),
              type = "warning",
              duration = 10
            )
            return()
          }
        }
        target_structure_use <- if (is.null(label_use) && !identical(target_type, "binary_mask")) {
          build_structure_signal_map(signal_use, valid_mask = is.finite(signal_use), positive_only = FALSE)
        } else {
          NULL
        }
        edge_cache <- if (!is.null(label_use)) {
          NULL
        } else if (identical(target_type, "binary_mask")) {
          build_mask_boundary_distance_cache(signal_mask)
        } else if (!is.null(target_structure_use) && !is.null(target_structure_use$signal)) {
          build_edge_distance_cache(target_structure_use$signal)
        } else {
          build_edge_distance_cache(signal_use)
        }
        message(sprintf(
          "[Histology Image Fit] Scoring window=%s target_type=%s mask_scoring=%s structure_scoring=%s",
          fit_window$label,
          if (!is.null(target_type) && nzchar(target_type)) target_type else "n/a",
          if (!is.null(binary_mask_cache)) "enabled" else "disabled",
          if (!is.null(target_structure_use) || identical(target_type, "binary_mask")) "enabled" else "disabled"
        ))

        score_histology_pair <- function(dx, dy) {
          tx_abs <- clamp(tx0 + as.numeric(dx), -1000, 1000)
          ty_abs <- clamp(ty0 + as.numeric(dy), -1000, 1000)

          hp <- try(
            histology_feature_payload(
              msi_obj = msi,
              tr_hist = tr_hist,
              tx = tx_abs,
              ty = ty_abs,
              feature_mode = feat_mode,
              gaussian_smooth = histology_use_smooth,
              smooth_sigma = histology_smooth_sigma
            ),
            silent = TRUE
          )
          if (inherits(hp, "try-error") || is.null(hp) || !is.matrix(hp$signal) || !is.matrix(hp$valid)) {
            return(list(
              dX = as.integer(dx),
              dY = as.integer(dy),
              raw_score = NA_real_,
              corr_score = NA_real_,
              nmi_score = NA_real_,
              edge_score = NA_real_,
              overlap = 0L,
              n_edge = 0L
            ))
          }

          hp_signal <- hp$signal[fit_window$rows, fit_window$cols, drop = FALSE]
          hp_valid <- hp$valid[fit_window$rows, fit_window$cols, drop = FALSE]
          keep <- if (!is.null(label_use)) {
            hp_valid & !is.na(label_use) & nzchar(trimws(as.character(label_use)))
          } else {
            hp_valid & is.finite(signal_use)
          }
          overlap_n <- suppressWarnings(sum(keep, na.rm = TRUE))
          if (!is.finite(overlap_n) || overlap_n < 25L) {
            return(list(
              dX = as.integer(dx),
              dY = as.integer(dy),
              raw_score = NA_real_,
              corr_score = NA_real_,
              nmi_score = NA_real_,
              edge_score = NA_real_,
              overlap = as.integer(ifelse(is.finite(overlap_n), overlap_n, 0L)),
              n_edge = 0L
            ))
          }

          hp_structure <- build_structure_signal_map(hp_signal, valid_mask = hp_valid, positive_only = TRUE)

          if (!is.null(label_use)) {
            sc_corr <- NA_real_
            sc_nmi <- NA_real_
            sc_intensity <- categorical_signal_association(hp_signal[keep], label_use[keep], min_group_n = 5L)
            sc_edge <- NA_real_
            n_edge <- 0L
          } else if (!is.null(binary_mask_cache)) {
            sc_corr <- signed_rank_correlation(hp_signal[keep], signal_use[keep])
            sc_nmi <- normalized_mutual_information(hp_signal[keep], signal_use[keep], nbins = 24L)
            sc_mask_raw <- combine_weighted_scores(
              c(
                score_mask_shift(binary_mask_cache, hp_signal, dx = 0L, dy = 0L),
                score_mask_centrality_shift(binary_mask_cache, hp_signal, dx = 0L, dy = 0L)
              ),
              c(0.35, 0.65)
            )
            sc_mask_structure <- if (!is.null(hp_structure) && !is.null(hp_structure$signal)) {
              combine_weighted_scores(
                c(
                  score_mask_shift(binary_mask_cache, hp_structure$signal, dx = 0L, dy = 0L),
                  score_mask_centrality_shift(binary_mask_cache, hp_structure$signal, dx = 0L, dy = 0L)
                ),
                c(0.25, 0.75)
              )
            } else {
              NA_real_
            }
            sc_mask <- combine_weighted_scores(
              c(sc_mask_raw, sc_mask_structure),
              c(0.85, 0.15)
            )
            sc_intensity <- switch(
              histology_fit_relation,
              direct = sc_mask,
              inverse = if (is.finite(sc_mask)) -sc_mask else NA_real_,
              either = if (is.finite(sc_mask)) abs(sc_mask) else NA_real_,
              if (is.finite(sc_mask)) abs(sc_mask) else NA_real_
            )
            edge_info <- if (!is.null(hp_structure) && !is.null(hp_structure$signal)) {
              build_edge_mask(hp_structure$signal, valid_mask = hp_structure$valid)
            } else {
              build_edge_mask(hp_signal, valid_mask = hp_valid)
            }
            sc_edge <- if (!is.null(edge_info)) score_edge_mask_against_distance(edge_info$edge, edge_cache) else NA_real_
            n_edge <- if (is.null(edge_info) || is.null(edge_info$n_edge)) 0L else as.integer(edge_info$n_edge)
          } else {
            sc_corr <- signed_rank_correlation(hp_signal[keep], signal_use[keep])
            sc_nmi <- normalized_mutual_information(hp_signal[keep], signal_use[keep], nbins = 24L)
            sc_relation <- switch(
              histology_fit_relation,
              direct = sc_corr,
              inverse = if (is.finite(sc_corr)) -sc_corr else NA_real_,
              either = if (is.finite(sc_corr)) abs(sc_corr) else NA_real_,
              if (is.finite(sc_corr)) abs(sc_corr) else NA_real_
            )
            sc_structure_corr <- NA_real_
            sc_structure_nmi <- NA_real_
            if (!is.null(hp_structure) &&
              !is.null(target_structure_use) &&
              !is.null(target_structure_use$signal) &&
              is.matrix(target_structure_use$signal)) {
              keep_structure <- hp_structure$valid & target_structure_use$valid
              keep_structure_n <- suppressWarnings(sum(keep_structure, na.rm = TRUE))
              if (is.finite(keep_structure_n) && keep_structure_n >= 25L) {
                sc_structure_corr <- signed_rank_correlation(
                  hp_structure$signal[keep_structure],
                  target_structure_use$signal[keep_structure]
                )
                sc_structure_nmi <- normalized_mutual_information(
                  hp_structure$signal[keep_structure],
                  target_structure_use$signal[keep_structure],
                  nbins = 24L
                )
              }
            }
            sc_structure <- combine_weighted_scores(
              c(sc_structure_corr, sc_structure_nmi),
              c(0.65, 0.35)
            )
            sc_intensity <- combine_weighted_scores(
              c(sc_relation, sc_structure),
              c(0.30, 0.70)
            )
            edge_info <- if (!is.null(hp_structure) && !is.null(hp_structure$signal)) {
              build_edge_mask(hp_structure$signal, valid_mask = hp_structure$valid)
            } else {
              build_edge_mask(hp_signal, valid_mask = hp_valid)
            }
            sc_edge <- if (!is.null(edge_info)) score_edge_mask_against_distance(edge_info$edge, edge_cache) else NA_real_
            n_edge <- if (is.null(edge_info) || is.null(edge_info$n_edge)) 0L else as.integer(edge_info$n_edge)
          }

          list(
            dX = as.integer(dx),
            dY = as.integer(dy),
            raw_score = as.numeric(sc_intensity),
            corr_score = as.numeric(sc_corr),
            nmi_score = as.numeric(sc_nmi),
            edge_score = as.numeric(sc_edge),
            overlap = as.integer(overlap_n),
            n_edge = as.integer(n_edge)
          )
        }

        score_histology_grid <- function(dx_set, dy_set, progress_weight = 0, detail = NULL, allow_parallel = TRUE) {
          empty_df <- data.frame(
            dX = integer(0),
            dY = integer(0),
            raw_score = numeric(0),
            corr_score = numeric(0),
            nmi_score = numeric(0),
            edge_score = numeric(0),
            overlap = integer(0),
            n_edge = integer(0),
            stringsAsFactors = FALSE
          )
          if (!is.null(detail)) incProgress(0, detail = detail)
          grid_pairs <- expand.grid(
            dX = as.integer(dx_set),
            dY = as.integer(dy_set),
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE
          )
          if (nrow(grid_pairs) == 0L) return(empty_df)

          res_list <- vector("list", nrow(grid_pairs))
          use_parallel <- isTRUE(allow_parallel) && isTRUE(can_parallel_histology) && nrow(grid_pairs) >= (2L * ncores_use)
          if (isTRUE(use_parallel)) {
            chunk_size <- max(64L, 16L * ncores_use)
            chunks <- split(seq_len(nrow(grid_pairs)), ceiling(seq_len(nrow(grid_pairs)) / chunk_size))
            parallel_failed <- FALSE
            message(sprintf(
              "[Histology Image Fit] Parallel grid launch pairs=%d chunks=%d chunk_size=%d cores=%d",
              nrow(grid_pairs),
              length(chunks),
              chunk_size,
              ncores_use
            ))
            for (chunk_idx in chunks) {
              chunk_df <- grid_pairs[chunk_idx, , drop = FALSE]
              out_chunk <- try(
                parallel::mclapply(
                  seq_len(nrow(chunk_df)),
                  function(k) score_histology_pair(dx = chunk_df$dX[k], dy = chunk_df$dY[k]),
                  mc.cores = ncores_use,
                  mc.preschedule = TRUE
                ),
                silent = TRUE
              )
              if (inherits(out_chunk, "try-error") || length(out_chunk) != length(chunk_idx)) {
                parallel_failed <- TRUE
                message("[Histology Image Fit] Parallel grid scoring failed; falling back to serial: ", if (inherits(out_chunk, "try-error")) as.character(out_chunk) else "unexpected result length")
                break
              }
              res_list[chunk_idx] <- out_chunk
              if (progress_weight > 0) incProgress(progress_weight * length(chunk_idx) / nrow(grid_pairs))
            }
            if (isTRUE(parallel_failed)) {
              rem <- which(vapply(res_list, is.null, logical(1)))
              step_update_local <- max(1L, as.integer(nrow(grid_pairs) / 40L))
              for (ii in rem) {
                res_list[[ii]] <- score_histology_pair(dx = grid_pairs$dX[ii], dy = grid_pairs$dY[ii])
                if (progress_weight > 0 && (ii %% step_update_local == 0L)) {
                  incProgress(progress_weight / ceiling(nrow(grid_pairs) / step_update_local))
                }
              }
            }
          } else {
            step_update_local <- max(1L, as.integer(nrow(grid_pairs) / 40L))
            for (ii in seq_len(nrow(grid_pairs))) {
              res_list[[ii]] <- score_histology_pair(dx = grid_pairs$dX[ii], dy = grid_pairs$dY[ii])
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
              corr_score = as.numeric(z$corr_score),
              nmi_score = as.numeric(z$nmi_score),
              edge_score = as.numeric(z$edge_score),
              overlap = as.integer(z$overlap),
              n_edge = as.integer(z$n_edge),
              stringsAsFactors = FALSE
            )
          }))
          if (is.null(out) || nrow(out) == 0L) out <- empty_df
          out
        }

        collapse_histology_candidates <- function(df, rng_local, stp_local) {
          if (is.null(df) || nrow(df) == 0L) return(data.frame())
          out <- df[is.finite(df$raw_score) | is.finite(df$edge_score), , drop = FALSE]
          if (nrow(out) == 0L) return(out)
          out <- out[order(out$dX, out$dY), , drop = FALSE]
          out <- out[!duplicated(out[, c("dX", "dY"), drop = FALSE]), , drop = FALSE]
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

          out$shift_px <- sqrt(out$dX^2 + out$dY^2)
          out$dist_norm <- sqrt((out$dX / rng_used)^2 + (out$dY / rng_used)^2)
          out$dist_norm <- pmax(0, pmin(out$dist_norm, 2))
          out$edge_margin <- pmin(rng_used - abs(out$dX), rng_used - abs(out$dY))
          out$edge_norm <- pmax(0, pmin(1, out$edge_margin / edge_band))
          # Histology-fit surfaces often have broad ridges. Use a stronger
          # displacement prior so small score fluctuations do not pull the fit
          # far from a plausible starting alignment.
          if (identical(target_type, "binary_mask")) {
            pen_shift <- 0.40 * out$dist_norm + 0.60 * (out$dist_norm^2)
            pen_edge <- 0.10 * ((1 - out$edge_norm)^2)
          } else {
            pen_shift <- 0.18 * out$dist_norm + 0.18 * (out$dist_norm^2)
            pen_edge <- 0.18 * ((1 - out$edge_norm)^2)
          }

          z_raw <- robust_z(out$raw_score)
          z_edge <- robust_z(out$edge_score)
          has_raw <- any(is.finite(out$raw_score))
          has_edge <- any(is.finite(out$edge_score))
          if (identical(target_type, "binary_mask")) {
            w_raw <- if (has_raw) 1 else 0
            w_edge <- 0
          } else {
            w_raw <- if (has_raw && has_edge) 0.45 else if (has_raw) 1 else 0
            w_edge <- if (has_raw && has_edge) 0.55 else if (has_edge) 1 else 0
          }
          out$score <- (w_raw * z_raw) + (w_edge * z_edge) - pen_shift - pen_edge
          out
        }

        pick_best_histology_candidate <- function(df, rng_local, stp_local) {
          if (is.null(df) || nrow(df) == 0L) return(NULL)
          ord <- order(df$score, decreasing = TRUE)
          df_ord <- df[ord, , drop = FALSE]
          best <- df_ord[1, , drop = FALSE]
          score_vals <- suppressWarnings(as.numeric(df_ord$score))
          score_vals <- score_vals[is.finite(score_vals)]
          score_mad <- if (length(score_vals) > 1L) {
            suppressWarnings(stats::mad(score_vals, center = stats::median(score_vals), constant = 1, na.rm = TRUE))
          } else {
            NA_real_
          }

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

          # Stability tie-break: if several candidates are effectively tied,
          # prefer the smallest movement from the current translation.
          tie_tol <- max(
            0.08,
            if (is.finite(score_mad) && score_mad > 0) min(0.50, 0.35 * score_mad) else 0
          )
          near_best <- df_ord[is.finite(df_ord$score) & df_ord$score >= (best$score[1] - tie_tol), , drop = FALSE]
          if (nrow(near_best) > 1L) {
            near_best$shift_px <- sqrt(near_best$dX^2 + near_best$dY^2)
            near_best <- near_best[order(near_best$shift_px, -near_best$score), , drop = FALSE]
            if (is.finite(near_best$score[1])) {
              best <- near_best[1, , drop = FALSE]
            }
          }

          if (identical(target_type, "binary_mask")) {
            local_radius <- max(4L, 2L * as.integer(stp_local))
            local_df <- df_ord[df_ord$shift_px <= local_radius, , drop = FALSE]
            if (nrow(local_df) > 0L) {
              local_best <- local_df[which.max(local_df$score), , drop = FALSE]
              shift_best <- sqrt(best$dX[1]^2 + best$dY[1]^2)
              local_tol <- max(
                0.45,
                if (is.finite(score_mad) && score_mad > 0) min(1.10, 0.90 * score_mad) else 0
              )
              if (is.finite(shift_best) && shift_best > local_radius) {
                local_tol <- local_tol + 0.02 * (shift_best - local_radius)
              }
              if (is.finite(local_best$score[1]) && local_best$score[1] >= (best$score[1] - local_tol)) {
                best <- local_best
              }
            }
          }
          best
        }

        dx_vals <- seq.int(-rng, rng, by = stp)
        dy_vals <- seq.int(-rng, rng, by = stp)
        message(sprintf(
          "[Histology Image Fit] Grid scoring mode=%s cores=%d",
          if (can_parallel_histology) "parallel" else "serial",
          ncores_use
        ))
        center_probe <- score_histology_pair(dx = 0L, dy = 0L)
        message(sprintf(
          "[Histology Image Fit] Center probe dX=0 dY=0 raw=%s corr=%s nmi=%s edge=%s overlap=%d n_edge=%d",
          if (is.finite(center_probe$raw_score)) sprintf("%.4f", center_probe$raw_score) else "NA",
          if (is.finite(center_probe$corr_score)) sprintf("%.4f", center_probe$corr_score) else "NA",
          if (is.finite(center_probe$nmi_score)) sprintf("%.4f", center_probe$nmi_score) else "NA",
          if (is.finite(center_probe$edge_score)) sprintf("%.4f", center_probe$edge_score) else "NA",
          as.integer(center_probe$overlap %||% 0L),
          as.integer(center_probe$n_edge %||% 0L)
        ))
        coarse_df <- score_histology_grid(dx_vals, dy_vals, progress_weight = 0.48, detail = "Scoring coarse XY grid", allow_parallel = TRUE)
        coarse_finite_n <- if (is.null(coarse_df) || nrow(coarse_df) == 0L) {
          0L
        } else {
          as.integer(sum(is.finite(coarse_df$raw_score) | is.finite(coarse_df$edge_score), na.rm = TRUE))
        }
        message(sprintf(
          "[Histology Image Fit] Coarse grid results total=%d finite_candidates=%d",
          if (is.null(coarse_df)) 0L else nrow(coarse_df),
          coarse_finite_n
        ))
        if (isTRUE(can_parallel_histology) && coarse_finite_n == 0L) {
          message("[Histology Image Fit] Parallel coarse grid produced no finite candidates; retrying serial.")
          coarse_df <- score_histology_grid(dx_vals, dy_vals, progress_weight = 0.48, detail = "Retrying coarse XY grid in serial", allow_parallel = FALSE)
          coarse_finite_n <- if (is.null(coarse_df) || nrow(coarse_df) == 0L) {
            0L
          } else {
            as.integer(sum(is.finite(coarse_df$raw_score) | is.finite(coarse_df$edge_score), na.rm = TRUE))
          }
          message(sprintf(
            "[Histology Image Fit] Serial coarse grid results total=%d finite_candidates=%d",
            if (is.null(coarse_df)) 0L else nrow(coarse_df),
            coarse_finite_n
          ))
        }
        cand_df <- collapse_histology_candidates(coarse_df, rng_local = rng, stp_local = stp)
        message(sprintf(
          "[Histology Image Fit] Collapsed coarse candidates=%d",
          if (is.null(cand_df)) 0L else nrow(cand_df)
        ))
        best_row <- pick_best_histology_candidate(cand_df, rng_local = rng, stp_local = stp)
        if (is.null(best_row) || nrow(best_row) == 0L || !is.finite(best_row$score[1])) {
          xh$histology_fit_grid <- NULL
          xh$histology_fit_candidates <- NULL
          xh$histology_fit_summary <- NULL
          xh$histology_fit_signature <- NULL
          xh$histology_fit_target_info <- NULL
          histology_fit_notify(
            "Histology Fit could not find a valid local translation score. Try a smaller range, different MSI display, or a different histology feature mode.",
            type = "warning",
            duration = 8
          )
          return()
        }

        best_dx <- as.integer(best_row$dX[1])
        best_dy <- as.integer(best_row$dY[1])
        fine_dx <- seq.int(best_dx - stp, best_dx + stp, by = 1L)
        fine_dy <- seq.int(best_dy - stp, best_dy + stp, by = 1L)
        fine_dx <- fine_dx[fine_dx >= -rng & fine_dx <= rng]
        fine_dy <- fine_dy[fine_dy >= -rng & fine_dy <= rng]
        fine_df <- score_histology_grid(fine_dx, fine_dy, progress_weight = 0.14, detail = "Scoring local XY refinement", allow_parallel = TRUE)
        fine_finite_n <- if (is.null(fine_df) || nrow(fine_df) == 0L) {
          0L
        } else {
          as.integer(sum(is.finite(fine_df$raw_score) | is.finite(fine_df$edge_score), na.rm = TRUE))
        }
        message(sprintf(
          "[Histology Image Fit] Fine grid results total=%d finite_candidates=%d",
          if (is.null(fine_df)) 0L else nrow(fine_df),
          fine_finite_n
        ))
        if (isTRUE(can_parallel_histology) && fine_finite_n == 0L && length(fine_dx) > 0L && length(fine_dy) > 0L) {
          message("[Histology Image Fit] Parallel fine grid produced no finite candidates; retrying serial.")
          fine_df <- score_histology_grid(fine_dx, fine_dy, progress_weight = 0.14, detail = "Retrying local XY refinement in serial", allow_parallel = FALSE)
          fine_finite_n <- if (is.null(fine_df) || nrow(fine_df) == 0L) {
            0L
          } else {
            as.integer(sum(is.finite(fine_df$raw_score) | is.finite(fine_df$edge_score), na.rm = TRUE))
          }
          message(sprintf(
            "[Histology Image Fit] Serial fine grid results total=%d finite_candidates=%d",
            if (is.null(fine_df)) 0L else nrow(fine_df),
            fine_finite_n
          ))
        }
        all_df <- if (nrow(fine_df) > 0L) rbind(coarse_df, fine_df) else coarse_df
        cand_df <- collapse_histology_candidates(all_df, rng_local = rng, stp_local = stp)
        message(sprintf(
          "[Histology Image Fit] Collapsed all candidates=%d",
          if (is.null(cand_df)) 0L else nrow(cand_df)
        ))
        best_row <- pick_best_histology_candidate(cand_df, rng_local = rng, stp_local = stp)
        if (is.null(best_row) || nrow(best_row) == 0L || !is.finite(best_row$score[1])) {
          xh$histology_fit_grid <- NULL
          xh$histology_fit_candidates <- NULL
          xh$histology_fit_summary <- NULL
          xh$histology_fit_signature <- NULL
          xh$histology_fit_target_info <- NULL
          histology_fit_notify(
            "Histology Fit could not refine a valid translation score. Try changing the signal source or search range.",
            type = "warning",
            duration = 8
          )
          return()
        }

        cand_df$translate_x <- clamp(tx0 + cand_df$dX, -1000, 1000)
        cand_df$translate_y <- clamp(ty0 + cand_df$dY, -1000, 1000)
        xh$histology_fit_grid <- cand_df
        start_row <- cand_df[cand_df$dX == 0L & cand_df$dY == 0L, , drop = FALSE]
        if (nrow(start_row) == 0L) {
          start_row <- data.frame(
            score = NA_real_,
            raw_score = NA_real_,
            corr_score = NA_real_,
            nmi_score = NA_real_,
            edge_score = NA_real_,
            overlap = NA_integer_,
            n_edge = NA_integer_,
            stringsAsFactors = FALSE
          )
        } else {
          start_row <- start_row[1, , drop = FALSE]
        }

        if (is.finite(start_row$score[1]) && is.finite(best_row$score[1])) {
          cand_scores <- suppressWarnings(as.numeric(cand_df$score))
          cand_scores <- cand_scores[is.finite(cand_scores)]
          score_mad <- if (length(cand_scores) > 1L) {
            suppressWarnings(stats::mad(cand_scores, center = stats::median(cand_scores), constant = 1, na.rm = TRUE))
          } else {
            NA_real_
          }
          gain_tol <- max(
            0.12,
            if (is.finite(score_mad) && score_mad > 0) min(0.60, 0.50 * score_mad) else 0
          )
          best_shift_px <- sqrt(best_row$dX[1]^2 + best_row$dY[1]^2)
          if (is.finite(best_shift_px) && best_shift_px > max(2, stp)) {
            gain_tol <- gain_tol + 0.015 * (best_shift_px - max(2, stp))
            if (identical(target_type, "binary_mask")) {
              gain_tol <- gain_tol + 0.030 * (best_shift_px - max(2, stp))
            }
          }
          score_gain_now <- best_row$score[1] - start_row$score[1]
          if (identical(target_type, "binary_mask")) {
            raw_gain_now <- best_row$raw_score[1] - start_row$raw_score[1]
            raw_gain_tol <- 0.008 + 0.0015 * max(0, best_shift_px - max(1, stp))
            if (is.finite(best_shift_px) &&
              best_shift_px > max(1, stp) &&
              is.finite(raw_gain_now) &&
              raw_gain_now < raw_gain_tol) {
              message(sprintf(
                "[Histology Image Fit] Retaining start translation for binary mask: raw gain %.4f below threshold %.4f for shift %.2f px",
                raw_gain_now, raw_gain_tol, best_shift_px
              ))
              best_row <- start_row
            }
          }
          if (is.finite(best_shift_px) &&
            best_shift_px > max(2, stp) &&
            is.finite(score_gain_now) &&
            score_gain_now < gain_tol) {
            message(sprintf(
              "[Histology Image Fit] Retaining start translation: gain %.4f below stability threshold %.4f for shift %.2f px",
              score_gain_now, gain_tol, best_shift_px
            ))
            best_row <- start_row
          }
        }

        best_idx <- which(cand_df$dX == best_row$dX[1] & cand_df$dY == best_row$dY[1])
        ord_hist <- order(cand_df$score, decreasing = TRUE)
        if (length(best_idx) > 0L) {
          ord_hist <- c(best_idx[1], ord_hist[ord_hist != best_idx[1]])
        }
        score_ord <- order(cand_df$score, decreasing = TRUE)
        start_idx <- which(cand_df$dX == 0L & cand_df$dY == 0L)
        display_idx <- unique(c(
          if (length(best_idx) > 0L) best_idx[1] else integer(0),
          score_ord[1],
          if (length(start_idx) > 0L) start_idx[1] else integer(0),
          ord_hist
        ))
        cand_top <- cand_df[display_idx, , drop = FALSE]
        if (nrow(cand_top) > 6L) cand_top <- cand_top[seq_len(6L), , drop = FALSE]
        cand_top$rank <- seq_len(nrow(cand_top))
        score_rank_map <- integer(nrow(cand_df))
        score_rank_map[score_ord] <- seq_along(score_ord)
        cand_top$score_rank <- score_rank_map[match(paste(cand_top$dX, cand_top$dY), paste(cand_df$dX, cand_df$dY))]
        cand_top$selected <- (cand_top$dX == best_row$dX[1]) & (cand_top$dY == best_row$dY[1])
        cand_top$peak_score <- cand_top$score_rank == 1L
        cand_top$start_position <- (cand_top$dX == 0L) & (cand_top$dY == 0L)
        cand_top_out <- cand_top[, c("rank", "score_rank", "selected", "peak_score", "start_position", "dX", "dY", "score", "raw_score", "corr_score", "nmi_score", "edge_score", "overlap", "n_edge", "translate_x", "translate_y"), drop = FALSE]
        names(cand_top_out)[names(cand_top_out) == "raw_score"] <- "intensity_score"
        names(cand_top_out)[names(cand_top_out) == "corr_score"] <- "signed_corr"
        names(cand_top_out)[names(cand_top_out) == "nmi_score"] <- "nmi"
        xh$histology_fit_candidates <- cand_top_out

        tx_new <- clamp(tx0 + best_row$dX[1], -1000, 1000)
        ty_new <- clamp(ty0 + best_row$dY[1], -1000, 1000)
        peak_row <- cand_df[score_ord[1], , drop = FALSE]

        xh$histology_fit_summary <- list(
          intensity_relation = histology_fit_relation,
          intensity_relation_label = histology_intensity_relation_label(histology_fit_relation),
          feature_mode = feat_mode,
          signal_source = if (!is.null(sig_info$source)) sig_info$source else "current",
          signal_field = if (!is.null(sig_info$field)) sig_info$field else NA_character_,
          signal_type = if (!is.null(sig_info$type)) sig_info$type else NA_character_,
          target_signal = if (!is.null(sig_info$label)) sig_info$label else if (!is.null(msi$display_label)) msi$display_label else NA_character_,
          target_signal_detail = if (!is.null(sig_info$detail)) sig_info$detail else NA_character_,
          target_feature_count = if (!is.null(sig_info$feature_count)) sig_info$feature_count else NA_integer_,
          msi_display_label = if (!is.null(msi$display_label)) msi$display_label else NA_character_,
          scoring_window = fit_window$label,
          search_range = rng,
          search_step = stp,
          start_translate_x = tx0,
          start_translate_y = ty0,
          start_score = start_row$score[1],
          start_intensity_score = start_row$raw_score[1],
          start_signed_correlation = start_row$corr_score[1],
          start_nmi = start_row$nmi_score[1],
          start_edge_score = start_row$edge_score[1],
          best_translate_x = tx_new,
          best_translate_y = ty_new,
          best_dX = best_row$dX[1],
          best_dY = best_row$dY[1],
          best_score = best_row$score[1],
          best_intensity_score = best_row$raw_score[1],
          best_signed_correlation = best_row$corr_score[1],
          best_nmi = best_row$nmi_score[1],
          best_edge_score = best_row$edge_score[1],
          score_gain = best_row$score[1] - start_row$score[1],
          peak_dX = peak_row$dX[1],
          peak_dY = peak_row$dY[1],
          peak_score = peak_row$score[1],
          overlap_pixels = best_row$overlap[1],
          histology_edge_pixels = best_row$n_edge[1],
          n_candidates = nrow(all_df)
        )
        xh$histology_fit_signature <- current_histology_fit_signature()
        message(sprintf(
          "[Histology Image Fit] Best candidate relation=%s target=%s dX=%.2f dY=%.2f -> X=%.2f Y=%.2f | intensity=%.4f corr=%.4f NMI=%.4f edge=%.4f score=%.4f",
          histology_intensity_relation_label(histology_fit_relation),
          if (!is.null(xh$histology_fit_summary$target_signal) && nzchar(as.character(xh$histology_fit_summary$target_signal))) {
            as.character(xh$histology_fit_summary$target_signal)
          } else {
            "n/a"
          },
          as.numeric(xh$histology_fit_summary$best_dX),
          as.numeric(xh$histology_fit_summary$best_dY),
          as.numeric(xh$histology_fit_summary$best_translate_x),
          as.numeric(xh$histology_fit_summary$best_translate_y),
          as.numeric(xh$histology_fit_summary$best_intensity_score),
          as.numeric(xh$histology_fit_summary$best_signed_correlation),
          as.numeric(xh$histology_fit_summary$best_nmi),
          as.numeric(xh$histology_fit_summary$best_edge_score),
          as.numeric(xh$histology_fit_summary$best_score)
        ))
        print(xh$histology_fit_summary)
        print(xh$histology_fit_candidates)

        showNotification(
          sprintf(
            "Histology Fit finished (%s relation) for target %s: best dX=%.2f, dY=%.2f (candidate X=%.2f, Y=%.2f; intensity=%.4f, corr=%.4f, edge=%.4f, score=%.4f). Top-5 candidates are available below; use Apply selected to update X/Y.",
            histology_intensity_relation_label(histology_fit_relation),
            if (!is.null(xh$histology_fit_summary$target_signal) && nzchar(as.character(xh$histology_fit_summary$target_signal))) {
              as.character(xh$histology_fit_summary$target_signal)
            } else {
              "n/a"
            },
            best_row$dX[1],
            best_row$dY[1],
            tx_new,
            ty_new,
            best_row$raw_score[1],
            best_row$corr_score[1],
            best_row$edge_score[1],
            best_row$score[1]
          ),
          type = "message",
          duration = 8
        )
      })
    }, ignoreInit = TRUE)

    output$stat_fit_preview_ui <- renderUI({
      cand <- xh$stat_fit_candidates
      if (is.null(cand) || nrow(cand) == 0) {
        return(tags$small("Stat-fit preview: run Stat Fit to generate up to 5 candidates."))
      }
      labels <- sprintf(
        "#%d dX=%d dY=%d | tx=%.1f ty=%.1f | n_sig=%d score=%.2f",
        cand$rank, cand$dX, cand$dY, as.numeric(cand$translate_x), as.numeric(cand$translate_y), cand$n_sig, cand$score
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
      if (length(pick_rank) != 1L || !is.finite(pick_rank)) pick_rank <- cand$rank[1]
      row <- cand[cand$rank == pick_rank, , drop = FALSE]
      if (nrow(row) == 0) row <- cand[1, , drop = FALSE]

      updateSliderInput(session, "translate_x", value = row$translate_x[1])
      updateSliderInput(session, "translate_y", value = row$translate_y[1])
      updateNumericInput(session, "translate_x_num", value = row$translate_x[1])
      updateNumericInput(session, "translate_y_num", value = row$translate_y[1])

      showNotification(
        sprintf(
          "Applied stat-fit candidate #%d: dX=%d, dY=%d => tx=%.1f, ty=%.1f, n_sig=%d, score=%.2f",
          row$rank[1], row$dX[1], row$dY[1], as.numeric(row$translate_x[1]), as.numeric(row$translate_y[1]), row$n_sig[1], row$score[1]
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

      txy <- current_translate_xy()
      tx <- as.numeric(txy$tx)
      ty <- as.numeric(txy$ty)

      outside_mode <- tolower(trimws(as.character(input$stat_fit_outside_mode)[1]))
      if (!outside_mode %in% c("local", "bbox", "global")) outside_mode <- "bbox"
      stat_metric_mode <- tolower(trimws(as.character(input$stat_fit_metric_mode)[1]))
      if (!stat_metric_mode %in% c("inside_outside", "polygon_cluster_groups")) stat_metric_mode <- "inside_outside"
      stat_group_spec <- parse_stat_fit_group_field()
      stat_objective <- tolower(trimws(as.character(input$stat_fit_objective)[1]))
      if (!stat_objective %in% c("min", "max")) stat_objective <- if (identical(stat_metric_mode, "polygon_cluster_groups")) "max" else "min"
      intensity_transform <- tolower(trimws(as.character(input$intensity_transform)[1]))
      if (!intensity_transform %in% c("none", "sqrt", "log1p", "asinh")) intensity_transform <- "none"

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

      message(sprintf(
        "[Histology Stat Fit] Intensity transform selected: %s",
        intensity_transform
      ))

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

      if (!identical(intensity_transform, "none")) {
        X <- transform_intensity(X, intensity_transform)
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
      cat(sprintf("intensity_transform: %s\n", intensity_transform))
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

      txy <- current_translate_xy()
      tx0 <- as.numeric(txy$tx)
      ty0 <- as.numeric(txy$ty)
      message(sprintf("[Histology Stat Fit] Anchor translation used: tx=%.3f ty=%.3f", tx0, ty0))

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
      intensity_transform <- tolower(trimws(as.character(input$intensity_transform)[1]))
      if (!intensity_transform %in% c("none", "sqrt", "log1p", "asinh")) intensity_transform <- "none"

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

      message(sprintf(
        "[Histology Stat Fit] Intensity transform selected: %s",
        intensity_transform
      ))

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
          if (identical(lab0$group_source, "pdata")) {
            showNotification(
              sprintf(
                "Stat Fit group source is pData (%s); labels are anchored to current translation (dX=0,dY=0).",
                if (!is.null(lab0$group_field)) lab0$group_field else "selected field"
              ),
              type = "message",
              duration = 6
            )
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
          inside_valid0 <- suppressWarnings(mask0[cbind(as.integer(msi$row_idx), as.integer(msi$x_norm))])
          inside_valid0[!is.finite(inside_valid0)] <- FALSE
          frac_in0 <- suppressWarnings(mean(as.logical(inside_valid0), na.rm = TRUE))
          if (!is.finite(frac_in0)) frac_in0 <- NA_real_
          rc <- which(mask0, arr.ind = TRUE)
          if (nrow(rc) > 0) {
            x_min <- min(rc[, "col"], na.rm = TRUE)
            x_max <- max(rc[, "col"], na.rm = TRUE)
            y_plot <- as.integer(msi$ny - rc[, "row"] + 1L)
            y_min <- min(y_plot, na.rm = TRUE)
            y_max <- max(y_plot, na.rm = TRUE)
            message(
              sprintf(
                "[Histology Stat Fit] Anchor mask bbox: x=[%d,%d] y=[%d,%d] inside=%d inside_frac=%.3f",
                x_min, x_max, y_min, y_max, as.integer(n_mask), ifelse(is.finite(frac_in0), frac_in0, NA_real_)
              )
            )
          }
          if (is.finite(frac_in0) && (frac_in0 > 0.90 || frac_in0 < 0.10)) {
            showNotification(
              sprintf(
                "Inside mask covers %.1f%% of MSI pixels at anchor. Stat-fit may be weakly sensitive; consider filtering polygons or using a different metric.",
                100 * frac_in0
              ),
              type = "warning",
              duration = 8
            )
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

        if (!identical(intensity_transform, "none")) {
          X <- transform_intensity(X, intensity_transform)
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
          anchor_translate_x = tx0,
          anchor_translate_y = ty0,
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
          adjusted_p = use_adjusted,
          intensity_transform = intensity_transform
        )

        showNotification(
          sprintf(
            "Stat Fit complete (%s; %s). Best dX=%d, dY=%d => tx=%.1f, ty=%.1f (n_sig=%d, score=%.2f).",
            if (identical(stat_objective, "min")) "min score" else "max score",
            if (identical(stat_metric_mode, "polygon_cluster_groups")) sprintf("%s groups (%s)", if (!is.null(lab0$group_field)) lab0$group_field else stat_group_spec$field, if (!is.null(lab0$group_source)) lab0$group_source else stat_group_spec$source) else "inside/outside",
            best$dX[1], best$dY[1], best$translate_x[1], best$translate_y[1], best$n_sig[1], best$score[1]
          ),
          type = "message",
          duration = 8
        )
      })
    }, ignoreInit = TRUE)

    observeEvent(
      list(
        input$translate_x,
        input$translate_y,
        input$rotate_deg,
        input$scale_x,
        input$scale_y,
        input$overlay_scale_mode,
        input$flip_histology_y,
        input$polygon_axis_mode
      ),
      {
        if (!is.null(xh$stat_fit_grid) || !is.null(xh$stat_fit_summary) || !is.null(xh$stat_fit_candidates)) {
          xh$stat_fit_grid <- NULL
          xh$stat_fit_candidates <- NULL
          xh$stat_fit_summary <- NULL
        }
      },
      ignoreInit = TRUE
    )

    observeEvent(current_histology_fit_signature(), {
      if (is.null(xh$histology_fit_grid) && is.null(xh$histology_fit_candidates) && is.null(xh$histology_fit_summary) && is.null(xh$histology_fit_target_info)) {
        return()
      }
      cur_sig <- current_histology_fit_signature()
      prev_sig <- xh$histology_fit_signature
      if (!is.null(prev_sig) && identical(cur_sig, prev_sig)) {
        return()
      }
      xh$histology_fit_grid <- NULL
      xh$histology_fit_candidates <- NULL
      xh$histology_fit_summary <- NULL
      xh$histology_fit_signature <- NULL
      xh$histology_fit_target_info <- NULL
    }, ignoreInit = TRUE)

    output$histology_fit_heatmap_ui <- renderUI({
      if (requireNamespace("plotly", quietly = TRUE)) {
        plotly::plotlyOutput(ns("histology_fit_heatmap"), height = "240px")
      } else {
        plotOutput(ns("histology_fit_heatmap"), height = "240px")
      }
    })

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

    get_histology_fit_intensity_relation <- function() {
      rel <- NULL
      if (!is.null(xh$histology_fit_summary) && !is.null(xh$histology_fit_summary$intensity_relation)) {
        rel <- tolower(trimws(as.character(xh$histology_fit_summary$intensity_relation)[1]))
      }
      if (!rel %in% c("direct", "inverse", "either")) {
        rel <- resolve_histology_intensity_relation(input$histology_fit_intensity_relation)
      }
      if (!rel %in% c("direct", "inverse", "either")) rel <- "either"
      rel
    }

    build_irregular_axis_bounds <- function(vals) {
      vals <- sort(unique(as.numeric(vals)))
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0L) {
        return(data.frame(value = numeric(0), vmin = numeric(0), vmax = numeric(0)))
      }
      if (length(vals) == 1L) {
        delta <- 0.5
        return(data.frame(value = vals, vmin = vals - delta, vmax = vals + delta))
      }
      mids <- (vals[-1] + vals[-length(vals)]) / 2
      lower_pad <- mids[1] - vals[1]
      upper_pad <- vals[length(vals)] - mids[length(mids)]
      data.frame(
        value = vals,
        vmin = c(vals[1] - lower_pad, mids),
        vmax = c(mids, vals[length(vals)] + upper_pad)
      )
    }

    regularize_histology_fit_plot_grid <- function(g2) {
      if (is.null(g2) || nrow(g2) == 0L) return(g2)
      gx <- sort(unique(as.integer(g2$dX)))
      gy <- sort(unique(as.integer(g2$dY)))
      if (length(gx) == 0L || length(gy) == 0L) return(g2)

      full <- expand.grid(
        dX = gx,
        dY = gy,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
      full <- merge(full, g2, by = c("dX", "dY"), all.x = TRUE, sort = FALSE)
      full$.filled <- !is.finite(full$score)

      have_score <- which(is.finite(full$score))
      need_fill <- which(full$.filled)
      if (length(have_score) > 0L && length(need_fill) > 0L) {
        fill_cols <- intersect(
          c("score", "raw_score", "corr_score", "nmi_score", "edge_score", "overlap", "n_edge", "translate_x", "translate_y"),
          names(full)
        )
        src_x <- full$dX[have_score]
        src_y <- full$dY[have_score]
        for (ii in need_fill) {
          dist2 <- ((src_x - full$dX[ii])^2) + ((src_y - full$dY[ii])^2)
          jj <- have_score[which.min(dist2)]
          full[ii, fill_cols] <- full[jj, fill_cols]
        }
      }

      xb <- build_irregular_axis_bounds(gx)
      yb <- build_irregular_axis_bounds(gy)
      names(xb) <- c("dX", "x_min", "x_max")
      names(yb) <- c("dY", "y_min", "y_max")
      full <- merge(full, xb, by = "dX", all.x = TRUE, sort = FALSE)
      full <- merge(full, yb, by = "dY", all.x = TRUE, sort = FALSE)
      full[order(full$dY, full$dX), , drop = FALSE]
    }

    build_stat_fit_heatmap <- function(g2, best, with_tooltip = FALSE, anchor_tx = NA_real_, anchor_ty = NA_real_) {
      origin <- data.frame(dx = 0, dy = 0)
      anchor_txt <- if (is.finite(anchor_tx) && is.finite(anchor_ty)) {
        sprintf(" | anchor tx=%.1f, ty=%.1f", as.numeric(anchor_tx), as.numeric(anchor_ty))
      } else {
        ""
      }
      if (isTRUE(with_tooltip)) {
        tx_abs <- if (is.finite(anchor_tx)) as.numeric(anchor_tx) + as.numeric(g2$dx) else NA_real_
        ty_abs <- if (is.finite(anchor_ty)) as.numeric(anchor_ty) + as.numeric(g2$dy) else NA_real_
        g2$.tip <- sprintf(
          "dX: %d<br>dY: %d<br>tx: %.1f<br>ty: %.1f<br>n_sig: %d<br>score: %.2f<br>n_inside: %d<br>n_outside: %d",
          as.integer(g2$dx), as.integer(g2$dy), tx_abs, ty_abs, as.integer(g2$n_sig), as.numeric(g2$score),
          as.integer(g2$n_inside), as.integer(g2$n_outside)
        )
        p <- ggplot2::ggplot(g2, ggplot2::aes(x = dx, y = dy, fill = n_sig, text = .tip))
      } else {
        p <- ggplot2::ggplot(g2, ggplot2::aes(x = dx, y = dy, fill = n_sig))
      }

      p +
        ggplot2::geom_tile() +
        ggplot2::geom_point(data = origin, ggplot2::aes(x = dx, y = dy), inherit.aes = FALSE, shape = 3, size = 3, stroke = 1.1, color = "deepskyblue3") +
        ggplot2::geom_point(data = best, ggplot2::aes(x = dx, y = dy), inherit.aes = FALSE, shape = 4, size = 4, stroke = 1.2, color = "red") +
        ggplot2::coord_equal() +
        ggplot2::scale_fill_viridis_c(option = "plasma", na.value = "grey85") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = "Stat-fit heatmap",
          subtitle = paste0("Fill = number of significant ions | blue + = current position (dX=0,dY=0)", anchor_txt),
          x = "dX",
          y = "dY",
          fill = "n_sig"
        )
    }

    if (requireNamespace("plotly", quietly = TRUE)) {
      output$stat_fit_heatmap <- plotly::renderPlotly({
        if (!isTRUE(input$show_fit_info)) return(NULL)
        grid <- xh$stat_fit_grid
        validate(need(!is.null(grid) && nrow(grid) > 0, "Run Stat Fit to view diagnostic heatmap."))
        g2 <- grid[is.finite(grid$n_sig), , drop = FALSE]
        validate(need(nrow(g2) > 0, "No valid Stat Fit points to plot."))
        best <- get_stat_fit_best_row(g2)
        s <- xh$stat_fit_summary
        p <- build_stat_fit_heatmap(
          g2,
          best,
          with_tooltip = TRUE,
          anchor_tx = if (!is.null(s$anchor_translate_x)) as.numeric(s$anchor_translate_x) else NA_real_,
          anchor_ty = if (!is.null(s$anchor_translate_y)) as.numeric(s$anchor_translate_y) else NA_real_
        )
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
        s <- xh$stat_fit_summary
        build_stat_fit_heatmap(
          g2,
          best,
          with_tooltip = FALSE,
          anchor_tx = if (!is.null(s$anchor_translate_x)) as.numeric(s$anchor_translate_x) else NA_real_,
          anchor_ty = if (!is.null(s$anchor_translate_y)) as.numeric(s$anchor_translate_y) else NA_real_
        )
      })
    }

    output$stat_fit_contour <- renderPlot({
      if (!isTRUE(input$show_fit_info)) return(invisible(NULL))
      grid <- xh$stat_fit_grid
      validate(need(!is.null(grid) && nrow(grid) > 0, "Run Stat Fit to view contour diagnostics."))
      g2 <- grid[is.finite(grid$score), , drop = FALSE]
      validate(need(nrow(g2) > 0, "No valid Stat Fit score values to plot."))
      best <- get_stat_fit_best_row(g2)
      origin <- data.frame(dx = 0, dy = 0)
      s <- xh$stat_fit_summary
      anchor_txt <- if (!is.null(s$anchor_translate_x) && !is.null(s$anchor_translate_y) &&
        is.finite(s$anchor_translate_x) && is.finite(s$anchor_translate_y)) {
        sprintf(" | anchor tx=%.1f, ty=%.1f", as.numeric(s$anchor_translate_x), as.numeric(s$anchor_translate_y))
      } else {
        ""
      }
      ggplot2::ggplot(g2, ggplot2::aes(x = dx, y = dy, z = score)) +
        ggplot2::geom_contour(bins = 12, linewidth = 0.5) +
        ggplot2::geom_point(data = origin, ggplot2::aes(x = dx, y = dy), inherit.aes = FALSE, shape = 3, size = 3, stroke = 1.1, color = "deepskyblue3") +
        ggplot2::geom_point(data = best, ggplot2::aes(x = dx, y = dy), shape = 4, size = 4, stroke = 1.2, color = "red") +
        ggplot2::coord_equal() +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = "Stat-fit contours",
          subtitle = paste0("Contours = composite statistical score | blue + = current position (dX=0,dY=0)", anchor_txt),
          x = "dX",
          y = "dY"
        )
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

    get_histology_fit_selected_row <- function(g2, sm = NULL) {
      if (!is.null(sm) &&
        !is.null(sm$best_dX) &&
        !is.null(sm$best_dY) &&
        is.finite(sm$best_dX) &&
        is.finite(sm$best_dY)) {
        hit <- g2[g2$dX == sm$best_dX & g2$dY == sm$best_dY, , drop = FALSE]
        if (nrow(hit) > 0L) return(hit[1, , drop = FALSE])
      }
      if ("selected" %in% names(g2)) {
        hit <- g2[!is.na(g2$selected) & g2$selected, , drop = FALSE]
        if (nrow(hit) > 0L) return(hit[1, , drop = FALSE])
      }
      idx <- which.max(g2$score)
      g2[idx, , drop = FALSE]
    }

    get_histology_fit_peak_row <- function(g2) {
      idx <- which.max(g2$score)
      g2[idx, , drop = FALSE]
    }

    build_histology_fit_heatmap <- function(g2, best, peak = NULL, with_tooltip = FALSE, anchor_tx = NA_real_, anchor_ty = NA_real_) {
      origin <- data.frame(dX = 0, dY = 0)
      rel <- get_histology_fit_intensity_relation()
      g_plot <- regularize_histology_fit_plot_grid(g2)
      filled_note <- if (isTRUE(any(g_plot$.filled))) {
        " | mixed coarse/fine grid completed for display"
      } else {
        ""
      }
      anchor_txt <- if (is.finite(anchor_tx) && is.finite(anchor_ty)) {
        sprintf(" | anchor tx=%.1f, ty=%.1f", as.numeric(anchor_tx), as.numeric(anchor_ty))
      } else {
        ""
      }
      if (isTRUE(with_tooltip)) {
        tx_abs <- if (is.finite(anchor_tx)) as.numeric(anchor_tx) + as.numeric(g_plot$dX) else NA_real_
        ty_abs <- if (is.finite(anchor_ty)) as.numeric(anchor_ty) + as.numeric(g_plot$dY) else NA_real_
        g_plot$.tip <- sprintf(
          "dX: %d<br>dY: %d<br>tx: %.1f<br>ty: %.1f<br>score: %.4f<br>intensity: %.4f<br>corr: %.4f<br>NMI: %.4f<br>edge: %.4f<br>overlap: %d<br>histology edges: %d<br>display cell: %s",
          as.integer(g_plot$dX), as.integer(g_plot$dY), tx_abs, ty_abs, as.numeric(g_plot$score),
          as.numeric(g_plot$raw_score), as.numeric(g_plot$corr_score), as.numeric(g_plot$nmi_score), as.numeric(g_plot$edge_score), as.integer(g_plot$overlap), as.integer(g_plot$n_edge),
          ifelse(g_plot$.filled, "nearest-filled", "scored")
        )
        p <- ggplot2::ggplot(g_plot, ggplot2::aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max, fill = score, text = .tip))
      } else {
        p <- ggplot2::ggplot(g_plot, ggplot2::aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max, fill = score))
      }

      p +
        ggplot2::geom_rect(color = NA) +
        ggplot2::geom_point(data = origin, ggplot2::aes(x = dX, y = dY), inherit.aes = FALSE, shape = 3, size = 3, stroke = 1.1, color = "deepskyblue3") +
        if (!is.null(peak) && nrow(peak) > 0L) ggplot2::geom_point(data = peak, ggplot2::aes(x = dX, y = dY), inherit.aes = FALSE, shape = 1, size = 4.2, stroke = 1.1, color = "orange2") else NULL +
        ggplot2::geom_point(data = best, ggplot2::aes(x = dX, y = dY), inherit.aes = FALSE, shape = 4, size = 4, stroke = 1.2, color = "red") +
        ggplot2::coord_equal() +
        ggplot2::scale_fill_viridis_c(option = "magma", na.value = "grey85") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = "Histology-fit heatmap",
          subtitle = paste0(
            "Fill = composite score | intensity relation: ",
            histology_intensity_relation_label(rel),
            " | blue + = current position (dX=0,dY=0)",
            " | red X = selected fit",
            " | orange circle = raw score peak",
            anchor_txt,
            filled_note
          ),
          x = "dX",
          y = "dY",
          fill = "Score"
        )
    }

    if (requireNamespace("plotly", quietly = TRUE)) {
      output$histology_fit_heatmap <- plotly::renderPlotly({
        if (!isTRUE(input$show_fit_info)) return(NULL)
        grid <- xh$histology_fit_grid
        validate(need(!is.null(grid) && nrow(grid) > 0, "Run Histology Fit to view diagnostic heatmap."))
        g2 <- grid[is.finite(grid$score), , drop = FALSE]
        validate(need(nrow(g2) > 0, "No valid Histology Fit points to plot."))
        s <- xh$histology_fit_summary
        best <- get_histology_fit_selected_row(g2, sm = s)
        peak <- get_histology_fit_peak_row(g2)
        if (nrow(best) > 0L && nrow(peak) > 0L && identical(best$dX[1], peak$dX[1]) && identical(best$dY[1], peak$dY[1])) {
          peak <- NULL
        }
        p <- build_histology_fit_heatmap(
          g2,
          best,
          peak = peak,
          with_tooltip = TRUE,
          anchor_tx = if (!is.null(s$start_translate_x)) as.numeric(s$start_translate_x) else NA_real_,
          anchor_ty = if (!is.null(s$start_translate_y)) as.numeric(s$start_translate_y) else NA_real_
        )
        plotly::ggplotly(p, tooltip = "text")
      })
    } else {
      output$histology_fit_heatmap <- renderPlot({
        if (!isTRUE(input$show_fit_info)) return(invisible(NULL))
        grid <- xh$histology_fit_grid
        validate(need(!is.null(grid) && nrow(grid) > 0, "Run Histology Fit to view diagnostic heatmap."))
        g2 <- grid[is.finite(grid$score), , drop = FALSE]
        validate(need(nrow(g2) > 0, "No valid Histology Fit points to plot."))
        s <- xh$histology_fit_summary
        best <- get_histology_fit_selected_row(g2, sm = s)
        peak <- get_histology_fit_peak_row(g2)
        if (nrow(best) > 0L && nrow(peak) > 0L && identical(best$dX[1], peak$dX[1]) && identical(best$dY[1], peak$dY[1])) {
          peak <- NULL
        }
        build_histology_fit_heatmap(
          g2,
          best,
          peak = peak,
          with_tooltip = FALSE,
          anchor_tx = if (!is.null(s$start_translate_x)) as.numeric(s$start_translate_x) else NA_real_,
          anchor_ty = if (!is.null(s$start_translate_y)) as.numeric(s$start_translate_y) else NA_real_
        )
      })
    }

      output$histology_fit_contour <- renderPlot({
        if (!isTRUE(input$show_fit_info)) return(invisible(NULL))
        grid <- xh$histology_fit_grid
      validate(need(!is.null(grid) && nrow(grid) > 0, "Run Histology Fit to view contour diagnostics."))
      g2 <- grid[is.finite(grid$score), , drop = FALSE]
      validate(need(nrow(g2) > 0, "No valid Histology Fit score values to plot."))
      g_plot <- regularize_histology_fit_plot_grid(g2)
      rel <- get_histology_fit_intensity_relation()
      origin <- data.frame(dX = 0, dY = 0)
        s <- xh$histology_fit_summary
        best <- get_histology_fit_selected_row(g2, sm = s)
        peak <- get_histology_fit_peak_row(g2)
        if (nrow(best) > 0L && nrow(peak) > 0L && identical(best$dX[1], peak$dX[1]) && identical(best$dY[1], peak$dY[1])) {
          peak <- NULL
        }
        anchor_txt <- if (!is.null(s$start_translate_x) && !is.null(s$start_translate_y) &&
          is.finite(s$start_translate_x) && is.finite(s$start_translate_y)) {
          sprintf(" | anchor tx=%.1f, ty=%.1f", as.numeric(s$start_translate_x), as.numeric(s$start_translate_y))
      } else {
        ""
        }
        ggplot2::ggplot(g_plot, ggplot2::aes(x = dX, y = dY, z = score)) +
          ggplot2::geom_contour(bins = 12, linewidth = 0.5) +
          ggplot2::geom_point(data = origin, ggplot2::aes(x = dX, y = dY), inherit.aes = FALSE, shape = 3, size = 3, stroke = 1.1, color = "deepskyblue3") +
          if (!is.null(peak) && nrow(peak) > 0L) ggplot2::geom_point(data = peak, ggplot2::aes(x = dX, y = dY), inherit.aes = FALSE, shape = 1, size = 4.2, stroke = 1.1, color = "orange2") else NULL +
          ggplot2::geom_point(data = best, ggplot2::aes(x = dX, y = dY), inherit.aes = FALSE, shape = 4, size = 4, stroke = 1.2, color = "red") +
          ggplot2::coord_equal() +
          ggplot2::theme_minimal() +
          ggplot2::labs(
          title = "Histology-fit contours",
          subtitle = paste0(
            "Contours = composite score | intensity relation: ",
            histology_intensity_relation_label(rel),
            " | blue + = current position (dX=0,dY=0)",
            " | red X = selected fit",
            " | orange circle = raw score peak",
            anchor_txt
          ),
            x = "dX",
            y = "dY"
          )
      })

    observeEvent(input$optimize_xy, {
      req(msi_for_pdata())
      if (is.null(input$polygon_file) || !nzchar(input$polygon_file$name)) {
        showNotification("Load a polygon file before running Auto-fit XY.", type = "warning", duration = 7)
        return()
      }

      txy <- current_translate_xy()
      tx0 <- as.numeric(txy$tx)
      ty0 <- as.numeric(txy$ty)

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

        msi <- make_msi_raster()
        obj <- msi_for_pdata()
        pd <- as.data.frame(Cardinal::pData(obj))
        validate(need(nrow(pd) == length(mapped_lab), "Length mismatch between pData and mapped clusters."))
        mapped_lab <- flip_values_vertically_to_match_display(mapped_lab, msi)
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
      poly_t <- transform_polygons_for_current_registration(poly, tx = input$translate_x, ty = input$translate_y, axis_mode = axis_mode)$display

      pts_df <- data.frame(
        x = msi$x_display,
        y = msi$y_display,
        idx = seq_along(msi$x_display)
      )
      poly_crs <- normalize_crs(try(sf::st_crs(poly_t), silent = TRUE))
      if (is.null(poly_crs)) {
        pts_sf <- sf::st_as_sf(pts_df, coords = c("x", "y"), remove = FALSE)
      } else {
        pts_sf <- sf::st_as_sf(pts_df, coords = c("x", "y"), crs = poly_crs, remove = FALSE)
      }
      hit <- sf::st_intersects(pts_sf, poly_t)

      outside_label <- "outside_polygon"
      poly_uid <- paste0("polygon_", sprintf("%05d", seq_len(nrow(poly_t))))
      has_hit <- lengths(hit) > 0
      raw_labels <- as.character(poly_t$map_label)
      is_cell_poly <- is_cell_like_polygon_label(raw_labels)
      # Guard against full-canvas/non-cell polygons dominating intersections.
      poly_area <- suppressWarnings(as.numeric(sf::st_area(poly_t)))
      poly_area[!is.finite(poly_area)] <- 0
      canvas_area <- as.numeric(msi$nx) * as.numeric(msi$ny)
      huge_poly <- is.finite(poly_area) & poly_area > (0.5 * canvas_area)
      usable_poly <- !huge_poly
      if (!any(usable_poly)) {
        usable_poly <- rep(TRUE, length(raw_labels))
      }
      if (any(huge_poly) && any(is_cell_poly & usable_poly)) {
        is_cell_poly <- is_cell_poly & usable_poly
      }
      if (!any(is_cell_poly) && any(usable_poly)) {
        is_cell_poly <- usable_poly
      }
      if (!any(is_cell_poly)) {
        # Fallback: if heuristics remove everything, use all polygons.
        is_cell_poly <- rep(TRUE, length(raw_labels))
      }

      map_polygon_values_from_hit <- function(value_vec_full, mask = rep(TRUE, length(value_vec_full)), outside = outside_label) {
        out <- rep(outside, nrow(pts_df))
        if (!any(has_hit)) return(out)
        mask <- rep_len(as.logical(mask), length(value_vec_full))
        mask[is.na(mask)] <- FALSE
        if (identical(input$polygon_overlap_rule, "all")) {
          out[has_hit] <- vapply(hit[has_hit], function(ix) {
            ix_use <- ix[mask[ix]]
            if (length(ix_use) == 0L) return(outside)
            paste(unique(as.character(value_vec_full[ix_use])), collapse = ";")
          }, character(1))
        } else {
          out[has_hit] <- vapply(hit[has_hit], function(ix) {
            ix_use <- ix[mask[ix]]
            if (length(ix_use) == 0L) return(outside)
            as.character(value_vec_full[ix_use[1]])
          }, character(1))
        }
        out
      }

      mapped_lab_all <- map_polygon_values_from_hit(poly_t$map_label, mask = usable_poly, outside = outside_label)
      mapped_uid_all <- map_polygon_values_from_hit(poly_uid, mask = usable_poly, outside = outside_label)
      mapped_lab <- map_polygon_values_from_hit(poly_t$map_label, mask = is_cell_poly, outside = outside_label)
      mapped_uid <- map_polygon_values_from_hit(poly_uid, mask = is_cell_poly, outside = outside_label)
      hit_cell <- mapped_uid != outside_label
      hit_non_cell <- mapped_uid_all != outside_label & !hit_cell

      mapped_lab_original <- NULL
      orig_col_name <- NULL
      if (!is.null(companion_base_field) && "map_label_original" %in% colnames(poly_t)) {
        mapped_lab_original <- map_polygon_values_from_hit(poly_t$map_label_original, mask = is_cell_poly, outside = outside_label)
        orig_col_name <- "polygon_region"
        if (identical(orig_col_name, col_name)) {
          orig_col_name <- make.names(paste0(col_name, "_original"))
        }
      }

      mapped_compartment_label <- NULL
      mapped_compartment_class <- NULL
      mapped_nucleus_uid <- NULL
      compartment_label_col <- make.names(paste0(col_name, "_compartment"))
      compartment_class_col <- make.names(paste0(col_name, "_compartment_class"))
      nucleus_uid_col <- "polygon_nucleus_id"
      n_nucleus_pixels <- 0L
      n_cytoplasm_pixels <- 0L
      n_nucleus_matched <- 0L

      if (!is.null(input$nucleus_polygon_file) && nzchar(input$nucleus_polygon_file$name)) {
        nucleus_poly_try <- try(nucleus_polygon_data(), silent = TRUE)
        if (inherits(nucleus_poly_try, "try-error") || is.null(nucleus_poly_try)) {
          showNotification("Nucleus polygon file could not be read. Remove it or provide a valid GeoJSON.", type = "error", duration = 8)
          return()
        }

        nucleus_poly <- nucleus_poly_try
        if (nrow(nucleus_poly) > 0) {
          nucleus_axis_mode <- resolve_polygon_axis_mode(nucleus_poly, input$polygon_axis_mode, get_histology_image_optional())
          nucleus_t <- transform_polygons_for_current_registration(
            nucleus_poly,
            tx = input$translate_x,
            ty = input$translate_y,
            axis_mode = nucleus_axis_mode
          )$display

          nucleus_uid <- paste0("nucleus_", sprintf("%05d", seq_len(nrow(nucleus_t))))
          nucleus_area <- suppressWarnings(as.numeric(sf::st_area(nucleus_t)))
          nucleus_keep <- !(is.finite(nucleus_area) & nucleus_area > (0.5 * canvas_area))
          if (!any(nucleus_keep)) nucleus_keep <- rep(TRUE, nrow(nucleus_t))
          nucleus_t <- nucleus_t[nucleus_keep, , drop = FALSE]
          nucleus_uid <- nucleus_uid[nucleus_keep]

          if (nrow(nucleus_t) > 0) {
            nucleus_parent_idx <- match_nucleus_polygons_to_cells(poly_t, nucleus_t, cell_keep = is_cell_poly)
            valid_nucleus <- is.finite(nucleus_parent_idx) & nucleus_parent_idx >= 1L & nucleus_parent_idx <= nrow(poly_t)
            n_nucleus_matched <- sum(valid_nucleus, na.rm = TRUE)

            nucleus_parent_uid <- rep(NA_character_, length(nucleus_parent_idx))
            nucleus_parent_lab <- rep(NA_character_, length(nucleus_parent_idx))
            if (any(valid_nucleus)) {
              nucleus_parent_uid[valid_nucleus] <- poly_uid[nucleus_parent_idx[valid_nucleus]]
              nucleus_parent_lab[valid_nucleus] <- as.character(poly_t$map_label[nucleus_parent_idx[valid_nucleus]])
            }

            nucleus_hit <- sf::st_intersects(pts_sf, nucleus_t)
            has_nucleus_hit <- lengths(nucleus_hit) > 0

            map_nucleus_values_from_hit <- function(value_vec_full, outside = outside_label) {
              out <- rep(outside, nrow(pts_df))
              if (!any(has_nucleus_hit)) return(out)
              if (identical(input$polygon_overlap_rule, "all")) {
                out[has_nucleus_hit] <- vapply(nucleus_hit[has_nucleus_hit], function(ix) {
                  ix_use <- ix[valid_nucleus[ix]]
                  if (length(ix_use) == 0L) return(outside)
                  vals <- unique(as.character(value_vec_full[ix_use]))
                  vals <- vals[!is.na(vals) & nzchar(trimws(vals))]
                  if (length(vals) == 0L) return(outside)
                  paste(vals, collapse = ";")
                }, character(1))
              } else {
                out[has_nucleus_hit] <- vapply(nucleus_hit[has_nucleus_hit], function(ix) {
                  ix_use <- ix[valid_nucleus[ix]]
                  if (length(ix_use) == 0L) return(outside)
                  as.character(value_vec_full[ix_use[1]])
                }, character(1))
              }
              out
            }

            mapped_nucleus_uid <- map_nucleus_values_from_hit(nucleus_uid, outside = outside_label)
            mapped_nucleus_parent_uid <- map_nucleus_values_from_hit(nucleus_parent_uid, outside = outside_label)
            mapped_nucleus_parent_lab <- map_nucleus_values_from_hit(nucleus_parent_lab, outside = outside_label)
            nucleus_pixel_hit <- mapped_nucleus_uid != outside_label

            if (any(nucleus_pixel_hit)) {
              mapped_uid[nucleus_pixel_hit] <- mapped_nucleus_parent_uid[nucleus_pixel_hit]
              mapped_lab[nucleus_pixel_hit] <- mapped_nucleus_parent_lab[nucleus_pixel_hit]
            }

            mapped_compartment_label <- mapped_lab_all
            mapped_compartment_class <- rep(outside_label, nrow(pts_df))
            mapped_compartment_class[hit_non_cell] <- "non_cell_polygon"
            mapped_compartment_class[hit_cell] <- "cytoplasm"
            mapped_compartment_label[hit_cell] <- paste0(mapped_lab[hit_cell], "__cytoplasm")
            if (any(nucleus_pixel_hit)) {
              mapped_compartment_class[nucleus_pixel_hit] <- "nucleus"
              mapped_compartment_label[nucleus_pixel_hit] <- paste0(mapped_nucleus_parent_lab[nucleus_pixel_hit], "__nucleus")
            }

            n_nucleus_pixels <- sum(mapped_compartment_class == "nucleus", na.rm = TRUE)
            n_cytoplasm_pixels <- sum(mapped_compartment_class == "cytoplasm", na.rm = TRUE)
          }
        }
      }

      hit_cell <- mapped_uid != outside_label

      obj <- msi_for_pdata()
      pd <- as.data.frame(Cardinal::pData(obj))
      validate(need(nrow(pd) == length(mapped_lab), "Length mismatch between pData and polygon mapping."))
      mapped_lab <- flip_values_vertically_to_match_display(mapped_lab, msi)
      mapped_uid <- flip_values_vertically_to_match_display(mapped_uid, msi)
      hit_cell <- mapped_uid != outside_label
      if (!is.null(mapped_lab_original) && length(mapped_lab_original) == length(mapped_lab)) {
        mapped_lab_original <- flip_values_vertically_to_match_display(mapped_lab_original, msi)
      }
      if (!is.null(mapped_compartment_label) && length(mapped_compartment_label) == length(mapped_lab)) {
        mapped_compartment_label <- flip_values_vertically_to_match_display(mapped_compartment_label, msi)
      }
      if (!is.null(mapped_compartment_class) && length(mapped_compartment_class) == length(mapped_lab)) {
        mapped_compartment_class <- flip_values_vertically_to_match_display(mapped_compartment_class, msi)
      }
      if (!is.null(mapped_nucleus_uid) && length(mapped_nucleus_uid) == length(mapped_lab)) {
        mapped_nucleus_uid <- flip_values_vertically_to_match_display(mapped_nucleus_uid, msi)
      }

      pd[[col_name]] <- mapped_lab
      if (!is.null(mapped_lab_original) && length(mapped_lab_original) == nrow(pd) && !is.null(orig_col_name)) {
        pd[[orig_col_name]] <- mapped_lab_original
      }
      uid_col <- "polygon_cell_id"
      pd[[uid_col]] <- mapped_uid
      binary_col <- "polygon_is_cell"
      pd[[binary_col]] <- as.integer(hit_cell)
      if (!is.null(mapped_compartment_label) && length(mapped_compartment_label) == nrow(pd)) {
        pd[[compartment_label_col]] <- mapped_compartment_label
      }
      if (!is.null(mapped_compartment_class) && length(mapped_compartment_class) == nrow(pd)) {
        pd[[compartment_class_col]] <- mapped_compartment_class
      }
      if (!is.null(mapped_nucleus_uid) && length(mapped_nucleus_uid) == nrow(pd)) {
        pd[[nucleus_uid_col]] <- mapped_nucleus_uid
      }

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
        if (!is.null(mapped_compartment_label) && !is.null(mapped_compartment_class) && !is.null(mapped_nucleus_uid)) {
          showNotification(
            sprintf(
              "Mapped polygons to '%s', original-label column '%s', '%s', '%s', '%s', '%s', and '%s' (cell=%d, nucleus=%d, cytoplasm=%d, matched nuclei=%d).",
              col_name,
              orig_col_name,
              uid_col,
              binary_col,
              compartment_label_col,
              compartment_class_col,
              nucleus_uid_col,
              sum(pd[[binary_col]] == 1, na.rm = TRUE),
              n_nucleus_pixels,
              n_cytoplasm_pixels,
              n_nucleus_matched
            ),
            type = "message",
            duration = 9
          )
        } else {
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
        }
      } else {
        if (!is.null(mapped_compartment_label) && !is.null(mapped_compartment_class) && !is.null(mapped_nucleus_uid)) {
          showNotification(
            sprintf(
              "Mapped polygons to '%s', '%s', '%s', '%s', '%s', and '%s' (cell=%d, nucleus=%d, cytoplasm=%d, matched nuclei=%d).",
              col_name,
              uid_col,
              binary_col,
              compartment_label_col,
              compartment_class_col,
              nucleus_uid_col,
              sum(pd[[binary_col]] == 1, na.rm = TRUE),
              n_nucleus_pixels,
              n_cytoplasm_pixels,
              n_nucleus_matched
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
      graphics::plot.window(xlim = c(1, msi$nx), ylim = c(msi$ny, 1), asp = 1, xaxs = "i", yaxs = "i")
      graphics::rasterImage(msi$raster, 1, 1, msi$nx, msi$ny, interpolate = FALSE)

      draw_polygons <- function(poly) {
        if (is.null(poly) || nrow(poly) == 0L) return(invisible(NULL))
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
        invisible(NULL)
      }

      if (identical(ov$layer, "polygon")) {
        draw_polygons(ov$polygons)
        graphics::title(main = sprintf("Polygon Overlay on MSI (%s)", msi_label))
      } else if (identical(ov$layer, "combined")) {
        if (!is.null(ov$histology)) {
          bbox_disp <- overlay_source_bbox_to_display(ov$histology)
          if (!is.null(bbox_disp)) {
            graphics::rasterImage(
              ov$histology$raster,
              bbox_disp[["xmin"]],
              bbox_disp[["ymin"]],
              bbox_disp[["xmax"]],
              bbox_disp[["ymax"]],
              interpolate = TRUE
            )
          }
        }
        draw_polygons(ov$polygons)
        graphics::title(main = sprintf("Combined Histology + Polygon Overlay on MSI (%s)", msi_label))
      } else {
        bbox_disp <- overlay_source_bbox_to_display(ov)
        validate(need(!is.null(bbox_disp), "Overlay display bounds are unavailable."))
        graphics::rasterImage(ov$raster, bbox_disp[["xmin"]], bbox_disp[["ymin"]], bbox_disp[["xmax"]], bbox_disp[["ymax"]], interpolate = TRUE)
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
      ov_hist <- if (identical(ov$layer, "combined")) ov$histology else ov
      reg_tr <- registration_transform()
      poly_corr <- list(fx = 1, fy = 1, source = "disabled")

      list(
        working_directory = setup_values()[["wd"]],
        overlay_shown = input$overlay_layer,
        frame_status = reg_tr$frame_status,
        slide_key = reg_tr$slide_key,
        overlay_source_frame_type = reg_tr$overlay_source_frame_type,
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
        source_image_dim = if (!is.null(ov_hist$source_width) && !is.null(ov_hist$source_height)) paste0(ov_hist$source_width, " x ", ov_hist$source_height) else NA_character_,
        transformed_overlay_dim = if (!is.null(ov_hist$width) && !is.null(ov_hist$height)) paste0(ov_hist$width, " x ", ov_hist$height) else NA_character_,
        suggested_export_factor_fit_x = if (!is.null(ov_hist$source_width) && is.finite(ov_hist$source_width) && ov_hist$source_width > 0 && is.finite(input$histology_um_per_px) && input$histology_um_per_px > 0 && is.finite(input$msi_um_per_px) && input$msi_um_per_px > 0) {
          msi$nx / (ov_hist$source_width * (input$histology_um_per_px / input$msi_um_per_px))
        } else NA_real_,
        suggested_export_factor_fit_y = if (!is.null(ov_hist$source_height) && is.finite(ov_hist$source_height) && ov_hist$source_height > 0 && is.finite(input$histology_um_per_px) && input$histology_um_per_px > 0 && is.finite(input$msi_um_per_px) && input$msi_um_per_px > 0) {
          msi$ny / (ov_hist$source_height * (input$histology_um_per_px / input$msi_um_per_px))
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

      if (!isTRUE(pdata_values_are_categorical(v))) {
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

    normalize_loaded_rotate_deg <- function(value, kv = NULL, upload_name = NULL) {
      deg <- suppressWarnings(as.numeric(value))
      deg <- if (length(deg) > 0L) deg[1] else NA_real_
      if (!is.finite(deg)) {
        return(list(value = NA_real_, converted = FALSE, reason = "invalid"))
      }

      kv <- if (is.null(kv)) list() else as.list(kv)
      conv <- ""
      if (!is.null(kv[["registration_rotation_convention"]])) {
        conv <- tolower(trimws(as.character(kv[["registration_rotation_convention"]])[1]))
      }
      if (identical(conv, "shared_overlay_rotation_v2")) {
        return(list(value = wrap_rotation_deg(deg), converted = FALSE, reason = "current"))
      }

      saved_layer <- ""
      if (!is.null(kv[["registration_saved_overlay_layer"]])) {
        saved_layer <- tolower(trimws(as.character(kv[["registration_saved_overlay_layer"]])[1]))
      } else if (!is.null(kv[["overlay_layer"]])) {
        saved_layer <- tolower(trimws(as.character(kv[["overlay_layer"]])[1]))
      }

      upload_label <- ""
      if (!is.null(upload_name) && length(upload_name) > 0) {
        upload_label <- tolower(trimws(as.character(upload_name)[1]))
      }

      looks_histology_legacy <- saved_layer %in% c("histology", "cluster") ||
        grepl("histology|hande|h&e", upload_label, perl = TRUE)

      if (isTRUE(looks_histology_legacy)) {
        return(list(value = wrap_rotation_deg(deg - 90), converted = TRUE, reason = "legacy_histology"))
      }

      list(value = wrap_rotation_deg(deg), converted = FALSE, reason = "legacy_default")
    }

    parse_saved_msi_coord_frame <- function(kv = NULL) {
      kv <- if (is.null(kv)) list() else as.list(kv)
      xmin <- suppressWarnings(as.numeric(kv[["msi_coord_xmin"]]))
      xmax <- suppressWarnings(as.numeric(kv[["msi_coord_xmax"]]))
      ymin <- suppressWarnings(as.numeric(kv[["msi_coord_ymin"]]))
      ymax <- suppressWarnings(as.numeric(kv[["msi_coord_ymax"]]))
      xmin <- if (length(xmin) > 0L) xmin[1] else NA_real_
      xmax <- if (length(xmax) > 0L) xmax[1] else NA_real_
      ymin <- if (length(ymin) > 0L) ymin[1] else NA_real_
      ymax <- if (length(ymax) > 0L) ymax[1] else NA_real_
      if (!all(is.finite(c(xmin, xmax, ymin, ymax)))) return(NULL)

      nx <- suppressWarnings(as.integer(kv[["msi_canvas_nx"]]))
      ny <- suppressWarnings(as.integer(kv[["msi_canvas_ny"]]))
      nx <- if (length(nx) > 0L) nx[1] else NA_integer_
      ny <- if (length(ny) > 0L) ny[1] else NA_integer_
      if (!is.finite(nx) || nx < 1L) nx <- suppressWarnings(as.integer(round(xmax - xmin + 1)))
      if (!is.finite(ny) || ny < 1L) ny <- suppressWarnings(as.integer(round(ymax - ymin + 1)))
      nx <- if (length(nx) > 0L) nx[1] else NA_integer_
      ny <- if (length(ny) > 0L) ny[1] else NA_integer_
      if (!is.finite(nx) || nx < 1L || !is.finite(ny) || ny < 1L) return(NULL)

      list(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        nx = nx,
        ny = ny,
        source = "saved_params"
      )
    }

    saved_translation_frame_mode <- function(kv = NULL) {
      kv <- if (is.null(kv)) list() else as.list(kv)
      mode <- ""
      if (!is.null(kv[["registration_translation_frame"]]) && length(kv[["registration_translation_frame"]]) > 0L) {
        mode <- tolower(trimws(as.character(kv[["registration_translation_frame"]])[1]))
      }
      if (!nzchar(mode)) mode <- "legacy_unknown"
      mode
    }

    current_frame_is_saved_subset <- function(current_frame, saved_frame, tol = 1e-8) {
      if (is.null(current_frame) || is.null(saved_frame)) return(FALSE)
      vals <- c(
        current_frame$xmin, current_frame$xmax, current_frame$ymin, current_frame$ymax,
        saved_frame$xmin, saved_frame$xmax, saved_frame$ymin, saved_frame$ymax
      )
      if (!all(is.finite(vals))) return(FALSE)
      inside_x <- current_frame$xmin >= (saved_frame$xmin - tol) &&
        current_frame$xmax <= (saved_frame$xmax + tol)
      inside_y <- current_frame$ymin >= (saved_frame$ymin - tol) &&
        current_frame$ymax <= (saved_frame$ymax + tol)
      smaller_extent <- (current_frame$nx <= saved_frame$nx) && (current_frame$ny <= saved_frame$ny) &&
        ((current_frame$nx < saved_frame$nx) || (current_frame$ny < saved_frame$ny))
      isTRUE(inside_x && inside_y && smaller_extent)
    }

    estimate_saved_overlay_canvas <- function(kv = NULL) {
      kv <- if (is.null(kv)) list() else as.list(kv)
      scale_mode <- if (!is.null(kv[["overlay_scale_mode"]]) && length(kv[["overlay_scale_mode"]]) > 0L) {
        tolower(trimws(as.character(kv[["overlay_scale_mode"]])[1]))
      } else {
        ""
      }
      if (!identical(scale_mode, "absolute")) return(NULL)

      src_dims <- get_overlay_source_dims()
      if (is.null(src_dims) || !all(is.finite(c(src_dims$width, src_dims$height))) ||
          src_dims$width <= 0 || src_dims$height <= 0) {
        return(NULL)
      }

      sx <- suppressWarnings(as.numeric(kv[["scale_x"]]))
      sy <- suppressWarnings(as.numeric(kv[["scale_y"]]))
      sx <- if (length(sx) > 0L) sx[1] else NA_real_
      sy <- if (length(sy) > 0L) sy[1] else NA_real_
      if (!is.finite(sx) || sx <= 0 || !is.finite(sy) || sy <= 0) return(NULL)

      rot_info <- normalize_loaded_rotate_deg(
        value = if (!is.null(kv[["rotate_deg"]])) kv[["rotate_deg"]] else 0,
        kv = kv,
        upload_name = NULL
      )
      rot_deg <- rot_info$value
      if (!is.finite(rot_deg)) rot_deg <- 0
      theta <- (rot_deg + 90) * pi / 180

      base_w <- as.numeric(src_dims$width) * sx
      base_h <- as.numeric(src_dims$height) * sy
      est_nx <- suppressWarnings(as.integer(round(abs(base_w * cos(theta)) + abs(base_h * sin(theta)))))
      est_ny <- suppressWarnings(as.integer(round(abs(base_w * sin(theta)) + abs(base_h * cos(theta)))))
      if (!is.finite(est_nx) || est_nx < 1L || !is.finite(est_ny) || est_ny < 1L) return(NULL)

      list(
        xmin = NA_real_,
        xmax = NA_real_,
        ymin = NA_real_,
        ymax = NA_real_,
        nx = est_nx,
        ny = est_ny,
        source = paste0("overlay_estimate:", src_dims$source)
      )
    }

    parse_saved_overlay_center_source <- function(kv = NULL) {
      kv <- if (is.null(kv)) list() else as.list(kv)
      cx <- suppressWarnings(as.numeric(kv[["overlay_center_source_x"]]))
      cy <- suppressWarnings(as.numeric(kv[["overlay_center_source_y"]]))
      cx <- if (length(cx) > 0L) cx[1] else NA_real_
      cy <- if (length(cy) > 0L) cy[1] else NA_real_
      if (all(is.finite(c(cx, cy)))) {
        return(list(x = cx, y = cy, source = "saved_center"))
      }

      saved_frame <- parse_saved_msi_coord_frame(kv)
      tx_num <- suppressWarnings(as.numeric(kv[["translate_x"]]))
      ty_num <- suppressWarnings(as.numeric(kv[["translate_y"]]))
      tx_num <- if (length(tx_num) > 0L) tx_num[1] else NA_real_
      ty_num <- if (length(ty_num) > 0L) ty_num[1] else NA_real_
      if (!is.null(saved_frame) && all(is.finite(c(saved_frame$xmin, saved_frame$xmax, saved_frame$ymin, saved_frame$ymax, tx_num, ty_num)))) {
        return(list(
          x = ((saved_frame$xmin + saved_frame$xmax) / 2) + tx_num,
          y = ((saved_frame$ymin + saved_frame$ymax) / 2) + ty_num,
          source = "saved_frame_plus_translate"
        ))
      }

      NULL
    }

    adjust_loaded_translate_xy <- function(tx, ty, kv = NULL, upload_name = NULL) {
      tx_num <- suppressWarnings(as.numeric(tx))
      ty_num <- suppressWarnings(as.numeric(ty))
      tx_num <- if (length(tx_num) > 0L) tx_num[1] else NA_real_
      ty_num <- if (length(ty_num) > 0L) ty_num[1] else NA_real_
      if (!is.finite(tx_num) && !is.finite(ty_num)) {
        return(list(tx = tx_num, ty = ty_num, adjusted = FALSE, mode = "none"))
      }

      current_frame <- current_msi_coord_frame()
      if (is.null(current_frame)) {
        return(list(tx = tx_num, ty = ty_num, adjusted = FALSE, mode = "no_current_frame"))
      }

      saved_frame <- parse_saved_msi_coord_frame(kv)
      translation_mode <- saved_translation_frame_mode(kv)
      if (current_frame_is_saved_subset(current_frame, saved_frame)) {
        return(list(
          tx = tx_num,
          ty = ty_num,
          adjusted = FALSE,
          mode = "preserved_subset_coords",
          current_frame = current_frame,
          saved_frame = saved_frame,
          saved_translation_mode = translation_mode
        ))
      }

      saved_center <- parse_saved_overlay_center_source(kv)
      if (!is.null(saved_center) && all(is.finite(c(saved_center$x, saved_center$y)))) {
        current_center_x <- (current_frame$xmin + current_frame$xmax) / 2
        current_center_y <- (current_frame$ymin + current_frame$ymax) / 2
        tx_adj <- saved_center$x - current_center_x
        ty_adj <- saved_center$y - current_center_y
        delta_tx <- if (is.finite(tx_num) && is.finite(tx_adj)) abs(tx_adj - tx_num) else 0
        delta_ty <- if (is.finite(ty_num) && is.finite(ty_adj)) abs(ty_adj - ty_num) else 0
        adjusted <- (delta_tx > 1e-8) || (delta_ty > 1e-8)
        return(list(
          tx = if (is.finite(tx_adj)) clamp(tx_adj, -1000, 1000) else tx_adj,
          ty = if (is.finite(ty_adj)) clamp(ty_adj, -1000, 1000) else ty_adj,
          adjusted = adjusted,
          mode = as.character(saved_center$source),
          current_frame = current_frame,
          saved_frame = parse_saved_msi_coord_frame(kv),
          saved_translation_mode = translation_mode
        ))
      }

      mode <- "saved_frame"
      if (is.null(saved_frame)) {
        saved_frame <- estimate_saved_overlay_canvas(kv)
        mode <- "heuristic_overlay_estimate"
        if (!is.null(saved_frame)) {
          guess_xmin <- if (is.finite(current_frame$xmin) && current_frame$xmin >= 1) 1 else 0
          guess_ymin <- if (is.finite(current_frame$ymin) && current_frame$ymin >= 1) 1 else 0
          saved_frame$xmin <- guess_xmin
          saved_frame$ymin <- guess_ymin
          saved_frame$xmax <- saved_frame$xmin + saved_frame$nx - 1
          saved_frame$ymax <- saved_frame$ymin + saved_frame$ny - 1

          # Only apply the heuristic if the current MSI frame looks like a crop
          # of a larger full-image canvas.
          bigger_canvas <- (saved_frame$nx >= current_frame$nx + 2L) || (saved_frame$ny >= current_frame$ny + 2L)
          cropped_extent <- (current_frame$xmin > saved_frame$xmin + 1) ||
            (current_frame$ymin > saved_frame$ymin + 1) ||
            (current_frame$xmax < saved_frame$xmax - 1) ||
            (current_frame$ymax < saved_frame$ymax - 1)
          if (!(bigger_canvas && cropped_extent)) {
            saved_frame <- NULL
          }
        }
      }

      if (is.null(saved_frame)) {
        return(list(tx = tx_num, ty = ty_num, adjusted = FALSE, mode = "no_saved_frame"))
      }

      tx_adj <- if (is.finite(tx_num)) {
        tx_num + (saved_frame$xmin - current_frame$xmin) + ((saved_frame$nx - current_frame$nx) / 2)
      } else {
        tx_num
      }
      ty_adj <- if (is.finite(ty_num)) {
        ty_num + (saved_frame$ymin - current_frame$ymin) + ((saved_frame$ny - current_frame$ny) / 2)
      } else {
        ty_num
      }

      delta_tx <- if (is.finite(tx_num) && is.finite(tx_adj)) abs(tx_adj - tx_num) else 0
      delta_ty <- if (is.finite(ty_num) && is.finite(ty_adj)) abs(ty_adj - ty_num) else 0
      adjusted <- (delta_tx > 1e-8) || (delta_ty > 1e-8)

      list(
        tx = if (is.finite(tx_adj)) clamp(tx_adj, -1000, 1000) else tx_adj,
        ty = if (is.finite(ty_adj)) clamp(ty_adj, -1000, 1000) else ty_adj,
        adjusted = adjusted,
        mode = mode,
        current_frame = current_frame,
        saved_frame = saved_frame,
        saved_translation_mode = translation_mode
      )
    }

    observeEvent(input$registration_params_upload, {
      req(input$registration_params_upload)
      xh$histology_fit_relation_auto <- TRUE
      xh$histology_fit_relation_last_default <- NULL
      xh$histology_fit_relation_updating <- FALSE
      kv <- read_registration_params(input$registration_params_upload$datapath)
      validate(need(!is.null(kv) && length(kv) > 0, "Could not parse registration parameter file."))
      kv <- as.list(kv)
      xh$loaded_registration_reference_frame <- parse_saved_msi_coord_frame(kv)

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
      } else if (!is.null(kv[["histology_downsample_factor"]])) {
        v <- suppressWarnings(as.numeric(kv[["histology_downsample_factor"]]))
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
        rot_info <- normalize_loaded_rotate_deg(
          value = kv[["rotate_deg"]],
          kv = kv,
          upload_name = input$registration_params_upload$name
        )
        v <- rot_info$value
        if (is.finite(v)) {
          updateSliderInput(session, "rotate_deg", value = v)
          updateNumericInput(session, "rotate_deg_num", value = v)
          if (isTRUE(rot_info$converted)) {
            showNotification(
              "Loaded legacy histology registration parameters with a 90-degree rotation correction for the shared overlay convention.",
              type = "message",
              duration = 7
            )
          }
        }
      }
      if (!is.null(kv[["translate_x"]]) || !is.null(kv[["translate_y"]])) {
        tr_info <- adjust_loaded_translate_xy(
          tx = kv[["translate_x"]],
          ty = kv[["translate_y"]],
          kv = kv,
          upload_name = input$registration_params_upload$name
        )
        if (is.finite(tr_info$tx)) {
          updateSliderInput(session, "translate_x", value = tr_info$tx)
          updateNumericInput(session, "translate_x_num", value = tr_info$tx)
        }
        if (is.finite(tr_info$ty)) {
          updateSliderInput(session, "translate_y", value = tr_info$ty)
          updateNumericInput(session, "translate_y_num", value = tr_info$ty)
        }
        msg <- if (identical(tr_info$mode, "saved_center")) {
          "Adjusted loaded translation to the current cropped MSI coordinate frame using the saved overlay center in MSI source coordinates."
        } else if (identical(tr_info$mode, "saved_frame_plus_translate")) {
          "Adjusted loaded translation to the current cropped MSI coordinate frame using the saved full-image frame and translation."
        } else if (identical(tr_info$mode, "saved_frame")) {
          "Adjusted loaded translation to the current cropped MSI coordinate frame using saved MSI extents."
        } else if (identical(tr_info$mode, "heuristic_overlay_estimate")) {
          "Adjusted loaded translation to the current cropped MSI coordinate frame using an estimated full-image canvas."
        } else if (identical(tr_info$mode, "preserved_subset_coords")) {
          "Preserved loaded translation because the current MSI frame appears to be a coordinate-preserving subset of the saved full-image frame."
        } else if (identical(tr_info$mode, "no_saved_frame")) {
          "Loaded translation without crop rebasing because no saved MSI frame metadata was available."
        } else if (identical(tr_info$mode, "no_current_frame")) {
          "Loaded translation without crop rebasing because the current MSI coordinate frame was unavailable."
        } else {
          sprintf("Loaded translation using mode '%s'.", tr_info$mode %||% "unknown")
        }
        if (!is.null(msg) && isTRUE(tr_info$adjusted)) {
          showNotification(msg, type = "message", duration = 7)
        }
        cur_frame_txt <- if (!is.null(tr_info$current_frame) && all(is.finite(unlist(tr_info$current_frame[c("xmin", "xmax", "ymin", "ymax")]))) ) {
          sprintf("[%.0f, %.0f] x [%.0f, %.0f]", tr_info$current_frame$xmin, tr_info$current_frame$xmax, tr_info$current_frame$ymin, tr_info$current_frame$ymax)
        } else {
          "n/a"
        }
        saved_frame_txt <- if (!is.null(tr_info$saved_frame) && all(is.finite(unlist(tr_info$saved_frame[c("xmin", "xmax", "ymin", "ymax")]))) ) {
          sprintf("[%.0f, %.0f] x [%.0f, %.0f]", tr_info$saved_frame$xmin, tr_info$saved_frame$xmax, tr_info$saved_frame$ymin, tr_info$saved_frame$ymax)
        } else {
          "n/a"
        }
        cat(sprintf(
          "[Registration Params] translate load mode=%s adjusted=%s tx=%.3f ty=%.3f current_frame=%s saved_frame=%s\n",
          tr_info$mode %||% "unknown",
          if (isTRUE(tr_info$adjusted)) "TRUE" else "FALSE",
          if (is.finite(tr_info$tx)) tr_info$tx else NA_real_,
          if (is.finite(tr_info$ty)) tr_info$ty else NA_real_,
          cur_frame_txt,
          saved_frame_txt
        ))
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
      if (!is.null(kv[["gaussian_smooth"]])) updateCheckboxInput(session, "gaussian_smooth", value = to_bool(kv[["gaussian_smooth"]], FALSE))
      if (!is.null(kv[["registration_saved_overlay_layer"]]) || !is.null(kv[["overlay_layer"]])) {
        overlay_val <- if (!is.null(kv[["registration_saved_overlay_layer"]])) kv[["registration_saved_overlay_layer"]] else kv[["overlay_layer"]]
        overlay_val <- tolower(trimws(as.character(overlay_val)[1]))
        if (overlay_val %in% c("combined", "polygon", "histology", "cluster")) {
          updateRadioButtons(session, "overlay_layer", selected = overlay_val)
        }
      }
      if (!is.null(kv[["orientation_override_active"]]) && to_bool(kv[["orientation_override_active"]], FALSE)) {
        updateCheckboxInput(session, "show_advanced_registration", value = TRUE)
      }

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
            if ("histology_fit_range" %in% names(kv) && !is.null(kv[["histology_fit_range"]])) {
              v <- suppressWarnings(as.numeric(kv[["histology_fit_range"]]))
              if (is.finite(v)) updateNumericInput(session, "histology_fit_range", value = clamp(v, 1, 500))
            }
            if ("histology_fit_step" %in% names(kv) && !is.null(kv[["histology_fit_step"]])) {
              v <- suppressWarnings(as.numeric(kv[["histology_fit_step"]]))
              if (is.finite(v)) updateNumericInput(session, "histology_fit_step", value = clamp(v, 1, 50))
            }
            if ("histology_feature_mode" %in% names(kv) && !is.null(kv[["histology_feature_mode"]])) {
              v <- tolower(trimws(as.character(kv[["histology_feature_mode"]])[1]))
              if (v %in% c("hematoxylin", "darkness", "purple")) updateSelectInput(session, "histology_feature_mode", selected = v)
            }
            if ("histology_fit_signal_source" %in% names(kv) && !is.null(kv[["histology_fit_signal_source"]])) {
              v <- tolower(trimws(as.character(kv[["histology_fit_signal_source"]])[1]))
              if (v %in% c("current", "pdata", "multi", "pca")) updateSelectInput(session, "histology_fit_signal_source", selected = v)
            }
            if ("histology_fit_pdata_field" %in% names(kv) && !is.null(kv[["histology_fit_pdata_field"]])) {
              v <- as.character(kv[["histology_fit_pdata_field"]])[1]
              if (!is.na(v) && nzchar(v)) updateSelectInput(session, "histology_fit_pdata_field", selected = v)
            }
            if ("histology_fit_multi_n" %in% names(kv) && !is.null(kv[["histology_fit_multi_n"]])) {
              v <- suppressWarnings(as.numeric(kv[["histology_fit_multi_n"]]))
              if (is.finite(v)) updateNumericInput(session, "histology_fit_multi_n", value = clamp(v, 2, 64))
            }
            if ("histology_fit_pca_component" %in% names(kv) && !is.null(kv[["histology_fit_pca_component"]])) {
              v <- suppressWarnings(as.numeric(kv[["histology_fit_pca_component"]]))
              if (is.finite(v)) updateNumericInput(session, "histology_fit_pca_component", value = clamp(v, 1, 64))
            }
            if ("histology_fit_pca_n_components" %in% names(kv) && !is.null(kv[["histology_fit_pca_n_components"]])) {
              v <- suppressWarnings(as.numeric(kv[["histology_fit_pca_n_components"]]))
              if (is.finite(v)) updateNumericInput(session, "histology_fit_pca_n_components", value = clamp(v, 1, 8))
            }
            if ("histology_fit_pca_n_features" %in% names(kv) && !is.null(kv[["histology_fit_pca_n_features"]])) {
              v <- suppressWarnings(as.numeric(kv[["histology_fit_pca_n_features"]]))
              if (is.finite(v)) updateNumericInput(session, "histology_fit_pca_n_features", value = clamp(v, 4, 128))
            }
            if ("histology_fit_intensity_relation" %in% names(kv) && !is.null(kv[["histology_fit_intensity_relation"]])) {
              v <- tolower(trimws(as.character(kv[["histology_fit_intensity_relation"]])[1]))
              if (v %in% c("direct", "inverse", "either")) updateSelectInput(session, "histology_fit_intensity_relation", selected = v)
            } else if ("histology_fit_objective" %in% names(kv) && !is.null(kv[["histology_fit_objective"]])) {
              v <- tolower(trimws(as.character(kv[["histology_fit_objective"]])[1]))
              v_map <- if (identical(v, "min")) "inverse" else "either"
              updateSelectInput(session, "histology_fit_intensity_relation", selected = v_map)
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
        if ("histology_fit_pdata_field" %in% names(kv) && !is.null(kv[["histology_fit_pdata_field"]])) {
          v <- as.character(kv[["histology_fit_pdata_field"]])[1]
          if (!is.na(v) && nzchar(v)) updateSelectInput(session, "histology_fit_pdata_field", selected = v)
        }
        if ("histology_fit_multi_n" %in% names(kv) && !is.null(kv[["histology_fit_multi_n"]])) {
          v <- suppressWarnings(as.numeric(kv[["histology_fit_multi_n"]]))
          if (is.finite(v)) updateNumericInput(session, "histology_fit_multi_n", value = clamp(v, 2, 64))
        }
        if ("histology_fit_pca_component" %in% names(kv) && !is.null(kv[["histology_fit_pca_component"]])) {
          v <- suppressWarnings(as.numeric(kv[["histology_fit_pca_component"]]))
          if (is.finite(v)) updateNumericInput(session, "histology_fit_pca_component", value = clamp(v, 1, 64))
        }
        if ("histology_fit_pca_n_components" %in% names(kv) && !is.null(kv[["histology_fit_pca_n_components"]])) {
          v <- suppressWarnings(as.numeric(kv[["histology_fit_pca_n_components"]]))
          if (is.finite(v)) updateNumericInput(session, "histology_fit_pca_n_components", value = clamp(v, 1, 8))
        }
        if ("histology_fit_pca_n_features" %in% names(kv) && !is.null(kv[["histology_fit_pca_n_features"]])) {
          v <- suppressWarnings(as.numeric(kv[["histology_fit_pca_n_features"]]))
          if (is.finite(v)) updateNumericInput(session, "histology_fit_pca_n_features", value = clamp(v, 4, 128))
        }
      }, once = TRUE)

      saved_slide_key <- normalize_slide_key(kv[["slide_key"]])
      current_slide_key <- normalize_slide_key(registration_transform()$slide_key)
      if (!is.na(saved_slide_key) && !is.na(current_slide_key) && !identical(saved_slide_key, current_slide_key)) {
        showNotification(
          sprintf("Loaded registration parameters from slide key '%s', but current assets resolve to '%s'. Treating this as best-effort transfer.", saved_slide_key, current_slide_key),
          type = "warning",
          duration = 8
        )
      } else {
        showNotification("Loaded registration parameters.", type = "message", duration = 5)
      }
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
      nudge_translation(dx = 0, dy = move_step())
    }, ignoreInit = TRUE)

    observeEvent(input$move_down, {
      nudge_translation(dx = 0, dy = -move_step())
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
      xh$loaded_registration_reference_frame <- NULL
      xh$histology_fit_relation_auto <- TRUE
      xh$histology_fit_relation_last_default <- NULL
      xh$histology_fit_relation_updating <- FALSE
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
      updateCheckboxInput(session, "show_advanced_registration", value = FALSE)
      updateSelectInput(session, "polygon_axis_mode", selected = "auto")
      updateSelectInput(session, "overlay_scale_mode", selected = "absolute")
      updateRadioButtons(session, "overlay_layer", selected = "combined")
      updateCheckboxInput(session, "fine_move", value = FALSE)
      updateCheckboxInput(session, "polygon_color_by_label", value = FALSE)
      updateCheckboxInput(session, "enhance_contrast", value = TRUE)
      updateCheckboxInput(session, "gaussian_smooth", value = FALSE)
      updateCheckboxInput(session, "show_fit_info", value = TRUE)
      updateSliderInput(session, "polygon_linewidth", value = 1)
      updateSelectInput(session, "cluster_palette", selected = "Alphabet")
            updateSelectInput(session, "rgb_render_mode", selected = "dominant")
            updateSliderInput(session, "rgb_bg_cutoff", value = 0.05)
            updateNumericInput(session, "optimize_edge_band", value = 4)
            updateNumericInput(session, "histology_fit_range", value = 30)
            updateNumericInput(session, "histology_fit_step", value = 2)
            updateSelectInput(session, "histology_feature_mode", selected = "purple")
            updateNumericInput(session, "histology_fit_pca_component", value = 1)
            updateNumericInput(session, "histology_fit_pca_n_components", value = 1)
            updateNumericInput(session, "histology_fit_pca_n_features", value = 32)
            updateSelectInput(session, "histology_fit_signal_source", selected = "multi")
            updateNumericInput(session, "histology_fit_multi_n", value = 12)
            updateSelectInput(session, "histology_fit_intensity_relation", selected = "inverse")
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
