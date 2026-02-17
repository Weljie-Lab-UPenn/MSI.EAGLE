HistologyIntegrationServer <- function(id, setup_values, preproc_values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    xh <- reactiveValues(
      mapped_obj = NULL,
      cluster_lookup = NULL,
      mapped_column = NULL,
      mapping_source = NULL
    )

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

    observeEvent(msi_data(), {
      mzv <- try(Cardinal::mz(msi_data()), silent = TRUE)
      if (inherits(mzv, "try-error") || length(mzv) == 0) {
        return()
      }

      idx_vals <- as.character(seq_along(mzv))
      labels <- format(mzv, digits = 10, scientific = FALSE, trim = TRUE)
      choices <- idx_vals
      names(choices) <- labels

      cur <- isolate(input$mz_select)
      sel <- if (!is.null(cur) && cur %in% idx_vals) cur else idx_vals[1]
      updateSelectizeInput(session, "mz_select", choices = choices, selected = sel, server = TRUE)
    }, ignoreInit = FALSE)

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
      default <- if ("classification" %in% choices) "classification" else choices[1]
      selectInput(ns("polygon_label_field"), "Polygon label field", choices = choices, selected = default)
    })

    transform_intensity <- function(x, trans) {
      if (identical(trans, "sqrt")) return(sqrt(pmax(x, 0)))
      if (identical(trans, "log1p")) return(log1p(pmax(x, 0)))
      if (identical(trans, "asinh")) return(asinh(x))
      x
    }

    clamp <- function(x, lo, hi) {
      max(lo, min(hi, x))
    }

    safe_color <- function(col, fallback = "#73FFFF") {
      col <- as.character(col)[1]
      if (is.null(col) || is.na(col) || !nzchar(trimws(col))) return(fallback)
      ok <- try(grDevices::col2rgb(col), silent = TRUE)
      if (inherits(ok, "try-error")) fallback else col
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
      list(
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
        intensity_transform = input$intensity_transform,
        msi_palette = input$msi_palette,
        cluster_palette = input$cluster_palette,
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
      req(msi_data())
      req(input$mz_select)

      obj <- msi_data()
      cd <- as.data.frame(Cardinal::coord(obj))
      validate(need(all(c("x", "y") %in% names(cd)), "MSI coordinates must contain x and y."))
      validate(need(nrow(cd) > 0, "MSI data has no coordinates."))

      mzv <- Cardinal::mz(obj)
      validate(need(length(mzv) > 0, "MSI data has no m/z axis."))
      idx <- suppressWarnings(as.integer(input$mz_select))
      validate(need(!is.na(idx) && idx >= 1L && idx <= length(mzv), "Please select a valid m/z value from the list."))

      vals <- try(as.numeric(Cardinal::spectra(obj)[idx, ]), silent = TRUE)
      validate(need(!inherits(vals, "try-error"), "Could not extract MSI intensities for selected m/z."))
      vals <- transform_intensity(vals, input$intensity_transform)

      x_norm <- as.integer(cd$x - min(cd$x, na.rm = TRUE) + 1L)
      y_norm <- as.integer(cd$y - min(cd$y, na.rm = TRUE) + 1L)
      nx <- max(x_norm, na.rm = TRUE)
      ny <- max(y_norm, na.rm = TRUE)

      mat <- matrix(NA_real_, nrow = ny, ncol = nx)
      row_idx <- ny - y_norm + 1L
      mat[cbind(row_idx, x_norm)] <- vals

      if (isTRUE(input$gaussian_smooth)) {
        sigma <- input$gaussian_sigma
        if (is.null(sigma) || !is.finite(sigma) || sigma <= 0) sigma <- 1
        mat <- gaussian_smooth_matrix(mat, sigma)
      }

      mat_scaled <- matrix(
        rescale01(as.vector(mat), enhance = isTRUE(input$enhance_contrast)),
        nrow = ny,
        ncol = nx
      )

      pal <- get_msi_palette(input$msi_palette)
      bin <- pmax(1L, pmin(256L, as.integer(mat_scaled * 255) + 1L))
      col_mat <- matrix(pal[bin], nrow = ny, ncol = nx)
      col_mat[is.na(mat_scaled)] <- "#00000000"

      list(
        raster = as.raster(col_mat),
        nx = nx,
        ny = ny,
        mz_selected = mzv[idx],
        mz_index = idx,
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
      if (is.null(label_field) || identical(label_field, "row_index") || !label_field %in% colnames(poly)) {
        return(paste0("polygon_", sprintf("%03d", seq_len(n))))
      }
      vals <- as.character(poly[[label_field]])
      vals[is.na(vals) | trimws(vals) == ""] <- paste0("polygon_", sprintf("%03d", which(is.na(vals) | trimws(vals) == "")))
      make.names(vals, unique = TRUE)
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

        x_src <- cc[, "X"]
        y_src <- cc[, "Y"]
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

        obj <- msi_data()
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

      if (is.null(input$polygon_file) || !nzchar(input$polygon_file)) {
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

      mapped_lab <- rep(NA_character_, nrow(pts_df))
      has_hit <- lengths(hit) > 0
      if (any(has_hit)) {
        if (identical(input$polygon_overlap_rule, "all")) {
          mapped_lab[has_hit] <- vapply(hit[has_hit], function(ix) {
            paste(unique(poly_t$map_label[ix]), collapse = ";")
          }, character(1))
        } else {
          mapped_lab[has_hit] <- vapply(hit[has_hit], function(ix) {
            as.character(poly_t$map_label[ix[1]])
          }, character(1))
        }
      }

      obj <- msi_data()
      pd <- as.data.frame(Cardinal::pData(obj))
      validate(need(nrow(pd) == length(mapped_lab), "Length mismatch between pData and polygon mapping."))
      pd[[col_name]] <- mapped_lab

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

      showNotification(
        sprintf("Mapped polygon annotations to pData column '%s' (%d non-NA assignments).", col_name, sum(!is.na(pd[[col_name]]))),
        type = "message",
        duration = 6
      )
    })

    output$overlay_plot <- renderPlot({
      req(make_msi_raster())
      req(transformed_overlay())

      msi <- make_msi_raster()
      ov <- transformed_overlay()

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
            lwd = 1,
            axes = FALSE,
            reset = FALSE
          )
        }
        graphics::title(main = sprintf("Polygon Overlay on MSI (m/z %.5f)", msi$mz_selected))
      } else {
        cx <- (msi$nx / 2) + input$translate_x
        cy <- (msi$ny / 2) - input$translate_y
        xleft <- cx - ov$width / 2
        xright <- cx + ov$width / 2
        ybottom <- cy - ov$height / 2
        ytop <- cy + ov$height / 2
        graphics::rasterImage(ov$raster, xleft, ybottom, xright, ytop, interpolate = TRUE)
        graphics::title(main = sprintf("Histology Overlay on MSI (m/z %.5f)", msi$mz_selected))
      }
      graphics::box()
    })

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
        mz_selected = round(msi$mz_selected, 6),
        mz_index = msi$mz_index,
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
        enhance_contrast = isTRUE(input$enhance_contrast),
        gaussian_smooth = isTRUE(input$gaussian_smooth),
        gaussian_sigma = input$gaussian_sigma,
        polygon_n = if (!is.null(ov$polygons)) nrow(ov$polygons) else NULL,
        polygon_axis_mode = if (!is.null(ov$axis_mode)) ov$axis_mode else NA_character_,
        polygon_source_dim = if (!is.null(ov$polygon_source_dim)) ov$polygon_source_dim else NA_character_,
        polygon_color = input$polygon_outline_color,
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

      if (is.numeric(v)) {
        ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = value)) +
          ggplot2::geom_raster() +
          ggplot2::scale_fill_gradientn(colors = get_msi_palette(input$msi_palette), na.value = "transparent") +
          ggplot2::coord_equal() +
          ggplot2::labs(fill = input$pdata_field, x = "x", y = "y", title = paste("pData:", input$pdata_field)) +
          ggplot2::theme_minimal()
      } else {
        vf <- as.factor(v)
        lev <- levels(vf)
        pal <- get_discrete_palette(length(lev), input$cluster_palette)
        names(pal) <- lev
        ggplot2::ggplot(data.frame(x = x, y = y_plot, value = vf), ggplot2::aes(x = x, y = y, fill = value)) +
          ggplot2::geom_raster() +
          ggplot2::scale_fill_manual(values = pal, na.value = "transparent") +
          ggplot2::coord_equal() +
          ggplot2::labs(fill = input$pdata_field, x = "x", y = "y", title = paste("pData:", input$pdata_field)) +
          ggplot2::theme_minimal()
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
        return(setNames(as.character(dat$value), as.character(dat$parameter)))
      }

      dat <- try(utils::read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE), silent = TRUE)
      if (!inherits(dat, "try-error") && all(c("parameter", "value") %in% colnames(dat))) {
        return(setNames(as.character(dat$value), as.character(dat$parameter)))
      }

      lines <- try(readLines(path, warn = FALSE), silent = TRUE)
      if (inherits(lines, "try-error")) return(NULL)
      kv <- strsplit(lines, "=", fixed = TRUE)
      ok <- lengths(kv) >= 2
      if (!any(ok)) return(NULL)
      keys <- trimws(vapply(kv[ok], `[`, character(1), 1))
      vals <- trimws(vapply(kv[ok], function(x) paste(x[-1], collapse = "="), character(1)))
      if (length(keys) == 0) return(NULL)
      setNames(vals, keys)
    }

    observeEvent(input$registration_params_upload, {
      req(input$registration_params_upload)
      kv <- read_registration_params(input$registration_params_upload$datapath)
      validate(need(!is.null(kv) && length(kv) > 0, "Could not parse registration parameter file."))

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
      if (!is.null(kv[["intensity_transform"]])) {
        v <- as.character(kv[["intensity_transform"]])
        if (v %in% c("none", "sqrt", "log1p", "asinh")) updateSelectInput(session, "intensity_transform", selected = v)
      }
      if (!is.null(kv[["msi_palette"]])) updateSelectInput(session, "msi_palette", selected = as.character(kv[["msi_palette"]]))
      if (!is.null(kv[["cluster_palette"]])) updateSelectInput(session, "cluster_palette", selected = as.character(kv[["cluster_palette"]]))

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
      if (!grepl("\\.imzML$", filen, ignore.case = TRUE)) filen <- paste0(filen, ".imzML")

      obj <- xh$mapped_obj
      writeImzML(obj, filen)
      saveRDS(obj, paste0(filen, ".rds"))
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
      updateNumericInput(session, "gaussian_sigma", value = 1)
      updateNumericInput(session, "histology_resample_factor", value = 1)
      apply_resolution_scale(show_message = FALSE)
    })
  })
}
