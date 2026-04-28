`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x)) || !nzchar(trimws(as.character(x)[1]))) y else x
}

ra_clamp <- function(x, lo, hi) {
  pmax(lo, pmin(hi, x))
}

ra_wrap_rotation_deg <- function(deg) {
  deg <- suppressWarnings(as.numeric(deg))
  if (!is.finite(deg)) return(NA_real_)
  out <- ((deg + 180) %% 360) - 180
  if (isTRUE(all.equal(out, -180))) out <- 180
  out
}

ra_to_bool <- function(x, default = FALSE) {
  if (is.logical(x) && length(x) > 0 && !is.na(x[1])) return(isTRUE(x[1]))
  sx <- tolower(trimws(as.character(x)[1]))
  if (sx %in% c("true", "t", "1", "yes", "y", "on")) return(TRUE)
  if (sx %in% c("false", "f", "0", "no", "n", "off")) return(FALSE)
  default
}

ra_parse_numeric_tokens <- function(x) {
  if (is.null(x) || length(x) == 0L) return(numeric(0))
  txt <- as.character(x)[1]
  if (is.na(txt) || !nzchar(trimws(txt))) return(numeric(0))
  toks <- trimws(unlist(strsplit(txt, "[,;[:space:]]+", perl = TRUE)))
  toks <- toks[nzchar(toks)]
  if (length(toks) == 0L) return(numeric(0))
  vals <- suppressWarnings(as.numeric(toks))
  vals[is.finite(vals)]
}

ra_safe_color <- function(col, fallback = "#73FFFF") {
  col <- as.character(col)[1]
  if (is.na(col) || !nzchar(trimws(col))) return(fallback)
  ok <- try(grDevices::col2rgb(col), silent = TRUE)
  if (inherits(ok, "try-error")) fallback else col
}

ra_normalize_labels <- function(x) {
  x <- as.character(x)
  x <- iconv(x, from = "", to = "UTF-8", sub = "")
  x <- gsub("[[:cntrl:]\u200B\u200C\u200D\uFEFF]", "", x, perl = TRUE)
  x <- gsub("\\s+", " ", x, perl = TRUE)
  x <- trimws(x)
  x[is.na(x) | !nzchar(x)] <- "NA"
  x
}

ra_are_valid_colors <- function(x) {
  x <- as.character(x)
  ok <- rep(FALSE, length(x))
  keep <- !is.na(x) & nzchar(x)
  if (any(keep)) {
    ok[keep] <- vapply(x[keep], function(val) {
      tryCatch({
        grDevices::col2rgb(val)
        TRUE
      }, error = function(e) FALSE)
    }, logical(1))
  }
  ok
}

ra_get_discrete_palette <- function(n, name = "Alphabet") {
  if (n <= 0L) return(character(0))
  name_chr <- trimws(as.character(name)[1] %||% "Alphabet")
  if (tolower(name_chr) == "alphabet") {
    base_n <- max(26L, n)
    hues <- seq(15, 375, length.out = base_n + 1L)[seq_len(base_n)]
    l_vals <- rep(c(62, 74), length.out = base_n)
    c_vals <- rep(c(85, 70), length.out = base_n)
    pal <- grDevices::hcl(h = hues, c = c_vals, l = l_vals)
  } else {
    pal <- try(grDevices::hcl.colors(max(n, 3L), name_chr), silent = TRUE)
  }
  if (inherits(pal, "try-error") || length(pal) < n) {
    pal <- grDevices::hcl.colors(max(n, 3L), "Set 2")
  }
  pal[seq_len(n)]
}

ra_get_msi_palette <- function(name = "Viridis") {
  pal <- try(grDevices::hcl.colors(256L, name), silent = TRUE)
  if (inherits(pal, "try-error") || length(pal) < 2L) {
    grDevices::hcl.colors(256L, "Inferno")
  } else {
    pal
  }
}

ra_read_registration_params <- function(path) {
  dat <- try(utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE), silent = TRUE)
  if (!inherits(dat, "try-error") && all(c("parameter", "value") %in% colnames(dat))) {
    return(as.list(stats::setNames(as.character(dat$value), as.character(dat$parameter))))
  }
  dat <- try(utils::read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE), silent = TRUE)
  if (!inherits(dat, "try-error") && all(c("parameter", "value") %in% colnames(dat))) {
    return(as.list(stats::setNames(as.character(dat$value), as.character(dat$parameter))))
  }
  lines <- try(readLines(path, warn = FALSE), silent = TRUE)
  if (inherits(lines, "try-error")) return(NULL)
  kv <- strsplit(lines, "=", fixed = TRUE)
  ok <- lengths(kv) >= 2L
  if (!any(ok)) return(NULL)
  keys <- trimws(vapply(kv[ok], `[`, character(1), 1))
  vals <- trimws(vapply(kv[ok], function(x) paste(x[-1], collapse = "="), character(1)))
  as.list(stats::setNames(vals, keys))
}

ra_normalize_loaded_rotate_deg <- function(value, kv = NULL, upload_name = NULL) {
  deg <- suppressWarnings(as.numeric(value))
  deg <- if (length(deg) > 0L) deg[1] else NA_real_
  if (!is.finite(deg)) return(list(value = NA_real_, converted = FALSE, reason = "invalid"))

  kv <- if (is.null(kv)) list() else as.list(kv)
  conv <- ""
  if (!is.null(kv[["registration_rotation_convention"]])) {
    conv <- tolower(trimws(as.character(kv[["registration_rotation_convention"]])[1]))
  }
  if (identical(conv, "shared_overlay_rotation_v2")) {
    return(list(value = ra_wrap_rotation_deg(deg), converted = FALSE, reason = "current"))
  }

  saved_layer <- ""
  if (!is.null(kv[["registration_saved_overlay_layer"]])) {
    saved_layer <- tolower(trimws(as.character(kv[["registration_saved_overlay_layer"]])[1]))
  } else if (!is.null(kv[["overlay_layer"]])) {
    saved_layer <- tolower(trimws(as.character(kv[["overlay_layer"]])[1]))
  }
  upload_label <- if (is.null(upload_name)) "" else tolower(trimws(as.character(upload_name)[1]))
  looks_histology_legacy <- saved_layer %in% c("histology", "cluster") || grepl("histology|hande|h&e", upload_label, perl = TRUE)
  if (isTRUE(looks_histology_legacy)) {
    return(list(value = ra_wrap_rotation_deg(deg - 90), converted = TRUE, reason = "legacy_histology"))
  }
  list(value = ra_wrap_rotation_deg(deg), converted = FALSE, reason = "legacy_default")
}

ra_parse_saved_msi_coord_frame <- function(kv = NULL) {
  kv <- if (is.null(kv)) list() else as.list(kv)
  xmin <- suppressWarnings(as.numeric(kv[["msi_coord_xmin"]]))[1]
  xmax <- suppressWarnings(as.numeric(kv[["msi_coord_xmax"]]))[1]
  ymin <- suppressWarnings(as.numeric(kv[["msi_coord_ymin"]]))[1]
  ymax <- suppressWarnings(as.numeric(kv[["msi_coord_ymax"]]))[1]
  if (!all(is.finite(c(xmin, xmax, ymin, ymax)))) return(NULL)
  nx <- suppressWarnings(as.integer(kv[["msi_canvas_nx"]]))[1]
  ny <- suppressWarnings(as.integer(kv[["msi_canvas_ny"]]))[1]
  if (!is.finite(nx) || nx < 1L) nx <- suppressWarnings(as.integer(round(xmax - xmin + 1)))
  if (!is.finite(ny) || ny < 1L) ny <- suppressWarnings(as.integer(round(ymax - ymin + 1)))
  if (!all(is.finite(c(nx, ny))) || nx < 1L || ny < 1L) return(NULL)
  list(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, nx = nx, ny = ny, source = "saved_params")
}

ra_current_frame_is_saved_subset <- function(current_frame, saved_frame, tol = 1e-8) {
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

ra_resolve_data_file <- function(path, exts) {
  if (file.exists(path) && !dir.exists(path)) return(path)
  if (!dir.exists(path)) stop(sprintf("Path not found: %s", path))
  files <- unlist(lapply(exts, function(ext) Sys.glob(file.path(path, paste0("*", ext)))))
  files <- unique(files[file.exists(files)])
  if (length(files) == 0L) stop(sprintf("No matching files (%s) found in %s", paste(exts, collapse = ", "), path))
  files[1]
}

ra_load_msi_dataset <- function(path) {
  path_use <- path
  if (dir.exists(path)) {
    rds_files <- Sys.glob(file.path(path, "*.rds"))
    if (length(rds_files) > 0L) {
      path_use <- rds_files[1]
    } else {
      imzml_files <- Sys.glob(file.path(path, "*.imzML"))
      if (length(imzml_files) == 0L) stop(sprintf("No .rds or .imzML found in %s", path))
      path_use <- imzml_files[1]
    }
  }
  ext <- tolower(tools::file_ext(path_use))
  if (identical(ext, "rds")) return(readRDS(path_use))
  if (identical(ext, "imzml")) return(Cardinal::readImzML(path_use))
  stop(sprintf("Unsupported MSI dataset type: %s", path_use))
}

ra_read_sf_optional <- function(path) {
  if (is.null(path) || !nzchar(trimws(path))) return(NULL)
  sf::st_read(path, quiet = TRUE, stringsAsFactors = FALSE)
}

ra_current_msi_coord_frame <- function(obj) {
  cd <- as.data.frame(Cardinal::coord(obj))
  if (!is.data.frame(cd) || nrow(cd) == 0L || !all(c("x", "y") %in% names(cd))) return(NULL)
  x <- suppressWarnings(as.numeric(cd$x))
  y <- suppressWarnings(as.numeric(cd$y))
  xmin <- suppressWarnings(min(x, na.rm = TRUE))
  xmax <- suppressWarnings(max(x, na.rm = TRUE))
  ymin <- suppressWarnings(min(y, na.rm = TRUE))
  ymax <- suppressWarnings(max(y, na.rm = TRUE))
  nx <- suppressWarnings(as.integer(round(xmax - xmin + 1)))
  ny <- suppressWarnings(as.integer(round(ymax - ymin + 1)))
  if (!all(is.finite(c(xmin, xmax, ymin, ymax, nx, ny)))) return(NULL)
  list(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, nx = nx, ny = ny)
}

ra_registration_reference_coord_frame <- function(current_frame, saved_frame = NULL) {
  if (ra_current_frame_is_saved_subset(current_frame, saved_frame)) saved_frame else current_frame
}

ra_registration_view_coord_frame <- function(current_frame, saved_frame = NULL) {
  if (!ra_current_frame_is_saved_subset(current_frame, saved_frame)) return(current_frame)
  list(
    xmin = current_frame$xmin,
    xmax = current_frame$xmax,
    ymin = saved_frame$ymax - current_frame$ymax + saved_frame$ymin,
    ymax = saved_frame$ymax - current_frame$ymin + saved_frame$ymin,
    nx = current_frame$nx,
    ny = current_frame$ny
  )
}

ra_make_msi_geometry <- function(obj, saved_frame = NULL) {
  cd <- as.data.frame(Cardinal::coord(obj))
  current_frame <- ra_current_msi_coord_frame(obj)
  reference_frame <- ra_registration_reference_coord_frame(current_frame, saved_frame)
  view_frame <- ra_registration_view_coord_frame(current_frame, saved_frame)
  x_current <- suppressWarnings(as.numeric(cd$x))
  y_current <- suppressWarnings(as.numeric(cd$y))
  x_source <- x_current
  y_source <- y_current
  if (ra_current_frame_is_saved_subset(current_frame, saved_frame)) {
    y_source <- saved_frame$ymax - y_source + saved_frame$ymin
  }
  x_norm <- as.integer(x_current - current_frame$xmin + 1L)
  y_norm <- as.integer(y_current - current_frame$ymin + 1L)
  row_idx <- current_frame$ny - y_norm + 1L
  pd <- try(as.data.frame(Cardinal::pData(obj)), silent = TRUE)
  if (inherits(pd, "try-error")) pd <- data.frame()
  list(
    obj = obj,
    pdata = pd,
    current_frame = current_frame,
    reference_frame = reference_frame,
    view_frame = view_frame,
    x_current = x_current,
    y_current = y_current,
    x_source = x_source,
    y_source = y_source,
    x_norm = x_norm,
    y_norm = y_norm,
    row_idx = row_idx,
    nx = current_frame$nx,
    ny = current_frame$ny
  )
}

ra_apply_subset_coord_mode <- function(msi_geom, saved_frame = NULL, mode = "saved_subset_view") {
  if (is.null(msi_geom) || !is.list(msi_geom) || is.null(msi_geom$current_frame)) return(msi_geom)
  mode_chr <- tolower(trimws(as.character(mode)[1] %||% "saved_subset_view"))
  if (!mode_chr %in% c("saved_subset_view", "preserve_current_frame", "saved_reference_current_view")) {
    mode_chr <- "saved_subset_view"
  }
  if (!identical(mode_chr, "preserve_current_frame")) {
    if (identical(mode_chr, "saved_reference_current_view") &&
        isTRUE(ra_current_frame_is_saved_subset(msi_geom$current_frame, saved_frame))) {
      msi_geom$reference_frame <- saved_frame
      # Match the Shiny app: keep the registration math in the saved full
      # frame, but remap the subset display window into that saved source
      # coordinate space before sampling the histology raster.
      msi_geom$view_frame <- ra_registration_view_coord_frame(msi_geom$current_frame, saved_frame)
      msi_geom$subset_coord_mode <- mode_chr
      return(msi_geom)
    }
    msi_geom$subset_coord_mode <- mode_chr
    return(msi_geom)
  }
  if (!isTRUE(ra_current_frame_is_saved_subset(msi_geom$current_frame, saved_frame))) {
    msi_geom$subset_coord_mode <- mode_chr
    return(msi_geom)
  }
  msi_geom$reference_frame <- msi_geom$current_frame
  msi_geom$view_frame <- msi_geom$current_frame
  msi_geom$x_source <- msi_geom$x_current
  msi_geom$y_source <- msi_geom$y_current
  msi_geom$subset_coord_mode <- mode_chr
  msi_geom
}

ra_transform_intensity <- function(x, trans = "none") {
  trans <- tolower(trimws(as.character(trans)[1] %||% "none"))
  if (identical(trans, "sqrt")) return(sqrt(pmax(x, 0)))
  if (identical(trans, "log1p")) return(log1p(pmax(x, 0)))
  if (identical(trans, "asinh")) return(asinh(x))
  x
}

ra_rescale01 <- function(x, enhance = FALSE) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(rep(0, length(x)))
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

ra_robust_rescale01 <- function(x, q_lo = 0.05, q_hi = 0.95) {
  out <- rep(0, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)
  vals <- as.numeric(x[ok])
  q <- suppressWarnings(stats::quantile(vals, probs = c(q_lo, q_hi), na.rm = TRUE, names = FALSE, type = 8))
  if (!all(is.finite(q)) || q[2] <= q[1]) {
    return(ra_rescale01(x, enhance = FALSE))
  }
  vals <- pmin(pmax(vals, q[1]), q[2])
  out[ok] <- (vals - q[1]) / (q[2] - q[1])
  pmin(1, pmax(0, out))
}

ra_gaussian_smooth_matrix <- function(mat, sigma) {
  sigma <- as.numeric(sigma)
  if (!is.finite(sigma) || sigma <= 0) return(mat)
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
  if (!is.matrix(out)) out <- matrix(out, nrow = nrow(mat), ncol = ncol(mat))
  out
}

ra_build_numeric_signal_matrix <- function(vals, msi_geom, intensity_transform = "none", apply_transform = TRUE, gaussian_smooth = FALSE, smooth_sigma = 1, enhance = FALSE) {
  vals <- as.numeric(vals)
  if (isTRUE(apply_transform)) {
    vals <- ra_transform_intensity(vals, intensity_transform)
  }
  mat <- matrix(NA_real_, nrow = msi_geom$ny, ncol = msi_geom$nx)
  mat[cbind(msi_geom$row_idx, msi_geom$x_norm)] <- vals
  valid_mat <- is.finite(mat)
  if (isTRUE(gaussian_smooth)) {
    mat <- ra_gaussian_smooth_matrix(mat, sigma = smooth_sigma)
  }
  sig <- matrix(ra_rescale01(as.vector(mat), enhance = isTRUE(enhance)), nrow = msi_geom$ny, ncol = msi_geom$nx)
  sig[!is.finite(sig)] <- 0
  list(signal = sig, valid = valid_mat)
}

ra_signal_matrix_to_raster <- function(signal, palette_name = "Viridis") {
  if (is.null(signal) || !is.matrix(signal) || nrow(signal) == 0L || ncol(signal) == 0L) return(NULL)
  mat <- suppressWarnings(matrix(as.numeric(signal), nrow = nrow(signal), ncol = ncol(signal)))
  valid <- is.finite(mat)
  if (!any(valid)) return(as.raster(matrix("#00000000", nrow = nrow(mat), ncol = ncol(mat))))
  mat01 <- matrix(ra_rescale01(as.vector(mat), enhance = FALSE), nrow = nrow(mat), ncol = ncol(mat))
  pal <- ra_get_msi_palette(palette_name)
  bin <- pmax(1L, pmin(256L, as.integer(mat01 * 255) + 1L))
  col_mat <- matrix(pal[bin], nrow = nrow(mat01), ncol = ncol(mat01))
  col_mat[!valid] <- "#00000000"
  as.raster(col_mat)
}

ra_decode_rgba_hex <- function(rgba_hex) {
  rgba_hex <- as.character(rgba_hex)
  n <- length(rgba_hex)
  out <- list(r = rep(0, n), g = rep(0, n), b = rep(0, n), a = rep(0, n))
  keep <- !is.na(rgba_hex) & nchar(rgba_hex) >= 7L
  if (!any(keep)) return(out)
  rr <- suppressWarnings(strtoi(substr(rgba_hex[keep], 2, 3), base = 16L))
  gg <- suppressWarnings(strtoi(substr(rgba_hex[keep], 4, 5), base = 16L))
  bb <- suppressWarnings(strtoi(substr(rgba_hex[keep], 6, 7), base = 16L))
  aa <- rep(255L, sum(keep))
  has_alpha <- nchar(rgba_hex[keep]) >= 9L
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

ra_image_to_raster_rgba <- function(img, alpha_scale = 1) {
  arr <- magick::image_data(img, channels = "rgba")
  h <- dim(arr)[2]
  w <- dim(arr)[3]
  r <- matrix(as.character(arr[1, , ]), nrow = h, ncol = w)
  g <- matrix(as.character(arr[2, , ]), nrow = h, ncol = w)
  b <- matrix(as.character(arr[3, , ]), nrow = h, ncol = w)
  a <- matrix(as.character(arr[4, , ]), nrow = h, ncol = w)
  a_byte <- suppressWarnings(strtoi(as.vector(a), base = 16L))
  a_byte[!is.finite(a_byte)] <- 255L
  alpha_scale <- suppressWarnings(as.numeric(alpha_scale))
  if (!is.finite(alpha_scale) || alpha_scale < 0) alpha_scale <- 1
  a_byte <- pmax(0L, pmin(255L, as.integer(round(a_byte * alpha_scale))))
  aa <- matrix(sprintf("%02X", a_byte), nrow = h, ncol = w)
  as.raster(matrix(paste0("#", as.vector(r), as.vector(g), as.vector(b), as.vector(aa)), nrow = h, ncol = w))
}

ra_overlay_uses_qupath_source_frame <- function(overlay_frame) {
  isTRUE(overlay_frame$roi_anchor_present) ||
    identical(overlay_frame$source_frame_type, "whole_slide") ||
    identical(overlay_frame$source_frame_type, "roi_relative") ||
    identical(overlay_frame$source_frame_type, "image_bbox")
}

ra_effective_overlay_rotate_deg <- function(rotate_deg) {
  deg <- suppressWarnings(as.numeric(rotate_deg))[1]
  if (!is.finite(deg)) deg <- 0
  deg + 90
}

ra_effective_polygon_flip_y <- function(overlay_frame, user_flip = FALSE) {
  xor(isTRUE(ra_overlay_uses_qupath_source_frame(overlay_frame)), isTRUE(user_flip))
}

ra_transform_overlay_image <- function(img, msi_geom, overlay_frame, transform, alpha_scale = 1) {
  img_work <- img
  if (isTRUE(transform$flip_histology_y)) img_work <- magick::image_flip(img_work)
  info0 <- magick::image_info(img_work)[1, ]
  scale_mode <- tolower(trimws(as.character(transform$overlay_scale_mode)[1] %||% "absolute"))
  if (!scale_mode %in% c("absolute", "fit")) scale_mode <- "absolute"
  if (identical(scale_mode, "fit")) {
    fit_scale <- min(msi_geom$nx / info0$width, msi_geom$ny / info0$height)
    target_w <- max(1L, as.integer(round(info0$width * fit_scale * transform$scale_x)))
    target_h <- max(1L, as.integer(round(info0$height * fit_scale * transform$scale_y)))
  } else {
    target_w <- max(1L, as.integer(round(info0$width * transform$scale_x)))
    target_h <- max(1L, as.integer(round(info0$height * transform$scale_y)))
  }
  img_work <- magick::image_resize(img_work, paste0(target_w, "x", target_h, "!"))
  img_work <- magick::image_background(img_work, "none")
  rot_deg <- ra_effective_overlay_rotate_deg(transform$rotate_deg)
  if (is.finite(rot_deg) && rot_deg != 0) {
    img_work <- magick::image_rotate(img_work, rot_deg)
  }
  ras <- ra_image_to_raster_rgba(img_work, alpha_scale = alpha_scale)
  coord_frame <- msi_geom$reference_frame
  src_center_x <- ((coord_frame$xmin + coord_frame$xmax) / 2) + transform$translate_x
  src_center_y <- ((coord_frame$ymin + coord_frame$ymax) / 2) + transform$translate_y
  list(
    raster = ras,
    width = ncol(ras),
    height = nrow(ras),
    translate_x_base = transform$translate_x,
    translate_y_base = transform$translate_y,
    xleft_source = src_center_x - (ncol(ras) / 2),
    xright_source = src_center_x + (ncol(ras) / 2),
    ybottom_source = src_center_y - (nrow(ras) / 2),
    ytop_source = src_center_y + (nrow(ras) / 2),
    source_width = as.numeric(info0$width),
    source_height = as.numeric(info0$height)
  )
}

ra_sample_transformed_overlay_to_msi <- function(msi_geom, tr, tx, ty) {
  w <- as.integer(tr$width)
  h <- as.integer(tr$height)
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
  u <- (msi_geom$x_source - xleft) / (xright - xleft)
  v <- (ytop - msi_geom$y_source) / (ytop - ybottom)
  inside <- is.finite(u) & is.finite(v) & u >= 0 & u <= 1 & v >= 0 & v <= 1
  col_idx <- pmin(w, pmax(1L, as.integer(floor(u * (w - 1L)) + 1L)))
  row_idx <- pmin(h, pmax(1L, as.integer(floor(v * (h - 1L)) + 1L)))
  col_mat <- as.matrix(tr$raster)
  sampled <- rep(NA_character_, length(msi_geom$x_source))
  sampled[inside] <- col_mat[cbind(row_idx[inside], col_idx[inside])]
  list(sampled_rgba = sampled, inside = inside, width = w, height = h)
}

ra_overlay_transform_source_bbox <- function(tr, tx = NULL, ty = NULL) {
  if (is.null(tr)) return(NULL)
  tx_base <- suppressWarnings(as.numeric(tr$translate_x_base))[1]
  ty_base <- suppressWarnings(as.numeric(tr$translate_y_base))[1]
  tx_now <- suppressWarnings(as.numeric(tx))[1]
  ty_now <- suppressWarnings(as.numeric(ty))[1]
  if (!is.finite(tx_base)) tx_base <- 0
  if (!is.finite(ty_base)) ty_base <- 0
  if (!is.finite(tx_now)) tx_now <- tx_base
  if (!is.finite(ty_now)) ty_now <- ty_base
  dx_shift <- tx_now - tx_base
  dy_shift <- ty_now - ty_base
  c(
    xmin = tr$xleft_source + dx_shift,
    xmax = tr$xright_source + dx_shift,
    ymin = tr$ybottom_source + dy_shift,
    ymax = tr$ytop_source + dy_shift
  )
}

ra_transform_polygon_sf <- function(poly_sf, nx, ny, scale_x, scale_y, translate_x, translate_y,
                                    rotate_deg, flip_y = FALSE, swap_xy = FALSE,
                                    scale_mode = "absolute",
                                    source_width = NA_real_, source_height = NA_real_,
                                    source_origin_x = NA_real_, source_origin_y = NA_real_,
                                    downsample_factor = 1,
                                    coord_frame = NULL) {
  geom_type <- as.character(sf::st_geometry_type(poly_sf))
  keep <- geom_type %in% c("POLYGON", "MULTIPOLYGON")
  poly <- poly_sf[keep, , drop = FALSE]
  if (nrow(poly) == 0L) stop("No polygon geometries found in input file.")

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
  if (!all(is.finite(c(bb_w, bb_h))) || bb_w <= 0 || bb_h <= 0) {
    stop("Polygon bounding box is invalid.")
  }

  has_source_dims <- is.finite(source_width) && is.finite(source_height) &&
    source_width > 0 && source_height > 0
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
  if (!all(is.finite(c(frame_w, frame_h))) || frame_w <= 0 || frame_h <= 0) {
    stop("Polygon reference frame is invalid.")
  }

  downsample_factor <- suppressWarnings(as.numeric(downsample_factor)[1])
  if (!is.finite(downsample_factor) || downsample_factor <= 0) downsample_factor <- 1
  frame_w_img <- frame_w / downsample_factor
  frame_h_img <- frame_h / downsample_factor

  mode <- tolower(trimws(as.character(scale_mode)[1]))
  if (!mode %in% c("absolute", "fit")) mode <- "absolute"
  if (identical(mode, "fit")) {
    fit_scale <- min(nx / frame_w_img, ny / frame_h_img)
    sx <- fit_scale * suppressWarnings(as.numeric(scale_x))[1]
    sy <- fit_scale * suppressWarnings(as.numeric(scale_y))[1]
  } else {
    sx <- suppressWarnings(as.numeric(scale_x))[1]
    sy <- suppressWarnings(as.numeric(scale_y))[1]
  }
  if (!all(is.finite(c(sx, sy))) || sx <= 0 || sy <= 0) {
    stop("Polygon scale must be positive.")
  }

  center_frame <- coord_frame
  if (is.null(center_frame)) {
    center_frame <- list(xmin = 1, xmax = nx, ymin = 1, ymax = ny, nx = nx, ny = ny)
  }
  cx <- ((center_frame$xmin + center_frame$xmax) / 2) + suppressWarnings(as.numeric(translate_x))[1]
  cy <- ((center_frame$ymin + center_frame$ymax) / 2) + suppressWarnings(as.numeric(translate_y))[1]
  theta <- suppressWarnings(as.numeric(rotate_deg))[1] * pi / 180
  if (!is.finite(theta)) theta <- 0
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
    if (nrow(cc) == 0L) stop("Polygon has empty coordinates after read.")
    cc_names <- colnames(cc)
    x_col <- if ("X" %in% cc_names) "X" else if ("x" %in% cc_names) "x" else NA_character_
    y_col <- if ("Y" %in% cc_names) "Y" else if ("y" %in% cc_names) "y" else NA_character_
    if (!is.character(x_col) || !nzchar(x_col) || !is.character(y_col) || !nzchar(y_col)) {
      xy_cols <- cc_names[cc_names %in% c("X", "Y", "x", "y")]
      if (length(xy_cols) >= 2L) {
        x_col <- xy_cols[1]
        y_col <- xy_cols[2]
      }
    }
    if (!is.character(x_col) || !nzchar(x_col) || !x_col %in% cc_names ||
        !is.character(y_col) || !nzchar(y_col) || !y_col %in% cc_names) {
      stop("Polygon coordinates missing X/Y columns after transformation.")
    }

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
          if (nrow(ring) > 2L && !all(ring[1, ] == ring[nrow(ring), ])) ring <- rbind(ring, ring[1, ])
          ring
        })
        sf::st_polygon(rings)
      })
      if (length(polys) == 1L) polys[[1]] else sf::st_multipolygon(polys)
    } else if ("L2" %in% colnames(cc_new)) {
      rings <- lapply(unique(cc_new[, "L2"]), function(l2) {
        ring <- as.matrix(cc_new[cc_new[, "L2"] == l2, c("x1", "y1"), drop = FALSE])
        if (nrow(ring) > 2L && !all(ring[1, ] == ring[nrow(ring), ])) ring <- rbind(ring, ring[1, ])
        ring
      })
      sf::st_polygon(rings)
    } else {
      ring <- as.matrix(cc_new[, c("x1", "y1"), drop = FALSE])
      if (nrow(ring) > 2L && !all(ring[1, ] == ring[nrow(ring), ])) ring <- rbind(ring, ring[1, ])
      sf::st_polygon(list(ring))
    }
  }

  geoms <- lapply(seq_len(nrow(poly)), function(i) try(rebuild_one(poly[i, , drop = FALSE]), silent = TRUE))
  good <- !vapply(geoms, inherits, logical(1), "try-error")
  poly <- poly[good, , drop = FALSE]
  geoms <- geoms[good]
  if (length(geoms) == 0L) stop("Could not transform polygon geometries.")

  sfc <- sf::st_sfc(geoms)
  poly <- sf::st_set_geometry(poly, sfc)
  poly <- sf::st_make_valid(poly)
  gtypes <- unique(as.character(sf::st_geometry_type(poly, by_geometry = TRUE)))
  if (any(!gtypes %in% c("POLYGON", "MULTIPOLYGON"))) {
    poly <- sf::st_collection_extract(poly, "POLYGON", warn = FALSE)
  }
  if (nrow(poly) == 0L) stop("No valid polygon geometries remain after transformation.")
  poly
}

ra_extract_msi_intensity_matrix <- function(obj) {
  sp <- try(as.matrix(Cardinal::spectra(obj)), silent = TRUE)
  if (!(inherits(sp, "try-error") || is.null(sp))) {
    if (ncol(sp) != ncol(obj) && nrow(sp) == ncol(obj)) sp <- t(sp)
    storage.mode(sp) <- "double"
    return(sp)
  }

  spectra_data <- try(methods::slot(obj, "spectraData"), silent = TRUE)
  if (!(inherits(spectra_data, "try-error") || is.null(spectra_data))) {
    intensity <- try(spectra_data[["intensity"]], silent = TRUE)
    if (inherits(intensity, "try-error") || is.null(intensity)) {
      intensity <- try(spectra_data$intensity, silent = TRUE)
    }
    if (!(inherits(intensity, "try-error") || is.null(intensity))) {
      sp <- try(as.matrix(intensity), silent = TRUE)
      if (!(inherits(sp, "try-error") || is.null(sp))) {
        if (ncol(sp) != ncol(obj) && nrow(sp) == ncol(obj)) sp <- t(sp)
        storage.mode(sp) <- "double"
        return(sp)
      }
    }
  }

  ra_extract_ibd_path_from_object <- function(obj) {
    spectra_data <- try(methods::slot(obj, "spectraData"), silent = TRUE)
    if (inherits(spectra_data, "try-error") || is.null(spectra_data)) return(NULL)
    intensity <- try(spectra_data[["intensity"]], silent = TRUE)
    if (inherits(intensity, "try-error") || is.null(intensity)) {
      intensity <- try(spectra_data[[1]], silent = TRUE)
    }
    if (inherits(intensity, "try-error") || is.null(intensity) || !isS4(intensity)) return(NULL)
    data_slot <- try(methods::slot(intensity, "data"), silent = TRUE)
    if (inherits(data_slot, "try-error") || is.null(data_slot) || !isS4(data_slot)) return(NULL)
    atoms <- try(methods::slot(data_slot, "data"), silent = TRUE)
    if (inherits(atoms, "try-error") || is.null(atoms) || !isS4(atoms)) return(NULL)
    src <- try(methods::slot(atoms, "source"), silent = TRUE)
    if (inherits(src, "try-error") || is.null(src)) return(NULL)
    lev <- try(levels(src), silent = TRUE)
    if (inherits(lev, "try-error") || length(lev) < 1L) return(NULL)
    ibd_path <- as.character(lev[1])
    if (!nzchar(ibd_path) || !file.exists(ibd_path)) return(NULL)
    ibd_path
  }

  ra_extract_msi_intensity_matrix_direct_imzml <- function(obj) {
    if (!requireNamespace("CardinalIO", quietly = TRUE)) return(NULL)
    ibd_path <- ra_extract_ibd_path_from_object(obj)
    if (is.null(ibd_path)) return(NULL)

    imzml_path <- sub("\\.[iI][bB][dD]$", ".imzML", ibd_path, perl = TRUE)
    if (!file.exists(imzml_path)) {
      imzml_path <- sub("\\.[iI][bB][dD]$", ".imzml", ibd_path, perl = TRUE)
    }
    if (!file.exists(imzml_path)) return(NULL)

    parsed <- try(CardinalIO::parseImzML(imzml_path), silent = TRUE)
    if (inherits(parsed, "try-error") || is.null(parsed)) return(NULL)
    spectrum_list <- try(methods::slot(parsed, "listData")$run$spectrumList, silent = TRUE)
    if (inherits(spectrum_list, "try-error") || is.null(spectrum_list)) return(NULL)
    intensity_arrays <- spectrum_list$intensityArrays
    if (!is.data.frame(intensity_arrays) || nrow(intensity_arrays) < 1L) return(NULL)

    offsets <- suppressWarnings(as.numeric(intensity_arrays[["external offset"]]))
    array_lengths <- suppressWarnings(as.integer(intensity_arrays[["external array length"]]))
    encoded_lengths <- suppressWarnings(as.integer(intensity_arrays[["external encoded length"]]))
    dtype <- tolower(trimws(as.character(intensity_arrays[["binary data type"]])))
    compression <- tolower(trimws(as.character(intensity_arrays[["binary data compression type"]])))

    if (!all(is.finite(offsets)) || !all(is.finite(array_lengths)) || !all(is.finite(encoded_lengths))) return(NULL)
    if (length(unique(compression)) != 1L || !compression[1] %in% c("no compression", "none", "")) return(NULL)
    if (length(unique(dtype)) != 1L) return(NULL)

    size_bytes <- switch(
      dtype[1],
      "32-bit float" = 4L,
      "64-bit float" = 8L,
      NULL
    )
    if (is.null(size_bytes)) return(NULL)

    n_features <- unique(array_lengths)
    block_bytes <- unique(encoded_lengths)
    if (length(n_features) != 1L || length(block_bytes) != 1L) return(NULL)
    n_features <- as.integer(n_features[1])
    block_bytes <- as.integer(block_bytes[1])
    if (!is.finite(n_features) || n_features < 1L || !is.finite(block_bytes) || block_bytes < 1L) return(NULL)
    if (!identical(block_bytes, as.integer(n_features * size_bytes))) return(NULL)

    n_pixels <- nrow(intensity_arrays)
    offsets_expected <- offsets[1] + block_bytes * seq.int(0, n_pixels - 1L)
    contiguous <- all(abs(offsets - offsets_expected) < 1e-8)

    con <- file(ibd_path, open = "rb")
    on.exit(try(close(con), silent = TRUE), add = TRUE)

    if (isTRUE(contiguous)) {
      seek(con, where = offsets[1], rw = "read")
      vals <- try(readBin(con, what = "numeric", n = n_features * n_pixels, size = size_bytes, endian = "little"), silent = TRUE)
      if (!(inherits(vals, "try-error") || length(vals) != (n_features * n_pixels))) {
        return(matrix(vals, nrow = n_features, ncol = n_pixels))
      }
    }

    mat <- matrix(NA_real_, nrow = n_features, ncol = n_pixels)
    for (jj in seq_len(n_pixels)) {
      seek(con, where = offsets[jj], rw = "read")
      vals <- try(readBin(con, what = "numeric", n = n_features, size = size_bytes, endian = "little"), silent = TRUE)
      if (inherits(vals, "try-error") || length(vals) != n_features) return(NULL)
      mat[, jj] <- vals
    }
    mat
  }

  sp <- try(ra_extract_msi_intensity_matrix_direct_imzml(obj), silent = TRUE)
  if (!(inherits(sp, "try-error") || is.null(sp))) {
    if (ncol(sp) != ncol(obj) && nrow(sp) == ncol(obj)) sp <- t(sp)
    storage.mode(sp) <- "double"
    return(sp)
  }

  stop("Could not extract MSI intensities.")
}

ra_build_histology_feature_from_sampled <- function(sampled_rgba, msi_geom, feature_mode = "hematoxylin", gaussian_smooth = FALSE, smooth_sigma = 1) {
  feat_mode <- tolower(trimws(as.character(feature_mode)[1]))
  if (!feat_mode %in% c("hematoxylin", "darkness", "purple")) feat_mode <- "hematoxylin"
  ch <- ra_decode_rgba_hex(sampled_rgba)
  valid_vec <- !is.na(sampled_rgba) & is.finite(ch$a) & ch$a > 0
  if (!any(valid_vec)) {
    return(list(signal = matrix(0, nrow = msi_geom$ny, ncol = msi_geom$nx), valid = matrix(FALSE, nrow = msi_geom$ny, ncol = msi_geom$nx)))
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
  mat <- matrix(0, nrow = msi_geom$ny, ncol = msi_geom$nx)
  mat[cbind(msi_geom$row_idx, msi_geom$x_norm)] <- feat
  valid_mat <- matrix(FALSE, nrow = msi_geom$ny, ncol = msi_geom$nx)
  valid_mat[cbind(msi_geom$row_idx, msi_geom$x_norm)] <- valid_vec
  if (isTRUE(gaussian_smooth)) mat <- ra_gaussian_smooth_matrix(mat, smooth_sigma)
  mat <- matrix(ra_rescale01(as.vector(mat), enhance = FALSE), nrow = msi_geom$ny, ncol = msi_geom$nx)
  mat[!is.finite(mat)] <- 0
  list(signal = mat, valid = valid_mat)
}

ra_build_binary_mask_signal_from_sampled <- function(sampled_rgba, msi_geom, gaussian_smooth = FALSE, smooth_sigma = 1) {
  ch <- ra_decode_rgba_hex(sampled_rgba)
  val <- as.numeric(ch$a)
  val[!is.finite(val)] <- 0
  mat <- matrix(0, nrow = msi_geom$ny, ncol = msi_geom$nx)
  mat[cbind(msi_geom$row_idx, msi_geom$x_norm)] <- val
  valid_mat <- matrix(FALSE, nrow = msi_geom$ny, ncol = msi_geom$nx)
  valid_mat[cbind(msi_geom$row_idx, msi_geom$x_norm)] <- !is.na(sampled_rgba)
  if (isTRUE(gaussian_smooth)) mat <- ra_gaussian_smooth_matrix(mat, smooth_sigma)
  mat <- matrix(ra_rescale01(as.vector(mat), enhance = FALSE), nrow = msi_geom$ny, ncol = msi_geom$nx)
  mat[!is.finite(mat)] <- 0
  list(signal = mat, valid = valid_mat)
}

ra_build_edge_mask <- function(signal, q = 0.85, valid_mask = NULL) {
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
  right <- cbind(sig[, -1, drop = FALSE], sig[, nx, drop = FALSE])
  down <- rbind(sig[-1, , drop = FALSE], sig[ny, , drop = FALSE])
  right_valid <- cbind(valid[, -1, drop = FALSE], FALSE)
  down_valid <- rbind(valid[-1, , drop = FALSE], rep(FALSE, nx))
  grad <- abs(right - sig) + abs(down - sig)
  grad[!(valid & right_valid & down_valid)] <- 0
  g_ok <- grad[valid & is.finite(grad) & grad > 0]
  if (length(g_ok) < 10L) return(NULL)
  thr <- suppressWarnings(stats::quantile(g_ok, probs = q, na.rm = TRUE, names = FALSE, type = 8))
  if (!is.finite(thr)) thr <- suppressWarnings(mean(g_ok, na.rm = TRUE))
  edge <- grad >= thr
  edge[!valid] <- FALSE
  n_edge <- suppressWarnings(sum(edge, na.rm = TRUE))
  if (!is.finite(n_edge) || n_edge < 5L) return(NULL)
  list(edge = edge, n_edge = as.integer(n_edge), threshold = thr)
}

ra_normalized_mutual_information <- function(x, y, nbins = 24L) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 20L) return(NA_real_)
  xk <- x[keep]
  yk <- y[keep]
  xr <- range(xk, na.rm = TRUE)
  yr <- range(yk, na.rm = TRUE)
  if (!all(is.finite(c(xr, yr))) || diff(xr) <= 1e-8 || diff(yr) <= 1e-8) return(NA_real_)
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

ra_signed_rank_correlation <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 20L) return(NA_real_)
  xk <- x[keep]
  yk <- y[keep]
  xr <- range(xk, na.rm = TRUE)
  yr <- range(yk, na.rm = TRUE)
  if (!all(is.finite(c(xr, yr))) || diff(xr) <= 1e-8 || diff(yr) <= 1e-8) return(NA_real_)
  sc <- suppressWarnings(stats::cor(xk, yk, method = "spearman"))
  if (!is.finite(sc)) return(NA_real_)
  max(-1, min(1, as.numeric(sc)))
}

ra_categorical_signal_association <- function(x, labels, min_group_n = 5L) {
  x <- as.numeric(x)
  labs <- trimws(as.character(labels))
  keep <- is.finite(x) & !is.na(labs) & nzchar(labs)
  if (sum(keep) < 25L) return(NA_real_)
  xk <- x[keep]
  lk <- labs[keep]
  group_n <- table(lk)
  keep_groups <- names(group_n[group_n >= min_group_n])
  if (length(keep_groups) < 2L) return(NA_real_)
  group_keep <- lk %in% keep_groups
  xk <- xk[group_keep]
  lk <- lk[group_keep]
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
  max(0, min(1, as.numeric(eta2)))
}

ra_binary_mask_label_association <- function(signal_mask, labels, min_group_n = 20L) {
  mask <- as.logical(signal_mask)
  labs <- trimws(as.character(labels))
  keep <- !is.na(mask) & !is.na(labs) & nzchar(labs)
  if (sum(keep) < 50L) {
    return(list(score = NA_real_, nmi = NA_real_, best_dice = NA_real_, best_label = NA_character_))
  }
  mask <- mask[keep]
  labs <- labs[keep]
  group_n <- table(labs)
  keep_groups <- names(group_n[group_n >= min_group_n])
  if (length(keep_groups) < 2L) {
    return(list(score = NA_real_, nmi = NA_real_, best_dice = NA_real_, best_label = NA_character_))
  }
  keep2 <- labs %in% keep_groups
  mask <- mask[keep2]
  labs <- labs[keep2]
  if (length(unique(mask)) < 2L) {
    return(list(score = NA_real_, nmi = NA_real_, best_dice = NA_real_, best_label = NA_character_))
  }

  tab <- table(mask = mask, label = labs)
  pxy <- prop.table(tab)
  px <- rowSums(pxy)
  py <- colSums(pxy)
  nz <- pxy > 0
  mi <- sum(pxy[nz] * log(pxy[nz] / (px[row(pxy)][nz] * py[col(pxy)][nz])))
  hx <- -sum(px[px > 0] * log(px[px > 0]))
  hy <- -sum(py[py > 0] * log(py[py > 0]))
  nmi <- if (is.finite(hx) && is.finite(hy) && (hx + hy) > 0) (2 * mi) / (hx + hy) else NA_real_

  best_dice <- -Inf
  best_label <- NA_character_
  for (lab in keep_groups) {
    lab_mask <- labs == lab
    inter <- sum(mask & lab_mask, na.rm = TRUE)
    denom <- sum(mask, na.rm = TRUE) + sum(lab_mask, na.rm = TRUE)
    if (denom <= 0) next
    dice <- (2 * inter) / denom
    if (is.finite(dice) && dice > best_dice) {
      best_dice <- dice
      best_label <- lab
    }
  }
  if (!is.finite(best_dice)) best_dice <- NA_real_
  score <- ra_combine_weighted_scores(c(nmi, best_dice), c(0.65, 0.35))
  list(score = score, nmi = nmi, best_dice = best_dice, best_label = best_label)
}

ra_score_binary_mask_against_categorical_target <- function(signal_mask, labels, min_group_n = 20L) {
  mask <- matrix(as.logical(signal_mask), nrow = nrow(signal_mask), ncol = ncol(signal_mask))
  labs <- matrix(trimws(as.character(labels)), nrow = nrow(labels), ncol = ncol(labels))
  if (!all(dim(mask) == dim(labs))) {
    return(list(
      score = NA_real_, best_label = NA_character_, best_label_score = NA_real_,
      best_dice = NA_real_, best_jaccard = NA_real_, best_boundary = NA_real_,
      n_labels = 0L
    ))
  }

  keep <- !is.na(labs) & nzchar(labs)
  if (sum(keep, na.rm = TRUE) < 50L || sum(mask, na.rm = TRUE) < 5L) {
    return(list(
      score = NA_real_, best_label = NA_character_, best_label_score = NA_real_,
      best_dice = NA_real_, best_jaccard = NA_real_, best_boundary = NA_real_,
      n_labels = 0L
    ))
  }

  label_counts <- sort(table(labs[keep]), decreasing = TRUE)
  keep_labels <- names(label_counts[label_counts >= min_group_n])
  if (length(keep_labels) < 1L) {
    return(list(
      score = NA_real_, best_label = NA_character_, best_label_score = NA_real_,
      best_dice = NA_real_, best_jaccard = NA_real_, best_boundary = NA_real_,
      n_labels = 0L
    ))
  }

  per_label <- lapply(keep_labels, function(lab) {
    lab_mask <- labs == lab
    bin <- ra_score_binary_masks(mask, lab_mask)
    inside_frac <- suppressWarnings(mean(lab_mask[mask], na.rm = TRUE))
    outside_frac <- suppressWarnings(mean(lab_mask[!mask], na.rm = TRUE))
    enrichment <- inside_frac - outside_frac
    score <- ra_combine_weighted_scores(
      c(bin$score, enrichment, bin$boundary),
      c(0.65, 0.20, 0.15)
    )
    list(
      label = lab,
      score = score,
      dice = bin$dice,
      jaccard = bin$jaccard,
      boundary = bin$boundary,
      enrichment = enrichment
    )
  })

  score_vec <- vapply(per_label, function(x) x$score, numeric(1))
  score_vec[!is.finite(score_vec)] <- -Inf
  best_idx <- which.max(score_vec)
  if (length(best_idx) == 0L || !is.finite(score_vec[best_idx])) {
    return(list(
      score = NA_real_, best_label = NA_character_, best_label_score = NA_real_,
      best_dice = NA_real_, best_jaccard = NA_real_, best_boundary = NA_real_,
      n_labels = length(keep_labels)
    ))
  }
  best <- per_label[[best_idx]]
  list(
    score = best$score,
    best_label = best$label,
    best_label_score = best$score,
    best_dice = best$dice,
    best_jaccard = best$jaccard,
    best_boundary = best$boundary,
    best_enrichment = best$enrichment,
    n_labels = length(keep_labels)
  )
}

ra_shift_bool <- function(mat, dr = 0L, dc = 0L) {
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

ra_dilate8 <- function(mat) {
  mat <- mat & TRUE
  out <- mat
  for (dr in -1L:1L) {
    for (dc in -1L:1L) {
      if (dr == 0L && dc == 0L) next
      out <- out | ra_shift_bool(mat, dr = dr, dc = dc)
    }
  }
  out
}

ra_erode8 <- function(mat) {
  mat <- mat & TRUE
  out <- mat
  for (dr in -1L:1L) {
    for (dc in -1L:1L) {
      if (dr == 0L && dc == 0L) next
      out <- out & ra_shift_bool(mat, dr = dr, dc = dc)
    }
  }
  out
}

ra_idx_to_coords <- function(idx, nr) {
  idx <- as.integer(idx)
  if (length(idx) == 0L) return(list(r = integer(0), c = integer(0), n_total = 0L))
  list(r = ((idx - 1L) %% nr) + 1L, c = ((idx - 1L) %/% nr) + 1L, n_total = length(idx))
}

ra_build_mask_bands <- function(mask, band_px = 4L) {
  mask <- (mask & TRUE)
  ny <- nrow(mask)
  band_px <- as.integer(min(50L, max(1L, band_px)))
  inside_rings <- vector("list", band_px)
  outside_rings <- vector("list", band_px)
  inside_work <- mask
  grown <- mask
  for (k in seq_len(band_px)) {
    eroded <- ra_erode8(inside_work)
    in_ring <- inside_work & !eroded
    inside_work <- eroded
    dilated <- ra_dilate8(grown)
    out_ring <- dilated & !grown
    grown <- dilated
    inside_rings[[k]] <- ra_idx_to_coords(which(in_ring), ny)
    outside_rings[[k]] <- ra_idx_to_coords(which(out_ring), ny)
  }
  list(
    ny = ny,
    nx = ncol(mask),
    band_px = band_px,
    inside_all = ra_idx_to_coords(which(mask), ny),
    inside_core = ra_idx_to_coords(which(inside_work), ny),
    boundary = ra_idx_to_coords(which(mask & !ra_erode8(mask)), ny),
    inside_rings = inside_rings,
    outside_rings = outside_rings
  )
}

ra_coord_obj_to_mask <- function(coord_obj, nr, nc) {
  out <- matrix(FALSE, nrow = nr, ncol = nc)
  if (is.null(coord_obj) || length(coord_obj$r) == 0L) return(out)
  keep <- coord_obj$r >= 1L & coord_obj$r <= nr & coord_obj$c >= 1L & coord_obj$c <= nc
  if (!any(keep)) return(out)
  out[cbind(coord_obj$r[keep], coord_obj$c[keep])] <- TRUE
  out
}

ra_connected_components <- function(mask) {
  mask <- matrix(as.logical(mask), nrow = nrow(mask), ncol = ncol(mask))
  nr <- nrow(mask)
  nc <- ncol(mask)
  labels <- matrix(0L, nrow = nr, ncol = nc)
  idx <- which(mask)
  if (length(idx) == 0L) {
    return(list(labels = labels, sizes = integer(0)))
  }
  next_id <- 0L
  sizes <- integer(0)
  neigh <- expand.grid(dr = -1L:1L, dc = -1L:1L, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  neigh <- neigh[!(neigh$dr == 0L & neigh$dc == 0L), , drop = FALSE]
  for (start in idx) {
    r0 <- ((start - 1L) %% nr) + 1L
    c0 <- ((start - 1L) %/% nr) + 1L
    if (!mask[r0, c0] || labels[r0, c0] != 0L) next
    next_id <- next_id + 1L
    queue_r <- integer(length(idx))
    queue_c <- integer(length(idx))
    head <- 1L
    tail <- 1L
    queue_r[tail] <- r0
    queue_c[tail] <- c0
    labels[r0, c0] <- next_id
    comp_n <- 0L
    while (head <= tail) {
      rr <- queue_r[head]
      cc <- queue_c[head]
      head <- head + 1L
      comp_n <- comp_n + 1L
      for (ii in seq_len(nrow(neigh))) {
        rn <- rr + neigh$dr[ii]
        cn <- cc + neigh$dc[ii]
        if (rn < 1L || rn > nr || cn < 1L || cn > nc) next
        if (!mask[rn, cn] || labels[rn, cn] != 0L) next
        tail <- tail + 1L
        queue_r[tail] <- rn
        queue_c[tail] <- cn
        labels[rn, cn] <- next_id
      }
    }
    sizes[next_id] <- comp_n
  }
  list(labels = labels, sizes = as.integer(sizes))
}

ra_mask_shape_metrics <- function(mask, valid_mask = NULL) {
  mask <- matrix(as.logical(mask), nrow = nrow(mask), ncol = ncol(mask))
  if (is.null(valid_mask)) valid_mask <- matrix(TRUE, nrow = nrow(mask), ncol = ncol(mask))
  valid_mask <- matrix(as.logical(valid_mask), nrow = nrow(mask), ncol = ncol(mask))
  total_valid <- suppressWarnings(sum(valid_mask, na.rm = TRUE))
  mask_valid <- mask & valid_mask
  n_mask <- suppressWarnings(sum(mask_valid, na.rm = TRUE))
  if (!is.finite(n_mask) || n_mask <= 0L || !is.finite(total_valid) || total_valid <= 0L) {
    return(list(
      support_fraction = 0,
      largest_component_fraction = 0,
      boundary_fraction = 0,
      n_components = 0L,
      n_mask = 0L
    ))
  }
  comps <- ra_connected_components(mask_valid)
  largest_component_fraction <- if (length(comps$sizes) > 0L) max(comps$sizes) / n_mask else 0
  boundary_mask <- mask_valid & !ra_erode8(mask_valid)
  boundary_fraction <- suppressWarnings(sum(boundary_mask, na.rm = TRUE)) / n_mask
  list(
    support_fraction = n_mask / total_valid,
    largest_component_fraction = largest_component_fraction,
    boundary_fraction = boundary_fraction,
    n_components = length(comps$sizes),
    n_mask = as.integer(n_mask)
  )
}

ra_sample_shifted_stats <- function(coord_obj, signal, dx = 0L, dy = 0L, min_n = 5L) {
  if (is.null(coord_obj) || length(coord_obj$r) == 0L) return(list(ok = FALSE, mean = NA_real_, q90 = NA_real_, n = 0L, cov = 0))
  ny <- nrow(signal)
  nx <- ncol(signal)
  rr <- coord_obj$r + as.integer(dy)
  cc <- coord_obj$c + as.integer(dx)
  keep <- rr >= 1L & rr <= ny & cc >= 1L & cc <= nx
  n_keep <- sum(keep)
  cov <- n_keep / max(1L, as.integer(coord_obj$n_total))
  if (!is.finite(n_keep) || n_keep < min_n) return(list(ok = FALSE, mean = NA_real_, q90 = NA_real_, n = as.integer(n_keep), cov = cov))
  vals <- as.numeric(signal[cbind(rr[keep], cc[keep])])
  vals <- vals[is.finite(vals)]
  n_vals <- length(vals)
  if (n_vals < min_n) return(list(ok = FALSE, mean = NA_real_, q90 = NA_real_, n = as.integer(n_vals), cov = cov))
  list(
    ok = TRUE,
    mean = suppressWarnings(mean(vals, na.rm = TRUE)),
    q90 = suppressWarnings(stats::quantile(vals, probs = 0.9, na.rm = TRUE, names = FALSE, type = 8)),
    n = as.integer(n_vals),
    cov = cov
  )
}

ra_sample_shifted_weighted_signal <- function(coord_obj, signal, weight, dx = 0L, dy = 0L, min_n = 5L) {
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
  if (!is.finite(n_keep) || n_keep < min_n) return(list(ok = FALSE, mean = NA_real_, n = as.integer(n_keep), cov = cov))
  vals <- as.numeric(signal[cbind(rr[keep], cc[keep])])
  ww <- as.numeric(weight[cbind(rr[keep], cc[keep])])
  good <- is.finite(vals) & is.finite(ww) & ww > 0
  if (sum(good) < min_n) return(list(ok = FALSE, mean = NA_real_, n = as.integer(sum(good)), cov = cov))
  ww <- ww[good]
  vals <- vals[good]
  wsum <- sum(ww)
  if (!is.finite(wsum) || wsum <= 0) return(list(ok = FALSE, mean = NA_real_, n = as.integer(length(vals)), cov = cov))
  list(ok = TRUE, mean = sum(vals * ww) / wsum, n = as.integer(length(vals)), cov = cov)
}

ra_score_mask_shift <- function(mask_cache, signal, dx = 0L, dy = 0L) {
  if (is.matrix(mask_cache)) mask_cache <- ra_build_mask_bands(mask_cache, band_px = 4L)
  if (is.null(mask_cache) || is.null(mask_cache$inside_all)) return(-Inf)
  dx <- as.integer(round(dx))
  dy <- as.integer(round(dy))
  in_all <- ra_sample_shifted_stats(mask_cache$inside_all, signal, dx = dx, dy = dy, min_n = 10L)
  in_relaxed <- if (isTRUE(in_all$ok)) in_all else ra_sample_shifted_stats(mask_cache$inside_all, signal, dx = dx, dy = dy, min_n = 1L)
  if (!isTRUE(in_relaxed$ok)) return(-Inf)
  coverage_all <- pmax(0, pmin(1, in_relaxed$cov))
  if (coverage_all < 0.01) return(-Inf)
  band_n <- min(length(mask_cache$inside_rings), length(mask_cache$outside_rings))
  w <- if (band_n > 0L) 1 / seq_len(band_n) else numeric(0)
  grad_num <- 0
  grad_den <- 0
  cov_num <- 0
  cov_den <- 0
  out_sum <- 0
  out_n <- 0
  if (band_n > 0L) {
    for (k in seq_len(band_n)) {
      in_k <- ra_sample_shifted_stats(mask_cache$inside_rings[[k]], signal, dx = dx, dy = dy, min_n = 3L)
      out_k <- ra_sample_shifted_stats(mask_cache$outside_rings[[k]], signal, dx = dx, dy = dy, min_n = 3L)
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
  global_contrast <- in_relaxed$mean - out_mean
  q_boost <- in_relaxed$q90 - out_mean
  if (!is.finite(q_boost)) q_boost <- 0
  core <- ra_sample_shifted_stats(mask_cache$inside_core, signal, dx = dx, dy = dy, min_n = 5L)
  core_term <- if (isTRUE(core$ok)) (core$mean - in_relaxed$mean) else 0
  coverage <- pmax(0.01, pmin(1, if (has_grad) pmin(coverage_all, pmax(0, cov_num / cov_den)) else coverage_all))
  base <- 0.60 * global_contrast + 0.25 * grad_contrast + 0.10 * q_boost + 0.05 * core_term
  score <- base - 1.20 * (1 - coverage)^2
  if (!is.finite(score)) score <- global_contrast - 1.50 * (1 - coverage_all)^2
  if (!is.finite(score)) return(-Inf)
  score
}

ra_build_histology_mask_fit_cache <- function(mask, band_px = 6L) {
  mask <- matrix(as.logical(mask), nrow = nrow(mask), ncol = ncol(mask))
  bands <- ra_build_mask_bands(mask, band_px = band_px)
  center_weight <- matrix(0, nrow = nrow(mask), ncol = ncol(mask))
  boundary_cache <- ra_build_mask_boundary_distance_cache(mask)
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

ra_score_mask_centrality_shift <- function(mask_cache, signal, dx = 0L, dy = 0L) {
  if (is.null(mask_cache) || is.null(mask_cache$inside_all)) return(-Inf)
  dx <- as.integer(round(dx))
  dy <- as.integer(round(dy))
  inside_all <- ra_sample_shifted_stats(mask_cache$inside_all, signal, dx = dx, dy = dy, min_n = 10L)
  if (!isTRUE(inside_all$ok)) inside_all <- ra_sample_shifted_stats(mask_cache$inside_all, signal, dx = dx, dy = dy, min_n = 1L)
  if (!isTRUE(inside_all$ok)) return(-Inf)
  center_term <- if (!is.null(mask_cache$center_weight)) {
    ra_sample_shifted_weighted_signal(mask_cache$inside_all, signal, mask_cache$center_weight, dx = dx, dy = dy, min_n = 5L)
  } else {
    list(ok = FALSE, mean = NA_real_, cov = inside_all$cov)
  }
  core_term <- ra_sample_shifted_stats(mask_cache$inside_core, signal, dx = dx, dy = dy, min_n = 5L)
  boundary_in <- if (length(mask_cache$inside_rings) >= 1L) ra_sample_shifted_stats(mask_cache$inside_rings[[1]], signal, dx = dx, dy = dy, min_n = 3L) else list(ok = FALSE, mean = NA_real_, cov = 0)
  boundary_out <- if (length(mask_cache$outside_rings) >= 1L) ra_sample_shifted_stats(mask_cache$outside_rings[[1]], signal, dx = dx, dy = dy, min_n = 3L) else list(ok = FALSE, mean = NA_real_, cov = 0)
  center_mean <- if (isTRUE(center_term$ok)) center_term$mean else if (isTRUE(core_term$ok)) core_term$mean else inside_all$mean
  boundary_in_mean <- if (isTRUE(boundary_in$ok)) boundary_in$mean else inside_all$mean
  boundary_out_mean <- if (isTRUE(boundary_out$ok)) boundary_out$mean else suppressWarnings(mean(as.numeric(signal), na.rm = TRUE))
  if (!is.finite(boundary_out_mean)) boundary_out_mean <- 0
  centeredness <- center_mean - boundary_in_mean
  inside_vs_outside <- center_mean - boundary_out_mean
  core_boost <- if (isTRUE(core_term$ok)) core_term$mean - inside_all$mean else 0
  global_contrast <- inside_all$mean - boundary_out_mean
  coverage <- pmax(0.01, pmin(1, min(
    if (is.finite(inside_all$cov)) inside_all$cov else 0,
    if (isTRUE(center_term$ok) && is.finite(center_term$cov)) center_term$cov else if (is.finite(inside_all$cov)) inside_all$cov else 0
  )))
  base <- 0.45 * centeredness + 0.30 * inside_vs_outside + 0.15 * core_boost + 0.10 * global_contrast
  score <- base - 1.10 * (1 - coverage)^2
  if (!is.finite(score)) return(-Inf)
  score
}

ra_build_distance_cache_from_edge_mask <- function(edge_mask, max_dist = 50L) {
  if (is.null(edge_mask) || !is.matrix(edge_mask) || nrow(edge_mask) < 2L || ncol(edge_mask) < 2L) return(NULL)
  edge <- matrix(as.logical(edge_mask), nrow = nrow(edge_mask), ncol = ncol(edge_mask))
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
  d <- pmin(d, max_dist)
  list(dist = d, edge = edge, n_edge = as.integer(n_edge))
}

ra_build_mask_boundary_distance_cache <- function(mask) {
  if (is.null(mask) || !is.matrix(mask)) return(NULL)
  mask <- matrix(as.logical(mask), nrow = nrow(mask), ncol = ncol(mask))
  if (sum(mask, na.rm = TRUE) < 10L) return(NULL)
  boundary <- mask & !ra_erode8(mask)
  ra_build_distance_cache_from_edge_mask(boundary)
}

ra_labels_to_boundary_mask <- function(labels, valid_mask = NULL) {
  labs <- matrix(as.character(labels), nrow = nrow(labels), ncol = ncol(labels))
  if (is.null(valid_mask)) valid_mask <- matrix(TRUE, nrow = nrow(labs), ncol = ncol(labs))
  valid_mask <- matrix(as.logical(valid_mask), nrow = nrow(labs), ncol = ncol(labs))
  nr <- nrow(labs)
  nc <- ncol(labs)
  boundary <- matrix(FALSE, nrow = nr, ncol = nc)
  right_lab <- cbind(labs[, -1, drop = FALSE], labs[, nc, drop = FALSE])
  down_lab <- rbind(labs[-1, , drop = FALSE], labs[nr, , drop = FALSE])
  right_valid <- cbind(valid_mask[, -1, drop = FALSE], FALSE)
  down_valid <- rbind(valid_mask[-1, , drop = FALSE], rep(FALSE, nc))
  diff_right <- valid_mask & right_valid & (labs != right_lab)
  diff_down <- valid_mask & down_valid & (labs != down_lab)
  boundary[diff_right | diff_down] <- TRUE
  boundary
}

ra_box_sum_matrix <- function(mat, radius = 1L) {
  mat <- matrix(as.numeric(mat), nrow = nrow(mat), ncol = ncol(mat))
  mat[!is.finite(mat)] <- 0
  radius <- max(0L, as.integer(radius))
  if (radius <= 0L) return(mat)
  nr <- nrow(mat)
  nc <- ncol(mat)
  padded <- matrix(0, nrow = nr + (2L * radius), ncol = nc + (2L * radius))
  padded[seq.int(radius + 1L, radius + nr), seq.int(radius + 1L, radius + nc)] <- mat
  kernel <- rep(1, (2L * radius) + 1L)
  sum_vec <- function(v) {
    out <- as.numeric(stats::filter(v, kernel, sides = 2))
    out[!is.finite(out)] <- 0
    out
  }
  tmp <- t(apply(padded, 1, sum_vec))
  out_pad <- apply(tmp, 2, sum_vec)
  if (!is.matrix(out_pad)) out_pad <- matrix(out_pad, nrow = nrow(padded), ncol = ncol(padded))
  out <- out_pad[seq.int(radius + 1L, radius + nr), seq.int(radius + 1L, radius + nc), drop = FALSE]
  out[!is.finite(out)] <- 0
  out
}

ra_label_map_diagnostics <- function(labels, valid_mask = NULL) {
  labs <- matrix(as.character(labels), nrow = nrow(labels), ncol = ncol(labels))
  if (is.null(valid_mask)) valid_mask <- matrix(TRUE, nrow = nrow(labs), ncol = ncol(labs))
  valid_mask <- matrix(as.logical(valid_mask), nrow = nrow(labs), ncol = ncol(labs))
  keep <- valid_mask & !is.na(labs) & nzchar(labs)
  n_valid <- suppressWarnings(sum(keep, na.rm = TRUE))
  if (!is.finite(n_valid) || n_valid <= 0L) {
    return(list(
      n_valid = 0L,
      n_labels = 0L,
      boundary_fraction = 0,
      component_count = 0L,
      component_density = 0,
      largest_label_fraction = 0
    ))
  }
  labs_valid <- labs[keep]
  uniq <- sort(unique(labs_valid))
  boundary_mask <- ra_labels_to_boundary_mask(labs, valid_mask = keep)
  component_count <- 0L
  label_frac <- numeric(0)
  for (lab in uniq) {
    mask <- keep & labs == lab
    comps <- ra_connected_components(mask)
    component_count <- component_count + length(comps$sizes)
    label_frac <- c(label_frac, suppressWarnings(sum(mask, na.rm = TRUE)) / n_valid)
  }
  list(
    n_valid = as.integer(n_valid),
    n_labels = length(uniq),
    boundary_fraction = suppressWarnings(sum(boundary_mask & keep, na.rm = TRUE)) / n_valid,
    component_count = as.integer(component_count),
    component_density = component_count / n_valid,
    largest_label_fraction = if (length(label_frac) > 0L) max(label_frac) else 0
  )
}

ra_majority_smooth_labels <- function(labels, valid_mask = NULL, radius = 2L, n_passes = 2L, min_fraction = 0.55) {
  labs <- matrix(as.character(labels), nrow = nrow(labels), ncol = ncol(labels))
  if (is.null(valid_mask)) valid_mask <- matrix(TRUE, nrow = nrow(labs), ncol = ncol(labs))
  valid_mask <- matrix(as.logical(valid_mask), nrow = nrow(labs), ncol = ncol(labs))
  radius <- max(0L, as.integer(radius))
  n_passes <- max(0L, as.integer(n_passes))
  min_fraction <- suppressWarnings(as.numeric(min_fraction)[1])
  if (!is.finite(min_fraction)) min_fraction <- 0.55
  if (radius <= 0L || n_passes <= 0L) return(labs)
  for (pass_idx in seq_len(n_passes)) {
    valid_count <- ra_box_sum_matrix(valid_mask, radius = radius)
    best_count <- matrix(-Inf, nrow = nrow(labs), ncol = ncol(labs))
    best_label <- labs
    uniq <- sort(unique(labs[valid_mask & !is.na(labs) & nzchar(labs)]))
    if (length(uniq) == 0L) break
    for (lab in uniq) {
      cnt <- ra_box_sum_matrix(valid_mask & labs == lab, radius = radius)
      replace <- cnt > best_count
      best_count[replace] <- cnt[replace]
      best_label[replace] <- lab
    }
    best_fraction <- ifelse(valid_count > 0, best_count / valid_count, 0)
    change <- valid_mask & is.finite(best_fraction) & (best_fraction >= min_fraction) & (best_label != labs)
    if (!any(change, na.rm = TRUE)) break
    labs[change] <- best_label[change]
  }
  labs
}

ra_merge_small_label_components <- function(labels, valid_mask = NULL, min_component_size = 64L, min_component_fraction = 0.005, max_passes = 2L, max_small_components = 256L) {
  labs <- matrix(as.character(labels), nrow = nrow(labels), ncol = ncol(labels))
  if (is.null(valid_mask)) valid_mask <- matrix(TRUE, nrow = nrow(labs), ncol = ncol(labs))
  valid_mask <- matrix(as.logical(valid_mask), nrow = nrow(labs), ncol = ncol(labs))
  min_component_size <- max(1L, as.integer(min_component_size))
  min_component_fraction <- suppressWarnings(as.numeric(min_component_fraction)[1])
  if (!is.finite(min_component_fraction) || min_component_fraction < 0) min_component_fraction <- 0
  max_passes <- max(1L, as.integer(max_passes))
  max_small_components <- max(1L, as.integer(max_small_components))
  for (pass_idx in seq_len(max_passes)) {
    changed_any <- FALSE
    uniq <- sort(unique(labs[valid_mask & !is.na(labs) & nzchar(labs)]))
    if (length(uniq) == 0L) break
    for (lab in uniq) {
      label_mask <- valid_mask & labs == lab
      label_n <- suppressWarnings(sum(label_mask, na.rm = TRUE))
      if (!is.finite(label_n) || label_n <= 0L) next
      comp_min <- max(min_component_size, as.integer(ceiling(label_n * min_component_fraction)))
      comps <- ra_connected_components(label_mask)
      if (length(comps$sizes) == 0L) next
      small_ids <- which(comps$sizes < comp_min)
      if (length(small_ids) == 0L) next
      if (length(small_ids) > max_small_components) next
      for (cid in small_ids) {
        comp_mask <- comps$labels == cid
        if (!any(comp_mask, na.rm = TRUE)) next
        neigh_mask <- ra_dilate8(comp_mask) & valid_mask & !comp_mask
        neigh_labels <- labs[neigh_mask]
        neigh_labels <- neigh_labels[!is.na(neigh_labels) & nzchar(neigh_labels) & neigh_labels != lab]
        if (length(neigh_labels) == 0L) next
        tab <- sort(table(neigh_labels), decreasing = TRUE)
        replacement <- names(tab)[1]
        if (is.na(replacement) || !nzchar(replacement)) next
        labs[comp_mask] <- replacement
        changed_any <- TRUE
      }
    }
    if (!changed_any) break
  }
  labs
}

ra_regularize_label_map <- function(labels, valid_mask = NULL, smooth_radius = 3L, smooth_passes = 3L, smooth_min_fraction = 0.50, min_component_size = 64L, min_component_fraction = 0.005, component_passes = 0L) {
  labs <- matrix(as.character(labels), nrow = nrow(labels), ncol = ncol(labels))
  if (is.null(valid_mask)) valid_mask <- matrix(TRUE, nrow = nrow(labs), ncol = ncol(labs))
  valid_mask <- matrix(as.logical(valid_mask), nrow = nrow(labs), ncol = ncol(labs))
  raw_diag <- ra_label_map_diagnostics(labs, valid_mask = valid_mask)
  labs_smoothed <- ra_majority_smooth_labels(
    labels = labs,
    valid_mask = valid_mask,
    radius = smooth_radius,
    n_passes = smooth_passes,
    min_fraction = smooth_min_fraction
  )
  labs_reg <- labs_smoothed
  if (is.finite(component_passes) && component_passes > 0L) {
    labs_reg <- ra_merge_small_label_components(
      labels = labs_smoothed,
      valid_mask = valid_mask,
      min_component_size = min_component_size,
      min_component_fraction = min_component_fraction,
      max_passes = component_passes
    )
  }
  reg_diag <- ra_label_map_diagnostics(labs_reg, valid_mask = valid_mask)
  changed_fraction <- suppressWarnings(mean((labs_reg != labs) & valid_mask, na.rm = TRUE))
  if (!is.finite(changed_fraction)) changed_fraction <- 0
  list(
    labels = labs_reg,
    valid = valid_mask,
    diagnostics = list(
      changed_fraction = changed_fraction,
      raw = raw_diag,
      regularized = reg_diag
    )
  )
}

ra_build_boundary_band_target <- function(mask, band_px = 3L, family = "pdata", field = "boundary_band", label = NULL) {
  mask <- matrix(as.logical(mask), nrow = nrow(mask), ncol = ncol(mask))
  if (sum(mask, na.rm = TRUE) < 5L) return(NULL)
  band_px <- max(1L, as.integer(band_px))
  bands <- ra_build_mask_bands(mask, band_px = max(2L * band_px, band_px + 1L))
  nr <- nrow(mask)
  nc <- ncol(mask)
  build_ring_mask <- function(ring_list, idx) {
    out <- matrix(FALSE, nrow = nr, ncol = nc)
    idx <- idx[idx >= 1L & idx <= length(ring_list)]
    if (length(idx) == 0L) return(out)
    for (ii in idx) out <- out | ra_coord_obj_to_mask(ring_list[[ii]], nr, nc)
    out
  }
  inner_mask <- build_ring_mask(bands$inside_rings, seq_len(min(band_px, length(bands$inside_rings))))
  outer_mask <- build_ring_mask(bands$outside_rings, seq_len(min(band_px, length(bands$outside_rings))))
  boundary_ring <- inner_mask | outer_mask
  context_idx <- seq.int(band_px + 1L, min(length(bands$inside_rings), 2L * band_px))
  inner_context <- build_ring_mask(bands$inside_rings, context_idx)
  outer_context <- build_ring_mask(bands$outside_rings, context_idx)
  context_mask <- inner_context | outer_context
  if (!any(context_mask, na.rm = TRUE)) context_mask <- boundary_ring
  edge_cache <- ra_build_distance_cache_from_edge_mask(mask & !ra_erode8(mask), max_dist = max(10L, 4L * band_px))
  list(
    family = family,
    field = field,
    label = label %||% field,
    type = "boundary_band",
    mask = mask,
    mask_cache = ra_build_histology_mask_fit_cache(mask, band_px = band_px),
    inner_mask = inner_mask,
    outer_mask = outer_mask,
    boundary_ring = boundary_ring,
    context_mask = context_mask,
    inner_context = inner_context,
    outer_context = outer_context,
    valid = boundary_ring | context_mask,
    edge_cache = edge_cache,
    band_px = band_px,
    local_shell_only = TRUE
  )
}

ra_build_contour_target <- function(mask, band_px = 1L, family = "pdata", field = "contour", label = NULL) {
  mask <- matrix(as.logical(mask), nrow = nrow(mask), ncol = ncol(mask))
  if (sum(mask, na.rm = TRUE) < 5L) return(NULL)
  band_px <- max(1L, as.integer(band_px))
  bands <- ra_build_mask_bands(mask, band_px = max(2L, band_px + 1L))
  nr <- nrow(mask)
  nc <- ncol(mask)
  build_ring_mask <- function(ring_list, idx) {
    out <- matrix(FALSE, nrow = nr, ncol = nc)
    idx <- idx[idx >= 1L & idx <= length(ring_list)]
    if (length(idx) == 0L) return(out)
    for (ii in idx) out <- out | ra_coord_obj_to_mask(ring_list[[ii]], nr, nc)
    out
  }
  boundary_line <- mask & !ra_erode8(mask)
  inner_mask <- build_ring_mask(bands$inside_rings, seq_len(min(band_px, length(bands$inside_rings))))
  outer_mask <- build_ring_mask(bands$outside_rings, seq_len(min(band_px, length(bands$outside_rings))))
  shell_mask <- inner_mask | outer_mask
  if (!any(shell_mask, na.rm = TRUE)) shell_mask <- boundary_line
  edge_cache <- ra_build_distance_cache_from_edge_mask(boundary_line, max_dist = max(10L, 4L * band_px))
  list(
    family = family,
    field = field,
    label = label %||% field,
    type = "contour_line",
    mask = mask,
    mask_cache = ra_build_histology_mask_fit_cache(mask, band_px = band_px),
    boundary_line = boundary_line,
    boundary_ring = boundary_line,
    inner_mask = inner_mask,
    outer_mask = outer_mask,
    shell_mask = shell_mask,
    context_mask = shell_mask,
    valid = boundary_line | shell_mask,
    edge_cache = edge_cache,
    band_px = band_px,
    local_shell_only = TRUE
  )
}

ra_build_edge_distance_cache <- function(signal, q = 0.85) {
  edge_info <- ra_build_edge_mask(signal, q = q)
  if (is.null(edge_info) || is.null(edge_info$edge)) return(NULL)
  out <- ra_build_distance_cache_from_edge_mask(edge_info$edge, max_dist = 50L)
  if (is.null(out)) return(NULL)
  out$threshold <- edge_info$threshold
  out
}

ra_score_edge_mask_against_distance <- function(edge_mask, edge_cache) {
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

ra_boundary_similarity <- function(mask_a, mask_b) {
  if (is.null(mask_a) || is.null(mask_b) || !is.matrix(mask_a) || !is.matrix(mask_b)) return(NA_real_)
  mask_a <- matrix(as.logical(mask_a), nrow = nrow(mask_a), ncol = ncol(mask_a))
  mask_b <- matrix(as.logical(mask_b), nrow = nrow(mask_b), ncol = ncol(mask_b))
  if (!all(dim(mask_a) == dim(mask_b))) return(NA_real_)
  edge_a <- mask_a & !ra_erode8(mask_a)
  edge_b <- mask_b & !ra_erode8(mask_b)
  if (sum(edge_a, na.rm = TRUE) < 5L || sum(edge_b, na.rm = TRUE) < 5L) return(NA_real_)
  cache_a <- ra_build_mask_boundary_distance_cache(mask_a)
  cache_b <- ra_build_mask_boundary_distance_cache(mask_b)
  if (is.null(cache_a) || is.null(cache_b)) return(NA_real_)
  ab <- cache_b$dist[edge_a]
  ba <- cache_a$dist[edge_b]
  ab <- ab[is.finite(ab)]
  ba <- ba[is.finite(ba)]
  if (length(ab) < 5L || length(ba) < 5L) return(NA_real_)
  sim_ab <- mean(1 / (1 + ab))
  sim_ba <- mean(1 / (1 + ba))
  mean(c(sim_ab, sim_ba), na.rm = TRUE)
}

ra_score_binary_masks <- function(signal_mask, target_mask) {
  if (is.null(signal_mask) || is.null(target_mask)) {
    return(list(score = NA_real_, dice = NA_real_, jaccard = NA_real_, boundary = NA_real_, overlap = 0L))
  }
  signal_mask <- matrix(as.logical(signal_mask), nrow = nrow(signal_mask), ncol = ncol(signal_mask))
  target_mask <- matrix(as.logical(target_mask), nrow = nrow(target_mask), ncol = ncol(target_mask))
  if (!all(dim(signal_mask) == dim(target_mask))) {
    return(list(score = NA_real_, dice = NA_real_, jaccard = NA_real_, boundary = NA_real_, overlap = 0L))
  }
  n_signal <- sum(signal_mask, na.rm = TRUE)
  n_target <- sum(target_mask, na.rm = TRUE)
  if (!is.finite(n_signal) || !is.finite(n_target) || n_signal < 5L || n_target < 5L) {
    return(list(score = NA_real_, dice = NA_real_, jaccard = NA_real_, boundary = NA_real_, overlap = 0L))
  }
  inter <- sum(signal_mask & target_mask, na.rm = TRUE)
  union <- sum(signal_mask | target_mask, na.rm = TRUE)
  dice <- (2 * inter) / (n_signal + n_target)
  jaccard <- inter / union
  boundary <- ra_boundary_similarity(signal_mask, target_mask)
  score <- ra_combine_weighted_scores(c(dice, jaccard, boundary), c(0.45, 0.20, 0.35))
  list(
    score = score,
    dice = dice,
    jaccard = jaccard,
    boundary = boundary,
    overlap = as.integer(inter)
  )
}

ra_combine_weighted_scores <- function(values, weights = NULL) {
  vals <- as.numeric(values)
  if (length(vals) == 0L) return(NA_real_)
  if (is.null(weights) || length(weights) != length(vals)) weights <- rep(1, length(vals))
  weights <- as.numeric(weights)
  keep <- is.finite(vals) & is.finite(weights) & weights > 0
  if (!any(keep)) return(NA_real_)
  sum(vals[keep] * weights[keep]) / sum(weights[keep])
}

ra_build_structure_signal_map <- function(signal, valid_mask = NULL, sigma = 4, positive_only = FALSE) {
  if (is.null(signal) || !is.matrix(signal) || nrow(signal) < 2L || ncol(signal) < 2L) return(NULL)
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
  sig_masked <- sig
  sig_masked[!valid] <- NA_real_
  smooth <- ra_gaussian_smooth_matrix(sig_masked, sigma = sigma)
  detail <- sig - smooth
  detail[!valid] <- NA_real_
  detail_work <- if (isTRUE(positive_only)) pmax(detail, 0) else abs(detail)
  right <- cbind(sig[, -1, drop = FALSE], sig[, nx, drop = FALSE])
  down <- rbind(sig[-1, , drop = FALSE], sig[ny, , drop = FALSE])
  right_valid <- cbind(valid[, -1, drop = FALSE], FALSE)
  down_valid <- rbind(valid[-1, , drop = FALSE], rep(FALSE, nx))
  grad <- abs(right - sig) + abs(down - sig)
  grad[!(valid & (right_valid | down_valid))] <- NA_real_
  detail_sig <- matrix(ra_rescale01(as.vector(detail_work), enhance = FALSE), nrow = ny, ncol = nx)
  grad_sig <- matrix(ra_rescale01(as.vector(grad), enhance = FALSE), nrow = ny, ncol = nx)
  out <- (0.70 * detail_sig) + (0.30 * grad_sig)
  out[!is.finite(out)] <- 0
  out[!valid] <- 0
  list(signal = out, valid = valid, detail = detail, gradient = grad)
}

ra_score_spatial_signal_map <- function(signal, valid_mask = NULL) {
  if (is.null(signal) || !is.matrix(signal)) return(-Inf)
  sig <- suppressWarnings(matrix(as.numeric(signal), nrow = nrow(signal), ncol = ncol(signal)))
  sig[!is.finite(sig)] <- 0
  valid <- if (is.null(valid_mask)) is.finite(signal) else matrix(as.logical(valid_mask), nrow = nrow(signal), ncol = ncol(signal))
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
  if (!all(is.finite(c(qv, qg)))) return(-Inf)
  contrast_score <- (qv[2] - qv[1]) + 0.5 * (qv[3] - qv[2])
  edge_score <- qg[2] + 0.5 * qg[1]
  0.60 * edge_score + 0.40 * contrast_score
}

ra_score_spatial_coherence_map <- function(signal, valid_mask = NULL, sigma = 2) {
  if (is.null(signal) || !is.matrix(signal)) return(-Inf)
  sig <- suppressWarnings(matrix(as.numeric(signal), nrow = nrow(signal), ncol = ncol(signal)))
  sig[!is.finite(sig)] <- 0
  valid <- if (is.null(valid_mask)) is.finite(signal) else matrix(as.logical(valid_mask), nrow = nrow(signal), ncol = ncol(signal))
  valid[!is.finite(valid)] <- FALSE
  if (sum(valid, na.rm = TRUE) < 25L) return(-Inf)
  sig_masked <- sig
  sig_masked[!valid] <- NA_real_
  smooth <- ra_gaussian_smooth_matrix(sig_masked, sigma = sigma)
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
  qv <- suppressWarnings(stats::quantile(raw_vals, probs = c(0.50, 0.90, 0.99), na.rm = TRUE, names = FALSE, type = 8))
  contrast <- if (all(is.finite(qv))) (qv[2] - qv[1]) + 0.35 * (qv[3] - qv[2]) else 0
  structure_score <- ra_score_spatial_signal_map(signal, valid_mask = valid)
  0.45 * corr + 0.30 * low_freq_frac + 0.15 * contrast + 0.10 * structure_score
}

ra_suggest_histology_fusion_indices <- function(obj, msi_geom, n_features = 12L, candidate_pool = NULL, intensity_transform = "none", max_cells = 1.5e7) {
  n_features <- as.integer(min(64L, max(2L, n_features)))
  mz_axis <- suppressWarnings(as.numeric(Cardinal::mz(obj)))
  n_feat <- length(mz_axis)
  n_pix <- ncol(obj)
  if (n_feat < 2L || n_pix < 2L) return(integer(0))
  if (!is.finite(candidate_pool) || candidate_pool < n_features) candidate_pool <- max(24L, 4L * n_features)
  candidate_pool <- as.integer(min(n_feat, max(n_features, candidate_pool)))
  sample_pix <- floor(max_cells / max(1L, n_feat))
  sample_pix <- max(200L, min(n_pix, sample_pix))
  pix_idx <- if (sample_pix < n_pix) sort(unique(as.integer(seq(1, n_pix, length.out = sample_pix)))) else seq_len(n_pix)
  sp <- try(as.matrix(Cardinal::spectra(obj)[, pix_idx, drop = FALSE]), silent = TRUE)
  if (inherits(sp, "try-error") || nrow(sp) < 2L || ncol(sp) < 2L) return(integer(0))
  storage.mode(sp) <- "double"
  sp[!is.finite(sp)] <- NA_real_
  if (!identical(intensity_transform, "none")) sp <- ra_transform_intensity(sp, intensity_transform)
  mean_intensity <- rowMeans(abs(sp), na.rm = TRUE)
  sd_intensity <- apply(sp, 1, stats::sd, na.rm = TRUE)
  nz <- rowMeans(is.finite(sp) & (abs(sp) > 0), na.rm = TRUE)
  mean_intensity[!is.finite(mean_intensity)] <- 0
  sd_intensity[!is.finite(sd_intensity)] <- 0
  nz[!is.finite(nz)] <- 0
  norm01 <- function(v) {
    v <- as.numeric(v)
    v[!is.finite(v)] <- 0
    rng <- range(v, na.rm = TRUE)
    if (!all(is.finite(rng)) || diff(rng) <= 1e-12) return(rep(0, length(v)))
    (v - rng[1]) / diff(rng)
  }
  mean_norm <- norm01(mean_intensity)
  sd_norm <- norm01(sd_intensity)
  nz_norm <- pmin(1, pmax(0, nz))
  score <- (0.45 * mean_norm + 0.45 * sd_norm + 0.10 * nz_norm)
  ord <- order(score, decreasing = TRUE, na.last = NA)
  ord <- ord[is.finite(score[ord])]
  if (length(ord) == 0L) return(integer(0))
  cand <- ord[seq_len(min(length(ord), max(candidate_pool, 96L)))]
  sp_eval <- try(as.matrix(Cardinal::spectra(obj)[cand, , drop = FALSE]), silent = TRUE)
  if (!inherits(sp_eval, "try-error") && !is.null(sp_eval)) {
    if (nrow(sp_eval) != length(cand) && ncol(sp_eval) == length(cand)) sp_eval <- t(sp_eval)
    if (nrow(sp_eval) == length(cand)) {
      storage.mode(sp_eval) <- "double"
      if (!identical(intensity_transform, "none")) sp_eval <- ra_transform_intensity(sp_eval, intensity_transform)
      sp_eval[!is.finite(sp_eval)] <- 0
      spatial_score <- rep(-Inf, length(cand))
      for (ii in seq_along(cand)) {
        sig_payload <- ra_build_numeric_signal_matrix(sp_eval[ii, ], msi_geom, intensity_transform = "none", apply_transform = FALSE, gaussian_smooth = FALSE, enhance = FALSE)
        spatial_score[ii] <- ra_combine_weighted_scores(
          c(
            ra_score_spatial_coherence_map(sig_payload$signal, sig_payload$valid, sigma = 2),
            ra_score_spatial_signal_map(sig_payload$signal, sig_payload$valid)
          ),
          c(0.75, 0.25)
        )
      }
      spatial_score[!is.finite(spatial_score)] <- min(spatial_score[is.finite(spatial_score)], 0)
      spatial_norm <- norm01(spatial_score)
      score[cand] <- 0.45 * score[cand] + 0.55 * spatial_norm
      ord <- order(score, decreasing = TRUE, na.last = NA)
      ord <- ord[is.finite(score[ord])]
    }
  }
  unique(ord[seq_len(min(length(ord), candidate_pool))])
}

ra_suggest_single_ion_indices <- function(obj, msi_geom, strategy = "most_intense_spatial", intensity_transform = "asinh", n_candidates = 5L) {
  strategy <- tolower(trimws(as.character(strategy)[1] %||% "most_intense_spatial"))
  n_candidates <- as.integer(max(1L, n_candidates))
  mz_n <- try(length(Cardinal::mz(obj)), silent = TRUE)
  if (inherits(mz_n, "try-error") || !is.finite(mz_n) || mz_n < 1L) return(integer(0))

  pix_n <- try(nrow(as.data.frame(Cardinal::coord(obj))), silent = TRUE)
  if (inherits(pix_n, "try-error") || !is.finite(pix_n) || pix_n < 1L) return(integer(0))
  pix_idx <- seq_len(min(5000L, pix_n))
  sp <- try(as.matrix(Cardinal::spectra(obj)[, pix_idx, drop = FALSE]), silent = TRUE)
  if (inherits(sp, "try-error") || is.null(sp)) return(integer(0))
  if (ncol(sp) != length(pix_idx) && nrow(sp) == length(pix_idx)) sp <- t(sp)
  storage.mode(sp) <- "double"
  sp[!is.finite(sp)] <- 0
  if (!identical(intensity_transform, "none")) sp <- ra_transform_intensity(sp, intensity_transform)

  mean_intensity <- rowMeans(sp, na.rm = TRUE)
  sd_intensity <- apply(sp, 1, stats::sd, na.rm = TRUE)
  nz_frac <- rowMeans(is.finite(sp) & abs(sp) > 0, na.rm = TRUE)
  mean_intensity[!is.finite(mean_intensity)] <- 0
  sd_intensity[!is.finite(sd_intensity)] <- 0
  nz_frac[!is.finite(nz_frac)] <- 0

  norm01 <- function(v) {
    v <- as.numeric(v)
    v[!is.finite(v)] <- 0
    rng <- range(v, na.rm = TRUE)
    if (!all(is.finite(rng)) || diff(rng) <= 1e-12) return(rep(0, length(v)))
    (v - rng[1]) / diff(rng)
  }

  if (identical(strategy, "most_intense")) {
    score <- norm01(mean_intensity) * (0.25 + 0.75 * pmin(1, nz_frac * 2))
    ord <- order(score, decreasing = TRUE, na.last = NA)
    return(ord[seq_len(min(length(ord), n_candidates))])
  }
  if (identical(strategy, "most_variable")) {
    score <- norm01(sd_intensity) * (0.25 + 0.75 * pmin(1, nz_frac * 2))
    ord <- order(score, decreasing = TRUE, na.last = NA)
    return(ord[seq_len(min(length(ord), n_candidates))])
  }

  base_score <- (0.45 * norm01(mean_intensity) + 0.45 * norm01(pmax(sd_intensity, 0)) + 0.10 * pmin(1, nz_frac)) *
    (0.25 + 0.75 * pmin(1, nz_frac * 2))
  score <- base_score
  score[!is.finite(score)] <- 0
  ord <- order(score, decreasing = TRUE, na.last = NA)
  ord <- ord[seq_len(min(length(ord), 48L))]
  if (length(ord) == 0L) return(integer(0))

  spatial_score <- rep(-Inf, length(ord))
  for (ii in seq_along(ord)) {
    sig_payload <- ra_build_numeric_signal_matrix(sp[ord[ii], ], msi_geom, intensity_transform = "none", apply_transform = FALSE, gaussian_smooth = FALSE, enhance = FALSE)
    spatial_score[ii] <- ra_combine_weighted_scores(
      c(
        ra_score_spatial_coherence_map(sig_payload$signal, sig_payload$valid, sigma = 2),
        ra_score_spatial_signal_map(sig_payload$signal, sig_payload$valid)
      ),
      c(0.75, 0.25)
    )
  }
  spatial_score[!is.finite(spatial_score)] <- 0
  final_score <- 0.45 * base_score[ord] + 0.55 * norm01(spatial_score)
  final_score[!is.finite(final_score)] <- -Inf
  ord2 <- ord[order(final_score, decreasing = TRUE, na.last = NA)]
  unique(ord2[seq_len(min(length(ord2), n_candidates))])
}

ra_suggest_single_ion_index <- function(obj, msi_geom, strategy = "most_intense_spatial", intensity_transform = "asinh") {
  out <- ra_suggest_single_ion_indices(
    obj = obj,
    msi_geom = msi_geom,
    strategy = strategy,
    intensity_transform = intensity_transform,
    n_candidates = 1L
  )
  if (length(out) < 1L) return(integer(0))
  out[1]
}

ra_histology_signal_from_fused_ions <- function(obj, msi_geom, n_ions = 12L, intensity_transform = "asinh", gaussian_smooth = FALSE, smooth_sigma = 1) {
  cand_idx <- ra_suggest_histology_fusion_indices(obj, msi_geom, n_features = n_ions, candidate_pool = max(24L, 4L * n_ions), intensity_transform = intensity_transform)
  if (length(cand_idx) == 0L) stop("Could not identify informative ions for fused target.")
  sp <- as.matrix(Cardinal::spectra(obj)[cand_idx, , drop = FALSE])
  if (nrow(sp) != length(cand_idx) && ncol(sp) == length(cand_idx)) sp <- t(sp)
  storage.mode(sp) <- "double"
  sig_list <- vector("list", length(cand_idx))
  sig_score <- rep(-Inf, length(cand_idx))
  for (ii in seq_along(cand_idx)) {
    sig_payload <- ra_build_numeric_signal_matrix(
      sp[ii, ], msi_geom,
      intensity_transform = intensity_transform,
      apply_transform = TRUE,
      gaussian_smooth = gaussian_smooth,
      smooth_sigma = smooth_sigma,
      enhance = FALSE
    )
    sig_list[[ii]] <- sig_payload$signal
    sig_score[ii] <- ra_combine_weighted_scores(
      c(
        ra_score_spatial_coherence_map(sig_payload$signal, sig_payload$valid, sigma = 2),
        ra_score_spatial_signal_map(sig_payload$signal, sig_payload$valid)
      ),
      c(0.75, 0.25)
    )
  }
  sig_score[!is.finite(sig_score)] <- 0
  ord <- order(sig_score, decreasing = TRUE, na.last = NA)
  ord <- ord[seq_len(min(length(ord), n_ions))]
  sel_idx <- cand_idx[ord]
  sel_scores <- sig_score[ord]
  sel_signals <- sig_list[ord]
  w <- sel_scores - min(sel_scores, na.rm = TRUE)
  if (!all(is.finite(w)) || sum(w) <= 1e-12) w <- rep(1, length(sel_scores))
  w <- (w + 1e-6) / sum(w + 1e-6)
  fused <- matrix(0, nrow = msi_geom$ny, ncol = msi_geom$nx)
  for (ii in seq_along(sel_signals)) fused <- fused + (w[ii] * sel_signals[[ii]])
  ok <- is.finite(fused) & fused > 0
  if (any(ok)) {
    q <- suppressWarnings(stats::quantile(fused[ok], probs = 0.75, na.rm = TRUE, names = FALSE, type = 8))
    if (is.finite(q) && q > 0) fused <- pmax(fused - q, 0)
  }
  mx <- suppressWarnings(max(fused, na.rm = TRUE))
  if (is.finite(mx) && mx > 0) fused <- fused / mx
  mz_axis <- suppressWarnings(as.numeric(Cardinal::mz(obj)))
  sel_mz <- mz_axis[sel_idx]
  out <- ra_build_numeric_target(fused, is.finite(fused), family = "msi", field = sprintf("Fused MSI (%d ions)", length(sel_idx)), type = "mz_fusion")
  out$label <- sprintf("Fused MSI (%d ions)", length(sel_idx))
  out$mz_values <- sel_mz
  out$feature_index <- sel_idx
  out
}

ra_histology_signal_from_pca <- function(obj, msi_geom, component = 1L, n_components = 1L, n_features = 32L, intensity_transform = "asinh", gaussian_smooth = FALSE, smooth_sigma = 1) {
  sel_idx <- ra_suggest_histology_fusion_indices(obj, msi_geom, n_features = n_features, candidate_pool = max(48L, 4L * n_features), intensity_transform = intensity_transform)
  if (length(sel_idx) < 2L) stop("Could not identify enough informative ions for PCA target.")
  sp <- as.matrix(Cardinal::spectra(obj)[sel_idx, , drop = FALSE])
  if (nrow(sp) != length(sel_idx) && ncol(sp) == length(sel_idx)) sp <- t(sp)
  storage.mode(sp) <- "double"
  if (!identical(intensity_transform, "none")) sp <- ra_transform_intensity(sp, intensity_transform)
  sp[!is.finite(sp)] <- 0
  feat_mean <- rowMeans(sp)
  sp_centered <- sp - feat_mean
  feat_scale <- sqrt(rowMeans(sp_centered^2))
  feat_scale[!is.finite(feat_scale) | feat_scale <= 1e-8] <- 1
  sp_scaled <- sp_centered / feat_scale
  eig <- eigen(tcrossprod(sp_scaled) / max(1, ncol(sp_scaled) - 1L), symmetric = TRUE)
  keep_comp <- is.finite(eig$values) & (eig$values > 1e-10)
  values <- as.numeric(eig$values[keep_comp])
  vectors <- eig$vectors[, keep_comp, drop = FALSE]
  comp_ids <- seq.int(component, min(ncol(vectors), component + n_components - 1L))
  signal_list <- vector("list", length(comp_ids))
  weight_vec <- rep(0, length(comp_ids))
  load_mat <- matrix(0, nrow = length(sel_idx), ncol = length(comp_ids))
  avg_signal <- suppressWarnings(colMeans(sp))
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
    sig_payload <- ra_build_numeric_signal_matrix(score_vec, msi_geom, intensity_transform = "none", apply_transform = FALSE, gaussian_smooth = gaussian_smooth, smooth_sigma = smooth_sigma, enhance = FALSE)
    signal_list[[ii]] <- sig_payload$signal
    var_frac <- values[comp_use] / sum(values)
    spatial_weight <- ra_combine_weighted_scores(
      c(
        ra_score_spatial_coherence_map(sig_payload$signal, sig_payload$valid, sigma = 2),
        ra_score_spatial_signal_map(sig_payload$signal, sig_payload$valid)
      ),
      c(0.75, 0.25)
    )
    weight_vec[ii] <- max(1e-6, ifelse(is.finite(var_frac), var_frac, 0)) * (1 + max(ifelse(is.finite(spatial_weight), spatial_weight, 0), 0))
  }
  if (!any(is.finite(weight_vec) & weight_vec > 0)) weight_vec <- rep(1, length(signal_list))
  weight_vec <- weight_vec / sum(weight_vec)
  pca_signal <- matrix(0, nrow = msi_geom$ny, ncol = msi_geom$nx)
  for (ii in seq_along(signal_list)) pca_signal <- pca_signal + (weight_vec[ii] * signal_list[[ii]])
  pca_signal <- matrix(ra_rescale01(as.vector(pca_signal), enhance = FALSE), nrow = msi_geom$ny, ncol = msi_geom$nx)
  comp_label <- if (length(comp_ids) == 1L) sprintf("PC%d", comp_ids[1]) else sprintf("PC%d-%d", min(comp_ids), max(comp_ids))
  out <- ra_build_numeric_target(pca_signal, is.finite(pca_signal), family = "msi", field = comp_label, type = "pca")
  out$label <- comp_label
  out$component_ids <- comp_ids
  out$feature_index <- sel_idx
  out
}

ra_build_fused_target_from_indices <- function(obj, msi_geom, sel_idx, intensity_transform = "asinh", gaussian_smooth = FALSE, smooth_sigma = 1, fusion_mode = "support_aware_mean") {
  sel_idx <- unique(as.integer(sel_idx[is.finite(sel_idx)]))
  if (length(sel_idx) == 0L) stop("No ion indices supplied for fused target.")
  sp <- as.matrix(Cardinal::spectra(obj)[sel_idx, , drop = FALSE])
  if (nrow(sp) != length(sel_idx) && ncol(sp) == length(sel_idx)) sp <- t(sp)
  storage.mode(sp) <- "double"
  mz_axis <- suppressWarnings(as.numeric(Cardinal::mz(obj)))
  ion_signals <- vector("list", length(sel_idx))
  ion_weights <- rep(0, length(sel_idx))
  ion_coverage <- rep(0, length(sel_idx))
  support_mask_list <- vector("list", length(sel_idx))
  for (ii in seq_along(sel_idx)) {
    vals <- as.numeric(sp[ii, ])
    vals <- ra_transform_intensity(vals, intensity_transform)
    vals[!is.finite(vals)] <- NA_real_
    vals_norm <- ra_robust_rescale01(vals, q_lo = 0.05, q_hi = 0.95)
    sig_payload <- ra_build_numeric_signal_matrix(vals_norm, msi_geom, intensity_transform = "none", apply_transform = FALSE, gaussian_smooth = gaussian_smooth, smooth_sigma = smooth_sigma, enhance = FALSE)
    ion_signals[[ii]] <- sig_payload$signal
    ion_coverage[ii] <- suppressWarnings(mean(sig_payload$signal[sig_payload$valid] >= 0.20, na.rm = TRUE))
    support_mask_list[[ii]] <- sig_payload$signal >= 0.25
    ion_weights[ii] <- sqrt(max(ion_coverage[ii], 0))
  }
  keep <- is.finite(ion_coverage) & ion_coverage >= 0.01
  if (!any(keep)) keep <- rep(TRUE, length(sel_idx))
  ion_signals <- ion_signals[keep]
  support_mask_list <- support_mask_list[keep]
  ion_weights <- ion_weights[keep]
  sel_idx <- sel_idx[keep]
  ion_coverage <- ion_coverage[keep]
  if (length(sel_idx) == 0L) stop("No informative ions remained after fusion coverage filtering.")
  if (!all(is.finite(ion_weights)) || sum(ion_weights) <= 1e-12) ion_weights <- rep(1, length(sel_idx))
  ion_weights <- ion_weights / sum(ion_weights)

  nr <- msi_geom$ny
  nc <- msi_geom$nx
  fused_mean <- matrix(0, nrow = nr, ncol = nc)
  support_map <- matrix(0, nrow = nr, ncol = nc)
  rank_sum <- matrix(0, nrow = nr, ncol = nc)
  for (ii in seq_along(sel_idx)) {
    fused_mean <- fused_mean + ion_weights[ii] * ion_signals[[ii]]
    support_map <- support_map + ion_weights[ii] * matrix(as.numeric(support_mask_list[[ii]]), nrow = nr, ncol = nc)
    rank_sum <- rank_sum + ion_weights[ii] * matrix(ra_robust_rescale01(as.vector(ion_signals[[ii]]), q_lo = 0.10, q_hi = 0.90), nrow = nr, ncol = nc)
  }
  fused <- switch(
    tolower(trimws(as.character(fusion_mode)[1] %||% "support_aware_mean")),
    robust_weighted_mean = fused_mean,
    rank_fusion = rank_sum,
    support_aware_mean = fused_mean * sqrt(pmax(support_map, 0)),
    fused_mean * sqrt(pmax(support_map, 0))
  )
  fused <- matrix(ra_robust_rescale01(as.vector(fused), q_lo = 0.05, q_hi = 0.95), nrow = nr, ncol = nc)
  fused_valid <- is.finite(fused)
  fused_valid[!is.finite(fused_valid)] <- FALSE
  support_binary <- support_map >= 0.25
  comps <- ra_connected_components(support_binary)
  largest_component_fraction <- if (length(comps$sizes) > 0L) max(comps$sizes) / max(1L, sum(support_binary, na.rm = TRUE)) else 0
  support_fraction <- suppressWarnings(mean(support_binary, na.rm = TRUE))
  out <- ra_build_numeric_target(fused, fused_valid, family = "msi", field = sprintf("Fused MSI (%d ions)", length(sel_idx)), type = "mz_fusion")
  out$label <- sprintf("Fused MSI (%d ions)", length(sel_idx))
  out$mz_values <- mz_axis[sel_idx]
  out$feature_index <- sel_idx
  out$fusion_mode <- fusion_mode
  out$fusion_diagnostics <- data.frame(
    feature_index = sel_idx,
    mz = mz_axis[sel_idx],
    weight = ion_weights,
    coverage_fraction = ion_coverage,
    stringsAsFactors = FALSE
  )
  out$support_fraction <- support_fraction
  out$largest_component_fraction <- largest_component_fraction
  out$support_map <- support_map
  out
}

ra_build_pca_target_from_indices <- function(obj, msi_geom, sel_idx, component = 1L, n_components = 1L, intensity_transform = "asinh", gaussian_smooth = FALSE, smooth_sigma = 1) {
  sel_idx <- unique(as.integer(sel_idx[is.finite(sel_idx)]))
  if (length(sel_idx) < 2L) stop("Need at least two ion indices for PCA target.")
  sp <- as.matrix(Cardinal::spectra(obj)[sel_idx, , drop = FALSE])
  if (nrow(sp) != length(sel_idx) && ncol(sp) == length(sel_idx)) sp <- t(sp)
  storage.mode(sp) <- "double"
  if (!identical(intensity_transform, "none")) sp <- ra_transform_intensity(sp, intensity_transform)
  sp[!is.finite(sp)] <- 0
  feat_mean <- rowMeans(sp)
  sp_centered <- sp - feat_mean
  feat_scale <- sqrt(rowMeans(sp_centered^2))
  feat_scale[!is.finite(feat_scale) | feat_scale <= 1e-8] <- 1
  sp_scaled <- sp_centered / feat_scale
  eig <- eigen(tcrossprod(sp_scaled) / max(1, ncol(sp_scaled) - 1L), symmetric = TRUE)
  keep_comp <- is.finite(eig$values) & (eig$values > 1e-10)
  values <- as.numeric(eig$values[keep_comp])
  vectors <- eig$vectors[, keep_comp, drop = FALSE]
  comp_ids <- seq.int(component, min(ncol(vectors), component + n_components - 1L))
  if (length(comp_ids) == 0L) stop("No PCA components available for strict target.")
  signal_list <- vector("list", length(comp_ids))
  weight_vec <- rep(0, length(comp_ids))
  avg_signal <- suppressWarnings(colMeans(sp))
  for (ii in seq_along(comp_ids)) {
    comp_use <- comp_ids[ii]
    load_vec <- as.numeric(vectors[, comp_use])
    score_vec <- as.numeric(crossprod(load_vec, sp_scaled))
    orient <- suppressWarnings(stats::cor(score_vec, avg_signal, use = "pairwise.complete.obs", method = "spearman"))
    if ((is.finite(orient) && orient < 0) || (!is.finite(orient) && sum(load_vec, na.rm = TRUE) < 0)) {
      load_vec <- -load_vec
      score_vec <- -score_vec
    }
    sig_payload <- ra_build_numeric_signal_matrix(score_vec, msi_geom, intensity_transform = "none", apply_transform = FALSE, gaussian_smooth = gaussian_smooth, smooth_sigma = smooth_sigma, enhance = FALSE)
    signal_list[[ii]] <- sig_payload$signal
    var_frac <- values[comp_use] / sum(values)
    spatial_weight <- ra_combine_weighted_scores(
      c(
        ra_score_spatial_coherence_map(sig_payload$signal, sig_payload$valid, sigma = 2),
        ra_score_spatial_signal_map(sig_payload$signal, sig_payload$valid)
      ),
      c(0.75, 0.25)
    )
    weight_vec[ii] <- max(1e-6, ifelse(is.finite(var_frac), var_frac, 0)) * (1 + max(ifelse(is.finite(spatial_weight), spatial_weight, 0), 0))
  }
  if (!any(is.finite(weight_vec) & weight_vec > 0)) weight_vec <- rep(1, length(signal_list))
  weight_vec <- weight_vec / sum(weight_vec)
  pca_signal <- matrix(0, nrow = msi_geom$ny, ncol = msi_geom$nx)
  for (ii in seq_along(signal_list)) pca_signal <- pca_signal + (weight_vec[ii] * signal_list[[ii]])
  pca_signal <- matrix(ra_rescale01(as.vector(pca_signal), enhance = FALSE), nrow = msi_geom$ny, ncol = msi_geom$nx)
  comp_label <- if (length(comp_ids) == 1L) sprintf("PC%d", comp_ids[1]) else sprintf("PC%d-%d", min(comp_ids), max(comp_ids))
  out <- ra_build_numeric_target(pca_signal, is.finite(pca_signal), family = "msi", field = comp_label, type = "pca")
  out$label <- comp_label
  out$component_ids <- comp_ids
  out$feature_index <- sel_idx
  out
}

ra_pdata_values_are_categorical <- function(vals, max_levels = 128L) {
  vals_num <- suppressWarnings(as.numeric(vals))
  is_num <- is.numeric(vals) || is.integer(vals)
  if (is.logical(vals)) return(TRUE)
  if (!is_num) return(TRUE)
  unique_n <- length(unique(vals[!is.na(vals)]))
  integer_like <- all(is.na(vals_num) | abs(vals_num - round(vals_num)) < 1e-8)
  isTRUE(integer_like) && isTRUE(unique_n <= max_levels)
}

ra_build_pdata_target <- function(msi_geom, field, gaussian_smooth = FALSE, smooth_sigma = 1) {
  pd <- msi_geom$pdata
  if (!field %in% names(pd)) stop(sprintf("pData field not found: %s", field))
  vals <- pd[[field]]
  vals_num <- suppressWarnings(as.numeric(vals))
  is_num <- is.numeric(vals) || is.integer(vals)
  unique_n <- length(unique(vals[!is.na(vals)]))
  integer_like <- is_num && all(is.na(vals_num) | abs(vals_num - round(vals_num)) < 1e-8)
  uniq_num <- sort(unique(vals_num[is.finite(vals_num)]))
  binary_like <- (is.logical(vals) || integer_like) && length(uniq_num) <= 2L && all(uniq_num %in% c(0, 1))
  as_categorical <- (!is_num) || (integer_like && unique_n <= 128L)
  if (isTRUE(binary_like)) {
    sig_payload <- ra_build_numeric_signal_matrix(vals_num, msi_geom, intensity_transform = "none", apply_transform = FALSE, gaussian_smooth = gaussian_smooth, smooth_sigma = smooth_sigma, enhance = FALSE)
    mask <- sig_payload$signal > 0.5
    return(list(
      family = "pdata",
      field = field,
      signal = sig_payload$signal,
      valid = sig_payload$valid,
      type = "binary_mask",
      mask = mask,
      mask_cache = ra_build_histology_mask_fit_cache(mask),
      edge_cache = ra_build_mask_boundary_distance_cache(mask)
    ))
  }
  if (!as_categorical) {
    sig_payload <- ra_build_numeric_signal_matrix(vals_num, msi_geom, intensity_transform = "none", apply_transform = FALSE, gaussian_smooth = gaussian_smooth, smooth_sigma = smooth_sigma, enhance = FALSE)
    return(list(
      family = "pdata",
      field = field,
      signal = sig_payload$signal,
      valid = sig_payload$valid,
      type = "numeric",
      structure = ra_build_structure_signal_map(sig_payload$signal, sig_payload$valid, positive_only = TRUE),
      edge_cache = ra_build_edge_distance_cache(sig_payload$signal)
    ))
  }
  labs <- ra_normalize_labels(vals)
  lab_mat <- matrix(NA_character_, nrow = msi_geom$ny, ncol = msi_geom$nx)
  lab_mat[cbind(msi_geom$row_idx, msi_geom$x_norm)] <- labs
  list(
    family = "pdata",
    field = field,
    type = "categorical_label_map",
    labels = lab_mat,
    valid = !is.na(lab_mat)
  )
}

ra_regularize_categorical_pdata_target <- function(target, smooth_radius = 3L, smooth_passes = 3L, smooth_min_fraction = 0.50, min_component_size = 96L, min_component_fraction = 0.01, component_passes = 0L) {
  if (is.null(target) || !identical(target$type, "categorical_label_map")) {
    return(list(target = target, diagnostics = NULL))
  }
  reg <- ra_regularize_label_map(
    labels = target$labels,
    valid_mask = target$valid,
    smooth_radius = smooth_radius,
    smooth_passes = smooth_passes,
    smooth_min_fraction = smooth_min_fraction,
    min_component_size = min_component_size,
    min_component_fraction = min_component_fraction,
    component_passes = component_passes
  )
  target$labels <- reg$labels
  target$valid <- reg$valid
  target$regularization <- reg$diagnostics
  list(target = target, diagnostics = reg$diagnostics)
}

ra_make_binary_mask_target <- function(mask, field, label, family = "pdata", band_px = 3L) {
  mask <- matrix(as.logical(mask), nrow = nrow(mask), ncol = ncol(mask))
  if (sum(mask, na.rm = TRUE) < 5L) return(NULL)
  list(
    family = family,
    field = field,
    label = label,
    type = "binary_mask",
    signal = matrix(as.numeric(mask), nrow = nrow(mask), ncol = ncol(mask)),
    valid = matrix(TRUE, nrow = nrow(mask), ncol = ncol(mask)),
    mask = mask,
    mask_cache = ra_build_histology_mask_fit_cache(mask, band_px = band_px),
    edge_cache = ra_build_mask_boundary_distance_cache(mask)
  )
}

ra_strict_pdata_priority <- function(metrics, target_type, route = "direct", consensus_score = 1) {
  if (is.null(metrics)) return(-Inf)
  support_fraction <- suppressWarnings(as.numeric(metrics$support_fraction)[1])
  largest_component_fraction <- suppressWarnings(as.numeric(metrics$largest_component_fraction)[1])
  boundary_fraction <- suppressWarnings(as.numeric(metrics$boundary_fraction)[1])
  n_components <- suppressWarnings(as.numeric(metrics$n_components)[1])
  n_mask <- suppressWarnings(as.numeric(metrics$n_mask)[1])
  component_density <- if (is.finite(n_components) && is.finite(n_mask) && n_mask > 0) n_components / n_mask else Inf
  if (!all(is.finite(c(support_fraction, largest_component_fraction, boundary_fraction, component_density)))) return(-Inf)
  target_type <- tolower(trimws(as.character(target_type)[1] %||% "unknown"))
  route <- tolower(trimws(as.character(route)[1] %||% "direct"))
  is_boundary_like <- target_type %in% c("boundary_band", "contour_line")
  consensus_score <- suppressWarnings(as.numeric(consensus_score)[1])
  if (!is.finite(consensus_score)) consensus_score <- 0
  consensus_score <- max(0, min(1, consensus_score))

  if (identical(target_type, "binary_mask")) {
    if (support_fraction < 0.03 || support_fraction > 0.25) return(-Inf)
    if (largest_component_fraction < 0.03 || component_density > 0.25) return(-Inf)
    if (boundary_fraction < 0.05 || boundary_fraction > 0.90) return(-Inf)
  } else if (is_boundary_like) {
    if (support_fraction < 0.02 || support_fraction > 0.65) return(-Inf)
    if (largest_component_fraction < 0.01 || component_density > if (identical(target_type, "contour_line")) 0.30 else 0.35) return(-Inf)
    if (boundary_fraction < 0.02 || boundary_fraction > if (identical(target_type, "contour_line")) 0.80 else 0.90) return(-Inf)
  } else {
    if (support_fraction < 0.02 || support_fraction > 0.75) return(-Inf)
  }

  support_opt <- if (is_boundary_like) 0.12 else 0.08
  support_span <- if (is_boundary_like) 0.18 else 0.10
  support_score <- max(0, min(1, 1 - abs(support_fraction - support_opt) / support_span))
  boundary_opt <- if (identical(target_type, "contour_line")) 0.22 else if (identical(target_type, "boundary_band")) 0.35 else 0.28
  boundary_span <- if (identical(target_type, "contour_line")) 0.24 else if (identical(target_type, "boundary_band")) 0.40 else 0.32
  boundary_score <- max(0, min(1, 1 - abs(boundary_fraction - boundary_opt) / boundary_span))
  connectedness_score <- max(0, min(1, (largest_component_fraction - 0.02) / 0.48))
  fragmentation_score <- max(0, min(1, 1 - component_density / 0.18))
  route_bonus <- switch(
    route,
    consensus_boundary = 0.12,
    auto_boundary = 0.08,
    companion_boundary = 0.05,
    direct_mask = 0.02,
    direct = 0.01,
    0
  )
  type_bonus <- if (identical(target_type, "contour_line")) 0.10 else if (identical(target_type, "boundary_band")) 0.06 else 0.02
  as.numeric(
    0.38 * support_score +
      0.18 * boundary_score +
      0.18 * connectedness_score +
      0.14 * fragmentation_score +
      0.12 * consensus_score +
      route_bonus +
      type_bonus
  )
}

ra_build_strict_pdata_targets <- function(msi_geom, fields, boundary_band_px = 3L, max_targets = 6L, max_per_field = 3L) {
  fields <- unique(as.character(fields))
  fields <- fields[fields %in% names(msi_geom$pdata)]
  if (length(fields) == 0L) {
    return(list(targets = list(), catalog = data.frame(), field_diagnostics = data.frame()))
  }
  max_targets <- max(1L, as.integer(max_targets[1] %||% 6L))
  max_per_field <- max(1L, as.integer(max_per_field[1] %||% max_targets))
  candidate_targets <- list()
  catalog_rows <- vector("list", 0L)
  field_diag_rows <- vector("list", 0L)
  add_target <- function(key, target, source_field, source_label, target_type, route, metrics, consensus_score = 1) {
    priority_score <- ra_strict_pdata_priority(metrics, target_type = target_type, route = route, consensus_score = consensus_score)
    target_group <- sprintf("%s::%s", source_field, source_label)
    if (!is.null(target)) {
      target$source_field <- source_field
      target$source_label <- source_label
      target$target_group <- target_group
      target$route <- route
      target$support_fraction <- metrics$support_fraction
      target$largest_component_fraction <- metrics$largest_component_fraction
      target$boundary_fraction <- metrics$boundary_fraction
      target$n_components <- metrics$n_components
      target$n_mask <- metrics$n_mask
      target$component_density <- if (is.finite(metrics$n_components) && is.finite(metrics$n_mask) && metrics$n_mask > 0) metrics$n_components / metrics$n_mask else NA_real_
      target$consensus_score <- consensus_score
      target$priority_score <- priority_score
      candidate_targets[[key]] <<- target
    }
    catalog_rows[[length(catalog_rows) + 1L]] <<- data.frame(
      target_key = key,
      field = source_field,
      label = source_label,
      target_group = target_group,
      target_type = target_type,
      route = route,
      support_fraction = metrics$support_fraction,
      largest_component_fraction = metrics$largest_component_fraction,
      boundary_fraction = metrics$boundary_fraction,
      n_components = metrics$n_components,
      n_mask = metrics$n_mask,
      component_density = if (is.finite(metrics$n_components) && is.finite(metrics$n_mask) && metrics$n_mask > 0) metrics$n_components / metrics$n_mask else NA_real_,
      consensus_score = consensus_score,
      priority_score = priority_score,
      strict_admitted = as.logical(FALSE),
      stringsAsFactors = FALSE
    )
    invisible(NULL)
  }

  categorical_fields <- character(0)
  for (field in fields) {
    target_base <- ra_build_pdata_target(msi_geom, field, gaussian_smooth = FALSE, smooth_sigma = 1)
    if (!identical(target_base$type, "categorical_label_map")) {
      metrics <- if (!is.null(target_base$mask) && is.matrix(target_base$mask)) {
        ra_mask_shape_metrics(target_base$mask, valid_mask = target_base$valid %||% matrix(TRUE, nrow = nrow(target_base$mask), ncol = ncol(target_base$mask)))
      } else {
        list(support_fraction = NA_real_, largest_component_fraction = NA_real_, boundary_fraction = NA_real_, n_components = NA_integer_)
      }
      key <- paste0("strict::", field)
      priority_score <- if (identical(target_base$type, "binary_mask")) {
        ra_strict_pdata_priority(metrics, target_type = target_base$type, route = "direct", consensus_score = 1)
      } else {
        -Inf
      }
      target_base$source_field <- field
      target_base$source_label <- field
      target_base$target_group <- sprintf("%s::%s", field, field)
      target_base$route <- "direct"
      target_base$support_fraction <- metrics$support_fraction
      target_base$largest_component_fraction <- metrics$largest_component_fraction
      target_base$boundary_fraction <- metrics$boundary_fraction
      target_base$n_components <- metrics$n_components
      target_base$n_mask <- metrics$n_mask
      target_base$component_density <- if (is.finite(metrics$n_components) && is.finite(metrics$n_mask) && metrics$n_mask > 0) metrics$n_components / metrics$n_mask else NA_real_
      target_base$consensus_score <- 1
      target_base$priority_score <- priority_score
      candidate_targets[[key]] <- target_base
      catalog_rows[[length(catalog_rows) + 1L]] <- data.frame(
        target_key = key,
        field = field,
        label = field,
        target_group = sprintf("%s::%s", field, field),
        target_type = target_base$type,
        route = "direct",
        support_fraction = metrics$support_fraction,
        largest_component_fraction = metrics$largest_component_fraction,
        boundary_fraction = metrics$boundary_fraction,
        n_components = metrics$n_components,
        n_mask = metrics$n_mask,
        component_density = if (is.finite(metrics$n_components) && is.finite(metrics$n_mask) && metrics$n_mask > 0) metrics$n_components / metrics$n_mask else NA_real_,
        consensus_score = 1,
        priority_score = priority_score,
        strict_admitted = identical(target_base$type, "binary_mask") && is.finite(priority_score),
        stringsAsFactors = FALSE
      )
      next
    }

    reg_out <- ra_regularize_categorical_pdata_target(target_base)
    target_base <- reg_out$target
    reg_diag <- reg_out$diagnostics
    if (!is.null(reg_diag)) {
      field_diag_rows[[length(field_diag_rows) + 1L]] <- data.frame(
        field = field,
        changed_fraction = reg_diag$changed_fraction %||% NA_real_,
        raw_boundary_fraction = reg_diag$raw$boundary_fraction %||% NA_real_,
        regularized_boundary_fraction = reg_diag$regularized$boundary_fraction %||% NA_real_,
        raw_component_density = reg_diag$raw$component_density %||% NA_real_,
        regularized_component_density = reg_diag$regularized$component_density %||% NA_real_,
        raw_component_count = reg_diag$raw$component_count %||% NA_integer_,
        regularized_component_count = reg_diag$regularized$component_count %||% NA_integer_,
        raw_largest_label_fraction = reg_diag$raw$largest_label_fraction %||% NA_real_,
        regularized_largest_label_fraction = reg_diag$regularized$largest_label_fraction %||% NA_real_,
        raw_n_labels = reg_diag$raw$n_labels %||% NA_integer_,
        regularized_n_labels = reg_diag$regularized$n_labels %||% NA_integer_,
        stringsAsFactors = FALSE
      )
    }

    categorical_fields <- c(categorical_fields, field)
    labs <- target_base$labels
    valid <- target_base$valid
    label_counts <- sort(table(labs[valid]), decreasing = TRUE)
    keep_labels <- names(label_counts)
    for (lab in keep_labels) {
      mask <- (labs == lab) & valid
      metrics <- ra_mask_shape_metrics(mask, valid_mask = valid)
      if (!is.finite(metrics$support_fraction) || metrics$support_fraction < 0.05 || metrics$support_fraction > 0.95) next
      broad_contiguous <- metrics$support_fraction >= 0.15 &&
        (metrics$largest_component_fraction >= 0.70 || metrics$boundary_fraction <= 0.30)
      if (isTRUE(broad_contiguous)) {
        contour_target <- ra_build_contour_target(mask, band_px = 1L, family = "pdata", field = field, label = lab)
        if (!is.null(contour_target)) {
          key <- sprintf("strict::%s::%s::contour", field, make.names(lab))
          add_target(key, contour_target, field, lab, "contour_line", "auto_boundary", metrics, consensus_score = 1)
        }
      } else {
        bin_target <- ra_make_binary_mask_target(mask, field = field, label = lab, family = "pdata", band_px = boundary_band_px)
        if (!is.null(bin_target)) {
          key <- sprintf("strict::%s::%s::mask", field, make.names(lab))
          add_target(key, bin_target, field, lab, "binary_mask", "direct_mask", metrics, consensus_score = 1)
        }
        contour_target <- ra_build_contour_target(mask, band_px = 1L, family = "pdata", field = field, label = lab)
        if (!is.null(contour_target)) {
          key <- sprintf("strict::%s::%s::contour", field, make.names(lab))
          add_target(key, contour_target, field, lab, "contour_line", "companion_boundary", metrics, consensus_score = 1)
        }
      }
    }
  }

  categorical_fields <- unique(categorical_fields)
  if (length(categorical_fields) >= 2L) {
    boundary_list <- lapply(categorical_fields, function(field) {
      target_base <- ra_build_pdata_target(msi_geom, field, gaussian_smooth = FALSE, smooth_sigma = 1)
      if (identical(target_base$type, "categorical_label_map")) {
        target_base <- ra_regularize_categorical_pdata_target(target_base)$target
      }
      ra_labels_to_boundary_mask(target_base$labels, valid_mask = target_base$valid)
    })
    boundary_arr <- simplify2array(boundary_list)
    if (length(dim(boundary_arr)) == 3L) {
      consensus_frac <- apply(boundary_arr, c(1, 2), mean, na.rm = TRUE)
      consensus_mask <- consensus_frac >= 0.5
      metrics <- ra_mask_shape_metrics(consensus_mask)
      if (is.finite(metrics$support_fraction) && metrics$support_fraction >= 0.02 && metrics$support_fraction <= 0.80) {
        target <- ra_build_contour_target(consensus_mask, band_px = 1L, family = "pdata", field = "cluster_boundary_consensus", label = "boundary_consensus")
        if (!is.null(target)) {
          add_target(
            "strict::cluster_boundary_consensus",
            target,
            paste(categorical_fields, collapse = "|"),
            "boundary_consensus",
            "contour_line",
            "consensus_boundary",
            metrics,
            consensus_score = suppressWarnings(mean(consensus_frac[consensus_mask], na.rm = TRUE))
          )
        }
      }
    }
  }

  catalog <- if (length(catalog_rows) > 0L) do.call(rbind, catalog_rows) else data.frame()
  field_diagnostics <- if (length(field_diag_rows) > 0L) do.call(rbind, field_diag_rows) else data.frame()
  admitted_keys <- character(0)
  if (nrow(catalog) > 0L) {
    eligible <- is.finite(catalog$priority_score)
    ord <- order(catalog$priority_score, decreasing = TRUE, na.last = TRUE)
    field_counts <- setNames(integer(0), character(0))
    seen_groups <- character(0)
    admitted_idx <- integer(0)
    for (ii in ord) {
      if (!eligible[ii]) next
      field_name <- as.character(catalog$field[ii])
      group_name <- as.character(catalog$target_group[ii])
      field_count <- if (field_name %in% names(field_counts)) field_counts[[field_name]] else 0L
      if (group_name %in% seen_groups || field_count >= max_per_field) next
      admitted_idx <- c(admitted_idx, ii)
      seen_groups <- c(seen_groups, group_name)
      field_counts[[field_name]] <- field_count + 1L
      if (length(admitted_idx) >= max_targets) break
    }
    if (length(admitted_idx) > 0L) {
      catalog$strict_admitted[admitted_idx] <- TRUE
      admitted_keys <- catalog$target_key[admitted_idx]
    }
  }
  targets <- candidate_targets[intersect(admitted_keys, names(candidate_targets))]
  list(targets = targets, catalog = catalog, field_diagnostics = field_diagnostics)
}

ra_score_signal_against_target <- function(signal, valid_mask, target, relation = "either") {
  relation <- tolower(trimws(as.character(relation)[1] %||% "either"))
  if (!relation %in% c("direct", "inverse", "either")) relation <- "either"
  if (identical(target$type, "categorical_label_map")) {
    keep <- valid_mask & target$valid
    raw <- ra_categorical_signal_association(signal[keep], target$labels[keep], min_group_n = 5L)
    return(list(score = raw, raw_score = raw, corr_score = NA_real_, nmi_score = NA_real_, edge_score = NA_real_, overlap = sum(keep, na.rm = TRUE), n_edge = 0L))
  }

  if (identical(target$type, "boundary_band")) {
    hp_structure <- ra_build_structure_signal_map(signal, valid_mask = valid_mask, positive_only = TRUE)
    sig_use <- if (!is.null(hp_structure)) hp_structure$signal else signal
    valid_use <- if (!is.null(hp_structure)) hp_structure$valid else valid_mask
    ring_keep <- target$boundary_ring & valid_use
    ctx_keep <- target$context_mask & valid_use
    ring_mean <- if (sum(ring_keep, na.rm = TRUE) >= 5L) suppressWarnings(mean(sig_use[ring_keep], na.rm = TRUE)) else NA_real_
    ctx_mean <- if (sum(ctx_keep, na.rm = TRUE) >= 5L) suppressWarnings(mean(sig_use[ctx_keep], na.rm = TRUE)) else NA_real_
    ring_enrich <- if (is.finite(ring_mean) && is.finite(ctx_mean)) ring_mean - ctx_mean else NA_real_

    inner_keep <- target$inner_mask & valid_mask
    outer_keep <- target$outer_mask & valid_mask
    inner_mean <- if (sum(inner_keep, na.rm = TRUE) >= 5L) suppressWarnings(mean(signal[inner_keep], na.rm = TRUE)) else NA_real_
    outer_mean <- if (sum(outer_keep, na.rm = TRUE) >= 5L) suppressWarnings(mean(signal[outer_keep], na.rm = TRUE)) else NA_real_
    shell_contrast <- if (is.finite(inner_mean) && is.finite(outer_mean)) abs(inner_mean - outer_mean) else NA_real_

    edge_info <- if (!is.null(hp_structure)) ra_build_edge_mask(hp_structure$signal, valid_mask = hp_structure$valid) else ra_build_edge_mask(signal, valid_mask = valid_mask)
    sc_edge <- if (!is.null(edge_info)) ra_score_edge_mask_against_distance(edge_info$edge, target$edge_cache) else NA_real_
    sc_ring <- switch(
      relation,
      direct = ring_enrich,
      inverse = if (is.finite(ring_enrich)) -ring_enrich else NA_real_,
      either = if (is.finite(ring_enrich)) abs(ring_enrich) else NA_real_,
      if (is.finite(ring_enrich)) abs(ring_enrich) else NA_real_
    )
    raw <- ra_combine_weighted_scores(c(sc_ring, shell_contrast), c(0.45, 0.55))
    score <- ra_combine_weighted_scores(c(raw, sc_edge), c(0.45, 0.55))
    return(list(
      score = score,
      raw_score = raw,
      corr_score = shell_contrast,
      nmi_score = NA_real_,
      edge_score = sc_edge,
      overlap = sum(ring_keep, na.rm = TRUE),
      n_edge = if (is.null(edge_info)) 0L else edge_info$n_edge
    ))
  }

  if (identical(target$type, "contour_line")) {
    hp_structure <- ra_build_structure_signal_map(signal, valid_mask = valid_mask, positive_only = TRUE)
    sig_use <- if (!is.null(hp_structure)) hp_structure$signal else signal
    valid_use <- if (!is.null(hp_structure)) hp_structure$valid else valid_mask
    line_keep <- target$boundary_line & valid_use
    shell_mask <- target$shell_mask %||% (target$inner_mask | target$outer_mask)
    shell_keep <- shell_mask & valid_use
    line_mean <- if (sum(line_keep, na.rm = TRUE) >= 5L) suppressWarnings(mean(sig_use[line_keep], na.rm = TRUE)) else NA_real_
    shell_mean <- if (sum(shell_keep, na.rm = TRUE) >= 5L) suppressWarnings(mean(sig_use[shell_keep], na.rm = TRUE)) else NA_real_
    line_enrich <- if (is.finite(line_mean) && is.finite(shell_mean)) line_mean - shell_mean else NA_real_

    inner_keep <- target$inner_mask & valid_mask
    outer_keep <- target$outer_mask & valid_mask
    inner_mean <- if (sum(inner_keep, na.rm = TRUE) >= 5L) suppressWarnings(mean(signal[inner_keep], na.rm = TRUE)) else NA_real_
    outer_mean <- if (sum(outer_keep, na.rm = TRUE) >= 5L) suppressWarnings(mean(signal[outer_keep], na.rm = TRUE)) else NA_real_
    shell_contrast_raw <- if (is.finite(inner_mean) && is.finite(outer_mean)) inner_mean - outer_mean else NA_real_
    sc_line <- switch(
      relation,
      direct = line_enrich,
      inverse = if (is.finite(line_enrich)) -line_enrich else NA_real_,
      either = if (is.finite(line_enrich)) abs(line_enrich) else NA_real_,
      if (is.finite(line_enrich)) abs(line_enrich) else NA_real_
    )
    sc_shell <- switch(
      relation,
      direct = shell_contrast_raw,
      inverse = if (is.finite(shell_contrast_raw)) -shell_contrast_raw else NA_real_,
      either = if (is.finite(shell_contrast_raw)) abs(shell_contrast_raw) else NA_real_,
      if (is.finite(shell_contrast_raw)) abs(shell_contrast_raw) else NA_real_
    )
    edge_info <- if (!is.null(hp_structure)) ra_build_edge_mask(hp_structure$signal, valid_mask = hp_structure$valid) else ra_build_edge_mask(signal, valid_mask = valid_mask)
    sc_edge <- if (!is.null(edge_info)) ra_score_edge_mask_against_distance(edge_info$edge, target$edge_cache) else NA_real_
    raw <- ra_combine_weighted_scores(c(sc_shell, sc_line), c(0.65, 0.35))
    score <- ra_combine_weighted_scores(c(sc_edge, sc_shell, sc_line), c(0.55, 0.30, 0.15))
    return(list(
      score = score,
      raw_score = raw,
      corr_score = sc_shell,
      nmi_score = NA_real_,
      edge_score = sc_edge,
      overlap = sum(line_keep, na.rm = TRUE),
      n_edge = if (is.null(edge_info)) 0L else edge_info$n_edge
    ))
  }

  if (identical(target$type, "binary_mask")) {
    sc_corr <- ra_signed_rank_correlation(signal[valid_mask], target$signal[valid_mask])
    sc_nmi <- ra_normalized_mutual_information(signal[valid_mask], target$signal[valid_mask], nbins = 24L)
    sc_mask_raw <- ra_combine_weighted_scores(
      c(
        ra_score_mask_shift(target$mask_cache, signal, dx = 0L, dy = 0L),
        ra_score_mask_centrality_shift(target$mask_cache, signal, dx = 0L, dy = 0L)
      ),
      c(0.35, 0.65)
    )
    hp_structure <- ra_build_structure_signal_map(signal, valid_mask = valid_mask, positive_only = TRUE)
    sc_mask_structure <- if (!is.null(hp_structure)) {
      ra_combine_weighted_scores(
        c(
          ra_score_mask_shift(target$mask_cache, hp_structure$signal, dx = 0L, dy = 0L),
          ra_score_mask_centrality_shift(target$mask_cache, hp_structure$signal, dx = 0L, dy = 0L)
        ),
        c(0.25, 0.75)
      )
    } else {
      NA_real_
    }
    sc_mask <- ra_combine_weighted_scores(c(sc_mask_raw, sc_mask_structure), c(0.85, 0.15))
    sc_intensity <- switch(
      relation,
      direct = sc_mask,
      inverse = if (is.finite(sc_mask)) -sc_mask else NA_real_,
      either = if (is.finite(sc_mask)) abs(sc_mask) else NA_real_,
      if (is.finite(sc_mask)) abs(sc_mask) else NA_real_
    )
    edge_info <- if (!is.null(hp_structure)) ra_build_edge_mask(hp_structure$signal, valid_mask = hp_structure$valid) else ra_build_edge_mask(signal, valid_mask = valid_mask)
    sc_edge <- if (!is.null(edge_info)) ra_score_edge_mask_against_distance(edge_info$edge, target$edge_cache) else NA_real_
    score <- ra_combine_weighted_scores(c(sc_intensity, sc_edge), c(0.85, 0.15))
    return(list(score = score, raw_score = sc_intensity, corr_score = sc_corr, nmi_score = sc_nmi, edge_score = sc_edge, overlap = sum(valid_mask, na.rm = TRUE), n_edge = if (is.null(edge_info)) 0L else edge_info$n_edge))
  }

  keep <- valid_mask & target$valid
  sc_corr <- ra_signed_rank_correlation(signal[keep], target$signal[keep])
  sc_nmi <- ra_normalized_mutual_information(signal[keep], target$signal[keep], nbins = 24L)
  sc_relation <- switch(
    relation,
    direct = sc_corr,
    inverse = if (is.finite(sc_corr)) -sc_corr else NA_real_,
    either = if (is.finite(sc_corr)) abs(sc_corr) else NA_real_,
    if (is.finite(sc_corr)) abs(sc_corr) else NA_real_
  )
  hp_structure <- ra_build_structure_signal_map(signal, valid_mask = valid_mask, positive_only = TRUE)
  sc_structure_corr <- NA_real_
  sc_structure_nmi <- NA_real_
  if (!is.null(hp_structure) && !is.null(target$structure) && is.matrix(target$structure$signal)) {
    keep_structure <- hp_structure$valid & target$structure$valid
    if (sum(keep_structure, na.rm = TRUE) >= 25L) {
      sc_structure_corr <- ra_signed_rank_correlation(hp_structure$signal[keep_structure], target$structure$signal[keep_structure])
      sc_structure_nmi <- ra_normalized_mutual_information(hp_structure$signal[keep_structure], target$structure$signal[keep_structure], nbins = 24L)
    }
  }
  sc_structure <- ra_combine_weighted_scores(c(sc_structure_corr, sc_structure_nmi), c(0.65, 0.35))
  sc_intensity <- ra_combine_weighted_scores(c(sc_relation, sc_structure), c(0.30, 0.70))
  edge_info <- if (!is.null(hp_structure)) ra_build_edge_mask(hp_structure$signal, valid_mask = hp_structure$valid) else ra_build_edge_mask(signal, valid_mask = valid_mask)
  sc_edge <- if (!is.null(edge_info)) ra_score_edge_mask_against_distance(edge_info$edge, target$edge_cache) else NA_real_
  score <- ra_combine_weighted_scores(c(sc_intensity, sc_edge), c(0.45, 0.55))
  list(score = score, raw_score = sc_intensity, corr_score = sc_corr, nmi_score = sc_nmi, edge_score = sc_edge, overlap = sum(keep, na.rm = TRUE), n_edge = if (is.null(edge_info)) 0L else edge_info$n_edge)
}

ra_build_numeric_target <- function(signal, valid = NULL, family, field, type = "numeric") {
  if (is.null(valid)) valid <- is.finite(signal)
  list(
    family = family,
    field = field,
    signal = signal,
    valid = valid,
    type = type,
    structure = ra_build_structure_signal_map(signal, valid, positive_only = TRUE),
    edge_cache = ra_build_edge_distance_cache(signal)
  )
}

ra_build_selected_ion_target <- function(obj, msi_geom, mz_value = NULL, mz_index = NULL, intensity_transform = "asinh", gaussian_smooth = FALSE, smooth_sigma = 1) {
  mz_axis <- suppressWarnings(as.numeric(Cardinal::mz(obj)))
  if (is.null(mz_index) || !is.finite(mz_index)) {
    if (is.null(mz_value) || !is.finite(mz_value)) stop("Need mz_value or mz_index for selected ion target.")
    mz_index <- which.min(abs(mz_axis - mz_value))
  }
  sp <- as.matrix(Cardinal::spectra(obj)[mz_index, , drop = FALSE])
  if (nrow(sp) != 1L && ncol(sp) == 1L) sp <- t(sp)
  sig_payload <- ra_build_numeric_signal_matrix(
    vals = as.numeric(sp[1, ]),
    msi_geom = msi_geom,
    intensity_transform = intensity_transform,
    apply_transform = TRUE,
    gaussian_smooth = gaussian_smooth,
    smooth_sigma = smooth_sigma,
    enhance = FALSE
  )
  out <- ra_build_numeric_target(sig_payload$signal, sig_payload$valid, family = "msi", field = sprintf("m/z %.5f", mz_axis[mz_index]), type = "mz")
  out$mz_index <- mz_index
  out$mz_value <- mz_axis[mz_index]
  out
}

ra_resolve_polygon_axis_mode <- function(poly_sf, mode = "auto", hist_img = NULL) {
  mode_chr <- tolower(trimws(as.character(mode)[1] %||% "auto"))
  if (mode_chr %in% c("xy", "yx")) return(mode_chr)
  if (is.null(hist_img) || is.null(poly_sf) || nrow(poly_sf) == 0L) return("xy")
  info <- try(magick::image_info(hist_img)[1, ], silent = TRUE)
  if (inherits(info, "try-error")) return("xy")
  img_w <- as.numeric(info$width)
  img_h <- as.numeric(info$height)
  bb <- try(sf::st_bbox(poly_sf), silent = TRUE)
  if (inherits(bb, "try-error")) return("xy")
  w <- as.numeric(bb$xmax - bb$xmin)
  h <- as.numeric(bb$ymax - bb$ymin)
  if (!all(is.finite(c(img_w, img_h, w, h))) || img_w <= 0 || img_h <= 0 || w <= 0 || h <= 0) return("xy")
  ar_img <- img_w / img_h
  d_xy <- abs(log((w / h) / ar_img))
  d_yx <- abs(log((h / w) / ar_img))
  if (is.finite(d_yx) && is.finite(d_xy) && d_yx < d_xy) "yx" else "xy"
}

ra_polygon_source_to_image_sf <- function(poly_sf, overlay_frame, image_width, image_height, axis_mode = "xy", base_flip_y = TRUE) {
  geom_type <- as.character(sf::st_geometry_type(poly_sf))
  poly <- poly_sf[geom_type %in% c("POLYGON", "MULTIPOLYGON"), , drop = FALSE]
  if (nrow(poly) == 0L) stop("No polygon geometries found.")
  swap_xy <- identical(axis_mode, "yx")
  source_origin_x <- overlay_frame$origin_x
  source_origin_y <- overlay_frame$origin_y
  downsample_factor <- overlay_frame$overlay_downsample_factor
  if (!is.finite(downsample_factor) || downsample_factor <= 0) downsample_factor <- 1

  rebuild_one <- function(one_row) {
    cc <- sf::st_coordinates(one_row)
    x_src <- cc[, if ("X" %in% colnames(cc)) "X" else "x"]
    y_src <- cc[, if ("Y" %in% colnames(cc)) "Y" else "y"]
    if (isTRUE(swap_xy)) {
      tmp <- x_src
      x_src <- y_src
      y_src <- tmp
    }
    x_frame <- (x_src - source_origin_x) / downsample_factor
    y_frame <- (y_src - source_origin_y) / downsample_factor
    x1 <- x_frame
    y1 <- if (isTRUE(base_flip_y)) (image_height - y_frame) else y_frame
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
      if (length(polys) == 1L) polys[[1]] else sf::st_multipolygon(polys)
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

  geoms <- lapply(seq_len(nrow(poly)), function(i) rebuild_one(poly[i, , drop = FALSE]))
  poly_out <- poly
  sf::st_geometry(poly_out) <- sf::st_sfc(geoms)
  poly_out
}

ra_render_polygon_mask_image <- function(poly_sf, overlay_frame, hist_img, axis_mode = "xy") {
  info <- magick::image_info(hist_img)[1, ]
  width <- as.integer(info$width)
  height <- as.integer(info$height)
  base_flip_y <- ra_overlay_uses_qupath_source_frame(overlay_frame)
  poly_img <- ra_polygon_source_to_image_sf(poly_sf, overlay_frame, image_width = width, image_height = height, axis_mode = axis_mode, base_flip_y = base_flip_y)
  tf <- tempfile(fileext = ".png")
  grDevices::png(tf, width = width, height = height, bg = "transparent")
  op <- graphics::par(mar = c(0, 0, 0, 0))
  dev_id <- grDevices::dev.cur()
  on.exit({
    graphics::par(op)
    if (dev_id %in% grDevices::dev.list()) grDevices::dev.off(dev_id)
  }, add = TRUE)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, width), ylim = c(0, height), xaxs = "i", yaxs = "i", asp = 1)
  graphics::rect(0, 0, width, height, col = "#00000000", border = NA)
  graphics::plot(sf::st_geometry(poly_img), add = TRUE, col = "white", border = NA, axes = FALSE, reset = FALSE)
  graphics::par(op)
  grDevices::dev.off(dev_id)
  img <- magick::image_read(tf)
  unlink(tf)
  img
}

ra_polygon_select_top_variable_features <- function(X, max_features = 300L) {
  max_features <- as.integer(max(20L, max_features))
  rr <- rowMeans(X, na.rm = TRUE)
  rr[!is.finite(rr)] <- 0
  vv <- apply(X, 1L, function(v) {
    v <- as.numeric(v)
    v <- v[is.finite(v)]
    if (length(v) < 2L) return(0)
    stats::var(v)
  })
  vv[!is.finite(vv)] <- 0
  detect <- rowMeans(is.finite(X) & X > 0, na.rm = TRUE)
  detect[!is.finite(detect)] <- 0
  score <- rank(rr, ties.method = "average") +
    rank(vv, ties.method = "average") +
    0.5 * rank(detect, ties.method = "average")
  ord <- order(score, decreasing = TRUE, na.last = NA)
  ord[seq_len(min(length(ord), max_features))]
}

ra_polygon_stat_ttest_count_sig <- function(X, idx_inside, idx_outside, alpha = 0.1, use_adjusted = TRUE, use_abs_lfc = TRUE) {
  if (length(idx_inside) < 2L || length(idx_outside) < 2L) {
    return(list(n_sig = 0L, score = -Inf, med_lfc = NA_real_, mean_logp = NA_real_, mean_abs_t_top = NA_real_))
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
  score <- as.numeric(n_sig) + 0.25 * med_lfc + 0.1 * mean_logp + 0.05 * mean_abs_t_top
  list(n_sig = n_sig, score = score, med_lfc = med_lfc, mean_logp = mean_logp, mean_abs_t_top = mean_abs_t_top)
}

ra_polygon_stat_feature_table <- function(X, idx_inside, idx_outside, mz_values = NULL, use_adjusted = TRUE) {
  if (length(idx_inside) < 2L || length(idx_outside) < 2L) return(data.frame())
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
  if (is.null(mz_values) || length(mz_values) != nrow(X)) mz_values <- rep(NA_real_, nrow(X))
  out <- data.frame(
    feature_index = seq_len(nrow(X)),
    mz = suppressWarnings(as.numeric(mz_values)),
    abs_log2FC = as.numeric(abs_lfc),
    t_stat = as.numeric(tstat),
    score_p = as.numeric(score_p),
    neglog10_score_p = as.numeric(neglogp),
    stringsAsFactors = FALSE
  )
  out <- out[is.finite(out$score_p) & is.finite(out$abs_log2FC), , drop = FALSE]
  if (nrow(out) > 0L) {
    out$fit_score <- out$neglog10_score_p + 0.25 * out$abs_log2FC
    out <- out[order(out$fit_score, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  }
  out
}

ra_bbox_outside_mask <- function(mask_inside, pad_px = 25L) {
  b <- suppressWarnings(as.integer(pad_px))
  if (!is.finite(b) || b < 0L) b <- 25L
  b <- as.integer(min(500L, max(0L, b)))
  idx <- which(mask_inside)
  if (length(idx) == 0L) return(matrix(FALSE, nrow = nrow(mask_inside), ncol = ncol(mask_inside)))
  rc <- arrayInd(idx, .dim = dim(mask_inside))
  rmin <- max(1L, min(rc[, 1]) - b)
  rmax <- min(nrow(mask_inside), max(rc[, 1]) + b)
  cmin <- max(1L, min(rc[, 2]) - b)
  cmax <- min(ncol(mask_inside), max(rc[, 2]) + b)
  out <- matrix(FALSE, nrow = nrow(mask_inside), ncol = ncol(mask_inside))
  out[rmin:rmax, cmin:cmax] <- TRUE
  out & !mask_inside
}

ra_prepare_polygon_stat_context <- function(obj, msi_geom, max_features = 200L, intensity_transform = "asinh") {
  sp <- ra_extract_msi_intensity_matrix(obj)
  if (!identical(intensity_transform, "none")) sp <- ra_transform_intensity(sp, intensity_transform)
  sp[!is.finite(sp)] <- 0
  feat_idx <- ra_polygon_select_top_variable_features(sp, max_features = max_features)
  valid_pix <- matrix(FALSE, nrow = msi_geom$ny, ncol = msi_geom$nx)
  valid_pix[cbind(msi_geom$row_idx, msi_geom$x_norm)] <- TRUE
  list(
    X = sp[feat_idx, , drop = FALSE],
    feature_index = feat_idx,
    mz_values = suppressWarnings(as.numeric(Cardinal::mz(obj)))[feat_idx],
    valid_pix = valid_pix,
    msi_geom = msi_geom
  )
}

ra_score_polygon_mask_discrimination <- function(context, mask_inside, outside_mode = "bbox", bbox_pad = 25L, buffer_px = 3L, min_pixels = 10L) {
  mask_i <- matrix(as.logical(mask_inside), nrow = nrow(mask_inside), ncol = ncol(mask_inside))
  outside_i <- switch(
    outside_mode,
    local = {
      grown <- mask_i
      for (ii in seq_len(max(1L, buffer_px))) grown <- ra_dilate8(grown)
      grown & !mask_i
    },
    global = !mask_i,
    bbox = ra_bbox_outside_mask(mask_i, pad_px = bbox_pad),
    ra_bbox_outside_mask(mask_i, pad_px = bbox_pad)
  )
  outside_i <- outside_i & !mask_i & context$valid_pix
  inside_vec <- mask_i[cbind(context$msi_geom$row_idx, context$msi_geom$x_norm)]
  outside_vec <- outside_i[cbind(context$msi_geom$row_idx, context$msi_geom$x_norm)]
  n_in <- as.integer(sum(inside_vec, na.rm = TRUE))
  n_out <- as.integer(sum(outside_vec, na.rm = TRUE))
  if (n_in < min_pixels || n_out < min_pixels) {
    return(list(
      n_inside = n_in, n_outside = n_out,
      n_sig = NA_integer_, score = NA_real_, med_lfc = NA_real_,
      mean_logp = NA_real_, mean_abs_t_top = NA_real_
    ))
  }
  out <- ra_polygon_stat_ttest_count_sig(context$X, which(inside_vec), which(outside_vec))
  out$n_inside <- n_in
  out$n_outside <- n_out
  out
}

ra_polygon_feature_ranking <- function(context, mask_inside, outside_mode = "bbox", bbox_pad = 25L, buffer_px = 3L, min_pixels = 10L) {
  mask_i <- matrix(as.logical(mask_inside), nrow = nrow(mask_inside), ncol = ncol(mask_inside))
  outside_i <- switch(
    outside_mode,
    local = {
      grown <- mask_i
      for (ii in seq_len(max(1L, buffer_px))) grown <- ra_dilate8(grown)
      grown & !mask_i
    },
    global = !mask_i,
    bbox = ra_bbox_outside_mask(mask_i, pad_px = bbox_pad),
    ra_bbox_outside_mask(mask_i, pad_px = bbox_pad)
  )
  outside_i <- outside_i & !mask_i & context$valid_pix
  inside_vec <- mask_i[cbind(context$msi_geom$row_idx, context$msi_geom$x_norm)]
  outside_vec <- outside_i[cbind(context$msi_geom$row_idx, context$msi_geom$x_norm)]
  if (sum(inside_vec, na.rm = TRUE) < min_pixels || sum(outside_vec, na.rm = TRUE) < min_pixels) return(data.frame())
  ra_polygon_stat_feature_table(context$X, which(inside_vec), which(outside_vec), mz_values = context$mz_values)
}

ra_select_strict_ion_indices <- function(obj, msi_geom, polygon_feature_table = NULL, pdata_targets = list(), n_candidates = 5L, candidate_pool = 64L, intensity_transform = "asinh") {
  n_candidates <- as.integer(max(1L, n_candidates))
  candidate_pool <- as.integer(max(candidate_pool, 4L * n_candidates, 16L))
  if (!is.null(polygon_feature_table) && nrow(polygon_feature_table) > 0L) {
    keep_n <- max(n_candidates, ceiling(0.20 * nrow(polygon_feature_table)))
    cand <- unique(as.integer(polygon_feature_table$feature_index[seq_len(min(nrow(polygon_feature_table), keep_n))]))
  } else {
    cand <- ra_suggest_histology_fusion_indices(obj, msi_geom, n_features = n_candidates, candidate_pool = candidate_pool, intensity_transform = intensity_transform)
  }
  if (length(cand) == 0L) return(integer(0))
  sp <- try(as.matrix(Cardinal::spectra(obj)[cand, , drop = FALSE]), silent = TRUE)
  if (inherits(sp, "try-error") || is.null(sp)) return(integer(0))
  if (nrow(sp) != length(cand) && ncol(sp) == length(cand)) sp <- t(sp)
  storage.mode(sp) <- "double"
  spatial_score <- rep(0, length(cand))
  pdata_score <- rep(0, length(cand))
  coverage <- rep(0, length(cand))
  for (ii in seq_along(cand)) {
    sig_payload <- ra_build_numeric_signal_matrix(sp[ii, ], msi_geom, intensity_transform = intensity_transform, apply_transform = TRUE, gaussian_smooth = FALSE, enhance = FALSE)
    spatial_score[ii] <- ra_combine_weighted_scores(
      c(
        ra_score_spatial_coherence_map(sig_payload$signal, sig_payload$valid, sigma = 2),
        ra_score_spatial_signal_map(sig_payload$signal, sig_payload$valid)
      ),
      c(0.75, 0.25)
    )
    coverage[ii] <- suppressWarnings(mean(sig_payload$signal[sig_payload$valid] >= 0.20, na.rm = TRUE))
    if (length(pdata_targets) > 0L) {
      pdata_hits <- vapply(pdata_targets, function(target) {
        sc <- ra_score_signal_against_target(sig_payload$signal, sig_payload$valid, target, relation = "either")
        if (is.finite(sc$score)) sc$score else NA_real_
      }, numeric(1))
      if (any(is.finite(pdata_hits))) pdata_score[ii] <- max(pdata_hits, na.rm = TRUE)
    }
  }
  norm01 <- function(v) {
    v <- as.numeric(v)
    v[!is.finite(v)] <- 0
    rr <- range(v, na.rm = TRUE)
    if (!all(is.finite(rr)) || diff(rr) <= 1e-12) return(rep(0, length(v)))
    (v - rr[1]) / diff(rr)
  }
  coverage[!is.finite(coverage)] <- 0
  keep <- coverage >= 0.01
  if (!any(keep)) keep <- rep(TRUE, length(cand))
  final_score <- (0.55 * norm01(pdata_score) + 0.25 * norm01(spatial_score) + 0.20 * norm01(coverage)) * sqrt(pmax(coverage, 0))
  final_score[!keep] <- -Inf
  ord <- order(final_score, decreasing = TRUE, na.last = NA)
  unique(cand[ord][seq_len(min(sum(is.finite(final_score[ord])), n_candidates))])
}

ra_rank_normalize_desc <- function(x) {
  x <- as.numeric(x)
  rank_raw <- rank(-x, ties.method = "average", na.last = "keep")
  out <- rep(NA_real_, length(rank_raw))
  keep <- is.finite(rank_raw)
  if (!any(keep)) return(out)
  if (sum(keep) == 1L) {
    out[keep] <- 1
  } else {
    out[keep] <- 1 - ((rank_raw[keep] - 1) / (sum(keep) - 1))
  }
  out
}

ra_confidence_label <- function(weighted_scores, agreement_counts, family_counts = NULL) {
  scores <- suppressWarnings(as.numeric(weighted_scores))
  scores <- scores[is.finite(scores)]
  if (length(scores) < 2L) return("confident")
  ord <- sort(scores, decreasing = TRUE)
  gap <- ord[1] - ord[2]
  agree_top <- if (length(agreement_counts) > 0L) suppressWarnings(max(as.numeric(agreement_counts), na.rm = TRUE)) else NA_real_
  family_top <- if (length(family_counts) > 0L) suppressWarnings(max(as.numeric(family_counts), na.rm = TRUE)) else NA_real_
  if (is.finite(family_top) && family_top <= 1) {
    if (is.finite(gap) && gap >= 0.03) return("confident")
    if (!is.finite(gap) || gap < 0.01) return("conflicted")
    return("ambiguous")
  }
  if (is.finite(gap) && gap >= 0.10 && is.finite(agree_top) && agree_top >= 2) return("confident")
  if (!is.finite(gap) || gap < 0.03 || (!is.na(agree_top) && agree_top <= 1)) return("conflicted")
  "ambiguous"
}
