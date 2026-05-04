normalize_slide_key <- function(x) {
  txt <- trimws(as.character(x)[1])
  if (is.na(txt) || !nzchar(txt)) return(NA_character_)
  key <- tolower(gsub("[^A-Za-z0-9]+", "_", txt))
  key <- gsub("^_+|_+$", "", key)
  if (!nzchar(key)) return(NA_character_)
  key
}

parse_histology_crop_from_name <- function(path_or_name) {
  txt <- trimws(as.character(path_or_name)[1])
  if (is.na(txt) || !nzchar(txt)) return(NULL)
  base <- basename(txt)
  m <- regexec(
    "(?i)\\((?:[^)]*?)x\\s*=\\s*([-0-9.]+)\\s*,\\s*y\\s*=\\s*([-0-9.]+)\\s*,\\s*w\\s*=\\s*([-0-9.]+)\\s*,\\s*h\\s*=\\s*([-0-9.]+)(?:[^)]*?)\\)",
    base,
    perl = TRUE
  )
  hit <- regmatches(base, m)[[1]]
  if (length(hit) != 5L) return(NULL)
  vals <- suppressWarnings(as.numeric(hit[2:5]))
  if (!all(is.finite(vals))) return(NULL)
  list(
    x = vals[1],
    y = vals[2],
    width = vals[3],
    height = vals[4],
    source = "filename_crop"
  )
}

polygon_geometry_bbox <- function(poly_sf) {
  if (is.null(poly_sf) || nrow(poly_sf) == 0L) return(NULL)
  bb <- try(sf::st_bbox(poly_sf), silent = TRUE)
  if (inherits(bb, "try-error")) return(NULL)
  out <- c(
    xmin = suppressWarnings(as.numeric(bb$xmin)),
    ymin = suppressWarnings(as.numeric(bb$ymin)),
    xmax = suppressWarnings(as.numeric(bb$xmax)),
    ymax = suppressWarnings(as.numeric(bb$ymax))
  )
  if (!all(is.finite(out))) return(NULL)
  out
}

detect_qupath_roi_anchor <- function(poly_sf) {
  if (is.null(poly_sf) || nrow(poly_sf) == 0L) return(NULL)
  props <- try(as.data.frame(poly_sf), silent = TRUE)
  if (inherits(props, "try-error") || is.null(props)) return(NULL)
  sf_col <- attr(poly_sf, "sf_column")
  if (!is.null(sf_col) && sf_col %in% colnames(props)) {
    props[[sf_col]] <- NULL
  }
  object_type <- rep("", nrow(poly_sf))
  if ("objectType" %in% colnames(props)) {
    object_type <- tolower(trimws(as.character(props$objectType)))
  }
  ann_idx <- which(object_type %in% c("annotation", "roi", "region"))
  if (length(ann_idx) == 0L) return(NULL)

  ann_poly <- poly_sf[ann_idx, , drop = FALSE]
  ann_poly_valid <- try(suppressWarnings(sf::st_make_valid(ann_poly)), silent = TRUE)
  if (!inherits(ann_poly_valid, "try-error") && !is.null(ann_poly_valid) && nrow(ann_poly_valid) == nrow(ann_poly)) {
    ann_poly <- ann_poly_valid
  }

  area <- try(suppressWarnings(as.numeric(sf::st_area(ann_poly))), silent = TRUE)
  if (inherits(area, "try-error") || length(area) != length(ann_idx)) {
    bb_area <- vapply(seq_len(nrow(ann_poly)), function(i) {
      bb <- try(sf::st_bbox(ann_poly[i, , drop = FALSE]), silent = TRUE)
      if (inherits(bb, "try-error")) return(NA_real_)
      w <- suppressWarnings(as.numeric(bb$xmax - bb$xmin))
      h <- suppressWarnings(as.numeric(bb$ymax - bb$ymin))
      if (!is.finite(w) || !is.finite(h) || w <= 0 || h <= 0) return(NA_real_)
      w * h
    }, numeric(1))
    area <- bb_area
  }
  area[!is.finite(area)] <- -Inf
  if (all(!is.finite(area))) return(NULL)
  ord <- order(area, decreasing = TRUE)
  best_idx <- ann_idx[ord[1]]
  list(
    polygon = poly_sf[best_idx, , drop = FALSE],
    index = best_idx,
    source = "embedded_annotation",
    area = area[ord[1]]
  )
}

overlay_frame_from_inputs <- function(
  hist_img = NULL,
  roi_poly = NULL,
  histology_resample_factor = 1,
  polygon_file_name = NULL,
  histology_file_name = NULL
) {
  downsample <- suppressWarnings(as.numeric(histology_resample_factor)[1])
  if (!is.finite(downsample) || downsample <= 0) downsample <- 1

  roi_bbox <- polygon_geometry_bbox(roi_poly)

  info <- NULL
  if (!is.null(hist_img)) {
    info <- try(magick::image_info(hist_img)[1, ], silent = TRUE)
    if (inherits(info, "try-error")) info <- NULL
  }

  slide_key <- normalize_slide_key(polygon_file_name %||% histology_file_name)

  out <- list(
    slide_key = slide_key,
    origin_x = NA_real_,
    origin_y = NA_real_,
    source_width = NA_real_,
    source_height = NA_real_,
    image_origin_x = NA_real_,
    image_origin_y = NA_real_,
    image_source_width = NA_real_,
    image_source_height = NA_real_,
    image_width = if (!is.null(info)) suppressWarnings(as.numeric(info$width)) else NA_real_,
    image_height = if (!is.null(info)) suppressWarnings(as.numeric(info$height)) else NA_real_,
    overlay_downsample_factor = downsample,
    source_frame_type = "legacy",
    frame_status = "legacy_best_effort",
    roi_anchor_present = FALSE,
    roi_anchor_bbox = roi_bbox,
    roi_source = "none"
  )

  crop_info <- parse_histology_crop_from_name(histology_file_name)
  if (!is.null(crop_info)) {
    out$image_origin_x <- crop_info$x
    out$image_origin_y <- crop_info$y
    out$image_source_width <- crop_info$width * downsample
    out$image_source_height <- crop_info$height * downsample
  }

  if (!is.null(roi_bbox)) {
    out$roi_anchor_present <- TRUE
    out$origin_x <- roi_bbox[["xmin"]]
    out$origin_y <- roi_bbox[["ymin"]]
    out$source_width <- roi_bbox[["xmax"]] - roi_bbox[["xmin"]]
    out$source_height <- roi_bbox[["ymax"]] - roi_bbox[["ymin"]]
    out$roi_source <- "roi"
    if (identical(out$frame_status, "legacy_best_effort")) out$frame_status <- "compatible"
    if (identical(out$source_frame_type, "legacy")) out$source_frame_type <- "roi_relative"
  }

  if (is.finite(out$image_width) && out$image_width > 0 &&
      is.finite(out$image_height) && out$image_height > 0 &&
      (!is.finite(out$source_width) || out$source_width <= 0 ||
       !is.finite(out$source_height) || out$source_height <= 0)) {
    out$source_width <- out$image_width * downsample
    out$source_height <- out$image_height * downsample
    if (!is.finite(out$origin_x)) out$origin_x <- 0
    if (!is.finite(out$origin_y)) out$origin_y <- 0
    if (identical(out$source_frame_type, "legacy")) out$source_frame_type <- "image_bbox"
  }

  if (!is.finite(out$image_source_width) || out$image_source_width <= 0) {
    out$image_source_width <- out$image_width * downsample
  }
  if (!is.finite(out$image_source_height) || out$image_source_height <= 0) {
    out$image_source_height <- out$image_height * downsample
  }
  if (!is.finite(out$image_origin_x)) out$image_origin_x <- out$origin_x %||% 0
  if (!is.finite(out$image_origin_y)) out$image_origin_y <- out$origin_y %||% 0

  out
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x)) || !nzchar(trimws(as.character(x)[1]))) y else x
}

source_xy_to_display_xy <- function(x, y, coord_frame, flip_y = TRUE) {
  if (is.null(coord_frame)) return(list(x = x, y = y))
  if (isTRUE(flip_y)) {
    list(
      x = as.numeric(x) - coord_frame$xmin + 1,
      y = coord_frame$ymax - as.numeric(y) + 1
    )
  } else {
    list(
      x = as.numeric(x) - coord_frame$xmin + 1,
      y = as.numeric(y) - coord_frame$ymin + 1
    )
  }
}

display_xy_to_source_xy <- function(x, y, coord_frame, flip_y = TRUE) {
  if (is.null(coord_frame)) return(list(x = x, y = y))
  if (isTRUE(flip_y)) {
    list(
      x = coord_frame$xmin + as.numeric(x) - 1,
      y = coord_frame$ymax - as.numeric(y) + 1
    )
  } else {
    list(
      x = coord_frame$xmin + as.numeric(x) - 1,
      y = coord_frame$ymin + as.numeric(y) - 1
    )
  }
}

source_bbox_to_display_bbox <- function(bbox, coord_frame, flip_y = TRUE) {
  if (is.null(bbox) || is.null(coord_frame)) return(NULL)
  if (isTRUE(flip_y)) {
    tl <- source_xy_to_display_xy(bbox[["xmin"]], bbox[["ymax"]], coord_frame, flip_y = TRUE)
    br <- source_xy_to_display_xy(bbox[["xmax"]], bbox[["ymin"]], coord_frame, flip_y = TRUE)
  } else {
    tl <- source_xy_to_display_xy(bbox[["xmin"]], bbox[["ymin"]], coord_frame, flip_y = FALSE)
    br <- source_xy_to_display_xy(bbox[["xmax"]], bbox[["ymax"]], coord_frame, flip_y = FALSE)
  }
  c(
    xmin = tl$x,
    xmax = br$x,
    ymin = tl$y,
    ymax = br$y
  )
}

source_sf_to_display_sf <- function(poly_sf, coord_frame, flip_y = TRUE) {
  if (is.null(poly_sf) || nrow(poly_sf) == 0L || is.null(coord_frame)) return(poly_sf)
  rebuild_one <- function(one_row) {
    cc <- sf::st_coordinates(one_row)
    if (nrow(cc) == 0L) return(NULL)

    x1 <- cc[, "X"] - coord_frame$xmin + 1
    if (isTRUE(flip_y)) {
      y1 <- coord_frame$ymax - cc[, "Y"] + 1
    } else {
      y1 <- cc[, "Y"] - coord_frame$ymin + 1
    }
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

  geoms <- lapply(seq_len(nrow(poly_sf)), function(i) rebuild_one(poly_sf[i, , drop = FALSE]))
  keep <- !vapply(geoms, is.null, logical(1))
  poly_out <- poly_sf[keep, , drop = FALSE]
  poly_out <- sf::st_set_geometry(poly_out, sf::st_sfc(geoms[keep]))
  poly_out
}
