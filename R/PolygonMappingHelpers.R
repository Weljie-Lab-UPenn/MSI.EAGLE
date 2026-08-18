rescue_zero_hit_polygons <- function(hit, pixel_xy, polygon_xy, polygon_index,
                                     polygon_area, eligible_mask,
                                     overlap_rule = "smallest",
                                     max_shift_steps = 2) {
  n_pixels <- length(hit)
  n_polygons <- length(polygon_area)
  if (n_pixels == 0L || n_polygons == 0L || length(polygon_index) == 0L) {
    return(list(hit = hit, assignments = data.frame()))
  }

  pixel_xy <- as.matrix(pixel_xy)
  polygon_xy <- as.matrix(polygon_xy)
  if (nrow(pixel_xy) != n_pixels || ncol(pixel_xy) < 2L ||
      nrow(polygon_xy) != length(polygon_index) || ncol(polygon_xy) < 2L) {
    stop("Polygon rescue coordinates do not match their pixel/polygon indices.")
  }

  polygon_index <- suppressWarnings(as.integer(polygon_index))
  polygon_area <- suppressWarnings(as.numeric(polygon_area))
  polygon_area[!is.finite(polygon_area) | polygon_area <= 0] <- Inf
  eligible_mask <- rep_len(as.logical(eligible_mask), n_polygons)
  eligible_mask[is.na(eligible_mask)] <- FALSE
  overlap_rule <- tolower(trimws(as.character(overlap_rule)[1]))
  if (!overlap_rule %in% c("smallest", "first", "all")) overlap_rule <- "smallest"

  x_step <- suppressWarnings(stats::median(diff(sort(unique(pixel_xy[, 1]))), na.rm = TRUE))
  y_step <- suppressWarnings(stats::median(diff(sort(unique(pixel_xy[, 2]))), na.rm = TRUE))
  if (!is.finite(x_step) || x_step <= 0) x_step <- 1
  if (!is.finite(y_step) || y_step <= 0) y_step <- 1
  max_step <- max(x_step, y_step)
  max_shift_steps <- suppressWarnings(as.numeric(max_shift_steps))
  if (!is.finite(max_shift_steps) || max_shift_steps < 1) max_shift_steps <- 2
  max_distance <- max_shift_steps * max_step

  rescue_area <- polygon_area[polygon_index]
  rescue_area[!is.finite(rescue_area) | rescue_area <= 0] <- Inf
  rescue_order <- order(rescue_area, polygon_index, na.last = TRUE)
  reserved_pixels <- integer(0)
  assignments <- vector("list", length(rescue_order))

  would_survive <- function(pixel_index, poly_index) {
    existing <- suppressWarnings(as.integer(hit[[pixel_index]]))
    existing <- existing[
      is.finite(existing) & existing >= 1L & existing <= n_polygons &
        eligible_mask[existing]
    ]
    if (length(existing) == 0L || identical(overlap_rule, "all")) return(TRUE)
    if (identical(overlap_rule, "first")) return(poly_index < min(existing))

    area_i <- polygon_area[existing]
    area_i[!is.finite(area_i) | area_i <= 0] <- Inf
    min_area <- min(area_i)
    if (polygon_area[poly_index] < min_area) return(TRUE)
    if (!isTRUE(all.equal(polygon_area[poly_index], min_area))) return(FALSE)
    poly_index < min(existing[area_i == min_area])
  }

  for (k in seq_along(rescue_order)) {
    j <- rescue_order[k]
    poly_index <- polygon_index[j]
    if (!is.finite(poly_index) || poly_index < 1L || poly_index > n_polygons) next

    dx <- pixel_xy[, 1] - polygon_xy[j, 1]
    dy <- pixel_xy[, 2] - polygon_xy[j, 2]
    distance2 <- dx * dx + dy * dy
    candidate_order <- order(distance2, seq_len(n_pixels), na.last = TRUE)
    candidate_order <- candidate_order[
      is.finite(distance2[candidate_order]) &
        sqrt(distance2[candidate_order]) <= (max_distance + sqrt(.Machine$double.eps))
    ]
    if (length(candidate_order) == 0L) next

    unreserved <- candidate_order[!candidate_order %in% reserved_pixels]
    viable <- unreserved[vapply(
      unreserved,
      function(px) would_survive(px, poly_index),
      logical(1)
    )]
    chosen <- if (length(viable) > 0L) viable[1] else if (length(unreserved) > 0L) unreserved[1] else candidate_order[1]

    survives_rule <- would_survive(chosen, poly_index)
    hit[[chosen]] <- sort(unique(c(suppressWarnings(as.integer(hit[[chosen]])), poly_index)))
    reserved_pixels <- c(reserved_pixels, chosen)
    assignments[[k]] <- data.frame(
      polygon_index = poly_index,
      pixel_index = chosen,
      distance = sqrt(distance2[chosen]),
      survives_rule = survives_rule,
      stringsAsFactors = FALSE
    )
  }

  assignments <- assignments[!vapply(assignments, is.null, logical(1))]
  assignments <- if (length(assignments) > 0L) do.call(rbind, assignments) else data.frame()
  list(hit = hit, assignments = assignments)
}
