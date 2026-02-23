#function to plot the results of the statistical tests
#used for plotting in the MSI.EAGLE app as well as outputting the plots
#most information is contained within the reactive x5 object, so this
#function is mostly a wrapper for the plot functions in the Cardinal package

#' Plot statistics results
plot_stats_results <- function(
    x5,
    stats_table_rows_selected,
    stats_test,
    phen_cols_stats,
    group_var,
    plot_choice,
    aov_vars1 = NULL,
    aov_vars2 = NULL,
    aov_vars3 = NULL,
    output_factors = NULL,
    anova_type = NULL,
    chunks = getCardinalNChunks(),
    # ADD THESE NEW PARAMETERS WITH DEFAULTS:
    export_contrast = "histogram",
    export_colorscale = "Inferno",
    export_smooth = "gaussian",
    export_scale = TRUE,
    export_dark_bg = FALSE
) {
  if (is.null(stats_table_rows_selected)) {
    stop("No rows selected")
  }
  
  # Helper function to get color palette
  cpal <- function(name) {
    if(name %in% hcl.pals()) {
      hcl.colors(255, name)
    } else {
      viridisLite::viridis(255, option = "A")
    }
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

  build_level_color_map <- function(vals, fallback_palette = "Dark 3") {
    labs <- normalize_labels(vals)
    lvls <- unique(labs)
    if (length(lvls) == 0) {
      lvls <- "NA"
    }

    map <- stats::setNames(rep("grey70", length(lvls)), lvls)
    is_col <- are_valid_colors(lvls)
    if (any(is_col)) {
      map[is_col] <- lvls[is_col]
    }
    if (any(!is_col)) {
      fallback <- grDevices::hcl.colors(sum(!is_col), palette = fallback_palette)
      map[!is_col] <- fallback
    }
    map
  }

  build_groupings_plot <- function(obj, group_var, color_map) {
    pdat <- as.data.frame(pData(obj))
    if (!group_var %in% colnames(pdat)) {
      stop(paste0("Grouping variable '", group_var, "' not found in pData"))
    }

    labels <- normalize_labels(pdat[[group_var]])
    coords_df <- as.data.frame(Cardinal::coord(obj))
    xv <- if ("x" %in% colnames(coords_df)) coords_df$x else coords_df[[1]]
    yv <- if ("y" %in% colnames(coords_df)) coords_df$y else if (ncol(coords_df) >= 2) coords_df[[2]] else rep(1, length(xv))
    runv <- as.character(Cardinal::run(obj))
    if (length(runv) != length(labels)) {
      runv <- rep("run", length(labels))
    }

    lvls <- names(color_map)
    plot_df <- data.frame(
      x = as.numeric(xv),
      y = as.numeric(yv),
      run = runv,
      label = factor(labels, levels = lvls),
      stringsAsFactors = FALSE
    )

    ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y, fill = label)) +
      ggplot2::geom_tile() +
      ggplot2::facet_wrap(. ~ run) +
      ggplot2::coord_fixed() +
      ggplot2::scale_y_continuous(trans = "reverse") +
      ggplot2::theme_minimal() +
      ggplot2::scale_fill_manual(
        values = unname(color_map[lvls]),
        name = group_var,
        drop = FALSE
      )
  }
  
  # Setup to ensure proper cleanup
  on.exit({
    if (exists("par_mode") && !is.null(par_mode())) {
      bpstop(par_mode())  # Ensure parallel workers are stopped
    }
  }, add = TRUE)
  
  #set chunk size
  setCardinalNChunks(chunks)
  
  dat <- x5$stats_table_filtered
  
  
  # Get row from original dataset corresponding to selected data
  m <- which(x5$stats_results$i %in% dat[stats_table_rows_selected, ]$i)
  
  
  
  if (stats_test == "ssctest") {
    m_model <- dat[stats_table_rows_selected, ]$model
    
    # Find which rows of x5$stats_results match both ion and model
    idx <- NULL
    for (i in seq_along(m)) {
      condition <- x5$stats_results$i %in% dat[stats_table_rows_selected[i], ]$i &
        x5$stats_results$model %in% m_model[i]
      if (sum(condition) > 0) {
        idx[i] <- which(condition)
      }
    }
  }
  
  if (stats_test %in% c("meanstest")) {
    if (is.null(x5$stats_results)) {
      stop("x5$stats_results is null")
    }
    
    phen_all <- as.data.frame(pData(x5$data_file_selected))[, phen_cols_stats]
    color_map <- build_level_color_map(phen_all)
    mycols <- unname(color_map)
    
    nplots <- length(stats_table_rows_selected)
    mplot <- FALSE
    
    if (plot_choice == "ggplot") {
      
      i_vals <- x5$stats_results$i[m]
      mz_vals <- x5$stats_results$mz[m]
      fdr_vals <- x5$stats_results$fdr[m]
      ion_names <- paste0("mz=", round(mz_vals, 4), " i=", i_vals, " FDR=", round(fdr_vals, 3))

      # Extract all selected ions in one subset call (much faster on large datasets)
      pdat_n <- nrow(as.data.frame(pData(x5$data_file_selected)))
      feat_idx <- rep(NA_integer_, length(m))
      fdat_df <- try(as.data.frame(fData(x5$data_file_selected)), silent = TRUE)

      if (!inherits(fdat_df, "try-error") && "i" %in% colnames(fdat_df)) {
        feat_idx <- match(i_vals, suppressWarnings(as.numeric(fdat_df$i)))
      }
      if (anyNA(feat_idx)) {
        mz_lookup <- suppressWarnings(as.numeric(mz(x5$data_file_selected)))
        feat_idx_mz <- match(mz_vals, mz_lookup)
        feat_idx[is.na(feat_idx)] <- feat_idx_mz[is.na(feat_idx)]
      }

      feat_idx_unique <- unique(feat_idx[is.finite(feat_idx)])
      if (length(feat_idx_unique) == 0) {
        stop("Could not match selected ions to dataset features for ggplot plotting.")
      }

      subset_data <- x5$data_file_selected[feat_idx_unique, ]
      spec_mat <- as.matrix(spectra(subset_data))
      if (!is.matrix(spec_mat)) {
        spec_mat <- matrix(spec_mat, nrow = length(feat_idx_unique))
      }

      # Normalize orientation to [features x pixels]
      if (ncol(spec_mat) == pdat_n) {
        spec_feat_pix <- spec_mat
      } else if (nrow(spec_mat) == pdat_n) {
        spec_feat_pix <- t(spec_mat)
      } else {
        stop(sprintf(
          "Unexpected spectra matrix shape for ggplot plotting: %s x %s (pixels=%s)",
          nrow(spec_mat), ncol(spec_mat), pdat_n
        ))
      }

      idx_map <- match(feat_idx, feat_idx_unique)
      dat_mat <- matrix(NA_real_, nrow = pdat_n, ncol = length(m))
      valid_cols <- which(!is.na(idx_map))
      if (length(valid_cols) > 0) {
        dat_mat[, valid_cols] <- t(spec_feat_pix[idx_map[valid_cols], , drop = FALSE])
      }
      dat_comb <- as.data.frame(dat_mat, check.names = FALSE, stringsAsFactors = FALSE)
      colnames(dat_comb) <- ion_names
      ion_cols <- colnames(dat_comb)
      if (anyNA(idx_map)) {
        warning(sprintf(
          "Could not match %d selected ion(s) to the current dataset; plotting with NA columns for unmatched ions.",
          sum(is.na(idx_map))
        ))
      }
      
      # Convert pData to proper data frame before subsetting.
      # Keep grouping columns unique to avoid duplicate-name artifacts (e.g., polygon_is_cell.1).
      pdat_sel <- as.data.frame(pData(x5$data_file_selected))
      group_cols <- unique(c(phen_cols_stats, group_var))
      group_cols <- group_cols[!is.na(group_cols) & nzchar(group_cols)]
      group_cols <- group_cols[group_cols %in% colnames(pdat_sel)]
      if (!phen_cols_stats %in% group_cols && phen_cols_stats %in% colnames(pdat_sel)) {
        group_cols <- c(phen_cols_stats, group_cols)
      }
      b <- pdat_sel[, group_cols, drop = FALSE]
      
      # Add phenotype data
      dat_comb <- cbind(dat_comb, b)
      
      # Check if this is pixel-level data (each pixel is a replicate)
      # Don't average if group_var makes each row unique
      # Now b is a proper data frame so interaction() will work
      is_pixel_level <- is.null(group_var) || 
        all(group_var %in% c("none", "1")) ||
        nrow(b) == length(unique(do.call(paste, c(b, sep = "_"))))
      
      if (is_pixel_level) {
        # Pixel-level: each row is an individual pixel/replicate
        message("Plotting pixel-level data - showing distribution within groups")

        long_format_df <- dat_comb %>%
          tidyr::pivot_longer(
            cols = dplyr::all_of(ion_cols),
            names_to = "ion",
            values_to = "value"
          ) %>% 
          dplyr::mutate(value = as.numeric(value)) %>%
          dplyr::filter(is.finite(value))
        
      } else {
        # Biological replicate level: average by grouping variable
        summarized_df <- dat_comb %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
          dplyr::summarize(dplyr::across(dplyr::all_of(ion_cols), ~mean(.x, na.rm = TRUE)), .groups = "drop")
        
        long_format_df <- summarized_df %>%
          tidyr::pivot_longer(
            cols = dplyr::all_of(ion_cols),
            names_to = "ion",
            values_to = "value"
          ) %>% 
          dplyr::mutate(value = as.numeric(value)) %>%
          dplyr::filter(is.finite(value))
      }
      
      n_samples <- long_format_df %>%
        dplyr::group_by(ion, !!rlang::sym(phen_cols_stats)) %>%
        dplyr::summarize(n = dplyr::n(),
                         min_value = min(value, na.rm = TRUE),
                         .groups = 'drop')

      long_format_df[[phen_cols_stats]] <- normalize_labels(long_format_df[[phen_cols_stats]])
      n_samples[[phen_cols_stats]] <- normalize_labels(n_samples[[phen_cols_stats]])
      
      p1 <- ggplot2::ggplot(long_format_df,
                            ggplot2::aes(x = !!ggplot2::sym(phen_cols_stats),
                                         y = value,
                                         fill = !!ggplot2::sym(phen_cols_stats))) +
        ggplot2::geom_jitter(alpha = 0.5, width = 0.2) +
        ggplot2::geom_boxplot(size = 1, alpha = 0.8) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          x = "",
          y = "Intensity (normalized, a.u.)"
        ) +
        ggprism::theme_prism() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::facet_wrap(~ion, scales = "free_y") +
        ggplot2::scale_fill_manual(values = color_map, drop = FALSE) +
        ggplot2::geom_text(
          data = n_samples,
          ggplot2::aes(x = !!ggplot2::sym(phen_cols_stats),
                       y = min_value - 0.05 * abs(min_value),
                       label = paste0("n=", n)),
          vjust = 1.5
        )
      
      return(p1)
      
    } else if (plot_choice == "cardinal"){
      
      # Cardinal plots require a valid Cardinal meansTest result object.
      # For pixel-level/fallback workflows, auto-fallback to ggplot instead of crashing.
      if (is.null(x5$test_result)) {
        message("Cardinal means model unavailable; switching to ggplot output.")
        return(plot_stats_results(
          x5 = x5,
          stats_table_rows_selected = stats_table_rows_selected,
          stats_test = stats_test,
          phen_cols_stats = phen_cols_stats,
          group_var = group_var,
          plot_choice = "ggplot",
          aov_vars1 = aov_vars1,
          aov_vars2 = aov_vars2,
          aov_vars3 = aov_vars3,
          output_factors = output_factors,
          anova_type = anova_type,
          chunks = chunks,
          export_contrast = export_contrast,
          export_colorscale = export_colorscale,
          export_smooth = export_smooth,
          export_scale = export_scale,
          export_dark_bg = export_dark_bg
        ))
      }
      
      m <- stats_table_rows_selected
      mz_vals <- x5$stats_table_filtered$mz[m]
      i_vals <- x5$stats_table_filtered$i[m]
      fdr_vals <- x5$stats_table_filtered$fdr[m]
      names(i_vals) <- paste0("mz= ", round(mz_vals, 4), " FDR= ", formatC(fdr_vals, digits = 3, format = "fg"))
      
      p1 <- try(
        plot(x5$test_result, i = i_vals, col = mycols, las = 0, fill = TRUE, free = "xy"),
        silent = TRUE
      )
      if (inherits(p1, "try-error")) {
        message("Cardinal means plot failed (likely fallback model). Switching to ggplot. Error: ", as.character(p1))
        return(plot_stats_results(
          x5 = x5,
          stats_table_rows_selected = stats_table_rows_selected,
          stats_test = stats_test,
          phen_cols_stats = phen_cols_stats,
          group_var = group_var,
          plot_choice = "ggplot",
          aov_vars1 = aov_vars1,
          aov_vars2 = aov_vars2,
          aov_vars3 = aov_vars3,
          output_factors = output_factors,
          anova_type = anova_type,
          chunks = chunks,
          export_contrast = export_contrast,
          export_colorscale = export_colorscale,
          export_smooth = export_smooth,
          export_scale = export_scale,
          export_dark_bg = export_dark_bg
        ))
      }
      print(p1)
      if (exists("par_mode")) bpstop(par_mode())
      return(p1)
      
    } else if (plot_choice=="msi_image"){
      
      # Build arguments list for image function
      img_args <- list(
        x5$data_file_selected,
        mz = (x5$stats_results$mz[m]),
        scale = export_scale,
        free = "xy"
      )
      
      # Add colorscale
      if(export_colorscale != "none") {
        img_args$col <- cpal(export_colorscale)
      }
      
      # Add contrast enhancement
      if(export_contrast != "none") {
        img_args$enhance <- export_contrast
      }
      
      # Add smoothing
      if(export_smooth != "none") {
        img_args$smooth <- export_smooth
      }
      
      # Add dark background
      if(export_dark_bg) {
        img_args$style <- "dark"
      }
      
      p2 <- do.call(Cardinal::image, img_args)
      return(p2)
      
    } else if (plot_choice=="groupings"){
      
      p3 <- build_groupings_plot(x5$data_file_selected, phen_cols_stats, color_map)
      return(p3)
    }
    
  } else if (stats_test %in% c("spatialDGMM")) {
    if (is.null(x5$stats_results)) {
      stop("x5$stats_results is null")
    }
    
    phen_all <- as.data.frame(pData(x5$data_file_selected))[, phen_cols_stats]
    color_map <- build_level_color_map(phen_all)
    mycols <- unname(color_map)
    
    mz_vals <- x5$stats_results$mz[m]
    i_vals <- x5$stats_results$i[m]
    fdr_vals <- x5$stats_results$fdr[m]
    names(i_vals) <- paste0("mz= ", round(mz_vals, 4), " FDR= ", round(fdr_vals, 3))
    
    if (length(m) > 1) {
      message("Multiple ions selected, only first ion will be plotted")
      idx <- m[1]
    } else {
      idx <- m
    }
    
    p1 <- plot(x5$test_result, i = i_vals, fill = TRUE, free = "xy")
    p2 <- plot(x5$test_result_feature_test, i = i_vals, col = mycols, fill = TRUE, free = "xy")
    p3 <- image(x5$test_result, i = i_vals, smooth = "bilateral", enhance = "adaptive", scale = TRUE)
    p4 <- image(
      x5$data_file_selected,
      mz = (x5$stats_results$mz[idx]),
      enhance = "histogram",
      scale = TRUE,
      free = "xy"
    )
    
    plot_choice <- switch(plot_choice,
                          "dgmm_means_test" = p2,
                          "dgmm_ranks" = p3,
                          "dgmm_params" = p1,
                          "msi_image" = p4
    )
    
    print(plot_choice)
    return(plot_choice)
    
  } else if (stats_test == "anova") {
    # Implement ANOVA plotting code here
    # ...
    message("ANOVA plotting not fully implemented in this function.")
    return(NULL)
    
  } else if (stats_test == "ssctest") {
    
    
    
    # if (length(idx) > 1) {
    #   message("Multiple ions selected, only first ion will be plotted for certain ssc plots")
    #   m <- idx[1]
    # }
    # 
    # 
    if(plot_choice=="means_plot") {
      
      #loop through unique mz values to create means plots
      mz_to_plot<-unique(x5$stats_results$mz[m])
      plot_list <- list()
      for(j in 1:length(mz_to_plot)) {
        
        a <- as.matrix(spectra(
          subsetFeatures(x5$data_file_selected, mz %in% mz_to_plot[j])
        ))
        b <- pData(subsetFeatures(x5$data_file_selected, mz %in% mz_to_plot[j]))
        dat <- data.frame(run = Cardinal::run(x5$data_file_selected), (b), value = t(a))
        
        dat_long_tech_avg <- as.data.frame(dat) %>%
          dplyr::group_by_at(unique(c("run", phen_cols_stats))) %>%
          dplyr::summarize(tech_avg = mean(value), .groups = "keep")
        
        fm <- as.formula(paste0("tech_avg~(", phen_cols_stats, ")"))
        
        dat_long_tech_avg[[phen_cols_stats]] <- normalize_labels(dat_long_tech_avg[[phen_cols_stats]])
        
        #create plot with means using ggplot2
        plot_list[[j]]<- ggplot2::ggplot(dat_long_tech_avg, ggplot2::aes_string(x=phen_cols_stats, y="tech_avg", fill=phen_cols_stats)) +
          ggplot2::geom_boxplot(alpha=0.5) +
          ggplot2::geom_point(ggplot2::aes(y=tech_avg), position=ggplot2::position_dodge(width=0.75), alpha=0.5) +
          ggplot2::theme_minimal() +
          ggplot2::labs(x="", y="Mean (norm. intensity, a.u.)") +
          ggprism::theme_prism() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
          ggplot2::theme(legend.position = "none") +
          ggplot2::scale_fill_manual(values = color_map, drop = FALSE)+
          #add title
          ggplot2::ggtitle(paste("mz =",round(mz_to_plot[j],4)))
        
      }
      
      p5<-gridExtra::grid.arrange(grobs=plot_list)
    }
    
    phen_all <- as.data.frame(pData(x5$data_file_selected))[, phen_cols_stats]
    color_map <- build_level_color_map(phen_all)
    mycols <- unname(color_map)
    
    
    
    
    if (is.null(names(x5$test_result))) {
      p1 <- image(x5$test_result[[1]], key = T, col = mycols)
    } else {
      p1 <- image(x5$test_result, i = unique(m_model), key = T, col = mycols)
    }
    
    mz_vals <- x5$stats_results$mz[m]
    i_vals <- x5$stats_results$i[m]
    names(i_vals) <- paste0("mz= ", round(mz_vals, 4))
    
    
    
    
    if (length(mz_vals) > 1) {
      message("Multiple ions selected, only first ion will be plotted")
      idx <- idx[1]
    }
    
    p2 <- image(
      x5$data_file_selected,
      mz = (x5$stats_results$mz[idx]),
      model = m_model,
      enhance = "histogram",
      scale = TRUE,
      free = "xy"
    )
    
    fm2 <- (paste0(phen_cols_stats))
    
    if (class(x5$test_result) == "SpatialShrunkenCentroids") {
      x5$test_result <- list(x5$test_result)
    }
    
    if (length(x5$test_result) > 1) {
      plot_list <- list()
      nmodels <- length(x5$test_result)
      for (i in seq_len(nmodels)) {
        plot_list[[i]] <- plot(x5$test_result[[i]], type = "statistic", linewidth = 2, col = mycols)
      }
      names(plot_list) <- names(x5$test_result)
    } else {
      plot_list <- plot(x5$test_result[[1]], type = "statistic", linewidth = 2, col = mycols)
    }
    p3 <- matter::as_facets(plot_list, ncol = 1)
    
    p4 <- build_groupings_plot(
      x5$data_file_selected[, Cardinal::run(x5$data_file_selected) %in% runNames(x5$data_file_selected)[]],
      fm2,
      color_map
    )
    
    
    plot_out <- switch(
      plot_choice,
      "ion_image" = p1,
      "means_plot" = p5,
      "t_statistic" = p3,
      "groupings" = p4,
      "msi_image" = p2
    )
    
    print(plot_out)
    return(plot_out)
    
  } else {
    message("Plot not supported yet")
    return(NULL)
  }
}
