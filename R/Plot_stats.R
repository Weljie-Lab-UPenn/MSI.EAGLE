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
    export_colorscale = "Spectral",
    export_smooth = "none",
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
    
    mycols <- ggsci::pal_npg()(length(levels(droplevels(factor(
      as.data.frame(pData(x5$data_file_selected))[, phen_cols_stats]
    )))))
    
    nplots <- length(stats_table_rows_selected)
    mplot <- FALSE
    
    if (plot_choice == "ggplot") {
      
      # Extract spectra for selected ions - properly handle Cardinal S4 objects
      a <- lapply(m, function(x) {
        # Get the subset for this specific m/z
        subset_data <- x5$data_file_selected[which(mz(x5$data_file_selected) %in% x5$stats_results$mz[x]), ]
        # Extract spectra as matrix and return as vector
        as.vector(as.matrix(spectra(subset_data)))
      })
      
      i_vals <- x5$stats_results$i[m]
      mz_vals <- x5$stats_results$mz[m]
      fdr_vals <- x5$stats_results$fdr[m]
      names(a) <- paste0("mz=", round(mz_vals, 4), " i=", i_vals, " FDR=", round(fdr_vals, 3))
      
      # Convert pData to proper data frame before subsetting
      b <- as.data.frame(pData(x5$data_file_selected))[, c(phen_cols_stats, group_var), drop = FALSE]
      
      # Combine data - ensure proper dimensions
      if (length(a) == 1) {
        # Single ion selected
        dat_comb <- as.data.frame(matrix(a[[1]], ncol = 1))
        colnames(dat_comb) <- names(a)[1]
      } else {
        # Multiple ions - combine as columns
        dat_comb <- as.data.frame(do.call(cbind, a))
        colnames(dat_comb) <- names(a)
      }
      
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
        
        # Get column names for ions (exclude phenotype columns)
        ion_cols <- setdiff(colnames(dat_comb), c(phen_cols_stats, group_var))
        
        long_format_df <- dat_comb %>%
          tidyr::pivot_longer(
            cols = dplyr::all_of(ion_cols),
            names_to = "ion",
            values_to = "value"
          ) %>% 
          na.omit() %>%
          dplyr::mutate(value = as.numeric(value))
        
      } else {
        # Biological replicate level: average by grouping variable
        summarized_df <- dat_comb %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(c(phen_cols_stats, group_var)))) %>%
          dplyr::summarize(dplyr::across(dplyr::everything(), mean, na.rm = TRUE), .groups = 'keep')
        
        # Get column names for ions (exclude phenotype columns)
        ion_cols <- setdiff(colnames(summarized_df), c(phen_cols_stats, group_var))
        
        long_format_df <- summarized_df %>%
          tidyr::pivot_longer(
            cols = dplyr::all_of(ion_cols),
            names_to = "ion",
            values_to = "value"
          ) %>% 
          na.omit() %>%
          dplyr::mutate(value = as.numeric(value))
      }
      
      n_samples <- long_format_df %>%
        dplyr::group_by(ion, !!rlang::sym(phen_cols_stats)) %>%
        dplyr::summarize(n = dplyr::n(),
                         min_value = min(value, na.rm = TRUE),
                         .groups = 'drop')
      
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
        ggplot2::scale_fill_manual(values = mycols) +
        ggplot2::geom_text(
          data = n_samples,
          ggplot2::aes(x = !!ggplot2::sym(phen_cols_stats),
                       y = min_value - 0.05 * abs(min_value),
                       label = paste0("n=", n)),
          vjust = 1.5
        )
      
      return(p1)
      
    } else if (plot_choice == "cardinal"){
      
      # Cardinal plots require x5$test_result which doesn't exist for pixel-level tests
      if (is.null(x5$test_result)) {
        message("Cardinal plotting not available for pixel-level means test. Please use 'ggplot' option instead.")
        stop("Cardinal plotting requires biological replicates with samples parameter. Use ggplot for pixel-level data.")
      }
      
      m <- stats_table_rows_selected
      mz_vals <- x5$stats_table_filtered$mz[m]
      i_vals <- x5$stats_table_filtered$i[m]
      fdr_vals <- x5$stats_table_filtered$fdr[m]
      names(i_vals) <- paste0("mz= ", round(mz_vals, 4), " FDR= ", formatC(fdr_vals, digits = 3, format = "fg"))
      
      p1 <- plot(x5$test_result, i = i_vals, col = mycols, las = 0, fill = TRUE, free = "xy")
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
      
      p3 <- image(
        x5$data_file_selected,
        phen_cols_stats,
        key = TRUE,
        col = mycols
      )
      return(p3)
    }
    
  } else if (stats_test %in% c("spatialDGMM")) {
    if (is.null(x5$stats_results)) {
      stop("x5$stats_results is null")
    }
    
    mycols <- ggsci::pal_npg()(length(levels(droplevels(factor(
      as.data.frame(pData(x5$data_file_selected))[, phen_cols_stats]
    )))))
    
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
        
        ncolors <- length(levels(as.factor(as.data.frame(
          pData(x5$data_file_selected)[, phen_cols_stats]
        )[,1])))
        
        mycols <- ggsci::pal_npg()(ncolors)
        
        #create plot with means using ggplot2
        plot_list[[j]]<- ggplot2::ggplot(dat_long_tech_avg, ggplot2::aes_string(x=phen_cols_stats, y="tech_avg", fill=phen_cols_stats)) +
          ggplot2::geom_boxplot(alpha=0.5) +
          ggplot2::geom_point(ggplot2::aes(y=tech_avg), position=ggplot2::position_dodge(width=0.75), alpha=0.5) +
          ggplot2::theme_minimal() +
          ggplot2::labs(x="", y="Mean (norm. intensity, a.u.)") +
          ggprism::theme_prism() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
          ggplot2::theme(legend.position = "none") +
          ggplot2::scale_fill_manual(values = mycols)+
          #add title
          ggplot2::ggtitle(paste("mz =",round(mz_to_plot[j],4)))
        
      }
      
      p5<-gridExtra::grid.arrange(grobs=plot_list)
    }
    
    ncolors <- length(levels(as.factor(as.data.frame(
      pData(x5$data_file_selected)[, phen_cols_stats]
    )[,1])))
    
    mycols <- ggsci::pal_npg()(ncolors)
    
    
    
    
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
    
    p4 <- image(
      x5$data_file_selected[, Cardinal::run(x5$data_file_selected) %in% runNames(x5$data_file_selected)[]],
      fm2,
      key = TRUE,
      col = mycols
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
