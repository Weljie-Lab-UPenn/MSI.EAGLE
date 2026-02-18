### R/UMAPServer.R
UMAPServer <- function(id, setup_values, preproc_values, preproc_values_umap = NULL) {
  moduleServer(id, function(input, output, session) {
    #for dynamic UI
    ns = session$ns
    
    output$variables <- renderPrint({
      setup_vals <- setup_values()
      paste("working directory=", setup_vals$wd)
    })
    
    #import parallel mode
    par_mode = reactive({
      setup_values()[["par_mode"]]
    })
    
    sanitize_tf_mask <- function(tf, n_pixels) {
      if (is.null(tf)) {
        return(NULL)
      }
      tf <- as.logical(tf)
      if (length(tf) != n_pixels) {
        return(NULL)
      }
      tf[is.na(tf)] <- FALSE
      tf
    }
    
    safe_subset_pixels <- function(dat, idx, context = "selection", notify = TRUE) {
      if (is.null(dat)) {
        return(NULL)
      }
      n_pixels <- ncol(dat)
      if (n_pixels == 0) {
        return(dat)
      }
      if (is.null(idx)) {
        return(dat)
      }
      
      if (is.logical(idx)) {
        idx <- sanitize_tf_mask(idx, n_pixels)
        if (is.null(idx)) {
          if (isTRUE(notify)) {
            showNotification(
              sprintf("%s mask length does not match number of pixels; skipping subsetting.", context),
              type = "warning",
              duration = 5
            )
          }
          return(dat)
        }
        n_sel <- sum(idx)
        if (!is.finite(n_sel) || n_sel <= 0) {
          if (isTRUE(notify)) {
            showNotification(sprintf("No pixels selected for %s.", context), type = "message", duration = 4)
          }
          return(NULL)
        }
        if (n_sel == n_pixels) {
          return(dat)
        }
      } else {
        idx <- suppressWarnings(as.integer(idx))
        idx <- unique(idx[is.finite(idx) & idx >= 1 & idx <= n_pixels])
        if (length(idx) == 0) {
          if (isTRUE(notify)) {
            showNotification(sprintf("No valid pixel indices for %s.", context), type = "message", duration = 4)
          }
          return(NULL)
        }
        if (length(idx) == n_pixels && identical(sort(idx), seq_len(n_pixels))) {
          return(dat)
        }
      }
      
      out <- try(subsetPixels(dat, idx), silent = TRUE)
      if (inherits(out, "try-error")) {
        if (isTRUE(notify)) {
          showNotification(
            sprintf("Unable to subset pixels for %s. Please adjust the selected filters.", context),
            type = "error",
            duration = 8
          )
        }
        return(NULL)
      }
      out
    }
    
    get_active_pdata <- function(x2) {
      if (is.null(x2$mytable_selected)) {
        return(NULL)
      }
      n_pixels <- ncol(x2$mytable_selected)
      if (!is.null(x2$pdat_anat) && nrow(x2$pdat_anat) == n_pixels) {
        return(x2$pdat_anat)
      }
      pData(x2$mytable_selected)
    }
    
    get_active_pdata_df <- function(x2) {
      pdat <- get_active_pdata(x2)
      if (is.null(pdat)) {
        return(NULL)
      }
      as.data.frame(pdat)
    }
    
    draw_empty_plot <- function(msg) {
      graphics::plot.new()
      graphics::text(0.5, 0.5, labels = msg, cex = 1.1)
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
              col2rgb(val)
              TRUE
            }, error = function(e) FALSE)
          },
          logical(1)
        )
      }
      ok
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
    
    build_fill_mapping <- function(labels, fallback_palette = "Dark 3") {
      labels <- normalize_labels(labels)
      uniq <- unique(labels)
      valid_color <- are_valid_colors(uniq)
      
      fill_values <- stats::setNames(rep("grey70", length(uniq)), uniq)
      if (any(valid_color)) {
        fill_values[valid_color] <- uniq[valid_color]
      }
      if (any(!valid_color)) {
        non_color <- uniq[!valid_color]
        pal <- grDevices::hcl.colors(max(length(non_color), 1), palette = fallback_palette)
        fill_values[non_color] <- pal[seq_along(non_color)]
      }
      
      list(labels = labels, values = fill_values)
    }
    
    log_fill_mapping <- function(tag, fill_values, n_show = 12) {
      if (length(fill_values) == 0) {
        return(invisible(NULL))
      }
      n <- min(length(fill_values), n_show)
      msg <- paste(
        sprintf("%s=%s", names(fill_values)[seq_len(n)], as.character(fill_values)[seq_len(n)]),
        collapse = ", "
      )
      message(sprintf("[%s] color map (%d labels): %s", tag, length(fill_values), msg))
      invisible(NULL)
    }
    
    
    
    
    ###START UMAP for background removal
    
    
    
    
    #list runs for selection
    output$mytable = DT::renderDataTable({
      x0 <- preproc_values()[["x0"]]
      req(x0$overview_peaks)
      ovps <- x0$overview_peaks
      #x2<-preproc_values()[["x2"]]
      #req(x2$overview_peaks_sel)
      #ovps<-x2$overview_peaks_sel
      DT::datatable(
        cbind(
          run = runNames(ovps),
          size = lapply(runNames(ovps), function(x)
            ncol(ovps[, Cardinal::run(ovps) %in% x]))
        ),
        #selection = list(mode = "multiple"),
        selection = 'none',
        caption = "Choose runs to process for UMAP clustering",
        extensions = c("Buttons", "Select"),
        options = list(
          dom = 'Bfrtip',
          select = TRUE,
          buttons = list('pageLength', "copy", "selectNone", "selectAll")
        )
      )
    },
    server = FALSE)
    
    #extract selected datasets from table (x2$mytable_selected)
    observeEvent(input$mytable_rows_selected, {
      x2 <- preproc_values()[["x2"]]
      
      if (is.null(x2$list_proc_img)) {
        #only do this if there is no table already created for anatomical segmentation
        req(x2$overview_peaks_sel)
        req(input$mytable_rows_selected)
        
        ids <- input$mytable_rows_selected
        x2$mytable_selected <-
          x2$overview_peaks_sel %>% subsetPixels(
            Cardinal::run(x2$overview_peaks_sel) %in% runNames(x2$overview_peaks_sel)[ids]
          )
      }
      #x2$tf_list<-coord(x2$mytable_selected) %in% coord(x2$mytable_selected)
      
      # Update tf_list accordingly TEMP, check carefully
      #x2$tf_list <- rep(TRUE, ncol(x2$mytable_selected))
    })
    
    #update the x2$mytable_selected variable if there is a change to the selected rows.
    #use processed data if available... need to think about tf list
    observeEvent(input$mytable_rows_selected, {
      
      x2 <- preproc_values()[["x2"]]
      req(x2$mytable_selected)
      
      if (input$seg_choice == "anat_seg") {
        if (!is.null(x2$list_proc_img)) {
          print("updating selected table from previously processed phenotype data")
          
          #store runNames from previous mytable_selected to help with tf_list tracking
          
          old_mytable_selected <- x2$mytable_selected
          x2$mytable_selected <-
            combine_card(x2$list_proc_img[input$mytable_rows_selected])
          
          #update tf list to reflect new dataset. create list of all tf
          
          tf_list_list <-
            lapply(x2$list_proc_img, function(x)
              rep(TRUE, ncol(x)))
          names(tf_list_list) <-
            unlist(lapply(x2$list_proc_img, runNames))
          
          x2$tf_list <-
            unlist(tf_list_list[runNames(x2$mytable_selected)])
          x2$pdat_anat <- (pData(x2$mytable_selected))
          
          
          
        }
      } else if (input$seg_choice == "bk_seg"){
        
        if (!is.null(x2$list_proc_img)) {
          print("updating selected table from previously processed phenotype data")
          
          
          x2$mytable_selected <-
            combine_card(x2$list_proc_img[input$mytable_rows_selected])
          
          x2$tf_list <-
            rep(TRUE, ncol(x2$mytable_selected))
        }
      }
    })
    
    
    #create reactive temporary variable to store dataset being worked on which matches the t/f list
    #start with the inital dataset selected from the table.
    tmp_umap_dat <- reactive({
      x2 <- preproc_values()[["x2"]]
      req(x2$mytable_selected)
      if (!is.null(x2$list_proc_img)) {
        img.dat <-
          combine_card(x2$list_proc_img[input$mytable_rows_selected])
        #do.call(cbind, x2$list_proc_img[input$mytable_rows_selected])
        #Cardinal::combine(x2$list_proc_img[input$mytable_rows_selected])
        #x2$tf_list<-rep(TRUE, ncol(img.dat))
        return(img.dat)
      } else {
        #if there is no processed data, default to the mytable selected original data
        
        x2$tf_list <- rep(TRUE, ncol(x2$mytable_selected))
        return(x2$mytable_selected)
        
      }
    })
    
    
    #create tf_list and tf_list_umap
    observe( {
      #true false list corresponsds to size of selected entries in table.
      x2 <- preproc_values()[["x2"]]
      req(x2$mytable_selected)
      
      #if there is a selection in the table, create a true/false list
      if (dim(x2$mytable_selected)[2] > 0) {
        if (is.null(x2$tf_list) | length(x2$tf_list) == 0) {
          x2$tf_list <-
            !is.na(prodlim::row.match(as.data.frame(coord(
              x2$mytable_selected
            )), as.data.frame(coord(
              x2$mytable_selected
            ))))
        }
        #x2$tf_list<-coord(x2$mytable_selected) %in% coord(x2$mytable_selected)
        x2$tf_list_umap <- NULL
      }
    })
    #
    #
    
    
    
    #create the plot using plot_card_server()
    observe( {
      req(input$mytable_rows_selected)
      x2 <- preproc_values()[["x2"]]
      req(x2$mytable_selected)
      #req(tmp_umap_dat())
      req(x2$tf_list)
      
      if (input$seg_choice == "anat_seg") {
        #look for data in processed image area first
        if (!is.null(x2$list_proc_img)) {
          img.dat <- combine_card(x2$list_proc_img[input$mytable_rows_selected])
        } else if (!is.null(x2$tf_list_umap) &&
                   !is.null(x2$umap_name)) {
          if (identical(x2$umap_name, runNames(x2$mytable_selected))) {
            img.dat <- safe_subset_pixels(x2$mytable_selected, x2$tf_list_umap, "UMAP subset", notify = FALSE)
          }
        } else {
          img.dat <- x2$mytable_selected #%>% subsetPixels(x2$tf_list)
        }
      } else {
        
        
        if(sum(x2$tf_list)==ncol(x2$mytable_selected)) {
          img.dat <- x2$mytable_selected
        } else if (length(x2$tf_list) == ncol(x2$mytable_selected)) { 
          
          img.dat <- safe_subset_pixels(x2$mytable_selected, x2$tf_list, "current segmentation mask", notify = FALSE)
        } else {
          
          message("tf_list is not the same size as the dataset, exiting")
          showNotification("tf_list is not the same size as the dataset, exiting", type="error")
          return(NULL)
        }
        
      }
      
      
      
      if (is.null(img.dat) || ncol(img.dat) == 0) {
        return(NULL)
      }
      
      if (sum(x2$tf_list) == 0) {
        print("nothing selected in tf_list, exiting")
        return(NULL)
        
      }
      
      
      #if there is already processed data, use that
      if (!is.null(x2$list_proc_img)) {
        img.dat <-
          combine_card(x2$list_proc_img[input$mytable_rows_selected])
      }
      
      plot_card_server("plot_card_umap",
                       overview_peaks_sel = img.dat,
                       spatialOnly = TRUE)
      
    })
    
    #create the plot when dataset updated
    observeEvent(req(tmp_umap_dat()), {
      plot_card_server("plot_card_umap",
                       overview_peaks_sel = tmp_umap_dat(),
                       spatialOnly = TRUE)
      
    })
    
    
    
    
    #run fix_pix to remove isolated pixels
    observeEvent(input$fix_pix, {
      x2 <- preproc_values()[["x2"]]
      if (is.null(x2$mytable_selected)) {
        showNotification("No dataset selected. Choose runs from the table before removing isolated pixels.", type = "warning", duration = 6)
        return()
      }
      tmp_dat <- fix_pix(
        x2$mytable_selected,
        remove = input$fix_pix_t_f,
        r = input$fix_pix_r,
        n_thresh = input$fix_pix_n_thresh
      )
      if (!is.null(tmp_dat)) {
        x2$mytable_selected <- tmp_dat
        
      }
      print("done fix_pix")
      
    })
    
    output$fm_params <- renderUI({
      if (input$fastmap == TRUE) {
        tagList(
          fluidRow(
            column(6, numericInput(ns("fm_r"), "Fastmap r", value = 1)),
            column(
              6,
              numericInput(
                ns("fm_ncomp"),
                "Fastmap components",
                value = 20,
                min = 3
              )
            )
          ),
          fluidRow(
            column(
              6,
              selectInput(
                ns("fm_method"),
                "Fastmap weights",
                choices = c("gaussian", "adaptive")
              )
            ),
            column(
              6,
              selectInput(
                ns("fm_metric"),
                "Fastmap metric",
                choices = c("euclidean", "manhattan", "cosine", "correlation"),
                selected = "euclidean"
              )
            )
          )
        )
      }
      
    })
    
    
    
    
    # Load necessary libraries at the beginning of your server code
    # library(dbscan)
    # library(mclust)
    # library(cluster)
    # library(e1071)
    # library(skmeans)
    # library(kohonen)
    # Add any other required libraries here
    
    # Run UMAP
    observeEvent(input$umap1, {
      if (is.null(input$mytable_rows_selected)) {
        showNotification("Select a dataset from the table first", type = "warning")
        print("Select a dataset from the table first")
        return()
      }
      print("Starting UMAP")
      # Reset colorscale
      x2 <- preproc_values()[["x2"]]
      
      gc()
      ptm <- proc.time()
      
      withProgress(message = "Performing UMAP analysis and clustering", value = 0.2, {
        if (input$seg_choice == "anat_seg") {
          img.dat <- safe_subset_pixels(x2$mytable_selected, x2$tf_list, "UMAP input", notify = TRUE)
          if (is.null(img.dat) || ncol(img.dat) == 0) {
            showNotification("No pixels selected for UMAP. Adjust filters and try again.", type = "warning", duration = 6)
            return()
          }
          x2$tf_list_anat <- rep(TRUE, ncol(img.dat))
          x2$tf_list_umap <- x2$tf_list  # Store T/F list to store UMAP pixels
          x2$umap_name <- runNames(img.dat)
          x2$pdat_anat <- pData(x2$mytable_selected)
        } else {
          img.dat <- tmp_umap_dat()
          x2$rcol_plot <- NULL
        }
        
        # Parse/validate UMAP parameters
        coalesce_input <- function(x, default) {
          if (is.null(x) || length(x) == 0) default else x
        }
        parse_pca_umap <- function(x) {
          txt <- trimws(as.character(coalesce_input(x, "")))
          if (identical(txt, "") || tolower(txt) == "null") {
            return(NULL)
          }
          val <- suppressWarnings(as.integer(txt))
          if (!is.finite(val) || val < 1) {
            showNotification("Invalid PCA component value. Using NULL (no PCA).", type = "warning", duration = 6)
            return(NULL)
          }
          val
        }
        pca_val <- parse_pca_umap(input$pca_umap)
        nn_val <- max(2L, as.integer(round(coalesce_input(input$nn, 8))))
        n_trees_val <- max(1L, as.integer(round(coalesce_input(input$n_trees, 100))))
        search_k_val <- as.integer(round(coalesce_input(input$search_k, 1000L)))
        if (!is.finite(search_k_val) || search_k_val < 1L) {
          search_k_val <- 1000L
        }
        min_dist_val <- suppressWarnings(as.numeric(input$min_dist))
        if (!is.finite(min_dist_val) || min_dist_val < 0) min_dist_val <- 0.02
        set_op_mix_ratio_val <- suppressWarnings(as.numeric(input$set_op_mix_ratio))
        if (!is.finite(set_op_mix_ratio_val)) set_op_mix_ratio_val <- 0.15
        set_op_mix_ratio_val <- min(1, max(0, set_op_mix_ratio_val))
        metric_val <- as.character(coalesce_input(input$umap_metric, "cosine"))
        spread_val <- suppressWarnings(as.numeric(input$umap_spread))
        if (!is.finite(spread_val) || spread_val <= 0) spread_val <- 1
        local_connectivity_val <- suppressWarnings(as.numeric(input$local_connectivity))
        if (!is.finite(local_connectivity_val) || local_connectivity_val < 0) local_connectivity_val <- 1
        repulsion_strength_val <- suppressWarnings(as.numeric(input$repulsion_strength))
        if (!is.finite(repulsion_strength_val) || repulsion_strength_val <= 0) repulsion_strength_val <- 1
        negative_sample_rate_val <- as.integer(round(coalesce_input(input$negative_sample_rate, 10)))
        if (!is.finite(negative_sample_rate_val) || negative_sample_rate_val < 1L) negative_sample_rate_val <- 10L
        learning_rate_val <- suppressWarnings(as.numeric(input$learning_rate))
        if (!is.finite(learning_rate_val) || learning_rate_val <= 0) learning_rate_val <- 1
        init_val <- as.character(coalesce_input(input$umap_init, "spectral"))
        n_epochs_val <- as.integer(round(coalesce_input(input$n_epochs, 0)))
        if (!is.finite(n_epochs_val) || n_epochs_val <= 0L) {
          n_epochs_val <- NULL
        }
        
        umap_args <- list(
          img.dat = img.dat,
          outliers = FALSE,
          nn = nn_val,
          search_k = search_k_val,
          min_dist = min_dist_val,
          n_trees = n_trees_val,
          set_op_mix_ratio = set_op_mix_ratio_val,
          pca = pca_val,
          metric = metric_val,
          spread = spread_val,
          local_connectivity = local_connectivity_val,
          repulsion_strength = repulsion_strength_val,
          negative_sample_rate = negative_sample_rate_val,
          learning_rate = learning_rate_val,
          init = init_val
        )
        if (!is.null(n_epochs_val)) {
          umap_args$n_epochs <- n_epochs_val
        }
        
        # Run UMAP
        if (input$fastmap == TRUE) {
          setCardinalBPPARAM(par_mode())
          umap_args$fastmap <- TRUE
          umap_args$fm_r <- input$fm_r
          umap_args$fm_method <- input$fm_method
          umap_args$fm_metric <- input$fm_metric
          umap_args$fm_ncomp <- input$fm_ncomp
          data_list <- try(do.call(get_umap, umap_args), silent = TRUE)
          if (class(data_list) %in% "try-error") {
            err_txt <- gsub("\\s+", " ", as.character(data_list))
            err_txt <- sub("^Error in .*?:\\s*", "", err_txt)
            err_txt <- substr(err_txt, 1, 180)
            showNotification(
              paste0("Fastmap UMAP failed: ", err_txt),
              type = "error",
              duration = 10
            )
            message("Fastmap UMAP failed: ", as.character(data_list))
            return()
          }
        } else {
          data_list <- try(do.call(get_umap, umap_args), silent = TRUE)
          if (class(data_list) %in% "try-error") {
            showNotification("UMAP failed, check parameters", type = "error")
            message("UMAP failed, check parameters")
            return()
          }
        }
        
        if (is.null(data_list$umap_separation$umap_out)) {
          print("UMAP matrix is null, exiting")
          return()
        }
        
        incProgress(amount = 0.5, message = "Clustering", detail = "Performing selected clustering methods")
        
        # Prepare embeddings
        embeddings <- as.matrix(data_list$umap_separation$umap_out)
        
        # Set number of clusters (if applicable)
        n_clusters <- input$k_clustering
        
        # Clustering methods based on user selection
        # Assume 'input$clustering_methods' is a character vector of selected methods
        # For example: input$clustering_methods <- c("kmeans", "hierarchical", "dbscan")
        # K-means Clustering
        if ("kmeans" %in% input$clustering_methods) {
          print("Starting K-means clustering")
          tryCatch({
            set.seed(123)  # For reproducibility
            data_list$kmeans_umap_separation <- kmeans(embeddings, centers = n_clusters)
          }, error = function(e) {
            print(paste("Error in K-means clustering:", e))
          })
        }
        
        # Hierarchical Clustering
        if ("hierarchical" %in% input$clustering_methods) {
          print("Starting hierarchical clustering")
          tryCatch({
            hc <- hclust(dist(embeddings), method = "ward.D2")
            data_list$hierarchical_umap_separation <- cutree(hc, k = n_clusters)
          }, error = function(e) {
            print(paste("Error in Hierarchical clustering:", e))
          })
        }
        
        # DBSCAN Clustering
        if ("dbscan" %in% input$clustering_methods) {
          print("Starting DBSCAN clustering")
          tryCatch({
            data_list$dbscan_umap_separation <- dbscan::dbscan(
              embeddings,
              eps = input$eps,
              minPts = input$minPts
            )
          }, error = function(e) {
            print(paste("Error in DBSCAN clustering:", e))
          })
        }
        
        # HDBSCAN Clustering
        if ("hdbscan" %in% input$clustering_methods) {
          print("Starting HDBSCAN clustering")
          tryCatch({
            data_list$hdbscan_umap_separation <- dbscan::hdbscan(
              embeddings,
              minPts = input$minPts
            )
          }, error = function(e) {
            print(paste("Error in HDBSCAN clustering:", e))
          })
        }
        
        # Spectral Clustering
        if ("spectral" %in% input$clustering_methods) {
          print("Starting spectral clustering")
          tryCatch({
            sigma <- 1
            similarity_matrix <- exp(-as.matrix(dist(embeddings))^2 / (2 * sigma^2))
            n_eigen <- n_clusters
            eig <- RSpectra::eigs(similarity_matrix, k = n_eigen)
            spectral_clusters <- kmeans(eig$vectors, centers = n_clusters)
            data_list$spectral_umap_separation <- spectral_clusters
          }, error = function(e) {
            print(paste("Error in Spectral clustering:", e))
          })
        }
        
        # K-medoids Clustering
        if ("kmedoids" %in% input$clustering_methods) {
          print("Starting K-medoids clustering")
          showNotification("Running K-medoids clustering. This method can be slow on large datasets.", type = "message", duration = 8)
          incProgress(amount = 0, detail = "Running K-medoids (can be slow)")
          tryCatch({
            data_list$kmedoids_umap_separation <- cluster::pam(embeddings, k = n_clusters)
          }, error = function(e) {
            print(paste("Error in K-medoids clustering:", e))
          })
        }
        
        # Fuzzy C-means Clustering
        if ("fuzzy" %in% input$clustering_methods) {
          print("Starting fuzzy C-means clustering")
          tryCatch({
            data_list$fuzzy_umap_separation <- e1071::cmeans(embeddings, centers = n_clusters)
          }, error = function(e) {
            print(paste("Error in Fuzzy C-means clustering:", e))
          })
        }
        
        
        # Model-based Clustering (Mclust)
        if ("mclust" %in% input$clustering_methods) {
          print("Starting model-based clustering (Mclust)")
          if (!requireNamespace("mclust", quietly = TRUE)) {
            showNotification("Skipping Mclust clustering: package 'mclust' is not installed.", type = "warning", duration = 8)
          } else {
            tryCatch({
              data_list$mclust_umap_separation <- mclust::Mclust(embeddings, G = 1:n_clusters)
            }, error = function(e) {
              print(paste("Error in Mclust clustering:", e))
            })
          }
        }
        
        # Self-Organizing Map (SOM) Clustering
        if ("som" %in% input$clustering_methods) {
          print("Starting self-organizing map clustering")
          if (!requireNamespace("kohonen", quietly = TRUE)) {
            showNotification("Skipping SOM clustering: package 'kohonen' is not installed.", type = "warning", duration = 8)
          } else {
            tryCatch({
              scaled_embeddings <- scale(embeddings)
              som_grid <- kohonen::somgrid(xdim = max(2, round(input$k_clustering/2)), 
                                           ydim = max(2, round(input$k_clustering/2)), 
                                           topo = "hexagonal")
              som_model <- kohonen::som(scaled_embeddings, grid = som_grid, rlen = 100)
              data_list$som_umap_separation <- som_model
            }, error = function(e) {
              print(paste("Error in SOM clustering:", e))
            })
          }
        }
        
        # Spherical K-means Clustering
        if ("skmeans" %in% input$clustering_methods) {
          print("Starting spherical K-means clustering")
          tryCatch({
            data_list$skmeans_umap_separation <- skmeans::skmeans(embeddings, k = n_clusters)
          }, error = function(e) {
            print(paste("Error in Spherical K-means clustering:", e))
          })
        }
        
        # Add any additional clustering methods here, following the same pattern
        
        print("Finished clustering methods")
        
        incProgress(amount = 0.9, detail = "Starting quantile color reduction")
        
        print("Started color reduction")
        
        # Color reduction logic (as in your original code)
        if (as.character(Sys.info()['sysname']) == "Windows") {
          # If on Windows, use a single core
          old_nc <- setup_values()[["ncores"]]
          ncores <- 1
        } else {
          ncores <- setup_values()[["ncores"]]
        }
        
        reduced <- unlist(
          parallel::mclapply(
            data_list$umap_separation$color_scheme,
            function(x) color.id(x)[1],
            mc.cores = ncores
          )
        )
        cols_red_list <- unique(names(table(reduced)[table(reduced) > quantile(table(reduced), input$quant)]))
        cols_red_90 <- unlist(
          parallel::mclapply(
            data_list$umap_separation$color_scheme,
            function(x) color.id(x, palette = cols_red_list)[1],
            mc.cores = ncores
          )
        )
        print("Finished color reduction")
        
        if (as.character(Sys.info()['sysname']) == "Windows") {
          # If on Windows, return original ncores
          ncores <- old_nc
        } else {
          ncores <- setup_values()[["ncores"]]
        }
        
        data_list$umap_separation$initial_reduced <- reduced
        data_list$umap_separation$col_reduced <- cols_red_90
        
        # Replace numeric clustering output with text colors to simplify plotting
        
        # Define all R colors
        r_colors <- data.frame(color = colors())
        r_colors <- cbind(r_colors, t(col2rgb(colors())))
        
        pc <- as.data.frame(t(col2rgb(pals::polychrome())))
        col_data <- cbind(Color = pals::polychrome(), pc)
        
        
        # Function to ensure that the max number of colors is not exceeded
        max_col <- function(num_vec_in, color_vec = pc_color_name) {
          a <- num_vec_in
          if (length(unique(a)) / length(color_vec) > 1) {
            ff <- floor(length(unique(a)) / length(color_vec)) + 1
            col_vector <- rep(as.vector(color_vec), ff)
            print("Warning: more than 36 colors, color replication performed")
          } else {
            col_vector <- as.vector(color_vec)
          }
          names(col_vector) <- unique(a)
          return(col_vector)
        }
        
        pc_color_name <- sapply(seq_along(col_data$Color),
                                function(i)
                                  r_colors$color[which.min(
                                    (r_colors$red - col_data$red[i])^2 +
                                      (r_colors$green - col_data$green[i])^2 +
                                      (r_colors$blue - col_data$blue[i])^2
                                  )])
        
        
        # Assign colors to clustering results
        assign_colors <- function(clusters) {
          if (is.null(clusters) || length(clusters) == 0) {
            return(NULL)
          }
          cols <- max_col(clusters, pc_color_name)
          colored_clusters <- cols[as.character(clusters)]
          return(colored_clusters)
        }
        
        # For each clustering method, assign colors
        if ("kmeans" %in% input$clustering_methods && !is.null(data_list$kmeans_umap_separation$cluster)) {
          data_list$kmeans_umap_separation$cluster <- assign_colors(data_list$kmeans_umap_separation$cluster)
        }
        
        if ("hierarchical" %in% input$clustering_methods && !is.null(data_list$hierarchical_umap_separation)) {
          data_list$hierarchical_umap_separation <- assign_colors(data_list$hierarchical_umap_separation)
        }
        
        if ("dbscan" %in% input$clustering_methods && !is.null(data_list$dbscan_umap_separation$cluster)) {
          data_list$dbscan_umap_separation$cluster <- assign_colors(data_list$dbscan_umap_separation$cluster)
        }
        
        if ("hdbscan" %in% input$clustering_methods && !is.null(data_list$hdbscan_umap_separation$cluster)) {
          data_list$hdbscan_umap_separation$cluster <- assign_colors(data_list$hdbscan_umap_separation$cluster)
        }
        
        if ("spectral" %in% input$clustering_methods && !is.null(data_list$spectral_umap_separation$cluster)) {
          data_list$spectral_umap_separation$cluster <- assign_colors(data_list$spectral_umap_separation$cluster)
        }
        
        if ("kmedoids" %in% input$clustering_methods && !is.null(data_list$kmedoids_umap_separation$clustering)) {
          data_list$kmedoids_umap_separation$clustering <- assign_colors(data_list$kmedoids_umap_separation$clustering)
        }
        
        if ("fuzzy" %in% input$clustering_methods && !is.null(data_list$fuzzy_umap_separation$membership)) {
          # For fuzzy clustering, get hard assignments
          fuzzy_clusters <- apply(data_list$fuzzy_umap_separation$membership, 1, which.max)
          data_list$fuzzy_umap_separation$cluster <- assign_colors(fuzzy_clusters)
        }
        
        if ("mclust" %in% input$clustering_methods && !is.null(data_list$mclust_umap_separation$classification)) {
          data_list$mclust_umap_separation$classification <- assign_colors(data_list$mclust_umap_separation$classification)
        }
        
        if ("som" %in% input$clustering_methods && !is.null(data_list$som_umap_separation$unit.classif)) {
          data_list$som_umap_separation$unit.classif <- assign_colors(data_list$som_umap_separation$unit.classif)
        }
        
        if ("skmeans" %in% input$clustering_methods && !is.null(data_list$skmeans_umap_separation$cluster)) {
          data_list$skmeans_umap_separation$cluster <- assign_colors(data_list$skmeans_umap_separation$cluster)
        }
        
        
        
        
        
        x2$data_list <- data_list
      })
      print(proc.time() - ptm)
    })
    
    #     observe({
    #       x2$data_list <- data_list_shiny()
    #     })
    
    #run clustering only
    observeEvent(input$action_dbscan, {
      x2 <- preproc_values()[["x2"]]
      if (is.null(x2$data_list)) {
        showNotification("No UMAP embedding found. Run 'Start UMAP clustering' first, then re-run clustering.", type = "warning", duration = 7)
        return()
      }
      if (is.null(x2$data_list$umap_separation) || is.null(x2$data_list$umap_separation$umap_out)) {
        showNotification("UMAP embedding is missing from current session. Re-run UMAP first.", type = "warning", duration = 7)
        return()
      }
      
      data_list <- x2$data_list
      
      print("Re-running selected clustering methods")
      
      withProgress(message = "Re-running clustering on existing UMAP", {
        embeddings <- as.matrix(data_list$umap_separation$umap_out)
        
        # Set the number of cores for parallel processing
        if (as.character(Sys.info()['sysname']) == "Windows") {
          # If on Windows, use a single core
          old_nc <- setup_values()[["ncores"]]
          ncores <- 1
        } else {
          ncores <- setup_values()[["ncores"]]
        }
        
        # Re-run clustering methods selected by the user
        # K-means Clustering
        if ("kmeans" %in% input$clustering_methods) {
          print("Re-running K-means clustering")
          tryCatch({
            n_clusters <- input$k_clustering  # Ensure 'input$k_clustering' is defined in your UI
            set.seed(123)
            data_list$kmeans_umap_separation <- kmeans(embeddings, centers = n_clusters)
          }, error = function(e) {
            print(paste("Error in K-means clustering:", e))
          })
        }
        
        # Hierarchical Clustering
        if ("hierarchical" %in% input$clustering_methods) {
          print("Re-running hierarchical clustering")
          tryCatch({
            hc <- hclust(dist(embeddings), method = "ward.D2")
            n_clusters <- input$k_clustering
            data_list$hierarchical_umap_separation <- cutree(hc, k = n_clusters)
          }, error = function(e) {
            print(paste("Error in Hierarchical clustering:", e))
          })
        }
        
        # DBSCAN Clustering
        if ("dbscan" %in% input$clustering_methods) {
          print("Re-running DBSCAN")
          tryCatch({
            data_list$dbscan_umap_separation <- dbscan::dbscan(
              embeddings,
              eps = input$eps,
              minPts = input$minPts
            )
          }, error = function(e) {
            print(paste("Error in DBSCAN clustering:", e))
          })
        }
        
        # HDBSCAN Clustering
        if ("hdbscan" %in% input$clustering_methods) {
          print("Re-running HDBSCAN")
          tryCatch({
            data_list$hdbscan_umap_separation <- dbscan::hdbscan(
              embeddings,
              minPts = input$minPts
            )
          }, error = function(e) {
            print(paste("Error in HDBSCAN clustering:", e))
          })
        }
        
        # Spectral Clustering
        if ("spectral" %in% input$clustering_methods) {
          print("Re-running spectral clustering")
          tryCatch({
            sigma <- 1
            similarity_matrix <- exp(-as.matrix(dist(embeddings))^2 / (2 * sigma^2))
            n_eigen <- input$k_clustering
            eig <- RSpectra::eigs(similarity_matrix, k = n_eigen)
            spectral_clusters <- kmeans(eig$vectors, centers = n_eigen)
            data_list$spectral_umap_separation <- spectral_clusters
          }, error = function(e) {
            print(paste("Error in Spectral clustering:", e))
          })
        }
        
        # K-medoids Clustering
        if ("kmedoids" %in% input$clustering_methods) {
          print("Re-running K-medoids clustering")
          showNotification("Re-running K-medoids clustering. This method can be slow on large datasets.", type = "message", duration = 8)
          incProgress(amount = 0, detail = "Running K-medoids (can be slow)")
          tryCatch({
            data_list$kmedoids_umap_separation <- cluster::pam(embeddings, k = input$k_clustering)
          }, error = function(e) {
            print(paste("Error in K-medoids clustering:", e))
          })
        }
        
        # Fuzzy C-means Clustering
        if ("fuzzy" %in% input$clustering_methods) {
          print("Re-running fuzzy C-means clustering")
          tryCatch({
            data_list$fuzzy_umap_separation <- e1071::cmeans(embeddings, centers = input$k_clustering)
          }, error = function(e) {
            print(paste("Error in Fuzzy C-means clustering:", e))
          })
        }
        
        # Model-based Clustering (Mclust)
        if ("mclust" %in% input$clustering_methods) {
          print("Re-running model-based clustering (Mclust)")
          if (!requireNamespace("mclust", quietly = TRUE)) {
            showNotification("Skipping Mclust clustering: package 'mclust' is not installed.", type = "warning", duration = 8)
          } else {
            tryCatch({
              data_list$mclust_umap_separation <- mclust::Mclust(embeddings, G = 1:input$k_clustering)
            }, error = function(e) {
              print(paste("Error in Mclust clustering:", e))
            })
          }
        }
        
        # Self-Organizing Map (SOM) Clustering
        if ("som" %in% input$clustering_methods) {
          print("Re-running self-organizing map clustering")
          if (!requireNamespace("kohonen", quietly = TRUE)) {
            showNotification("Skipping SOM clustering: package 'kohonen' is not installed.", type = "warning", duration = 8)
          } else {
            tryCatch({
              scaled_embeddings <- scale(embeddings)
              som_grid <- kohonen::somgrid(xdim = max(2, round(input$k_clustering/2)), 
                                           ydim = max(2, round(input$k_clustering/2)), 
                                           topo = "hexagonal")
              som_model <- kohonen::som(scaled_embeddings, grid = som_grid, rlen = 100)
              data_list$som_umap_separation <- som_model
            }, error = function(e) {
              print(paste("Error in SOM clustering:", e))
            })
          }
        }
        
        # Spherical K-means Clustering
        if ("skmeans" %in% input$clustering_methods) {
          print("Re-running spherical K-means clustering")
          tryCatch({
            data_list$skmeans_umap_separation <- skmeans::skmeans(embeddings, k = input$k_clustering)
          }, error = function(e) {
            print(paste("Error in Spherical K-means clustering:", e))
          })
        }
        
        # Color reduction logic (if needed)
        reduced <- unlist(
          parallel::mclapply(
            data_list$umap_separation$color_scheme,
            function(x) color.id(x)[1],
            mc.cores = ncores
          )
        )
        cols_red_list <- unique(names(table(reduced)[table(reduced) > quantile(table(reduced), input$quant)]))
        cols_red_90 <- unlist(
          parallel::mclapply(
            data_list$umap_separation$color_scheme,
            function(x) color.id(x, palette = cols_red_list)[1],
            mc.cores = ncores
          )
        )
        
        if (as.character(Sys.info()['sysname']) == "Windows") {
          # If on Windows, restore the original number of cores
          ncores <- old_nc
        }
        
        data_list$umap_separation$initial_reduced <- reduced
        data_list$umap_separation$col_reduced <- cols_red_90
        
        # Replace numeric clustering output with text colors to simplify plotting
        
        # Define all R colors
        r_colors <- data.frame(color = colors())
        r_colors <- cbind(r_colors, t(col2rgb(colors())))
        
        pc <- as.data.frame(t(col2rgb(pals::polychrome())))
        col_data <- cbind(Color = pals::polychrome(), pc)
        
        pc_color_name <- sapply(seq_along(col_data$Color),
                                function(i)
                                  r_colors$color[which.min(
                                    (r_colors$red - col_data$red[i])^2 +
                                      (r_colors$green - col_data$green[i])^2 +
                                      (r_colors$blue - col_data$blue[i])^2
                                  )])
        
        # Function to ensure that the max number of colors is not exceeded
        max_col <- function(num_vec_in, color_vec = pc_color_name) {
          a <- num_vec_in
          if (length(unique(a)) / length(color_vec) > 1) {
            ff <- floor(length(unique(a)) / length(color_vec)) + 1
            col_vector <- rep(as.vector(color_vec), ff)
            print("Warning: more than 36 colors, color replication performed")
          } else {
            col_vector <- as.vector(color_vec)
          }
          names(col_vector) <- unique(a)
          return(col_vector)
        }
        
        # Assign colors to clustering results
        assign_colors <- function(clusters) {
          if (is.null(clusters) || length(clusters) == 0) {
            return(NULL)
          }
          cols <- max_col(clusters, pc_color_name)
          colored_clusters <- cols[as.character(clusters)]
          return(colored_clusters)
        }
        
        # For each clustering method, assign colors
        if ("kmeans" %in% input$clustering_methods && !is.null(data_list$kmeans_umap_separation$cluster)) {
          data_list$kmeans_umap_separation$cluster <- assign_colors(data_list$kmeans_umap_separation$cluster)
        }
        
        if ("hierarchical" %in% input$clustering_methods && !is.null(data_list$hierarchical_umap_separation)) {
          data_list$hierarchical_umap_separation <- assign_colors(data_list$hierarchical_umap_separation)
        }
        
        if ("dbscan" %in% input$clustering_methods && !is.null(data_list$dbscan_umap_separation$cluster)) {
          data_list$dbscan_umap_separation$cluster <- assign_colors(data_list$dbscan_umap_separation$cluster)
        }
        
        if ("hdbscan" %in% input$clustering_methods && !is.null(data_list$hdbscan_umap_separation$cluster)) {
          data_list$hdbscan_umap_separation$cluster <- assign_colors(data_list$hdbscan_umap_separation$cluster)
        }
        
        if ("spectral" %in% input$clustering_methods && !is.null(data_list$spectral_umap_separation$cluster)) {
          data_list$spectral_umap_separation$cluster <- assign_colors(data_list$spectral_umap_separation$cluster)
        }
        
        if ("kmedoids" %in% input$clustering_methods && !is.null(data_list$kmedoids_umap_separation$clustering)) {
          data_list$kmedoids_umap_separation$clustering <- assign_colors(data_list$kmedoids_umap_separation$clustering)
        }
        
        if ("fuzzy" %in% input$clustering_methods && !is.null(data_list$fuzzy_umap_separation$membership)) {
          # For fuzzy clustering, get hard assignments
          fuzzy_clusters <- apply(data_list$fuzzy_umap_separation$membership, 1, which.max)
          data_list$fuzzy_umap_separation$cluster <- assign_colors(fuzzy_clusters)
        }
        
        if ("mclust" %in% input$clustering_methods && !is.null(data_list$mclust_umap_separation$classification)) {
          data_list$mclust_umap_separation$classification <- assign_colors(data_list$mclust_umap_separation$classification)
        }
        
        if ("som" %in% input$clustering_methods && !is.null(data_list$som_umap_separation$unit.classif)) {
          data_list$som_umap_separation$unit.classif <- assign_colors(data_list$som_umap_separation$unit.classif)
        }
        
        if ("skmeans" %in% input$clustering_methods && !is.null(data_list$skmeans_umap_separation$cluster)) {
          data_list$skmeans_umap_separation$cluster <- assign_colors(data_list$skmeans_umap_separation$cluster)
        }
        
        print("Finished re-running selected clustering methods")
        
        x2$data_list <- data_list
      })
    })
    
    
    
    #clust<- reactive({
    
    observe({
      x2 <- preproc_values()[["x2"]]
      if (is.null(x2$data_list)) {
        return()
      }
      x2$bkcols <- switch(
        input$umap_cols,
        # Get color info
        #"RGB" = data_list_shiny()$umap_separation$color_scheme,
        "Reduced" = x2$data_list$umap_separation$initial_reduced,
        "Reduced2" = x2$data_list$umap_separation$col_reduced,
        "kmeans" = x2$data_list$kmeans_umap_separation$cluster,
        "hierarchical" = x2$data_list$hierarchical_umap_separation,
        "dbscan" = x2$data_list$dbscan_umap_separation$cluster,
        "hdbscan" = x2$data_list$hdbscan_umap_separation$cluster,
        "spectral" = x2$data_list$spectral_umap_separation$cluster,
        "kmedoids" = x2$data_list$kmedoids_umap_separation$clustering,
        "fuzzy" = {
          if (!is.null(x2$data_list$fuzzy_umap_separation$cluster) &&
              length(x2$data_list$fuzzy_umap_separation$cluster) > 0) {
            x2$data_list$fuzzy_umap_separation$cluster
          } else if (!is.null(x2$data_list$fuzzy_umap_separation$membership)) {
            # Fallback to hard assignments if colorized labels are not present.
            apply(x2$data_list$fuzzy_umap_separation$membership, 1, which.max)
          } else {
            rep("grey70", nrow(x2$data_list$umap_separation$umap_out))
          }
        },
        "mclust" = x2$data_list$mclust_umap_separation$classification,
        "som" = x2$data_list$som_umap_separation$unit.classif,
        "skmeans" = x2$data_list$skmeans_umap_separation$cluster
      )
      if (!is.null(x2$bkcols)) {
        x2$bkcols <- normalize_labels(x2$bkcols)
      }
    })
    
    output$data_list_table <- renderPrint({
      x2 <- preproc_values()[["x2"]]
      paste(table(x2$bkcols), table(unique(x2$tf_list)))
    })
    
    
    #UMAP embeddings pairwise plot
    output$plot6 <- renderImage({
      x2 <- preproc_values()[["x2"]]
      req(x2$data_list)
      req(x2$tf_list)
      print("plot6 continuing")
      # A temp file to save the output.
      # This file will be removed later by renderImage
      outfile <- tempfile(fileext = '.png')
      
      png(outfile, width = 800, height = 600)
      
      
      if (input$seg_choice == "anat_seg") {
        tf_list = x2$tf_list_anat
      } else {
        tf_list <- x2$tf_list
      }
      
      tf_list <- sanitize_tf_mask(tf_list, nrow(x2$data_list$umap_separation$umap_out))
      if (is.null(tf_list) || sum(tf_list) == 0) {
        draw_empty_plot("No pixels selected for UMAP display.")
        dev.off()
        return(list(
          src = outfile,
          contentType = 'image/png',
          width = 800,
          height = 600,
          alt = "No pixels selected"
        ))
      }
      
      if (!is.null(x2$bkcols) && length(x2$bkcols) == length(tf_list)) {
        cols <- normalize_labels(x2$bkcols[tf_list])
      } else {
        cols <- rep("grey70", sum(tf_list))
      }
      emb_fill <- build_fill_mapping(cols)
      cols_plot <- unname(emb_fill$values[emb_fill$labels])
      
      updateSelectizeInput(session,
                           'Color_choices',
                           choices = unique(emb_fill$labels),
                           server = TRUE)
      
      print(
        pairs(
          x2$data_list$umap_separation$umap_out[tf_list,],
          col = cols_plot,
          pch = ".",
          cex = input$cex,
          main = runNames(x2$overview_peaks_sel)
        )
      )
      
      #
      #
      # }
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
    
    
    output$segmentation_plots <- renderUI({
      if (input$seg_choice == "bk_seg") {
        list(
          imageOutput(ns("plot7"), width = "800px", height = "600px"),
          imageOutput(ns("plot6"), width = "800px", height = "600px"),
          imageOutput(ns("colormaps"), width = "800px", height = "300px")
        )
      } else if (input$seg_choice == "anat_seg") {
        list(
          imageOutput(
            ns("plot7_tissue_umap"),
            width = "800px",
            height = "600px"
          ),
          imageOutput(ns("colormaps"), width = "800px", height = "300px"),
          imageOutput(ns("plot7"), width = "800px", height = "600px"),
          imageOutput(ns("plot6"), width = "800px", height = "600px")
          
        )
      }
    })
    
    
    
    output$plot7_tissue_umap <- renderImage({
      
      
      req(input$mytable_rows_selected)
      req(input$seg_pdat_col)
      x2 <- preproc_values()[["x2"]]
      
      x2$rcol_plot <- input$seg_pdat_col
      
      req(x2$rcol_plot)
      
      print("plot7_tissue_umap_continuting")
      outfile <- tempfile(fileext = '.png')
      
      png(outfile, width = 800, height = 600)
      if (is.null(input$show_dat) || length(input$show_dat) == 0) {
        draw_empty_plot("No pData values selected to display.")
        dev.off()
        return(list(
          src = outfile,
          contentType = 'image/png',
          width = 800,
          height = 600,
          alt = "No pData values selected"
        ))
      }
      
      
      
      
      if (input$seg_choice == "anat_seg") {
        #test for existing proc data w/o umap, then test for umap to visualize
        if (!is.null(x2$list_proc_img) && is.null(x2$umap_name)) {
          img.dat <- combine_card(x2$list_proc_img[input$mytable_rows_selected])
        } else if (!is.null(x2$tf_list_umap) &&
                   !is.null(x2$umap_name)) {
          if (identical(x2$umap_name, runNames(x2$mytable_selected))) {
            img.dat <- safe_subset_pixels(x2$mytable_selected, x2$tf_list_umap, "UMAP pData plot", notify = FALSE)
            if (!is.null(x2$pdat_anat)) {
              tf_umap <- sanitize_tf_mask(x2$tf_list_umap, nrow(x2$pdat_anat))
              if (!is.null(tf_umap) && !is.null(img.dat) && ncol(img.dat) == sum(tf_umap)) {
                pData(img.dat) <- x2$pdat_anat[tf_umap,]
              }
            }
          }
        } else {
          #test to see if x2$pdat_anat matches and use that for viz
          img.dat <-
            x2$mytable_selected #%>% subsetPixels(x2$tf_list)
          if (!is.null(x2$pdat_anat) &&
              nrow(x2$pdat_anat) == ncol(img.dat)) {
            pData(img.dat) <- x2$pdat_anat
          }
        }
      } else {
        img.dat <- safe_subset_pixels(x2$mytable_selected, x2$tf_list, "Segmentation image", notify = FALSE)
      }
      
      
      
      
      #img.dat<-x2$mytable_selected #x2$overview_peaks_sel  #subset somehow??
      
      #browser()
      
      #if(x2$rcol_plot=="Rcol_reduced"){
      if (is.null(img.dat) || ncol(img.dat) == 0) {
        draw_empty_plot("No image data available for this selection.")
        dev.off()
        return(list(
          src = outfile,
          contentType = 'image/png',
          width = 800,
          height = 600,
          alt = "No image data available"
        ))
      }
      
      pdat_df <- as.data.frame(pData(img.dat))
      if (!x2$rcol_plot %in% colnames(pdat_df)) {
        draw_empty_plot(sprintf("Field '%s' not found in pData.", x2$rcol_plot))
        dev.off()
        return(list(
          src = outfile,
          contentType = 'image/png',
          width = 800,
          height = 600,
          alt = "pData field not found"
        ))
      }
      cols <- normalize_labels(pdat_df[, x2$rcol_plot, drop = TRUE])
      show_vals <- normalize_labels(input$show_dat)
      tf_list <- cols %in% show_vals
      tf_list[is.na(tf_list)] <- FALSE
      
      img.dat <- safe_subset_pixels(img.dat, tf_list, "pData display filter", notify = FALSE)
      if (is.null(img.dat) || ncol(img.dat) == 0) {
        draw_empty_plot("No pixels left after applying display filters.")
        dev.off()
        return(list(
          src = outfile,
          contentType = 'image/png',
          width = 800,
          height = 600,
          alt = "No pixels left after filtering"
        ))
      }
      
      x <-
        coord(img.dat)[, 1]
      y <- coord(img.dat)[, 2]
      runs = Cardinal::run(img.dat)
      aa <- data.frame(x, y, runs)
      cols <- normalize_labels((as.data.frame(pData(img.dat))[, x2$rcol_plot, drop = TRUE]))
      
      # } else {
      #   
      #   cols=(as.factor(as.data.frame(pData(img.dat))[,x2$rcol_plot]))
      #   tf_list<-cols%in%input$show_dat
      #
      #   img.dat <- img.dat %>% subsetPixels(tf_list)
      #
      #   x<-coord(img.dat)[,1] ; y<-coord(img.dat)[,2]; runs=Cardinal::run(img.dat)
      #   aa<-data.frame(x,y,runs=req(runs))
      #   cols=as.numeric(as.factor(as.data.frame(pData(img.dat))[,x2$rcol_plot]))+1
      #
      #
      #
      # }
      
      
      
      if (dim(aa)[1] == 0) {
        draw_empty_plot("No pixels available to plot.")
        dev.off()
        return(list(
          src = outfile,
          contentType = 'image/png',
          width = 800,
          height = 600,
          alt = "No pixels available"
        ))
      }
      
      fill_map <- build_fill_mapping(cols)
      fill_factor <- factor(fill_map$labels, levels = names(fill_map$values))
      log_fill_mapping("plot7_tissue_umap", fill_map$values)
      p1 <- ggplot2::ggplot(aa, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_tile(ggplot2::aes(fill = fill_factor)) +
        ggplot2::facet_wrap(. ~ runs) +
        ggplot2::theme_minimal() +
        ggplot2::scale_y_continuous(trans = "reverse") +
        ggplot2::theme(legend.position = "bottom") +
        ggplot2::scale_fill_manual(values = fill_map$values, name = x2$rcol_plot, drop = FALSE)
      
      print(p1)
      
      dev.off()
      graphics.off()
      
      # Return a list containing the filename
      list(
        src = outfile,
        contentType = 'image/png',
        width = 800,
        height = 600,
        alt = "This is alternate text"
      )
    }, deleteFile = TRUE)
    
    
    
    #plot to make sure it looks okay
    output$plot7 <- renderImage({
      x2 <- preproc_values()[["x2"]]
      req(x2$mytable_selected)
      req(x2$tf_list)
      
      outfile <- tempfile(fileext = '.png')
      
      png(outfile, width = 800, height = 600)
      
      
      if (input$seg_choice == "anat_seg") {
        tf_list = x2$tf_list_anat
      } else {
        tf_list <- x2$tf_list
      }
      
      
      
      #require(ggplot2)
      
      req(x2$data_list)
      req(x2$bkcols)
      
      # A temp file to save the output.
      # This file will be removed later by renderImage
      #outfile <- tempfile(fileext = '.png')
      
      #png(outfile, width = 800, height = 600)
      
      
      img.dat <-
        safe_subset_pixels(x2$data_list$clean_img.dat, tf_list, "segmentation map", notify = FALSE)
      if (is.null(img.dat) || ncol(img.dat) == 0) {
        draw_empty_plot("No image data left after filtering.")
        dev.off()
        return(list(
          src = outfile,
          contentType = 'image/png',
          width = 800,
          height = 600,
          alt = "No image data"
        ))
      }
      #previous line will revert to original image after saving to x2$overview_peaks_sel
      #img.dat<-x2$mytable_selected %>% subsetPixels(x2$tf_list)
      
      x <-
        coord(img.dat)[, 1]
      y <- coord(img.dat)[, 2]
      runs = Cardinal::run(img.dat)
      req(runs)
      a <- data.frame(x, y, runs)
      if (!is.null(x2$bkcols) && length(x2$bkcols) == length(tf_list)) {
        cols <- normalize_labels(x2$bkcols[tf_list])
      } else {
        cols <- rep("grey70", nrow(a))
      }
      
      cat(names(table(cols)))
      
      
      
      
      
      #if(input$umap_cols=="Reduced2" | input$umap_cols=="Reduced" ){
      fill_map <- build_fill_mapping(cols)
      fill_factor <- factor(fill_map$labels, levels = names(fill_map$values))
      log_fill_mapping("plot7", fill_map$values)
      p1 <- ggplot2::ggplot(a, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_tile(ggplot2::aes(fill = fill_factor)) +
        ggplot2::facet_wrap(. ~ runs) +
        ggplot2::theme_minimal() +
        ggplot2::scale_y_continuous(trans = "reverse") +
        ggplot2::theme(legend.position = "bottom") +
        ggplot2::scale_fill_manual(values = fill_map$values, name = input$umap_cols, drop = FALSE)
      # } else { #for h/dbscan
      #   if(length(table(cols))/length(pals::polychrome()) >1) {
      #     ff<-floor(length(table(cols))/length(pals::polychrome()))+1
      #     col_vector=rep(as.vector(pals::polychrome()), ff)
      #   }else {
      #     col_vector=as.vector(pals::polychrome())
      #   }
      #
      #   #h/dbscan plot
      #   p1<-ggplot2::ggplot(a, ggplot2::aes(x=x, y=y))+
      #     ggplot2::geom_tile(ggplot2::aes( fill=as.factor(cols)))+
      #     ggplot2::facet_wrap(.~runs)+
      #     ggplot2::theme_minimal()+
      #     ggplot2::scale_y_continuous(trans = "reverse")+
      #     ggplot2::scale_fill_manual(values=col_vector) +
      #     ggplot2::guides(colour=FALSE, fill=FALSE)
      # }
      #
      
      print(p1)
      
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
    
    #lets have the user choose which values to keep from clust, then regenerate the dataset based on those values.
    
    output$colormaps <- renderImage({
      if (input$seg_choice == "anat_seg") {
        #only output one colormap for anatomical segmentatin
        
        outfile <- tempfile(fileext = '.png')
        
        png(outfile, width = 800, height = 600)
        
        x2 <- preproc_values()[["x2"]]
        req(x2$mytable_selected)
        req(input$seg_pdat_col)
        req(input$show_dat)
        #browser()
        
        #TODO - deal with segmentation from exisitng info, without running UMAP
        img.dat <- x2$mytable_selected
        if (!is.null(x2$pdat_anat) &&
            identical(x2$umap_name, runNames(x2$mytable_selected))) {
          pData(img.dat) <- x2$pdat_anat
        } else if (!is.null(x2$pdat_anat) &&
                   dim(x2$pdat_anat)[1] == ncol(img.dat)) {
          pData(img.dat) <- x2$pdat_anat
        }
        dat <- as.data.frame(pData(img.dat))
        
        cols <-
          dat[, input$seg_pdat_col][dat[, input$seg_pdat_col] %in% input$show_dat]
        a <- as.data.frame(table(cols))
        
        if (input$seg_pdat_col != "Rcol_reduced") {
          a$cols2 <-
            as.numeric(as.factor(as.data.frame(table(cols))$cols)) + 1
        } else {
          a$cols2 <- a$cols
        }
        
        
        
        req(a$cols2)
        
        p1 <-
          ggplot2::ggplot(a, ggplot2::aes(x = cols2, y = log(Freq))) +
          ggplot2::geom_col(fill = a$cols2) +
          ggplot2::geom_text(
            ggplot2::aes(label = cols),
            position = ggplot2::position_dodge(width = 0.9),
            angle = 90,
            size = 10,
            hjust = 1
          ) +
          ggplot2::theme_minimal() +
          ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 1)) +
          ggplot2::ggtitle(input$show_dat)
        #browser() TODO figure out plotting errors during segmentation
        
        print(p1)
        
        
      } else {
        x2 <- preproc_values()[["x2"]]
        req(x2$data_list)
        # A temp file to save the output.
        # This file will be removed later by renderImage
        outfile <- tempfile(fileext = '.png')
        
        png(outfile, width = 800, height = 600)
        
        
        #   
        #   
        #   col_red <- x2$data_list$umap_separation$col_reduced
        #   col_dbscan <- x2$data_list$dbscan_umap_separation$cluster
        #   col_hdbscan <- x2$data_list$hdbscan_umap_separation$cluster
        #   
        #   #helper function to generate color map plots; a is the color vector
        #   plot_col_map <- function(a, titles = "plot") {
        #     #require(ggplot2)
        #     
        #     b <- data.frame(table(a))
        #     colnames(b) <- c("color",  "Freq")
        #     #if(input$umap_cols=="Reduced2" | input$umap_cols=="Reduced"){
        #     
        #     p1 <-
        #       ggplot2::ggplot(b, ggplot2::aes(x = color, y = log(Freq))) +
        #       ggplot2::geom_col(fill = b$color) +
        #       ggplot2::geom_text(
        #         ggplot2::aes(label = color),
        #         position = ggplot2::position_dodge(width = 0.9),
        #         angle = 90,
        #         size = 10,
        #         hjust = 1
        #       ) +
        #       ggplot2::theme_minimal() +
        #       ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2)) +
        #       ggplot2::ggtitle(titles)
        #     # } else {
        #     #
        #     #   browser
        #     #   if(length(table(a))/length(pals::polychrome()) >1) {
        #     #     ff<-floor(length(table(a))/length(pals::polychrome()))+1
        #     #     col_vector=rep(as.vector(pals::polychrome()), ff)
        #     #     print("Warning: more than 36 colors, color replication performed")
        #     #   }else {
        #     #     col_vector=as.vector(pals::polychrome())
        #     #   }
        #     #
        #     #   p1<-ggplot2::ggplot(b, ggplot2::aes(x=color,y=log(Freq)))+
        #     #     ggplot2::geom_col(ggplot2::aes(fill=color))+
        #     #     ggplot2::geom_text(ggplot2::aes(label=color), position=ggplot2::position_dodge(width=0.9), angle=90, size=10)+
        #     #     ggplot2::theme_minimal()+
        #     #     ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2))+
        #     #     ggplot2::ggtitle(titles)+
        #     #     ggplot2::scale_fill_manual(values=col_vector)+
        #     #     ggplot2::guides(fill="none")
        #     # }
        #     
        #     return(p1)
        #   }
        #   
        #   p1 <- plot_col_map(col_red, titles = "Reduced R colors")
        #   p2 <- plot_col_map(col_dbscan, titles = "DBSCAN colors")
        #   
        #   
        #   if (sum(is.na(x2$data_list$hdbscan_umap_separation$cluster)) <
        #       1) {
        #     p3 <- plot_col_map(col_hdbscan, titles = "HDBSCAN colors")
        #     print(gridExtra::grid.arrange(p1, p2, p3, nrow = 1))
        #   } else {
        #     print(gridExtra::grid.arrange(p1, p2, nrow = 1))
        #   }
        # }
        # 
        # 
        
        #    x2 <- preproc_values()[["x2"]]
        #    req(x2$data_list)
        # #New plottign code:
        #Helper function to generate color map plots; 'a' is the color vector
        plot_col_map <- function(a, titles = "plot") {
          if (is.null(a)) {
            return(NULL)
          }
          if (is.matrix(a) || is.data.frame(a)) a <- as.vector(a)
          if (is.list(a)) a <- unlist(a, use.names = FALSE)
          a <- as.character(a)
          a <- a[!is.na(a) & nzchar(a)]
          if (length(a) < 2 || length(unique(a)) < 2) {
            # If 'a' is empty or has only one unique value, skip this plot.
            return(NULL)
          }
          
          freq_tab <- sort(table(a), decreasing = TRUE)
          b <- data.frame(
            color = names(freq_tab),
            Freq = as.numeric(freq_tab),
            stringsAsFactors = FALSE
          )
          
          p1 <- ggplot2::ggplot(b, ggplot2::aes(x = color, y = log(Freq))) +
            ggplot2::geom_col(fill = b$color) +
            ggplot2::geom_text(
              ggplot2::aes(label = color),
              position = ggplot2::position_dodge(width = 0.9),
              angle = 90,
              size = 4,
              hjust = 1
            ) +
            ggplot2::theme_minimal() +
            ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2)) +
            ggplot2::ggtitle(titles) +
            ggplot2::theme(axis.title.x = ggplot2::element_blank())
          
          return(p1)
        }
        
        
        # Initialize a list to hold the plots
        plot_list <- list()
        
        # First, add the Reduced R colors plot
        col_red <- x2$data_list$umap_separation$col_reduced
        p_red <- plot_col_map(col_red, titles = "Reduced R colors")
        if (!is.null(p_red)) {
          plot_list <- c(plot_list, list(p_red))
        }
        
        # Define a list of clustering methods and their corresponding data paths and titles
        clustering_methods <- list(
          "dbscan" = list(
            data = x2$data_list$dbscan_umap_separation$cluster,
            title = "DBSCAN Clusters"
          ),
          "hdbscan" = list(
            data = x2$data_list$hdbscan_umap_separation$cluster,
            title = "HDBSCAN Clusters"
          ),
          "kmeans" = list(
            data = if (!is.null(x2$data_list$kmeans_umap_separation)) x2$data_list$kmeans_umap_separation$cluster else NULL,
            title = "K-means Clusters"
          ),
          "hierarchical" = list(
            data = if (!is.null(x2$data_list$hierarchical_umap_separation)) x2$data_list$hierarchical_umap_separation else NULL,
            title = "Hierarchical Clusters"
          ),
          "spectral" = list(
            data = if (!is.null(x2$data_list$spectral_umap_separation)) x2$data_list$spectral_umap_separation$cluster else NULL,
            title = "Spectral Clusters"
          ),
          "kmedoids" = list(
            data = if (!is.null(x2$data_list$kmedoids_umap_separation)) x2$data_list$kmedoids_umap_separation$clustering else NULL,
            title = "K-medoids Clusters"
          ),
          "fuzzy" = list(
            data = if (!is.null(x2$data_list$fuzzy_umap_separation$cluster)) {
              x2$data_list$fuzzy_umap_separation$cluster
            } else if (!is.null(x2$data_list$fuzzy_umap_separation$membership)) {
              apply(x2$data_list$fuzzy_umap_separation$membership, 1, which.max)
            } else NULL,
            title = "Fuzzy C-means Clusters"
          ),
          "mclust" = list(
            data = if (!is.null(x2$data_list$mclust_umap_separation)) x2$data_list$mclust_umap_separation$classification else NULL,
            title = "Mclust Clusters"
          ),
          "som" = list(
            data = if (!is.null(x2$data_list$som_umap_separation)) x2$data_list$som_umap_separation$unit.classif else NULL,
            title = "Self-Organizing Map Clusters"
          ),
          "skmeans" = list(
            data = if (!is.null(x2$data_list$skmeans_umap_separation)) x2$data_list$skmeans_umap_separation$cluster else NULL,
            title = "Spherical K-means Clusters"
          )
        )
        
        
        # Loop over the clustering methods and create plots for those that have been run
        for (method_name in names(clustering_methods)) {
          method <- clustering_methods[[method_name]]
          col_cluster <- method$data
          if (!is.null(col_cluster) && length(col_cluster) > 0) {
            # Create the plot only if the result is valid (not NULL)
            p <- plot_col_map(col_cluster, titles = method$title)
            
            if (!is.null(p)) {
              plot_list <- c(plot_list, list(p))  # Add plot only if it's valid
            }
          }
        }
        
        # Arrange the plots in a grid
        n_plots <- length(plot_list)
        if (n_plots == 0) {
          draw_empty_plot("No clustering color maps are available yet.")
        } else {
          n_cols <- 3  # Number of columns in the grid
          print(do.call(gridExtra::grid.arrange, c(plot_list, ncol = n_cols)))
        }
      }
      
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
    
    #color choices for umap selection
    output$Color_choices <-
      renderUI({
        #https://gist.github.com/wch/4211337
        # If missing input, return to avoid error later in function
        #if(is.null(x2$bkcols))
        #   return()
        
        # Get the data set with the appropriate name
        x2 <- preproc_values()[["x2"]]
        dat <- (x2$bkcols)
        cols <- unique(dat)
        
        # Create the checkboxes and select them all by default
        checkboxGroupInput(ns("cols"),
                           "Choose colors",
                           choices  = cols,
                           selected = cols)
      })
    
    
    
    # adjust T/F list from pData value filter (anat mode)
    observeEvent({
      input$seg_choice
      input$seg_pdat_col
      input$show_dat
      1
    }, {
      x2 <- preproc_values()[["x2"]]
      req(x2$mytable_selected)
      if (!identical(input$seg_choice, "anat_seg")) {
        return()
      }
      
      dat <- get_active_pdata_df(x2)
      if (is.null(dat) || nrow(dat) != ncol(x2$mytable_selected)) {
        x2$tf_list <- rep(TRUE, ncol(x2$mytable_selected))
        x2$tf_list_anat <- x2$tf_list
        return()
      }
      
      if (!is.null(input$show_dat) && input$seg_pdat_col %in% colnames(dat)) {
        tf_tmp <- dat[, input$seg_pdat_col] %in% input$show_dat
        tf_tmp[is.na(tf_tmp)] <- FALSE
        x2$tf_list <- tf_tmp
      } else {
        x2$tf_list <- rep(TRUE, nrow(dat))
      }
      
      # Keep anat mask aligned with displayed/editable pData subset.
      x2$tf_list_anat <- x2$tf_list
    }, ignoreInit = FALSE)
    
    # adjust T/F list from UMAP color filter (background mode)
    observeEvent({
      input$seg_choice
      input$cols
      1
    }, {
      x2 <- preproc_values()[["x2"]]
      if (!identical(input$seg_choice, "bk_seg")) {
        return()
      }
      req(x2$data_list)
      tf_tmp <- x2$bkcols %in% input$cols
      tf_tmp[is.na(tf_tmp)] <- FALSE
      x2$tf_list <- tf_tmp
    }, ignoreInit = FALSE)
    
    #save processed data in a list when save button pressed and plot results
    
    # observeEvent(input$store_proc, {
    #   x2 <- preproc_values()[["x2"]]
    #   
    #   req(x2$tf_list) #this will tell us  which pixels to use / store??
    #   print("")
    #   print("Storing processed data")
    #   showNotification("Storing processed Data")
    #   
    #   
    #   if (input$seg_choice == "anat_seg") {
    #     #find colnames only in processed data, and add those to the overview data
    #     
    #     if (!is.null(x2$tf_list_anat) 
    #          &&
    #          (sum(x2$tf_list_anat) != sum(x2$tf_list))) {
    #         
    #       print("Storing UMAP results as processed data")
    #       
    #       
    #       x2$tf_list <- x2$tf_list_anat
    #       aa <- x2$tf_list_anat
    #       x2$tf_list_anat <- NULL
    #       
    #       req(x2$data_list)
    #       #x2<-preproc_values()[["x2"]]
    #       
    #       if (is.null(x2$overview_peaks_sel$Rcol_reduced)) {
    #         x2$overview_peaks_sel$Rcol_reduced <- "NA"
    #       }
    #       if (is.null(x2$overview_peaks_sel$col_dbscan)) {
    #         x2$overview_peaks_sel$col_dbscan <- "NA"
    #       }
    #       
    #       if (is.null(x2$overview_peaks_sel$col_hdbscan)) {
    #         x2$overview_peaks_sel$col_hdbscan <- "NA"
    #       }
    #       
    #       tmp.img <- x2$mytable_selected[x2$tf_list]
    #       x2$pdat_anat <- x2$pdata_anat[x2$tf_list,]
    #       
    #       
    #       
    #       #store umap embeddings
    #       #tmp.img$x_umaps <-
    #       #  x2$data_list$umap_separation$[x2$tf_list]
    #       
    #       
    #       tmp.img$Rcol_reduced <-
    #         x2$data_list$umap_separation$col_reduced[x2$tf_list]
    #       tmp.img$col_dbscan <-
    #         x2$data_list$dbscan_umap_separation$cluster[x2$tf_list]
    #       tmp.img$col_hdbscan <-
    #         x2$data_list$hdbscan_umap_separation$cluster[x2$tf_list]
    #       
    #       mytable_selected = tmp.img
    #       x2$tf_list <- rep(TRUE, ncol(mytable_selected))
    #     } else {
    #       mytable_selected <- x2$mytable_selected
    #       if (dim(pData(mytable_selected))[1] != dim(x2$pdat_anat)[1]) {
    #         print(paste0(
    #           "processed data size, ",
    #           dim(x2$pdat_anat)[1],
    #           " does not match image size, ",
    #           dim(pData(mytable_selected))[1]
    #         ))
    #         print("choose image matching processed data size before storing")
    #         print(paste0("possible name: ", x2$umap_name))
    #         return()
    #       } else {
    #         pData(mytable_selected) <- x2$pdat_anat
    #       }
    #     }
    #     
    #     if (!setequal(colnames(pData(x2$overview_peaks_sel)), colnames(pData(mytable_selected)))) {
    #       cols <-
    #         union(colnames(pData(x2$overview_peaks_sel)), colnames(pData(mytable_selected)))
    #       cols_overview <-
    #         cols[!cols %in% (colnames(pData(x2$overview_peaks_sel)))]
    #       cols_selected <-
    #         cols[!cols %in% (colnames(pData(mytable_selected)))]
    #       
    #       pData(x2$overview_peaks_sel)[, cols_overview] <- "NA"
    #       pData(mytable_selected)[, cols_selected] <- "NA"
    #     }
    #     
    #     
    #     
    #     tmp.img <- mytable_selected %>% subsetPixels(x2$tf_list)
    #     
    #   } else if (input$seg_choice == "fix_pix") {
    #     tmp.img <- x2$mytable_selected
    #     
    #     all_runs <-
    #       runNames(x2$overview_peaks_sel) #from the total set being analyzed
    #     runs <- runNames(tmp.img)
    #     idx <-
    #       which(all_runs %in% runs)  #are we assuming that the runs are stored in the same order?
    #     
    #     #put all runs into a separate list from the originally imported / peakpicked data
    #     for (i in 1:length(all_runs)) {
    #       x2$list_proc_img[[i]] <-
    #         x2$overview_peaks_sel %>% subsetPixels(Cardinal::run(x2$overview_peaks_sel) ==
    #                                                  all_runs[i])
    #       fData(x2$list_proc_img[[i]]) <-
    #         fData(x2$list_proc_img[[i]])[unique(colnames(fData(x2$list_proc_img[[i]])))] #clean up featuredata
    #     }
    #     
    #     #replace list elements with updated processed data
    #     for (i in 1:length(runs)) {
    #       x2$list_proc_img[[idx[i]]] <-
    #         tmp.img %>% subsetPixels(Cardinal::run(tmp.img) %in% runs[i])
    #       fData(x2$list_proc_img[[idx[i]]]) <-
    #         fData(x2$list_proc_img[[idx[i]]])[unique(colnames(fData(x2$list_proc_img[[idx[i]]])))] #clean up featuredata
    #     }
    #     
    #     
    #     
    #     
    #   } else{
    #     
    #     
    #     req(x2$data_list)
    #     #x2<-preproc_values()[["x2"]]
    # 
    #     if (is.null(x2$overview_peaks_sel$Rcol_reduced)) {
    #       x2$overview_peaks_sel$Rcol_reduced <- "NA"
    #     }
    #     if (is.null(x2$overview_peaks_sel$col_dbscan)) {
    #       x2$overview_peaks_sel$col_dbscan <- "NA"
    #     }
    #     
    #     if (is.null(x2$overview_peaks_sel$col_hdbscan)) {
    #       x2$overview_peaks_sel$col_hdbscan <- "NA"
    #     }
    #     
    #     #input$seg_field add this column to the rest of the data!
    #     
    #     #tmp.img<-x2$mytable_selected %>% subsetPixels(x2$tf_list)
    #     
    #     
    #     
    #     
    #     tmp.img <- tmp_umap_dat()[x2$tf_list]
    #     
    #     tmp.img$Rcol_reduced <-
    #       x2$data_list$umap_separation$col_reduced[x2$tf_list]
    #     tmp.img$col_dbscan <-
    #       x2$data_list$dbscan_umap_separation$cluster[x2$tf_list]
    #     tmp.img$col_hdbscan <-
    #       x2$data_list$hdbscan_umap_separation$cluster[x2$tf_list]
    #     
    #     #save umap embeddings
    #     embeddings<-x2$data_list$umap_separation$umap_out[x2$tf_list,]
    #     tmp.img$x_umap <- embeddings[,1]
    #     tmp.img$y_umap <- embeddings[,2]
    #     tmp.img$z_umap <- embeddings[,3]
    #     
    #   }
    #   
    #   
    #   tmp.img <- tmp_umap_dat()[x2$tf_list]
    #   
    #   tmp.img$Rcol_reduced <- x2$data_list$umap_separation$col_reduced[x2$tf_list]
    #   tmp.img$col_dbscan <- x2$data_list$dbscan_umap_separation$cluster[x2$tf_list]
    #   tmp.img$col_hdbscan <- x2$data_list$hdbscan_umap_separation$cluster[x2$tf_list]
    #   
    #   #browser()
    #   
    #   # Add assignments for new clustering methods
    #   if ("kmeans" %in% input$clustering_methods) {
    #     tmp.img$col_kmeans <- x2$data_list$kmeans_umap_separation$cluster[x2$tf_list]
    #   }
    #   if ("hierarchical" %in% input$clustering_methods) {
    #     tmp.img$col_hierarchical <- x2$data_list$hierarchical_umap_separation[x2$tf_list]
    #   }
    #   if ("spectral" %in% input$clustering_methods) {
    #     tmp.img$col_spectral <- x2$data_list$spectral_umap_separation$cluster[x2$tf_list]
    #   }
    #   if ("kmedoids" %in% input$clustering_methods) {
    #     tmp.img$col_kmedoids <- x2$data_list$kmedoids_umap_separation$clustering[x2$tf_list]
    #   }
    #   if ("fuzzy" %in% input$clustering_methods) {
    #     # For fuzzy clustering, get the hard assignments
    #     tmp.img$col_fuzzy <- apply(x2$data_list$fuzzy_umap_separation$membership[x2$tf_list, ], 1, which.max)
    #   }
    #   if ("mclust" %in% input$clustering_methods) {
    #     tmp.img$col_mclust <- x2$data_list$mclust_umap_separation$classification[x2$tf_list]
    #   }
    #   if ("som" %in% input$clustering_methods) {
    #     tmp.img$col_som <- x2$data_list$som_umap_separation$unit.classif[x2$tf_list]
    #   }
    #   if ("skmeans" %in% input$clustering_methods) {
    #     tmp.img$col_skmeans <- x2$data_list$skmeans_umap_separation$cluster[x2$tf_list]
    #   }
    #   
    #   # Save UMAP embeddings
    #   embeddings <- x2$data_list$umap_separation$umap_out[x2$tf_list, ]
    #   tmp.img$x_umap <- embeddings[, 1]
    #   tmp.img$y_umap <- embeddings[, 2]
    #   tmp.img$z_umap <- embeddings[, 3]
    #   
    #   
    #   
    #   all_runs <-
    #     runNames(x2$overview_peaks_sel) #from the total set being analyzed
    #   runs <- runNames(tmp.img)
    #   idx <-
    #     which(all_runs %in% runs)  #are we assuming that the runs are stored in the same order?
    #   
    #   # if x2$list_proc_img is NULL, then we need to create it
    #   if (is.null(x2$list_proc_img)) {
    #     x2$list_proc_img <- vector("list", length = length(all_runs))
    #            
    #     #put all runs into a separate list from the originally imported / peakpicked data
    #     for (i in 1:length(all_runs)) {
    #       x2$list_proc_img[[i]] <-
    #         x2$overview_peaks_sel %>% subsetPixels(Cardinal::run(x2$overview_peaks_sel) ==
    #                                                  all_runs[i])
    #       
    #         a<-(fData(x2$list_proc_img[[i]])[unique(colnames(fData(x2$list_proc_img[[i]])))]) #clean up featuredata
    #         
    #         #need to use massdatafram explicitly
    #         fData(x2$list_proc_img[[i]]) <- MassDataFrame(mz=a$mz, a %>% as.data.frame() %>% dplyr::select(-mz))
    #         
    #     }
    #   }
    #   
    #   #replace list elements with updated processed data ###TODO
    #   for (i in 1:length(runs)) {
    #     x2$list_proc_img[[idx[i]]] <-
    #       tmp.img %>% subsetPixels(Cardinal::run(tmp.img) %in% runs[i])
    #      a<- fData(x2$list_proc_img[[idx[i]]])[unique(colnames(fData(x2$list_proc_img[[idx[i]]])))] #clean up featuredata
    #      
    #      fData(x2$list_proc_img[[idx[i]]]) <- MassDataFrame(mz=a$mz, a %>% as.data.frame() %>% dplyr::select(-mz))
    #   }
    #   
    #   
    # 
    #   x2$data_list <- NULL
    # 
    #   x2$mytable_selected <-
    #     combine_card(x2$list_proc_img[input$mytable_rows_selected])
    #   #browser()
    #     #do.call(cbind, x2$list_proc_img[input$mytable_rows_selected]) #TODO check carefully!!
    #   # tmp_umap_dat<-reactive({return(tmp.img)})
    #   x2$tf_list<-rep(TRUE, ncol(combine_card(x2$list_proc_img[input$mytable_rows_selected])))
    #   
    #   
    # })
    # 
    
    observeEvent(input$store_proc, {
      x2 <- preproc_values()[["x2"]]
      if (is.null(x2$mytable_selected)) {
        showNotification("No active dataset to store. Select runs and run UMAP/filters first.", type = "warning", duration = 6)
        return()
      }
      if (is.null(x2$tf_list)) {
        showNotification("No pixel-selection mask available to store. Run clustering or filtering first.", type = "warning", duration = 6)
        return()
      }
      tf_valid <- sanitize_tf_mask(x2$tf_list, ncol(x2$mytable_selected))
      if (is.null(tf_valid) || sum(tf_valid) == 0) {
        showNotification("No pixels selected to store.", type = "warning", duration = 6)
        return()
      }
      print("")
      print("Storing processed data")
      showNotification("Storing processed Data")
      
      if (input$seg_choice == "anat_seg") {
        # Find colnames only in processed data, and add those to the overview data
        
        if (!is.null(x2$tf_list_anat) &&
            (sum(x2$tf_list_anat) != sum(x2$tf_list))) {
          
          print("Storing UMAP results as processed data")
          
          x2$tf_list <- x2$tf_list_anat
          x2$tf_list_anat <- NULL
          tf_valid <- sanitize_tf_mask(x2$tf_list, ncol(x2$mytable_selected))
          if (is.null(tf_valid) || sum(tf_valid) == 0) {
            showNotification("No pixels selected after applying anatomical filter.", type = "warning", duration = 6)
            return()
          }
          
          req(x2$data_list)
          
          # Initialize columns if they don't exist
          if (is.null(x2$overview_peaks_sel$Rcol_reduced)) {
            x2$overview_peaks_sel$Rcol_reduced <- "NA"
          }
          if (is.null(x2$overview_peaks_sel$col_dbscan)) {
            x2$overview_peaks_sel$col_dbscan <- "NA"
          }
          if (is.null(x2$overview_peaks_sel$col_hdbscan)) {
            x2$overview_peaks_sel$col_hdbscan <- "NA"
          }
          # Add initializations for new clustering methods
          if (is.null(x2$overview_peaks_sel$col_kmeans)) {
            x2$overview_peaks_sel$col_kmeans <- "NA"
          }
          if (is.null(x2$overview_peaks_sel$col_hierarchical)) {
            x2$overview_peaks_sel$col_hierarchical <- "NA"
          }
          if (is.null(x2$overview_peaks_sel$col_spectral)) {
            x2$overview_peaks_sel$col_spectral <- "NA"
          }
          if (is.null(x2$overview_peaks_sel$col_kmedoids)) {
            x2$overview_peaks_sel$col_kmedoids <- "NA"
          }
          if (is.null(x2$overview_peaks_sel$col_fuzzy)) {
            x2$overview_peaks_sel$col_fuzzy <- "NA"
          }
          if (is.null(x2$overview_peaks_sel$col_mclust)) {
            x2$overview_peaks_sel$col_mclust <- "NA"
          }
          if (is.null(x2$overview_peaks_sel$col_som)) {
            x2$overview_peaks_sel$col_som <- "NA"
          }
          if (is.null(x2$overview_peaks_sel$col_skmeans)) {
            x2$overview_peaks_sel$col_skmeans <- "NA"
          }
          
          tmp.img <- safe_subset_pixels(x2$mytable_selected, tf_valid, "store processed data", notify = FALSE)
          if (is.null(tmp.img) || ncol(tmp.img) == 0) {
            showNotification("No image data available to store after filtering.", type = "warning", duration = 6)
            return()
          }
          if (!is.null(x2$pdat_anat) && nrow(x2$pdat_anat) == length(tf_valid)) {
            x2$pdat_anat <- x2$pdat_anat[tf_valid, ]
          }
          
          # Store UMAP embeddings
          # tmp.img$x_umaps <- x2$data_list$umap_separation$[x2$tf_list]
          
          tmp.img$Rcol_reduced <- x2$data_list$umap_separation$col_reduced[tf_valid]
          tmp.img$col_dbscan <- x2$data_list$dbscan_umap_separation$cluster[tf_valid]
          tmp.img$col_hdbscan <- x2$data_list$hdbscan_umap_separation$cluster[tf_valid]
          
          # Add assignments for new clustering methods
          if ("kmeans" %in% input$clustering_methods) {
            tmp.img$col_kmeans <- x2$data_list$kmeans_umap_separation$cluster[tf_valid]
          }
          if ("hierarchical" %in% input$clustering_methods) {
            tmp.img$col_hierarchical <- x2$data_list$hierarchical_umap_separation[tf_valid]
          }
          if ("spectral" %in% input$clustering_methods) {
            tmp.img$col_spectral <- x2$data_list$spectral_umap_separation$cluster[tf_valid]
          }
          if ("kmedoids" %in% input$clustering_methods) {
            tmp.img$col_kmedoids <- x2$data_list$kmedoids_umap_separation$clustering[tf_valid]
          }
          if ("fuzzy" %in% input$clustering_methods) {
            # For fuzzy clustering, get the hard assignments
            tmp.img$col_fuzzy <- apply(x2$data_list$fuzzy_umap_separation$membership[tf_valid, ], 1, which.max)
          }
          if ("mclust" %in% input$clustering_methods) {
            tmp.img$col_mclust <- x2$data_list$mclust_umap_separation$classification[tf_valid]
          }
          if ("som" %in% input$clustering_methods) {
            tmp.img$col_som <- x2$data_list$som_umap_separation$unit.classif[tf_valid]
          }
          if ("skmeans" %in% input$clustering_methods) {
            tmp.img$col_skmeans <- x2$data_list$skmeans_umap_separation$cluster[tf_valid]
          }
          
          mytable_selected <- tmp.img
          x2$tf_list <- rep(TRUE, ncol(mytable_selected))
        } else {
          mytable_selected <- x2$mytable_selected
          if (dim(pData(mytable_selected))[1] != dim(x2$pdat_anat)[1]) {
            print(paste0(
              "Processed data size, ",
              dim(x2$pdat_anat)[1],
              " does not match image size, ",
              dim(pData(mytable_selected))[1]
            ))
            print("Choose image matching processed data before storing")
            print(paste0("Possible name: ", x2$umap_name))
            return()
          } else {
            pData(mytable_selected) <- x2$pdat_anat
          }
        }
        
        if (!setequal(colnames(pData(x2$overview_peaks_sel)), colnames(pData(mytable_selected)))) {
          cols <- union(colnames(pData(x2$overview_peaks_sel)), colnames(pData(mytable_selected)))
          cols_overview <- cols[!cols %in% (colnames(pData(x2$overview_peaks_sel)))]
          cols_selected <- cols[!cols %in% (colnames(pData(mytable_selected)))]
          
          pData(x2$overview_peaks_sel)[, cols_overview] <- "NA"
          pData(mytable_selected)[, cols_selected] <- "NA"
        }
        
        tmp.img <- safe_subset_pixels(mytable_selected, x2$tf_list, "store processed data", notify = FALSE)
        if (is.null(tmp.img) || ncol(tmp.img) == 0) {
          showNotification("No image data available to store after filtering.", type = "warning", duration = 6)
          return()
        }
        
      } else if (input$seg_choice == "fix_pix") {
        tmp.img <- x2$mytable_selected
        
        all_runs <- runNames(x2$overview_peaks_sel)  # From the total set being analyzed
        runs <- runNames(tmp.img)
        idx <- which(all_runs %in% runs)  # Are we assuming that the runs are stored in the same order?
        
        # Put all runs into a separate list from the originally imported/peakpicked data
        for (i in 1:length(all_runs)) {
          x2$list_proc_img[[i]] <-
            x2$overview_peaks_sel %>% subsetPixels(Cardinal::run(x2$overview_peaks_sel) == all_runs[i])
          fData(x2$list_proc_img[[i]]) <-
            fData(x2$list_proc_img[[i]])[unique(colnames(fData(x2$list_proc_img[[i]])))]  # Clean up feature data
        }
        
        # Replace list elements with updated processed data
        for (i in 1:length(runs)) {
          x2$list_proc_img[[idx[i]]] <-
            tmp.img %>% subsetPixels(Cardinal::run(tmp.img) %in% runs[i])
          fData(x2$list_proc_img[[idx[i]]]) <-
            fData(x2$list_proc_img[[idx[i]]])[unique(colnames(fData(x2$list_proc_img[[idx[i]]])))]  # Clean up feature data
        }
        
      } else {
        req(x2$data_list)
        # x2 <- preproc_values()[["x2"]]
        # browser()
        
        if (is.null(x2$overview_peaks_sel$Rcol_reduced)) {
          x2$overview_peaks_sel$Rcol_reduced <- "NA"
        }
        if (is.null(x2$overview_peaks_sel$col_dbscan)) {
          x2$overview_peaks_sel$col_dbscan <- "NA"
        }
        if (is.null(x2$overview_peaks_sel$col_hdbscan)) {
          x2$overview_peaks_sel$col_hdbscan <- "NA"
        }
        # Initialize columns for new clustering methods
        if (is.null(x2$overview_peaks_sel$col_kmeans)) {
          x2$overview_peaks_sel$col_kmeans <- "NA"
        }
        if (is.null(x2$overview_peaks_sel$col_hierarchical)) {
          x2$overview_peaks_sel$col_hierarchical <- "NA"
        }
        if (is.null(x2$overview_peaks_sel$col_spectral)) {
          x2$overview_peaks_sel$col_spectral <- "NA"
        }
        if (is.null(x2$overview_peaks_sel$col_kmedoids)) {
          x2$overview_peaks_sel$col_kmedoids <- "NA"
        }
        if (is.null(x2$overview_peaks_sel$col_fuzzy)) {
          x2$overview_peaks_sel$col_fuzzy <- "NA"
        }
        if (is.null(x2$overview_peaks_sel$col_mclust)) {
          x2$overview_peaks_sel$col_mclust <- "NA"
        }
        if (is.null(x2$overview_peaks_sel$col_som)) {
          x2$overview_peaks_sel$col_som <- "NA"
        }
        if (is.null(x2$overview_peaks_sel$col_skmeans)) {
          x2$overview_peaks_sel$col_skmeans <- "NA"
        }
        
        # input$seg_field: add this column to the rest of the data!
        
        base_dat <- tmp_umap_dat()
        tf_valid <- sanitize_tf_mask(x2$tf_list, ncol(base_dat))
        if (is.null(tf_valid) || sum(tf_valid) == 0) {
          showNotification("No pixels selected to store from UMAP results.", type = "warning", duration = 6)
          return()
        }
        tmp.img <- safe_subset_pixels(base_dat, tf_valid, "store UMAP processed data", notify = FALSE)
        if (is.null(tmp.img) || ncol(tmp.img) == 0) {
          showNotification("No image data available after UMAP filtering.", type = "warning", duration = 6)
          return()
        }
        
        tmp.img$Rcol_reduced <- x2$data_list$umap_separation$col_reduced[tf_valid]
        tmp.img$col_dbscan <- x2$data_list$dbscan_umap_separation$cluster[tf_valid]
        tmp.img$col_hdbscan <- x2$data_list$hdbscan_umap_separation$cluster[tf_valid]
        
        
        # Add assignments for new clustering methods
        if ("kmeans" %in% input$clustering_methods) {
          tmp.img$col_kmeans <- x2$data_list$kmeans_umap_separation$cluster[tf_valid]
        }
        if ("hierarchical" %in% input$clustering_methods) {
          tmp.img$col_hierarchical <- x2$data_list$hierarchical_umap_separation[tf_valid]
        }
        if ("spectral" %in% input$clustering_methods) {
          tmp.img$col_spectral <- x2$data_list$spectral_umap_separation$cluster[tf_valid]
        }
        if ("kmedoids" %in% input$clustering_methods) {
          tmp.img$col_kmedoids <- x2$data_list$kmedoids_umap_separation$clustering[tf_valid]
        }
        if ("fuzzy" %in% input$clustering_methods) {
          # For fuzzy clustering, get the hard assignments
          tmp.img$col_fuzzy <- apply(x2$data_list$fuzzy_umap_separation$membership[tf_valid, ], 1, which.max)
        }
        if ("mclust" %in% input$clustering_methods) {
          tmp.img$col_mclust <- x2$data_list$mclust_umap_separation$classification[tf_valid]
        }
        if ("som" %in% input$clustering_methods) {
          tmp.img$col_som <- x2$data_list$som_umap_separation$unit.classif[tf_valid]
        }
        if ("skmeans" %in% input$clustering_methods) {
          tmp.img$col_skmeans <- x2$data_list$skmeans_umap_separation$cluster[tf_valid]
        }
        
        # Save UMAP embeddings
        embeddings <- x2$data_list$umap_separation$umap_out[tf_valid, ]
        tmp.img$x_umap <- embeddings[, 1]
        tmp.img$y_umap <- embeddings[, 2]
        tmp.img$z_umap <- embeddings[, 3]
        
      }
      
      all_runs <- runNames(x2$overview_peaks_sel)  # From the total set being analyzed
      runs <- runNames(tmp.img)
      idx <- which(all_runs %in% runs)  # Are we assuming that the runs are stored in the same order?
      
      # If x2$list_proc_img is NULL, then we need to create it
      if (is.null(x2$list_proc_img)) {
        x2$list_proc_img <- vector("list", length = length(all_runs))
        
        # Put all runs into a separate list from the originally imported/peakpicked data
        for (i in 1:length(all_runs)) {
          x2$list_proc_img[[i]] <-
            x2$overview_peaks_sel %>% subsetPixels(Cardinal::run(x2$overview_peaks_sel) == all_runs[i])
          a <- fData(x2$list_proc_img[[i]])[unique(colnames(fData(x2$list_proc_img[[i]])))]  # Clean up feature data
          
          # Need to use MassDataFrame explicitly
          fData(x2$list_proc_img[[i]]) <- MassDataFrame(mz = a$mz, a %>% as.data.frame() %>% dplyr::select(-mz))
        }
      }
      # Replace list elements with updated processed data
      for (i in 1:length(runs)) {
        x2$list_proc_img[[idx[i]]] <-
          tmp.img %>% subsetPixels(Cardinal::run(tmp.img) %in% runs[i])
        a <- fData(x2$list_proc_img[[idx[i]]])[unique(colnames(fData(x2$list_proc_img[[idx[i]]])))]  # Clean up feature data
        
        fData(x2$list_proc_img[[idx[i]]]) <- MassDataFrame(mz = a$mz, a %>% as.data.frame() %>% dplyr::select(-mz))
      }
      
      x2$data_list <- NULL
      
      x2$mytable_selected <- combine_card(x2$list_proc_img[input$mytable_rows_selected])
      # browser()
      # do.call(cbind, x2$list_proc_img[input$mytable_rows_selected])  # TODO check carefully!!
      # tmp_umap_dat <- reactive({return(tmp.img)})
      x2$tf_list <- rep(TRUE, ncol(combine_card(x2$list_proc_img[input$mytable_rows_selected])))
    })
    
    
    output$umap_choice <- renderUI({
      make_umap_controls <- function(show_anat_editor = FALSE) {
        tagList(
          tags$div(style = "font-weight:600; margin-bottom:6px;", "Common parameters"),
          fluidRow(
            column(
              6,
              numericInput(
                ns("nn"),
                "Nearest neighbors",
                min = 2,
                max = 100,
                value = 8
              )
            ),
            column(
              6,
              numericInput(ns("min_dist"), "Min distance", min = 0, value = 0.02, step = 0.01)
            )
          ),
          fluidRow(
            column(6, textInput(ns("pca_umap"), "PCA components (NULL = none)", value = "")),
            column(
              6,
              selectInput(
                ns("umap_metric"),
                "Metric",
                choices = c("cosine", "euclidean", "manhattan", "correlation"),
                selected = "cosine"
              )
            )
          ),
          fluidRow(
            column(
              6,
              numericInput(
                ns("set_op_mix_ratio"),
                "set_op_mix_ratio",
                min = 0,
                max = 1,
                value = 0.15,
                step = 0.01
              )
            ),
            column(6, numericInput(ns("search_k"), "search_k", value = 1000, min = 1))
          ),
          fluidRow(
            column(6, numericInput(ns("n_trees"), "n_trees", value = 100, min = 1))
          ),
          tags$div(style = "font-weight:600; margin-top:6px; margin-bottom:6px;", "Advanced parameters"),
          fluidRow(
            column(6, numericInput(ns("umap_spread"), "Spread", value = 1, min = 0.01, step = 0.05)),
            column(6, numericInput(ns("local_connectivity"), "Local connectivity", value = 1, min = 0, step = 0.5))
          ),
          fluidRow(
            column(6, numericInput(ns("repulsion_strength"), "Repulsion", value = 1.2, min = 0.01, step = 0.1)),
            column(6, numericInput(ns("negative_sample_rate"), "Negative sample rate", value = 10, min = 1))
          ),
          fluidRow(
            column(6, numericInput(ns("learning_rate"), "Learning rate", value = 1, min = 0.001, step = 0.1)),
            column(6, numericInput(ns("n_epochs"), "n_epochs (0 = auto)", value = 600, min = 0))
          ),
          selectInput(
            ns("umap_init"),
            "Initialization",
            choices = c("spectral", "random", "pca", "spca"),
            selected = "spectral"
          ),
          checkboxInput(
            ns("fastmap"),
            "Use spatialFastmap (experimental)",
            value = FALSE
          ),
          uiOutput(ns("fm_params")),
          actionButton(ns("umap1"), label = "Start UMAP"),
          if (show_anat_editor) {
            tagList(
              tags$hr(style = "margin-top:8px;margin-bottom:8px;"),
              uiOutput(ns("anat_pdata"))
            )
          }
        )
      }
      
      switch(
        input$seg_choice,
        "bk_seg" = make_umap_controls(show_anat_editor = FALSE),
        "anat_seg" = make_umap_controls(show_anat_editor = TRUE),
        NULL
      )
    })
    
    output$use_python <- renderUI({
      # if (input$use_python) {
      #   list(
      #     numericInput(ns("minSize"), "Python HDBSCAN min size", value = 15L),
      #     numericInput(
      #       ns("epsilon"),
      #       "Python HDBSCAN cluster_selection_epsilon",
      #       value = 0.2
      #     ),
      #     selectInput(
      #       ns("hdbscan_method"),
      #       "Method for HDBSCAN calculation",
      #       choices = list(
      #         "R DBSCAN package" = "R_hdbscan",
      #         "python HDBSCAN" = "py_hdbscan",
      #         None = "hdbscan_none"
      #       )
      #     )
      #
      #   )
      # } else {
      # list(
      #   selectInput(
      #   ns("hdbscan_method"),
      #   "Method for HDBSCAN calculation",
      #   choices = list("R DBSCAN package" = "R_hdbscan", None =
      #                    "hdbscan_none")
      # ))
      #}
    })
    
    output$anat_pdata <- renderUI({
      x2 <- preproc_values()[["x2"]]
      
      #extract pdata
      req(x2$mytable_selected)
      
      pdat <- get_active_pdata_df(x2)
      req(pdat)
      pdat_cols <- colnames(pdat)
      if (length(pdat_cols) == 0) {
        return(tags$div("No pData fields available for editing."))
      }
      default_col <- if ("Rcol_reduced" %in% pdat_cols) "Rcol_reduced" else pdat_cols[1]
      
      tagList(
        selectInput(
          ns("seg_pdat_col"),
          "pData field for visualization/editing",
          choices = pdat_cols,
          selected = default_col
        ),
        uiOutput(ns("seg_pdat_show")),
        radioButtons(
          ns("seg_rename"),
          "Edit action",
          choices = c(
            "Edit values in existing field" = "replace",
            "Create new field from selected pixels" = "add",
            "Rename pData field" = "rename_field"
          ),
          selected = "replace"
        ),
        uiOutput(ns("seg_field")),
        actionButton(ns("apply_annotation"), "Apply pData Edit")
      )
      
    })
    
    output$seg_field <- renderUI ({
      x2 <- preproc_values()[["x2"]]
      req(x2$mytable_selected)
      pdat <- get_active_pdata_df(x2)
      req(pdat)
      
      choices <- colnames(pdat)
      if (length(choices) == 0) {
        return(NULL)
      }
      editable_choices <- choices[!choices %in% c("x", "y", "run")]
      if (length(editable_choices) == 0) {
        editable_choices <- choices
      }
      selected_field <- if (!is.null(input$seg_pdat_col) && input$seg_pdat_col %in% editable_choices) {
        input$seg_pdat_col
      } else {
        editable_choices[1]
      }
      
      switch(
        req(input$seg_rename),
        "replace" = list(
          selectInput(
            ns("seg_field"),
            "Field to edit",
            choices = editable_choices,
            selected = selected_field
          ),
          uiOutput(ns("current_seg_choices")),
          textInput(ns("new_seg_value"), "New value for selected pixels")
        ),
        "add" = list(
          textInput(ns("seg_field"), "Name of new field to create"),
          textInput(ns("new_seg_value"), "Value for selected pixels")
        ),
        "rename_field" = list(
          selectInput(ns("seg_field"), "Field to rename", choices = editable_choices, selected = selected_field),
          textInput(ns("new_seg_value"), "New field name")
        )
      )
      
    })
    
    output$current_seg_choices <- renderUI ({
      x2 <- preproc_values()[["x2"]]
      req(x2$mytable_selected, input$seg_field)
      pdat <- get_active_pdata_df(x2)
      req(pdat)
      if (!input$seg_field %in% colnames(pdat)) {
        return(NULL)
      }
      vals <- unique(pdat[, input$seg_field])
      vals <- vals[!is.na(vals)]
      vals <- unique(as.character(vals))
      vals <- vals[order(vals)]
      
      selectInput(
        ns("current_seg_value"),
        "Current value filter (optional)",
        choices = c("[Any current value]" = "__ANY__", vals),
        selected = "__ANY__"
      )
      
      
    })
    
    output$seg_pdat_show <- renderUI ({
      x2 <- preproc_values()[["x2"]]
      req(x2$mytable_selected)
      
      pdat <- get_active_pdata_df(x2)
      req(pdat, input$seg_pdat_col)
      if (!input$seg_pdat_col %in% colnames(pdat)) {
        return(selectizeInput(
          ns("show_dat"),
          "Values to display / edit",
          choices = character(0),
          selected = character(0),
          multiple = TRUE
        ))
      }
      vals <- unique(as.character(pdat[, input$seg_pdat_col]))
      vals <- vals[!is.na(vals)]
      vals <- vals[order(vals)]
      current_sel <- isolate(input$show_dat)
      selected_vals <- if (!is.null(current_sel)) intersect(as.character(current_sel), vals) else vals
      if (length(selected_vals) == 0 && length(vals) > 0) {
        selected_vals <- vals
      }
      
      selectizeInput(
        ns("show_dat"),
        "Values to display / edit",
        choices = vals,
        selected = selected_vals,
        multiple = TRUE
      )
    })
    
    output$fix_pix_opts <- renderUI({
      if (input$seg_choice == "fix_pix") {
        list(
          numericInput(ns("fix_pix_r"), "radius to search for neighbors", 2),
          numericInput(ns("fix_pix_n_thresh"), "min # of neighbors allowed", 3),
          checkboxInput(ns("fix_pix_t_f"), "Actually remove pixels?", value =
                          TRUE),
          actionButton(ns("fix_pix"), HTML("Remove isolated pixels"))
        )
      }
    })
    
    #Apply the results of annotations for specific tissue regions
    observeEvent(input$apply_annotation, {
      x2 <- preproc_values()[["x2"]]
      if (is.null(x2$mytable_selected)) {
        showNotification("No dataset selected for annotation. Choose runs first.", type = "warning", duration = 7)
        return()
      }
      if (is.null(input$seg_rename) || !nzchar(input$seg_rename)) {
        showNotification("Choose an annotation edit mode before applying changes.", type = "warning", duration = 7)
        return()
      }
      
      pdat <- get_active_pdata(x2)
      if (is.null(pdat)) {
        showNotification("No pData available to edit.", type = "error")
        return()
      }
      pdat_df <- as.data.frame(pdat)
      n_pixels <- nrow(pdat_df)
      if (n_pixels != ncol(x2$mytable_selected)) {
        showNotification("pData length does not match selected image. Reload data selection and retry.", type = "error", duration = 8)
        return()
      }
      
      field_name <- trimws(input$seg_field)
      new_value <- trimws(input$new_seg_value)
      
      selected_mask <- rep(TRUE, n_pixels)
      if (!is.null(input$show_dat) && length(input$show_dat) > 0 && input$seg_pdat_col %in% colnames(pdat_df)) {
        selected_mask <- pdat_df[, input$seg_pdat_col] %in% input$show_dat
        selected_mask[is.na(selected_mask)] <- FALSE
      }
      tf_current <- sanitize_tf_mask(x2$tf_list, n_pixels)
      if (!is.null(tf_current)) {
        selected_mask <- selected_mask & tf_current
      }
      
      if (identical(input$seg_rename, "rename_field")) {
        if (!nzchar(field_name) || !field_name %in% colnames(pdat_df)) {
          showNotification("Choose an existing pData field to rename.", type = "error")
          return()
        }
        if (!nzchar(new_value)) {
          showNotification("Enter a new field name.", type = "error")
          return()
        }
        if (new_value %in% colnames(pdat_df) && !identical(new_value, field_name)) {
          showNotification("That field name already exists. Choose a different name.", type = "error")
          return()
        }
        
        colnames(pdat)[colnames(pdat) == field_name] <- new_value
        x2$pdat_anat <- pdat
        
        updated_cols <- colnames(as.data.frame(x2$pdat_anat))
        updateSelectInput(session, "seg_pdat_col", choices = updated_cols, selected = new_value)
        updateSelectInput(session, "seg_field", choices = updated_cols, selected = new_value)
        showNotification(sprintf("Renamed pData field '%s' to '%s'.", field_name, new_value), type = "message", duration = 5)
        return()
      }
      
      selected_idx <- which(selected_mask)
      if (length(selected_idx) == 0) {
        showNotification("No pixels selected for editing. Adjust 'Values to display / edit' first.", type = "warning", duration = 6)
        return()
      }
      
      if (identical(input$seg_rename, "replace")) {
        if (!nzchar(field_name) || !field_name %in% colnames(pdat_df)) {
          showNotification("Choose an existing field to edit.", type = "error")
          return()
        }
        if (!nzchar(new_value)) {
          showNotification("Enter a new value to apply.", type = "error")
          return()
        }
        
        current_filter <- if (is.null(input$current_seg_value)) "__ANY__" else input$current_seg_value
        if (!identical(current_filter, "__ANY__")) {
          existing_vals <- as.character(pdat_df[, field_name])
          keep <- existing_vals[selected_idx] %in% as.character(current_filter)
          selected_idx <- selected_idx[which(keep)]
        }
        
        if (length(selected_idx) == 0) {
          showNotification("No pixels matched the selected current-value filter.", type = "warning", duration = 6)
          return()
        }
        
        updated_vals <- as.character(pdat_df[, field_name])
        updated_vals[selected_idx] <- new_value
        pdat[[field_name]] <- updated_vals
        x2$pdat_anat <- pdat
        showNotification(sprintf("Updated %d pixels in '%s'.", length(selected_idx), field_name), type = "message", duration = 5)
        
      } else if (identical(input$seg_rename, "add")) {
        if (!nzchar(field_name)) {
          showNotification("Enter a name for the new field.", type = "error")
          return()
        }
        if (field_name %in% colnames(pdat_df)) {
          showNotification("Field already exists. Use 'Edit values' to modify it.", type = "error")
          return()
        }
        if (!nzchar(new_value)) {
          showNotification("Enter a value to assign in the new field.", type = "error")
          return()
        }
        
        new_col <- rep(NA_character_, n_pixels)
        new_col[selected_idx] <- new_value
        pdat[[field_name]] <- new_col
        x2$pdat_anat <- pdat
        
        updated_cols <- colnames(as.data.frame(x2$pdat_anat))
        updateSelectInput(session, "seg_pdat_col", choices = updated_cols, selected = field_name)
        showNotification(sprintf("Created field '%s' with %d annotated pixels.", field_name, length(selected_idx)), type = "message", duration = 5)
      }
      
    })
    
    observeEvent(input$save_imzml, {
      
      req(input$save_imzml)
      x2 <- preproc_values()[["x2"]]
      if (is.null(x2$list_proc_img) || length(x2$list_proc_img) == 0) {
        showNotification("No processed data to save. Click 'Store processed data' first.", type = "error", duration = 8)
        return()
      }
      #create variables required for saving
      volumes <- c(wd = setup_values()[["wd"]], home = fs::path_home())
      shinyFiles::shinyFileSave(input, "save_imzml", roots = volumes, session = session)
      
      
      
      save_path <- shinyFiles::parseSavePath(volumes, input$save_imzml)
      if (nrow(save_path) == 0) return(NULL)
      
      
      filen <- as.character(save_path$datapath)
      #   
      # 
      # output$save_proc <- downloadHandler(
      #   filename = function() {
      #     if (input$seg_choice == "anat_seg") {
      #       fn = "/MSI-proc-anat_proc-"
      #     } else {
      #       fn = "/MSI-proc-bk_proc-"
      #     }
      #     
      #     paste0(getwd(), fn, Sys.Date(), ".rds")
      #   },
      #   content = function(file) {
      #pk_img <- x0$overview_peaks
      #pk_img <- Cardinal::combine(x2$list_proc_img)
      
      #issue with combining data, need to reorder with largest first
      pk_img_reorder <- x2$list_proc_img[order(sapply(x2$list_proc_img, ncol), decreasing = TRUE)]
      
      
      #use do.call to cbind the list
      pk_img<-combine_card(pk_img_reorder)
      #do.call(cbind, pk_img_reorder)
      
      
      #may need to restore the order based on the original data
      
      
      
      #cleanup feature data
      
      a<-fData(pk_img)[unique(colnames(fData(pk_img)))]
      
      fData(pk_img) <- MassDataFrame(mz=a$mz, a %>% as.data.frame() %>% dplyr::select(-mz))
      
      writeImzML(pk_img, filen)
      
      
      
      
    }
    )
    
    
    return(preproc_values)
    
  })
}
