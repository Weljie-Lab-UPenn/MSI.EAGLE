
  plot_card_UI<-function(id) {
    
    ns <- NS(id)
    
    col_choices<-c(hcl.pals())
    initial_cols<-c("Inferno", "Spectral", "Cividis", "Viridis", "Plasma", 
                    "Zissou 1", "Purple-Green", "Berlin", "PiYG", "Grays", 
                    "Batlow", "turku", "YlOrRd", "Terrain", 
                    "PrGn", "Green-Brown", "Hawaii", "Cork", "Rocket", "RdYlBu")
    
    col_choices<-c(initial_cols, setdiff(col_choices, initial_cols))
  
    tagList(  
      
        fluidRow(
          tags$head(
            #tags$style("label{font-family: BentonSans Book;}")
            #tags$style("label{font-size: 11px;} ")
          ),
          
          column(3, offset = 0,
                 radioButtons(ns("ion_viz3"), "Image visualization ions",
                                 c("First ion" = "viz_first", 
                                   "All ions (TIC)"="viz_all", 
                                    "Custom single / multiple"="custom")),
                 #numericInput(ns("mz_viz3"), "mz value for visualization",255.2),
                 uiOutput(ns("mz_viz3a")),
                 
                 fluidRow(
                   column(6, selectInput(ns("contrast3"), "Contrast enhancement", c( "none", "histogram",  "adaptive"))),
                   column(6, selectInput(ns("color3"), "Colorscale", col_choices, selected="Inferno"))
                 ),
                 fluidRow(
                   column(6, selectInput(ns("smooth3"), "Smoothing options", c("none", "mean", "gaussian", "bilateral", "adaptive", "diffusion", "guided"))),
                   column(6, checkboxInput(ns("normalize3"), "Scale multiple images?", value = FALSE))
                   ),
                 
                 fluidRow(
                   column(6, checkboxInput(ns("colorkey3"), "Draw color key?", value = TRUE)),
                   column(6, checkboxInput(ns("dark_bg"), "Dark background?", value = TRUE))
                   ),
                   
                
                 
                 fluidRow(
                   column(6, numericInput(ns("width_im"), "Image plot width (px)", value = 800, step = 50)),
                   column(6, numericInput(ns("height_im"), "Image plot height (px)", value=600, step = 50))
                 ),
                 checkboxInput(ns("plot_pdata"), "Plot Phenotype data?", value=FALSE),
                 uiOutput(ns("plotpdata")),
                 checkboxInput(ns("expand_fonts"), "Extended font options?", value=FALSE),
                 uiOutput(ns("fonts")),                 
                 checkboxInput(ns("expand_runs"), "Select individual runs for plotting only?", value=FALSE),
                 uiOutput(ns("select_runs")),
                 checkboxInput(ns("apply_mask"), "Apply pixel mask?", value=FALSE),
                 uiOutput(ns("mask_options")),
                 p("___________________________________________________________________")
                 

                 
          ),
          column(9, uiOutput(ns("plot.window"))
          )  
        ),
        
        uiOutput(ns("spectrum"))
    )
          
  }
  

  plot_card_server <- function(id, overview_peaks_sel, spatialOnly=FALSE) {
    
    moduleServer(id, function(input, output, session){
      
      #for dynamic UI
      ns = session$ns
      
      graphics.off()
      
      # NEW: RGB Overlay Function
      create_rgb_overlay <- function(msi_data, ions, tolerance, 
                                     channel_mapping = c("red", "green", "blue"),
                                     smooth = NULL, enhance = NULL, 
                                     normalize = FALSE, dark_bg = FALSE) {
        
        # Validate inputs
        if(length(ions) < 2) {
          showNotification("RGB overlay requires at least 2 ions", type = "warning")
          return(NULL)
        }
        
        if(length(ions) > 3) {
          showNotification("RGB overlay limited to 3 ions (using first 3)", type = "warning")
          ions <- ions[1:3]
        }
        
        # Get spatial coordinates - IMPORTANT: get from the actual (possibly masked) data
        coords <- coord(msi_data)
        
        # FIXED: Get the FULL coordinate space (not just unique from masked data)
        # This ensures we create a matrix that covers the entire spatial range
        x_min <- min(coords$x)
        x_max <- max(coords$x)
        y_min <- min(coords$y)
        y_max <- max(coords$y)
        
        # Create full coordinate grids
        all_x_coords <- x_min:x_max
        all_y_coords <- y_min:y_max
        
        width <- length(all_x_coords)
        height <- length(all_y_coords)
        
        # Initialize RGB array
        rgb_array <- array(0, dim = c(height, width, 3))
        
        # Process each ion
        for(i in seq_along(ions)) {
          # Extract ion image
          ion_data <- subsetFeatures(msi_data, mz = ions[i], tolerance = tolerance, units = "mz")
          
          if(nrow(ion_data) == 0) {
            showNotification(paste("No data found for m/z =", round(ions[i], 4)), type = "warning")
            next
          }
          
          # Get intensities
          intensities <- spectra(ion_data)[1,]
          
          # Get coordinates for this ion data
          ion_coords <- coord(ion_data)
          
          # Create 2D matrix matching spatial layout
          ion_matrix <- matrix(0, nrow = height, ncol = width)
          
          # FIXED: Map each pixel to its correct position in the matrix
          for(j in seq_along(intensities)) {
            # Get actual x,y coordinates for this pixel
            pixel_x <- ion_coords$x[j]
            pixel_y <- ion_coords$y[j]
            
            # Convert to matrix indices (1-based)
            x_idx <- pixel_x - x_min + 1
            y_idx <- pixel_y - y_min + 1
            
            # Store intensity at correct position
            ion_matrix[y_idx, x_idx] <- intensities[j]
          }
          
          # Apply enhancement FIRST (works better on raw intensities)
          if(!is.null(enhance) && enhance != "none") {
            ion_matrix <- apply_enhancement(ion_matrix, enhance)
          }
          
          # Apply smoothing SECOND
          if(!is.null(smooth) && smooth != "none") {
            ion_matrix <- apply_smoothing(ion_matrix, smooth)
          }
          
          # Normalize intensities
          if(normalize && max(ion_matrix) > 0) {
            # Normalize across all ions together for relative comparison
            ion_matrix <- ion_matrix / max(ion_matrix)
          } else if(max(ion_matrix) > 0) {
            # Individual normalization per ion
            ion_matrix <- (ion_matrix - min(ion_matrix)) / (max(ion_matrix) - min(ion_matrix))
          }
          
          # Map to RGB channel
          channel <- switch(channel_mapping[i],
                            "red" = 1,
                            "green" = 2, 
                            "blue" = 3,
                            "cyan" = 2:3,  # Cyan = Green + Blue
                            "magenta" = c(1,3),  # Magenta = Red + Blue
                            "yellow" = 1:2,  # Yellow = Red + Green
                            1)  # Default to red
          
          # Handle multi-channel colors (cyan, magenta, yellow)
          if(length(channel) > 1) {
            for(ch in channel) {
              rgb_array[,,ch] <- rgb_array[,,ch] + ion_matrix
            }
          } else {
            rgb_array[,,channel] <- ion_matrix
          }
        }
        
        # Clip values to [0,1]
        rgb_array[rgb_array > 1] <- 1
        rgb_array[rgb_array < 0] <- 0
        
        # Convert to raster
        raster_img <- as.raster(rgb_array)
        
        # Set plot parameters
        par(mar = c(4, 4, 3, 2))
        if(dark_bg) {
          par(bg = "black", col.main = "white", col.lab = "white", col.axis = "white")
        }
        
        # FIXED: Create empty plot with proper coordinate system
        plot(x = NULL, y = NULL, type = "n", 
             xlim = c(x_min, x_max), 
             ylim = c(y_min, y_max),
             xlab = "x", ylab = "y",
             main = paste("RGB Overlay:", paste(round(ions, 4), collapse = ", ")),
             asp = 1)
        
        # FIXED: Add the raster image using actual coordinate ranges with proper centering
        rasterImage(raster_img, 
                    xleft = x_min - 0.5, 
                    ybottom = y_min - 0.5, 
                    xright = x_max + 0.5, 
                    ytop = y_max + 0.5, 
                    interpolate = FALSE)
        
        # Add color legend
        if(length(ions) <= 3) {
          legend_colors <- channel_mapping[1:length(ions)]
          legend_text <- paste0("m/z ", round(ions, 4))
          legend("topright", legend = legend_text, 
                 fill = legend_colors, 
                 bg = if(dark_bg) "black" else "white",
                 text.col = if(dark_bg) "white" else "black",
                 cex = 1.5,           # Increased from 0.8 to 1.5
                 pt.cex = 2,          # Makes the color boxes larger
                 box.lwd = 2)         # Makes the legend box border thicker
        }
        
        return(invisible(rgb_array))
      }
      
      
      # Helper function for smoothing
      apply_smoothing <- function(mat, method) {
        # Ensure matrix is numeric
        mat <- as.matrix(mat)
        
        if(method == "gaussian") {
          # Gaussian smoothing using a 3x3 kernel
          kernel <- matrix(c(1, 2, 1, 2, 4, 2, 1, 2, 1), nrow = 3) / 16
          mat <- apply_kernel(mat, kernel)
          
        } else if(method == "mean") {
          # Mean filter (3x3)
          kernel <- matrix(1, nrow = 3, ncol = 3) / 9
          mat <- apply_kernel(mat, kernel)
          
        } else if(method == "bilateral") {
          # Simplified bilateral filter
          # This preserves edges while smoothing
          mat <- bilateral_filter(mat, sigma_spatial = 2, sigma_range = 0.1)
          
        } else if(method == "adaptive") {
          # Adaptive smoothing based on local variance
          mat <- adaptive_smooth(mat)
          
        } else if(method == "diffusion") {
          # Anisotropic diffusion (edge-preserving)
          mat <- anisotropic_diffusion(mat, iterations = 5)
          
        } else if(method == "guided") {
          # Guided filter (edge-preserving)
          mat <- guided_filter(mat, radius = 3, epsilon = 0.01)
        }
        
        return(mat)
      }
      
      # Helper function for enhancement
      apply_enhancement <- function(mat, method) {
        # Ensure matrix is numeric
        mat <- as.matrix(mat)
        
        if(method == "histogram") {
          # Histogram equalization
          mat <- histogram_equalization(mat)
          
        } else if(method == "adaptive") {
          # Adaptive histogram equalization (CLAHE)
          mat <- adaptive_histogram_equalization(mat, clip_limit = 0.03)
        }
        
        return(mat)
      }
      
      # Kernel convolution helper
      apply_kernel <- function(mat, kernel) {
        nrow_mat <- nrow(mat)
        ncol_mat <- ncol(mat)
        nrow_ker <- nrow(kernel)
        ncol_ker <- ncol(kernel)
        
        # Pad the matrix
        pad_row <- floor(nrow_ker / 2)
        pad_col <- floor(ncol_ker / 2)
        
        padded <- matrix(0, nrow = nrow_mat + 2 * pad_row, ncol = ncol_mat + 2 * pad_col)
        padded[(pad_row + 1):(pad_row + nrow_mat), (pad_col + 1):(pad_col + ncol_mat)] <- mat
        
        # Apply convolution
        result <- mat
        for(i in 1:nrow_mat) {
          for(j in 1:ncol_mat) {
            window <- padded[i:(i + nrow_ker - 1), j:(j + ncol_ker - 1)]
            result[i, j] <- sum(window * kernel)
          }
        }
        
        return(result)
      }
      
      # Bilateral filter
      bilateral_filter <- function(mat, sigma_spatial = 2, sigma_range = 0.1, window_size = 5) {
        nrow_mat <- nrow(mat)
        ncol_mat <- ncol(mat)
        result <- mat
        half_win <- floor(window_size / 2)
        
        for(i in 1:nrow_mat) {
          for(j in 1:ncol_mat) {
            # Get window bounds
            i_min <- max(1, i - half_win)
            i_max <- min(nrow_mat, i + half_win)
            j_min <- max(1, j - half_win)
            j_max <- min(ncol_mat, j + half_win)
            
            # Extract window
            window <- mat[i_min:i_max, j_min:j_max]
            center_val <- mat[i, j]
            
            # Calculate spatial weights
            rows <- i_min:i_max - i
            cols <- j_min:j_max - j
            spatial_weights <- outer(rows^2, cols^2, "+")
            spatial_weights <- exp(-spatial_weights / (2 * sigma_spatial^2))
            
            # Calculate range weights
            range_weights <- exp(-(window - center_val)^2 / (2 * sigma_range^2))
            
            # Combined weights
            weights <- spatial_weights * range_weights
            weights <- weights / sum(weights)
            
            # Apply filter
            result[i, j] <- sum(window * weights)
          }
        }
        
        return(result)
      }
      
      # Adaptive smoothing
      adaptive_smooth <- function(mat, window_size = 5) {
        nrow_mat <- nrow(mat)
        ncol_mat <- ncol(mat)
        result <- mat
        half_win <- floor(window_size / 2)
        
        for(i in 1:nrow_mat) {
          for(j in 1:ncol_mat) {
            i_min <- max(1, i - half_win)
            i_max <- min(nrow_mat, i + half_win)
            j_min <- max(1, j - half_win)
            j_max <- min(ncol_mat, j + half_win)
            
            window <- mat[i_min:i_max, j_min:j_max]
            
            # Adapt based on local variance
            local_var <- var(as.vector(window))
            if(local_var < 0.01) {
              # Low variance - smooth more
              result[i, j] <- mean(window)
            } else {
              # High variance - preserve detail
              result[i, j] <- mat[i, j]
            }
          }
        }
        
        return(result)
      }
      
      # Anisotropic diffusion
      anisotropic_diffusion <- function(mat, iterations = 5, kappa = 0.1, lambda = 0.25) {
        result <- mat
        
        for(iter in 1:iterations) {
          # Calculate gradients
          grad_n <- rbind(result[1,], result[1:(nrow(result)-1),]) - result
          grad_s <- rbind(result[2:nrow(result),], result[nrow(result),]) - result
          grad_e <- cbind(result[,2:ncol(result)], result[,ncol(result)]) - result
          grad_w <- cbind(result[,1], result[,1:(ncol(result)-1)]) - result
          
          # Calculate diffusion coefficients (edge-preserving)
          c_n <- exp(-(grad_n / kappa)^2)
          c_s <- exp(-(grad_s / kappa)^2)
          c_e <- exp(-(grad_e / kappa)^2)
          c_w <- exp(-(grad_w / kappa)^2)
          
          # Update
          result <- result + lambda * (c_n * grad_n + c_s * grad_s + c_e * grad_e + c_w * grad_w)
        }
        
        return(result)
      }
      
      # Guided filter
      guided_filter <- function(mat, radius = 3, epsilon = 0.01) {
        # Use the input as both guide and input
        guide <- mat
        input <- mat
        
        nrow_mat <- nrow(mat)
        ncol_mat <- ncol(mat)
        result <- matrix(0, nrow = nrow_mat, ncol = ncol_mat)
        
        for(i in 1:nrow_mat) {
          for(j in 1:ncol_mat) {
            i_min <- max(1, i - radius)
            i_max <- min(nrow_mat, i + radius)
            j_min <- max(1, j - radius)
            j_max <- min(ncol_mat, j + radius)
            
            # Local window
            guide_window <- guide[i_min:i_max, j_min:j_max]
            input_window <- input[i_min:i_max, j_min:j_max]
            
            # Calculate local statistics
            mean_guide <- mean(guide_window)
            mean_input <- mean(input_window)
            var_guide <- var(as.vector(guide_window))
            cov_guide_input <- mean((guide_window - mean_guide) * (input_window - mean_input))
            
            # Calculate a and b
            a <- cov_guide_input / (var_guide + epsilon)
            b <- mean_input - a * mean_guide
            
            # Apply
            result[i, j] <- a * guide[i, j] + b
          }
        }
        
        return(result)
      }
      
      # Histogram equalization
      histogram_equalization <- function(mat) {
        # Normalize to 0-1 range
        mat_norm <- (mat - min(mat)) / (max(mat) - min(mat))
        
        # Calculate histogram
        breaks <- seq(0, 1, length.out = 256)
        hist_data <- hist(mat_norm, breaks = breaks, plot = FALSE)
        
        # Calculate CDF
        cdf <- cumsum(hist_data$counts) / sum(hist_data$counts)
        
        # Map values
        result <- mat_norm
        for(i in 1:length(breaks[-1])) {
          mask <- mat_norm >= breaks[i] & mat_norm < breaks[i + 1]
          result[mask] <- cdf[i]
        }
        
        return(result)
      }
      
      # Adaptive histogram equalization (CLAHE)
      adaptive_histogram_equalization <- function(mat, clip_limit = 0.03, tile_size = 8) {
        nrow_mat <- nrow(mat)
        ncol_mat <- ncol(mat)
        
        # Determine number of tiles
        n_tiles_row <- ceiling(nrow_mat / tile_size)
        n_tiles_col <- ceiling(ncol_mat / tile_size)
        
        result <- mat
        
        # Process each tile
        for(tile_i in 1:n_tiles_row) {
          for(tile_j in 1:n_tiles_col) {
            # Get tile boundaries
            row_start <- (tile_i - 1) * tile_size + 1
            row_end <- min(tile_i * tile_size, nrow_mat)
            col_start <- (tile_j - 1) * tile_size + 1
            col_end <- min(tile_j * tile_size, ncol_mat)
            
            # Extract tile
            tile <- mat[row_start:row_end, col_start:col_end]
            
            # Normalize
            tile_norm <- (tile - min(tile)) / (max(tile) - min(tile) + 1e-10)
            
            # Calculate histogram with clipping
            hist_data <- hist(tile_norm, breaks = seq(0, 1, length.out = 50), plot = FALSE)
            counts <- hist_data$counts
            
            # Clip histogram
            clip_value <- clip_limit * length(tile_norm) / length(counts)
            clipped <- pmin(counts, clip_value)
            
            # Redistribute clipped values
            excess <- sum(counts - clipped)
            redistribution <- excess / length(counts)
            clipped <- clipped + redistribution
            
            # Calculate CDF
            cdf <- cumsum(clipped) / sum(clipped)
            
            # Map values
            tile_result <- tile_norm
            breaks <- hist_data$breaks
            for(i in 1:(length(breaks) - 1)) {
              mask <- tile_norm >= breaks[i] & tile_norm < breaks[i + 1]
              tile_result[mask] <- cdf[i]
            }
            
            # Store result
            result[row_start:row_end, col_start:col_end] <- tile_result
          }
        }
        
        return(result)
      }

      #create new overview_peaks_sel object with mean values
      if(is.null(fData(overview_peaks_sel)$mean)) {
        overview_peaks_sel<-Cardinal::summarizeFeatures(overview_peaks_sel, verbose=F)
        if(class(overview_peaks_sel)=="try-error") {
          showNotification("No data available, please check your parameters or dataset", type="error")
          return()
        }
      }

      # NEW: Add UI for pixel masking options
      observe({
        output$mask_options <- renderUI({
          req(overview_peaks_sel)
          
          if(input$apply_mask) {
            # Get available pData variables (exclude coordinate columns)
            pdata_vars <- colnames(as.data.frame(pData(overview_peaks_sel)))
            pdata_vars <- pdata_vars[!pdata_vars %in% c("x", "y", "run")]
            
            tagList(
              fluidRow(
                column(6, selectInput(ns("mask_variable"), 
                                      label = "pData variable for masking", 
                                      choices = pdata_vars)),
                column(6, uiOutput(ns("mask_values")))
              ),
              fluidRow(
                column(6, checkboxInput(ns("invert_mask"), "Invert mask?", value = FALSE)),
                column(6, checkboxInput(ns("show_mask_info"), "Show mask info?", value = TRUE))
              )
            )
          }
        })
      })
      
      # NEW: Dynamic UI for selecting mask values
      observe({
        output$mask_values <- renderUI({
          req(input$mask_variable)
          req(overview_peaks_sel)
          
          # Get unique values for the selected pData variable
          pdata_df <- as.data.frame(pData(overview_peaks_sel))
          unique_values <- unique(pdata_df[[input$mask_variable]])
          unique_values <- unique_values[!is.na(unique_values)]
          
          if(is.numeric(unique_values)) {
            # For numeric variables, provide a range slider
            tagList(
              sliderInput(ns("mask_range"), 
                          label = "Value range to include",
                          min = min(unique_values, na.rm = TRUE),
                          max = max(unique_values, na.rm = TRUE),
                          value = c(min(unique_values, na.rm = TRUE), 
                                    max(unique_values, na.rm = TRUE)),
                          step = (max(unique_values, na.rm = TRUE) - min(unique_values, na.rm = TRUE))/20)
            )
          } else {
            # For categorical variables, provide checkboxes
            selectInput(ns("mask_categories"), 
                        label = "Categories to include",
                        choices = unique_values,
                        selected = unique_values,
                        multiple = TRUE)
          }
        })
      })
      
      # NEW: Function to apply pixel mask
      apply_pixel_mask <- function(data) {
        if(!input$apply_mask) {
          return(data)
        }
        
        req(input$mask_variable)
        
        pdata_df <- as.data.frame(pData(data))
        mask_var <- pdata_df[[input$mask_variable]]
        
        if(is.numeric(mask_var)) {
          req(input$mask_range)
          # Create mask for numeric variables
          pixel_mask <- mask_var >= input$mask_range[1] & mask_var <= input$mask_range[2]
        } else {
          req(input$mask_categories)
          # Create mask for categorical variables
          pixel_mask <- mask_var %in% input$mask_categories
        }
        
        # Remove NA values from mask
        pixel_mask[is.na(pixel_mask)] <- FALSE
        
        # Invert mask if requested
        if(input$invert_mask) {
          pixel_mask <- !pixel_mask
        }
        
        # Show mask information if requested
        if(input$show_mask_info) {
          n_total <- length(pixel_mask)
          n_selected <- sum(pixel_mask, na.rm = TRUE)
          showNotification(
            paste("Pixel mask applied:", n_selected, "of", n_total, "pixels selected"),
            type = "message", duration = 3
          )
        }
        
        # Apply mask using subsetPixels
        if(sum(pixel_mask, na.rm = TRUE) == 0) {
          showNotification("No pixels selected with current mask settings!", type = "warning")
          return(data)
        }
        
        # Create logical vector for pixel selection
        pixel_indices <- which(pixel_mask)
        masked_data <- subsetPixels(data, pixel_indices)
        
        return(masked_data)
      }
 
      observe({
        
        output$plot.window <- renderUI({
          
          if(input$ion_viz3=="custom") {
            fluidRow(
              column(12, imageOutput(ns("plot3_pk"))),
              fluidRow(
                column(12, br()),  # Adding space between the image and the table
                column(12, DT::dataTableOutput(ns('tbl')))
              )
            )
          } else {
            fluidRow(
              column(12, imageOutput(ns("plot3_pk")))
            )
          }
          
        })
        
      })
      
      

      
      output$mz_viz3a <- renderUI ({
        req(overview_peaks_sel)
        
        if(input$ion_viz3=="custom") {
          
          mz_list<- mz(overview_peaks_sel)
          tbl<-as.data.frame(fData(overview_peaks_sel))
          
          # Round existing columns as before
          if(!is.null(tbl$freq)) {
            tbl$freq<-round(tbl$freq, 2)
          }
          if(!is.null(tbl$mean)) {
            tbl$mean<-round(tbl$mean, 1)
          }
          tbl$mz<-round(tbl$mz, 4)
          
          # Add mass error in ppm if available (for targeted analysis)
          if(!is.null(tbl$mass_diff_ppm)) {
            tbl$mass_error_ppm <- round(as.numeric(tbl$mass_diff_ppm), 2)
            tbl$mass_diff_ppm <- NULL
            
            col_order <- names(tbl)
            if("mass_error_ppm" %in% col_order) {
              mz_pos <- which(col_order == "mz")
              if(length(mz_pos) > 0) {
                new_order <- c(col_order[1:mz_pos], 
                               "mass_error_ppm", 
                               col_order[!col_order %in% c(names(tbl)[1:mz_pos], "mass_error_ppm")])
                tbl <- tbl[, new_order]
              }
            }
          }
          
          updateNumericInput(session, ("width_im"), value=800, step=50)
          updateNumericInput(session, ("height_im"), value=450, step=50)
          
          output$tbl <-DT::renderDataTable({
            tbl
          }, selection = "multiple")
          
          list(
            checkboxInput(ns("superpose"), "Superpose images?", value=FALSE),
            
            # NEW: Conditional UI for RGB overlay
            conditionalPanel(
              condition = "input.superpose == true",
              ns = ns,
              radioButtons(ns("overlay_method"), "Overlay method:",
                           choices = c("Cardinal default (single palette)" = "cardinal",
                                       "RGB overlay (multi-color)" = "rgb"),
                           selected = "cardinal"),
              
              conditionalPanel(
                condition = "input.overlay_method == 'rgb'",
                ns = ns,
                selectInput(ns("rgb_mapping"), "Channel assignment:",
                            choices = list(
                              "Ion 1: Red, Ion 2: Green, Ion 3: Blue" = "RGB",
                              "Ion 1: Green, Ion 2: Red, Ion 3: Blue" = "GRB",
                              "Ion 1: Blue, Ion 2: Green, Ion 3: Red" = "BGR",
                              "Ion 1: Red, Ion 2: Green, Ion 3: Blue (norm together)" = "RGB_norm",
                              "Ion 1: Cyan, Ion 2: Magenta, Ion 3: Yellow" = "CMY"
                            ),
                            selected = "RGB"),
                helpText("Note: RGB overlay uses custom rendering. Smoothing/enhancement applied per channel.")
              )
            ),
            
            selectInput(ns("display_mode"), "Ion math?", 
                        c("none", "sum", "ratio", "subtract", "min", "max", "mean", "sd", "var", "multiply")),
            numericInput(ns("plusminus_viz3"), "+/- m/z for visualization", 0.05)
          )
        }
      })
         
      observe({
        req(input$ion_viz3 == "custom")
        req(input$superpose)
        
        if(!is.null(input$overlay_method) && input$overlay_method == "rgb") {
          if(is.null(mz_viz3()) || length(mz_viz3()) < 2) {
            showNotification(
              "Please select 2-3 ions from the table for RGB overlay",
              type = "message",
              duration = 5
            )
          } else if(length(mz_viz3()) > 3) {
            showNotification(
              "RGB overlay uses first 3 selected ions only",
              type = "warning",
              duration = 5
            )
          }
        }
      })   
    
      
      mz_viz3 <- reactive({
        req(overview_peaks_sel)
        req(input$tbl_rows_selected)
        tbl<-as.data.frame(fData(overview_peaks_sel))
        return(tbl[input$tbl_rows_selected, "mz"])
        
      })

      

      observe({
        
        output$select_runs <- renderUI ({
          req(overview_peaks_sel)
          
          if(input$expand_runs) {
            
            run_list<- unique(run((overview_peaks_sel)))
            
            list(
              selectInput(ns("select_runs"),
                          label= "run selection (plot only)",
                          multiple=TRUE,
                          choices = run_list,
                          selected = run_list)
            )
          }
        })
        
        
      })
      
      observe({
        output$fonts <- renderUI ({
          
          if(input$expand_fonts){
            list(
              numericInput(ns("axis_size"), label = "Axis font scaling (%)", min = 0, value=100),
              numericInput(ns("title_size"), label = "Title font scaling (%)", min = 0, value=100),
              numericInput(ns("label_size"), label = "Label font scaling (%)", min = 0, value=100),
              numericInput(ns("subtitle_size"), label = "Subtitle font scaling (%)", min = 0, value=100)
            )
          }
          
          
        })
      })
      
      observe({
        output$plotpdata <- renderUI ({
          
          if(input$plot_pdata){
            list(
              selectInput(ns("pdata_var_plot"), label="pData variable to plot", 
                          choices = colnames(as.data.frame(pData(overview_peaks_sel)))[-c(1:2)]
                          )
            )
          }
          
          
        })
      })

       #add spectrum or not
      output$spectrum<-renderUI({
        
        if(input$ion_viz3=="custom") {
          DT::dataTableOutput(ns('tbl'))
        }

        
        if(spatialOnly==FALSE){
          tagList(
            fluidRow(
              column(6,uiOutput(ns("plot_ranges"))),
              column(6,uiOutput(ns("plot_ranges2"))),
            ),
           
            fluidRow(
              
              column(2,  uiOutput(ns("x_target"))),
              column(2,  numericInput(ns("x_tol"), "+/- m/z window", value=5)),
              column(2,  uiOutput(ns("slider_y_max"))),
              column(3,numericInput(ns("width_ms"), "Plot width (px)", value = 800)),
              column(3,numericInput(ns("height_ms"), "Plot height (px)", value=300))
            ),
            fluidRow(
              column(10, imageOutput(ns("plot4")), style = "margin-bottom: 0px; padding-bottom: 0px;")
              ),
            fluidRow(
              column(3, checkboxInput(ns("calc_ppm"), label = "Show ppm error?", width = '100%')),
              column(3, checkboxInput(ns("show_int"), label = "Show intensity?", width = '100%')),
              column(3, numericInput(ns("show_mz"), label = "# mz values to show?", width = '100%', value=0))
            ),
            checkboxInput(ns("spectrum_expand_fonts"), "Extended font options?", value=FALSE),
            uiOutput(ns("spectrum_fonts")),
          )
        } 
      })
      
      observe({
        output$spectrum_fonts <- renderUI ({
          
          if(input$spectrum_expand_fonts) {
            fluidRow(
              column(
                3,
                numericInput(
                  ns("spectrum_axis_size"),
                  label = "Axis font scaling (%)",
                  min = 0,
                  value = 100
                )
              ),
              # column(
              #   2,
              #   numericInput(
              #     ns("spectrum_title_size"),
              #     label = "Title font scaling (%)",
              #     min = 0,
              #     value = 100
              #   )
              # ),
              column(
                3,
                numericInput(
                  ns("spectrum_label_size"),
                  label = "Label font scaling (%)",
                  min = 0,
                  value = 100
                )
              ),
              # column(
              #   2,
              #   numericInput(
              #     ns("spectrum_subtitle_size"),
              #     label = "Subtitle font scaling (%)",
              #     min = 0,
              #     value = 100
              #   )
              # ),
              column(
                3,
                numericInput(
                  ns("linewidth"),
                  label="Linewidth",
                  min=0,
                  value=100
                )
              )
            )
            
          }
          
          
        })
      })
          
       observe({
        output$plot3_pk <- renderImage( {  #plot image in overview after peakpicking / reading file
          
          #req(overview_peaks_sel)
          
          # A temp file to save the output.
          # This file will be removed later by renderImage
          outfile <- tempfile(fileext = '.png')
          
          png(outfile, width = input$width_im, height = input$height_im)
          
          #ion <- switch(input$mode,
          #              "p"=786,
          #             "n"=255.2)
          
          
          if(!is.null(input$select_runs)) {
            overview_peaks_sel <- subsetPixels(overview_peaks_sel, run %in% input$select_runs)
          }
          
          # Apply pixel mask if enabled
          overview_peaks_sel <- apply_pixel_mask(overview_peaks_sel)
          
          vp_orig<-vizi_par()
          
          #set sizes
          if(input$expand_fonts) {
            req(input$axis_size)
            
            vizi_par(
              cex.axis=input$axis_size/100,
              cex.lab=input$label_size/100,
              cex.main=input$title_size/100,
              cex.sub=input$subtitle_size/100
            )
            #get margins
            # cur_mar<-par()$mar
            # 
            # new_mar<-c(cur_mar[1]+cex.labp/8, cur_mar[2]+cex.labp/2, cur_mar[3]+cex.mainp/2, cur_mar[4])
            # 
            # #if drawing colorkey, add a little on the left
            # if(input$colorkey3) {
            #   new_mar[4]=new_mar[4]+cex.axisp
            # }
            # 
            # cur_mgp<-par()$mgp
            # 
            # #new_mgp<-c(cur_mgp[1]+cex.labp/8, cur_mgp[2], 0)
            # new_mgp<-c(3+max(new_mar[1:2])/20, 1,0)
            # 
            
            
            
            
           } else {
              vizi_par(
                cex.axis=1,
                cex.lab=1,
                cex.main=1,
                cex.sub=1
              )

           }
           
          if(input$plot_pdata){
            req(input$pdata_var_plot)
            
            #create list of arguments for image
            arg_list<-list(overview_peaks_sel, 
                       input$pdata_var_plot,
                        key=(input$colorkey3),
                        col=pals::alphabet())
                        
            if(input$dark_bg) {
              arg_list$style <- "dark"
            }
            
            plt_tmp<-do.call(Cardinal::image, arg_list)
            
            
            # Cardinal::image(overview_peaks_sel,
            #                         input$pdata_var_plot,
            #                         key=(input$colorkey3),
            #                         #superpose=input$superpose,
            #                         col=pals::alphabet())
            print(plt_tmp,
                                  #cex.axis=req(cex.axisp),
                                  #cex.lab=cex.labp,
                                  #cex.main=cex.mainp,
                                  #cex.sub=cex.subp,
                                  #mar=new_mar,
                                  #mgp=new_mgp
                  )
            
            vizi_par(vp_orig)
          } else if (input$ion_viz3=="viz_all") {
            
            
            mz_range=range(mz(overview_peaks_sel))
            
            #find closest mz value to middle of range
            test_value <- mean(mz_range)
            
            
            # Calculate the absolute differences
            differences <- abs(mz(overview_peaks_sel) - test_value)
            
            # Find the index of the minimum difference
            closest_index <- which.min(differences)
            
            mz_set=mz(overview_peaks_sel)[closest_index]
            tol=max(differences) + differences[closest_index]+1
            
            plusminus=tol
            
            #old way-- may be more memory efficient?
            # smoothing_option <- if (input$smooth3 != "none") paste0(", smooth ='", input$smooth3,"'") else ""
            # enhance_option <- if (input$contrast3 != "none") paste0(", enhance ='", input$contrast3,"'") else ""
            # 
            # image_command <- paste("Cardinal::image(overview_peaks_sel, 
            #                       mz=mz_set,
            #                       tolerance=round(plusminus,3), 
            #                       units='mz',
            #                       col=cpal(input$color3)",
            #                       enhance_option,
            #                       smoothing_option,",
            #                       scale=input$normalize3,
            #                       #superpose=input$superpose,
            #                       key=(input$colorkey3),
            #                       #cex.axis=req(cex.axisp),
            #                       #cex.lab=cex.labp,
            #                       #cex.main=cex.mainp,
            #                       #cex.sub=cex.subp,
            #                       #mar=new_mar,
            #                       #mgp=new_mgp
            # )")
            # 
            
            #print(eval(parse(text = image_command)))
            
            
            
            smoothing_option <- if (input$smooth3 != "none")  input$smooth3 else NULL
            enhance_option <- if (input$contrast3 != "none")  input$contrast3 else NULL
            
            
            arg_list<-list(overview_peaks_sel, 
                           mz=mz_set,
                           tolerance=round(plusminus,3), 
                           units='mz',
                           col=cpal(input$color3),
                            enhance=enhance_option,
                           smooth=smoothing_option,
                           scale=input$normalize3,
                           #superpose=input$superpose,
                           key=(input$colorkey3))
            
            if(input$dark_bg) {
              arg_list$style <- "dark"
            }
            
            print(do.call(Cardinal::image, arg_list))
            
            
            
            

            vizi_par(vp_orig)
            
            
          } else if (input$ion_viz3=="custom"){
            
            if(is.null(mz_viz3())){
              ion=mz(overview_peaks_sel[1,])
            } else {
              ion=as.numeric(mz_viz3())
            }
            
            # NEW: Check if using RGB overlay
            if(input$superpose && !is.null(input$overlay_method) && 
               input$overlay_method == "rgb" && length(ion) >= 2) {
              
              # Parse RGB mapping
              mapping_choice <- if(is.null(input$rgb_mapping)) "RGB" else input$rgb_mapping
              
              channel_map <- switch(mapping_choice,
                                    "RGB" = c("red", "green", "blue"),
                                    "GRB" = c("green", "red", "blue"),
                                    "BGR" = c("blue", "green", "red"),
                                    "RGB_norm" = c("red", "green", "blue"),
                                    "CMY" = c("cyan", "magenta", "yellow"),
                                    c("red", "green", "blue")  # default
              )
              
              # Determine if using collective normalization
              normalize_together <- grepl("norm", mapping_choice)
              
              smoothing_option <- if (input$smooth3 != "none") input$smooth3 else NULL
              enhance_option <- if (input$contrast3 != "none") input$contrast3 else NULL
              
              # Create RGB overlay
              create_rgb_overlay(
                msi_data = overview_peaks_sel,
                ions = ion,
                tolerance = input$plusminus_viz3,
                channel_mapping = channel_map,
                smooth = smoothing_option,
                enhance = enhance_option,
                normalize = normalize_together,
                dark_bg = input$dark_bg
              )
              
              vizi_par(vp_orig)
              
            } else if(!is.null(input$display_mode) && input$display_mode!="none"){
              if(input$display_mode%in%c("min", "max", "min", "mean", "sum", "sd", "var")) { 
          
                
                
                select_vec<-as.character(mz(overview_peaks_sel)) %in% as.character(ion)
                #test to make sure there are 2 or more elements
                if(sum(select_vec)<2){
                  showNotification("At least two ions required for this calculation", type="error")
                  message("At least two ions required for this calculation")
                  return()
                }
                
                sm<-summarizePixels(overview_peaks_sel[select_vec,], stat=c(xic=input$display_mode), as="DataFrame")
                pData(overview_peaks_sel)$xic<-sm$xic
                
                label_txt=paste(input$display_mode, "mz(s)=", paste(ion, collapse=", "))
            
              }else if(input$display_mode=="ratio"){
                if(length(ion)!=2){
                  showNotification("Exactly two ions required for ratio (mz1/mz2)", type="error")
                  message("Exactly two ions required for ratio (mz1/mz2)")
                  return()
                } else {
                  mz1 <- spectra(subsetFeatures(overview_peaks_sel, mz=ion[1]))[1,]
                  mz2 <- spectra(subsetFeatures(overview_peaks_sel, mz=ion[2]))[1,]
                  overview_peaks_sel$xic <- (1 + mz1) / (1 + mz2)
                  
                  ion=round(ion, 4)
                  label_txt=paste("ratio of",ion[1],"/",ion[2])
                }
                
              } else if(input$display_mode=="subtract") {
                if(length(ion)!=2){
                  showNotification("Exactly two ions required for subtraction (mz1-mz2)", type="error")
                  message("Exactly two ions required for subtraction (mz1-mz2)")
                  return()
                } else {
                  mz1 <- spectra(subsetFeatures(overview_peaks_sel, mz=ion[1]))[1,]
                  mz2 <- spectra(subsetFeatures(overview_peaks_sel, mz=ion[2]))[1,]
                  overview_peaks_sel$xic <- (mz1) - (mz2)
                  ion=round(ion, 4)
                  label_txt=paste("difference of",ion[1],"-",ion[2])
                  
                }
                
                
              }else if(input$display_mode=="multiply"){
                nelements=length(ion)
                xic <- 1+spectra(subsetFeatures(overview_peaks_sel, mz=ion[1]))[1,]
                
                for(i in 2:nelements){
                  mz2= 1+spectra(subsetFeatures(overview_peaks_sel, mz=ion[i]))[1,]
                  overview_peaks_sel$xic=(xic) * (mz2)
                  ion=round(ion, 4)
                  label_txt=paste(ion[1],"*",ion[2])
                  }
                }
            
  
              
              
              if(sum(is.na(pData(overview_peaks_sel)$xic))==length(overview_peaks_sel)) {
                showNotification("This calculation does not work!")
                return()
              }
              
              
              
              # smoothing_option <- if (input$smooth3 != "none") paste0(", smooth ='", input$smooth3,"'") else ""
              # enhance_option <- if (input$contrast3 != "none") paste0(", enhance ='", input$contrast3,"'") else ""
              # 
              plusminus=input$plusminus_viz3
              
            #   image_command <- paste("Cardinal::image(overview_peaks_sel, 'xic',
            #                       tolerance=round(plusminus,3), 
            #                       units='mz',
            #                       col=cpal(input$color3)",
            #                          enhance_option,
            #                          smoothing_option,",
            #                       scale=input$normalize3,
            #                       #superpose=input$superpose,
            #                       key=(input$colorkey3),
            #                       #cex.axis=req(cex.axisp),
            #                       #cex.lab=cex.labp,
            #                       #cex.main=cex.mainp,
            #                       #cex.sub=cex.subp,
            #                       #mar=new_mar,
            #                       #mgp=new_mgp
            # )")
            #   
              smoothing_option <- if (input$smooth3 != "none")  input$smooth3 else NULL
              enhance_option <- if (input$contrast3 != "none")  input$contrast3 else NULL
              
              
              arg_list<-list(overview_peaks_sel,
                             'xic',
                             tolerance=round(plusminus,3), 
                             units='mz',
                             col=cpal(input$color3),
                             enhance=enhance_option,
                             smooth=smoothing_option,
                             scale=input$normalize3,
                             #superpose=input$superpose,
                             key=(input$colorkey3))
              
              if(input$dark_bg) {
                arg_list$style <- "dark"
              }
              
              print(matter::as_facets(do.call(Cardinal::image, arg_list), labels=label_txt))
              
              vizi_par(vp_orig)

            } else {
              
              
              
              mz_set=ion
              
              mz_range=range(mz(overview_peaks_sel))
              
              
              # smoothing_option <- if (input$smooth3 != "none") paste0(", smoothing ='", input$smooth3,"'") else ""
              # enhance_option <- if (input$contrast3 != "none") paste0(", enhance ='", input$contrast3,"'") else ""
              # 
              plusminus=input$plusminus_viz3
              
            #   image_command <- paste("Cardinal::image(overview_peaks_sel, 
            #                       mz=mz_set,
            #                       tolerance=round(plusminus,3), 
            #                       units='mz',
            #                       col=cpal(input$color3)",
            #                          enhance_option,
            #                          smoothing_option,",
            #                       scale=input$normalize3,
            #                       superpose=input$superpose,
            #                       key=(input$colorkey3),
            #                       #cex.axis=req(cex.axisp),
            #                       #cex.lab=cex.labp,
            #                       #cex.main=cex.mainp,
            #                       #cex.sub=cex.subp,
            #                       #mar=new_mar,
            #                       #mgp=new_mgp
            # )")
            #   
              
              
              smoothing_option <- if (input$smooth3 != "none")  input$smooth3 else NULL
              enhance_option <- if (input$contrast3 != "none")  input$contrast3 else NULL
              
              
              arg_list<-list(overview_peaks_sel,
                             #'xic',
                             mz=mz_set,
                             tolerance=round(plusminus,3), 
                             units='mz',
                             col=cpal(input$color3),
                             enhance=enhance_option,
                             smooth=smoothing_option,
                             scale=input$normalize3,
                             superpose=input$superpose,
                             key=(input$colorkey3))
              
              if(input$dark_bg) {
                arg_list$style <- "dark"
              }
              
              
              
              
              print(do.call(Cardinal::image, arg_list))
              
              vizi_par(vp_orig)

              
            }
          } else if (input$ion_viz3=="viz_first") {
            
            
            tol=0.05
            #browser()
            image_command <-Cardinal::image(overview_peaks_sel, 
                                  col=cpal(input$color3),
                                  #enhance=input$contrast3,
                                  #smooth=input$smooth3,
                                  scale=input$normalize3,
                                  #superpose=input$superpose,
                                  key=(input$colorkey3),
                                  #cex.axis=req(cex.axisp),
                                  #cex.lab=cex.labp,
                                  #cex.main=cex.mainp,
                                  #cex.sub=cex.subp,
                                  #mar=new_mar,
                                  #mgp=new_mgp
            )
            
            print(image_command)
            
            vizi_par(vp_orig)
            
          }
          
          
          dev.off()
          
          
          # Return a list containing the filename
          list(src = outfile,
               contentType = 'image/png',
               width = input$width_im,
               height = input$height_im,
               alt = "This is alternate text")
        }, deleteFile = TRUE)
       }) 
       
       
       if(spatialOnly==FALSE) {
         observe({  
          output$plot_ranges<- renderUI( {
            
            req(overview_peaks_sel)
            
            a<-Cardinal::plot(overview_peaks_sel)
            
  
                  sliderInput(ns("mass_range_plot"), 
                              label = p("m/z range for MS plot (X)"), 
                              min = round(a$channels$x$limits[1]),
                              max = round(a$channels$x$limits[2]), 
                              value = round(a$channels$x$limits), 
                              step = NULL,
                              round = TRUE)
  
            
          })
         })
         
         observe({
          output$plot_ranges2<- renderUI( {
            req(overview_peaks_sel)
            #browser()
            #overview_peaks_sel<-Cardinal::summarizeFeatures(overview_peaks_sel)
            a<-Cardinal::plot(overview_peaks_sel, "mean")
            
            
           
            
            sliderInput (ns("int_range_plot"), 
                        label = p("intensity range for MS plot (Y)"), 
                        min = round(a$channels$y$limits[1],0),
                        max = round(a$channels$y$limits[2],0)*1.05, 
                        value = req(input$param_numeric), #round(a$par$ylim,0), 
                        step = a$channels$y$limits[2]/20,
                        round = TRUE
            )
            
                        
                        
            
          })
         })
         
         
         #set center of observed spectrum if using custom ion visualization
         observe({
           output$x_target <- renderUI({
             
             if(input$ion_viz3!="custom") {
               numericInput(ns("x_target"), "Center m/z value", value = NULL)
             } else {
               numericInput(ns("x_target"), "Center m/z value", value = mz_viz3())
             }
             
           })
         })
         
        
         
         
         
         observe( {
           
           output$slider_y_max <- renderUI({
             
             req(overview_peaks_sel)
            
             
             a<-Cardinal::plot(overview_peaks_sel, "mean")
             
             
             numericInput(ns("param_numeric"),
                          "Manual y-axis intensity max value",
                          min = round(a$channels$y$limits[1],0),
                          max = round(a$channels$y$limits[2],0),
                          value = round(a$channels$y$limits[2],0)
             )
           })
         })
               
          
         if(!is.null(overview_peaks_sel)) {
           
           updateSliderInput(session, ns("int_range_plot"), value = c(round(Cardinal::plot(req(overview_peaks_sel))$channels$y$limits,0) 
                                                                      ))
           updateNumericInput(session, ns("param_numeric"), value = round(Cardinal::plot(req(overview_peaks_sel))$channels$y$limits[2],0))
          }
        
        
        
         
         observe({
           output$plot4 <- renderImage( {
             
             
             req(overview_peaks_sel)
             req(input$mass_range_plot)
             
             # A temp file to save the output.
             # This file will be removed later by renderImage
             outfile <- tempfile(fileext = '.png')
             
             png(outfile, width = input$width_ms, height = input$height_ms)
             
            
             
             if(length(input$int_range_plot)==1) {
               ylim=c(0, input$int_range_plot)
             } else {
               ylim = input$int_range_plot
             }
             
             #change xlimits based on custom ion or not
             if(is.null(input$x_target) || is.na(input$x_target)){
               xlim=input$mass_range_plot
               
               overview_peaks_sel_plot<-overview_peaks_sel
               
             } else {
               xlim=c(input$x_target-input$x_tol, input$x_target+input$x_tol)
               
               #subsetFeatures to only include mz values within range
               overview_peaks_sel_plot<-subsetFeatures(overview_peaks_sel, mz > xlim[1] & mz < xlim[2])
               
             }
             
             if(!is.finite(xlim[1])){
               xlim=input$mass_range_plot
             }
             
             vp_orig<-vizi_par()
             if(input$spectrum_expand_fonts) {
               req(input$spectrum_axis_size)
               browser()
               vizi_par(
                 cex.axis=req(input$spectrum_axis_size)/100,
                 cex.lab=input$spectrum_label_size/100,
                 cex.main=input$spectrum_label_size/100,
                 #cex.subp=input$spectrum_subtitle_size/100
                 lwd=input$linewidth/100,
                 mar = c(0, 0, 1, 1)
               )
               
               
               #get margins
               #cur_mar<-par()$mar
               
               #new_mar<-c(cur_mar[1]+cex.labp/8, cur_mar[2]+cex.labp/2, cur_mar[3], cur_mar[4])
               
               #cur_mgp<-par()$mgp
               
               #new_mgp<-c(cur_mgp[1]+cex.labp/8, cur_mgp[2], 0)
               #new_mgp<-c(3+max(new_mar[1:2])/20, 1,0)
               
               #lwd=lwdp
                          
               
               
               
               
             } else {
               vizi_par(
                 cex.axis=1,
                 cex.lab=1,
                 cex.main=1,
                 cex.sub=1
               )
             }
             
             #browser()
             #overview_peaks_sel<-Cardinal::summarizeFeatures(overview_peaks_sel)
             
             p1<-Cardinal::plot(overview_peaks_sel_plot,
                                xlim=xlim,
                                ylim =ylim,
                                #cex.axis=req(cex.axisp),
                                #cex.lab=cex.labp,
                                #cex.main=cex.mainp,
                                #cex.sub=cex.subp,
                                #lwd=lwdp,
                                # mar=new_mar, 
                                # mgp=new_mgp, 
                                "mean",
                                annPeaks=input$show_mz,
                                free="y")
             print(p1)
             vizi_par(vp_orig)
             
             #check for ppm calc
             if(input$calc_ppm) {
               
               dat=overview_peaks_sel_plot
               x=mz(dat)
               targ_mz<-req(input$x_target)
               x_sel<-subset(x, x>=xlim[1] & x<= xlim[2])
               
               ppm_error<- round(1e6*(x_sel-targ_mz)/targ_mz, 2)
               
               p1_coord<-p1[[1]][[1]]$marks$peaks$encoding
               
               y_labs<-p1_coord$y[p1_coord$x %in% x_sel]
               
               if(length(ppm_error)==0){
                 showNotification("No ppm error calculated, are there any peaks?", type="warning")
                 return()
               } else {
                print(text(x=x_sel, y=y_labs+ylim[2]*.25, req(ppm_error)))
               }
               
             }
             
             
             
             if(input$show_int) {
              
               
               ###
               p1_coord<-p1[[1]][[1]]$marks$peaks$encoding
               
               
               
               dat=overview_peaks_sel_plot
               x=mz(dat)
               
               x_sel<-subset(x, x>=xlim[1] & x<= xlim[2])
               
               
               
               y_labs<-p1_coord$y[p1_coord$x %in% x_sel]
               
               if(length(y_labs)==0){
                 showNotification("No intensities found, are there any peaks?", type="warning")
                 return()
               } else {
                 print(text(x=x_sel, y=y_labs+ylim[2]*.15, req(round(y_labs, 0))))
               }
               
             }
             
             

             
             dev.off()
             
             # Return a list containing the filename
             list(src = outfile,
                  contentType = 'image/png',
                  width = input$width,
                  height = input$height,
                  alt = "This is alternate text")
           }, deleteFile = TRUE)
           
         })
       }
  })
        
  
}

