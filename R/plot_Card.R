plot_card_UI<-function(id) {
  
  ns <- NS(id)
  
  col_choices<-c(hcl.pals())
  initial_cols<-c("Inferno","Plasma", "Spectral", "Cividis", "Viridis", 
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
               column(6, selectInput(ns("color3"), "Colorscale", col_choices, selected= "Plasma"))
             ),
             fluidRow(
               column(6, selectInput(ns("smooth3"), "Smoothing options", c("none", "mean", "gaussian", "bilateral", "adaptive", "diffusion", "guided"))),
               column(6, checkboxInput(ns("normalize3"), "Scale multiple images?", value = TRUE))
             ),
             
             fluidRow(
               column(6, checkboxInput(ns("colorkey3"), "Draw color key?", value = TRUE)),
               column(6, checkboxInput(ns("dark_bg"), "Dark background?", value = FALSE))
             ),
             
             # NEW: Add pixel masking options
             fluidRow(
               column(12, checkboxInput(ns("apply_mask"), "Apply pData pixel mask?", value = FALSE))
             ),
             uiOutput(ns("mask_options")),
             
             
             fluidRow(
               column(6, numericInput(ns("width_im"), "Image plot width (px)", value = 800, step = 50)),
               column(6, numericInput(ns("height_im"), "Image plot height (px)", value=600, step = 50))
             ),
             downloadButton(ns("download_plot3_pdf"), "Download image (PDF)"),
             checkboxInput(ns("plot_pdata"), "Plot Phenotype data?", value=FALSE),
             uiOutput(ns("plotpdata")),
             checkboxInput(ns("expand_fonts"), "Extended font options?", value=FALSE),
             uiOutput(ns("fonts")),                 
             checkboxInput(ns("expand_runs"), "Select individual runs for plotting only?", value=FALSE),
             uiOutput(ns("select_runs")),
             p("___________________________________________________________________")
             
             
             
      ),
      column(9, uiOutput(ns("plot.window"))
      )  
    ),
    
    uiOutput(ns("spectrum"))
  )
  
}


plot_card_server <- function(id, overview_peaks_sel, spatialOnly=FALSE, allInputs=NULL) {
  
  moduleServer(id, function(input, output, session){
    
    #for dynamic UI
    ns = session$ns
    
    graphics.off()
    plot3_last_png <- reactiveVal(NULL)
    plot3_last_recorded <- reactiveVal(NULL)
    
    cache_plot3_png <- function(src_png) {
      if (is.null(src_png) || !file.exists(src_png)) {
        return(invisible(FALSE))
      }
      prev_png <- isolate(plot3_last_png())
      if (!is.null(prev_png) && file.exists(prev_png)) {
        unlink(prev_png)
      }
      cached_png <- tempfile(pattern = paste0("plot3_cache_", session$token, "_"), fileext = ".png")
      copied <- isTRUE(file.copy(src_png, cached_png, overwrite = TRUE))
      if (copied) {
        plot3_last_png(cached_png)
      }
      invisible(copied)
    }

    capture_plot3_recording <- function() {
      rec <- try(grDevices::recordPlot(), silent = TRUE)
      if (!inherits(rec, "try-error")) {
        plot3_last_recorded(rec)
      }
      invisible(NULL)
    }
    
    make_plot3_image_result <- function(src_png, alt_text) {
      capture_plot3_recording()
      cache_plot3_png(src_png)
      list(
        src = src_png,
        contentType = 'image/png',
        width = input$width_im,
        height = input$height_im,
        alt = alt_text
      )
    }
    
    output$download_plot3_pdf <- downloadHandler(
      filename = function() {
        paste0("plot_card_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
      },
      content = function(file) {
        width_px <- suppressWarnings(as.numeric(input$width_im))
        height_px <- suppressWarnings(as.numeric(input$height_im))
        if (!is.finite(width_px) || width_px <= 0) width_px <- 800
        if (!is.finite(height_px) || height_px <= 0) height_px <- 600
        
        grDevices::pdf(file = file, width = width_px / 72, height = height_px / 72, onefile = TRUE)
        on.exit(grDevices::dev.off(), add = TRUE)
        
        rec <- plot3_last_recorded()
        if (!is.null(rec)) {
          try(grDevices::replayPlot(rec), silent = TRUE)
          return(invisible(NULL))
        }
        
        src_png <- plot3_last_png()
        if (is.null(src_png) || !file.exists(src_png)) {
          showNotification("No image available yet. Render a plot first, then download PDF.", type = "warning", duration = 5)
          stop("No rendered image available for PDF export.")
        }
        
        img <- try(png::readPNG(src_png), silent = TRUE)
        if (inherits(img, "try-error")) {
          stop("Could not read cached image for PDF export.")
        }
        
        grid::grid.newpage()
        grid::grid.raster(img, x = 0.5, y = 0.5, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"), interpolate = TRUE)
      }
    )
    
    #create new overview_peaks_sel object with mean values
    if(is.null(fData(overview_peaks_sel)$mean)) {
      overview_peaks_sel <- try(Cardinal::summarizeFeatures(overview_peaks_sel, verbose = FALSE), silent = TRUE)
      if (inherits(overview_peaks_sel, "try-error")) {
        showNotification("No data available, please check your parameters or dataset", type="error")
        message("plot_card_server summarizeFeatures failed: ", as.character(overview_peaks_sel))
        return()
      }
    }
    
    sanitize_limits <- function(x, default = c(0, 1), clamp_nonnegative = FALSE) {
      x <- suppressWarnings(as.numeric(x))
      x <- x[is.finite(x)]
      if (length(x) < 2) {
        x <- default
      } else {
        x <- range(x)
      }
      if (clamp_nonnegative) {
        x[1] <- max(0, x[1])
      }
      if (!is.finite(x[1]) || !is.finite(x[2]) || x[2] <= x[1]) {
        x <- default
      }
      x
    }
    
    get_safe_plot_limits <- function(obj) {
      mz_vals <- try(suppressWarnings(as.numeric(Cardinal::mz(obj))), silent = TRUE)
      mz_limits <- if (inherits(mz_vals, "try-error")) c(0, 1) else sanitize_limits(mz_vals, default = c(0, 1))
      
      mean_vals <- try(suppressWarnings(as.numeric(fData(obj)$mean)), silent = TRUE)
      if (inherits(mean_vals, "try-error")) {
        y_limits <- c(0, 1)
      } else {
        y_limits <- sanitize_limits(mean_vals, default = c(0, 1), clamp_nonnegative = TRUE)
      }
      
      list(x = mz_limits, y = y_limits)
    }

    normalize_pdata_labels <- function(x) {
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
        if(!is.null(tbl$freq)) {
          tbl$freq<-round(tbl$freq, 2)
        }
        if(!is.null(tbl$mean)) {
          tbl$mean<-round(tbl$mean, 1)
        }
        tbl$mz<-round(tbl$mz, 4)
        
        
        updateNumericInput(session, ("width_im"), value=800, step=50)
        updateNumericInput(session, ("height_im"), value=450, step=50)
        
        
        
        #add something here at some point to only select mz, freq, count, and mean columns
        
        
        output$tbl <-DT::renderDataTable({
          tbl
          
        }, selection = "multiple")
        
        list(
          
          checkboxInput(ns("superpose"), "Superpose images?", value=FALSE),
          selectInput(ns("display_mode"), "Ion math?", c("none", "RGB (2-3 ions)" = "rgb", "sum", "ratio", "subtract", "min", "max", "mean", "sd", "var", "multiply")),
          numericInput(ns("plusminus_viz3"), "+/- m/z for visualization",0.05)
          
        )
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
            column(
              3,
              numericInput(
                ns("spectrum_label_size"),
                label = "Label font scaling (%)",
                min = 0,
                value = 100
              )
            ),
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
        
        # MODIFIED: Apply pixel mask to the data before plotting
        overview_peaks_sel_masked <- apply_pixel_mask(overview_peaks_sel)
        
        if(!is.null(input$select_runs)) {
          overview_peaks_sel_masked <- subsetPixels(overview_peaks_sel_masked, run %in% input$select_runs)
        }

        if (is.null(overview_peaks_sel_masked) || ncol(overview_peaks_sel_masked) == 0) {
          graphics::plot.new()
          graphics::text(0.5, 0.5, "No pixels available after mask/run filtering.")
          dev.off()
          return(make_plot3_image_result(outfile, "No pixels available"))
        }
        
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
          pdata_df <- as.data.frame(pData(overview_peaks_sel_masked))
          if (!input$pdata_var_plot %in% colnames(pdata_df)) {
            showNotification(
              paste0("pData field '", input$pdata_var_plot, "' is not available in the current selection."),
              type = "warning",
              duration = 5
            )
            graphics::plot.new()
            graphics::text(0.5, 0.5, "Selected pData field is not available.")
            dev.off()
            return(make_plot3_image_result(outfile, "Invalid pData field"))
          }

          pdata_vals <- pdata_df[[input$pdata_var_plot]]
          if (length(pdata_vals) != ncol(overview_peaks_sel_masked)) {
            showNotification(
              "pData length mismatch for selected field; skipping phenotype plot.",
              type = "error",
              duration = 6
            )
            graphics::plot.new()
            graphics::text(0.5, 0.5, "pData length mismatch.")
            dev.off()
            return(make_plot3_image_result(outfile, "pData mismatch"))
          }

          is_numeric_pdata <- is.numeric(pdata_vals) || is.integer(pdata_vals)
          if (is_numeric_pdata) {
            #create list of arguments for image
            arg_list <- list(
              overview_peaks_sel_masked,
              input$pdata_var_plot,
              key = (input$colorkey3),
              col = cpal(input$color3),
              grid = FALSE
            )

            smoothing_option <- if (input$smooth3 != "none") input$smooth3 else NULL
            if (!is.null(smoothing_option)) {
              arg_list$smooth <- smoothing_option
            }

            if (input$dark_bg) {
              arg_list$style <- "dark"
            }

            plt_tmp <- do.call(Cardinal::image, arg_list)
            print(plt_tmp)
          } else {
            # Categorical pData plotted via ggplot to avoid pData mutation issues.
            labels <- normalize_pdata_labels(pdata_vals)
            if (length(labels) == 0) {
              labels <- rep("NA", ncol(overview_peaks_sel_masked))
            }
            levels_lab <- unique(labels)
            valid_col <- are_valid_colors(levels_lab)
            col_map <- stats::setNames(rep("grey70", length(levels_lab)), levels_lab)

            if (any(valid_col)) {
              col_map[valid_col] <- levels_lab[valid_col]
            }
            if (any(!valid_col)) {
              fallback_cols <- grDevices::hcl.colors(sum(!valid_col), palette = "Dark 3")
              col_map[!valid_col] <- fallback_cols
            }

            coords_df <- as.data.frame(Cardinal::coord(overview_peaks_sel_masked))
            xv <- if ("x" %in% colnames(coords_df)) coords_df$x else coords_df[[1]]
            yv <- if ("y" %in% colnames(coords_df)) coords_df$y else if (ncol(coords_df) >= 2) coords_df[[2]] else rep(1, length(xv))
            run_vec <- as.character(Cardinal::run(overview_peaks_sel_masked))
            if (length(run_vec) != length(labels)) {
              run_vec <- rep("run", length(labels))
            }

            plot_df <- data.frame(
              x = as.numeric(xv),
              y = as.numeric(yv),
              run = run_vec,
              label = factor(labels, levels = levels_lab),
              stringsAsFactors = FALSE
            )

            p_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y, fill = label)) +
              ggplot2::geom_raster() +
              ggplot2::facet_wrap(. ~ run) +
              ggplot2::coord_fixed() +
              ggplot2::scale_y_continuous(trans = "reverse") +
              ggplot2::theme_minimal() +
              ggplot2::scale_fill_manual(
                values = unname(col_map[levels_lab]),
                name = input$pdata_var_plot,
                drop = FALSE
              )

            if (!isTRUE(input$colorkey3)) {
              p_plot <- p_plot + ggplot2::guides(fill = "none")
            }

            if (isTRUE(input$dark_bg)) {
              p_plot <- p_plot +
                ggplot2::theme(
                  panel.background = ggplot2::element_rect(fill = "black", colour = NA),
                  plot.background = ggplot2::element_rect(fill = "black", colour = NA),
                  strip.background = ggplot2::element_rect(fill = "black", colour = NA),
                  strip.text = ggplot2::element_text(color = "white"),
                  text = ggplot2::element_text(color = "white"),
                  axis.text = ggplot2::element_text(color = "white"),
                  axis.title = ggplot2::element_text(color = "white"),
                  legend.background = ggplot2::element_rect(fill = "black", colour = NA),
                  legend.key = ggplot2::element_rect(fill = "black", colour = NA),
                  legend.text = ggplot2::element_text(color = "white"),
                  legend.title = ggplot2::element_text(color = "white")
                )
            }

            print(p_plot)

            if (input$smooth3 != "none") {
              showNotification(
                paste0("Smoothing is only applied to numeric pData. Skipping for '", input$pdata_var_plot, "'."),
                type = "message",
                duration = 3
              )
            }
          }
          
          vizi_par(vp_orig)
        } else if (input$ion_viz3=="viz_all") {
          
          
          mz_range=range(mz(overview_peaks_sel_masked))
          
          #find closest mz value to middle of range
          test_value <- mean(mz_range)
          
          
          # Calculate the absolute differences
          differences <- abs(mz(overview_peaks_sel_masked) - test_value)
          
          # Find the index of the minimum difference
          closest_index <- which.min(differences)
          
          mz_set=mz(overview_peaks_sel_masked)[closest_index]
          tol=max(differences) + differences[closest_index]+1
          
          plusminus=tol
          
          smoothing_option <- if (input$smooth3 != "none")  input$smooth3 else NULL
          enhance_option <- if (input$contrast3 != "none")  input$contrast3 else NULL
          
          
          arg_list<-list(overview_peaks_sel_masked, 
                         mz=mz_set,
                         tolerance=round(plusminus,3), 
                         units='mz',
                         col=cpal(input$color3),
                         enhance=enhance_option,
                         smooth=smoothing_option,
                         scale=input$normalize3,
                         key=(input$colorkey3),
                         grid=FALSE)
          
          if(input$dark_bg) {
            arg_list$style <- "dark"
          }
          
          print(do.call(Cardinal::image, arg_list))
          
          vizi_par(vp_orig)
          
          
        } else if (input$ion_viz3=="custom"){
          
          if(is.null(mz_viz3())){
            ion=mz(overview_peaks_sel_masked[1,])
          } else {
            ion=as.numeric(mz_viz3())
          }

          mz_axis_all <- suppressWarnings(as.numeric(Cardinal::mz(overview_peaks_sel_masked)))
          n_pix_masked <- ncol(overview_peaks_sel_masked)

          coerce_spectra_matrix <- function(spec_obj, n_feat_expected, n_pix_expected) {
            mat <- try(as.matrix(spec_obj), silent = TRUE)
            if (inherits(mat, "try-error")) return(NULL)
            if (!is.matrix(mat)) {
              vec <- suppressWarnings(as.numeric(mat))
              if (length(vec) != (n_feat_expected * n_pix_expected)) return(NULL)
              mat <- matrix(vec, nrow = n_feat_expected, ncol = n_pix_expected, byrow = TRUE)
            }
            if (nrow(mat) == n_pix_expected && ncol(mat) == n_feat_expected) {
              mat <- t(mat)
            }
            if (nrow(mat) != n_feat_expected || ncol(mat) != n_pix_expected) return(NULL)
            mat
          }

          resolve_feature_index_exact <- function(target_mz) {
            if (!is.finite(target_mz) || length(mz_axis_all) == 0) return(NA_integer_)
            idx <- tryCatch(Cardinal::features(overview_peaks_sel_masked, mz = target_mz), error = function(e) integer(0))
            idx <- suppressWarnings(as.integer(idx))
            idx <- idx[is.finite(idx) & idx >= 1 & idx <= length(mz_axis_all)]
            if (length(idx) > 0) return(idx[[1]])
            diffs <- abs(mz_axis_all - target_mz)
            diffs[!is.finite(diffs)] <- Inf
            if (!any(is.finite(diffs))) return(NA_integer_)
            which.min(diffs)
          }

          extract_exact_ions_matrix_once <- function(target_ions) {
            target_ions <- suppressWarnings(as.numeric(target_ions))
            if (length(target_ions) == 0 || any(!is.finite(target_ions))) return(NULL)
            idx_map <- vapply(target_ions, resolve_feature_index_exact, integer(1))
            if (any(!is.finite(idx_map)) || any(idx_map < 1)) return(NULL)
            uniq_idx <- unique(idx_map)
            spec_obj <- try(Cardinal::spectra(overview_peaks_sel_masked[uniq_idx, ]), silent = TRUE)
            if (inherits(spec_obj, "try-error")) return(NULL)
            mat_all <- coerce_spectra_matrix(spec_obj, length(uniq_idx), n_pix_masked)
            if (is.null(mat_all)) return(NULL)
            row_map <- match(idx_map, uniq_idx)
            list(
              mat = mat_all[row_map, , drop = FALSE],
              feature_idx = idx_map,
              ions = target_ions
            )
          }

          extract_rgb_channels_once <- function(target_ions, tol) {
            target_ions <- suppressWarnings(as.numeric(target_ions))
            if (length(target_ions) == 0 || any(!is.finite(target_ions))) return(NULL)
            sel_list <- lapply(target_ions, function(target_mz) {
              which(is.finite(mz_axis_all) & abs(mz_axis_all - target_mz) <= tol)
            })
            if (any(lengths(sel_list) == 0)) return(NULL)
            uniq_idx <- sort(unique(unlist(sel_list, use.names = FALSE)))
            if (length(uniq_idx) == 0) return(NULL)
            spec_obj <- try(Cardinal::spectra(overview_peaks_sel_masked[uniq_idx, ]), silent = TRUE)
            if (inherits(spec_obj, "try-error")) return(NULL)
            mat_all <- coerce_spectra_matrix(spec_obj, length(uniq_idx), n_pix_masked)
            if (is.null(mat_all)) return(NULL)
            channels <- lapply(sel_list, function(idx) {
              rows <- match(idx, uniq_idx)
              rows <- rows[is.finite(rows)]
              if (length(rows) == 0) return(NULL)
              vals <- mat_all[rows, , drop = FALSE]
              if (nrow(vals) == 1) {
                as.numeric(vals[1, ])
              } else {
                colMeans(vals, na.rm = TRUE)
              }
            })
            if (any(vapply(channels, is.null, logical(1)))) return(NULL)
            channels
          }
          
          
          if(!is.null(input$display_mode) && input$display_mode!="none"){
            if (identical(input$display_mode, "rgb")) {
              if (!(length(ion) %in% c(2, 3))) {
                showNotification("RGB mode requires 2 or 3 selected ions.", type = "error")
                message("RGB mode requires 2 or 3 selected ions.")
                return()
              }
              
              coords_df <- as.data.frame(Cardinal::coord(overview_peaks_sel_masked))
              validate(need(all(c("x", "y") %in% colnames(coords_df)), "MSI coordinates must contain x and y."))
              run_vec <- as.character(Cardinal::run(overview_peaks_sel_masked))
              
              tol_rgb <- suppressWarnings(as.numeric(input$plusminus_viz3))
              if (!is.finite(tol_rgb) || tol_rgb <= 0) {
                tol_rgb <- 0.05
              }
              
              smooth_mode <- tolower(trimws(as.character(input$smooth3)[1]))
              if (!smooth_mode %in% c("none", "mean", "gaussian", "bilateral", "adaptive", "diffusion", "guided")) {
                smooth_mode <- "none"
              }
              
              contrast_mode <- tolower(trimws(as.character(input$contrast3)[1]))
              if (!contrast_mode %in% c("none", "histogram", "adaptive")) {
                contrast_mode <- "none"
              }
              
              scale01 <- function(x, shared_range = NULL) {
                x <- as.numeric(x)
                x[!is.finite(x)] <- NA_real_
                if (is.null(shared_range)) {
                  ok <- is.finite(x)
                  if (!any(ok)) return(rep(0, length(x)))
                  lo <- min(x[ok], na.rm = TRUE)
                  hi <- max(x[ok], na.rm = TRUE)
                } else {
                  lo <- shared_range[1]
                  hi <- shared_range[2]
                }
                if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
                  return(rep(0, length(x)))
                }
                out <- (x - lo) / (hi - lo)
                out[!is.finite(out)] <- 0
                pmin(1, pmax(0, out))
              }
              
              to_01 <- function(x, mode = "none", shared_vals = NULL) {
                x <- as.numeric(x)
                x[!is.finite(x)] <- NA_real_
                ref <- if (is.null(shared_vals)) x else as.numeric(shared_vals)
                ref <- ref[is.finite(ref)]
                if (length(ref) == 0) {
                  return(rep(0, length(x)))
                }
                
                if (identical(mode, "adaptive")) {
                  lo_ref <- min(ref, na.rm = TRUE)
                  hi_ref <- max(ref, na.rm = TRUE)
                  if (!is.finite(lo_ref) || !is.finite(hi_ref) || hi_ref <= lo_ref) {
                    return(rep(0, length(x)))
                  }
                  ec <- stats::ecdf(ref)
                  base <- ec(lo_ref)
                  den <- 1 - base
                  if (!is.finite(den) || den <= 0) {
                    return(rep(0, length(x)))
                  }
                  out <- (ec(x) - base) / den
                  out[!is.finite(out)] <- 0
                  return(pmin(1, pmax(0, out)))
                }
                
                if (identical(mode, "histogram")) {
                  q <- stats::quantile(ref, probs = c(0.01, 0.99), na.rm = TRUE, names = FALSE, type = 8)
                  lo <- q[1]
                  hi <- q[2]
                  x2 <- pmin(pmax(x, lo), hi)
                } else {
                  lo <- min(ref, na.rm = TRUE)
                  hi <- max(ref, na.rm = TRUE)
                  x2 <- x
                }
                
                if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
                  return(rep(0, length(x)))
                }
                
                out <- (x2 - lo) / (hi - lo)
                out[!is.finite(out)] <- 0
                pmin(1, pmax(0, out))
              }
              
              smooth_vec <- function(v, k) {
                w <- as.numeric(!is.na(v))
                v0 <- v
                v0[is.na(v0)] <- 0
                num <- as.numeric(stats::filter(v0, k, sides = 2))
                den <- as.numeric(stats::filter(w, k, sides = 2))
                out <- num / den
                out[!is.finite(out) | den <= 0] <- NA_real_
                out
              }
              
              smooth_matrix <- function(mat, method) {
                if (identical(method, "none")) {
                  return(mat)
                }
                
                if (method %in% c("mean", "adaptive", "guided")) {
                  k <- rep(1 / 3, 3)
                } else {
                  xk <- -3:3
                  sigma <- 1
                  k <- exp(-(xk * xk) / (2 * sigma * sigma))
                  k <- k / sum(k)
                }
                
                tmp <- t(apply(mat, 1, smooth_vec, k = k))
                out <- apply(tmp, 2, smooth_vec, k = k)
                if (!is.matrix(out)) {
                  out <- matrix(out, nrow = nrow(mat), ncol = ncol(mat))
                }
                out
              }
              
              smooth_channel_by_run <- function(vals, coords, runs, method) {
                if (identical(method, "none")) {
                  return(vals)
                }
                
                out <- vals
                run_split <- split(seq_along(runs), ifelse(is.na(runs), "__NA_RUN__", runs))
                
                for (idx in run_split) {
                  crd <- coords[idx, , drop = FALSE]
                  x_norm <- as.integer(crd$x - min(crd$x, na.rm = TRUE) + 1L)
                  y_norm <- as.integer(crd$y - min(crd$y, na.rm = TRUE) + 1L)
                  nx <- max(x_norm, na.rm = TRUE)
                  ny <- max(y_norm, na.rm = TRUE)
                  row_idx <- ny - y_norm + 1L
                  
                  mat <- matrix(NA_real_, nrow = ny, ncol = nx)
                  mat[cbind(row_idx, x_norm)] <- vals[idx]
                  
                  mat_s <- smooth_matrix(mat, method)
                  out[idx] <- as.numeric(mat_s[cbind(row_idx, x_norm)])
                }
                out
              }
              
              rgb_channels <- extract_rgb_channels_once(ion[seq_len(min(length(ion), 3))], tol_rgb)
              ch_r <- if (!is.null(rgb_channels) && length(rgb_channels) >= 1) rgb_channels[[1]] else NULL
              ch_g <- if (!is.null(rgb_channels) && length(rgb_channels) >= 2) rgb_channels[[2]] else NULL
              ch_b <- if (length(ion) == 3) {
                if (!is.null(rgb_channels) && length(rgb_channels) >= 3) rgb_channels[[3]] else NULL
              } else {
                rep(0, ncol(overview_peaks_sel_masked))
              }
              
              if (is.null(ch_r) || is.null(ch_g) || is.null(ch_b)) {
                showNotification("Could not resolve one or more RGB ions within tolerance. Increase +/- m/z or adjust ion choices.", type = "error")
                return()
              }
              
              if (!identical(smooth_mode, "none")) {
                ch_r <- smooth_channel_by_run(ch_r, coords_df, run_vec, smooth_mode)
                ch_g <- smooth_channel_by_run(ch_g, coords_df, run_vec, smooth_mode)
                ch_b <- smooth_channel_by_run(ch_b, coords_df, run_vec, smooth_mode)
              }
              
              if (isTRUE(input$normalize3)) {
                r01 <- to_01(ch_r, mode = contrast_mode)
                g01 <- to_01(ch_g, mode = contrast_mode)
                b01 <- to_01(ch_b, mode = contrast_mode)
              } else {
                pooled <- c(ch_r, ch_g, ch_b)
                pooled <- pooled[is.finite(pooled)]
                if (length(pooled) == 0) {
                  r01 <- g01 <- b01 <- rep(0, length(ch_r))
                } else {
                  if (identical(contrast_mode, "adaptive")) {
                    r01 <- to_01(ch_r, mode = contrast_mode, shared_vals = pooled)
                    g01 <- to_01(ch_g, mode = contrast_mode, shared_vals = pooled)
                    b01 <- to_01(ch_b, mode = contrast_mode, shared_vals = pooled)
                  } else {
                    shared_rng <- range(pooled, na.rm = TRUE)
                    if (identical(contrast_mode, "histogram")) {
                      shared_rng <- stats::quantile(pooled, probs = c(0.01, 0.99), na.rm = TRUE, names = FALSE, type = 8)
                    }
                    r01 <- scale01(if (identical(contrast_mode, "histogram")) pmin(pmax(ch_r, shared_rng[1]), shared_rng[2]) else ch_r, shared_rng)
                    g01 <- scale01(if (identical(contrast_mode, "histogram")) pmin(pmax(ch_g, shared_rng[1]), shared_rng[2]) else ch_g, shared_rng)
                    b01 <- scale01(if (identical(contrast_mode, "histogram")) pmin(pmax(ch_b, shared_rng[1]), shared_rng[2]) else ch_b, shared_rng)
                  }
                }
              }
              
              # In 2-ion RGB mode, keep blue strictly zero after any smoothing/contrast transform.
              if (length(ion) == 2) {
                b01 <- rep(0, length(r01))
              }
              
              rgb_df <- data.frame(
                x = coords_df$x,
                y = coords_df$y,
                runs = run_vec,
                r = r01,
                g = g01,
                b = b01,
                stringsAsFactors = FALSE
              )
              rgb_df$fill_col <- grDevices::rgb(rgb_df$r, rgb_df$g, rgb_df$b)
              
              title_txt <- if (length(ion) == 3) {
                paste0("RGB ion image (R=", round(ion[1], 4), ", G=", round(ion[2], 4), ", B=", round(ion[3], 4), ")")
              } else {
                paste0("RGB ion image (R=", round(ion[1], 4), ", G=", round(ion[2], 4), ", B=0)")
              }
              
              p_rgb <- ggplot2::ggplot(rgb_df, ggplot2::aes(x = x, y = y, fill = fill_col)) +
                ggplot2::geom_raster() +
                ggplot2::facet_wrap(. ~ runs) +
                ggplot2::scale_fill_identity() +
                ggplot2::coord_equal() +
                ggplot2::scale_y_continuous(trans = "reverse") +
                ggplot2::theme_minimal() +
                ggplot2::labs(
                  title = title_txt,
                  subtitle = paste0("+/- ", round(tol_rgb, 4), " m/z | contrast=", contrast_mode, " | smooth=", smooth_mode)
                )
              
              if (isTRUE(input$dark_bg)) {
                p_rgb <- p_rgb +
                  ggplot2::theme(
                    panel.background = ggplot2::element_rect(fill = "black", colour = NA),
                    plot.background = ggplot2::element_rect(fill = "black", colour = NA),
                    panel.grid = ggplot2::element_blank(),
                    strip.background = ggplot2::element_rect(fill = "grey20", colour = NA),
                    strip.text = ggplot2::element_text(color = "white"),
                    axis.text = ggplot2::element_text(color = "white"),
                    axis.title = ggplot2::element_text(color = "white"),
                    plot.title = ggplot2::element_text(color = "white"),
                    plot.subtitle = ggplot2::element_text(color = "grey80")
                  )
              }
              
              print(p_rgb)
              vizi_par(vp_orig)
            } else {
              if(input$display_mode%in%c("min", "max", "min", "mean", "sum", "sd", "var")) { 
              
              
              
              select_vec<-as.character(mz(overview_peaks_sel_masked)) %in% as.character(ion)
              #test to make sure there are 2 or more elements
              if(sum(select_vec)<2){
                showNotification("At least two ions required for this calculation", type="error")
                message("At least two ions required for this calculation")
                return()
              }
              
              sm<-summarizePixels(overview_peaks_sel_masked[select_vec,], stat=c(xic=input$display_mode), as="DataFrame")
              pData(overview_peaks_sel_masked)$xic<-sm$xic
              
              label_txt=paste(input$display_mode, "mz(s)=", paste(ion, collapse=", "))
              
            }else if(input$display_mode=="ratio"){
              if(length(ion)!=2){
                showNotification("Exactly two ions required for ratio (mz1/mz2)", type="error")
                message("Exactly two ions required for ratio (mz1/mz2)")
                return()
              } else {
                exact_ion_vals <- extract_exact_ions_matrix_once(ion[1:2])
                if (is.null(exact_ion_vals) || nrow(exact_ion_vals$mat) < 2) {
                  showNotification("Unable to extract one or more ions for ratio calculation.", type="error")
                  message("Unable to extract one or more ions for ratio calculation.")
                  return()
                }
                mz1 <- as.numeric(exact_ion_vals$mat[1, ])
                mz2 <- as.numeric(exact_ion_vals$mat[2, ])
                overview_peaks_sel_masked$xic <- (1 + mz1) / (1 + mz2)
                
                ion=round(ion, 4)
                label_txt=paste("ratio of",ion[1],"/",ion[2])
              }
              
            } else if(input$display_mode=="subtract") {
              if(length(ion)!=2){
                showNotification("Exactly two ions required for subtraction (mz1-mz2)", type="error")
                message("Exactly two ions required for subtraction (mz1-mz2)")
                return()
              } else {
                exact_ion_vals <- extract_exact_ions_matrix_once(ion[1:2])
                if (is.null(exact_ion_vals) || nrow(exact_ion_vals$mat) < 2) {
                  showNotification("Unable to extract one or more ions for subtraction calculation.", type="error")
                  message("Unable to extract one or more ions for subtraction calculation.")
                  return()
                }
                mz1 <- as.numeric(exact_ion_vals$mat[1, ])
                mz2 <- as.numeric(exact_ion_vals$mat[2, ])
                overview_peaks_sel_masked$xic <- (mz1) - (mz2)
                ion=round(ion, 4)
                label_txt=paste("difference of",ion[1],"-",ion[2])
                
              }
              
              
            }else if(input$display_mode=="multiply"){
              nelements=length(ion)
              exact_ion_vals <- extract_exact_ions_matrix_once(ion)
              if (is.null(exact_ion_vals) || nrow(exact_ion_vals$mat) < 1) {
                showNotification("Unable to extract ions for multiply calculation.", type="error")
                message("Unable to extract ions for multiply calculation.")
                return()
              }
              mat_mul <- 1 + exact_ion_vals$mat
              overview_peaks_sel_masked$xic <- apply(mat_mul, 2, prod)
              ion=round(ion, 4)
              label_txt=paste(ion, collapse=" * ")
            }
            
            
            
            if(sum(is.na(pData(overview_peaks_sel_masked)$xic))==length(overview_peaks_sel_masked)) {
              showNotification("This calculation does not work!")
              return()
            }
            
            plusminus=input$plusminus_viz3
            
            smoothing_option <- if (input$smooth3 != "none")  input$smooth3 else NULL
            enhance_option <- if (input$contrast3 != "none")  input$contrast3 else NULL
            
            
            arg_list<-list(overview_peaks_sel_masked,
                           'xic',
                           tolerance=round(plusminus,3), 
                           units='mz',
                           col=cpal(input$color3),
                           enhance=enhance_option,
                           smooth=smoothing_option,
                           scale=input$normalize3,
                           key=(input$colorkey3),
                           grid=FALSE)
            
            if(input$dark_bg) {
              arg_list$style <- "dark"
            }
            
            print(matter::as_facets(do.call(Cardinal::image, arg_list), labels=label_txt))
            
            vizi_par(vp_orig)
            }
            
          } else {
            
            mz_set=ion
            
            mz_range=range(mz(overview_peaks_sel_masked))
            
            plusminus=input$plusminus_viz3
            
            smoothing_option <- if (input$smooth3 != "none")  input$smooth3 else NULL
            enhance_option <- if (input$contrast3 != "none")  input$contrast3 else NULL
            
            
            arg_list<-list(overview_peaks_sel_masked,
                           mz=mz_set,
                           tolerance=round(plusminus,3), 
                           units='mz',
                           col=cpal(input$color3),
                           enhance=enhance_option,
                           smooth=smoothing_option,
                           scale=input$normalize3,
                           superpose=input$superpose,
                           key=(input$colorkey3),
                           grid=FALSE)
            
            if(input$dark_bg) {
              arg_list$style <- "dark"
            }
            
            print(do.call(Cardinal::image, arg_list))
            
            vizi_par(vp_orig)
            
            
          }
        } else if (input$ion_viz3=="viz_first") {
          
          tol=0.05
          image_command <-Cardinal::image(overview_peaks_sel_masked, 
                                          col=cpal(input$color3),
                                          scale=input$normalize3,
                                          key=(input$colorkey3),
                                          grid=FALSE
          )
          
          print(image_command)
          
          vizi_par(vp_orig)
          
        }
        
        
        #browser()
        if (!is.null(allInputs) && !is.null(allInputs$histology_upload) && !is.null(allInputs$msi_upload)) {
          
          
          ensure_3d <- function(image) {
            # Check if image is grayscale (2 dimensions)
            if (length(dim(image)) == 2) {
              # Convert grayscale to RGB by replicating the matrix across 3 layers
              image <- abind::abind(image, image, image, along = 3)
            }
            # Check if image has more than 3 dimensions and reduce if necessary
            else if (length(dim(image)) > 3) {
              image <- image[,,,1]  # Take the first layer if multiple are present
            }
            return(image)
          }
          
          #browser()
          # Load the uploaded histology image
          histology_path <- allInputs$histology_upload$datapath
          file_type <- tools::file_ext(histology_path)
          histology_image <- switch(file_type,
                                    "png" = png::readPNG(histology_path),
                                    "jpg" = jpeg::readJPEG(histology_path),
                                    "jpeg" = jpeg::readJPEG(histology_path),
                                    stop("Unsupported file type"))
          
          histology_image <- ensure_3d(histology_image)
          
          # Assuming histology_image has dimensions in pixels
          image_width <- dim(histology_image)[2]
          image_height <- dim(histology_image)[1]
          
          if(allInputs$debug==T){
            showNotification("Debug mode enabled for overlay.", type = "message")
          }
          
          # Create a grob (graphics object) from the histology image with the specified alpha
          
          # For PNG images with RGBA channels, keep the alpha as is.
          if (file_type == "png") {
            histology_image[,,4] <- histology_image[,,4] * allInputs$alpha
          }
          
          # For JPEG images, create an RGBA image where the alpha layer is set to the transparency value
          if (file_type %in% c("jpg", "jpeg")) {
            # Create an empty alpha channel with the specified transparency
            alpha_channel <- array(allInputs$alpha, dim = c(dim(histology_image)[1], dim(histology_image)[2]))
            # Combine the RGB layers with the new alpha channel
            histology_image <- abind::abind(histology_image, alpha_channel, along = 3)
          }
          
          # Now create the rasterGrob with the RGBA image data
          histology_grob <- grid::rasterGrob(histology_image, interpolate = TRUE)
          
          histology_grob <- grid::editGrob(
            histology_grob,
            vp = grid::viewport(
              x = grid::unit(0.5, "npc") + grid::unit(allInputs$translate_x, "mm"),
              y = grid::unit(0.5, "npc") + grid::unit(allInputs$translate_y, "mm"),
              angle = allInputs$rotate,
              width = grid::unit(input$width_im * allInputs$scalex, "points"),
              height = grid::unit(input$height_im * allInputs$scaley, "points"),
              just = c("center", "center")
            )
          )
          
          (grid::pushViewport(grid::viewport(x = 0.5, y = 0.5, width = 1, height = 1)))
          (grid::grid.draw(histology_grob))
          (grid::popViewport())
          
        }
        
        dev.off()
        
        # Return a list containing the filename
        make_plot3_image_result(outfile, "This is alternate text")
      }, deleteFile = TRUE)
    }) 
    
    
    if(spatialOnly==FALSE) {
      observe({  
        output$plot_ranges<- renderUI( {
          
          req(overview_peaks_sel)
          
          lims <- get_safe_plot_limits(overview_peaks_sel)
          
          
          sliderInput(ns("mass_range_plot"), 
                      label = p("m/z range for MS plot (X)"), 
                      min = round(lims$x[1], 4),
                      max = round(lims$x[2], 4), 
                      value = round(lims$x, 4), 
                      step = NULL,
                      round = TRUE)
          
          
        })
      })
      
      observe({
        output$plot_ranges2<- renderUI( {
          req(overview_peaks_sel)
          lims <- get_safe_plot_limits(overview_peaks_sel)
          y_step <- max((lims$y[2] - lims$y[1]) / 20, .Machine$double.eps)
          y_top <- if (is.null(input$param_numeric) || !is.finite(input$param_numeric)) lims$y[2] else input$param_numeric
          y_top <- min(max(y_top, lims$y[1]), lims$y[2] * 1.05)
          
          sliderInput (ns("int_range_plot"), 
                       label = p("intensity range for MS plot (Y)"), 
                       min = round(lims$y[1], 0),
                       max = round(lims$y[2] * 1.05, 0), 
                       value = c(round(lims$y[1], 0), round(y_top, 0)),
                       step = y_step,
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
          
          
          lims <- get_safe_plot_limits(overview_peaks_sel)
          
          
          numericInput(ns("param_numeric"),
                       "Manual y-axis intensity max value",
                       min = round(lims$y[1], 0),
                       max = round(lims$y[2], 0),
                       value = round(lims$y[2], 0)
          )
        })
      })
      
      
      if(!is.null(overview_peaks_sel)) {
        lims <- get_safe_plot_limits(overview_peaks_sel)
        
        updateSliderInput(session, ns("int_range_plot"), value = c(round(lims$y, 0)))
        updateNumericInput(session, ns("param_numeric"), value = round(lims$y[2], 0))
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
            vizi_par(
              cex.axis=req(input$spectrum_axis_size)/100,
              cex.lab=input$spectrum_label_size/100,
              cex.main=input$spectrum_label_size/100,
              lwd=input$linewidth/100,
              mar = c(0, 0, 1, 1)
            )
            
          } else {
            vizi_par(
              cex.axis=1,
              cex.lab=1,
              cex.main=1,
              cex.sub=1
            )
          }
          
          # MODIFIED: Apply pixel mask to spectrum plot data as well
          overview_peaks_sel_plot_masked <- apply_pixel_mask(overview_peaks_sel_plot)
          
          p1 <- try(
            Cardinal::plot(overview_peaks_sel_plot_masked,
                           xlim=xlim,
                           ylim =ylim,
                           grid=FALSE,
                           "mean",
                           annPeaks=input$show_mz,
                           free="y"),
            silent = TRUE
          )
          if (inherits(p1, "try-error")) {
            message("plot_card_server spectrum plotting failed: ", as.character(p1))
            showNotification(
              "Spectrum plotting failed for this file. Ion image plotting is still available.",
              type = "warning",
              duration = 8
            )
            plot.new()
            text(0.5, 0.6, "Spectrum plot unavailable for this dataset.", cex = 1.1)
            text(0.5, 0.45, "Try a source .imzML export if this came from a converted file.", cex = 0.9)
            vizi_par(vp_orig)
            dev.off()
            return(list(src = outfile,
                        contentType = 'image/png',
                        width = input$width,
                        height = input$height,
                        alt = "This is alternate text"))
          }
          print(p1)
          vizi_par(vp_orig)
          
          #check for ppm calc
          if(input$calc_ppm) {
            
            dat=overview_peaks_sel_plot_masked
            x=mz(dat)
            targ_mz<-req(input$x_target)
            x_sel<-subset(x, x>=xlim[1] & x<= xlim[2])
            
            ppm_error<- round(1e6*(x_sel-targ_mz)/targ_mz, 2)
            
            p1_coord <- try(p1[[1]][[1]]$marks$peaks$encoding, silent = TRUE)
            if (inherits(p1_coord, "try-error") || is.null(p1_coord)) {
              showNotification("Could not access plotted peaks for ppm labels.", type = "warning")
              p1_coord <- NULL
            }
            if (is.null(p1_coord)) {
              dev.off()
              return(list(src = outfile,
                          contentType = 'image/png',
                          width = input$width,
                          height = input$height,
                          alt = "This is alternate text"))
            }
            
            y_labs<-p1_coord$y[p1_coord$x %in% x_sel]
            
            if(length(ppm_error)==0){
              showNotification("No ppm error calculated, are there any peaks?", type="warning")
              return()
            } else {
              print(text(x=x_sel, y=y_labs+ylim[2]*.25, req(ppm_error)))
            }
            
          }
          
          
          
          if(input$show_int) {
            
            p1_coord <- try(p1[[1]][[1]]$marks$peaks$encoding, silent = TRUE)
            if (inherits(p1_coord, "try-error") || is.null(p1_coord)) {
              showNotification("Could not access plotted peaks for intensity labels.", type = "warning")
              p1_coord <- NULL
            }
            if (is.null(p1_coord)) {
              dev.off()
              return(list(src = outfile,
                          contentType = 'image/png',
                          width = input$width,
                          height = input$height,
                          alt = "This is alternate text"))
            }
            
            dat=overview_peaks_sel_plot_masked
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
