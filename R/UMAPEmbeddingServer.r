### R/UMAPEmbeddingServer.R 

UMAPEmbeddingServer <- function(id, setup_values, preproc_values, preproc_values_umap = NULL) {
  moduleServer(id, function(input, output, session) {
    ns = session$ns
    
    # Add reactive values to store color selections
    color_state <- reactiveValues(
      colors = NULL,
      selected_colors = NULL
    )

    get_x2 <- function() {
      if (is.function(preproc_values_umap)) {
        x_umap <- try(preproc_values_umap(), silent = TRUE)
        if (!inherits(x_umap, "try-error") &&
            !is.null(x_umap) &&
            !is.null(x_umap[["x2"]])) {
          return(x_umap[["x2"]])
        }
      }
      preproc_values()[["x2"]]
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

    build_display_mapping <- function(labels, match_original = FALSE, palette = "Dark 3") {
      labels <- normalize_labels(labels)
      uniq <- unique(labels)
      fill_values <- stats::setNames(rep("grey70", length(uniq)), uniq)

      # Always honor literal color labels (e.g. "darkseagreen3", "#33AA88"),
      # regardless of the "match original" toggle.
      valid <- are_valid_colors(uniq)
      if (any(valid)) {
        fill_values[valid] <- uniq[valid]
      }
      remaining <- uniq[!valid]

      if (length(remaining) > 0) {
        base_pal <- grDevices::hcl.colors(12, palette = palette)
        gen_pal <- grDevices::colorRampPalette(base_pal)(length(remaining))
        fill_values[remaining] <- gen_pal
      }

      list(labels = labels, values = fill_values)
    }
    
    # List runs for selection
    output$mytable = DT::renderDataTable({
      x2 <- get_x2()
      
      req(x2$overview_peaks_sel)
      ovps <- x2$overview_peaks_sel
      
      DT::datatable(
        cbind(
          run = runNames(ovps),
          size = lapply(runNames(ovps), function(x)
            ncol(ovps[, Cardinal::run(ovps) %in% x]))
        ),
        selection = 'none',
        caption = "Choose runs to process for UMAP visualization",
        extensions = c("Buttons", "Select"),
        options = list(
          dom = 'Bfrtip',
          select = TRUE,
          buttons = list('pageLength', "copy", "selectNone", "selectAll")
        )
      )
    }, server = FALSE)
    
    # Get active dataset based on table selection
    data_img <- reactive({
      x2 <- get_x2()
      req(x2$overview_peaks_sel)
      req(input$mytable_rows_selected)
      
      ids <- input$mytable_rows_selected
      x2$overview_peaks_sel %>% 
        subsetPixels(Cardinal::run(x2$overview_peaks_sel) %in% 
                       runNames(x2$overview_peaks_sel)[ids])
    })
    
    # Simple check for UMAP coordinates 
    has_umap <- reactive({
      #browser()
      req(data_img())
      has_umap<-all(c("x_umap", "y_umap", "z_umap") %in% colnames(pData(data_img())))
      if(!has_umap){
        showNotification("No UMAP coordinates found in data. Run UMAP first.", type = "error")
      }
      return(has_umap)
      
    })
    
    # Get UMAP coordinates
    umap_coords <- reactive({
      req(data_img(), has_umap())
      if (!has_umap()) {
        showNotification("No UMAP coordinates found in data. Run UMAP first.", type = "error")
        return(NULL)
      }
      
      #browser()
      as.matrix(pData(data_img())[, c("x_umap", "y_umap", "z_umap")])
    })
    
    # UI elements
    output$mz_selector <- renderUI({
      req(data_img())
      
      mz_values <- fData(data_img())$mz
      names(mz_values) <- paste("m/z:", round(mz_values, 4))
      
      selectizeInput(ns("mz_ion"), 
                     "Select m/z Ion:",
                     choices = mz_values, 
                     selected = mz_values[1],
                     options = list(maxOptions = 10000))
    })
    
    output$color_selector <- renderUI({
      req(data_img())
      
      color_vars <- colnames(pData(data_img()))
      
      selectInput(ns("color_var"),
                  "Select Color Variable:", 
                  choices = color_vars,
                  selected = "Rcol_reduced")
    })
    
    # Selected m/z index
    selected_index <- reactive({
      req(data_img(), input$mz_ion)
      which(fData(data_img())$mz == as.numeric(input$mz_ion))
    })
    
    # UMAP plot with m/z intensity
    output$umapPlot <- renderImage({
      req(data_img(), selected_index(), umap_coords(), color_state$selected_colors)
      
      outfile <- tempfile(fileext = '.png')
      png(outfile, width = 800, height = 600)
      
      # Get metadata values for masking
      color_values <- normalize_labels(pData(data_img())[[input$color_var]])
      selected_vals <- normalize_labels(color_state$selected_colors)
      mask <- color_values %in% selected_vals
      
      # Get intensity values
      intensity_values <- spectra(data_img())[selected_index(), ]
      
      if (input$log_scale) {
        intensity_values <- log(intensity_values + 1)
      }
      
      # Create color vector starting with all grey
      plot_colors <- rep("grey80", length(intensity_values))
      
      # Only color the selected points based on intensity
      if(any(mask)) {
        # Normalize intensities for selected points only
        selected_intensities <- intensity_values[mask]
        rng <- range(selected_intensities, na.rm = TRUE)
        if (!all(is.finite(rng))) {
          normalized_intensities <- rep(0.5, length(selected_intensities))
        } else if ((rng[2] - rng[1]) <= 0) {
          normalized_intensities <- rep(0.5, length(selected_intensities))
        } else {
          normalized_intensities <- (selected_intensities - rng[1]) / (rng[2] - rng[1])
        }
        
        # Create color palette and map values
        color_palette <- pals::viridis(100)
        midpoint <- input$color_midpoint
        rescaled_numbers <- (normalized_intensities - midpoint) / (1 - midpoint) + 0.5
        rescaled_numbers <- pmax(0, pmin(1, rescaled_numbers))
        color_breaks <- seq(0, 1, length.out = length(color_palette) + 1)
        
        # Assign colors only to selected points
        plot_colors[mask] <- color_palette[findInterval(rescaled_numbers, color_breaks)]
      }
      
      print(pairs(umap_coords(),
                  col = plot_colors,
                  pch = ".",
                  cex = input$point_size))
      
      dev.off()
      
      list(src = outfile,
           contentType = 'image/png',
           width = 800,
           height = 600,
           alt = "UMAP plot colored by m/z intensity")
    }, deleteFile = TRUE)
    
    # Color choices UI
    output$color_choices <- renderUI({
      req(data_img(), input$color_var)
      
      # Get unique values from the selected variable
      unique_values <- unique(normalize_labels(pData(data_img())[[input$color_var]]))
      unique_values <- sort(unique_values)
      
      # Check if the number of unique values exceeds the limit
      if (length(unique_values) > 20) {
        return(tags$div(
          style = "color: red;",
          "The selected variable has more than 20 unique values. Please select a different variable."
        ))
        showNotification("The selected variable has more than 20 unique values. Please select a different variable.", type = "error")
      }
      
      # If within the limit, use unique values as color choices
      colors <- unique_values
      
      # Update color_state with the available colors
      color_state$colors <- colors
      
      # Automatically select all colors if a new variable is chosen
      if (is.null(color_state$selected_colors) || input$color_var != color_state$last_selected_var) {
        color_state$selected_colors <- colors
        color_state$last_selected_var <- input$color_var  # Track the last selected variable
      } else {
        color_state$selected_colors <- intersect(normalize_labels(color_state$selected_colors), colors)
      }
      
      # Display the color choices with checkboxes, defaulting to all selected
      checkboxGroupInput(ns("selected_colors"),
                         "Select colors to display:",
                         choices = colors,
                         selected = color_state$selected_colors)
    })
    
    
    # Update selected colors when changed
    observeEvent(input$selected_colors, {
      color_state$selected_colors <- normalize_labels(input$selected_colors)
    })
    
    # Modified UMAP plot with metadata colors and selection
    output$umapPlot2 <- renderImage({
      req(data_img(), input$color_var, umap_coords(), color_state$selected_colors)
      
      outfile <- tempfile(fileext = '.png')
      png(outfile, width = 800, height = 600)
      
      # Get metadata values
      color_values <- normalize_labels(pData(data_img())[[input$color_var]])
      selected_vals <- normalize_labels(color_state$selected_colors)
      
      # Create mask for selected colors
      mask <- color_values %in% selected_vals
      
      # Use a stable label->color map shared with spatial rendering.
      plot_colors <- rep("grey80", length(color_values))
      if (any(mask)) {
        mapping <- build_display_mapping(color_values[mask], match_original = isTRUE(input$col_match))
        plot_colors[mask] <- unname(mapping$values[mapping$labels])
      }

      print(pairs(umap_coords(),
                  col = plot_colors,
                  pch = ".",
                  cex = input$point_size))
      
      dev.off()
      
      list(src = outfile,
           contentType = 'image/png',
           width = 800,
           height = 600,
           alt = "UMAP plot colored by metadata")
    }, deleteFile = TRUE)
    
    # Modified phenotype plot with color selection
    output$phenoplot <- renderImage({
      req(data_img(), color_state$selected_colors)
      
      outfile <- tempfile(fileext = '.png')
      png(outfile, width = 800, height = 600)
      
      if(input$color_var != "" || input$plot_mz) {
        if(input$plot_mz) {
          print(Cardinal::image(data_img(),
                                mz = as.numeric(input$mz_ion),
                                tol = 0.005,
                                units = "mz",
                                enhance = "hist",
                                smooth = "gaussian"))
        } else {
          img_vec <- normalize_labels(pData(data_img())[[input$color_var]])
          select_vec <- normalize_labels(color_state$selected_colors)
          keep <- img_vec %in% select_vec

          subset_data <- data_img() %>% subsetPixels(keep)
          if (ncol(subset_data) > 0) {
            subset_labels <- normalize_labels(pData(subset_data)[[input$color_var]])
            mapping <- build_display_mapping(subset_labels, match_original = isTRUE(input$col_match))
            pData(subset_data)$.__embed_color_label <- factor(
              mapping$labels,
              levels = names(mapping$values)
            )
            print(Cardinal::image(
              subset_data,
              ".__embed_color_label",
              col = unname(mapping$values)
            ))
          } else {
            graphics::plot.new()
            graphics::text(0.5, 0.5, "No selected labels to display.")
          }
        }
      }
      
      dev.off()
      
      list(src = outfile,
           contentType = 'image/png',
           width = 800,
           height = 600,
           alt = "Phenotype/spatial plot")
    }, deleteFile = TRUE)
    
  })
}
