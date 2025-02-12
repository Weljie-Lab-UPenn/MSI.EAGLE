### R/UMAPEmbeddingServer.R 

UMAPEmbeddingServer <- function(id, setup_values, preproc_values, preproc_values_umap) {
  moduleServer(id, function(input, output, session) {
    ns = session$ns
    
    # Add reactive values to store color selections
    color_state <- reactiveValues(
      colors = NULL,
      selected_colors = NULL
    )
    
    # List runs for selection
    output$mytable = DT::renderDataTable({

      
      
      x2 <- try(preproc_values_umap()[["x2"]])
      if(inherits(x2, "try-error")) {
        x2<-preproc_values()[["x2"]]
      }
      #browser()
      
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
      
      x2 <- try(preproc_values_umap()[["x2"]])
      if(inherits(x2, "try-error")) {
        x2<-preproc_values()[["x2"]]
      }
      req(x2$overview_peaks_sel)
      req(input$mytable_rows_selected)
      #browser()
      
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
      color_values <- pData(data_img())[[input$color_var]]
      mask <- color_values %in% color_state$selected_colors
      
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
        normalized_intensities <- (selected_intensities - min(selected_intensities)) / 
          (max(selected_intensities) - min(selected_intensities))
        
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
      unique_values <- unique(pData(data_img())[[input$color_var]])
      
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
      }
      
      # Display the color choices with checkboxes, defaulting to all selected
      checkboxGroupInput(ns("selected_colors"),
                         "Select colors to display:",
                         choices = colors,
                         selected = color_state$selected_colors)
    })
    
    
    # Update selected colors when changed
    observeEvent(input$selected_colors, {
      color_state$selected_colors <- input$selected_colors
    })
    
    # Modified UMAP plot with metadata colors and selection
    output$umapPlot2 <- renderImage({
      req(data_img(), input$color_var, umap_coords(), color_state$selected_colors)
      
      outfile <- tempfile(fileext = '.png')
      png(outfile, width = 800, height = 600)
      
      # Get metadata values
      color_values <- pData(data_img())[[input$color_var]]
      
      # Create mask for selected colors
      mask <- color_values %in% color_state$selected_colors
      
      if(input$col_match) {
        # Use original colors for selected values, grey for unselected
        plot_colors <- rep("grey80", length(color_values))
        plot_colors[mask] <- color_values[mask]
      } else {
        unique_values <- unique(color_values[mask])
        tol_palette <- pals::tol(max(12, length(unique_values)))
        color_ramp <- colorRampPalette(tol_palette)
        
        # Create factor levels only for selected values
        plot_colors <- rep("grey80", length(color_values))
        if(length(unique_values) > 0) {
          factor_levels <- as.numeric(factor(color_values[mask], levels = unique_values))
          plot_colors[mask] <- color_ramp(length(unique_values))[factor_levels]
        }
      }
      #browser()
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
          # Create a subset of data containing only selected colors
          
          #browser()
          img_vec <- pData(data_img())[[input$color_var]]
          select_vec <- color_state$selected_colors
          
          # Check if select_vec only contains "TRUE" and "FALSE"
          if (all(select_vec %in% c("TRUE", "FALSE"))) {
            # Convert to logical
            select_vec <- as.logical(select_vec)
          } else {
            select_vec <- as.character(select_vec)
          }
          
          if (is.logical(img_vec) && is.logical(select_vec) ||
              (is.character(img_vec) && is.character(select_vec) )) {
            
              subset_data <- data_img() %>%
                subsetPixels(img_vec %in% select_vec)
          } else {
            return()
          }
          
          if(ncol(subset_data) > 0) {
            print(Cardinal::image(subset_data,
                                  input$color_var,
                                  col = pals::tol()))
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
