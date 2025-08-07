### UMAPEmbeddingServer.R 
UMAPEmbeddingServer <- function(id, setup_values, preproc_values) {
  moduleServer(id, function(input, output, session) {
    ns = session$ns
    
    output$variables <- renderPrint({
      setup_vals <- setup_values()
      paste("working directory=", setup_vals$wd)
    })
    
    # List runs for selection - matches UMAPServer.r exactly
    output$mytable = DT::renderDataTable({
      x2 <- preproc_values()[["x2"]]
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
    data <- reactive({
      x2 <- preproc_values()[["x2"]]
      req(x2$overview_peaks_sel)
      req(input$mytable_rows_selected)
      
      ids <- input$mytable_rows_selected
      x2$overview_peaks_sel %>% 
        subsetPixels(Cardinal::run(x2$overview_peaks_sel) %in% 
                       runNames(x2$overview_peaks_sel)[ids])
    })
    
    # Simple check for UMAP coordinates 
    has_umap <- reactive({
      req(data())
      all(c("UMAP1", "UMAP2", "UMAP3") %in% colnames(pData(data())))
    })
    
    # Get UMAP coordinates
    umap_coords <- reactive({
      req(data(), has_umap())
      if (!has_umap()) {
        showNotification("No UMAP coordinates found in data. Run UMAP first.", type = "error")
        return(NULL)
      }
      as.matrix(pData(data())[, c("UMAP1", "UMAP2", "UMAP3")])
    })
    
    # UI elements
    output$mz_selector <- renderUI({
      req(data())
      
      mz_values <- fData(data())$mz
      names(mz_values) <- paste("m/z:", round(mz_values, 4))
      
      selectizeInput(ns("mz_ion"), 
                     "Select m/z Ion:",
                     choices = mz_values, 
                     selected = mz_values[1],
                     options = list(maxOptions = 10000))
    })
    
    output$color_selector <- renderUI({
      req(data())
      
      color_vars <- colnames(pData(data()))
      
      selectInput(ns("color_var"),
                  "Select Color Variable:", 
                  choices = color_vars,
                  selected = "Rcol_reduced")
    })
    
    # Selected m/z index
    selected_index <- reactive({
      req(data(), input$mz_ion)
      which(fData(data())$mz == as.numeric(input$mz_ion))
    })
    
    # UMAP plot with m/z intensity
    output$umapPlot <- renderImage({
      req(data(), selected_index(), umap_coords())
      
      outfile <- tempfile(fileext = '.png')
      png(outfile, width = 800, height = 600)
      
      # Get intensity values
      intensity_values <- spectra(data())[selected_index(), ]
      
      if (input$log_scale) {
        intensity_values <- log(intensity_values + 1)
      }
      
      normalized_intensities <- (intensity_values - min(intensity_values)) / 
        (max(intensity_values) - min(intensity_values))
      
      # Create color palette and map values
      color_palette <- pals::viridis(100)
      midpoint <- input$color_midpoint
      rescaled_numbers <- (normalized_intensities - midpoint) / (1 - midpoint) + 0.5
      rescaled_numbers <- pmax(0, pmin(1, rescaled_numbers))
      color_breaks <- seq(0, 1, length.out = length(color_palette) + 1)
      color_values <- color_palette[findInterval(rescaled_numbers, color_breaks)]
      
      print(pairs(umap_coords(),
                  col = color_values,
                  pch = ".",
                  cex = input$point_size))
      
      dev.off()
      
      list(src = outfile,
           contentType = 'image/png',
           width = 800,
           height = 600,
           alt = "UMAP plot colored by m/z intensity")
    }, deleteFile = TRUE)
    
    # UMAP plot with metadata colors
    output$umapPlot2 <- renderImage({
      req(data(), input$color_var, umap_coords())
      
      outfile <- tempfile(fileext = '.png')
      png(outfile, width = 800, height = 600)
      
      if(input$col_match) {
        color_values <- pData(data())[[input$color_var]]
      } else {
        unique_values <- unique(pData(data())[[input$color_var]])
        tol_palette <- pals::tol(max(12, length(unique_values)))
        color_ramp <- colorRampPalette(tol_palette)
        factor_levels <- as.numeric(factor(pData(data())[[input$color_var]]))
        color_values <- color_ramp(length(unique_values))[factor_levels]
      }
      
      print(pairs(umap_coords(),
                  col = color_values,
                  pch = ".",
                  cex = input$point_size))
      
      dev.off()
      
      list(src = outfile,
           contentType = 'image/png',
           width = 800,
           height = 600,
           alt = "UMAP plot colored by metadata")
    }, deleteFile = TRUE)
    
    # Phenotype plot
    output$phenoplot <- renderImage({
      req(data())
      
      outfile <- tempfile(fileext = '.png')
      png(outfile, width = 800, height = 600)
      
      if(input$color_var != "" || input$plot_mz) {
        if(input$plot_mz) {
          print(Cardinal::image(data(),
                                mz = as.numeric(input$mz_ion),
                                tol = 0.005,
                                units = "mz",
                                enhance = "hist",
                                smooth = "gaussian"))
        } else {
          print(Cardinal::image(data(),
                                input$color_var,
                                col = pals::tol()))
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
