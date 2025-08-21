### R/CorrelationServer.R
CorrelationServer <- function(id, proc_values, setup_values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns # For dynamic variable namespace
    
    # Reactive expression to access x5 from proc_values
    x5 <- reactive({
      req(proc_values()[["x5"]])
      proc_values()[["x5"]]
    })
    
    # Functions to check and get files from the working directory
    has_new_files <- function() {
      unique(list.files(req(setup_values())[["wd"]], recursive = TRUE))
    }
    
    get_files <- function() {
      list.files(req(setup_values())[["wd"]], recursive = TRUE)
    }
    
    # Reactive polling to detect new files every 10 seconds (10000 ms)
    my_files <- reactivePoll(10000, session, checkFunc = has_new_files, valueFunc = get_files)
    
    # Update the 'correlation_file' selectInput when new files are detected
    observeEvent(my_files(), ignoreInit = TRUE, ignoreNULL = TRUE, {
      file_choices <- grep("\\.imzML$|\\.rds$", my_files(), ignore.case = TRUE, value = TRUE)
      updateSelectInput(session, ns('correlation_file'), choices = file_choices)
    })
    
    # ReactiveValues to store data_file and hmap_choices
    x6 <- reactiveValues(
      data_file = NULL,
      hmap_choices = NULL
    )
    
    ### Render UI for Source Selection ###
    output$source_ui <- renderUI({
      if (input$corr_source == "from_file") {
        tagList(
          selectInput(
            ns("correlation_file"),
            "Select File for Correlation Analysis",
            choices = grep("\\.imzML$|\\.rds$", my_files(), ignore.case = TRUE, value = TRUE),
            selected = NULL
          ),
          actionButton(ns("action_seg"), label = HTML("Read Selected File"))
        )
      } else {
        # If 'from_stats' is selected, no additional UI elements are needed here
        NULL
      }
    })
    
    ### Observe Event for Reading File ###
    observeEvent(input$action_seg, {
      req(input$correlation_file) # Ensure a file is selected
      
      # Construct full file path
      file_path <- file.path(setup_values()[["wd"]], input$correlation_file)
      
      # Read the selected file based on its extension
      data_file <- NULL
      if (grepl("\\.rds$", input$correlation_file, ignore.case = TRUE)) {
        data_file <- tryCatch({
          readRDS(file_path)
        }, error = function(e) {
          showNotification(paste("Error reading .rds file:", e$message), type = "error")
          NULL
        })
      } else if (grepl("\\.imzML$", input$correlation_file, ignore.case = TRUE)) {
        data_file <- tryCatch({
          readMSIData(file_path)
        }, error = function(e) {
          showNotification(paste("Error reading .imzML file:", e$message), type = "error")
          NULL
        })
      } else {
        showNotification("Unsupported file type selected.", type = "error")
      }
      
      if (!is.null(data_file)) {
        x6$data_file <- data_file
        showNotification("Data file loaded successfully.", type = "message")
      }
    })
    
    ### Observe Event for Source Selection ###
    observeEvent(input$corr_source, {
      if (input$corr_source == "from_stats") {
        req(x5()$data_file_selected)
        x6$data_file <- x5()$data_file_selected
        
        # Extract hmap_choices from stats_results, excluding 'mz' and 'feature'
        x6$hmap_choices <- setdiff(colnames(x5()$stats_results), c("mz", "feature"))
      } else {
        # If 'from_file' is selected, clear hmap_choices
        x6$hmap_choices <- NULL
      }
    }, ignoreNULL = FALSE) # Trigger on initial load as well
    
    ### Render UI for Correlation Variables ###
    output$corr_vars <- renderUI({
      if (input$corr_source == "from_stats") {
        req(x6$hmap_choices)
        tagList(
          selectInput(
            ns("corr_sig_select"),
            "Choose Filtering Column for Correlation Choices",
            choices = x6$hmap_choices,
            selected = "fdr"
          ),
          radioButtons(
            ns("corr_sig_direction"),
            label = "Filtering Statistic Direction:",
            choices = list("Ascending (≤)" = "ascending", "Descending (≥)" = "descending"),
            selected = "descending"
          ),
          numericInput(ns("corr_sig"), "Significance Cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
          checkboxInput(ns("key_on"), "Display Color Key", value = TRUE)
        )
      } else {
        # If 'from_file', no additional correlation variables are needed here
        NULL
      }
    })
    
    ### Render UI for m/z Selection ###
    output$corr_mz <- renderUI({
      if (input$corr_source == "from_stats") {
        req(x5()$stats_results)
        req(input$corr_sig_select)
        req(input$corr_sig_direction)
        req(input$corr_sig)
        
        # Filter m/z based on statistical test results
        if (input$corr_sig_direction == "descending") {
          sig_mz <- x5()$stats_results %>%
            dplyr::filter(.data[[input$corr_sig_select]] >= input$corr_sig) %>%
            dplyr::pull(mz)
        } else {
          sig_mz <- x5()$stats_results %>%
            dplyr::filter(.data[[input$corr_sig_select]] <= input$corr_sig) %>%
            dplyr::pull(mz)
        }
        
        # Remove duplicates and sort
        sig_mz <- unique(sort(sig_mz))
        
        selectInput(
          ns("corr_mz"),
          "Choose m/z Value for Correlation Analysis",
          choices = sig_mz,
          selected = NULL
        )
      } else if (input$corr_source == "from_file") {
        req(x6$data_file)
        
          mz_values <- mz(x6$data_file)
        
        selectInput(
          ns("corr_mz"),
          "Choose m/z Value for Correlation Analysis",
          choices = mz_values,
          selected = NULL
        )
      } else {
        NULL
      }
    })
    
    ### Render Correlation Plot ###
    observeEvent(input$action_corr, {
      # Validate inputs based on source
      if (input$corr_source == "from_stats") {
        req(input$corr_mz)
      } else if (input$corr_source == "from_file") {
        req(input$corr_mz)
      } else {
        showNotification("Invalid data source selected.", type = "error")
        return(NULL)
      }
      
      output$plot15 <- renderImage({
        req(x6$data_file)
        req(input$corr_mz)
        req(input$corr_n)
        req(input$corr_sort)
        req(input$corr_plot_n)
        
        # Get the selected m/z or i value
        selected_mz <- input$corr_mz
        
        # Check if selected_mz is valid
        if (is.null(selected_mz) || selected_mz == "") {
          showNotification("Please select a valid m/z or i value for correlation.", type = "error")
          return(NULL)
        }
      
        # Generate the correlation data using colocalized
        coloc <- tryCatch({
          colocalized(
            x6$data_file,
            mz = as.numeric(selected_mz),
            n = input$corr_n,
            sort.by = input$corr_sort
          )
        }, error = function(e) {
          showNotification(paste("Error in correlation computation:", e$message), type = "error")
          NULL
        })
        
        if (is.null(coloc)) {
          return(NULL)
        }
        
        print(coloc)
        # Ensure 'coloc$mz' is available
        if (!"mz" %in% names(coloc)) {
          showNotification("Correlation result does not contain 'mz' values.", type = "error")
          return(NULL)
        }
        
        # Generate the plot
        outfile <- tempfile(fileext = '.png')
        png(outfile, width = 800, height = 600)
        
        # Determine m/z values to plot based on 'corr_plot_n'
        plot_n <- min(input$corr_plot_n, length(coloc$mz))
        plot_mz <- coloc$mz[1:plot_n]
        
        # Generate the image using Cardinal's image function
        p1 <- tryCatch({
          image(
            x6$data_file,
            mz = plot_mz,
            layout = c(length(runNames(x6$data_file)), plot_n),
            enhance = "histogram",
            scale=T,
            col = cpal("Spectral"),
            colorkey = input$key_on
          )
        }, error = function(e) {
          showNotification(paste("Error in plotting:", e$message), type = "error")
          NULL
        })
        
        if (!is.null(p1)) {
          print(p1)
        }
        
        dev.off()
        
        # Return a list containing the filename
        list(
          src = outfile,
          contentType = 'image/png',
          width = 800,
          height = 600,
          alt = "Correlation Heatmap"
        )
      }, deleteFile = TRUE)
    })
  })
}
