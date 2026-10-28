MzIdentificationServer <- function(id, setup_values, calibration_values) {
  moduleServer(id, function(input, output, session) {
    matches <- reactiveVal(data.frame()); accepted <- reactiveVal(data.frame()); library_revision <- reactiveVal(0L); saved_data <- reactiveVal(NULL); status <- reactiveVal("Apply a Mass Calibration result or load a saved calibrated dataset to begin.")
    current_data <- reactive({
      x <- saved_data()
      if (is.null(x)) x <- calibration_values()
      validate(need(!is.null(x), "Apply Mass Calibration or load a saved calibrated dataset before m/z Identification."))
      x
    })
    directory_files <- function() {
      root <- setup_values()[["wd"]]
      list(root = root, files = list.files(root, recursive = TRUE))
    }
    working_files <- reactivePoll(10000, session, checkFunc = directory_files, valueFunc = directory_files)
    dataset_files <- reactive(msi_dataset_choices(working_files()$files))
    observeEvent(dataset_files(), {
      choices <- dataset_files(); selected <- isolate(input$dataset_file)
      if (is.null(selected) || !selected %in% choices) selected <- if (length(choices)) choices[1] else character(0)
      updateSelectInput(session, "dataset_file", choices = choices, selected = selected)
    }, ignoreNULL = FALSE)
    observeEvent(input$load_dataset, {
      tryCatch({
        selected <- input$dataset_file
        if (is.null(selected) || !nzchar(selected) || !selected %in% dataset_files()) stop("Choose a saved dataset in the working directory first")
        x <- mzid_load_dataset(file.path(setup_values()[["wd"]], selected))
        saved_data(x); matches(data.frame()); accepted(mzid_accepted_from_dataset(x))
        status(sprintf("Loaded saved calibrated dataset: %s (%d features, %d pixels).", selected, nrow(x), ncol(x)))
      }, error = function(e) { status(conditionMessage(e)); showNotification(conditionMessage(e), type = "error", duration = 10) })
    })
    observeEvent(input$use_calibration_result, {
      saved_data(NULL); matches(data.frame()); accepted(if (is.null(calibration_values())) data.frame() else mzid_accepted_from_dataset(calibration_values()))
      if (is.null(calibration_values())) status("No active Mass Calibration result is available.") else status("Using the current Mass Calibration result.")
    })
    library_state <- reactive({ library_revision(); mzid_available(setup_values()[["wd"]]) })
    output$library_status <- renderUI({ x <- library_state(); tags$ul(lapply(seq_len(nrow(x)), function(i) tags$li(paste(x$label[i], if (x$installed[i]) "installed" else "not installed")))) })
    output$library_select <- renderUI({ x <- library_state(); choices <- setNames(x$id[x$installed], x$label[x$installed]); checkboxGroupInput(session$ns("libraries"), "Databases to include", choices = choices) })
    output$adducts <- renderUI({ a <- mzid_adducts(input$polarity); checkboxGroupInput(session$ns("selected_adducts"), "Adducts", choices = a$adduct, selected = a$adduct) })
    install <- function(name) observeEvent(input[[paste0("install_", name)]], {
      installed <- tryCatch(withProgress(message = paste("Installing", name), value = 0.1, { m <- mzid_install(setup_values()[["wd"]], name, function(detail) setProgress(0.6, detail = detail)); setProgress(1); m }), error = function(e) { status(conditionMessage(e)); showNotification(conditionMessage(e), type = "error", duration = 12); NULL })
      if (!is.null(installed)) { library_revision(library_revision() + 1L); status(paste("Installed", mzid_specs()[[name]]$label, "locally.")) }
    })
    install("hmdb")
    observeEvent(input$open_hmdb, {
      utils::browseURL("https://hmdb.ca/downloads")
      status("Opened the official HMDB download page. After its browser challenge/download completes, click 'Install HMDB from Downloads'.")
    })
    observeEvent(input$install_hmdb_downloaded, {
      installed <- tryCatch(withProgress(message = "Installing downloaded HMDB", value = 0.1, {
        archive <- mzid_downloaded_hmdb_archive(setup_values()[["wd"]])
        m <- mzid_install_archive(setup_values()[["wd"]], "hmdb", archive, function(detail) setProgress(0.6, detail = detail)); setProgress(1); m
      }), error = function(e) { status(conditionMessage(e)); showNotification(conditionMessage(e), type = "error", duration = 12); NULL })
      if (!is.null(installed)) { library_revision(library_revision() + 1L); status("Installed HMDB from the official downloaded archive.") }
    })
    observeEvent(input$open_lmsd, {
      utils::browseURL("https://www.lipidmaps.org/databases/lmsd/download")
      status("Opened the official LIPID MAPS LMSD download page. Download the LMSD Structure-data SDF ZIP, then click 'Install LIPID MAPS from Downloads'.")
    })
    observeEvent(input$install_lmsd_downloaded, {
      installed <- tryCatch(withProgress(message = "Installing downloaded LIPID MAPS LMSD", value = 0.1, {
        archive <- mzid_downloaded_lmsd_archive(setup_values()[["wd"]])
        m <- mzid_install_archive(setup_values()[["wd"]], "lmsd", archive, function(detail) setProgress(0.6, detail = detail)); setProgress(1); m
      }), error = function(e) { status(conditionMessage(e)); showNotification(conditionMessage(e), type = "error", duration = 12); NULL })
      if (!is.null(installed)) { library_revision(library_revision() + 1L); status("Installed LIPID MAPS LMSD from the official downloaded archive.") }
    })
    observeEvent(input$match, {
      tryCatch({ d <- current_data(); ids <- input$libraries; if (!length(ids)) stop("Install and select at least one database")
        x <- mzid_match(as.numeric(Cardinal::mz(d)), ids, setup_values()[["wd"]], input$polarity, input$selected_adducts, input$ppm)
        matches(x); status(sprintf("Matched %d candidate records across all %d selected database(s).", nrow(x), length(ids)))
      }, error = function(e) { status(conditionMessage(e)); showNotification(conditionMessage(e), type = "error", duration = 10) })
    })
    features <- reactive({ mzid_feature_table(current_data(), matches(), accepted()) })
    output$features <- DT::renderDataTable(DT::datatable(features(), selection = "single", options = list(scrollX = TRUE, pageLength = 12), caption = "Corrected m/z feature summary. All selected databases contribute candidates."))
    selected_feature <- reactive({ rows <- input$features_rows_selected; if (!length(rows)) return(NULL); features()$feature_index[rows[1]] })
    candidate_view <- reactive({ f <- selected_feature(); x <- matches(); if (is.null(f) || !nrow(x)) return(data.frame()); x[x$feature_index == f, c("library_id", "accession", "name", "formula", "adduct", "theoretical_mz", "ppm_error", "compound_class", "candidate_rank"), drop = FALSE] })
    output$candidates <- DT::renderDataTable(DT::datatable(candidate_view(), selection = "single", options = list(scrollX = TRUE, pageLength = 10), caption = "Candidates are MS1 mass matches; source-library records remain distinct."))
    observeEvent(input$accept, {
      f <- selected_feature(); r <- input$candidates_rows_selected; x <- matches(); if (is.null(f) || !length(r)) { showNotification("Select a feature and one candidate.", type = "warning"); return() }
      pool <- x[x$feature_index == f, , drop = FALSE]; chosen <- pool[r[1], , drop = FALSE]; old <- accepted(); old <- old[old$feature_index != f, , drop = FALSE]; accepted(rbind(old, chosen)); status(paste("Accepted tentative annotation:", chosen$name));
    })
    save_roots <- function() isolate(c(wd = setup_values()[["wd"]], home = path.expand("~")))
    shinyFiles::shinyFileSave(input, "save_annotated", roots = save_roots(), session = session)
    observeEvent(input$save_annotated, {
      tryCatch({
        if (!nrow(accepted())) stop("Accept at least one annotation before saving an annotated dataset")
        save_path <- shinyFiles::parseSavePath(save_roots(), input$save_annotated)
        if (!nrow(save_path)) return(invisible(NULL))
        folder <- withProgress(message = "Saving annotated bundle", value = 0.1, {
          out <- mzid_save_bundle(current_data(), accepted(), matches(), input$libraries %||% character(), setup_values()[["wd"]], as.character(save_path$datapath)); setProgress(1); out
        })
        status(sprintf("Saved %d curator-accepted annotations in %s.", nrow(accepted()), folder))
      }, error = function(e) { status(conditionMessage(e)); showNotification(conditionMessage(e), type = "error", duration = 12) })
    })
    output$status <- renderText(status())
    reactive({ d <- current_data(); list(data = mzid_attach_accepted(d, accepted()), candidates = matches(), accepted = accepted(), feature_table = features()) })
  })
}
