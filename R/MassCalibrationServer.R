MassCalibrationServer <- function(id, setup_values, preproc_values) {
  moduleServer(id,function(input,output,session) {
    source <- reactiveVal(NULL); result <- reactiveVal(NULL); anchors <- reactiveVal(NULL); fit <- reactiveVal(NULL)
    candidate_data <- reactiveVal(NULL)
    stamp <- reactiveVal(NULL); status <- reactiveVal("Load a processed dataset to begin.")
    guard <- function(expr) tryCatch(force(expr),error=function(e) {status(conditionMessage(e)); showNotification(conditionMessage(e),type="error",duration=10); NULL})
    invalidate <- function() {fit(NULL); result(NULL)}
    observeEvent(list(input$polarity,input$extended,input$references,input$anchors_rows_selected,input$show_all_matches,input$adducts,input$model),{invalidate()},ignoreInit=TRUE)
    observeEvent(list(input$polarity,input$extended,input$references,input$adducts), {
      candidate_data(NULL); anchors(NULL)
    },ignoreInit=TRUE)
    output$adduct_controls <- renderUI({
      req(input$polarity != "unknown")
      refs <- cal_references(input$polarity)
      choices <- sort(unique(refs$adduct))
      defaults <- setdiff(choices,c("[M+H-H2O]+","[M-H-H2O]-","[M+HCOO]-"))
      checkboxGroupInput(session$ns("adducts"),"Ion adducts to search",choices=choices,selected=defaults)
    })
    anchor_view <- reactive({
      a <- anchors(); req(a)
      if(!isTRUE(input$show_all_matches)) a <- a[a$use,,drop=FALSE]
      a
    })
    # Match File Restore: discover saved files recursively and poll every ten seconds.
    directory_files <- function() {
      root <- setup_values()[["wd"]]
      list(root=root, files=list.files(root, recursive=TRUE))
    }
    my_files <- reactivePoll(10000,session,checkFunc=directory_files,valueFunc=directory_files)
    dataset_files <- reactive({
      files <- my_files()$files
      files[grepl("\\.(imzML|rds|RData)$",basename(files),ignore.case=TRUE) &
            !grepl("_calibration[.]rds$",basename(files),ignore.case=TRUE)]
    })
    observeEvent(dataset_files(), {
      files <- dataset_files()
      selected <- isolate(input$dataset_file)
      if(is.null(selected) || !selected %in% files) selected <- if(length(files)) files[1] else character(0)
      updateSelectInput(session,"dataset_file",choices=files,selected=selected)
    },ignoreNULL=FALSE)
    observeEvent(input$load,guard({
      selected <- input$dataset_file
      if(is.null(selected) || !nzchar(selected)) stop("Choose a saved dataset first")
      if(!selected %in% dataset_files()) stop("Selected dataset is no longer in the working directory")
      p <- file.path(setup_values()[["wd"]],selected)
      x <- if(grepl("\\.(rds|RData)$",p,ignore.case=TRUE)) readRDS(p) else Cardinal::readMSIData(p)
      x <- msi_clean_calibration_metadata(x)
      if(is.null(x)) stop("No dataset available from this source")
      m <- as.numeric(Cardinal::mz(x))
      if(!length(m) || any(!is.finite(m)) || any(diff(m)<=0)) stop("A processed dataset with a unique ordered feature axis is required")
      source(x); invalidate(); anchors(NULL); stamp(cal_signature(x))
      status(sprintf("Loaded %d features, %d pixels, %d runs.",nrow(x),ncol(x),length(Cardinal::runNames(x))))
    }))
    observeEvent(input$find,guard({
      req(source()); invalidate()
      if(input$polarity=="unknown") stop("Select the acquisition polarity")
      if(msi_is_calibrated(source())) stop("Restore original masses before fitting again")
      refs <- if(is.null(input$references)) cal_references(input$polarity) else read.delim(input$references$datapath,stringsAsFactors=FALSE)
      if(is.null(input$adducts) || !length(input$adducts)) stop("Select at least one ion adduct")
      if(!"adduct" %in% names(refs)) stop("Reference table is missing adducts")
      if(any(grepl(if(input$polarity=="positive") "-$" else "\\+$",refs$adduct))) stop("Reference polarity conflicts with acquisition polarity")
      refs <- cal_prepare_references(refs,input$adducts)
      withProgress(message="Summarizing feature intensities",{
        s <- Cardinal::summarizeFeatures(source(),stat=c(mean="mean"),verbose=FALSE)
        a <- cal_discover(as.numeric(Cardinal::mz(source())),as.numeric(Cardinal::fData(s)$mean),refs,input$extended)
        candidate_data(list(mass=as.numeric(Cardinal::mz(source())),abundance=as.numeric(Cardinal::fData(s)$mean),refs=refs,extended=input$extended))
        anchors(a); status(sprintf("%d eligible / %d matched ion compositions. Recommended rows are selected across the mass range. Enable weaker/ambiguous matches to inspect the rest.",sum(a$use),nrow(a)))
      })
    }))
    observeEvent(input$fit,guard({
      req(anchors(),source()); result(NULL)
      d <- candidate_data(); req(d)
      a <- cal_selected_candidates(anchor_view(),input$anchors_rows_selected)
      f <- cal_fit(a,input$extended)
      if(input$model != "auto") {
        if(!input$model %in% f$models$model || !is.finite(f$models$loo_rmse_ppm[match(input$model,f$models$model)])) stop("Selected model is unsupported by these anchors/settings")
        aa <- f$anchors[f$anchors$use,,drop=FALSE]
        f$manual_model <- TRUE; f$model <- input$model; f$coef <- cal_fit_one(aa$obs,aa$ppm_raw,f$model)
      }
       f$source_signature <- stamp(); f$polarity <- input$polarity
      f$created <- Sys.time()
      f$library_version <- unique(d$refs$library_version); f$adducts <- input$adducts
      f$library_references <- d$refs
      fit(f); status(paste("Recommended model:",f$model,". Review LOO scores and mass coverage; Apply creates a calibrated copy."))
    }))
    active_fit <- reactive({
      f <- fit(); req(f)
      rows <- input$models_rows_selected
      model <- f$model
      if(length(rows)) {
        validate(need(length(rows)==1L && !is.na(rows) && rows>=1L && rows<=nrow(f$models),"Select one model row"))
        model <- f$models$model[rows]
      }
      validate(need(is.finite(f$models$loo_rmse_ppm[match(model,f$models$model)]),"This model is unavailable: insufficient independent anchors or a singular fit. Select an eligible model."))
      cal_choose_model(f,model)
    })
    observeEvent(input$models_rows_selected, {
      result(NULL)
      if(!is.null(fit())) status("Model selection changed. The plot reflects the selected model; Apply to new dataset uses that model. Reapply before saving.")
    },ignoreInit=TRUE)
    observeEvent(input$apply,guard({
      req(active_fit(),source()); if(!identical(cal_signature(source()),active_fit()$source_signature)) stop("Source changed; refit required")
      result(cal_apply(source(),active_fit())); plot_card_server("preview",result())
      status(paste("Applied",active_fit()$model,"to a new dataset. Save or open it in Peak Picking."))
    }))
    observeEvent(input$undo,guard({
      req(source()); x <- if(!is.null(result())) result() else source()
      x <- cal_restore(x); source(x); result(NULL); anchors(NULL); fit(NULL); stamp(cal_signature(x))
      plot_card_server("preview",x); status("Original masses restored; refit before applying another correction.")
    }))
    observeEvent(input$open,guard({req(result()); preproc_values()$replace_dataset(result()); status("Calibrated result opened in Peak Picking. Existing mass-based identifications require review.")}))
    save_roots <- reactive(c(wd=setup_values()[["wd"]],home=fs::path_home()))
    observe({
      shinyFiles::shinyFileSave(input,"save_imzml",roots=save_roots(),session=session)
    })
    observeEvent(input$save_imzml,guard({
      req(result(),active_fit())
      save_path <- shinyFiles::parseSavePath(save_roots(),input$save_imzml)
      if(!nrow(save_path)) return(invisible(NULL))
      folder <- withProgress(message="Saving calibrated bundle",{
        cal_save_bundle(result(),active_fit(),as.character(save_path$datapath))
      })
      status(paste("Saved calibrated dataset and reports in",folder))
    }))
    output$status <- renderText(status())
    output$anchors <- DT::renderDataTable({
      a <- anchor_view(); display <- cal_anchor_display(a)
      dominance_column <- match("dominance",names(display)) # DT includes row names as column zero.
      DT::datatable(display,selection=list(mode="multiple",selected=which(a$recommended & a$use),target="row"),
        options=list(scrollX=TRUE,pageLength=10,columnDefs=list(list(targets=dominance_column,
          render=DT::JS("function(data,type,row,meta) { if(type === 'sort' || type === 'type') { if(data === 'No competing peak') return Infinity; if(data === 'Unavailable') return -Infinity; return parseFloat(data); } return data; }")))),
        caption="Mean intensity is averaged across the loaded dataset's pixels, in its existing intensity units. Dominance is the ratio to the next strongest peak in the search window; No competing peak means that ratio is unbounded. Click rows to change the recommended selection.")
    })
    output$models <- DT::renderDataTable({req(fit()); DT::datatable(fit()$models,selection=list(mode="single",selected=match(fit()$model,fit()$models$model),target="row"),caption="Click a model to display its correction and residuals. Apply uses the selected model.")})
    output$errors <- renderPlot({
      req(active_fit()); f <- active_fit(); a <- f$anchors[f$anchors$use,,drop=FALSE]
      plot(a$obs,a$ppm_raw,xlab="Observed m/z",ylab="Mass error (ppm)",pch=16,main=paste("Selected model:",f$model))
      abline(h=0,lty=2); m <- seq(min(a$obs),max(a$obs),length.out=300)
      lines(m,cal_predict(m,f$model,f$coef),col="steelblue",lwd=2)
      points(a$obs,a$ppm_raw-cal_predict(a$obs,f$model,f$coef),col="darkgreen",pch=1)
      legend("topright",c("Before","After (fit residual)","Correction"),col=c("black","darkgreen","steelblue"),pch=c(16,1,NA),lty=c(NA,NA,1))
    })
    reactive(result())
  })
}
