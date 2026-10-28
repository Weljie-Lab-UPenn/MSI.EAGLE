plot_card_UI<-function(id) {
  
  ns <- NS(id)
  
  col_choices<-c(hcl.pals())
  initial_cols<-c("Inferno","Plasma", "Spectral", "Cividis", "Viridis", 
                  "Zissou 1", "Purple-Green", "Berlin", "PiYG", "Grays", 
                  "Batlow", "turku", "YlOrRd", "Terrain", 
                  "PrGn", "Green-Brown", "Hawaii", "Cork", "Rocket", "RdYlBu")
  
  col_choices<-c(initial_cols, setdiff(col_choices, initial_cols))
  
  tagList(
    
    div(
      class = "msie-panel",
      div(class = "msie-panel-title", "Ion image"),
      
      fluidRow(
        column(
          3, offset = 0,
          class = "msie-opts",
          
          # --- render gate ---------------------------------------------------
          # With "Auto-update" ticked the plot redraws on every control change
          # (original behaviour). Untick it to set up several options -- ions,
          # colorscale, masks -- and redraw once with "Update plot".
          div(
            class = "msie-gate",
            fluidRow(
              column(7, actionButton(ns("update_plot"), "Update plot",
                                     class = "btn-primary", width = "100%")),
              column(5, checkboxInput(ns("auto_update"), "Auto-update", value = TRUE))
            )
          ),
          
          tabsetPanel(
            id = ns("img_opt_tabs"),
            type = "pills",
            
            tabPanel(
              "Ions",
              div(
                class = "msie-opt-body",
                radioButtons(ns("ion_viz3"), "Image visualization ions",
                             c("First ion" = "viz_first",
                               "All ions (TIC)" = "viz_all",
                               "Custom single / multiple" = "custom")),
                # audit PLOTCARD-07: these controls used to be built by a renderUI
                # (mz_viz3a) that re-created them -- resetting them to their defaults --
                # on every dataset change and ion-mode switch. They are now static,
                # with the same ids, labels, choices and defaults (FALSE / "none" /
                # 0.01), shown only in Custom mode. The first draw is unchanged;
                # the user's settings are simply no longer silently reset.
                conditionalPanel(
                  condition = "input.ion_viz3 == 'custom'",
                  ns = ns,
                  checkboxInput(ns("superpose"), "Superpose images?", value=FALSE),
                  selectInput(ns("display_mode"), "Ion math?", c("none", "RGB (2-3 ions)" = "rgb", "sum", "ratio", "subtract", "min", "max", "mean", "sd", "var", "multiply")),
                  numericInput(ns("plusminus_viz3"), "+/- m/z for visualization",0.01)
                )
              )
            ),
            
            tabPanel(
              "Display",
              div(
                class = "msie-opt-body",
                selectInput(ns("color3"), "Colorscale", col_choices, selected = "Plasma"),
                checkboxInput(ns("clip_outliers"), "Clip bright outliers (applied before contrast)", value = FALSE),
                numericInput(ns("clip_pct"), "Clip above percentile", value = 99.9,
                             min = 50, max = 100, step = 0.1),
                checkboxInput(ns("floor_clip"), "Raise floor (clip below percentile)", value = FALSE),
                numericInput(ns("floor_pct"), "Clip below percentile", value = 1,
                             min = 0, max = 50, step = 0.5),
                tags$small(class = "msie-hint",
                           "Percentiles use the pixels on screen -- apply the tissue mask (Pixels tab) first so the background does not count. Most useful with contrast = none."),
                selectInput(ns("contrast3"), "Contrast enhancement", c("none", "histogram", "adaptive")),
                checkboxInput(ns("fixed_scale"), "Fixed intensity scale (optional)", value = FALSE),
                conditionalPanel(
                  condition = "input.fixed_scale",
                  ns = ns,
                  fluidRow(
                    column(6, numericInput(ns("zmin"), "Scale min", value = 0)),
                    column(6, numericInput(ns("zmax"), "Scale max", value = NA))
                  ),
                  actionButton(ns("set_zmax"), "Set max from current ion(s)"),
                  tags$small(class = "msie-hint",
                             "Same intensity = same colour in every image (custom ions and ion math). Per-image scaling and contrast enhancement are turned off while this is on; values outside min/max are shown at the ends of the scale.")
                ),
                checkboxInput(ns("outlier_check"), "Report outlier pixels (custom ions)", value = FALSE),
                selectInput(ns("smooth3"), "Smoothing options",
                            c("none", "mean", "gaussian", "bilateral", "adaptive", "diffusion", "guided")),
                fluidRow(
                  column(6, checkboxInput(ns("normalize3"), "Scale multiple images?", value = TRUE)),
                  column(6, checkboxInput(ns("colorkey3"), "Draw color key?", value = TRUE))
                ),
                checkboxInput(ns("dark_bg"), "Dark background?", value = FALSE),
                tags$hr(class = "msie-rule"),
                checkboxInput(ns("plot_pdata"), "Plot Phenotype data?", value = FALSE),
                uiOutput(ns("plotpdata"))
              )
            ),
            
            tabPanel(
              "Pixels",
              div(
                class = "msie-opt-body",
                checkboxInput(ns("apply_mask"), "Apply pData pixel mask?", value = FALSE),
                uiOutput(ns("mask_options")),
                tags$hr(class = "msie-rule"),
                checkboxInput(ns("expand_runs"), "Select individual runs for plotting only?", value = FALSE),
                uiOutput(ns("select_runs"))
              )
            ),
            
            tabPanel(
              "Size",
              div(
                class = "msie-opt-body",
                checkboxInput(ns("auto_aspect"), "Auto height (match image aspect)", value = TRUE),
                checkboxInput(ns("fit_width"), "Fit image to window width", value = FALSE),
                fluidRow(
                  column(6, numericInput(ns("width_im"), "Width (px)", value = 800, step = 50)),
                  column(6, numericInput(ns("height_im"), "Height (px)", value = 600, step = 50))
                ),
                tags$small(class = "msie-hint",
                           "Height is ignored while Auto height is ticked."),
                tags$hr(class = "msie-rule"),
                checkboxInput(ns("expand_fonts"), "Extended font options?", value = FALSE),
                uiOutput(ns("fonts")),
                tags$hr(class = "msie-rule"),
                textInput(ns("png_basename"), "PNG base name", value = "ionmap"),
                tags$small(class = "msie-hint",
                           "Saved as <base>_<metabolite>_<m/z>.png in the working directory."),
                checkboxInput(ns("png_transparent"), "Transparent background", value = FALSE),
                actionButton(ns("save_png"), "Save PNG to working directory"),
                actionButton(ns("save_png_batch"), "Save PNG for each selected ion"),
                tags$small(class = "msie-hint",
                           "Batch: one image per selected m/z row, current display settings, plain ion images (Ion math = none). Tick Fixed intensity scale for comparable colours."),
                downloadButton(ns("download_plot3_png"), "Download PNG"),
                tags$hr(class = "msie-rule"),
                downloadButton(ns("download_plot3_pdf"), "Download image (PDF)"),
                tags$hr(class = "msie-rule"),
                actionButton(ns("reset_display"), "Reset display settings"),
                tags$small(class = "msie-hint",
                           "Display and size settings are remembered between datasets and sessions.")
              )
            )
          )
        ),
        
        column(9, uiOutput(ns("plot.window")))
      )
    ),
    
    uiOutput(ns("spectrum"))
  )
  
}


# Callers invoke plot_card_server() from inside observers, once per dataset
# change. Originally every call built a brand-new module instance whose
# observers and data were never released (memory grew with every file
# restore, and old instances still reacted to buttons -- e.g. Save PNG wrote
# one copy per instance). Now the module is created once per id and later
# calls only hand it the new dataset.
plot_card_server <- function(id, overview_peaks_sel, spatialOnly=FALSE, allInputs=NULL) {
  session <- shiny::getDefaultReactiveDomain()
  if (is.null(session)) {
    return(plot_card_server_impl(id, shiny::reactiveVal(overview_peaks_sel),
                                 shiny::reactiveVal(1L), spatialOnly, allInputs))
  }
  reg <- get0("msie_plot_cards", envir = session$userData, inherits = FALSE)
  if (is.null(reg)) {
    reg <- new.env(parent = emptyenv())
    assign("msie_plot_cards", reg, envir = session$userData)
  }
  key <- session$ns(id)
  entry <- get0(key, envir = reg, inherits = FALSE)
  if (!is.null(entry)) {
    entry$data(overview_peaks_sel)
    entry$version(shiny::isolate(entry$version()) + 1L)
    return(invisible(NULL))
  }
  data_rv <- shiny::reactiveVal(overview_peaks_sel)
  version_rv <- shiny::reactiveVal(1L)
  out <- plot_card_server_impl(id, data_rv, version_rv, spatialOnly, allInputs)
  # Register only once the module exists, so a failed build is retried.
  assign(key, list(data = data_rv, version = version_rv), envir = reg)
  invisible(out)
}

plot_card_server_impl <- function(id, data_rv, data_version, spatialOnly=FALSE, allInputs=NULL) {
  
  moduleServer(id, function(input, output, session){
    
    #for dynamic UI
    ns = session$ns
    
    # Current dataset (with per-feature means, which the m/z table and the
    # spectrum limits use). Replaces the fixed object the module used to be
    # built around.
    overview_peaks_sel_r <- reactive({
      d <- data_rv()
      req(!is.null(d))
      if (is.null(fData(d)$mean)) {
        d2 <- try(Cardinal::summarizeFeatures(d, verbose = FALSE), silent = TRUE)
        if (inherits(d2, "try-error")) {
          showNotification("No data available, please check your parameters or dataset", type="error")
          message("plot_card_server summarizeFeatures failed: ", as.character(d2))
          req(FALSE)
        }
        d <- d2
      }
      d
    })
    plot3_last_png <- reactiveVal(NULL)
    # audit PLOTCARD-09: plot3_last_recorded / capture_plot3_recording removed -- the
    # helper was never called, so the recording was always NULL (export-only dead code).

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

    # Height:width of one image panel, from the pixel coordinates of the data.
    # Cardinal draws ion images with square pixels (asp = 1), but the colour key
    # is drawn at 2/3 of the *device* height regardless. So when the device is
    # much wider than the data, the image is letterboxed into a thin strip while
    # the key still spans the full height -- the "huge colour bar, lots of white
    # space" default. Sizing the device to the data's own aspect fixes it.
    image_panel_aspect <- function(obj) {
      co <- try(as.data.frame(Cardinal::coord(obj)), silent = TRUE)
      if (inherits(co, "try-error") || is.null(co) || ncol(co) < 2L) return(NULL)
      xv <- suppressWarnings(as.numeric(co[[1]]))
      yv <- suppressWarnings(as.numeric(co[[2]]))
      xv <- xv[is.finite(xv)]
      yv <- yv[is.finite(yv)]
      if (length(xv) == 0L || length(yv) == 0L) return(NULL)
      nx <- diff(range(xv)) + 1
      ny <- diff(range(yv)) + 1
      if (!is.finite(nx) || !is.finite(ny) || nx <= 0 || ny <= 0) return(NULL)
      ny / nx
    }

    # How many panels Cardinal will facet the image into (one per run, and one
    # per ion unless the ions are superposed).
    n_image_panels <- function(obj, n_ions = NULL) {
      n_runs <- try(length(unique(as.character(Cardinal::run(obj)))), silent = TRUE)
      if (inherits(n_runs, "try-error") || !is.finite(n_runs) || n_runs < 1L) n_runs <- 1L
      if (!is.null(n_ions)) {
        return(max(1L, as.integer(n_runs) * as.integer(n_ions)))
      }
      n_ions <- 1L
      if (identical(input$ion_viz3, "custom") && !isTRUE(input$superpose)) {
        sel <- try(mz_viz3(), silent = TRUE)
        if (!inherits(sel, "try-error") && length(sel) > 0L) n_ions <- length(sel)
      }
      max(1L, as.integer(n_runs) * as.integer(n_ions))
    }

    # Image dimensions. "Fit image to window width" draws at the width the browser
    # actually gives the image element (which grows when the options sidebar is
    # collapsed); "Auto height" then derives the height from the data's aspect.
    plot_dims <- function(obj = NULL, n_ions = NULL) {
      w0 <- suppressWarnings(as.numeric(input$width_im)[1])
      h0 <- suppressWarnings(as.numeric(input$height_im)[1])
      if (!is.finite(w0) || w0 <= 0) w0 <- 800
      if (!is.finite(h0) || h0 <= 0) h0 <- 600
      w <- w0
      if (isTRUE(input$fit_width)) {
        key <- paste0("output_", session$ns("plot3_pk"), "_width")
        cw <- suppressWarnings(as.numeric(session$clientData[[key]])[1])
        if (is.finite(cw) && cw >= 200) w <- cw
      }
      h <- h0 * (w / w0)
      if (isTRUE(input$auto_aspect) && !is.null(obj)) {
        asp <- image_panel_aspect(obj)
        if (!is.null(asp)) {
          n <- n_image_panels(obj, n_ions)
          nc <- max(1L, as.integer(ceiling(sqrt(n))))
          nr <- max(1L, as.integer(ceiling(n / nc)))
          # ~28% of the device width goes to axes, margins and the colour key.
          h_try <- (w * 0.72) * asp * (nr / nc) + 70
          if (is.finite(h_try) && h_try > 0) h <- h_try
        }
      }
      w <- min(max(w, 200), 8000)
      h <- min(max(h, 200), 8000)
      list(w = as.integer(round(w)), h = as.integer(round(h)))
    }

    # Every input the image render reads. In auto mode the render depends on all
    # of them (as it always did); otherwise only on the Update plot button.
    PLOT3_INPUT_IDS <- c(
      "ion_viz3", "tbl_rows_selected", "display_mode", "superpose", "plusminus_viz3",
      "contrast3", "color3", "smooth3", "normalize3", "colorkey3", "dark_bg",
      "apply_mask", "mask_variable", "mask_range", "mask_categories", "invert_mask",
      "show_mask_info", "width_im", "height_im", "fit_width", "auto_aspect",
      "outlier_check", "png_transparent", "clip_outliers", "clip_pct",
      "floor_clip", "floor_pct", "fixed_scale", "zmin", "zmax",
      "plot_pdata", "pdata_var_plot", "expand_fonts", "axis_size", "label_size",
      "title_size", "subtitle_size", "expand_runs", "select_runs"
    )

    plot3_trigger <- reactive({
      # A new dataset always redraws (as when each dataset got its own module).
      dv <- data_version()
      if (isTRUE(input$auto_update) || is.null(input$auto_update)) {
        # With 'Fit image to window width', a change of the image area's width
        # (window resize, options panel collapsed) also redraws.
        cw <- if (isTRUE(input$fit_width)) {
          session$clientData[[paste0("output_", session$ns("plot3_pk"), "_width")]]
        } else NULL
        c(list(data_version = dv, fit_width_px = cw), lapply(PLOT3_INPUT_IDS, function(nm) input[[nm]]))
      } else {
        list(data_version = dv, update = input$update_plot)
      }
    })

    # --- Remembered display settings -------------------------------------------
    # Colours, contrast, clipping, fixed scale and size are saved (per plot
    # panel, plus a shared "last used" copy for panels with nothing saved yet)
    # in the user config folder and restored when the panel is created.
    PERSIST_DEFAULTS <- list(
      color3 = "Plasma", contrast3 = "none", smooth3 = "none",
      normalize3 = TRUE, colorkey3 = TRUE, dark_bg = FALSE,
      clip_outliers = FALSE, clip_pct = 99.9, floor_clip = FALSE, floor_pct = 1,
      fixed_scale = FALSE, zmin = 0, zmax = NA_real_, outlier_check = FALSE,
      auto_aspect = TRUE, fit_width = FALSE, width_im = 800, height_im = 600,
      png_basename = "ionmap", png_transparent = FALSE, auto_update = TRUE
    )
    PERSIST_SELECT <- c("color3", "contrast3", "smooth3")
    PERSIST_CHECK <- c("normalize3", "colorkey3", "dark_bg", "clip_outliers", "floor_clip",
                       "fixed_scale", "outlier_check", "auto_aspect", "fit_width",
                       "png_transparent", "auto_update")
    PERSIST_TEXT <- c("png_basename")
    settings_file <- tryCatch(
      file.path(tools::R_user_dir("MSI.EAGLE", which = "config"), "plot_settings.rds"),
      error = function(e) NULL)
    settings_key <- session$ns("settings")
    read_settings_store <- function() {
      if (is.null(settings_file) || !file.exists(settings_file)) return(list())
      st <- tryCatch(readRDS(settings_file), error = function(e) list())
      if (is.list(st)) st else list()
    }
    write_settings_store <- function(st) {
      if (is.null(settings_file)) return(invisible(FALSE))
      tryCatch({
        dir.create(dirname(settings_file), recursive = TRUE, showWarnings = FALSE)
        saveRDS(st, settings_file)
        TRUE
      }, error = function(e) FALSE)
    }
    apply_display_settings <- function(vals) {
      for (nm in intersect(names(vals), names(PERSIST_DEFAULTS))) {
        v <- vals[[nm]]
        if (length(v) != 1L) next
        if (nm %in% PERSIST_SELECT) {
          updateSelectInput(session, nm, selected = as.character(v))
        } else if (nm %in% PERSIST_CHECK) {
          updateCheckboxInput(session, nm, value = isTRUE(as.logical(v)))
        } else if (nm %in% PERSIST_TEXT) {
          updateTextInput(session, nm, value = as.character(v))
        } else {
          updateNumericInput(session, nm, value = suppressWarnings(as.numeric(v)))
        }
      }
    }
    local({
      st <- read_settings_store()
      saved <- st[[settings_key]]
      if (is.null(saved)) saved <- st[["_last"]]
      if (is.list(saved) && length(saved)) apply_display_settings(saved)
    })
    settings_first_run <- TRUE
    observe({
      vals <- lapply(names(PERSIST_DEFAULTS), function(nm) input[[nm]])
      names(vals) <- names(PERSIST_DEFAULTS)
      # The first run sees the page defaults, before the restored values
      # reach the browser -- never save those.
      if (settings_first_run) {
        settings_first_run <<- FALSE
        return(invisible(NULL))
      }
      vals <- vals[!vapply(vals, is.null, logical(1))]
      isolate({
        st <- read_settings_store()
        st[[settings_key]] <- vals
        st[["_last"]] <- vals
        write_settings_store(st)
      })
    })
    observeEvent(input$reset_display, {
      apply_display_settings(PERSIST_DEFAULTS)
      showNotification("Display settings reset to defaults.", type = "message", duration = 5)
    }, ignoreInit = TRUE)

    # --- PNG export ------------------------------------------------------------
    plot3_last_label <- reactiveVal("image")
    outlier_report <- reactiveVal(NULL)

    # A file-name label for what is currently displayed: target/ID name plus m/z
    # for each selected ion (windowed picks can have several sub-peaks per
    # target, so the m/z is always kept).
    plot3_label <- function(rows = NULL) {
      mode <- input$ion_viz3
      batch_rows <- !is.null(rows)
      if (!batch_rows) {
        if (isTRUE(input$plot_pdata) && !is.null(input$pdata_var_plot)) {
          return(paste0("pdata_", input$pdata_var_plot))
        }
        if (identical(mode, "viz_all")) return("TIC")
      }
      fd <- try(as.data.frame(fData(overview_peaks_sel_r())), silent = TRUE)
      if (inherits(fd, "try-error") || is.null(fd) || nrow(fd) == 0L) return("image")
      if (!batch_rows) rows <- if (identical(mode, "custom")) input$tbl_rows_selected else 1L
      rows <- suppressWarnings(as.integer(rows))
      rows <- rows[is.finite(rows) & rows >= 1L & rows <= nrow(fd)]
      if (length(rows) == 0L) return("image")
      labs <- vapply(rows, function(r) {
        mzv <- sprintf("%.4f", as.numeric(fd[["mz"]][r]))
        nm <- ""
        for (col in c("target", "ID")) {
          if (col %in% names(fd)) {
            v <- as.character(fd[[col]][r])
            if (length(v) == 1L && !is.na(v) && nzchar(trimws(v))) { nm <- v; break }
          }
        }
        if (nzchar(nm)) paste0(nm, "_", mzv) else paste0("mz", mzv)
      }, character(1))
      if (!batch_rows && identical(mode, "custom") && !is.null(input$display_mode) &&
          !identical(input$display_mode, "none")) {
        labs <- c(input$display_mode, labs)
      }
      paste(labs, collapse = "__")
    }

    # Cap values above a percentile (computed over acquired pixels). Used to
    # stop a few very bright pixels setting the top of the colour scale; any
    # contrast enhancement is then applied to the capped image.
    clip_pct_value <- function() {
      pct <- suppressWarnings(as.numeric(input$clip_pct)[1])
      if (length(pct) != 1L || !is.finite(pct)) pct <- 99.9
      min(max(pct, 50), 100)
    }
    cap_values <- function(v, pct) {
      v <- as.numeric(v)
      ok <- is.finite(v)
      if (!any(ok) || pct >= 100) return(v)
      cap <- as.numeric(stats::quantile(v[ok], pct / 100, names = FALSE))
      v[ok] <- pmin(v[ok], cap)
      v
    }

    floor_pct_value <- function() {
      pct <- suppressWarnings(as.numeric(input$floor_pct)[1])
      if (length(pct) != 1L || !is.finite(pct)) pct <- 1
      min(max(pct, 0), 50)
    }
    floor_values <- function(v, pct) {
      v <- as.numeric(v)
      ok <- is.finite(v)
      if (!any(ok) || pct <= 0) return(v)
      lo <- as.numeric(stats::quantile(v[ok], pct / 100, names = FALSE))
      v[ok] <- pmax(v[ok], lo)
      v
    }
    # Clip top / raise floor, as ticked. Applied to each ion image before any
    # contrast enhancement or fixed scale.
    transform_ion_values <- function(v) {
      if (isTRUE(input$clip_outliers)) v <- cap_values(v, clip_pct_value())
      if (isTRUE(input$floor_clip)) v <- floor_values(v, floor_pct_value())
      v
    }
    # Fixed colour-scale limits. Blank max falls back to the current image(s)
    # maximum, with a reminder that the scale is then not fixed.
    fixed_zlim <- function(all_vals) {
      zmin <- suppressWarnings(as.numeric(input$zmin)[1])
      zmax <- suppressWarnings(as.numeric(input$zmax)[1])
      if (length(zmin) != 1L || !is.finite(zmin)) zmin <- 0
      if (length(zmax) != 1L || !is.finite(zmax)) {
        fin <- all_vals[is.finite(all_vals)]
        zmax <- if (length(fin)) max(fin) else zmin + 1
        if (!plot3_batch_mode) {
          showNotification("Fixed scale: no max set, so this image uses its own maximum. Enter a max (or press 'Set max from current ion(s)') to compare images.",
                           type = "warning", duration = 8)
        }
      }
      if (zmax <= zmin) zmax <- zmin + abs(zmin) * 1e-6 + 1e-12
      c(zmin, zmax)
    }

    safe_file_part <- function(x) {
      x <- gsub("[^A-Za-z0-9._+-]+", "_", as.character(x))
      x <- gsub("_+", "_", x)
      gsub("^_+|_+$", "", x)
    }

    png_filename <- function(label = NULL) {
      base <- input$png_basename
      if (is.null(base) || !nzchar(trimws(base))) base <- "ionmap"
      if (is.null(label)) label <- isolate(plot3_last_label())
      paste0(safe_file_part(base), "_", safe_file_part(label), ".png")
    }

    # Never overwrite: name.png, name_2.png, name_3.png, ...
    unique_png_dest <- function(fn) {
      dest <- file.path(getwd(), fn)
      if (file.exists(dest)) {
        stem <- sub("\\.png$", "", fn)
        k <- 2L
        repeat {
          dest <- file.path(getwd(), paste0(stem, "_", k, ".png"))
          if (!file.exists(dest)) break
          k <- k + 1L
        }
      }
      dest
    }

    observeEvent(input$save_png, {
      src <- isolate(plot3_last_png())
      if (is.null(src) || !file.exists(src)) {
        showNotification("Render an image first, then save it.", type = "warning", duration = 6)
        return(invisible(NULL))
      }
      dest <- unique_png_dest(png_filename())
      ok <- isTRUE(file.copy(src, dest, overwrite = FALSE))
      if (ok) {
        message(sprintf("[plot_card:%s] saved %s", id, dest))
        showNotification(paste("Saved", dest), type = "message", duration = 8)
      } else {
        showNotification(paste("Could not save", dest), type = "error", duration = 10)
      }
    }, ignoreInit = TRUE)

    output$download_plot3_png <- downloadHandler(
      filename = function() png_filename(),
      content = function(file) {
        src <- isolate(plot3_last_png())
        if (is.null(src) || !file.exists(src)) {
          showNotification("Render an image first, then download it.", type = "warning", duration = 6)
          stop("No rendered image available.")
        }
        file.copy(src, file, overwrite = TRUE)
      }
    )

    output$outlier_report <- renderText({
      txt <- outlier_report()
      if (is.null(txt)) "" else txt
    })

    make_plot3_image_result <- function(src_png, alt_text, dims = NULL) {
      if (!plot3_batch_mode) {
        # (recordPlot() used to run here, after the PNG device was closed, so
        # it recorded nothing -- or, under RStudio, whatever the RStudio
        # device showed. The PDF download now always embeds the cached PNG.)
        cache_plot3_png(src_png)
      }
      if (is.null(dims)) dims <- plot_dims()
      list(
        src = src_png,
        contentType = 'image/png',
        width = dims$w,
        height = dims$h,
        alt = alt_text
      )
    }
    
    output$download_plot3_pdf <- downloadHandler(
      filename = function() {
        paste0("plot_card_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
      },
      content = function(file) {
        # audit PLOTCARD-09: the page used to be Width x Height (px/72) while the PNG was
        # drawn at plot_dims() (auto height / fit width), so the raster was stretched. The
        # page is now the cached PNG's own size (px/72), so the PDF holds the same pixels
        # at their true aspect ratio; interpolate = FALSE keeps ion pixels sharp. The dead
        # recordPlot replay branch (recording was never captured) is removed. Export only:
        # no analysis value or on-screen image changes.
        src_png <- isolate(plot3_last_png())
        if (is.null(src_png) || !file.exists(src_png)) {
          showNotification("No image available yet. Render a plot first, then download PDF.", type = "warning", duration = 5)
          stop("No rendered image available for PDF export.")
        }

        img <- try(png::readPNG(src_png), silent = TRUE)
        if (inherits(img, "try-error")) {
          stop("Could not read cached image for PDF export.")
        }

        grDevices::pdf(file = file, width = dim(img)[2] / 72, height = dim(img)[1] / 72, onefile = TRUE)
        on.exit(grDevices::dev.off(), add = TRUE)

        grid::grid.newpage()
        grid::grid.raster(img, x = 0.5, y = 0.5, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"), interpolate = FALSE)
      }
    )
    
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

    discrete_fallback_colors <- function(n) {
      n <- suppressWarnings(as.integer(n))
      if (!is.finite(n) || n <= 0) return(character(0))
      # Prefer Alphabet family for categorical pData defaults.
      # On some R builds, Alphabet is available via palette.colors() but not hcl.colors().
      pal_try <- function(name) {
        tryCatch(grDevices::palette.colors(n, palette = name), error = function(e) character(0))
      }
      cols <- pal_try("Alphabet")
      if (length(cols) != n) {
        cols <- pal_try("Alphabet2")
      }
      if (length(cols) != n) {
        cols <- tryCatch(grDevices::hcl.colors(n, palette = "Dark 3"), error = function(e) character(0))
      }
      if (length(cols) != n) {
        cols <- grDevices::rainbow(n)
      }
      as.character(cols)
    }
    
    # NEW: Add UI for pixel masking options
    observe({
      output$mask_options <- renderUI({
        req(overview_peaks_sel_r())
        
        if(input$apply_mask) {
          # Get available pData variables (exclude coordinate columns)
          pdata_vars <- colnames(as.data.frame(pData(overview_peaks_sel_r())))
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
        req(overview_peaks_sel_r())
        
        # Get unique values for the selected pData variable
        pdata_df <- as.data.frame(pData(overview_peaks_sel_r()))
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
    
    # audit PLOTCARD-11: the pixel-selection logic is factored out of
    # apply_pixel_mask() unchanged (same statements, same order), so the mask
    # notification observer below can count pixels without subsetting.
    pixel_mask_vector <- function(data) {
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
      pixel_mask
    }

    # NEW: Function to apply pixel mask
    # audit PLOTCARD-10/11: `notify = FALSE` (used by plot_obj_r) and batch export skip
    # only the "Pixel mask applied" notification; the pixels selected and the
    # returned object are exactly as before.
    apply_pixel_mask <- function(data, notify = TRUE) {
      if(!input$apply_mask) {
        return(data)
      }

      pixel_mask <- pixel_mask_vector(data)

      # Show mask information if requested
      if(notify && !plot3_batch_mode && input$show_mask_info) {
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

    # audit PLOTCARD-11: the object the ion image is drawn from (mask, then run
    # subset) is computed once here and shared by the image render, batch export
    # and the slot-height observer. It is the same expression plot3_draw used to
    # evaluate inline (apply_pixel_mask + subsetPixels(run %in% select_runs)) on
    # the same inputs, and both steps are deterministic, so the drawn object is
    # identical; only the reserved slot height now uses the masked, run-subset
    # extent that the image actually has.
    plot_obj_r <- reactive({
      d <- apply_pixel_mask(overview_peaks_sel_r(), notify = FALSE)
      if(!is.null(input$select_runs)) {
        d <- subsetPixels(d, run %in% input$select_runs)
      }
      d
    })

    # audit PLOTCARD-10/11: the "Pixel mask applied" message is shown once per
    # change of the mask settings (or dataset) instead of once per redraw / per
    # batch image. Notification only; nothing drawn or computed changes.
    observe({
      if (!isTRUE(input$apply_mask) || !isTRUE(input$show_mask_info)) return(invisible(NULL))
      pm <- try(pixel_mask_vector(overview_peaks_sel_r()), silent = TRUE)
      if (inherits(pm, "try-error") || length(pm) == 0L) return(invisible(NULL))
      showNotification(
        paste("Pixel mask applied:", sum(pm, na.rm = TRUE), "of", length(pm), "pixels selected"),
        type = "message", duration = 3
      )
    })

    # audit PLOTCARD-07: the plot window is now built once (it no longer depends on
    # input$ion_viz3). The Custom-only parts -- the slot hint, the outlier/report
    # text and the m/z table -- sit in conditionalPanels, so switching the ion mode
    # hides them instead of destroying the DT (which cleared the row selection).
    # Same elements, ids and classes as the old Custom layout; the image output is
    # the same single imageOutput in every mode. Layout only.
    observe({
      
      output$plot.window <- renderUI({
        
        fluidRow(
          # height = "auto" lets the container grow with the rendered PNG (the
          # default 400px made a taller image overlap the list below), and the
          # slot div keeps the space reserved so the list does not jump when an
          # image appears or is cleared.
          column(12,
                 div(class = "msie-image-slot", id = ns("image_slot"),
                     conditionalPanel(
                       condition = "input.ion_viz3 == 'custom'",
                       ns = ns,
                       div(class = "msie-slot-hint",
                           "Select one or more m/z below -- arrow keys step through ions, Shift+arrow adds")
                     ),
                     imageOutput(ns("plot3_pk"), height = "auto"))),
          conditionalPanel(
            condition = "input.ion_viz3 == 'custom'",
            ns = ns,
            column(12, tags$pre(class = "msie-hint", style = "white-space: pre-wrap; border: none; background: none; padding: 0;",
                                textOutput(ns("outlier_report"), inline = TRUE))),
            fluidRow(
              column(12, tags$hr(class = "msie-rule")),
              column(12, div(class = "msie-subtitle", "m/z list -- click rows to choose ions")),
              column(12, div(class = "msie-ion-table", `data-ns` = ns(""),
                             DT::dataTableOutput(ns('tbl'))))
            )
          )
        )
        
      })
      
    })
    
    
    
    
    # audit PLOTCARD-07: the m/z table is (re)assigned only when the dataset
    # changes -- the same table, rounding and DT call as before, formerly inside
    # the mz_viz3a renderUI, which also re-ran on every ion-mode switch and so
    # cleared the row selection. A new dataset still re-renders (and clears) it.
    # The superpose / Ion math / +/- m/z controls moved to the static UI.
    observe({
      req(overview_peaks_sel_r())
      
      tbl<-as.data.frame(fData(overview_peaks_sel_r()))
      if(!is.null(tbl$freq)) {
        tbl$freq<-round(tbl$freq, 2)
      }
      if(!is.null(tbl$mean)) {
        tbl$mean<-round(tbl$mean, 1)
      }
      for (column in intersect(c("mz", "mz_original", "correction_ppm"), names(tbl))) {
        tbl[[column]] <- round(tbl[[column]], 4)
      }
      
      
      # (No longer resets the image width/height here: sizes now persist
      # across datasets and sessions.)
      
      #add something here at some point to only select mz, freq, count, and mean columns
      
      
      # Windowed peak picking adds ppm_from_target. What matters there is how
      # far a sub-peak is from the target, not its sign, so that column sorts
      # by |ppm| (via a hidden helper column) while still showing the sign.
      tbl_opts <- list()
      if ("ppm_from_target" %in% names(tbl)) {
        tbl$abs_ppm <- abs(as.numeric(tbl$ppm_from_target))
        ppm_col <- which(names(tbl) == "ppm_from_target")
        abs_col <- which(names(tbl) == "abs_ppm")
        # DataTables indices are 0-based with the row-name column at 0, so a
        # data column's R position equals its DataTables index.
        tbl_opts <- list(columnDefs = list(
          list(orderData = abs_col, targets = ppm_col),
          list(visible = FALSE, targets = abs_col)
        ))
      }

      output$tbl <-DT::renderDataTable({
        tbl
        
      }, selection = "multiple", options = tbl_opts)
    })
    
    
    
    
    # Tell the browser how tall the next image will be, so the slot reserves
    # exactly that much and the m/z list never shifts. This runs whenever a size
    # control or the ion selection changes -- i.e. before the image is redrawn.
    observe({
      # audit PLOTCARD-11: size from the masked, run-subset object the image is
      # drawn from (was the unmasked, all-runs data). Layout only.
      dims <- try(plot_dims(plot_obj_r()), silent = TRUE)
      if (inherits(dims, "try-error") || is.null(dims$h)) return(invisible(NULL))
      h <- min(max(as.integer(dims$h) + 16L, 160L), 8016L)
      try(session$sendCustomMessage("msieSlotHeight",
                                    list(id = session$ns("image_slot"), h = h)),
          silent = TRUE)
    })

    # Arrow-key navigation of the m/z table (see the keydown handler in
    # msi_ui_assets). Down/Up move to the neighbouring ion; Shift+Down/Up extend
    # the selection. Selection is driven through the DT proxy, so everything
    # downstream sees an ordinary tbl_rows_selected change.
    observeEvent(input$mz_key_nav, {
      info <- input$mz_key_nav
      if (is.null(info) || is.null(info$dir)) return(invisible(NULL))
      if (!identical(input$ion_viz3, "custom")) return(invisible(NULL))
      n_rows <- try(nrow(as.data.frame(fData(overview_peaks_sel_r()))), silent = TRUE)
      if (inherits(n_rows, "try-error") || length(n_rows) != 1L ||
          is.na(n_rows) || n_rows < 1L) {
        return(invisible(NULL))
      }
      step <- if (identical(info$dir, "down")) 1L else -1L
      # Step through the table as it is currently DISPLAYED -- tbl_rows_all is
      # the data-row order after the active sort and search filter. Using raw
      # data indices would walk in m/z order even when the table is sorted by
      # intensity, or jump outside a search filter.
      order_idx <- suppressWarnings(as.integer(input$tbl_rows_all))
      order_idx <- order_idx[is.finite(order_idx) & order_idx >= 1L & order_idx <= n_rows]
      if (length(order_idx) == 0L) order_idx <- seq_len(as.integer(n_rows))

      cur <- input$tbl_rows_selected
      cur <- suppressWarnings(as.integer(cur))
      cur <- cur[is.finite(cur) & cur >= 1L & cur <= n_rows]

      cur_pos <- match(cur, order_idx)
      cur_pos <- cur_pos[!is.na(cur_pos)]
      if (length(cur_pos) == 0L) {
        new_pos <- if (step > 0L) 1L else length(order_idx)
      } else {
        anchor_pos <- if (step > 0L) max(cur_pos) else min(cur_pos)
        new_pos <- anchor_pos + step
      }
      if (new_pos < 1L || new_pos > length(order_idx)) return(invisible(NULL))
      target <- order_idx[new_pos]
      new_sel <- if (isTRUE(info$add) && length(cur) > 0L) {
        sort(unique(c(cur, target)))
      } else {
        target
      }
      proxy <- try(DT::dataTableProxy("tbl", session = session), silent = TRUE)
      if (inherits(proxy, "try-error")) return(invisible(NULL))
      # Keep the target row on screen when stepping past a page boundary.
      page_rows <- suppressWarnings(as.integer(input$tbl_rows_current))
      page_rows <- page_rows[is.finite(page_rows)]
      if (length(page_rows) > 0L && !(target %in% page_rows)) {
        try(DT::selectPage(proxy, ceiling(new_pos / length(page_rows))), silent = TRUE)
      }
      try(DT::selectRows(proxy, new_sel), silent = TRUE)
    }, ignoreInit = TRUE)

    # Ion images (summed within +/- tolerance, clip/floor applied) for a set of
    # m/z values on the pixels currently displayed. Used by "Set max".
    current_ion_values <- function(ions) {
      dat <- apply_pixel_mask(overview_peaks_sel_r())
      if (!is.null(input$select_runs)) {
        dat <- subsetPixels(dat, run %in% input$select_runs)
      }
      tol <- suppressWarnings(as.numeric(input$plusminus_viz3)[1])
      if (length(tol) != 1L || !is.finite(tol) || tol < 0) tol <- 0.01
      tol <- round(tol, 3)
      mz_axis <- suppressWarnings(as.numeric(Cardinal::mz(dat)))
      lapply(ions, function(m) {
        idx <- which(is.finite(mz_axis) & abs(mz_axis - m) <= tol)
        if (length(idx) == 0L) return(NULL)
        mat <- try(as.matrix(Cardinal::spectra(dat[idx, ])), silent = TRUE)
        if (inherits(mat, "try-error")) return(NULL)
        if (nrow(mat) != length(idx) && ncol(mat) == length(idx)) mat <- t(mat)
        transform_ion_values(colSums(mat, na.rm = TRUE))
      })
    }

    observeEvent(input$set_zmax, {
      if (!identical(input$ion_viz3, "custom")) {
        showNotification("Choose 'Custom single / multiple' ions first.", type = "warning", duration = 6)
        return(invisible(NULL))
      }
      ions <- tryCatch(as.numeric(mz_viz3()), error = function(e) numeric(0))
      if (length(ions) == 0L) {
        showNotification("Select one or more ions in the m/z table first.", type = "warning", duration = 6)
        return(invisible(NULL))
      }
      vals <- tryCatch(unlist(current_ion_values(ions), use.names = FALSE), error = function(e) NULL)
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0L) {
        showNotification("Could not compute the selected ion image(s).", type = "error", duration = 6)
        return(invisible(NULL))
      }
      zmax <- signif(max(vals), 4)
      updateNumericInput(session, "zmax", value = zmax)
      msg <- paste0("Scale max set to ", format(zmax), " (maximum of the selected ion image(s)",
                    if (isTRUE(input$clip_outliers)) " after clipping" else "", ").")
      if (!is.null(input$display_mode) && !identical(input$display_mode, "none")) {
        msg <- paste(msg, "Note: this uses the plain ion images, not the ion-math result.")
      }
      showNotification(msg, type = "message", duration = 8)
    }, ignoreInit = TRUE)

    mz_viz3 <- reactive({
      req(overview_peaks_sel_r())
      req(input$tbl_rows_selected)
      tbl<-as.data.frame(fData(overview_peaks_sel_r()))
      rows <- suppressWarnings(as.integer(input$tbl_rows_selected))
      rows <- rows[is.finite(rows) & rows >= 1L & rows <= nrow(tbl)]
      req(length(rows) > 0L)
      return(tbl[rows, "mz"])
      
    })
    
    
    
    observe({
      
      output$select_runs <- renderUI ({
        req(overview_peaks_sel_r())
        
        if(input$expand_runs) {
          
          run_list<- unique(run((overview_peaks_sel_r())))
          
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
                        choices = colnames(as.data.frame(pData(overview_peaks_sel_r())))[-c(1:2)]
            )
          )
        }
        
        
      })
    })
    
    # The image options now live in tab panels, and Shiny suspends outputs that
    # are not visible -- which would mean the mask / font / run / phenotype
    # controls did not exist until their tab was opened. Keep them live so the
    # inputs behave exactly as they did when every control was in one column.
    # (audit PLOTCARD-07: "mz_viz3a" dropped from this list -- that output no longer exists.)
    for (.nm in c("mask_options", "plotpdata", "fonts", "select_runs", "plot.window")) {
      try(outputOptions(output, .nm, suspendWhenHidden = FALSE), silent = TRUE)
    }

    #add spectrum or not
    output$spectrum<-renderUI({

      # audit PLOTCARD-06: removed a dead `if (input$ion_viz3 == "custom")
      # DT::dataTableOutput(...)` whose value was discarded (renderUI returns only
      # the last expression). The rendered UI is identical; the panel just no longer
      # rebuilds -- and resets its spectrum controls -- on every ion-mode switch.

      if(spatialOnly==FALSE){
        div(
          class = "msie-panel",
          div(class = "msie-panel-title", "Mass spectrum"),
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
            column(10, imageOutput(ns("plot4"), height = "auto"), style = "margin-bottom: 0px; padding-bottom: 0px;")
          ),
          fluidRow(
            column(3, checkboxInput(ns("calc_ppm"), label = "Show ppm error?", width = '100%')),
            column(3, checkboxInput(ns("show_int"), label = "Show intensity?", width = '100%')),
            column(3, numericInput(ns("show_mz"), label = "# mz values to show?", width = '100%', value=0))
          ),
          checkboxInput(ns("spectrum_expand_fonts"), "Extended font options?", value=FALSE),
          uiOutput(ns("spectrum_fonts"))
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
    
    # Draws the ion image into `outfile`. Used by the on-screen render and by
    # batch export (which passes one ion at a time via ion_override).
    plot3_batch_mode <- FALSE
    plot3_draw <- function(outfile = tempfile(fileext = '.png'), ion_override = NULL) {
        
        #req(overview_peaks_sel_r())
        
        .plot3_t0 <- proc.time()[["elapsed"]]
        # Only report redraws slow enough to be worth knowing about.
        on.exit({
          .el <- proc.time()[["elapsed"]] - .plot3_t0
          if (is.finite(.el) && .el >= 0.5) {
            message(sprintf("[plot_card:%s] image redraw total %.2fs", id, .el))
          }
        }, add = TRUE)
        
        # MODIFIED: Apply pixel mask to the data before plotting
        # audit PLOTCARD-10/11: same mask + run subset, now taken from the shared
        # plot_obj_r() (identical expression and inputs), so it is computed once per
        # mask/run change instead of once per redraw or per batch image.
        overview_peaks_sel_masked <- plot_obj_r()
        
        .dims <- plot_dims(overview_peaks_sel_masked,
                           n_ions = if (!is.null(ion_override)) length(ion_override) else NULL)
        png(outfile, width = .dims$w, height = .dims$h,
            bg = if (isTRUE(input$png_transparent)) "transparent" else "white")
        # Several early exits below return() without dev.off(); make sure this
        # device is always closed so devices do not pile up (this is what the
        # old graphics.off() at module start was papering over).
        .plot3_dev <- grDevices::dev.cur()
        on.exit({
          if (.plot3_dev %in% grDevices::dev.list()) try(grDevices::dev.off(.plot3_dev), silent = TRUE)
        }, add = TRUE)
        if (!plot3_batch_mode) {
          # Remember what this image shows, for the PNG file name.
          try(plot3_last_label(plot3_label()), silent = TRUE)
          if (!isTRUE(input$outlier_check) || !identical(input$ion_viz3, "custom")) {
            outlier_report(NULL)
          }
        }

        if (is.null(overview_peaks_sel_masked) || ncol(overview_peaks_sel_masked) == 0) {
          graphics::plot.new()
          graphics::text(0.5, 0.5, "No pixels available after mask/run filtering.")
          dev.off()
          return(make_plot3_image_result(outfile, "No pixels available", .dims))
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
            return(make_plot3_image_result(outfile, "Invalid pData field", .dims))
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
            return(make_plot3_image_result(outfile, "pData mismatch", .dims))
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
              fallback_cols <- discrete_fallback_colors(sum(!valid_col))
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
          
          if (!is.null(ion_override)) {
            ion <- as.numeric(ion_override)
          } else if(is.null(mz_viz3())){
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
            uniq_idx <- sort(unique(idx_map))
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
              if (!contrast_mode %in% c("none", "histogram", "adaptive", "adjust")) {
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
                
                if (mode %in% c("histogram", "adjust")) {
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
            
            overview_peaks_sel_masked$xic <- transform_ion_values(overview_peaks_sel_masked$xic)
            xic_zlim <- NULL
            if (isTRUE(input$fixed_scale)) {
              xic_zlim <- fixed_zlim(overview_peaks_sel_masked$xic)
              overview_peaks_sel_masked$xic <- pmin(pmax(overview_peaks_sel_masked$xic, xic_zlim[1]), xic_zlim[2])
            }
            
            plusminus=input$plusminus_viz3
            
            smoothing_option <- if (input$smooth3 != "none")  input$smooth3 else NULL
            enhance_option <- if (input$contrast3 != "none")  input$contrast3 else NULL
            
            
            arg_list<-list(overview_peaks_sel_masked,
                           'xic',
                           tolerance=round(plusminus,3), 
                           units='mz',
                           col=cpal(input$color3),
                           enhance=if (is.null(xic_zlim)) enhance_option else NULL,
                           smooth=smoothing_option,
                           scale=if (is.null(xic_zlim)) input$normalize3 else FALSE,
                           key=(input$colorkey3),
                           grid=FALSE)
            if (!is.null(xic_zlim)) arg_list$zlim <- xic_zlim
            
            if(input$dark_bg) {
              arg_list$style <- "dark"
            }
            
            print(matter::as_facets(do.call(Cardinal::image, arg_list), labels=label_txt))
            
            vizi_par(vp_orig)
            }
            
          } else {
            
            mz_set=ion
            
            # Optional outlier report: a single very bright pixel sets the top
            # of the colour scale and squeezes everything else. Compare the
            # maximum with the 99.9th percentile for each selected ion.
            if (isTRUE(input$outlier_check) && !plot3_batch_mode) {
              ov <- try(extract_exact_ions_matrix_once(ion), silent = TRUE)
              if (!inherits(ov, "try-error") && !is.null(ov)) {
                lines_out <- vapply(seq_len(nrow(ov$mat)), function(r) {
                  v <- as.numeric(ov$mat[r, ])
                  v <- v[is.finite(v)]
                  if (length(v) < 10L) return(sprintf("m/z %.4f: too few pixels", ov$ions[r]))
                  q999 <- as.numeric(stats::quantile(v, 0.999, names = FALSE))
                  q99 <- as.numeric(stats::quantile(v, 0.99, names = FALSE))
                  mx <- max(v)
                  n_hot <- if (q999 > 0) sum(v > 2 * q999) else 0L
                  ratio <- if (q999 > 0) mx / q999 else NA_real_
                  flag <- if (is.finite(ratio) && ratio > 2) "  (outliers likely: tick Clip bright outliers)" else ""
                  sprintf("m/z %.4f: max %.4g | 99.9th pct %.4g | 99th pct %.4g | max = %.1fx the 99.9th pct | %d pixel(s) > 2x the 99.9th pct%s",
                          ov$ions[r], mx, q999, q99, ratio, as.integer(n_hot), flag)
                }, character(1))
                outlier_report(paste(lines_out, collapse = "\n"))
              } else {
                outlier_report("Outlier check: could not extract the selected ion(s).")
              }
            }
            
            mz_range=range(mz(overview_peaks_sel_masked))
            
            plusminus=input$plusminus_viz3
            
            smoothing_option <- if (input$smooth3 != "none")  input$smooth3 else NULL
            enhance_option <- if (input$contrast3 != "none")  input$contrast3 else NULL

            # audit PLOTCARD-V02: note under the image when a panel sums several
            # features. It only counts features of mz(x) within +/- round(tol, 3) --
            # the same |m/z - target| <= tol rule as Cardinal's kdsearch(tol.ref =
            # "abs") (matter search.h:598-604) used for the image -- and appends text
            # to the report line. The image, the tolerance and the summation are
            # untouched; any failure here is ignored.
            if (!plot3_batch_mode) {
              try({
                tol_hint <- round(as.numeric(plusminus)[1], 3)
                if (length(tol_hint) == 1L && is.finite(tol_hint)) {
                  n_in <- vapply(mz_set, function(m) {
                    as.integer(sum(is.finite(mz_axis_all) & abs(mz_axis_all - m) <= tol_hint))
                  }, integer(1))
                  multi <- which(n_in > 1L)
                  if (length(multi) > 0L) {
                    hint_lines <- sprintf("m/z %.4f: image is the sum of %d features within +/- %s",
                                          mz_set[multi], n_in[multi], format(tol_hint))
                    prev_txt <- isolate(outlier_report())
                    outlier_report(paste(c(prev_txt, hint_lines), collapse = "\n"))
                  }
                }
              }, silent = TRUE)
            }

            if (isTRUE(input$clip_outliers) || isTRUE(input$floor_clip) || isTRUE(input$fixed_scale)) {
              # Same ion image Cardinal would draw -- every feature within
              # +/- tolerance of each m/z, summed per pixel (Cardinal's
              # reduce = "+") -- built here so it can be clipped and/or put on
              # a fixed scale first; any contrast enhancement then works on
              # the clipped image.
              tol_abs <- round(plusminus, 3)
              notes <- character(0)
              if (isTRUE(input$clip_outliers)) notes <- c(notes, sprintf("clipped at %s%%", format(clip_pct_value())))
              if (isTRUE(input$floor_clip)) notes <- c(notes, sprintf("floor at %s%%", format(floor_pct_value())))
              note_txt <- if (length(notes)) paste0(" (", paste(notes, collapse = ", "), ")") else ""
              vals_list <- list()
              for (k in seq_along(mz_set)) {
                idx <- which(is.finite(mz_axis_all) & abs(mz_axis_all - mz_set[k]) <= tol_abs)
                if (length(idx) == 0L) {
                  showNotification(sprintf("No features within +/- %s of m/z %.4f", tol_abs, mz_set[k]),
                                   type = "warning", duration = 6)
                  next
                }
                spec_k <- try(Cardinal::spectra(overview_peaks_sel_masked[idx, ]), silent = TRUE)
                if (inherits(spec_k, "try-error")) next
                mat_k <- coerce_spectra_matrix(spec_k, length(idx), n_pix_masked)
                if (is.null(mat_k)) next
                col_k <- sprintf("m/z = %.4f +/- %s%s", mz_set[k], tol_abs, note_txt)
                vals_list[[col_k]] <- transform_ion_values(colSums(mat_k, na.rm = TRUE))
              }
              if (length(vals_list) == 0L) {
                graphics::plot.new()
                graphics::text(0.5, 0.5, "No ion image could be built for these settings.")
                dev.off()
                return(make_plot3_image_result(outfile, "No ion image", .dims))
              }
              ion_zlim <- NULL
              if (isTRUE(input$fixed_scale)) {
                ion_zlim <- fixed_zlim(unlist(vals_list, use.names = FALSE))
                vals_list <- lapply(vals_list, function(v) pmin(pmax(v, ion_zlim[1]), ion_zlim[2]))
              }
              for (nm in names(vals_list)) {
                pData(overview_peaks_sel_masked)[[nm]] <- vals_list[[nm]]
              }
              arg_list <- list(overview_peaks_sel_masked,
                               names(vals_list),
                               col = cpal(input$color3),
                               enhance = if (is.null(ion_zlim)) enhance_option else NULL,
                               smooth = smoothing_option,
                               scale = if (is.null(ion_zlim)) input$normalize3 else FALSE,
                               superpose = input$superpose,
                               key = (input$colorkey3),
                               grid = FALSE)
              if (!is.null(ion_zlim)) arg_list$zlim <- ion_zlim
            } else {
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
            }
            
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
              width = grid::unit(.dims$w * allInputs$scalex, "points"),
              height = grid::unit(.dims$h * allInputs$scaley, "points"),
              just = c("center", "center")
            )
          )
          
          (grid::pushViewport(grid::viewport(x = 0.5, y = 0.5, width = 1, height = 1)))
          (grid::grid.draw(histology_grob))
          (grid::popViewport())
          
        }
        
        dev.off()
        
        # Return a list containing the filename
        make_plot3_image_result(outfile, "This is alternate text", .dims)
    }

    observe({
      output$plot3_pk <- shiny::bindEvent(renderImage({
        plot3_draw()
      }, deleteFile = TRUE), plot3_trigger())
    })

    # Batch export: one PNG per selected m/z row, drawn with the current
    # display settings (mask, runs, colours, clip/floor, fixed scale, size).
    observeEvent(input$save_png_batch, {
      if (!identical(input$ion_viz3, "custom")) {
        showNotification("Batch export needs 'Custom single / multiple' ions.", type = "warning", duration = 6)
        return(invisible(NULL))
      }
      if (!is.null(input$display_mode) && !identical(input$display_mode, "none")) {
        showNotification("Batch export saves plain ion images -- set 'Ion math?' to none first.", type = "warning", duration = 8)
        return(invisible(NULL))
      }
      if (isTRUE(input$plot_pdata)) {
        showNotification("Untick 'Plot Phenotype data?' to batch-export ion images.", type = "warning", duration = 6)
        return(invisible(NULL))
      }
      fd <- tryCatch(as.data.frame(fData(overview_peaks_sel_r())), error = function(e) NULL)
      if (is.null(fd) || nrow(fd) == 0L) {
        showNotification("No dataset loaded.", type = "warning", duration = 6)
        return(invisible(NULL))
      }
      rows <- suppressWarnings(as.integer(input$tbl_rows_selected))
      rows <- rows[is.finite(rows) & rows >= 1L & rows <= nrow(fd)]
      if (length(rows) == 0L) {
        showNotification("Select the m/z rows to export in the table first.", type = "warning", duration = 6)
        return(invisible(NULL))
      }
      if (isTRUE(input$fixed_scale) && !is.finite(suppressWarnings(as.numeric(input$zmax)[1]))) {
        showNotification("Fixed scale has no max, so each image uses its own maximum.", type = "warning", duration = 8)
      }
      saved <- character(0)
      failed <- character(0)
      t0 <- proc.time()[["elapsed"]]
      plot3_batch_mode <<- TRUE
      tryCatch({
        withProgress(message = "Saving ion images", value = 0, {
          for (i in seq_along(rows)) {
            r <- rows[i]
            mzv <- suppressWarnings(as.numeric(fd[["mz"]][r]))
            lab <- plot3_label(rows = r)
            setProgress(value = (i - 1) / length(rows), detail = sprintf("%d of %d: %s", i, length(rows), lab))
            tmp <- tempfile(fileext = ".png")
            res <- tryCatch(plot3_draw(outfile = tmp, ion_override = mzv), error = function(e) {
              message(sprintf("[plot_card:%s] batch export failed for %s: %s", id, lab, conditionMessage(e)))
              NULL
            })
            # Placeholder images ("No ion image", "No pixels ...") are not saved.
            if (!is.null(res) && identical(res$alt, "This is alternate text") && file.exists(tmp)) {
              dest <- unique_png_dest(png_filename(lab))
              if (isTRUE(file.copy(tmp, dest, overwrite = FALSE))) saved <- c(saved, dest) else failed <- c(failed, lab)
            } else {
              failed <- c(failed, lab)
            }
            unlink(tmp)
          }
          setProgress(value = 1, detail = "done")
        })
      }, finally = {
        plot3_batch_mode <<- FALSE
      })
      el <- proc.time()[["elapsed"]] - t0
      message(sprintf("[plot_card:%s] batch export: %d saved, %d failed in %.1fs -> %s",
                      id, length(saved), length(failed), el, getwd()))
      if (length(saved)) {
        showNotification(sprintf("Saved %d image(s) to %s", length(saved), getwd()), type = "message", duration = 10)
      }
      if (length(failed)) {
        showNotification(paste("Could not export:", paste(failed, collapse = ", ")), type = "error", duration = 12)
      }
    }, ignoreInit = TRUE)
    
    
    if(spatialOnly==FALSE) {
      observe({  
        output$plot_ranges<- renderUI( {
          
          req(overview_peaks_sel_r())
          
          lims <- get_safe_plot_limits(overview_peaks_sel_r())
          
          
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
          req(overview_peaks_sel_r())
          lims <- get_safe_plot_limits(overview_peaks_sel_r())
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
          
          req(overview_peaks_sel_r())
          
          
          lims <- get_safe_plot_limits(overview_peaks_sel_r())
          
          
          numericInput(ns("param_numeric"),
                       "Manual y-axis intensity max value",
                       min = round(lims$y[1], 0),
                       max = round(lims$y[2], 0),
                       value = round(lims$y[2], 0)
          )
        })
      })
      
      
      observe({
        output$plot4 <- renderImage( {
          
          
          req(overview_peaks_sel_r())
          req(input$mass_range_plot)
          
          # A temp file to save the output.
          # This file will be removed later by renderImage
          outfile <- tempfile(fileext = '.png')
          
          # audit PLOTCARD-12: an empty/non-numeric/non-positive size field used to make png() error
          # (red box); fall back to the UI defaults (800 x 300) only in that case. Any valid
          # positive size is passed through unchanged, so every plot that rendered before is identical.
          w_ms <- suppressWarnings(as.numeric(input$width_ms)[1])
          if (length(w_ms) == 0 || !is.finite(w_ms) || w_ms <= 0) w_ms <- 800
          h_ms <- suppressWarnings(as.numeric(input$height_ms)[1])
          if (length(h_ms) == 0 || !is.finite(h_ms) || h_ms <= 0) h_ms <- 300
          png(outfile, width = w_ms, height = h_ms)
          .plot4_dev <- grDevices::dev.cur()
          on.exit({
            if (.plot4_dev %in% grDevices::dev.list()) try(grDevices::dev.off(.plot4_dev), silent = TRUE)
          }, add = TRUE)
          
          
          
          if(length(input$int_range_plot)==1) {
            ylim=c(0, input$int_range_plot)
          } else {
            ylim = input$int_range_plot
          }
          
          #change xlimits based on custom ion or not
          if(is.null(input$x_target) || is.na(input$x_target)){
            xlim=input$mass_range_plot
            
            overview_peaks_sel_plot<-overview_peaks_sel_r()
            
          } else {
            xlim=c(input$x_target-input$x_tol, input$x_target+input$x_tol)
            
            #subsetFeatures to only include mz values within range
            overview_peaks_sel_plot<-subsetFeatures(overview_peaks_sel_r(), mz > xlim[1] & mz < xlim[2])
            
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
              # audit PLOTCARD-12: previously a bare return() handed renderImage NULL (blank/error
              # output); now close the device and return the already-drawn spectrum, exactly as the
              # other early-exit branches here do. Nothing drawn or computed changes.
              dev.off()
              return(list(src = outfile,
                          contentType = 'image/png',
                          width = input$width,
                          height = input$height,
                          alt = "This is alternate text"))
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
              # audit PLOTCARD-12: previously a bare return() handed renderImage NULL (blank/error
              # output); now close the device and return the already-drawn spectrum, exactly as the
              # other early-exit branches here do. Nothing drawn or computed changes.
              dev.off()
              return(list(src = outfile,
                          contentType = 'image/png',
                          width = input$width,
                          height = input$height,
                          alt = "This is alternate text"))
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
