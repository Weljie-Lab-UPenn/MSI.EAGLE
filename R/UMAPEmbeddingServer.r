### R/UMAPEmbeddingServer.R 

UMAPEmbeddingServer <- function(id, setup_values, preproc_values, preproc_values_umap = NULL) {
  moduleServer(id, function(input, output, session) {
    ns = session$ns
    
    # Store color choices/selection for embedding plots.
    color_state <- reactiveValues(
      choices = character(0),
      selected_colors = character(0),
      last_selected_var = NULL
    )
    selected_rows_state <- reactiveVal(NULL)
    embed_data_state <- reactiveVal(NULL)
    embed_source_sig <- reactiveVal(NULL)

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

    same_labels <- function(a, b) {
      identical(normalize_labels(a), normalize_labels(b))
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
      uniq <- sort(unique(labels))
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

    cluster_to_ids <- function(cluster_vec) {
      out <- rep(NA_integer_, length(cluster_vec))
      keep <- !is.na(cluster_vec)
      if (any(keep)) {
        out[keep] <- suppressWarnings(as.integer(cluster_vec[keep]))
      }
      out
    }

    ensure_pdata_col <- function(dat, col_name, default = NA_character_) {
      pdat <- pData(dat)
      if (!(col_name %in% colnames(pdat))) {
        pdat[[col_name]] <- rep(default, ncol(dat))
      }
      pData(dat) <- pdat
      dat
    }

    dataset_signature <- function(dat) {
      rn <- paste(runNames(dat), collapse = "|")
      mz <- suppressWarnings(as.numeric(fData(dat)$mz))
      mz <- mz[is.finite(mz)]
      mz_min <- if (length(mz) > 0) round(min(mz), 6) else NA_real_
      mz_max <- if (length(mz) > 0) round(max(mz), 6) else NA_real_
      paste(nrow(dat), ncol(dat), rn, mz_min, mz_max, sep = "::")
    }

    store_overview_to_proc_list <- function(x2, dat = NULL) {
      if (is.null(dat)) dat <- x2$overview_peaks_sel
      req(dat)
      all_runs <- runNames(dat)
      if (length(all_runs) == 0) return(FALSE)
      proc_list <- lapply(all_runs, function(rn) {
        img <- subsetPixels(dat, Cardinal::run(dat) == rn)
        a <- fData(img)[unique(colnames(fData(img)))]
        a_df <- as.data.frame(a)
        if ("mz" %in% colnames(a_df)) {
          other_cols <- setdiff(colnames(a_df), "mz")
          fData(img) <- MassDataFrame(mz = a_df$mz, a_df[, other_cols, drop = FALSE])
        }
        img
      })
      names(proc_list) <- all_runs
      x2$list_proc_img <- proc_list
      TRUE
    }

    active_overview <- reactive({
      dat <- embed_data_state()
      if (!is.null(dat)) return(dat)
      x2 <- get_x2()
      req(x2$overview_peaks_sel)
      x2$overview_peaks_sel
    })

    resolve_selected_ids <- function(n_runs) {
      ids <- suppressWarnings(as.integer(input$mytable_rows_selected))
      ids <- ids[is.finite(ids)]
      if (length(ids) == 0) {
        ids <- suppressWarnings(as.integer(isolate(selected_rows_state())))
        ids <- ids[is.finite(ids)]
      }
      if (length(ids) == 0 && n_runs > 0) {
        ids <- seq_len(n_runs)
      }
      ids <- ids[ids >= 1 & ids <= n_runs]
      unique(ids)
    }

    observe({
      x2 <- get_x2()
      src <- x2$overview_peaks_sel
      req(src)
      sig <- dataset_signature(src)
      prev_sig <- isolate(embed_source_sig())
      if (is.null(isolate(embed_data_state())) || is.null(prev_sig) || !identical(prev_sig, sig)) {
        embed_data_state(src)
        embed_source_sig(sig)
        if (is.null(isolate(selected_rows_state())) || length(isolate(selected_rows_state())) == 0) {
          n_runs <- length(runNames(src))
          if (n_runs > 0) selected_rows_state(seq_len(n_runs))
        }
      }
    })

    available_color_labels <- reactive({
      req(data_img(), input$color_var)
      vals <- unique(normalize_labels(pData(data_img())[[input$color_var]]))
      sort(vals)
    })

    selected_color_labels <- reactive({
      choices <- available_color_labels()
      if (length(choices) == 0) return(character(0))
      if (length(choices) > 20) return(character(0))
      sel <- normalize_labels(input$selected_colors)
      if (is.null(input$selected_colors) || length(sel) == 0) {
        return(choices)
      }
      sel <- choices[choices %in% sel]
      if (length(sel) == 0) choices else sel
    })

    observeEvent({
      data_img()
      input$color_var
      1
    }, {
      choices <- available_color_labels()
      color_state$choices <- choices
      if (length(choices) == 0 || length(choices) > 20) {
        if (!same_labels(color_state$selected_colors, character(0))) {
          color_state$selected_colors <- character(0)
        }
        color_state$last_selected_var <- input$color_var
        return()
      }
      prev <- normalize_labels(isolate(color_state$selected_colors))
      next_sel <- choices
      if (is.null(isolate(color_state$last_selected_var)) ||
          !identical(isolate(color_state$last_selected_var), input$color_var)) {
        next_sel <- choices
      } else {
        keep <- intersect(prev, choices)
        next_sel <- if (length(keep) > 0) choices[choices %in% keep] else choices
      }
      if (!same_labels(color_state$selected_colors, next_sel)) {
        color_state$selected_colors <- next_sel
      }
      color_state$last_selected_var <- input$color_var
    }, ignoreInit = FALSE)
    
    # List runs for selection
    output$mytable = DT::renderDataTable({
      ovps <- active_overview()
      rn <- runNames(ovps)
      n_runs <- length(rn)
      sel <- isolate(selected_rows_state())
      if (is.null(sel) || length(sel) == 0) {
        sel <- seq_len(n_runs)
      }
      sel <- as.integer(sel)
      sel <- sel[is.finite(sel) & sel >= 1 & sel <= n_runs]
      if (length(sel) == 0 && n_runs > 0) sel <- 1L
      
      DT::datatable(
        data.frame(
          run = rn,
          size = vapply(rn, function(x) sum(Cardinal::run(ovps) %in% x), integer(1)),
          stringsAsFactors = FALSE
        ),
        selection = list(mode = "multiple", selected = sel),
        caption = "Choose runs to process for UMAP visualization",
        extensions = c("Buttons"),
        options = list(
          dom = 'Bfrtip',
          buttons = list('pageLength', "copy")
        )
      )
    }, server = FALSE)

    observeEvent(input$mytable_rows_selected, {
      rows <- suppressWarnings(as.integer(input$mytable_rows_selected))
      rows <- rows[is.finite(rows)]
      if (length(rows) > 0) {
        selected_rows_state(rows)
      }
    }, ignoreInit = FALSE)
    
    # Get active dataset based on table selection
    data_img <- reactive({
      ovps <- active_overview()
      ids <- resolve_selected_ids(length(runNames(ovps)))
      req(length(ids) > 0)
      selected_rows_state(ids)
      ovps %>% 
        subsetPixels(Cardinal::run(ovps) %in% runNames(ovps)[ids])
    })

    observeEvent(input$rerun_embed_clustering, {
      x2 <- get_x2()
      full_data <- active_overview()
      req(full_data)

      methods <- input$embed_clustering_methods
      if (is.null(methods) || length(methods) == 0) {
        showNotification("Select at least one clustering method.", type = "warning", duration = 5)
        return()
      }

      all_runs <- runNames(full_data)
      run_ids <- resolve_selected_ids(length(all_runs))
      run_ids <- run_ids[run_ids >= 1 & run_ids <= length(all_runs)]
      if (length(run_ids) == 0) {
        showNotification("Select at least one run in the table first.", type = "warning", duration = 6)
        return()
      }
      selected_runs <- all_runs[run_ids]

      pix_mask <- Cardinal::run(full_data) %in% selected_runs
      if (!any(pix_mask)) {
        showNotification("No pixels found for selected run(s).", type = "warning", duration = 6)
        return()
      }

      pdat_full <- pData(full_data)
      umap_cols <- c("x_umap", "y_umap", "z_umap")
      if (!all(umap_cols %in% colnames(pdat_full))) {
        showNotification("Stored data has no x_umap/y_umap/z_umap columns. Run UMAP first.", type = "error", duration = 8)
        return()
      }

      embed_all <- as.matrix(pdat_full[pix_mask, umap_cols, drop = FALSE])
      valid_embed <- apply(embed_all, 1, function(z) all(is.finite(z)))
      n_valid <- sum(valid_embed)
      if (n_valid < 3) {
        showNotification("Not enough valid UMAP points for clustering.", type = "error", duration = 8)
        return()
      }
      embeddings <- embed_all[valid_embed, , drop = FALSE]

      k_clusters <- suppressWarnings(as.integer(input$embed_k_clusters))
      if (!is.finite(k_clusters) || k_clusters < 2) k_clusters <- 2L
      min_pts <- suppressWarnings(as.integer(input$embed_minPts))
      if (!is.finite(min_pts) || min_pts < 2) min_pts <- 2L
      eps_val <- suppressWarnings(as.numeric(input$embed_eps))
      if (!is.finite(eps_val) || eps_val <= 0) eps_val <- 0.5
      allow_highmem <- isTRUE(input$embed_allow_highmem)
      pairwise_limit <- suppressWarnings(as.integer(input$embed_pairwise_limit))
      if (!is.finite(pairwise_limit) || pairwise_limit < 2000) pairwise_limit <- 12000L
      spectral_limit <- suppressWarnings(as.integer(input$embed_spectral_limit))
      if (!is.finite(spectral_limit) || spectral_limit < 1000) spectral_limit <- 8000L
      hdbscan_limit <- max(25000L, pairwise_limit)

      format_bytes <- function(bytes) {
        if (!is.finite(bytes) || bytes < 0) return("NA")
        units <- c("B", "KB", "MB", "GB", "TB")
        idx <- if (bytes > 0) floor(log(bytes, 1024)) else 0
        idx <- max(0, min(idx, length(units) - 1))
        paste0(round(bytes / (1024^idx), 2), " ", units[idx + 1])
      }

      estimate_method_mem <- function(method_name, n_pts, dims, k) {
        n <- as.double(n_pts)
        d <- as.double(dims)
        k <- as.double(max(2, k))
        if (method_name %in% c("hierarchical", "kmedoids")) {
          return((n * (n - 1) / 2) * 8)
        }
        if (identical(method_name, "spectral")) {
          return((n * n * 8) * 2) # dense similarity + overhead
        }
        if (identical(method_name, "fuzzy")) {
          return((n * d * 8) + (n * k * 8))
        }
        if (identical(method_name, "som")) {
          return((n * d * 8) * 2)
        }
        if (identical(method_name, "mclust")) {
          return((n * d * 8) * 3)
        }
        (n * d * 8)
      }

      guard_reason <- function(method_name, n_pts) {
        if (allow_highmem) return(NULL)
        if (method_name %in% c("hierarchical", "kmedoids") && n_pts > pairwise_limit) {
          return(paste0(method_name, " skipped by safety guard (n=", n_pts, " > ", pairwise_limit, "). Enable high-memory override or subset runs."))
        }
        if (identical(method_name, "spectral") && n_pts > spectral_limit) {
          return(paste0("spectral skipped by safety guard (n=", n_pts, " > ", spectral_limit, "). Enable high-memory override or subset runs."))
        }
        if (identical(method_name, "hdbscan") && n_pts > hdbscan_limit) {
          return(paste0("hdbscan skipped by safety guard (n=", n_pts, " > ", hdbscan_limit, "). Enable high-memory override or subset runs."))
        }
        if (identical(method_name, "mclust") && n_pts > pairwise_limit) {
          return(paste0("mclust skipped by safety guard (n=", n_pts, " > ", pairwise_limit, "). Enable high-memory override or subset runs."))
        }
        NULL
      }

      write_diag <- isTRUE(input$embed_write_diag)
      diag_file <- NULL
      if (write_diag) {
        wd <- tryCatch(as.character(setup_values()[["wd"]]), error = function(e) "")
        if (!nzchar(wd) || !dir.exists(wd)) wd <- getwd()
        diag_file <- file.path(wd, paste0("umap_embedding_cluster_diag_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
      }
      log_diag <- function(...) {
        msg <- paste0(...)
        line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", msg)
        message("[UMAPEmbedding] ", line)
        if (!is.null(diag_file)) {
          try(write(line, file = diag_file, append = TRUE), silent = TRUE)
        }
      }

      if (!is.null(diag_file)) {
        try(write("# UMAP Embedding clustering diagnostics", file = diag_file), silent = TRUE)
      }
      log_diag("selected_runs=", paste(selected_runs, collapse = ","), " methods=", paste(methods, collapse = ","), " n_selected_pixels=", sum(pix_mask), " n_valid_umap=", n_valid, " dims=", ncol(embeddings))
      log_diag("params: k=", k_clusters, " minPts=", min_pts, " eps=", eps_val)
      log_diag("guards: allow_highmem=", allow_highmem, " pairwise_limit=", pairwise_limit, " spectral_limit=", spectral_limit, " hdbscan_limit=", hdbscan_limit)
      log_diag("gc_pre=", paste(capture.output(gc()), collapse = " | "))

      run_one_method <- function(method_name) {
        if (identical(method_name, "kmeans")) {
          set.seed(123)
          km <- stats::kmeans(embeddings, centers = k_clusters)
          return(cluster_to_ids(km$cluster))
        }
        if (identical(method_name, "hierarchical")) {
          hc <- stats::hclust(stats::dist(embeddings), method = "ward.D2")
          return(cluster_to_ids(stats::cutree(hc, k = k_clusters)))
        }
        if (identical(method_name, "dbscan")) {
          db <- dbscan::dbscan(embeddings, eps = eps_val, minPts = min_pts)
          return(cluster_to_ids(db$cluster))
        }
        if (identical(method_name, "hdbscan")) {
          hdb <- dbscan::hdbscan(embeddings, minPts = min_pts)
          return(cluster_to_ids(hdb$cluster))
        }
        if (identical(method_name, "spectral")) {
          if (!requireNamespace("RSpectra", quietly = TRUE)) {
            stop("Spectral clustering requires package 'RSpectra'.")
          }
          sim <- exp(-as.matrix(stats::dist(embeddings))^2 / 2)
          eig <- RSpectra::eigs(sim, k = k_clusters)
          set.seed(123)
          sp <- stats::kmeans(eig$vectors, centers = k_clusters)
          return(cluster_to_ids(sp$cluster))
        }
        if (identical(method_name, "kmedoids")) {
          if (!requireNamespace("cluster", quietly = TRUE)) {
            stop("K-medoids clustering requires package 'cluster'.")
          }
          showNotification("Running K-medoids clustering. This method can be slow on large datasets.", type = "message", duration = 8)
          km <- cluster::pam(embeddings, k = k_clusters)
          return(cluster_to_ids(km$clustering))
        }
        if (identical(method_name, "fuzzy")) {
          if (!requireNamespace("e1071", quietly = TRUE)) {
            stop("Fuzzy C-means requires package 'e1071'.")
          }
          fz <- e1071::cmeans(embeddings, centers = k_clusters)
          return(cluster_to_ids(apply(fz$membership, 1, which.max)))
        }
        if (identical(method_name, "mclust")) {
          if (!requireNamespace("mclust", quietly = TRUE)) {
            stop("Mclust requires package 'mclust'.")
          }
          mc <- mclust::Mclust(embeddings, G = 1:k_clusters)
          return(cluster_to_ids(mc$classification))
        }
        if (identical(method_name, "som")) {
          if (!requireNamespace("kohonen", quietly = TRUE)) {
            stop("SOM requires package 'kohonen'.")
          }
          scaled_embeddings <- scale(embeddings)
          som_grid <- kohonen::somgrid(
            xdim = max(2, round(k_clusters / 2)),
            ydim = max(2, round(k_clusters / 2)),
            topo = "hexagonal"
          )
          sm <- kohonen::som(scaled_embeddings, grid = som_grid, rlen = 100)
          return(cluster_to_ids(sm$unit.classif))
        }
        if (identical(method_name, "skmeans")) {
          if (!requireNamespace("skmeans", quietly = TRUE)) {
            stop("Spherical K-means requires package 'skmeans'.")
          }
          sk <- skmeans::skmeans(embeddings, k = k_clusters)
          return(cluster_to_ids(sk$cluster))
        }
        stop(paste0("Unsupported method: ", method_name))
      }

      # Remove legacy cluster_* columns from older embedding rerun logic.
      legacy_cluster_cols <- intersect(paste0("cluster_", methods), colnames(pData(full_data)))
      if (length(legacy_cluster_cols) > 0) {
        pdat <- pData(full_data)
        keep_cols <- setdiff(colnames(pdat), legacy_cluster_cols)
        pData(full_data) <- pdat[, keep_cols, drop = FALSE]
      }

      # Ensure output columns exist before assigning values.
      for (m in methods) {
        full_data <- ensure_pdata_col(full_data, paste0("col_", m))
      }

      method_status <- setNames(rep("pending", length(methods)), methods)

      withProgress(message = "Re-running clustering on existing UMAP", value = 0, {
        n_methods <- length(methods)
        for (ii in seq_along(methods)) {
          m <- methods[[ii]]
          incProgress(1 / n_methods, detail = paste("Method:", m))

          est_mem <- estimate_method_mem(m, n_valid, ncol(embeddings), k_clusters)
          log_diag("method=", m, " est_mem=", format_bytes(est_mem))
          reason <- guard_reason(m, n_valid)
          if (!is.null(reason)) {
            method_status[m] <- "skipped_guard"
            showNotification(reason, type = "warning", duration = 8)
            log_diag("method=", m, " status=skipped_guard reason=", reason)
            next
          }

          invisible(gc())
          t0 <- proc.time()
          lbl <- tryCatch(
            run_one_method(m),
            error = function(e) {
              showNotification(paste0("Failed ", m, ": ", conditionMessage(e)), type = "error", duration = 8)
              log_diag("method=", m, " status=failed error=", conditionMessage(e))
              NULL
            }
          )
          elapsed <- unname((proc.time() - t0)[["elapsed"]])
          if (is.null(lbl) || length(lbl) != n_valid) {
            method_status[m] <- "failed"
            log_diag("method=", m, " status=failed elapsed_s=", round(elapsed, 3))
            next
          }

          full_lbl <- rep(NA_integer_, sum(pix_mask))
          full_lbl[valid_embed] <- suppressWarnings(as.integer(lbl))

          col_col <- paste0("col_", m)
          pdat <- pData(full_data)
          if (!(col_col %in% colnames(pdat))) pdat[[col_col]] <- rep(NA_integer_, ncol(full_data))
          pdat[[col_col]] <- suppressWarnings(as.integer(pdat[[col_col]]))
          pdat[[col_col]][pix_mask] <- full_lbl
          pData(full_data) <- pdat
          method_status[m] <- "ok"
          log_diag("method=", m, " status=ok elapsed_s=", round(elapsed, 3), " unique_labels=", length(unique(full_lbl[!is.na(full_lbl)])))
          invisible(gc())
        }
      })

      embed_data_state(full_data)
      embed_source_sig(dataset_signature(full_data))
      x2$overview_peaks_sel <- full_data
      selected_rows_state(run_ids)
      try(store_overview_to_proc_list(x2, full_data), silent = TRUE)
      selected_row_idx <- run_ids
      session$onFlushed(function() {
        proxy <- DT::dataTableProxy("mytable", session = session)
        try(DT::selectRows(proxy, selected_row_idx), silent = TRUE)
      }, once = TRUE)
      summary_txt <- paste(paste(names(method_status), method_status, sep = ":"), collapse = ", ")
      log_diag("summary=", summary_txt)
      log_diag("gc_post=", paste(capture.output(gc()), collapse = " | "))

      showNotification(
        paste0(
          "Clustering re-run complete (", summary_txt, "). Added/updated columns: ",
          paste(unlist(lapply(methods, function(m) c(paste0("col_", m)))), collapse = ", ")
        ),
        type = "message",
        duration = 9
      )
      if (!is.null(diag_file)) {
        showNotification(paste0("Diagnostics saved: ", diag_file), type = "message", duration = 10)
      }
    })

    observeEvent(input$embed_store_proc, {
      x2 <- get_x2()
      dat <- active_overview()
      req(dat)
      x2$overview_peaks_sel <- dat
      ok <- try(store_overview_to_proc_list(x2, dat), silent = TRUE)
      if (inherits(ok, "try-error") || !isTRUE(ok)) {
        showNotification("Unable to store embedding results into processed data list.", type = "error", duration = 8)
        return()
      }
      showNotification(
        paste0("Stored embedding results for ", length(x2$list_proc_img), " run(s) into processed data list."),
        type = "message",
        duration = 7
      )
    })

    observeEvent(input$embed_save_imzml, {
      x2 <- get_x2()
      dat <- active_overview()
      req(dat)
      x2$overview_peaks_sel <- dat
      if (is.null(x2$list_proc_img) || length(x2$list_proc_img) == 0) {
        ok <- try(store_overview_to_proc_list(x2, dat), silent = TRUE)
        if (inherits(ok, "try-error") || !isTRUE(ok)) {
          showNotification("No processed data available to save. Store failed.", type = "error", duration = 8)
          return()
        }
      }
      wd_root <- tryCatch(as.character(setup_values()[["wd"]]), error = function(e) "")
      if (!nzchar(wd_root) || !dir.exists(wd_root)) wd_root <- getwd()
      volumes <- c(wd = wd_root, home = fs::path_home())
      shinyFiles::shinyFileSave(input, "embed_save_imzml", roots = volumes, session = session)
      save_path <- shinyFiles::parseSavePath(volumes, input$embed_save_imzml)
      if (nrow(save_path) == 0) return(NULL)
      filen <- as.character(save_path$datapath)

      pk_img_reorder <- x2$list_proc_img[order(sapply(x2$list_proc_img, ncol), decreasing = TRUE)]
      pk_img <- combine_card(pk_img_reorder)
      a <- fData(pk_img)[unique(colnames(fData(pk_img)))]
      a_df <- as.data.frame(a)
      if ("mz" %in% colnames(a_df)) {
        other_cols <- setdiff(colnames(a_df), "mz")
        fData(pk_img) <- MassDataFrame(mz = a_df$mz, a_df[, other_cols, drop = FALSE])
      }
      writeImzML(pk_img, filen)
      showNotification(paste0("Saved imzML: ", filen), type = "message", duration = 7)
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
      
      mz_values <- suppressWarnings(as.numeric(fData(data_img())$mz))
      mz_values <- mz_values[is.finite(mz_values)]
      if (length(mz_values) == 0) {
        return(tags$div(style = "color:#a94442;", "No m/z values available in current selection."))
      }
      names(mz_values) <- paste("m/z:", round(mz_values, 4))
      selected_mz <- suppressWarnings(as.numeric(isolate(input$mz_ion)))
      if (length(selected_mz) != 1 || is.na(selected_mz) || !is.finite(selected_mz)) {
        selected_mz <- mz_values[1]
      } else if (!isTRUE(any(abs(mz_values - selected_mz) < 1e-9))) {
        selected_mz <- mz_values[which.min(abs(mz_values - selected_mz))]
      }
      
      selectizeInput(ns("mz_ion"), 
                     "Select m/z Ion:",
                     choices = mz_values, 
                     selected = selected_mz,
                     options = list(maxOptions = 10000))
    })
    
    output$color_selector <- renderUI({
      req(data_img())
      
      color_vars <- colnames(pData(data_img()))
      if (length(color_vars) == 0) {
        return(tags$div(style = "color:#a94442;", "No pData fields available in current selection."))
      }
      selected_color <- isolate(input$color_var)
      if (!is.character(selected_color) || length(selected_color) != 1 || !(selected_color %in% color_vars)) {
        selected_color <- if ("Rcol_reduced" %in% color_vars) "Rcol_reduced" else color_vars[1]
      }
      
      selectInput(ns("color_var"),
                  "Select Color Variable:", 
                  choices = color_vars,
                  selected = selected_color)
    })
    
    # Selected m/z index
    selected_index <- reactive({
      req(data_img(), input$mz_ion)
      mz_all <- suppressWarnings(as.numeric(fData(data_img())$mz))
      req(length(mz_all) > 0)
      target <- suppressWarnings(as.numeric(input$mz_ion))
      if (length(target) != 1 || is.na(target) || !is.finite(target)) {
        finite_idx <- which(is.finite(mz_all))
        req(length(finite_idx) > 0)
        target <- mz_all[finite_idx[1]]
      }
      dif <- abs(mz_all - target)
      dif[!is.finite(dif)] <- Inf
      idx <- which.min(dif)
      req(length(idx) == 1, is.finite(idx), idx >= 1)
      idx
    })
    
    # UMAP plot with m/z intensity
    output$umapPlot <- renderImage({
      req(data_img(), selected_index(), umap_coords())
      
      outfile <- tempfile(fileext = '.png')
      png(outfile, width = 800, height = 600)
      
      # Get metadata values for masking
      color_values <- normalize_labels(pData(data_img())[[input$color_var]])
      selected_vals <- selected_color_labels()
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
      
      pairs(umap_coords(),
            col = plot_colors,
            pch = ".",
            cex = input$point_size)
      
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
      choices <- color_state$choices
      if (length(choices) > 20) {
        return(tags$div(
          style = "color: red;",
          "The selected variable has more than 20 unique values. Please select a different variable."
        ))
      }

      checkboxGroupInput(ns("selected_colors"),
                         "Select colors to display:",
                         choices = choices,
                         selected = color_state$selected_colors)
    })
    
    
    # Update selected colors when changed
    observeEvent(input$selected_colors, {
      choices <- isolate(color_state$choices)
      if (length(choices) == 0 || length(choices) > 20) {
        if (!same_labels(color_state$selected_colors, character(0))) {
          color_state$selected_colors <- character(0)
        }
        return()
      }
      sel <- normalize_labels(input$selected_colors)
      next_sel <- choices
      if (is.null(input$selected_colors) || length(sel) == 0) {
        next_sel <- choices
      } else {
        keep <- intersect(sel, choices)
        next_sel <- if (length(keep) > 0) choices[choices %in% keep] else choices
      }
      if (!same_labels(color_state$selected_colors, next_sel)) {
        color_state$selected_colors <- next_sel
      }
    }, ignoreInit = TRUE)
    
    # Modified UMAP plot with metadata colors and selection
    output$umapPlot2 <- renderImage({
      req(data_img(), input$color_var, umap_coords())
      
      outfile <- tempfile(fileext = '.png')
      png(outfile, width = 800, height = 600)
      
      # Get metadata values
      color_values <- normalize_labels(pData(data_img())[[input$color_var]])
      selected_vals <- selected_color_labels()
      
      # Create mask for selected colors
      mask <- color_values %in% selected_vals
      
      # Use a stable label->color map shared with spatial rendering.
      plot_colors <- rep("grey80", length(color_values))
      if (any(mask)) {
        mapping <- build_display_mapping(color_values[mask], match_original = isTRUE(input$col_match))
        plot_colors[mask] <- unname(mapping$values[mapping$labels])
      }

      pairs(umap_coords(),
            col = plot_colors,
            pch = ".",
            cex = input$point_size)
      
      dev.off()
      
      list(src = outfile,
           contentType = 'image/png',
           width = 800,
           height = 600,
           alt = "UMAP plot colored by metadata")
    }, deleteFile = TRUE)
    
    # Modified phenotype plot with color selection
    output$phenoplot <- renderPlot({
      tryCatch({
        dat <- try(data_img(), silent = TRUE)
        if (inherits(dat, "try-error") || is.null(dat)) {
          graphics::plot.new()
          graphics::text(0.5, 0.5, "Select at least one run to render the phenotype/spatial plot.")
        } else {
          pdat_cols <- colnames(pData(dat))
          color_var <- input$color_var
          has_color_var <- !is.null(color_var) &&
            length(color_var) == 1 &&
            !is.na(color_var) &&
            nzchar(color_var) &&
            color_var %in% pdat_cols
          plot_mz <- isTRUE(input$plot_mz)

          if (plot_mz || has_color_var) {
            if (plot_mz) {
              mz_val <- suppressWarnings(as.numeric(input$mz_ion))
              if (!is.finite(mz_val)) {
                graphics::plot.new()
                graphics::text(0.5, 0.5, "Invalid m/z value selected.")
              } else {
                Cardinal::image(
                  dat,
                  mz = mz_val,
                  tol = 0.005,
                  units = "mz",
                  enhance = "hist",
                  smooth = "gaussian"
                )
              }
            } else {
              img_vec <- normalize_labels(pData(dat)[[color_var]])
              select_vec <- selected_color_labels()
              keep <- img_vec %in% select_vec
              keep[is.na(keep)] <- FALSE

              subset_data <- try(subsetPixels(dat, keep), silent = TRUE)
              if (inherits(subset_data, "try-error")) {
                stop(paste0("subsetPixels failed: ", as.character(subset_data)))
              }

              if (ncol(subset_data) > 0) {
                subset_labels <- normalize_labels(pData(subset_data)[[color_var]])
                mapping <- build_display_mapping(subset_labels, match_original = isTRUE(input$col_match))
                pData(subset_data)$.__embed_color_label <- factor(
                  mapping$labels,
                  levels = names(mapping$values)
                )
                Cardinal::image(
                  subset_data,
                  ".__embed_color_label",
                  col = unname(mapping$values)
                )
              } else {
                graphics::plot.new()
                graphics::text(0.5, 0.5, "No selected labels to display.")
              }
            }
          } else {
            graphics::plot.new()
            graphics::text(0.5, 0.5, "Choose a valid color field or enable m/z plot.")
          }
        }
      }, error = function(e) {
        err <- paste0("Phenotype plot error: ", conditionMessage(e))
        message("[UMAPEmbedding] ", err)
        graphics::plot.new()
        graphics::text(0.5, 0.5, err, cex = 0.9)
      })
    }, width = 800, height = 600, res = 96)
    
  })
}
