### R/HeatmapServer.R
HeatmapServer <- function(id,  proc_values) {
  moduleServer(id, function(input, output, session) {
    ns = session$ns #for dyanamic variable namespace
    x5 =  observe({
      proc_values()[["x5"]]
    }) # if more than one value in list; not sure how to do it otherwise
    #x5<- reactive({proc_values()$x5})
    
    
    ###HEATMAP####
    
    output$hmap_vars <-
      renderUI({
        #https://gist.github.com/wch/4211337
        # If missing input, return to avoid error later in function
        req(proc_values()$x5$stats_results)
        req(proc_values()$x5$phen_options)
        
        
        
        selectInput(
          ns("hmap_vars"),
          "Choose variables for heatmap",
          choices  = proc_values()$x5$phen_options,
          multiple = T,
          selected = proc_values()$x5$phen_cols_stats
        )
      })
    
    
    #    output$hm_sig_select<-renderUI(req("hm_sigOnly" %in% input$hmap_params), {  #https://gist.github.com/wch/4211337
    output$hm_sig_select <-
      renderUI({
        #https://gist.github.com/wch/4211337
        # If missing input, return to avoid error later in function
        req(proc_values()$x5$stats_results)
        
        
        hmap_choices <-
          colnames(proc_values()$x5$stats_results)[!colnames(proc_values()$x5$stats_results) %in%
                                                     c("mz", "feature")]
        
        list(
          selectInput(
            ns("hm_sig_select"),
            "Choose filtering column for heatmap",
            choices  = hmap_choices,
            selected = "fdr"
          ),
          radioButtons(
            ns("sig_direction"),
            label = "Filtering stat direction?",
            choices = list("Ascending" = "ascending",
                           "Descending" = "descending"),
            selected = "descending"
          ),
          numericInput(ns("hm_sig"), "Significance cutoff", value = .05)
        )
      })
    
    
    output$hm_labels <- renderUI({
      selectInput(
        ns("hm_labels"),
        label = "Choose labels",
        choices = list("m/z" = "mz",
                       "ID" = "ID"),
        selected = "mz"
      )
    })
    
    x6 <- reactiveValues(
      dat_long_tech_avg = NULL,
      hmap_plot = NULL,
      hmap_dat = NULL,
      annotation_row = NULL
    )
    
    observeEvent(input$action_hmap, {
      if (is.null(proc_values()$x5$data_file_selected)) {
        message("No stats data, please run a stats analysis first")
        return()
      }
      withProgress(message = "Generating heatmap...", detail = "summarizing data", {
        message("Generating heatmap!")
        
        
        
        #data_file_selected<-x5$data_file_selected
        data_file_selected <- proc_values()$x5$data_file_selected
        
        #library(data.table)
        numdat <- as.matrix(spectra(data_file_selected))
        rownames(numdat) <- mz(data_file_selected)
        pdat <- as.data.frame(pData(data_file_selected))
        
        # Convert to data.table for faster binding
        pdat_dt <- data.table::as.data.table(pdat)
        numdat_dt <- data.table::as.data.table(t(numdat))
        
        # Ensure column names are unique and compatible
        data.table::setnames(numdat_dt, make.unique(colnames(numdat_dt)))
        
        # Bind columns using data.table's `cbind`
        dat_dt <- cbind(pdat_dt, numdat_dt)
        
        # Convert back to data.frame if necessary
        dat <- as.data.frame(dat_dt)
        
        
        
        
        
        dat_long <-
          (
            tidyr::pivot_longer(
              dat,
              cols = (ncol(pdat) + 1):ncol(dat),
              names_to = "mz",
              values_to = "response"
            )
          )
        
        
        #dat_long_tech_avg<- dat_long %>% group_by(eval(input$phen_interaction_stats), mz) %>%
        dat_long_tech_avg <-
          dat_long %>% dplyr::group_by_at(c(unique(
            c(
              input$grouping_variables,
              input$output_factors,
              input$hmap_vars
            )
          ), "mz")) %>%
          dplyr::summarize(tech_avg = mean(response), .groups = "keep")
        
        print("data summarized!")
        
        if (length(colnames(dat_long_tech_avg)) == 3 &&
            colnames(dat_long_tech_avg)[!colnames(dat_long_tech_avg) %in% c("mz", "tech_avg")] ==
            "Plate") {
          message(
            "Plate is the only variable in the data-- are you sure? Check Stats tab for modeling and export"
          )
        }
        # FIX THIS!!
        # dat<-proc_values()$x5$data_file_selected
        # #create interaction of hmap_vars for grouping from pData of dat
        # group_vars <- unique(c(input$grouping_variables, input$output_factors))
        # pData(dat)$interaction <- interaction(pData(dat)[, group_vars])
        # 
        # summarizeFeatures(proc_values()$x5$data_file_selected, groups = interaction(input$hmap_vars))
        # 
        x6$dat_long_tech_avg <- dat_long_tech_avg
        
      })
    })
    
    output$plot14 <- renderImage({
      req(x6$dat_long_tech_avg)
      req(input$hmap_vars)
      # A temp file to save the output.
      # This file will be removed later by renderImage
      outfile <- tempfile(fileext = '.png')
      
      png(outfile, width = 800, height = 600)
      
      #library(pheatmap)
      
      
      
      
      if ("hm_sigOnly" %in% input$hmap_params) {
        req(input$hm_sig_select)
        
        stats_results <- proc_values()$x5$stats_results
        if (is.null(stats_results$ID)) {
          #try to get IDs from featureData
          stats_results <- stats_results[order(stats_results$mz),]
          
          IDs <- fData(proc_values()$x5$data_file_selected)$ID
          
          if (length(IDs) == dim(stats_results)[1]) {
            stats_results$ID <- IDs
          } else{
            message("could not find IDs in stats_file or featureData")
            message("check stats_table length to make sure it matches number of mzs")
            stats_results$ID <- NA
          }
          
          
        }
        
        if (input$sig_direction == "descending") {
          sig_mz <-
            stats_results %>% dplyr::filter(.data[[input$hm_sig_select]] <= input$hm_sig) %>% dplyr::select(mz, ID)
        } else {
          sig_mz <-
            stats_results %>% dplyr::filter(.data[[input$hm_sig_select]] >= input$hm_sig) %>% dplyr::select(mz, ID)
        }
        x6$dat_long_tech_avg <-
          x6$dat_long_tech_avg[round(as.numeric(x6$dat_long_tech_avg$mz), 4) %in% round(sig_mz[, 1], 4),]
      } else {
        sig_mz <- proc_values()$x5$stats_results[, "mz", drop = F]
        if (is.null(proc_values()$x5$stats_results$ID)) {
          sig_mz$ID <- NA
        } else {
          sig_mz$ID <- proc_values()$x5$stats_results$ID
          sig_mz$ID[sig_mz$ID == ""] <- NA
        }
        stats_results <- proc_values()$x5$stats_results
      }
      
      if (all(is.na(stats_results$ID))) {
        stats_results <-
          stats_results[,!colnames(stats_results) %in% c("ID")]
        sig_mz <- data.frame(mz = sig_mz[, "mz"])
      }
      
      if (dim(sig_mz)[1] == 0) {
        message ("No m/z values selected. Check FDR cutoff and other selection parameters")
        return()
      }
      
      
      
      
      dat_short <-
        tidyr::pivot_wider(x6$dat_long_tech_avg,
                           names_from = mz,
                           values_from = tech_avg)
      
      if (length(complete.cases(dat_short)) != sum(complete.cases(dat_short))) {
        print(
          "NAs detected in dataset, only complete cases will be used!
              Be sure to check the input data carefully"
        )
        dat_short <- dat_short[complete.cases(dat_short),]
        
      }
      
      
      #dat_short_num<-as.matrix(dat_short[,!colnames(dat_short)%in%c( unique(c(proc_values()$x5$group_var, input$hmap_vars)))])
      dat_short_num <-
        as.matrix(dat_short[,!colnames(dat_short) %in% c(input$hmap_vars, unique(c(proc_values()$x5$group_var)))])
      
      #look at first set
      #as.data.frame(dat_short_num) %>% dplyr::select(1:10)
      
      if (dim(dat_short_num)[2] < 2) {
        print("not enough peaks for heatmaps, check thresholds and variables")
        return()
      }
      
      rownames(dat_short_num) <- 1:dim(dat_short_num)[1]
      colnames(dat_short_num) <-
        as.character(round(as.numeric(colnames(
          dat_short_num
        )), 4))
      
      #need to make sure this is ordered by numeric ions...
      order(as.numeric(colnames(dat_short_num)))
      dat_short_num <-
        dat_short_num[, order(as.numeric(colnames(dat_short_num)))]
      
      if (sum(input$hmap_params %in% c("rowScale", "colScale")) == 0) {
        scale_val = "none"
      } else if (sum(input$hmap_params %in% c("rowScale", "colScale")) ==
                 2) {
        print("must choose row or column scaling!")
        scale_val = "none"
        return()
      } else if (c("rowScale") %in% input$hmap_params) {
        scale_val = "row"
      } else if (c("colScale") %in% input$hmap_params) {
        scale_val = "column"
      } else {
        scale_val = "none"
      }
      
      
      
      if ("hm_plotCols" %in% input$hmap_params) {
        annotation_col =  as.data.frame(dat_short[, colnames(dat_short) %in% input$hmap_vars])
        rownames(annotation_col) = 1:dim(dat_short_num)[1]
        
        
        if (!proc_values()$x5$phen_cols_stats %in% colnames(annotation_col)) {
          print(
            "column annotation requires the Test variable to be set in the Stats tab and match a heatmap variable"
          )
          return()
        }
        
        
        if (length(levels(as.factor(annotation_col[, 1]))) < 11) {
          annotation_colors = list()
          annotation_colors[[1]] = ggsci::pal_npg()(length(levels(as.factor(
            annotation_col[, proc_values()$x5$phen_cols_stats]
          ))))
          names(annotation_colors[[1]]) <-
            levels(as.factor(annotation_col[, proc_values()$x5$phen_cols_stats]))
          names(annotation_colors) <-
            proc_values()$x5$phen_cols_stats
          
        } else {
          print("too many colors to match columns with npg palette, using alphabet")
          #annotation_col=NA
          annotation_colors = list(pals::alphabet2())
        }
        
        
      } else {
        annotation_col = NA
        annotation_colors = NA
      }
      
      if (!"rowclust" %in% input$hmap_params) {
        cluster_rows = F
      } else{
        cluster_rows = T
      }
      
      if (!"colclust" %in% input$hmap_params) {
        cluster_cols = F
      } else{
        cluster_cols = T
      }
      
      if (!"hmap_legend" %in% input$hmap_params) {
        mylegend = T
      } else{
        mylegend = F
      }
      #browser()
      
      if ("hm_plotRows" %in% input$hmap_params) {
        annotation_row = (sig_mz[order(sig_mz$mz), , drop = F])
        rownames(annotation_row) <- annotation_row[, "mz"]
        
        if (dim(annotation_row)[1] < 2) {
          print("No IDs, check targeted mass file and be sure ID column is included")
        }
        
        annotation_row = annotation_row[,!colnames(annotation_row) %in% "mz", drop =
                                          F]
        annotation_row[annotation_row == ""] <- NA
        
        
        
      } else {
        annotation_row = NA
      }
      
      if ("hm_rowNames" %in% input$hmap_params) {
        show_rownames = F
      } else {
        show_rownames = T
      }
      
      if ("hm_colNames" %in% input$hmap_params) {
        show_colnames = F
      } else {
        show_colnames = T
      }
      
      #palettes...
      
      if(input$heatmap_colors!="default") {
        
        
        col_choices<-c(hcl.pals())
        
        hmap_colors<-
          switch(input$heatmap_colors,
                 viridis = cpal("viridis")(100),
                 spectral=cpal("Spectral")(100),
                 plasma=cpal("plasma")(100),
                 inferno=cpal("inferno")(100),
                 cividis=cpal("cividis")(100)
          )
        #browser()
      } else (
        
        hmap_colors=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                                    "RdYlBu")))(100)
      )
      
      
      #browser()
      if (is.numeric(sig_mz[, input$hm_labels])) {
        labels_row = round(sig_mz[order(sig_mz$mz), input$hm_labels], 4)
      } else {
        labels_row = sig_mz[order(sig_mz$mz), input$hm_labels]
      }
      
      
      x6$annotation_row <- labels_row
      
      #browser()
      
      x6$hmap_dat <- as.matrix(t(dat_short_num))
      
      hmap <- pheatmap::pheatmap(
        x6$hmap_dat,
        scale = scale_val,
        annotation_col = annotation_col,
        annotation_colors = annotation_colors,
        cluster_rows = cluster_rows,
        cluster_cols = cluster_cols,
        clustering_method = input$clust_method,
        clustering_distance_columns = input$clustering_distance_columns,
        clustering_distance_rows = input$clustering_distance_rows,
        cutree_rows = input$hm_cut_row,
        cutree_cols = input$hm_cut_col,
        kmeans_k = input$hm_k,
        border_color = NA,
        legend = mylegend,
        drop_levels = T,
        annotation_row = annotation_row,
        show_colnames = show_colnames,
        show_rownames = show_rownames,
        labels_row = labels_row,
        color=hmap_colors
      )
      
      x6$hmap_plot <- hmap
      print(hmap)
      #print(pheatmap((t(dat_short_num)), scale=scale_val ))# , annotation_row = annotation_row))
      #browser()
      
      dev.off()
      
      # Return a list containing the filename
      list(
        src = outfile,
        contentType = 'image/png',
        width = 800,
        height = 600,
        alt = "This is alternate text"
      )
    }, deleteFile = TRUE)
    
    output$Export = downloadHandler(
      filename = function() {
        "heatmap_plot.pdf"
      },
      content = function(file) {
        pdf(file, onefile = TRUE)
        print(x6$hmap_plot)
        message(
          "writing data to 'heatmap_dat_raw.tsv and heatmap_dat_scaled.tsv; will overwrite existing file!"
        )
        #browser()
        
        dat_raw <- x6$hmap_dat
        dat_scaled <- t(scale(t(x6$hmap_dat)))
        
        mz <- x6$hmap_plot$tree_row$labels
        mz_ord <- x6$hmap_plot$tree_row$order
        
        sample_ord <- x6$hmap_plot$tree_col$order
        sample_lab <- x6$hmap_plot$tree_col$labels
        
        dat_raw <- dat_raw[mz_ord, sample_ord]
        dat_scaled <- dat_scaled[mz_ord, sample_ord]
        
        
        readr::write_tsv(as.data.frame(dat_raw) %>% tibble::add_column(ID =
                                                                         x6$annotation_row[mz_ord]),
                         file = "heatmap_dat_raw.tsv")
        readr::write_tsv(
          as.data.frame(dat_scaled) %>% tibble::add_column(ID = x6$annotation_row[mz_ord]),
          file = "heatmap_dat_scaled.tsv"
        )
        
        dev.off()
      }
    )
    
  })
}
