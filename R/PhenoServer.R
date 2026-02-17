### R/PhenoServer.R
PhenoServer <- function(id,  setup_values, preproc_values) {
  moduleServer(id, function(input, output, session){
    # output$table <- DT::renderDataTable(mtcars)
    
    ns = session$ns #for dyanamic variable namespace
    
    output$variables <- renderPrint({
      setup_vals <- setup_values()
      paste("working directory=", setup_vals$wd)
    })
    
    # Integer-safe tiling from coord(msi)
    # Modes:
    #   - Fixed tile size:   w, h (pixels). nx, ny derived from data range.
    #   - Fixed tile count:  nx, ny (number of tiles). Ignores w, h if both nx,ny given.
    # Notes:
    #   * Uses inclusive integer ranges [min, max]; last break is max+1 so findInterval catches the max.
    #   * Distributes remainders (if any) to the first intervals so widths are all integers.
    #   * Returns a factor like "t003_007". Also returns tile_x/tile_y if requested.
    
    .compute_tile_factor <- function(msi_obj,
                                     w = NULL, h = NULL,
                                     nx = NULL, ny = NULL,
                                     add_indices = FALSE) {
      coords <- as.data.frame(Cardinal::coord(msi_obj))
      if (!all(c("x", "y") %in% names(coords))) {
        stop("coord(msi) does not contain columns 'x' and 'y'.")
      }
      if (nrow(coords) == 0L) {
        stop("No coordinates found.")
      }
      
      # Inclusive integer span
      x_min <- min(coords$x, na.rm = TRUE)
      x_max <- max(coords$x, na.rm = TRUE)
      y_min <- min(coords$y, na.rm = TRUE)
      y_max <- max(coords$y, na.rm = TRUE)
      
      # Defensive: require integer-like coordinates
      if (any(is.na(coords$x)) || any(is.na(coords$y))) {
        stop("Missing values in coord(msi).")
      }
      if (!all(coords$x == as.integer(coords$x)) ||
          !all(coords$y == as.integer(coords$y))) {
        warning("Non-integer coords detected; coercing to integer via as.integer().")
        coords$x <- as.integer(coords$x)
        coords$y <- as.integer(coords$y)
        # Recompute bounds just in case
        x_min <- min(coords$x); x_max <- max(coords$x)
        y_min <- min(coords$y); y_max <- max(coords$y)
      }
      
      # Helper: inclusive integer breaks producing n_intervals tiles
      # Returns a length (n_intervals + 1) vector of edges where edges[1]=start, edges[end]=end+1
      integer_breaks_inclusive <- function(start, end, n_intervals) {
        total <- (end - start + 1L)                    # inclusive width
        base  <- total %/% n_intervals
        extra <- total %%  n_intervals
        widths <- rep.int(base, n_intervals)
        if (extra > 0L) widths[seq_len(extra)] <- widths[seq_len(extra)] + 1L
        c(start, start + cumsum(widths))               # last edge is end+1
      }
      
      # Determine nx, ny if user supplied w/h (tile size)
      if (!is.null(nx) && !is.null(ny)) {
        # use fixed tile counts
        stopifnot(nx >= 1L, ny >= 1L)
      } else {
        # fall back to fixed tile size; derive nx, ny from ranges
        if (is.null(w) || is.null(h)) {
          # sensible defaults if nothing provided
          w <- 10L; h <- 10L
        }
        w <- as.integer(w); h <- as.integer(h)
        if (w < 1L || h < 1L) stop("Tile sizes w and h must be >= 1.")
        xrange <- (x_max - x_min + 1L)
        yrange <- (y_max - y_min + 1L)
        nx <- ceiling(xrange / w)
        ny <- ceiling(yrange / h)
      }
      
      # Build integer edges; remainders are distributed across the first tiles
      x_edges <- integer_breaks_inclusive(x_min, x_max, nx)  # length nx+1, last is x_max+1
      y_edges <- integer_breaks_inclusive(y_min, y_max, ny)  # length ny+1, last is y_max+1
      
      # Map coordinates to 1..nx and 1..ny with inclusive right edge
      # findInterval returns indices in 1..length(edges)-1 when using our edges (last is end+1)
      tx <- findInterval(coords$x, x_edges, rightmost.closed = TRUE)
      ty <- findInterval(coords$y, y_edges, rightmost.closed = TRUE)
      
      # Safety: clamp any 0 (below first edge) or >n (above last) due to numeric quirks
      tx[tx < 1L] <- 1L;      tx[tx > nx] <- nx
      ty[ty < 1L] <- 1L;      ty[ty > ny] <- ny
      
      tile_lab <- sprintf("t%03d_%03d", tx, ty)
      out <- factor(tile_lab)
      
      if (isTRUE(add_indices)) {
        attr(out, "tile_x") <- tx
        attr(out, "tile_y") <- ty
        attr(out, "x_edges") <- x_edges
        attr(out, "y_edges") <- y_edges
      }
      out
    }
    
    
    
    ####PHENOTYPING####
    #create new reactive object for phenotyping section since
    #this is largely independent now.
    
    x3=reactiveValues(img.dat=NULL,   #image to be phenotyped
                      img_p.dat=NULL,  #slot for image file with phenotype data
                      pdata=NULL,     #currentphenotype data  
                      txt_data=NULL, #input phenotype info
                      fdata=NULL)  #featuredata slot
    
    
    output$phenotype_table <- DT::renderDataTable({
      file <- input$phenotype_file
      ext <- tools::file_ext(file$datapath)
      
      req(file)
      validate(need(ext == "txt", "Please upload a tab delimited .txt file"))
      
      x3$txt_data<-read_samples(file$datapath, method="auto")
      
      #browser()
      
      if(which(names(x3$txt_data)%in%"Spot") !=3) {
        message("Spot is not in first column, phenotyping will likely not work!")
        showNotification("Spot is not in first column, phenotyping will likely not work!", type="error")
      }
      
      #code for select button selecting only filtered rows from here:
      #https://community.rstudio.com/t/select-only-filtered-rows-using-select-all-button-that-comes-with-select-extension-in-shinys-dt-package/66749/2
      DT::datatable(x3$txt_data,escape=F,
                    rownames=F,
                    filter = 'top',
                    #  colnames = c("Data Type","Variable","Description", "Filename"),
                    class = "compact hover row-border",
                    extensions = c('Scroller','Select', 'Buttons'),
                    
                    options = list(
                      select = list(style = "multi", items = "row"),
                      columnDefs = list(list(className = 'dt-center', targets = "_all")),
                      language = list(
                        info = 'Showing _START_ to _END_ of _TOTAL_ variables'),
                      deferRender = TRUE,
                      scrollY = 500,
                      scroller = TRUE,
                      dom = "Blfrtip",
                      buttons = list(list(extend='selectAll',className='selectAll',
                                          text="select all rows",
                                          action=DT::JS("function () {
                                var table = $('#DataTables_Table_0').DataTable();
                                table.rows({ search: 'applied'}).deselect();
                                table.rows({ search: 'applied'}).select();
                }")
                      ), list(extend='selectNone',
                              text="DeselectAll",
                              action=DT::JS("function () {
                                var table = $('#DataTables_Table_0').DataTable();
                                table.rows({ search: 'applied'}).select();
                                table.rows({ search: 'applied'}).deselect();
                }")
                      ))
                      
                    ),
                    selection="none"
      ) }, server = F
      )
    
    
    observeEvent(input$read_phenotype, {
      
      x2<-preproc_values()[["x2"]]
      
      #browser()
      if(is.null(x2$list_proc_img) && is.null(x2$overview_peaks_sel)) {
        print("open a dataset or used stored data to start")
        showNotification("open a dataset or used stored data to start", type="error")
        return()
      }
        
      
      
      x3$img.dat<-switch(input$img_phen, 
                         "phen_no_bk"="proc",
                         "phen_pp"=x2$overview_peaks_sel)
      
      
      if(is.character(x3$img.dat) && x3$img.dat=="proc") {
        if(is.null(x2$list_proc_img)) {
          print("No processed data, be sure to store data first.")
          showNotification("No processed data found. Store processed data first, then click 'Read dataset to be phenotyped'.", type = "warning", duration = 8)
          return()
        } else {
          x3$img.dat<-combine_card(x2$list_proc_img)
        }
      }
      
      x3$pdata<-pData(x3$img.dat)
      x3$fdata<-fData(x3$img.dat)
      
    })
    
    output$phen_method_opts<- renderUI({
      switch(input$phen_method,
             "period"=NULL,
             "breaks"=list(
               numericInput(ns("breaks_thresh"), "Threshold for breaks method",4,  min=2 )
              )
      )
      
    })
    
    
    
    output$phen_cols <- renderUI({  #https://gist.github.com/wch/4211337
      # If missing input, return to avoid error later in function
      #if(is.null(x3$pdata))
      #  return()
      req(x3$txt_data)
      # Get the data set with the appropriate name
      
      dat <- (as.data.frame(x3$txt_data))
      cols <- colnames(dat)
      print(head(cols))
      
      # Create the checkboxes and select them all by default
      checkboxGroupInput(ns("phen_cols_out"), "Choose columns to add to pData", 
                         choices  = cols,
                         selected = "Plate")
    }) 
    
    observeEvent(input$start_phenotype, {
      
      #browser()
      if(is.null(x3$img.dat) | is.null(x3$txt_data)) {
        showNotification("Load both MSI data and phenotype table first, then click 'Start phenotyping'.", type = "warning", duration = 8)
        return()
      }
      #try(if(x3$img.dat=="proc") #no processed data available
      #  return())
    
      if(input$phen_method=="breaks"){
        req(input$breaks_thresh)
      }
      
      if(is.null(x3$txt_data$Plate)){
        print("No 'Plate' column detected in sample list, please check!")
        showNotification("No 'Plate' column detected in phenotype table.", type = "error", duration = 8)
        return()
      }
      
      
      
      #reorder input file to match datafile
      plates=unique(x3$txt_data$Plate)
      
      #browser()
      #reorder input file so plates are in the sample order
      #first only keep overlapping set of runs?
      
      plates<-try(plates[sapply(plates, function(x) grep(x, runNames(x3$img.dat)))])
      
      if(class(plates)=="try-error") {
        message("Plate names / run names do not match, please check!")
        showNotification("Plate names / run names do not match, please check!", type="error")
        
        print(plates)
        print(runNames(x3$img.dat))
        return()
      }
      
      x3$img.dat<-combine_card(lapply(1:length(plates), function(x) x3$img.dat[,grep(plates[x], Cardinal::run(x3$img.dat))]))
      
      #ADD code to remove all other ID columns?
      a<-as.data.frame(fData(x3$img.dat)) %>% dplyr::select(!contains("ID."))
      
      if(!is.null(a$ID)) {
      
        Cardinal::fData(x3$img.dat)<-Cardinal::MassDataFrame(mz(x3$img.dat), ID=a$ID)
      }
      
      
      
      if(input$phen_method %in% c("spec_density", "period", "breaks")) {
        #browser()
        x3$img_p.dat <- pixDatFill_mult(datas = x3$img.dat, 
                                      sample_list = x3$txt_data, 
                                      variables = input$phen_cols_out, 
                                      method=input$phen_method, 
                                      inflect_thresh=input$breaks_thresh,
                                      lsp_plot=FALSE)
      #tmp_pdata<-pixDatFill_mult(x3$img.dat, x3$txt_data, input$phen_cols_out, method=input$phen_method, inflect_thresh=input$breaks_thresh)
      
      
      } else if (input$phen_method %in% "manual"){
        #browser()
        
        x3$img_p.dat<-pixDatFill_manual(datas = x3$img.dat, 
                                     sample_list = x3$txt_data, 
                                     variables = input$phen_cols_out)
        
      }
      
      x3$pdata<-x3$img_p.dat
      x3$fdata<-fData(x3$img.dat)
      
      
    })
    
    
    output$phenotype_table_new <- DT::renderDataTable({
      
      DT::datatable(
        as.data.frame(x3$pdata)
      )
    })
    
    output$plot_cols <- renderUI({  #https://gist.github.com/wch/4211337
      # If missing input, return to avoid error later in function
      #if(is.null(x3$pdata))
      #  return()
      req(x3$pdata)
      # Get the data set with the appropriate name
      
      dat <- (as.data.frame(x3$pdata))
      cols <- colnames(dat)
      print(head(cols))
      
      # Create the checkboxes and select them all by default
      radioButtons(ns("phen_plot"), "Choose pdata variable to plot", 
                   choices  = cols,
                   selected = "Plate")
    }) 
    
    
    output$phen_interaction <- renderUI({  #https://gist.github.com/wch/4211337
      # If missing input, return to avoid error later in function
      #if(is.null(x3$pdata))
      #  return()
      req(x3$pdata)
      
      # Get the data set with the appropriate name
      
      dat <- (as.data.frame(x3$pdata))
      cols <- colnames(dat)
      print(head(cols))
      
      # Create the checkboxes and select them all by default
      checkboxGroupInput(ns("int_cols_out"), "Choose 2 columns for an interaction term", 
                         choices  = cols,
                         selected = "")
    }) 
    
    observeEvent(input$interaction, {
      if(is.null(input$int_cols_out)) {
        showNotification("Select two pData fields to build an interaction term.", type = "warning", duration = 7)
        return()
      }
      
      if(length(input$int_cols_out)!=2) {
        print("Need exactly 2 fields selected for interaction mapping")
        showNotification("Need exactly 2 fields selected for interaction mapping.", type = "warning", duration = 7)
        return()
      }
      
      int=interaction(as.data.frame(x3$pdata)[,input$int_cols_out[1]], as.data.frame(x3$pdata)[,input$int_cols_out[2]], drop=T)
      x3$pdata[, paste0(input$int_cols_out, collapse = ".")]<-int
    })
    
    observeEvent(input$make_tiles, {
      if (is.null(x3$img.dat)) {
        showNotification("No dataset loaded for tiling. Read and phenotype data first.", type = "warning", duration = 8)
        return()
      }
      
      # safely pull user inputs with defaults
      w <- input$tile_w; if (is.null(w) || is.na(w) || w < 1) w <- 10
      h <- input$tile_h; if (is.null(h) || is.na(h) || h < 1) h <- 10
      
      # compute factor from coord(msi)
      tile_fac <- try(.compute_tile_factor(x3$img.dat, w = w, h = h, add_indices = TRUE), silent = TRUE)
      if (inherits(tile_fac, "try-error")) {
        showNotification("Could not compute tiles from coord(msi).", type = "error")
        return()
      }
      
      # ensure pdata exists and is aligned; x3$pdata is set after read_phenotype
      if (is.null(x3$pdata)) {
        # initialize from current pData if needed
        x3$pdata <- (Cardinal::pData(x3$img.dat))
      }
      
      # length check
      if (nrow(as.data.frame(x3$pdata)) != length(tile_fac)) {
        showNotification("Length mismatch between pData and coord(msi).", type = "error")
        return()
      }
      
      # add/replace column 'tile' (factor)
      tmp <- (x3$pdata)
      tmp$tile <- tile_fac
      tx <- attr(tile_fac, "tile_x")
      ty <- attr(tile_fac, "tile_y")
      if (!is.null(tx) && length(tx) == nrow(as.data.frame(tmp))) tmp$tile_x <- as.integer(tx)
      if (!is.null(ty) && length(ty) == nrow(as.data.frame(tmp))) tmp$tile_y <- as.integer(ty)
      x3$pdata <- tmp
      
      showNotification(sprintf("Added 'tile' to pData (%dx%d px, %d unique tiles).", w, h, length(unique(as.character(tile_fac)))), type = "message")
    })
    
    observeEvent(input$save_imzml, {
      
      req(input$save_imzml)
      if (is.null(x3$img.dat) || is.null(x3$pdata)) {
        showNotification("No phenotyped dataset to save. Read and phenotype data first.", type = "error", duration = 8)
        return()
      }
      
      
      volumes <- c(wd = setup_values()[["wd"]], home = fs::path_home())
      shinyFiles::shinyFileSave(input, "save_imzml", roots = volumes, session = session)
      
      
      save_path <- shinyFiles::parseSavePath(volumes, input$save_imzml)
      if (nrow(save_path) == 0) return(NULL)
      
      #browser()
      filen <- as.character(save_path$datapath)
      pk_img <- x3$img.dat
      pData(pk_img)<-x3$pdata
      
      
      fdata_tmp<-fData(x3$img.dat)
      if(dim(fdata_tmp)[2]>1) {
        coln<-unique(colnames(fdata_tmp))
        #browser()
        #only keep one ID column to prevent all sort of issues. would be a problem if additional featuredata is kept longer term
        
        fdata_unique<-fdata_tmp[,unique(colnames(fdata_tmp))]
        
        fData(pk_img)<-MassDataFrame(mz=fdata_unique$mz, fdata_unique %>% as.data.frame() %>% dplyr::select(!mz))
      }
      
        
        writeImzML(pk_img, filen)
      }
    )
    
    #plot to make sure it looks okay
    output$plot8 <- renderImage({
      req(x3$pdata)
      
      req(input$phen_plot)
        
      
      #browser()
      # A temp file to save the output.
      # This file will be removed later by renderImage
      outfile <- tempfile(fileext = '.png')
      
      png(outfile, width = 800, height = 600)
      
      # #img.dat<-x3$img_p.dat  #should the be pdata?
      # img.dat<-x3$img.dat
      # p<-as.data.frame(x3$pdata)[,input$phen_plot]
      # 
      
      pk_img <- x3$img.dat
      pData(pk_img)<-x3$pdata
      
      #browser()
      
      if(isFALSE(input$key_tf) || is.null(input$key_tf)) {
        keyval=TRUE
      } else if(isTRUE(input$key_tf)) {
        keyval=FALSE
      }
      
      #browser()
      print(image(pk_img, input$phen_plot))
      
      
      #print(Cardinal::image(mytable_selected(), mz=ion, plusminus=input$plusminus_viz))
      
      dev.off()
      
      # Return a list containing the filename
      list(src = outfile,
           contentType = 'image/png',
           width = 800,
           height = 600,
           alt = "This is alternate text")
    }, deleteFile = TRUE)
    
    
    
  })
}
