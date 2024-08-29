### R/PhenoServer.R
PhenoServer <- function(id,  setup_values, preproc_values) {
  moduleServer(id, function(input, output, session){
    # output$table <- DT::renderDataTable(mtcars)
    
    ns = session$ns #for dyanamic variable namespace
    
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
      
      
      if(is.null(x3$img.dat) | is.null(x3$txt_data))
        return()
      #try(if(x3$img.dat=="proc") #no processed data available
      #  return())
    
      if(input$phen_method=="breaks"){
        req(input$breaks_thresh)
      }
      
      if(is.null(x3$txt_data$Plate)){
        print("No 'Plate' column detected in sample list, please check!")
        return()
      }
      
      
      
      #reorder input file to match datafile
      plates=unique(x3$txt_data$Plate)
      
      #browser()
      #reorder input file so plates are in the sample order
      #first only keep overlapping set of runs?
      
      plates<-plates[sapply(plates, function(x) grep(x, runNames(x3$img.dat)))]
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
      if(is.null(input$int_cols_out))
        return()
      
      if(length(input$int_cols_out)!=2) {
        print("Need exactly 2 fields selected for interaction mapping")
      }
      
      int=interaction(as.data.frame(x3$pdata)[,input$int_cols_out[1]], as.data.frame(x3$pdata)[,input$int_cols_out[2]], drop=T)
      x3$pdata[, paste0(input$int_cols_out, collapse = ".")]<-int
    })
    
    
    
    observeEvent(input$save_imzml, {
      
      req(input$save_imzml)
      
      
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
