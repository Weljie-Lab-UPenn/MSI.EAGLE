
  plot_card_UI<-function(id) {
    
    ns <- NS(id)
    
  
    tagList(  
      
        fluidRow(
          tags$head(
            #tags$style("label{font-family: BentonSans Book;}")
            #tags$style("label{font-size: 11px;} ")
          ),
          
          column(3, offset = 0,
                 radioButtons(ns("ion_viz3"), "Image visualization ions",
                                 c( "All"="viz_all", "Custom"="custom")),
                 #numericInput(ns("mz_viz3"), "mz value for visualization",255.2),
                 uiOutput(ns("mz_viz3a")),
                 
                 fluidRow(
                   column(6, selectInput(ns("contrast3"), "Contrast enhancement", c( "none", "histogram",  "adaptive"))),
                   column(6, selectInput(ns("color3"), "Colorscale", c(hcl.pals()), selected="Spectral"))
                 ),
                 fluidRow(
                   column(6, selectInput(ns("smooth3"), "Smoothing options", c("none", "mean", "gaussian", "bilateral", "adaptive", "diffusion", "guided"))),
                   column(6, checkboxInput(ns("normalize3"), "Scale multiple images?", value = TRUE))
                   ),
                   
                 
                 checkboxInput(ns("colorkey3"), "Draw color key?", value = TRUE),
                 fluidRow(
                   column(6, numericInput(ns("width_im"), "Image plot width (px)", value = 800)),
                   column(6, numericInput(ns("height_im"), "Image plot height (px)", value=600))
                 ),
                 checkboxInput(ns("plot_pdata"), "Plot Phenotype data?", value=FALSE),
                 uiOutput(ns("plotpdata")),
                 checkboxInput(ns("expand_fonts"), "Extended font options?", value=FALSE),
                 uiOutput(ns("fonts")),                 
                 checkboxInput(ns("expand_runs"), "Select individual runs for plotting only?", value=FALSE),
                 uiOutput(ns("select_runs")),
                 p("___________________________________________________________________")
                 

                 
          ),
          column(9, uiOutput(ns("plot.window"))
          )  
        ),
        
        uiOutput(ns("spectrum"))
    )
          
  }
  

  plot_card_server <- function(id, overview_peaks_sel, spatialOnly=FALSE) {

    moduleServer(id, function(input, output, session){
      
      #for dynamic UI
      ns = session$ns
      
      graphics.off()
      
      #create new overview_peaks_sel object with mean values
      overview_peaks_sel<-Cardinal::summarizeFeatures(overview_peaks_sel)
      
      observe({
        
        output$plot.window <- renderUI({
          
          if(input$ion_viz3=="custom") {
            list(
              imageOutput(ns("plot3_pk") ), #, width="800px", height="600px"),
              DT::dataTableOutput(ns('tbl'))
              
            )
          } else {
            list(
              imageOutput(ns("plot3_pk") )) #, width="800px", height="600px")
            
          }
          
        })
        
      })
      
      
      observe({
        
         output$mz_viz3a <- renderUI ({
           req(overview_peaks_sel)
           
           if(input$ion_viz3=="custom") {
             
            
             mz_list<- mz(overview_peaks_sel)
             tbl<-as.data.frame(fData(overview_peaks_sel))
             if(!is.null(tbl$freq)) {
               tbl$freq<-round(tbl$freq, 2)
             }
             if(!is.null(tbl$mean)) {
               tbl$mean<-round(tbl$mean, 1)
             }
             tbl$mz<-round(tbl$mz, 4)
             
             #add something here at some point to only select mz, freq, count, and mean columns
             
             
             output$tbl <-DT::renderDataTable({
               tbl
               
             }, selection = "multiple")
  
             list(
               # selectInput(ns("mz_viz3"),
               #           label= "m/z selection",
               #           multiple=TRUE,
               #           choices = mz_list
               #           ),
               
               
               checkboxInput(ns("superpose"), "Superpose images?", value=FALSE),
               selectInput(ns("display_mode"), "Ion math?", c("none", "sum", "ratio", "subtract", "min", "max", "mean", "sd", "var",  "multiply")),
               #checkboxInput(ns("display_mode"), 
              #               "Sum m/z TICs?",
              #                FALSE),
               numericInput(ns("plusminus_viz3"), "plusminus value for visualization",0.05)
             )
           }
         
          })
         
        
      })
      
      mz_viz3 <- reactive({
        req(overview_peaks_sel)
        req(input$tbl_rows_selected)
        tbl<-as.data.frame(fData(overview_peaks_sel))
        return(tbl[input$tbl_rows_selected, "mz"])
        
      })

      

      observe({
        
        output$select_runs <- renderUI ({
          req(overview_peaks_sel)
          
          if(input$expand_runs) {
            
            run_list<- unique(run((overview_peaks_sel)))
            
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
              selectInput(ns("pdata_var_plot"), label="pData variable to plot", choices = colnames(as.data.frame(pData(overview_peaks_sel))))
            )
          }
          
          
        })
      })

       #add spectrum or not
      output$spectrum<-renderUI({
        
        if(input$ion_viz3=="custom") {
          DT::dataTableOutput(ns('tbl'))
        }

        
        if(spatialOnly==FALSE){
          tagList(
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
            imageOutput(ns("plot4")), #, width="800px", height="600px")
            fluidRow(
              column(3, checkboxInput(ns("calc_ppm"), label = "Show ppm error?", width = '100%')),
              column(3, checkboxInput(ns("show_int"), label = "Show intensity?", width = '100%')),
              column(3, checkboxInput(ns("show_mz"), label = "Show mz value?", width = '100%'))
            ),
            checkboxInput(ns("spectrum_expand_fonts"), "Extended font options?", value=FALSE),
            uiOutput(ns("spectrum_fonts")),
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
              # column(
              #   2,
              #   numericInput(
              #     ns("spectrum_title_size"),
              #     label = "Title font scaling (%)",
              #     min = 0,
              #     value = 100
              #   )
              # ),
              column(
                3,
                numericInput(
                  ns("spectrum_label_size"),
                  label = "Label font scaling (%)",
                  min = 0,
                  value = 100
                )
              ),
              # column(
              #   2,
              #   numericInput(
              #     ns("spectrum_subtitle_size"),
              #     label = "Subtitle font scaling (%)",
              #     min = 0,
              #     value = 100
              #   )
              # ),
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
          
       observe({
        output$plot3_pk <- renderImage( {  #plot image in overview after peakpicking / reading file
          
          #req(overview_peaks_sel)
          
          # A temp file to save the output.
          # This file will be removed later by renderImage
          outfile <- tempfile(fileext = '.png')
          
          png(outfile, width = input$width_im, height = input$height_im)
          
          #ion <- switch(input$mode,
          #              "p"=786,
          #             "n"=255.2)
          #browser()
          
          if(!is.null(input$select_runs)) {
            overview_peaks_sel <- subsetPixels(overview_peaks_sel, run %in% input$select_runs)
          }
          
          #set sizes
          if(input$expand_fonts) {
            cex.axisp=input$axis_size/100
            cex.labp=input$label_size/100
            cex.mainp=input$title_size/100
            cex.subp=input$subtitle_size/100
            
            #get margins
            cur_mar<-par()$mar
            
            new_mar<-c(cur_mar[1]+cex.labp/8, cur_mar[2]+cex.labp/2, cur_mar[3]+cex.mainp/2, cur_mar[4])
            
            #if drawing colorkey, add a little on the left
            if(input$colorkey3) {
              new_mar[4]=new_mar[4]+cex.axisp
            }
            
            cur_mgp<-par()$mgp
            
            #new_mgp<-c(cur_mgp[1]+cex.labp/8, cur_mgp[2], 0)
            new_mgp<-c(3+max(new_mar[1:2])/20, 1,0)
            
            
            
            
            
          } else {
            cex.axisp=1
            cex.labp=1
            cex.mainp=1
            cex.subp=1
            new_mar<-par()$mar
            new_mgp<-par()$mgp
            
          }
          
          if(input$plot_pdata){
            req(input$pdata_var_plot)
            print(Cardinal::image(overview_peaks_sel,
                                  as.data.frame(pData(overview_peaks_sel))[,input$pdata_var_plot]~x*y,
                                  key=(input$colorkey3),
                                  superpose=input$superpose,
                                  col=pals::alphabet()),
                                  cex.axis=req(cex.axisp),
                                  cex.lab=cex.labp,
                                  cex.main=cex.mainp,
                                  cex.sub=cex.subp,
                                  mar=new_mar,
                                  mgp=new_mgp
                  )
          } else if (input$ion_viz3=="viz_all") {
            
            
            mz_range=range(mz(overview_peaks_sel))
            
            #find closest mz value to middle of range
            test_value <- mean(mz_range)
            
            
            # Calculate the absolute differences
            differences <- abs(mz(overview_peaks_sel) - test_value)
            
            # Find the index of the minimum difference
            closest_index <- which.min(differences)
            
            mz_set=mz(overview_peaks_sel)[closest_index]
            tol=max(differences) + differences[closest_index]+1
            
            smoothing_option <- if (input$smooth3 != "none") paste0(", smooth ='", input$smooth3,"'") else ""
            enhance_option <- if (input$contrast3 != "none") paste0(", enhance ='", input$contrast3,"'") else ""
            
            plusminus=tol
            
            image_command <- paste("Cardinal::image(overview_peaks_sel, 
                                  mz=mz_set,
                                  tolerance=round(plusminus,3), 
                                  units='mz',
                                  col=cpal(input$color3)",
                                  enhance_option,
                                  smoothing_option,",
                                  scale=input$normalize3,
                                  #superpose=input$superpose,
                                  key=(input$colorkey3),
                                  cex.axis=req(cex.axisp),
                                  cex.lab=cex.labp,
                                  cex.main=cex.mainp,
                                  cex.sub=cex.subp,
                                  mar=new_mar,
                                  mgp=new_mgp
            )")
            
            
            print(eval(parse(text = image_command)))

            
            
            
          } else if (input$ion_viz3=="custom"){
            
            if(is.null(mz_viz3())){
              ion=mz(overview_peaks_sel[1,])
            } else {
              
              # observeEvent(input$mz_viz3,{
                ion=as.numeric(mz_viz3())
              # })
            }
            
            
            if(!is.null(input$display_mode) && input$display_mode!="none"){
              if(input$display_mode%in%c("min", "max", "min", "mean", "sum", "sd", "var")) { 
          
                
                browser()
                select_vec<-as.character(mz(overview_peaks_sel)) %in% as.character(ion)
                #test to make sure there are 2 or more elements
                if(sum(select_vec)<2){
                  showNotification("At least two ions required for this calculation", type="error")
                  message("At least two ions required for this calculation")
                  return()
                }
                
                sm<-summarizePixels(overview_peaks_sel[select_vec,], stat=c(xic=input$display_mode), as="DataFrame")
                pData(overview_peaks_sel)$xic<-sm$xic
            
              }else if(input$display_mode=="ratio"){
                if(length(ion)!=2){
                  showNotification("Exactly two ions required for ratio (mz1/mz2)", type="error")
                  message("Exactly two ions required for ratio (mz1/mz2)")
                  return()
                } else {
                  mz1 <- spectra(subsetFeatures(overview_peaks_sel, mz=ion[1]))[1,]
                  mz2 <- spectra(subsetFeatures(overview_peaks_sel, mz=ion[2]))[1,]
                  xic <- (1 + mz1) / (1 + mz2)
                }
                
              } else if(input$display_mode=="subtract") {
                if(length(ion)!=2){
                  showNotification("Exactly two ions required for subtraction (mz1-mz2)", type="error")
                  message("Exactly two ions required for subtraction (mz1-mz2)")
                  return()
                } else {
                  mz1 <- spectra(subsetFeatures(overview_peaks_sel, mz=ion[1]))[1,]
                  mz2 <- spectra(subsetFeatures(overview_peaks_sel, mz=ion[2]))[1,]
                  xic <- (mz1) - (mz2)
                }
                
                
              }else if(input$display_mode=="multiply"){
                nelements=length(ion)
                xic <- 1+spectra(subsetFeatures(overview_peaks_sel, mz=ion[1]))[1,]
                
                for(i in 2:nelements){
                  mz2= 1+spectra(subsetFeatures(overview_peaks_sel, mz=ion[i]))[1,]
                  xic=(xic) * (mz2)
                }
                }else if(input$display_mode=="subtract") {
                
                mz1 <- spectra(subsetFeatures(overview_peaks_sel, mz=ion[1]))[1,]
                mz2 <- spectra(subsetFeatures(overview_peaks_sel, mz=ion[2]))[1,]
                xic <- (mz1) - (mz2)
                
              }
            
  
              browser()
              
              if(sum(is.na(pData(overview_peaks_sel)$xic))==length(overview_peaks_sel)) {
                showNotification("This calculation does not work!")
                return()
              }
              
              print(Cardinal::image(overview_peaks_sel, xic~x*y, 
                    enhance=input$contrast3,
                    colorscale=matter::cpal(input$color3),
                    smooth=input$smooth3,
                    scale=input$normalize3,
                    superpose=input$superpose,
                    key=input$colorkey3,
                    cex.axis=cex.axisp,
                    cex.lab=cex.labp,
                    cex.main=cex.mainp,
                    cex.sub=cex.subp,
                    mar=new_mar,
                    mgp=new_mgp
                    ))
              
            } else {
              
              
              
              mz_set=ion
              
              mz_range=range(mz(overview_peaks_sel))
              
              
              smoothing_option <- if (input$smooth3 != "none") paste0(", smoothing ='", input$smooth3,"'") else ""
              enhance_option <- if (input$contrast3 != "none") paste0(", enhance ='", input$contrast3,"'") else ""
              
              plusminus=input$plusminus_viz3
              
              image_command <- paste("Cardinal::image(overview_peaks_sel, 
                                  mz=mz_set,
                                  tolerance=round(plusminus,3), 
                                  units='mz',
                                  col=cpal(input$color3)",
                                     enhance_option,
                                     smoothing_option,",
                                  scale=input$normalize3,
                                  #superpose=input$superpose,
                                  key=(input$colorkey3),
                                  cex.axis=req(cex.axisp),
                                  cex.lab=cex.labp,
                                  cex.main=cex.mainp,
                                  cex.sub=cex.subp,
                                  mar=new_mar,
                                  mgp=new_mgp
            )")
              
              
              print(eval(parse(text = image_command)))

              
            }
          }
          
          
          dev.off()
          
          
          # Return a list containing the filename
          list(src = outfile,
               contentType = 'image/png',
               width = input$width_im,
               height = input$height_im,
               alt = "This is alternate text")
        }, deleteFile = TRUE)
       }) 
       
       
       if(spatialOnly==FALSE) {
         observe({  
          output$plot_ranges<- renderUI( {
            
            req(overview_peaks_sel)
            
            a<-Cardinal::plot(overview_peaks_sel)
            
  
                  sliderInput(ns("mass_range_plot"), 
                              label = p("m/z range for MS plot (X)"), 
                              min = round(a$channels$x$limits[1]),
                              max = round(a$channels$x$limits[2]), 
                              value = round(a$channels$x$limits), 
                              step = NULL,
                              round = TRUE)
  
            
          })
         })
         
         observe({
          output$plot_ranges2<- renderUI( {
            req(overview_peaks_sel)
            #browser()
            #overview_peaks_sel<-Cardinal::summarizeFeatures(overview_peaks_sel)
            a<-Cardinal::plot(overview_peaks_sel, "mean")
            
            
           
            
            sliderInput (ns("int_range_plot"), 
                        label = p("intensity range for MS plot (Y)"), 
                        min = round(a$channels$y$limits[1],0),
                        max = round(a$channels$y$limits[2],0)*1.05, 
                        value = req(input$param_numeric), #round(a$par$ylim,0), 
                        step = a$channels$y$limits[2]/20,
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
             
             req(overview_peaks_sel)
            
             
             a<-Cardinal::plot(overview_peaks_sel, "mean")
             
             
             numericInput(ns("param_numeric"),
                          "Manual y-axis intensity max value",
                          min = round(a$channels$y$limits[1],0),
                          max = round(a$channels$y$limits[2],0),
                          value = round(a$channels$y$limits[2],0)
             )
           })
         })
               
          
         if(!is.null(overview_peaks_sel)) {
           
           updateSliderInput(session, ns("int_range_plot"), value = c(round(Cardinal::plot(req(overview_peaks_sel))$channels$y$limits[1],0), 
                                                                      round(Cardinal::plot(req(overview_peaks_sel))$channels$y$limits[1],0)))
           updateNumericInput(session, ns("param_numeric"), value = round(Cardinal::plot(req(overview_peaks_sel))$channels$y$limits[2],0))
          }
        
        
        
         
         observe({
           output$plot4 <- renderImage( {
             
             
             req(overview_peaks_sel)
             req(input$mass_range_plot)
             
             # A temp file to save the output.
             # This file will be removed later by renderImage
             outfile <- tempfile(fileext = '.png')
             
             png(outfile, width = input$width_ms, height = input$height_ms)
             
            
             
             if(length(input$int_range_plot)==1) {
               ylim=c(0, input$int_range_plot)
             } else {
               ylim = input$int_range_plot
             }
             
             #change xlimits based on custom ion or not
             if(is.null(input$x_target) || is.na(input$x_target)){
               xlim=input$mass_range_plot
               
               
               
             } else {
               xlim=c(input$x_target-input$x_tol, input$x_target+input$x_tol)
               
             }
             
             if(!is.finite(xlim[1])){
               xlim=input$mass_range_plot
             }
             
             
             if(input$spectrum_expand_fonts) {
               cex.axisp=input$spectrum_axis_size/100
               cex.labp=input$spectrum_label_size/100
               #cex.mainp=input$spectrum_title_size/100
               #cex.subp=input$spectrum_subtitle_size/100
               lwdp=input$linewidth/100
               
               
               #get margins
               cur_mar<-par()$mar
               
               new_mar<-c(cur_mar[1]+cex.labp/8, cur_mar[2]+cex.labp/2, cur_mar[3], cur_mar[4])
               
               cur_mgp<-par()$mgp
               
               #new_mgp<-c(cur_mgp[1]+cex.labp/8, cur_mgp[2], 0)
               new_mgp<-c(3+max(new_mar[1:2])/20, 1,0)
                          
               
               
               
               
             } else {
               cex.axisp=1
               cex.labp=1
               #cex.mainp=1
               #cex.subp=1
               lwdp=1
               new_mar<-par()$mar
               new_mgp<-par()$mgp
             }
             
             #browser()
             #overview_peaks_sel<-Cardinal::summarizeFeatures(overview_peaks_sel)
             
             p1<-Cardinal::plot(overview_peaks_sel,
                                xlim=xlim,
                                ylim =ylim,
                                cex.axis=req(cex.axisp),
                                cex.lab=cex.labp,
                                #cex.main=cex.mainp,
                                #cex.sub=cex.subp,
                                lwd=lwdp,
                                mar=new_mar, mgp=new_mgp, "mean")
             print(p1)
             
             #check for ppm calc
             if(input$calc_ppm) {
               #browser()
               dat=overview_peaks_sel
               x=mz(dat)
               targ_mz<-req(input$x_target)
               x_sel<-subset(x, x>xlim[1] & x< xlim[2])
               
               ppm_error<- round(1e6*(x_sel-targ_mz)/targ_mz, 2)
               
               y_labs<-p1$facets[[1]][[1]]$y[p1$facets[[1]][[1]]$x %in% x_sel]
               
               if(length(ppm_error)==0){
                 showNotification("No ppm error calculated, are there any peaks?", type="warning")
                 return()
               } else {
                print(text(x=x_sel, y=y_labs+ylim[2]*.05, req(ppm_error)))
               }
               
             }
             
             if(input$show_mz) {
               #browser()
               dat=overview_peaks_sel
               x=mz(dat)
               x_sel<-subset(x, x>xlim[1] & x< xlim[2])
               
               
               y_labs<-p1$facets[[1]][[1]]$y[p1$facets[[1]][[1]]$x %in% x_sel]
               
               if(length(x_sel)==0){
                 showNotification("No mz values in range, are there any peaks?", type="warning")
                 return()
               } else {
                 print(text(x=x_sel, y=y_labs+ylim[2]*.15, req(round(x_sel, 4))))
               }
               
             }
             
             if(input$show_int) {
               #browser()
               dat=overview_peaks_sel
               x=mz(dat)
               
               x_sel<-subset(x, x>xlim[1] & x< xlim[2])
               
               
               
               y_labs<-p1$facets[[1]][[1]]$y[p1$facets[[1]][[1]]$x %in% x_sel]
               
               if(length(y_labs)==0){
                 showNotification("No intensities found, are there any peaks?", type="warning")
                 return()
               } else {
                 print(text(x=x_sel, y=y_labs+ylim[2]*.10, req(round(y_labs, 0))))
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

