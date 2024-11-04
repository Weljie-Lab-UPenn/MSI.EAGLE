#require(Cardinal)
library(magrittr)

#function to process data
HTS_reproc<-function(data1, SN, res=10, align_tol=15, pix_include=NULL,mzrange=NULL, freq.min=0.01, method="adaptive") {
  
  #remove z-axis as it messes with the segmentation algorithms
  coord(data1)$z<-NULL
  
  
  #remove non-relevant peaks for reprocessing-- does this make a difference?
  if(!is.null(pix_include)) {
    data1<-data1[coord(data1)%in%pix_include]
  }
  
  #setCardinalBPPARAM(MulticoreParam())
  
  
  data2 <- data1 %>% #reduceBaseline(method="median", blocks=100) %>%
    normalize(method="tic") %>% 
    peakPick(
      #pixel=seq(1, ncol(data1), by=1), 
      method=method, SNR=SN) %>%
    process()

  #peakAlign no longer available in Cardinal 3.6, but peakAlign adds info
  data3<- data2 %>% peakAlign(units="ppm", tolerance=align_tol) %>%
    process()
  
  
  #if no peaks meet filtering criteria, return original data3
  if(sum(fData(data3)$freq>freq.min)>0) {
    data4 <- subsetFeatures(data3, freq> freq.min)
  } else {
    message("No peaks meet filtering criteria, returning peak aligned data")
    data4 <- data3
  }

  return(data4)
}



#Use a set of reference peaks to normalize and peak bin rawdata####
#inputs are raw data file and peak list (eg mz_ref=mz(data))
ref_to_peaks <- function(raw_data, mz_ref, tol=15, units="ppm") {
  print(runNames(raw_data))
  data_peaks <- raw_data %>%
    normalize(method="tic") %>%
    peakBin(ref=mz_ref,
            tolerance=tol,
            units=units) %>%
    process()
  coord(data_peaks)$z <-NULL
  return(data_peaks)
}

#Use a set of reference peaks to normalize and peak bin rawdata####
#inputs are raw data file and peak list (eg mz_ref=mz(data)). 
#This version is for the new Cardinal version 3.6
ref_to_peaks_3_6 <- function(raw_data, mz_ref, tol=15, units_mz="ppm") {
  print(runNames(raw_data))
  data_peaks <- raw_data %>%
    #normalize(method="tic") %>%
    peakPick(ref=mz_ref[1:3],
             type="area",
            tolerance=tol,
            units=units_mz) %>%
    process()
  coord(data_peaks)$z <-NULL
  return(data_peaks)
}


# helper function to run umap on imaging dataset, visualize and look for outliers
#the thresh parameter can be optimzied as a percent of the outlier gap to the max distance value
# the umap parameters intially are optimized for outlier detection and so may cause artificial
# separation of otherwise desired clusters. The second umap run is performed with slighly different paramters for this reason
get_umap<-function(img.dat, outliers=T, thresh=0.15, min_dist=0, set_op_mix_ratio=0.25,ncenters=3,...) {
  
  if(outliers==T){
    test_umap<-NULL
    test_umap<-umap_run(img.dat = img.dat, 
                        min_dist=min_dist, 
                        set_op_mix_ratio=set_op_mix_ratio, 
                        ...)
    
    kmeans.result <- kmeans(test_umap$umap_out, centers=ncenters)
    centers <- kmeans.result$centers[kmeans.result$cluster, ] # "centers" is a data frame of 3 centers but the length of the dataset so we can canlculate distance difference easily.
    distances <- sqrt(rowSums((test_umap$umap_out - centers)^2))
    dist_ord<-distances[order(distances, decreasing=T)[1:50]]
    plot(dist_ord, main=runNames(img.dat))
    #thresh<-as.numeric(readline("Threshold for outlier removal? (input value)\n "))
    #look for the maximum difference in distances across the first 50 points
    max_diff<-min(diff(dist_ord))
    #calculate a threshold which is some percentage of the max distance
    percent_diff<-abs(max_diff)/max(distances)
    if(percent_diff>thresh) {
      outliers<- which(distances>=dist_ord[which(diff(dist_ord)==max_diff)])
      print(outliers)
      
    } else (outliers<-NULL)
    
    
    labels <- 1:nrow(test_umap$umap_out)
    
    pch <- rep(".", nrow(test_umap$umap_out))
    
    col <- rep("black", nrow(test_umap$umap_out))
    
    
    if(length(outliers)>0) {
      cat("Outliers in red being removed\n")
      cat("If this is not desired, increase the threshold value 'thresh'\n")
      labels[-outliers] <- "."
      pch[outliers] <- "+"
      col[outliers] <- "red"
      img.dat_clean<-img.dat[,-outliers]
      
    } else {
      img.dat_clean<-img.dat
    }
    
    pairs(test_umap$umap_out, pch=pch, col=col, main=runNames(img.dat))
    if(length(outliers>0)) {print(image(img.dat, col~x*y, main=paste(runNames(img.dat))))}
  } else if(outliers==F){
    img.dat_clean<-img.dat
    test_umap<-NULL
  }
  
  
  separation_umap<-umap_run(img.dat_clean, ...)
  return(list(clean_img.dat=img.dat_clean, umap_separation=separation_umap, outliers=outliers, umap_outliers=test_umap))
}

#function to run umap on imaging dataset
#hardcoded for 3 components, number of nearest neighbors changed with nn parameter
#choice of cosine metric is based on https://pubs.acs.org/doi/10.1021/acs.analchem.8b05827
#optimize output function with a, b, and min_dist, spread as shown here
#https://jlmelville.github.io/uwot/abparams.html
umap_run<-function(img.dat, nn=5, metric="cosine", n_components = 3, fastmap=FALSE, fm_r=1, fm_method="adaptive", fm_metric="average", fm_ncomp=3,...) {
  #library(uwot)
  #library(Cardinal)
  #library(scales)
  if(fastmap==TRUE){
    fm<-spatialFastmap(img.dat, r=fm_r, method=fm_method, metric=fm_metric, ncomp=fm_ncomp)
    print(summary(fm))
    fm2<-resultData(fm)[[1]] #Note, selecting only the first model
    fmdata<-fm2$scores
    if(sum(is.na(fmdata))>0){
      print("NA or NaN detected in fastmap output. Try 'gaussian' weights")
      return()
    }
    aa_no_ssc_umap_trim2<-uwot::umap((fmdata), metric = metric, n_neighbors = nn, n_components=n_components,...)
    
  } else {
    #browser()
    idata<-spectra(img.dat)
    rownames(idata)<-mz(img.dat)
    nn=nn
    aa_no_ssc_umap_trim2<-uwot::umap(t(as.matrix(idata)), metric = metric, n_neighbors = nn, n_components=n_components,...)
    
  }
  
  #plot(aa)
  colnames(aa_no_ssc_umap_trim2) <- c("x_umap", "y_umap", "z_umap")
  cols_no_ssc_umap_trim2<-as.data.frame(aa_no_ssc_umap_trim2) %>% dplyr::mutate( one=scales::rescale(x_umap),two=scales::rescale(y_umap), three=scales::rescale(z_umap)) %>% dplyr::select(4:6) %>% rgb()
  return(list(umap_out=aa_no_ssc_umap_trim2, color_scheme=cols_no_ssc_umap_trim2))
}


#function to run umap on spatial fastmap result 
#hardcoded for 3 components, number of nearest neighbors changed with nn parameter
#choice of cosine metric is based on https://pubs.acs.org/doi/10.1021/acs.analchem.8b05827
#optimize output function with a, b, and min_dist, spread as shown here
#https://jlmelville.github.io/uwot/abparams.html
umap_run_fm<-function(fm, model=1, img.dat, nn=5, metric="cosine", n_components = 3, ...) {
  #library(uwot)
  #library(Cardinal)
  #library(scales)
  fm2<-resultData(fm)[[model]]
  fmdata<-fm2$scores
  #rownames(fmdata)<-mz(img.dat)
  nn=nn
  aa_no_ssc_umap_trim2<-uwot::umap((fmdata), metric = metric, n_neighbors = nn, n_components=n_components,...)
  #plot(aa)
  colnames(aa_no_ssc_umap_trim2) <- c("x_umap", "y_umap", "z_umap")
  cols_no_ssc_umap_trim2<-as.data.frame(aa_no_ssc_umap_trim2) %>% mutate( one=scales::rescale(x_umap),two=scales::rescale(y_umap), three=scales::rescale(z_umap)) %>% select(4:6) %>% rgb()
  return(list(umap_out=aa_no_ssc_umap_trim2, color_scheme=cols_no_ssc_umap_trim2))
}


#function to find closest match between rgb color and palette
#palette defaults to R colors()
color.id<-function (col, palette=colors()) {
  c2 <- col2rgb(col)
  coltab <- col2rgb(palette)
  cdist <- apply(coltab, 2, function(z) sum((z - c2)^2))
  palette[which(cdist == min(cdist))]
}


#function to look for isolated pixels with no adjacent neighbors 
#and optionally remove them
#the n_thresh variable remove pixels with the specified number of neighbors or less
fix_pix<-function(dat, remove=F, r=1, n_thresh=1) {
  #dat<-brains_neg_selected[[1]]
  
  pdata <- pData(dat)
  
  neighbors=findNeighbors(pdata, r=as.numeric(r))
  
  sing_pix<-which(sapply(neighbors, length)<=n_thresh)
  cat("Pixels ", sing_pix, " have fewer than", n_thresh, " neighbors.\nUse the remove=T option to delete and return a modified dataset\n")
  
  if(remove==F){
    cat("nothing done")
    return()
  } else if (remove==T) {
    dat2<-dat %>% subsetPixels(!1:ncol(dat)%in%sing_pix  )
    cat("returning modified dataset with ",length(sing_pix)," pixels removed")
    return(dat2)
    print(dat2)
  } else cat("remove value doesn't make sense, must be T/F")
  
}


##read samples spotted plates####
#read sample list and format in a dataframe for other functions
#assumes that samples are ordered first by row, then column
#Eg if columns are 1,2,3..n, and rows are A,B,C..Z, will read
#A1, A2, A3,  etc
#can specify type as 'auto' based on first column in sample list, or manual to ignore first column

read_samples<-function(name, method=c('none', 'man', 'auto'), rownum=NULL, colnum=NULL) {
  #require(tidyr)
  #read list
  plate_dat<-read.csv(name, sep="\t")
  
  if(method=="man") {
    print("method is man")
    
    #create well IDs from rownum and column
    if(is.null(rownum)||is.null(colnum)){
      print("need rownum and colnum values")
      return()
    }
    
    id<-NULL
    for(i in 1:rownum){
      id<-c(id, paste(LETTERS[i], 1:colnum, sep=""))
    }
    id=data.frame(samples=id)
    
    id<-id %>%
      tidyr::separate(samples, 
               into = c("row_id", "col_id"), 
               sep = "(?<=[A-Za-z])(?=[0-9])"
      )
    
    #make sure length matches
    if(dim(id)[1]!=dim(plate_dat)[1]) {
      print("id length and sample list length do not match!")
      return()
    }
    
  } else if (method=="auto") {
    print("method is auto")
    
    id=data.frame(samples=as.character(plate_dat[,1]))
    
    id<-id %>%
      tidyr::separate(samples, 
               into = c("row_id", "col_id"), 
               sep = "(?<=[A-Za-z])(?=[0-9])"
      )
    
    
  } else {
    print("'method' not specified")
    return()
  }
  plate_dat<-cbind(id, plate_dat)
  return(plate_dat)
}


#Add phenotype data to multiple plates.
#Easiest way to use this function is to have a second column in the
#sample file with the plate string in addition to the regular spot ID in the first column.
#the plate string must be in the runname and unique for each run
pixDatFill_mult<-function(datas, sample_list, variables, method=c("spec_density"), inflect_thresh=7, lsp_plot=TRUE) {
  
  #require(Cardinal)
  
  #browser()
  
  plate_dat=sample_list
  data1_samples<-datas
  
  #clean up sample list to remove empty rows
  sample_list<-sample_list[rowSums(is.na(sample_list)) != ncol(sample_list),]
  
  #establish plates to work with
  plates=unique(plate_dat$Plate)
  #make sure there is data there!
 
  
  pdat<-NULL
  for(p in plates) {
    
    #check to make sure p is not NA or empty string, can happen with extra lines in data
    if(p %in% c("", NA)) {
      message("skipping plate with blank or NA value.")
      showNotification("skipping plate with blank or NA value. Check sample list.",type="warning")
      next
    }
    
    #create dataset of only selected plate
    cat("Working on plate: ", p,"\n")
    
    
    
    select_pix<-grep(p, Cardinal::run(data1_samples))
    if(length(select_pix)==0) {
      showNotification(paste0("Plate: ", p, " from sample list not found! skipping"))
      next
    } else {
      data1_plate<-subsetPixels(data1_samples, grep(p, Cardinal::run(data1_samples)))
      if(length(data1_plate)==0)
        next
      if(length(runNames(data1_plate))>1){
        cat("More than one plate selected for",p,", check naming \n")
        showNotification(paste0("More than one plate selected for",p," , check naming \n"))
        break
      }
    }
    
    #create data table for single plate
    plate_dat_single<-plate_dat[plate_dat$Plate %in% p,]
    
    #figure out row/column numbers from sample list
    nrows=length(unique(plate_dat_single$row_id))
    ncols=length(unique(plate_dat_single$col_id))
    
    
    d<-coord(data1_plate)
    
    #calculate x and y breaks and store in vector
    
    if(method=="breaks") {
      
      print("calculating breaks from spot positions")
      
      #function to calculate breaks within data
      calc_breaks<-function(x){
        a<-hist(x, breaks=min(x):max(x))
        
        #function to calc min/max from Evan Friedland
        #https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima/6836924
        inflect <- function(x, threshold = 1){
          up   <- sapply(1:threshold, function(n) c(x[-(seq(n))], rep(NA, n)))
          down <-  sapply(-1:-threshold, function(n) c(rep(NA,abs(n)), x[-seq(length(x), length(x) - abs(n) + 1)]))
          a    <- cbind(x,up,down)
          list(minima = which(apply(a, 1, min) == a[,1]), maxima = which(apply(a, 1, max) == a[,1]))
        }
        scatter.smooth(x=1:length(a$counts), a$counts, degree=2, family="gaussian", evaluation=length(a$counts), span=.2)
        aa<-loess.smooth(x=1:length(a$counts), a$counts, degree=2, family="gaussian", evaluation=length(a$counts), span=.2, col="blue")
        
        #inflect_thresh=7 #need this for debugging if not assigned earlier
        b<-inflect(aa$y, threshold = inflect_thresh)$minima
        if(length(b)==0) {
          idx<-c(0, max(x))
        } else {
          idx<-c(0,c(b[1],b[-1][diff(b)>1])+min(x), max(x))
        }
        abline(v=idx, col="red")
        return(idx)
      }
      idx<-calc_breaks(d$x)
      idy<-calc_breaks(d$y)
    } else if (method=="period") {
      
      
      #estimate spot size based on number of spots and pixel length, approx 50%
      period_est<-function(x) {
        
        a<-hist(x, breaks=min(x):max(x))
        
        #require(lomb)
        f.data<-lomb::lsp(a$counts,type='period',from=min(x),to=max(x),ofac=5, plot=lsp_plot)
        period=max(f.data$peak.at)
        
        if(f.data$p.value<0.05){
          print(cat("Note: LS p-value is not significant!\n"))
        }
        
        xc<-cos(2*pi*a$breaks[-length(a$breaks)]/period)
        xs<-sin(2*pi*a$breaks[-length(a$breaks)]/period)
        fit.lm <- lm(a$counts~xc+xs)
        
        # access the fitted series (for plotting)
        fit <- fitted(fit.lm)  
        
        # find predictions for original time series
        pred <- predict(fit.lm, newdata=data.frame(Time=a$breaks[-length(a$breaks)]))    
        
        plot(a$counts ~ a$breaks[-length(a$breaks)], xlim=c(1, max(a$breaks)), type="l")
        lines(fit, col="red")
        lines(a$breaks[-length(a$breaks)], pred, col="blue")
        
        #function to calc min/max from Evan Friedland
        #https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima/6836924
        inflect <- function(x, threshold = 1){
          up   <- sapply(1:threshold, function(n) c(x[-(seq(n))], rep(NA, n)))
          down <-  sapply(-1:-threshold, function(n) c(rep(NA,abs(n)), x[-seq(length(x), length(x) - abs(n) + 1)]))
          a    <- cbind(x,up,down)
          list(minima = which(apply(a, 1, min) == a[,1]), maxima = which(apply(a, 1, max) == a[,1]))
        }
        
        b<-inflect(pred, threshold = 1)$minima
        idx<-c(min(x),c(b[1],b[-1][diff(b)>1])+min(x), max(x))
        return(idx)
      }
      
      
      #if only one row, do not calculate breaks but use range
      if(dim(table(plate_dat_single$row_id))>1) {
        idy<-period_est(d$y)
      } else {
        idy<-range(d$y)
      }
      
      #if only one column, do not calculate breaks but use range
      if(dim(table(plate_dat_single$col_id))>1) {
        idx<-period_est(d$x)
      } else {
        idx<-range(d$x)
      }
      
      
    } else if(method=="manual") {
      #idx<- 
      
      print("manual method not implemented yet")
      
      
    }  else if (method=="spec_density") {
      period_est<-function(x) {
        
        a<-hist(x, breaks=min(x):max(x))
        
        #using this spec.ar function from Rob Hyndman
        #https://stats.stackexchange.com/questions/1207/period-detection-of-a-generic-time-series
        find.freq <- function(x)
        {
          n <- length(x)
          spec <- spec.ar(c(x),plot=FALSE)
          if(max(spec$spec)>10) # Arbitrary threshold chosen by trial and error.
          {
            period <- round(1/spec$freq[which.max(spec$spec)])
            if(period==Inf) # Find next local maximum
            {
              j <- which(diff(spec$spec)>0)
              if(length(j)>0)
              {
                nextmax <- j[1] + which.max(spec$spec[j[1]:500])
                period <- round(1/spec$freq[nextmax])
              }
              else
                period <- 1
            }
          }
          else
            period <- 1
          return(period)
        }
        
        
        period=find.freq(a$counts)
        
        if(period==1){
          print(cat("No period detected using spectra density, trying LS fit!!\n"))
          showNotification("No period detected using spectra density, trying LS fit!!", type="warning")
          
          f.data<-lomb::lsp(a$counts,type='period',from=min(x),to=max(x),ofac=5, plot=lsp_plot)
          period=max(f.data$peak.at)
          
          #browser()
          
        }
        
        xc<-cos(2*pi*a$breaks[-length(a$breaks)]/period)
        xs<-sin(2*pi*a$breaks[-length(a$breaks)]/period)
        fit.lm <- lm(a$counts~xc+xs)
        
        # access the fitted series (for plotting)
        fit <- fitted(fit.lm)  
        
        # find predictions for original time series
        pred <- predict(fit.lm, newdata=data.frame(Time=a$breaks[-length(a$breaks)]))    
        
        plot(a$counts ~ a$breaks[-length(a$breaks)], xlim=c(1, max(a$breaks)), type="l")
        lines(fit, col="red")
        lines(a$breaks[-length(a$breaks)], pred, col="blue")
        
        #function to calc min/max from Evan Friedland
        #https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima/6836924
        inflect <- function(x, threshold = 1){
          up   <- sapply(1:threshold, function(n) c(x[-(seq(n))], rep(NA, n)))
          down <-  sapply(-1:-threshold, function(n) c(rep(NA,abs(n)), x[-seq(length(x), length(x) - abs(n) + 1)]))
          a    <- cbind(x,up,down)
          list(minima = which(apply(a, 1, min) == a[,1]), maxima = which(apply(a, 1, max) == a[,1]))
        }
        
        b<-inflect(pred, threshold = 1)$minima
        idx<-c(min(x),c(b[1],b[-1][diff(b)>1])+min(x), max(x))
        return(idx)
      }
      
      #if only one row, do not calculate breaks but use range
      if(dim(table(plate_dat_single$row_id))>1) {
        idy<-period_est(d$y)
      } else {
        idy<-range(d$y)
      }
      
      #if only one column, do not calculate breaks but use range
      if(dim(table(plate_dat_single$col_id))>1) {
        idx<-period_est(d$x)
      } else {
        idx<-range(d$x)
      }
      
        
    }else {
      print(" method not recognized, use 'breaks', 'period' or 'manual'")
    }
    
    
    
    
    #give an message if the predicted index doesn't match the
    #number of rows/columns
    if(length(idx)!=ncols+1) {
      cat(runNames(data1_plate), "Dimension mismatch in spot columns between",
          "sample list and detected spots. Check for spurious spots.",
          "Consider using pix_fix() function or remove specific rows/columns.", 
          paste("idx= ",idx), paste("ncols=",ncols), sep="\n")
      
      cat("Median x size length = ", median(diff(idx)), "\n")
      bad_spot<-which(diff(idx)<(median(diff(idx))/2)+1)
      if(length(bad_spot)==1) {
        cat("Looks like row/column ", bad_spot, " is less than 50% of the median + 1, removing\n")
        idx<-idx[-c(bad_spot+1)]
      } else {
        print("\n")
        message("\nSample list mismatch with auto detected phenotype list. Please check blank spots or extra pixels and sample list.")
        message(cat("Plate is:", p, " and detected x breaks at: ", idx, " with y breaks at: ", idy))
        
        showNotification("Sample list mismatch with auto detected phenotype list. Please check blank spots or extra pixels and sample list. Check console for details", duration=10, type="error")
        
        
        #browser()
        return()
        
      }
      
    }
    
    if(length(idy)!=nrows+1) {
      #browser()
      
      cat(runNames(data1_plate), "Dimension mismatch in spot rows between",
          "sample list and detected spots. Check for spurious spots.",
          "Consider using pix_fix() function or removing specific rows/columns.", 
          paste("idy= ",idy), paste("nrows=",nrows), sep="\n")
      break_y_lengths=diff(idy)
      cat("Median y size length = ", median(break_y_lengths),"\n")
      bad_spot<-which(break_y_lengths<median(break_y_lengths)/2)
      if(length(bad_spot)==1) {
        cat("Looks like row/column ", bad_spot, " is less than 50% of the median, removing\n")
        idy<-idy[-c(bad_spot+1)]
      } else if(sum(break_y_lengths[bad_spot]<(median(break_y_lengths)/2))==length(bad_spot)) { #bit of a hack-- removing breaks less
        cat("Looks like detected break ", bad_spot, " is less than 50% of the median, removing\n")
        idy<-idy[-c(bad_spot+1)]
      } else {
        print("\n")
        message("\nSample list mismatch with auto detected phenotype list. Please check blank spots or extra pixels and sample list.")
        message(cat("Plate is:", p, " and detected x breaks at: ", idx, " with y breaks at: ", idy))
        
        showNotification("Sample list mismatch with auto detected phenotype list. Please check blank spots or extra pixels and sample list. Check console for details", duration=10, type="error")
        
        browser()
        return()
      }
        
      
      
    }
    
    #browser()
    
    
    #for each variable in the variable list, lets loop through
    for(v in 1:length(variables)) {
      #lets write some loops to provide the group information
      
      
      for(i in 1:nrows){
        
        for(j in 1:ncols) {
          
          #DEBUG
          #print(paste("v = ", v))
          #print(paste("i = ", i))
          #print(paste("j = ", j))

          #test of zero length and skip instead of failing
          if (length(as.character(
            dplyr::filter(
              plate_dat_single,
              row_id == unique(plate_dat_single$row_id)[i],
              col_id == unique(plate_dat_single$col_id)[j]
            )[, variables[v]]
          )) == 0)
            next
          
          #extract value to use for spot to be applied to pixel data based on row and column ID
          pdat_val = as.character(dplyr::filter(
            plate_dat_single,
            row_id == unique(plate_dat_single$row_id)[i],
            col_id == unique(plate_dat_single$col_id)[j]
          )[, variables[v]])
          
          try(pixelData(data1_plate)[[variables[v]]][coord(data1_plate)$x %in% (idx[j]):(idx[j +
                                                                                               1]) &
                                                       coord(data1_plate)$y %in% (idy[i]):(idy[i +
                                                                                                 1])] <- pdat_val
              
          )
          
          #if(is.na(pdat_val)) next
        }
      }
      
    }
    if(is.null(pdat)){
      pdat<-pixelData(data1_plate)
    } else {
      pdat<-rbind(pdat, pixelData(data1_plate))      
    }
    
    
  }
  if(is.null(pdat)){
    print("No new pdata created! Ensure the Plate column of the sample list is matched in run data.")
    return(pixelData(data1_samples))
  } else {
    return(pdat)
  }
}

#when the tissue / spot positions are not regular, this function will assign the requested values
#the sample_list information must contain xstart, xend, ystart, and yend columns in addition to other sample info
pixDatFill_manual<-function(datas, sample_list, variables) {
  
  
  plate_dat=sample_list
  data1_samples<-datas
  
  #establish plates to work with
  plates=unique(plate_dat$Plate)
  
  
  
  pdat<-NULL
  for(p in plates) {
    #create dataset of only selected plate
    data1_plate<-subsetPixels(data1_samples, grep(p, run(data1_samples)))
    if(length(data1_plate)==0)
      next
    
    #create data table for single plate
    plate_dat_single<-plate_dat %>% dplyr::filter(Plate==p)
    
    
    
    
    #for each variable in the variable list, lets loop through
    for(v in 1:length(variables)) {
      #lets write some loops to provide the group information
      #because each spot is defined, simply loop through spots
      
      for(i in 1:nrow(plate_dat_single)) {       
        
        #extract value to use for spot to be applied to pixel data
        pdat_val=as.character(plate_dat_single[i,variables[v]])
        
        xstart=plate_dat_single[i,"xstart"]
        xend=plate_dat_single[i,"xend"]
        ystart=plate_dat_single[i,"ystart"]
        yend=plate_dat_single[i,"yend"]
        
        try(
          pixelData(data1_plate)[[variables[v]]][coord(data1_plate)$x %in% (xstart:xend) & 
                                                   coord(data1_plate)$y %in% (ystart:yend)] <- pdat_val
        )
        
        #if(is.na(pdat_val)) next
      }
      
      
    }
    if(is.null(pdat)){
      pdat<-pixelData(data1_plate)
    } else {
      pdat<-rbind(pdat, pixelData(data1_plate))      
    }
    
    
  }
  if(is.null(pdat)){
    print("No new pdata created! Ensure the Plate column of the sample list is matched in run data.")
    return(pixelData(data1_samples))
  } else {
    #redorder rows of pdat based on original order
    #browser()
    
    
    
    # Step 1: Create an index in 'original_df' to preserve the original order
    original_df <- pData(datas) %>% as.data.frame() %>%
      dplyr::mutate(original_order = dplyr::row_number())
    
    
    #Step 2: Use left_join() to join 'new_df' to 'original_df' and reorder based on original_order
    reordered_df <- as.data.frame(pdat) %>%
      dplyr::left_join(original_df %>% dplyr::select(x, y, run, original_order), by = c("x", "y", "run")) %>%
      dplyr::arrange(original_order) %>%
      dplyr::select(-original_order)  # Drop the original_order column if not needed
    

    
    pixelData(data1_samples)<-PositionDataFrame(run=reordered_df$run, coord=coord(datas), reordered_df[,!colnames(reordered_df)%in%c("run", "x", "y")]) 
    
    return(pixelData(data1_samples))
  }
}


#for Cardinal 3.6+ the combine function doesn't work the same way
#in many case can use do.call(cbind) to combine from a list, but doesn't work with some datasets
#this function is a workaround to combine datasets

combine_card <- function(x) {
  
  if (length(x) == 1) {
    return(x[[1]])
  }
  if (length(x) == 0) {
    return(NULL)
  }
  if (length(x) > 1) {
    
    #get run order to preserve order of runs
    run_order <- unlist(lapply(x, runNames))
    
    #reorder list based on ncol
    x <- x[order(unlist(lapply(x, ncol)), decreasing = TRUE)]
    
    tmp <- cbind(x[[1]])
    for (i in 2:length(x)) {
      tmp <- cbind(tmp, x[[i]])
    }
    
    #reorder pixels based on orginal order
    pData(tmp)$run <- factor(run(tmp), levels = run_order)
    
    return(tmp)
  }
}

#convert old cardinal object to new one
convert_card<-function(obj) {
  
  #obj<-readRDS("_Neg_761_oxid_targeted_MSI-depth_picked_SN-peak_picked-2022-10-19.rds")
  spectra <- obj@imageData$data[[1L]]
  
  coord <- obj@elementMetadata@coord
  run <- obj@elementMetadata@run
  pdat2 <- obj@elementMetadata@listData
  
  mz <- obj@featureData@mz
  ID <- obj@featureData@listData[[1]]
  
  #create new cardinal object
  
  pdata <- PositionDataFrame(run=run, coord=coord, pdat2)
  fdata <- MassDataFrame(mz=mz, ID=ID)
  
  out <- MSImagingExperiment(spectraData=spectra,
                             featureData=fdata,
                             pixelData=pdata)
  
  
  
}

#round to higher integer, ie 2.5 becomes 3, not 2 as in round()
round2 = function(x, digits) {
  posneg = sign(x)
  z = abs(x)*10^digits
  z = z + 0.5 + sqrt(.Machine$double.eps)
  z = trunc(z)
  z = z/10^digits
  z*posneg
}
