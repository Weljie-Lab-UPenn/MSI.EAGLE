#require(Cardinal)
library(magrittr)

#### Centroided spectra ####
# Centroided data (e.g. Waters .raw centroided and lock-mass corrected with AFAMM) must be aligned directly. Cardinal
# only flags a file as centroided when the imzML fileContent lists "centroid spectrum" (MS:1000127); AFAMM writes it
# per spectrum instead, so the data load as unflagged. Unflagged processed data are resampled onto a profile m/z grid
# (readImzML with mass.range/resolution, or convertMSImagingArrays2Experiment) and re-picked, which moves every peak
# 1-2 grid steps (~20-38 ppm with the default ~17.5 ppm grid) below its true m/z.

# spectrum type declared in an imzML file ("centroid", "profile" or NA), from the header and the first spectra
imzml_spectrum_type <- function(path, n_bytes = 4e6) {
  if (is.null(path) || !file.exists(path) || !grepl("\\.imzml$", path, ignore.case = TRUE)) return(NA_character_)
  con <- file(path, "rb"); on.exit(close(con))
  txt <- readChar(con, n_bytes, useBytes = TRUE)
  centroid <- grepl("MS:1000127", txt, fixed = TRUE)
  profile <- grepl("MS:1000128", txt, fixed = TRUE)
  if (centroid && !profile) "centroid" else if (profile && !centroid) "profile" else NA_character_
}

is_centroided_msi <- function(x) isTRUE(tryCatch(Cardinal::isCentroided(x), error = function(e) FALSE))

all_centroided <- function(lst) length(lst) > 0L && all(vapply(lst, is_centroided_msi, logical(1)))

mark_centroided <- function(x) {
  Cardinal::centroided(x) <- TRUE
  x
}

# combine raw files for peak processing. Centroided: kept as centroids (MSImagingArrays, flag re-applied); profile:
# converted to an MSImagingExperiment on a profile m/z grid, as before
prepare_for_peaks <- function(raw_list, mass.range = NULL, resolution = NA, units = "ppm") {
  arr <- Cardinal::combine(lapply(raw_list, convertMSImagingExperiment2Arrays))
  if (all_centroided(raw_list)) {
    message("[Peaks] centroided spectra: kept as centroids (no profile resampling, no re-picking)")
    return(mark_centroided(arr))
  }
  convertMSImagingArrays2Experiment(arr, mass.range = mass.range, resolution = resolution, units = units)
}

# Reference peak list for centroided data. Aligning single-pixel centroids directly splits weak ions over several
# low-frequency bins and has no noise filter, so instead: the centroids of (a sample of) all spectra are summed into a
# fine log-spaced histogram (res_ppm), lightly smoothed, and peaks are picked on that pooled spectrum with the usual
# SNR / noise method. Each peak's m/z is the intensity-weighted mean of the centroids within its bounds (at most
# +/- tolerance/2), so no m/z is snapped to a coarse grid.
# Pooling follows the Data Setup parallel settings: the backend from getCardinalBPPARAM() (Serial, Multicore, Snow or
# Future) and the number of chunks from getCardinalNChunks(). Chunks run in waves of one chunk per worker (bplapply,
# which every backend supports); each chunk's sums are added as soon as its wave returns, with a progress message.

# one chunk of spectra -> non-empty log-spaced bins: index k, summed intensity s, summed intensity x m/z sm.
# 'spec' is either a list(mz, intensity) holding just this chunk's spectra (in-memory data), or the full file-backed
# matter lists plus the chunk's indices (only a reference is sent to the worker, which reads its own spectra)
.pool_centroid_chunk <- function(spec, lo, step, nb) {
  if (is.null(spec$idx)) { m <- unlist(spec$mz, use.names = FALSE); v <- unlist(spec$intensity, use.names = FALSE) }
  else {
    suppressPackageStartupMessages(requireNamespace("matter", quietly = TRUE))   # fresh Snow / Future workers
    m <- unlist(as.list(spec$mz[spec$idx]), use.names = FALSE); v <- unlist(as.list(spec$intensity[spec$idx]), use.names = FALSE)
  }
  b <- floor(log(m / lo) / step) + 1L; ok <- b >= 1L & b <= nb & v > 0
  if (!any(ok)) return(list(k = integer(0), s = numeric(0), sm = numeric(0)))
  r <- rowsum(cbind(v[ok], v[ok] * m[ok]), b[ok])
  list(k = as.integer(rownames(r)), s = r[, 1], sm = r[, 2])     # sparse: only non-empty bins go back to the main process
}

centroid_reference_peaks <- function(x, SNR = 20, method = "diff", tolerance = 30, res_ppm = 1, smooth_ppm = 3,
                                     sampleSize = NA, mass.range = NULL, nchunks = NULL,
                                     BPPARAM = Cardinal::getCardinalBPPARAM()) {
  t0 <- proc.time()[["elapsed"]]
  mzl <- Cardinal::mz(x); inl <- Cardinal::intensity(x); n <- length(mzl)
  idx <- if (is.na(sampleSize)) seq_len(n) else
    unique(round(seq(1, n, length.out = if (sampleSize < 1) max(1, ceiling(sampleSize * n)) else min(n, sampleSize))))
  if (is.null(mass.range) || length(mass.range) < 2L || any(!is.finite(mass.range))) {
    probe <- idx[unique(round(seq(1, length(idx), length.out = min(500L, length(idx)))))]
    mass.range <- range(unlist(lapply(mzl[probe], range)))
  }
  lo <- min(mass.range); hi <- max(mass.range); step <- log1p(res_ppm * 1e-6)
  nb <- ceiling(log(hi / lo) / step) + 1L
  if (is.null(nchunks)) nchunks <- Cardinal::getCardinalNChunks()             # Data Setup "# of chunks"
  nchunks <- max(1L, min(as.integer(nchunks), length(idx)))
  blocks <- unname(split(idx, cut(seq_along(idx), nchunks, labels = FALSE)))
  if (is.null(BPPARAM)) BPPARAM <- BiocParallel::SerialParam()                 # Data Setup "Parallelization mode"
  nw <- max(1L, tryCatch(as.integer(BiocParallel::bpnworkers(BPPARAM)), error = function(e) 1L))
  file_backed <- inherits(mzl, "matter") && inherits(inl, "matter")
  message(sprintf("[Peaks] pooled centroid spectrum: %d spectra in %d chunks of ~%d (%s, %d worker%s), m/z %.1f-%.1f at %g ppm",
                  length(idx), nchunks, ceiling(length(idx) / nchunks), class(BPPARAM)[1], nw, if (nw == 1L) "" else "s",
                  lo, hi, res_ppm))
  S <- numeric(nb); SM <- numeric(nb); done <- 0L               # summed intensity and intensity x m/z per bin
  if (!BiocParallel::bpisup(BPPARAM)) {                         # start the workers once for all waves
    BiocParallel::bpstart(BPPARAM); on.exit(BiocParallel::bpstop(BPPARAM), add = TRUE)
  }
  for (wave in split(seq_along(blocks), ceiling(seq_along(blocks) / nw))) {
    spec <- lapply(blocks[wave], function(b) if (file_backed) list(mz = mzl, intensity = inl, idx = b)
                                             else list(mz = mzl[b], intensity = inl[b]))
    res <- BiocParallel::bplapply(spec, .pool_centroid_chunk, lo = lo, step = step, nb = nb, BPPARAM = BPPARAM)
    for (r in res) {
      S[r$k] <- S[r$k] + r$s; SM[r$k] <- SM[r$k] + r$sm; done <- done + 1L
      el <- proc.time()[["elapsed"]] - t0
      message(sprintf("[Peaks] pooled centroid spectrum: chunk %d/%d done (%.0f s elapsed, ~%.0f s left)",
                      done, nchunks, el, el / done * (nchunks - done)))
    }
  }
  sd_bins <- smooth_ppm / res_ppm; kern <- stats::dnorm(seq(-ceiling(3 * sd_bins), ceiling(3 * sd_bins)), sd = sd_bins)
  Ss <- as.numeric(stats::filter(S, kern / sum(kern), sides = 2)); Ss[is.na(Ss)] <- 0
  noise <- switch(method, quantile = "quant", cwt = , adaptive = "diff", method)
  pk <- matter::findpeaks(Ss, snr = SNR, noise = noise, relheight = NULL)   # as Cardinal peakPick (no relative-height cut)
  if (!length(pk)) return(numeric(0))
  lb <- attr(pk, "left_bounds"); rb <- attr(pk, "right_bounds"); half <- ceiling(0.5 * tolerance / res_ppm)
  ref <- vapply(seq_along(pk), function(j) {
    a <- max(lb[j], pk[j] - half, 1L); z <- min(rb[j], pk[j] + half, nb); sum(SM[a:z]) / sum(S[a:z])
  }, numeric(1))
  ref <- sort(unique(ref[is.finite(ref)]))
  message(sprintf("[Peaks] centroided spectra: %d reference peaks from the pooled centroid spectrum of %d spectra (SNR %s; %.0f s)",
                  length(ref), length(idx), SNR, proc.time()[["elapsed"]] - t0))
  ref
}

# Spectrum type of the loaded raw data ("centroid", "profile" or NA), for UI defaults. A reactiveVal only invalidates
# when its value changes, so panels that use it re-render when newly imported data change type, not on every Data
# Setup edit
spectrum_mode_tracker <- function(setup_values) {
  mode <- shiny::reactiveVal(NA_character_)
  shiny::observe({
    rl <- setup_values()[["x1"]]$raw_list
    mode(if (is.null(rl) || !length(rl)) NA_character_ else if (all_centroided(rl)) "centroid" else "profile")
  })
  mode
}

# S/N box that switches with the spectrum type: for centroided data it sets the S/N on the pooled centroid spectrum
# (ref_SNR, default 3); for profile data it is the per-spectrum peak-picking S/N with the tab's usual default
snr_input <- function(id, mode, label, profile_default) {
  cen <- identical(mode, "centroid")
  shiny::numericInput(id, if (cen) paste(label, "- centroided data: S/N on the pooled centroid spectrum") else label,
                      if (cen) 3 else profile_default, min = 0)
}

# keep features inside a mass range (applied after alignment for centroided data, which are not resampled)
restrict_mass_range <- function(x, mass.range = NULL) {
  if (is.null(mass.range) || length(mass.range) < 2L || any(!is.finite(mass.range))) return(x)
  subsetFeatures(x, mz >= min(mass.range) & mz <= max(mass.range))
}

#function to process data
# ref_SNR: SNR for picking reference peaks on the pooled centroid spectrum (centroided data only; SN is the per-spectrum
# SNR used for profile data; in the app both come from the same S/N box, which defaults to 3 for centroided data). On
# the 2 x 2 um DESI test data, ref_SNR = 3 recovered 93 % of the features that the
# profile route (SN = 20) found in >= 5 % of pixels.
HTS_reproc<-function(data1, SN, res=10, align_tol=15, pix_include=NULL,mzrange=NULL, freq.min=0.01, method="adaptive", ref_SNR=3) {
  
  #remove z-axis as it messes with the segmentation algorithms
  coord(data1)$z<-NULL
  
  
  #remove non-relevant peaks for reprocessing-- does this make a difference?
  if(!is.null(pix_include)) {
    data1<-data1[coord(data1)%in%pix_include]
  }
  
  #setCardinalBPPARAM(MulticoreParam())
  
  
  if (is_centroided_msi(data1)) {
    # centroided input: no profile resampling and no re-picking. Peaks are picked once on the pooled centroid spectrum
    # (ref_SNR applies there) and every pixel's centroid intensities are extracted at those reference m/z
    message("[HTS_reproc] centroided spectra: reference peaks from the pooled centroids, then per-pixel extraction")
    ref <- centroid_reference_peaks(data1, SNR = ref_SNR, method = method, tolerance = align_tol, mass.range = mzrange)
    t_ex <- proc.time()[["elapsed"]]
    message(sprintf("[HTS_reproc] extracting intensities at %d reference peaks from %d spectra (Cardinal chunks follow)",
                    length(ref), length(data1)))
    data3 <- data1 %>% normalize(method="tic") %>%
      peakProcess(ref = ref, tolerance = align_tol, units = "ppm", verbose = TRUE) %>% process(verbose = TRUE)
    message(sprintf("[HTS_reproc] extraction done (%.0f s)", proc.time()[["elapsed"]] - t_ex))
    # peakProcess(ref = ...) does not count detections, so add the per-feature pixel frequency used by the filter below
    nz <- Cardinal::summarizeFeatures(data3, stat = "nnzero")
    fData(data3)$count <- fData(nz)[[grep("nnzero", names(fData(nz)))[1]]]
    fData(data3)$freq <- fData(data3)$count / ncol(data3)
  } else {
  data2 <- data1 %>% #reduceBaseline(method="median", blocks=100) %>%
    normalize(method="tic") %>% 
    peakPick(
      #pixel=seq(1, ncol(data1), by=1), 
      method=method, SNR=SN) %>%
    process()

  #peakAlign no longer available in Cardinal 3.6, but peakAlign adds info
  data3<- data2 %>% peakAlign(units="ppm", tolerance=align_tol) %>%
    process()
  }
  
  
  #if no peaks meet filtering criteria, return original data3
  if(sum(fData(data3)$freq>freq.min)>0) {
    data4 <- subsetFeatures(data3, freq> freq.min)
  } else {
    message("No peaks meet filtering criteria, returning peak aligned data")
    data4 <- data3
  }

  # centroided data are not resampled to a mass range, so restrict the aligned features instead
  data4 <- restrict_mass_range(data4, mzrange)

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
umap_run<-function(img.dat, nn=5, metric="cosine", n_components = 3, fastmap=FALSE, fm_r=1, fm_method="adaptive", fm_metric="euclidean", fm_ncomp=3,...) {
  #library(uwot)
  #library(Cardinal)
  #library(scales)
  if(fastmap==TRUE){
    fm_r <- suppressWarnings(as.numeric(fm_r))
    if (!is.finite(fm_r) || fm_r <= 0) {
      fm_r <- 1
    }
    fm_weights <- tolower(as.character(fm_method)[1])
    if (!fm_weights %in% c("gaussian", "adaptive")) {
      fm_weights <- "adaptive"
    }
    fm_metric <- tolower(as.character(fm_metric)[1])
    valid_fm_metrics <- c("euclidean", "manhattan", "cosine", "correlation")
    if (!fm_metric %in% valid_fm_metrics) {
      warning(sprintf("Invalid fastmap metric '%s'; using 'euclidean'.", fm_metric))
      fm_metric <- "euclidean"
    }
    
    # Force coordinate storage to double before neighbor search to avoid
    # integer-only kdsearch failures in some Cardinal builds.
    coord_df <- try(as.data.frame(Cardinal::coord(img.dat)), silent = TRUE)
    if (inherits(coord_df, "try-error") || nrow(coord_df) == 0) {
      stop("Could not extract coordinates for spatialFastmap.")
    }
    coord_num <- coord_df[vapply(coord_df, is.numeric, logical(1))]
    if (ncol(coord_num) < 2) {
      stop("spatialFastmap requires at least two numeric coordinate columns.")
    }
    coord_mat <- as.matrix(coord_num)
    storage.mode(coord_mat) <- "double"
    
    groups <- try(as.character(Cardinal::run(img.dat)), silent = TRUE)
    if (inherits(groups, "try-error") || length(groups) != nrow(coord_mat)) {
      groups <- NULL
    }
    
    neighbors <- try(
      Cardinal::findNeighbors(coord_mat, r = fm_r, groups = groups),
      silent = TRUE
    )
    if (inherits(neighbors, "try-error")) {
      stop(sprintf("neighbor search failed: %s", as.character(neighbors)))
    }
    
    fm <- try(
      spatialFastmap(
        img.dat,
        r = fm_r,
        neighbors = neighbors,
        weights = fm_weights,
        metric = fm_metric,
        ncomp = fm_ncomp
      ),
      silent = TRUE
    )
    if (inherits(fm, "try-error")) {
      stop(sprintf("spatialFastmap failed: %s", as.character(fm)))
    }
    
    print(summary(fm))
    fmdata <- NULL
    
    if (inherits(fm, "fastmap")) {
      fmdata <- fm$x
    } else if (inherits(fm, "SpatialFastmap")) {
      fm_x <- try(fm$x, silent = TRUE)
      if (!inherits(fm_x, "try-error") && !is.null(fm_x)) {
        fmdata <- fm_x
      }
      if (is.null(fmdata)) {
        fm_model <- try(modelData(fm), silent = TRUE)
        if (!inherits(fm_model, "try-error") && is.list(fm_model) && !is.null(fm_model$x)) {
          fmdata <- fm_model$x
        }
      }
      if (is.null(fmdata)) {
        fm_pred <- try(stats::predict(fm, newdata = img.dat), silent = TRUE)
        if (!inherits(fm_pred, "try-error") && !is.null(fm_pred)) {
          fmdata <- fm_pred
        }
      }
    } else {
      fm_result <- try(resultData(fm), silent = TRUE)
      if (!inherits(fm_result, "try-error") && length(fm_result) >= 1) {
        fm2 <- fm_result[[1]]
        if (!is.null(fm2$scores)) {
          fmdata <- fm2$scores
        }
      }
    }
    
    if (is.null(fmdata)) {
      stop(sprintf("spatialFastmap produced no usable scores object (class=%s).", paste(class(fm), collapse = ",")))
    }
    fmdata <- as.matrix(fmdata)
    if (nrow(fmdata) < 2 || ncol(fmdata) < 2) {
      stop(sprintf("spatialFastmap scores have invalid shape: %s x %s", nrow(fmdata), ncol(fmdata)))
    }
    if (any(!is.finite(fmdata))) {
      stop("NA/NaN/Inf detected in fastmap output. Try 'gaussian' weights or lower fastmap r.")
    }
    aa_no_ssc_umap_trim2 <- uwot::umap((fmdata), metric = metric, n_neighbors = nn, n_components=n_components,...)
    
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
    
    #browser()
    
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

  # Preserve pixel identity across subset/rbind operations. Coordinates are not
  # guaranteed to be unique in stitched or polygon-derived datasets.
  pdat_cols <- colnames(as.data.frame(pixelData(data1_samples)))
  pixel_order_col <- ".msi_eagle_pixel_order"
  while (pixel_order_col %in% pdat_cols) {
    pixel_order_col <- paste0(pixel_order_col, "_")
  }
  pixelData(data1_samples)[[pixel_order_col]] <- seq_len(nrow(as.data.frame(pixelData(data1_samples))))
  for (v in variables) {
    if (!v %in% colnames(as.data.frame(pixelData(data1_samples)))) {
      pixelData(data1_samples)[[v]] <- NA
    }
  }
  
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
    return(pixelData(datas))
  } else {
    original_df <- as.data.frame(pixelData(data1_samples))
    reordered_df <- original_df
    pdat_df <- as.data.frame(pdat)

    if (!pixel_order_col %in% colnames(pdat_df)) {
      message("Could not preserve pixel order during manual phenotyping.")
      showNotification("Could not preserve pixel order during manual phenotyping.", type = "error")
      return(pixelData(datas))
    }

    pdat_df <- pdat_df[order(pdat_df[[pixel_order_col]]), , drop = FALSE]
    if (any(duplicated(pdat_df[[pixel_order_col]]))) {
      message("Duplicate pixel assignments detected during manual phenotyping; keeping the first assignment for each pixel.")
      showNotification("Duplicate pixel assignments detected during manual phenotyping; keeping the first assignment for each pixel.", type = "warning")
      pdat_df <- pdat_df[!duplicated(pdat_df[[pixel_order_col]]), , drop = FALSE]
    }

    update_idx <- as.integer(pdat_df[[pixel_order_col]])
    valid_idx <- !is.na(update_idx) & update_idx >= 1L & update_idx <= nrow(reordered_df)
    if (!all(valid_idx)) {
      message("Some manual phenotype pixel assignments were outside the source pixel range and were skipped.")
      showNotification("Some manual phenotype pixel assignments were outside the source pixel range and were skipped.", type = "warning")
      pdat_df <- pdat_df[valid_idx, , drop = FALSE]
      update_idx <- update_idx[valid_idx]
    }

    reordered_df[update_idx, variables] <- pdat_df[, variables, drop = FALSE]
    reordered_df <- reordered_df[, !colnames(reordered_df) %in% pixel_order_col, drop = FALSE]

    pixelData(data1_samples)<-PositionDataFrame(run=reordered_df$run, coord=coord(datas), reordered_df[,!colnames(reordered_df)%in%c("run", "x", "y")]) 
    
    return(pixelData(data1_samples))
  }
}


#for Cardinal 3.6+ the combine function doesn't work the same way
#in many case can use do.call(cbind) to combine from a list, but doesn't work with some datasets
#this function is a workaround to combine datasets

combine_card <- function(x) {
  
  # Drop empty entries to avoid method dispatch errors during combine.
  x <- Filter(Negate(is.null), x)
  
  if (length(x) == 1) {
    return(x[[1]])
  }
  if (length(x) == 0) {
    return(NULL)
  }
  if (length(x) > 1) {
    msi_check_feature_compatibility(x)
    
    #get run order to preserve order of runs
    run_order <- unique(as.character(unlist(lapply(x, runNames))))
    
    #reorder list based on ncol
    x <- x[order(unlist(lapply(x, ncol)), decreasing = TRUE)]
    
    # Avoid one-object cbind() since it can dispatch through dplyr::combine()
    # for certain metadata payloads.
    tmp <- x[[1]]
    for (i in 2:length(x)) {
      tmp <- cbind(tmp, x[[i]])
    }
    
    #reorder pixels based on orginal order
    run_vals <- as.character(run(tmp))
    run_levels <- unique(c(run_order, run_vals))
    pData(tmp)$run <- factor(run_vals, levels = run_levels)
    
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

# Pixel-level means test function
# Performs t-test or ANOVA at pixel level when pixel coordinates are used as grouping
pixel_level_means_test <- function(data, test_var, grouping_var = NULL) {
  
  # Extract spectra and metadata
  spectra_mat <- as.matrix(spectra(data))
  pdat <- as.data.frame(pData(data))
  mz_vals <- mz(data)
  feature_ids <- features(data)
  
  # Get test variable factor
  test_groups <- factor(pdat[[test_var]])
  
  # Check if this is a pixel-level comparison.
  # `grouping_var` can be a vector of field names (e.g. c("histo_cluster","x","y")).
  grouping_chr <- as.character(grouping_var)
  grouping_chr <- grouping_chr[!is.na(grouping_chr) & nzchar(grouping_chr)]
  is_pixel_level <- FALSE
  if (length(grouping_chr) == 0) {
    is_pixel_level <- TRUE
  } else if (length(grouping_chr) == 1 && grouping_chr %in% c("none", "1")) {
    is_pixel_level <- TRUE
  } else {
    is_pixel_level <- all(c("x", "y") %in% grouping_chr) || any(grouping_chr %in% c("x.y", "x_y"))
  }
  
  if (!is_pixel_level) {
    stop("This function is for pixel-level comparisons only. Use regular meansTest for biological replicates.")
  }
  
  message("Performing pixel-level means test - each pixel treated as independent replicate")
  
  # Perform t-test for each m/z
  n_features <- nrow(spectra_mat)
  results_list <- vector("list", n_features)
  
  for (i in 1:n_features) {
    intensities <- spectra_mat[i, ]
    
    # Remove NA values
    valid_idx <- !is.na(intensities)
    y <- intensities[valid_idx]
    groups <- test_groups[valid_idx]
    
    # Perform test based on number of groups
    if (nlevels(droplevels(groups)) == 2) {
      # Two-sample t-test
      test_result <- t.test(y ~ groups)
      statistic <- test_result$statistic
      pvalue <- test_result$p.value
      
    } else {
      # ANOVA for >2 groups
      test_result <- oneway.test(y ~ groups)
      statistic <- test_result$statistic
      pvalue <- test_result$p.value
    }
    
    # Store results
    results_list[[i]] <- list(
      i = i,
      mz = mz_vals[i],
      feature = feature_ids[i],
      statistic = statistic,
      pvalue = pvalue
    )
  }
  
  # Convert to data frame matching topFeatures format
  stats_df <- data.frame(
    i = sapply(results_list, `[[`, "i"),
    mz = sapply(results_list, `[[`, "mz"),
    feature = sapply(results_list, `[[`, "feature"),
    statistic = sapply(results_list, `[[`, "statistic"),
    pvalue = sapply(results_list, `[[`, "pvalue"),
    stringsAsFactors = FALSE
  )
  
  # Calculate FDR
  stats_df$fdr <- p.adjust(stats_df$pvalue, method = "BH")
  stats_df <- stats_df[, c("i", "mz", "feature", "statistic", "pvalue", "fdr")]
  
  # Sort by p-value
  stats_df <- stats_df[order(stats_df$pvalue), ]
  
  return(stats_df)
}

# --- Collapsible options panel -------------------------------------------------
# Every analysis tab is a sidebarLayout: options on the left, results on the
# right. msi_options_toggle() emits a small button that hides the options column
# and lets the results column expand to the full window width; clicking again
# restores it. The toggle is done entirely in the browser (no Shiny round trip,
# no server code per tab), by finding the sidebarLayout row that follows the
# button and swapping the Bootstrap column classes.
msi_ui_assets <- function() {
  shiny::tags$head(
    shiny::tags$style(shiny::HTML("
      .msie-toggle-bar { margin: 4px 0 6px 0; }
      .msie-toggle-btn { font-size: 12px; }

      /* Plot panels: one bordered card per plot, with a title bar, so the ion
         image controls and the mass spectrum controls are visibly separate. */
      .msie-panel {
        border: 1px solid #dcd6cd;
        border-radius: 4px;
        background: #fff;
        padding: 10px 12px 6px 12px;
        margin: 0 0 14px 0;
      }
      .msie-panel-title {
        font-size: 12px;
        font-weight: 700;
        text-transform: uppercase;
        letter-spacing: 0.05em;
        color: #8a7f72;
        border-bottom: 1px solid #eee7dd;
        padding-bottom: 6px;
        margin-bottom: 10px;
      }
      .msie-subtitle {
        font-size: 12px;
        font-weight: 600;
        color: #8a7f72;
        margin: 4px 0 6px 0;
      }
      .msie-gate {
        background: #f6f3ee;
        border: 1px solid #e6e1da;
        border-radius: 4px;
        padding: 8px 8px 0 8px;
        margin-bottom: 10px;
      }
      .msie-opts .form-group { margin-bottom: 8px; }
      .msie-opts .control-label { font-size: 12px; margin-bottom: 2px; }
      .msie-opts .checkbox { margin-top: 2px; margin-bottom: 2px; }
      .msie-opts .nav-pills > li > a { padding: 4px 10px; font-size: 12px; }
      .msie-opt-body {
        padding-top: 10px;
        max-height: 72vh;
        overflow-y: auto;
        overflow-x: hidden;
      }
      /* While the server is busy for more than half a second, dim the tab
         content as a visual cue. The content stays interactive, so the user can
         keep selecting and clicking; Shiny processes those inputs once R is
         free. */
      body.msie-busy .tab-content {
        opacity: 0.40;
        transition: opacity 150ms ease;
      }
      body.msie-busy::after {
        content: 'Working...';
        position: fixed;
        top: 6px;
        left: 50%;
        transform: translateX(-50%);
        background: #3b3b3b;
        color: #fff;
        padding: 4px 14px;
        border-radius: 12px;
        font-size: 12px;
        letter-spacing: 0.03em;
        z-index: 10000;
        pointer-events: none;
      }

      .msie-rule { margin: 10px 0; border-top: 1px solid #eee7dd; }
      .msie-hint { color: #9a9a9a; display: block; margin-top: -4px; }

      /* Reserve the image area so the m/z list does not jump up and down as
         images appear and disappear. */
      .msie-image-slot {
        /* Starting reservation only -- the server sends the real height of the
           image it is about to draw, so the slot matches the plot instead of
           leaving a block of empty space. */
        min-height: 240px;
        /* Draw the image at exactly the requested size; if that is larger than
           the window, scroll rather than shrink (shrinking only the width is
           what distorted the aspect ratio). */
        overflow: auto;
        border: 1px dashed #e6e1da;
        border-radius: 4px;
        padding: 4px;
        margin-bottom: 8px;
        position: relative;
      }
      .msie-image-slot img { max-width: none; }
      .msie-image-slot .msie-slot-hint {
        position: absolute;
        top: 50%;
        left: 0;
        right: 0;
        transform: translateY(-50%);
        text-align: center;
        color: #b8b0a6;
        font-size: 13px;
        pointer-events: none;
        z-index: 0;
      }
      .msie-image-slot .shiny-image-output { position: relative; z-index: 1; }
    ")),
    shiny::tags$script(shiny::HTML("
      function msieToggleOptions(btn) {
        var bar = btn.parentNode;
        while (bar && !(bar.className && bar.className.indexOf('msie-toggle-bar') >= 0)) {
          bar = bar.parentNode;
        }
        if (!bar) { return; }
        var row = bar.nextElementSibling;
        while (row && !(row.className && / row /.test(' ' + row.className + ' '))) {
          row = row.nextElementSibling;
        }
        if (!row) { return; }
        var cols = [];
        for (var i = 0; i < row.children.length; i++) {
          var c = row.children[i];
          if (/\\bcol-(xs|sm|md|lg)-\\d+\\b/.test(c.className)) { cols.push(c); }
        }
        if (cols.length < 2) { return; }
        var sb = cols[0], mn = cols[1];
        if (sb.getAttribute('data-msie-hidden') === '1') {
          sb.style.display = '';
          sb.setAttribute('data-msie-hidden', '0');
          var restore = mn.getAttribute('data-msie-cls');
          if (restore) { mn.className = restore; }
          btn.innerHTML = '&#9664; Hide options';
        } else {
          if (!mn.getAttribute('data-msie-cls')) {
            mn.setAttribute('data-msie-cls', mn.className);
          }
          sb.style.display = 'none';
          sb.setAttribute('data-msie-hidden', '1');
          var keep = mn.className.split(/\\s+/).filter(function (c) {
            return c && !/^col-(xs|sm|md|lg)-\\d+$/.test(c);
          });
          keep.push('col-sm-12');
          mn.className = keep.join(' ');
          btn.innerHTML = '&#9654; Show options';
        }
        if (window.jQuery) { jQuery(window).trigger('resize'); }
        else { window.dispatchEvent(new Event('resize')); }
      }

      // Busy state. Shiny reports when the server is working; after a short
      // grace period (so quick updates do not flicker) the tab content is dimmed
      // and a 'Working...' label is shown. Inputs remain clickable.
      document.addEventListener('DOMContentLoaded', function () {
        if (!window.jQuery) { return; }
        var busyTimer = null;
        jQuery(document).on('shiny:busy', function () {
          if (busyTimer !== null) { return; }
          busyTimer = setTimeout(function () {
            busyTimer = null;
            document.body.classList.add('msie-busy');
          }, 500);
        });
        jQuery(document).on('shiny:idle', function () {
          if (busyTimer !== null) { clearTimeout(busyTimer); busyTimer = null; }
          document.body.classList.remove('msie-busy');
        });
      });

      // The plot card reports the height of the image it is about to draw, so
      // the reserved slot matches the plot exactly.
      if (window.Shiny && Shiny.addCustomMessageHandler) {
        Shiny.addCustomMessageHandler('msieSlotHeight', function (msg) {
          if (!msg || !msg.id) { return; }
          var el = document.getElementById(msg.id);
          if (el && msg.h) { el.style.minHeight = msg.h + 'px'; }
        });
      }

      // Arrow-key navigation of the m/z list. Up/Down move the selection to the
      // neighbouring ion; Shift+Up/Down extend the selection instead of
      // replacing it. The handler finds the visible ion table and reports to
      // that plot card's namespace, so it works on every tab that plots.
      // Only after the user last clicked inside an ion table, so the arrow keys
      // still scroll the page and work in other widgets otherwise.
      document.addEventListener('mousedown', function (e) {
        var t = e.target;
        window.msieLastIonTable = (t && t.closest) ? t.closest('.msie-ion-table') : null;
      }, true);
      document.addEventListener('keydown', function (e) {
        if (e.key !== 'ArrowDown' && e.key !== 'ArrowUp') { return; }
        if (e.metaKey || e.ctrlKey || e.altKey) { return; }
        var t = e.target;
        if (t && (t.tagName === 'INPUT' || t.tagName === 'TEXTAREA' ||
                  t.tagName === 'SELECT' || t.isContentEditable)) { return; }
        var slot = window.msieLastIonTable || null;
        if (!slot || slot.offsetParent === null || !document.body.contains(slot)) { return; }
        var nsPrefix = slot.getAttribute('data-ns') || '';
        if (!window.Shiny || !Shiny.setInputValue) { return; }
        e.preventDefault();
        Shiny.setInputValue(nsPrefix + 'mz_key_nav', {
          dir: (e.key === 'ArrowDown') ? 'down' : 'up',
          add: e.shiftKey === true,
          nonce: Date.now() + Math.random()
        }, { priority: 'event' });
      });
    "))
  )
}

msi_options_toggle <- function(label = "options") {
  shiny::div(
    class = "msie-toggle-bar",
    shiny::tags$button(
      type = "button",
      class = "btn btn-default btn-sm msie-toggle-btn",
      onclick = "msieToggleOptions(this)",
      shiny::HTML(paste0("&#9664; Hide ", label))
    )
  )
}

# Physical RAM in MB, or NA when it can't be determined. Used to cap R's
# vector-memory limit (see configure_cardinal_parallel) below physical RAM.

# Canonical metadata at calibrated imzML export boundaries. Summaries follow saved pixels.
msi_prepare_export <- function(x) {
  fd <- as.data.frame(Cardinal::fData(x))
  if(!msi_is_calibrated(x)) return(x)
  x <- msi_clean_calibration_metadata(x)
  fd <- as.data.frame(Cardinal::fData(x))
  if(!nrow(x) || !ncol(x)) stop("Cannot export an empty calibrated dataset")
  sf <- Cardinal::summarizeFeatures(x,stat=c(mean="mean",count="nnzero"),verbose=FALSE)
  fd$mean <- as.numeric(Cardinal::fData(sf)$mean)
  fd$count <- as.numeric(Cardinal::fData(sf)$count)
  fd$freq <- fd$count/ncol(x)
  fd$mz <- NULL
  Cardinal::fData(x) <- Cardinal::MassDataFrame(mz=Cardinal::mz(x),fd)
  x
}
msi_write_bundle <- function(object, file, ...) {
  object <- msi_prepare_export(object)
  Cardinal::writeImzML(object,file,...)
  invisible(object)
}
msi_dataset_choices <- function(files) {
  files[!is.na(files) & grepl("\\.(imzML|rds|RData)$",basename(files),ignore.case=TRUE) &
          !grepl("_calibration[.]rds$",basename(files),ignore.case=TRUE)]
}
msi_read_bundle <- function(path) {
  if(!file.exists(path)) stop("The selected imzML file is unavailable. Select its current bundle location.")
  if(!grepl("[.]imzML$",path,ignore.case=TRUE)) stop("Select a saved imzML dataset")
  stem <- sub("[.]imzML$","",path,ignore.case=TRUE)
  if(any(file.exists(paste0(stem,c("_masses.tsv","_calibration.txt"))))) {
    missing <- paste0(stem,c(".fdata",".pdata"))
    missing <- missing[!file.exists(missing)]
    if(length(missing)) stop("Calibrated bundle is incomplete: missing ",paste(basename(missing),collapse=", "))
  }
  obj <- msi_clean_calibration_metadata(Cardinal::readMSIData(path))
  if(!inherits(obj,"MSImagingExperiment") || !nrow(obj) || !ncol(obj)) stop("A non-empty processed MSI dataset is required")
  if(nrow(Cardinal::fData(obj))!=nrow(obj) || nrow(Cardinal::pData(obj))!=ncol(obj)) stop("Bundle metadata dimensions do not match the spectra")
  probe <- try(as.matrix(Cardinal::spectra(obj)[unique(c(1L,nrow(obj))),unique(c(1L,ncol(obj))),drop=FALSE]),silent=TRUE)
  if(inherits(probe,"try-error")) stop("Cannot read bundle intensities. Keep the imzML and its linked IBD together. ",as.character(probe))
  obj
}
msi_check_feature_compatibility <- function(objects) {
  if(length(objects)<2) return(invisible(TRUE))
  first <- objects[[1]]; fd <- as.data.frame(Cardinal::fData(first))
  for(other in objects[-1]) {
    od <- as.data.frame(Cardinal::fData(other))
    if(!msi_is_calibrated(first) && !msi_is_calibrated(other)) next
    if(!isTRUE(all.equal(as.numeric(Cardinal::mz(first)),as.numeric(Cardinal::mz(other)),tolerance=0))) stop("Calibrated feature axes differ")
    for(n in c("mz_original","correction_ppm")) {
      if(!n %in% names(fd) || !n %in% names(od) || !identical(fd[[n]],od[[n]])) stop("Incompatible calibrated feature metadata: ",n)
    }
  }
  invisible(TRUE)
}

msi_feature_result_index <- function(results,features) {
  if("feature_id" %in% intersect(names(results),names(features))) return(match(results$feature_id,features$feature_id))
  if(!"mz" %in% names(results) || !"mz" %in% names(features)) stop("Feature results require feature identity or exact m/z")
  match(results$mz,features$mz)
}

# Calibration-specific details belong in reports; only these two fields travel with features.
msi_is_calibrated <- function(x) {
  fd <- Cardinal::fData(x)
  all(c("mz_original","correction_ppm") %in% names(fd)) || "calibration_id" %in% names(fd)
}
msi_clean_calibration_metadata <- function(x) {
  if(!inherits(x,"MSImagingExperiment")) return(x)
  fd <- as.data.frame(Cardinal::fData(x))
  if(!msi_is_calibrated(x)) return(x)
  legacy <- "calibration_id" %in% names(fd)
  drop <- c("mz_observed","outside_calibration_range","calibration_id","calibration_model",
            "calibration_library_version","calibration_schema","calibration_annotation_status")
  # Earlier calibration generated these synthetic IDs. Keep IDs supplied by source datasets.
  if(legacy && "feature_id" %in% names(fd) && identical(as.character(fd$feature_id),paste0("feature_",seq_len(nrow(fd))))) drop <- c(drop,"feature_id")
  if(!any(drop %in% names(fd))) return(x)
  fd <- fd[,!names(fd)%in%c("mz",drop),drop=FALSE]
  Cardinal::fData(x) <- Cardinal::MassDataFrame(mz=Cardinal::mz(x),fd)
  x
}
