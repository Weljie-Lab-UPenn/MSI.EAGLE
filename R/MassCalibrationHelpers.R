# Pure calibration functions. All calculations use unrounded masses.
cal_signature <- function(x) paste(nrow(x), ncol(x), paste(format(Cardinal::mz(x), digits=17), collapse=","), sep=":")
cal_references <- function(polarity) {
  path <- system.file("calibration", paste0(polarity, ".tsv"), package="MSI.EAGLE")
  if (!nzchar(path)) path <- file.path("inst", "calibration", paste0(polarity, ".tsv"))
  read.delim(path, stringsAsFactors=FALSE, check.names=FALSE)
}
cal_candidates <- function(mass, abundance, refs, extended=FALSE) {
  stopifnot(length(mass)==length(abundance), all(is.finite(mass)), all(mass>0))
  q <- rank(abundance, ties.method="average") / length(abundance)
  out <- refs[, !names(refs) %in% "evidence_complete", drop=FALSE]
  out$feature <- rep(NA_integer_,nrow(out)); out$obs <- out$mean_intensity <- out$dominance <- out$int_q <- out$iso_gap <- out$m2_interf <- rep(NA_real_,nrow(out))
  for (k in seq_len(nrow(out))) {
    t <- out$theo[k]; win <- if (extended && t<300) 60 else 20
    w <- which(abs(mass/t-1)*1e6 <= win); w <- w[order(-abundance[w])]
    if (!length(w)) next
    j <- w[1]; out$feature[k] <- j; out$obs[k] <- mass[j]; out$mean_intensity[k] <- abundance[j]; out$int_q[k] <- q[j]
    out$dominance[k] <- if(length(w)>1) abundance[j]/max(abundance[w[2]], .Machine$double.eps) else Inf
    z <- which(abs(mass-(mass[j]+1.0033548))/(t+1.0033548)*1e6<=8)
    if(length(z)) out$iso_gap[k] <- (mass[j]/t-1)*1e6 - (mass[z[which.max(abundance[z])]]/(t+1.0033548)-1)*1e6
    z <- which(abs(mass/(t-2.0156501)-1)*1e6<=20)
    out$m2_interf[k] <- if(length(z)) max(abundance[z])*(0.0107*(t-2.0156501)/14)^2/2/max(abundance[j], .Machine$double.eps) else 0
  }
  out$ppm_raw <- (out$obs/out$theo-1)*1e6
  out$reason <- vapply(seq_len(nrow(out)), function(i) {
    if(is.na(out$obs[i])) return("No feature")
    if(out$n_alt_5ppm[i]>0) return("Plausible alternative within 5 ppm")
    if(out$dominance[i]<5) return("Competing peak")
    if(out$int_q[i]<0.75) return("Below top intensity quartile")
    if(!is.na(out$iso_gap[i]) && abs(out$iso_gap[i])>5) return("Isotope disagreement")
    if(out$m2_interf[i]>0.05) return("M+2 interference")
    "Accepted"
  }, character(1))
  out$use <- out$reason=="Accepted"
  # One observation cannot serve as several independent anchors.
  dup <- rep(FALSE,nrow(out)); accepted <- which(out$use)
  accepted <- accepted[order(abs(out$ppm_raw[accepted]),-out$int_q[accepted])]
  dup[accepted] <- duplicated(out$feature[accepted])
  out$use[dup] <- FALSE; out$reason[dup] <- "Duplicate observed feature"
  out$isotope_evidence <- ifelse(is.na(out$iso_gap), "Unavailable", "Observed")
  out
}
cal_design <- function(m, model) {
  x <- m/1000
  switch(model, unchanged=matrix(numeric(), length(m), 0), constant=cbind(intercept=1+0*x),
    linear=cbind(1,x), quadratic=cbind(1,x,x^2), cubic=cbind(1,x,x^2,x^3),
    inverse_linear=cbind(1,1/x,x), inverse_quadratic=cbind(1,1/x,x,x^2),
    log_quadratic=cbind(1,log(x),log(x)^2), log_cubic=cbind(1,log(x),log(x)^2,log(x)^3))
}
cal_fit_one <- function(m,y,model) {
  X <- cal_design(m,model)
  if(model=="unchanged") return(numeric())
  if(length(m)<ncol(X)+1 || qr(X)$rank<ncol(X)) stop("Insufficient independent anchors")
  fit <- lm.fit(X,y)
  if(any(!is.finite(fit$coefficients))) stop("Singular model")
  fit$coefficients
}
cal_predict <- function(m,model,coef) if(model=="unchanged") rep(0,length(m)) else drop(cal_design(m,model)%*%coef)
cal_fit <- function(anchors, extended=FALSE, robust=TRUE) {
  a <- anchors[anchors$use & is.finite(anchors$obs), , drop=FALSE]
  if(nrow(a)<3) stop("At least three independent accepted anchors are required")
  # Bounded rejection, recorded in the returned table.
  pilot <- if(extended) "inverse_quadratic" else "quadratic"
  if(robust) repeat {
    b <- tryCatch(cal_fit_one(a$obs,a$ppm_raw,pilot),error=function(e) NULL)
    if(is.null(b)) break
    r <- a$ppm_raw-cal_predict(a$obs,pilot,b)
    bad <- abs(r-median(r))>3*max(mad(r),0.5)
    if(!any(bad)) break
    j <- which.max(abs(r-median(r)))
    ix <- which(anchors$feature==a$feature[j] & anchors$use)
    anchors$use[ix] <- FALSE; anchors$reason[ix] <- "Robust fit outlier (>3 MAD)"
    a <- a[-j,,drop=FALSE]
  }
  models <- c("unchanged","constant","linear","quadratic","cubic",if(extended)c("inverse_linear","inverse_quadratic","log_quadratic","log_cubic"))
  scores <- lapply(models,function(model) {
    tryCatch({
      if(model != "unchanged" && nrow(a)<ncol(cal_design(a$obs,model))+2) stop("Insufficient independent anchors")
      b <- cal_fit_one(a$obs,a$ppm_raw,model)
      e <- vapply(seq_len(nrow(a)),function(i) {
        bb <- cal_fit_one(a$obs[-i],a$ppm_raw[-i],model)
        a$ppm_raw[i]-cal_predict(a$obs[i],model,bb)
      },numeric(1))
      list(model=model, rmse=sqrt(mean(e^2)), coef=b, reason="Eligible")
    },error=function(e) list(model=model,rmse=Inf,coef=NULL,reason=conditionMessage(e)))
  })
  tab <- data.frame(model=models,loo_rmse_ppm=vapply(scores,`[[`,numeric(1),"rmse"),status=vapply(scores,`[[`,character(1),"reason"))
  best <- which.min(tab$loo_rmse_ppm)
  list(model=models[best],coef=scores[[best]]$coef,domain=range(a$obs),anchors=anchors,models=tab,
       validation="Internal LOO after anchor selection; not an independent accuracy estimate",schema=1L)
}
cal_axis <- function(m,fit) {
  ppm <- cal_predict(pmin(pmax(m,fit$domain[1]),fit$domain[2]),fit$model,fit$coef)
  z <- m/(1+ppm*1e-6)
  if(any(!is.finite(z)) || any(z<=0) || any(diff(z)<=0)) stop("Correction produces invalid or unordered masses")
  data.frame(mz_observed=m,mz_calibrated=z,correction_ppm=ppm,outside_range=m<fit$domain[1]|m>fit$domain[2])
}
cal_apply <- function(x,fit) {
  fd <- as.data.frame(Cardinal::fData(x))
  if(msi_is_calibrated(x)) stop("Restore the original axis before replacing a calibration")
  m <- as.numeric(Cardinal::mz(x)); axis <- cal_axis(m,fit)
  fd$mz <- NULL
  fd$mz_original <- m; fd$correction_ppm <- axis$correction_ppm
  Cardinal::fData(x) <- Cardinal::MassDataFrame(mz=axis$mz_calibrated,fd)
  x
}
cal_restore <- function(x) {
  fd <- as.data.frame(Cardinal::fData(x)); if(!"mz_original" %in% names(fd)) stop("Original axis unavailable")
  m <- fd$mz_original
  fd <- fd[,!names(fd)%in%c("mz","mz_original","mz_observed","correction_ppm","outside_calibration_range","calibration_id","calibration_model","calibration_library_version","calibration_schema","calibration_annotation_status"),drop=FALSE]
  Cardinal::fData(x) <- Cardinal::MassDataFrame(mz=m,fd); x
}

# Row selections are source row indices supplied by DT, including after sorting/filtering.
cal_selected_anchors <- function(mass, abundance, refs, selected, extended=FALSE) {
  if(!length(selected)) stop("Select candidate anchor rows in the table first")
  if(anyNA(selected) || any(selected!=as.integer(selected)) || any(selected<1 | selected>nrow(refs))) stop("Invalid candidate selection; find anchors again")
  selected <- unique(as.integer(selected))
  a <- cal_candidates(mass,abundance,refs,extended)
  omitted <- !seq_len(nrow(a)) %in% selected
  a$use[omitted] <- FALSE; a$reason[omitted] <- "Not selected"
  a
}

# Refit an eligible comparison model on the same accepted anchor set.
cal_choose_model <- function(fit, model) {
  k <- match(model,fit$models$model)
  if(is.na(k) || !is.finite(fit$models$loo_rmse_ppm[k])) stop("This model cannot be fitted with the accepted anchors")
  a <- fit$anchors[fit$anchors$use,,drop=FALSE]
  fit$manual_model <- !identical(model,fit$model)
  fit$model <- model
  fit$coef <- cal_fit_one(a$obs,a$ppm_raw,model)
  fit
}

# Normalize custom/reference ions before mass matching. Same ion formula/charge is one mass target.
cal_prepare_references <- function(refs, adducts=NULL) {
  needed <- c("name","adduct","formula","theo")
  if(!all(needed %in% names(refs))) stop("References require name, adduct, formula and theo")
  if(any(!is.finite(refs$theo)) || any(refs$theo<=0)) stop("Reference masses must be finite and positive")
  if(!is.null(adducts)) refs <- refs[refs$adduct %in% adducts,,drop=FALSE]
  if(!nrow(refs)) stop("No reference ions for the selected adducts")
  if(any(grepl("[2-9][+-]$",refs$adduct))) stop("Only singly charged references are supported")
  if(!"charge" %in% names(refs)) refs$charge <- ifelse(grepl("-$",refs$adduct),-1L,1L)
  if(anyNA(refs$charge) || any(refs$charge != ifelse(grepl("-$",refs$adduct),-1L,1L))) stop("Reference charge conflicts with its adduct")
  if(any(!refs$charge %in% c(-1L,1L))) stop("Only singly charged references are supported")
  if(!"ion_formula" %in% names(refs)) refs$ion_formula <- vapply(seq_len(nrow(refs)),function(i) cal_ion_formula(refs$formula[i],refs$adduct[i]),character(1))
  if(anyNA(refs$ion_formula) || any(!nzchar(refs$ion_formula))) stop("Missing ion formulas")
  if(!"reference_family" %in% names(refs)) refs$reference_family <- refs$name
  if(!"prior" %in% names(refs)) refs$prior <- 0
  groups <- split(seq_len(nrow(refs)),paste(refs$ion_formula,refs$charge,sep="|"))
  refs <- do.call(rbind,lapply(groups,function(ii) {
    x <- refs[ii[1],,drop=FALSE]
    if(diff(range(refs$theo[ii]))>2e-5) stop("Conflicting reference masses for one ion composition")
    for(n in intersect(c("name","adduct","formula","reference_family","source_ids","source_scope","class_name"),names(refs))) x[[n]] <- paste(sort(unique(refs[[n]][ii])),collapse=" / ")
    x$prior <- min(refs$prior[ii]); x
  }))
  rownames(refs) <- NULL
  refs <- refs[order(refs$theo,refs$ion_formula),,drop=FALSE]
  # Search every distinct ion composition in the enabled library, not historical alternatives.
  m <- refs$theo; refs$n_alt_5ppm <- integer(nrow(refs)); refs$alternatives <- character(nrow(refs))
  for(i in seq_along(m)) {
    lo <- findInterval(m[i]*(1-5e-6),m)+1L; hi <- findInterval(m[i]*(1+5e-6),m)
    ix <- if(lo<=hi) seq.int(lo,hi) else integer(0); ix <- setdiff(ix,i)
    refs$n_alt_5ppm[i] <- length(ix)
    if(length(ix)) refs$alternatives[i] <- paste(sprintf("%s %s (%+.3f ppm; %s)",refs$name[ix],refs$adduct[ix],(m[ix]/m[i]-1)*1e6,refs$ion_formula[ix]),collapse="; ")
  }
  refs
}
cal_ion_formula <- function(formula,adduct) {
  hits <- regmatches(formula,gregexpr("[A-Z][a-z]?[0-9]*",formula))[[1]]
  if(paste(hits,collapse="")!=formula) stop("Unsupported elemental formula")
  elements <- sub("[0-9]+$","",hits); counts <- sub("^[A-Z][a-z]?","",hits)
  counts <- ifelse(nzchar(counts),as.integer(counts),1L)
  atoms <- tapply(counts,elements,sum)
  delta <- switch(adduct,"[M+H]+"=c(H=1),"[M+Na]+"=c(Na=1),"[M+K]+"=c(K=1),"[M+NH4]+"=c(N=1,H=4),"[M+H-H2O]+"=c(H=-1,O=-1),"[M]+"=numeric(),"[M-H]-"=c(H=-1),"[M+Cl]-"=c(Cl=1),"[M+HCOO]-"=c(C=1,H=1,O=2),"[M-CH3]-"=c(C=-1,H=-3),"[M-H-H2O]-"=c(H=-3,O=-1),stop("Unsupported adduct"))
  for(el in names(delta)) {if(!el %in% names(atoms)) atoms[el] <- 0; atoms[el] <- atoms[el]+delta[el]}
  if(any(atoms<0)) stop("Adduct gives negative atom counts")
  atoms <- atoms[atoms>0]; atoms <- atoms[order(names(atoms)!="C",names(atoms)!="H",names(atoms))]
  paste0(names(atoms),ifelse(atoms==1,"",atoms),collapse="")
}
cal_discover <- function(mass, abundance, refs, extended=FALSE) {
  # Skip reference rows outside all observed feature windows; broad libraries remain available offline.
  radius <- if(extended) 60 else 20
  refs <- refs[refs$theo>=min(mass)/(1+radius*1e-6) & refs$theo<=max(mass)/(1-radius*1e-6),,drop=FALSE]
  a <- cal_candidates(mass,abundance,refs,extended)
  a <- a[!is.na(a$obs),,drop=FALSE]
  a <- a[order(!a$use,-a$int_q,abs(a$ppm_raw),a$theo),,drop=FALSE]
  a$recommended <- rep(FALSE,nrow(a))
  if(nrow(a)) {
    bins <- if(diff(range(a$obs))>0) cut(a$obs,breaks=seq(min(a$obs)-1e-8,max(a$obs)+1e-8,length.out=9),labels=FALSE) else rep(1L,nrow(a))
    families <- character(); counts <- integer(8)
    for(i in which(a$use)) {
      fam <- strsplit(a$reference_family[i]," / ",fixed=TRUE)[[1]]
      if(counts[bins[i]]>=2L || any(fam %in% families)) next
      a$recommended[i] <- TRUE; counts[bins[i]] <- counts[bins[i]]+1L; families <- c(families,fam)
    }
  }
  rownames(a) <- NULL; a
}
cal_selected_candidates <- function(a,selected) {
  if(!length(selected)) stop("Select candidate anchor rows in the table first")
  if(anyNA(selected)||any(selected!=as.integer(selected))||any(selected<1|selected>nrow(a))) stop("Invalid candidate selection; find anchors again")
  omitted <- !seq_len(nrow(a)) %in% selected
  a$use[omitted] <- FALSE; a$reason[omitted] <- "Not selected"
  a
}

# Format non-finite dominance explicitly for display; leave fitting/report data numeric.
cal_anchor_display <- function(a) {
  ratio <- a$dominance
  a$dominance <- ifelse(is.na(ratio),"Unavailable",ifelse(is.infinite(ratio),"No competing peak",format(signif(ratio,4),trim=TRUE)))
  first <- intersect(c("name","adduct","obs","mean_intensity","dominance","recommended","reason","use"),names(a))
  columns <- c(first,setdiff(names(a),first))
  if(all(c("theo","ppm_raw") %in% columns)) {
    columns <- setdiff(columns,"ppm_raw")
    columns <- append(columns,"ppm_raw",after=match("theo",columns))
  }
  a[,columns,drop=FALSE]
}

# Match other tabs: Cardinal creates a named bundle directory. Keep every companion inside it.
cal_save_bundle <- function(dataset, fit, path) {
  path <- sub("\\.(imzML|rds)$","",path.expand(trimws(path)),ignore.case=TRUE)
  if(length(path)!=1L || !nzchar(path) || basename(path) %in% c(".","..")) stop("Choose a new bundle name")
  if(file.exists(path) || dir.exists(path)) stop("Output already exists; choose a new bundle name")
  if(!dir.exists(dirname(path))) stop("Output directory does not exist")
  dataset <- msi_prepare_export(dataset)
  Cardinal::writeImzML(dataset,path,bundle=TRUE)
  stem <- file.path(path,basename(path))
  exported <- Cardinal::readMSIData(paste0(stem,".imzML"))
  if(!isTRUE(all.equal(as.numeric(Cardinal::mz(exported)),as.numeric(Cardinal::mz(dataset)),tolerance=1e-7))) stop("Exported mass-axis validation failed")
  for(column in names(Cardinal::fData(dataset))) {
    if(!column %in% names(Cardinal::fData(exported))) stop("Feature metadata missing after export: ",column)
    if(!isTRUE(all.equal(Cardinal::fData(dataset)[[column]],Cardinal::fData(exported)[[column]],check.attributes=FALSE,tolerance=1e-7))) stop("Feature metadata changed after export: ",column)
  }
  original_pd <- as.data.frame(Cardinal::pData(dataset)); saved_pd <- as.data.frame(Cardinal::pData(exported))
  rownames(original_pd) <- NULL; rownames(saved_pd) <- NULL
  if(!isTRUE(all.equal(original_pd,saved_pd,check.attributes=FALSE,tolerance=1e-7))) stop("Pixel metadata changed after export")
  dput(fit,file=paste0(stem,"_calibration.txt"))
  write.table(fit$anchors,paste0(stem,"_anchors.tsv"),sep="\t",row.names=FALSE)
  write.table(fit$models,paste0(stem,"_models.tsv"),sep="\t",row.names=FALSE)
  report <- as.data.frame(Cardinal::fData(dataset))
  report$feature_index <- seq_len(nrow(report))
  axis <- cal_axis(report$mz_original,fit)
  report$mz_observed <- axis$mz_observed
  report$outside_calibration_range <- axis$outside_range
  report$calibration_model <- fit$model
  report$calibration_library_version <- if(length(fit$library_version)) paste(fit$library_version,collapse=";") else "unspecified"
  write.table(report,paste0(stem,"_masses.tsv"),sep="\t",row.names=FALSE)
  normalizePath(path)
}
