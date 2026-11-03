pkgload::load_all('.',quiet=TRUE)
r <- cal_prepare_references(cal_references('positive'),c('[M+H]+','[M+Na]+','[M]+'))
ix <- which(r$n_alt_5ppm==0 & r$theo>100 & r$theo<1000);ix <- ix[unique(round(seq(1,length(ix),length.out=12)))]
m <- c(r$theo[ix]*(1+2e-6),seq(50,95,length.out=50));ab <- c(rep(100,12),rep(1,50));o <- order(m);m <- m[o];ab <- ab[o]
a <- cal_discover(m,ab,r)
x <- Cardinal::MSImagingExperiment(spectraData=S4Vectors::SimpleList(intensity=matrix(ab,length(ab),1)),featureData=Cardinal::MassDataFrame(mz=m),pixelData=Cardinal::PositionDataFrame(coord=data.frame(x=1,y=1),run='run1'),centroided=TRUE)
shiny::testServer(MassCalibrationServer,args=list(setup_values=function()list(wd=tempdir(),rawd=tempdir()),preproc_values=function()NULL),{
 session$setInputs(polarity='positive',extended=FALSE,model='auto',adducts=c('[M+H]+','[M+Na]+','[M]+'),show_all_matches=FALSE)
 source(x);stamp(cal_signature(x));anchors(a);candidate_data(list(refs=r))
 session$flushReact()
 view <- anchor_view();stopifnot(all(view$use),nrow(view)<=nrow(a))
 session$setInputs(anchors_rows_selected=which(view$recommended),fit=1)
 stopifnot(!is.null(fit()),fit()$adducts[1]=='[M+H]+',max(abs(cal_axis(m,fit())$correction_ppm-2))<1e-5)
 session$setInputs(show_all_matches=TRUE)
 stopifnot(nrow(anchor_view())==nrow(a),is.null(fit()))
})
cat('GUI candidate view mapping, recommended row fitting, and view-change invalidation passed\n')
