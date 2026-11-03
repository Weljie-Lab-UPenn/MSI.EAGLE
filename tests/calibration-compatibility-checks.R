source('tests/calibration-save-checks.R')
path <- file.path(moved,'calibrated.imzML')
obj <- msi_read_bundle(path)
fd <- as.data.frame(Cardinal::fData(obj));mat <- as.matrix(Cardinal::spectra(obj))
stopifnot(isTRUE(all.equal(fd$mean,rowMeans(mat))),all(fd$count==rowSums(mat!=0)),all(fd$freq==rowMeans(mat!=0)))
sub <- obj[,1:2];dest <- tempfile();msi_write_bundle(sub,dest)
back <- msi_read_bundle(file.path(dest,paste0(basename(dest),'.imzML')))
stopifnot(isTRUE(all.equal(Cardinal::fData(back)$mean,rowMeans(mat[,1:2]))),identical(as.character(Cardinal::fData(back)$feature_id),as.character(fd$feature_id)))
stopifnot(!any(grepl('_calibration.rds',msi_dataset_choices(c('a.imzML','a.rds','a_calibration.rds','a.tsv')),fixed=TRUE)))
bad <- obj;Cardinal::fData(bad)$correction_ppm <- Cardinal::fData(bad)$correction_ppm+1
stopifnot(inherits(try(msi_check_feature_compatibility(list(obj,bad)),silent=TRUE),'try-error'))
# Use the real Histology module and its actual bundle-load event/ion table.
current <- shiny::reactiveValues(overview_peaks_sel=obj)
shiny::testServer(HistologyIntegrationServer,args=list(setup_values=function()list(wd=root,chunks=1,par_mode=BiocParallel::SerialParam()),preproc_values=function()list(x2=current)),{
 session$setInputs(translate_x=0,translate_y=0,rotate_deg=0,scale_x=1,scale_y=1,msi_source='bundle',bundle_file='calibrated_moved/calibrated.imzML',msi_plot_mode='mz',mz_sort_mode='m/z')
 session$setInputs(load_bundle=1)
 loaded <- msi_data();stopifnot(isTRUE(all.equal(Cardinal::mz(loaded),Cardinal::mz(obj))))
 tab <- build_mz_ion_table(loaded)
 big_n <- 250001L
 big_fd <- Cardinal::fData(obj); big_fd$count <- rep(big_n,nrow(obj)); big_fd$freq <- rep(1,nrow(obj))
 big <- Cardinal::MSImagingExperiment(spectraData=S4Vectors::SimpleList(intensity=matrix(rep(fd$mean,big_n),nrow(obj),big_n)),featureData=big_fd,
   pixelData=Cardinal::PositionDataFrame(coord=data.frame(x=seq_len(big_n),y=1),run=rep('large',big_n)),centroided=TRUE)
 big_tab <- build_mz_ion_table(big)
 stopifnot(nrow(big)*ncol(big)>5e6,identical(big_tab$mean_intensity,as.numeric(fd$mean)),all(big_tab$detect_fraction==1))
 rm(big)

 stopifnot(isTRUE(all.equal(tab$mean_intensity,fd$mean)),all(tab$detect_fraction==fd$freq),identical(tab$feature_id,as.character(seq_len(nrow(obj)))))
 xh$mapped_obj <- obj;xh$stat_fit_summary <- list(stale=TRUE)
 session$setInputs(load_bundle=2)
 stopifnot(is.null(xh$mapped_obj),is.null(xh$stat_fit_summary),identical(current$overview_peaks_sel,obj))
 session$setInputs(msi_source='current')
 stopifnot(isTRUE(all.equal(build_mz_ion_table(msi_data())$mz,tab$mz)))
 mapped <- obj[,1:2]; Cardinal::pData(mapped)$histo_cluster <- c('A','B'); xh$mapped_obj <- mapped
 session$setInputs(save_mapped_imzml=list(root='wd',path=list(),name='mapped_again',type=list()))
 mapped_back <- msi_read_bundle(file.path(root,'mapped_again','mapped_again.imzML'))
 stopifnot(identical(as.character(Cardinal::pData(mapped_back)$histo_cluster),c('A','B')),
   isTRUE(all.equal(Cardinal::fData(mapped_back)$mean,rowMeans(mat[,1:2]))),
   identical(as.character(Cardinal::fData(mapped_back)$feature_id),as.character(fd$feature_id)))
})
cat('Export summaries, subset re-export, companion filtering, calibration compatibility and actual Histology direct/current loading passed\n')
# Missing companion data must fail instead of silently dropping provenance.
fdata_path <- file.path(moved,'calibrated.fdata');temp_fdata <- paste0(fdata_path,'.hold')
stopifnot(file.rename(fdata_path,temp_fdata))
stopifnot(inherits(try(msi_read_bundle(path),silent=TRUE),'try-error'))
stopifnot(file.rename(temp_fdata,fdata_path))
ibd <- file.path(moved,'calibrated.ibd');stopifnot(file.rename(ibd,paste0(ibd,'.hold')))
stopifnot(inherits(try(msi_read_bundle(path),silent=TRUE),'try-error'))
stopifnot(file.rename(paste0(ibd,'.hold'),ibd))
# Recombining split runs retains metadata; full-pixel statistics are refreshed on export.
run2 <- obj[,1:2];Cardinal::pData(run2)$run <- factor(rep('run2',ncol(run2)))
combined <- combine_card(list(obj,run2))
stopifnot(all(c('mz_original','correction_ppm')%in%names(Cardinal::fData(combined))))
idx <- msi_feature_result_index(data.frame(mz=rev(Cardinal::mz(obj))),as.data.frame(Cardinal::fData(obj)))
stopifnot(identical(idx,rev(seq_len(nrow(obj)))))
cat('Missing companions, calibrated run recombination and reordered feature matching passed\n')
