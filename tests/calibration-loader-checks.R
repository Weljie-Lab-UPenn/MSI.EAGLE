source('tests/calibration-save-checks.R')
rel <- 'calibrated_moved/calibrated.imzML'
setup <- function() list(wd=root,rawd=root,x1=list(raw_list=NULL),tol=30,res=10,mz_min=50,mz_max=1200,par_mode=BiocParallel::SerialParam(),units='ppm')
restored <- NULL
shiny::testServer(PeakPickServer,args=list(setup_values=setup),{
 session$setInputs(peak_pick_status='pp_y',peakPickfile=rel,action=1)
 restored <<- overview_peaks()
 stopifnot('correction_ppm'%in%names(Cardinal::fData(restored)))
})
shiny::testServer(StatsPrepServer,args=list(setup_values=setup),{
 choices <- stats_file_choices()
 stopifnot(rel %in% choices,all(file.exists(file.path(root,choices))),
   !any(grepl('__msieagle_|m/z Identification',choices)))
 selector <- output$final_data$html
 stopifnot(grepl(rel,selector,fixed=TRUE),!grepl('m/z Identification',selector,fixed=TRUE))
 session$setInputs(stats_input_file=rel,debug_n=NA_real_,action_read_stats=1)
 stopifnot('correction_ppm'%in%names(Cardinal::fData(x5$data_file)),all(Cardinal::fData(x5$data_file)$freq==1),
   identical(x5$source_file,rel),isTRUE(all.equal(Cardinal::mz(x5$data_file),Cardinal::mz(z))),
   isTRUE(all.equal(as.matrix(Cardinal::spectra(x5$data_file)),as.matrix(Cardinal::spectra(z)))))
})
cat('Actual File Restore and Statistics load events retain calibration metadata\n')

# Previously annotated files remain ordinary datasets after module removal.
annotated <- z
Cardinal::fData(annotated)$annotation <- rep('Existing annotation',nrow(annotated))
Cardinal::fData(annotated)$annotation_accession <- rep('existing-accession',nrow(annotated))
msi_write_bundle(annotated,file.path(root,'previous_annotations'))
shiny::testServer(StatsPrepServer,args=list(setup_values=setup),{
 session$setInputs(stats_input_file='previous_annotations/previous_annotations.imzML',debug_n=NA_real_,action_read_stats=1)
 stopifnot(all(Cardinal::fData(x5$data_file)$annotation=='Existing annotation'),
   all(Cardinal::fData(x5$data_file)$annotation_accession=='existing-accession'),msi_is_calibrated(x5$data_file))
})
cat('Statistics uses saved files and preserves existing annotations\n')

# Render the actual application UI without starting a browser or server.
app <- MSI.EAGLE()
response <- app$httpHandler(list(REQUEST_METHOD='GET',PATH_INFO='/',QUERY_STRING=''))
html <- paste(as.character(response$content),collapse='\n')
stopifnot(response$status==200L,grepl('Mass Calibration',html,fixed=TRUE),
  !grepl('m/z Identification',html,fixed=TRUE),!grepl('tab_identification',html,fixed=TRUE))
calibration_ui <- as.character(MassCalibrationUI('check'))
for(control in c('load','find','fit','apply','undo','open','save_imzml')) {
 stopifnot(grepl(paste0('check-',control,'"'),calibration_ui,fixed=TRUE))
}
cat('Application renders with calibration controls and without identification\n')
