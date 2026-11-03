pkgload::load_all('.',quiet=TRUE)
for(p in c('positive','negative')) {
 r <- cal_references(p); stopifnot(nrow(r)>10000,all(is.finite(r$theo)),all(r$charge==if(p=='positive')1 else -1))
 z <- cal_prepare_references(r)
 stopifnot(!anyDuplicated(paste(z$ion_formula,z$charge)),nrow(z)<nrow(r))
 # Independently check each ambiguity count against a direct full-library comparison.
 for(i in unique(round(seq(1,nrow(z),length.out=30)))) stopifnot(z$n_alt_5ppm[i]==sum(abs(z$theo/z$theo[i]-1)*1e6<=5)-1)
 x <- cal_prepare_references(r,if(p=='positive')'[M+H]+' else '[M-H]-')
 stopifnot(all(x$adduct==if(p=='positive')'[M+H]+' else '[M-H]-'))
}
stopifnot(cal_ion_formula('C2H6O','[M+H]+')=='C2H7O',cal_ion_formula('C2H6O','[M-H]-')=='C2H5O')
# Same ion composition, different names/adduct pathways: one target, not competitors.
r <- data.frame(name=c('ethanol','ethyl cation','nearby'),formula=c('C2H6O','C2H7O','CH4'),adduct=c('[M+H]+','[M]+','[M+H]+'),theo=c(47.0491413,47.0491413,47.0492413))
z <- cal_prepare_references(r)
stopifnot(nrow(z)==2,all(z$n_alt_5ppm==1),any(grepl('ethanol / ethyl cation',z$name,fixed=TRUE)))
# Expanded-library discovery: use theoretical peaks spread over the range and low-abundance background.
r <- cal_prepare_references(cal_references('positive'),c('[M+H]+','[M+Na]+','[M]+'))
ix <- which(r$n_alt_5ppm==0 & r$theo>100 & r$theo<1000)
ix <- ix[unique(round(seq(1,length(ix),length.out=12)))]
m <- c(r$theo[ix]*(1+2e-6),seq(50,95,length.out=50));ab <- c(rep(100,12),rep(1,50));o <- order(m)
a <- cal_discover(m[o],ab[o],r)
stopifnot(nrow(a)>0,any(a$use),sum(a$recommended)>2,!anyDuplicated(a$feature[a$use]))
selected <- which(a$recommended)
f <- cal_fit(cal_selected_candidates(a,selected),robust=FALSE)
stopifnot(max(abs(cal_axis(m[o],f)$correction_ppm-2))<1e-5)
cat('Library formula grouping, adduct filters, expanded ambiguity search, unique feature anchors, ranking and correction recovery passed\n')

stopifnot(nrow(cal_discover(c(1,2),c(1,1),r))==0)
bad <- data.frame(name="bad",formula="C2H6O",adduct="[M+2H]2+",theo=24.0)
stopifnot(inherits(try(cal_prepare_references(bad),silent=TRUE),"try-error"))
