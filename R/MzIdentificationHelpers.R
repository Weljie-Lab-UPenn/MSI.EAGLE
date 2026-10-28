# Local exact-mass libraries and MS1 matching.  Database files are deliberately
# separate from the package: LMSD is CC BY and HMDB is downloaded by the user.

mzid_specs <- function() list(
  hmdb = list(id = "hmdb", label = "HMDB", license = "CC BY-NC 4.0", url = "https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip", format = "hmdb_xml"),
  lmsd = list(id = "lmsd", label = "LIPID MAPS LMSD", license = "CC BY 4.0", url = "https://www.lipidmaps.org/files/?ext=sdf.zip&file=LMSD", format = "sdf")
)

mzid_library_dir <- function(wd) {
  out <- file.path(normalizePath(wd, mustWork = FALSE), "MSI.EAGLE_libraries")
  if (!dir.exists(out)) dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out
}
mzid_db_path <- function(wd, id) file.path(mzid_library_dir(wd), paste0(id, ".sqlite"))
mzid_manifest_path <- function(wd) file.path(mzid_library_dir(wd), "manifest.json")
mzid_available <- function(wd) {
  specs <- mzid_specs(); ids <- names(specs)
  data.frame(id = ids, label = vapply(specs, `[[`, "", "label"), installed = file.exists(vapply(ids, function(id) mzid_db_path(wd, id), "")), stringsAsFactors = FALSE)
}
mzid_connect <- function(path) {
  if (!requireNamespace("DBI", quietly = TRUE) || !requireNamespace("RSQLite", quietly = TRUE)) stop("DBI and RSQLite are required for local library matching")
  DBI::dbConnect(RSQLite::SQLite(), path)
}
mzid_adducts <- function(polarity) {
  if (identical(polarity, "positive")) {
    data.frame(adduct = c("[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+"), mass_add = c(1.007276466621, 22.989218, 38.963158, 18.033823), charge = 1L)
  } else {
    data.frame(adduct = c("[M-H]-", "[M+Cl]-", "[M+HCOO]-", "[M+CH3COO]-"), mass_add = c(-1.007276466621, 34.969402, 44.998201, 59.013851), charge = -1L)
  }
}
mzid_feature_ids <- function(dataset) {
  fd <- as.data.frame(Cardinal::fData(dataset))
  if ("msieagle_feature_id" %in% names(fd) && all(!is.na(fd$msieagle_feature_id) & nzchar(fd$msieagle_feature_id))) return(as.character(fd$msieagle_feature_id))
  original <- if ("mz_original" %in% names(fd)) as.numeric(fd$mz_original) else as.numeric(Cardinal::mz(dataset))
  paste0("msieagle:", seq_along(original), ":", format(original, digits = 17, scientific = FALSE, trim = TRUE))
}
mzid_load_dataset <- function(path) {
  if (!file.exists(path)) stop("The selected dataset is unavailable")
  x <- if (grepl("\\.(rds|RData)$", path, ignore.case = TRUE)) readRDS(path) else msi_read_bundle(path)
  x <- msi_clean_calibration_metadata(x)
  if (!inherits(x, "MSImagingExperiment") || !nrow(x) || !ncol(x)) stop("A non-empty processed MSI dataset is required")
  mz <- as.numeric(Cardinal::mz(x))
  if (!length(mz) || any(!is.finite(mz)) || any(diff(mz) <= 0)) stop("The dataset must have a finite, unique, ordered m/z axis")
  if (!msi_is_calibrated(x)) stop("m/z Identification requires a calibrated dataset with corrected m/z metadata")
  x
}
mzid_create_db <- function(path, compounds, metadata) {
  if (file.exists(path)) unlink(path)
  con <- mzid_connect(path); on.exit(DBI::dbDisconnect(con), add = TRUE)
  DBI::dbWriteTable(con, "metadata", as.data.frame(metadata), overwrite = TRUE)
  DBI::dbWriteTable(con, "compounds", compounds, overwrite = TRUE)
  DBI::dbExecute(con, "CREATE INDEX compounds_mass_idx ON compounds(exact_mass)")
  invisible(path)
}
mzid_write_manifest <- function(wd) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) return(invisible(NULL))
  installed <- lapply(mzid_available(wd)$id, function(id) mzid_metadata(wd, id))
  installed <- Filter(Negate(is.null), installed)
  jsonlite::write_json(installed, mzid_manifest_path(wd), pretty = TRUE, auto_unbox = TRUE)
  invisible(mzid_manifest_path(wd))
}
mzid_normalize_compounds <- function(x, library_id) {
  needed <- c("accession", "name", "formula", "exact_mass", "inchikey", "compound_class")
  for (n in needed) if (!n %in% names(x)) x[[n]] <- rep(NA_character_, nrow(x))
  x <- x[is.finite(suppressWarnings(as.numeric(x$exact_mass))) & as.numeric(x$exact_mass) > 0, , drop = FALSE]
  if (!nrow(x)) return(data.frame(library_id = character(), accession = character(), name = character(), formula = character(), exact_mass = numeric(), inchikey = character(), compound_class = character(), stringsAsFactors = FALSE))
  data.frame(library_id = library_id, accession = as.character(x$accession), name = as.character(x$name), formula = as.character(x$formula), exact_mass = as.numeric(x$exact_mass),
    inchikey = as.character(x$inchikey), compound_class = as.character(x$compound_class), stringsAsFactors = FALSE)
}
mzid_lipid_sum_composition <- function(name) {
  if (is.na(name) || !nzchar(trimws(name))) return(name)
  name <- trimws(name)
  parts <- regmatches(name, regexec("^([A-Za-z][A-Za-z0-9_-]*)\\s*\\(?\\s*((?:[OP]-)?[0-9]{1,2}:[0-9]{1,2})", name))[[1]]
  if (length(parts) != 3L) return(name)
  lipid_class <- parts[2]
  chains <- regmatches(name, gregexpr("(?:[OP]-)?[0-9]{1,2}:[0-9]{1,2}", name, perl = TRUE))[[1]]
  if (!length(chains)) return(name)
  chain_parts <- regexec("^((?:[OP]-)?)([0-9]{1,2}):([0-9]{1,2})$", chains, perl = TRUE)
  chain_parts <- regmatches(chains, chain_parts)
  if (any(lengths(chain_parts) != 4L)) return(name)
  prefixes <- vapply(chain_parts, `[[`, "", 2L)
  carbons <- as.integer(vapply(chain_parts, `[[`, "", 3L))
  double_bonds <- as.integer(vapply(chain_parts, `[[`, "", 4L))
  prefix <- unique(prefixes[nzchar(prefixes)])
  if (length(prefix) > 1L) return(name)
  paste0(lipid_class, "(", if (length(prefix)) prefix else "", sum(carbons), ":", sum(double_bonds), ")")
}
mzid_sdf_records <- function(path, library_id) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  blocks <- split(lines, cumsum(lines == "$$$$"))
  rows <- lapply(blocks, function(b) {
    get_field <- function(field) { i <- grep(paste0("^> <", field, ">"), b); if (!length(i) || i[1] == length(b)) NA_character_ else b[i[1] + 1L] }
    raw_name <- get_field("NAME")
    data.frame(accession = get_field("LM_ID"), name = mzid_lipid_sum_composition(raw_name), formula = get_field("FORMULA"), exact_mass = suppressWarnings(as.numeric(get_field("EXACT_MASS"))),
      inchikey = get_field("INCHI_KEY"), compound_class = paste(na.omit(c(get_field("CATEGORY"), get_field("MAIN_CLASS"), get_field("SUB_CLASS"))), collapse = "; "), stringsAsFactors = FALSE)
  })
  mzid_normalize_compounds(do.call(rbind, rows), library_id)
}
mzid_hmdb_records <- function(path, library_id = "hmdb") {
  # Stream the multi-gigabyte HMDB XML; do not construct a full XML document in memory.
  if (!requireNamespace("XML", quietly = TRUE)) stop("XML is required to install HMDB")
  rows <- list(); state <- new.env(parent = emptyenv()); state$stack <- character(); state$record <- NULL
  wanted <- c("accession", "name", "chemical_formula", "monisotopic_molecular_weight", "inchikey")
  start <- function(name, attrs, ...) {
    state$stack <- c(state$stack, name)
    if (identical(name, "metabolite")) state$record <- stats::setNames(as.list(rep(NA_character_, length(wanted))), wanted)
  }
  text <- function(value, ...) {
    if (is.null(state$record) || length(state$stack) != 3L) return()
    field <- state$stack[3]
    if (field %in% wanted) {
      old <- state$record[[field]]
      state$record[[field]] <- paste0(if (is.na(old)) "" else old, value)
    }
  }
  finish <- function(name, ...) {
    if (identical(name, "metabolite") && !is.null(state$record)) {
      rows[[length(rows) + 1L]] <<- data.frame(accession = state$record$accession, name = mzid_lipid_sum_composition(state$record$name), formula = state$record$chemical_formula,
        exact_mass = suppressWarnings(as.numeric(state$record$monisotopic_molecular_weight)), inchikey = state$record$inchikey, compound_class = NA_character_, stringsAsFactors = FALSE)
      state$record <- NULL
    }
    state$stack <- head(state$stack, -1L)
  }
  XML::xmlEventParse(path, handlers = list(startElement = start, endElement = finish, text = text), addContext = FALSE, replaceEntities = FALSE)
  if (!length(rows)) stop("No HMDB metabolite records were parsed")
  mzid_normalize_compounds(do.call(rbind, rows), library_id)
}
mzid_install_archive <- function(wd, id, archive, progress = NULL) {
  specs <- mzid_specs(); spec <- specs[[id]]; if (is.null(spec)) stop("Unknown library: ", id)
  if (!file.exists(archive) || file.info(archive)$size < 100) stop("Official archive is unavailable or empty")
  lib_dir <- mzid_library_dir(wd); unpack <- file.path(lib_dir, paste0(id, "_source_", as.integer(Sys.time())))
  on.exit(unlink(unpack, recursive = TRUE), add = TRUE)
  dir.create(unpack, recursive = TRUE); utils::unzip(archive, exdir = unpack)
  files <- list.files(unpack, recursive = TRUE, full.names = TRUE)
  source_file <- if (identical(spec$format, "sdf")) files[grepl("\\.sdf$", files, ignore.case = TRUE)][1] else files[grepl("hmdb_metabolites\\.xml$", files, ignore.case = TRUE)][1]
  if (is.na(source_file) || !length(source_file)) stop("Official archive did not contain the expected source file")
  if (!is.null(progress)) progress("Building local exact-mass index")
  compounds <- if (identical(spec$format, "sdf")) mzid_sdf_records(source_file, id) else mzid_hmdb_records(source_file, id)
  if (!nrow(compounds)) stop("No usable exact-mass records were found")
  metadata <- data.frame(library_id = id, label = spec$label, release = basename(source_file), installed_at = as.character(Sys.time()), source_url = spec$url, license = spec$license, n_compounds = nrow(compounds), stringsAsFactors = FALSE)
  mzid_create_db(mzid_db_path(wd, id), compounds, metadata)
  mzid_write_manifest(wd)
  metadata
}
mzid_install <- function(wd, id, progress = NULL) {
  specs <- mzid_specs(); spec <- specs[[id]]; if (is.null(spec)) stop("Unknown library: ", id)
  archive <- file.path(mzid_library_dir(wd), paste0(id, "_download.zip")); on.exit(unlink(archive), add = TRUE)
  if (!is.null(progress)) progress("Downloading official release")
  if (!requireNamespace("curl", quietly = TRUE)) stop("The curl package is required to download local libraries")
  response <- tryCatch(curl::curl_fetch_disk(spec$url, archive, handle = curl::new_handle(useragent = "MSI.EAGLE local annotation installer")), error = function(e) e)
  if (inherits(response, "error")) stop("Official download failed: ", conditionMessage(response))
  if (identical(response$status_code, 403L) && identical(id, "hmdb")) {
    stop("HMDB blocked the automated download with a Cloudflare challenge. Click 'Open HMDB download page', download hmdb_metabolites.zip in your browser, then click 'Install HMDB from Downloads'.")
  }
  if (response$status_code < 200L || response$status_code >= 300L) stop("Official download returned HTTP status ", response$status_code)
  mzid_install_archive(wd, id, archive, progress)
}
mzid_downloaded_hmdb_archive <- function(wd) {
  names <- c("hmdb_metabolites.zip", "hmdb_metabolites.xml.zip")
  roots <- unique(c(normalizePath(wd, mustWork = FALSE), path.expand("~/Downloads")))
  candidates <- unlist(lapply(roots, function(root) file.path(root, names)), use.names = FALSE)
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) stop("HMDB archive not found in the working directory or Downloads. Download hmdb_metabolites.zip from the official HMDB page first.")
  candidates[which.max(file.info(candidates)$mtime)]
}
mzid_downloaded_lmsd_archive <- function(wd) {
  roots <- unique(c(normalizePath(wd, mustWork = FALSE), path.expand("~/Downloads")))
  candidates <- unlist(lapply(roots, function(root) {
    if (!dir.exists(root)) return(character())
    list.files(root, pattern = "^LMSD.*[.]zip$", ignore.case = TRUE, full.names = TRUE)
  }), use.names = FALSE)
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) stop("LIPID MAPS LMSD archive not found in the working directory or Downloads. Download the LMSD Structure-data SDF ZIP from the official LIPID MAPS page first.")
  candidates[which.max(file.info(candidates)$mtime)]
}
mzid_metadata <- function(wd, id) {
  path <- mzid_db_path(wd, id); if (!file.exists(path)) return(NULL)
  con <- mzid_connect(path); on.exit(DBI::dbDisconnect(con), add = TRUE); DBI::dbReadTable(con, "metadata")
}
mzid_match_one_library <- function(mz, id, wd, adducts, ppm) {
  path <- mzid_db_path(wd, id); if (!file.exists(path)) return(data.frame())
  con <- mzid_connect(path); on.exit(DBI::dbDisconnect(con), add = TRUE)
  adducts <- adducts[adducts$adduct %in% adducts$adduct, , drop = FALSE]
  pieces <- lapply(seq_len(nrow(adducts)), function(i) {
    a <- adducts[i, ]; low <- (mz * (1 - ppm * 1e-6) - a$mass_add) * abs(a$charge); high <- (mz * (1 + ppm * 1e-6) - a$mass_add) * abs(a$charge)
    q <- DBI::dbGetQuery(con, "SELECT library_id, accession, name, formula, exact_mass, inchikey, compound_class FROM compounds WHERE exact_mass BETWEEN ? AND ?", params = list(low, high))
    if (!nrow(q)) return(NULL)
    q$adduct <- a$adduct; q$charge <- a$charge; q$theoretical_mz <- (q$exact_mass / abs(a$charge) + a$mass_add / abs(a$charge)); q$ppm_error <- 1e6 * (mz - q$theoretical_mz) / q$theoretical_mz; q
  })
  out <- do.call(rbind, Filter(Negate(is.null), pieces)); if (is.null(out)) data.frame() else out[abs(out$ppm_error) <= ppm, , drop = FALSE]
}
mzid_match <- function(mz, ids, wd, polarity, selected_adducts, ppm) {
  all_adducts <- mzid_adducts(polarity); adducts <- all_adducts[all_adducts$adduct %in% selected_adducts, , drop = FALSE]
  if (!length(ids) || !nrow(adducts)) return(data.frame())
  out <- lapply(seq_along(mz), function(i) {
    x <- do.call(rbind, lapply(ids, mzid_match_one_library, mz = mz[i], wd = wd, adducts = adducts, ppm = ppm))
    if (is.null(x) || !nrow(x)) return(NULL)
    x$feature_index <- i; x$mz_corrected <- mz[i]; x
  })
  out <- do.call(rbind, Filter(Negate(is.null), out)); if (is.null(out)) return(data.frame())
  out$abs_ppm_error <- abs(out$ppm_error); out <- out[order(out$feature_index, out$abs_ppm_error, out$library_id, out$accession), , drop = FALSE]
  out$candidate_rank <- ave(out$abs_ppm_error, out$feature_index, FUN = rank, ties.method = "first")
  rownames(out) <- NULL; out
}
mzid_feature_table <- function(dataset, candidates, accepted = data.frame()) {
  mz <- as.numeric(Cardinal::mz(dataset)); out <- data.frame(feature_index = seq_along(mz), mz_original = if ("mz_original" %in% names(Cardinal::fData(dataset))) as.numeric(Cardinal::fData(dataset)$mz_original) else mz, mz_corrected = mz, stringsAsFactors = FALSE)
  fd <- as.data.frame(Cardinal::fData(dataset))
  out$annotation <- if ("annotation" %in% names(fd)) as.character(fd$annotation) else NA_character_
  out$database_count <- 0L; out$candidate_count <- 0L; out$top_candidate <- NA_character_; out$top_adduct <- NA_character_; out$top_ppm_error <- NA_real_
  if (nrow(candidates)) {
    counts <- aggregate(library_id ~ feature_index, candidates, function(x) length(unique(x)))
    names(counts)[2] <- "database_count"; hits <- aggregate(accession ~ feature_index, candidates, length); names(hits)[2] <- "candidate_count"
    top <- candidates[order(candidates$feature_index, abs(candidates$ppm_error), candidates$library_id, candidates$accession), c("feature_index", "name", "adduct", "ppm_error")]
    top <- top[!duplicated(top$feature_index), , drop = FALSE]
    names(top)[2:4] <- c("top_candidate", "top_adduct", "top_ppm_error")
    ix <- match(counts$feature_index, out$feature_index); out$database_count[ix] <- counts$database_count
    ix <- match(hits$feature_index, out$feature_index); out$candidate_count[ix] <- hits$candidate_count
    ix <- match(top$feature_index, out$feature_index); out$top_candidate[ix] <- top$top_candidate; out$top_adduct[ix] <- top$top_adduct; out$top_ppm_error[ix] <- top$top_ppm_error
  }
  out$annotation_status <- "Unreviewed"
  if (nrow(accepted)) { ix <- match(accepted$feature_index, out$feature_index); out$annotation_status[ix] <- "Curator-accepted tentative annotation"; out$annotation[ix] <- accepted$name }
  out
}
mzid_accepted_from_dataset <- function(dataset) {
  fd <- as.data.frame(Cardinal::fData(dataset))
  if (!"annotation" %in% names(fd)) return(data.frame())
  keep <- !is.na(fd$annotation) & nzchar(as.character(fd$annotation))
  if (!any(keep)) return(data.frame())
  pull <- function(name, default = NA) if (name %in% names(fd)) fd[[name]][keep] else rep(default, sum(keep))
  data.frame(feature_index = which(keep), name = as.character(fd$annotation[keep]), formula = as.character(pull("annotation_formula")),
    adduct = as.character(pull("annotation_adduct")), library_id = as.character(pull("annotation_database")), accession = as.character(pull("annotation_accession")),
    ppm_error = as.numeric(pull("annotation_ppm_error")), stringsAsFactors = FALSE)
}
mzid_attach_accepted <- function(dataset, accepted) {
  fd <- as.data.frame(Cardinal::fData(dataset)); n <- nrow(fd)
  if (!"msieagle_feature_id" %in% names(fd)) fd$msieagle_feature_id <- mzid_feature_ids(dataset)
  annotation_fields <- c("annotation", "annotation_status", "annotation_formula", "annotation_adduct", "annotation_database", "annotation_accession", "annotation_ppm_error")
  retained <- fd[, !names(fd) %in% c("mz", annotation_fields, "accepted_annotation", "accepted_formula", "accepted_adduct", "accepted_database", "accepted_accession", "accepted_ppm_error"), drop = FALSE]
  ann <- data.frame(annotation = rep(NA_character_, n), annotation_status = rep("Unreviewed", n), annotation_formula = rep(NA_character_, n), annotation_adduct = rep(NA_character_, n), annotation_database = rep(NA_character_, n), annotation_accession = rep(NA_character_, n), annotation_ppm_error = rep(NA_real_, n), stringsAsFactors = FALSE)
  if (!is.null(accepted) && nrow(accepted)) {
    a <- accepted[!duplicated(accepted$feature_index), , drop = FALSE]; ix <- a$feature_index[a$feature_index >= 1 & a$feature_index <= n]
    aa <- a[match(ix, a$feature_index), , drop = FALSE]
    ann$annotation_status[ix] <- "Curator-accepted tentative annotation"; ann$annotation[ix] <- aa$name; ann$annotation_formula[ix] <- aa$formula; ann$annotation_adduct[ix] <- aa$adduct; ann$annotation_database[ix] <- aa$library_id; ann$annotation_accession[ix] <- aa$accession; ann$annotation_ppm_error[ix] <- aa$ppm_error
  }
  Cardinal::fData(dataset) <- Cardinal::MassDataFrame(mz = Cardinal::mz(dataset), cbind(ann, retained)); dataset
}
mzid_save_bundle <- function(dataset, accepted, candidates, selected_libraries, wd, path) {
  path <- sub("\\.(imzML|rds)$", "", path.expand(trimws(path)), ignore.case = TRUE)
  if (length(path) != 1L || !nzchar(path) || basename(path) %in% c(".", "..")) stop("Choose a new bundle name")
  if (file.exists(path) || dir.exists(path)) stop("Output already exists; choose a new bundle name")
  if (!dir.exists(dirname(path))) stop("Output directory does not exist")
  annotated <- msi_prepare_export(mzid_attach_accepted(dataset, accepted))
  Cardinal::writeImzML(annotated, path, bundle = TRUE)
  stem <- file.path(path, basename(path)); exported <- Cardinal::readMSIData(paste0(stem, ".imzML"))
  expected <- c("mz", "annotation", "annotation_status", "annotation_formula", "annotation_adduct", "annotation_database", "annotation_accession", "annotation_ppm_error")
  if (!identical(names(Cardinal::fData(exported))[seq_along(expected)], expected)) stop("Exported feature metadata does not retain annotation columns beside m/z")
  for (field in expected[-1]) if (!isTRUE(all.equal(Cardinal::fData(annotated)[[field]], Cardinal::fData(exported)[[field]], check.attributes = FALSE, tolerance = 1e-7))) stop("Exported annotation metadata changed: ", field)
  feature_ids <- mzid_feature_ids(annotated)
  accepted_out <- accepted
  if (!is.null(accepted_out) && nrow(accepted_out)) accepted_out$feature_id <- feature_ids[accepted_out$feature_index]
  candidates_out <- candidates
  if (!is.null(candidates_out) && nrow(candidates_out)) candidates_out$feature_id <- feature_ids[candidates_out$feature_index]
  utils::write.table(accepted_out, paste0(stem, "_annotations.tsv"), sep = "\t", row.names = FALSE)
  utils::write.table(candidates_out, paste0(stem, "_annotation_candidates.tsv"), sep = "\t", row.names = FALSE)
  manifest <- list(saved_at = as.character(Sys.time()), selected_libraries = selected_libraries, libraries = lapply(selected_libraries, function(id) mzid_metadata(wd, id)))
  jsonlite::write_json(manifest, paste0(stem, "_annotation_manifest.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
  normalizePath(path)
}
