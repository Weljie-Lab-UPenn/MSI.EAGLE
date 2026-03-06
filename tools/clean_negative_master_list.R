#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) args[1] else "/Users/aalim/Library/CloudStorage/Box-Box/DESI (Aalim Weljie)/Peak_list_sources/Negative_master_list_2025-02-01-TOH.txt"
out_dir <- if (length(args) >= 2) args[2] else "inst/peak_lists"
ppm_tol <- if (length(args) >= 3) suppressWarnings(as.numeric(args[3])) else 5
if (!is.finite(ppm_tol) || ppm_tol <= 0) ppm_tol <- 5
resolve_mode <- if (length(args) >= 4) tolower(trimws(args[4])) else "permissive"
if (!resolve_mode %in% c("permissive", "conservative")) resolve_mode <- "permissive"

normalize_spaces <- function(x) {
  gsub("\\s+", " ", trimws(as.character(x)))
}

clean_non_lipid_name <- function(x) {
  y <- normalize_spaces(x)
  y <- gsub('^"|"$', "", y)
  y <- gsub("(?i)^nicotine\\s+amide$", "Nicotinamide", y, perl = TRUE)
  y <- gsub("(?i)^alpha-ketoglutarate\\s*/\\s*oxoglutaric acid$", "Alpha-ketoglutarate", y, perl = TRUE)
  y <- gsub("(?i)^acetyl-glycine$", "N-Acetyl-glycine", y, perl = TRUE)
  y <- gsub("(?i)^malic acid$", "Malate", y, perl = TRUE)
  y <- gsub("(?i)^pyroglutamate$", "Pyroglutamic acid", y, perl = TRUE)
  y <- gsub("(?i)^unknown\\s*-\\s*", "Unknown: ", y, perl = TRUE)
  y <- gsub("\\s*\\((acrolein metab|or related isomer)\\)", "", y, ignore.case = TRUE)
  normalize_spaces(y)
}

detect_lipid_classes <- function(up) {
  out <- character(0)
  if (grepl("\\bCL\\b|CL\\(", up)) out <- c(out, "Cardiolipin")
  if (grepl("\\bLPC\\b|\\[LPC", up)) out <- c(out, "LPC")
  if (grepl("\\bLPE\\b|\\[LPE", up)) out <- c(out, "LPE")
  if (grepl("\\bPC\\b|\\[PC", up)) out <- c(out, "PC")
  if (grepl("\\bPE\\b|\\[PE", up)) out <- c(out, "PE")
  if (grepl("\\bPG\\b|\\[PG", up)) out <- c(out, "PG")
  if (grepl("\\bPI\\b|\\[PI", up)) out <- c(out, "PI")
  if (grepl("\\bPS\\b|\\[PS", up)) out <- c(out, "PS")
  if (grepl("\\bPA\\b|\\[PA", up)) out <- c(out, "PA")
  if (grepl("\\bSM\\b", up)) out <- c(out, "SM")
  if (grepl("GLCCER", up)) out <- c(out, "GlcCer")
  if (grepl("LACCER", up)) out <- c(out, "LacCer")
  if (grepl("\\bCER\\b|CER\\(", up) && !grepl("GLCCER|LACCER", up)) out <- c(out, "Cer")
  if (grepl("\\bTG\\b|^TG", up)) out <- c(out, "TG")
  if (grepl("\\bDG\\b|^DG", up)) out <- c(out, "DG")
  if (grepl("CHOE|CHE\\(|\\bCE\\b", up)) out <- c(out, "CE")
  if (grepl("\\bST\\(|SULFATIDE", up)) out <- c(out, "ST")
  unique(out)
}

extract_comp <- function(up) {
  pair_match <- regexec("([0-9]{1,2})\\s*:\\s*([0-9]{1,2})\\s*/\\s*([0-9]{1,2})\\s*:\\s*([0-9]{1,2})", up, perl = TRUE)
  pair_vals <- regmatches(up, pair_match)[[1]]
  if (length(pair_vals) == 5) {
    tc <- as.integer(pair_vals[2]) + as.integer(pair_vals[4])
    td <- as.integer(pair_vals[3]) + as.integer(pair_vals[5])
    return(sprintf("%d:%d", tc, td))
  }

  sum_match <- regexec("([0-9]{1,2})\\s*:\\s*([0-9]{1,2})", up, perl = TRUE)
  sum_vals <- regmatches(up, sum_match)[[1]]
  if (length(sum_vals) == 3) {
    return(sprintf("%d:%d", as.integer(sum_vals[2]), as.integer(sum_vals[3])))
  }

  NA_character_
}

extract_carbon_only <- function(up) {
  carbon_match <- regexec(
    "(?i)(?:LPC|LPE|PC|PE|PG|PI|PS|PA|SM|CER|TG|DG|CHOE|CHE|CE|ST)\\s*(?:\\(|-|\\s)*(?:O-|P-)?([0-9]{1,2})\\b",
    up,
    perl = TRUE
  )
  carbon_vals <- regmatches(up, carbon_match)[[1]]
  carbon <- if (length(carbon_vals) >= 2) as.integer(carbon_vals[2]) else NA_integer_
  if (is.finite(carbon) && carbon >= 10L) carbon else NA_integer_
}

mono_mass <- c(
  C = 12,
  H = 1.00782503223,
  N = 14.00307400443,
  O = 15.99491461957,
  P = 30.97376199842
)
proton_mass <- 1.007276466812
chloride_mass <- 34.968852682
formate_mass <- mono_mass[["C"]] + mono_mass[["H"]] + 2 * mono_mass[["O"]]
h2_mass <- 2 * mono_mass[["H"]]
max_resolve_err_da <- 0.04

make_name_result <- function(name, rank = 0L, note = "") {
  list(name = name, rank = as.integer(rank), note = normalize_spaces(note))
}

format_lipid_label <- function(class_text, carbon, db, subtype = "") {
  if (identical(class_text, "Cardiolipin")) {
    return(sprintf("Cardiolipin CL(%d:%d) (putative)", carbon, db))
  }
  if (nzchar(subtype)) {
    return(sprintf("%s %s%d:%d (putative)", class_text, subtype, carbon, db))
  }
  sprintf("%s %d:%d (putative)", class_text, carbon, db)
}

format_lipid_placeholder <- function(class_text, carbon, subtype = "") {
  if (identical(class_text, "Cardiolipin")) {
    return("Cardiolipin (putative)")
  }
  if (nzchar(subtype)) {
    return(sprintf("%s %s%d:x (putative)", class_text, subtype, carbon))
  }
  sprintf("%s %d:x (putative)", class_text, carbon)
}

formula_mz <- function(C, H, N = 0, O = 0, P = 0, adduct = 0) {
  C * mono_mass[["C"]] +
    H * mono_mass[["H"]] +
    N * mono_mass[["N"]] +
    O * mono_mass[["O"]] +
    P * mono_mass[["P"]] +
    adduct
}

theoretical_shorthand_mz <- function(label, carbon, db) {
  if (!is.finite(carbon) || carbon <= 0 || !is.finite(db) || db < 0) {
    return(NA_real_)
  }

  switch(
    label,
    LPC_FORMATE = formula_mz(8 + carbon, 18 + 2 * carbon - 2 * db, N = 1, O = 7, P = 1, adduct = formate_mass),
    LPC_O_FORMATE = formula_mz(8 + carbon, 20 + 2 * carbon - 2 * db, N = 1, O = 6, P = 1, adduct = formate_mass),
    LPE_HLOSS = formula_mz(5 + carbon, 12 + 2 * carbon - 2 * db, N = 1, O = 7, P = 1, adduct = -proton_mass),
    LPE_O_HLOSS = formula_mz(5 + carbon, 14 + 2 * carbon - 2 * db, N = 1, O = 6, P = 1, adduct = -proton_mass),
    PC_CL = formula_mz(8 + carbon, 16 + 2 * carbon - 2 * db, N = 1, O = 8, P = 1, adduct = chloride_mass),
    PC_FORMATE = formula_mz(8 + carbon, 16 + 2 * carbon - 2 * db, N = 1, O = 8, P = 1, adduct = formate_mass),
    PC_O_CL = formula_mz(8 + carbon, 18 + 2 * carbon - 2 * db, N = 1, O = 7, P = 1, adduct = chloride_mass),
    PC_O_FORMATE = formula_mz(8 + carbon, 18 + 2 * carbon - 2 * db, N = 1, O = 7, P = 1, adduct = formate_mass),
    PE_HLOSS = formula_mz(5 + carbon, 10 + 2 * carbon - 2 * db, N = 1, O = 8, P = 1, adduct = -proton_mass),
    PE_O_HLOSS = formula_mz(5 + carbon, 12 + 2 * carbon - 2 * db, N = 1, O = 7, P = 1, adduct = -proton_mass),
    PG_HLOSS = formula_mz(6 + carbon, 11 + 2 * carbon - 2 * db, O = 10, P = 1, adduct = -proton_mass),
    PI_HLOSS = formula_mz(9 + carbon, 15 + 2 * carbon - 2 * db, O = 13, P = 1, adduct = -proton_mass),
    PS_HLOSS = formula_mz(6 + carbon, 12 + 2 * carbon - 2 * db, N = 1, O = 10, P = 1, adduct = -proton_mass),
    PA_HLOSS = formula_mz(3 + carbon, 5 + 2 * carbon - 2 * db, O = 8, P = 1, adduct = -proton_mass),
    DG_CL = formula_mz(3 + carbon, 4 + 2 * carbon - 2 * db, O = 5, adduct = chloride_mass),
    SM_CL = formula_mz(5 + carbon, 13 + 2 * carbon - 2 * db, N = 2, O = 6, P = 1, adduct = chloride_mass),
    SM_FORMATE = formula_mz(5 + carbon, 13 + 2 * carbon - 2 * db, N = 2, O = 6, P = 1, adduct = formate_mass),
    CER_CL = formula_mz(carbon, 2 * carbon - 2 * db + 1, N = 1, O = 3, adduct = chloride_mass),
    NA_real_
  )
}

resolve_shorthand_candidate <- function(label, carbon, mz, prefix, subtype = "") {
  sat_mz <- theoretical_shorthand_mz(label, carbon, 0L)
  if (!is.finite(sat_mz)) {
    return(NULL)
  }

  db <- round((sat_mz - mz) / h2_mass)
  if (!is.finite(db) || db < 0L || db > carbon + 10L) {
    return(NULL)
  }

  pred_mz <- theoretical_shorthand_mz(label, carbon, db)
  err <- abs(pred_mz - mz)
  if (!is.finite(err) || err > max_resolve_err_da) {
    return(NULL)
  }

  data.frame(
    prefix = prefix,
    subtype = subtype,
    carbon = carbon,
    db = db,
    err = err,
    stringsAsFactors = FALSE
  )
}

resolve_shorthand_lipid_name <- function(class_text, carbon, mz, subtype = "") {
  candidates <- vector("list", 0)

  add_candidate <- function(label, carbon_adj) {
    if (!is.finite(carbon_adj) || carbon_adj <= 0) {
      return()
    }
    cand <- resolve_shorthand_candidate(label, carbon_adj, mz, class_text, subtype = subtype)
    if (!is.null(cand)) {
      candidates[[length(candidates) + 1L]] <<- cand
    }
  }

  if (identical(class_text, "LPC") && !nzchar(subtype)) {
    add_candidate("LPC_FORMATE", carbon)
    add_candidate("LPC_FORMATE", carbon - 1L)
    add_candidate("LPC_FORMATE", carbon + 1L)
  } else if (identical(class_text, "LPC") && identical(subtype, "O-")) {
    add_candidate("LPC_O_FORMATE", carbon)
    add_candidate("LPC_O_FORMATE", carbon - 1L)
    add_candidate("LPC_O_FORMATE", carbon + 1L)
  } else if (identical(class_text, "LPE") && !nzchar(subtype)) {
    add_candidate("LPE_HLOSS", carbon)
    add_candidate("LPE_HLOSS", carbon - 1L)
    add_candidate("LPE_HLOSS", carbon + 1L)
  } else if (identical(class_text, "LPE") && identical(subtype, "O-")) {
    add_candidate("LPE_O_HLOSS", carbon)
    add_candidate("LPE_O_HLOSS", carbon - 1L)
    add_candidate("LPE_O_HLOSS", carbon + 1L)
  } else if (identical(class_text, "PC") && !nzchar(subtype)) {
    add_candidate("PC_CL", carbon)
    add_candidate("PC_FORMATE", carbon)
    add_candidate("PC_CL", carbon - 1L)
    add_candidate("PC_FORMATE", carbon - 1L)
    add_candidate("PC_CL", carbon + 1L)
    add_candidate("PC_FORMATE", carbon + 1L)
  } else if (identical(class_text, "PC") && identical(subtype, "O-")) {
    add_candidate("PC_O_CL", carbon)
    add_candidate("PC_O_FORMATE", carbon)
    add_candidate("PC_O_CL", carbon - 1L)
    add_candidate("PC_O_FORMATE", carbon - 1L)
    add_candidate("PC_O_CL", carbon + 1L)
    add_candidate("PC_O_FORMATE", carbon + 1L)
  } else if (identical(class_text, "PE") && !nzchar(subtype)) {
    add_candidate("PE_HLOSS", carbon)
  } else if (identical(class_text, "PE") && identical(subtype, "O-")) {
    add_candidate("PE_O_HLOSS", carbon)
  } else if (identical(class_text, "PG") && !nzchar(subtype)) {
    add_candidate("PG_HLOSS", carbon)
  } else if (identical(class_text, "PI") && !nzchar(subtype)) {
    add_candidate("PI_HLOSS", carbon)
  } else if (identical(class_text, "PS") && !nzchar(subtype)) {
    add_candidate("PS_HLOSS", carbon)
  } else if (identical(class_text, "PA") && !nzchar(subtype)) {
    add_candidate("PA_HLOSS", carbon)
  } else if (identical(class_text, "DG")) {
    add_candidate("DG_CL", carbon)
  } else if (identical(class_text, "SM")) {
    add_candidate("SM_CL", carbon)
    add_candidate("SM_FORMATE", carbon)
    add_candidate("SM_CL", carbon + 18L)
    add_candidate("SM_FORMATE", carbon + 18L)
  } else if (identical(class_text, "Cer")) {
    add_candidate("CER_CL", carbon)
    add_candidate("CER_CL", carbon + 2L)
  }

  if (length(candidates) == 0L) {
    return(make_name_result(NA_character_, 0L, ""))
  }

  cand <- do.call(rbind, candidates)
  cand$carbon_shift <- abs(cand$carbon - carbon)
  cand <- cand[order(cand$err, cand$carbon_shift, cand$db, cand$carbon), , drop = FALSE]
  resolved_name <- format_lipid_label(cand$prefix[1], cand$carbon[1], cand$db[1], cand$subtype[1])
  if (cand$carbon_shift[1] > 0L) {
    shift_note <- sprintf(
      "carbon_shift_candidate %d->%d: %s",
      carbon,
      cand$carbon[1],
      resolved_name
    )
    if (identical(resolve_mode, "conservative")) {
      return(make_name_result(format_lipid_placeholder(class_text, carbon, subtype), 100L, shift_note))
    }
    return(make_name_result(resolved_name, 50L, shift_note))
  }
  make_name_result(resolved_name, 50L, "")
}

extract_fa_comp <- function(up) {
  m <- regexec("(?i)(?:^|\\b)(?:FA\\s*)?([0-9]{1,2})\\s*:\\s*([0-9]{1,2})(?:\\b|$)", up, perl = TRUE)
  r <- regmatches(up, m)[[1]]
  if (length(r) >= 3) {
    c_num <- suppressWarnings(as.integer(r[2]))
    db_num <- suppressWarnings(as.integer(r[3]))
    if (is.finite(c_num) && is.finite(db_num) && c_num >= 2 && c_num <= 36) {
      return(sprintf("%d:%d", c_num, db_num))
    }
  }
  NA_character_
}

clean_fa_name <- function(x) {
  y <- normalize_spaces(x)
  up <- toupper(y)

  if (grepl("PALMITIC ACID", up)) return("FA 16:0 (putative)")
  if (grepl("STEARIC ACID", up)) return("FA 18:0 (putative)")
  if (grepl("OLEIC ACID", up)) return("FA 18:1 (putative)")

  fa_comp <- extract_fa_comp(up)
  if (is.na(fa_comp)) return(NA_character_)

  is_hydroxy <- grepl("OH|HYDROXY", up)
  if (is_hydroxy) {
    sprintf("FA %s;O (putative)", fa_comp)
  } else {
    sprintf("FA %s (putative)", fa_comp)
  }
}

clean_lipid_name_info <- function(x, mz = NA_real_) {
  up <- toupper(normalize_spaces(x))
  classes <- detect_lipid_classes(up)
  if (length(classes) == 0) return(make_name_result(NA_character_, 0L, ""))

  class_text <- paste(classes, collapse = "/")
  subtype <- if (grepl("\\bP-", up)) "P-" else if (grepl("\\bO-", up)) "O-" else ""
  comp <- extract_comp(up)

  if (!is.na(comp)) {
    if (identical(class_text, "Cardiolipin")) {
      return(make_name_result(sprintf("Cardiolipin CL(%s) (putative)", comp), 0L, ""))
    }
    comp_parts <- strsplit(comp, ":", fixed = TRUE)[[1]]
    if (length(comp_parts) == 2 && all(nzchar(comp_parts))) {
      return(make_name_result(
        format_lipid_label(class_text, as.integer(comp_parts[1]), as.integer(comp_parts[2]), subtype),
        0L,
        ""
      ))
    }
    if (nzchar(subtype) && length(classes) == 1) {
      return(make_name_result(sprintf("%s %s%s (putative)", class_text, subtype, comp), 0L, ""))
    }
    return(make_name_result(sprintf("%s %s (putative)", class_text, comp), 0L, ""))
  }

  carbon <- extract_carbon_only(up)
  if (is.finite(carbon)) {
    if (length(classes) == 1 && !identical(class_text, "Cardiolipin") && is.finite(mz)) {
      resolved <- resolve_shorthand_lipid_name(class_text, carbon, mz, subtype = subtype)
      if (!is.na(resolved$name)) return(resolved)
    }
    return(make_name_result(format_lipid_placeholder(class_text, carbon, subtype), 100L, ""))
  }

  if (identical(class_text, "Cardiolipin")) {
    return(make_name_result("Cardiolipin (putative)", 0L, ""))
  }
  make_name_result(sprintf("%s (putative)", class_text), 0L, "")
}

clean_name_info <- function(x, mz = NA_real_) {
  lipid <- clean_lipid_name_info(x, mz = mz)
  if (!is.na(lipid$name)) return(lipid)

  fa <- clean_fa_name(x)
  if (!is.na(fa)) return(make_name_result(fa, 0L, ""))

  make_name_result(clean_non_lipid_name(x), 0L, "")
}

canonical_key <- function(x) {
  y <- tolower(normalize_spaces(x))
  y <- gsub("\\(putative\\)", "", y)
  y <- gsub("choe|che", "ce", y)
  y <- gsub("[^a-z0-9]+", "", y)
  y
}

collapse_by_ppm <- function(d, ppm = 5) {
  d <- d[order(d$exact.mass), , drop = FALSE]
  if (!"clean_rank" %in% names(d)) {
    d$clean_rank <- 0L
  }
  if (!"clean_note" %in% names(d)) {
    d$clean_note <- ""
  }
  out <- vector("list", 0)
  i <- 1L

  while (i <= nrow(d)) {
    rep_mass <- d$exact.mass[i]
    use <- i
    j <- i + 1L

    while (j <= nrow(d)) {
      ppm_err <- abs(d$exact.mass[j] - rep_mass) / rep_mass * 1e6
      if (is.finite(ppm_err) && ppm_err <= ppm) {
        use <- c(use, j)
        j <- j + 1L
      } else {
        break
      }
    }

    dd <- d[use, , drop = FALSE]
    pick <- order(dd$clean_rank, nchar(dd$clean_name), dd$clean_name)[1]
    cluster_notes <- unique(dd$clean_note[nzchar(dd$clean_note)])
    out[[length(out) + 1L]] <- data.frame(
      ID = dd$clean_name[pick],
      exact.mass = median(dd$exact.mass),
      name_rank = dd$clean_rank[pick],
      resolve_note = if (length(cluster_notes) > 0) paste(cluster_notes, collapse = " | ") else "",
      n_merged = nrow(dd),
      source_ids = paste(unique(dd$ID), collapse = " | "),
      source_clean_names = paste(unique(dd$clean_name), collapse = " | "),
      mass_min = min(dd$exact.mass),
      mass_max = max(dd$exact.mass),
      stringsAsFactors = FALSE
    )

    i <- max(use) + 1L
  }

  do.call(rbind, out)
}

name_specificity_score <- function(x) {
  x <- as.character(x)
  ifelse(grepl(":x\\b", x, ignore.case = TRUE), 100L, 0L) +
    ifelse(grepl("\\bunknown\\b", x, ignore.case = TRUE), 50L, 0L) +
    ifelse(grepl("\\(putative\\)", x, ignore.case = TRUE), 5L, 0L) +
    nchar(x)
}

cluster_indices_by_ppm <- function(mass, ppm = 5) {
  n <- length(mass)
  if (n == 0) return(integer(0))
  ord <- order(mass)
  cl_ord <- integer(n)
  cid <- 0L
  i <- 1L
  while (i <= n) {
    cid <- cid + 1L
    cl_ord[i] <- cid
    j <- i
    while (j < n) {
      a <- mass[ord[j]]
      b <- mass[ord[j + 1L]]
      ppm_err <- abs(b - a) / max(a, .Machine$double.eps) * 1e6
      if (is.finite(ppm_err) && ppm_err <= ppm) {
        j <- j + 1L
        cl_ord[j] <- cid
      } else {
        break
      }
    }
    i <- j + 1L
  }
  cl <- integer(n)
  cl[ord] <- cl_ord
  cl
}

collapse_global_ppm <- function(d, ppm = 5) {
  keep_cols <- intersect(c("ID", "exact.mass", "name_rank", "resolve_note"), names(d))
  d <- d[is.finite(d$exact.mass) & nzchar(d$ID), keep_cols, drop = FALSE]
  if (!"name_rank" %in% names(d)) {
    d$name_rank <- 0L
  }
  if (!"resolve_note" %in% names(d)) {
    d$resolve_note <- ""
  }
  if (nrow(d) == 0) {
    return(list(
      main = data.frame(ID = character(0), exact.mass = numeric(0), name_rank = integer(0), resolve_note = character(0), stringsAsFactors = FALSE),
      audit = data.frame(ID = character(0), exact.mass = numeric(0), name_rank = integer(0), resolve_note = character(0), n_merged = integer(0),
                         source_ids = character(0), mass_min = numeric(0), mass_max = numeric(0),
                         ppm_span = numeric(0), stringsAsFactors = FALSE)
    ))
  }
  d <- d[order(d$exact.mass, d$ID), , drop = FALSE]
  cl <- cluster_indices_by_ppm(d$exact.mass, ppm = ppm)
  parts <- split(d, cl)

  summarize_cluster <- function(dd) {
    score <- name_specificity_score(dd$ID)
    pick <- order(dd$name_rank, score, dd$ID)[1]
    rep_name <- as.character(dd$ID[pick])
    rep_mass <- stats::median(dd$exact.mass)
    mmin <- min(dd$exact.mass)
    mmax <- max(dd$exact.mass)
    ppm_span <- if (is.finite(rep_mass) && rep_mass > 0) 1e6 * (mmax - mmin) / rep_mass else NA_real_
    cluster_notes <- unique(dd$resolve_note[nzchar(dd$resolve_note)])
    data.frame(
      ID = rep_name,
      exact.mass = rep_mass,
      name_rank = dd$name_rank[pick],
      resolve_note = if (length(cluster_notes) > 0) paste(cluster_notes, collapse = " | ") else "",
      n_merged = nrow(dd),
      source_ids = paste(unique(dd$ID), collapse = " | "),
      mass_min = mmin,
      mass_max = mmax,
      ppm_span = ppm_span,
      stringsAsFactors = FALSE
    )
  }

  audit <- do.call(rbind, lapply(parts, summarize_cluster))
  rownames(audit) <- NULL
  audit <- audit[order(audit$exact.mass, audit$ID), , drop = FALSE]
  main <- audit[, c("ID", "exact.mass", "name_rank", "resolve_note"), drop = FALSE]
  list(main = main, audit = audit)
}

cardiolipin_negative_table <- function() {
  cl <- data.frame(
    species = c(
      "CL(70:5)", "CL(70:4)",
      "CL(72:8)", "CL(72:7)", "CL(72:6)", "CL(72:5)", "CL(72:4)",
      "CL(74:9)", "CL(74:8)",
      "CL(76:11)", "CL(76:10)", "CL(76:9)", "CL(76:12)",
      "CL(78:12)", "CL(80:8)"
    ),
    mz_mh = c(
      1425.9806, 1427.9963,
      1447.9650, 1449.9806, 1451.9963, 1454.0119, 1456.0276,
      1473.9806, 1475.9963,
      1497.9806, 1499.9963, 1502.0119, 1495.9650,
      1523.9963, 1560.0902
    ),
    mz_m2h2 = c(
      712.4867, 713.4945,
      723.4788, 724.4867, 725.4945, 726.5023, 727.5101,
      736.4867, 737.4945,
      748.4867, 749.4945, 750.5023, 747.4788,
      761.4945, 779.5414
    ),
    stringsAsFactors = FALSE
  )

  neg1 <- data.frame(
    ID = sprintf("Cardiolipin %s [M-H]- (putative)", cl$species),
    exact.mass = cl$mz_mh,
    name_rank = 0L,
    resolve_note = "",
    stringsAsFactors = FALSE
  )
  neg2 <- data.frame(
    ID = sprintf("Cardiolipin %s [M-2H]2- (putative)", cl$species),
    exact.mass = cl$mz_m2h2,
    name_rank = 0L,
    resolve_note = "",
    stringsAsFactors = FALSE
  )
  rbind(neg1, neg2)
}

if (!file.exists(input_path)) {
  stop(sprintf("Input file not found: %s", input_path))
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

df <- read.delim(input_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("ID", "exact.mass") %in% colnames(df))) {
  stop("Input file must contain columns: ID, exact.mass")
}

df$ID <- normalize_spaces(df$ID)
df$exact.mass <- suppressWarnings(as.numeric(df$exact.mass))
df <- df[is.finite(df$exact.mass) & nzchar(df$ID), c("ID", "exact.mass")]
name_info <- Map(clean_name_info, df$ID, df$exact.mass)
df$clean_name <- vapply(name_info, function(x) x$name, character(1))
df$clean_rank <- vapply(name_info, function(x) x$rank, integer(1))
df$clean_note <- vapply(name_info, function(x) x$note, character(1))
df$key <- vapply(df$clean_name, canonical_key, character(1))

parts <- lapply(split(df, df$key), collapse_by_ppm, ppm = ppm_tol)
cleaned <- do.call(rbind, parts)
cleaned <- cleaned[order(cleaned$exact.mass), , drop = FALSE]
rownames(cleaned) <- NULL

clean_main <- cleaned[, c("ID", "exact.mass", "name_rank", "resolve_note")]
cl_block <- cardiolipin_negative_table()

clean_main <- rbind(clean_main, cl_block)
pre_global_n <- nrow(clean_main)
global <- collapse_global_ppm(clean_main, ppm = ppm_tol)
clean_main <- global$main[, c("ID", "exact.mass"), drop = FALSE]
cleaned <- global$audit[, setdiff(colnames(global$audit), "name_rank"), drop = FALSE]
removed_ppm_dups <- pre_global_n - nrow(clean_main)

clean_main$exact.mass <- round(clean_main$exact.mass, 4)
cleaned$exact.mass <- round(cleaned$exact.mass, 4)
cleaned$mass_min <- round(cleaned$mass_min, 6)
cleaned$mass_max <- round(cleaned$mass_max, 6)
cleaned$ppm_span <- round(cleaned$ppm_span, 4)

base_name <- "Negative_master_list_2025-02-01-TOH_cleaned"
out_main <- file.path(out_dir, paste0(base_name, ".tsv"))
out_audit <- file.path(out_dir, paste0(base_name, "_audit.tsv"))

write.table(clean_main, file = out_main, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(cleaned, file = out_audit, sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("Input rows: %d\n", nrow(df)))
cat(sprintf("Resolver mode: %s\n", resolve_mode))
cat(sprintf("Output rows: %d\n", nrow(clean_main)))
cat(sprintf("Merged clusters: %d\n", sum(cleaned$n_merged > 1)))
cat(sprintf("Audit rows with resolver notes: %d\n", sum(nzchar(cleaned$resolve_note))))
cat(sprintf("Added cardiolipin entries: %d\n", nrow(cl_block)))
cat(sprintf("Removed final m/z clusters within %.2f ppm: %d\n", ppm_tol, removed_ppm_dups))
cat(sprintf("Wrote cleaned list: %s\n", out_main))
cat(sprintf("Wrote audit list:   %s\n", out_audit))
