## Consolidate results into Parquet outputs
#' Merge all *_performance.tsv into one table (plus metadata) and write Parquet inside results path
#'
#' - Uses createMLResultDir() to find the ML_performance directory under `path`
#' - Parses filenames using the same semantics as createMLinputList()
#' - Binds rows from all TSVs, adds parsed columns
#' - Writes a single Parquet file **inside** the ML_performance directory
#'
#' @param path Root results path (the same 'path' you pass to createMLResultDir)
#' @param stratify_by NULL | "year" | "country"
#' @param LOO logical; default FALSE
#' @param MDR logical; default FALSE
#' @param cross_test logical; default FALSE
#' @param out_parquet optional filename (no directories). If NULL, defaults to "all_performance.parquet".
#'                    If a path is given, only its basename is used; it is written in ML_performance/.
#' @param compression parquet compression ("zstd" or "snappy"); default "zstd"
#' @param verbose logical; print progress messages
#' @return A tibble with all performance rows + parsed metadata columns
buildPerfPq <- function(
  path,
  stratify_by = NULL,
  LOO = FALSE,
  MDR = FALSE,
  cross_test = FALSE,
  out_parquet = NULL,     # only filename; will be written under ML_performance/
  compression = "zstd",
  verbose = TRUE
) {
  # -----------------------
  # Validate inputs
  # -----------------------
  if (!is.character(path) || length(path) != 1 || is.na(path) || nchar(path) == 0) {
    stop("`path` must be a non-empty character scalar.")
  }
  path <- normalizePath(path)

  if (!is.null(stratify_by) && !stratify_by %in% c("year", "country")) {
    stop("`stratify_by` must be NULL, 'year', or 'country'.")
  }
  if (isTRUE(LOO) && is.null(stratify_by)) {
    stop("With LOO=TRUE, stratify_by must be 'year' or 'country'.")
  }
  if (isTRUE(MDR) && (!is.null(stratify_by) || isTRUE(LOO) || isTRUE(cross_test))) {
    stop("MDR can only run when stratify_by = NULL, LOO = FALSE, cross_test = FALSE.")
  }

  # -----------------------
  # Resolve directories from your function (ensures they exist)
  # -----------------------
  paths <- createMLResultDir(path, stratify_by = stratify_by, LOO = LOO,
                             cross_test = cross_test, MDR = MDR)
  perf_dir <- paths$ML_performance

  # -----------------------
  # Locate all performance TSVs
  # -----------------------
  perf_files <- list.files(perf_dir, pattern = "performance\\.tsv$", full.names = TRUE, recursive = TRUE)
  if (length(perf_files) == 0) {
    if (verbose) message("No *_performance.tsv files found under: ", perf_dir)
    return(tibble::tibble())
  }

  # -----------------------
  # Helpers
  # -----------------------
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  .NA_chr <- function() NA_character_

  .find_drug_label_value <- function(tokens) {
    # Finds first "drug" and determines if it's "drug" or "drug_class"
    idx <- which(tokens == "drug")
    if (length(idx) == 0) {
      return(list(drug_label = .NA_chr(), drug_value = .NA_chr(), label_end = NA_integer_))
    }
    i <- idx[1]
    if (i < length(tokens) && identical(tokens[i + 1], "class")) {
      list(drug_label = "drug_class",
           drug_value = if (i + 2 <= length(tokens)) tokens[i + 2] else .NA_chr(),
           label_end  = i + 1)
    } else {
      list(drug_label = "drug",
           drug_value = if (i + 1 <= length(tokens)) tokens[i + 1] else .NA_chr(),
           label_end  = i)
    }
  }

  .parse_base <- function(base_no_suffix) {
    xs <- strsplit(base_no_suffix, "_", fixed = TRUE)[[1]]
    n  <- length(xs)

    # Initialize
    species         <- .NA_chr()
    mdr_tag         <- .NA_chr()
    phenotype       <- .NA_chr()
    drug_label      <- .NA_chr()
    drug_or_class   <- .NA_chr()
    strat_label     <- .NA_chr()
    strat_value     <- .NA_chr()
    strat_value_test<- .NA_chr()
    leaveout        <- FALSE
    is_cross        <- FALSE
    ref_drug        <- .NA_chr()
    test_drug       <- .NA_chr()
    prefix_key      <- .NA_chr()
    feature         <- .NA_chr()
    feature_type    <- .NA_chr()
    feature_subtype <- .NA_chr()

    # Require at least 2 tokens for feature
    if (n >= 2) {
      feature         <- paste(xs[(n - 1):n], collapse = "_")
      feature_type    <- xs[n - 1]
      feature_subtype <- xs[n]
    }
    core <- if ((n - 2) >= 1) xs[1:(n - 2)] else character(0)
    core_str <- paste(core, collapse = "_")

    # MDR output prefixes are "MDR_<phenotype>_<feature>"
    if (length(core) > 0 && identical(core[1], "MDR")) {
      mdr_tag <- "MDR"
      if (length(core) >= 2) phenotype <- paste(core[-1], collapse = "_")
      prefix_key <- "MDR"
      return(list(
        species=species, mdr_tag=mdr_tag, phenotype=phenotype,
        drug_label=drug_label, drug_or_class=drug_or_class,
        strat_label=strat_label, strat_value=strat_value, strat_value_test=strat_value_test,
        leaveout=leaveout, is_cross=is_cross,
        ref_drug=ref_drug, test_drug=test_drug,
        prefix_key=prefix_key,
        feature=feature, feature_type=feature_type, feature_subtype=feature_subtype
      ))
    }

    # Cross-test variants
    if (grepl("_tested_on_", core_str, fixed = TRUE)) {
      is_cross <- TRUE
      # 1) LOO cross-test: ... <prefix_key>_<drug_value>_leaveout_tested_on_<strat_value>_<feature>
      if (grepl("_leaveout_tested_on_", core_str, fixed = TRUE)) {
        leaveout <- TRUE
        di <- .find_drug_label_value(core)
        drug_label    <- di$drug_label
        drug_or_class <- di$drug_value
        label_end     <- di$label_end
        if (!is.na(label_end)) {
          prefix_key <- paste(core[1:label_end], collapse = "_")
          i_val <- label_end + 1
          # expect: value, leaveout, tested, on, strat_value
          if ((i_val + 4) <= length(core) &&
              core[i_val + 1] == "leaveout" && core[i_val + 2] == "tested" && core[i_val + 3] == "on") {
            strat_label <- stratify_by %||% .NA_chr()
            strat_value <- core[i_val + 4]
          }
        }
        species <- if (length(core) >= 1) core[1] else .NA_chr()
        return(list(
          species=species, mdr_tag=mdr_tag, phenotype=phenotype,
          drug_label=drug_label, drug_or_class=drug_or_class,
          strat_label=strat_label, strat_value=strat_value, strat_value_test=strat_value_test,
          leaveout=leaveout, is_cross=is_cross,
          ref_drug=ref_drug, test_drug=test_drug,
          prefix_key=prefix_key,
          feature=feature, feature_type=feature_type, feature_subtype=feature_subtype
        ))
      }

      # 2) Cross by strat group: ... <prefix_key>_<drug_value>_cross_<strat_value>_tested_on_<strat_value_test>_<feature>
      if (grepl("_cross_", core_str, fixed = TRUE)) {
        di <- .find_drug_label_value(core)
        drug_label    <- di$drug_label
        drug_or_class <- di$drug_value
        label_end     <- di$label_end
        if (!is.na(label_end)) {
          prefix_key <- paste(core[1:label_end], collapse = "_")
          i_val <- label_end + 1
          # expect: value, cross, strat_value, tested, on, strat_value_test
          if ((i_val + 5) <= length(core) &&
              core[i_val + 1] == "cross" && core[i_val + 3] == "tested" && core[i_val + 4] == "on") {
            strat_label     <- stratify_by %||% .NA_chr()
            strat_value     <- core[i_val + 2]
            strat_value_test<- core[i_val + 5]
          }
        }
        species <- if (length(core) >= 1) core[1] else .NA_chr()
        return(list(
          species=species, mdr_tag=mdr_tag, phenotype=phenotype,
          drug_label=drug_label, drug_or_class=drug_or_class,
          strat_label=strat_label, strat_value=strat_value, strat_value_test=strat_value_test,
          leaveout=leaveout, is_cross=is_cross,
          ref_drug=ref_drug, test_drug=test_drug,
          prefix_key=prefix_key,
          feature=feature, feature_type=feature_type, feature_subtype=feature_subtype
        ))
      }

      # 3) Cross-test (drug vs drug): ... <prefix_key>_<ref_drug>_tested_on_<test_drug>_<feature>
      di <- .find_drug_label_value(core)
      label_end <- di$label_end
      if (!is.na(label_end)) {
        prefix_key <- paste(core[1:label_end], collapse = "_")
        if ((label_end + 4) <= length(core) &&
            core[label_end + 2] == "tested" && core[label_end + 3] == "on") {
          ref_drug  <- core[label_end + 1]
          test_drug <- core[label_end + 4]
          drug_label <- if (endsWith(prefix_key, "drug_class")) "drug_class" else "drug"
          drug_or_class <- ref_drug
        }
      }
      species <- if (length(core) >= 1) core[1] else .NA_chr()
      return(list(
        species=species, mdr_tag=mdr_tag, phenotype=phenotype,
        drug_label=drug_label, drug_or_class=drug_or_class,
        strat_label=strat_label, strat_value=strat_value, strat_value_test=strat_value_test,
        leaveout=leaveout, is_cross=is_cross,
        ref_drug=ref_drug, test_drug=test_drug,
        prefix_key=prefix_key,
        feature=feature, feature_type=feature_type, feature_subtype=feature_subtype
      ))
    }

# Non-cross variants (may be stratified)
species <- if (length(core) >= 1) core[1] else .NA_chr()
di <- .find_drug_label_value(core)
drug_label    <- di$drug_label
label_end     <- di$label_end

if (!is.na(label_end)) {

  # -------- STRATIFIED CASE --------
  # pattern: species drug[_class] strat_label drug_value strat_value
  if (label_end + 3 <= length(core) &&
      core[label_end + 1] %in% c("year", "country")) {

    strat_label <- core[label_end + 1]
    drug_or_class <- core[label_end + 2]     # <-- fixed: this is the FLQ/MAC/etc
    strat_value <- core[label_end + 3]       # <-- fixed: this is "2015-2019"
    prefix_key <- paste(core[1:label_end], collapse = "_")

  # -------- UNSTRATIFIED CASE --------
  } else {
    drug_or_class <- di$drug_value
    prefix_key <- paste(core[1:label_end], collapse = "_")
  }
}

    list(
      species=species, mdr_tag=mdr_tag, phenotype=phenotype,
      drug_label=drug_label, drug_or_class=drug_or_class,
      strat_label=strat_label, strat_value=strat_value, strat_value_test=strat_value_test,
      leaveout=leaveout, is_cross=is_cross,
      ref_drug=ref_drug, test_drug=test_drug,
      prefix_key=prefix_key,
      feature=feature, feature_type=feature_type, feature_subtype=feature_subtype
    )
  }

  # -----------------------
  # Read, parse and bind
  # -----------------------
  out <- purrr::map_dfr(perf_files, function(f) {
    base <- basename(f)
    base_no_suffix <- sub("_performance\\.tsv$", "", base)
    if (identical(base_no_suffix, base)) {
      base_no_suffix <- sub("\\.tsv$", "", base)  # fallback
    }

    meta <- .parse_base(base_no_suffix)

    df <- tryCatch(
      readr::read_tsv(f, show_col_types = FALSE, progress = FALSE),
      error = function(e) {
        if (verbose) message("Failed to read TSV: ", f, " (", conditionMessage(e), ") — using metadata only.")
        tibble()
      }
    )

    md_cols <- tibble(

      output_prefix   = base_no_suffix,
      species         = meta$species,
      mdr_tag         = meta$mdr_tag,
      phenotype       = meta$phenotype,
      drug_label      = meta$drug_label,
      drug_or_class   = meta$drug_or_class,
      strat_label     = meta$strat_label,
      strat_value     = meta$strat_value,
      strat_value_test= meta$strat_value_test,
      leaveout        = meta$leaveout,
      cross_test      = meta$is_cross,
      ref_drug        = meta$ref_drug,
      test_drug       = meta$test_drug,
      prefix_key      = meta$prefix_key,
      feature         = meta$feature,
      feature_type    = meta$feature_type,
      feature_subtype = meta$feature_subtype
    )

    if (nrow(df) == 0) {
      md_cols
    } else {
      dplyr::bind_cols(md_cols[rep(1, nrow(df)), ], df)
    }
  })

  # -----------------------
  # Compute the output Parquet path (ALWAYS inside perf_dir)
  # -----------------------
  # If user gives a name/path, we keep only the basename and write under perf_dir.
  parquet_name <- (out_parquet %||% "all_performance.parquet")
  parquet_name <- basename(parquet_name)
  out_path <- file.path(perf_dir, parquet_name)

  # -----------------------
  # Write Parquet
  # -----------------------
  suppressPackageStartupMessages(library(arrow))
  arrow::write_parquet(out, out_path, compression = compression)
  if (verbose) message("Wrote merged Parquet: ", out_path, "  [", nrow(out), " rows]")

  out
}


#' Merge all *_top_features.tsv into one table + metadata, write Parquet inside results path
#'
#' - Uses createMLResultDir() to find the ML_top_features directory under `path`
#' - Parses filenames to derive metadata (aligned with createMLinputList() semantics)
#' - Binds rows from all top-features TSVs (keeps all original columns)
#' - Writes a single Parquet file **inside** the ML_top_features directory
#'
#' @param path Root results path (same `path` used for createMLResultDir)
#' @param stratify_by NULL | "year" | "country"
#' @param LOO logical; default FALSE
#' @param MDR logical; default FALSE
#' @param cross_test logical; default FALSE
#' @param out_parquet optional filename (no directories). If NULL, defaults to "all_top_features.parquet".
#'                    If a path is given, only its basename is used; the file is written in ML_top_features/.
#' @param compression parquet compression ("zstd" or "snappy"); default "zstd"
#' @param verbose logical; print progress messages
#' @return A tibble with all top-features rows + parsed metadata columns
buildTopFeatsPq <- function(
  path,
  stratify_by = NULL,
  LOO = FALSE,
  MDR = FALSE,
  cross_test = FALSE,
  out_parquet = NULL,   # only filename; will be written under ML_top_features/
  compression = "zstd",
  verbose = TRUE
) {
  # -----------------------
  # Validate inputs
  # -----------------------
  if (!is.character(path) || length(path) != 1 || is.na(path) || nchar(path) == 0) {
    stop("`path` must be a non-empty character scalar.")
  }
  path <- normalizePath(path)

  if (!is.null(stratify_by) && !stratify_by %in% c("year", "country")) {
    stop("`stratify_by` must be NULL, 'year', or 'country'.")
  }
  if (isTRUE(LOO) && is.null(stratify_by)) {
    stop("With LOO=TRUE, stratify_by must be 'year' or 'country'.")
  }
  if (isTRUE(MDR) && (!is.null(stratify_by) || isTRUE(LOO) || isTRUE(cross_test))) {
    stop("MDR can only run when stratify_by = NULL, LOO = FALSE, cross_test = FALSE.")
  }

  # -----------------------
  # Resolve directories (ensures existence)
  # -----------------------
  paths <- createMLResultDir(path, stratify_by = stratify_by, LOO = LOO,
                             cross_test = cross_test, MDR = MDR)
  top_dir <- paths$ML_top_features

  # -----------------------
  # Locate all top-features TSVs
  # -----------------------
  top_files <- list.files(top_dir, pattern = "top_features\\.tsv$", full.names = TRUE, recursive = TRUE)
  if (length(top_files) == 0) {
    if (verbose) message("No *_top_features.tsv files found under: ", top_dir)
    return(tibble::tibble())
  }

  # -----------------------
  # Helpers (shared with performance parser; includes stratified fix)
  # -----------------------
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  .NA_chr <- function() NA_character_

  .find_drug_label_value <- function(tokens) {
    # Finds first "drug" and determines if it's "drug" or "drug_class"
    idx <- which(tokens == "drug")
    if (length(idx) == 0) {
      return(list(drug_label = .NA_chr(), drug_value = .NA_chr(), label_end = NA_integer_))
    }
    i <- idx[1]
    if (i < length(tokens) && identical(tokens[i + 1], "class")) {
      list(drug_label = "drug_class",
           drug_value = if (i + 2 <= length(tokens)) tokens[i + 2] else .NA_chr(),
           label_end  = i + 1)
    } else {
      list(drug_label = "drug",
           drug_value = if (i + 1 <= length(tokens)) tokens[i + 1] else .NA_chr(),
           label_end  = i)
    }
  }

  .parse_base <- function(base_no_suffix) {
    xs <- strsplit(base_no_suffix, "_", fixed = TRUE)[[1]]
    n  <- length(xs)

    # Initialize
    species         <- .NA_chr()
    mdr_tag         <- .NA_chr()
    phenotype       <- .NA_chr()
    drug_label      <- .NA_chr()
    drug_or_class   <- .NA_chr()
    strat_label     <- .NA_chr()
    strat_value     <- .NA_chr()
    strat_value_test<- .NA_chr()
    leaveout        <- FALSE
    is_cross        <- FALSE
    ref_drug        <- .NA_chr()
    test_drug       <- .NA_chr()
    prefix_key      <- .NA_chr()
    feature         <- .NA_chr()
    feature_type    <- .NA_chr()
    feature_subtype <- .NA_chr()

    # Feature from last 2 tokens
    if (n >= 2) {
      feature         <- paste(xs[(n - 1):n], collapse = "_")
      feature_type    <- xs[n - 1]
      feature_subtype <- xs[n]
    }
    core <- if ((n - 2) >= 1) xs[1:(n - 2)] else character(0)
    core_str <- paste(core, collapse = "_")

    # MDR: "MDR_<phenotype>_<feature>"
    if (length(core) > 0 && identical(core[1], "MDR")) {
      mdr_tag <- "MDR"
      if (length(core) >= 2) phenotype <- paste(core[-1], collapse = "_")
      prefix_key <- "MDR"
      return(list(
        species=species, mdr_tag=mdr_tag, phenotype=phenotype,
        drug_label=drug_label, drug_or_class=drug_or_class,
        strat_label=strat_label, strat_value=strat_value, strat_value_test=strat_value_test,
        leaveout=leaveout, is_cross=is_cross,
        ref_drug=ref_drug, test_drug=test_drug,
        prefix_key=prefix_key,
        feature=feature, feature_type=feature_type, feature_subtype=feature_subtype
      ))
    }

    # Cross-test variants
    if (grepl("_tested_on_", core_str, fixed = TRUE)) {
      is_cross <- TRUE

      # LOO cross-test: ... <prefix_key>_<drug_value>_leaveout_tested_on_<strat_value>_<feature>
      if (grepl("_leaveout_tested_on_", core_str, fixed = TRUE)) {
        leaveout <- TRUE
        di <- .find_drug_label_value(core)
        drug_label    <- di$drug_label
        drug_or_class <- di$drug_value
        label_end     <- di$label_end
        if (!is.na(label_end)) {
          prefix_key <- paste(core[1:label_end], collapse = "_")
          i_val <- label_end + 1
          if ((i_val + 4) <= length(core) &&
              core[i_val + 1] == "leaveout" && core[i_val + 2] == "tested" && core[i_val + 3] == "on") {
            strat_label <- stratify_by %||% .NA_chr()
            strat_value <- core[i_val + 4]
          }
        }
        species <- if (length(core) >= 1) core[1] else .NA_chr()
        return(list(
          species=species, mdr_tag=mdr_tag, phenotype=phenotype,
          drug_label=drug_label, drug_or_class=drug_or_class,
          strat_label=strat_label, strat_value=strat_value, strat_value_test=strat_value_test,
          leaveout=leaveout, is_cross=is_cross,
          ref_drug=ref_drug, test_drug=test_drug,
          prefix_key=prefix_key,
          feature=feature, feature_type=feature_type, feature_subtype=feature_subtype
        ))
      }

      # Cross by strat group: ... <prefix_key>_<drug_value>_cross_<strat_value>_tested_on_<strat_value_test>_<feature>
      if (grepl("_cross_", core_str, fixed = TRUE)) {
        di <- .find_drug_label_value(core)
        drug_label    <- di$drug_label
        drug_or_class <- di$drug_value
        label_end     <- di$label_end
        if (!is.na(label_end)) {
          prefix_key <- paste(core[1:label_end], collapse = "_")
          i_val <- label_end + 1
          if ((i_val + 5) <= length(core) &&
              core[i_val + 1] == "cross" && core[i_val + 3] == "tested" && core[i_val + 4] == "on") {
            strat_label     <- stratify_by %||% .NA_chr()
            strat_value     <- core[i_val + 2]
            strat_value_test<- core[i_val + 5]
          }
        }
        species <- if (length(core) >= 1) core[1] else .NA_chr()
        return(list(
          species=species, mdr_tag=mdr_tag, phenotype=phenotype,
          drug_label=drug_label, drug_or_class=drug_or_class,
          strat_label=strat_label, strat_value=strat_value, strat_value_test=strat_value_test,
          leaveout=leaveout, is_cross=is_cross,
          ref_drug=ref_drug, test_drug=test_drug,
          prefix_key=prefix_key,
          feature=feature, feature_type=feature_type, feature_subtype=feature_subtype
        ))
      }

      # Cross-test (drug vs drug): ... <prefix_key>_<ref_drug>_tested_on_<test_drug>_<feature>
      di <- .find_drug_label_value(core)
      label_end <- di$label_end
      if (!is.na(label_end)) {
        prefix_key <- paste(core[1:label_end], collapse = "_")
        if ((label_end + 4) <= length(core) &&
            core[label_end + 2] == "tested" && core[label_end + 3] == "on") {
          ref_drug  <- core[label_end + 1]
          test_drug <- core[label_end + 4]
          drug_label <- if (endsWith(prefix_key, "drug_class")) "drug_class" else "drug"
          drug_or_class <- ref_drug
        }
      }
      species <- if (length(core) >= 1) core[1] else .NA_chr()
      return(list(
        species=species, mdr_tag=mdr_tag, phenotype=phenotype,
        drug_label=drug_label, drug_or_class=drug_or_class,
        strat_label=strat_label, strat_value=strat_value, strat_value_test=strat_value_test,
        leaveout=leaveout, is_cross=is_cross,
        ref_drug=ref_drug, test_drug=test_drug,
        prefix_key=prefix_key,
        feature=feature, feature_type=feature_type, feature_subtype=feature_subtype
      ))
    }

    # Non-cross variants (un/stratified), with stratified FIX
    species <- if (length(core) >= 1) core[1] else .NA_chr()
    di <- .find_drug_label_value(core)
    drug_label <- di$drug_label
    label_end  <- di$label_end

    if (!is.na(label_end)) {
      # STRATIFIED pattern (fix): species drug[_class] strat_label drug_value strat_value
      if (label_end + 3 <= length(core) && core[label_end + 1] %in% c("year", "country")) {
        strat_label  <- core[label_end + 1]
        drug_or_class<- core[label_end + 2]      # e.g., FLQ, MAC, CIP, etc.
        strat_value  <- core[label_end + 3]      # e.g., 2015-2019 or country name
        prefix_key   <- paste(core[1:label_end], collapse = "_")
      } else {
        # UNSTRATIFIED pattern: species drug[_class] drug_value
        drug_or_class <- di$drug_value
        prefix_key    <- paste(core[1:label_end], collapse = "_")
      }
    }

    list(
      species=species, mdr_tag=mdr_tag, phenotype=phenotype,
      drug_label=drug_label, drug_or_class=drug_or_class,
      strat_label=strat_label, strat_value=strat_value, strat_value_test=strat_value_test,
      leaveout=leaveout, is_cross=is_cross,
      ref_drug=ref_drug, test_drug=test_drug,
      prefix_key=prefix_key,
      feature=feature, feature_type=feature_type, feature_subtype=feature_subtype
    )
  }

  # -----------------------
  # Read, parse and bind
  # -----------------------
  out <- purrr::map_dfr(top_files, function(f) {
    base <- basename(f)
    # Accept either "..._top_features.tsv" (preferred) or "<prefix>.tsv" fallback
    base_no_suffix <- sub("_top_features\\.tsv$", "", base)
    if (identical(base_no_suffix, base)) {
      base_no_suffix <- sub("\\.tsv$", "", base)
    }

    meta <- .parse_base(base_no_suffix)

    df <- tryCatch(
      readr::read_tsv(f, show_col_types = FALSE, progress = FALSE),
      error = function(e) {
        if (verbose) message("Failed to read TSV: ", f, " (", conditionMessage(e), ") — using metadata only.")
        tibble()
      }
    )

    # Attach metadata columns
    md_cols <- tibble(

      output_prefix   = base_no_suffix,
      species         = meta$species,
      mdr_tag         = meta$mdr_tag,
      phenotype       = meta$phenotype,
      drug_label      = meta$drug_label,
      drug_or_class   = meta$drug_or_class,
      strat_label     = meta$strat_label,
      strat_value     = meta$strat_value,
      strat_value_test= meta$strat_value_test,
      leaveout        = meta$leaveout,
      cross_test      = meta$is_cross,
      ref_drug        = meta$ref_drug,
      test_drug       = meta$test_drug,
      prefix_key      = meta$prefix_key,
      feature         = meta$feature,
      feature_type    = meta$feature_type,
      feature_subtype = meta$feature_subtype
    )

    if (nrow(df) == 0) {
      md_cols
    } else {
      dplyr::bind_cols(md_cols[rep(1, nrow(df)), ], df)
    }
  })

  # -----------------------
  # Compute the output Parquet path (ALWAYS inside top_dir)
  # -----------------------
  parquet_name <- (out_parquet %||% "all_top_features.parquet")
  parquet_name <- basename(parquet_name)
  out_path <- file.path(top_dir, parquet_name)

  # -----------------------
  # Write Parquet
  # -----------------------
  suppressPackageStartupMessages(library(arrow))
  arrow::write_parquet(out, out_path, compression = compression)
  if (verbose) message("Wrote merged Parquet: ", out_path, "  [", nrow(out), " rows]")

  out
}
