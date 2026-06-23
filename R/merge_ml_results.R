#' Parse ML result filename into structured metadata
#'
#' Extracts metadata fields encoded in machine learning result filenames.
#' Supports shuffled runs, drug vs drug_class, optional stratification,
#' feature types, and seed information.
#'
#' @param filename Character. Filename (not full path) to parse.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{shuffled}{Logical; whether dataset was shuffled}
#'   \item{species}{Species name}
#'   \item{drug_label}{"drug" or "drug_class"}
#'   \item{drug_or_class}{Drug name or class}
#'   \item{strat_label}{Stratification type (year/country) or NA}
#'   \item{strat_value}{Stratification value or NA}
#'   \item{feature_type}{Feature category}
#'   \item{feature_subtype}{binary/counts}
#'   \item{seed}{Integer seed}
#' }
#'
#' @examples
#' parse_ml_filename("Csp_drug_AMX_genes_binary_42_top_features.tsv")
#'
#' @export
parse_ml_filename <- function(filename) {
  # Strip known suffix
  base <- sub("_(top_features|performance)\\.tsv$", "", filename)
  xs <- strsplit(base, "_", fixed = TRUE)[[1]]

  out <- list(
    shuffled = FALSE,
    species = NA_character_,
    drug_label = NA_character_,
    drug_or_class = NA_character_,
    strat_label = NA_character_,
    strat_value = NA_character_,
    feature_type = NA_character_,
    feature_subtype = NA_character_,
    seed = NA_integer_
  )

  i <- 1

  # ---------------------------
  # 1. Shuffled detection
  # ---------------------------
  if (xs[i] == "shuffled") {
    out$shuffled <- TRUE
    i <- i + 1
  }

  # ---------------------------
  # 2. Species
  # ---------------------------
  out$species <- xs[i]
  i <- i + 1

  # ---------------------------
  # 3. Drug or drug_class
  # ---------------------------
  if (xs[i] == "drug") {
    # Case A: drug_class
    if (xs[i + 1] == "class") {
      out$drug_label <- "drug_class"
      out$drug_or_class <- xs[i + 2]
      i <- i + 3
    }
    # Case B: simple drug
    else {
      out$drug_label <- "drug"
      out$drug_or_class <- xs[i + 1]
      i <- i + 2
    }
  } else {
    stop("ERROR: expected 'drug' token after species")
  }

  # ---------------------------
  # 4. Stratified?
  # ---------------------------
  if (i <= length(xs) && xs[i] %in% c("year", "country")) {
    out$strat_label <- xs[i]
    i <- i + 1
    out$strat_value <- xs[i]
    i <- i + 1
  }

  # ---------------------------
  # 5. Feature type + subtype
  # ---------------------------
  out$feature_type <- xs[i]
  i <- i + 1
  out$feature_subtype <- xs[i]
  i <- i + 1

  # ---------------------------
  # 6. Trailing strat label (mirror)
  # ---------------------------
  if (!is.na(out$strat_label)) {
    stopifnot(xs[i] == out$strat_label)
    i <- i + 1
  }

  # ---------------------------
  # 7. Seed
  # ---------------------------
  out$seed <- as.integer(xs[i])

  return(out)
}

#' Build Parquet file of top ML features
#'
#' Reads all `_top_features.tsv` files, parses metadata from filenames,
#' combines them into a single table, and writes a Parquet file.
#'
#' @param path Base directory containing ML results
#' @param stratify_by NULL, "year", or "country"
#' @param LOO Logical; leave-one-out analysis
#' @param MDR Logical; MDR mode
#' @param cross_test Logical; cross-testing mode
#' @param out_parquet Output file name
#' @param compression Compression method (default: "zstd")
#' @param verbose Logical; print progress
#'
#' @return A tibble with metadata columns + feature importance data
#'
#' @examples
#' buildTopFeatsPq("data/Campylobacter")
#'
#' @export
buildTopFeatsPq <- function(
  path,
  stratify_by = NULL,
  LOO = FALSE,
  MDR = FALSE,
  cross_test = FALSE,
  out_parquet = "all_top_features.parquet",
  compression = "zstd",
  verbose = TRUE
) {
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
  paths <- createMLResultDir(path,
    stratify_by = stratify_by, LOO = LOO,
    cross_test = cross_test, MDR = MDR
  )
  top_dir <- paths$ML_top_features

  files <- list.files(
    top_dir,
    pattern = "_top_features\\.tsv$",
    full.names = TRUE,
    recursive = TRUE
  )

  if (!length(files)) {
    return(tibble::tibble())
  }

  out <- purrr::map_dfr(files, function(f) {
    meta <- parse_ml_filename(basename(f))
    df <- readr::read_tsv(f,
      col_types = readr::cols(
        Variable   = readr::col_character(),
        Importance = readr::col_double(),
        Sign       = readr::col_character()
      ),
      show_col_types = FALSE,
      na = c("NA", "", "NaN")
    )

    meta_tbl <- tibble::as_tibble(meta)
    dplyr::bind_cols(meta_tbl[rep(1, nrow(df)), ], df)
  })

  out_path <- file.path(top_dir, basename(out_parquet))
  arrow::write_parquet(out, out_path, compression = compression)

  if (verbose) message("Wrote ", out_path)
  out
}

#' Build Parquet file of ML performance results
#'
#' Reads all `_performance.tsv` files, parses metadata from filenames,
#' combines them into a single table, and writes a Parquet output.
#'
#' @inheritParams buildTopFeatsPq
#'
#' @return A tibble with metadata columns + performance metrics
#'
#' @examples
#' buildPerfPq("data/Campylobacter")
#'
#' @export
buildPerfPq <- function(
  path,
  stratify_by = NULL,
  LOO = FALSE,
  MDR = FALSE,
  cross_test = FALSE,
  out_parquet = "all_perf.parquet",
  compression = "zstd",
  verbose = TRUE
) {
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
  paths <- createMLResultDir(path,
    stratify_by = stratify_by, LOO = LOO,
    cross_test = cross_test, MDR = MDR
  )
  perf_dir <- paths$ML_performance

  files <- list.files(
    perf_dir,
    pattern = "_performance\\.tsv$",
    full.names = TRUE,
    recursive = TRUE
  )

  if (!length(files)) {
    return(tibble::tibble())
  }

  out <- purrr::map_dfr(files, function(f) {
    meta <- parse_ml_filename(basename(f))
    df <- readr::read_tsv(f, show_col_types = FALSE)

    # ---- SAFETY FIX ----
    # If TSV already has seed, drop parsed seed
    if ("seed" %in% names(df)) {
      meta$seed <- NULL
    }
    # --------------------

    meta_tbl <- tibble::as_tibble(meta)
    dplyr::bind_cols(meta_tbl[rep(1, nrow(df)), ], df)
  })

  out_path <- file.path(perf_dir, basename(out_parquet))
  arrow::write_parquet(out, out_path, compression = compression)

  if (verbose) message("Wrote ", out_path)
  out
}

#' Build Parquet file from year/country ML performances
#'
#' @param perf_dir_path Directory containing performance TSV files
#'
#' @return Writes a Parquet file to the same directory
#'
#' # Year stratified
#' buildPerfPqYearCountry(perf_dir_path = "data/Campylobacter/ML_year_performance/")
#' # Country stratified
#' buildPerfPqYearCountry(perf_dir_path = "data/Campylobacter/ML_country_performance/")
#' @export
buildPerfPqYearCountry <- function(
  perf_dir_path
) {
  files <- list.files(perf_dir_path, pattern = "\\.tsv$", full.names = TRUE)

  # Read and combine
  merged_df <- files |>
    purrr::set_names() |> # keeps file names attached
    purrr::map_dfr(~ readr::read_tsv(.x, show_col_types = FALSE) |>
      dplyr::mutate(filename = basename(.x))) |>
    tidyr::extract(
      filename,
      into = c(
        "species",
        "drug_label",
        "strat_label",
        "drug_or_class",
        "strat_value",
        "feature_type",
        "feature_subtype",
        "strat_label2",
        "seed_from_name"
      ),
      regex = "^([^_]+)_((?:[^_]+(?:_[^_]+)?))_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_performance\\.tsv$",
      remove = FALSE
    ) |>
    dplyr::select(-c("strat_label2", "seed_from_name"))

  strat_label <- merged_df |>
    dplyr::distinct(strat_label) |>
    dplyr::pull()

  arrow::write_parquet(merged_df, file.path(perf_dir_path, paste0(strat_label, "_perf.parquet")))
}

#' Build Parquet file from cross drug testing ML performances
#'
#' @param perf_dir_path Directory containing performance TSV files
#' @return Writes a Parquet file to the same directory
#'
#' @examples
#' # Cross testing
#' buildPerfPqCrossDrug(perf_dir_path = "data/Campylobacter/cross_test_ML_performance/")
#'
#' @export
buildPerfPqCrossDrug <- function(
  perf_dir_path
) {
  files <- list.files(perf_dir_path, pattern = "\\.tsv$", full.names = TRUE)

  # Read and combine
  merged_df <- files |>
    purrr::set_names() |> # keeps file names attached
    purrr::map_dfr(~ readr::read_tsv(.x, show_col_types = FALSE) |>
      dplyr::mutate(filename = basename(.x))) |>
    tidyr::extract(
      filename,
      into = c("species", "drug_or_class", "tested_on", "feature_type", "feature_subtype"),
      regex = "cross_test_([A-Za-z0-9-]+)_drug_([A-Z]+)_cross_tested_on_([A-Z-]+)_([a-z]+)_(binary|counts)"
    )

  arrow::write_parquet(merged_df, file.path(perf_dir_path, "cross_drug_perf.parquet"))
}

#' Build Parquet file from cross year testing ML performances
#'
#' @param perf_dir_path Directory containing performance TSV files
#'
#' @return Writes a Parquet file to the same directory
#'
#' @examples
#' buildPerfPqCrossYear(perf_dir_path = "data/Campylobacter/cross_test_ML_year_performance/")
#'
#' @export
buildPerfPqCrossYear <- function(
  perf_dir_path
) {
  files <- list.files(perf_dir_path, pattern = "\\.tsv$", full.names = TRUE)

  # Read and combine
  merged_df <- files |>
    purrr::set_names() |> # keeps file names attached
    purrr::map_dfr(~ readr::read_tsv(.x, show_col_types = FALSE) |>
      dplyr::mutate(filename = basename(.x))) |>
    tidyr::extract(
      filename,
      into = c("species", "drug_label", "drug_or_class", "train_year", "test_year", "feature_type", "feature_subtype"),
      regex = "cross_test_([A-Za-z0-9-]+)_(drug|drug_class)_([A-Za-z0-9-]+)_cross_([0-9]{4}-[0-9]{4})_tested_on_([0-9]{4}-[0-9]{4})_([a-z]+)_(binary|counts)"
    )

  arrow::write_parquet(merged_df, file.path(perf_dir_path, "cross_year_perf.parquet"))
}

#' Build Parquet file from cross country testing ML performances
#'
#' @param perf_dir_path Directory containing performance TSV files
#'
#' @return Writes a Parquet file to the same directory
#'
#' @examples
#' buildPerfPqCrossCountry(perf_dir_path = "data/Campylobacter/cross_test_ML_country_performance/")
#'
#' @export
buildPerfPqCrossCountry <- function(
  perf_dir_path
) {
  files <- list.files(perf_dir_path, pattern = "\\.tsv$", full.names = TRUE)

  # Read and combine
  merged_df <- files |>
    purrr::set_names() |> # keeps file names attached
    purrr::map_dfr(~ readr::read_tsv(.x, show_col_types = FALSE) |>
      dplyr::mutate(filename = basename(.x))) |>
    tidyr::extract(
      filename,
      into = c("species", "drug_label", "drug_or_class", "train_country", "test_country", "feature_type", "feature_subtype"),
      regex = "cross_test_([A-Za-z0-9-]+)_(drug|drug_class)_([A-Za-z0-9-]+)_cross_([A-Z]{3})_tested_on_([A-Z]{3})_([a-z]+)_(binary|counts)"
    )

  arrow::write_parquet(merged_df, file.path(perf_dir_path, "cross_country_perf.parquet"))
}

#' Build Parquet file from LOO drug ML performances
#'
#' @param perf_dir_path Directory containing performance TSV files
#'
#' @return Writes a Parquet file to the same directory
#'
#' @examples
#' # LOO drug
#' buildPerfPqLOODrug(perf_dir_path = "data/Campylobacter/LOO_ML_performance/")
#'
#' @export
buildPerfPqLOODrug <- function(
  perf_dir_path
) {
  files <- list.files(perf_dir_path, pattern = "\\.tsv$", full.names = TRUE)

  # Read and combine
  merged_df <- files |>
    purrr::set_names() |> # keeps file names attached
    purrr::map_dfr(~ readr::read_tsv(.x, show_col_types = FALSE) |>
      dplyr::mutate(filename = basename(.x))) |>
    tidyr::extract(
      filename,
      into = c("species", "drug_or_class", "feature_type", "feature_subtype"),
      regex = "LOO_([A-Za-z0-9-]+)_drug_leaveout_([A-Z]+)_([a-z]+)_(binary|counts)"
    )

  arrow::write_parquet(merged_df, file.path(perf_dir_path, "LOO_drug_perf.parquet"))
}

#' Build Parquet file from year/country ML top features
#'
#' @param top_feat_dir_path Directory containing top feature TSV files
#'
#' @return Writes a Parquet file to the same directory
#'
#' @examples
#' # Year stratified
#' buildTopFeatsPqYearCountry(top_feat_dir_path = "data/Campylobacter/ML_year_top_features/")
#'
#' # Country stratified
#' buildTopFeatsPqYearCountry(top_feat_dir_path = "data/Campylobacter/ML_country_top_features/")
#'
#' @export
buildTopFeatsPqYearCountry <- function(
  top_feat_dir_path
) {
  files <- list.files(top_feat_dir_path, pattern = "\\.tsv$", full.names = TRUE)

  # Read and combine
  merged_df <- files |>
    purrr::set_names() |>
    purrr::map_dfr(~ readr::read_tsv(.x, show_col_types = FALSE) |>
      dplyr::mutate(
        Importance = as.numeric(Importance),
        filename = basename(.x)
      )) |>
    tidyr::extract(
      filename,
      into = c(
        "species",
        "drug_label",
        "strat_label",
        "drug_or_class",
        "strat_value",
        "feature_type",
        "feature_subtype",
        "strat_label2",
        "seed_from_name"
      ),
      regex = "^([^_]+)_((?:[^_]+(?:_[^_]+)?))_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_top_features\\.tsv$",
      remove = FALSE
    ) |>
    dplyr::select(-c("strat_label2"))

  strat_label <- merged_df |>
    dplyr::distinct(strat_label) |>
    dplyr::pull()

  arrow::write_parquet(merged_df, file.path(top_feat_dir_path, paste0(strat_label, "_top_features.parquet")))
}

#' Build Parquet file from LOO drug ML top features
#'
#' @param top_feat_dir_path Directory containing top feature TSV files
#'
#' @return Writes a Parquet file to the same directory
#'
#' @examples
#' # LOO drug
#' buildTopFeatsPqLOODrug(top_feat_dir_path = "data/Campylobacter/LOO_ML_top_features/")
#'
#' @export
buildTopFeatsPqLOODrug <- function(
  top_feat_dir_path
) {
  files <- list.files(top_feat_dir_path, pattern = "\\.tsv$", full.names = TRUE)

  # Read and combine
  merged_df <- files |>
    purrr::set_names() |> # keeps file names attached
    purrr::map_dfr(~ readr::read_tsv(.x, show_col_types = FALSE) |>
      dplyr::mutate(filename = basename(.x))) |>
    tidyr::extract(
      filename,
      into = c("species", "drug_or_class", "feature_type", "feature_subtype"),
      regex = "LOO_([A-Za-z0-9-]+)_drug_leaveout_([A-Z]+)_([a-z]+)_(binary|counts)"
    )

  arrow::write_parquet(merged_df, file.path(top_feat_dir_path, "LOO_drug_top_features.parquet"))
}

#' Build Parquet file from MDR ML top features
#'
#' @param top_feat_dir_path Directory containing top feature TSV files
#'
#' @return Writes a Parquet file to the same directory
#'
#' @examples
#' # MDR
#' buildTopFeatsPqMDR(top_feat_dir_path = "data/Campylobacter/MDR_ML_top_features/")
#'
#' @export
buildTopFeatsPqMDR <- function(
  top_feat_dir_path
) {
  files <- list.files(top_feat_dir_path, pattern = "\\.tsv$", full.names = TRUE)

  # Read and combine
  merged_df <- files |>
    purrr::set_names() |> # keeps file names attached
    purrr::map_dfr(~ readr::read_tsv(.x, show_col_types = FALSE) |>
      dplyr::mutate(filename = basename(.x))) |>
    dplyr::mutate(
      feature_type = stringr::str_extract(filename, "args|cogs|genes|domains|proteins|struct"),
      feature_subtype = stringr::str_extract(filename, "binary|counts")
    )

  arrow::write_parquet(merged_df, file.path(top_feat_dir_path, "MDR_top_features.parquet"))
}

#' Build Parquet file from MDR ML performances
#'
#'
#' @param perf_dir_path Directory containing performance TSV files
#'
#' @return Writes a Parquet file to the same directory
#'
#' @examples
#' # MDR
#' buildPerfPqMDR(perf_dir_path = "data/Campylobacter/MDR_ML_performance/")
#'
#' @export
buildPerfPqMDR <- function(
  perf_dir_path
) {
  files <- list.files(perf_dir_path, pattern = "\\.tsv$", full.names = TRUE)

  # Read and combine
  merged_df <- files |>
    purrr::set_names() |> # keeps file names attached
    purrr::map_dfr(~ readr::read_tsv(.x, show_col_types = FALSE) |>
      dplyr::mutate(filename = basename(.x))) |>
    dplyr::mutate(
      feature_type = stringr::str_extract(filename, "args|cogs|genes|domains|proteins|struct"),
      feature_subtype = stringr::str_extract(filename, "binary|counts")
    )

  arrow::write_parquet(merged_df, file.path(perf_dir_path, "MDR_perf.parquet"))
}

#' Build Parquet file from LOO country ML performances
#'
#' @param perf_dir_path Directory containing performance TSV files
#'
#' @return Writes a Parquet file to the same directory
#'
#' @examples
#' # LOO country
#' buildPerfPqLOOCountry(perf_dir_path = "data/Campylobacter/LOO_ML_country_performance/")
#'
#' @export
buildPerfPqLOOCountry <- function(
  perf_dir_path
) {
  files <- list.files(perf_dir_path, pattern = "\\.tsv$", full.names = TRUE)

  # Read and combine
  merged_df <- files |>
    purrr::set_names() |> # keeps file names attached
    purrr::map_dfr(~ readr::read_tsv(.x, show_col_types = FALSE) |>
      dplyr::mutate(filename = basename(.x))) |>
    tidyr::extract(
      filename,
      into = c("species", "drug_label", "drug_or_class", "leaveout_country", "feature_type", "feature_subtype"),
      regex = "LOO_([A-Za-z0-9-]+)_(drug|drug_class)_country_([A-Za-z0-9-]+)_leaveout_([A-Z]{3})_([a-z]+)_(binary|counts)"
    )

  arrow::write_parquet(merged_df, file.path(perf_dir_path, "LOO_country_perf.parquet"))
}

#' Build Parquet file from LOO year ML performances
#' @param perf_dir_path Directory containing performance TSV files
#' @param top_feat_dir_path Directory containing top feature TSV files
#'
#' @return Writes a Parquet file to the same directory
#'
#' @examples
#' # LOO year
#' buildPerfPqLOOYear(perf_dir_path = "data/Campylobacter/LOO_ML_year_performance/")
#'
#' @export
buildPerfPqLOOYear <- function(
  perf_dir_path
) {
  files <- list.files(perf_dir_path, pattern = "\\.tsv$", full.names = TRUE)

  # Read and combine
  merged_df <- files |>
    purrr::set_names() |> # keeps file names attached
    purrr::map_dfr(~ readr::read_tsv(.x, show_col_types = FALSE) |>
      dplyr::mutate(filename = basename(.x))) |>
    tidyr::extract(
      filename,
      into = c("species", "drug_label", "drug_or_class", "leaveout_year", "feature_type", "feature_subtype"),
      regex = "LOO_([A-Za-z0-9-]+)_(drug|drug_class)_year_([A-Za-z0-9-]+)_leaveout_([0-9]{4}-[0-9]{4})_([a-z]+)_(binary|counts)"
    )

  arrow::write_parquet(merged_df, file.path(perf_dir_path, "LOO_year_perf.parquet"))
}

#' Build Parquet file from LOO country ML top features
#'
#' @param top_feat_dir_path Directory containing top feature TSV files
#'
#' @return Writes a Parquet file to the same directory
#'
#' @examples
#'
#' # LOO country
#' buildTopFeatsPqLOOCountry(top_feat_dir_path = "data/Campylobacter/LOO_ML_country_top_features/")
#'
#' @export
buildTopFeatsPqLOOCountry <- function(
  top_feat_dir_path
) {
  files <- list.files(top_feat_dir_path, pattern = "\\.tsv$", full.names = TRUE)

  # Read and combine
  merged_df <- files |>
    purrr::set_names() |> # keeps file names attached
    purrr::map_dfr(~ readr::read_tsv(.x, show_col_types = FALSE) |>
      dplyr::mutate(filename = basename(.x))) |>
    tidyr::extract(
      filename,
      into = c("species", "drug_label", "drug_or_class", "leaveout_country", "feature_type", "feature_subtype"),
      regex = "LOO_([A-Za-z0-9-]+)_(drug|drug_class)_country_([A-Za-z0-9-]+)_leaveout_([A-Z]{3})_([a-z]+)_(binary|counts)"
    )

  arrow::write_parquet(merged_df, file.path(top_feat_dir_path, "LOO_country_top_features.parquet"))
}

#' Build Parquet file from LOO year ML top features
#'
#' @param top_feat_dir_path Directory containing top feature TSV files
#'
#' @return Writes a Parquet file to the same directory
#'
#' @examples
#'
#' # LOO year
#' buildTopFeatsPqLOOYear(top_feat_dir_path = "data/Campylobacter/LOO_ML_year_top_features/")
#'
#' @export
buildTopFeatsPqLOOYear <- function(
  top_feat_dir_path
) {
  files <- list.files(top_feat_dir_path, pattern = "\\.tsv$", full.names = TRUE)

  # Read and combine
  merged_df <- files |>
    purrr::set_names() |> # keeps file names attached
    purrr::map_dfr(~ readr::read_tsv(.x, show_col_types = FALSE) |>
      dplyr::mutate(filename = basename(.x))) |>
    tidyr::extract(
      filename,
      into = c("species", "drug_label", "drug_or_class", "leaveout_year", "feature_type", "feature_subtype"),
      regex = "LOO_([A-Za-z0-9-]+)_(drug|drug_class)_year_([A-Za-z0-9-]+)_leaveout_([0-9]{4}-[0-9]{4})_([a-z]+)_(binary|counts)"
    )

  arrow::write_parquet(merged_df, file.path(top_feat_dir_path, "LOO_year_top_features.parquet"))
}
