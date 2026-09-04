#' Parse ML result filename into structured metadata
#'
#' Extracts metadata fields encoded in machine learning result filenames.
#' Supports shuffled runs, drug vs drug_class, feature types, and seed
#' information.
#'
#' This helper handles only unstratified result filenames. Year/country
#' stratified filenames carry extra tokens and are aggregated by
#' [buildPerfPqYearCountry()] instead; passing one here raises an error
#' rather than returning silently wrong metadata.
#'
#' @param filename Character. Filename (not full path) to parse.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{shuffled}{Logical; whether dataset was shuffled}
#'   \item{species}{Species name}
#'   \item{drug_label}{"drug" or "drug_class"}
#'   \item{drug_or_class}{Drug name or class}
#'   \item{feature_type}{Feature category}
#'   \item{feature_subtype}{binary/counts}
#'   \item{seed}{Integer seed}
#' }
#'
#' @examples
#' parse_ml_filename("Csp_drug_AMX_genes_binary_42_top_features.parquet")
#' parse_ml_filename("Csp_drug_class_AMG_genes_binary_42_performance.parquet")
#'
#' @export
parse_ml_filename <- function(filename) {
  # Strip known suffix
  base <- sub("_(top_features|performance)\\.parquet$", "", filename)
  xs <- strsplit(base, "_", fixed = TRUE)[[1]]

  out <- list(
    shuffled = FALSE,
    species = NA_character_,
    drug_label = NA_character_,
    drug_or_class = NA_character_,
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
      i <- i + 2
    }
    # Case B: simple drug
    else {
      out$drug_label <- "drug"
      i <- i + 1
    }
  } else {
    stop("ERROR: expected 'drug' token after species")
  }

  # ---------------------------
  # 4. Reject stratified filenames
  #
  # In stratified filenames the strat label ("year"/"country") sits right
  # after "drug"/"drug_class", e.g. "..._drug_year_AMX_2010-2015_...".
  # Those carry extra tokens this parser does not account for, so bail out
  # with a clear message instead of silently mislabelling drug and seed.
  # ---------------------------
  if (i <= length(xs) && xs[i] %in% c("year", "country")) {
    stop(
      "parse_ml_filename() does not support the '", xs[i],
      "' stratified filename '", filename,
      "'; use buildPerfPqYearCountry() for stratified results."
    )
  }

  # ---------------------------
  # 5. Drug or class value
  # ---------------------------
  out$drug_or_class <- xs[i]
  i <- i + 1

  # ---------------------------
  # 6. Feature type + subtype
  # ---------------------------
  out$feature_type <- xs[i]
  i <- i + 1
  out$feature_subtype <- xs[i]
  i <- i + 1

  # ---------------------------
  # 7. Seed
  # ---------------------------
  out$seed <- as.integer(xs[i])

  return(out)
}

#' Build Parquet file of top ML features
#'
#' Reads all `_top_features.parquet` files, parses metadata from filenames,
#' combines them into a single table, and writes a Parquet file.
#'
#' @param top_feat_dir_path The directory containing ML top features files
#' @param out_parquet Output file name
#' @param compression Compression method (default: "zstd")
#' @param verbose Logical; print progress
#'
#' @return A tibble with metadata columns + feature importance data
#'
#' @examples
#' \dontrun{
#' buildTopFeatsPq("data/Campylobacter/ML_top_features/")
#' }
#'
#' @export
buildTopFeatsPq <- function(
  top_feat_dir_path,
  # stratify_by = NULL,
  # LOO = FALSE,
  # MDR = FALSE,
  # cross_test = FALSE,
  out_parquet = "all_top_features.parquet",
  compression = "zstd",
  verbose = TRUE
) {
  if (!is.character(top_feat_dir_path) || length(top_feat_dir_path) != 1 || is.na(top_feat_dir_path) || nchar(top_feat_dir_path) == 0) {
    stop("`path` must be a non-empty character scalar.")
  }
  top_feat_dir_path <- normalizePath(top_feat_dir_path)

  # if (!is.null(stratify_by) && !stratify_by %in% c("year", "country")) {
  #   stop("`stratify_by` must be NULL, 'year', or 'country'.")
  # }
  # if (isTRUE(LOO) && is.null(stratify_by)) {
  #   stop("With LOO=TRUE, stratify_by must be 'year' or 'country'.")
  # }
  # if (isTRUE(MDR) && (!is.null(stratify_by) || isTRUE(LOO) || isTRUE(cross_test))) {
  #   stop("MDR can only run when stratify_by = NULL, LOO = FALSE, cross_test = FALSE.")
  # }

  # -----------------------
  # Resolve directories (ensures existence)
  # -----------------------
  # paths <- createMLResultDir(path,
  #   stratify_by = NULL,
  # LOO = FALSE,
  # MDR = FALSE,
  # cross_test = FALSE
  # )
  # top_dir <- paths$ML_top_features

  files <- list.files(
  top_feat_dir_path,
  pattern = "_top_features\\.parquet$",
  full.names = TRUE,
  recursive = TRUE
)
files <- files[basename(files) != "all_top_features.parquet"]
  
if (length(files) == 0) {
  stop(
    "No files matching '_top_features.parquet' were found in: ",
    normalizePath(top_feat_dir_path, mustWork = FALSE)
  )
}

  if (!length(files)) {
    return(tibble::tibble())
  }

  out <- purrr::map_dfr(files, function(f) {
  meta <- parse_ml_filename(basename(f))

  df <- arrow::read_parquet(f)

  meta_tbl <- tibble::as_tibble(meta)
  dplyr::bind_cols(meta_tbl[rep(1, nrow(df)), ], df)
})

  out_path <- file.path(top_feat_dir_path, basename(out_parquet))
  arrow::write_parquet(out, out_path, compression = compression)

  if (verbose) message("Wrote ", out_path)
  out
}

#' Build Parquet file of ML performance results
#'
#' Reads all `_performance.parquet` files, parses metadata from filenames,
#' combines them into a single table, and writes a Parquet output.
#'
#' @param perf_dir_path The directory containing ML performance files
#' @param out_parquet Output file name
#' @param compression Compression method (default: "zstd")
#' @param verbose Logical; print progress
#'
#' @return A tibble with metadata columns + performance metrics
#'
#' @examples
#' \dontrun{
#' buildPerfPq("data/Campylobacter/ML_performance/")
#' }
#'
#' @export
buildPerfPq <- function(
  perf_dir_path,
  # stratify_by = NULL,
  # LOO = FALSE,
  # MDR = FALSE,
  # cross_test = FALSE,
  out_parquet = "all_perf.parquet",
  compression = "zstd",
  verbose = TRUE
) {
  if (!is.character(perf_dir_path) || length(perf_dir_path) != 1 || is.na(perf_dir_path) || nchar(perf_dir_path) == 0) {
    stop("`perf_dir_path` must be a non-empty character scalar.")
  }
  perf_dir_path <- normalizePath(perf_dir_path)

  # if (!is.null(stratify_by) && !stratify_by %in% c("year", "country")) {
  #   stop("`stratify_by` must be NULL, 'year', or 'country'.")
  # }
  # if (isTRUE(LOO) && is.null(stratify_by)) {
  #   stop("With LOO=TRUE, stratify_by must be 'year' or 'country'.")
  # }
  # if (isTRUE(MDR) && (!is.null(stratify_by) || isTRUE(LOO) || isTRUE(cross_test))) {
  #   stop("MDR can only run when stratify_by = NULL, LOO = FALSE, cross_test = FALSE.")
  # }

  # -----------------------
  # Resolve directories from your function (ensures they exist)
  # -----------------------
  # paths <- createMLResultDir(path,
  #   stratify_by = NULL, LOO = FALSE,
  #   cross_test = FALSE, MDR = FALSE
  # )
  # perf_dir <- paths$ML_performance

  files <- list.files(
    perf_dir_path,
    pattern = "_performance\\.parquet$",
    full.names = TRUE,
    recursive = TRUE
  )

if (length(files) == 0) {
  stop(
    "No files matching '_performance.parquet' were found in: ",
    normalizePath(perf_dir_path, mustWork = FALSE)
  )
}
  
  if (!length(files)) {
    return(tibble::tibble())
  }

  out <- purrr::map_dfr(files, function(f) {
    meta <- parse_ml_filename(basename(f))
    df <- arrow::read_parquet(f)

    # ---- SAFETY FIX ----
    # If TSV already has seed, drop parsed seed
    if ("seed" %in% names(df)) {
      meta$seed <- NULL
    }
    # --------------------

    meta_tbl <- tibble::as_tibble(meta)
    dplyr::bind_cols(meta_tbl[rep(1, nrow(df)), ], df)
  })

  out_path <- file.path(perf_dir_path, basename(out_parquet))
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
#' @examples
#' \dontrun{
#' # Year stratified
#' buildPerfPqYearCountry(perf_dir_path = "data/Campylobacter/ML_year_performance/")
#' # Country stratified
#' buildPerfPqYearCountry(perf_dir_path = "data/Campylobacter/ML_country_performance/")
#' }
#' @export
buildPerfPqYearCountry <- function(
  perf_dir_path
) {
  files <- list.files(perf_dir_path, pattern = "\\.parquet$", full.names = TRUE)

  # Read and combine
  merged_df <- files |>
    purrr::set_names() |> # keeps file names attached
    purrr::map_dfr(~ arrow::read_parquet(.x) |>
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
        # "strat_label2",
        "seed_from_name"
      ),
      regex = "^([^_]+)_((?:[^_]+(?:_[^_]+)?))_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_performance\\.parquet$",
      remove = FALSE
    ) |>
    dplyr::select(-c("seed_from_name"))

  strat_label <- merged_df |>
    dplyr::distinct(strat_label) |>
    dplyr::pull()

  arrow::write_parquet(merged_df, file.path(perf_dir_path, paste0(strat_label, "_perf.parquet")))
  merged_df
}

#' Build Parquet file from cross drug testing ML performances
#'
#' @param perf_dir_path Directory containing performance TSV files
#' @return Writes a Parquet file to the same directory
#'
#' @examples
#' \dontrun{
#' # Cross testing
#' buildPerfPqCrossDrug(perf_dir_path = "data/Campylobacter/cross_test_ML_performance/")
#' }
#'
#' @export
buildPerfPqCrossDrug <- function(
  perf_dir_path
) {
  files <- list.files(perf_dir_path, pattern = "\\.parquet$", full.names = TRUE)

  # Read and combine
  merged_df <- files |>
    purrr::set_names() |> # keeps file names attached
    purrr::map_dfr(~ arrow::read_parquet(.x) |>
      dplyr::mutate(filename = basename(.x))) |>
    tidyr::extract(
      filename,
      into = c("species", "drug", "test_drug", "feature_type", "feature_subtype"),
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
#' \dontrun{
#' buildPerfPqCrossYear(perf_dir_path = "data/Campylobacter/cross_test_ML_year_performance/")
#' }
#'
#' @export
buildPerfPqCrossYear <- function(
  perf_dir_path
) {
  files <- list.files(perf_dir_path, pattern = "\\.parquet$", full.names = TRUE)

  # Read and combine
  merged_df <- files |>
    purrr::set_names() |> # keeps file names attached
    purrr::map_dfr(~ arrow::read_parquet(.x) |>
      dplyr::mutate(filename = basename(.x))) |>
    tidyr::extract(
      filename,
      into = c("species", "drug_label", "drug_or_class", "strat_value", "strat_value_test", "feature_type", "feature_subtype"),
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
#' \dontrun{
#' buildPerfPqCrossCountry(perf_dir_path = "data/Campylobacter/cross_test_ML_country_performance/")
#' }
#'
#' @export
buildPerfPqCrossCountry <- function(
  perf_dir_path
) {
  files <- list.files(perf_dir_path, pattern = "\\.parquet$", full.names = TRUE)

  # Read and combine
  merged_df <- files |>
    purrr::set_names() |> # keeps file names attached
    purrr::map_dfr(~ arrow::read_parquet(.x) |>
      dplyr::mutate(filename = basename(.x))) |>
    tidyr::extract(
      filename,
      into = c("species", "drug_label", "drug_or_class", "strat_value", "strat_value_test", "feature_type", "feature_subtype"),
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
#' \dontrun{
#' # LOO drug
#' buildPerfPqLOODrug(perf_dir_path = "data/Campylobacter/LOO_ML_performance/")
#' }
#'
#' @export
buildPerfPqLOODrug <- function(
  perf_dir_path
) {
  files <- list.files(perf_dir_path, pattern = "\\.parquet$", full.names = TRUE)

  # Read and combine
  merged_df <- files |>
    purrr::set_names() |> # keeps file names attached
    purrr::map_dfr(~ arrow::read_parquet(.x) |>
      dplyr::mutate(filename = basename(.x))) |>
    tidyr::extract(
      filename,
      into = c("species", "LOO_drug", "feature_type", "feature_subtype"),
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
#' \dontrun{
#' # Year stratified
#' buildTopFeatsPqYearCountry(top_feat_dir_path = "data/Campylobacter/ML_year_top_features/")
#'
#' # Country stratified
#' buildTopFeatsPqYearCountry(top_feat_dir_path = "data/Campylobacter/ML_country_top_features/")
#' }
#'
#' @export
buildTopFeatsPqYearCountry <- function(
  top_feat_dir_path
) {
  files <- list.files(top_feat_dir_path, pattern = "\\.parquet$", full.names = TRUE)
files <- files[!basename(files) %in% c("year_top_features.parquet", "country_top_features.parquet")]
  
  # Read and combine
  merged_df <- files |>
    purrr::set_names() |>
    purrr::map_dfr(~ arrow::read_parquet(.x) |>
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
        # "strat_label2",
        "seed_from_name"
      ),
      regex = "^([^_]+)_((?:[^_]+(?:_[^_]+)?))_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_top_features\\.parquet$",
      remove = FALSE
    ) 
    

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
#' \dontrun{
#' # LOO drug
#' buildTopFeatsPqLOODrug(top_feat_dir_path = "data/Campylobacter/LOO_ML_top_features/")
#' }
#'
#' @export
buildTopFeatsPqLOODrug <- function(
  top_feat_dir_path
) {
  files <- list.files(top_feat_dir_path, pattern = "\\.parquet$", full.names = TRUE)
files <- files[!basename(files) %in% c("LOO_drug_top_features.parquet")]
  
  # Read and combine
  merged_df <- files |>
    purrr::set_names() |> # keeps file names attached
    purrr::map_dfr(~ arrow::read_parquet(.x) |>
      dplyr::mutate(filename = basename(.x))) |>
    tidyr::extract(
      filename,
      into = c("species", "LOO_drug", "feature_type", "feature_subtype"),
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
#' \dontrun{
#' # MDR
#' buildTopFeatsPqMDR(top_feat_dir_path = "data/Campylobacter/MDR_ML_top_features/")
#' }
#'
#' @export
buildTopFeatsPqMDR <- function(
    top_feat_dir_path) {
  
  files <- list.files(top_feat_dir_path, pattern = "\\.parquet$", full.names = TRUE)
  files <- files[!basename(files) %in% c("MDR_top_features.parquet")]
  
  # Read and combine
  merged_df <- files|>
    purrr::set_names() |>   # keeps file names attached
    purrr::map_dfr(~ arrow::read_parquet(.x) |>
                     dplyr::mutate(filename = basename(.x))) |>
     tidyr::extract(
      filename,
      into = c("feature_type", "feature_subtype", "seed"),
      regex = "classes_([a-z]+)_(binary|counts)_([0-9]+)"
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
#' \dontrun{
#' # MDR
#' buildPerfPqMDR(perf_dir_path = "data/Campylobacter/MDR_ML_performance/")
#' }
#'
#' @export
buildPerfPqMDR <- function(
    perf_dir_path) {
  
  files <- list.files(perf_dir_path, pattern = "\\.parquet$", full.names = TRUE)
  
  # Read and combine
  merged_df <- files|>
    purrr::set_names() |>   # keeps file names attached
    purrr::map_dfr(~ arrow::read_parquet(.x) |>
                     dplyr::mutate(filename = basename(.x))) |>
    tidyr::extract(
      filename,
      into = c("feature_type", "feature_subtype", "seed"),
      regex = "classes_([a-z]+)_(binary|counts)_([0-9]+)"
    )
  
  arrow::write_parquet(merged_df, file.path(perf_dir_path, "MDR_perf.parquet"))
  
}

#' Build Parquet file from MDR ML predictions
#'
#'
#' @param pred_dir_path Directory containing prediction TSV files
#'
#' @return Writes a Parquet file to the same directory
#'
#' @examples
#' \dontrun{
#' # MDR
#' buildPredPqMDR(pred_dir_path = "data/Campylobacter/MDR_ML_pred/")
#' }
#'
#' @export
buildPredPqMDR <- function(
    pred_dir_path) {
  
  files <- list.files(pred_dir_path, pattern = "\\.parquet$", full.names = TRUE)
  
  # Read and combine
  merged_df <- files|>
    purrr::set_names() |>   # keeps file names attached
    purrr::map_dfr(~ arrow::read_parquet(.x) |>
                  dplyr::select(1:resistant_classes)|>
                     dplyr::mutate(filename = basename(.x)) ) |>
   tidyr::extract(
      filename,
      into = c("feature_type", "feature_subtype", "seed"),
      regex = "classes_([a-z]+)_(binary|counts)_([0-9]+)"
    )
  
  arrow::write_parquet(merged_df, file.path(pred_dir_path, "MDR_pred.parquet"))
  
}
#' Build Parquet file from LOO country ML performances
#'
#' @param perf_dir_path Directory containing performance TSV files
#'
#' @return Writes a Parquet file to the same directory
#'
#' @examples
#' \dontrun{
#' # LOO country
#' buildPerfPqLOOCountry(perf_dir_path = "data/Campylobacter/LOO_ML_country_performance/")
#' }
#'
#' @export
buildPerfPqLOOCountry <- function(
  perf_dir_path
) {
  files <- list.files(perf_dir_path, pattern = "\\.parquet$", full.names = TRUE)

  # Read and combine
  merged_df <- files |>
    purrr::set_names() |> # keeps file names attached
    purrr::map_dfr(~ arrow::read_parquet(.x) |>
      dplyr::mutate(filename = basename(.x))) |>
    tidyr::extract(
      filename,
      into = c("species", "drug_label", "drug_or_class", "leaveout_strat_value", "feature_type", "feature_subtype"),
      regex = "LOO_([A-Za-z0-9-]+)_(drug|drug_class)_country_([A-Za-z0-9-]+)_leaveout_([A-Z]{3})_([a-z]+)_(binary|counts)"
    )

  arrow::write_parquet(merged_df, file.path(perf_dir_path, "LOO_country_perf.parquet"))
}

#' Build Parquet file from LOO year ML performances
#' @param perf_dir_path Directory containing performance TSV files
#'
#' @return Writes a Parquet file to the same directory
#'
#' @examples
#' \dontrun{
#' # LOO year
#' buildPerfPqLOOYear(perf_dir_path = "data/Campylobacter/LOO_ML_year_performance/")
#' }
#'
#' @export
buildPerfPqLOOYear <- function(
  perf_dir_path
) {
  files <- list.files(perf_dir_path, pattern = "\\.parquet$", full.names = TRUE)

  # Read and combine
  merged_df <- files |>
    purrr::set_names() |> # keeps file names attached
    purrr::map_dfr(~ arrow::read_parquet(.x) |>
      dplyr::mutate(filename = basename(.x))) |>
    tidyr::extract(
      filename,
      into = c("species", "drug_label", "drug_or_class", "leaveout_strat_value", "feature_type", "feature_subtype"),
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
#' \dontrun{
#'
#' # LOO country
#' buildTopFeatsPqLOOCountry(top_feat_dir_path = "data/Campylobacter/LOO_ML_country_top_features/")
#' }
#'
#' @export
buildTopFeatsPqLOOCountry <- function(
  top_feat_dir_path
) {
  files <- list.files(top_feat_dir_path, pattern = "\\.parquet$", full.names = TRUE)
  files <- files[files != file.path(top_feat_dir_path, "LOO_country_top_features.parquet")]
  
  # Read and combine
  merged_df <- files |>
    purrr::set_names() |> # keeps file names attached
    purrr::map_dfr(~ arrow::read_parquet(.x) |>
      dplyr::mutate(filename = basename(.x))) |>
    tidyr::extract(
      filename,
      into = c("species", "drug_label", "drug_or_class", "leaveout_strat_value", "feature_type", "feature_subtype"),
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
#' \dontrun{
#'
#' # LOO year
#' buildTopFeatsPqLOOYear(top_feat_dir_path = "data/Campylobacter/LOO_ML_year_top_features/")
#' }
#'
#' @export
buildTopFeatsPqLOOYear <- function(
  top_feat_dir_path
) {
  files <- list.files(top_feat_dir_path, pattern = "\\.parquet$", full.names = TRUE)
  files <- files[files != file.path(top_feat_dir_path, "LOO_year_top_features.parquet")]
  
  # Read and combine
  merged_df <- files |>
    purrr::set_names() |> # keeps file names attached
    purrr::map_dfr(~ arrow::read_parquet(.x) |>
      dplyr::mutate(filename = basename(.x))) |>
    tidyr::extract(
      filename,
      into = c("species", "drug_label", "drug_or_class", "leaveout_strat_value", "feature_type", "feature_subtype"),
      regex = "LOO_([A-Za-z0-9-]+)_(drug|drug_class)_year_([A-Za-z0-9-]+)_leaveout_([0-9]{4}-[0-9]{4})_([a-z]+)_(binary|counts)"
    )

  arrow::write_parquet(merged_df, file.path(top_feat_dir_path, "LOO_year_top_features.parquet"))
}


#' Merge all ML result files into consolidated parquet files
#'
#' Automatically detects supported ML result directories and
#' builds merged parquet summaries for each.
#'
#' @param path Base analysis directory.
#'
#' @return Invisibly returns NULL.
#'
#' @examples
#' \dontrun{
#' mergeMLresults("data/Campylobacter")
#' }
#'
#' @export
mergeMLresults <- function(path) {

  path <- normalizePath(path)

  # -----------------------
  # Standard ML
  # -----------------------
  if (dir.exists(file.path(path, "ML_performance"))) {
    buildPerfPq(
      perf_dir_path = file.path(path, "ML_performance")
    )
  }

  if (dir.exists(file.path(path, "ML_top_features"))) {
    buildTopFeatsPq(
      top_feat_dir_path = file.path(path, "ML_top_features")
    )
  }

  # -----------------------
  # Year stratified
  # -----------------------
  if (dir.exists(file.path(path, "ML_year_performance"))) {
    buildPerfPqYearCountry(
      perf_dir_path = file.path(path, "ML_year_performance")
    )
  }

  if (dir.exists(file.path(path, "ML_year_top_features"))) {
    buildTopFeatsPqYearCountry(
      top_feat_dir_path = file.path(path, "ML_year_top_features")
    )
  }

  # -----------------------
  # Country stratified
  # -----------------------
  if (dir.exists(file.path(path, "ML_country_performance"))) {
    buildPerfPqYearCountry(
      perf_dir_path = file.path(path, "ML_country_performance")
    )
  }

  if (dir.exists(file.path(path, "ML_country_top_features"))) {
    buildTopFeatsPqYearCountry(
      top_feat_dir_path = file.path(path, "ML_country_top_features")
    )
  }

  # -----------------------
  # Cross-testing
  # -----------------------
  if (dir.exists(file.path(path, "cross_test_ML_performance"))) {
    buildPerfPqCrossDrug(
      perf_dir_path = file.path(path, "cross_test_ML_performance")
    )
  }

  if (dir.exists(file.path(path, "cross_test_ML_year_performance"))) {
    buildPerfPqCrossYear(
      perf_dir_path = file.path(path, "cross_test_ML_year_performance")
    )
  }

  if (dir.exists(file.path(path, "cross_test_ML_country_performance"))) {
    buildPerfPqCrossCountry(
      perf_dir_path = file.path(path, "cross_test_ML_country_performance")
    )
  }

  # -----------------------
  # LOO
  # -----------------------
  if (dir.exists(file.path(path, "LOO_ML_performance"))) {
    buildPerfPqLOODrug(
      perf_dir_path = file.path(path, "LOO_ML_performance")
    )
  }

  if (dir.exists(file.path(path, "LOO_ML_top_features"))) {
    buildTopFeatsPqLOODrug(
      top_feat_dir_path = file.path(path, "LOO_ML_top_features")
    )
  }

  if (dir.exists(file.path(path, "LOO_ML_country_performance"))) {
    buildPerfPqLOOCountry(
      perf_dir_path = file.path(path, "LOO_ML_country_performance")
    )
  }

  if (dir.exists(file.path(path, "LOO_ML_country_top_features"))) {
    buildTopFeatsPqLOOCountry(
      top_feat_dir_path = file.path(path, "LOO_ML_country_top_features")
    )
  }

  if (dir.exists(file.path(path, "LOO_ML_year_performance"))) {
    buildPerfPqLOOYear(
      perf_dir_path = file.path(path, "LOO_ML_year_performance")
    )
  }

  if (dir.exists(file.path(path, "LOO_ML_year_top_features"))) {
    buildTopFeatsPqLOOYear(
      top_feat_dir_path = file.path(path, "LOO_ML_year_top_features")
    )
  }

  # -----------------------
  # MDR
  # -----------------------
  if (dir.exists(file.path(path, "MDR_ML_performance"))) {
    buildPerfPqMDR(
      perf_dir_path = file.path(path, "MDR_ML_performance")
    )
  }

  if (dir.exists(file.path(path, "MDR_ML_top_features"))) {
    buildTopFeatsPqMDR(
      top_feat_dir_path = file.path(path, "MDR_ML_top_features")
    )
  }

  if (dir.exists(file.path(path, "MDR_ML_pred"))) {
    buildPredPqMDR(
      pred_dir_path = file.path(path, "MDR_ML_pred")
    )
  }

  invisible(NULL)
}