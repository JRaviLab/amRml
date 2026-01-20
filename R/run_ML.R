#' .resolveSplitParams()
#' 
#' Pulls the ML parameters .json and reads what model split parameters are to be
#' used. These defaults can be overridden if you so choose, but consider regenerating
#' the ML matrices with these new split/CV values instead.
#'
#' @param parquet_path Path to the output directory created by amR_data::cleanData()
#' @param defaults
#'
#' @returns ML splits read from the parameters .json at the entered output directory
#'
#' @keywords internal
#' @examples NULL
.resolveSplitParams <- function(parquet_path,
                                defaults = list(split = c(0.8, 0),
                                                seed  = 5280,
                                                n_fold = 5)) {
  # matrix_dir is the directory that contains the parquet files
  matrix_dir <- normalizePath(dirname(parquet_path))
  params_json <- .readMLParameters(matrix_dir)

  if (is.null(params_json)) {
    return(defaults)
  }

  list(
    split  = if (!is.null(params_json$split))  params_json$split  else defaults$split,
    seed   = if (!is.null(params_json$seed))   params_json$seed   else defaults$seed,
    n_fold = if (!is.null(params_json$n_fold)) params_json$n_fold else defaults$n_fold
  )
}

#' Create machine learning result directories
#'
#' Creates a structured directory hierarchy for storing machine learning results
#' including matrices, performance metrics, feature importance, models, and predictions.
#'
#' @param path Character scalar. Base directory path where subdirectories will be created.
#' @param stratify_by Character scalar or NULL. Stratification method: \code{"country"},
#'   \code{"year"}, or \code{NULL}. Default is \code{NULL} (no stratification).
#' @param LOO Logical. Whether to create directories for Leave-One-Out analysis.
#'   Default is \code{FALSE}. Only valid when \code{stratify_by} is not \code{NULL}.
#' @param cross_test Logical. Whether to create directories for cross-testing.
#'   Default is \code{FALSE}.
#' @param MDR Logical. Whether to create directories for Multi-Drug Resistance (MDR) analysis.
#'   Default is \code{FALSE}. When \code{TRUE}, all other parameters must be \code{FALSE}/\code{NULL}.
#'
#' @return A named list containing paths to:
#'   \describe{
#'     \item{matrix_path}{Directory for input matrices}
#'     \item{ML_performance}{Directory for performance metrics}
#'     \item{ML_top_features}{Directory for top feature rankings}
#'     \item{ML_models}{Directory for saved model objects}
#'     \item{ML_prediction}{Directory for prediction results}
#'   }
#'
#' @examples
#' \dontrun{
#' # Basic directory structure
#' paths <- createMLResultDir("/path/to/results")
#'
#' # LOO analysis stratified by year
#' paths_loo <- createMLResultDir("/path/to/results",
#'                                 stratify_by = "year",
#'                                 LOO = TRUE)
#'
#' # MDR analysis
#' paths_mdr <- createMLResultDir("/path/to/results", MDR = TRUE)
#' }
#'
#' @export
createMLResultDir <- function(path,
                              stratify_by = NULL,
                              LOO = FALSE,
                              cross_test = FALSE,
                              MDR = FALSE) {
  # Basic input validation
  if (!is.character(path) || length(path) != 1 || is.na(path) || nchar(path) == 0) {
    stop("`path` must be a non-empty character scalar.")
  }
  path <- normalizePath(path)

  if (!is.null(stratify_by) && length(stratify_by) != 1) {
    stop("`stratify_by` must be NULL or a single value: 'country' or 'year'.")
  }

  # MDR mode: mutually exclusive with LOO/cross_test/stratification
  if (MDR) {
    if (!is.null(stratify_by) || isTRUE(LOO) || isTRUE(cross_test)) {
      stop("MDR can only run when stratify_by is NULL and both LOO and cross_test are FALSE.")
    }
    paths <- list(
      matrix_path     = file.path(path, "MDR_matrix"),
      ML_performance  = file.path(path, "MDR_ML_performance"),
      ML_top_features = file.path(path, "MDR_ML_top_features"),
      ML_models       = file.path(path, "MDR_ML_models"),
      ML_prediction   = file.path(path, "MDR_ML_pred")
    )
  } else {
    # Determine prefixes (only in non-MDR mode)
    full_prefix <- paste0(ifelse(isTRUE(LOO), "LOO_", ""),
                          ifelse(isTRUE(cross_test), "cross_test_", ""))
    half_prefix <- ifelse(isTRUE(LOO), "LOO_", "")

    # Determine suffix
    suffix <- if (is.null(stratify_by) || identical(stratify_by, "")) {
      ""
    } else {
      switch(
        stratify_by,
        "country" = "_country",
        "year"    = "_year",
        stop("`stratify_by` must be NULL, 'country', or 'year'.")
      )
    }

    # Build paths
    paths <- list(
      matrix_path     = file.path(path, paste0(half_prefix, "matrix", suffix)),
      ML_performance  = file.path(path, paste0(full_prefix, "ML", suffix, "_performance")),
      ML_top_features = file.path(path, paste0(full_prefix, "ML", suffix, "_top_features")),
      ML_models       = file.path(path, paste0(full_prefix, "ML", suffix, "_models")),
      ML_prediction   = file.path(path, paste0(full_prefix, "ML", suffix, "_pred"))
    )
  }

  # Create directories (ensure matrix_path is created too)
  dirs <- unlist(paths, use.names = FALSE)
  for (d in dirs) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE)
    }
  }

  return(paths)
}

 # createAllMLResultDir <- function(path) {
 #    createMLResultDir(path, stratify_by = NULL, LOO = FALSE, cross_test = FALSE, MDR = FALSE)
 #    createMLResultDir(path, stratify_by = NULL, LOO = FALSE, cross_test = TRUE, MDR = FALSE)
 #    createMLResultDir(path, stratify_by = NULL, LOO = FALSE, cross_test = FALSE, MDR = TRUE)
 #    createMLResultDir(path, stratify_by = "year", LOO = FALSE, cross_test = FALSE, MDR = FALSE)
 #    createMLResultDir(path, stratify_by = "year", LOO = FALSE, cross_test = TRUE, MDR = FALSE)
 #    createMLResultDir(path, stratify_by = "year", LOO = TRUE, cross_test = FALSE, MDR = FALSE)
 #    createMLResultDir(path, stratify_by = "year", LOO = TRUE, cross_test = TRUE, MDR = FALSE)
 #    createMLResultDir(path, stratify_by = "country", LOO = FALSE, cross_test = FALSE, MDR = FALSE)
 #    createMLResultDir(path, stratify_by = "country", LOO = FALSE, cross_test = TRUE, MDR = FALSE)
 #    createMLResultDir(path, stratify_by = "country", LOO = TRUE, cross_test = FALSE, MDR = FALSE)
 #    createMLResultDir(path, stratify_by = "country", LOO = TRUE, cross_test = TRUE, MDR = FALSE)
 #    }
 #

#' Create machine learning input list
#'
#' Parses parquet file names in the matrix directory and generates a tibble
#' mapping input files to output paths. Handles multiple analysis modes including
#' standard, cross-test, Leave-One-Out (LOO), and Multi-Drug Resistance (MDR).
#'
#' @param path Character scalar. Base directory path containing matrix subdirectories.
#' @param stratify_by Character scalar or NULL. Stratification method: \code{"country"},
#'   \code{"year"}, or \code{NULL}.
#' @param LOO Logical. Whether to perform Leave-One-Out analysis. Requires \code{stratify_by}.
#' @param MDR Logical. Whether to perform Multi-Drug Resistance analysis.
#' @param cross_test Logical. Whether to perform cross-testing between groups.
#'
#' @return A tibble with columns:
#'   \describe{
#'     \item{ref_file}{Path to reference/training parquet file}
#'     \item{test_file}{Path to test parquet file (NA for non-cross-test)}
#'     \item{output_prefix}{Prefix for output files}
#'     \item{matrix_path}{Directory containing matrix files}
#'     \item{out_perf}{Directory for performance output}
#'     \item{out_top}{Directory for top features output}
#'     \item{out_models}{Directory for model objects}
#'     \item{out_pred}{Directory for predictions}
#'   }
#'
#' @examples
#' \dontrun{
#' # Standard ML input list
#' inputs <- createMLinputList("/path/to/results")
#'
#' # Cross-test with year stratification
#' inputs_ct <- createMLinputList("/path/to/results",
#'                                stratify_by = "year",
#'                                cross_test = TRUE)
#'
#' # MDR analysis
#' inputs_mdr <- createMLinputList("/path/to/results", MDR = TRUE)
#' }
#'
#' @export
createMLinputList <- function(path,
                              stratify_by = NULL,
                              LOO = FALSE,
                              MDR = FALSE,
                              cross_test = FALSE) {

  # Validate inputs
  if (!is.character(path) || length(path) != 1 || is.na(path))
    stop("`path` must be a valid file path string.")

  path <- normalizePath(path)

  if (isTRUE(LOO) && (is.null(stratify_by) || !(stratify_by %in% c("year", "country")))) {
    stop("For Leave-One-Out (LOO) models, stratify_by must be 'year' or 'country'.")
  }

  if (isTRUE(MDR) && (!is.null(stratify_by) || LOO || cross_test)) {
    stop("MDR can only run when stratify_by = NULL, LOO = FALSE, cross_test = FALSE.")
  }

  # Create results directories
  paths <- createMLResultDir(path, stratify_by, LOO, cross_test, MDR)

  # List ML matrix parquet files
  files_vec <- list.files(paths$matrix_path, pattern = "\\.parquet$", full.names = TRUE)

  if (length(files_vec) == 0) {
    message("No parquet files found in matrix path: ", paths$matrix_path)
    return(tibble::tibble())
  }

  # Get index of token
  .get_idx <- function(v, token) {
    w <- which(v == token)
    if (length(w) == 0) NA_integer_ else w[1]
  }

  # Parsing tokens from filename

  # ============================
  # Multi-drug resistance models
  # ============================
  if (MDR) {

    parsed <- tibble::tibble(ref_file = files_vec) |>
      dplyr::mutate(
        parts = stringr::str_split(basename(ref_file), "_"),

        species = purrr::map_chr(parts, ~ .x[1]),
        mdr_tag = purrr::map_chr(parts, ~ .x[2]),   # always "MDR"
        phenotype = purrr::map_chr(parts, ~ paste(.x[3:4], collapse = "_")),

        # Feature is 5th + 6th tokens
        feature_type = purrr::map_chr(parts, ~ .x[5]),
        feature_subtype = purrr::map_chr(parts, ~ stringr::str_remove(.x[6], "_sparse.parquet")),

        feature = purrr::map2_chr(feature_type, feature_subtype, paste, sep = "_"),

        output_prefix = paste0("MDR_", phenotype, "_", feature)
      )

    out <- parsed |>
      dplyr::mutate(
        test_file = NA_character_,
        matrix_path = paths$matrix_path,
        out_perf  = paths$ML_performance,
        out_top   = paths$ML_top_features,
        out_models= paths$ML_models,
        out_pred  = paths$ML_prediction
      )

    return(out)

  # ============================
  # For all other modeling types
  # ============================
  } else {

    parsed <- tibble::tibble(ref_file = files_vec) |>
      dplyr::mutate(
        parts    = stringr::str_split(basename(ref_file), "_"),
        i_sparse = purrr::map_int(parts, ~ .get_idx(.x, "sparse.parquet")),
        i_strat  = purrr::map_int(parts, ~ {
          if (is.null(stratify_by)) return(NA_integer_)
          .get_idx(.x, stratify_by)
        }),

        # Feature = last two tokens before sparse.parquet
        feature = purrr::map2_chr(parts, i_sparse, ~ {
          i <- .y; x <- .x
          if (is.na(i) || i < 3) return(NA_character_)
          paste(x[(i - 2):(i - 1)], collapse = "_")
        }),

        # Drug or drug class extraction
        drug_or_class = purrr::map2_chr(parts, i_strat, ~ {
          i <- .y; x <- .x

          # Stratified models
          if (!is.na(i)) {
            return(x[i + 1])
          }

          # Unstratified models
          # Find where "drug" or "drug_class" ends

          if (x[2] == "drug" && x[3] != "class") {
            # Case A: Cje_drug_X
            return(x[3])
          }

          if (x[2] == "drug" && x[3] == "class") {
            # Case B: Cje_drug_class_X
            return(x[4])
          }

          # Fallback for unexpected structures
          NA_character_
        }),

        # Stratification value (if present)
        strat_value = purrr::map2_chr(parts, i_strat, ~ {
          i <- .y; x <- .x
          if (is.na(i)) return("")
          # default position is two tokens after the strat label
          j <- i + 2
          # if there's an intervening 'leaveout', skip over it
          if (j <= length(x) && identical(x[j], "leaveout")) j <- j + 1
          if (j <= length(x)) return(x[j])
          ""   # no stratification
        }),

        # Prefix key for grouping
        prefix_key = purrr::map2_chr(parts, i_strat, ~ {
          i <- .y; x <- .x

          # Case A: stratified -> prefix before the stratify label
          if (!is.na(i)) {
            if (i - 1 >= 1) return(paste(x[1:(i - 1)], collapse = "_"))
            return("")
          }

          # Case B: unstratified -> prefix is first two tokens
          if (x[2] == "drug" && x[3] != "class"){
            # Case A: Cje_drug_X
            return(paste(x[1:2], collapse = "_"))
          }
          if (x[2] == "drug" && x[3] == "class"){
            # Case A: Cje_drug_X
            return(paste(x[1:3], collapse = "_"))
          }
        })
      )

    # ============================
    # Non-cross-test modeling
    # ============================
    if (!MDR && !cross_test) {
      out <- parsed |>
        dplyr::mutate(
          test_file = NA_character_,
          output_prefix = gsub("_sparse\\.parquet$", "", basename(ref_file)),
          matrix_path = paths$matrix_path,
          out_perf  = paths$ML_performance,
          out_top   = paths$ML_top_features,
          out_models= paths$ML_models,
          out_pred  = paths$ML_prediction
        )
      return(out)

    # ============================
    # Cross-test modeling, no LOO
    # ============================
    } else if (cross_test && !LOO) {

      if (is.null(stratify_by)) {
        # Case A: stratify_by = NULL, pair across abx within same feature + prefix
        pairs <- parsed |>
          dplyr::select(ref_file, feature, prefix_key, strat_value, ref_drug = drug_or_class) |>
          dplyr::inner_join(
            parsed |>
              dplyr::select(test_file = ref_file, feature, prefix_key, strat_value, test_drug = drug_or_class),
            by = c("feature", "prefix_key", "strat_value")
          ) |>
          dplyr::filter(ref_file != test_file,
                        ref_drug != test_drug) |>
          dplyr::distinct() |>
          dplyr::mutate(
            output_prefix = paste0(
              prefix_key, "_",
              ref_drug, "_tested_on_", test_drug, "_",
              feature
            )
          )

        out <- pairs |>
          dplyr::mutate(
            matrix_path = paths$matrix_path,
            out_perf  = paths$ML_performance,
            out_top   = paths$ML_top_features,
            out_models= paths$ML_models,
            out_pred  = paths$ML_prediction
          )

        return(out)
      }

      # Case B: stratify_by != NULL, pair same drug/class, prefix, feature,
      # but across different stratification groups
      pairs <- parsed |>
        dplyr::select(ref_file, feature, prefix_key, strat_value,
                      drug_or_class) |>

        # self-join ONLY on prefix_key, drug/class, feature
        dplyr::inner_join(
          parsed |>
            dplyr::select(test_file = ref_file,
                          feature, prefix_key, strat_value_test = strat_value,
                          drug_or_class),
          by = c("prefix_key", "feature", "drug_or_class")
        ) |>

        # do NOT test file against itself
        dplyr::filter(ref_file != test_file) |>

        # enforce different stratification group
        dplyr::filter(strat_value != strat_value_test) |>

        # remove symmetric duplicates (A,B == B,A)
        dplyr::rowwise() |>
        dplyr::mutate(pair_id = paste(sort(c(ref_file, test_file)), collapse = "||")) |>
        dplyr::ungroup() |>
        dplyr::distinct(pair_id, .keep_all = TRUE) |>

        dplyr::mutate(
          output_prefix = paste0(
            prefix_key, "_",
            drug_or_class, "_cross_",
            strat_value, "_tested_on_", strat_value_test, "_",
            feature
          )
        ) |>
        dplyr::select(ref_file, test_file, output_prefix)

      out <- pairs |>
        dplyr::mutate(
          matrix_path = paths$matrix_path,
          out_perf  = paths$ML_performance,
          out_top   = paths$ML_top_features,
          out_models= paths$ML_models,
          out_pred  = paths$ML_prediction
        )

      return(out)

    # ============================
    # Cross-test + LOO modeling
    # ============================
    } else if (cross_test && LOO) {

      # LOO requires special directory structure resolution
      test_path <- file.path(path, stringr::str_remove(basename(paths$matrix_path), "^LOO_"))
      test_path <- normalizePath(test_path)

      loo_pairs <- parsed |>
        dplyr::transmute(
          ref_file,
          test_file = file.path(test_path, stringr::str_remove(basename(ref_file), "_leaveout")),
          output_prefix = paste0(
            prefix_key, "_",
            drug_or_class, "_leaveout_tested_on_", strat_value, "_",
            feature
          )
        ) |>
        dplyr::filter(file.exists(test_file))

      out <- loo_pairs |>
        dplyr::mutate(
          matrix_path = paths$matrix_path,
          out_perf  = paths$ML_performance,
          out_top   = paths$ML_top_features,
          out_models= paths$ML_models,
          out_pred  = paths$ML_prediction
        )

      return(out)
    }
  }

  # If we ever get here, something wasn't covered
  stop("Unhandled combination of arguments: ",
       "MDR=", MDR, ", cross_test=", cross_test, ", LOO=", LOO,
       ", stratify_by=", if (is.null(stratify_by)) "NULL" else stratify_by)
}


#  createAllMLInputList <- function(path) {
#
#     createMLinputList(path, stratify_by = NULL, LOO = FALSE, cross_test = FALSE, MDR = FALSE) |>
#     readr::write_tsv(file.path(path, "ML.tsv"))
#
#     createMLinputList(path, stratify_by = NULL, LOO = FALSE, cross_test = TRUE, MDR = FALSE)|>
#     readr::write_tsv(file.path(path, "cross_test_ML.tsv"))
#
#     createMLinputList(path, stratify_by = NULL, LOO = FALSE, cross_test = FALSE, MDR = TRUE)|>
#     readr::write_tsv(file.path(path, "MDR_ML.tsv"))
#
#     createMLinputList(path, stratify_by = "year", LOO = FALSE, cross_test = FALSE, MDR = FALSE)|>
#     readr::write_tsv(file.path(path, "ML_year.tsv"))
#
#     createMLinputList(path, stratify_by = "year", LOO = FALSE, cross_test = TRUE, MDR = FALSE)|>
#     readr::write_tsv(file.path(path, "cross_test_ML_year.tsv"))
#
#     createMLinputList(path, stratify_by = "year", LOO = TRUE, cross_test = FALSE, MDR = FALSE)|>
#     readr::write_tsv(file.path(path, "LOO_ML_year.tsv"))
#
#     createMLinputList(path, stratify_by = "year", LOO = TRUE, cross_test = TRUE, MDR = FALSE)|>
#     readr::write_tsv(file.path(path, "LOO_cross_test_ML_year.tsv"))
#
#     createMLinputList(path, stratify_by = "country", LOO = FALSE, cross_test = FALSE, MDR = FALSE)|>
#     readr::write_tsv(file.path(path, "ML_country.tsv"))
#
#     createMLinputList(path, stratify_by = "country", LOO = FALSE, cross_test = TRUE, MDR = FALSE)|>
#     readr::write_tsv(file.path(path, "cross_test_ML_country.tsv"))
#
#     createMLinputList(path, stratify_by = "country", LOO = TRUE, cross_test = FALSE, MDR = FALSE)|>
#     readr::write_tsv(file.path(path, "LOO_ML_country.tsv"))
#
#     createMLinputList(path, stratify_by = "country", LOO = TRUE, cross_test = TRUE, MDR = FALSE)|>
#     readr::write_tsv(file.path(path, "LOO_cross_test_ML_country.tsv"))
#
# }

#' Run MDR (multi-drug resistance) machine learning models
#'
#' Executes machine learning pipeline for MDR analysis using logistic regression
#' with parallel processing via the future backend. Trains models on all MDR
#' parquet files and saves results to designated output directories.
#'
#' @param path Character scalar. Base directory containing MDR matrix files.
#' @param threads Integer. Number of parallel workers for model training. Default is 16.
#' @param split Numeric vector of length 2. Train/validation split proportions.
#' @param n_fold Integer. Number of cross-validation folds. Default 5.
#' @param prop_vi_top_feats Numeric vector of length 2. Proportion range for variable-importance selection.
#' @param pca_threshold Numeric. PCA variance threshold. Default 0.99.
#' @param verbose Logical. Print progress messages during model training. Default TRUE.
#' @param return_tune_res Logical. Return tuning results from cross-validation. Default TRUE.
#' @param return_fit Logical. Return fitted model objects. Default TRUE.
#' @param return_pred Logical. Return prediction results. Default TRUE.
#' @param use_saved_split Logical. Whether to inherit split/seed/n_fold from ml_parameters.json. Default TRUE.
#' @param shuffle_labels Logical. Randomly shuffle labels for baseline runs. Default FALSE.
#' @param use_pca Logical. Use PCA on predictors. Default FALSE.
#'
#' @return NULL (invisible). Called for side effects (model training and result saving).
#'
#' @examples
#' \dontrun{
#' # Run MDR models with default settings
#' runMDRmodels("/path/to/results")
#'
#' # Run with more threads and minimal output
#' runMDRmodels("/path/to/results",
#'              threads = 32,
#'              verbose = FALSE)
#'
#' # Run without saving model fits (save disk space)
#' runMDRmodels("/path/to/results",
#'              threads = 16,
#'              return_fit = FALSE)
#' }
#'
#' @seealso
#' \code{\link{createMLinputList}} for generating input file lists,
#' \code{\link{runMLmodels}} for non-MDR model execution
#
#' @export
runMDRmodels <- function(path,
                         threads = 16,
                         split = c(0.8, 0),
                         n_fold = 5,
                         prop_vi_top_feats = c(0, 1),
                         pca_threshold = 0.99,
                         verbose = TRUE,
                         return_tune_res = TRUE,
                         return_fit = TRUE,
                         return_pred = TRUE,
                         use_saved_split = TRUE,
                         shuffle_labels = FALSE,
                         use_pca = FALSE) {

  files <- createMLinputList(path,
                             stratify_by = NULL,
                             LOO         = FALSE,
                             cross_test  = FALSE,
                             MDR         = TRUE)

  if (nrow(files) == 0) {
    message("No MDR files found to process. Exiting.")
    return(invisible(NULL))
  }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = threads)

  if (isTRUE(verbose)) {
    nw <- tryCatch(future::nbrOfWorkers(), error = function(e) NA_integer_)
    message("runMDRmodels(): enabling multisession with workers = ", nw)
  }

  # Auto tags for shuffled and PCA
  shuffle_tag <- if (isTRUE(shuffle_labels)) "shuffled_" else ""
  pca_tag     <- if (isTRUE(use_pca)) paste0("_pca", as.character(pca_threshold)) else ""

  results_list <- future.apply::future_lapply(
    seq_len(nrow(files)),
    FUN = function(i) {

      ref_parquet   <- files$ref_file[i]
      output_prefix <- files$output_prefix[i]

      if (interactive()) {
        message(sprintf("[runMDRmodels] %d/%d: %s",
                        i, nrow(files), basename(ref_parquet)))
      }

      ml_input <- loadMLInputTibble(ref_parquet)

      sp <- if (isTRUE(use_saved_split)) {
        .resolveSplitParams(
          parquet_path = ref_parquet,
          defaults = list(split = split, seed = 5280, n_fold = n_fold)
        )
      } else {
        list(split = split, seed = 5280, n_fold = n_fold)
      }

      res <- try({
        runMLPipeline(
          ml_input_tibble = ml_input,
          test_data       = NA,
          model           = "LR",
          split           = sp$split,
          n_fold          = sp$n_fold,
          prop_vi_top_feats = prop_vi_top_feats,
          n_top_feats       = NA,
          use_pca           = use_pca,
          pca_threshold     = pca_threshold,
          shuffle_labels    = shuffle_labels,
          penalty_vec       = 10^seq(-4, -1, length.out = 10),
          mix_vec           = 0:5 / 5,
          select_best_metric = "mcc",
          seed              = sp$seed,
          verbose           = verbose,
          return_tune_res   = return_tune_res,
          return_fit        = return_fit,
          return_pred       = return_pred
        )
      }, silent = TRUE)

      if (inherits(res, "try-error")) {
        warning("Model failed for: ", output_prefix,
                "\n  Error: ", attr(res, "condition")$message)
        return(NULL)
      }

      # Final base filename: shuffled_ + <matrix prefix> + _pcaXX
      base <- paste0(shuffle_tag, output_prefix, pca_tag)

      if (!is.null(res$performance_tibble)) {
        readr::write_tsv(res$performance_tibble,
                         file.path(files$out_perf[i], paste0(base, "_performance.tsv")))
      }
      if (!is.null(res$top_feat_tibble)) {
        readr::write_tsv(res$top_feat_tibble,
                         file.path(files$out_top[i], paste0(base, "_top_features.tsv")))
      }
      if (!is.null(res$fit)) {
        saveRDS(res$fit, file.path(files$out_models[i], paste0(base, "_model_fit.rds")))
      }
      if (!is.null(res$pred)) {
        readr::write_tsv(res$pred,
                         file.path(files$out_pred[i], paste0(base, "_prediction.tsv")))
      }

      NULL
    },
    future.seed = TRUE
  )

  if (verbose) {
    message("\nML modeling complete.")
    message("Results saved under:")
    message("  ", normalizePath(path))
  }

  invisible(NULL)
}

#' Run machine learning models with multiple configurations
#'
#' Executes machine learning pipeline with support for stratification,
#' leave-one-out (LOO), and cross-testing configurations using logistic regression
#' with parallel processing. Provides flexible model training across different
#' experimental designs. Split/seed/n_fold are resolved from `ml_parameters.json`
#' when available via `.resolveSplitParams()`.
#'
#' @param path Character scalar. Base directory containing matrix files.
#' @param stratify_by Character scalar or NULL. One of "year", "country", or NULL (no stratification).
#' @param LOO Logical. Perform Leave-One-Out analysis. Default is FALSE.
#' @param cross_test Logical. Perform cross-testing between groups. Default is FALSE.
#' @param threads Integer. Number of parallel workers for model training. Default is 16.
#' @param split Numeric vector of length 2. Train/validation split proportions.
#' @param n_fold Integer. Number of cross-validation folds. Default is 5.
#' @param prop_vi_top_feats Numeric vector of length 2. Proportion range for variable-importance selection.
#' @param pca_threshold Numeric. PCA variance threshold. Default 0.99.
#' @param verbose Logical. Print progress messages during model training. Default TRUE.
#' @param return_tune_res Logical. Return tuning results from cross-validation. Default TRUE.
#' @param return_fit Logical. Return fitted model objects. Default TRUE.
#' @param return_pred Logical. Return prediction results. Default TRUE.
#' @param use_saved_split Logical. Whether to inherit split/seed/n_fold from ml_parameters.json. Default TRUE.
#' @param shuffle_labels Logical. Randomly shuffle labels for baseline runs. Default FALSE.
#' @param use_pca Logical. Use PCA on predictors. Default FALSE.
#'
#' @return NULL (invisible). Called for side effects (model training and result saving).
#'
#' @details
#' This function supports multiple analysis configurations:
#'
#' \strong{Standard mode} (\code{stratify_by = NULL, LOO = FALSE, cross_test = FALSE}):
#' \itemize{
#'   \item Trains models using train/test split from the same dataset
#'   \item Saves results to \code{ML_*} directories
#' }
#'
#' \strong{Cross-test without stratification} (\code{stratify_by = NULL, cross_test = TRUE}):
#' \itemize{
#'   \item Trains on one drug/class, tests on another drug/class
#'   \item Pairs different drugs within same feature type
#'   \item Saves results to \code{cross_test_ML_*} directories
#' }
#'
#' \strong{Cross-test with stratification} (\code{stratify_by != NULL, cross_test = TRUE}):
#' \itemize{
#'   \item Trains on one stratum (year/country), tests on another stratum
#'   \item Same drug/class across different stratification groups
#'   \item Saves results to \code{cross_test_ML_year_*} or \code{cross_test_ML_country_*} directories
#' }
#'
#' \strong{LOO with cross-test} (\code{LOO = TRUE, cross_test = TRUE}):
#' \itemize{
#'   \item Trains on leave-out dataset (one stratum excluded)
#'   \item Tests on the full dataset including the left-out stratum
#'   \item Saves results to \code{LOO_cross_test_ML_year_*} or \code{LOO_cross_test_ML_country_*} directories
#' }
#'
#' \strong{Model configuration}:
#' \itemize{
#'   \item Algorithm: Logistic Regression with elastic net regularization
#'   \item Penalty values: \code{10^seq(-4, -1, length.out = 10)}
#'   \item Mixture (alpha): 0, 0.2, 0.4, 0.6, 0.8, 1.0 (ridge to lasso)
#'   \item Selection metric: Matthews Correlation Coefficient (MCC)
#'   \item Random seed: 5280 (for reproducibility)
#'   \item PCA: Disabled
#' }
#'
#' \strong{Output file naming}:
#' Files are saved with prefixes and suffixes indicating the configuration:
#' \itemize{
#'   \item LOO: Prefixed with \code{"LOO_"}
#'   \item Cross-test: Prefixed with \code{"cross_test_"}
#'   \item Stratification: Suffixed with \code{"_country"} or \code{"_year"}
#' }
#'
#' For example: \code{"LOO_cross_test_ML_year_performance.tsv"}
#'
#' @note
#' This function requires the following packages:
#' \itemize{
#'   \item \pkg{future} - for parallel processing backend
#'   \item \pkg{future.apply} - for parallel lapply
#'   \item \pkg{readr} - for reading/writing TSV files
#'   \item \pkg{dplyr}, \pkg{purrr}, \pkg{stringr}, \pkg{tibble} - for data manipulation
#' }
#'
#' Ensure that \code{loadMLInputTibble()}, \code{runMLPipeline()}, and
#' \code{createMLinputList()} are available in your environment before calling this function.
#'
#' @examples
#' \dontrun{
#' # Standard ML models (no stratification)
#' runMLmodels("/path/to/results")
#'
#' # Cross-test between drugs (no stratification)
#' runMLmodels("/path/to/results", cross_test = TRUE)
#'
#' # Stratified by year
#' runMLmodels("/path/to/results", stratify_by = "year")
#'
#' # Cross-test with year stratification
#' runMLmodels("/path/to/results",
#'             stratify_by = "year",
#'             cross_test = TRUE,
#'             threads = 32)
#'
#' # LOO analysis stratified by country with cross-testing
#' runMLmodels("/path/to/results",
#'             stratify_by = "country",
#'             LOO = TRUE,
#'             cross_test = TRUE,
#'             verbose = TRUE)
#'
#' # Run without saving model fits (save disk space)
#' runMLmodels("/path/to/results",
#'             stratify_by = "year",
#'             return_fit = FALSE)
#' }
#'
#' @seealso
#' \code{\link{createMLinputList}} for generating input file lists,
#' \code{\link{runMDRmodels}} for MDR-specific model execution,
#' \code{\link{createMLResultDir}} for directory structure creation
#'
#' @return NULL (invisible). Called for side effects (writes results).
#' @export
runMLmodels <- function(path,
                        stratify_by = NULL,
                        LOO = FALSE,
                        cross_test = FALSE,
                        threads = 16,
                        split = c(0.8, 0),
                        n_fold = 5,
                        prop_vi_top_feats = c(0, 1),
                        pca_threshold = 0.99,
                        verbose = TRUE,
                        return_tune_res = TRUE,
                        return_fit = TRUE,
                        return_pred = TRUE,
                        use_saved_split = TRUE,
                        shuffle_labels = FALSE,
                        use_pca = FALSE) {

  if (!is.null(stratify_by)) {
    if (!is.character(stratify_by) || length(stratify_by) != 1L)
      stop("`stratify_by` must be NULL or a single string: 'year' or 'country'.")
    if (!stratify_by %in% c("year", "country"))
      stop("`stratify_by` must be NULL, 'year', or 'country'.")
  }

  files <- createMLinputList(path,
                             stratify_by = stratify_by,
                             LOO         = LOO,
                             MDR         = FALSE,
                             cross_test  = cross_test)

  if (nrow(files) == 0) {
    message("No files found to process. Exiting.")
    return(invisible(NULL))
  }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = threads)

  if (isTRUE(verbose)) {
    nw <- tryCatch(future::nbrOfWorkers(), error = function(e) NA_integer_)
    message("runMLmodels(): enabling multisession with workers = ", nw)
  }

  # LOO / cross-test configuration prefix
  config_prefix <- paste0(
    paste0(c(
      if (isTRUE(LOO)) "LOO",
      if (isTRUE(cross_test)) "cross_test"
    ), collapse = "_"),
    if (isTRUE(LOO) || isTRUE(cross_test)) "_" else ""
  )

  # Stratification suffix
  strat_suffix <- if (is.null(stratify_by) || identical(stratify_by, "")) {
    ""
  } else {
    switch(
      stratify_by,
      "country" = "_country",
      "year"    = "_year",
      stop("`stratify_by` must be NULL, 'year', or 'country'.")
    )
  }

  # Auto naming for shuffled and PCA
  shuffle_tag <- if (isTRUE(shuffle_labels)) "shuffled_" else ""
  pca_tag     <- if (isTRUE(use_pca)) paste0("_pca", as.character(pca_threshold)) else ""

  results_list <- future.apply::future_lapply(
    seq_len(nrow(files)),
    FUN = function(i) {

      ref_parquet   <- files$ref_file[i]
      output_prefix <- files$output_prefix[i]

      if (interactive()) {
        message(sprintf("[runMLmodels] %d/%d: %s",
                        i, nrow(files), basename(ref_parquet)))
      }

      ml_input <- loadMLInputTibble(ref_parquet)

      if ("test_file" %in% names(files)) {
        val <- files$test_file[i]
        test_data <- if (!is.na(val) && !is.null(val)) {
          loadMLInputTibble(as.character(val))
        } else {
          NULL
        }
      } else {
        test_data <- NULL
      }

      sp <- if (isTRUE(use_saved_split)) {
        .resolveSplitParams(
          parquet_path = ref_parquet,
          defaults = list(split = split, seed = 5280, n_fold = n_fold)
        )
      } else {
        list(split = split, seed = 5280, n_fold = n_fold)
      }

      res <- try({
        runMLPipeline(
          ml_input_tibble = ml_input,
          test_data       = test_data,
          model           = "LR",
          split           = sp$split,
          n_fold          = sp$n_fold,
          prop_vi_top_feats = prop_vi_top_feats,
          n_top_feats       = NA,
          use_pca           = use_pca,
          pca_threshold     = pca_threshold,
          shuffle_labels    = shuffle_labels,
          penalty_vec       = 10^seq(-4, -1, length.out = 10),
          mix_vec           = 0:5 / 5,
          select_best_metric = "mcc",
          seed              = sp$seed,
          verbose           = verbose,
          return_tune_res   = return_tune_res,
          return_fit        = return_fit,
          return_pred       = return_pred
        )
      }, silent = TRUE)

      if (inherits(res, "try-error")) {
        warning("Model failed for: ", output_prefix,
                "\n  Error: ", attr(res, "condition")$message)
        return(NULL)
      }

      # Final base filename: shuffled_ + [LOO_/cross_test_] + <matrix prefix> + _pcaXX + _year/_country
      base <- paste0(shuffle_tag, config_prefix, output_prefix, pca_tag, strat_suffix)

      if (!is.null(res$performance_tibble)) {
        readr::write_tsv(res$performance_tibble,
                         file.path(files$out_perf[i], paste0(base, "_performance.tsv")))
      }
      if (!is.null(res$top_feat_tibble)) {
        readr::write_tsv(res$top_feat_tibble,
                         file.path(files$out_top[i], paste0(base, "_top_features.tsv")))
      }
      if (!is.null(res$fit)) {
        saveRDS(res$fit, file.path(files$out_models[i], paste0(base, "_model_fit.rds")))
      }
      if (!is.null(res$pred)) {
        readr::write_tsv(res$pred,
                         file.path(files$out_pred[i], paste0(base, "_prediction.tsv")))
      }

      NULL
    },
    future.seed = TRUE
  )

  if (verbose) {
    message("\nML modeling complete.")
    message("Results saved under:")
    message("  ", normalizePath(path))
  }

  invisible(NULL)
}



#' Run the entire AMR ML pipeline from a parquet-backed DuckDB
#'
#' This function provides a complete end-to-end AMR machine learning workflow.
#' Given a DuckDB file produced by `runDataProcessing()`, it:
#'   1. Generates all ML feature matrices (drug, class, year, country, MDR, LOO)
#'   2. Creates all ML directory structures
#'   3. Prepares ML input lists for every mode
#'   4. Runs logistic regression ML models (standard + stratified + cross-test + MDR)
#'   5. Saves performance metrics, fitted models, predictions, and top feature rankings
#'
#' @param parquet_duckdb_path Path to a `<bug>_parquet.duckdb` produced by data_processing.R
#' @param threads Number of parallel workers. Default: 16
#' @param n_fold Cross-validation folds (default: 5). Use 0 or NULL for classical splits.
#' @param split Training/validation split (default: c(1,0) for CV mode)
#' @param min_n Minimum samples per drug class for MDR matrices (default: 25)
#' @param prop_vi_top_feats Proportion of variable importance for top features (default: c(0,1))
#' @param pca_threshold PCA variance threshold (not used unless `use_pca = TRUE`)
#' @param verbose Print progress updates? Default: TRUE
#' @param use_saved_split Whether to inherit split/seed/n_fold from ml_parameters.json
#'
#' @return Invisibly returns the output directory used for ML results.
#'
#' @export
runModelingPipeline <- function(parquet_duckdb_path,
                                threads = 16,
                                n_fold = 5,
                                split = c(1, 0),
                                min_n = 25,
                                prop_vi_top_feats = c(0, 1),
                                pca_threshold = 0.99,
                                verbose = TRUE,
                                use_saved_split = TRUE) {

  parquet_duckdb_path <- normalizePath(parquet_duckdb_path)
  if (!file.exists(parquet_duckdb_path)) {
    stop("Parquet-backed DuckDB at ", parquet_duckdb_path, " not found.\n",
         "Are you using `{Bug}.duckdb` instead of `{Bug}_parquet.duckdb?`")
  }

  out_root <- dirname(parquet_duckdb_path)

  if (verbose) {
    message("\n=== amR_ml: Full pipeline runner ===")
    message("Using parquet-backed DB:\n  ", parquet_duckdb_path)
  }

  if (verbose) message("\n[1/4] Generating ML feature matrices.")
  generateMLInputs(
    parquet_duckdb_path = parquet_duckdb_path,
    out_path = out_root,
    n_fold   = n_fold,
    split    = split,
    min_n    = min_n,
    verbosity = if (verbose) "minimal" else "debug"
  )

  if (verbose) message("\n[2/4] Running standard ML models.")
  runMLmodels(
    path = out_root,
    stratify_by = NULL,
    LOO = FALSE,
    cross_test = FALSE,
    threads = threads,
    split = split,
    n_fold = n_fold,
    prop_vi_top_feats = prop_vi_top_feats,
    pca_threshold = pca_threshold,
    verbose = verbose,
    use_saved_split = use_saved_split
  )

  if (verbose) message("\n[3/4] Running stratified (year) ML models.")
  runMLmodels(
    path = out_root,
    stratify_by = "year",
    LOO = FALSE,
    cross_test = FALSE,
    threads = threads,
    split = split,
    n_fold = n_fold,
    prop_vi_top_feats = prop_vi_top_feats,
    pca_threshold = pca_threshold,
    verbose = verbose,
    use_saved_split = use_saved_split
  )

  if (verbose) message("\n[3/4] Running stratified (country) ML models.")
  runMLmodels(
    path = out_root,
    stratify_by = "country",
    LOO = FALSE,
    cross_test = FALSE,
    threads = threads,
    split = split,
    n_fold = n_fold,
    prop_vi_top_feats = prop_vi_top_feats,
    pca_threshold = pca_threshold,
    verbose = verbose,
    use_saved_split = use_saved_split
  )

  if (verbose) message("\n[4/4] Running MDR ML models.")
  runMDRmodels(
    path = out_root,
    threads = threads,
    split = split,
    n_fold = n_fold,
    prop_vi_top_feats = prop_vi_top_feats,
    pca_threshold = pca_threshold,
    verbose = verbose,
    use_saved_split = use_saved_split
  )
  # All done!
  if (verbose) {
    message("\n=== AMR-ML Pipeline Complete ===")
    message("All matrices, models, top feature lists, and performance metrics saved under:\n  ",
            out_root)
    message("\nTo inspect model outputs, see directories such as:")
    message("  ML_performance/, ML_models/, ML_prediction/, ML_top_features/")
  }

  invisible(out_root)
}

