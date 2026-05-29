#' @importFrom DBI dbConnect dbDisconnect dbExecute dbGetQuery dbWriteTable
#' @importFrom duckdb duckdb
#' @importFrom glue glue
#' @importFrom readr write_lines
#' @importFrom stringr str_remove str_split
#' @importFrom tidyr drop_na
#' @importFrom dplyr bind_rows count desc distinct filter group_by mutate
#' @importFrom dplyr n_distinct pull select ungroup arrange
#' @importFrom purrr map map_int pmap_chr walk
#' @importFrom arrow read_parquet write_parquet
#' @importFrom tibble tibble
NULL

# A logging verbosity function to add debugging info
.make_logger <- function(verbosity = c("minimal", "debug")) {
  v <- match.arg(verbosity)
  function(level = c("info", "debug"), ...) {
    lvl <- match.arg(level)
    # Print basic info messages in either mode, print full info only in debug mode
    if (lvl == "info" || (lvl == "debug" && v == "debug")) {
      message(...)
    }
  }
}

#' .calculateMinSamples()
#'
#' Returns the minimum number of total observations needed (one bug-drug combo)
#' to have at least one observation of the rarer class in each CV fold
#' or, in validation mode, across train/val/test.
#'
#' @param n_fold [numeric] Number of CV folds; must be >= 2 in CV mode
#' @param split [numeric] Vector length 2: c(train_prop, val_prop).
#'        For pure CV we enforce c(1, 0). For classical splits: val > 0, test > 0.
#' @param res_prop [numeric] Proportion of resistant genomes in the combo
#' @return Minimum total observations required so each fold (CV) or each partition (val mode) contains the rarer class at least once
#' @keywords internal
.calculateMinSamples <- function(n_fold, split, res_prop) {
  .checkArgSplit(split) # validates split per new policy
  .checkArgResProp(res_prop) # validates [0,1]

  # If one class is absent entirely, we cannot model it, so we skip using Inf
  if (res_prop == 0 || res_prop == 1) {
    return(Inf)
  }

  # Rarer class proportion (min of resistant vs. susceptible)
  lowest_prop_rs <- if (res_prop <= 0.5) res_prop else (1 - res_prop)

  train_prop <- split[1]
  val_prop <- split[2]
  test_prop <- 1 - train_prop - val_prop

  if (train_prop == 1 && val_prop == 0) {
    # Cross-validation mode (CV)
    .checkArgNFold(n_fold)
    # Req.: rarer_class_count / n_fold >= 1  => total_obs >= n_fold / lowest_prop_rs
    min_samples <- n_fold / lowest_prop_rs
  } else if (val_prop > 0 && test_prop > 0) {
    # Train-validate-test mode (splits)
    # All three partitions must be positive
    lowest_prop_tvt <- min(train_prop, val_prop, test_prop)
    min_samples <- 1 / (lowest_prop_rs * lowest_prop_tvt)
  } else {
    # Unsupported combinations (like a CV + separate test holdout)
    stop("Unsupported split configuration in calculateMinSamples().")
  }

  return(ceiling(min_samples))
}


# Drop where class balance prevents CV from adequately showing all classes
#' skipImbalancedMatrix()
#'
#' Evaluates whether a drug-bug combination should be skipped based on
#' insufficient samples or class imbalance for cross-validation or train/val/test splits.
#'
#' @param genome_ids The subsetted genome IDs for a given matrix combination
#' @param phenotype_counts The count of AMR observations for that combination
#' @param n_fold Number of CV folds
#' @param split For CV: c(1, 0); for classical splits: c(train_prop, val_prop)
#' @param min_total_obs Minimum total observations required before attempting CV (default 40)
#' @param verbosity "minimal" or "debug"; "debug" offers more info per processed matrix
#' @return Logical; TRUE if matrix should be skipped, FALSE if matrix should be generated
#' @keywords internal
#' @examples
#' \dontrun{
#' # Check if a drug-bug combination should be skipped
#' genome_ids <- c("genome1", "genome2", "genome3")
#' phenotype_counts <- data.frame(
#'   phenotype = c("Resistant", "Susceptible"),
#'   count = c(2, 1)
#' )
#'
#' # For 5-fold CV
#' skip <- skipImbalancedMatrix(
#'   genome_ids = genome_ids,
#'   phenotype_counts = phenotype_counts,
#'   n_fold = 5,
#'   split = c(1, 0),
#'   min_total_obs = 40,
#'   verbosity = "minimal"
#' )
#' }
skipImbalancedMatrix <- function(genome_ids,
                                 phenotype_counts,
                                 n_fold = 5,
                                 split = c(1, 0),
                                 min_total_obs = 40,
                                 verbosity = c("minimal", "debug")) {
  log <- .make_logger(verbosity)

  total_obs <- sum(phenotype_counts$count)
  res_count <- phenotype_counts$count[phenotype_counts$phenotype == "Resistant"]
  res_prop <- if (length(res_count) == 0) 0 else res_count / total_obs

  min_samples <- .calculateMinSamples(n_fold = n_fold, split = split, res_prop = res_prop)

  # Debug logging block (printed if verbosity == "debug")
  log("debug", "DEBUG: Checking matrix skip conditions")
  log("debug", "  total_obs: ", total_obs)
  log("debug", "  genome_ids length: ", length(genome_ids))
  log("debug", "  resistant count: ", ifelse(length(res_count) == 0, 0, res_count))
  log("debug", "  resistant proportion: ", round(res_prop, 4))
  log("debug", "  split: ", paste(split, collapse = ", "))
  log("debug", "  n_fold: ", n_fold)
  log("debug", "  min_samples required: ", min_samples)
  log("debug", "  min_total_obs threshold: ", min_total_obs)

  skip <- (total_obs < min_total_obs) || (length(genome_ids) < min_samples)

  if (skip) {
    log("debug", "  --> SKIPPED! This matrix has too few samples for CV/split.")
  } else {
    log("debug", "  --> KEPT! This matrix has sufficient sample counts and class balance.")
  }

  return(skip)
}


#' Metadata and Feature parquets to bug-drug-feature parquet as ML input
#'
#' @param parquet_duckdb_path [character] The path to the DuckDB that contains the view of metadata and feature parquets
#' @param path [character] The path to the working directory
#' @param n_fold [numeric] the number of cross-validation folds
#' @param split if folds has been set to NULL, indicating classical splits instead
#' @param stratify_by [character] default is NULL; other options are "year" or "country"
#' @param verbosity [character] "minimal" or "debug"
#' @return invisible tibble of created files
#' @keywords internal
.parquet2Matrix <- function(parquet_duckdb_path,
                            path,
                            n_fold,
                            split,
                            stratify_by = NULL,
                            verbosity = c("minimal", "debug")) {
  verbosity <- match.arg(verbosity)
  log <- .make_logger(verbosity)

  parquet_duckdb_path <- normalizePath(parquet_duckdb_path)
  path <- normalizePath(path)

  # Choose stratification column
  if (is.null(stratify_by) || !(stratify_by %in% c("year", "country"))) {
    stratify_column <- NULL
  } else if (stratify_by == "year") {
    stratify_column <- "year_bin"
  } else {
    stratify_column <- "country_abbr"
  }

  # Determine matrix output directory
  if (!is.null(stratify_column) && stratify_by == "country") {
    matrix_path <- file.path(path, "matrix_country")
    grouping_variables <- list(
      drug_country       = c("drug_abbr", "country_abbr"),
      drug_class_country = c("class_abbr", "country_abbr")
    )
  } else if (!is.null(stratify_column) && stratify_by == "year") {
    matrix_path <- file.path(path, "matrix_year")
    grouping_variables <- list(
      drug_year       = c("drug_abbr", "year_bin"),
      drug_class_year = c("class_abbr", "year_bin")
    )
  } else {
    matrix_path <- file.path(path, "matrix")
    grouping_variables <- list(
      drug       = "drug_abbr",
      drug_class = "class_abbr"
    )
    stratify_column <- NULL
  }

  if (!dir.exists(matrix_path)) dir.create(matrix_path, recursive = TRUE)

  log("info", paste0("Matrix output directory: ", matrix_path))
  log("debug", paste0(
    "Stratification: ",
    ifelse(is.null(stratify_column), "None", stratify_column)
  ))

  # Feature and matrix types
  feature_types <- list(
    genes    = list(view = "gene_count", id_col = "gene"),
    proteins = list(view = "protein_count", id_col = "protein"),
    domains  = list(view = "domain_count", id_col = "domain"),
    struct   = list(view = "struct", id_col = "struct")
  )

  matrix_types <- list(
    counts        = list(value_col = "value", filter = "variance > 0", binary_only = FALSE, label = "counts"),
    binary        = list(value_col = "present", filter = "variance > 0", binary_only = FALSE, label = "binary"),
    struct_binary = list(value_col = "present", filter = "variance > 0", binary_only = TRUE, label = "binary")
  )

  # Initialize skipped log
  log_path <- file.path(matrix_path, "skipped_drugs.log")
  readr::write_lines("Drug\tReason\tTotal_obs\tPhenotype_distribution", log_path)

  # Safe DBI-quoting
  quote_condition <- function(group_cols, group_values, con) {
    ids <- vapply(
      group_cols,
      function(col) DBI::dbQuoteIdentifier(con, col),
      character(1)
    )
    vals <- vapply(
      group_cols,
      function(col) {
        val <- group_values[[col]][1]
        DBI::dbQuoteString(con, as.character(val))
      },
      character(1)
    )

    paste(sprintf("%s = %s", ids, vals), collapse = " AND ")
  }

  con0 <- DBI::dbConnect(duckdb::duckdb(), parquet_duckdb_path)
  on.exit(DBI::dbDisconnect(con0, shutdown = FALSE), add = TRUE)
  DBI::dbExecute(con0, sprintf(
    "SET file_search_path='%s'",
    dirname(parquet_duckdb_path)
  ))

  bug <- stringr::str_remove(basename(parquet_duckdb_path), "_parquet\\.duckdb$")
  log("info", paste0("Connected to DuckDB for bug: ", bug))

  # Iterate over grouping variables
  all_group_definitions <- lapply(grouping_variables, function(cols) {
    idcols <- paste0('"', cols, '"', collapse = ", ")
    DBI::dbGetQuery(con0, sprintf("SELECT DISTINCT %s FROM metadata", idcols)) |>
      tidyr::drop_na()
  })

  # Each group gets processed with a fresh DB connection
  for (group_type in names(grouping_variables)) {
    group_cols <- grouping_variables[[group_type]]
    all_groups <- all_group_definitions[[group_type]]

    log("debug", paste0("Found ", nrow(all_groups), " groups for type: ", group_type))

    for (i in seq_len(nrow(all_groups))) {
      # New connection for this group
      con <- DBI::dbConnect(duckdb::duckdb(), parquet_duckdb_path)
      DBI::dbExecute(con, sprintf(
        "SET file_search_path='%s'",
        dirname(parquet_duckdb_path)
      ))

      group_values <- all_groups[i, , drop = FALSE]
      group_label <- paste(group_values, collapse = "_")

      log("info", paste0("Building ML matrices for ", group_type, ": ", group_label))

      condition_string <- quote_condition(group_cols, group_values, con)

      # Strat filter
      strat_filter <- if (!is.null(stratify_column)) {
        sprintf("AND \"%s\" IS NOT NULL AND \"%s\" != ''", stratify_column, stratify_column)
      } else {
        ""
      }

      # Genome selection logic
      if (group_type %in% c("drug_class", "drug_class_year", "drug_class_country")) {
        genome_ids <- DBI::dbGetQuery(con, sprintf("
          WITH class_phenotypes AS (
            SELECT \"genome.genome_id\" AS genome_id,
                   MAX(CASE WHEN \"genome_drug.resistant_phenotype\" = 'Resistant'
                            THEN 1 ELSE 0 END) AS any_resistant,
                   MIN(CASE WHEN \"genome_drug.resistant_phenotype\" = 'Susceptible'
                            THEN 1 ELSE 0 END) AS all_susceptible
            FROM metadata
            WHERE %s
            GROUP BY \"genome.genome_id\"
          )
          SELECT genome_id
          FROM class_phenotypes
          WHERE any_resistant = 1 OR all_susceptible = 1
        ", condition_string))[[1]]
      } else {
        genome_ids <- DBI::dbGetQuery(con, sprintf("
          SELECT DISTINCT \"genome.genome_id\"
          FROM metadata
          WHERE %s
            AND \"genome_drug.resistant_phenotype\" IN ('Resistant','Susceptible')
            %s
        ", condition_string, strat_filter))[[1]]
      }

      # Phenotype counts
      phenotype_counts_all <- DBI::dbGetQuery(con, sprintf("
        SELECT \"genome_drug.resistant_phenotype\" AS phenotype, COUNT(*) AS count
        FROM metadata
        WHERE %s
        GROUP BY \"genome_drug.resistant_phenotype\"
      ", condition_string))

      phenotype_summary <- paste(
        apply(
          phenotype_counts_all, 1,
          function(row) paste0(row["phenotype"], "=", row["count"])
        ),
        collapse = "; "
      )

      # Apply skip logic
      if (skipImbalancedMatrix(genome_ids, phenotype_counts_all, n_fold, split,
        verbosity = verbosity
      )) {
        readr::write_lines(
          sprintf(
            "%s\tToo few samples for CV/split\t%d\t%s",
            group_label, length(genome_ids), phenotype_summary
          ),
          log_path,
          append = TRUE
        )

        DBI::dbDisconnect(con, shutdown = FALSE)
        next
      }

      if (length(genome_ids) < 40) {
        readr::write_lines(
          sprintf(
            "%s\tToo few observations\t%d\t%s",
            group_label, length(genome_ids), phenotype_summary
          ),
          log_path,
          append = TRUE
        )

        DBI::dbDisconnect(con, shutdown = FALSE)
        next
      }

      phen2 <- DBI::dbGetQuery(con, sprintf("
        SELECT \"genome_drug.resistant_phenotype\", COUNT(*) AS count
        FROM metadata
        WHERE %s
          AND \"genome_drug.resistant_phenotype\" IN ('Resistant','Susceptible')
          %s
        GROUP BY \"genome_drug.resistant_phenotype\"
      ", condition_string, strat_filter))

      if (nrow(phen2) < 2) {
        readr::write_lines(
          sprintf(
            "%s\tOnly one phenotype class\t%d\t%s",
            group_label, length(genome_ids), phenotype_summary
          ),
          log_path,
          append = TRUE
        )

        DBI::dbDisconnect(con, shutdown = FALSE)
        next
      }

      # Create selected_genomes
      DBI::dbExecute(con, "CREATE OR REPLACE TEMP TABLE selected_genomes (genome_id VARCHAR)")
      DBI::dbWriteTable(con, "selected_genomes",
        data.frame(genome_id = genome_ids),
        append = TRUE
      )

      # Feature and matrix generation steps
      for (ftype in names(feature_types)) {
        fview <- feature_types[[ftype]]$view
        fid <- feature_types[[ftype]]$id_col

        # binary view
        DBI::dbExecute(con, sprintf("
          CREATE OR REPLACE VIEW %s_binary AS
          SELECT genome_id, %s,
                 CASE WHEN value > 0 THEN 1 ELSE 0 END AS present
          FROM %s
        ", ftype, fid, fview))

        # counts view
        if (ftype != "struct") {
          DBI::dbExecute(con, sprintf("
            CREATE OR REPLACE VIEW %s_counts AS
            SELECT genome_id, %s, value
            FROM %s
          ", ftype, fid, fview))
        }

        for (mtype in names(matrix_types)) {
          binary_only <- matrix_types[[mtype]]$binary_only
          if (ftype == "struct" && !binary_only) next

          mview <- sprintf(
            "%s_%s", ftype,
            ifelse(grepl("binary", mtype), "binary", "counts")
          )
          value_col <- matrix_types[[mtype]]$value_col
          filter_clause <- matrix_types[[mtype]]$filter

          # select features with non-zero variance
          keep_query <- sprintf("
            SELECT %s AS feature_id, VAR_POP(%s) AS variance
            FROM %s
            JOIN selected_genomes USING (genome_id)
            GROUP BY %s
            HAVING %s
          ", fid, value_col, mview, fid, filter_clause)

          keep_features <- DBI::dbGetQuery(con, keep_query)[["feature_id"]]
          if (length(keep_features) == 0) {
            log("info", paste0(
              "All features filtered for ",
              ftype, " - ", mtype, " - ", group_label
            ))
            next
          }

          DBI::dbExecute(
            con,
            "CREATE OR REPLACE TEMP TABLE keep_features (feature_id VARCHAR)"
          )
          DBI::dbWriteTable(con,
            "keep_features",
            data.frame(feature_id = keep_features),
            append = TRUE
          )

          mtype_label <- matrix_types[[mtype]]$label

          long_out_path <- file.path(
            matrix_path,
            sprintf(
              "%s_%s_%s_%s_%s_sparse.parquet",
              bug, group_type, group_label, ftype, mtype_label
            )
          )

          long_out_path_sql <- gsub("\\\\", "/", long_out_path)

          # phenotype case
          phenotype_case <- if (group_type %in%
            c("drug_class", "drug_class_year", "drug_class_country")) {
            "
               CASE
                  WHEN MAX(CASE WHEN f.\"genome_drug.resistant_phenotype\"='Resistant'
                                THEN 1 ELSE 0 END)=1 THEN 'Resistant'
                  WHEN MIN(CASE WHEN f.\"genome_drug.resistant_phenotype\"='Susceptible'
                                THEN 1 ELSE 0 END)=1 THEN 'Susceptible'
               END
            "
          } else {
            "
              MAX(
                CASE f.\"genome_drug.resistant_phenotype\"
                   WHEN 'Resistant' THEN 'Resistant'
                   WHEN 'Susceptible' THEN 'Susceptible'
                END
              )
            "
          }

          strat_col_select <- if (!is.null(stratify_by)) {
            sprintf(", f.\"%s\"", stratify_column)
          } else {
            ""
          }

          strat_col_group <- if (!is.null(stratify_by)) {
            sprintf(", f.\"%s\"", stratify_column)
          } else {
            ""
          }

          copy_sql <- sprintf(
            "
            COPY (
              SELECT
                f.\"genome.genome_id\" AS genome_id,
                %s AS feature_id,
                MAX(CAST(%s AS DOUBLE)) AS value,
                %s AS \"genome_drug.resistant_phenotype\"
                %s
              FROM %s
              JOIN selected_genomes USING (genome_id)
              JOIN keep_features kf ON %s = kf.feature_id
              JOIN metadata f ON genome_id = f.\"genome.genome_id\"
              WHERE %s
                AND f.\"genome_drug.resistant_phenotype\" IN ('Resistant','Susceptible')
                %s
              GROUP BY f.\"genome.genome_id\", %s %s
              ORDER BY f.\"genome.genome_id\", %s
            )
            TO '%s'
            (FORMAT 'parquet', COMPRESSION 'zstd')
          ",
            fid, value_col, phenotype_case, strat_col_select,
            mview, fid, condition_string,
            strat_filter, fid, strat_col_group, fid,
            long_out_path_sql
          )

          ok <- try(DBI::dbExecute(con, copy_sql), silent = TRUE)

          # On copy failure, log + continue without stopping entire pipeline
          if (inherits(ok, "try-error")) {
            readr::write_lines(
              sprintf(
                "%s\tCOPY_failed\t%d\t%s",
                group_label, length(genome_ids), phenotype_summary
              ),
              log_path,
              append = TRUE
            )
            log("debug", paste0("COPY failed for ", group_label, ". Continuing."))
            next
          }

          log("info", paste0("Exported matrix: ", long_out_path))
        }
      }

      # Close group connection to avoid stack overflow behaviors
      DBI::dbDisconnect(con, shutdown = FALSE)
    }
  }
}

#' Build leave-one-out (LOO) merged parquet matrices from stratified parquet files.
#'
#' @param path Character. Base directory containing stratified parquet matrices.
#'             Expected subdirs: matrix_year/ or matrix_country/.
#' @param stratify_by Character or NULL. One of "year", "country".
#' @param verbosity Character. "minimal" or "debug"; when "debug", prints detailed steps.
#' @return Invisibly returns a tibble with paths of created LOO parquet files.
#' @import arrow dplyr purrr stringr tibble
.parquet2LOOMatrix <- function(path, stratify_by = "year", verbosity = c("minimal", "debug")) {
  verbosity <- match.arg(verbosity)
  log <- .make_logger(verbosity)

  # Validate intended stratification strategy
  if (!is.null(stratify_by) && !(stratify_by %in% c("year", "country"))) {
    stop("`stratify_by` must be 'year', 'country', or NULL.")
  }
  if (is.null(stratify_by)) {
    log("info", "No stratification provided (stratify_by = NULL). Leave-one-out is undefined. Exiting.")
    return(invisible(tibble::tibble(created_file = character())))
  }

  # Normalize paths to forward slashes for consistency
  matrix_path <- gsub("\\\\", "/", file.path(path, paste0("matrix_", stratify_by)))
  LOO_path <- gsub("\\\\", "/", file.path(path, paste0("LOO_matrix_", stratify_by)))

  if (!dir.exists(matrix_path)) {
    log("info", paste0("The matrix directory ", matrix_path, " does not exist."))
    return(invisible(tibble::tibble(created_file = character())))
  }
  if (!dir.exists(LOO_path)) dir.create(LOO_path, recursive = TRUE)

  files <- list.files(matrix_path, pattern = "\\.parquet$", full.names = TRUE)
  if (length(files) == 0) {
    log("info", paste0("No parquet files found in ", matrix_path))
    return(invisible(tibble::tibble(created_file = character())))
  }
  files <- gsub("\\\\", "/", files)

  log("debug", paste0("Found ", length(files), " parquet files."))

  # Parse filename structure
  parsed_info <- tibble::tibble(
    file  = files,
    parts = stringr::str_split(basename(files), "_")
  ) |>
    dplyr::mutate(
      idx_sparse = purrr::map_int(parts, function(x) {
        w <- which(x == "sparse.parquet")
        if (length(w) == 0) NA_integer_ else w[1]
      }),
      idx_strat = purrr::map_int(parts, function(x) {
        w <- which(x == stratify_by)
        if (length(w) == 0) NA_integer_ else w[1]
      }),
      feature = purrr::pmap_chr(list(parts, idx_sparse), function(x, i) {
        if (is.na(i) || i < 8) {
          return(NA_character_)
        }
        paste(x[(i - 2):(i - 1)], collapse = "_")
      }),
      prefix = purrr::pmap_chr(list(parts, idx_strat), function(x, i) {
        if (is.na(i) || i < 3) {
          return(NA_character_)
        }
        paste(x[(i - 3):(i - 1)], collapse = "_")
      }),
      drug_or_class = purrr::pmap_chr(list(parts, idx_strat), function(x, i) {
        if (is.na(i) || (i + 1) > length(x)) {
          return(NA_character_)
        }
        x[i + 1]
      }),
      stratification = purrr::pmap_chr(list(parts, idx_strat), function(x, i) {
        if (is.na(i) || (i + 2) > length(x)) {
          return(NA_character_)
        }
        x[i + 2]
      })
    ) |>
    dplyr::filter(!is.na(prefix), !is.na(feature), !is.na(drug_or_class), !is.na(stratification))

  if (nrow(parsed_info) == 0) {
    log("info", paste0("No files matched the expected naming pattern for stratify_by = '", stratify_by, "'."))
    return(invisible(tibble::tibble(created_file = character())))
  }

  log("debug", paste0("Parsed ", nrow(parsed_info), " files."))

  prefix_feature <- parsed_info |> dplyr::distinct(prefix, feature)
  created <- character(0)

  for (i in seq_len(nrow(prefix_feature))) {
    sub_prefix <- prefix_feature$prefix[i]
    sub_feature <- prefix_feature$feature[i]

    group_df <- parsed_info |>
      dplyr::filter(prefix == sub_prefix, feature == sub_feature) |>
      dplyr::group_by(drug_or_class) |>
      dplyr::filter(dplyr::n_distinct(stratification) >= 3) |>
      dplyr::ungroup()

    if (nrow(group_df) == 0) next

    grouped <- split(group_df, group_df$drug_or_class)

    purrr::walk(names(grouped), function(drug_class) {
      group <- grouped[[drug_class]]

      purrr::walk(unique(group$stratification), function(leave_one_out) {
        subset <- dplyr::filter(group, stratification != leave_one_out)
        if (nrow(subset) == 0) {
          return(NULL)
        }

        combined <- purrr::map(subset$file, arrow::read_parquet) |>
          dplyr::bind_rows()

        out_file <- gsub("\\\\", "/", file.path(
          LOO_path,
          paste0(
            sub_prefix, "_", stratify_by, "_",
            drug_class, "_leaveout_", leave_one_out, "_",
            sub_feature, "_sparse.parquet"
          )
        ))
        arrow::write_parquet(combined, out_file)
        created <<- c(created, out_file)

        log("debug", paste0("Created LOO file: ", out_file))
      })
    })
  }

  log("info", "All LOO matrices generated and saved.")
  invisible(tibble::tibble(created_file = created))
}


#' Generate ML input matrix for multidrug resistance analysis
#'
#' @param parquet_duckdb_path [character] Path to the DuckDB that contains the view of metadata and feature parquets
#' @param path [character] Working directory
#' @param min_n [numeric] Minimum number of samples for each combination of drug classes
#' @param verbosity [character] "minimal" or "debug"; when "debug", prints detailed steps
#' @return invisible(TRUE); side effect: writes MDR matrices to disk
#' @keywords internal
.parquet2MDRMatrix <- function(parquet_duckdb_path,
                               path,
                               min_n,
                               verbosity = c("minimal", "debug")) {
  verbosity <- match.arg(verbosity)
  log <- .make_logger(verbosity)

  parquet_duckdb_path <- normalizePath(parquet_duckdb_path)
  path <- normalizePath(path)
  MDR_matrix_path <- gsub("\\\\", "/", file.path(path, "MDR_matrix"))
  if (!dir.exists(MDR_matrix_path)) dir.create(MDR_matrix_path, recursive = TRUE)

  # Listing feature types and DuckDB views
  feature_types <- list(
    genes    = list(view = "gene_count", id_col = "gene"),
    proteins = list(view = "protein_count", id_col = "protein"),
    domains  = list(view = "domain_count", id_col = "domain"),
    struct   = list(view = "struct", id_col = "struct")
  )
  # Listing matrix types and filters
  matrix_types <- list(
    counts        = list(value_col = "value", filter = "variance > 0", binary_only = FALSE, label = "counts"),
    binary        = list(value_col = "present", filter = "variance > 0", binary_only = FALSE, label = "binary"),
    struct_binary = list(value_col = "present", filter = "variance > 0", binary_only = TRUE, label = "binary")
  )

  # Connect to DuckDB
  con0 <- DBI::dbConnect(duckdb::duckdb(), parquet_duckdb_path)
  DBI::dbExecute(con0, sprintf(
    "SET file_search_path='%s'",
    dirname(parquet_duckdb_path)
  ))
  # Discover bug name and which genomes to include
  bug <- stringr::str_remove(basename(parquet_duckdb_path), "_parquet\\.duckdb$")
  metadata_all <- DBI::dbGetQuery(con0, "SELECT * FROM metadata")
  DBI::dbDisconnect(con0, shutdown = FALSE)

  classes <- metadata_all |>
    dplyr::select(genome.genome_id, resistant_classes) |>
    dplyr::distinct() |>
    dplyr::group_by(resistant_classes) |>
    dplyr::count() |>
    dplyr::arrange(dplyr::desc(n)) |>
    dplyr::filter(n >= min_n, resistant_classes != "Intermediate") |>
    dplyr::pull(resistant_classes)

  log("info", paste0("Building MDR matrices for classes: ", paste(classes, collapse = ", ")))

  genomes_to_keep <- metadata_all |>
    dplyr::filter(resistant_classes %in% classes) |>
    dplyr::pull(genome.genome_id)

  # Build one matrix per feature type and matrix type
  for (ftype in names(feature_types)) {
    fview <- feature_types[[ftype]]$view
    fid <- feature_types[[ftype]]$id_col

    for (mtype in names(matrix_types)) {
      binary_only <- matrix_types[[mtype]]$binary_only
      if (ftype == "struct" && !binary_only) next

      mtype_label <- matrix_types[[mtype]]$label

      out_file <- file.path(
        MDR_matrix_path,
        sprintf("%s_MDR_resistant_classes_%s_%s_sparse.parquet", bug, ftype, mtype_label)
      )
      out_file_sql <- gsub("\\\\", "/", out_file)

      # Fresh connection for this matrix
      con <- DBI::dbConnect(duckdb::duckdb(), parquet_duckdb_path)
      DBI::dbExecute(con, sprintf(
        "SET file_search_path='%s'",
        dirname(parquet_duckdb_path)
      ))

      # Selected genomes
      DBI::dbExecute(con, "CREATE OR REPLACE TEMP TABLE selected_genomes (genome_id VARCHAR)")
      DBI::dbWriteTable(con, "selected_genomes",
        data.frame(genome_id = genomes_to_keep),
        append = TRUE
      )

      # Binary view
      DBI::dbExecute(con, sprintf("
        CREATE OR REPLACE VIEW %s_binary AS
        SELECT genome_id, %s,
               CASE WHEN value > 0 THEN 1 ELSE 0 END AS present
        FROM %s
      ", ftype, fid, fview))

      # Counts view
      if (ftype != "struct") {
        DBI::dbExecute(con, sprintf("
          CREATE OR REPLACE VIEW %s_counts AS
          SELECT genome_id, %s, value
          FROM %s
        ", ftype, fid, fview))
      }

      mview <- sprintf("%s_%s", ftype, ifelse(grepl("binary", mtype), "binary", "counts"))
      value_col <- matrix_types[[mtype]]$value_col
      filter_clause <- matrix_types[[mtype]]$filter

      keep_query <- sprintf("
        SELECT %s AS feature_id, VAR_POP(%s) AS variance
        FROM %s
        JOIN selected_genomes USING (genome_id)
        GROUP BY %s
        HAVING %s
      ", fid, value_col, mview, fid, filter_clause)

      keep_features <- DBI::dbGetQuery(con, keep_query)[["feature_id"]]
      if (length(keep_features) == 0) {
        log("info", paste0("All features filtered for MDR: ", ftype, " - ", mtype))
        DBI::dbDisconnect(con, shutdown = FALSE)
        next
      }

      DBI::dbExecute(con, "CREATE OR REPLACE TEMP TABLE keep_features (feature_id VARCHAR)")
      DBI::dbWriteTable(con, "keep_features",
        data.frame(feature_id = keep_features),
        append = TRUE
      )


      copy_sql <- sprintf(
        "
        COPY (
          SELECT
            f.\"genome.genome_id\" AS genome_id,
            %s AS feature_id,
            MAX(CAST(%s AS DOUBLE)) AS value,
            resistant_classes
          FROM %s
          JOIN selected_genomes USING (genome_id)
          JOIN keep_features kf ON %s = kf.feature_id
          JOIN metadata f ON genome_id = f.\"genome.genome_id\"
          WHERE resistant_classes <> 'Intermediate'
          GROUP BY
            f.\"genome.genome_id\",
            %s,
            resistant_classes
          ORDER BY
            f.\"genome.genome_id\",
            %s
        )
        TO '%s'
        (FORMAT 'parquet', COMPRESSION 'zstd')
      ",
        fid, # %s -> feature_id expression column name
        value_col, # %s -> value column to CAST
        mview, # %s -> source view (binary or counts)
        fid, # %s -> join to keep_features
        fid, # %s -> group by feature id
        fid, # %s -> order by feature id
        out_file_sql # %s -> destination parquet file
      )

      ok <- try(DBI::dbExecute(con, copy_sql), silent = TRUE)
      if (inherits(ok, "try-error")) {
        log("info", paste0("COPY failed for MDR matrix: ", out_file))
      } else {
        log("info", paste0("Exported MDR matrix: ", out_file))
      }

      DBI::dbDisconnect(con, shutdown = FALSE)
    }
  }

  log("info", "The MDR matrix has been generated and saved.")
  invisible(TRUE)
}

# Normalize split (support split = 0 for pure CV)
.normalize_split <- function(split) {
  if (length(split) == 1) {
    if (is.numeric(split) && split == 0) {
      return(c(1, 0))
    }
    stop("If 'split' is not a vector, only 0 is allowed (this will enable CV).")
  }
  split
}

#' Generate all ML input matrices
#'
#' a) each drug/drugclass matrices for all data and feature types
#' b) year bin holdouts for drug/class
#' c) country holdouts for drug/class
#' d) Leave-one-out from the years and countries for drug/class
#' e) MDR based on classes
#'
#' @param parquet_duckdb_path [character] path to the DuckDB that contains the view of metadata and feature parquets
#' @param out_path [character] path to the directory where the results files (matrices) will be written
#' @param n_fold [numeric] number of cross-validation folds; default is 5
#' @param split [numeric] training/validation split specification. Two formats accepted:
#'   - Shorthand for CV: `split = 0` (converted internally to `c(1, 0)`)
#'   - Vector form: `c(train_prop, val_prop)` where test_prop = 1 - train - val
#'     * For CV: `c(1, 0)` means 80% training data with k-fold CV, 20% stratified testing
#'     * For classical splits: all three partitions must be > 0
#'       Example: `c(0.7, 0.15)` = 70% train, 15% val, 15% test
#' @param min_n [numeric] minimum number of samples for each combination of drug classes for MDR matrix; default is 25
#' @param verbosity [character] "minimal" or "debug"; when "debug", prints full diagnostics per matrix
#' @return invisible(TRUE) on success; side effects: writes matrices/files
#' @examples
#' \dontrun{
#' # Generate ML input matrices with 5-fold cross-validation (using shorthand)
#' generateMLInputs(
#'   parquet_duckdb_path = "results/Cje_parquet.duckdb",
#'   out_path = "results/",
#'   n_fold = 5,
#'   split = 0, # shorthand for CV mode
#'   min_n = 25,
#'   verbosity = "minimal"
#' )
#'
#' # Same as above but using vector notation
#' generateMLInputs(
#'   parquet_duckdb_path = "results/Cje_parquet.duckdb",
#'   out_path = "results/",
#'   n_fold = 5,
#'   split = c(1, 0), # explicit CV mode
#'   verbosity = "minimal"
#' )
#'
#' # Generate with classical train/val/test split (70/15/15)
#' generateMLInputs(
#'   parquet_duckdb_path = "results/Cje_parquet.duckdb",
#'   out_path = "results/",
#'   n_fold = NULL,
#'   split = c(0.7, 0.15), # 70% train, 15% val, 15% test
#'   min_n = 25,
#'   verbosity = "debug"
#' )
#' }
#' @export
generateMLInputs <- function(parquet_duckdb_path = "results/Cje_parquet.duckdb",
                             out_path = "results/",
                             n_fold = 5,
                             split = c(1, 0), # Default: CV
                             min_n = 25,
                             verbosity = c("minimal", "debug")) {
  verbosity <- match.arg(verbosity)
  log <- .make_logger(verbosity)

  # Validate inputs before processing
  if (!file.exists(parquet_duckdb_path)) {
    stop("DuckDB file not found: ", parquet_duckdb_path)
  }

  if (!dir.exists(dirname(out_path))) {
    stop("Output directory does not exist: ", dirname(out_path))
  }

  # Normalize input paths
  parquet_duckdb_path <- normalizePath(parquet_duckdb_path)
  path <- normalizePath(out_path)

  .checkArgPath(parquet_duckdb_path)

  split <- .normalize_split(split)

  # Enforce mode with upstream arg checks (CV vs. train/val/test)
  mode <- .checkArgSplitAndMode(split, n_fold) # returns "cv" or "splits"
  train_prop <- split[1]
  val_prop <- split[2]
  test_prop <- 1 - train_prop - val_prop

  # Logging how we split with a wee .json file so ML scripts can read it in
  params <- list(
    mode        = mode,
    n_fold      = if (mode == "cv") n_fold else NULL,
    split       = split,
    train_prop  = train_prop,
    val_prop    = val_prop,
    test_prop   = test_prop,
    seed        = 5280,
    timestamp <- as.character(Sys.time())
  )

  jsonlite::write_json(
    params,
    file.path(out_path, "ml_parameters.json"),
    pretty = TRUE, auto_unbox = TRUE
  )

  # Log resolved mode
  if (identical(mode, "cv")) {
    log("info", paste0("Selected mode: CV (n_fold = ", n_fold, "), train = 1, val = 0, test = 0"))
  } else if (identical(mode, "splits")) {
    log("info", paste0(
      "Selected mode: splits with validation, train = ", train_prop,
      ", val = ", val_prop, ", test = ", round(test_prop, 3)
    ))
  }

  # Proceed (propagate verbosity to downstream steps)
  .parquet2Matrix(parquet_duckdb_path, out_path, n_fold, split,
    stratify_by = "year",
    verbosity = verbosity
  )
  .parquet2LOOMatrix(out_path,
    stratify_by = "year",
    verbosity = verbosity
  )

  .parquet2Matrix(parquet_duckdb_path, out_path, n_fold, split,
    stratify_by = "country",
    verbosity = verbosity
  )
  .parquet2LOOMatrix(out_path,
    stratify_by = "country",
    verbosity = verbosity
  )

  .parquet2Matrix(parquet_duckdb_path, out_path, n_fold, split,
    stratify_by = NULL,
    verbosity = verbosity
  )
  .parquet2MDRMatrix(parquet_duckdb_path, out_path, min_n,
    verbosity = verbosity
  )

  log("info", "All matrices generated and saved.")
  invisible(TRUE)
}
