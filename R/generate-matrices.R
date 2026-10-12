# Planning and building the ML input matrices of eligible scopes.


# Encodings per feature table: binary, plus counts if any value in the table is above 1.
.encodings <- function(feature_tables) {
  has_counts <- vapply(feature_tables$path, function(path) {
    values <- arrow::read_parquet(path, col_select = "value", as_data_frame = FALSE)$value
    isTRUE(as.vector(max(values, na.rm = TRUE)) > 1)
  }, logical(1), USE.NAMES = FALSE)

  tibble::tibble(
    feature_type = rep(feature_tables$feature_type, 1L + has_counts),
    encoding = unlist(lapply(has_counts, function(counts) c("binary", if (counts) "counts")))
  )
}

# Readable, file-safe IDs, numbered because file-safe names can collide; nothing parses them.
.matrixIds <- function(identity) {
  parts <- lapply(identity, function(x) ifelse(is.na(x), "-", x))
  names <- gsub("[^A-Za-z0-9._-]", "_", do.call(paste, c(parts, sep = "__")))
  sprintf("%04d__%s", seq_along(names), names)
}

# A feature table as `genome_id`, `feature_id` and `value`, for the genomes asked for.
.readFeatureTable <- function(path, feature_type, genomes, call) {
  columns <- c("genome_id", feature_type, "value")
  table <- arrow::read_parquet(path, col_select = dplyr::all_of(columns))
  names(table) <- c("genome_id", "feature_id", "value")
  table <- table[table$genome_id %in% genomes, , drop = FALSE]

  # Encoding and the zero-fill treat each value as a count: present, and 0 or more.
  if (any(!is.finite(table$value)) || any(table$value < 0)) {
    .amrAbort(
      "feature_values_invalid",
      paste0("The ", feature_type, " table has missing, nonfinite or negative values."),
      observed = list(path = path, n_missing = sum(is.na(table$value))),
      call = call
    )
  }

  if (anyNA(table$genome_id) || anyNA(table$feature_id)) {
    .amrAbort("feature_values_invalid", "Feature rows need genome and feature IDs.", call = call)
  }

  # One value per genome and feature; duplicated() is slow, so it only runs to name the pairs.
  if (dplyr::n_distinct(table$genome_id, table$feature_id) < nrow(table)) {
    pairs <- table[duplicated(table[c("genome_id", "feature_id")]), c("genome_id", "feature_id")]
    .amrAbort(
      "feature_rows_repeated",
      c(
        paste0("The ", feature_type, " table has more than one row for a genome and feature."),
        .bullets(utils::head(paste(pairs$genome_id, pairs$feature_id), 5))
      ),
      observed = list(path = path, pairs = pairs),
      call = call
    )
  }

  table
}

# A shared sparse feature table, include genome rows without observed features.
# The stored matrix is counts. Binary is derived only when scope needs it
.featureStore <- function(table, genomes) {
  features <- sort(unique(table$feature_id), method = "radix")
  values <- Matrix::sparseMatrix(
    i = match(table$genome_id, genomes),
    j = match(table$feature_id, features),
    x = as.numeric(table$value),
    dims = c(length(genomes), length(features)),
    giveCsparse = TRUE
  )

  list(
    values = Matrix::drop0(values),
    genomes = genomes,
    features = features
  )
}

# Adding indexing ~magic~ to .variableFeatures(), avoids full scan feature table.
# Features vary if present in some, but not all, training genomes, or (for counts)
# has different (positive) counts in every training genome
.variableFeaturesIndexed <- function(store, train_ids, encoding) {
  train <- store$values[match(train_ids, store$genomes), , drop = FALSE]
  n_present <- Matrix::colSums(train != 0)
  variable <- n_present > 0 & n_present < nrow(train)

  if (encoding == "counts") {
    # Only fully present columns need a count-variability check
    full <- which(n_present == nrow(train))
    p <- train@p
    x <- train@x
    for (j in full) {
      positions <- seq.int(p[[j]] + 1L, p[[j + 1L]])
      variable[[j]] <- any(x[positions] != x[positions[[1L]]])
    }
  }
  which(variable)
}

# Scope membership computed just once
.scopeIndex <- function(matrices, members) {
  keys <- c("mode_id", .SCOPE_KEYS)
  scopes <- dplyr::distinct(matrices[keys])
  scopes$scope_index <- seq_len(nrow(scopes))
  member_scopes <- dplyr::left_join(members, scopes, by = keys)$scope_index
  if (anyNA(member_scopes)) stop("A member has no planned scope.")
  groups <- split(seq_len(nrow(members)), factor(member_scopes, levels = scopes$scope_index))
  rows <- lapply(groups, function(i) {
    selected <- members[i, c("genome_id", "role", "phenotype"), drop = FALSE]
    selected[order(selected$genome_id, method = "radix"), , drop = FALSE]
  })
  matrix_scopes <- dplyr::left_join(matrices[keys], scopes, by = keys)$scope_index
  list(rows = rows, matrix_scopes = matrix_scopes)
}

#' Plan the matrices to build
#'
#' Lists every matrix the eligible scopes need, before any is built: one per
#' eligible scope, feature type and encoding. Each feature type is encoded as
#' binary, and also as counts when its table has values above 1.
#'
#' @param eligibility An `amr_eligibility` from [eligibleScopes()].
#'
#' @return An `amr_matrix_plan` list:
#'   * `eligibility`: the eligibility used.
#'   * `matrices`: one row per matrix: `matrix_id`, the scope (`mode_id`, `unit`,
#'     `target`, `train_group`, `test_group`), `feature_type`, `encoding`, and
#'     the scope's `n_train` and `n_test` genomes (`NA` without a test set).
#'
#' @section Errors:
#' All errors have class `amrml_error`, plus `amrml_invalid_argument` if
#' `eligibility` is not from [eligibleScopes()].
#'
#' @examples
#' \dontrun{
#' plan <- planMatrices(eligibleScopes(profileORB(readORB("path/to/amRdata/output"))))
#' plan$matrices
#' }
#' @export
planMatrices <- function(eligibility) {
  call <- rlang::current_env()
  .checkArgClass(eligibility, "eligibility", "amr_eligibility", "eligibleScopes()", call = call)

  scopes <- eligibility$scopes[eligibility$scopes$eligible, , drop = FALSE]
  scopes <- tibble::tibble(
    scopes[c("mode_id", .SCOPE_KEYS)],
    n_train = scopes$n_genomes,
    n_test = scopes$test_n_genomes
  )
  encodings <- .encodings(eligibility$profile$orb$feature_tables)

  keys <- c("mode_id", .SCOPE_KEYS, "feature_type", "encoding")
  matrices <- dplyr::cross_join(scopes, encodings)
  matrices <- tibble::tibble(
    matrix_id = .matrixIds(matrices[keys]),
    matrices[c(keys, "n_train", "n_test")]
  )

  structure(
    list(eligibility = eligibility, matrices = matrices),
    class = "amr_matrix_plan"
  )
}

#' Build the planned matrices
#'
#' Writes each sparse count matrix per feature type in `features/`.
#' Binary encoding and scope-specific feature selection are indexed without
#' materializing matrices for speed. amRdata's files in the ORB are only read.
#'
#' @param plan An `amr_matrix_plan` from [planMatrices()].
#' @param out_dir Folder to write to, created if needed.
#' @param overwrite Replace matrices already in `out_dir`.
#'
#' @return An `amr_matrices` list: `eligibility`, `out_dir`, `matrices`, and
#'   `storage = "shared_sparse"`. The `file` column points to a shared feature
#'   store rather than a scope-specific Parquet file. `matrix_index.rds` records
#'   membership rows and selected feature column indices per matrix
#'
#' @section Errors:
#' All errors have class `amrml_error`, plus one of:
#'
#' * `amrml_invalid_argument`: invalid `plan`/`out_dir`, or output exists without overwrite.
#' * `amrml_feature_rows_repeated`: duplicate genome/feature pairs.
#' * `amrml_feature_values_invalid`: missing, nonfinite, or negative feature values.
#'
#' @export
buildMatrices <- function(plan, out_dir = NULL, overwrite = FALSE) {
  call <- rlang::current_env()
  .checkArgClass(plan, "plan", "amr_matrix_plan", "planMatrices()", call = call)
  orb <- plan$eligibility$profile$orb
  out_dir <- out_dir %||% .defaultOutDir(orb$directory)
  .checkArgOutDir(out_dir, orb$directory, overwrite, call = call)

  # Preserve outputs until the new index and feature stores exist
  staging <- tempfile(".matrices-building-", tmpdir = out_dir)
  dir.create(file.path(staging, "features"), recursive = TRUE)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)

  matrices <- plan$matrices
  matrices$status <- "built"
  matrices$rule_id <- NA_character_
  matrices$n_features <- NA_integer_
  matrices$file <- NA_character_
  members <- plan$eligibility$members
  scopes <- .scopeIndex(matrices, members)
  columns <- vector("list", nrow(matrices))
  genomes <- sort(unique(members$genome_id), method = "radix")

  # Construction happens once per feature type
  for (j in seq_along(unique(matrices$feature_type))) {
    feature_type <- unique(matrices$feature_type)[[j]]
    path <- orb$feature_tables$path[orb$feature_tables$feature_type == feature_type]
    table <- .readFeatureTable(path, feature_type, genomes, call = call)
    store <- .featureStore(table, genomes)
    rm(table)

    file <- file.path("features", sprintf("%04d.rds", j))
    saveRDS(store, file.path(staging, file), compress = FALSE)
    for (i in which(matrices$feature_type == feature_type)) {
      rows <- scopes$rows[[scopes$matrix_scopes[[i]]]]
      selected <- .variableFeaturesIndexed(store,
                                           rows$genome_id[rows$role == "train"],
                                           matrices$encoding[[i]]
                                           )
      columns[[i]] <- selected
      if (!length(selected)) {
        matrices$status[[i]] <- "skipped"
        matrices$rule_id[[i]] <- "no_variable_features"
        next
      }
      matrices$n_features[[i]] <- length(selected)
      matrices$file[[i]] <- file
    }
    rm(store)
  }

  saveRDS(
    list(rows = scopes$rows, matrix_scopes = scopes$matrix_scopes, columns = columns),
    file.path(staging, "matrix_index.rds"), compress = FALSE
  )
  arrow::write_parquet(matrices, file.path(staging, "matrices.parquet"))

  # Move the previous output
  previous <- c("matrices", "features", "matrices.parquet", "matrix_index.rds")
  produced <- c("features", "matrices.parquet", "matrix_index.rds")
  if (!overwrite && any(file.exists(file.path(out_dir, previous)))) {
    .amrAbort("invalid_argument", "Matrix outputs already exist; set `overwrite = TRUE`.", call = call)
  }
  backup <- tempfile(".matrices-previous-", tmpdir = out_dir)
  dir.create(backup)
  old_moved <- character()
  new_moved <- character()
  tryCatch({
    for (name in previous) {
      if (file.exists(file.path(out_dir, name))) {
        if (!file.rename(file.path(out_dir, name), file.path(backup, name))) {
          stop("Could not move previous matrix output: ", name)
        }
        old_moved <- c(old_moved, name)
      }
    }
    for (name in produced) {
      if (!file.rename(file.path(staging, name), file.path(out_dir, name))) {
        stop("Could not install new matrix output: ", name)
      }
      new_moved <- c(new_moved, name)
    }
  }, error = function(e) {
    for (name in rev(new_moved)) unlink(file.path(out_dir, name), recursive = TRUE)
    restored <- vapply(rev(old_moved), function(name) {
      file.rename(file.path(backup, name), file.path(out_dir, name))
    }, logical(1))
    if (all(restored)) unlink(backup, recursive = TRUE)
    stop(conditionMessage(e), if (!all(restored)) paste0(
      " Previous outputs may be recoverable from: ", backup
    ), call. = FALSE)
  })
  unlink(backup, recursive = TRUE)

  structure(
    list(eligibility = plan$eligibility, out_dir = out_dir,
         matrices = matrices, storage = "shared_sparse"),
    class = "amr_matrices"
  )
}


# How many of `matrices` were built and skipped.
.builtAndSkipped <- function(matrices) {
  paste0(sum(matrices$status == "built"), " built, ", sum(matrices$status == "skipped"), " skipped")
}

#' @export
print.amr_matrix_plan <- function(x, ...) {
  orb <- x$eligibility$profile$orb
  m <- x$matrices
  n_scopes <- nrow(dplyr::distinct(m[c("mode_id", .SCOPE_KEYS)]))

  cat("<amr_matrix_plan>", orb$dataset_id, "-", orb$dataset_label, "\n")
  cat("  ", nrow(m), " matrices from ", n_scopes, " eligible scopes\n", sep = "")
  n <- table(factor(m$mode_id, unique(m$mode_id)))
  cat(.modeLines(names(n), sprintf("%d matrices", n)), sep = "\n")
  invisible(x)
}

#' @export
print.amr_matrices <- function(x, ...) {
  orb <- x$eligibility$profile$orb
  m <- x$matrices

  cat("<amr_matrices>", orb$dataset_id, "-", orb$dataset_label, "\n")
  cat("  in", x$out_dir, "\n")
  cat("  ", .builtAndSkipped(m), "\n", sep = "")
  by_mode <- split(m, factor(m$mode_id, unique(m$mode_id)))
  cat(.modeLines(names(by_mode), vapply(by_mode, .builtAndSkipped, character(1))), sep = "\n")
  invisible(x)
}
