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
  if (anyNA(table$value) || any(table$value < 0)) {
    .amrAbort(
      "feature_values_invalid",
      paste0("The ", feature_type, " table has missing or negative values."),
      observed = list(path = path, n_missing = sum(is.na(table$value))),
      call = call
    )
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

# Values in an encoding: binary is 1 if the value is above 0; counts are the values.
.encode <- function(value, encoding) {
  if (encoding == "binary") as.numeric(value > 0) else as.numeric(value)
}

# Features whose encoded value differs among the training genomes; a missing row counts as 0.
.variableFeatures <- function(table, train_ids, encoding) {
  present <- table[table$genome_id %in% train_ids & table$value != 0, , drop = FALSE]
  present$value <- .encode(present$value, encoding)
  stats <- dplyr::summarise(
    dplyr::group_by(present, .data$feature_id),
    n_present = dplyr::n(),
    n_values = dplyr::n_distinct(.data$value),
    .groups = "drop"
  )

  # A feature absent from every training genome has no row here, so is never kept.
  constant <- stats$n_present == length(train_ids) & stats$n_values == 1
  sort(stats$feature_id[!constant], method = "radix")
}

# One long matrix: `rows` genomes × features varying among its training genomes. NULL if none.
.buildMatrix <- function(table, rows, encoding) {
  features <- .variableFeatures(table, rows$genome_id[rows$role == "train"], encoding)
  if (!length(features)) {
    return(NULL)
  }

  rows <- rows[order(rows$genome_id, method = "radix"), , drop = FALSE]
  n_features <- length(features)

  # Put each non-zero value in a grid of zeros, indexing the shared table without copying it.
  feature <- match(table$feature_id, features)
  genome <- match(table$genome_id, rows$genome_id)
  kept <- !is.na(feature) & !is.na(genome) & table$value != 0
  value <- numeric(nrow(rows) * n_features)
  value[(genome[kept] - 1) * n_features + feature[kept]] <- .encode(table$value[kept], encoding)

  tibble::tibble(
    genome_id = rep(rows$genome_id, each = n_features),
    role = rep(rows$role, each = n_features),
    phenotype = rep(rows$phenotype, each = n_features),
    feature_id = rep(features, times = nrow(rows)),
    value = value
  )
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
#' Writes each planned matrix to `out_dir/matrices/`, and the list of matrices
#' to `out_dir/matrices.parquet`. A matrix has one row per genome and per
#' feature that varies among its training genomes (`genome_id`, `role`,
#' `phenotype`, `feature_id`, `value`); a matrix with no varying feature is
#' skipped. amRdata's files in the ORB are only read.
#'
#' @param plan An `amr_matrix_plan` from [planMatrices()].
#' @param out_dir Folder to write to, created if needed. By default a new run
#'   folder in the ORB's `amRml/` folder, named by date and time. Elsewhere in
#'   the ORB is refused.
#' @param overwrite Replace matrices already in `out_dir`.
#'
#' @return An `amr_matrices` list: `eligibility`, `out_dir`, and `matrices`, the
#'   plan's rows plus `status` (`"built"`/`"skipped"`), `rule_id`, `n_features`
#'   and `file` (relative to `out_dir`).
#'
#' @section Errors:
#' All errors have class `amrml_error`, plus one of:
#'
#' * `amrml_invalid_argument`: `plan` is not from [planMatrices()], `out_dir`
#'   is not one writable path outside the ORB or in its `amRml/` folder, or it
#'   already has output (eligibility, matrices, tasks or results) and
#'   `overwrite` is `FALSE`.
#' * `amrml_feature_rows_repeated`: a feature table has two rows for one genome
#'   and feature.
#' * `amrml_feature_values_invalid`: a feature table has missing or negative
#'   values.
#'
#' @examples
#' \dontrun{
#' built <- buildMatrices(plan, "path/to/output")
#' built$matrices
#' }
#' @export
buildMatrices <- function(plan, out_dir = NULL, overwrite = FALSE) {
  call <- rlang::current_env()
  .checkArgClass(plan, "plan", "amr_matrix_plan", "planMatrices()", call = call)
  orb <- plan$eligibility$profile$orb
  out_dir <- out_dir %||% .defaultOutDir(orb$directory)
  .checkArgOutDir(out_dir, orb$directory, overwrite, call = call)

  # Build in a new staging folder and swap it in at the end; a failure while building keeps
  # earlier output, though the final delete-then-rename is not atomic.
  staging <- tempfile(".matrices-building-", tmpdir = out_dir)
  dir.create(file.path(staging, "matrices"), recursive = TRUE)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)

  matrices <- plan$matrices
  matrices$status <- "built"
  matrices$rule_id <- NA_character_
  matrices$n_features <- NA_integer_
  matrices$file <- NA_character_
  members <- plan$eligibility$members

  # Each feature table is read once, for all its matrices.
  for (feature_type in unique(matrices$feature_type)) {
    path <- orb$feature_tables$path[orb$feature_tables$feature_type == feature_type]
    table <- .readFeatureTable(path, feature_type, unique(members$genome_id), call = call)

    for (i in which(matrices$feature_type == feature_type)) {
      scope <- matrices[i, c("mode_id", .SCOPE_KEYS)]
      rows <- dplyr::semi_join(members, scope, by = c("mode_id", .SCOPE_KEYS))
      matrix <- .buildMatrix(table, rows, matrices$encoding[[i]])

      if (is.null(matrix)) {
        matrices$status[[i]] <- "skipped"
        matrices$rule_id[[i]] <- "no_variable_features"
        next
      }

      file <- file.path("matrices", paste0(matrices$matrix_id[[i]], ".parquet"))
      arrow::write_parquet(matrix, file.path(staging, file), compression = "zstd")
      matrices$n_features[[i]] <- nrow(matrix) %/% nrow(rows)
      matrices$file[[i]] <- file
    }
  }

  arrow::write_parquet(matrices, file.path(staging, "matrices.parquet"))

  # Only what this function writes is replaced.
  outputs <- c("matrices", "matrices.parquet")
  unlink(file.path(out_dir, outputs), recursive = TRUE)
  if (!all(file.rename(file.path(staging, outputs), file.path(out_dir, outputs)))) {
    rlang::abort(
      paste0("Could not move the built matrices into '", out_dir, "'."),
      class = c("amrml_internal_error", "amrml_error")
    )
  }

  structure(
    list(eligibility = plan$eligibility, out_dir = out_dir, matrices = matrices),
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
