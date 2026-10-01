# Reading an amRdata ORB.

# Names and formats set by amRdata, kept in one place.
# The only manifest schema this code can read; other versions are rejected.
.AMRDATA_SCHEMA_VERSION <- 1L
# Tells an amRdata dataset manifest apart from other JSON files.
.AMRDATA_MANIFEST_TYPE <- "amR_dataset"
# The producer-run stage whose outputs list every file in the ORB.
.AMRDATA_EXPORT_STAGE <- "clean_metadata_and_export"
# Name of the manifest entry amRdata sets to "ready" when its output can be used by amRml.
.AMRDATA_ML_ARTIFACT <- "amRml_input"
# File names amRdata gives its manifests, used to find them in the ORB directory.
.MANIFEST_PATTERN <- "^manifest_.*[.]json$"
# Columns every feature table has, alongside one feature ID column.
.FEATURE_TABLE_COLUMNS <- c("genome_id", "value")


# Tibble of the files the producer run's export stage lists. Manifest fields are
# read with `[[`, since `$` also matches a longer name that starts the same way.
.declaredFiles <- function(run) {
  stage <- Find(function(s) identical(s[["name"]], .AMRDATA_EXPORT_STAGE), run[["stages"]])

  if (is.null(stage)) {
    rlang::abort(paste0("The producer run has no '", .AMRDATA_EXPORT_STAGE, "' stage."))
  }

  outputs <- stage[["outputs"]]
  path <- vapply(outputs, function(o) o[["path"]] %||% NA_character_, character(1))
  size_bytes <- vapply(outputs, function(o) {
    size <- o[["size_bytes"]]
    if (is.numeric(size) && length(size) == 1L) size else NA_real_
  }, numeric(1))

  invalid <- is.na(path) | !nzchar(path) | !is.finite(size_bytes) | size_bytes < 0
  if (any(invalid)) {
    rlang::abort(paste0(
      "Declared file ", which(invalid)[[1]], " has a missing or invalid path or size."
    ))
  }

  tibble::tibble(
    path = path,
    name = basename(path),
    size_bytes = size_bytes,
    modified_at = vapply(outputs, function(o) o[["modified_at"]] %||% NA_character_, character(1))
  )
}

# The fields readORB() uses, checked; they become the ORB object's fields.
# Returns an error (not thrown) if the ORB isn't ready, and stops with a message
# if a field is missing or malformed.
.manifestRecord <- function(manifest_path, json) {
  string <- function(x, name) {
    if (!rlang::is_string(x) || !nzchar(x)) {
      rlang::abort(paste0("`", name, "` is missing or not a single string."))
    }
    x
  }

  supported <- is.list(json) &&
    identical(json[["schema_version"]], .AMRDATA_SCHEMA_VERSION) &&
    identical(json[["manifest_type"]], .AMRDATA_MANIFEST_TYPE)

  if (!supported) {
    rlang::abort(paste0(
      "Not an amRdata v", .AMRDATA_SCHEMA_VERSION, " '", .AMRDATA_MANIFEST_TYPE, "' manifest."
    ))
  }

  artifact <- json[["artifacts"]][[.AMRDATA_ML_ARTIFACT]]
  status <- artifact[["status"]] %||% "absent"

  if (!identical(status, "ready")) {
    return(.amrError(
      "orb_not_ready",
      paste0("'", .AMRDATA_ML_ARTIFACT, "' is ", status, ", not ready."),
      observed = list(status = status)
    ))
  }

  producer_run_id <- string(artifact[["producer_run_id"]], "producer_run_id")
  run <- Find(function(r) identical(r[["run_id"]], producer_run_id), json[["runs"]])
  run_status <- run[["status"]] %||% "missing"

  if (!identical(run_status, "success")) {
    return(.amrError(
      "producer_run_invalid",
      paste0("Producer run '", producer_run_id, "' is ", run_status, "."),
      observed = list(producer_run_id = producer_run_id, status = run_status)
    ))
  }

  finished_at <- string(run[["finished_at"]], "finished_at")

  if (!grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}", finished_at)) {
    rlang::abort(paste0("`finished_at` is not a timestamp: ", finished_at))
  }

  # The metadata file must be declared, so the file checks cover it.
  files <- .declaredFiles(run)
  metadata_parquet <- string(artifact[["metadata_parquet"]], "metadata_parquet")

  if (!metadata_parquet %in% files$path) {
    rlang::abort(paste0(
      "`metadata_parquet` is not among the declared files: ", basename(metadata_parquet)
    ))
  }

  dataset_id <- string(json[["dataset_id"]], "dataset_id")
  # What amRdata was asked for: species names or taxon IDs.
  requested <- unlist(json[["dataset"]][["selection"]][["user_bacs"]], use.names = FALSE)

  list(
    dataset_id = dataset_id,
    label = if (length(requested)) paste(requested, collapse = ", ") else dataset_id,
    directory = string(artifact[["directory"]], "directory"),
    manifest_path = manifest_path,
    metadata_parquet = metadata_parquet,
    producer = string(artifact[["producer"]], "producer"),
    producer_run_id = producer_run_id,
    finished_at = finished_at,
    files = files
  )
}

# Check one manifest file. Returns an error (not thrown) if it can't be used,
# otherwise its checked record from `.manifestRecord()`.
.manifestReadiness <- function(manifest_path) {
  json <- tryCatch(
    jsonlite::read_json(manifest_path, simplifyVector = FALSE),
    error = function(e) e
  )

  if (inherits(json, "error")) {
    return(.amrError(
      "manifest_unreadable",
      "The manifest is not readable JSON.",
      observed = list(error = conditionMessage(json))
    ))
  }

  # Anything wrong with the manifest's contents is reported as manifest_invalid,
  # not an R error.
  tryCatch(
    .manifestRecord(manifest_path, json),
    error = function(e) {
      .amrError(
        "manifest_invalid",
        c("The manifest can't be used.", x = conditionMessage(e)),
        observed = list(error = conditionMessage(e))
      )
    }
  )
}

# Of the usable manifests, the one whose producer run finished last; ties go to
# the first manifest file name. Timestamps are year first and zero-padded, so
# they sort correctly as text.
.newestReady <- function(ready) {
  finished <- vapply(ready, function(x) x$finished_at, character(1))
  files <- vapply(ready, function(x) basename(x$manifest_path), character(1))

  ready[[order(finished, files, decreasing = c(TRUE, FALSE), method = "radix")[[1]]]]
}

# Throw an error unless the ORB is where its manifest says, and its files exist
# at their recorded size.
.checkOrbFiles <- function(record, dir, call) {
  recorded_dir <- record$directory

  # Check this first, since in a moved ORB every declared file would look missing.
  if (!identical(normalizePath(recorded_dir, mustWork = FALSE), dir)) {
    .amrAbort(
      "orb_moved",
      c(
        "This ORB is not where its manifest says it was written.",
        x = paste0("Manifest records: ", recorded_dir),
        i = "amRdata records absolute paths, so a moved or copied ORB can't be read yet."
      ),
      observed = list(recorded = recorded_dir, actual = dir),
      call = call
    )
  }

  # A file that can't be read is as unusable as a missing one.
  files <- record$files
  absent <- files$path[file.access(files$path, mode = 4) != 0]

  if (length(absent)) {
    .amrAbort(
      "orb_files_missing",
      c("Files the manifest declares are missing or unreadable.", .bullets(basename(absent))),
      observed = list(missing = absent),
      call = call
    )
  }

  files$actual_bytes <- file.info(files$path)$size
  changed <- files[files$actual_bytes != files$size_bytes, , drop = FALSE]

  if (nrow(changed)) {
    listed <- paste0(
      changed$name, ": ", changed$size_bytes, " -> ", changed$actual_bytes, " bytes"
    )

    .amrAbort(
      "orb_file_changed",
      c("Declared files have changed size since amRdata wrote them.", .bullets(listed)),
      observed = list(
        file = changed$path,
        size_bytes = changed$size_bytes,
        actual_bytes = changed$actual_bytes
      ),
      call = call
    )
  }

  invisible(NULL)
}

# Tibble of `feature_type` and `path` for declared parquets with the columns of a
# feature table.
.featureTables <- function(files, call) {
  parquets <- files$path[grepl("[.]parquet$", files$name)]

  schemas <- lapply(parquets, function(f) {
    tryCatch(
      arrow::open_dataset(f, format = "parquet")$schema$names,
      error = function(e) {
        .amrAbort(
          "parquet_unreadable",
          c("A declared parquet file could not be read.", x = basename(f)),
          observed = list(file = f, error = conditionMessage(e)),
          call = call
        )
      }
    )
  })

  # Each table's feature ID column, or NA if it isn't a feature table.
  feature_type <- vapply(schemas, function(cols) {
    if (length(cols) == 3L && all(.FEATURE_TABLE_COLUMNS %in% cols)) {
      setdiff(cols, .FEATURE_TABLE_COLUMNS)
    } else {
      NA_character_
    }
  }, character(1))

  tables <- tibble::tibble(feature_type = feature_type, path = parquets)
  tables <- tables[!is.na(tables$feature_type), , drop = FALSE]

  if (!nrow(tables)) {
    .amrAbort(
      "orb_no_feature_tables",
      c(
        "The ORB declares no feature tables.",
        i = "Expected columns: genome_id, value, and one feature ID."
      ),
      call = call
    )
  }

  # Later stages use feature_type as a key, so two tables can't share one.
  duplicated_types <- unique(tables$feature_type[duplicated(tables$feature_type)])

  if (length(duplicated_types)) {
    .amrAbort(
      "orb_duplicate_feature_type",
      paste0("Feature tables share a feature ID column: ", toString(duplicated_types), "."),
      observed = list(feature_type = duplicated_types),
      call = call
    )
  }

  # "radix" sorts the same way in every locale.
  tables[order(tables$feature_type, method = "radix"), , drop = FALSE]
}

#' Read an amRdata ORB
#'
#' Validates the ORB in `path` and returns an object describing it.
#'
#' @details
#' Uses the manifest whose `amRml_input` is ready and whose producer run
#' finished last, with a message if more than one is ready. The ORB must be
#' where that manifest says, with every declared file present at its recorded
#' size. If it is not, `readORB()` stops with an error; it never falls back to
#' an older manifest, which would read older data. Feature tables are the
#' declared parquet files with the columns `genome_id`, `value` and one feature
#' ID column such as `gene` or `Pfam`.
#'
#' @section Errors:
#' All errors have class `amrml_error`, plus one of:
#'
#' * `amrml_invalid_argument`: `path` is not a single directory path.
#' * `amrml_orb_dir_not_found`: `path` does not exist.
#' * `amrml_manifest_not_found`: the directory holds no `manifest_*.json`.
#' * `amrml_orb_no_usable_manifest`: no manifest is usable. `observed` holds
#'   each manifest's error.
#' * `amrml_orb_moved`: the ORB is not where its manifest says it was written.
#' * `amrml_orb_files_missing`: declared files are missing or unreadable.
#' * `amrml_orb_file_changed`: a declared file's size has changed.
#' * `amrml_parquet_unreadable`: a declared parquet cannot be read.
#' * `amrml_orb_no_feature_tables`: no declared parquet is a feature table.
#' * `amrml_orb_duplicate_feature_type`: two feature tables share a feature ID
#'   column.
#'
#' @param path Character. The ORB directory: the folder amRdata wrote its
#'   output to, containing `manifest_*.json` (e.g.
#'   `data/Staphylococcus_argenteus`), not its parent.
#'
#' @return An `amr_orb` list: `dataset_id`, `label` (what amRdata was asked for,
#'   species names or taxon IDs, otherwise `dataset_id`), `directory`,
#'   `manifest_path`, `metadata_parquet`, `producer`, `producer_run_id`,
#'   `finished_at`, `files` (the producer run's declared outputs) and
#'   `feature_tables` (`feature_type` and `path`).
#'
#' @examples
#' \dontrun{
#' orb <- readORB("path/to/amRdata/output")
#' orb$feature_tables
#' }
#' @export
readORB <- function(path) {
  call <- rlang::current_env()

  # Check the path is a single existing directory.
  if (missing(path) || !rlang::is_string(path) || !nzchar(path)) {
    .amrAbort("invalid_argument", "`path` must be a single ORB directory path.", call = call)
  }

  if (!dir.exists(path)) {
    .amrAbort(
      "orb_dir_not_found",
      c("ORB directory not found.", x = path),
      call = call
    )
  }

  # Find the manifests in the ORB directory.
  dir <- normalizePath(path, mustWork = TRUE)
  manifest_paths <- list.files(dir, pattern = .MANIFEST_PATTERN, full.names = TRUE)

  if (!length(manifest_paths)) {
    .amrAbort(
      "manifest_not_found",
      c("No amRdata manifest in the ORB directory.", x = dir),
      call = call
    )
  }

  # Check each manifest; stop with every manifest's error if none is usable.
  manifests <- lapply(manifest_paths, .manifestReadiness)
  names(manifests) <- basename(manifest_paths)
  usable <- !vapply(manifests, inherits, logical(1), "amrml_error")

  if (!any(usable)) {
    reasons <- paste0(
      names(manifests), ": ",
      vapply(manifests, function(e) paste(e$message, collapse = " "), character(1))
    )

    .amrAbort(
      "orb_no_usable_manifest",
      c("No manifest in the ORB declares a usable amRml input.", .bullets(reasons)),
      observed = manifests,
      call = call
    )
  }

  # Use the newest usable manifest.
  ready <- manifests[usable]
  chosen <- .newestReady(ready)

  if (length(ready) > 1L) {
    rlang::inform(paste0(
      "Using manifest ", basename(chosen$manifest_path), ", the newest of ",
      length(ready), " ready manifests."
    ))
  }

  # Check the ORB hasn't moved and its files exist at their recorded size.
  .checkOrbFiles(chosen, dir, call = call)

  # Find the feature tables among the declared files.
  feature_tables <- .featureTables(chosen$files, call = call)

  # Build the ORB object from the checked record.
  orb <- chosen
  orb$directory <- dir
  orb$feature_tables <- feature_tables
  structure(orb, class = "amr_orb")
}


#' @export
print.amr_orb <- function(x, ...) {
  cat("<amr_orb>", x$dataset_id, "-", x$label, "\n")
  cat("  directory :", x$directory, "\n")
  cat("  manifest  :", basename(x$manifest_path), "\n")
  cat("  finished  :", x$finished_at, "\n")
  cat("  declares  :", nrow(x$files), "files\n")
  cat("  features  :", paste(x$feature_tables$feature_type, collapse = ", "), "\n")
  invisible(x)
}
