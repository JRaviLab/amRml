#########################
# BiocFileCache helpers #
#########################

#' Shared BFC used across the amR package suite
#'
#' Makes sure BiocFileCache exists, and if it does, opens the cache and
#' sets the behavior to create directories silently without interrupting a user
#' 's lunch to ask if they should create each new directory.
#'
#' @return A `BiocFileCache` object
#' @keywords internal
.amr_bfc <- function() {
  if (!requireNamespace("BiocFileCache", quietly = TRUE)) {
    stop(
      "Package 'BiocFileCache' is required for amR dataset discovery."
    )
  }

  BiocFileCache::BiocFileCache(ask = FALSE)
}

#' Find dataset manifests registered across the amR suite
#'
#' @return BFC-registered amRdata manifests on the system
#' @keywords internal
.amr_registered_manifests <- function() {
  bfc <- .amr_bfc()

  BiocFileCache::bfcquery(
    bfc,
    query = "^amR_dataset_manifest_",
    field = "rname",
    exact = FALSE
  )
}

#' Discover completed amRdata datasets available for amRml
#'
#' @return Existing amRdata datasets that are ready for modeling
#' @keywords internal
.discoverAmrDatasets <- function() {
  bfc <- .amr_bfc()
  hits <- .amr_registered_manifests()

  if (!nrow(hits)) {
    return(tibble::tibble())
  }

  datasets <- purrr::map_dfr(
    seq_len(nrow(hits)),
    function(i) {
      manifest_path <- tryCatch(
        BiocFileCache::bfcrpath(
          bfc,
          rids = hits$rid[[i]],
          exact = TRUE
        ),
        error = function(e) NA_character_
      )

      if (
        length(manifest_path) != 1L ||
        is.na(manifest_path) ||
        !file.exists(manifest_path)
      ) {
        return(NULL)
      }

      manifest <- tryCatch(
        jsonlite::read_json(
          manifest_path,
          simplifyVector = FALSE
        ),
        error = function(e) NULL
      )

      if (is.null(manifest)) {
        return(NULL)
      }

      ml_input <- tryCatch(
        .manifest_ml_input(manifest),
        error = function(e) NULL
      )

      if (is.null(ml_input)) {
        return(NULL)
      }

      artifact <- ml_input$artifact
      producer_run <- ml_input$producer_run

      selection <- manifest$dataset$selection$user_bacs

      label <- if (
        is.null(selection) ||
        !length(selection)
      ) {
        manifest$dataset_id
      } else {
        paste(
          unlist(selection, use.names = FALSE),
          collapse = ", "
        )
      }

      tibble::tibble(
        dataset_id = manifest$dataset_id,
        manifest_id = manifest$manifest_id,
        label = label,
        completed_at = producer_run$finished_at,
        parquet_dir = normalizePath(
          artifact$directory,
          mustWork = TRUE
        ),
        parquet_duckdb = normalizePath(
          artifact$parquet_duckdb,
          mustWork = TRUE
        ),
        metadata_parquet = normalizePath(
          artifact$metadata_parquet,
          mustWork = TRUE
        ),
        manifest_path = normalizePath(
          manifest_path,
          mustWork = TRUE
        ),
        bfc_rid = hits$rid[[i]]
      )
    }
  )

  if (!nrow(datasets)) {
    return(datasets)
  }

  # Multiple manifests can describe the same physical dataset, so keep the newest
  # usable one for each dataset directory
  datasets |>
    dplyr::arrange(
      dplyr::desc(.data$completed_at)
    ) |>
    dplyr::distinct(
      .data$parquet_dir,
      .keep_all = TRUE
    )
}

#' Choose an amRdata dataset if no explicit path was supplied
#'
#' @return A user-selected dataset choice for modeling.
#' @keywords internal
.selectAmrDataset <- function() {
  datasets <- .discoverAmrDatasets()

  if (!nrow(datasets)) {
    stop(
      "No completed amRdata datasets were found in BiocFileCache.\n",
      "Run the amRdata workflow first or provide `parquet_dir` explicitly."
    )
  }

  if (nrow(datasets) == 1L) {
    return(datasets[1, , drop = FALSE])
  }

  if (!interactive()) {
    stop(
      "Multiple completed amRdata datasets were found. ",
      "You can provide `parquet_dir` explicitly for non-interactive use."
    )
  }

  choices <- paste0(
    datasets$label,
    " | completed ",
    datasets$completed_at,
    " | ",
    datasets$parquet_dir
  )

  choice <- utils::menu(
    choices,
    title = "Please select an amRdata dataset for modeling:"
  )

  if (choice == 0L) {
    stop("No dataset selected.")
  }

  datasets[choice, , drop = FALSE]
}



#########################
#   Manifest helpers    #
#########################

#' Returns the basics about a file for manifest logging
#'
#' @param path Character vector of file paths.
#' @param hash Logical. If TRUE, calculate SHA-256 checksums.
#'
#' @return A list of file records.
#' @keywords internal
.manifest_file_info <- function(path, hash = FALSE) {
  path <- unique(as.character(path))
  path <- path[nzchar(path)]

  if (!length(path)) {
    return(list())
  }

  # See what exists
  purrr::map(path, function(x) {
    exists <- file.exists(x)

    out <- list(
      path = x,
      exists = exists,
      size_bytes = if (exists) file.info(x)$size else NA_real_,
      modified_at = if (exists) as.character(file.info(x)$mtime) else NA_character_
    )

    # Hash what exists, if desired
    if (isTRUE(hash) && exists && !dir.exists(x)) {
      out$sha256 <- unname(tools::sha256(x))
    }

    out
  })
}


#' Capture basic GitHub repo state for manifest provenance
#'
#' @param base_dir Character. Project root.
#'
#' @return A named list.
#' @keywords internal
.manifest_git_info <- function(base_dir = ".") {
  base_dir <- normalizePath(base_dir, mustWork = FALSE)

  # Find Git
  git <- Sys.which("git")

  if (!nzchar(git)) {
    return(list(
      available = FALSE
    ))
  }

  # Run Git through system commands
  run_git <- function(args) {
    tryCatch(
      system2(
        git,
        args = args,
        stdout = TRUE,
        stderr = FALSE
      ),
      error = function(e) character()
    )
  }

  inside <- run_git(c("-C", shQuote(base_dir), "rev-parse", "--is-inside-work-tree"))

  if (!length(inside) || !identical(trimws(inside[[1]]), "true")) {
    return(list(
      available = TRUE,
      repository = FALSE
    ))
  }

  commit <- run_git(c("-C", shQuote(base_dir), "rev-parse", "HEAD"))
  branch <- run_git(c("-C", shQuote(base_dir), "rev-parse", "--abbrev-ref", "HEAD"))
  dirty <- run_git(c("-C", shQuote(base_dir), "status", "--porcelain"))

  list(
    available = TRUE,
    repository = TRUE,
    commit = if (length(commit)) trimws(commit[[1]]) else NA_character_,
    branch = if (length(branch)) trimws(branch[[1]]) else NA_character_,
    dirty = length(dirty) > 0L
  )
}


#' Capture package versions currently loaded in the R session
#'
#' @return Named character vector of package versions.
#' @keywords internal
.manifest_package_versions <- function() {
  pkgs <- sort(loadedNamespaces())

  stats::setNames(
    as.list(
      purrr::map_chr(
        pkgs,
        function(pkg) {
          tryCatch(
            as.character(utils::packageVersion(pkg)),
            error = function(e) NA_character_
          )
        }
      )
    ),
    pkgs
  )
}


#' Generate a unique manifest run identifier
#'
#' @return Character scalar.
#' @keywords internal
.manifest_run_id <- function() {
  paste0(
    "run_",
    format(Sys.time(), "%Y%m%dT%H%M%OS3", tz = "UTC"),
    "_pid",
    Sys.getpid()
  ) |>
    gsub("[^A-Za-z0-9_]", "", x = _)
}


#' Start or load a dataset provenance manifest
#'
#' @param manifest_path Character. Path to the JSON manifest.
#' @param dataset_id Character scalar.
#' @param duckdb_path Character scalar.
#' @param base_dir Character scalar.
#' @param selection Optional named list describing the dataset selection.
#' @param hash_files Logical. Calculate SHA-256 for manifest-recorded files.
#'
#' @return A manifest object with `path` and `run_index`.
#' @keywords internal
.manifest_start <- function(
    manifest_path,
    dataset_id,
    duckdb_path,
    base_dir = ".",
    selection = list(),
    hash_files = FALSE
) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for manifest generation.")
  }

  manifest_path <- normalizePath(
    manifest_path,
    mustWork = FALSE
  )

  dir.create(
    dirname(manifest_path),
    recursive = TRUE,
    showWarnings = FALSE
  )

  manifest_id <- tools::file_path_sans_ext(
    basename(manifest_path)
  )

  manifest <- list(
    schema_version = 1L,
    manifest_type = "amR_dataset",
    manifest_id = manifest_id,
    manifest_created_at = as.character(Sys.time()),
    manifest_updated_at = as.character(Sys.time()),
    dataset_id = dataset_id,
    dataset = list(
      duckdb = duckdb_path,
      selection = selection
    ),
    artifacts = list(),
    runs = list()
  )

  run <- list(
    run_id = .manifest_run_id(),
    status = "running",
    started_at = as.character(Sys.time()),
    finished_at = NA_character_,
    command = commandArgs(trailingOnly = FALSE),
    working_directory = getwd(),
    host = as.list(Sys.info()),
    r = list(
      version = R.version.string,
      platform = R.version$platform
    ),
    git = .manifest_git_info(base_dir),
    packages = .manifest_package_versions(),
    stages = list(),
    events = list()
  )

  if (is.null(manifest$runs)) {
    manifest$runs <- list()
  }

  manifest$runs[[length(manifest$runs) + 1L]] <- run
  manifest$manifest_updated_at <- as.character(Sys.time())

  run_index <- length(manifest$runs)

  jsonlite::write_json(
    manifest,
    manifest_path,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )

  structure(
    list(
      manifest = manifest,
      path = manifest_path,
      run_index = run_index,
      hash_files = isTRUE(hash_files)
    ),
    class = "amr_manifest"
  )
}


#' Update a manifest stage
#'
#' @param manifest_state Manifest state returned by [.manifest_start()].
#' @param name Character stage name.
#' @param status Character stage status.
#' @param parameters Optional named list.
#' @param inputs Optional character vector of input paths.
#' @param outputs Optional character vector of output paths.
#' @param tool Optional named list describing the tool.
#' @param metrics Optional named list of metrics.
#' @param message Optional log message.
#'
#' @return Updated manifest state.
#' @keywords internal
.manifest_stage <- function(
    manifest_state,
    name,
    status = "success",
    parameters = list(),
    inputs = character(),
    outputs = character(),
    tool = list(),
    metrics = list(),
    message = NULL
) {
  if (!inherits(manifest_state, "amr_manifest")) {
    stop("Invalid manifest state.")
  }

  stage_index <- which(
    purrr::map_lgl(
      manifest_state$manifest$runs[[manifest_state$run_index]]$stages,
      ~ identical(.x$name, name) && identical(.x$status, "running")
    )
  )

  stage <- list(
    name = name,
    status = status,
    started_at = as.character(Sys.time()),
    parameters = parameters,
    inputs = .manifest_file_info(inputs, hash = manifest_state$hash_files),
    outputs = .manifest_file_info(outputs, hash = manifest_state$hash_files),
    tool = tool,
    metrics = metrics
  )

  if (!is.null(message)) {
    stage$message <- as.character(message)
  }

  if (length(stage_index) == 1L) {
    existing <- manifest_state$manifest$runs[[manifest_state$run_index]]$stages[[stage_index]]

    stage$started_at <- existing$started_at
    stage$finished_at <- if (status != "running") {
      as.character(Sys.time())
    } else {
      NULL
    }

    manifest_state$manifest$runs[[manifest_state$run_index]]$stages[[stage_index]] <- stage
  } else {
    if (status != "running") {
      stage$finished_at <- as.character(Sys.time())
    }

    manifest_state$manifest$runs[[manifest_state$run_index]]$stages <-
      append(
        manifest_state$manifest$runs[[manifest_state$run_index]]$stages,
        list(stage)
      )
  }

  manifest_state$manifest$manifest_updated_at <- as.character(Sys.time())

  jsonlite::write_json(
    manifest_state$manifest,
    manifest_state$path,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )

  manifest_state
}


#' Validate an amR dataset manifest
#'
#' @param manifest The previously parsed manifest list.
#' @return `TRUE` invisibly if valid, else throws error.
#' @keywords internal
.manifest_validate <- function(manifest) {
  if (!is.list(manifest)) {
    stop("Manifest must be a list.")
  }

  schema_version <- if (is.null(manifest$schema_version)) {
    "missing"
  } else {
    manifest$schema_version
  }

  if (
    is.null(manifest$schema_version) ||
    !identical(as.integer(manifest$schema_version), 1L)
  ) {
    stop(
      "Unsupported amR manifest schema version: ",
      schema_version,
      ". Expected schema version 1."
    )
  }

  if (!identical(manifest$manifest_type, "amR_dataset")) {
    stop("Manifest is not an amR dataset manifest.")
  }

  required <- c(
    "manifest_id",
    "dataset_id",
    "dataset",
    "artifacts",
    "runs"
  )

  missing <- setdiff(required, names(manifest))

  if (length(missing)) {
    stop(
      "Manifest is missing required field(s): ",
      paste(missing, collapse = ", ")
    )
  }

  invisible(TRUE)
}


#' Append a provenance event to the active manifest run
#'
#' @param manifest_state Manifest state returned by [.manifest_start()].
#' @param level Character event level.
#' @param message Character message.
#' @param details Optional named list.
#'
#' @return Updated manifest state.
#' @keywords internal
.manifest_event <- function(
    manifest_state,
    level = "info",
    message,
    details = list()
) {
  manifest_state$manifest$runs[[manifest_state$run_index]]$events <-
    append(
      manifest_state$manifest$runs[[manifest_state$run_index]]$events,
      list(
        list(
          timestamp = as.character(Sys.time()),
          level = level,
          message = message,
          details = details
        )
      )
    )

  manifest_state$manifest$manifest_updated_at <- as.character(Sys.time())

  jsonlite::write_json(
    manifest_state$manifest,
    manifest_state$path,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )

  manifest_state
}

#' Finish an active provenance manifest run
#'
#' @param manifest_state Manifest state returned by [.manifest_start()].
#' @param status Final run status.
#' @param error Optional error message.
#'
#' @return Invisibly returns the final manifest state.
#' @keywords internal
.manifest_finish <- function(
    manifest_state,
    status = "success",
    error = NULL
) {
  manifest_state$manifest$runs[[manifest_state$run_index]]$status <- status
  manifest_state$manifest$runs[[manifest_state$run_index]]$finished_at <-
    as.character(Sys.time())

  # Patching to resolve an indefinite `running` failure state in the manifest
  if (identical(status, "failed")) {
    stages <- manifest_state$manifest$runs[[manifest_state$run_index]]$stages
    running_stage <- which(purrr::map_lgl(stages, ~ identical(.x$status, "running")))

    if (length(running_stage)) {
      stage_error <- if (!is.null(error)) {
        as.character(error)
      } else {
        "Parent run failed before this stage completed."
      }

      for (i in running_stage) {
        stages[[i]]$status <- "failed"
        stages[[i]]$finished_at <- as.character(Sys.time())
        stages[[i]]$error <- stage_error
      }

      manifest_state$manifest$runs[[manifest_state$run_index]]$stages <- stages
    }
  }

  if (!is.null(error)) {
    manifest_state$manifest$runs[[manifest_state$run_index]]$error <- as.character(error)
  }

  manifest_state$manifest$manifest_updated_at <- as.character(Sys.time())

  jsonlite::write_json(
    manifest_state$manifest,
    manifest_state$path,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )

  invisible(manifest_state)
}

#' Resume provenance logging in an existing manifest
#'
#' Loads an existing manifest and appends a new run.
#'
#' @param manifest_path Character. Path to an existing JSON manifest.
#' @param base_dir Character. Project root.
#' @param hash_files Logical. Calculate SHA-256 checksums for manifest-recorded files.
#'
#' @return A manifest object with `path` and `run_index`.
#' @keywords internal
.manifest_resume <- function(
    manifest_path,
    base_dir = ".",
    hash_files = FALSE
) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for manifest generation.")
  }

  manifest_path <- normalizePath(
    manifest_path,
    mustWork = TRUE
  )

  manifest <- jsonlite::read_json(
    manifest_path,
    simplifyVector = FALSE
  )

  .manifest_validate(manifest)

  if (is.null(manifest$runs)) {
    manifest$runs <- list()
  }

  run <- list(
    run_id = .manifest_run_id(),
    status = "running",
    started_at = as.character(Sys.time()),
    finished_at = NA_character_,
    command = commandArgs(trailingOnly = FALSE),
    working_directory = getwd(),
    host = as.list(Sys.info()),
    r = list(
      version = R.version.string,
      platform = R.version$platform
    ),
    git = .manifest_git_info(base_dir),
    packages = .manifest_package_versions(),
    stages = list(),
    events = list()
  )

  manifest$runs[[length(manifest$runs) + 1L]] <- run
  manifest$manifest_updated_at <- as.character(Sys.time())

  run_index <- length(manifest$runs)

  jsonlite::write_json(
    manifest,
    manifest_path,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )

  structure(
    list(
      manifest = manifest,
      path = manifest_path,
      run_index = run_index,
      hash_files = isTRUE(hash_files)
    ),
    class = "amr_manifest"
  )
}

#' Find the newest amRml-ready manifest associated with a parquet directory
#'
#' You could have multiple runs and manifests in a single data directory, so find
#' the most recently completed manifest.
#'
#' @return Path to the most recent successful manifest for a dataset, or yells `NULL`
#' @keywords internal
.manifest_find_latest_ml <- function(parquet_dir) {
  parquet_dir <- normalizePath(
    parquet_dir,
    mustWork = TRUE
  )

  manifests <- list.files(
    parquet_dir,
    pattern = "^manifest_.*\\.json$",
    full.names = TRUE
  )

  if (!length(manifests)) {
    return(NULL)
  }

  manifests <- manifests[
    order(
      file.info(manifests)$mtime,
      decreasing = TRUE
    )
  ]

  for (manifest_path in manifests) {
    manifest <- tryCatch(
      jsonlite::read_json(
        manifest_path,
        simplifyVector = FALSE
      ),
      error = function(e) NULL
    )

    if (is.null(manifest)) {
      next
    }

    ml_input <- tryCatch(
      .manifest_ml_input(manifest),
      error = function(e) NULL
    )

    if (is.null(ml_input)) {
      next
    }

    artifact_dir <- tryCatch(
      normalizePath(
        ml_input$artifact$directory,
        mustWork = TRUE
      ),
      error = function(e) NA_character_
    )

    if (
      length(artifact_dir) == 1L &&
      !is.na(artifact_dir) &&
      identical(artifact_dir, parquet_dir)
    ) {
      return(manifest_path)
    }
  }

  NULL
}

#' Validate and extract an amRml-ready manifest artifact
#'
#' @param manifest Parsed amR dataset manifest
#' @return The ready artifact and run it came from, or `NULL` if unusable
#' @keywords internal
.manifest_ml_input <- function(manifest) {
  .manifest_validate(manifest)

  artifact <- manifest$artifacts$amRml_input

  if (
    is.null(artifact) ||
    !identical(artifact$status, "ready")
  ) {
    return(NULL)
  }

  producer_run_id <- artifact$producer_run_id

  if (
    is.null(producer_run_id) ||
    !length(producer_run_id) ||
    !nzchar(producer_run_id)
  ) {
    return(NULL)
  }

  producer_runs <- purrr::keep(
    manifest$runs,
    ~ identical(.x$run_id, producer_run_id)
  )

  if (
    length(producer_runs) != 1L ||
    !identical(producer_runs[[1]]$status, "success")
  ) {
    return(NULL)
  }

  parquet_dir <- artifact$directory
  parquet_duckdb <- artifact$parquet_duckdb
  metadata_parquet <- artifact$metadata_parquet

  if (
    is.null(parquet_dir) ||
    !dir.exists(parquet_dir) ||
    is.null(parquet_duckdb) ||
    !file.exists(parquet_duckdb) ||
    is.null(metadata_parquet) ||
    !file.exists(metadata_parquet)
  ) {
    return(NULL)
  }

  list(
    artifact = artifact,
    producer_run = producer_runs[[1]]
  )
}
