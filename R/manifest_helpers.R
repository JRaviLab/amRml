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

  manifest <- list(
    schema_version = 1L,
    manifest_created_at = as.character(Sys.time()),
    manifest_updated_at = as.character(Sys.time()),
    dataset_id = dataset_id,
    dataset = list(
      duckdb = duckdb_path,
      selection = selection
    ),
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

# To distinguish multiple manifests in the same bug directory
.manifest_find_latest <- function(duckdb_path) {
  manifest_dir <- dirname(normalizePath(
    duckdb_path,
    mustWork = FALSE
  ))

  manifests <- list.files(
    manifest_dir,
    pattern = "^manifest_.*\\.json$",
    full.names = TRUE
  )

  if (!length(manifests)) {
    return(NULL)
  }

  manifests[which.max(file.info(manifests)$mtime)]
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
