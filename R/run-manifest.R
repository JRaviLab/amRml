# Paths within an amRml run ORB
.runPaths <- function(out_dir) {
  list(
    manifest = file.path(out_dir, "manifest.json"),
    metadata = file.path(out_dir, "metadata"),
    eligibility = file.path(out_dir, "metadata", "eligibility.parquet"),
    matrices = file.path(out_dir, "metadata", "matrices.parquet"),
    state = file.path(out_dir, "state"),
    features = file.path(out_dir, "state", "features"),
    matrix_index = file.path(out_dir, "state", "matrix_index.rds"),
    workflow = file.path(out_dir, "state", "workflow.rds"),
    checkpoints = file.path(out_dir, "state", "fit-checkpoints"),
    completion = file.path(out_dir, "state", "fit-completion.rds"),
    results = file.path(out_dir, "results"),
    tasks = file.path(out_dir, "results", "tasks.parquet"),
    performance = file.path(out_dir, "results", "performance.parquet"),
    importance = file.path(out_dir, "results", "importance.parquet"),
    predictions = file.path(out_dir, "results", "predictions.parquet"),
    magic_ready_performance = file.path(out_dir, "results", "magic_ready_performance.parquet"),
    magic_ready_importance = file.path(out_dir, "results", "magic_ready_importance.parquet"),
    rescoring = file.path(out_dir, "results", "rescoring"),
    notes = file.path(out_dir, "results", "notes.parquet"),
    models = file.path(out_dir, "results", "models"),
    merged = file.path(out_dir, "results", "merged")
  )
}

.writeParquet <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  arrow::write_parquet(x, path, compression = "zstd")
  invisible(path)
}

.manifestRunId <- function() {
  gsub("[^A-Za-z0-9_]", "", paste0(
    "run_", format(Sys.time(), "%Y%m%dT%H%M%OS3", tz = "UTC"), "_pid", Sys.getpid()
  ))
}

.manifestGitInfo <- function(base_dir = ".") {
  git <- Sys.which("git")
  if (!nzchar(git)) return(list(available = FALSE))
  run <- function(args) tryCatch(
    suppressWarnings(system2(
      git, c("-C", shQuote(base_dir), args), stdout = TRUE, stderr = FALSE
    )),
    error = function(e) character()
  )
  inside <- run(c("rev-parse", "--is-inside-work-tree"))
  if (!length(inside) || trimws(inside[[1]]) != "true") {
    return(list(available = TRUE, repository = FALSE))
  }
  commit <- run(c("rev-parse", "HEAD"))
  branch <- run(c("rev-parse", "--abbrev-ref", "HEAD"))
  list(
    available = TRUE, repository = TRUE,
    commit = if (length(commit)) trimws(commit[[1]]) else NA_character_,
    branch = if (length(branch)) trimws(branch[[1]]) else NA_character_,
    dirty = length(run(c("status", "--porcelain"))) > 0L
  )
}

.manifestPackageVersions <- function() {
  packages <- sort(loadedNamespaces())
  stats::setNames(lapply(packages, function(package) {
    tryCatch(as.character(utils::packageVersion(package)), error = function(e) NA_character_)
  }), packages)
}

.manifestWrite <- function(manifest, path) {
  manifest$manifest_updated_at <- as.character(Sys.time())
  jsonlite::write_json(manifest, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  invisible(manifest)
}

.manifestFileInfo <- function(paths, out_dir, hash = FALSE) {
  paths <- unique(as.character(paths))
  paths <- paths[!is.na(paths) & nzchar(paths)]
  lapply(paths, function(path) {
    exists <- file.exists(path)
    relative <- if (startsWith(path, paste0(out_dir, "/"))) {
      substring(path, nchar(out_dir) + 2L)
    } else path
    record <- list(
      path = relative,
      exists = exists,
      size_bytes = if (exists) unname(file.info(path)$size) else NA_real_,
      modified_at = if (exists) as.character(file.info(path)$mtime) else NA_character_
    )
    if (hash && exists && !dir.exists(path)) record$md5 <- unname(tools::md5sum(path))
    record
  })
}

.manifestNewRun <- function(parameters) {
  list(
    run_id = .manifestRunId(), status = "running",
    started_at = as.character(Sys.time()), finished_at = NA_character_,
    command = commandArgs(trailingOnly = FALSE), working_directory = getwd(),
    host = as.list(Sys.info()),
    r = list(version = R.version.string, platform = R.version$platform),
    git = .manifestGitInfo(), packages = .manifestPackageVersions(),
    parameters = parameters, stages = list(), events = list()
  )
}

.manifestStart <- function(out_dir, orb, parameters, resume = FALSE) {
  path <- .runPaths(out_dir)$manifest
  if (file.exists(path)) {
    manifest <- jsonlite::read_json(path, simplifyVector = FALSE)
    if (!identical(manifest$manifest_type, "amR_ml") ||
        !identical(as.integer(manifest$schema_version), 1L)) {
      stop("The saved amRml manifest is invalid or unsupported.", call. = FALSE)
    }
  } else {
    if (resume) stop("No amRml manifest was found; this run cannot be resumed.", call. = FALSE)
    created <- as.character(Sys.time())
    manifest <- list(
      schema_version = 1L, manifest_type = "amR_ml",
      manifest_id = paste0("manifest_", .manifestRunId()),
      manifest_created_at = created, manifest_updated_at = created,
      dataset_id = orb$dataset_id,
      dataset = list(
        label = orb$dataset_label,
        input_orb = list(
          directory = orb$directory,
          manifest = orb$manifest_path,
          producer_run_id = orb$producer_run_id,
          finished_at = orb$finished_at
        )
      ),
      artifacts = list(), runs = list()
    )
  }
  manifest$runs[[length(manifest$runs) + 1L]] <- .manifestNewRun(parameters)
  .manifestWrite(manifest, path)
  path
}

.manifestCurrentRunId <- function(path) {
  manifest <- jsonlite::read_json(path, simplifyVector = FALSE)
  manifest$runs[[length(manifest$runs)]]$run_id
}

.manifestStage <- function(path, name, status = "success", parameters = list(),
                           inputs = character(), outputs = character(), metrics = list(),
                           message = NULL, hash = FALSE) {
  manifest <- jsonlite::read_json(path, simplifyVector = FALSE)
  i <- length(manifest$runs)
  stages <- manifest$runs[[i]]$stages
  running <- which(vapply(stages, function(x) {
    identical(x$name, name) && identical(x$status, "running")
  }, logical(1)))
  stage <- list(
    name = name, status = status, started_at = as.character(Sys.time()),
    parameters = parameters,
    inputs = .manifestFileInfo(inputs, dirname(path), hash),
    outputs = .manifestFileInfo(outputs, dirname(path), hash),
    metrics = metrics
  )
  if (!is.null(message)) stage$message <- as.character(message)
  if (length(running) == 1L) {
    stage$started_at <- stages[[running]]$started_at
    if (status != "running") stage$finished_at <- as.character(Sys.time())
    stages[[running]] <- stage
  } else {
    if (status != "running") stage$finished_at <- as.character(Sys.time())
    stages[[length(stages) + 1L]] <- stage
  }
  manifest$runs[[i]]$stages <- stages
  .manifestWrite(manifest, path)
  invisible(path)
}

.manifestEvent <- function(path, message, level = "info", details = list()) {
  manifest <- jsonlite::read_json(path, simplifyVector = FALSE)
  i <- length(manifest$runs)
  manifest$runs[[i]]$events[[length(manifest$runs[[i]]$events) + 1L]] <- list(
    timestamp = as.character(Sys.time()), level = level,
    message = as.character(message), details = details
  )
  .manifestWrite(manifest, path)
  invisible(path)
}

.manifestArtifact <- function(path, name, status = "ready", details = list()) {
  manifest <- jsonlite::read_json(path, simplifyVector = FALSE)
  manifest$artifacts[[name]] <- c(list(
    status = status, updated_at = as.character(Sys.time())
  ), details)
  .manifestWrite(manifest, path)
  invisible(path)
}

.manifestFinish <- function(path, status = "success", error = NULL) {
  if (!file.exists(path)) return(invisible(path))
  manifest <- jsonlite::read_json(path, simplifyVector = FALSE)
  i <- length(manifest$runs)
  manifest$runs[[i]]$status <- status
  manifest$runs[[i]]$finished_at <- as.character(Sys.time())
  if (!is.null(error)) manifest$runs[[i]]$error <- as.character(error)
  if (status == "failed") {
    for (j in seq_along(manifest$runs[[i]]$stages)) {
      if (identical(manifest$runs[[i]]$stages[[j]]$status, "running")) {
        manifest$runs[[i]]$stages[[j]]$status <- "failed"
        manifest$runs[[i]]$stages[[j]]$finished_at <- as.character(Sys.time())
        manifest$runs[[i]]$stages[[j]]$error <- error %||% "Parent run failed."
      }
    }
  }
  .manifestWrite(manifest, path)
  invisible(path)
}

.registerManifest <- function(path) {
  manifest <- jsonlite::read_json(path, simplifyVector = FALSE)
  rname <- paste("amR_ml_manifest", manifest$dataset_id, manifest$manifest_id, sep = "_")
  bfc <- BiocFileCache::BiocFileCache(ask = FALSE)
  hit <- BiocFileCache::bfcquery(bfc, query = rname, field = "rname", exact = TRUE)
  if (!nrow(hit)) {
    BiocFileCache::bfcadd(
      bfc, rname = rname, fpath = normalizePath(path), rtype = "local",
      action = "asis", progress = FALSE
    )
  } else if (!identical(normalizePath(hit$rpath[[1]], mustWork = FALSE), normalizePath(path))) {
    BiocFileCache::bfcupdate(bfc, hit$rid[[1]], rpath = normalizePath(path), rname = rname)
  }
  invisible(rname)
}
