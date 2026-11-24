# Internal SLURM accounting and OOM error recovery
.slurmPause <- function(seconds) Sys.sleep(seconds)

.slurmCommand <- function(command, args) {
  result <- system2(command, vapply(args, shQuote, character(1)), stdout = TRUE)
  status <- attr(result, "status")
  if (!is.null(status) && status != 0L) stop(command, " exited with status ", status, call. = FALSE)
  result
}

.slurmAudit <- function(state, event, ...) {
  record <- c(list(time = format(Sys.time(), tz = "UTC", usetz = TRUE), event = event), list(...))
  directory <- file.path(state$run_dir, "diagnostics")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  cat(jsonlite::toJSON(record, auto_unbox = TRUE, null = "null"), "\n",
      file = file.path(directory, "slurm-events.jsonl"), append = TRUE)
}

.slurmAccounting <- function(text, job_id, chunks) {
  terminal <- c("COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY",
                "NODE_FAIL", "PREEMPTED", "BOOT_FAIL", "DEADLINE", "REVOKED")
  states <- stats::setNames(rep(NA_character_, length(chunks)), as.character(chunks))
  oom <- integer()
  prefix <- paste0(job_id, "_")
  for (line in text) {
    fields <- strsplit(line, "|", fixed = TRUE)[[1L]]
    if (length(fields) < 2L || !startsWith(trimws(fields[[1L]]), prefix)) next
    suffix <- substring(trimws(fields[[1L]]), nchar(prefix) + 1L)
    if (!grepl("^[0-9]+([.].+)?$", suffix)) next
    chunk <- as.integer(sub("[.].*$", "", suffix))
    if (!chunk %in% chunks) next
    status <- sub("[+ ]+$", "", sub(" .*$", "", trimws(fields[[2L]])))
    if (!grepl(".", suffix, fixed = TRUE)) states[[as.character(chunk)]] <- status
    if (identical(status, "OUT_OF_MEMORY")) oom <- c(oom, chunk)
  }
  if (anyNA(states) || any(!states %in% terminal)) {
    stop("Accounting is incomplete or jobs are still active; no locks removed.", call. = FALSE)
  }
  list(states = states, oom = sort(unique(oom)))
}

.slurmCleanLocks <- function(state, job_id, chunks) {
  root <- file.path(state$run_dir, "state", "fit-checkpoints")
  locks <- list.files(root, pattern = "^(chunk|group)-.*[.]lock$", full.names = TRUE)
  for (lock in locks) {
    owner <- file.path(lock, "owner.json")
    link <- Sys.readlink(lock)
    if (!dir.exists(lock) || (!is.na(link) && nzchar(link)) || !file.exists(owner)) next
    metadata <- tryCatch(jsonlite::read_json(owner), error = function(e) NULL)
    if (!is.list(metadata)) next
    if (length(metadata$array_job_id) != 1L || length(metadata$array_task_id) != 1L ||
        !identical(as.character(metadata$array_job_id), as.character(job_id)) ||
        !as.character(metadata$array_task_id) %in% as.character(chunks)) next
    .slurmAudit(state, "remove_terminated_lock", job_id = job_id, lock = lock, owner = metadata)
    if (unlink(lock, recursive = TRUE) != 0L) stop("Could not remove terminated lock: ", lock)
  }
}

.slurmSubmit <- function(state, path, chunks, attempt) {
  memory <- state$memory_mb[[attempt + 1L]]
  wave <- list(chunks = chunks, memory_mb = memory, attempt = attempt, status = "submitting")
  index <- length(state$waves) + 1L
  state$waves[[index]] <- wave
  .saveFitCheckpoint(state, path) # Ambiguous failures never auto-repeated
  args <- c("--parsable", state$worker_args, "--nodes=1", "--ntasks=1", "--no-requeue",
            paste0("--cpus-per-task=", state$cpus_per_task),
            paste0("--output=", state$run_dir, "/hpc/worker-%A_%a.out"),
            paste0("--error=", state$run_dir, "/hpc/worker-%A_%a.err"),
            paste0("--mem=", memory, "M"),
            paste0("--array=", paste(chunks, collapse = ","), "%", state$array_limit),
            state$runner, "worker")
  reply <- .slurmCommand("sbatch", args)
  if (length(reply) != 1L || !grepl("^[0-9]+(;[A-Za-z0-9_.-]+)?$", reply)) {
    stop("Unexpected sbatch reply; inspect SLURM before resubmitting.")
  }
  wave$job_id <- sub(";.*$", "", reply)
  wave$status <- "submitted"
  state$waves[[index]] <- wave
  .saveFitCheckpoint(state, path)
  .slurmAudit(state, "array_submitted", wave = wave)
  monitor <- .slurmCommand("sbatch", c(
    "--parsable", state$worker_args, "--nodes=1", "--ntasks=1", "--cpus-per-task=1", "--no-requeue",
    "--mem=512M", "--time=00:10:00", paste0("--dependency=afterany:", wave$job_id),
    "--job-name=amRml-oom-monitor", paste0("--output=", state$run_dir, "/hpc/monitor-%j.out"),
    state$runner, "watch", as.character(attempt)
  ))
  if (length(monitor) != 1L || !grepl("^[0-9]+(;[A-Za-z0-9_.-]+)?$", monitor)) {
    stop("Unexpected monitor sbatch reply; array ", wave$job_id,
         " is already submitted. Inspect SLURM before resubmitting.")
  }
  wave$monitor_job_id <- sub(";.*$", "", monitor)
  state$waves[[index]] <- wave
  .saveFitCheckpoint(state, path)
  message("Array ", wave$job_id, "; monitor ", paste(monitor, collapse = " "), "; state: ", path)
  invisible(list(job_id = wave$job_id, monitor_job_id = monitor, state_file = path))
}

# Submission holds the same guard as monitors until monitor job ID saved
.slurmGuard <- function(path) {
  guard <- paste0(path, ".watch.lock")
  for (attempt in seq_len(6L)) {
    if (dir.create(guard, showWarnings = FALSE)) return(guard)
    if (attempt < 6L) .slurmPause(1)
  }
  stop("Monitor guard exists: ", guard)
}

.slurmWatch <- function(path, attempt) {
  guard <- .slurmGuard(path)
  on.exit(unlink(guard, recursive = TRUE), add = TRUE)
  state <- readRDS(path)
  index <- length(state$waves)
  wave <- state$waves[[index]]
  if (wave$attempt != attempt || !identical(wave$status, "submitted")) {
    stop("This attempt has already been handled or may require manual review.")
  }
  for (poll in seq_len(6L)) {
    checked <- tryCatch({
      text <- .slurmCommand("sacct", c(
        "--array", "-n", "-P", "-j", wave$job_id,
        "--format=JobID%100,State%40,ExitCode,ReqMem,MaxRSS"
      ))
      result <- .slurmAccounting(text, wave$job_id, wave$chunks)
      queued <- trimws(.slurmCommand("squeue", c("-h", "-u", Sys.info()[["user"]], "-o", "%F")))
      if (wave$job_id %in% queued) stop("Array still appears in squeue")
      c(result, list(raw = text))
    }, error = identity)
    if (!inherits(checked, "error")) break
    if (poll == 6L) stop(checked)
    .slurmPause(10)
  }
  .slurmAudit(state, "accounting", job_id = wave$job_id, states = as.list(checked$states),
              oom_chunks = checked$oom, raw = checked$raw)
  wave$status <- "accounted"
  wave$states <- checked$states
  wave$oom_chunks <- checked$oom
  wave$failed_chunks <- as.integer(names(checked$states)[checked$states != "COMPLETED"])
  state$waves[[index]] <- wave
  .saveFitCheckpoint(state, path)
  non_oom <- setdiff(wave$failed_chunks, checked$oom)
  if (length(non_oom)) {
    message("Chunks failed without confirmed OOM error: ", paste(non_oom, collapse = ", "),
            ". These are not resubmitted. Please inspect scheduler and worker logs.")
  }
  if (!length(checked$oom)) {
    if (!length(wave$failed_chunks)) {
      message("All submitted chunks completed. Inspect task-level failures during collection.")
    } else message("No OOM error retries; non-OOM failures require inspection.")
    return(invisible(NULL))
  }
  if (attempt + 1L >= length(state$memory_mb)) {
    .slurmAudit(state, "retry_limit", job_id = wave$job_id, oom_chunks = checked$oom)
    message("OOM error retry limit reached; checkpoints and locks have been retained.")
    return(invisible(NULL))
  }
  .slurmCleanLocks(state, wave$job_id, checked$oom)
  .slurmSubmit(state, path, checked$oom, attempt + 1L)
}
