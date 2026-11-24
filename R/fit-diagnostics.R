# Diagnostics inhabit a realm beyond state so they don't get wiped out when a run
# is cleaned up after finishing
.fitOwner <- function() {
  list(host = unname(Sys.info()[["nodename"]]), pid = Sys.getpid(),
       created = format(Sys.time(), tz = "UTC", usetz = TRUE),
       job_id = Sys.getenv("SLURM_JOB_ID"),
       array_job_id = Sys.getenv("SLURM_ARRAY_JOB_ID"),
       array_task_id = Sys.getenv("SLURM_ARRAY_TASK_ID"))
}

# When models crash (from OOM, for example) they leave behind a lock file, which
# needs to be dealt with before resuming and otherwise keeps things safe
.acquireFitLock <- function(path, label) {
  if (!dir.create(path, showWarnings = FALSE)) {
    owner_path <- file.path(path, "owner.json")
    owner <- if (file.exists(owner_path)) paste(readLines(owner_path, warn = FALSE), collapse = "") else
      "ownership unavailable (legacy or interrupted lock creation)"
    stop(label, " is already running or has a stale lock: ", path,
         "\nLock owner: ", owner,
         "\nConfirm the owning job has exited before removing this lock.", call. = FALSE)
  }
  tryCatch(jsonlite::write_json(.fitOwner(), file.path(path, "owner.json"), auto_unbox = TRUE),
           error = function(e) {
             unlink(path, recursive = TRUE)
             stop(e)
           })
  invisible(path)
}

.isMemoryError <- function(message) {
  grepl("cannot allocate|out[ -]of[ -]memory|std::bad_alloc|vector memory (exhausted|limit.*reached)|memory exhausted",
        message, ignore.case = TRUE)
}

.fitDiagnostic <- function(directory, event, details = list()) {
  if (is.null(directory)) return(invisible(NULL))
  # One stream per process avoids unrelated workers spamming things
  tryCatch({
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    owner <- .fitOwner()
    host <- gsub("[^A-Za-z0-9_.-]", "_", owner$host)
    path <- file.path(directory, paste0("worker-", host, "-", owner$pid, ".jsonl"))
    record <- c(list(event = event), owner, details)
    cat(jsonlite::toJSON(record, auto_unbox = TRUE, null = "null"), "\n", file = path, append = TRUE)
  }, error = function(e) message("Could not write fit diagnostics: ", conditionMessage(e)))
  invisible(NULL)
}
