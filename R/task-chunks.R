# Split without cut(): empty and singleton plans are valid.
.splitGroupIds <- function(n, n_chunks) {
  if (!n) return(list())
  k <- min(n, n_chunks)
  unname(split(seq_len(n), (seq_len(n) - 1) %% k + 1))
}

#' Split fit groups into chunks
#'
#' Keeps tasks that share the same fitted model together, how beautiful. Return 
#' contains group IDs. For array execution, use the chunk map saved by 
#' `fitModels(execution = "prepare")` or the workflow preparation phase.
#' @param tasks An `amr_tasks` object from [expandTasks()].
#' @param n_chunks Maximum number of chunks. Fewer chunks assign more fit groups 
#' to each array element. Groups sharing a fit stay together.
#' @param reuse_fits Respect fit-group reuse. Must match fit settings.
#' @return A list of integer group index vectors. Empty plan returns `list()`.
#' @export
splitTaskGroups <- function(tasks, n_chunks = 100L, reuse_fits = TRUE) {
  .checkArgClass(tasks, "tasks", "amr_tasks", "expandTasks()")
  .checkArgCount(n_chunks, "n_chunks")
  .checkArgFlag(reuse_fits, "reuse_fits")
  index <- readRDS(.runPaths(tasks$built$out_dir)$matrix_index)
  groups <- if (reuse_fits) .fitTaskGroups(tasks$tasks, tasks$built$matrices, index) else
    as.list(seq_len(nrow(tasks$tasks)))
  .splitGroupIds(length(groups), n_chunks)
}

#' Get one chunk from a chunk map
#' @param chunks Output from [splitTaskGroups()].
#' @param chunk_id One-based array task ID. A numeric string is accepted.
#' @return Integer vector of group IDs.
#' @export
taskChunk <- function(chunks, chunk_id) {
  if (is.character(chunk_id) && length(chunk_id) == 1L && !is.na(chunk_id) &&
      grepl("^[0-9]+$", chunk_id)) chunk_id <- as.numeric(chunk_id)
  .checkArgCount(chunk_id, "chunk_id")
  if (!is.list(chunks) || chunk_id > length(chunks)) {
    .amrAbort("invalid_argument", "`chunk_id` is outside the saved chunk map.")
  }
  chunks[[chunk_id]]
}

.readTaskArray <- function(out_dir) {
  path <- file.path(.runPaths(out_dir)$checkpoints, "array.rds")
  if (!file.exists(path)) stop("No prepared array plan found in: ", out_dir, call. = FALSE)
  readRDS(path)
}

#' Run one prepared array chunk
#'
#' Runs a saved chunk from `fitModels(execution = "prepare")` or
#' [runModelingWorkflowHPC()]. Completed fit groups are skipped on retries.
#' Workers write only their own checkpoints and models; shared results are
#' written during collection. Do not change the plan or inputs while jobs run.
#'
#' Directory locks stop workers from simultaneously fiddling with a chunk. A lock's
#' `owner.json` records hostname, PID, time and SLURM IDs. Confirm a job has exited 
#' before removing a stale lock. Group diagnostics are retained in
#' `diagnostics/worker-*.jsonl`, including task/matrix IDs and R heap use.
#' R memory errors are recorded as `fit_oom`.
#'
#' @seealso [runModelingWorkflowHPC()] for preparation, generated SLURM runners,
#'   bounded OOM error retries and collection.
#' @param out_dir Prepared run directory on a filesystem shared by all jobs.
#' @param chunk_id One-based chunk ID; defaults to `SLURM_ARRAY_TASK_ID`.
#' @param threads Number of R workers in this array element. Default `1` keeps
#'   the full job memory allocation available to one fit. Must not exceed
#'   `SLURM_CPUS_PER_TASK` when set. Limit BLAS/OpenMP threads in the batch script.
#'   Currently, modeling logic is single-threaded and will not benefit from extra
#'   CPU resources.
#' @return Invisibly, a list of the completed group IDs.
#' @export
runTaskChunk <- function(out_dir, chunk_id = Sys.getenv("SLURM_ARRAY_TASK_ID"), threads = 1L) {
  .checkArgCount(threads, "threads")
  allocated <- Sys.getenv("SLURM_CPUS_PER_TASK")
  if (nzchar(allocated)) {
    cpus <- suppressWarnings(as.numeric(allocated))
    .checkArgCount(cpus, "SLURM_CPUS_PER_TASK")
    if (threads > cpus) .amrAbort("invalid_argument", "`threads` exceeds SLURM_CPUS_PER_TASK.")
  }
  message("Array worker: loading saved task plan.")
  plan <- .readTaskArray(out_dir)
  groups <- taskChunk(plan$chunks, chunk_id)
  id <- match(TRUE, vapply(plan$chunks, identical, logical(1), groups))
  lock <- file.path(.runPaths(out_dir)$checkpoints, sprintf("chunk-%06d.lock", id))
  .acquireFitLock(lock, "Chunk")
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  message("Array chunk ", id, ": ", length(groups), " fit groups assigned.")
  backend <- if (threads == 1L) BiocParallel::SerialParam() else if (.Platform$OS.type == "windows") {
    BiocParallel::SnowParam(workers = threads)
  } else BiocParallel::MulticoreParam(workers = threads)
  fitModels(plan$tasks, save_fits = plan$save_fits, save_state = plan$save_state,
            reuse_fits = plan$reuse_fits, resume = TRUE, BPPARAM = backend,
            execution = "worker", group_subset = groups)
}

#' Collect a completed task array
#'
#' Requires all group checkpoints. Run after all workers finish. Missing or corrupt
#' checkpoints stop collection and preserve state. Individual model failures are
#' retained in the task table, as in [fitModels()]. Workflow arrays also export
#' merged results and finish the workflow manifest using saved settings.
#' @param out_dir Prepared run directory.
#' @return Invisibly, an `amr_fit`.
#' @export
collectTaskArray <- function(out_dir) {
  snapshot <- .runPaths(out_dir)$workflow
  if (file.exists(snapshot)) {
    settings <- readRDS(snapshot)$settings
    return(invisible(do.call(runModelingWorkflow, c(settings, list(
      out_dir = out_dir, resume = TRUE, execution = "collect"
    )))))
  }
  plan <- .readTaskArray(out_dir)
  invisible(fitModels(
    plan$tasks, save_fits = plan$save_fits, save_state = plan$save_state,
    reuse_fits = plan$reuse_fits, resume = TRUE,
    BPPARAM = BiocParallel::SerialParam(), execution = "collect"
  ))
}
