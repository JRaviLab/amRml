#' Prepare, launch and collect a modeling workflow on an HPC cluster
#'
#' Uses the same models, eligibility rules and checkpoints as
#' [runModelingWorkflow()]. Currently supports SLURM clusters.
#'
#' @param stage `"prep"` prepares matrices and chunks; `"run"` writes a runner and
#'   returns its submission command; `"finish"` collects model results.
#' @param out_dir Run directory shared by the submission host and compute nodes.
#' @param ... Arguments to [runModelingWorkflow()] during `"prep"`, including
#'   `orb_dir`, `feature_types`, `model_types`, `seeds`, and `n_chunks`.
#'   Do not supply `execution` or `resume`. Other stages accept no extra arguments.
#' @param scheduler Submission backend, currently `"slurm"`.
#' @param cpus_per_task CPUs and parallel R workers per array element. Default `1`.
#'   Modeling logic does not currently benefit from multiple CPUs per model!
#' @param max_concurrent Maximum simultaneously running array elements.
#' @param memory_mb For `"run"`, a required, strictly increasing vector of total
#'   per-node memory allocations in SLURM M units. The first is the initial
#'   allocation per task. Subsequent values are for retries of scheduler-confirmed 
#'   OOM kills. i.e., If a task fails due to OOM errors, it will be retried
#'   automatically with the next memory value in the vector, and so on.
#'   A single value disables retries. This does not apply to prep or finish.
#' @param scheduler_options Named list of SLURM `account`, `partition`, `qos`,
#'   `time`, `constraint`, or `reservation` options. No cluster-specific defaults.
#'   Check with your HPC cluster for which of these are required for job submission.
#' @param preamble Vector of trusted Bash commands executed before R, such as for
#'   module loading and environment activation. Compute nodes must have `Rscript` 
#'   and the same installed amRml package and R libraries available. This can be
#'   commonly done using Anaconda environments. Environment configuration is up 
#'   to individual users depending on their cluster and its limitations!
#'
#' @details
#' None of these stages submit jobs. Execute the returned `commands` on your cluster
#' to submit the array and its `afterany` monitor. The generated `hpc/run.sh`
#' handles submission, workers and monitors. Resources and job IDs are saved in
#' `hpc/submission.rds`; logs are written under `hpc/` and `diagnostics/`.
#' Repeating the submission command is refused once a submission was attempted.
#' Please use a new run directory to change the plan or generated resource settings.
#'
#' The monitor requires `sbatch`, `sacct` and `squeue` on compute nodes. It checks
#' that every submitted job has terminated and removes only locks owned by the 
#' terminated OOM attempt. Unknown locks are retained for user inspection. 
#' Other failures and exit 137 alone are not retried. R allocation errors saved 
#' as `fit_oom` are task failures and are not refitted. We recommend that you 
#' inspect collected task status even when scheduler jobs completed successfully.
#'
#' Monitor jobs use 1 CPU, 512MB RAM and 10 minutes, inherit all scheduler options, 
#' and override the worker walltime. After all waves and monitors end, run `"finish"` 
#' on a compute node for result collection and downstream analyses. A failed monitor 
#' or ambiguous submission requires inspection of saved state and scheduler logs 
#' before recovery. Please (PLEASE) do not modify files during a run.
#'
#' @return For prep and run, a list with `out_dir`, `n_chunks`, `n_models`, and a
#'   `models` summary by mode, target and feature type. Prep also returns the
#'   eligibility scope table. Run adds `runner`, `state_file`, `commands`, and
#'   resource settings. Finish returns an invisible `amr_fit`, as [collectTaskArray()].
#' @examples
#' \dontrun{
#' runModelingWorkflowHPC("prep", out_dir = "results/run", orb_dir = "data/orb",
#'                        feature_types = "standard", n_chunks = 100L)
#' job <- runModelingWorkflowHPC("run", out_dir = "results/run",
#'   memory_mb = c(2000L, 8000L, 16000L), max_concurrent = 50L,
#'   scheduler_options = list(account = "my_account", partition = "my_partition"))
#' job$commands # Review, then execute on the cluster.
#' fit <- runModelingWorkflowHPC("finish", out_dir = "results/run")
#' }
#' @export
runModelingWorkflowHPC <- function(stage = c("prep", "run", "finish"), out_dir,
                                   ..., scheduler = "slurm", cpus_per_task = 1L,
                                   max_concurrent = 50L, memory_mb = NULL,
                                   scheduler_options = list(), preamble = character()) {
  stage <- match.arg(stage)
  dots <- list(...)
  if (stage != "prep" && length(dots)) stop("Extra workflow arguments are accepted only for stage = 'prep'.")
  if (stage == "prep") {
    if (any(c("execution", "resume") %in% names(dots))) stop("HPC preparation sets execution and resume.")
    plan <- do.call(runModelingWorkflow, c(list(out_dir = out_dir, execution = "prepare"), dots))
    return(c(.hpcSummary(plan), list(eligibility = plan$tasks$built$eligibility$scopes)))
  }
  out_dir <- normalizePath(out_dir, mustWork = TRUE)
  path <- file.path(out_dir, "hpc", "submission.rds")
  if (stage == "finish") {
    if (file.exists(path)) {
      waves <- readRDS(path)$waves
      if (dir.exists(paste0(path, ".watch.lock")) || !length(waves) ||
          !identical(waves[[length(waves)]]$status, "accounted")) {
        stop("HPC monitoring is unfinished or has a stale guard; inspect submission state and scheduler logs.")
      }
    }
    return(invisible(collectTaskArray(out_dir)))
  }
  if (!identical(scheduler, "slurm")) stop("The supported submission scheduler is 'slurm'.")
  .checkArgCount(cpus_per_task, "cpus_per_task")
  .checkArgCount(max_concurrent, "max_concurrent")
  if (!is.numeric(memory_mb) || !length(memory_mb) || anyNA(memory_mb) ||
      any(!is.finite(memory_mb) | memory_mb <= 0 | memory_mb > .Machine$integer.max |
          memory_mb != floor(memory_mb)) || any(diff(memory_mb) <= 0)) {
    stop("memory_mb must be a strictly increasing vector of positive integer M allocations.")
  }
  options <- scheduler_options
  if (!is.list(options) || (length(options) && (is.null(names(options)) ||
      anyDuplicated(names(options)) || any(!names(options) %in%
        c("account", "partition", "qos", "time", "constraint", "reservation")))) ||
      any(!vapply(options, function(x) is.character(x) && length(x) == 1L && !is.na(x) &&
        nzchar(x) && !grepl("[\r\n]", x), logical(1)))) stop("Invalid named scheduler_options.")
  if (!is.character(preamble) || anyNA(preamble)) stop("preamble must be Bash command lines.")
  plan <- .readTaskArray(out_dir)
  if (!length(plan$chunks)) stop("The prepared array has no chunks; finish without submitting.")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  guard <- .slurmGuard(path)
  on.exit(unlink(guard, recursive = TRUE), add = TRUE)
  if (file.exists(path)) stop("A submission is already prepared at ", path, "; use its saved runner command.")
  runner <- file.path(out_dir, "hpc", "run.sh")
  writeLines(c("#!/bin/bash", "set -euo pipefail", preamble,
    "export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1",
    "export VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1",
    paste("exec Rscript -e", shQuote("amRml:::.runHPC(commandArgs(TRUE))"),
          shQuote(path), '"$@"')), runner)
  Sys.chmod(runner, "0755")
  state <- list(run_dir = out_dir, runner = runner, memory_mb = as.integer(memory_mb),
                cpus_per_task = cpus_per_task, array_limit = max_concurrent,
                worker_args = if (length(options)) paste0("--", names(options), "=", unlist(options)) else character(),
                chunks = seq_along(plan$chunks), waves = list())
  .saveFitCheckpoint(state, path)
  c(.hpcSummary(plan), list(runner = runner, state_file = path,
    commands = paste("bash", shQuote(runner), "submit"),
    memory_mb = state$memory_mb, cpus_per_task = cpus_per_task,
    max_concurrent = max_concurrent, scheduler_options = options))
}

.hpcSummary <- function(plan) {
  tasks <- plan$tasks$tasks
  matrices <- plan$tasks$built$matrices
  models <- matrices[match(tasks$matrix_id, matrices$matrix_id),
                      c("mode_id", "unit", "target", "feature_type"), drop = FALSE]
  if (nrow(models)) models <- stats::aggregate(list(n_models = rep(1L, nrow(models))), models, sum)
  list(out_dir = plan$tasks$built$out_dir, n_chunks = length(plan$chunks),
       n_models = nrow(tasks), models = models)
}

# Invoked by the single generated runner, no included cluster configuration 
.runHPC <- function(args) {
  if (length(args) < 2L || !args[[2L]] %in% c("submit", "worker", "watch") ||
      length(args) != (if (args[[2L]] == "watch") 3L else 2L)) stop("Invalid HPC runner arguments.")
  path <- args[[1L]]
  if (args[[2L]] == "watch") {
    attempt <- suppressWarnings(as.numeric(args[[3L]]))
    .checkArgCount(attempt, "attempt", min = 0)
    return(.slurmWatch(path, attempt))
  }
  state <- readRDS(path)
  if (args[[2L]] == "worker") return(runTaskChunk(state$run_dir, threads = state$cpus_per_task))
  guard <- .slurmGuard(path)
  on.exit(unlink(guard, recursive = TRUE), add = TRUE)
  state <- readRDS(path)
  if (length(state$waves)) stop("Submission was already attempted; inspect saved state before recovery.")
  .slurmSubmit(state, path, state$chunks, 0L)
}
