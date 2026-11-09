# The single entry point: every stage, from the ORB to fitted models.

#' Run the modeling workflow
#'
#' Runs every stage on an ORB: [readORB()], [profileORB()], [eligibleScopes()],
#' [planMatrices()], [buildMatrices()], [expandTasks()] and [fitModels()]. Each
#' can also be called on its own. By default every task is fitted on both the real
#' and the shuffled labels; `labels = "real"` skips the shuffled controls.
#'
#' @details
#' The default `execution = "fit"` runs the complete workflow locally using
#' `threads` or `BPPARAM`. SLURM arrays are an opt-in advanced workflow: explicitly
#' prepare a run, invoke [runTaskChunk()] from user-managed batch jobs, then call
#' [collectTaskArray()]. SLURM environment variables do not change the default
#' execution path, and the package does not submit scheduler jobs.
#'
#' Everything goes to `out_dir`: `manifest.json`, planning tables in `metadata/`,
#' and durable model outputs in `results/` (including `models/` when
#' `save_fits = TRUE` and `merged/` from [exportMergedResults()]). Sparse stores
#' and resume files are temporary in `state/`; successful workflows remove them
#' unless `save_state = TRUE`. amRdata's files are only read.
#' Live progress reports completed task counts and elapsed minutes:seconds.
#' Progress and completion counts are also recorded as manifest events.
#'
#' @param orb_dir An amRdata ORB directory. If `NULL`, discover registered ORBs
#'   interactively: use the only usable ORB, or prompt to select among several.
#'   Required in non-interactive sessions.
#' @param out_dir Folder to write to. By default a new run folder,
#'   `<ORB>/amRml/<date-time>/`.
#' @param model_types Mode IDs from [modelingModes()], or `NULL` for all of them.
#' @param threads Number of parallel workers. Default `1` runs serially. Larger
#'   values use multicore workers on Unix and socket workers on Windows.
#' @param feature_types Feature types to model: `"quick"`, `"standard"` (default),
#'   `"full"`, or an explicit character vector of ORB feature types.
#'   Selected before matrices are planned or built; see [planMatrices()].
#' @param seeds Seeds passed to [expandTasks()]. Normal mode defaults to
#'   `c(123L, 456L, 789L)`; easy mode defaults to `123L`.
#' @param labels `"real"`, `"shuffled"` or both (the default), passed to
#'   [expandTasks()]. 
#' @param save_fits Retain fitted model RDS files. Default `FALSE`.
#' @param save_state Retain internal RDS matrix stores and resume state after a
#'   successful workflow. Default `FALSE`; failed or interrupted workflows keep
#'   state so they can be resumed. Both flags must match when resuming.
#' @param min_genomes,n_fold,min_test_genomes,min_test_minority,min_fold_minority,min_analysis_minority Passed to
#'   [eligibleScopes()].
#' @param overwrite Replace an earlier run's matrices and results in `out_dir`;
#'   if `FALSE`, an `out_dir` holding earlier output is refused. The default run
#'   folder is always new.
#' @param resume Resume interrupted model fitting in an existing `out_dir`,
#'   using its saved matrices and completed task groups, or reopen unchanged
#'   completed results. Supply the same modeling arguments and the exact run
#'   folder (e.g. `<ORB>/amRml/<date-time>/`); `BPPARAM` may be changed.
#'   If `orb_dir` is omitted, it is taken from the saved workflow.
#'   `overwrite` must be `FALSE`.
#' @param BPPARAM Optional advanced BiocParallel backend, passed to [fitModels()].
#'   Supply either `threads` or `BPPARAM`. May be changed when resuming.
#' @param execution Advanced options for HPC job arrays. `"fit"` runs the full 
#'   workflow normally (default). `"prepare"` builds shared matrices and an array 
#'   plan without fitting. Run [runTaskChunk()] for each saved chunk, then 
#'   [collectTaskArray()] to finish exports. `"collect"` with `resume = TRUE` 
#'   gathers only completed checkpoints.
#' @param n_chunks Maximum number of nonempty array chunks when preparing.
#'   Fewer chunks assign more fit groups to each worker; see [runTaskChunk()].
#' @param mdr_min_classes Minimum resistant classes for MDR, passed to [profileORB()].
#' @param min_mdr_label_genomes Minimum training genomes per MDR combination label, passed to [eligibleScopes()].
#' @param reuse_fits Passed to [fitModels()]; reuse compatible cross-drug,
#'   cross-year, and cross-country fits.
#'
#' @param easy_mode Run a preview with vanilla models, the `"quick"` feature
#'   preset (available gene, protein, Pfam tables), real labels, and one seed.
#'   Default `FALSE`. Conflicting explicit arguments are rejected. Model fitting
#'   and cross-validation are unchanged; automatic rescoring is skipped.
#' @param rescore Automatically run [rescoreFeatures()] after merged export when
#'   at least three fitted real/shuffled seed pairs exist in a vanilla scope.
#'   Default `TRUE`. Expected empty results are recorded without failing the run.
#'
#' @return Preparation returns a plan with tasks, chunks and fit settings.
#'   Otherwise the `amr_fit` from [fitModels()]; every earlier stage's result can be
#'   reached from it (`$built`, `$built$eligibility`, ...). `$rescoring` contains
#'   the discovery result when available; `$rescoring_status` records the outcome.
#'   Rescoring tables and `discovery.rds` are saved in `results/rescoring/`, with
#'   provenance and skip/empty reasons in the manifest. Binary MDR and individual
#'   resistance-combination labels are included. Cross-target tables compare
#'   shared MDR, drug and class feature/dyad memberships. This uses returned-feature
#'   ranking, not the separate candidate-aware scoring helper or network export.
#'
#' @section Errors:
#' Errors from any stage, including `amrml_invalid_argument` for mismatched resume
#' settings, propagate to the caller.
#'
#' @examples
#' \dontrun{
#' fit <- runModelingWorkflow("path/to/amRdata/output/orb", feature_types = "standard",
#'                            threads = 8, out_dir = "path/to/output")
#' fit$performance
#' }
#' @export
runModelingWorkflow <- function(orb_dir = NULL,
                                feature_types = "standard",
                                model_types = NULL,
                                threads = 1L,
                                out_dir = NULL,
                                seeds = c(123L, 456L, 789L),
                                labels = c("real", "shuffled"),
                                save_fits = FALSE,
                                save_state = FALSE,
                                min_genomes = 40,
                                n_fold = 5,
                                min_test_genomes = 5,
                                min_test_minority = 2,
                                overwrite = FALSE,
                                reuse_fits = TRUE,
                                resume = FALSE,
                                BPPARAM = NULL,
                                min_fold_minority = 2,
                                min_analysis_minority = 8,
                                min_mdr_label_genomes = 25,
                                mdr_min_classes = 2,
                                easy_mode = FALSE,
                                rescore = TRUE,
                                execution = c("fit", "prepare", "collect"),
                                n_chunks = 100L) {
  call <- rlang::current_env()
  started <- Sys.time()
  .checkArgFlag(easy_mode, "easy_mode", call = call)
  .checkArgFlag(rescore, "rescore", call = call)
  if (easy_mode) {
    if ((!missing(feature_types) && !identical(feature_types, "quick")) ||
        (!missing(model_types) && !identical(model_types, "vanilla")) ||
        (!missing(labels) && !identical(labels, "real")) ||
        (!missing(seeds) && length(seeds) != 1L)) {
      .amrAbort("invalid_argument", paste0(
        "Easy mode requires feature_types = \"quick\", model_types = \"vanilla\", ",
        "labels = \"real\", and a single seed. Use easy_mode = FALSE to customize these settings."
      ), call = call)
    }
    feature_types <- "quick"
    model_types <- "vanilla"
    labels <- "real"
    if (missing(seeds)) seeds <- 123L
    message("Easy mode fits vanilla gene, protein, and Pfam domain models using one seed, ",
            "without shuffled-label controls or MAGIC rescoring. ",
            "We recommend following these runs with more rigorous analyses, ",
            "using 3 or more seeds, shuffled label controls, and additional features.")
  }
  .checkArgCount(threads, "threads", call = call)
  if (!is.null(BPPARAM) && !missing(threads)) {
    .amrAbort("invalid_argument", "Supply either `threads` or `BPPARAM`, not both.", call = call)
  }
  if (is.null(BPPARAM)) {
    BPPARAM <- if (threads == 1L) BiocParallel::SerialParam() else if (.Platform$OS.type == "windows") {
      BiocParallel::SnowParam(workers = threads)
    } else BiocParallel::MulticoreParam(workers = threads)
  }
  # Checked first, so a bad value doesn't surface after the matrices are built.
  .checkArgCount(mdr_min_classes, "mdr_min_classes", min = 2, call = call)
  .checkArgSeeds(seeds, call = call)
  .checkArgLabels(labels, call = call)
  .checkArgFlag(overwrite, "overwrite", call = call)
  .checkArgClass(BPPARAM, "BPPARAM", "BiocParallelParam", "BiocParallel", call = call)
  .checkArgFlag(reuse_fits, "reuse_fits", call = call)
  .checkArgFlag(resume, "resume", call = call)
  .checkArgFlag(save_fits, "save_fits", call = call)
  .checkArgFlag(save_state, "save_state", call = call)
  execution <- match.arg(execution)
  .checkArgCount(n_chunks, "n_chunks", call = call)
  if ((execution == "collect" && !resume) || (execution == "prepare" && resume)) {
    .amrAbort("invalid_argument", "Collect requires resume = TRUE; prepare requires resume = FALSE.", call = call)
  }
  if (resume && (is.null(out_dir) || overwrite)) {
    .amrAbort("invalid_argument", "Resuming requires `out_dir` and `overwrite = FALSE`.", call = call)
  }

  if (resume) {
    out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
    snapshot <- .runPaths(out_dir)$workflow
    if (!file.exists(snapshot)) {
      .amrAbort("invalid_argument", c(
        "No saved workflow state was found; this run cannot be resumed.",
        i = paste0("Expected: ", snapshot),
        i = "Completed workflows remove state unless save_state = TRUE; supply the interrupted run directory."
      ), call = call)
    }
    saved <- readRDS(snapshot)
    # The normal default changed from one seed to three. Reuse the saved seeds
    # when the caller omits them, including for legacy interrupted workflows.
    if (missing(seeds) && !easy_mode) seeds <- saved$settings$seeds
    if (is.null(orb_dir)) orb_dir <- saved$settings$orb_dir
  }
  if (is.null(orb_dir)) orb_dir <- .selectORB(call, argument = "orb_dir")

  settings <- list(
    orb_dir = normalizePath(orb_dir, winslash = "/", mustWork = TRUE),
    model_types = model_types, seeds = as.integer(seeds), labels = labels,
    min_genomes = min_genomes, n_fold = n_fold,
    min_test_genomes = min_test_genomes, min_test_minority = min_test_minority,
    min_fold_minority = min_fold_minority, min_analysis_minority = min_analysis_minority,
    reuse_fits = reuse_fits, feature_types = feature_types, save_fits = save_fits,
        min_mdr_label_genomes = min_mdr_label_genomes, mdr_min_classes = mdr_min_classes,
    save_state = save_state, easy_mode = easy_mode, rescore = rescore
  )
  # Legacy snapshots already have an eligibility plan. Preserve their settings
  # when resuming without an explicit request to change the new CV thresholds.
  if (resume && is.null(saved$settings$min_fold_minority) && missing(min_fold_minority)) {
    settings$min_fold_minority <- NULL
  }
  if (resume && is.null(saved$settings$min_analysis_minority) && missing(min_analysis_minority)) {
    settings$min_analysis_minority <- NULL
  }
  if (resume && is.null(saved$settings$min_mdr_label_genomes) && missing(min_mdr_label_genomes)) settings$min_mdr_label_genomes <- NULL
  if (resume && is.null(saved$settings$mdr_min_classes) && missing(mdr_min_classes)) settings$mdr_min_classes <- NULL
  if (resume && is.null(saved$settings$easy_mode) && missing(easy_mode)) settings$easy_mode <- NULL
  if (resume && is.null(saved$settings$rescore) && missing(rescore)) settings$rescore <- NULL
  manifest_path <- NULL
  run_failed <- TRUE
  on.exit(
    if (run_failed && !is.null(manifest_path)) {
      .manifestFinish(
        manifest_path, status = "failed",
        error = "runModelingWorkflow() exited before successful completion."
      )
    },
    add = TRUE
  )

  if (resume) {
    out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
    paths <- .runPaths(out_dir)
    if (!identical(saved$settings, settings) ||
        !identical(normalizePath(saved$built$out_dir, winslash = "/", mustWork = TRUE),
                   normalizePath(out_dir, winslash = "/", mustWork = TRUE))) {
      .amrAbort("invalid_argument", "The saved workflow does not match these arguments or output directory.", call = call)
    }
    built <- saved$built
    orb <- built$eligibility$profile$orb
    manifest_path <- .manifestStart(out_dir, orb, settings, resume = TRUE)
    .manifestArtifact(
      manifest_path, "amRviz_input", status = "building",
      details = list(producer = "amRml",
                     producer_run_id = .manifestCurrentRunId(manifest_path),
                     directory = ".", path_mode = "relative_to_manifest")
    )
    .manifestEvent(manifest_path, "Started resumed modeling workflow.")
    reused <- list(
      read_orb = orb$manifest_path,
      profile_orb = orb$metadata_parquet,
      determine_eligibility = paths$eligibility,
      build_matrices = c(paths$matrices, paths$matrix_index,
                         list.files(paths$features, full.names = TRUE)),
      expand_tasks = paths$workflow
    )
    for (name in names(reused)) {
      .manifestStage(
        manifest_path, name, outputs = reused[[name]],
        message = "Resumed: reused successful output from a previous run."
      )
    }
  } else {
    orb <- readORB(orb_dir)
    if (identical(feature_types, "full") && "struct" %in% orb$feature_tables$feature_type) {
      warning("The 'full' feature preset includes struct models, which are very slow and can significantly increase runtime.",
              call. = FALSE, immediate. = TRUE)
    }
    out_dir <- out_dir %||% .defaultOutDir(orb$directory)
    .checkArgOutDir(out_dir, orb$directory, overwrite, call = call)
    out_dir <- .resolvePath(out_dir, call = call)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    paths <- .runPaths(out_dir)
    manifest_path <- .manifestStart(out_dir, orb, settings, resume = FALSE)
    producer_run_id <- .manifestCurrentRunId(manifest_path)
    .manifestArtifact(
      manifest_path, "amRviz_input", status = "building",
      details = list(producer = "amRml", producer_run_id = producer_run_id,
                     directory = ".", path_mode = "relative_to_manifest")
    )
    .manifestStage(
      manifest_path, "read_orb", inputs = orb$manifest_path,
      metrics = list(dataset_id = orb$dataset_id, producer_run_id = orb$producer_run_id)
    )

    .manifestStage(manifest_path, "profile_orb", status = "running", inputs = orb$metadata_parquet)
    profile <- profileORB(orb, mdr_min_classes = mdr_min_classes)
    .manifestStage(
      manifest_path, "profile_orb", inputs = orb$metadata_parquet,
      metrics = list(n_genomes = nrow(profile$genomes), n_targets = nrow(profile$targets))
    )

    .manifestStage(manifest_path, "determine_eligibility", status = "running")
    eligibility <- eligibleScopes(
      profile, modes = model_types, min_genomes = min_genomes, n_fold = n_fold,
      min_test_genomes = min_test_genomes, min_test_minority = min_test_minority,
      min_fold_minority = min_fold_minority, min_analysis_minority = min_analysis_minority,
      min_mdr_label_genomes = min_mdr_label_genomes
    )
    .manifestStage(manifest_path, "build_matrices", status = "running",
                   parameters = list(feature_types = feature_types))
    built <- buildMatrices(planMatrices(eligibility, feature_types = feature_types),
                           out_dir, overwrite = overwrite)
    out_dir <- built$out_dir
    paths <- .runPaths(out_dir)
    .writeParquet(eligibility$scopes, paths$eligibility)
    .manifestStage(
      manifest_path, "determine_eligibility", outputs = paths$eligibility,
      parameters = eligibility$settings,
      metrics = list(n_scopes = nrow(eligibility$scopes),
                     n_eligible = sum(eligibility$scopes$eligible))
    )
    .manifestStage(
      manifest_path, "build_matrices",
      inputs = orb$feature_tables$path[
        orb$feature_tables$feature_type %in% built$matrices$feature_type
      ],
      outputs = paths$matrices,
      parameters = list(feature_types = feature_types, storage = built$storage),
      metrics = list(n_built = sum(built$matrices$status == "built"),
                     n_skipped = sum(built$matrices$status == "skipped"))
    )
    # Save the original task inputs so an interrupted run can be resumed without
    # rebuilding matrices or reassigning its genome splits.
    .saveFitCheckpoint(
      list(settings = settings, built = built), paths$workflow
    )
  }

  producer_run_id <- .manifestCurrentRunId(manifest_path)
  eligibility <- built$eligibility
  tasks <- expandTasks(built, seeds = seeds, labels = labels)
  if (!resume) {
    .manifestStage(
      manifest_path, "expand_tasks",
      parameters = list(seeds = I(as.integer(seeds)), labels = I(labels)),
      metrics = list(n_tasks = nrow(tasks$tasks))
    )
  }

  fit_parameters <- list(
    save_fits = save_fits, save_state = save_state, reuse_fits = reuse_fits,
    backend = class(BPPARAM)[[1]], workers = BiocParallel::bpnworkers(BPPARAM),
    resume = resume
  )
  .manifestStage(
    manifest_path, "fit_models", status = "running",
    parameters = fit_parameters
  )
  fit_started <- Sys.time()
  fit <- fitModels(
    tasks,
    save_fits = save_fits, save_state = TRUE,
    overwrite = overwrite, BPPARAM = BPPARAM, reuse_fits = reuse_fits, resume = resume,
    execution = execution, n_chunks = n_chunks
  )
  if (execution == "prepare") {
    .manifestEvent(manifest_path, "Array prepared; waiting for workers and collection.",
                   details = list(n_chunks = length(fit$chunks)))
    .manifestFinish(manifest_path, status = "prepared")
    run_failed <- FALSE
    message("Prepared ", length(fit$chunks), " array chunks in ", out_dir)
    return(invisible(fit))
  }
  paths <- .runPaths(out_dir)
  fit_outputs <- c(paths$tasks, paths$performance, paths$importance, paths$predictions, paths$notes)
  if (save_fits) fit_outputs <- c(fit_outputs, list.files(paths$models, full.names = TRUE))
  .manifestStage(
    manifest_path, "fit_models", parameters = fit_parameters, outputs = fit_outputs,
    metrics = list(
      n_tasks = nrow(fit$tasks),
      n_fitted = sum(fit$tasks$status == "fitted"),
      n_failed = sum(fit$tasks$status == "failed"),
      n_notes = nrow(fit$notes),
      
      n_collapsed = sum(
        fit$performance$collapsed_model,
        na.rm = TRUE
      ),
      
      n_all_susceptible = sum(
        fit$performance$prediction_behavior ==
          "all_susceptible",
        na.rm = TRUE
      ),
      
      n_all_resistant = sum(
        fit$performance$prediction_behavior ==
          "all_resistant",
        na.rm = TRUE
      ),
      
      elapsed_seconds = as.numeric(difftime(Sys.time(), fit_started, units = "secs")))
  )

  .manifestStage(manifest_path, "export_merged_results", status = "running")
  exportMergedResults(fit, overwrite = overwrite || resume)
  merged <- list.files(paths$merged, recursive = TRUE, full.names = TRUE)
  .manifestStage(manifest_path, "export_merged_results", outputs = merged,
                 metrics = list(n_files = length(merged)))

  fit <- .rescoreWorkflow(fit, manifest_path, enabled = rescore && !easy_mode,
                          reason = if (easy_mode) "Easy mode does not run MAGIC rescoring." else
                            "Feature rescoring was disabled with rescore = FALSE.")

  if (!save_state) unlink(paths$state, recursive = TRUE)

  .manifestArtifact(
    manifest_path, "amRviz_input", status = "ready",
    details = list(
      producer = "amRml", directory = ".", path_mode = "relative_to_manifest",
      producer_run_id = producer_run_id,
      source_orb_manifest = orb$manifest_path,
      metadata = list(eligibility = "metadata/eligibility.parquet",
                      matrices = "metadata/matrices.parquet"),
      results = list(tasks = "results/tasks.parquet",
                     performance = "results/performance.parquet",
                     importance = "results/importance.parquet",
                     predictions = "results/predictions.parquet",
                     notes = "results/notes.parquet",
                     merged = "results/merged"),
      models = if (save_fits) "results/models" else NULL,
      modes = I(eligibility$modes$mode_id), seeds = I(as.integer(seeds)), labels = I(labels),
      feature_types = I(unique(built$matrices$feature_type)),
      encodings = I(unique(built$matrices$encoding)),
      parquet_compression = "zstd"
    )
  )
  .manifestArtifact(
    manifest_path, "resume_state", status = if (save_state) "ready" else "removed",
    details = list(
      producer = "amRml", producer_run_id = producer_run_id,
      directory = "state", workflow = "state/workflow.rds",
      matrix_index = "state/matrix_index.rds", feature_stores = "state/features",
      completion_receipt = "state/fit-completion.rds", rds_compression = "none",
      retained = save_state
    )
  )
  summary <- sprintf("Modeling workflow complete: %d of %d models fitted; %d failed. %s elapsed.",
                     sum(fit$tasks$status == "fitted"), nrow(fit$tasks),
                     sum(fit$tasks$status == "failed"), .elapsedTime(started))
  
  perf <- fit$performance
  
  n_collapsed <- sum(perf$collapsed_model, na.rm = TRUE)
  
  n_all_resistant <- sum(
    perf$prediction_behavior == "all_resistant",
    na.rm = TRUE
  )
  
  n_all_susceptible <- sum(
    perf$prediction_behavior == "all_susceptible",
    na.rm = TRUE
  )
  
  diagnostic_summary <- paste(
    "Model diagnostics:",
    sprintf("Collapsed models: %d", n_collapsed),
    sprintf(
      "  All susceptible: %d | All resistant: %d",
      n_all_susceptible,
      n_all_resistant
    ),
    sep = "\n"
  )
  
  .manifestEvent(
    manifest_path,
    sprintf(
      "Model diagnostics: %d collapsed models (%d all-susceptible, %d all-resistant).",
      n_collapsed,
      n_all_susceptible,
      n_all_resistant
    ),
    details = list(
      collapsed_models = n_collapsed,
      all_susceptible = n_all_susceptible,
      all_resistant = n_all_resistant
    )
  )
  
  .manifestEvent(manifest_path, summary, details = list(
    n_tasks = nrow(fit$tasks), n_failed = sum(fit$tasks$status == "failed"),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
  ))
  .manifestFinish(manifest_path, status = "success")
  
  run_failed <- FALSE
  
  message(
    summary,
    "\n",
    "\nOutputs: ",
    out_dir
  )
  
  invisible(fit)
}
