# The single entry point: every stage, from the ORB to fitted models.

#' Run the modeling workflow
#'
#' Runs every stage on an ORB: [readORB()], [profileORB()], [eligibleScopes()],
#' [planMatrices()], [buildMatrices()], [expandTasks()] and [fitModels()]. Each
#' can also be called on its own. By default every task is fitted on both the real
#' and the shuffled labels; `labels = "real"` skips the shuffled controls.
#'
#' @details
#' Everything goes to `out_dir`: `eligibility.parquet`, `matrices.parquet`,
#' `features/`, `matrix_index.rds`, `workflow.rds`, `tasks.parquet` and `results/`.
#' amRdata's files are only read.
#'
#' @param orb_dir An amRdata ORB directory.
#' @param out_dir Folder to write to. By default a new run folder,
#'   `<ORB>/amRml/<date-time>/`.
#' @param modes Mode IDs from [modelingModes()], or `NULL` for all of them.
#' @param feature_types Feature types to model: `"quick"`, `"standard"` (default),
#'   `"full"`, or an explicit character vector of ORB feature types.
#'   Selected before matrices are planned or built; see [planMatrices()].
#' @param seeds Passed to [expandTasks()].
#' @param labels `"real"`, `"shuffled"` or both (the default), passed to
#'   [expandTasks()]. 
#' @param min_genomes,n_fold,min_test_genomes,min_test_minority,min_label_genomes Passed to
#'   [eligibleScopes()].
#' @param mdr_min_classes Passed to [profileORB()]: fewest drug classes a genome
#'   must be Resistant in to be MDR; at least 2.
#' @param overwrite Replace an earlier run's matrices and results in `out_dir`;
#'   if `FALSE`, an `out_dir` holding earlier output is refused. The default run
#'   folder is always new.
#' @param resume Resume interrupted model fitting in an existing `out_dir`,
#'   using its saved matrices and completed task groups, or reopen unchanged
#'   completed results. Supply the same modeling arguments and the exact run
#'   folder (e.g. `<ORB>/amRml/<date-time>/`); `BPPARAM` may be changed.
#'   `overwrite` must be `FALSE`.
#' @param BPPARAM Passed to [fitModels()]. May be changed when resuming.
#' @param reuse_fits Passed to [fitModels()]; reuse compatible cross-drug fits.
#'
#' @return The `amr_fit` from [fitModels()]; every earlier stage's result can be
#'   reached from it (`$built`, `$built$eligibility`, ...).
#'
#' @section Errors:
#' Errors from any stage, including `amrml_invalid_argument` for mismatched resume
#' settings, propagate to the caller.
#'
#' @examples
#' \dontrun{
#' fit <- runModelingWorkflow("path/to/amRdata/output", "path/to/output")
#' fit$performance
#' }
#' @export
runModelingWorkflow <- function(orb_dir,
                                out_dir = NULL,
                                modes = NULL,
                                seeds = 123L,
                                labels = c("real", "shuffled"),
                                min_genomes = 40,
                                n_fold = 5,
                                min_test_genomes = 5,
                                min_test_minority = 2,
                                min_label_genomes = 25,
                                mdr_min_classes = 2,
                                overwrite = FALSE,
                                BPPARAM = BiocParallel::SnowParam(workers = 8L, tasks = 128L),
                                reuse_fits = TRUE,
                                feature_types = "standard",
                                resume = FALSE) {
  call <- rlang::current_env()
  # Checked first, so a bad value doesn't surface after the matrices are built.
  .checkArgSeeds(seeds, call = call)
  .checkArgLabels(labels, call = call)
  .checkArgFlag(overwrite, "overwrite", call = call)
  .checkArgClass(BPPARAM, "BPPARAM", "BiocParallelParam", "BiocParallel", call = call)
  .checkArgCount(mdr_min_classes, "mdr_min_classes", min = 2, call = call)
  .checkArgFlag(reuse_fits, "reuse_fits", call = call)
  .checkArgFlag(resume, "resume", call = call)
  if (resume && (is.null(out_dir) || overwrite)) {
    .amrAbort("invalid_argument", "Resuming requires `out_dir` and `overwrite = FALSE`.", call = call)
  }

  settings <- list(
    orb_dir = normalizePath(orb_dir, winslash = "/", mustWork = TRUE),
    modes = modes, seeds = as.integer(seeds), labels = labels,
    min_genomes = min_genomes, n_fold = n_fold,
    min_test_genomes = min_test_genomes, min_test_minority = min_test_minority,
    min_label_genomes = min_label_genomes, mdr_min_classes = mdr_min_classes,
    reuse_fits = reuse_fits, feature_types = feature_types
  )
  if (resume) {
    snapshot <- file.path(out_dir, "workflow.rds")
    if (!file.exists(snapshot)) {
      .amrAbort("invalid_argument", "No saved workflow state was found; this run cannot be resumed.", call = call)
    }
    saved <- readRDS(snapshot)
    if (!identical(saved$settings, settings) ||
        !identical(normalizePath(saved$built$out_dir, winslash = "/", mustWork = TRUE),
                   normalizePath(out_dir, winslash = "/", mustWork = TRUE))) {
      .amrAbort("invalid_argument", "The saved workflow does not match these arguments or output directory.", call = call)
    }
    built <- saved$built
  } else {
    orb <- readORB(orb_dir)
    out_dir <- out_dir %||% .defaultOutDir(orb$directory)
    eligibility <- eligibleScopes(
      profileORB(orb, mdr_min_classes = mdr_min_classes),
      modes = modes, min_genomes = min_genomes, n_fold = n_fold,
      min_test_genomes = min_test_genomes, min_test_minority = min_test_minority,
      min_label_genomes = min_label_genomes
    )
    built <- buildMatrices(planMatrices(eligibility, feature_types = feature_types),
                           out_dir, overwrite = overwrite)
    out_dir <- built$out_dir
    arrow::write_parquet(eligibility$scopes, file.path(out_dir, "eligibility.parquet"))
    # Save the original task inputs so an interrupted run can be resumed without
    # rebuilding matrices or reassigning its genome splits.
    .saveFitCheckpoint(list(settings = settings, built = built), file.path(out_dir, "workflow.rds"))
  }

  fit <- fitModels(
    expandTasks(built, seeds = seeds, labels = labels),
    overwrite = overwrite, BPPARAM = BPPARAM, reuse_fits = reuse_fits, resume = resume
  )
  fit
}
