# The single entry point: every stage, from the ORB to fitted models.

#' Run the modeling workflow
#'
#' Runs every stage on an ORB: [readORB()], [profileORB()], [eligibleScopes()],
#' [planMatrices()], [buildMatrices()], [expandTasks()] and [fitModels()]. Each
#' can also be called on its own. Old amRml's standard run is a choice of
#' `modes`; its intense run adds `seeds` and `labels`.
#'
#' @details
#' Everything goes to `out_dir`: `eligibility.parquet`, `matrices.parquet`,
#' `matrices/`, `tasks.parquet` and `results/`. amRdata's files are only read.
#'
#' @param orb_dir An amRdata ORB directory.
#' @param out_dir Folder to write to. By default a new run folder,
#'   `<ORB>/amRml/<date-time>/`.
#' @param modes Mode IDs from [modelingModes()], or `NULL` for all of them.
#' @param seeds,labels Passed to [expandTasks()].
#' @param min_genomes,n_fold,min_test_genomes,min_test_minority,min_label_genomes Passed to
#'   [eligibleScopes()].
#' @param overwrite Replace an earlier run's matrices and results in `out_dir`;
#'   if `FALSE`, an `out_dir` holding earlier output is refused. The default run
#'   folder is always new.
#' @param BPPARAM BiocParallel backend passed to [fitModels()]. Defaults to a
#'   eight-worker cluster with 128 chunks of work to chew through.
#'
#' @return The `amr_fit` from [fitModels()]; every earlier stage's result can be
#'   reached from it (`$built`, `$built$eligibility`, ...).
#'
#' @section Errors:
#' Any error from the stages it runs, all with class `amrml_error`.
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
                                labels = "real",
                                min_genomes = 40,
                                n_fold = 5,
                                min_test_genomes = 5,
                                min_test_minority = 2,
                                min_label_genomes = 25,
                                overwrite = FALSE,
                                BPPARAM = BiocParallel::SnowParam(workers = 8L, tasks = 128L)) {
  call <- rlang::current_env()

  # Checked first, so a bad value doesn't surface after the matrices are built.
  .checkArgSeeds(seeds, call = call)
  .checkArgLabels(labels, call = call)
  .checkArgClass(BPPARAM, "BPPARAM", "BiocParallelParam", "BiocParallel::SnowParam()", call = call)

  orb <- readORB(orb_dir)
  out_dir <- out_dir %||% .defaultOutDir(orb$directory)
  eligibility <- eligibleScopes(
    profileORB(orb),
    modes = modes, min_genomes = min_genomes, n_fold = n_fold,
    min_test_genomes = min_test_genomes, min_test_minority = min_test_minority,
    min_label_genomes = min_label_genomes
  )
  built <- buildMatrices(planMatrices(eligibility), out_dir, overwrite = overwrite)
  arrow::write_parquet(eligibility$scopes, file.path(out_dir, "eligibility.parquet"))

  fitModels(expandTasks(built, seeds = seeds, labels = labels),
            overwrite = overwrite, BPPARAM = BPPARAM)
}
