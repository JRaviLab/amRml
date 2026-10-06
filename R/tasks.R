# Expanding built matrices into the tasks to fit.

#' List the tasks to fit
#'
#' Lists every model to fit, before any is fitted: one task per built matrix,
#' label kind and seed. Shuffled tasks are the no-signal baseline: asked for,
#' they're added for every model. As in old amRml, their labels are shuffled
#' before the split, so their test genomes differ from their real twin's. More
#' seeds repeat every model with a different split.
#'
#' @param built An `amr_matrices` from [buildMatrices()].
#' @param seeds Seeds; each gives every model its own split and folds.
#' @param labels `"real"` (the default, as in old amRml's standard run),
#'   `"shuffled"`, or both.
#'
#' @return An `amr_tasks` list:
#'   * `built`: the built matrices the tasks use.
#'   * `tasks`: one row per task: `task_id`, `matrix_id`, `labels`, `seed`, and
#'     the `n_fold` and `holdout` from eligibility's settings.
#'
#' @section Errors:
#' All errors have class `amrml_error`, plus `amrml_invalid_argument`: `built`
#' is not from [buildMatrices()], `seeds` are not distinct whole numbers, or
#' `labels` are not `"real"` and/or `"shuffled"`.
#'
#' @examples
#' \dontrun{
#' tasks <- expandTasks(buildMatrices(plan, "path/to/output"), seeds = c(31, 79, 51))
#' tasks$tasks
#' }
#' @export
expandTasks <- function(built, seeds = 123L, labels = "real") {
  call <- rlang::current_env()
  .checkArgClass(built, "built", "amr_matrices", "buildMatrices()", call = call)
  .checkArgSeeds(seeds, call = call)
  .checkArgLabels(labels, call = call)

  matrix_ids <- built$matrices$matrix_id[built$matrices$status == "built"]
  tasks <- expand.grid(
    labels = labels, seed = as.integer(seeds), matrix_id = matrix_ids,
    stringsAsFactors = FALSE
  )
  settings <- built$eligibility$settings

  tasks <- tibble::tibble(
    task_id = paste(tasks$matrix_id, tasks$labels, tasks$seed, sep = "__"),
    matrix_id = tasks$matrix_id,
    labels = tasks$labels,
    seed = tasks$seed,
    n_fold = as.integer(settings$n_fold),
    holdout = settings$holdout
  )

  structure(list(built = built, tasks = tasks), class = "amr_tasks")
}


#' @export
print.amr_tasks <- function(x, ...) {
  orb <- x$built$eligibility$profile$orb
  t <- x$tasks
  modes <- .taskModes(t, x$built$matrices)
  n <- table(factor(modes, unique(modes)))

  cat("<amr_tasks>", orb$dataset_id, "-", orb$dataset_label, "\n")
  cat(
    "  ", nrow(t), " tasks: ", dplyr::n_distinct(t$matrix_id), " matrices x ",
    dplyr::n_distinct(t$labels), " label kinds x ", dplyr::n_distinct(t$seed), " seeds\n",
    sep = ""
  )
  cat(.modeLines(names(n), sprintf("%d tasks", n)), sep = "\n")
  invisible(x)
}
