# step_scale() and step_normalize(), applied to all columns at once instead of one by one.

.stepScaleNew <- function(terms, center, trained, means, sds, id) {
  recipes::step(subclass = "amr_scale", terms = terms, center = center, role = NA,
                trained = trained, means = means, sds = sds, skip = FALSE, id = id)
}

# Scales the selected columns; center = TRUE also centers them, like step_normalize().
.stepScale <- function(recipe, ..., center = FALSE) {
  # Same random id draw as recipes' steps, so the folds drawn after it don't change.
  id <- recipes::rand_id(if (center) "normalize" else "scale")
  recipes::add_step(recipe, .stepScaleNew(rlang::enquos(...), center, FALSE, NULL, NULL, id))
}

# Computed as recipes computes them; a constant column keeps sd 1.
#' @exportS3Method recipes::prep
prep.step_amr_scale <- function(x, training, info = NULL, ...) {
  columns <- recipes::recipes_eval_select(x$terms, training, info)
  data <- training[, columns]
  sds <- sqrt(vapply(data, stats::var, numeric(1), na.rm = TRUE))
  sds[sds < .Machine$double.eps] <- 1
  means <- if (x$center) colMeans(data, na.rm = TRUE)
  .stepScaleNew(x$terms, x$center, TRUE, means, sds, x$id)
}

#' @exportS3Method recipes::bake
bake.step_amr_scale <- function(object, new_data, ...) {
  columns <- names(object$sds)
  recipes::check_new_data(columns, object, new_data)
  if (!length(columns)) return(new_data)
  sparse <- !object$center &&
    all(vapply(new_data[columns], sparsevctrs::is_sparse_vector, logical(1)))
  if (sparse) {
    # One sparse-matrix division: dividing 30,000 sparse vectors one by one takes seconds.
    values <- sparsevctrs::coerce_to_sparse_matrix(new_data[columns])
    values@x <- values@x / rep(unname(object$sds), diff(values@p))
    new_data[columns] <- sparsevctrs::coerce_to_sparse_tibble(values)
    return(new_data)
  }
  new_data[columns] <- lapply(columns, function(column) {
    x <- new_data[[column]]
    if (object$center) x <- x - object$means[column]
    x / object$sds[column]
  })
  new_data
}

#' @exportS3Method
print.step_amr_scale <- function(x, width = max(20, options()$width - 30), ...) {
  recipes::print_step(names(x$sds), x$terms, x$trained, "Scaling for", width = width)
  invisible(x)
}
