# Argument checks. Each throws `amrml_invalid_argument`.

# Throw an error unless `value` has class `expected`, as returned by `from`.
.checkArgClass <- function(value, name, expected, from, call = rlang::caller_env()) {
  if (!inherits(value, expected)) {
    .amrAbort(
      "invalid_argument",
      c(
        paste0("`", name, "` must be an `", expected, "` from ", from, "."),
        x = paste0("Got class ", class(value)[[1]], ".")
      ),
      observed = list(class = class(value)),
      call = call
    )
  }

  invisible(TRUE)
}

# Throw an error unless `value` is a single finite whole number of at least `min`.
.checkArgCount <- function(value, name, min = 1, call = rlang::caller_env()) {
  ok <- is.numeric(value) && length(value) == 1L && is.finite(value) &&
    value == round(value) && value >= min

  if (!ok) {
    .amrAbort(
      "invalid_argument",
      paste0("`", name, "` must be a single finite whole number of at least ", min, "."),
      observed = list(value = value),
      call = call
    )
  }

  invisible(TRUE)
}

# Throw an error unless `modes` holds mode IDs from modelingModes().
.checkArgModes <- function(modes, call = rlang::caller_env()) {
  unknown <- if (is.character(modes)) setdiff(modes, .MODES$mode_id) else character()

  if (!is.character(modes) || anyNA(modes) || length(unknown)) {
    .amrAbort(
      "invalid_argument",
      c("`modes` must be mode IDs from modelingModes().", .bullets(unknown)),
      observed = list(unknown = unknown),
      call = call
    )
  }

  invisible(TRUE)
}
