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

# Throw an error unless `value` is TRUE or FALSE.
.checkArgFlag <- function(value, name, call = rlang::caller_env()) {
  if (!rlang::is_bool(value)) {
    .amrAbort("invalid_argument", paste0("`", name, "` must be TRUE or FALSE."), call = call)
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

# `path` with symlinks resolved, including when it, or its parents, don't exist yet. The
# part that doesn't exist can't hold a symlink, unless ".." steps back out of it; so ".."
# there is refused.
.resolvePath <- function(path, call = rlang::caller_env()) {
  path <- path.expand(path)
  missing <- character()
  while (!file.exists(path)) {
    missing <- c(basename(path), missing)
    path <- dirname(path)
  }

  if (".." %in% missing) {
    .amrAbort(
      "invalid_argument",
      "`out_dir` can't use \"..\" after a folder that doesn't exist yet.",
      call = call
    )
  }
  do.call(file.path, as.list(c(normalizePath(path), missing[missing != "."])))
}

# Throw an error unless `out_dir` is one writable path outside the ORB at `orb_dir`, or in
# its amRml/ folder, and holds no earlier output unless `overwrite` is TRUE.
.checkArgOutDir <- function(out_dir, orb_dir, overwrite, call = rlang::caller_env()) {
  if (!rlang::is_string(out_dir) || !nzchar(out_dir)) {
    .amrAbort("invalid_argument", "`out_dir` must be a single directory path.", call = call)
  }
  .checkArgFlag(overwrite, "overwrite", call = call)

  # The nearest folder that exists must be writable, or nothing can be created in it.
  existing <- path.expand(out_dir)
  while (!file.exists(existing)) {
    existing <- dirname(existing)
  }
  if (!dir.exists(existing) || file.access(existing, mode = 2) != 0) {
    .amrAbort(
      "invalid_argument",
      c("`out_dir` must be a folder you can write to.", x = paste0("Not writable: ", existing)),
      observed = list(out_dir = out_dir),
      call = call
    )
  }

  # Inside the ORB, only a run folder under its amRml/ folder; amRdata's files are read-only.
  out <- .resolvePath(out_dir, call = call)
  orb <- .resolvePath(orb_dir, call = call)
  in_orb <- out == orb || startsWith(out, paste0(orb, "/"))
  in_runs <- startsWith(out, paste0(file.path(orb, "amRml"), "/"))
  if (in_orb && !in_runs) {
    .amrAbort(
      "invalid_argument",
      c(
        "`out_dir` must be outside the ORB, or a run folder in its `amRml/` folder.",
        i = paste0("ORB: ", orb)
      ),
      observed = list(out_dir = out),
      call = call
    )
  }

  # Every stage's output, so a run never mixes new matrices with earlier results.
  outputs <- c(
    "eligibility.parquet", "matrices", "matrices.parquet", "tasks.parquet", "run.json",
    "results", "merged"
  )
  existing <- file.exists(file.path(out_dir, outputs))
  if (any(existing) && !overwrite) {
    .amrAbort(
      "invalid_argument",
      "`out_dir` already has output; set `overwrite = TRUE` to replace it.",
      observed = list(out_dir = out),
      call = call
    )
  }

  invisible(TRUE)
}

# Throw an error unless `seeds` is one or more distinct whole numbers.
.checkArgSeeds <- function(seeds, call = rlang::caller_env()) {
  ok <- is.numeric(seeds) && length(seeds) > 0 && all(is.finite(seeds)) &&
    all(seeds == round(seeds)) && !anyDuplicated(seeds)

  if (!ok) {
    .amrAbort(
      "invalid_argument",
      "`seeds` must be one or more distinct whole numbers.",
      observed = list(seeds = seeds),
      call = call
    )
  }

  invisible(TRUE)
}

# Throw an error unless `labels` is one or both of "real" and "shuffled".
.checkArgLabels <- function(labels, call = rlang::caller_env()) {
  ok <- is.character(labels) && length(labels) > 0 && !anyNA(labels) &&
    all(labels %in% c("real", "shuffled")) && !anyDuplicated(labels)

  if (!ok) {
    .amrAbort(
      "invalid_argument",
      "`labels` must be one or both of \"real\" and \"shuffled\".",
      observed = list(labels = labels),
      call = call
    )
  }

  invisible(TRUE)
}
