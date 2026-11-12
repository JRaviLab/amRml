# Helpers shared by more than one file.

# One "mode : label" line per mode, with the modes padded to one width.
.modeLines <- function(modes, labels) {
  paste0("  ", format(modes), " : ", labels, recycle0 = TRUE)
}

# A new run folder in the ORB's amRml/ folder, named by the time, and announced.
.defaultOutDir <- function(orb_dir) {
  out_dir <- file.path(orb_dir, "amRml", format(Sys.time(), "%Y-%m-%d_%H%M%S"))
  rlang::inform(paste0("Writing to ", out_dir))
  out_dir
}

# Each task's mode, looked up through its matrix.
.taskModes <- function(tasks, matrices) {
  matrices$mode_id[match(tasks$matrix_id, matrices$matrix_id)]
}

# BiocParallel's default backend, except where forking the R session is unsafe
# (Positron and RStudio refuse it): there, fresh R processes with as many workers.
.defaultBPPARAM <- function() {
  param <- BiocParallel::bpparam()
  no_fork <- nzchar(Sys.getenv("POSITRON")) || nzchar(Sys.getenv("RSTUDIO"))
  if (inherits(param, "MulticoreParam") && no_fork) {
    return(BiocParallel::SnowParam(workers = BiocParallel::bpnworkers(param)))
  }
  param
}
