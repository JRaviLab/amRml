# Argument checks. Each throws `amrml_invalid_argument`.

# Throw an error unless `orb` is an ORB object from readORB().
.checkArgORB <- function(orb, call = rlang::caller_env()) {
  if (!inherits(orb, "amr_orb")) {
    .amrAbort(
      "invalid_argument",
      c(
        "`orb` must be an `amr_orb` from readORB().",
        x = paste0("Got class ", class(orb)[[1]], ".")
      ),
      observed = list(class = class(orb)),
      call = call
    )
  }

  invisible(TRUE)
}
