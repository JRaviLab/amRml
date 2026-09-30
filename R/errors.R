# amRml's error codes, and the functions that build and throw its errors.

# `stage` is the pipeline stage; `grain` says whether the error is about the
# whole call or one manifest file.
.RULES <- tibble::tribble(
  ~rule_id,                     ~stage,    ~grain,     ~description,

  # General
  "invalid_argument",           "input",   "call",     "Wrong type, length or class of argument.",

  # Reading the ORB: readORB()
  # Locating the ORB
  "orb_dir_not_found",          "orb",     "call",     "The ORB directory does not exist.",
  "manifest_not_found",         "orb",     "call",     "No amRdata manifest in the ORB directory.",

  # Manifest selection
  "manifest_unreadable",        "orb",     "manifest", "Manifest is not readable JSON.",
  "manifest_invalid",           "orb",     "manifest", "Manifest is not a valid amRdata manifest.",
  "orb_not_ready",              "orb",     "manifest", "No 'ready' amRml_input artifact.",
  "producer_run_invalid",       "orb",     "manifest", "Producer run is missing or failed.",
  "orb_no_usable_manifest",     "orb",     "call",     "No manifest declares a usable amRml input.",

  # Checking the chosen ORB's files
  "orb_moved",                  "orb",     "call",     "ORB is not where its manifest was written.",
  "orb_files_missing",          "orb",     "call",     "Declared files are missing or unreadable.",
  "orb_file_changed",           "orb",     "call",     "A declared file's size has changed.",

  # Feature tables
  "parquet_unreadable",         "orb",     "call",     "A declared parquet could not be read.",
  "orb_no_feature_tables",      "orb",     "call",     "No declared parquet is a feature table.",
  "orb_duplicate_feature_type", "orb",     "call",     "Feature tables share a feature ID column.",

  # Profiling the ORB: profileORB()
  "metadata_columns_missing",   "profile", "call",     "Metadata lacks a column profiling needs.",
  "stratum_inconsistent",       "profile", "call",     "A genome has two values for a stratum.",
  "no_labelled_genomes",        "profile", "call",     "No genome with features has a phenotype.",
)

# Build an error of class `amrml_<rule_id>` without throwing it. `message` is a
# headline plus rlang-style bullets; `call = NULL` means the error names no function.
.amrError <- function(rule_id,
                      message,
                      observed = NULL,
                      threshold = NULL,
                      call = NULL) {
  # An unknown rule is a bug in amRml, so it gets its own class.
  if (!rlang::is_string(rule_id) || !rule_id %in% .RULES$rule_id) {
    rlang::abort(
      c(
        "Internal error: unknown `rule_id`.",
        x = paste0("Got: ", paste(format(rule_id), collapse = ", "))
      ),
      class = c("amrml_internal_error", "amrml_error")
    )
  }

  rlang::error_cnd(
    class = c(paste0("amrml_", rule_id), "amrml_error"),
    message = message,
    rule_id = rule_id,
    stage = .RULES$stage[.RULES$rule_id == rule_id],
    observed = observed,
    threshold = threshold,
    call = call,
    use_cli_format = TRUE
  )
}

# Build an error and throw it, blaming the function whose environment is `call`.
.amrAbort <- function(..., call = rlang::caller_env()) {
  rlang::cnd_signal(.amrError(..., call = call))
}

# Name every element "x", so rlang shows each as a bullet in an error message.
# For example, listing all missing files in an ORB.
.bullets <- function(x) {
  rlang::set_names(x, rep("x", length(x)))
}
