# Profiling an amRdata ORB: each genome's phenotypes, and the counts per target.

# amRdata's metadata.parquet columns that profiling needs, named as they are used here.
.METADATA_COLUMNS <- c(
  genome_id = "genome.genome_id",
  phenotype = "genome_drug.resistant_phenotype",
  drug_abbr = "drug_abbr",
  class_abbr = "class_abbr"
)


# The metadata as character columns with empty strings as NA. Strata columns get
# their registry names.
.readMetadata <- function(path, call) {
  metadata <- arrow::read_parquet(path)
  missing_columns <- setdiff(.METADATA_COLUMNS, names(metadata))

  if (length(missing_columns)) {
    .amrAbort(
      "metadata_columns_missing",
      c("The metadata is missing columns profiling needs.", .bullets(missing_columns)),
      observed = list(missing = missing_columns),
      call = call
    )
  }

  columns <- c(.METADATA_COLUMNS, .STRATA)
  columns <- columns[columns %in% names(metadata)]
  metadata <- metadata[columns]
  names(metadata) <- names(columns)

  # as.character() also drops the unused factor levels of amRdata's year bins.
  metadata[] <- lapply(metadata, function(x) {
    x <- as.character(x)
    x[x %in% ""] <- NA_character_
    x
  })

  metadata
}

# The genomes in `metadata` that a feature table lacks, with the feature types lacking them.
.missingFeatures <- function(metadata, feature_tables) {
  genomes <- unique(metadata$genome_id)

  lacking <- lapply(feature_tables$path, function(f) {
    table <- arrow::read_parquet(f, col_select = "genome_id", as_data_frame = FALSE)
    setdiff(genomes, as.vector(unique(table$genome_id)))
  })

  missing <- tibble::tibble(
    genome_id = unlist(lacking),
    feature_type = rep(feature_tables$feature_type, lengths(lacking))
  )

  dplyr::summarise(
    dplyr::group_by(missing, .data$genome_id),
    missing_from = toString(.data$feature_type)
  )
}

# One row per genome and drug; repeated rows that agree count once. amRdata gives
# each drug one class, from a lookup table.
.drugPhenotypes <- function(metadata, call) {
  phenotypes <- dplyr::distinct(metadata[c("genome_id", "drug_abbr", "phenotype")])
  pairs <- phenotypes[c("genome_id", "drug_abbr")]
  conflicts <- unique(pairs[duplicated(pairs), , drop = FALSE])

  if (nrow(conflicts)) {
    .amrAbort(
      "phenotype_conflict",
      c(
        "Genomes are recorded as both Resistant and Susceptible for a drug.",
        .bullets(paste(conflicts$genome_id, conflicts$drug_abbr))
      ),
      observed = list(conflicts = conflicts),
      call = call
    )
  }

  phenotypes$class_abbr <- metadata$class_abbr[match(phenotypes$drug_abbr, metadata$drug_abbr)]
  phenotypes
}

# One row per genome and drug class: Resistant if any member drug is Resistant,
# otherwise Susceptible.
.classPhenotypes <- function(drug_phenotypes) {
  dplyr::summarise(
    dplyr::group_by(
      drug_phenotypes[!is.na(drug_phenotypes$class_abbr), , drop = FALSE],
      .data$genome_id, .data$class_abbr
    ),
    phenotype = if (any(.data$phenotype == "Resistant")) "Resistant" else "Susceptible",
    .groups = "drop"
  )
}

# One row per genome with its strata, which must have one value per genome.
.genomeTable <- function(metadata, call) {
  columns <- setdiff(names(metadata), c("drug_abbr", "class_abbr", "phenotype"))
  genomes <- dplyr::distinct(metadata[columns])
  repeated <- unique(genomes$genome_id[duplicated(genomes$genome_id)])

  if (length(repeated)) {
    .amrAbort(
      "stratum_inconsistent",
      c(
        "Strata have more than one value within a genome.",
        x = paste0("Genomes: ", toString(repeated))
      ),
      observed = list(genomes = repeated),
      call = call
    )
  }

  genomes[order(genomes$genome_id), , drop = FALSE]
}

# Counts per target over the whole dataset (stratum NA) and within each observed
# value of each stratum.
.targetCounts <- function(phenotypes, genomes) {
  scopes <- dplyr::bind_rows(
    tibble::tibble(
      genome_id = genomes$genome_id,
      stratum = NA_character_,
      stratum_value = NA_character_
    ),
    lapply(intersect(names(.STRATA), names(genomes)), function(stratum_name) {
      values <- genomes[[stratum_name]]
      tibble::tibble(
        genome_id = genomes$genome_id[!is.na(values)],
        stratum = stratum_name,
        stratum_value = values[!is.na(values)]
      )
    })
  )

  counts <- dplyr::summarise(
    dplyr::group_by(
      dplyr::inner_join(phenotypes, scopes, by = "genome_id", relationship = "many-to-many"),
      .data$unit, .data$target, .data$stratum, .data$stratum_value
    ),
    n_genomes = dplyr::n(),
    n_resistant = sum(.data$phenotype == "Resistant"),
    n_susceptible = sum(.data$phenotype == "Susceptible"),
    .groups = "drop"
  )

  dplyr::arrange(
    counts,
    !is.na(.data$stratum), .data$stratum, .data$stratum_value, .data$unit, .data$target
  )
}

#' Profile an amRdata ORB
#'
#' Labels each genome Resistant or Susceptible per drug and per drug class, and
#' counts the genomes behind each target, over the whole dataset and within each
#' value of each stratum (year and country).
#'
#' @details
#' A genome is excluded, and listed with the reason, if it has no Resistant or
#' Susceptible row or is missing from a feature table. A class is Resistant if
#' any member drug is, and otherwise Susceptible. Repeated rows for a genome and
#' drug count once if they agree; if they disagree, profiling stops.
#'
#' @section Errors:
#' All errors have class `amrml_error`, plus one of:
#'
#' * `amrml_invalid_argument`: `orb` is not from [readORB()].
#' * `amrml_metadata_columns_missing`: the metadata lacks a needed column.
#' * `amrml_phenotype_conflict`: a genome is recorded as both Resistant and
#'   Susceptible for a drug.
#' * `amrml_no_labelled_genomes`: no genome in every feature table has a
#'   usable phenotype.
#' * `amrml_stratum_inconsistent`: a genome has two values for a stratum.
#'
#' @param orb An `amr_orb` from [readORB()].
#'
#' @return An `amr_orb_profile` list:
#'   * `orb`: the ORB profiled.
#'   * `genomes`: each profiled genome and its strata values.
#'   * `phenotypes`: each genome's phenotype per `unit` (drug or drug class) and
#'     `target`.
#'   * `targets`: genome counts per target and stratum value; `stratum` is NA
#'     for the whole dataset.
#'   * `excluded_genomes`: each excluded genome and its `reason`
#'     (`"no_usable_phenotype"` or `"missing_features"`), with the feature types
#'     it is `missing_from`.
#'
#' @examples
#' \dontrun{
#' profile <- profileORB(readORB("path/to/amRdata/output"))
#' profile$targets
#' }
#' @export
profileORB <- function(orb) {
  call <- rlang::current_env()

  # Check the argument is an ORB from readORB().
  .checkArgClass(orb, "orb", "amr_orb", "readORB()", call = call)

  # Read the metadata and keep its Resistant and Susceptible rows.
  metadata <- .readMetadata(orb$metadata_parquet, call = call)
  usable <- metadata$phenotype %in% c("Resistant", "Susceptible") &
    !is.na(metadata$genome_id) & !is.na(metadata$drug_abbr)
  no_usable_phenotype <- setdiff(metadata$genome_id, c(NA, metadata$genome_id[usable]))
  metadata <- metadata[usable, , drop = FALSE]

  # Each genome's phenotype per drug.
  drug_phenotypes <- .drugPhenotypes(metadata, call = call)

  # Drop genomes that aren't in every feature table.
  missing <- .missingFeatures(metadata, orb$feature_tables)
  drug_phenotypes <- drug_phenotypes[!drug_phenotypes$genome_id %in% missing$genome_id, ]

  # Record every genome left out, and why.
  excluded <- dplyr::bind_rows(
    tibble::tibble(genome_id = no_usable_phenotype, reason = "no_usable_phenotype"),
    tibble::tibble(
      genome_id = missing$genome_id,
      reason = "missing_features",
      missing_from = as.character(missing$missing_from)
    )
  )
  excluded <- excluded[order(excluded$genome_id), , drop = FALSE]
  reasons <- table(excluded$reason)
  reasons <- paste0(names(reasons), ": ", reasons, collapse = ", ")

  if (!nrow(drug_phenotypes)) {
    .amrAbort(
      "no_labelled_genomes",
      c(
        "No genome in every feature table has a usable phenotype.",
        i = paste0("Excluded: ", reasons)
      ),
      observed = list(excluded = excluded),
      call = call
    )
  }

  if (nrow(excluded)) {
    rlang::inform(paste0(
      nrow(excluded), " genomes are excluded (", reasons, "); see `$excluded_genomes`."
    ))
  }

  # Each genome's phenotype per drug class.
  class_phenotypes <- .classPhenotypes(drug_phenotypes)
  phenotypes <- dplyr::bind_rows(
    tibble::tibble(
      genome_id = drug_phenotypes$genome_id,
      unit = "drug",
      target = drug_phenotypes$drug_abbr,
      phenotype = drug_phenotypes$phenotype
    ),
    tibble::tibble(
      genome_id = class_phenotypes$genome_id,
      unit = "drug_class",
      target = class_phenotypes$class_abbr,
      phenotype = class_phenotypes$phenotype
    )
  )

  # Collect the strata of each genome with a phenotype.
  profiled <- metadata[metadata$genome_id %in% phenotypes$genome_id, , drop = FALSE]
  genomes <- .genomeTable(profiled, call = call)

  # Count genomes per target, over the whole dataset and within each stratum value.
  targets <- .targetCounts(phenotypes, genomes)

  structure(
    list(
      orb = orb,
      genomes = genomes,
      phenotypes = phenotypes,
      targets = targets,
      excluded_genomes = excluded
    ),
    class = "amr_orb_profile"
  )
}


#' @export
print.amr_orb_profile <- function(x, ...) {
  whole <- x$targets[is.na(x$targets$stratum), , drop = FALSE]
  strata <- vapply(intersect(names(.STRATA), names(x$genomes)), function(s) {
    values <- x$genomes[[s]]
    paste0(
      s, " (values: ", dplyr::n_distinct(values, na.rm = TRUE),
      ", missing: ", sum(is.na(values)), ")"
    )
  }, character(1))

  # Wraps "target:resistant/susceptible" items to the console width.
  counts <- function(unit) {
    rows <- whole[whole$unit == unit, , drop = FALSE]
    items <- paste0(rows$target, ":", rows$n_resistant, "/", rows$n_susceptible)
    strwrap(paste(items, collapse = ", "), indent = 4, exdent = 4)
  }

  cat("<amr_orb_profile>", x$orb$dataset_id, "-", x$orb$dataset_label, "\n")
  cat("  genomes   :", nrow(x$genomes), "profiled,", nrow(x$excluded_genomes), "excluded\n")
  cat("  strata    :", paste(strata, collapse = "; "), "\n")
  cat("  Resistant/Susceptible genomes per drug:\n")
  cat(counts("drug"), sep = "\n")
  cat("  per drug class:\n")
  cat(counts("drug_class"), sep = "\n")
  invisible(x)
}
