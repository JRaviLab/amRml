#' Filter model fits to those passing quality thresholds
#'
#' @param all_performance_parquet The path to the 'all performance parquet' file
#' @param MCC_threshold The minimum non-shuffled MCC required to keep a model fit (default is \code{NULL}, no filtering)
#' @param compare_to_shuffled Logical indicating whether to require the non-shuffled MCC to exceed (outperform) the shuffled-label MCC for the same fit (default is \code{TRUE}). A fit with no shuffled counterpart (\code{shuffled_MCC} is \code{NA} after pivoting) always passes this check regardless.
#'
#' @returns a tibble of model fits (one row per species, drug_label, drug_or_class, seed, feature_type, feature_subtype, model, fit_penalty, fit_mixture) that pass the requested quality thresholds, with \code{nonshuffled_MCC}, \code{shuffled_MCC}, and \code{MCC_diff} columns added. This does not pick a single "best" fit per group; it filters out fits that fail the MCC/shuffled-comparison criteria, so multiple passing fits per group can remain.
#'
#' @keywords internal
#' @examples
#' filterOptimalModel(all_performance_parquet = "inst/extdata/all_perf.parquet")
filterOptimalModel <- function(all_performance_parquet,
  MCC_threshold = NULL,
  compare_to_shuffled = TRUE)
  {
  stopifnot(file.exists(all_performance_parquet))
  all_perf <- arrow::read_parquet(normalizePath(all_performance_parquet))

  if (!all(c(TRUE, FALSE) %in% unique(all_perf$shuffled))) {
    stop("The 'shuffled' column must contain both TRUE and FALSE values. Run runModelingPipelineIntense() ")
  }

  all_perf |>
    dplyr::select(
      species, drug_label, drug_or_class,
      seed, feature_type, feature_subtype,
      model, fit_penalty, fit_mixture,
      shuffled, mcc
    ) |>
    tidyr::pivot_wider(
      names_from = shuffled,
      values_from = mcc,
      names_prefix = "shuffled_"
    ) |>
    dplyr::rename(
      nonshuffled_MCC = shuffled_FALSE,
      shuffled_MCC = shuffled_TRUE
    ) |>
    dplyr::filter(!is.na(nonshuffled_MCC)) |>
    dplyr::mutate(
      MCC_diff = nonshuffled_MCC - shuffled_MCC
    ) |>
    dplyr::filter(
      if (!is.null(MCC_threshold)) (nonshuffled_MCC >= MCC_threshold) else TRUE,
      if (compare_to_shuffled) (MCC_diff > 0 | is.na(shuffled_MCC)) else TRUE
    )
}

#' Score features within each seed
#'
#' @param all_top_features_parquet The path to the Parquet file containing all top features with their importance scores.
#' @param core_contribution_threshold The cumulative-contribution cutoff, in \[0, 1\] (default is \code{0.75}, i.e. 75%), used to flag whether a feature falls within the "core" set of features that jointly account for that share of a seed's total importance.
#' @param exclude_feature_types Feature types to drop before any scoring happens (default is \code{NULL}, i.e. no feature types are excluded and struct is currently included). struct variables are composite IDs (e.g. \code{polA.group_211.group_2176}, three dot-joined gene/domain identifiers) representing a co-occurrence/structural motif rather than a single molecular entity like the other five scales, and its candidate-variable count (tens of thousands per group) dwarfs the other scales by orders of magnitude — pooling it into this rank_score/contribution machinery would compare a compound signal against five primary ones on an incomparable scale. Pass \code{"struct"} here to exclude it once struct is meant to be reserved for post-hoc biological annotation rather than scoring/ranking/thresholding.
#' @param filter_model Logical indicating whether to restrict scoring to (species, drug_label, drug_or_class, feature_type, feature_subtype, seed) groups that have at least one model fit passing \code{filterOptimalModel()}'s MCC/shuffled-comparison quality thresholds (default is \code{TRUE}). The join is not keyed on \code{model}/\code{fit_penalty}/\code{fit_mixture}; this is only safe because there is currently never more than one fit per (species, drug_label, drug_or_class, feature_type, feature_subtype, seed) group. If that assumption changes, this join (and the \code{all_perf} join below) would need a finer-grained key to avoid silently keeping/duplicating rows from multiple fits.
#' @param all_performance_parquet The path to the all performance parquet file. Always required — it is read unconditionally 
#' @param MCC_threshold The minimum MCC threshold passed through to \code{filterOptimalModel()} (default is \code{NULL}, no filtering)
#' @param compare_to_shuffled Logical indicating whether to compare the model to shuffled data, passed through to \code{filterOptimalModel()} (default is \code{TRUE})
#'
#' @returns a tibble of scored top features, one row per feature within each species/drug_label/drug_or_class/feature_type/feature_subtype/seed group, with the following columns added:
#' \itemize{
#'   \item \code{contribution}: the feature's importance divided by the sum of importance across all features in the group.
#'   \item \code{feat_return_ratio}: \code{n_feats_returned / n_feat} for the fit that produced this group. This assumes there is never more than one surviving fit per (species, drug_label, drug_or_class, seed, feature_type, feature_subtype) group; if that ever stops being true, this join would need a finer-grained key (\code{all_top_features_parquet} carries no \code{model}/\code{fit_penalty}/\code{fit_mixture} column to join on directly) or the ratio would need to be aggregated across fits before joining, to avoid fanning out and double-counting feature rows.
#'   \item \code{rank}: descending rank of \code{contribution} within the group; ties receive the average of the ranks they span.
#'   \item \code{n_features}: the number of rows (features) in the group.
#'   \item \code{rank_score}: \code{(n_features - rank) / (n_features - 1)}; ranges 0-1 with higher values indicating higher importance (a single-feature group scores 1).
#'   \item \code{cum_contrib}: the cumulative sum of \code{contribution} in descending order; features tied on \code{contribution} share the same \code{cum_contrib}, equal to the cumulative sum through the end of their tied block, so a tie is never split by arbitrary sort order.
#'   \item \code{in_core}: TRUE when \code{cum_contrib <= core_contribution_threshold}; a tied block that would push the cumulative total past the threshold is excluded in its entirety (conservative: stays at-or-under the threshold rather than overshooting it).
#' }
#'
#' @keywords internal
#' @examples
#' scoreFeaturesWithinSeed(all_top_features_parquet = "inst/extdata/all_top_features.parquet",
#' all_performance_parquet = "inst/extdata/all_perf.parquet")
#'
scoreFeaturesWithinSeed <- function(all_top_features_parquet,
  core_contribution_threshold = 0.75,
  exclude_feature_types = NULL,
  filter_model = TRUE,
  all_performance_parquet,
  MCC_threshold = NULL,
  compare_to_shuffled = TRUE)
  {
  # check for the all_perf.parquet and all_top_features.parquet files
  stopifnot(file.exists(all_top_features_parquet))
  stopifnot(file.exists(all_performance_parquet))

  all_top_features <- arrow::read_parquet(normalizePath(all_top_features_parquet)) |>
    dplyr::filter(!shuffled, !feature_type %in% exclude_feature_types)

  if (filter_model) {
    filtered_model <- filterOptimalModel(
      normalizePath(all_performance_parquet),
      MCC_threshold = MCC_threshold,
      compare_to_shuffled = compare_to_shuffled
    )

    all_top_features <- all_top_features |>
      dplyr::semi_join(
        filtered_model,
        by = dplyr::join_by(
          species, drug_label, drug_or_class,
          feature_type, feature_subtype, seed
        )
      )
  }

  all_perf <- arrow::read_parquet(normalizePath(all_performance_parquet)) |>
    dplyr::filter(!shuffled) |>
    dplyr::select(species, drug_label, drug_or_class, seed,
      feature_type, feature_subtype, fit_penalty, fit_mixture,
      mcc, n_feat, n_feats_returned) |>
    dplyr::mutate(
      feat_return_ratio = n_feats_returned / n_feat
    )

  # add different layers of scoring to the features within each seed
  scored_top_features <- all_top_features |>
    dplyr::select(
      species, drug_label, drug_or_class, seed,
      feature_type, feature_subtype, variable = Variable,
      importance = Importance, sign = Sign
    ) |>
    dplyr::mutate(
      variable = dplyr::case_when(
        feature_type == "protein" ~ sub("fig.", "fig|", variable, fixed = TRUE),
        feature_type == "AMRFinder" ~ sub(
          "^X", "",
          gsub("\\.NCBIFAM", "", variable)
        ),
        TRUE ~ variable
      )
    ) |>
    dplyr::left_join(all_perf, by = dplyr::join_by(species, drug_label, drug_or_class, seed,
      feature_type, feature_subtype)) |>
    dplyr::group_by(
      species,
      drug_label,
      drug_or_class,
      feature_type,
      feature_subtype,
      seed
    ) |>
    dplyr::mutate(
      contribution = importance / sum(importance, na.rm = TRUE)
    ) |>
    dplyr::arrange(dplyr::desc(contribution), .by_group = TRUE) |>
    dplyr::mutate(
      rank = rank(dplyr::desc(contribution), ties.method = "average"),
      n_features = dplyr::n(),
      rank_score = dplyr::if_else(
        n_features > 1,
        (n_features - rank) / (n_features - 1),
        1
      ),
      running_contrib = cumsum(contribution)
    ) |>
    dplyr::group_by(contribution, .add = TRUE) |>
    dplyr::mutate(
      cum_contrib = max(running_contrib),
      in_core = cum_contrib <= core_contribution_threshold
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-running_contrib)

  return(scored_top_features)
}

#' Summarize a feature's scoring across seeds
#'
#' Answers "which molecular features are consistently important?" by
#' collapsing the per-seed rows from `scoreFeaturesWithinSeed()` down to one
#' row per feature.
#'
#' @param scored_top_features The tibble of scored top features with their contribution, rank, and rank score within each seed generated from `scoreFeaturesWithinSeed()`
#'
#' @returns a tibble with one row per species/drug_label/drug_or_class/feature_type/feature_subtype/variable, with:
#' \itemize{
#'   \item \code{seed_ratio}: the number of seeds the feature appears in, divided by the total number of distinct \code{seed} values present anywhere in \code{scored_top_features} (a single count computed once for the whole call, not per feature or per group — so this assumes every group was fit with the same set of seeds).
#'   \item \code{mean_rank_score}, \code{median_rank_score}: mean/median of \code{rank_score} across seeds; ranges 0-1 with higher values indicating higher importance.
#'   \item \code{median_rank}: median of \code{rank} across seeds. (\code{mean_rank} is not currently computed.)
#'   \item \code{median_contribution}: median of \code{contribution} across seeds. (\code{mean_contribution} is not currently computed.)
#'   \item \code{median_cum_contrib}: median of \code{cum_contrib} across seeds.
#'   \item \code{best_rank}: the best (lowest) rank seen across seeds.
#'   \item \code{rank_consistent}: TRUE if the feature's rank is identical in every seed.
#'   \item \code{rank_score_sd}, \code{rank_score_cv}: standard deviation and coefficient of variation of \code{rank_score} across seeds, i.e. how much the normalized rank_score varies — computed on \code{rank_score} rather than raw rank so it is comparable across groups with different numbers of features.
#'   \item \code{in_core_consistent}, \code{in_core}: whether the feature's \code{in_core} flag is identical across every seed; \code{in_core} is that shared value if consistent, else \code{FALSE}.
#'   \item \code{sign_consistent}, \code{sign}: whether the feature's \code{sign} is identical across every seed; \code{sign} is that shared value if consistent, else \code{"MIXED"}.
#' }
#'
#' @keywords internal
#' @examples
#' summariseFeaturesAcrossSeeds(scoreFeaturesWithinSeed(all_top_features.parquet))
summariseFeaturesAcrossSeeds <- function(scored_top_features) {

  # find max number of seeds
  max_seeds <- scored_top_features |>
    dplyr::summarise(n_seeds = dplyr::n_distinct(seed)) |>
    dplyr::pull(n_seeds)

  feature_summary <- scored_top_features |>
    dplyr::group_by(
      species,
      drug_label,
      drug_or_class,
      feature_type,
      feature_subtype,
      variable
    ) |>
    dplyr::summarise(
      seed_ratio = dplyr::n_distinct(seed) / max_seeds,

      mean_rank_score = mean(rank_score, na.rm = TRUE),
      median_rank_score = median(rank_score, na.rm = TRUE),
      rank_score_sd = sd(rank_score, na.rm = TRUE),
      # coefficient of variation: how large rank_score_sd is relative to mean_rank_score
      rank_score_cv = dplyr::if_else(
        mean_rank_score != 0,
        rank_score_sd / mean_rank_score,
        NA_real_
      ),

      median_rank = median(rank, na.rm = TRUE),
      best_rank = min(rank, na.rm = TRUE),
      rank_consistent = dplyr::n_distinct(rank) == 1,

      median_contribution = median(contribution, na.rm = TRUE),
      median_cum_contrib = median(cum_contrib, na.rm = TRUE),

      in_core_consistent = dplyr::n_distinct(in_core) == 1,
      in_core = if (in_core_consistent) dplyr::first(in_core) else FALSE,

      sign_consistent = dplyr::n_distinct(sign) == 1,
      sign = if (sign_consistent) dplyr::first(sign) else "MIXED",
      .groups = "drop"
    )

  return(feature_summary)
}

#' Build a per-drug top-feature table with cutoffs
#'
#' Runs `scoreFeaturesWithinSeed()` and `summariseFeaturesAcrossSeeds()` on
#' `all_top_features_parquet`, then filters the resulting per-feature summary
#' down to the top features for each drug/class.
#'
#' @inheritParams scoreFeaturesWithinSeed
#' @param rank_score_quantile A value in \[0, 1\] (default is \code{0.95}). Keep only features whose \code{median_rank_score} is at or above this quantile of \code{median_rank_score}, where the quantile is computed separately within each (species, drug_label, drug_or_class) group (with \code{na.rm = TRUE}) — so a feature is compared only against other features from the same drug/class, not against the whole panel. This means a drug/class with weaker or more diffuse signal overall can still contribute its own top features, rather than being crowded out by drugs/classes with stronger or more concentrated signal. As with the other conditions listed below, this is not restricted to rows that already pass them: \code{dplyr::filter()} evaluates every condition passed to a single call against the same original data, so this threshold does not narrow as other conditions are applied
#' @param cv_threshold The maximum allowed coefficient of variation (default is \code{1}). Keep only features with \code{rank_score_cv <= cv_threshold}, i.e. drop features whose rank_score is inconsistent across seeds relative to its mean. A feature present in every seed (\code{seed_ratio == 1}) always passes this check regardless of its \code{rank_score_cv} — including when \code{rank_score_cv} is \code{NA}, which happens whenever a group has only a single seed
#' @param cumulative_contribution_threshold The cumulative-contribution cutoff, in \[0, 1\] (default is \code{0.75}, i.e. 75%). Keep only features with \code{median_cum_contrib <= cumulative_contribution_threshold}
#' @param seed_ratio_threshold If not \code{NULL} (the default), keep only features whose \code{seed_ratio} exactly equals this value
#' @param found_in_both_subtypes Logical indicating whether to additionally restrict to features that survive the filters above in both the binary and counts \code{feature_subtype} (default is \code{FALSE})
#' @param compare_median_to_sd_rank_score Logical indicating whether to additionally require \code{median_rank_score > rank_score_sd} (default is \code{FALSE})
#'
#' @returns a tibble of top features for each drug/class: the `summariseFeaturesAcrossSeeds()` output (species, drug label, drug or class, feature type, feature subtype, variable, seed_ratio, mean/median rank score, median rank, median contribution, median cumulative contribution, best rank, rank/sign/in_core consistency flags, sign), filtered in two \code{dplyr::filter()} passes.
#' The first pass keeps rows where all of the following hold, evaluated together against \code{feature_summary} (see \code{rank_score_quantile} for how the last condition is grouped):
#' \itemize{
#'   \item \code{seed_ratio == seed_ratio_threshold}, only applied when \code{seed_ratio_threshold} is not \code{NULL},
#'   \item \code{seed_ratio == 1} OR \code{rank_score_cv <= cv_threshold} (see \code{cv_threshold}),
#'   \item \code{in_core} is TRUE,
#'   \item \code{sign_consistent} is TRUE (sign is the same in every seed; this does not require the sign to be negative),
#'   \item \code{median_cum_contrib <= cumulative_contribution_threshold}, and
#'   \item \code{median_rank_score} is at or above the \code{rank_score_quantile} quantile of \code{median_rank_score}, computed separately per (species, drug_label, drug_or_class) group.
#' }
#' Two columns are then added, grouped by (species, drug_label, drug_or_class, feature_type, variable): \code{n_subtype} and \code{subtype_csv}, recording how many/which \code{feature_subtype} values each combination has among the rows that survived the first pass.
#' A second \code{dplyr::filter()} pass then optionally keeps only rows where \code{subtype_csv == "binary,counts"} (when \code{found_in_both_subtypes = TRUE}) and/or \code{median_rank_score > rank_score_sd} (when \code{compare_median_to_sd_rank_score = TRUE}).
#' Every drug/class may not have variables from all feature types. Because the \code{rank_score_quantile} cutoff is now per drug/class rather than global, every drug/class with at least one feature surviving the other conditions will generally contribute something here, rather than potentially being excluded entirely by comparison to stronger-signal drugs/classes elsewhere in the panel.
#'
#' @export
topFeaturesPerDrugOrClass <- function(
  all_top_features_parquet,
  core_contribution_threshold = 0.75,
  exclude_feature_types = NULL,
  filter_model = TRUE,
  all_performance_parquet,
  MCC_threshold = NULL,
  compare_to_shuffled = TRUE,
  rank_score_quantile = 0.95,
  cv_threshold = 1,
  cumulative_contribution_threshold = 0.75,
  # additional filters
  seed_ratio_threshold = NULL,
  found_in_both_subtypes = FALSE,
  compare_median_to_sd_rank_score = FALSE
  )
  {

  scored_features <- scoreFeaturesWithinSeed(
    all_top_features_parquet,
    core_contribution_threshold = core_contribution_threshold,
    exclude_feature_types = exclude_feature_types,
    filter_model = filter_model,
    all_performance_parquet,
    MCC_threshold = MCC_threshold,
    compare_to_shuffled = compare_to_shuffled
  )

  feature_summary <- summariseFeaturesAcrossSeeds(scored_features)

  top_filtered_features <- feature_summary |>
    dplyr::group_by(species, drug_label, drug_or_class) |>
    dplyr::mutate(
      rank_score_cutoff = quantile(median_rank_score, rank_score_quantile, na.rm = TRUE)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(
      if (!is.null(seed_ratio_threshold)) seed_ratio == seed_ratio_threshold else TRUE,
      seed_ratio == 1 | rank_score_cv <= cv_threshold, # rank_score_cv can be NA if there is only one seed.
      in_core,
      sign_consistent,
      median_cum_contrib <= cumulative_contribution_threshold,
      median_rank_score >= rank_score_cutoff
    ) |>
    dplyr::select(-rank_score_cutoff) |>
    dplyr::group_by(species, drug_label, drug_or_class, feature_type, variable) |>
    dplyr::mutate(
      n_subtype = dplyr::n_distinct(feature_subtype),
      subtype_csv = paste(sort(unique(feature_subtype)), collapse = ",")
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(
      if (found_in_both_subtypes) (subtype_csv == "binary,counts") else TRUE,
      if (compare_median_to_sd_rank_score) (median_rank_score > rank_score_sd) else TRUE
    )

  return(top_filtered_features)
}

#' Aggregate mapped features to protein dyads
#'
#' Internal helper called on the output of `topFeaturesPerDrugOrClass()`.
#'
#' @param top_filtered_features The tibble of top features for each drug/class generated from `topFeaturesPerDrugOrClass()`.
#' @param dyad_feature_parquet The path to the Parquet file containing the mapping of features to protein dyads. Must have a \code{feature} column formatted as \code{"<feature_type>:<variable>"} (using the short feature-type codes \code{amr}/\code{cog}/\code{defense}/\code{pfam}/\code{protein}) and a \code{dyad} column giving the dyad id.
#'
#' @returns a tibble with one row per species/drug_label/drug_or_class/dyad (\code{dyad} is the protein-dyad id), with:
#' \itemize{
#'   \item \code{frequency}: the number of top-feature rows mapped to this dyad.
#'   \item \code{n_variables}, \code{variables_csv}: number of, and comma-separated list of, distinct \code{feature} (\code{"<feature_type>:<variable>"}) values mapped to this dyad.
#'   \item \code{n_feature_types}, \code{feature_types_csv}: number of, and comma-separated list of, distinct \code{feature_type_subtype} (\code{"<feature_type>:<feature_subtype>"}) values mapped to this dyad.
#'   \item \code{dyad_median_rank_score}: median of \code{median_rank_score} across the top-feature rows mapped to this dyad.
#'   \item \code{dyad_median_contribution}: median of \code{median_contribution} across the top-feature rows mapped to this dyad.
#'   \item \code{sign_consistent}, \code{sign}: whether \code{sign} is identical across every top-feature row mapped to this dyad; \code{sign} is that shared value if consistent, else \code{"MIXED"}.
#' }
#'
#' @keywords internal
summariseDyads <- function(top_filtered_features,
                            dyad_feature_parquet
                            ) {
  stopifnot(is.data.frame(top_filtered_features))
  stopifnot(file.exists(dyad_feature_parquet))

  dyad_feature <- arrow::read_parquet(normalizePath(dyad_feature_parquet)) |> 
    dplyr::distinct()

  top_dyads <- top_filtered_features |>
    dplyr::mutate(
      feature_type = dplyr::case_when(
        feature_type == "AMRFinder"  ~ "amr",
        feature_type == "COG"        ~ "cog",
        feature_type == "DefenseCas" ~ "defense",
        feature_type == "Pfam"       ~ "pfam",
        TRUE                         ~ feature_type
      )
    ) |>
    tidyr::unite("feature", feature_type, variable, sep = ":", remove = FALSE) |>
    tidyr::unite("feature_type_subtype", feature_type, feature_subtype, sep = ":", remove = FALSE) |>
    dplyr::left_join(dyad_feature, by = "feature", relationship = "many-to-many") |>
    dplyr::filter(!is.na(dyad)) |>
    dplyr::select(species, drug_label, drug_or_class, feature_type_subtype, dyad, feature, sign,
      median_rank_score, median_contribution, subtype_csv) |>
    dplyr::group_by(species, drug_label, drug_or_class, dyad) |>
    dplyr::summarise(
      frequency = dplyr::n(),
      n_variables = dplyr::n_distinct(feature),
      variables_csv = paste(sort(unique(feature)), collapse = ","),
      n_feature_types = dplyr::n_distinct(feature_type_subtype),
      feature_types_csv = paste(sort(unique(feature_type_subtype)), collapse = ","),
      dyad_median_rank_score = median(median_rank_score, na.rm = TRUE),
      dyad_median_contribution = median(median_contribution, na.rm = TRUE),
      sign_consistent = dplyr::n_distinct(sign) == 1,
      sign = if (sign_consistent) dplyr::first(sign) else "MIXED",
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(dyad_median_rank_score), dplyr::desc(n_feature_types))

  return(top_dyads)
}

#' Build a feature network from selected top features and top dyads
#'
#' @param top_features Output of \code{topFeaturesPerDrugOrClass()}.
#' @param top_dyads Output of \code{summariseDyads()}.
#' @param dyad_feature_parquet Path to the same Parquet file passed to \code{summariseDyads()}: a \code{feature} column formatted as \code{"<feature_type>:<variable>"} (using the short feature-type codes \code{amr}/\code{cog}/\code{defense}/\code{pfam}/\code{protein}) and a \code{dyad} column giving the dyad id. \code{feature} is reconstructed here from \code{top_features$feature_type}/\code{variable} the same way \code{summariseDyads()} builds it, so the feature-to-dyad edges use the same mapping as the dyad table itself.
#' @param protein_names_parquet Path to the Parquet file with dyad name annotations.
#'
#' @returns A list with \code{feature_table}, \code{dyad_table}, \code{nodes}, \code{edges}, and \code{graph}.
#' Node/edge weights are built from \code{median_rank_score} (features) and \code{dyad_median_rank_score} (dyads, via the \code{dyad_score} column, itself median-based per \code{summariseDyads()}), consistent with the seed-noise-robust selection made in \code{topFeaturesPerDrugOrClass()}.
#' @export
buildFeatureNetwork <- function(top_features,
                               top_dyads,
                               dyad_feature_parquet,
                               protein_names_parquet
                              ) {
  stopifnot(is.data.frame(top_features))
  stopifnot(is.data.frame(top_dyads))
  stopifnot(file.exists(dyad_feature_parquet))
  stopifnot(file.exists(protein_names_parquet))

  required_feature_cols <- c(
    "species", "drug_label", "drug_or_class",
    "feature_type", "variable",
    "median_rank_score"
  )
  required_dyad_cols <- c(
    "species", "drug_label", "drug_or_class",
    "dyad", "dyad_median_rank_score"
  )

  missing_feature_cols <- setdiff(required_feature_cols, names(top_features))
  missing_dyad_cols <- setdiff(required_dyad_cols, names(top_dyads))

  if (length(missing_feature_cols) > 0) {
    stop("top_features is missing required columns: ",
         paste(missing_feature_cols, collapse = ", "))
  }
  if (length(missing_dyad_cols) > 0) {
    stop("top_dyads is missing required columns: ",
         paste(missing_dyad_cols, collapse = ", "))
  }

  dyad_feature <- arrow::read_parquet(normalizePath(dyad_feature_parquet)) |>
    dplyr::distinct()

  if (!all(c("dyad", "feature") %in% names(dyad_feature))) {
    stop(
      "dyad_feature_parquet is expected to have 'dyad' (dyad id) and ",
      "'feature' ('<feature_type>:<variable>') columns, matching what ",
      "summariseDyads() expects; found: ", paste(names(dyad_feature), collapse = ", ")
    )
  }

  protein_names <- arrow::read_parquet(normalizePath(protein_names_parquet)) |>
    dplyr::distinct()

  make_model_id <- function(drug_label, drug_or_class) {
    paste(drug_label, drug_or_class, sep = ".")
  }

  # Same short feature-type codes used by summariseDyads() to build `feature`.
  shorten_feature_type <- function(feature_type) {
    dplyr::case_when(
      feature_type == "AMRFinder"  ~ "amr",
      feature_type == "COG"        ~ "cog",
      feature_type == "DefenseCas" ~ "defense",
      feature_type == "Pfam"       ~ "pfam",
      TRUE                         ~ feature_type
    )
  }

  feature_table <- top_features |>
    dplyr::mutate(model_id = make_model_id(drug_label, drug_or_class)) |>
    dplyr::group_by(species, model_id, feature_type, variable) |>
    dplyr::summarise(
      # mean of median_rank_score across subtype (bin/count) rows for this variable --
      # this is where bin/count reconciliation currently happens (implicitly)
      feature_score = mean(median_rank_score, na.rm = TRUE),
      .groups = "drop"
    )

  dyad_table <- top_dyads |>
    dplyr::mutate(model_id = make_model_id(drug_label, drug_or_class)) |>
    dplyr::group_by(species, model_id, dyad) |>
    dplyr::summarise(
      dyad_score = mean(dyad_median_rank_score, na.rm = TRUE),
      .groups = "drop"
    )

  model_nodes <- dplyr::bind_rows(
    feature_table |>
      dplyr::distinct(species, model_id),
    dyad_table |>
      dplyr::distinct(species, model_id)
  ) |>
    dplyr::distinct(species, model_id) |>
    dplyr::transmute(
      name = model_id,
      label = model_id,
      node_type = "model",
      species = species,
      score = NA_real_,
      breadth = NA_real_,
      node_size = 4
    )

  feature_nodes <- feature_table |>
    dplyr::group_by(species, variable) |>
    dplyr::summarise(
      score = mean(feature_score, na.rm = TRUE),
      breadth = dplyr::n_distinct(model_id),
      .groups = "drop"
    ) |>
    dplyr::transmute(
      name = variable,
      label = variable,
      node_type = "feature",
      species = species,
      score = score,
      breadth = breadth,
      node_size = pmax(3, pmin(10, breadth + 2))
    )

  dyad_nodes <- dyad_table |>
    dplyr::group_by(species, dyad) |>
    dplyr::summarise(
      score = median(dyad_score, na.rm = TRUE),
      breadth = dplyr::n_distinct(model_id),
      .groups = "drop"
    ) |>
    dplyr::transmute(
      name = dyad,
      label = dyad,
      node_type = "dyad",
      species = species,
      score = score,
      breadth = breadth,
      node_size = pmax(3, pmin(10, breadth + 2))
    )

  nodes <- dplyr::bind_rows(model_nodes, feature_nodes, dyad_nodes) |>
    dplyr::distinct(name, .keep_all = TRUE)

  feature_edges <- feature_table |>
    dplyr::transmute(
      from = model_id,
      to = variable,
      weight = feature_score,
      edge_type = "model_feature"
    ) |>
    dplyr::distinct(from, to, edge_type, .keep_all = TRUE)

  dyad_edges <- dyad_table |>
    dplyr::transmute(
      from = model_id,
      to = dyad,
      weight = dyad_score,
      edge_type = "model_dyad"
    ) |>
    dplyr::distinct(from, to, edge_type, .keep_all = TRUE)

  feature_dyad_edges <- feature_table |>
    dplyr::mutate(short_feature_type = shorten_feature_type(feature_type)) |>
    tidyr::unite("feature", short_feature_type, variable, sep = ":", remove = FALSE) |>
    dplyr::left_join(
      dyad_feature |> dplyr::add_count(feature, name = "n_dyads"),
      by = "feature",
      relationship = "many-to-many"
    ) |>
    dplyr::filter(!is.na(dyad)) |>
    dplyr::transmute(
      from = variable,
      to = dyad,
      weight = 1 / n_dyads,
      edge_type = "feature_dyad"
    ) |>
    dplyr::distinct(from, to, edge_type, .keep_all = TRUE)

  edges <- dplyr::bind_rows(feature_edges, feature_dyad_edges, dyad_edges)

  missing_vertices <- setdiff(unique(c(edges$from, edges$to)), nodes$name)
  if (length(missing_vertices) > 0) {
    extra_nodes <- tibble::tibble(name = missing_vertices) |>
      dplyr::mutate(
        label = name,
        node_type = dplyr::case_when(
          grepl("^drug\\.|^drug_class\\.", name) ~ "model",
          grepl("^fig\\||^dyad", name) ~ "dyad",
          TRUE ~ "feature"
        ),
        species = NA_character_,
        score = NA_real_,
        breadth = NA_real_,
        node_size = 4
      )

    nodes <- dplyr::bind_rows(nodes, extra_nodes) |>
      dplyr::distinct(name, .keep_all = TRUE)
  }

  graph <- if (nrow(edges) > 0) {
    igraph::graph_from_data_frame(
      d = edges,
      directed = FALSE,
      vertices = nodes
    )
  } else {
    igraph::graph_from_data_frame(
      d = data.frame(from = character(), to = character()),
      directed = FALSE,
      vertices = nodes
    )
  }

  feature_network <- list(
    feature_table = feature_table,
    dyad_table = dyad_table,
    nodes = nodes,
    edges = edges,
    graph = graph
  )

  return(feature_network)
}

#' Plot the feature network with networkD3
#'
#' @param feature_network Output of \code{buildFeatureNetwork()}.
#' @param height Widget height in pixels (default is \code{800}).
#' @param width Widget width (default is \code{"100\%"}).
#'
#' @returns A \code{networkD3} widget.
#' @export
plotFeatureNetworkD3 <- function(feature_network,
                                height = 800,
                                width = "100%"
                              ) {

  stopifnot(is.list(feature_network))
  stopifnot(!is.null(feature_network$nodes))
  stopifnot(!is.null(feature_network$edges))

  nodes <- feature_network$nodes |>
    dplyr::distinct(name, .keep_all = TRUE) |>
    dplyr::mutate(
      id = dplyr::row_number() - 1L,
      group = node_type,
      title = paste0(
        "<b>", label, "</b>",
        ifelse(is.na(species), "", paste0("<br>Species: ", species)),
        ifelse(is.na(score), "", paste0("<br>Score: ", signif(score, 3))),
        ifelse(is.na(breadth), "", paste0("<br>Breadth: ", breadth))
      )
    )

  links <- feature_network$edges |>
    dplyr::filter(!is.na(from), !is.na(to)) |>
    dplyr::left_join(
      nodes |> dplyr::select(name, id),
      by = c("from" = "name")
    ) |>
    dplyr::rename(dyad = id) |>
    dplyr::left_join(
      nodes |> dplyr::select(name, id),
      by = c("to" = "name")
    ) |>
    dplyr::rename(feature = id) |>
    dplyr::filter(!is.na(dyad), !is.na(feature)) |>
    dplyr::mutate(
      value = dplyr::if_else(is.na(weight), 1, weight)
    ) |>
    dplyr::select(dyad, feature, value, edge_type)

  stopifnot(nrow(nodes) > 0)
  stopifnot(nrow(links) > 0)

  colour_scale <- networkD3::JS(
    "d3.scaleOrdinal()
      .domain(['model', 'feature', 'dyad'])
      .range(['#4C78A8', '#F58518', '#54A24B'])"
  )

  networkD3::forceNetwork(
    Links = links,
    Nodes = nodes,
    Source = "dyad",
    Target = "feature",
    Value = "value",
    NodeID = "label",
    Group = "group",
    opacity = 0.9,
    zoom = TRUE,
    fontSize = 14,
    height = height,
    width = width,
    colourScale = colour_scale,
    linkDistance = networkD3::JS(
      "function(d) {
         if (d.edge_type === 'feature_dyad') return 60;
         if (d.edge_type === 'model_feature') return 120;
         return 90;
       }"
    )
  )
}

#' Find dyads that appear across multiple drugs/classes
#'
#' @param top_dyads The tibble of summarized dyads generated from `summarisedyads()`
#' @param label The \code{drug_label} value to filter dyads by, either \code{"drug"} or \code{"drug_class"} (default is \code{"drug"})
#' @param min_drugs_or_classes The minimum number of distinct drugs or classes required for a dyad to be considered shared (default is \code{2})
#'
#' @returns a tibble with one row per shared \code{dyad}, with \code{n_drug_or_class} (the number of distinct \code{drug_or_class} values the dyad appears in) and \code{drug_or_class_csv} (a comma-separated string of those values), sorted by \code{n_drug_or_class} descending.
#'
#' @export
findSharedDyads <- function(top_dyads = summariseDyads(top_features, dyad_feature_parquet),
                                label = "drug",
                                 min_drugs_or_classes = 2
                                ) {
  shared_dyads <- top_dyads |>
    dplyr::filter(!is.na(dyad), drug_label == label) |>
    dplyr::group_by(dyad) |>
    dplyr::mutate(
      n_drug_or_class = dplyr::n_distinct(drug_or_class),
      drug_or_class_csv = paste(sort(unique(drug_or_class)), collapse = ", ")
    ) |>
    dplyr::filter(n_drug_or_class >= min_drugs_or_classes) |>
    dplyr::ungroup() |>
    dplyr::select(dyad, n_drug_or_class, drug_or_class_csv) |>
    dplyr::arrange(dplyr::desc(n_drug_or_class))

  return(shared_dyads)
}

#' Find the dyads that are unique to a single drug/class
#'
#' @param top_dyads The tibble of summarized dyads generated from `summarisedyads()`
#' @param label The \code{drug_label} value to filter dyads by, either \code{"drug"} or \code{"drug_class"} (default is \code{"drug"})
#' @param protein_names_parquet The path to the Parquet file containing the annotations to protein dyad names
#'
#' @returns a tibble with one row per dyad unique to a single drug/class, with \code{drug_or_class}, \code{dyad}, \code{dyad_name} (from the protein name annotations), and \code{dyad_mean_rank_score}, sorted by \code{dyad_mean_rank_score} descending.
#'
#' @export
#' @examples
#' findUniquedyads(summarisedyads(top_features, dyad_feature_parquet), label = "drug", protein_names_parquet)
findUniqueDyads <- function(top_dyads = summariseDyads(top_features, dyad_feature_parquet),
                            label = "drug",
                            protein_names_parquet
) {

  protein_names <- arrow::read_parquet(normalizePath(protein_names_parquet)) |>
    dplyr::distinct()

  unique_dyads <- top_dyads |>
    dplyr::filter(!is.na(dyad), drug_label == label) |>
    dplyr::group_by(dyad) |>
    dplyr::mutate(
      n_drug_or_class = dplyr::n_distinct(drug_or_class),
      drug_or_class_csv = paste(sort(unique(drug_or_class)), collapse = ", ")
    ) |>
    dplyr::filter(n_drug_or_class == 1) |>
    dplyr::ungroup() |>
    dplyr::select(dyad, drug_or_class_csv, dyad_median_rank_score) |>
    dplyr::arrange(dplyr::desc(dyad_median_rank_score)) |>
dplyr::rename(drug_or_class = drug_or_class_csv) 

  return(unique_dyads)
}

# final run would be:
# top_features <- topFeaturesPerDrugOrClass(rank_score_quantile = 0.75)
# top_dyads <- summarisedyads(top_features, dyad_feature_parquet = dyad_feature_parquet)
# feature_network <- buildFeatureNetwork(top_features = top_features, top_dyads = top_dyads,
#   dyad_feature_parquet = dyad_feature_parquet, protein_names_parquet = protein_names_parquet)
  # plotFeatureNetworkD3(feature_network)

#' Discover high-performing models, stable features, and protein dyads
#'
#' Runs the model-quality, feature-ranking, and protein-dyad aggregation
#' workflow and organizes the results into question-oriented tables.
#'
#' The returned object can be used to identify:
#' \itemize{
#'   \item model fits with good non-shuffled performance;
#'   \item model fits with good separation from shuffled-label models;
#'   \item model fits satisfying both performance conditions;
#'   \item top features with high seed coverage and median rank score;
#'   \item top features with high median contribution;
#'   \item top features present across every analyzed seed;
#'   \item top features shared across drugs or drug classes;
#'   \item top features unique to one drug or drug class;
#'   \item high-scoring, shared, and unique protein dyads.
#' }
#'
#' This is a high-level wrapper around \code{filterOptimalModel()},
#' \code{topFeaturesPerDrugOrClass()}, and \code{summariseDyads()}.
#'
#' @param all_top_features_parquet Path to the Parquet file containing
#' top features and their importance scores.
#' @param all_performance_parquet Path to the Parquet file containing
#' model-performance results.
#' @param dyad_feature_parquet Path to the Parquet file mapping features
#' to protein dyads.
#'
#' @param MCC_threshold Minimum non-shuffled MCC used to define good model
#' performance. Default is \code{0.5}, meaning that only models with a non-shuffled MCC of 0.5 or higher are considered.
#' @param compare_to_shuffled Logical indicating whether the models used for
#' feature selection must outperform their shuffled-label counterparts.
#' Default is \code{TRUE}.
#'
#' @param core_contribution_threshold Cumulative-contribution threshold
#' passed to \code{scoreFeaturesWithinSeed()}. Default is \code{0.9}.
#' @param exclude_feature_types Feature types to exclude before feature
#' scoring. Default is \code{NULL}.
#' @param filter_model Logical indicating whether feature scoring should be
#' restricted to model groups passing the requested model-quality criteria.
#' Default is \code{TRUE}.
#' @param rank_score_quantile Quantile cutoff for median rank score, computed
#' separately within each species/drug-label/drug-or-class group.
#' Default is \code{0.95}.
#' @param cv_threshold Maximum allowed rank-score coefficient of variation.
#' Default is \code{1}.
#' @param cumulative_contribution_threshold Maximum allowed median cumulative
#' contribution. Default is \code{0.75}.
#' @param seed_ratio_threshold Optional exact seed-ratio requirement passed
#' to \code{topFeaturesPerDrugOrClass()}. Default is \code{1}.
#' @param found_in_both_subtypes Logical indicating whether selected features
#' must occur in both binary and counts feature subtypes.
#' Default is \code{FALSE}.
#' @param compare_median_to_sd_rank_score Logical indicating whether selected
#' features must have \code{median_rank_score > rank_score_sd}.
#' Default is \code{FALSE}.
#'
#' @param consistent_seed_ratio Minimum seed ratio used to define a
#' consistently selected feature. Default is \code{0.8}.
#' @param consistent_rank_score Minimum median rank score used to define a
#' consistently highly ranked feature. Default is \code{0.75}.
#' @param high_contribution_quantile Quantile of \code{median_contribution}
#' used to define high-contribution features. It is computed separately
#' within each species/drug-label/drug-or-class group. Default is
#' \code{0.75}.
#' @param high_dyad_score_quantile Quantile of
#' \code{dyad_median_rank_score} used to define high-scoring dyads. It is
#' computed separately within each species/drug-label/drug-or-class group.
#' Default is \code{0.75}.
#' @param min_models_shared Minimum number of distinct drugs or drug classes
#' required for a feature or dyad to be considered shared.
#' Default is \code{2}.
#'
#' @returns An object of class \code{feature_dyad_discovery}, implemented as
#' a named list containing:
#' \itemize{
#'   \item \code{models$good_performance}: fit-level models satisfying
#'   \code{MCC_threshold}, without requiring shuffled-model separation.
#'   \item \code{models$good_shuffled_separation}: fit-level models whose
#'   non-shuffled MCC exceeds the shuffled-label MCC.
#'   \item \code{models$qualified}: fit-level models satisfying the model
#'   criteria used for feature selection.
#'   \item \code{models$summary}: model-level summary across seeds and feature
#'   scales.
#'   \item \code{features$top}: all selected top features.
#'   \item \code{features$consistent}: selected features with high seed ratio
#'   and high median rank score.
#'   \item \code{features$high_contribution}: selected features with high
#'   median contribution.
#'   \item \code{features$all_seeds}: selected features present across all
#'   analyzed seeds.
#'   \item \code{features$shared}: selected features associated with at least
#'   \code{min_models_shared} drugs or drug classes.
#'   \item \code{features$unique}: selected features associated with exactly
#'   one drug or drug class.
#'   \item \code{dyads$top}: dyads mapped from the selected top features.
#'   \item \code{dyads$high_score}: dyads with high median rank scores.
#'   \item \code{dyads$shared}: dyads associated with at least
#'   \code{min_models_shared} drugs or drug classes.
#'   \item \code{dyads$unique}: dyads associated with exactly one drug or
#'   drug class.
#'   \item \code{parameters}: parameter values used for the analysis.
#' }
#'
#' @export
runFeatureDyadDiscovery <- function(
    all_top_features_parquet,
    all_performance_parquet,
    dyad_feature_parquet,
    MCC_threshold = 0.5,
    compare_to_shuffled = TRUE,
    core_contribution_threshold = 0.9,
    exclude_feature_types = NULL,
    filter_model = TRUE,
    rank_score_quantile = 0.95,
    cv_threshold = 1,
    cumulative_contribution_threshold = 0.75,
    seed_ratio_threshold = 1,
    found_in_both_subtypes = FALSE,
    compare_median_to_sd_rank_score = FALSE,
    consistent_seed_ratio = 0.8,
    consistent_rank_score = 0.75,
    high_contribution_quantile = 0.75,
    high_dyad_score_quantile = 0.75,
    min_models_shared = 2
) {

  # -------------------------------------------------------------------------
  # Validate files
  # -------------------------------------------------------------------------

  input_files <- c(
    all_top_features_parquet = all_top_features_parquet,
    all_performance_parquet = all_performance_parquet,
    dyad_feature_parquet = dyad_feature_parquet
  )

  missing_files <- input_files[!file.exists(input_files)]

  if (length(missing_files) > 0) {
    stop(
      "The following input file(s) do not exist: ",
      paste(names(missing_files), collapse = ", ")
    )
  }

  # -------------------------------------------------------------------------
  # Validate analysis thresholds
  # -------------------------------------------------------------------------

  unit_interval_parameters <- c(
    core_contribution_threshold = core_contribution_threshold,
    rank_score_quantile = rank_score_quantile,
    cumulative_contribution_threshold =
      cumulative_contribution_threshold,
    consistent_seed_ratio = consistent_seed_ratio,
    consistent_rank_score = consistent_rank_score,
    high_contribution_quantile = high_contribution_quantile,
    high_dyad_score_quantile = high_dyad_score_quantile
  )

  invalid_parameters <- names(unit_interval_parameters)[
    is.na(unit_interval_parameters) |
      unit_interval_parameters < 0 |
      unit_interval_parameters > 1
  ]

  if (length(invalid_parameters) > 0) {
    stop(
      "The following parameter(s) must be between 0 and 1: ",
      paste(invalid_parameters, collapse = ", ")
    )
  }

  if (!is.null(seed_ratio_threshold)) {
    if (
      length(seed_ratio_threshold) != 1 ||
      is.na(seed_ratio_threshold) ||
      seed_ratio_threshold < 0 ||
      seed_ratio_threshold > 1
    ) {
      stop(
        "seed_ratio_threshold must be NULL or a single value ",
        "between 0 and 1."
      )
    }
  }

  if (
    length(cv_threshold) != 1 ||
    is.na(cv_threshold) ||
    cv_threshold < 0
  ) {
    stop("cv_threshold must be a single non-negative value.")
  }

  if (
    length(min_models_shared) != 1 ||
    is.na(min_models_shared) ||
    min_models_shared < 2 ||
    min_models_shared != as.integer(min_models_shared)
  ) {
    stop(
      "min_models_shared must be a single integer greater than ",
      "or equal to 2."
    )
  }

  # -------------------------------------------------------------------------
  # 1. Find fits with good absolute performance
  #
  # This table answers:
  # "Which model fits have sufficiently high non-shuffled MCC?"
  # -------------------------------------------------------------------------

  good_performance_models <- filterOptimalModel(
    all_performance_parquet = all_performance_parquet,
    MCC_threshold = MCC_threshold,
    compare_to_shuffled = FALSE
  )

  # -------------------------------------------------------------------------
  # 2. Find fits with good separation from shuffled models
  #
  # This table answers:
  # "Which model fits outperform their shuffled-label counterparts?"
  #
  # MCC_threshold is deliberately NULL here so that separation can be
  # examined independently of absolute MCC.
  # -------------------------------------------------------------------------

  good_shuffled_separation_models <- filterOptimalModel(
    all_performance_parquet = all_performance_parquet,
    MCC_threshold = NULL,
    compare_to_shuffled = TRUE
  )

  # -------------------------------------------------------------------------
  # 3. Find models satisfying the criteria used for feature selection
  #
  # If compare_to_shuffled is TRUE, these satisfy both the requested
  # MCC threshold and shuffled-label comparison.
  # -------------------------------------------------------------------------

  qualified_models <- filterOptimalModel(
    all_performance_parquet = all_performance_parquet,
    MCC_threshold = MCC_threshold,
    compare_to_shuffled = compare_to_shuffled
  )

  if (nrow(qualified_models) == 0) {
    stop(
      "No model fits passed the requested model-quality criteria."
    )
  }

  # Summarize fit-level results to the species/drug or drug-class level.
  model_summary <- qualified_models |>
    dplyr::group_by(
      species,
      drug_label,
      drug_or_class
    ) |>
    dplyr::summarise(
      n_passing_fits = dplyr::n(),
      n_seeds = dplyr::n_distinct(seed),
      n_feature_types = dplyr::n_distinct(feature_type),
      n_feature_scales = dplyr::n_distinct(
        paste(feature_type, feature_subtype, sep = ":")
      ),
      median_nonshuffled_MCC = median(
        nonshuffled_MCC,
        na.rm = TRUE
      ),
      best_nonshuffled_MCC = max(
        nonshuffled_MCC,
        na.rm = TRUE
      ),
      median_shuffled_MCC = median(
        shuffled_MCC,
        na.rm = TRUE
      ),
      median_MCC_diff = median(
        MCC_diff,
        na.rm = TRUE
      ),
      if(!is.na(median_MCC_diff)) {
        minimum_MCC_diff = min(
          MCC_diff,
          na.rm = TRUE
        )
      } else {
        minimum_MCC_diff = NULL
      },
      if(!is.na(median_MCC_diff)) {
        maximum_MCC_diff = max(
          MCC_diff,
          na.rm = TRUE
        )
      } else {
        maximum_MCC_diff = NULL
      },
      proportion_outperforming_shuffled = median(
        MCC_diff > 0 | is.na(shuffled_MCC),
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      dplyr::desc(median_nonshuffled_MCC),
      dplyr::desc(median_MCC_diff)
    )

  # -------------------------------------------------------------------------
  # 4. Select top features
  # -------------------------------------------------------------------------

  top_features <- topFeaturesPerDrugOrClass(
    all_top_features_parquet = all_top_features_parquet,
    core_contribution_threshold = core_contribution_threshold,
    exclude_feature_types = exclude_feature_types,
    filter_model = filter_model,
    all_performance_parquet = all_performance_parquet,
    MCC_threshold = MCC_threshold,
    compare_to_shuffled = compare_to_shuffled,
    rank_score_quantile = rank_score_quantile,
    cv_threshold = cv_threshold,
    cumulative_contribution_threshold =
      cumulative_contribution_threshold,
    seed_ratio_threshold = seed_ratio_threshold,
    found_in_both_subtypes = found_in_both_subtypes,
    compare_median_to_sd_rank_score =
      compare_median_to_sd_rank_score
  )

  if (nrow(top_features) == 0) {
    stop(
      "No features passed the requested feature-selection criteria."
    )
  }

  # -------------------------------------------------------------------------
  # 5. Consistently top-ranked features
  #
  # Answers:
  # "Which features have both high seed coverage and high median rank?"
  # -------------------------------------------------------------------------

  consistent_features <- top_features |>
    dplyr::filter(
      seed_ratio >= consistent_seed_ratio,
      median_rank_score >= consistent_rank_score
    ) |>
    dplyr::arrange(
      species,
      drug_label,
      drug_or_class,
      dplyr::desc(seed_ratio),
      dplyr::desc(median_rank_score)
    )

  # -------------------------------------------------------------------------
  # 6. High-contribution features
  #
  # The contribution cutoff is calculated separately for each species and
  # drug/drug-class model so each model is evaluated relative to its own
  # selected feature set.
  # -------------------------------------------------------------------------

  high_contribution_features <- top_features |>
    dplyr::group_by(
      species,
      drug_label,
      drug_or_class
    ) |>
    dplyr::mutate(
      contribution_cutoff = stats::quantile(
        median_contribution,
        probs = high_contribution_quantile,
        na.rm = TRUE
      )
    ) |>
    dplyr::filter(
      median_contribution >= contribution_cutoff
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-contribution_cutoff) |>
    dplyr::arrange(
      species,
      drug_label,
      drug_or_class,
      dplyr::desc(median_contribution)
    )

  # -------------------------------------------------------------------------
  # 7. Features present across all seeds
  #
  # near() avoids problems caused by floating-point representation.
  # -------------------------------------------------------------------------

  all_seed_features <- top_features |>
    dplyr::filter(dplyr::near(seed_ratio, 1)) |>
    dplyr::arrange(
      species,
      drug_label,
      drug_or_class,
      dplyr::desc(median_rank_score),
      dplyr::desc(median_contribution)
    )

  # -------------------------------------------------------------------------
  # 8. Summarize sharedness of top features
  #
  # Sharedness is evaluated independently for drugs and drug classes because
  # drug_label is part of the grouping.
  #
  # Feature subtype is not part of the feature identity here. Thus, binary
  # and counts representations of the same variable are treated as the same
  # molecular feature.
  # -------------------------------------------------------------------------

  feature_sharedness <- top_features |>
    dplyr::group_by(
      species,
      drug_label,
      feature_type,
      variable
    ) |>
    dplyr::summarise(
      n_drug_or_class = dplyr::n_distinct(drug_or_class),
      drug_or_class_csv = paste(
        sort(unique(drug_or_class)),
        collapse = ", "
      ),
      maximum_seed_ratio = max(
        seed_ratio,
        na.rm = TRUE
      ),
      median_seed_ratio = median(
        seed_ratio,
        na.rm = TRUE
      ),
      maximum_rank_score = max(
        median_rank_score,
        na.rm = TRUE
      ),
      overall_median_rank_score = median(
        median_rank_score,
        na.rm = TRUE
      ),
      overall_median_contribution = median(
        median_contribution,
        na.rm = TRUE
      ),
      sign_consistent_across_models =
        dplyr::n_distinct(sign) == 1,
      sign = if (sign_consistent_across_models) {
        dplyr::first(sign)
      } else {
        "MIXED"
      },
      .groups = "drop"
    )

  # -------------------------------------------------------------------------
  # 9. Features shared across drugs or drug classes
  # -------------------------------------------------------------------------

  shared_features <- feature_sharedness |>
    dplyr::filter(
      n_drug_or_class >= min_models_shared
    ) |>
    dplyr::arrange(
      species,
      drug_label,
      dplyr::desc(n_drug_or_class),
      dplyr::desc(overall_median_rank_score),
      dplyr::desc(overall_median_contribution)
    )

  # -------------------------------------------------------------------------
  # 10. Features unique to one drug or drug class
  #
  # Join back to top_features so the output retains the associated
  # drug_or_class and feature-level statistics.
  # -------------------------------------------------------------------------

  unique_feature_ids <- feature_sharedness |>
    dplyr::filter(n_drug_or_class == 1) |>
    dplyr::select(
      species,
      drug_label,
      feature_type,
      variable
    )

  unique_features <- top_features |>
    dplyr::semi_join(
      unique_feature_ids,
      by = c(
        "species",
        "drug_label",
        "feature_type",
        "variable"
      )
    ) |>
    dplyr::arrange(
      species,
      drug_label,
      drug_or_class,
      dplyr::desc(median_rank_score),
      dplyr::desc(median_contribution)
    )

  # -------------------------------------------------------------------------
  # 11. Aggregate top features to protein dyads
  # -------------------------------------------------------------------------

  top_dyads <- summariseDyads(
    top_filtered_features = top_features,
    dyad_feature_parquet = dyad_feature_parquet
  )

  # -------------------------------------------------------------------------
  # 12. High-scoring dyads
  #
  # A dyad score is high relative to the other selected dyads for the same
  # species and drug/drug-class model.
  # -------------------------------------------------------------------------

  high_score_dyads <- top_dyads |>
    dplyr::group_by(
      species,
      drug_label,
      drug_or_class
    ) |>
    dplyr::mutate(
      dyad_score_cutoff = stats::quantile(
        dyad_median_rank_score,
        probs = high_dyad_score_quantile,
        na.rm = TRUE
      )
    ) |>
    dplyr::filter(
      dyad_median_rank_score >= dyad_score_cutoff
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-dyad_score_cutoff) |>
    dplyr::arrange(
      species,
      drug_label,
      drug_or_class,
      dplyr::desc(dyad_median_rank_score),
      dplyr::desc(dyad_median_contribution)
    )

  # -------------------------------------------------------------------------
  # 13. Summarize dyad sharedness
  #
  # A dyad is identified independently within each species and drug-label
  # level.
  # -------------------------------------------------------------------------

  dyad_sharedness <- top_dyads |>
    dplyr::group_by(
      species,
      drug_label,
      dyad
    ) |>
    dplyr::summarise(
      n_drug_or_class = dplyr::n_distinct(drug_or_class),
      drug_or_class_csv = paste(
        sort(unique(drug_or_class)),
        collapse = ", "
      ),
      total_frequency = sum(
        frequency,
        na.rm = TRUE
      ),
      n_distinct_variables = dplyr::n_distinct(
        unlist(strsplit(variables_csv, ",", fixed = TRUE))
      ),
      maximum_dyad_rank_score = max(
        dyad_median_rank_score,
        na.rm = TRUE
      ),
      overall_dyad_median_rank_score = median(
        dyad_median_rank_score,
        na.rm = TRUE
      ),
      overall_dyad_median_contribution = median(
        dyad_median_contribution,
        na.rm = TRUE
      ),
      sign_consistent_across_models =
        dplyr::n_distinct(sign) == 1,
      sign = if (sign_consistent_across_models) {
        dplyr::first(sign)
      } else {
        "MIXED"
      },
      .groups = "drop"
    )

  # -------------------------------------------------------------------------
  # 14. Dyads shared across drugs or drug classes
  # -------------------------------------------------------------------------

  shared_dyads <- dyad_sharedness |>
    dplyr::filter(
      n_drug_or_class >= min_models_shared
    ) |>
    dplyr::arrange(
      species,
      drug_label,
      dplyr::desc(n_drug_or_class),
      dplyr::desc(overall_dyad_median_rank_score),
      dplyr::desc(total_frequency)
    )

  # -------------------------------------------------------------------------
  # 15. Dyads unique to one drug or drug class
  #
  # Join back to top_dyads to retain the drug/class identity and the
  # original dyad statistics.
  # -------------------------------------------------------------------------

  unique_dyad_ids <- dyad_sharedness |>
    dplyr::filter(n_drug_or_class == 1) |>
    dplyr::select(
      species,
      drug_label,
      dyad
    )

  unique_dyads <- top_dyads |>
    dplyr::semi_join(
      unique_dyad_ids,
      by = c(
        "species",
        "drug_label",
        "dyad"
      )
    ) |>
    dplyr::arrange(
      species,
      drug_label,
      drug_or_class,
      dplyr::desc(dyad_median_rank_score),
      dplyr::desc(dyad_median_contribution)
    )

  # -------------------------------------------------------------------------
  # 16. Return an organized analysis object
  # -------------------------------------------------------------------------

  result <- list(
    models = list(
      good_performance = good_performance_models,
      good_shuffled_separation =
        good_shuffled_separation_models,
      qualified = qualified_models,
      summary = model_summary
    ),
    features = list(
      top = top_features,
      consistent = consistent_features,
      high_contribution = high_contribution_features,
      all_seeds = all_seed_features,
      sharedness = feature_sharedness,
      shared = shared_features,
      unique = unique_features
    ),
    dyads = list(
      top = top_dyads,
      high_score = high_score_dyads,
      sharedness = dyad_sharedness,
      shared = shared_dyads,
      unique = unique_dyads
    ),
    parameters = list(
      MCC_threshold = MCC_threshold,
      compare_to_shuffled = compare_to_shuffled,
      core_contribution_threshold =
        core_contribution_threshold,
      exclude_feature_types = exclude_feature_types,
      filter_model = filter_model,
      rank_score_quantile = rank_score_quantile,
      cv_threshold = cv_threshold,
      cumulative_contribution_threshold =
        cumulative_contribution_threshold,
      seed_ratio_threshold = seed_ratio_threshold,
      found_in_both_subtypes = found_in_both_subtypes,
      compare_median_to_sd_rank_score =
        compare_median_to_sd_rank_score,
      consistent_seed_ratio = consistent_seed_ratio,
      consistent_rank_score = consistent_rank_score,
      high_contribution_quantile =
        high_contribution_quantile,
      high_dyad_score_quantile =
        high_dyad_score_quantile,
      min_models_shared = min_models_shared
    )
  )

  class(result) <- c(
    "feature_dyad_discovery",
    class(result)
  )

  return(result)
}