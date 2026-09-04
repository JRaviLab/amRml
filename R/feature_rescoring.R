
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
  dplyr::mutate(
    MCC_diff = nonshuffled_MCC - shuffled_MCC
  ) |>
  dplyr::filter(
  if (!is.null(MCC_threshold))(nonshuffled_MCC >= MCC_threshold) else TRUE,
  if (compare_to_shuffled) (MCC_diff > 0 | is.na(shuffled_MCC)) else TRUE
)
}

#' Score features within each seed
#'
#' @param all_top_features_parquet The path to the Parquet file containing all top features with their importance scores.
#' @param core_contribution_threshold The cumulative-contribution cutoff, in \[0, 1\] (default is \code{0.75}, i.e. 75%), used to flag whether a feature falls within the "core" set of features that jointly account for that share of a seed's total importance.
#' @param exclude_feature_types Feature types to drop before any scoring happens (default is \code{"struct"}). struct variables are composite IDs (e.g. \code{polA.group_211.group_2176}, three dot-joined gene/domain identifiers) representing a co-occurrence/structural motif rather than a single molecular entity like the other five scales, and its candidate-variable count (tens of thousands per group) dwarfs the other scales by orders of magnitude — pooling it into this rank_score/contribution machinery would compare a compound signal against five primary ones on an incomparable scale. struct is reserved for post-hoc biological annotation once top clusters are identified, not for scoring/ranking/thresholding here.
#' @param filter_model Logical indicating whether to restrict scoring to (species, drug_label, drug_or_class, feature_type, feature_subtype, seed) groups that have at least one model fit passing \code{filterOptimalModel()}'s MCC/shuffled-comparison quality thresholds (default is \code{TRUE}). The join is not keyed on \code{model}/\code{fit_penalty}/\code{fit_mixture}, so if \code{all_top_features_parquet} contains multiple fits per group, all of that group's rows survive as soon as any one fit passes.
#' @param all_performance_parquet The path to the all performance parquet file. Always required — it is read unconditionally (regardless of \code{filter_model} or \code{add_lasso_advtg}) to compute the per-fit sparsity score.
#' @param MCC_threshold The minimum MCC threshold passed through to \code{filterOptimalModel()} (default is \code{NULL}, no filtering)
#' @param compare_to_shuffled Logical indicating whether to compare the model to shuffled data, passed through to \code{filterOptimalModel()} (default is \code{TRUE})
#' @param add_lasso_advtg Logical indicating whether to weight each feature's contribution by a sparsity score that rewards model fits returning fewer features relative to the candidate feature space (default is \code{TRUE})
#'
#' @returns a tibble of scored top features, one row per feature within each species/drug_label/drug_or_class/feature_type/feature_subtype/seed group, with the following columns added:
#' \itemize{
#'   \item \code{contribution}: the feature's importance divided by the sum of importance across all features in the group.
#'   \item \code{feat_return_ratio}: \code{n_feats_returned / n_feat} for the fit that produced this group (from \code{all_performance_parquet}).
#'   \item \code{sparsity_score}: \code{1 - feat_return_ratio} when \code{add_lasso_advtg = TRUE} (fits returning fewer features relative to the candidate space score closer to 1), else 1 for every row.
#'   \item \code{adjusted_contribution}: \code{contribution * sparsity_score}. This, not raw \code{contribution}, is what every downstream column below is actually computed from.
#'   \item \code{rank}: descending rank of \code{adjusted_contribution} within the group; ties receive the average of the ranks they span.
#'   \item \code{n_features}: the number of rows (features) in the group.
#'   \item \code{rank_score}: \code{(n_features - rank) / (n_features - 1)}; ranges 0-1 with higher values indicating higher importance (a single-feature group scores 1).
#'   \item \code{cum_contrib}: the cumulative sum of \code{adjusted_contribution} in descending order; features tied on \code{adjusted_contribution} share the same \code{cum_contrib}, equal to the cumulative sum through the end of their tied block, so a tie is never split by arbitrary sort order.
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
  exclude_feature_types = "struct",
  filter_model = TRUE, 
  all_performance_parquet, 
  MCC_threshold = NULL, 
  compare_to_shuffled = TRUE, 
add_lasso_advtg = TRUE) 
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

  # The models are elastic net, fit with a mix of penalties from lasso to ridge.
  # Lasso shrinks the feature space to fewer selected variables, while
  # ridge retains variables (keeps correlated variables together).
  # Calculate the sparsity score to give an advantage to lasso-like fits.
all_perf <- arrow::read_parquet(normalizePath(all_performance_parquet)) |>
  dplyr::filter(!shuffled) |>
  dplyr::select(species, drug_label, drug_or_class, seed, 
    feature_type, feature_subtype, fit_penalty, fit_mixture, 
    mcc, n_feat, n_feats_returned) |> 
  dplyr::mutate(feat_return_ratio = n_feats_returned / n_feat,
  if(add_lasso_advtg) (sparsity_score = 1 - feat_return_ratio) else (sparsity_score = 1))
  
  # add different layers of scoring to the features within each seed 
  scored_top_features <- all_top_features |>
    dplyr::select(
      species, drug_label, drug_or_class, seed,
      feature_type, feature_subtype, variable=Variable,
      importance=Importance, sign=Sign
    ) |>
    dplyr::mutate(
      variable = dplyr::case_when(
        feature_type == "domains" ~ sub("_.+$", "", variable),
        feature_type == "proteins" ~ sub("fig.", "fig|", variable, fixed = TRUE),
        feature_type == "args" ~ sub(
          "^X", "",
          gsub("\\.NCBIFAM", "", variable)
        ),
        TRUE ~ variable
      )
    )|>
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
      contribution = importance / sum(importance, na.rm = TRUE),
      adjusted_contribution = contribution * sparsity_score
    ) |>
    dplyr::arrange(dplyr::desc(adjusted_contribution), .by_group = TRUE) |>
    dplyr::mutate(
      rank = rank(dplyr::desc(adjusted_contribution), ties.method = "average"),
      n_features = dplyr::n(),
      rank_score = dplyr::if_else(
        n_features > 1,
        (n_features - rank) / (n_features - 1),
        1
      ),
      running_contrib = cumsum(adjusted_contribution)
    ) |>
    dplyr::group_by(adjusted_contribution, .add = TRUE) |>
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
#' @param scored_features The tibble of scored top features with their contribution, rank, and rank score within each seed generated from `scoreFeaturesWithinSeed()`
#'
#' @returns a tibble with one row per species/drug_label/drug_or_class/feature_type/feature_subtype/variable, with:
#' \itemize{
#'   \item \code{seed_ratio}: the number of seeds the feature appears in, divided by the total number of distinct \code{seed} values present anywhere in \code{scored_features} (a single count computed once for the whole call, not per feature or per group — so this assumes every group was fit with the same set of seeds).
#'   \item \code{mean_rank_score}, \code{median_rank_score}: mean/median of \code{rank_score} across seeds; ranges 0-1 with higher values indicating higher importance.
#'   \item \code{median_rank}: median of \code{rank} across seeds. (\code{mean_rank} is not currently computed.)
#'   \item \code{median_contribution}: median of \code{adjusted_contribution} across seeds. (\code{mean_contribution} is not currently computed.)
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
summariseFeaturesAcrossSeeds <- function(scored_features) {
  
  # find max number of seeds
  max_seeds <- scored_features |>
  dplyr::summarise(n_seeds = dplyr::n_distinct(seed)) |>
  dplyr::pull(n_seeds)

 feature_summary <- scored_features |>
    dplyr::group_by(
      species,
      drug_label,
      drug_or_class,
      feature_type,
      feature_subtype,
      variable
    ) |>
    dplyr::summarise(
    seed_ratio = dplyr::n_distinct(seed)/max_seeds,
    # mean_rank = mean(rank, na.rm = TRUE),
    # mean_contribution = mean(contribution, na.rm = TRUE),
    # mean_cum_contrib = mean(cum_contrib, na.rm = TRUE),
    mean_rank_score = mean(rank_score, na.rm = TRUE),
    median_rank = median(rank, na.rm = TRUE),
    median_contribution = median(adjusted_contribution, na.rm = TRUE),
    median_cum_contrib = median(cum_contrib, na.rm = TRUE),
    median_rank_score = median(rank_score, na.rm = TRUE),
    best_rank = min(rank, na.rm = TRUE),
    rank_consistent = dplyr::n_distinct(rank) == 1,
    rank_score_sd = sd(rank_score, na.rm = TRUE),
    # coefficient of variation: how large rank_score_sd is relative to mean_rank_score
    rank_score_cv = dplyr::if_else(
      mean_rank_score != 0,
      rank_score_sd / mean_rank_score,
      NA_real_
    ),
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
#' @param rank_score_quantile A value in \[0, 1\] (default is \code{0.95}). Keep only features whose \code{median_rank_score} is at or above this quantile of \code{median_rank_score}. The quantile is computed once over every row of \code{feature_summary} — globally across all species/drugs/drug classes/feature types/subtypes, not per group — and is not restricted to rows that already pass the other conditions listed below: \code{dplyr::filter()} evaluates every condition passed to a single call against the same original, ungrouped data, so this threshold does not narrow as other conditions are applied
#' @param cv_threshold The maximum allowed coefficient of variation (default is \code{1}). Keep only features with \code{rank_score_cv <= cv_threshold}, i.e. drop features whose rank_score is inconsistent across seeds relative to its mean
#' @param cumulative_contribution_threshold The cumulative-contribution cutoff, in \[0, 1\] (default is \code{0.75}, i.e. 75%). Keep only features with \code{median_cum_contrib <= cumulative_contribution_threshold}
#' @param seed_ratio_threshold If not \code{NULL} (the default), keep only features whose \code{seed_ratio} exactly equals this value
#' @param both_subtypes Logical indicating whether to additionally restrict to features that survive the filters above in both the binary and counts \code{feature_subtype} (default is \code{FALSE})
#' @param compare_median_to_sd_rank_score Logical indicating whether to additionally require \code{median_rank_score > rank_score_sd} (default is \code{FALSE})
#'
#' @returns a tibble of top features for each drug/class: the `summariseFeaturesAcrossSeeds()` output (species, drug label, drug or class, feature type, feature subtype, variable, seed_ratio, mean/median rank score, median rank, median contribution, median cumulative contribution, best rank, rank/sign/in_core consistency flags, sign), filtered in two \code{dplyr::filter()} passes.
#' The first pass keeps rows where all of the following hold, evaluated together against the full, ungrouped \code{feature_summary} (see \code{rank_score_quantile} for what that means for the last condition):
#' \itemize{
#'   \item \code{seed_ratio == seed_ratio_threshold}, only applied when \code{seed_ratio_threshold} is not \code{NULL},
#'   \item \code{rank_score_cv <= cv_threshold},
#'   \item \code{in_core} is TRUE,
#'   \item \code{sign_consistent} is TRUE (sign is the same in every seed; this does not require the sign to be negative),
#'   \item \code{median_cum_contrib <= cumulative_contribution_threshold}, and
#'   \item \code{median_rank_score} is at or above the \code{rank_score_quantile} quantile of \code{median_rank_score}.
#' }
#' Two columns are then added, grouped by (species, drug_label, drug_or_class, feature_type, variable): \code{n_subtype} and \code{subtype_csv}, recording how many/which \code{feature_subtype} values each combination has among the rows that survived the first pass.
#' A second \code{dplyr::filter()} pass then optionally keeps only rows where \code{subtype_csv == "binary,counts"} (when \code{both_subtypes = TRUE}) and/or \code{median_rank_score > rank_score_sd} (when \code{compare_median_to_sd_rank_score = TRUE}).
#' Every drug/class may not have variables from all feature types.
#'
#' @export
topFeaturesPerDrugOrClass <- function( 
  all_top_features_parquet, 
  core_contribution_threshold = 0.75, 
  exclude_feature_types = "struct",
  filter_model = TRUE, 
  all_performance_parquet, 
  MCC_threshold = NULL, 
  compare_to_shuffled = TRUE,
  add_lasso_advtg = TRUE,
  # feature_summary = summariseFeaturesAcrossSeeds(scored_features),
                                        rank_score_quantile = 0.95,
                                        cv_threshold = 1,
                                        cumulative_contribution_threshold = 0.75,
                                        #additional filters
                                        seed_ratio_threshold = NULL,
                                        both_subtypes = FALSE,
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
  compare_to_shuffled = compare_to_shuffled,
  add_lasso_advtg = add_lasso_advtg
  )

feature_summary <- summariseFeaturesAcrossSeeds(scored_features)

top_features <- feature_summary |>
  dplyr::filter(
    if (!is.null(seed_ratio_threshold)) seed_ratio == seed_ratio_threshold else TRUE,
    rank_score_cv <= cv_threshold,
    in_core,
    sign_consistent,
    median_cum_contrib <= cumulative_contribution_threshold,
    median_rank_score >= quantile(median_rank_score, rank_score_quantile)
  ) |>
  dplyr::group_by(species, drug_label, drug_or_class, feature_type, variable) |>
    dplyr::mutate( n_subtype = dplyr::n_distinct(feature_subtype),
  subtype_csv = paste(sort(unique(feature_subtype)), collapse = ",") ) |>
  dplyr::ungroup() |>
dplyr::filter(
  if(both_subtypes) (subtype_csv == "binary,counts") else TRUE,
  if(compare_median_to_sd_rank_score) (median_rank_score > rank_score_sd) else TRUE
    )
  
  return(top_features)
}

#' Aggregate mapped features to protein clusters
#'
#' @param top_features The tibble of top features for each drug/class generated from `topFeaturesPerDrugOrClass()`
#' @param cluster_feature_parquet The path to the Parquet file containing the mapping of features to protein clusters
#'
#' @returns a tibble with one row per species/drug_label/drug_or_class/cluster, with:
#' \itemize{
#'   \item \code{frequency}: the number of top-feature rows mapped to this cluster.
#'   \item \code{n_variables}, \code{variables_csv}: number of, and comma-separated list of, distinct \code{variable} values mapped to this cluster.
#'   \item \code{n_feature_types}, \code{feature_types_csv}: number of, and comma-separated list of, distinct \code{feature_type} values mapped to this cluster.
#'   \item \code{cluster_mean_rank_score}, \code{cluster_median_rank_score}: mean/median of \code{median_rank_score}, down-weighted by \code{1 / n_clusters} for features that map to multiple clusters so promiscuous domains/COGs don't inflate every cluster they touch.
#'   \item \code{cluster_max_rank_score}: the max of the unweighted \code{median_rank_score}, i.e. the single strongest feature backing this cluster regardless of its fan-out to other clusters.
#'   \item \code{cluster_rank_score_sd}: standard deviation of the same down-weighted score used for \code{cluster_mean_rank_score}/\code{cluster_median_rank_score}.
#'   \item \code{cluster_best_rank}: the best (lowest) \code{best_rank} among the features mapped to this cluster.
#' }
#' \code{median_rank_score} is used throughout (rather than \code{mean_rank_score}) for consistency with the seed-noise-robust selection made in \code{topFeaturesPerDrugOrClass()}.
#' \code{top_features} rows whose \code{variable} has no match in \code{cluster_feature_parquet} are not dropped: they collapse into one \code{cluster = NA} row per species/drug_label/drug_or_class, aggregating every unmapped feature for that group. Callers that want only real clusters must filter this out explicitly (as \code{findSharedClusters()} and \code{findUniqueClusters()} do).
#'
#' @export
summariseClusters <- function(top_features = topFeaturesPerDrugOrClass(feature_summary),
                               cluster_feature_parquet
                              ) {
  stopifnot(file.exists(cluster_feature_parquet))
  cluster_feature <- arrow::read_parquet(normalizePath(cluster_feature_parquet)) |>
    dplyr::add_count(feature, name = "n_clusters")

  top_clusters <- top_features |>
    dplyr::left_join(
      cluster_feature,
      by = dplyr::join_by(variable == feature),
      relationship = "many-to-many"
    ) |>
    dplyr::mutate(
      n_clusters = dplyr::coalesce(n_clusters, 1L),
      weighted_rank_score = median_rank_score / n_clusters
    ) |>
    dplyr::group_by(
      species,
      drug_label,
      drug_or_class,
      cluster
    ) |>
    dplyr::summarise(
      frequency = dplyr::n(),
      n_variables = dplyr::n_distinct(variable),
      variables_csv = paste(sort(unique(variable)), collapse = ","),
      n_feature_types = dplyr::n_distinct(feature_type),
      feature_types_csv = paste(sort(unique(feature_type)), collapse = ","),
      cluster_mean_rank_score = mean(weighted_rank_score, na.rm = TRUE),
      cluster_median_rank_score = median(weighted_rank_score, na.rm = TRUE),
      cluster_max_rank_score = max(median_rank_score, na.rm = TRUE),
      cluster_rank_score_sd = sd(weighted_rank_score, na.rm = TRUE),
      cluster_best_rank = min(best_rank, na.rm = TRUE),
      .groups = "drop"
    ) |>
  dplyr::arrange(dplyr::desc(frequency), dplyr::desc(n_feature_types))

  return(top_clusters)
}

#' Find clusters that appear across multiple drugs/classes
#'
#' @param top_clusters The tibble of summarized clusters generated from `summariseClusters()`
#' @param label The \code{drug_label} value to filter clusters by, either \code{"drug"} or \code{"drug_class"} (default is \code{"drug"})
#' @param min_drugs_or_classes The minimum number of distinct drugs or classes required for a cluster to be considered shared (default is \code{2})
#'
#' @returns a tibble with one row per shared \code{cluster}, with \code{n_drug_or_class} (the number of distinct \code{drug_or_class} values the cluster appears in) and \code{drug_or_class_csv} (a comma-separated string of those values), sorted by \code{n_drug_or_class} descending.
#'
#' @export
findSharedClusters <- function(top_clusters = summariseClusters(top_features, cluster_feature_parquet),
                                label = "drug",
                                 min_drugs_or_classes = 2
                                ) {
  shared_clusters <- top_clusters |>
    dplyr::filter(!is.na(cluster), drug_label == label) |>
    dplyr::group_by(cluster) |>
    dplyr::mutate(
      n_drug_or_class = dplyr::n_distinct(drug_or_class),
      drug_or_class_csv = paste(sort(unique(drug_or_class)), collapse = ", ")
    ) |>
    dplyr::filter(n_drug_or_class >= min_drugs_or_classes) |>
    dplyr::ungroup() |>
    dplyr::select(cluster, n_drug_or_class, drug_or_class_csv) |>
    dplyr::arrange(dplyr::desc(n_drug_or_class))

  return(shared_clusters)
}

#' Find the clusters that are unique to a single drug/class
#'
#' @param top_clusters The tibble of summarized clusters generated from `summariseClusters()`
#' @param label The \code{drug_label} value to filter clusters by, either \code{"drug"} or \code{"drug_class"} (default is \code{"drug"})
#' @param protein_names_parquet The path to the Parquet file containing the annotations to protein cluster names
#'
#' @returns a tibble with one row per cluster unique to a single drug/class, with \code{drug_or_class}, \code{cluster}, \code{cluster_name} (from the protein name annotations), and \code{cluster_mean_rank_score}, sorted by \code{cluster_mean_rank_score} descending.
#'
#' @export
#' @examples
#' findUniqueClusters(summariseClusters(top_features, cluster_feature_parquet), label = "drug", protein_names_parquet)
findUniqueClusters <- function(top_clusters = summariseClusters(top_features, cluster_feature_parquet),
                            label = "drug",
                            protein_names_parquet
) {

  protein_names <- arrow::read_parquet(normalizePath(protein_names_parquet)) |>
    dplyr::distinct()

  unique_clusters <- top_clusters |>
    dplyr::filter(!is.na(cluster), drug_label == label) |>
    dplyr::group_by(cluster) |>
    dplyr::mutate(
      n_drug_or_class = dplyr::n_distinct(drug_or_class),
      drug_or_class_csv = paste(sort(unique(drug_or_class)), collapse = ", ")
    ) |>
    dplyr::filter(n_drug_or_class == 1) |>
    dplyr::ungroup() |>
    dplyr::select(cluster, drug_or_class_csv, cluster_mean_rank_score) |>
    dplyr::arrange(dplyr::desc(cluster_mean_rank_score)) |>
dplyr::rename(drug_or_class = drug_or_class_csv) |>
    dplyr::left_join(protein_names, by = dplyr::join_by(cluster == proteinID)) |>
    dplyr::rename(cluster_name = proteinName) |>
    dplyr::select(drug_or_class, cluster, cluster_name, cluster_mean_rank_score)

  return(unique_clusters)
}

#' Build a feature network from selected top features and top clusters
#'
#' @param top_features Output of \code{topFeaturesPerDrugOrClass()}.
#' @param top_clusters Output of \code{summariseClusters()}.
#' @param cluster_feature_parquet Path to the Parquet file mapping variables to clusters.
#' @param protein_names_parquet Path to the Parquet file with cluster name annotations.
#'
#' @returns A list with \code{feature_table}, \code{cluster_table}, \code{nodes}, \code{edges}, and \code{graph}.
#' Node/edge weights are built from \code{median_rank_score} (features) and \code{cluster_mean_rank_score} (clusters, itself median-based per \code{summariseClusters()}), consistent with the seed-noise-robust selection made in \code{topFeaturesPerDrugOrClass()}.
#' @export
buildFeatureNetwork <- function(top_features,
                               top_clusters,
                               cluster_feature_parquet,
                               protein_names_parquet
                              ) {
  stopifnot(is.data.frame(top_features))
  stopifnot(is.data.frame(top_clusters))
  stopifnot(file.exists(cluster_feature_parquet))
  stopifnot(file.exists(protein_names_parquet))

  required_feature_cols <- c(
    "species", "drug_label", "drug_or_class",
    "feature_type", "variable",
    "median_rank_score"
  )
  required_cluster_cols <- c(
    "species", "drug_label", "drug_or_class",
    "cluster", "cluster_mean_rank_score"
  )

  missing_feature_cols <- setdiff(required_feature_cols, names(top_features))
  missing_cluster_cols <- setdiff(required_cluster_cols, names(top_clusters))

  if (length(missing_feature_cols) > 0) {
    stop("top_features is missing required columns: ",
         paste(missing_feature_cols, collapse = ", "))
  }
  if (length(missing_cluster_cols) > 0) {
    stop("top_clusters is missing required columns: ",
         paste(missing_cluster_cols, collapse = ", "))
  }

  cluster_feature <- arrow::read_parquet(normalizePath(cluster_feature_parquet)) |>
    dplyr::distinct()

  protein_names <- arrow::read_parquet(normalizePath(protein_names_parquet)) |>
    dplyr::distinct()

  make_model_id <- function(drug_label, drug_or_class) {
    paste(drug_label, drug_or_class, sep = ".")
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

  cluster_table <- top_clusters |>
    dplyr::mutate(model_id = make_model_id(drug_label, drug_or_class)) |>
    dplyr::left_join(
      protein_names,
      by = dplyr::join_by(cluster == proteinID)
    ) |>
    dplyr::group_by(species, model_id, cluster, proteinName) |>
    dplyr::summarise(
      cluster_score = mean(cluster_mean_rank_score, na.rm = TRUE),
      .groups = "drop"
    )

  model_nodes <- dplyr::bind_rows(
    feature_table |>
      dplyr::distinct(species, model_id),
    cluster_table |>
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

  cluster_nodes <- cluster_table |>
    dplyr::group_by(species, cluster, proteinName) |>
    dplyr::summarise(
      score = mean(cluster_score, na.rm = TRUE),
      breadth = dplyr::n_distinct(model_id),
      .groups = "drop"
    ) |>
    dplyr::transmute(
      name = cluster,
      label = dplyr::if_else(
        is.na(proteinName) | proteinName == "",
        cluster,
        proteinName
      ),
      node_type = "cluster",
      species = species,
      score = score,
      breadth = breadth,
      node_size = pmax(3, pmin(10, breadth + 2))
    )

  nodes <- dplyr::bind_rows(model_nodes, feature_nodes, cluster_nodes) |>
    dplyr::distinct(name, .keep_all = TRUE)

  feature_edges <- feature_table |>
    dplyr::transmute(
      from = model_id,
      to = variable,
      weight = feature_score,
      edge_type = "model_feature"
    ) |>
    dplyr::distinct(from, to, edge_type, .keep_all = TRUE)

  cluster_edges <- cluster_table |>
    dplyr::transmute(
      from = model_id,
      to = cluster,
      weight = cluster_score,
      edge_type = "model_cluster"
    ) |>
    dplyr::distinct(from, to, edge_type, .keep_all = TRUE)

  feature_cluster_edges <- feature_table |>
    dplyr::left_join(
      cluster_feature |> dplyr::add_count(feature, name = "n_clusters"),
      by = dplyr::join_by(variable == feature),
      relationship = "many-to-many"
    ) |>
    dplyr::filter(!is.na(cluster)) |>
    dplyr::transmute(
      from = variable,
      to = cluster,
      weight = 1 / n_clusters,
      edge_type = "feature_cluster"
    ) |>
    dplyr::distinct(from, to, edge_type, .keep_all = TRUE)

  edges <- dplyr::bind_rows(feature_edges, feature_cluster_edges)

    edges <- dplyr::bind_rows(edges, cluster_edges)


  missing_vertices <- setdiff(unique(c(edges$from, edges$to)), nodes$name)
  if (length(missing_vertices) > 0) {
    extra_nodes <- tibble::tibble(name = missing_vertices) |>
      dplyr::mutate(
        label = name,
        node_type = dplyr::case_when(
          grepl("^drug\\.|^drug_class\\.", name) ~ "model",
          grepl("^fig\\||^cluster", name) ~ "cluster",
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
    cluster_table = cluster_table,
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
    # dplyr::filter(
    #   show_direct_model_cluster | edge_type != "model_cluster"
    # ) |>
    dplyr::left_join(
      nodes |> dplyr::select(name, id),
      by = c("from" = "name")
    ) |>
    dplyr::rename(source = id) |>
    dplyr::left_join(
      nodes |> dplyr::select(name, id),
      by = c("to" = "name")
    ) |>
    dplyr::rename(target = id) |>
    dplyr::filter(!is.na(source), !is.na(target)) |>
    dplyr::mutate(
      value = dplyr::if_else(is.na(weight), 1, weight)
    ) |>
    dplyr::select(source, target, value, edge_type)

  stopifnot(nrow(nodes) > 0)
  stopifnot(nrow(links) > 0)

  colour_scale <- networkD3::JS(
    "d3.scaleOrdinal()
      .domain(['model', 'feature', 'cluster'])
      .range(['#4C78A8', '#F58518', '#54A24B'])"
  )

  networkD3::forceNetwork(
    Links = links,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "label",
    Group = "group",
    opacity = 0.9,
    zoom = TRUE,
    fontSize = 14,
    # nodeWidth = 24,
    height = height,
    width = width,
    colourScale = colour_scale,
    linkDistance = networkD3::JS(
      "function(d) {
         if (d.edge_type === 'feature_cluster') return 60;
         if (d.edge_type === 'model_feature') return 120;
         return 90;
       }"
    )
  )
}


# final run would be:
# top_features <- topFeaturesPerDrugOrClass(rank_score_quantile = 0.75)
# top_clusters <- summariseClusters(top_features, cluster_feature_parquet = cluster_feature_parquet)
# feature_network <- buildFeatureNetwork(top_features = top_features, top_clusters = top_clusters,
#   cluster_feature_parquet = cluster_feature_parquet, protein_names_parquet = protein_names_parquet)
# plotFeatureNetworkD3(feature_network)
