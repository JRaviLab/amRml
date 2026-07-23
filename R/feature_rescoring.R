# computeFeatureImprovement <- function(
#   all_top_features_parquet,
#   cluster_feature_parquet
# ) {
#   stopifnot(file.exists(all_top_features_parquet))
#   cluster_feature <- arrow::read_parquet(normalizePath(cluster_feature_parquet))


#   features_rescored <- arrow::read_parquet(normalizePath(all_top_features_parquet)) |>
#     dplyr::select(
#       species, drug_label, drug_or_class, shuffled,
#       feature_type, feature_subtype, Variable,
#       Importance, Sign
#     ) |>
#     dplyr::filter(!pca) |>
#     dplyr::mutate(
#       Variable = dplyr::case_when(
#         feature_type == "domains" ~ sub("_.+$", "", Variable),
#         feature_type == "proteins" ~ sub("fig.", "fig|", Variable, fixed = TRUE),
#         TRUE ~ Variable
#       )
#     ) |>
#     dplyr::group_by(output_prefix) |>
#     dplyr::mutate(
#       rank = dplyr::dense_rank(dplyr::desc(Importance)),
#       denom = sum(Importance, na.rm = TRUE),
#       contribution = dplyr::if_else(denom > 0, Importance / denom, 0),

#       # Safer rescaling within each output_prefix
#       min_imp = min(Importance, na.rm = TRUE),
#       max_imp = max(Importance, na.rm = TRUE),
#       range_imp = max_imp - min_imp,
#       rescaled = dplyr::if_else(range_imp > 0, (Importance - min_imp) / range_imp, 0)
#     ) |>
#     dplyr::ungroup() |>
#     dplyr::group_by(drug_label, drug_or_class, feature_type, feature_subtype, shuffled, Variable) |>
#     dplyr::mutate(median_datatype = stats::median(rank, na.rm = TRUE)) |>
#     dplyr::ungroup() |>
#     dplyr::group_by(drug_label, drug_or_class, feature_type, shuffled, Variable) |>
#     dplyr::mutate(median_scale = stats::median(median_datatype, na.rm = TRUE)) |>
#     dplyr::ungroup() |>
#     # NOTE: correct join syntax (no quotes in join_by)
#     dplyr::left_join(feature_cluster, by = dplyr::join_by(Variable == feature)) |>
#     dplyr::group_by(drug_label, drug_or_class, shuffled, cluster) |>
#     dplyr::mutate(
#       median_drug_or_class       = stats::median(rank, na.rm = TRUE),
#       count_scales_for_cluster   = dplyr::n_distinct(feature_type),
#       feature_types_csv          = paste(sort(unique(feature_type)), collapse = ","),
#       lowest_contri              = min(contribution, na.rm = TRUE),
#       highest_contri             = max(contribution, na.rm = TRUE)
#     ) |>
#     dplyr::ungroup() |>
#     dplyr::group_by(drug_label, drug_or_class, cluster) 

#   return(features_rescored)
# }

# #' Score and filter top features based on importance, contribution, and rank across seeds.
# #'
# #' @param all_top_features_parquet The path to the Parquet file containing all top features with their importance scores.
# #'
# #' @returns A data frame of top features with their scores.
# #'
# #' @export
# #' @examples
# scoreTopFeatures <- function(
#   all_top_features_parquet
# ) {
#   stopifnot(file.exists(all_top_features_parquet))
#   all_top_features <- arrow::read_parquet(normalizePath(all_top_features_parquet)) |>
#     dplyr::filter(!shuffled) |>
#     dplyr::select(
#       species, drug_label, drug_or_class, seed,
#       feature_type, feature_subtype, Variable,
#       Importance, Sign
#     ) |>
#     dplyr::mutate(
#       Variable = dplyr::case_when(
#         feature_type == "domains" ~ sub("_.+$", "", Variable),
#         feature_type == "proteins" ~ sub("fig.", "fig|", Variable, fixed = TRUE),
#         feature_type == "args" ~ sub(
#           "^X", "",
#           gsub("\\.NCBIFAM", "", Variable)
#         ),
#         TRUE ~ Variable
#       )
#     )

#   ## --------------------------------
#   ## Importance → contribution → rank
#   ## --------------------------------

#   ranked_features <- all_top_features |>
#   dplyr::group_by(
#     species,
#     drug_label,
#     drug_or_class,
#     feature_type,
#     feature_subtype,
#     seed
#   ) |>
#   dplyr::mutate(
#     contribution = Importance / sum(Importance, na.rm = TRUE),
#     rank = dplyr::dense_rank(dplyr::desc(contribution)),
#     n_features = dplyr::n(),
#     rank_score = ifelse(n_features > 1, (n_features - rank) / (n_features - 1), 1)
#   ) |>
#   dplyr::ungroup() |>
#   dplyr::group_by(
#     species,
#     drug_label,
#     drug_or_class,
#     feature_type,
#     feature_subtype,
#     Variable
#   ) |>
#   dplyr::summarise(
#     n_seeds = dplyr::n_distinct(seed),
#     mean_rank = mean(rank, na.rm = TRUE),
#     mean_contribution = mean(contribution, na.rm = TRUE),
#     median_rank = median(rank, na.rm = TRUE),
#     median_contribution = median(contribution, na.rm = TRUE),
#     mean_rank_score = mean(rank_score, na.rm = TRUE),
#     best_rank = min(rank, na.rm = TRUE),
#     rank_consistent = dplyr::n_distinct(rank) == 1,
#     rank_sd = sd(rank, na.rm = TRUE),
#     sign_consistent = dplyr::n_distinct(Sign) == 1,
#     sign = if (sign_consistent) dplyr::first(Sign) else "mixed",
#     .groups = "drop"
#   ) |>
#   dplyr::arrange(dplyr::desc(mean_rank_score), mean_rank, best_rank)
  
#   # filter for top features: negative sign, low rank SD, high mean rank score, and max seeds
#   top_features <- ranked_features |>
#   dplyr::filter(
#     sign == "NEG",
#     rank_sd <= 1
#   ) |>
#   dplyr::group_by(drug_label, drug_or_class) |>
#   dplyr::filter(
#     n_seeds == max(n_seeds),
#     mean_rank_score >= quantile(mean_rank_score, 0.95)
#   ) |>
#   dplyr::ungroup() |>
#     dplyr::select(
#       species, drug_label, drug_or_class,
#       feature_type, feature_subtype, Variable,
#       mean_rank_score, mean_rank, best_rank,
#       median_rank, mean_contribution, median_contribution,
#       n_seeds, rank_sd,
#       sign_consistent, sign
#     )
  
#   return(top_features)
# }

# findTopClusters <- function(
#   all_top_features_parquet,
#   cluster_feature_parquet
# ) {
#   stopifnot(file.exists(all_top_features_parquet))
#   stopifnot(file.exists(cluster_feature_parquet))
#   all_top_features <- arrow::read_parquet(normalizePath(all_top_features_parquet)) |>
#     dplyr::filter(!shuffled) |>
#     dplyr::select(
#       species, drug_label, drug_or_class, seed,
#       feature_type, feature_subtype, Variable,
#       Importance, Sign
#     ) |>
#     dplyr::mutate(
#       Variable = dplyr::case_when(
#         feature_type == "domains" ~ sub("_.+$", "", Variable),
#         feature_type == "proteins" ~ sub("fig.", "fig|", Variable, fixed = TRUE),
#         feature_type == "args" ~ sub(
#           "^X", "",
#           gsub("\\.NCBIFAM", "", Variable)
#         ),
#         TRUE ~ Variable
#       )
#     )
  
#   cluster_feature <- arrow::read_parquet(normalizePath(cluster_feature_parquet))

#   ## --------------------------------
#   ##  Scale capability (data-driven)
#   ## --------------------------------
#   scale_capability <- all_top_features |>
#     dplyr::distinct(feature_type, feature_subtype) |>
#     dplyr::group_by(feature_type) |>
#     dplyr::summarise(
#       expected_types = dplyr::n_distinct(feature_subtype),
#       expected_types_csv = paste(sort(feature_subtype), collapse = ","),
#       .groups = "drop"
#     )

#   ## --------------------------------
#   ##  Collapse to protein level
#   ## --------------------------------

#   protein_scale_realization <- all_top_features |>
#   dplyr::group_by(
#     species,
#     drug_label,
#     drug_or_class,
#     feature_type,
#     feature_subtype,
#     seed
#   ) |>
#   dplyr::mutate(
#     contribution = Importance / sum(Importance, na.rm = TRUE),
#     rank = dplyr::dense_rank(dplyr::desc(contribution)),
#     n_features = dplyr::n(),
#     rank_score = ifelse(n_features > 1, (n_features - rank) / (n_features - 1), 1)
#   ) |>
#  dplyr::ungroup() |>  
#     dplyr::filter(rank_score >= 0.95, Sign == "NEG") |>
#     dplyr::left_join(cluster_feature, by = dplyr::join_by(Variable == feature)) |>
#     dplyr::filter(!is.na(cluster)) |>
#     dplyr::group_by(
#       species,
#       drug_label,
#       drug_or_class,
#       cluster
#     ) |>
#     dplyr::mutate(
#       mean_rank_cluster = mean(rank, na.rm = TRUE),
#       sd_rank_cluster = sd(rank, na.rm = TRUE),
#       mean_contribution_cluster = mean(contribution, na.rm = TRUE)
#     )
#     dplyr::group_by(
#       species,
#       drug_label,
#       drug_or_class,
#       cluster
#     ) |>
#     dplyr::mutate(
#       observed_scales = dplyr::n_distinct(feature_type),
#       observed_scales_csv = paste(sort(unique(feature_type)), collapse = ",")
#     ) |>
#     dplyr::group_by(
#       feature_type) |>
#     dplyr::mutate(
#       observed_types = dplyr::n_distinct(feature_subtype),
#       observed_types_csv = paste(sort(unique(feature_subtype)), collapse = ",")
#     ) |>
#     dplyr::left_join(scale_capability, by = "feature_type") |>
#     dplyr::mutate(
#       scale_realization = dplyr::case_when(
#         observed_types == expected_types ~ "full_realization",
#         observed_types < expected_types ~ "partial_realization"
#       )
#     ) |>
#     dplyr::ungroup() |>
#     dplyr::ungroup()

#   protein_scale_summary <- protein_scale_realization |>
#     dplyr::group_by(species, drug_label, drug_or_class, cluster) |>
#     dplyr::summarise(
#       n_scales = dplyr::n_distinct(feature_type),
#       fully_realized_scales =
#         sum(scale_realization == "full_realization"),
#       partially_realized_scales =
#         sum(scale_realization == "partial_realization"),
#       scale_support_csv =
#         paste(
#           feature_type,
#           "(", observed_types_csv, "/", expected_types, ")",
#           collapse = "; "
#         ),
#       .groups = "drop"
#     ) |>
#     dplyr::mutate(
#       realization_score =
#         (fully_realized_scales +
#           0.5 * partially_realized_scales) /
#           (fully_realized_scales + partially_realized_scales)
#     ) |>
#     dplyr::mutate(
#       coverage_boost =
#         n_scales / max(n_scales)
#     ) |>
#     dplyr::mutate(
#       scale_factor = realization_score * coverage_boost
#     )

#   ## --------------------------------
#   ## 4. Shuffle vs non-shuffle
#   ## --------------------------------

#   shuffle_delta <- ranked_features |>
#     dplyr::filter(!is.na(cluster)) |>
#     dplyr::group_by(
#       species, drug_label,
#       drug_or_class,
#       cluster
#     ) |>
#     dplyr::summarise(
#       mean_rank_nonshuffle = mean(mean_rank[shuffled == FALSE], na.rm = TRUE),
#       mean_rank_shuffle = mean(mean_rank[shuffled == TRUE], na.rm = TRUE),
#       delta_rank = mean_rank_shuffle - mean_rank_nonshuffle,
#       # If non-shuffled is missing -> NA (no evidence). If shuffled missing -> +Inf improvement.
#       improvement = dplyr::case_when(
#         !is.na(mean_rank_nonshuffle) ~ tidyr::replace_na(mean_rank_shuffle, Inf) - mean_rank_nonshuffle,
#         TRUE ~ NA_real_
#       ),

#       # "Good" if non-shuffled exists AND (non-shuffled < shuffled OR shuffled is missing)
#       good_feature = !is.na(mean_rank_nonshuffle) & improvement > 0,
#       .groups = "drop"
#     )

#   ## --------------------------------
#   ## 5. Cluster scoring
#   ## --------------------------------
#   robustness <- ranked_features |>
#     dplyr::filter(!is.na(cluster), shuffled == FALSE) |>
#     dplyr::group_by(species, drug_label, drug_or_class, cluster) |>
#     dplyr::summarize(
#       median_contribution = median(contribution),
#       median_rank = median(mean_rank), .groups = "drop"
#     ) |>
#     dplyr::left_join(
#       protein_scale_summary |>
#         dplyr::distinct(species, drug_label, drug_or_class, cluster, n_scales, scale_factor),
#       by = c("species", "drug_label", "drug_or_class", "cluster")
#     ) |>
#     dplyr::left_join(
#       shuffle_delta |>
#         dplyr::distinct(species, drug_label, drug_or_class, cluster, delta_rank, good_feature),
#       by = c("species", "drug_label", "drug_or_class", "cluster")
#     ) |>
#     dplyr::group_by(species, drug_label, drug_or_class) |>
#     dplyr::mutate(
#       # higher = better
#       contrib_score = dplyr::percent_rank(median_contribution),

#       # lower rank = better → invert
#       stability_score = 1 - dplyr::percent_rank(median_rank),

#       # robustness: handle NaN = strongest case (missing in shuffle)
#       delta_score = dplyr::case_when(
#         is.nan(delta_rank) ~ 1, # best possible signal
#         delta_rank > 0 ~ dplyr::percent_rank(delta_rank), # reward
#         delta_rank == 0 ~ 0, # neutral
#         delta_rank < 0 ~ -dplyr::percent_rank(abs(delta_rank)) # penalize
#       ),

#       # already bounded [0,1]
#       scale_score = scale_factor,
#       robustness_score =
#         contrib_score *
#           stability_score *
#           scale_factor *
#           delta_score
#     )

#   return(robustness)
# }

#' Score features within each seed
#'
#' @param all_top_features_parquet The path to the Parquet file containing all top features with their importance scores.
#'
#' @returns a tibble of scored top features with their contribution, rank, and rank score within each seed.
#' within a seed
#' contribution is calculated as the importance of a feature divided by the sum of importance scores for all features. 
#' Rank is assigned based on the descending order of contribution, and 
#' rank score is calculated as (n_features - rank) / (n_features - 1), where n_features is the total number of features within the same seed.
#' rank score is ranged between 0 and 1, with higher values indicating higher importance.
#' 
#' @export
#' @examples
scoreFeaturesWithinSeed <- function(all_top_features_parquet) {

  stopifnot(file.exists(all_top_features_parquet))
  scored_top_features <- arrow::read_parquet(normalizePath(all_top_features_parquet)) |>
    dplyr::filter(!shuffled) |>
    dplyr::select(
      species, drug_label, drug_or_class, seed,
      feature_type, feature_subtype, Variable,
      Importance, Sign
    ) |>
    dplyr::mutate(
      Variable = dplyr::case_when(
        feature_type == "domains" ~ sub("_.+$", "", Variable),
        feature_type == "proteins" ~ sub("fig.", "fig|", Variable, fixed = TRUE),
        feature_type == "args" ~ sub(
          "^X", "",
          gsub("\\.NCBIFAM", "", Variable)
        ),
        TRUE ~ Variable
      )
    )|>
    dplyr::group_by(
      species,
      drug_label,
      drug_or_class,
      feature_type,
      feature_subtype,
      seed
    ) |>
    dplyr::mutate(
      contribution = Importance / sum(Importance, na.rm = TRUE),
      rank = dplyr::dense_rank(dplyr::desc(contribution)),
      n_features = dplyr::n(),
      rank_score = dplyr::if_else(
        n_features > 1,
        (n_features - rank) / (n_features - 1),
        1
      )
    ) |>
    dplyr::ungroup()

  return(scored_top_features)
}

#' Summarize a variable across seeds
#    "which molecular features are consistent?"
#'
#' @param scored_features The tibble of scored top features with their contribution, rank, and rank score within each seed generated from `scoreFeaturesWithinSeed()`
#'
#' @returns a tibble of summarized features across seeds, including the number of seeds, 
#' mean and median rank, 
#' mean and median contribution, 
#' mean rank score (scaled between 0 and 1 with higher values indicating higher importance), 
#' best rank, 
#' rank consistency (TRUE if the rank is same across seeds), rank standard deviation (how much the rank varies across seeds), 
#' sign consistency (TRUE if the sign is same across seeds), and sign (the sign of the feature; POS, NEG or MIXED).
#'
#' @export
#' @examples
summariseFeatureAcrossSeeds <- function(scored_features = scoreFeaturesWithinSeed(all_top_features_parquet)) {
 feature_summary <- scored_features |>
    dplyr::group_by(
      species,
      drug_label,
      drug_or_class,
      feature_type,
      feature_subtype,
      Variable
    ) |>
    dplyr::summarise(
    n_seeds = dplyr::n_distinct(seed),
    mean_rank = mean(rank, na.rm = TRUE),
    mean_contribution = mean(contribution, na.rm = TRUE),
    median_rank = median(rank, na.rm = TRUE),
    median_contribution = median(contribution, na.rm = TRUE),
    mean_rank_score = mean(rank_score, na.rm = TRUE),
    best_rank = min(rank, na.rm = TRUE),
    rank_consistent = dplyr::n_distinct(rank) == 1,
    rank_sd = sd(rank, na.rm = TRUE),
    sign_consistent = dplyr::n_distinct(Sign) == 1,
    sign = if (sign_consistent) dplyr::first(Sign) else "MIXED",
    .groups = "drop"
  ) 

  return(feature_summary)
}

#' Build a per-drug top-feature table with cutoffs
#' This returns the top features for each drug/class.
#' 
#' @param feature_summary The tibble of summarized features across seeds generated from `summariseFeatureAcrossSeeds()`
#' @param rank_score_quantile The quantile threshold for filtering features based on their mean rank scores
#' @param threshold_sd_rank The threshold for filtering features based on the standard deviation of ranks across seeds
#'
#' @returns a tibble of top features for each drug/class, 
#' filtered based on the specified rank score quantile and standard deviation threshold. 
#' The resulting filtered table includes species, drug label, drug or class, feature type, feature subtype, variable, 
#' mean rank score, mean rank, best rank, median rank, mean contribution, median contribution, number of seeds, 
#' rank standard deviation, sign consistency, and sign.
#' The features are filtered to include only those with a negative sign, low rank standard deviation,
#' and are grouped by drug label and drug or class, retaining only the features with the
#' maximum number of seeds and mean rank scores above the specified quantile threshold.
#' Every drug/class may not have Variables from all feature types. 
#'
#' @export
#' @examples
topFeaturesPerDrugOrClass <- function(feature_summary = summariseFeatureAcrossSeeds(scored_features),
                                        rank_score_quantile = 0.95,
                                        threshold_sd_rank = 1
                                        ) {
  
top_features <- feature_summary |>
  dplyr::filter(
    sign == "NEG",
    rank_sd <= threshold_sd_rank
  ) |> 
  dplyr::group_by(species, drug_label, drug_or_class, feature_type, Variable) |>
    dplyr::mutate( n_subtype = dplyr::n_distinct(feature_subtype), 
  subtype_csv = paste(sort(unique(feature_subtype)), collapse = ",") ) |>
  dplyr::ungroup() |>
  dplyr::group_by(drug_label, drug_or_class) |>
  dplyr::filter(
    n_seeds == max(n_seeds),
    mean_rank_score >= quantile(mean_rank_score, rank_score_quantile)
  ) |>
  dplyr::ungroup() |>
    dplyr::distinct(
      species, drug_label, drug_or_class,
      feature_type, feature_subtype, Variable, n_subtype, subtype_csv,
      mean_rank_score, mean_rank, best_rank,
      median_rank, mean_contribution, median_contribution,
      n_seeds, rank_sd, 
      sign_consistent, sign
    ) |>
dplyr::filter(sign_consistent)

  return(top_features)
}

#' Aggregate mapped features to protein clusters
#'
#' @param top_features The tibble of top features for each drug/class generated from `topFeaturesPerDrugOrClass()`
#' @param cluster_feature_parquet The path to the Parquet file containing the mapping of features to protein clusters
#'
#' @returns a tibble of summarized clusters, including species, drug label, drug or class, cluster, frequency (number of features in the cluster),
#' number of distinct variables, number of distinct feature types, feature types as a comma-separated string,
#' cluster mean rank score (mean of mean_rank_score), cluster median rank score (median of mean_rank_score), 
#' cluster max rank score (max of mean_rank_score), cluster rank score standard deviation, and cluster best rank.
#'
#' @export
#' @examples
summariseClusters <- function(top_features = topFeaturesPerDrugOrClass(feature_summary),
                               cluster_feature_parquet
                              ) {
  stopifnot(file.exists(cluster_feature_parquet))
  cluster_feature <- arrow::read_parquet(normalizePath(cluster_feature_parquet)) 

  top_clusters <- top_features |>
    dplyr::left_join(cluster_feature, by = dplyr::join_by(Variable == feature)) |>
    dplyr::group_by(
      species,
      drug_label,
      drug_or_class,
      cluster 
    ) |>
    dplyr::summarise(
      frequency = dplyr::n(),
      n_variables = dplyr::n_distinct(Variable),
      variables_csv = paste(sort(unique(Variable)), collapse = ","),
      n_feature_types = dplyr::n_distinct(feature_type),
      feature_types_csv = paste(sort(unique(feature_type)), collapse = ","),
      cluster_mean_rank_score = mean(mean_rank_score, na.rm = TRUE),
      cluster_median_rank_score = median(mean_rank_score, na.rm = TRUE),
      cluster_max_rank_score = max(mean_rank_score, na.rm = TRUE),
      cluster_rank_score_sd = sd(mean_rank_score, na.rm = TRUE),
      cluster_best_rank = min(best_rank, na.rm = TRUE),
      .groups = "drop"
    ) |>
  dplyr::arrange(dplyr::desc(frequency), dplyr::desc(n_feature_types))

  return(top_clusters)
}

#' Find clusters that appear across multiple drugs/classes
#'
#' @param top_clusters The tibble of summarized clusters generated from `summariseClusters()`
#' @param label The label to filter clusters by, either "drug" or "drug_class"
#' @param  min_drugs_or_classes The minimum number of distinct drugs or classes required for a cluster to be considered shared
#'
#' @returns a tibble of shared clusters, including cluster, number of distinct drugs or classes, and a comma-separated string of the distinct drugs or classes.
#'
#' @export
#' @examples
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

#' find the clusters that are unique to a single drug/class
#'
#' @param top_clusters The tibble of summarized clusters generated from `summariseClusters()`
#' @param label The label to filter clusters by, either "drug" or "drug_class"
#' @param protein_names_parquet The path to the Parquet file containing the annotations to protein cluster names
#'
#' @returns
#'
#' @export
#' @examples
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

#' build the wide table for features while calculating the global score and breadth
#' 
#' First aggregate the feature subtype to feature type level, then calculate the global score and breadth for each feature type across all drugs/classes.
#' Then create a wide table representation 
#' 
#' @param feature_summary The tibble of summarized features across seeds generated from `summariseFeatureAcrossSeeds()`
#'
#' @returns a wide tibble with each drug/class score and global score for individual features from different scales. 
#'
#' @export
#' @examples
buildFeatureWideTable <- function(feature_summary
                              ) {
 
  id_cols = c("drug_label", "drug_or_class")
 
  adv_feat_summary <- feature_summary |>
  dplyr::group_by(
    species, drug_label, drug_or_class, feature_type, Variable
  ) |>
  dplyr::summarise(
    n_subtype = dplyr::n_distinct(feature_subtype),
    subtype_csv = paste(sort(unique(feature_subtype)), collapse = ","),
    type_mean_score = mean(mean_rank_score, na.rm = TRUE),
    type_median_rank = median(median_rank, na.rm = TRUE),
    type_rank_sd = sd(median_rank, na.rm = TRUE),
    frequency = sum(n_seeds),
    sign = if (dplyr::n_distinct(sign) == 1) dplyr::first(sign) else "MIXED",
    .groups = "drop"
  ) |>
    tidyr::unite(
      col = "model_id",
      dplyr::all_of(id_cols),
      sep = ".",
      remove = FALSE
    )
  
   row_cols = c("species", "feature_type", "Variable")

  global_summary <- adv_feat_summary |>
    dplyr::group_by(dplyr::across(dplyr::all_of(row_cols))) |>
    dplyr::summarise(
      global_breadth = dplyr::n_distinct(model_id),
      global_score = mean(type_mean_score, na.rm = TRUE),
      global_sd = sd(type_rank_sd, na.rm = TRUE),
      .groups = "drop"
    )

  wide_cols = c("sign", "frequency", "type_mean_score")

  wide_part <- adv_feat_summary |>
    dplyr::select(
      dplyr::all_of(row_cols),
      model_id,
      dplyr::all_of(wide_cols)
    ) |>
    tidyr::pivot_wider(
      names_from = model_id,
      values_from = dplyr::all_of(wide_cols),
      names_sep = "."
    )

  wide_table <- dplyr::left_join(wide_part, global_summary, by = row_cols) |>
    dplyr::arrange(dplyr::desc(global_breadth), dplyr::desc(global_score))

  return(wide_table)
}


#------------------------------------------------------------
# Cluster-wide table:
# one row per cluster, with per-model columns
# and global_cluster_score / global_cluster_breadth
#------------------------------------------------------------
#' Cluster wide table with global scores. 
#'
#' @param cluster_summary
#' @param id_cols
#' @param row_cols
#' @param score_col
#'
#' @returns
#'
#' @export
#' @examples
buildClusterWideTable <- function(feature_summary, cluster_feature_parquet) {
  
  build <- buildFeatureWideTable(
    feature_summary
  )

  cluster_feature <- arrow::read_parquet(cluster_feature_parquet)
  build |> 
    dplyr::rename(
      global_cluster_score = global_score,
      global_cluster_sd = global_sd
    ) |>
    dplyr::arrange(dplyr::desc(global_cluster_score), dplyr::desc(global_breadth))
}