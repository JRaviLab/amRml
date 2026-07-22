computeFeatureImprovement <- function(
  all_top_features_parquet,
  cluster_feature_parquet
) {
  stopifnot(file.exists(all_top_features_parquet))
  cluster_feature <- arrow::read_parquet(normalizePath(cluster_feature_parquet))


  features_rescored <- arrow::read_parquet(normalizePath(all_top_features_parquet)) |>
    dplyr::select(
      species, drug_label, drug_or_class, shuffled,
      feature_type, feature_subtype, Variable,
      Importance, Sign
    ) |>
    dplyr::filter(!pca) |>
    dplyr::mutate(
      Variable = dplyr::case_when(
        feature_type == "domains" ~ sub("_.+$", "", Variable),
        feature_type == "proteins" ~ sub("fig.", "fig|", Variable, fixed = TRUE),
        TRUE ~ Variable
      )
    ) |>
    dplyr::group_by(output_prefix) |>
    dplyr::mutate(
      rank = dplyr::dense_rank(dplyr::desc(Importance)),
      denom = sum(Importance, na.rm = TRUE),
      contribution = dplyr::if_else(denom > 0, Importance / denom, 0),

      # Safer rescaling within each output_prefix
      min_imp = min(Importance, na.rm = TRUE),
      max_imp = max(Importance, na.rm = TRUE),
      range_imp = max_imp - min_imp,
      rescaled = dplyr::if_else(range_imp > 0, (Importance - min_imp) / range_imp, 0)
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(drug_label, drug_or_class, feature_type, feature_subtype, shuffled, Variable) |>
    dplyr::mutate(median_datatype = stats::median(rank, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::group_by(drug_label, drug_or_class, feature_type, shuffled, Variable) |>
    dplyr::mutate(median_scale = stats::median(median_datatype, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    # NOTE: correct join syntax (no quotes in join_by)
    dplyr::left_join(feature_cluster, by = dplyr::join_by(Variable == feature)) |>
    dplyr::group_by(drug_label, drug_or_class, shuffled, cluster) |>
    dplyr::mutate(
      median_drug_or_class       = stats::median(rank, na.rm = TRUE),
      count_scales_for_cluster   = dplyr::n_distinct(feature_type),
      feature_types_csv          = paste(sort(unique(feature_type)), collapse = ","),
      lowest_contri              = min(contribution, na.rm = TRUE),
      highest_contri             = max(contribution, na.rm = TRUE)
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(drug_label, drug_or_class, cluster) |>
    dplyr::mutate(
      shuffled_rank = dplyr::if_else(shuffled, median_drug_or_class, NA_real_),
      non_shuffled_rank = dplyr::if_else(!shuffled, median_drug_or_class, NA_real_),

      # If non-shuffled is missing -> NA (no evidence). If shuffled missing -> +Inf improvement.
      improvement = dplyr::case_when(
        !is.na(non_shuffled_rank) ~ tidyr::replace_na(shuffled_rank, Inf) - non_shuffled_rank,
        TRUE ~ NA_real_
      ),

      # "Good" if non-shuffled exists AND (non-shuffled < shuffled OR shuffled is missing)
      good_feature = !is.na(non_shuffled_rank) & improvement > 0
    ) |>
    dplyr::ungroup()

  return(features_rescored)
}

computeFeatureScore <- function(
  all_feature_parquet,
  feature_cluster_parquet
) {
  top_features <- arrow::read_parquet(normalizePath(all_feature_parquet))
  top_features <- top_features |>
    dplyr::mutate(
      shuffled = stringr::str_detect(prefix_key, "^shuffled_"),
      species = prefix_key |>
        stringr::str_remove("^shuffled_") |>
        stringr::str_remove("(_drug_class|_drug).*")
    )

  ## --------------------------------
  ## 1. Scale capability (data-driven)
  ## --------------------------------

  scale_capability <- top_features |>
    dplyr::distinct(feature_type, feature_subtype) |>
    dplyr::group_by(feature_type) |>
    dplyr::summarise(
      expected_types = dplyr::n_distinct(feature_subtype),
      expected_types_csv = paste(sort(feature_subtype), collapse = ","),
      .groups = "drop"
    )

  ## --------------------------------
  ## 2. Importance → contribution → rank
  ## --------------------------------

  ranked_features <- top_features |>
    dplyr::group_by(
      species,
      drug_label,
      drug_or_class,
      feature_type,
      feature_subtype,
      shuffled
    ) |>
    dplyr::mutate(
      contribution = Importance / sum(Importance, na.rm = TRUE),
      rank = dplyr::dense_rank(dplyr::desc(contribution)),
      max_rank = max(rank),
      min_rank = min(rank),
      mean_rank = mean(rank),
      rank_score = 1 - (mean_rank - 1) / max_rank,
    ) |>
    dplyr::ungroup() |>
    dplyr::select(
      species,
      drug_label,
      drug_or_class,
      feature_type,
      feature_subtype,
      shuffled,
      Variable,
      contribution,
      rank,
      max_rank,
      min_rank,
      mean_rank
    )

  ## --------------------------------
  ## 3. Collapse to protein level
  ## --------------------------------

  feature_cluster <- arrow::read_parquet(normalizePath(feature_cluster_parquet))

  ranked_features <- ranked_features |>
    dplyr::mutate(
      Variable = dplyr::case_when(
        feature_type == "domains" ~ sub("_.+$", "", Variable),
        feature_type == "proteins" ~ sub("fig.", "fig|", Variable, fixed = TRUE),
        TRUE ~ Variable
      )
    ) |>
    dplyr::left_join(feature_cluster, by = dplyr::join_by(Variable == feature))

  protein_scale_realization <- ranked_features |>
    dplyr::filter(!is.na(cluster), shuffled == FALSE) |>
    dplyr::group_by(
      species,
      drug_label,
      drug_or_class,
      cluster,
      feature_type
    ) |>
    dplyr::summarise(
      observed_types = dplyr::n_distinct(feature_subtype),
      observed_types_csv = paste(sort(unique(feature_subtype)), collapse = ","),
      .groups = "drop"
    ) |>
    dplyr::left_join(scale_capability, by = "feature_type") |>
    dplyr::mutate(
      scale_realization = dplyr::case_when(
        observed_types == expected_types ~ "full_realization",
        observed_types < expected_types ~ "partial_realization"
      )
    )

  protein_scale_summary <- protein_scale_realization |>
    dplyr::group_by(species, drug_label, drug_or_class, cluster) |>
    dplyr::summarise(
      n_scales = dplyr::n_distinct(feature_type),
      fully_realized_scales =
        sum(scale_realization == "full_realization"),
      partially_realized_scales =
        sum(scale_realization == "partial_realization"),
      scale_support_csv =
        paste(
          feature_type,
          "(", observed_types_csv, "/", expected_types, ")",
          collapse = "; "
        ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      realization_score =
        (fully_realized_scales +
          0.5 * partially_realized_scales) /
          (fully_realized_scales + partially_realized_scales)
    ) |>
    dplyr::mutate(
      coverage_boost =
        n_scales / max(n_scales)
    ) |>
    dplyr::mutate(
      scale_factor = realization_score * coverage_boost
    )

  ## --------------------------------
  ## 4. Shuffle vs non-shuffle
  ## --------------------------------

  shuffle_delta <- ranked_features |>
    dplyr::filter(!is.na(cluster)) |>
    dplyr::group_by(
      species, drug_label,
      drug_or_class,
      cluster
    ) |>
    dplyr::summarise(
      mean_rank_nonshuffle = mean(mean_rank[shuffled == FALSE], na.rm = TRUE),
      mean_rank_shuffle = mean(mean_rank[shuffled == TRUE], na.rm = TRUE),
      delta_rank = mean_rank_shuffle - mean_rank_nonshuffle,
      # If non-shuffled is missing -> NA (no evidence). If shuffled missing -> +Inf improvement.
      improvement = dplyr::case_when(
        !is.na(mean_rank_nonshuffle) ~ tidyr::replace_na(mean_rank_shuffle, Inf) - mean_rank_nonshuffle,
        TRUE ~ NA_real_
      ),

      # "Good" if non-shuffled exists AND (non-shuffled < shuffled OR shuffled is missing)
      good_feature = !is.na(mean_rank_nonshuffle) & improvement > 0,
      .groups = "drop"
    )

  ## --------------------------------
  ## 5. Cluster scoring
  ## --------------------------------
  robustness <- ranked_features |>
    dplyr::filter(!is.na(cluster), shuffled == FALSE) |>
    dplyr::group_by(species, drug_label, drug_or_class, cluster) |>
    dplyr::summarize(
      median_contribution = median(contribution),
      median_rank = median(mean_rank), .groups = "drop"
    ) |>
    dplyr::left_join(
      protein_scale_summary |>
        dplyr::distinct(species, drug_label, drug_or_class, cluster, n_scales, scale_factor),
      by = c("species", "drug_label", "drug_or_class", "cluster")
    ) |>
    dplyr::left_join(
      shuffle_delta |>
        dplyr::distinct(species, drug_label, drug_or_class, cluster, delta_rank, good_feature),
      by = c("species", "drug_label", "drug_or_class", "cluster")
    ) |>
    dplyr::group_by(species, drug_label, drug_or_class) |>
    dplyr::mutate(
      # higher = better
      contrib_score = dplyr::percent_rank(median_contribution),

      # lower rank = better → invert
      stability_score = 1 - dplyr::percent_rank(median_rank),

      # robustness: handle NaN = strongest case (missing in shuffle)
      delta_score = dplyr::case_when(
        is.nan(delta_rank) ~ 1, # best possible signal
        delta_rank > 0 ~ dplyr::percent_rank(delta_rank), # reward
        delta_rank == 0 ~ 0, # neutral
        delta_rank < 0 ~ -dplyr::percent_rank(abs(delta_rank)) # penalize
      ),

      # already bounded [0,1]
      scale_score = scale_factor,
      robustness_score =
        contrib_score *
          stability_score *
          scale_factor *
          delta_score
    )

  return(robustness)
}
