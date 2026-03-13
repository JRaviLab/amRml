computeFeatureImprovement <- function(
  all_feature_parquet,
  feature_cluster_parquet
) {
  stopifnot(file.exists(all_feature_parquet))
feature_cluster <- arrow::read_parquet(normalizePath(feature_cluster_parquet))
    
    
 features_rescored <- arrow::read_parquet(normalizePath(all_feature_parquet)) |>
    dplyr::select(
      output_prefix, 
      drug_label, drug_or_class, shuffled, pca,
      feature_type, feature_subtype, Variable,
      Importance, Sign
    ) |>
    dplyr::filter(!pca) |>
    dplyr::mutate(
      Variable = dplyr::case_when(
        feature_type == "domains"  ~ sub("_.+$", "", Variable),
        feature_type == "proteins" ~ sub("fig.", "fig|", Variable, fixed = TRUE),
        TRUE ~ Variable
      )
    ) |>
    dplyr::group_by(output_prefix) |>
    dplyr::mutate(
      rank  = dplyr::dense_rank(dplyr::desc(Importance)),
      denom = sum(Importance, na.rm = TRUE),
      contribution = dplyr::if_else(denom > 0, Importance / denom, 0),

      # Safer rescaling within each output_prefix
      min_imp   = min(Importance, na.rm = TRUE),
      max_imp   = max(Importance, na.rm = TRUE),
      range_imp = max_imp - min_imp,
      rescaled  = dplyr::if_else(range_imp > 0, (Importance - min_imp) / range_imp, 0)
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
      shuffled_rank     = dplyr::if_else(shuffled,  median_drug_or_class, NA_real_),
      non_shuffled_rank = dplyr::if_else(!shuffled, median_drug_or_class, NA_real_),

      # If non-shuffled is missing -> NA (no evidence). If shuffled missing -> +Inf improvement.
      improvement = dplyr::case_when(
        !is.na(non_shuffled_rank) ~ tidyr::replace_na(shuffled_rank, Inf) - non_shuffled_rank,
        TRUE                      ~ NA_real_
      ),

      # "Good" if non-shuffled exists AND (non-shuffled < shuffled OR shuffled is missing)
      good_feature = !is.na(non_shuffled_rank) & improvement > 0
    ) |>
    dplyr::ungroup() 
    
    return(features_rescored)
    
    }
  
