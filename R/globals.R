# Package-level documentation.

#' @keywords internal
#' @importFrom rlang .data
#' @importFrom glmnet glmnet
#' @importFrom stats median quantile sd
#' @importFrom utils head
"_PACKAGE"

# Column names used in feature rescoring's dplyr code (R/feature-rescoring.R).
utils::globalVariables(c(
  "contribution", "contribution_cutoff", "cum_contrib", "drug_label", "drug_or_class", "dyad",
  "dyad_median_contribution", "dyad_median_rank_score", "dyad_score_cutoff", "feature",
  "feature_subtype", "feature_type", "feature_type_subtype", "fit_mixture", "fit_penalty",
  "frequency",
  "importance", "Importance", "in_core", "in_core_consistent", "mcc", "MCC_diff",
  "mean_rank_score", "median_contribution", "median_cum_contrib", "median_MCC_diff",
  "median_nonshuffled_MCC", "median_rank_score", "model", "n_drug_or_class", "n_feat",
  "n_feats_returned", "n_feature_types", "n_features", "n_seeds", "nonshuffled_MCC",
  "overall_dyad_median_rank_score", "overall_median_contribution", "overall_median_rank_score",
  "rank_score", "rank_score_cutoff", "rank_score_cv", "rank_score_sd", "running_contrib", "seed",
  "seed_ratio", "shuffled", "shuffled_FALSE", "shuffled_MCC", "shuffled_TRUE", "Sign",
  "sign_consistent", "sign_consistent_across_models", "species", "subtype_csv",
  "total_frequency", "variable", "Variable", "variables_csv"
))

# Additional columns used by model reliability and network evidence summaries.
utils::globalVariables(c(
  "above_background", "annotation_key", "assessable_pair", "direct_score", "encoding",
  "evidence_status", "gene_cluster", "gene_representation_supported", "gene_robust75",
  "has_real", "has_shuffled", "high_priority", "legacy_key", "matrix_id", "mcc_adjusted",
  "mcc_for_reliability", "mcc_real", "mcc_shuffled", "mean_returned_contribution",
  "mean_shuffled_mcc", "mode_id", "n_core_seeds", "n_eligible_seeds", "n_fitted_seeds",
  "n_missing", "n_paired_assessable", "n_predicted", "n_returned_features",
  "n_returned_seeds", "n_shuffled_defined", "n_supported_seeds", "n_true", "node_id",
  "pair_reliability_supported", "pair_robust75", "predicted", "priority_under_reliability",
  "protein_cluster", "protein_representation_supported", "protein_robust75",
  "reproducibly_important", "returned", "returned_nonzero", "sd_shuffled_mcc",
  "shuffled_pair_reliability_supported", "shuffled_pair_robust75", "shuffled_upper", "status",
  "supported", "target", "task_id", "truth", "unit"
))
