# Package-level documentation.

#' @keywords internal
#' @importFrom rlang .data
#' @importFrom glmnet glmnet
#' @importFrom stats median quantile sd
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
