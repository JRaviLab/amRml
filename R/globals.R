# Global variables declaration to avoid R CMD check NOTEs for NSE
# (non-standard evaluation) variables used with dplyr/tidyr/ggplot2

# Package-level imports for packages used with :: notation
#' @importFrom jsonlite fromJSON write_json
#' @importFrom glmnet glmnet
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c(
  # Prediction columns from tidymodels
  ".estimate",
  ".pred_Resistant",
  ".pred_Susceptible",
  ".pred_class",

  # Variable importance columns
  "Importance",
  "Sign",
  "Variable",
  "coefficient",
  "cum_imp",

  # Data columns
  "adj_p_value",
  "antibiotic",
  "bal_acc",
  "drug_or_class",
  "feature",
  "feature_id",
  "feature_subtype",
  "feature_type",
  "fit_min_n",
  "fit_mixture",
  "fit_mtry",
  "fit_trees",
  "gene",
  "genome_drug.genome_id",
  "genome_drug.resistant_phenotype",
  "genome_id",
  "i_sparse",
  "i_strat",
  "idx_sparse",
  "idx_strat",
  "model",
  "neg_log10_adj_p",
  "nmcc",
  "num_obs",
  "output_prefix",
  "p_value",
  "pair_id",
  "parts",
  "phenotype",
  "precision",
  "prefix",
  "prefix_key",
  "prop",
  "recall",
  "ref_drug",
  "ref_file",
  "res_prop",
  "resistant_classes",
  "run_time_sec",
  "sig_after_bh",
  "significance",
  "strat_value",
  "strat_value_test",
  "stratification",
  "test_drug",
  "test_file",
  "train_prop",
  "value",

  # Base R functions that need explicit import
  "barplot",
  "reformulate"
))
