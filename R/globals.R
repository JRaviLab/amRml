# Global variables declaration to avoid R CMD check NOTEs for NSE
# (non-standard evaluation) variables used with dplyr/tidyr/ggplot2

# Package-level imports for packages used with :: notation
#' @importFrom jsonlite fromJSON write_json
#' @importFrom glmnet glmnet
#' @importFrom grDevices colorRampPalette
#' @importFrom stats median
#' @importFrom stats reorder
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
  "abs_imp",
  "adj_p_value",
  "antibiotic",
  "bal_acc",
  "category",
  "class_abbr",
  "cluster",
  "contribution",
  "diff_top2",
  "drug",
  "drug_abbr",
  "drug_class",
  "drug_label",
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
  "genome_drug.antibiotic",
  "genome_drug.genome_id",
  "genome_drug.resistant_phenotype",
  "genome_id",
  "i_sparse",
  "i_strat",
  "idx_sparse",
  "idx_strat",
  "label",
  "model",
  "neg_log10_adj_p",
  "mcc",
  "mean_margin",
  "median_mcc",
  "n_feat_types",
  "nmcc",
  "num_obs",
  "seed",
  "sens",
  "spec",
  "output_prefix",
  "p_value",
  "pair_id",
  "parts",
  "phenotype",
  "proteinName",
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
  "shuffled",
  "shuffled_label",
  "sig_after_bh",
  "significance",
  "strat_value",
  "strat_value_test",
  "stratification",
  "test_country",
  "test_drug",
  "test_file",
  "test_year",
  "tested_on",
  "total",
  "train_country",
  "train_prop",
  "train_year",
  "value",

  # Base R functions that need explicit import
  "barplot",
  "reformulate"
))
