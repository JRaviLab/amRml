# This script contains functions to check the arguments of ML functions.

#' @importFrom methods is
#' @importFrom tibble is_tibble
NULL

#' Checks the validity of file path arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param path [chr] A file path
#'
.checkArgPath <- function(path) {
  if (!is.character(path)) {
    stop("The provided file path must be a character value.")
  }

  if (!file.exists(path)) {
    stop("No file exists at the provided path.")
  }
}

#' Checks the validity of tibble arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param tibble A tibble object
#' @param ml [bool] Whether the tibble to check is for ML
#'
.checkArgTibble <- function(tibble, ml = FALSE) {
  if (!tibble::is_tibble(tibble)) {
    stop("A tibble was expected but not received.")
  }

  if (nrow(tibble) == 0 || ncol(tibble) == 0) {
    stop("The provided tibble has 0 rows and/or 0 columns.")
  }

  if (ml) {
    if (!xor(
      "genome_drug.resistant_phenotype" %in% colnames(tibble),
      "resistant_classes" %in% colnames(tibble)
    )) {
      stop(paste(
        "The input tibble must have a target variable column named",
        "either `genome_drug.resistant_phenotype` or `resistant_classes`, but",
        "not both."
      ))
    }
  }
}

#' Checks the validity of `seed` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param seed [num] For reproducible analysis
#'
.checkArgSeed <- function(seed) {
  if (!is.numeric(seed)) {
    stop("The `seed` argument can only take numeric values.")
  }
}

#' Checks the validity of `n_fold` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param n_fold [num] Number of folds of cross-validation
#'
.checkArgNFold <- function(n_fold) {
  if (!is.numeric(n_fold)) {
    stop("The `n_fold` argument can only take numeric values.")
  }

  if (n_fold < 2 || n_fold %% 1 != 0) {
    stop(
      "The `n_fold` argument must be a whole, positive number greater than 1."
    )
  }
}


#' Checks the validity of `split` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param split [num] Vector of length 2 indicating the proportion of data to
#' be designated as training and validation, respectively.
.checkArgSplit <- function(split) {
  if (!is.numeric(split) || length(split) != 2) {
    stop("The `split` argument must be a vector of length 2 with numeric values.")
  }
  if (any(is.na(split))) {
    stop("The `split` argument must not contain NA values.")
  }
  if (any(split < 0) || any(split > 1)) {
    stop("The `split` values must be in the range [0, 1].")
  }

  train_prop <- split[1]
  val_prop <- split[2]
  test_prop <- 1 - train_prop - val_prop

  # Allow either CV (c(1,0)) OR classical train/val/test with all > 0
  is_pure_cv <- (train_prop == 1 && val_prop == 0)
  is_train_val_test <- (train_prop > 0 && val_prop > 0 && test_prop > 0)

  if (!(is_pure_cv || is_train_val_test)) {
    stop(paste0(
      "Unsupported split proportions. Use either:\n",
      "  - pure CV: split = c(1, 0)\n",
      "  - classical splits (train/val/test): c(train_prop, val_prop) with train > 0, val > 0, and test = 1 - train - val > 0."
    ))
  }

  # Sum must be <= 1 for all valid cases
  if ((train_prop + val_prop) > 1) {
    stop("The sum of train and validation proportions must be <= 1 (test_prop = 1 - train - val).")
  }

  invisible(TRUE)
}


#' Checks consistency of `split` and `n_fold` to enforce matrix generation modes
#'
#' @noRd
#' @keywords internal
#' @param split \[num\] Vector length 2 c(train_prop, val_prop)
#' @param n_fold \[numeric or NULL\] Number of CV folds; NULL for classical splits
.checkArgSplitAndMode <- function(split, n_fold) {
  .checkArgSplit(split) # validates proportions and allowed patterns

  train_prop <- split[1]
  val_prop <- split[2]
  test_prop <- 1 - train_prop - val_prop

  is_pure_cv <- (train_prop == 1 && val_prop == 0)
  is_train_val_test <- (train_prop > 0 && val_prop > 0 && test_prop > 0)

  if (is_pure_cv) {
    # CV requires n_fold >= 2
    if (is.null(n_fold) || n_fold < 2) {
      stop("Pure CV requires 'n_fold' to be provided and >= 2.")
    }
    .checkArgNFold(n_fold)
    return(invisible("cv"))
  }

  if (is_train_val_test) {
    # Classical splits require NO CV
    if (!is.null(n_fold)) {
      stop("Classical splits (train/val/test) require 'n_fold = NULL' (no CV).")
    }
    return(invisible("splits"))
  }

  # Should not reach here due to .checkArgSplit(), but just in case:
  stop("Unsupported split/n_fold combination. Use pure CV (split=c(1,0), n_fold>=2) OR classical splits (val>0, test>0, n_fold=NULL).")
}


#' Checks the validity of `res_prop` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param res_prop [num] The proportion of antibiotic resistant genomes
#' phenotype labels only.
#'
.checkArgResProp <- function(res_prop) {
  if (!is.numeric(res_prop)) {
    stop("The `res_prop` argument can only take numeric values.")
  }

  if (res_prop < 0 || res_prop > 1) {
    stop("The `res_prop` argument must be [0, 1].")
  }
}

#' Checks the validity of `smallest_n_obs_rs` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param smallest_n_obs_rs [num] The smallest number of "Resistant" or
#' "Susceptible" observations allowable in the smallest subset of data used
#' during ML.
#'
.checkArgSmallestNObsRS <- function(smallest_n_obs_rs) {
  if (!is.numeric(smallest_n_obs_rs)) {
    stop("The `smallest_n_obs_rs` argument can only take numeric values.")
  }

  if (smallest_n_obs_rs < 1 || smallest_n_obs_rs %% 1 != 0) {
    stop(paste(
      "The `smallest_n_obs_rs` argument must be a whole, positive",
      "number greater than 0."
    ))
  }
}

#' Checks the validity of `use_pca` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param use_pca [bool] Set to `TRUE` to use PCA instead of all features.
#'
.checkArgUsePCA <- function(use_pca) {
  if (!is.logical(use_pca)) {
    stop(paste(
      "The `use_pca` argument can only take logical values:",
      "`TRUE` or `FALSE`."
    ))
  }
}

#' Checks the validity of `pca_threshold` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param pca_threshold [num] The proportion of total variance for which the
#' principle components account
#'
.checkArgPCAThreshold <- function(pca_threshold) {
  if (!is.numeric(pca_threshold)) {
    stop("The `pca_threshold` argument can only take numeric values.")
  }

  if (pca_threshold < 0 || pca_threshold > 1) {
    stop("The `pca_threshold` argument must be [0, 1].")
  }
}

#' Checks the validity of `multi_class` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param multi_class [bool] Whether to construct a model for multi-class
#' classification
#'
.checkArgMultiClass <- function(multi_class) {
  if (!is.logical(multi_class)) {
    stop(paste(
      "The `multi_class` argument can only take logical values:",
      "`TRUE` or `FALSE`."
    ))
  }
}

#' Checks the validity of `parsnip_mod` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param parsnip_mod A parsnip model object, such as the output of
#' `buildLRModel()` (random forest and boosted tree support planned)
#'
.checkArgParsnipMod <- function(parsnip_mod) {
  if (!is(parsnip_mod, "model_spec")) {
    stop("A `parsnip` model was expected but not received.")
  }
}

#' Checks the validity of `recipe` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param recipe A `recipe` object, such as the output of `buildRecipe()`
#'
.checkArgRecipe <- function(recipe) {
  if (!inherits(recipe, "recipe")) {
    stop("A `recipe` object was expected but not received.")
  }
}

#' Checks the validity of `model` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param model [chr]  Logistic regression ("LR"), random forest ("RF"), or
#' boosted tree ("BT")
#'
.checkArgModel <- function(model) {
  if (!(model %in% c("LR", "RF", "BT"))) {
    stop(paste(
      "The `model` argument must be one of:",
      "'LR' (logistic regression), 'RF' (random forest), 'BT' (boosted tree)."
    ))
  }
}

#' Checks the validity of `penalty_vec` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param penalty_vec \[num\] A vector containing `penalty` (regularization
#' strength) values to try (for logistic regression). It is recommended to
#' choose values in the range 10^-4 to 10^4.
#'
.checkArgPenaltyVec <- function(penalty_vec) {
  if (!is.numeric(penalty_vec)) {
    stop("The `penalty_vec` argument can only take numeric values.")
  }

  if (any(penalty_vec < 0)) {
    stop("No penalty value may be less than zero.")
  }
}

#' Checks the validity of `mix_vec` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param mix_vec [num] A vector containing `mixture` values to try for logistic
#' regression. 0 corresponds to L2 regularization; 1 corresponds to L1;
#' intermediate values correspond to elastic net.
#'
.checkArgMixVec <- function(mix_vec) {
  if (!is.numeric(mix_vec)) {
    stop("The `mix_vec` argument can only take numeric values.")
  }

  for (mix in mix_vec) {
    if (mix > 1 || mix < 0) {
      stop("All elements of `mix_vec` must be [0, 1].")
    }
  }
}

#' Checks the validity of `n_feat` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param n_feat [num] Number of features in the pangenome, used
#' to calculate `mtry` values for a subsequent grid search (for random forest or
#' boosted tree). Output of `getNumFeat()`.
#'
.checkArgNFeat <- function(n_feat) {
  if (is.na(n_feat)) {
    stop(
      "`n_feat` must be specified for random forest and boosted tree models."
    )
  }

  if (!is.numeric(n_feat)) {
    stop("The `n_feat` argument can only take numeric values.")
  }

  if (n_feat < 1) {
    stop("The `n_feat` argument must be at least 1.")
  }
}

#' Checks the validity of `min_n_vec` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param min_n_vec \[num\] A vector containing `min_n` values (the number of data
#' points in a node required for the node to be split) to try for random forest
#' or boosted tree. It is recommended to choose values in the range 1 to 100.
#'
.checkArgMinNVec <- function(min_n_vec) {
  if (!is.numeric(min_n_vec)) {
    stop("The `min_n_vec` argument can only take numeric values.")
  }

  if (any(min_n_vec < 1)) {
    stop("No `min_n` value may be less than one.")
  }
}

#' Checks the validity of `tree_vec` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param tree_vec \[num\] A vector containing values to try for the number of
#' `trees` in random forest or boosted tree. It is recommended to choose values
#' in the range 100 to 1000.
#'
.checkArgTreeVec <- function(tree_vec) {
  if (!is.numeric(tree_vec)) {
    stop("The `tree_vec` argument can only take numeric values.")
  }

  if (any(tree_vec < 1)) {
    stop("No element of `tree_vec` may be less than one.")
  }
}

#' Checks the validity of `wflow` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param wflow A `workflow` object
#'
.checkArgWflow <- function(wflow) {
  if (!inherits(wflow, "workflow")) {
    stop("The `wflow` argument can only take `workflow` objects.")
  }
}

#' Checks the validity of `data_split` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param data_split An `rsplit` object, such as the output of
#' `splitMLInputTibble()`. Provide if not using cross-validation.
#'
.checkArgDataSplit <- function(data_split) {
  if (!(class(data_split)[1] %in%
    c("initial_split", "initial_validation_split"))) {
    stop("The `data_split` argument must be an `rsplit` object.")
  }
}

#' Checks the validity of `tune_res` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param tune_res Results of grid tuning, such as the output of `tuneGrid()`
#'
.checkArgTuneRes <- function(tune_res) {
  if (!is(tune_res, "tune_results")) {
    stop("The `tune_res` argument can only take `tune_results` objects.")
  }
}

#' Checks the validity of `select_best_metric` arguments in ML
#' functions
#'
#' @noRd
#' @keywords internal
#' @param select_best_metric [chr] Metric to select best model: "f_meas",
#' "pr_auc", "mcc", or "bal_accuracy"
#'
.checkArgSelectBestMetric <- function(select_best_metric) {
  if (!(select_best_metric %in% c("f_meas", "pr_auc", "mcc", "bal_accuracy"))) {
    stop(paste(
      "The `select_best_metric` argument must be one of: 'f_meas',",
      "'pr_auc', 'mcc', 'bal_accuracy'."
    ))
  }
}

#' Checks the validity of `test_data_plus_predictions` arguments in
#' ML functions
#'
#' @noRd
#' @keywords internal
#' @param test_data_plus_predictions Test data (tibble) with an added column for
#' predicted phenotype labels
#'
.checkArgTestDataPlusPredictions <- function(test_data_plus_predictions) {
  .checkArgTibble(test_data_plus_predictions, ml = TRUE)

  if (!(".pred_class" %in% colnames(test_data_plus_predictions))) {
    stop(paste(
      "`test_data_plus_predictions` is missing a `.pred_class`",
      "column. Try using the output of `predictML()`."
    ))
  }
}

#' Checks the validity of `n_top_feats arguments` in ML functions
#'
#' @noRd
#' @keywords internal
#' @param n_top_feats [num] Number of top features to plot (`plotTopFeatsVI()`)
#' or extract (`extractTopFeats()`)
#'
.checkArgNTopFeats <- function(n_top_feats) {
  if (!is.numeric(n_top_feats)) {
    stop("The `n_top_feats` argument can only take numeric values.")
  }

  if (n_top_feats < 1) {
    stop("The `n_top_feats` argument must be at least 1.")
  }
}

#' Checks the validity of `prop_vi_top_feats` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param prop_vi_top_feats [num] A vector of length 2 with elements together
#' indicating the proportion of total variable importance the top features
#' should comprise. To get the features that contribute to the top 10 to 20% of
#' total variable importance, for example, set
#' `prop_vi_top_feats = c(0.1, 0.2)`. Returns all features by default.
#'
.checkArgPropVITopFeats <- function(prop_vi_top_feats) {
  if (!is.numeric(prop_vi_top_feats) || length(prop_vi_top_feats) != 2) {
    stop(paste(
      "The `prop_vi_top_feats` argument must be a vector of length",
      "2 with numeric values."
    ))
  }

  if (any(prop_vi_top_feats > 1) || any(prop_vi_top_feats < 0)) {
    stop("Both elements of `prop_vi_top_feats` must be [0, 1].")
  }

  if (prop_vi_top_feats[1] >= prop_vi_top_feats[2]) {
    stop(
      "The first element of `prop_vi_top_feats` must be less than the second."
    )
  }
}

#' Checks the validity of `shuffle_labels` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param shuffle_labels [bool] Set to `TRUE` to randomly shuffle AMR phenotype
#' labels for baseline comparisons.
#'
.checkArgShuffleLabels <- function(shuffle_labels) {
  if (!is.logical(shuffle_labels)) {
    stop(paste(
      "The `shuffle_labels` argument can only take logical values:",
      "`TRUE` or `FALSE`."
    ))
  }
}

#' Checks the validity of `return_tune_res` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param return_tune_res [bool] Set to `TRUE` to return tuning results.
#'
.checkArgReturnTuneRes <- function(return_tune_res) {
  if (!is.logical(return_tune_res)) {
    stop(paste(
      "The `return_tune_res` argument can only take logical values:",
      "`TRUE` or `FALSE`."
    ))
  }
}

#' Checks the validity of `return_fit` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param return_fit [bool] Set to `TRUE` to return the model `fit`.
#'
.checkArgReturnFit <- function(return_fit) {
  if (!is.logical(return_fit)) {
    stop(paste(
      "The `return_fit` argument can only take logical values:",
      "`TRUE` or `FALSE`."
    ))
  }
}

#' Checks the validity of `return_pred` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param return_pred [bool] Set to `TRUE` to return the predicted and actual
#' AMR phenotypes.
#'
.checkArgReturnPred <- function(return_pred) {
  if (!is.logical(return_pred)) {
    stop(paste(
      "The `return_pred` argument can only take logical values:",
      "`TRUE` or `FALSE`."
    ))
  }
}

#' Checks the validity of `verbose` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param verbose [bool] The function will stay quiet if set to `FALSE`.
#'
.checkArgVerbose <- function(verbose) {
  if (!is.logical(verbose)) {
    stop(paste(
      "The `verbose` argument can only take logical values:",
      "`TRUE` or `FALSE`."
    ))
  }
}

#' Checks the validity of `n_reps` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param n_reps [num] Number of iterations (and thus, random seeds) over
#' which average performance metrics and runtime will be calculated for a
#' given combination of model, train data proportion, and number of folds of
#' cross-validation
#'
.checkArgNReps <- function(n_reps) {
  if (!is.numeric(n_reps)) {
    stop("The `n_reps` argument can only take numeric values.")
  }

  if (n_reps < 1 || n_reps %% 1 != 0) {
    stop(
      "The `n_reps` argument must be a whole, positive number greater than 0."
    )
  }
}

#' Checks the validity of `by_vi` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param by_vi [bool] Set to `TRUE` if removing top features by their
#' contribution to the total variable importance.
#'
.checkArgByVI <- function(by_vi) {
  if (!is.logical(by_vi)) {
    stop(paste(
      "The `by_vi` argument can only take logical values:",
      "`TRUE` or `FALSE`."
    ))
  }
}

#' Checks the validity of `by_num` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param by_num [bool] Set to `TRUE` if removing top features as a percentage
#' of the total number of features.
#'
.checkArgByNum <- function(by_num) {
  if (!is.logical(by_num)) {
    stop(paste(
      "The `by_num` argument can only take logical values:",
      "`TRUE` or `FALSE`."
    ))
  }
}

#' Checks the validity of `percent_removal_vec` args in ML functions
#'
#' @noRd
#' @keywords internal
#' @param percent_removal_vec [num] A vector of percentages of total top
#' features removed (if `by_num = TRUE` and `by_vi = FALSE`) or percent
#' contributions of top features removed to the total variable importance (if
#' `by_vi = TRUE` and `by_num = FALSE`) at each iteration. The function will
#' automatically train with all features to start, so there is no need to
#' specify 0 in this vector.
#'
.checkArgPercentRemovalVec <- function(percent_removal_vec) {
  if (any(!is.numeric(percent_removal_vec))) {
    stop("The `percent_removal_vec` argument can only take numeric values.")
  }

  if (any(percent_removal_vec > 100) || any(percent_removal_vec < 0)) {
    stop("The `percent_removal_vec` argument must contain elements [0, 100].")
  }

  if (any(diff(percent_removal_vec) <= 0)) {
    stop(paste(
      "The elements of `percent_removal_vec` must be in ascending",
      "order with no repeats."
    ))
  }
}

#' Checks the validity of `return_feats` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param return_feats [bool] Whether to return top features from each iteration
#' of IFE.
#'
.checkArgReturnFeats <- function(return_feats) {
  if (!is.logical(return_feats)) {
    stop(paste(
      "The `return_feats` argument can only take logical values:",
      "`TRUE` or `FALSE`."
    ))
  }
}

#' Checks the validity of `x_default_eval` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param x_default_eval [chr] x value of default evaluation plot: "train_prop"
#' or "n_fold"
#'
.checkArgXDefaultEval <- function(x_default_eval) {
  if (!is.character(x_default_eval)) {
    stop("The `x_default_eval` argument can only take character values.")
  }

  if (!(x_default_eval %in% c("train_prop", "n_fold"))) {
    stop("`x_default_eval` must be one of: 'train_prop', 'n_fold'.")
  }
}

#' Checks the validity of `y_default_eval` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param y_default_eval [chr] y value of default evaluation plot. It can be
#' "avg_runtime_sec" or one of the following performance metrics:
#' "avg_f1_score", "avg_log2_apop", "avg_bal_acc", "avg_mcc", or "avg_nmcc"
#'
.checkArgYDefaultEval <- function(y_default_eval) {
  if (!is.character(y_default_eval)) {
    stop("The `y_default_eval` argument can only take character values.")
  }

  if (!(y_default_eval %in%
    c("avg_f1_score", "avg_log2_apop", "avg_bal_acc", "avg_mcc", "avg_nmcc"))
  ) {
    stop(paste(
      "`y_default_eval` must be one of:",
      "'avg_f1_score', 'avg_log2_apop', 'avg_bal_acc', 'avg_mcc', 'avg_nmcc'."
    ))
  }
}

#' Checks the validity of `xlab` and `ylab` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param xlab [chr] Label for x axis
#' @param ylab [chr] Label for y axis
#'
.checkArgXYLabs <- function(xlab, ylab) {
  if (!is.character(xlab) || !is.character(ylab)) {
    stop("The `xlab` and `ylab` arguments must be of class character.")
  }
}

#' Checks the validity of `n_drugs_per_bug` arguments in ML functions
#'
#' @noRd
#' @keywords internal
#' @param n_drugs_per_bug [num] Number of drugs to plot per bug. These will be
#' chosen based on which have the most observations in their smaller class
#' ("Resistant" or "Susceptible") to reduce noise when plotting.
#'
.checkArgNDrugsPerBug <- function(n_drugs_per_bug) {
  if (!is.numeric(n_drugs_per_bug)) {
    stop("The `n_drugs_per_bug` argument can only take numeric values.")
  }
}
