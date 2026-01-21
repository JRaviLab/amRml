## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)

## ----setup--------------------------------------------------------------------
# library(amRml)

## ----load-data----------------------------------------------------------------
# # Load Parquet file and convert to wide format
# ml_data <- loadMLInputTibble("features.parquet")
# 
# # Get number of features
# n_features <- getNumFeat(ml_data)
# 
# # Get target variable name
# target_var <- getTargetVarName(ml_data)

## ----split-data---------------------------------------------------------------
# # Split data into train/validation/test sets
# # split = c(0.6, 0.2) means 60% train, 20% validation, 20% test
# data_split <- splitMLInputTibble(ml_data, split = c(0.6, 0.2), seed = 123)
# 
# # Access train and test sets
# library(rsample)
# train_data <- training(data_split)
# test_data <- testing(data_split)

## ----utils--------------------------------------------------------------------
# # Calculate minimum samples needed for valid ML
# min_samples <- calculateMinSamples(
#   n_fold = 5,
#   split = c(0.6, 0.2),
#   res_prop = 0.3
# )
# 
# # Shuffle labels for baseline comparison
# shuffled_data <- shuffleLabels(ml_data, seed = 123)

## ----run-pipeline-------------------------------------------------------------
# results <- runMLPipeline(
#   ml_input_tibble = ml_data,
#   model = "LR",                              # "LR", "RF", or "BT"
#   split = c(0.6, 0.2),                       # Train/validation proportions
#   n_fold = 2,                                # CV folds (when validation = 0)
#   prop_vi_top_feats = c(0, 1),               # Feature importance range
#   n_top_feats = NA,                          # Or specify exact number
#   use_pca = FALSE,                           # Use PCA dimensionality reduction
#   penalty_vec = 10^seq(-4, -1, length.out = 10),  # LR penalty values
#   mix_vec = 0:5 / 5,                         # LR mixture values (0=L2, 1=L1)
#   min_n_vec = c(2, 6, 12),                   # RF/BT min_n values
#   tree_vec = c(100, 500, 1000),              # RF/BT tree counts
#   select_best_metric = "mcc",                # "mcc", "f_meas", "pr_auc", "bal_accuracy"
#   seed = 123,
#   shuffle_labels = FALSE,                    # For baseline comparisons
#   return_tune_res = FALSE,
#   return_fit = FALSE,
#   return_pred = FALSE,
#   verbose = TRUE
# )
# 
# # Results
# results$performance_tibble
# results$top_feat_tibble

## ----build-models-------------------------------------------------------------
# # Build recipe (preprocessing)
# recipe <- buildRecipe(train_data, use_pca = FALSE, pca_threshold = 0.95)
# 
# # Build model specification
# lr_model <- buildLRModel(multi_class = FALSE)
# rf_model <- buildRFModel()
# bt_model <- buildBTModel()
# 
# # Build workflow
# wflow <- buildWflow(lr_model, recipe)

## ----tune-fit-----------------------------------------------------------------
# # Build tuning grid
# grid <- buildTuningGrid(
#   model = "LR",
#   penalty_vec = 10^seq(-4, -1, length.out = 10),
#   mix_vec = 0:5 / 5
# )
# 
# # For random forest or boosted tree:
# # grid <- buildTuningGrid(
# #   model = "RF",
# #   n_feat = getNumFeat(train_data),
# #   min_n_vec = c(2, 6, 12),
# #   tree_vec = c(100, 500, 1000)
# # )
# 
# # Tune hyperparameters
# tune_res <- tuneGrid(wflow, data_split, grid, n_fold = 5)
# 
# # Select and fit best model
# best_wflow <- selectBestModel(tune_res, wflow, select_best_metric = "mcc")
# fit <- fitBestModel(best_wflow, train_data)
# 
# # Get fitted hyperparameters
# fit_hps <- getFitHps(fit)

## ----predict------------------------------------------------------------------
# # Make predictions on test data
# predictions <- predict(fit, test_data)
# 
# # Get confusion matrix
# conf_mat <- getConfusionMatrix(predictions)

## ----metrics------------------------------------------------------------------
# # Individual metrics
# nmcc <- calculatenMCC(predictions)       # Normalized MCC (0-1 scale)
# f1 <- calculateF1(predictions)           # F1 score
# bal_acc <- calculateBalAcc(predictions)  # Balanced accuracy
# auprc <- calculateAUPRC(predictions)     # Area under PR curve
# log2_apop <- calculateLog2APOP(predictions)  # log2(AUPRC/prior)
# 
# # All metrics at once
# all_metrics <- calculateEvalMets(predictions)

## ----feature-importance-------------------------------------------------------
# # Extract top features from fitted model
# top_features <- extractTopFeats(
#   fit,
#   n_top_feats = 20              # Exact number
# )
# 
# # Or by proportion of variable importance
# top_features <- extractTopFeats(
#   fit,
#   prop_vi_top_feats = c(0, 0.2)  # Top 20% of variable importance
# )
# 
# # Feature importance columns: Variable, Importance, Sign
# # Positive Sign = associated with resistance
# # Negative Sign = associated with susceptibility

## ----ife----------------------------------------------------------------------
# # Run IFE by removing percentage of features
# ife_results <- runIFE(
#   ml_data,
#   by_num = TRUE,                         # Remove by number of features
#   by_vi = FALSE,                         # Or remove by variable importance
#   percent_removal_vec = 10 * 1:9,        # Remove 10%, 20%, ..., 90%
#   mix_vec = 0,                           # L2 regularization for IFE
#   return_feats = TRUE,                   # Return features removed at each step
#   verbose = TRUE
# )
# 
# # Results include nMCC at each iteration
# ife_results$ife_performance_tibble
# ife_results$feats_removed  # If return_feats = TRUE
# 
# # Remove specific features manually
# ml_data_reduced <- removeTopFeats(ml_data, top_features)

## ----plot-prc-----------------------------------------------------------------
# # Plot precision-recall curve
# prc_plot <- plotPRC(predictions)
# print(prc_plot)

## ----plot-vi------------------------------------------------------------------
# # Plot top features by variable importance
# vi_plot <- plotTopFeatsVI(fit, n_top_feats = 10)
# print(vi_plot)

## ----plot-default-eval--------------------------------------------------------
# # Plot performance metric vs training proportion or CV folds
# eval_plot <- plotDefaultEval(
#   default_eval_tibble,
#   x_default_eval = "train_prop",    # or "n_fold"
#   y_default_eval = "avg_f1_score",  # or other metrics
#   xlab = "Training data proportion",
#   ylab = "Average F1 score"
# )
# print(eval_plot)

## ----plot-baseline------------------------------------------------------------
# # Compare model performance with and without shuffled labels
# baseline_plot <- getBaselineComparisonBarplot(
#   non_shuffled_label_results = results$performance_tibble,
#   shuffled_label_results = baseline_results$performance_tibble
# )

## ----fisher-------------------------------------------------------------------
# # Complete Fisher pipeline
# fisher_results <- runFishers(
#   matrix_path = "features.parquet",
#   Q = 0.05,                              # FDR threshold
#   alternative = "two.sided",
#   susceptible_label = "Susceptible",
#   resistant_label = "Resistant"
# )
# 
# # Or step-by-step:
# ml_data <- loadMLInputTibble("features.parquet")
# encoded <- encodePhenotype(ml_data)
# fisher_res <- runFisherTests(encoded$df, encoded$target, alternative = "two.sided")
# fisher_res <- applyBenjaminiHochberg(fisher_res, Q = 0.05)
# fisher_res <- computeFeatureFreq(encoded$df, fisher_res, encoded$target)

## ----workflow-----------------------------------------------------------------
# library(amRml)
# 
# # 1. Load and prepare data
# ml_data <- loadMLInputTibble("species_drug_features.parquet")
# message("Loaded ", nrow(ml_data), " samples with ", getNumFeat(ml_data), " features")
# 
# # 2. Run ML pipeline with logistic regression
# results <- runMLPipeline(
#   ml_input_tibble = ml_data,
#   model = "LR",
#   split = c(0.6, 0.2),
#   n_top_feats = 50,
#   select_best_metric = "mcc",
#   return_fit = TRUE,
#   return_pred = TRUE
# )
# 
# # 3. View results
# print(results$performance_tibble)
# print(head(results$top_feat_tibble, 20))
# 
# # 4. Baseline comparison with shuffled labels
# baseline_results <- runMLPipeline(
#   ml_input_tibble = ml_data,
#   model = "LR",
#   split = c(0.6, 0.2),
#   shuffle_labels = TRUE
# )
# 
# # 5. Compare real vs baseline performance
# cat("Real nMCC:", results$performance_tibble$nmcc, "\n")
# cat("Baseline nMCC:", baseline_results$performance_tibble$nmcc, "\n")
# 
# # 6. Run iterative feature elimination
# ife_results <- runIFE(
#   ml_data,
#   by_num = TRUE,
#   percent_removal_vec = c(10, 25, 50, 75, 90),
#   return_feats = TRUE
# )
# 
# # 7. Generate plots
# vi_plot <- plotTopFeatsVI(results$fit, n_top_feats = 15)
# print(vi_plot)

## ----session-info, eval = TRUE------------------------------------------------
sessionInfo()

