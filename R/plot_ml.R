# This script contains functions to generate ML-related plots for the amR
# project.

#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom dplyr select
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_path
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylim
#' @importFrom tune extract_fit_parsnip
#' @importFrom vip vip
#' @importFrom yardstick pr_curve
NULL

#' Plot a Precision-Recall Curve
#' 
#' Generates a precision-recall curve (PRC) for AMR phenotype prediction results.
#' @param test_data_plus_predictions A tibble containing test data with added
#' prediction columns, typically the output of `runMLmodels()`.
#' 
#' @return A `ggplot2` object showing the precision-recall curve.
#' 
#' @details  
#' The function uses `yardstick::pr_curve()` to compute the PR curve and then
#' visualizes it using `ggplot2`.
#' 
#' @examples
#' \dontrun{
#' test_data_plus_predictions <- readr::read_tsv(results/ML_pred/Sfl_drug_AMP_domains_binary_prediction.tsv)
#' plotPRC(test_data_plus_predictions)
#'  }
#'  
#'  @export
plotPRC <- function(test_data_plus_predictions) {
  .checkArgTestDataPlusPredictions(test_data_plus_predictions)
test_data_plus_predictions <- test_data_plus_predictions |>
dplyr::mutate(
genome_drug.resistant_phenotype = factor(
genome_drug.resistant_phenotype,
levels = c("Resistant", "Susceptible")
)
)

  prc <- yardstick::pr_curve(
    test_data_plus_predictions,
    genome_drug.resistant_phenotype, .pred_Resistant
  ) |>
    ggplot2::ggplot(ggplot2::aes(x = recall, y = precision)) +
    ggplot2::geom_path() +
    ggplot2::ylim(0, 1) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())

  return(prc)
}

#' Plot a Receiver Operating Characteristic (ROC) Curve
#' 
#' Generates a ROC curve for AMR phenotype prediction results.
#' 
#' @param test_data_plus_predictions A tibble with test data and prediction
#' columns (output of `runMLmodels()`).
#' 
#' @return A ROC curve plotted using `ggplot2::autoplot()`.
#' 
#' @export
plotROC <- function(test_data_plus_predictions) {
  .checkArgTestDataPlusPredictions(test_data_plus_predictions)
  test_data_plus_predictions <- test_data_plus_predictions |>
dplyr::mutate(
genome_drug.resistant_phenotype = factor(
genome_drug.resistant_phenotype,
levels = c("Resistant", "Susceptible")
)
)

  roc <- yardstick::roc_curve(
    test_data_plus_predictions,
    genome_drug.resistant_phenotype, .pred_Resistant
  ) |>
    ggplot2::autoplot(type = "se") +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  
  return(roc)
}

#' Plot a Confusion Matrix Heatmap
#' 
#' Produces a heatmap visualization of the confusion matrix for AMR predictions.
#' 
#' @param test_data_plus_predictions A tibble containing true and predicted
#' phenotype labels.
#' 
#' @return A heatmap (`ggplot2` object) showing the confusion matrix.
#' 
#' @export
plotCM <- function(test_data_plus_predictions) {
  .checkArgTestDataPlusPredictions(test_data_plus_predictions)
  test_data_plus_predictions <- test_data_plus_predictions |>
    dplyr::mutate(
      genome_drug.resistant_phenotype = factor(
        genome_drug.resistant_phenotype,
        levels = c("Resistant", "Susceptible")
      ),
      .pred_class = factor(
        .pred_class,
        levels = c("Resistant", "Susceptible")
      )
    )
  test_data_plus_predictions |>
yardstick::conf_mat(truth = genome_drug.resistant_phenotype,
 estimate = .pred_class) |>
ggplot2::autoplot(type = "heatmap")
}

#' Plot Density of Predicted Class Probabilities
#' 
#' Visualizes how predicted class probabilities differ between resistant and
#' susceptible genome-drug combinations.
#' 
#' @param test_data_plus_predictions Tibble with prediction probabilities and
#' true labels.
#' 
#' @return A ggplot2 density plot.
#' 
#' @export
plotDensity <- function(test_data_plus_predictions) {
  test_data_plus_predictions |>
ggplot2::ggplot(ggplot2::aes(x = .pred_Resistant,
fill = genome_drug.resistant_phenotype)) +
ggplot2::geom_density(alpha = 0.5)
}

#' Plot Top Feature Importances
#' 
#' Creates a bar plot showing the most important features affecting
#' AMR phenotype predictions.
#' 
#' @param topfeat A tibble containing feature importance scores
#' (output of `runMLmodels()`).
#' @param n_top_feats Number of top features to display (default: 10).
#' 
#' @return A bar plot of variable importance (`ggplot2` object).
#' 
#' @examples
#' \dontrun{
#' topfeat <- readr::read_tsv(results/ML_top_features/Sfl_drug_AMP_domains_binary_top_features.tsv)
#' plotTopFeatsVI(topfeat)
#'  }
#' 
#' @export
plotTopFeatsVI <- function(topfeat, n_top_feats = 10) {
  .checkArgNTopFeats(n_top_feats)

  vip <- topfeat |> 
    dplyr::slice_max(order_by = Importance, n = n_top_feats) |>
    dplyr::mutate(
      Variable = factor(Variable, levels = rev(Variable)),   # preserve order as shown in table
      Sign = factor(Sign, levels = c("POS", "NEG"))
    ) |> 
    ggplot2::ggplot(ggplot2::aes(x = Importance, y = Variable, fill = Sign)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(
      values = c(
        "POS" = "#c6d8d3",  
        "NEG" = "#f6c9a1"   
      )
    ) +
    ggplot2::labs(
      x = "Importance",
      y = "Features"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 10)
    )
  
 return(vip) 
}

#' Compare Baseline Performance With and Without Shuffled Labels
#' 
#' Produces a bar plot comparing balanced accuracy for each antibiotic using
#' true AMR labels vs. randomly shuffled labels.
#' 
#' @param non_shuffled_label_results Output of `runMLPipeline(shuffle_labels = FALSE)`
#' @param shuffled_label_results Output of `runMLPipeline(shuffle_labels = TRUE)`
#' 
#' @return A base R barplot comparing balanced accuracy across models.
#' 
#' @export
getBaselineComparisonBarplot <- function(
  non_shuffled_label_results,
  shuffled_label_results
) {
  .checkArgTibble(non_shuffled_label_results$performance_tibble)
  .checkArgTibble(shuffled_label_results$performance_tibble)

  non_shuffled_bal_acc <- non_shuffled_label_results$performance_tibble |>
    dplyr::select(bal_acc) |>
    dplyr::pull()

  shuffled_bal_acc <- shuffled_label_results$performance_tibble |>
    dplyr::select(bal_acc) |>
    dplyr::pull()

  bal_acc_matrix <- matrix(c(non_shuffled_bal_acc, shuffled_bal_acc),
    nrow = 2, byrow = TRUE
  )

  rownames(bal_acc_matrix) <- c("Non-Shuffled Labels", "Shuffled Labels")

  baseline_comparison_barplot <- barplot(bal_acc_matrix,
    beside = TRUE,
    legend.text = TRUE, col = c("skyblue", "lightpink"),
    ylab = "Balanced Accuracy"
  )

  return(baseline_comparison_barplot)
}
