### Incorporating Fisher tests with BH correction as a
### baseline comparison against our classical ML models.

#' @importFrom dplyr arrange
#' @importFrom dplyr case_when
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom dplyr row_number
#' @importFrom dplyr rowwise
#' @importFrom dplyr select
#' @importFrom dplyr summarise
#' @importFrom dplyr ungroup
#' @importFrom stats fisher.test
#' @importFrom tibble tibble
#' @importFrom tidyr pivot_longer
NULL

#' encodePhenotype
#' @param df A tibble with AMR column `genome_drug.resistant_phenotype`
#' @param susceptible_label Label for "Susceptible" phenotype
#' @param resistant_label Label for "Resistant" phenotype
#'
#' @return A list with:
#'   - `df`: the input data with AMR phenotype now encoded as integer
#'   - `target`: a binarized vector of the encoded AMR phenotypes
#' @export
encodePhenotype <- function(df, susceptible_label = "Susceptible", resistant_label = "Resistant") {
  stopifnot("genome_drug.resistant_phenotype" %in% colnames(df)) # Throw a fit if it's the wrong matrix

  df <- df |>
    dplyr::filter(genome_drug.resistant_phenotype %in% c(susceptible_label, resistant_label)) |>
    dplyr::mutate(
      genome_drug.resistant_phenotype = dplyr::case_when(
        genome_drug.resistant_phenotype == susceptible_label ~ 0L,
        genome_drug.resistant_phenotype == resistant_label ~ 1L
      )
    )

  target <- df$genome_drug.resistant_phenotype
  return(list(df = df, target = target))
}

### I think we want two-sided here since genes could predispose in either direction
### but open to feedback from the team on that!
#' runFisherTests
#' @param df Binary feature columns and a binary target column
#' @param target The above-mentioned binarized vector of AMR phenotypes
#' @param alternative Type of Fisher test: "two.sided", "greater", or "less"
#'
#' @return A tibble with columns: gene, p_value
#' @export
runFisherTests <- function(df, target, alternative = "two.sided") {
  feature_cols <- df |>
    dplyr::select(-genome_id, -genome_drug.resistant_phenotype)

  fisher_test_on_column <- function(column) {
    column <- factor(column, levels = c(0, 1))
    result <- stats::fisher.test(x = column, y = target, alternative = alternative)
    return(result$p.value)
  }

  df_fisher <- feature_cols |>
    dplyr::summarise(across(everything(), fisher_test_on_column)) |>
    tidyr::pivot_longer(cols = everything(), names_to = "gene", values_to = "p_value") |>
    dplyr::arrange(p_value)

  df_fisher <- df_fisher |>
    dplyr::mutate(alternative = alternative)

  return(df_fisher)
}


#' applyBenjaminiHochberg
#' @param df_fisher The prior Fisher result data with gene & p_value cols
#' @param Q The FDR threshold (default = 0.05)
#'
#' @return Returns Fisher results with added cols for adj_p_value, sig_after_bh, Q
#' @export
applyBenjaminiHochberg <- function(df_fisher, Q = 0.05) {
  stopifnot("p_value" %in% colnames(df_fisher))

  # Slather on the BH corrections
  bh_result <- sgof::BH(u = df_fisher$p_value, alpha = Q)

  df_fisher <- df_fisher |>
    dplyr::mutate(
      adj_p_value = bh_result$Adjusted.pvalues,
      sig_after_bh = dplyr::row_number() <= bh_result$Rejections,
      Q = Q
    ) |>
    dplyr::relocate(adj_p_value, .after = p_value) |> # Order the new cols
    dplyr::relocate(sig_after_bh, .after = adj_p_value)

  return(df_fisher)
}

#' computeFeatureFreq
#' @param df The input matrix with binary features and the phenotype column
#' @param df_fisher The post-BH fisher results with gene names to annotate
#' @param target The binarized AMR phenotype "target" output from before
#'
#' @return The BH Fisher data table with added feature frequency columns
#' @export
computeFeatureFreq <- function(df, df_fisher, target) {
  df_features <- df |>
    dplyr::select(-genome_id, -genome_drug.resistant_phenotype)

  gene_pres_given_target <- function(gene, target_val) {
    gene_col <- df_features[[gene]]
    df_gene_target <- tibble::tibble(gene = gene_col, target = target)
    n_target <- sum(target == target_val)
    n_gene_present <- sum(df_gene_target$gene == 1 & target == target_val)
    return(n_gene_present / n_target)
  }

  df_fisher <- df_fisher |>
    dplyr::rowwise() |>
    dplyr::mutate(
      freq_susceptible_gene_pres = gene_pres_given_target(gene, 0L),
      freq_resistant_gene_pres = gene_pres_given_target(gene, 1L)
    ) |>
    dplyr::ungroup()

  return(df_fisher)
}

### The function that runs it all
#' runFishers
#' @param matrix_path Path to the long sparse Parquet matrices
#' @param Q FDR cutoff for Benjamini Hochberg correction (default = 0.05)
#' @param alternative Type of Fisher test: "two.sided", "greater", or "less"
#' @param susceptible_label "Susceptible" phenotype label in the data
#' @param resistant_label "Resistant" phenotype label in the data
#'
#' @return The Fisher test results, BH correction, and feature frequencies
#' @export
runFishers <- function(
  matrix_path,
  Q = 0.05,
  alternative = "two.sided",
  susceptible_label = "Susceptible",
  resistant_label = "Resistant"
) {
  # Step 1: Load and pivot the long matrix to wide
  df_wide <- loadMLInputTibble(matrix_path)

  # Step 2: Encode phenotype from factorized characters to binary
  encoded <- encodePhenotype(df_wide, susceptible_label, resistant_label)
  df <- encoded$df
  target <- encoded$target

  # Step 3: Run basic Fisher tests
  df_fisher <- runFisherTests(df, target, alternative)

  # Step 4: Apply BH correction to the Fisher results
  df_fisher <- applyBenjaminiHochberg(df_fisher, Q)

  # Step 5: Compute the feature frequencies across the phenotype classes
  df_fisher <- computeFeatureFreq(df, df_fisher, target)

  return(df_fisher)
}
