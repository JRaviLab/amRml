# Shared fixture helpers for unit tests.

make_binary_ml_tibble <- function(n_genomes = 10, n_feats = 4, seed = 1) {
  set.seed(seed)
  feats <- matrix(
    sample(c(0L, 1L), n_genomes * n_feats, replace = TRUE),
    nrow = n_genomes,
    ncol = n_feats
  )
  colnames(feats) <- paste0("feat_", seq_len(n_feats))
  tibble::as_tibble(
    cbind(
      data.frame(
        genome_id = paste0("g", seq_len(n_genomes)),
        genome_drug.resistant_phenotype = factor(
          rep(c("Resistant", "Susceptible"), length.out = n_genomes),
          levels = c("Resistant", "Susceptible")
        ),
        stringsAsFactors = FALSE
      ),
      feats
    )
  )
}

make_multiclass_ml_tibble <- function() {
  tibble::tibble(
    genome_id = c("g1", "g2", "g3"),
    resistant_classes = factor(c("A", "B", "A")),
    feat_1 = c(1L, 0L, 1L),
    feat_2 = c(0L, 1L, 1L)
  )
}
