# Builds a small, signal-rich ML tibble used by pipeline tests. Keeping the
# dataset tiny keeps glmnet-based tuning fast while retaining enough samples
# per class for stratified splits + vfold CV.
make_pipeline_fixture <- function(n_per_class = 20, n_informative = 3,
                                  n_noise = 2, seed = 123) {
  set.seed(seed)
  n <- 2 * n_per_class
  label <- factor(
    rep(c("Resistant", "Susceptible"), each = n_per_class),
    levels = c("Resistant", "Susceptible")
  )

  informative <- lapply(seq_len(n_informative), function(i) {
    # 85% chance the feature aligns with the class label.
    ifelse(
      label == "Resistant",
      rbinom(n, 1, 0.85),
      rbinom(n, 1, 0.15)
    )
  })
  names(informative) <- paste0("feat_inf_", seq_len(n_informative))

  noise <- lapply(seq_len(n_noise), function(i) rbinom(n, 1, 0.5))
  names(noise) <- paste0("feat_noise_", seq_len(n_noise))

  tibble::as_tibble(
    c(
      list(
        genome_id = paste0("g", seq_len(n)),
        genome_drug.resistant_phenotype = label
      ),
      informative,
      noise
    )
  )
}
