
<!-- README.md is generated from README.Rmd. Please edit that file -->

# amRml: machine learning for antimicrobial resistance prediction

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/amRml)](https://CRAN.R-project.org/package=amRml)
<!-- badges: end -->

amRml is part of the [amR suite](https://github.com/JRaviLab/amR) for
antimicrobial resistance prediction. This package trains interpretable
machine learning models using genomic features from bacterial isolates.

## Installation

``` r
# Install from GitHub
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("JRaviLab/amRml")
```

## Quick start

``` r
library(amRml)

# Load feature data from a Parquet file
ml_data <- loadMLInputTibble("path/to/features.parquet")

# Run the ML pipeline
results <- runMLPipeline(
 ml_input_tibble = ml_data,
  model = "LR",
  split = c(0.6, 0.2),
  n_top_feats = 20,
  return_fit = TRUE
)

# View results
results$performance_tibble
results$top_feat_tibble
```

## Features

- **Data preparation**: Load Parquet files and prepare ML-ready datasets
- **Model training**: Logistic regression, random forest, and boosted
  trees via tidymodels
- **Evaluation**: nMCC, F1, balanced accuracy, AUPRC, and confusion
  matrices
- **Feature importance**: Extract and rank predictive features
- **Iterative feature elimination**: Identify minimal predictive feature
  sets
- **Baseline comparisons**: Fisher’s exact tests with Benjamini-Hochberg
  correction

See the [package
vignette](https://jravilab.github.io/amRml/articles/intro.html) for
detailed usage.

## Related packages

- [amR](https://github.com/JRaviLab/amR): Suite metapackage
- [amRdata](https://github.com/JRaviLab/amRdata): Data curation
- [amRshiny](https://github.com/JRaviLab/amRshiny): Interactive
  dashboard

## Citation

    Wolfe EP^, Brenner EP^, Ravi J. (2025). amRml: Machine learning for antimicrobial
    resistance prediction. R package version 0.99.0.
    https://github.com/JRaviLab/amRml

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md)
for guidelines.

## Reporting issues

Report bugs and request features at:
<https://github.com/JRaviLab/amRml/issues>

## License

BSD 3-Clause License. See [LICENSE](LICENSE) for details.

## Contact

**Corresponding author**: Janani Ravi (<janani.ravi@cuanschutz.edu>)

**Lab website**: <https://jravilab.github.io>

## Code of conduct

Please note that `amRml` is released with a [Contributor Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
