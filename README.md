
<!-- README.md is generated from README.Rmd. Please edit that file -->

## amRml: Machine Learning for Antimicrobial Resistance Prediction

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

amRml is part of the [amR suite](https://github.com/JRaviLab/amR) for
antimicrobial resistance prediction. amRml is the machine learning (ML)
engine of the <https://github.com/JRaviLab/amR>. It trains interpretable
ML models to predict antimicrobial resistance (AMR) from the multi‑scale
genomic features prepared by amRdata. The package includes streamlined
pipelines for:

- Single‑drug AMR prediction
- Multi‑drug resistance (MDR) prediction
- Stratified and cross‑testing workflows (year, country)
- Leave‑one‑out (LOO) generalization tests
- Baseline shuffled‑label controls
- PCA and variable‑importance‑based feature reduction

All with reproducible, parallelized model execution.

amRml produces ML‑ready tibbles, trains logistic regression models, and
exports clean performance summaries, feature rankings, predictions, and
tuned model objects.

## Overview

amRml provides functions to:

- Convert multi‑scale genomic feature Parquets into ML‑ready tibbles
- Train interpretable AMR classifiers with tidymodels
- Perform robust evaluation using MCC, F1, AuPRC, balanced accuracy, and
  confusion matrices
- Extract top predictive features to connect models’ performance with
  biological underpinnings

The models are optimized for reproducibility, interpretability, and
benchmarking across species and antibiotic classes. See the
<https://jravilab.github.io/amRml/articles/intro.html> for full
examples.

## Installation

``` r
# Install from GitHub
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("JRaviLab/amRml")
```

## Quick start

### 1. Generate ML‑ready data and metadata inputs

`generateMLInputs()` converts the Parquet‑backed DuckDB file created by
amRdata into a standard set of ready-to-model matrices.

``` r
library(amRml)

generateMLInputs(
  parquet_duckdb_path = "data/Shigella_flexneri/Sfl_parquet.duckdb",
  out_path            = "results/Shigella_flexneri/",
  n_fold              = 5,
  verbosity           = "minimal"
)
```

This will analyze the Parquet tables attached to the
{Bug}\_parquet.duckdb database created by the amRdata::cleanData() step
at the end of the `amRdata::runDataProcessing()` workflow. For each drug
or drug class, it evaluates how many isolates have paired genotype and
AMR phenotype data, and produced subsetted ML-ready matrices for each
feature scale.

The result: All modeling matrices that can be created for a given
dataset, ready for input into the ML modeling workflow next.

### 2. Run standard AMR prediction models

`runMLmodels()` executes logistic regression ML models across all
prepared matrices.

``` r
runMLmodels(
  path        = "results/Shigella_flexneri/",
  threads     = 16
)
```

This produces 4 output directories: - `performance`: Containing
performance metrics per bug-drug-feature_scale model - `pred`:
Containing the specific predictions made per model - `top_features`:
Containing which specific features drove model predictions per model -
`models`: Containing the model fits themselves in .Rds format

These models can also be adjusted with parameters like
`shuffle_labels = TRUE` (creates random baseline models to compare
against), `use_pca = TRUE` (to reduce feature space by using PCs as
features), `stratify_by = "country"` (to stratify samples into countries
of origin to identify regional trends in AMR and how well models
generalize), and many other options.

### 3. Multi‑drug resistance (MDR) modeling

`runMDRmodels()` trains models that predict aggregated multi‑drug
resistance phenotypes.

``` r
runMDRmodels(
  path        = "results/Shigella_flexneri/",
  threads     = 16
)
```

This uses specific matrices to test whether ML models can predict
resistance against multiple drug classes, and identify any features
associated with MDR.

## Features

- **Data preparation**: Load Parquet files and prepare ML-ready datasets
- **Model training**: User-customizable logistic regression via
  tidymodels
- **Evaluation**: MCC, nMCC, F1, balanced accuracy, AuPRC, and confusion
  matrices
- **Feature importance**: Extract and rank predictive features

See the [package
vignette](https://jravilab.github.io/amRml/articles/intro.html) for
detailed usage.

## Integration with amR suite

amRml is designed to work seamlessly with other amR packages:

``` r
library(amRdata)
library(amRml)
library(amRviz)

# 1. Curate data
prepareGenomes("Shigella flexneri")
runDataProcessing("amRdata/data/Shigella_flexneri/Sfl.duckdb")

# 2. Train models
runMLmodels("amRdata/data/Shigella_flexneri/Sfl_parquet.duckdb")

# 3. Visualize
launchAMRDashboard()
```

## Related packages

- [amR](https://github.com/JRaviLab/amR): Suite metapackage
- [amRdata](https://github.com/JRaviLab/amRdata): Data preparation for
  AMR prediction
- [amRviz](https://github.com/JRaviLab/amRviz): Interactive dashboard

## Citation

If you use `amRml` in your research, please cite:

> Ghosh A^, Brenner EP^, Boyer EA, McKim AP, Vang CK, Wolfe EP, Mayer D,
> Lesiyon RL, Ravi J.
>
> amR: an R package suite to predict antimicrobial resistance in
> bacterial pathogens.
>
> bioRxiv. 2026. DOI:
> [10.64898/2026.07.10.734579](https://doi.org/10.64898/2026.07.10.734579).

^ Co-first authors

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
