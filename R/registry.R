# The registry of what amRml can model: the strata, and the modelling modes.

# Strata, as amRml names them, and the genome-level amRdata column holding each.
.STRATA <- c(year = "year_bin", country = "country_abbr")

# The modelling modes of the amRml v1 design. `grouping` is what a mode splits genomes
# by: a stratum, drugs, or nothing.
.MODES <- tibble::tribble(
  ~mode_id,             ~grouping, ~LOO,  ~cross_test, ~supported,
  ~description,
  "vanilla",            "none",    FALSE, FALSE,       TRUE,
  "Each drug or class, with internal cross-validation.",
  "stratified_year",    "year",    FALSE, FALSE,       TRUE,
  "Each drug or class within each year bin.",
  "stratified_country", "country", FALSE, FALSE,       TRUE,
  "Each drug or class within each country.",
  "cross_year",         "year",    FALSE, TRUE,        TRUE,
  "Trained on one year bin, tested on another.",
  "cross_country",      "country", FALSE, TRUE,        TRUE,
  "Trained on one country, tested on another.",
  "cross_drug",         "drug",    FALSE, TRUE,        TRUE,
  "Trained on one drug, tested on another drug's other genomes.",
  "loto",               "year",    TRUE,  FALSE,       TRUE,
  "Trained on all year bins but one, tested on the held-out bin.",
  "logo",               "country", TRUE,  FALSE,       TRUE,
  "Trained on all countries but one, tested on the held-out country.",
  "lodo",               "drug",    TRUE,  FALSE,       FALSE,
  "Trained on all drugs but one, tested on the held-out drug."
)

#' Modelling modes
#'
#' The kinds of model amRml can build, from the amRml v1 design. `LOO` modes
#' train on every group but one and test on the held-out group; `cross_test`
#' modes train on one group and test on another. Modes with `supported = FALSE`
#' can't run yet.
#'
#' @return A tibble with one row per mode: `mode_id`, `grouping` (`"none"`,
#'   `"year"`, `"country"` or `"drug"`), `LOO`, `cross_test`, `supported` and
#'   `description`.
#'
#' @examples
#' modelingModes()
#' @export
modelingModes <- function() {
  .MODES
}
