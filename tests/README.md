# Tests

Run `testthat::test_local()`.

Each file tests one file in `R/`, in the order that code runs:

| test file | tests |
|---|---|
| `test-errors.R` | `R/errors.R`: the error codes, and building and throwing errors |
| `test-orb-read.R` | `R/orb-read.R`: `readORB()`, from the path argument to the ORB object |
| `test-orb-profile.R` | `R/orb-profile.R` and `R/registry.R`: `profileORB()`, from the argument to the profile object |
| `test-orb-integration.R` | reading and profiling a real ORB (see below) |

## Fixtures

`testthat/helper-orb.R` has two main ones:

- `fxManifest()` returns an amRdata manifest as an R list, without writing
  anything. Use it for tests of the manifest record and manifest selection.
- `fxOrb()` writes that manifest plus small Parquet files and an empty DuckDB
  file to a temporary directory. Use it for `readORB()` and `profileORB()`
  tests. Pass `metadata = fxMetadataRows(...)` to profile specific genomes,
  drugs and phenotypes.


The data are synthetic and too small for modeling, however.

## A real ORB

Set `AMRML_TEST_ORB` to an amRdata output directory to also run
`test-orb-integration.R`, which reads and profiles it and checks that no
top-level file changed. CI doesn't set it, so this is a manual check.
