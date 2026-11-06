# MDR: modeling combinations of resistant classes

Status: both MDR targets are implemented. See "Checking against old amRml" for results.

## Goal

Two MDR models, both going through the same stages as every other target
(profile, eligibility, matrices, fitting) and running in the rewrite's modes,
not only on the whole dataset:

1. **Combination model** (`target = "resistant_classes"`): port of old amRml's
   MDR model (`.parquet2MDRMatrix()` in `R/generate_matrices_ml.R` on `main`).
   Each combination of resistant classes is one label of a multi-class model.
2. **MDR yes/no model** (`target = "mdr"`): Resistant if the genome is resistant
   to **2 or more** drug classes, Susceptible otherwise. The lab chose 2, not
   the clinical standard of 3 (Magiorakos et al., 2012).

Multi-label models (one yes/no per drug, nested in classes), random forests and
neural networks are later work and out of scope here.

## What old amRml did

- **Label:** amRdata's `resistant_classes`, the genome's resistant classes
  joined with `_` (e.g. `MAC_PEN_TMD`), or `Susceptible`. Each combination is
  one class of a multi-class model.
- **Filter:** keep combinations with at least `min_n` genomes (default 25) over
  the whole dataset, and drop `Intermediate`.
- **Features:** variance > 0 across all kept genomes, test genomes included.
  Counts and binary per feature type, binary only for struct (13 matrices).
- **Checks:** none. On S. argenteus it builds a model with a single label (`PEN`).
- **Fitting:** `multinom_reg()` with glmnet, the same tuning grid, chosen by MCC.
  Only MCC and nMCC are reported. Importance is one coefficient table per label.
- **Modes:** whole dataset only (`runMLmodels()` stops if MDR is combined with
  stratification, LOO or cross-testing).

## Data check (amRdata main, October 2026)

`resistant_classes` equals the combination recomputed from the drug-level
phenotypes (Resistant to any drug in a class means that class is resistant) for
every genome: S. flexneri 379/379, S. argenteus 73/73, E. faecium 2,907/2,907.
No genome has more than one value, and none is `Intermediate`.

With `min_n = 25`:

| Dataset | Labels kept | Genomes kept | Years / countries among kept genomes |
|---|---|---|---|
| S. flexneri | 4 | 288 of 379 | one year bin (2015–2019), one country (AUS) |
| S. argenteus | 1 (`PEN`) | 65 of 73 | one year bin, one country (THA) |
| E. faecium | 19 | 2,537 of 2,907 | several; many genomes have no year |

MDR yes/no (2+ classes). "Can't tell" means the genome was tested in fewer
than 2 classes, so it can't be called not MDR:

| Dataset | MDR | Not MDR | Can't tell |
|---|---|---|---|
| S. flexneri | 300 | 75 | 4 |
| S. argenteus | 5 | 68 | 0 |
| E. faecium | 1,413 | 563 | 931 |

So only E. faecium exercises the stratified, cross and leave-one-out modes. Its
label mix differs strongly by country (e.g. `GLY_MAC_PEN`: 90 in CAN, 0 in GBR),
which matters for the cross modes below.

## Design, by stage

The combination target is declared once, in `.MULTICLASS_TARGETS` in
`R/registry.R`, next to `.MODES`. Eligibility and fitting read it through
`.isMulticlass(unit, target)` instead of checking names or guessing from
labels. A new multi-class target only needs a row there and a profile builder.

### 1. Profile (`profileORB()`)

Add two kinds of `unit = "mdr"` rows to `profile$phenotypes`:

- `target = "mdr"` (**implemented**): each genome gets a status in the new
  `profile$mdr` table, with `n_classes_tested` and `n_classes_resistant`:
  - `mdr`: Resistant in at least `mdr_min_classes` classes (new
    `profileORB()` argument, default 2). Phenotype `Resistant`.
  - `not_mdr`: tested in at least `mdr_min_classes` classes but Resistant in
    fewer. Phenotype `Susceptible`.
  - `undetermined`: tested in fewer than `mdr_min_classes` classes, so MDR
    can't be ruled out (931 genomes in E. faecium). Kept as its own label in
    `profile$mdr`, but given no `mdr` phenotype, so it never enters the yes/no
    model. (Agreed.)

  This is a plain two-phenotype target, so every later stage handles it
  unchanged.
- `target = "resistant_classes"`: `phenotype` set to the genome's combination
  (or `Susceptible`).

For the combination target, Build it from
the class phenotypes the profile already computes, sorting class names with
radix order and joining them with `_`. That matches amRdata's string exactly (see the data check)
and keeps a single source of labels. The genome exclusions and contradiction
checks then apply to MDR too.

`targets` counts stay per target. For the combination target, `n_resistant`/`n_susceptible` don't
apply. Use `n_resistant` = genomes with any resistant class, and
`n_susceptible` = `Susceptible` genomes, so the counts still add up.

### 2. Eligibility (`eligibleScopes()`)

The `mdr` yes/no target needs no change: it is judged like a drug or class.
The rest of this section is about the combination target.

MDR uses the existing scope builders, so it gets scopes in every mode with a
grouping of `none`, `year` or `country`. The drug-grouped modes (`cross_drug`,
`lodo`) only take `unit == "drug"` rows, so they leave MDR out with no change.

New argument `min_mdr_label_genomes = 25` (old `min_n`). Per MDR scope:

1. Keep the labels with at least `min_mdr_label_genomes` **training** genomes, and
   drop the scope's genomes (train and test) with any other label. On the whole
   dataset this gives exactly old amRml's genome set. Doing it per scope keeps a
   label with 3 genomes in one year from breaking that year's folds.
2. Then apply the usual rules, with "phenotype" meaning label:
   - `single_phenotype`: fewer than 2 labels kept.
   - `single_class_labels`: the kept labels involve fewer than 2 drug classes
     (e.g. only `CLX` and `Susceptible`), so the model repeats that class's model.
     This covers a dataset tested on one class only, where the yes/no target
     also has no model (every genome is `undetermined`). (Suggested by Abhirupa.)
   - `too_few_genomes`: fewer than `min_genomes` genomes kept.
   - `too_few_for_cv`: the rarest kept label can't fill `n_fold` folds after the
     holdout. With 25 genomes per label this always passes.
   - `test_too_small`: the test set has fewer than `min_test_genomes` genomes
     with kept labels.
3. Record the dropped labels and genomes per scope, so the eligibility figure
   can show them.

The rest of eligibility assumes two phenotypes (`pmin(n_resistant, n_susceptible)`).
MDR needs its own count of the rarest label. The cleanest fix is a per-scope
`n_labels` and `n_rarest` computed for every unit; for drugs, `n_rarest` is the
existing `pmin()`.

In cross and leave-one-out modes, test genomes whose label never appears in
training can't be predicted correctly, so they are dropped from the test set
(rule 1 already does this) and the number dropped is recorded. (Agreed.)

The combination label has the same blind spot as the yes/no one: a genome
tested in only one class is labelled from that class alone. Old amRml did the
same, so the port keeps these genomes on purpose, while the yes/no target leaves
them out as `undetermined`. **Open:** should the combination target also leave
out `undetermined` genomes? That would change the results compared with old amRml.

### 3. Matrices (`planMatrices()`, `buildMatrices()`)

No change expected. The indexed store already keeps each scope's rows with
their `phenotype`, and features are chosen from training genomes only. That is a
deliberate difference from old MDR (it used all genomes), as in #51.

### 4. Fitting (`fitModels()`)

The `mdr` yes/no target fits exactly like a drug. For the combination target,
use `multinom_reg()` instead of `logistic_reg()`.
Same engine, grid and MCC selection. Changes needed in `R/fit-models.R`:

- `.modelData()`: factor levels from the scope's labels, not `.PHENOTYPES`.
- `.stratifiedFolds()`: loop over the labels present, not `.PHENOTYPES`.
  Otherwise MDR rows are never assigned a fold.
- `.splitTask()`: already splits per phenotype; no change.
- `.taskMetrics()`: for MDR, MCC and nMCC only, as old amRml did. Columns that
  only make sense for Resistant vs Susceptible are `NA`.
- `.taskImportance()`: one coefficient table per label, with a new `label`
  column (`NA` for binary models).
- Predictions: `prob_resistant` is `NA` for MDR; add `prob_predicted`.
- Performance: `res_prop` and `res_prop_test` are `NA` for MDR; add `n_labels`.
- Tuning metric set: drop `pr_auc` for MDR, since it needs one probability
  column per label.

The new columns (`label`, `prob_predicted`, `n_labels`) go in the shared
results tables, not in separate MDR tables as old amRml had. (Agreed.) They
should be included in the global rename.

### 5. Checking against old amRml

Done on S. flexneri (amRdata main, October 2026), with old amRml from #42:

1. **Matrices:** all 13 old MDR matrices (`.parquet2MDRMatrix()`, `min_n = 25`)
   are identical to the new whole-dataset combination matrices: the same 288
   genomes and labels, the same features, and the same non-zero values. In the
   whole-dataset mode every scope genome is a training genome at matrix time
   (the holdout happens at fitting), so selecting features from training genomes
   gives the old all-genome filter.
2. **Fit, AMRFinder binary, seed 123:** old `runMLPipeline()` and new
   `fitModels()` choose the same penalty (1e-4) and mixture (0.2), test on the
   same 59 genomes, predict the same label for every one, and score MCC 0.85.
   All 521 non-zero coefficients match exactly.

Still to check:

3. **S. argenteus:** old builds a 1-label MDR model (`PEN`, 65 genomes); new
   records `single_phenotype` and drops the other 8 genomes. Its yes/no target
   (5 MDR, 68 not) fails the fold check.

Still to check: E. faecium in every mode, once its ORB is ready.

## Order of work

1. Profile: `mdr` yes/no rows, tests. **Done.**
2. Profile: combination rows, tests. **Done.**
3. Eligibility: label filter, rules, tests. **Done.**
4. Fitting: multinomial path, metrics, importance, tests. **Done.**
5. Comparison with old amRml on S. flexneri. **Done.**
