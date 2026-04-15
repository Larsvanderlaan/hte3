# hte3

`hte3` provides causal machine learning tools for heterogeneous treatment
effects with a high-level wrapper API on top of the `sl3` and
`tmle3` ecosystems.

The package has two supported workflows:

- High-level workflow: use `hte_task()`, `fit_cate()`, and `fit_crr()`.
- Advanced workflow: use the lower-level `Lrnr_*` classes and
  `make_hte3_Task_tx()` directly.

Paper reproducibility is preserved separately through the frozen GitHub ref
`legacy-paper-repro`. The research scripts remain in this repository, but
they are no longer shipped in the package tarball.

The production API has three main entry points:

- `hte_task()` builds an `hte3_Task` and estimates nuisance functions.
- `fit_cate()` fits conditional average treatment effect models.
- `fit_crr()` fits conditional risk-ratio models with non-negative outcomes.

## Installation

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("Larsvanderlaan/hte3")
```

For the frozen paper workflow:

```r
source("paper_EPlearner_experiments/install_legacy_paper_repro.R")
```

## Quickstart

```r
library(hte3)
library(sl3)

data <- hte3_example_data(n = 150, seed = 1)

task <- hte_task(
  data = data,
  modifiers = c("W1", "W2", "W3"),
  confounders = c("W1", "W2", "W3"),
  treatment = "A",
  outcome = "Y",
  propensity_learner = Lrnr_mean$new(),
  outcome_learner = Lrnr_mean$new(),
  mean_learner = Lrnr_mean$new(),
  cross_fit = FALSE
)

model <- fit_cate(
  task,
  method = "dr",
  base_learner = Lrnr_mean$new(),
  cross_validate = FALSE
)

pred <- predict(model, data)
head(pred)
summary(model)
```

For real analyses, the most important tuning concepts are:

- `cross_fit` in `hte_task()`: cross-fits nuisance learners such as the
  propensity score and outcome regression.
- `cross_validate` in `fit_cate()` or `fit_crr()`: cross-validates across
  candidate HTE learners.
- `cv_control = list(V = ...)`: sets the outer HTE cross-validation folds.

`cross_fit` and `cross_validate` control different stages of the pipeline.

## Cross-Validation

The wrapper API supports three common cross-validation patterns.

### 1. Cross-validate across learner families

Pass multiple methods to `fit_cate()` or `fit_crr()`:

```r
cate_portfolio <- fit_cate(
  task,
  method = c("dr", "r", "ep"),
  base_learner = Lrnr_mean$new(),
  cross_validate = TRUE,
  cv_control = list(V = 5)
)

crr_task <- hte_task(
  data = data,
  modifiers = c("W1", "W2", "W3"),
  confounders = c("W1", "W2", "W3"),
  treatment = "A",
  outcome = "Y_binary",
  propensity_learner = Lrnr_mean$new(),
  outcome_learner = Lrnr_mean$new(),
  mean_learner = Lrnr_mean$new(),
  cross_fit = FALSE
)

crr_portfolio <- fit_crr(
  crr_task,
  method = c("ipw", "t", "ep"),
  base_learner = Lrnr_mean$new(),
  cross_validate = TRUE,
  cv_control = list(V = 5)
)
```

### 2. Cross-validate across EP basis sizes

If `method = "ep"` and `cross_validate = TRUE`, the wrapper also evaluates a
portfolio of EP learners with different sieve basis sizes and selects among
them with the appropriate selector loss.

### 3. Cross-validate from the low-level API

If you want full control over the candidate learner list, use
`cross_validate_cate()` or `cross_validate_crr()` directly:

```r
cv_fit <- cross_validate_cate(
  list(
    Lrnr_cate_DR$new(base_learner = Lrnr_mean$new()),
    Lrnr_cate_R$new(base_learner = Lrnr_mean$new()),
    Lrnr_cate_EP$new(base_learner = Lrnr_mean$new(), sieve_num_basis = 6)
  ),
  task,
  cv_control = list(V = 5)
)
```

## Notes

- Continuous-treatment CATE tasks currently support `method = "r"` only.
- CRR workflows require a non-negative outcome.
- In the examples above, `cross_fit = FALSE` keeps the code lightweight. For
  analyses beyond simple examples, nuisance cross-fitting is generally preferred.

## What Changed

- Added a high-level API for common CATE and CRR workflows.
- Hardened package metadata and dependency declarations for GitHub-first
  releases.
- Added CI, pkgdown scaffolding, tests, and a legacy reproducibility
  smoke-test entrypoint.
- Kept the research code in-repo while separating it from the package build.

## Documentation

The documentation site is organized around:

- Quickstart
- CATE workflow
- CRR workflow
- Advanced `sl3` integration
- Legacy paper reproduction

The vignettes in `vignettes/` mirror that documentation. In particular:

- [`quickstart`](vignettes/quickstart.Rmd) explains the high-level workflow and the difference between cross-fitting and cross-validation.
- [`cate`](vignettes/cate.Rmd) shows single-learner and portfolio CV CATE workflows.
- [`crr`](vignettes/crr.Rmd) shows the analogous CRR workflows.
- [`advanced-sl3`](vignettes/advanced-sl3.Rmd) shows low-level learner portfolios with `cross_validate_cate()` and `cross_validate_crr()`.
