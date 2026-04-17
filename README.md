# hte3

`hte3` provides R workflows for estimating conditional average treatment
effects (CATE) and conditional risk ratios (CRR). Use the high-level wrappers
for standard analyses, the GRF wrappers for forest-based fits, and the
lower-level `sl3` interface when you need direct control over nuisance
learners or candidate HTE models.

## Installation

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("Larsvanderlaan/hte3", dependencies = TRUE)
```

Install `grf` separately if you want the causal-forest wrappers:

```r
install.packages("grf")
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

fit <- fit_cate(
  task,
  method = "dr",
  base_learner = Lrnr_mean$new(),
  cross_validate = FALSE
)

pred <- predict(fit, data)
head(pred)
summary(fit)
```

## Choose Your Workflow

### High-level wrappers

Use `hte_task()`, `fit_cate()`, and `fit_crr()` when you want the standard
package interface for difference-scale or risk-ratio targets.

### GRF / causal-forest path

Use `grf_cate()` or `grf_crr()` when you want the optional `grf` engine and a
wrapper that accepts raw data plus optional nuisance estimates.

### Advanced / low-level sl3

Use `make_hte3_Task_tx()`, `cross_validate_cate()`, `cross_validate_crr()`,
and the `Lrnr_*` classes when you want to assemble learner portfolios
directly.

## Where to Find Documentation

- [Quickstart](articles/quickstart.html): one end-to-end example with the
  high-level wrappers.
- [CATE Workflow](articles/cate.html): choosing among `dr`, `r`, `t`, and
  `ep` for conditional mean differences.
- [CRR Workflow](articles/crr.html): fitting conditional risk ratios with
  `fit_crr()`.
- [Causal Forests with GRF](articles/grf.html): `grf_cate()`, `grf_crr()`,
  direct nuisance inputs, and GRF-specific arguments.
- [Advanced sl3 Integration](articles/advanced-sl3.html): low-level task
  construction, learner control, and manual cross-validation.
- [Reference](reference/index.html): exported functions, helper constructors,
  and learner classes.

## Important Constraints

- `fit_crr()` and `grf_crr()` require a non-negative outcome.
- `grf_cate()` and `grf_crr()` currently support binary treatment only.
- Continuous-treatment CATE currently goes through the R-learner path via
  `fit_cate(..., method = "r")`.
