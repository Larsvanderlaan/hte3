# hte3

`hte3` provides R workflows for estimating conditional average treatment
effects (CATE) and conditional risk ratios (CRR). Use the high-level
wrappers for standard analyses, the GRF wrappers for forest-based fits,
and the lower-level `sl3` interface when you need direct control over
nuisance learners or candidate HTE models.

## Installation

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("Larsvanderlaan/hte3", dependencies = TRUE)
```

Install `grf` separately if you want the causal-forest wrappers:

``` r
install.packages("grf")
```

## Quickstart

``` r
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

Use
[`hte_task()`](https://larsvanderlaan.github.io/hte3/reference/hte_task.md),
[`fit_cate()`](https://larsvanderlaan.github.io/hte3/reference/fit_cate.md),
and
[`fit_crr()`](https://larsvanderlaan.github.io/hte3/reference/fit_crr.md)
when you want the standard package interface for difference-scale or
risk-ratio targets.

### GRF / causal-forest path

Use
[`grf_cate()`](https://larsvanderlaan.github.io/hte3/reference/grf_cate.md)
or
[`grf_crr()`](https://larsvanderlaan.github.io/hte3/reference/grf_crr.md)
when you want the optional `grf` engine and a wrapper that accepts raw
data plus optional nuisance estimates.

### Advanced / low-level sl3

Use
[`make_hte3_Task_tx()`](https://larsvanderlaan.github.io/hte3/reference/make_hte3_Task_tx.md),
[`cross_validate_cate()`](https://larsvanderlaan.github.io/hte3/reference/cross_validate_cate.md),
[`cross_validate_crr()`](https://larsvanderlaan.github.io/hte3/reference/cross_validate_crr.md),
and the `Lrnr_*` classes when you want to assemble learner portfolios
directly.

## Where to Find Documentation

- [Quickstart](https://larsvanderlaan.github.io/hte3/articles/quickstart.md):
  one end-to-end example with the high-level wrappers.
- [CATE
  Workflow](https://larsvanderlaan.github.io/hte3/articles/cate.md):
  choosing among `dr`, `r`, `t`, and `ep` for conditional mean
  differences.
- [CRR Workflow](https://larsvanderlaan.github.io/hte3/articles/crr.md):
  fitting conditional risk ratios with
  [`fit_crr()`](https://larsvanderlaan.github.io/hte3/reference/fit_crr.md).
- [Causal Forests with
  GRF](https://larsvanderlaan.github.io/hte3/articles/grf.md):
  [`grf_cate()`](https://larsvanderlaan.github.io/hte3/reference/grf_cate.md),
  [`grf_crr()`](https://larsvanderlaan.github.io/hte3/reference/grf_crr.md),
  direct nuisance inputs, and GRF-specific arguments.
- [Advanced sl3
  Integration](https://larsvanderlaan.github.io/hte3/articles/advanced-sl3.md):
  low-level task construction, learner control, and manual
  cross-validation.
- [Reference](https://larsvanderlaan.github.io/hte3/reference/index.md):
  exported functions, helper constructors, and learner classes.

## Important Constraints

- [`fit_crr()`](https://larsvanderlaan.github.io/hte3/reference/fit_crr.md)
  and
  [`grf_crr()`](https://larsvanderlaan.github.io/hte3/reference/grf_crr.md)
  require a non-negative outcome.
- [`grf_cate()`](https://larsvanderlaan.github.io/hte3/reference/grf_cate.md)
  and
  [`grf_crr()`](https://larsvanderlaan.github.io/hte3/reference/grf_crr.md)
  currently support binary treatment only.
- Continuous-treatment CATE currently goes through the R-learner path
  via `fit_cate(..., method = "r")`.
