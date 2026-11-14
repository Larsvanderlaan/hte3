# hte3

`hte3` estimates conditional average treatment effects (CATE) and conditional
risk ratios (CRR) in R. Fit DR-, R-, T-, and efficient plug-in (EP) learners,
compare models by cross-validation, and predict effects for new observations.
Use learners from `sl3` or fit forests with the optional `grf` package.

## Installation

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("Larsvanderlaan/hte3", dependencies = TRUE)
```

To use `grf_cate()` or `grf_crr()`, also install `grf`:

```r
install.packages("grf")
```

## Quickstart

Create a task with your data and covariates, fit a model, and predict treatment
effects. This example uses simulated data and simple mean learners.

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

## Guides

| Guide | What it covers |
| --- | --- |
| [Quickstart](articles/quickstart.html) | Fit and predict with `hte_task()` and `fit_cate()`. |
| [Treatment effects](articles/cate.html) | Choose among DR-, R-, T-, and EP-learners for conditional mean differences. |
| [Risk ratios](articles/crr.html) | Estimate conditional risk ratios with `fit_crr()`. |
| [Causal forests](articles/grf.html) | Fit models with `grf_cate()` and `grf_crr()`, with optional nuisance estimates supplied by you. |
| [EP-learner](articles/ep-learner.html) | Choose a sieve basis size and compare EP-DR with EP-R. |
| [Using sl3 directly](articles/advanced-sl3.html) | Construct tasks, train learners, and compare custom sets of models. |

See the [function reference](reference/index.html) for arguments and return values.

## Supported outcomes and treatments

- `fit_crr()` and `grf_crr()` require a non-negative outcome.
- `grf_cate()` and `grf_crr()` currently support binary treatment only.
- For continuous treatments, use `fit_cate(..., method = "r")`, which assumes
  a partially linear effect model `A * tau(X)`.
