# Quickstart

Create a task, fit a DR-learner, and predict conditional average
treatment effects. The example uses simulated data and simple mean
learners, with cross-fitting and model selection disabled.

``` r
library(hte3)
library(sl3)

data <- hte3_example_data(n = 120, seed = 1)

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

head(predict(model, data))
#> [1] 1.027376 1.027376 1.027376 1.027376 1.027376 1.027376
summary(model)
#> <summary.hte3_model>
#>   Target: CATE
#>   Engine: sl3
#>   Method: dr
#>   Cross-validated: no
#>   Rows: 120
#>   Modifiers: W1, W2, W3
#>   Treatment variable: A
```

## Cross-fitting and model selection

Use these arguments to control nuisance estimation and model selection:

- `cross_fit` in
  [`hte_task()`](https://larsvanderlaan.github.io/hte3/reference/hte_task.md)
  controls nuisance estimation.
- `cross_validate` in
  [`fit_cate()`](https://larsvanderlaan.github.io/hte3/reference/fit_cate.md)
  or
  [`fit_crr()`](https://larsvanderlaan.github.io/hte3/reference/fit_crr.md)
  controls selection across candidate HTE learners.

Set `cross_fit = TRUE` to estimate each observation’s nuisance values
using other folds. Set `cross_validate = TRUE` to compare candidate HTE
learners.

## Next steps

- [Treatment
  effects](https://larsvanderlaan.github.io/hte3/articles/cate.md):
  choose among DR-, R-, T-, and EP-learners.
- [Risk ratios](https://larsvanderlaan.github.io/hte3/articles/crr.md):
  estimate relative effects for non-negative outcomes.
- [Causal
  forests](https://larsvanderlaan.github.io/hte3/articles/grf.md): use
  [`grf_cate()`](https://larsvanderlaan.github.io/hte3/reference/grf_cate.md)
  and
  [`grf_crr()`](https://larsvanderlaan.github.io/hte3/reference/grf_crr.md).
- [Using sl3
  directly](https://larsvanderlaan.github.io/hte3/articles/advanced-sl3.md):
  construct and train custom learners.
