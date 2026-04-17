# Quickstart

Read this article if you want one minimal end-to-end example with the
high-level wrappers before choosing a more specific workflow.

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

## Cross-Fit vs Cross-Validate

Two arguments control different stages of the workflow:

- `cross_fit` in
  [`hte_task()`](https://larsvanderlaan.github.io/hte3/reference/hte_task.md)
  controls nuisance estimation.
- `cross_validate` in
  [`fit_cate()`](https://larsvanderlaan.github.io/hte3/reference/fit_cate.md)
  or
  [`fit_crr()`](https://larsvanderlaan.github.io/hte3/reference/fit_crr.md)
  controls selection across candidate HTE learners.

These are separate decisions. A common analysis cross-fits nuisances and
only turns on wrapper-level cross-validation when comparing learner
families.

## What to Read Next

- For CATE method choice, see [CATE
  Workflow](https://larsvanderlaan.github.io/hte3/articles/cate.md).
- For relative effects with non-negative outcomes, see [CRR
  Workflow](https://larsvanderlaan.github.io/hte3/articles/crr.md).
- For
  [`grf_cate()`](https://larsvanderlaan.github.io/hte3/reference/grf_cate.md)
  and
  [`grf_crr()`](https://larsvanderlaan.github.io/hte3/reference/grf_crr.md),
  see [Causal Forests with
  GRF](https://larsvanderlaan.github.io/hte3/articles/grf.md).
- For low-level `sl3` control, see [Advanced sl3
  Integration](https://larsvanderlaan.github.io/hte3/articles/advanced-sl3.md).
