# Risk ratios

Use
[`hte_task()`](https://larsvanderlaan.github.io/hte3/reference/hte_task.md)
and
[`fit_crr()`](https://larsvanderlaan.github.io/hte3/reference/fit_crr.md)
to estimate conditional risk ratios for non-negative outcomes. Choose
among EP-, IPW-, and T-learners, or compare them by cross-validation.

``` r
library(hte3)
library(sl3)

data <- hte3_example_data(n = 150, seed = 3)

task <- hte_task(
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

crr_model <- fit_crr(
  task,
  method = "ipw",
  base_learner = Lrnr_mean$new(),
  cross_validate = FALSE
)

head(predict(crr_model, data))
#> [1] 0.1677935 0.1677935 0.1677935 0.1677935 0.1677935 0.1677935
summary(crr_model)
#> <summary.hte3_model>
#>   Target: CRR
#>   Engine: sl3
#>   Method: ipw
#>   Cross-validated: no
#>   Rows: 150
#>   Modifiers: W1, W2, W3
#>   Treatment variable: A
```

## Choose a method

[`fit_crr()`](https://larsvanderlaan.github.io/hte3/reference/fit_crr.md)
supports `method = c("ep", "ipw", "t")`.

- `method = "ep"` is the default. It updates initial outcome regressions
  using a sieve basis before fitting the risk-ratio model.
- `method = "ipw"` uses inverse treatment probabilities to weight the
  fit.
- `method = "t"` estimates outcome means in each treatment arm and fits
  a regression of their log ratio on the modifiers.

As with CATE, create a task, choose a method, and compare candidates
with `cross_validate = TRUE`.

## Compare methods

``` r
cv_model <- fit_crr(
  task,
  method = c("ipw", "t", "ep"),
  base_learner = Lrnr_mean$new(),
  cross_validate = TRUE,
  cv_control = list(V = 2)
)

cv_model
#> <hte3_model>
#>   Target: CRR
#>   Engine: sl3
#>   Method: portfolio(ipw, t, ep)
#>   Cross-validated: yes
#>   Selected candidate: ep[basis=24, interaction=3]
#>   Selected method: ep
#>   Selected EP basis size: 24
#>   EP basis grid: 3, 6, 12, 18, 24
#>   Modifiers: W1, W2, W3
head(predict(cv_model, data))
#> [1] 0.151008 0.151008 0.151008 0.151008 0.151008 0.151008
```

If `method = "ep"` and `cross_validate = TRUE`, the wrapper also
compares EP fits across basis sizes. Use `sieve_basis_grid` when you
want to control that grid explicitly, and use `sieve_num_basis` for one
fixed EP fit.

## Outcomes and cross-fitting

- [`fit_crr()`](https://larsvanderlaan.github.io/hte3/reference/fit_crr.md)
  requires a non-negative outcome.
- `cross_fit` in
  [`hte_task()`](https://larsvanderlaan.github.io/hte3/reference/hte_task.md)
  controls nuisance estimation and is separate from model selection with
  `cross_validate`.
- To construct and compare your own learners, see [Using sl3
  directly](https://larsvanderlaan.github.io/hte3/articles/advanced-sl3.md).
