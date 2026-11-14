# Treatment effects

Use
[`hte_task()`](https://larsvanderlaan.github.io/hte3/reference/hte_task.md)
and
[`fit_cate()`](https://larsvanderlaan.github.io/hte3/reference/fit_cate.md)
to estimate conditional average treatment effects. Choose the covariates
to adjust for with `confounders` and the covariates to condition the
effect on with `modifiers`.

``` r
library(hte3)
library(sl3)

data <- hte3_example_data(n = 150, seed = 2)

task <- hte_task(
  data = data,
  modifiers = c("W1", "W2"),
  confounders = c("W1", "W2", "W3"),
  treatment = "A",
  outcome = "Y",
  propensity_learner = Lrnr_mean$new(),
  outcome_learner = Lrnr_mean$new(),
  mean_learner = Lrnr_mean$new(),
  cross_fit = FALSE
)

dr_model <- fit_cate(
  task,
  method = "dr",
  base_learner = Lrnr_mean$new(),
  cross_validate = FALSE
)

head(predict(dr_model, data))
#> [1] 1.305974 1.305974 1.305974 1.305974 1.305974 1.305974
summary(dr_model)
#> <summary.hte3_model>
#>   Target: CATE
#>   Engine: sl3
#>   Method: dr
#>   Cross-validated: no
#>   Rows: 150
#>   Modifiers: W1, W2
#>   Treatment variable: A
```

## Choose a method

[`fit_cate()`](https://larsvanderlaan.github.io/hte3/reference/fit_cate.md)
supports `method = c("dr", "r", "t", "ep")`.

- `method = "dr"` is the default for binary or categorical treatments.
  It uses estimated treatment probabilities and outcome regressions to
  construct a doubly robust pseudo-outcome, then regresses it on the
  modifiers.
- `method = "r"` regresses outcome residuals on treatment residuals and
  also supports continuous treatments. When the modifiers are a subset
  of the adjustment covariates, it targets an overlap-weighted
  projection rather than the unweighted `E[Y(1) - Y(0) | V]`.
- `method = "t"` fits outcome models by treatment arm, takes their
  difference, and regresses that difference on the modifiers.
- `method = "ep"` updates the initial outcome estimates using a sieve
  basis before fitting the effect regression. See
  [EP-learner](https://larsvanderlaan.github.io/hte3/articles/ep-learner.md)
  for basis-size selection and the differences between EP-DR and EP-R.

## Conditioning on fewer covariates

When `modifiers = V` and `confounders = W` with `V` a strict subset of
`W`, the target of interest is `E[Y(1) - Y(0) | V] = E[tau(W) | V]`.

In the supported binary/categorical-treatment setting:

- DR, EP-DR, and the default two-stage T-learner estimate this target
- the current R-learner instead targets an overlap-weighted projection
  onto functions of `V`

EP-R generally also targets an overlap-weighted projection in this
setting. Both R and EP-R issue a warning when the modifiers are a strict
subset of the adjustment covariates.

## Compare methods

Set `cross_validate = TRUE` to compare methods and use `cv_control` to
set the number of validation folds.

``` r
cv_model <- fit_cate(
  task,
  method = c("dr", "r", "ep"),
  base_learner = Lrnr_mean$new(),
  cross_validate = TRUE,
  cv_control = list(V = 2)
)

cv_model
#> <hte3_model>
#>   Target: CATE
#>   Engine: sl3
#>   Method: portfolio(dr, r, ep)
#>   Cross-validated: yes
#>   Selected candidate: r
#>   Selected method: r
#>   EP basis grid: 2, 4, 8, 12, 16
#>   Modifiers: W1, W2
head(predict(cv_model, data))
#> [1] 1.305974 1.305974 1.305974 1.305974 1.305974 1.305974
```

`summary(cv_model)` reports the selected learner. If the selected
learner is an EP candidate, the summary also reports its targeting
style, basis size, and the grid of basis sizes compared.

## Treatment types and tuning

- Continuous-treatment CATE currently supports `method = "r"` only,
  using the partially linear effect model `A * tau(X)`.
- Use `sieve_num_basis` when you want one fixed EP basis size and
  `sieve_basis_grid` to compare EP basis sizes by cross-validation.
- For direct control over candidate learner objects, see [Using sl3
  directly](https://larsvanderlaan.github.io/hte3/articles/advanced-sl3.md).
