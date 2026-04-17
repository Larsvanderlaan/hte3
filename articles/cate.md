# CATE Workflow

Read this article if you are fitting conditional average treatment
effects with
[`hte_task()`](https://larsvanderlaan.github.io/hte3/reference/hte_task.md)
and
[`fit_cate()`](https://larsvanderlaan.github.io/hte3/reference/fit_cate.md)
and need to choose among the supported CATE learner families.

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

## Choose the Wrapper Method

[`fit_cate()`](https://larsvanderlaan.github.io/hte3/reference/fit_cate.md)
supports `method = c("dr", "r", "t", "ep")`.

- `method = "dr"` is the main difference-scale default when treatment is
  binary or categorical and you want a target-aligned CATE fit that uses
  both propensity and outcome nuisance estimates.
- `method = "r"` is the residual-on-residual path and the current
  continuous-treatment option. In reduced-modifier settings, it targets
  an overlap-weighted projection rather than the unweighted
  `E[Y(1) - Y(0) | V]`.
- `method = "t"` is the baseline that fits outcome models by treatment
  arm and subtracts them.
- `method = "ep"` is the sieve-based plug-in path. Use it when you want
  the EP family; see [EP-Learner Sieve
  Tuning](https://larsvanderlaan.github.io/hte3/articles/ep-learner.md)
  for basis-size tuning and EP-DR versus EP-R details.

## Reduced-Modifier Targets

When `modifiers = V` and `confounders = W` with `V` a strict subset of
`W`, the target of interest is `E[Y(1) - Y(0) | V] = E[tau(W) | V]`.

In the supported binary/categorical-treatment setting:

- DR, EP, and the default two-stage T-learner align with that target
- the current R-learner instead targets an overlap-weighted projection
  onto functions of `V`

This is why `fit_cate(..., method = "r")` warns in reduced-modifier
settings.

## Compare CATE Learner Families

Wrapper-level CATE cross-validation is controlled by `cross_validate`
and `cv_control`.

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
learner is an EP candidate, the summary also reports the selected basis
size and any EP selection metadata.

## Important Constraints

- Continuous-treatment CATE currently supports `method = "r"` only,
  using the partially linear effect model `A * tau(X)`.
- Use `sieve_num_basis` when you want one fixed EP basis size and
  `sieve_basis_grid` when you want wrapper-level EP basis-size
  selection.
- For direct control over candidate learner objects, see [Advanced sl3
  Integration](https://larsvanderlaan.github.io/hte3/articles/advanced-sl3.md).
