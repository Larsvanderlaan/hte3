# Causal forests

Use
[`grf_cate()`](https://larsvanderlaan.github.io/hte3/reference/grf_cate.md)
and
[`grf_crr()`](https://larsvanderlaan.github.io/hte3/reference/grf_crr.md)
to fit treatment-effect models with the `grf` package. Pass your data
directly, with optional nuisance estimates. Both functions return an
`hte3_model` with [`predict()`](https://rdrr.io/r/stats/predict.html)
and [`summary()`](https://rdrr.io/r/base/summary.html) methods.

``` r

library(hte3)
```

## Choose the target

- Use
  [`grf_cate()`](https://larsvanderlaan.github.io/hte3/reference/grf_cate.md)
  when the target is a conditional mean difference.
- Use
  [`grf_crr()`](https://larsvanderlaan.github.io/hte3/reference/grf_crr.md)
  when the target is a conditional risk ratio and the outcome is
  non-negative.

## Choose a CATE method

[`grf_cate()`](https://larsvanderlaan.github.io/hte3/reference/grf_cate.md)
supports `method = c("dr", "r", "ep")`.

- `method = "r"` trains
  [`grf::causal_forest()`](https://rdrr.io/pkg/grf/man/causal_forest.html).
- `method = "dr"` uses the DR-learner with GRF nuisance and base
  learners.
- `method = "ep"` uses the EP-learner with GRF nuisance and base
  learners.

Within the EP family, `ep_targeting_style = "dr"` selects EP-DR and
`ep_targeting_style = "r"` selects EP-R. If you pass
`ep_targeting_style = c("dr", "r")`, cross-validation can compare both
EP variants.

If you omit both `method` and `cross_validate`,
[`grf_cate()`](https://larsvanderlaan.github.io/hte3/reference/grf_cate.md)
compares the default EP-R candidate and the GRF causal forest.

``` r

data <- hte3_example_data(n = 200, seed = 1)

fit <- grf_cate(
  data = data,
  modifiers = c("W1", "W2"),
  confounders = c("W1", "W2", "W3"),
  treatment = "A",
  outcome = "Y",
  method = "r",
  cross_fit = FALSE
)

head(predict(fit, data))
summary(fit)
```

## `grf_crr()`

[`grf_crr()`](https://larsvanderlaan.github.io/hte3/reference/grf_crr.md)
supports `method = c("ep", "ipw", "t")`. If you omit `method`, the
wrapper fits one EP model.

``` r

fit <- grf_crr(
  data = data,
  modifiers = c("W1", "W2"),
  confounders = c("W1", "W2", "W3"),
  treatment = "A",
  outcome = "Y_binary",
  method = "ep",
  cross_fit = FALSE
)

head(predict(fit, data))
summary(fit)
```

## Supply nuisance estimates

The GRF wrappers accept direct nuisance estimates:

- `mu.hat`: an `n x 2` matrix ordered as `(control, treatment)`
- `pi.hat`: a length-`n` vector for `P(A = treatment_level | X)`
- `m.hat`: an optional length-`n` vector estimating `E[Y | X]`

If `m.hat` is omitted and both `mu.hat` and `pi.hat` are supplied,
`hte3` computes `m.hat` internally.

``` r

mu.hat <- cbind(data$mu0, data$mu1)
pi.hat <- data$pi1

fit <- grf_cate(
  data = data,
  modifiers = c("W1", "W2"),
  confounders = c("W1", "W2", "W3"),
  treatment = "A",
  outcome = "Y",
  mu.hat = mu.hat,
  pi.hat = pi.hat,
  cross_fit = FALSE
)
```

## Tuning and cross-validation

Use these arguments to control the fit:

- `cross_fit` and `folds` for nuisance estimation
- `cross_validate` and `cv_control` for outer learner selection
- `ep_targeting_style` and `ep_r_targeting_basis` for the EP variants
- `treatment_level` and `control_level` when treatment labels are not
  the default `1` and `0`
- `grf_params` for GRF hyperparameters such as `num.trees` or
  `min.node.size`
- `tune` for the amount of GRF tuning: `"light"`, `"none"`, or `"all"`

In the GRF learners, `tune = "light"` tunes `sample.fraction`, `mtry`,
and `min.node.size` with a modest search budget. `tune = "none"`
disables tuning, and `tune = "all"` tunes all parameters supported by
GRF.

## Supported treatments and outcomes

- [`grf_cate()`](https://larsvanderlaan.github.io/hte3/reference/grf_cate.md)
  and
  [`grf_crr()`](https://larsvanderlaan.github.io/hte3/reference/grf_crr.md)
  currently support binary treatment only.
- [`grf_crr()`](https://larsvanderlaan.github.io/hte3/reference/grf_crr.md)
  requires a non-negative outcome.
- For continuous-treatment CATE, use
  [`hte_task()`](https://larsvanderlaan.github.io/hte3/reference/hte_task.md)
  and `fit_cate(..., method = "r")`.

## Construct learners directly

Use these helpers to construct GRF learners for a custom `sl3` model
library or manual cross-validation:

- `make_grf_base_learner()`
- [`make_grf_nuisance_learners()`](https://larsvanderlaan.github.io/hte3/reference/make_grf_nuisance_learners.md)
- [`make_grf_cate_learners()`](https://larsvanderlaan.github.io/hte3/reference/make_grf_cate_learners.md)
- [`make_grf_crr_learners()`](https://larsvanderlaan.github.io/hte3/reference/make_grf_crr_learners.md)
