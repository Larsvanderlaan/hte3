# Causal Forests with GRF

Read this article if you want the optional `grf` wrappers rather than
the standard
[`hte_task()`](https://larsvanderlaan.github.io/hte3/reference/hte_task.md)
plus `fit_*()` workflow. The GRF path keeps the same `hte3_model` output
interface, but it uses `grf` as the forest engine and accepts raw data
plus optional nuisance estimates directly.

``` r
library(hte3)
```

## When to Use This Workflow

Use the GRF wrappers when all of the following are true:

- you want a forest-based workflow built around the optional `grf`
  package
- you want to call a wrapper directly on raw data instead of building
  [`hte_task()`](https://larsvanderlaan.github.io/hte3/reference/hte_task.md)
  first
- you still want an `hte3_model` that works with
  [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html)

## Choose the Wrapper from the Estimand

- Use
  [`grf_cate()`](https://larsvanderlaan.github.io/hte3/reference/grf_cate.md)
  when the target is a conditional mean difference.
- Use
  [`grf_crr()`](https://larsvanderlaan.github.io/hte3/reference/grf_crr.md)
  when the target is a conditional risk ratio and the outcome is
  non-negative.

## Choose the CATE Method

[`grf_cate()`](https://larsvanderlaan.github.io/hte3/reference/grf_cate.md)
supports `method = c("dr", "r", "ep")`.

- `method = "r"` is the direct causal-forest path and trains
  [`grf::causal_forest()`](https://rdrr.io/pkg/grf/man/causal_forest.html).
- `method = "dr"` uses the DR-learner with GRF nuisance and base
  learners.
- `method = "ep"` uses the EP-learner with GRF nuisance and base
  learners.

Within the EP family, `ep_targeting_style = "dr"` selects EP-DR and
`ep_targeting_style = "r"` selects EP-R. If you pass
`ep_targeting_style = c("dr", "r")`, wrapper-level cross-validation can
compare both EP variants.

If you omit both `method` and `cross_validate`,
[`grf_cate()`](https://larsvanderlaan.github.io/hte3/reference/grf_cate.md)
compares the default EP-R candidate and the direct GRF R-learner.

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

## Direct Nuisance Inputs

The GRF wrappers accept direct nuisance estimates:

- `mu.hat`: an `n x 2` matrix ordered as `(control, treatment)`
- `pi.hat`: a length-`n` vector for `P(A = treatment_level | X)`
- `m.hat`: an optional length-`n` vector of marginal outcome means

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

## Important Arguments

The GRF-specific arguments worth understanding are:

- `cross_fit` and `folds` for nuisance estimation
- `cross_validate` and `cv_control` for outer learner selection
- `ep_targeting_style` and `ep_r_targeting_basis` for the GRF-backed EP
  variants
- `treatment_level` and `control_level` when treatment labels are not
  the default `1` and `0`
- `grf_params` for GRF hyperparameters such as `num.trees` or
  `min.node.size`
- `tune` for the amount of GRF tuning: `"light"`, `"none"`, or `"all"`

In the GRF learners, `tune = "light"` tunes `sample.fraction`, `mtry`,
and `min.node.size` with a modest search budget. `tune = "none"`
disables tuning, and `tune = "all"` requests GRF’s full tuning path.

## Important Constraints

- [`grf_cate()`](https://larsvanderlaan.github.io/hte3/reference/grf_cate.md)
  and
  [`grf_crr()`](https://larsvanderlaan.github.io/hte3/reference/grf_crr.md)
  currently support binary treatment only.
- [`grf_crr()`](https://larsvanderlaan.github.io/hte3/reference/grf_crr.md)
  requires a non-negative outcome.
- If you want the standard
  [`hte_task()`](https://larsvanderlaan.github.io/hte3/reference/hte_task.md)
  plus `fit_*()` interface or continuous-treatment CATE, use the non-GRF
  workflow instead.

## Helper Constructors

When the wrappers are too narrow, these helpers expose the forest-backed
building blocks directly:

- `make_grf_base_learner()`
- [`make_grf_nuisance_learners()`](https://larsvanderlaan.github.io/hte3/reference/make_grf_nuisance_learners.md)
- [`make_grf_cate_learners()`](https://larsvanderlaan.github.io/hte3/reference/make_grf_cate_learners.md)
- [`make_grf_crr_learners()`](https://larsvanderlaan.github.io/hte3/reference/make_grf_crr_learners.md)

Use these when you want manual cross-validation or a custom `sl3`
learner library rather than the default wrapper behavior.
