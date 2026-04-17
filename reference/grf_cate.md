# Fit GRF-Backed CATE Models

Fit GRF-Backed CATE Models

## Usage

``` r
grf_cate(
  data,
  modifiers,
  confounders = modifiers,
  treatment,
  outcome,
  ...,
  method = NULL,
  cross_validate = NULL,
  mu.hat = NULL,
  pi.hat = NULL,
  m.hat = NULL,
  treatment_level = 1,
  control_level = 0,
  cross_fit = TRUE,
  folds = 10,
  cv_control = NULL,
  ep_targeting_style = "r",
  ep_r_targeting_basis = "v_plus_propensity",
  tune = c("light", "none", "all"),
  grf_params = list()
)
```

## Arguments

- data:

  A data frame or data.table.

- modifiers:

  Effect modifiers.

- confounders:

  Adjustment covariates.

- treatment:

  Treatment column.

- outcome:

  Outcome column.

- ...:

  Additional arguments passed to
  [`hte_task()`](https://larsvanderlaan.github.io/hte3/reference/hte_task.md).

- method:

  Optional GRF CATE method specification.

- cross_validate:

  Optional logical flag for outer learner selection.

- mu.hat:

  Optional `n x 2` matrix of nuisance outcome regressions ordered as
  `(control, treatment)`.

- pi.hat:

  Optional length-`n` vector of treatment propensities
  `P(A = treatment_level | X)`.

- m.hat:

  Optional length-`n` vector of marginal outcome means.

- treatment_level:

  Treated level for the target contrast.

- control_level:

  Control level for the target contrast.

- cross_fit:

  Whether nuisance learners should be cross-fitted.

- folds:

  Number of folds used for nuisance cross-fitting.

- cv_control:

  Optional outer cross-validation control list.

- ep_targeting_style:

  EP targeting variant used when `method` includes `"ep"`.

- ep_r_targeting_basis:

  First-stage basis construction used for EP-R.

- tune:

  Tuning mode for the final GRF learner(s).

- grf_params:

  Optional named list of GRF arguments passed to the learners.

## Value

An `hte3_model`.
