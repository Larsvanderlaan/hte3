# Helper function for estimation of outcome regression using `sl3`.

Helper function for estimation of outcome regression using `sl3`.

## Usage

``` r
estimate_mu(
  W,
  A,
  Y,
  treatment_level,
  learner = Lrnr_gam$new(family = "gaussian"),
  weights = NULL,
  cross_fit_and_cv = TRUE,
  stratified_by_trt = TRUE,
  return_learner = FALSE,
  folds = 10,
  ...
)
```

## Arguments

- W:

  Covariate matrix or data frame used to estimate the outcome
  regression.

- A:

  Observed treatment assignments.

- Y:

  Observed outcomes.

- treatment_level:

  Treatment level for which the outcome regression should be estimated.

- learner:

  Base `sl3` learner used for the regression fit.

- weights:

  Optional observation weights.

- cross_fit_and_cv:

  Whether to wrap the nuisance learner in cross-fitting and learner
  selection.

- stratified_by_trt:

  Whether to fit the outcome regression separately within the specified
  treatment arm.

- return_learner:

  Whether to include the trained learner in the returned list.

- folds:

  Cross-validation fold specification passed to the constructed
  `sl3_Task`.

- ...:

  Additional arguments forwarded to
  [`sl3_Task`](https://tlverse.org/sl3/reference/sl3_Task.html).
