# Helper function for estimation of treatment-averaged outcome regression using `sl3`.

Helper function for estimation of treatment-averaged outcome regression
using `sl3`.

## Usage

``` r
estimate_m(
  W,
  Y,
  learner = Lrnr_gam$new(family = "gaussian"),
  weights = NULL,
  cross_fit_and_cv = TRUE,
  return_learner = FALSE,
  folds = 10,
  ...
)
```

## Arguments

- W:

  Covariate matrix or data frame used to estimate the conditional mean
  outcome.

- Y:

  Observed outcomes.

- learner:

  Base `sl3` learner used for the regression fit.

- weights:

  Optional observation weights.

- cross_fit_and_cv:

  Whether to wrap the nuisance learner in cross-fitting and learner
  selection.

- return_learner:

  Whether to include the trained learner in the returned list.

- folds:

  Cross-validation fold specification passed to the constructed
  `sl3_Task`.

- ...:

  Additional arguments forwarded to
  [`sl3_Task`](https://tlverse.org/sl3/reference/sl3_Task.html).
