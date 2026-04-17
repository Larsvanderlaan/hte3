# Helper function for estimation of propensity score using `sl3`.

Helper function for estimation of propensity score using `sl3`.

## Usage

``` r
estimate_pi(
  W,
  A,
  binomial_learner = Lrnr_gam$new(family = "binomial"),
  weights = NULL,
  treatment_level = max(A),
  cross_fit_and_cv = TRUE,
  return_learner = FALSE,
  folds = 10,
  ...
)
```

## Arguments

- W:

  Covariate matrix or data frame used to estimate the propensity score.

- A:

  Observed treatment assignments.

- binomial_learner:

  Base `sl3` learner used for the binomial treatment regression.

- weights:

  Optional observation weights.

- treatment_level:

  Treatment level defining the positive class for the returned
  propensity estimates.

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
