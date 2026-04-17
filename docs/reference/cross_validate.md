# Cross-Validate Heterogeneous Treatment Effect Models

Cross-validates a collection of heterogeneous treatment effect learners.

## Usage

``` r
cross_validate(
  hte_learners,
  hte3_task,
  cv_metalearner = Lrnr_cv_selector$new(loss_squared_error),
  cv_control = NULL,
  ...
)
```

## Arguments

- hte_learners:

  A single `Lrnr_hte` learner, a list of `Lrnr_hte` learners, or a
  `Stack` of `Lrnr_hte` learners to cross-validate.

- hte3_task:

  An `hte3_Task` object containing the data and necessary information
  for heterogeneous treatment effect estimation.

- cv_metalearner:

  An optional metalearner (`Lrnr_base` object) used to combine the
  cross-validated learners. Default is
  `Lrnr_cv_selector$new(loss_squared_error)`.

- cv_control:

  A list of control parameters for cross-validation passed to
  [`Lrnr_sl`](https://tlverse.org/sl3/reference/Lrnr_sl.html). Default
  is `NULL`.

- ...:

  Additional arguments to pass to the loss function and other functions.

## Value

A trained `Lrnr_sl` object containing the cross-validated ensemble of
heterogeneous treatment effect learners.
