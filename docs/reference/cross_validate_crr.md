# Cross-Validate CRR Models

Cross-validates a collection of CRR learners using the CRR selector
loss.

## Usage

``` r
cross_validate_crr(
  hte_learners,
  hte3_task,
  cv_control = NULL,
  treatment_level = NULL,
  control_level = NULL,
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

- cv_control:

  A list of control parameters for cross-validation passed to
  [`Lrnr_sl`](https://tlverse.org/sl3/reference/Lrnr_sl.html). Default
  is `NULL`.

- treatment_level:

  The treatment level to be considered the treated group in the contrast
  used for selection.

- control_level:

  The treatment level to be considered the control group in the contrast
  used for selection.

- ...:

  Additional arguments to pass to the loss function and other functions.

## Value

A list containing the following elements:

- `lrnr_sl`: A `Lrnr_sl` object that represents the cross-validated
  ensemble of CRR learners.

- `cv_risk`: The cross-validation risk associated with the ensemble,
  which serves as a measure of the ensemble's predictive performance.

- `coefficients`: The coefficients derived from the cross-validation
  process, providing insights into the relative importance of different
  learners within the ensemble.

- `selection_summary`: A compact summary of candidate labels, selected
  candidate, selected method, and any EP basis-size metadata.

## Examples

``` r
if (FALSE) { # \dontrun{
library(sl3)

data <- hte3_example_data(n = 80, seed = 1)
task <- hte_task(
  data = data,
  modifiers = c("W1", "W2"),
  confounders = c("W1", "W2", "W3"),
  treatment = "A",
  outcome = "Y_binary",
  propensity_learner = Lrnr_mean$new(),
  outcome_learner = Lrnr_mean$new(),
  mean_learner = Lrnr_mean$new(),
  cross_fit = FALSE
)

cv_fit <- cross_validate_crr(
  list(
    Lrnr_crr_IPW$new(base_learner = Lrnr_mean$new()),
    Lrnr_crr_EP$new(base_learner = Lrnr_mean$new(), sieve_num_basis = 4)
  ),
  task,
  cv_control = list(V = 2)
)

cv_fit$selection_summary
} # }
```
