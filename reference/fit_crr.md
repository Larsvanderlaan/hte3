# Fit a CRR Model

High-level wrapper for fitting conditional risk-ratio models with
non-negative outcomes. The wrapper supports both single-method fits and
cross-validated portfolios across multiple learner families.

## Usage

``` r
fit_crr(
  task,
  method = supported_crr_methods(),
  base_learner = get_autoML(),
  treatment_level = NULL,
  control_level = NULL,
  cross_validate = NULL,
  cv_control = NULL,
  sieve_num_basis = NULL,
  sieve_basis_grid = NULL,
  sieve_interaction_order = 3
)
```

## Arguments

- task:

  An `hte3_Task`.

- method:

  Either a single method or a character vector of methods drawn from
  `"ep"`, `"ipw"`, and `"t"`. A vector creates a learner portfolio that
  is selected by cross-validation.

- base_learner:

  Base supervised learner used by the meta-learner.

- treatment_level:

  Optional treated level used for the contrast.

- control_level:

  Optional control level used for the contrast.

- cross_validate:

  Whether to cross-validate across candidate learners. Defaults to
  `TRUE` for stacks and is automatically enabled for learner portfolios.

- cv_control:

  Optional `sl3` cross-validation control list.

- sieve_num_basis:

  Optional sieve basis size for EP-learner fits.

- sieve_basis_grid:

  Optional EP basis-size grid used when `method` includes `"ep"` and
  wrapper-level cross-validation is enabled. If `NULL`, the wrapper uses
  the heuristic default grid `c(d, 2d, 4d, 6d, 8d)`, where `d` is the
  number of modifiers.

- sieve_interaction_order:

  Interaction order for EP-learner sieve construction.

## Value

An object of class `hte3_model`.

## Details

CRR models are supported for binary or categorical treatment tasks and
require a non-negative outcome.

For EP-learner fits, `sieve_num_basis` controls a single fixed sieve
size when `cross_validate = FALSE`. When `cross_validate = TRUE`,
`sieve_basis_grid` controls which EP basis sizes are compared by the
wrapper. After fitting,
[`summary()`](https://rdrr.io/r/base/summary.html) reports the selected
candidate and the selected EP basis size when applicable.

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

fit <- fit_crr(
  task,
  method = "ep",
  base_learner = Lrnr_mean$new(),
  cross_validate = TRUE,
  cv_control = list(V = 2),
  sieve_basis_grid = c(3, 6, 9)
)

summary(fit)
} # }
```
