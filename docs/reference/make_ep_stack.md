# Create an Ensemble of CATE EP-learners with Varying Sieve Basis Sizes

This function creates an ensemble of CATE EP-learners of varying sieve
dimensions for use with cross-validation.

## Usage

``` r
make_ep_stack(
  base_learner,
  hte3_task,
  treatment_level = 1,
  control_level = 0,
  sieve_basis_grid = NULL,
  targeting_style = "dr",
  r_targeting_basis = "v_plus_propensity"
)
```

## Arguments

- base_learner:

  The base learner of `Lrnr_cate_EP`.

- hte3_task:

  A `hte3_Task` object containing the data and necessary information for
  the heterogeneous treatment effect estimation.

- treatment_level:

  Treatment level used for the treated arm in each EP learner in the
  stack.

- control_level:

  Reference treatment level used for the control arm in each EP learner
  in the stack.

- sieve_basis_grid:

  Optional vector of sieve basis sizes to include in the ensemble. If
  `NULL`, the default heuristic grid uses the first-stage EP sieve
  dimension for the chosen targeting style.

- targeting_style:

  EP targeting style or styles to include in the stack. Use
  `c("dr", "r")` to build both EP variants across the requested basis
  grid.

- r_targeting_basis:

  First-stage basis construction used for EP-R fits in the stack.

## Value

A stack of CATE EP meta-learners with varying sieve basis sizes.

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
  outcome = "Y",
  propensity_learner = Lrnr_mean$new(),
  outcome_learner = Lrnr_mean$new(),
  mean_learner = Lrnr_mean$new(),
  cross_fit = FALSE
)

ep_stack <- make_ep_stack(
  base_learner = Lrnr_mean$new(),
  hte3_task = task,
  sieve_basis_grid = c(2, 4, 6),
  targeting_style = c("dr", "r")
)

ep_stack
} # }
```
