# Build GRF-Backed CATE Learners

Build GRF-Backed CATE Learners

## Usage

``` r
make_grf_cate_learners(
  method = c("ep", "r"),
  hte3_task = NULL,
  treatment_level = 1,
  control_level = 0,
  tune = c("light", "none", "all"),
  grf_params = list(),
  sieve_num_basis = NULL,
  sieve_basis_grid = NULL,
  sieve_interaction_order = 3,
  screen_basis_with_lasso = FALSE,
  ep_targeting_style = "r",
  ep_r_targeting_basis = "v_plus_propensity"
)
```

## Arguments

- method:

  Character vector of GRF-supported CATE methods.

- hte3_task:

  Optional task used when expanding EP basis-size portfolios.

- treatment_level:

  Treated level for the target contrast.

- control_level:

  Control level for the target contrast.

- tune:

  Tuning mode for the final GRF learner(s).

- grf_params:

  Optional named list of GRF arguments passed to the learners.

- sieve_num_basis:

  Optional EP basis size for single-fit EP learners.

- sieve_basis_grid:

  Optional EP basis grid.

- sieve_interaction_order:

  EP interaction order.

- screen_basis_with_lasso:

  Whether EP basis screening is enabled.

- ep_targeting_style:

  EP targeting variant used when `method` includes `"ep"`.

- ep_r_targeting_basis:

  First-stage basis construction used for EP-R.

## Value

A list of learner objects.
