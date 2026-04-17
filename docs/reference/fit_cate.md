# Fit a CATE Model

High-level wrapper for fitting conditional average treatment effect
models. The wrapper supports both single-method fits and cross-validated
portfolios across multiple learner families.

## Usage

``` r
fit_cate(
  task,
  method = supported_cate_methods(),
  base_learner = get_autoML(),
  treatment_level = NULL,
  control_level = NULL,
  cross_validate = NULL,
  cv_control = NULL,
  sieve_num_basis = NULL,
  sieve_basis_grid = NULL,
  sieve_interaction_order = 3,
  screen_basis_with_lasso = FALSE,
  ep_targeting_style = "dr",
  ep_r_targeting_basis = "v_plus_propensity"
)
```

## Arguments

- task:

  An `hte3_Task`.

- method:

  Either a single method or a character vector of methods drawn from
  `"dr"`, `"r"`, `"t"`, and `"ep"`. A vector creates a learner portfolio
  that is selected by cross-validation.

- base_learner:

  Base supervised learner used by the meta-learner.

- treatment_level:

  Optional treated level used for binary or categorical contrasts.

- control_level:

  Optional control level used for binary or categorical contrasts.

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

- screen_basis_with_lasso:

  Whether to screen EP sieve terms with lasso.

- ep_targeting_style:

  EP targeting variant or variants used when `method` includes `"ep"`.
  Use `"dr"` for EP-DR, `"r"` for EP-R, or `c("dr", "r")` to let
  wrapper-level cross-validation compare both EP variants.

- ep_r_targeting_basis:

  First-stage basis construction used for `ep_targeting_style = "r"`.
  The default `"v_plus_propensity"` uses the modifier set plus the
  treated-arm propensity score; `"full_w"` uses the full confounder set.

## Value

An object of class `hte3_model`.

## Details

For continuous-treatment CATE tasks, the supported high-level method is
`"r"` only. That path uses the partially linear R-learner effect model
`A * tau(X)`.

If `modifiers = V` and `confounders = W` with `V` a strict subset of
`W`, the natural target is `E[Y(1) - Y(0) | V]`. In the supported
binary/categorical-treatment setting, DR-, EP-, and the default
two-stage T-learner target that modifier-conditional CATE surface. The
current R-learner does not generally target that object in the
reduced-modifier setting and warns at fit time.

For EP-learner fits, `sieve_num_basis` controls a single fixed sieve
size when `cross_validate = FALSE`. When `cross_validate = TRUE`,
`sieve_basis_grid` controls which EP basis sizes are compared by the
wrapper. After fitting,
[`summary()`](https://rdrr.io/r/base/summary.html) reports the selected
candidate, the selected EP targeting style, and the selected EP basis
size when applicable. Wrapper-level EP tuning continues to use the DR
selector loss from
[`cross_validate_cate()`](https://larsvanderlaan.github.io/hte3/reference/cross_validate_cate.md),
including for the EP-R variant.

The CATE EP family now has two variants. EP-DR
(`ep_targeting_style = "dr"`) is the standard plug-in DR-style EP
variant. EP-R (`ep_targeting_style = "r"`) uses an overlap-weighted
targeting objective, can be more stable under severe overlap violations,
and currently supports binary-treatment CATE tasks only. For EP-R, the
default first-stage basis is
`ep_r_targeting_basis = "v_plus_propensity"`, which uses the modifier
set plus the treated-arm propensity score as a dimension-reduced summary
of the adjustment set. To compare EP-DR and EP-R directly at the wrapper
level, set `ep_targeting_style = c("dr", "r")`.

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

fit <- fit_cate(
  task,
  method = "ep",
  base_learner = Lrnr_mean$new(),
  cross_validate = TRUE,
  cv_control = list(V = 2),
  ep_targeting_style = c("dr", "r"),
  sieve_basis_grid = c(2, 4, 6, 8)
)

summary(fit)
} # }
```
