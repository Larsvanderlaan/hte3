# EP-Learner Sieve Tuning

This vignette focuses on the EP-learner. It explains what the `Sieve`
package is doing inside `hte3`, how the sieve basis is constructed and
ordered, how the two CATE EP variants differ, how wrapper-level EP
basis-size selection works, and how to reproduce the paper-style sieve
cross-validation grids.

``` r
library(hte3)
library(sl3)

data <- hte3_example_data(n = 120, seed = 4)

task <- make_hte3_Task_tx(
  data = data,
  modifiers = c("W1", "W2"),
  confounders = c("W1", "W2", "W3"),
  treatment = "A",
  outcome = "Y",
  learner_pi = Lrnr_mean$new(),
  learner_mu = Lrnr_mean$new(),
  learner_m = Lrnr_mean$new(),
  cross_fit_and_cv = FALSE
)

crr_task <- make_hte3_Task_tx(
  data = data,
  modifiers = c("W1", "W2", "W3"),
  confounders = c("W1", "W2", "W3"),
  treatment = "A",
  outcome = "Y_binary",
  learner_pi = Lrnr_mean$new(),
  learner_mu = Lrnr_mean$new(),
  learner_m = Lrnr_mean$new(),
  cross_fit_and_cv = FALSE
)
```

## What `Sieve` Does in `hte3`

For CATE, `hte3` now exposes two EP variants through
`Lrnr_cate_EP(..., targeting_style = ...)`:

- `targeting_style = "dr"` is the original EP update and acts as a
  plug-in analogue of the DR-learner.
- `targeting_style = "r"` uses an overlap-weighted targeting objective
  and acts as a plug-in analogue of the R-learner.

In both cases, the `Sieve` package is used only to build the basis
matrix. The EP update itself is then fit by `glmnet` inside
`Lrnr_cate_EP` and `Lrnr_crr_EP`.

Before evaluating the basis, `Sieve` rescales each modifier to `[0, 1]`.
In the current `hte3` implementation, the sieve family is fixed to
cosine basis functions, so the univariate building blocks are
`1, cos(pi x), cos(2 pi x), ...`. Multivariate basis functions are
tensor products of those univariate terms across modifiers.

The package does not currently expose alternative `Sieve` basis families
at the wrapper or learner level. The practical tuning knobs are:

- `sieve_num_basis`: sieve size for one fixed EP fit
- `sieve_basis_grid`: basis-size grid for wrapper-level EP
  cross-validation
- `sieve_interaction_order`: maximum interaction order
- `screen_basis_with_lasso`: experimental CATE-only alternative to
  explicit sieve-grid cross-validation
- `targeting_style`: choose between EP-DR and EP-R for CATE
- `r_targeting_basis`: choose the EP-R first-stage basis, either full
  `W` or the reduced `V + e(W)` construction

If `sieve_num_basis = NULL`, `hte3` uses the default
`ceiling(n^(1/3) * d)`, where `n` is the sample size and `d` is the
dimension of the first-stage EP sieve covariates.

## EP-DR versus EP-R

The important difference is the first-stage targeting objective:

- EP-DR uses the original inverse-propensity-weighted targeting step on
  a sieve built from the modifier set `V`.
- EP-R starts from the plug-in effect
  `tau0.hat(W) = mu.hat(1, W) - mu.hat(0, W)` and then fits an
  overlap-weighted residual update.

For EP-R, you can choose between two first-stage basis constructions:

- `r_targeting_basis = "v_plus_propensity"`: the default; build the
  first-stage sieve on the modifier set `V` plus the treated-arm
  propensity score `e(W)`
- `r_targeting_basis = "full_w"`: build the first-stage sieve on the
  full confounder set `W`

The `v_plus_propensity` option is the dimension-reduced EP-R path: it
uses the scalar treated-arm propensity score as a 1D summary of `W`
while still allowing the first-stage sieve to interact that score with
the modifier set. For this mode, the effective first-stage interaction
order is always at least 2.

The practical target distinction is:

- EP-DR keeps the usual `V`-conditional CATE target in the supported
  binary/categorical setting.
- EP-R can be more stable under severe overlap violations, but in
  reduced-modifier settings it generally targets an overlap-weighted
  projection rather than the unweighted `E[Y(1) - Y(0) | V]`.

EP-R currently supports binary-treatment CATE tasks only. `Lrnr_crr_EP`
remains unchanged.

## Wrapper-Level EP Tuning

At the high-level wrapper, EP tuning now follows two simple rules:

- Use `sieve_num_basis` when `cross_validate = FALSE` and you want one
  fixed sieve size.
- Use `sieve_basis_grid` when `cross_validate = TRUE` and you want the
  wrapper to compare multiple EP candidates.

For CATE EP learners, wrapper-level tuning still uses the DR selector
loss from
[`cross_validate_cate()`](https://larsvanderlaan.github.io/hte3/reference/cross_validate_cate.md),
even when `targeting_style = "r"`.

For example:

``` r
ep_cv <- fit_cate(
  hte_task(
    data = data,
    modifiers = c("W1", "W2"),
    confounders = c("W1", "W2", "W3"),
    treatment = "A",
    outcome = "Y",
    propensity_learner = Lrnr_mean$new(),
    outcome_learner = Lrnr_mean$new(),
    mean_learner = Lrnr_mean$new(),
    cross_fit = FALSE
  ),
  method = "ep",
  base_learner = Lrnr_mean$new(),
  cross_validate = TRUE,
  cv_control = list(V = 2),
  ep_targeting_style = "dr",
  sieve_basis_grid = c(2, 4, 6, 8)
)

summary(ep_cv)
#> <summary.hte3_model>
#>   Target: CATE
#>   Engine: sl3
#>   Method: ep
#>   Cross-validated: yes
#>   Rows: 120
#>   Selected candidate: ep[style=dr, basis=8, interaction=3]
#>   Selected method: ep
#>   Selected EP basis size: 8
#>   EP targeting style: dr
#>   EP basis grid: 2, 4, 6, 8
#>   Modifiers: W1, W2
#>   Treatment variable: A
```

You can also ask the wrapper to compare EP-DR and EP-R directly:

``` r
ep_style_cv <- suppressWarnings(fit_cate(
  hte_task(
    data = data,
    modifiers = c("W1", "W2"),
    confounders = c("W1", "W2", "W3"),
    treatment = "A",
    outcome = "Y",
    propensity_learner = Lrnr_mean$new(),
    outcome_learner = Lrnr_mean$new(),
    mean_learner = Lrnr_mean$new(),
    cross_fit = FALSE
  ),
  method = "ep",
  base_learner = Lrnr_mean$new(),
  cross_validate = TRUE,
  cv_control = list(V = 2),
  sieve_basis_grid = c(2, 4, 6, 8),
  ep_targeting_style = c("dr", "r")
))
#> Error : x has missing values; consider using makeX() to impute them
#> Error : x has missing values; consider using makeX() to impute them
#> Error : x has missing values; consider using makeX() to impute them
#> Error : x has missing values; consider using makeX() to impute them
#> Error : x has missing values; consider using makeX() to impute them
#> Error : x has missing values; consider using makeX() to impute them
#> Error : x has missing values; consider using makeX() to impute them
#> Error : x has missing values; consider using makeX() to impute them
#> Error : x has missing values; consider using makeX() to impute them

summary(ep_style_cv)
#> <summary.hte3_model>
#>   Target: CATE
#>   Engine: sl3
#>   Method: ep
#>   Cross-validated: yes
#>   Rows: 120
#>   Selected candidate: candidate_2
#>   Selected method: candidate_2
#>   Modifiers: W1, W2
#>   Treatment variable: A
```

If `sieve_basis_grid = NULL`, the wrapper uses the heuristic default
grid `c(d, 2d, 4d, 6d, 8d)`, where `d` is the first-stage EP sieve
dimension for the chosen variant. The fitted model summary reports the
selected candidate, the selected EP targeting style, the selected EP
basis size, and the EP grid that was compared.

## How to Choose the Grid

Use the grid to control the flexibility-compute tradeoff:

- Smaller grids are faster and are usually safer in moderate sample
  sizes.
- Larger grids give the EP correction more flexibility but can become
  noisier in higher-dimensional modifier spaces.
- Lower `sieve_interaction_order` values reduce complexity when the
  modifier dimension is large.

In practice:

- start with the default heuristic grid when you want a simple first
  pass
- use `sieve_basis_grid = (3:8) * d` to match the low-dimensional paper
  grid
- use smaller additive grids such as `(1:6) * d` with
  `sieve_interaction_order = 1` when the modifier dimension is high

## How the Basis Is Ordered

`Sieve` orders the multivariate basis by increasing index-product
complexity. In practice, this means the constant term comes first, then
lower-frequency univariate terms, then higher-frequency and interaction
terms.

The argument `sieve_interaction_order` limits how many coordinates can
vary inside a single tensor-product basis term:

- `sieve_interaction_order = 1` gives an additive sieve
- `sieve_interaction_order = 2` allows pairwise interactions
- `sieve_interaction_order = 3` allows up to three-way interactions

So with two modifiers, the early cosine basis terms look like the
constant, then terms like `cos(pi x1)` and `cos(pi x2)`, then
higher-frequency terms, then interaction terms such as
`cos(pi x1) * cos(pi x2)`.

## Fitting a Single EP-Learner

This is the most direct low-level EP fit:

``` r
ep_fit <- Lrnr_cate_EP$new(
  base_learner = Lrnr_mean$new(),
  sieve_num_basis = 4,
  sieve_interaction_order = 2
)$train(task)
#> Warning in learner$predict_fold(learner_task, fold_number):
#> Lrnr_independent_binomial is not cv-aware: self$predict_fold reverts to
#> self$predict
#> Warning in learner$predict_fold(learner_task, fold_number): strat_A_Lrnr_mean
#> is not cv-aware: self$predict_fold reverts to self$predict
#> Warning in learner$predict_fold(learner_task, fold_number): Lrnr_mean is not
#> cv-aware: self$predict_fold reverts to self$predict

head(predict_hte3(ep_fit, data))
#> [1] 1.429171 1.429171 1.429171 1.429171 1.429171 1.429171
```

This is the corresponding EP-R fit with the dimension-reduced
first-stage basis:

``` r
ep_r_fit <- Lrnr_cate_EP$new(
  base_learner = Lrnr_mean$new(),
  sieve_num_basis = 4,
  sieve_interaction_order = 1,
  targeting_style = "r",
  r_targeting_basis = "v_plus_propensity"
)$train(task)

head(predict_hte3(ep_r_fit, data))
```

In practice, the paper often paired EP-learner with richer base-learner
stacks. The examples below keep `Lrnr_mean$new()` so the sieve tuning
pattern is easy to read; you can replace `base_learner` with any `sl3`
learner or stack.

## Paper-Faithful CATE Sieve CV

The main low-dimensional CATE simulations in the paper used a sieve grid
of `3d, 4d, ..., 8d` with `sieve_interaction_order = 3`.

``` r
d <- length(task$npsem$modifiers$variables)
paper_cate_basis_grid <- (3:8) * d

paper_cate_ep_stack <- lapply(paper_cate_basis_grid, function(basis_size) {
  Lrnr_cate_EP$new(
    base_learner = Lrnr_mean$new(),
    sieve_num_basis = basis_size,
    sieve_interaction_order = 3
  )
})

paper_cate_cv <- cross_validate_cate(
  paper_cate_ep_stack,
  task,
  cv_control = list(V = 5)
)
#> Warning in learner$predict_fold(learner_task, fold_number):
#> Lrnr_independent_binomial is not cv-aware: self$predict_fold reverts to
#> self$predict
#> Warning in learner$predict_fold(learner_task, fold_number): strat_A_Lrnr_mean
#> is not cv-aware: self$predict_fold reverts to self$predict
#> Warning in learner$predict_fold(learner_task, fold_number): Lrnr_mean is not
#> cv-aware: self$predict_fold reverts to self$predict
#> Warning in learner$predict_fold(learner_task, fold_number):
#> Lrnr_independent_binomial is not cv-aware: self$predict_fold reverts to
#> self$predict
#> Warning in learner$predict_fold(learner_task, fold_number): strat_A_Lrnr_mean
#> is not cv-aware: self$predict_fold reverts to self$predict
#> Warning in learner$predict_fold(learner_task, fold_number): Lrnr_mean is not
#> cv-aware: self$predict_fold reverts to self$predict
#> Warning in learner$predict_fold(learner_task, fold_number):
#> Lrnr_independent_binomial is not cv-aware: self$predict_fold reverts to
#> self$predict
#> Warning in learner$predict_fold(learner_task, fold_number): strat_A_Lrnr_mean
#> is not cv-aware: self$predict_fold reverts to self$predict
#> Warning in learner$predict_fold(learner_task, fold_number): Lrnr_mean is not
#> cv-aware: self$predict_fold reverts to self$predict
#> Warning in learner$predict_fold(learner_task, fold_number):
#> Lrnr_independent_binomial is not cv-aware: self$predict_fold reverts to
#> self$predict
#> Warning in learner$predict_fold(learner_task, fold_number): strat_A_Lrnr_mean
#> is not cv-aware: self$predict_fold reverts to self$predict
#> Warning in learner$predict_fold(learner_task, fold_number): Lrnr_mean is not
#> cv-aware: self$predict_fold reverts to self$predict
#> Warning in learner$predict_fold(learner_task, fold_number):
#> Lrnr_independent_binomial is not cv-aware: self$predict_fold reverts to
#> self$predict
#> Warning in learner$predict_fold(learner_task, fold_number): strat_A_Lrnr_mean
#> is not cv-aware: self$predict_fold reverts to self$predict
#> Warning in learner$predict_fold(learner_task, fold_number): Lrnr_mean is not
#> cv-aware: self$predict_fold reverts to self$predict
#> Warning in learner$predict_fold(learner_task, fold_number):
#> Lrnr_independent_binomial is not cv-aware: self$predict_fold reverts to
#> self$predict
#> Warning in learner$predict_fold(learner_task, fold_number): strat_A_Lrnr_mean
#> is not cv-aware: self$predict_fold reverts to self$predict
#> Warning in learner$predict_fold(learner_task, fold_number): Lrnr_mean is not
#> cv-aware: self$predict_fold reverts to self$predict

paper_cate_basis_grid
#> [1]  6  8 10 12 14 16
paper_cate_cv$coefficients
#>  ep[style=dr, basis=6, interaction=3]  ep[style=dr, basis=8, interaction=3] 
#>                                     0                                     1 
#> ep[style=dr, basis=10, interaction=3] ep[style=dr, basis=12, interaction=3] 
#>                                     0                                     0 
#> ep[style=dr, basis=14, interaction=3] ep[style=dr, basis=16, interaction=3] 
#>                                     0                                     0
```

The high-level wrapper
`fit_cate(..., method = "ep", cross_validate = TRUE)` uses the heuristic
default grid `c(d, 2d, 4d, 6d, 8d)` unless you provide
`sieve_basis_grid` yourself. If you want paper-faithful tuning, set
`sieve_basis_grid = (3:8) * d` or use an explicit EP stack like the one
above. The same wrapper can compare EP-DR and EP-R directly with
`ep_targeting_style = c("dr", "r")`.

For the high-dimensional CATE simulations, the paper used a more
conservative additive sieve:

``` r
d <- length(task$npsem$modifiers$variables)
paper_cate_highdim_grid <- (1:6) * d

paper_cate_highdim_ep_stack <- lapply(
  paper_cate_highdim_grid,
  function(basis_size) {
    Lrnr_cate_EP$new(
      base_learner = Lrnr_mean$new(),
      sieve_num_basis = basis_size,
      sieve_interaction_order = 1
    )
  }
)

paper_cate_highdim_cv <- cross_validate_cate(
  paper_cate_highdim_ep_stack,
  task,
  cv_control = list(V = 5)
)
```

## Paper-Faithful CRR Sieve CV

The main CRR simulations used the same low-dimensional sieve grid
`3d, 4d, ..., 8d` with `sieve_interaction_order = 3`.

``` r
d_crr <- length(crr_task$npsem$modifiers$variables)
paper_crr_basis_grid <- (3:8) * d_crr

paper_crr_ep_stack <- lapply(paper_crr_basis_grid, function(basis_size) {
  Lrnr_crr_EP$new(
    base_learner = Lrnr_mean$new(),
    sieve_num_basis = basis_size,
    sieve_interaction_order = 3
  )
})

paper_crr_cv <- cross_validate_crr(
  paper_crr_ep_stack,
  crr_task,
  cv_control = list(V = 5)
)
#> Warning in learner$predict_fold(learner_task, fold_number):
#> Lrnr_independent_binomial is not cv-aware: self$predict_fold reverts to
#> self$predict
#> Warning in learner$predict_fold(learner_task, fold_number): strat_A_Lrnr_mean
#> is not cv-aware: self$predict_fold reverts to self$predict
#> Warning in learner$predict_fold(learner_task, fold_number): Lrnr_mean is not
#> cv-aware: self$predict_fold reverts to self$predict
#> Warning in learner$predict_fold(learner_task, fold_number):
#> Lrnr_independent_binomial is not cv-aware: self$predict_fold reverts to
#> self$predict
#> Warning in learner$predict_fold(learner_task, fold_number): strat_A_Lrnr_mean
#> is not cv-aware: self$predict_fold reverts to self$predict
#> Warning in learner$predict_fold(learner_task, fold_number): Lrnr_mean is not
#> cv-aware: self$predict_fold reverts to self$predict
#> Warning: from glmnet C++ code (error code -1); Convergence for 1th lambda value
#> not reached after maxit=100000 iterations; solutions for larger lambdas
#> returned
#> Warning in getcoef(fit, nvars, nx, vnames): an empty model has been returned;
#> probably a convergence issue
#> Warning in learner$predict_fold(learner_task, fold_number):
#> Lrnr_independent_binomial is not cv-aware: self$predict_fold reverts to
#> self$predict
#> Warning in learner$predict_fold(learner_task, fold_number): strat_A_Lrnr_mean
#> is not cv-aware: self$predict_fold reverts to self$predict
#> Warning in learner$predict_fold(learner_task, fold_number): Lrnr_mean is not
#> cv-aware: self$predict_fold reverts to self$predict
#> Warning in learner$predict_fold(learner_task, fold_number):
#> Lrnr_independent_binomial is not cv-aware: self$predict_fold reverts to
#> self$predict
#> Warning in learner$predict_fold(learner_task, fold_number): strat_A_Lrnr_mean
#> is not cv-aware: self$predict_fold reverts to self$predict
#> Warning in learner$predict_fold(learner_task, fold_number): Lrnr_mean is not
#> cv-aware: self$predict_fold reverts to self$predict
#> Warning in learner$predict_fold(learner_task, fold_number):
#> Lrnr_independent_binomial is not cv-aware: self$predict_fold reverts to
#> self$predict
#> Warning in learner$predict_fold(learner_task, fold_number): strat_A_Lrnr_mean
#> is not cv-aware: self$predict_fold reverts to self$predict
#> Warning in learner$predict_fold(learner_task, fold_number): Lrnr_mean is not
#> cv-aware: self$predict_fold reverts to self$predict
#> Warning: from glmnet C++ code (error code -1); Convergence for 1th lambda value
#> not reached after maxit=100000 iterations; solutions for larger lambdas
#> returned
#> Warning in getcoef(fit, nvars, nx, vnames): an empty model has been returned;
#> probably a convergence issue
#> Warning in learner$predict_fold(learner_task, fold_number):
#> Lrnr_independent_binomial is not cv-aware: self$predict_fold reverts to
#> self$predict
#> Warning in learner$predict_fold(learner_task, fold_number): strat_A_Lrnr_mean
#> is not cv-aware: self$predict_fold reverts to self$predict
#> Warning in learner$predict_fold(learner_task, fold_number): Lrnr_mean is not
#> cv-aware: self$predict_fold reverts to self$predict
#> Warning in learner$predict_fold(learner_task, fold_number):
#> Lrnr_independent_binomial is not cv-aware: self$predict_fold reverts to
#> self$predict
#> Warning in learner$predict_fold(learner_task, fold_number): strat_A_Lrnr_mean
#> is not cv-aware: self$predict_fold reverts to self$predict
#> Warning in learner$predict_fold(learner_task, fold_number): Lrnr_mean is not
#> cv-aware: self$predict_fold reverts to self$predict
#> Warning in initialize(...): The supplied outcome_type (continuous) does not
#> correspond to the detected type (categorical). Ensure outcome_type is specified
#> appropriately.

paper_crr_basis_grid
#> [1]  9 12 15 18 21 24
paper_crr_cv$coefficients
#>  ep[basis=9, interaction=3] ep[basis=12, interaction=3] 
#>                           0                           0 
#> ep[basis=15, interaction=3] ep[basis=18, interaction=3] 
#>                           0                           0 
#> ep[basis=21, interaction=3] ep[basis=24, interaction=3] 
#>                           0                           1
```

## Notes

- For wrapper-level EP selection, `sieve_basis_grid` is the user-visible
  knob for basis-size cross-validation and `sieve_num_basis` is the
  user-visible knob for a single fixed EP fit.
- The paper-faithful path is explicit sieve-grid selection with
  [`cross_validate_cate()`](https://larsvanderlaan.github.io/hte3/reference/cross_validate_cate.md)
  or
  [`cross_validate_crr()`](https://larsvanderlaan.github.io/hte3/reference/cross_validate_crr.md).
- For CATE, `screen_basis_with_lasso = TRUE` is an experimental
  alternative when explicit sieve CV is too expensive.
- For broader low-level `sl3` portfolio control beyond EP-learner, see
  the `advanced-sl3` vignette.
