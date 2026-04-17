# Task object for meta-learners in the point-treatment setting.

Constructs a `hte3_Task` object for meta-learners in the point-treatment
setting containing relevant data and nuisance function estimators.

## Usage

``` r
make_hte3_Task_tx(
  data,
  modifiers,
  confounders,
  treatment,
  outcome,
  id = NULL,
  weights = NULL,
  treatment_type = c("default", "binomial", "categorical", "continuous"),
  pi.hat = NULL,
  mu.hat = NULL,
  m.hat = NULL,
  learner_pi = get_autoML(),
  learner_mu = get_autoML(),
  learner_m = NULL,
  multinomial_learner = Lrnr_independent_binomial,
  cross_fit_and_cv = TRUE,
  folds = 10,
  warn = TRUE,
  for_prediction = FALSE,
  ...
)
```

## Arguments

- data:

  A named data frame or data.table containing treatment effect
  modifiers, potential confounders, treatment, and outcome. Optionally,
  the dataset can contain weights and subject IDs. See the `data`
  argument of
  [`sl3_Task`](https://tlverse.org/sl3/reference/sl3_Task.html) for
  further options.

- modifiers:

  A character vector of variable names in `colnames(data)` for treatment
  effect moderators.

- confounders:

  A character vector of variable names in `colnames(data)` for potential
  confounders `W` for which to adjust.

- treatment:

  A character specifying the variable name in `colnames(data)` for the
  numeric treatment assignment `A`.

- outcome:

  A character specifying the variable name in `colnames(data)` for the
  outcome variable `Y`.

- id:

  An (optional) character specifying the variable name in
  `colnames(data)` for observation IDs.

- weights:

  An (optional) character specifying the variable name in
  `colnames(data)` for observation weights.

- pi.hat:

  An (optional) numeric matrix of dimension `n by nlevels(A)` containing
  estimates of the propensity score `a -> pi(a | W_i)` at each treatment
  level, where column `j` corresponds to treatment level
  `sort(unique(A))[j]`. This argument can be used by DR-type learners.
  Alternatively, the estimates can be learned internally by passing an
  `sl3_Learner` to the `learner_pi` argument.

- mu.hat:

  An (optional) numeric matrix of dimension `n by nlevels(A)` containing
  estimates of the outcome regression `a -> mu(a, W_i)` at each
  treatment level, where column `j` corresponds to treatment level
  `sort(unique(A))[j]`. This argument can be used by DR-type learners.
  Alternatively, the estimates can be learned internally by passing an
  `sl3_Learner` to the `learner_m` argument.

- m.hat:

  An (optional) numeric vector of size `n` containing estimates of the
  conditional mean outcome `E[Y | W_i]`. This argument can be used by
  R-type learners. Alternatively, the estimates can be learned
  internally by passing an `sl3_Learner` to the `learner_m` argument.

- learner_pi:

  A binomial `sl3_Learner` or `Stack` object specifying the learning
  algorithm to estimate the propensity score `w -> pi(a_star | w)` at a
  given treatment level `a_star`. During training, this learner is
  passed an `sl3_Task` object `task` that contains a feature matrix
  `task$X` of `confounders` and an outcome vector `task$Y` corresponding
  to the binary treatment indicator `1(A = a_star)`. By default,
  `cross_fit_and_cv = TRUE`, and `learner_m` is fit with 10-fold
  cross-fitting using
  [`make_cross_fitted`](https://larsvanderlaan.github.io/hte3/reference/make_cross_fitted.md).
  If a `Stack` object, the best model is selected using
  cross-validation.

- learner_mu:

  A `sl3_Learner` or `sl3_Task` object specifying the learning algorithm
  to estimate the outcome regression `w -> mu(a_star, w)` at a given
  treatment level `a_star`. The outcome regression is estimated via a
  treatment-stratified regression. During training, this learner is
  passed an `sl3_Task` object `task` that is subsetted to variables with
  treatment `A = a_star`, and contains a feature matrix `task$X` of
  `confounders` and an outcome vector `task$Y` corresponding to
  `outcome`. By default, `cross_fit_and_cv = TRUE`, and `learner_m` is
  fit with 10-fold cross-fitting using
  [`make_cross_fitted`](https://larsvanderlaan.github.io/hte3/reference/make_cross_fitted.md).
  If a `Stack` object, the best model is selected using
  cross-validation.

- learner_m:

  A `sl3_Learner` or `sl3_Task` object specifying the learning algorithm
  to estimate the conditional mean outcome `m`. During training, this
  learner is passed an `sl3_Task` object `task` that contains a feature
  matrix `task$X` of `confounders` and an outcome vector `task$Y`
  corresponding to `outcome`. By default, `cross_fit_and_cv = TRUE`, and
  `learner_m` is fit with 10-fold cross-fitting using
  [`make_cross_fitted`](https://larsvanderlaan.github.io/hte3/reference/make_cross_fitted.md).
  If a `Stack` object, the best model is selected using
  cross-validation.

- multinomial_learner:

  A `multinomial` `Lrnr_base` object to convert `learner_pi` from a
  `binomial` learner to a `multinomial` learner.

- cross_fit_and_cv:

  Whether to cross-fit the specified nuisance learners by applying
  `learner <- make_cross_fitted(learner)`.

- for_prediction:

  A `boolean` of whether to return an hte3_Task without the nuisance
  training and estimates. This can be useful when you wish to construct
  a task for predicting on a new dataset.

- ...:

  Additional arguments to pass to the
  [`sl3_Task`](https://tlverse.org/sl3/reference/sl3_Task.html) and
  `tmle_Task` constructors.

## Value

A `hte3_Task` object for the point-treatment data-structure
