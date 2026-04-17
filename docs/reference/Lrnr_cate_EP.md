# Lrnr_cate_EP Class

This class defines the EP-learner of van der Laan et al. (2023) for
estimation of the conditional average treatment effect.

## Format

An R6 class with public methods to initialize the learner, create a
regression task, and access the base learner.

## Details

The EP-learner is a robust and doubly-robust meta-learner that inherits
desirable properties of both T-learner and DR-learner. The current
implementation exposes two CATE EP variants through `targeting_style`:
`"dr"` is the plug-in analogue of the DR-learner and `"r"` is the
plug-in analogue of the R-learner. The `"r"` variant uses an
overlap-weighted targeting objective that can be more stable under
severe overlap violations, but it generally targets an overlap-weighted
projection in reduced-modifier settings.

## Super classes

[`sl3::Lrnr_base`](https://tlverse.org/sl3/reference/Lrnr_base.html) -\>
[`hte3::Lrnr_hte`](https://larsvanderlaan.github.io/hte3/reference/Lrnr_hte.md)
-\> `Lrnr_cate_EP`

## Methods

### Public methods

- [`Lrnr_cate_EP$new()`](#method-Lrnr_cate_EP-new)

- [`Lrnr_cate_EP$get_pseudo_data()`](#method-Lrnr_cate_EP-get_pseudo_data)

- [`Lrnr_cate_EP$clone()`](#method-Lrnr_cate_EP-clone)

Inherited methods

- [`sl3::Lrnr_base$assert_trained()`](https://tlverse.org/sl3/reference/Lrnr_base.html#method-assert_trained)
- [`sl3::Lrnr_base$base_chain()`](https://tlverse.org/sl3/reference/Lrnr_base.html#method-base_chain)
- [`sl3::Lrnr_base$base_predict()`](https://tlverse.org/sl3/reference/Lrnr_base.html#method-base_predict)
- [`sl3::Lrnr_base$base_train()`](https://tlverse.org/sl3/reference/Lrnr_base.html#method-base_train)
- [`sl3::Lrnr_base$chain()`](https://tlverse.org/sl3/reference/Lrnr_base.html#method-chain)
- [`sl3::Lrnr_base$custom_chain()`](https://tlverse.org/sl3/reference/Lrnr_base.html#method-custom_chain)
- [`sl3::Lrnr_base$get_outcome_range()`](https://tlverse.org/sl3/reference/Lrnr_base.html#method-get_outcome_range)
- [`sl3::Lrnr_base$get_outcome_type()`](https://tlverse.org/sl3/reference/Lrnr_base.html#method-get_outcome_type)
- [`sl3::Lrnr_base$predict()`](https://tlverse.org/sl3/reference/Lrnr_base.html#method-predict)
- [`sl3::Lrnr_base$predict_fold()`](https://tlverse.org/sl3/reference/Lrnr_base.html#method-predict_fold)
- [`sl3::Lrnr_base$print()`](https://tlverse.org/sl3/reference/Lrnr_base.html#method-print)
- [`sl3::Lrnr_base$process_formula()`](https://tlverse.org/sl3/reference/Lrnr_base.html#method-process_formula)
- [`sl3::Lrnr_base$reparameterize()`](https://tlverse.org/sl3/reference/Lrnr_base.html#method-reparameterize)
- [`sl3::Lrnr_base$retrain()`](https://tlverse.org/sl3/reference/Lrnr_base.html#method-retrain)
- [`sl3::Lrnr_base$sample()`](https://tlverse.org/sl3/reference/Lrnr_base.html#method-sample)
- [`sl3::Lrnr_base$set_train()`](https://tlverse.org/sl3/reference/Lrnr_base.html#method-set_train)
- [`sl3::Lrnr_base$subset_covariates()`](https://tlverse.org/sl3/reference/Lrnr_base.html#method-subset_covariates)
- [`sl3::Lrnr_base$train()`](https://tlverse.org/sl3/reference/Lrnr_base.html#method-train)
- [`sl3::Lrnr_base$train_sublearners()`](https://tlverse.org/sl3/reference/Lrnr_base.html#method-train_sublearners)
- [`hte3::Lrnr_hte$check_treatment_type()`](https://larsvanderlaan.github.io/hte3/reference/Lrnr_hte.html#method-check_treatment_type)
- [`hte3::Lrnr_hte$get_modifiers()`](https://larsvanderlaan.github.io/hte3/reference/Lrnr_hte.html#method-get_modifiers)
- [`hte3::Lrnr_hte$make_metalearner_task()`](https://larsvanderlaan.github.io/hte3/reference/Lrnr_hte.html#method-make_metalearner_task)

------------------------------------------------------------------------

### Method `new()`

#### Usage

    Lrnr_cate_EP$new(
      base_learner,
      sieve_num_basis = NULL,
      sieve_interaction_order = 3,
      screen_basis_with_lasso = FALSE,
      targeting_style = c("dr", "r"),
      r_targeting_basis = c("v_plus_propensity", "full_w"),
      treatment_level = NULL,
      control_level = NULL,
      ...
    )

#### Arguments

- `base_learner`:

  A `sl3` learner object inheriting from
  [`Lrnr_base`](https://tlverse.org/sl3/reference/Lrnr_base.html) that
  specifies the base supervised learning algorithm used by the
  meta-learner.

- `sieve_num_basis`:

  The number of trigonometric basis functions to construct the
  EP-learner sieve space. This argument is passed as the argument
  `basisN` to the function
  [`sieve_preprocess`](https://rdrr.io/pkg/Sieve/man/sieve_preprocess.html).
  By default, `sieve_num_basis = ceiling((n^(1/3)*d)` where `n` denotes
  the sample size and `d` is the dimension of the first-stage EP sieve
  covariates. For wrapper-level EP cross-validation, use
  `sieve_basis_grid` in
  [`fit_cate()`](https://larsvanderlaan.github.io/hte3/reference/fit_cate.md)
  or
  [`fit_crr()`](https://larsvanderlaan.github.io/hte3/reference/fit_crr.md)
  to compare multiple basis sizes.

- `sieve_interaction_order`:

  The maximum interaction degree of tensor-product basis functions in
  the EP-learner sieve basis. Default is 3.
  `sieve_interaction_order = 1` corresponds to an additive sieve model
  and `sieve_interaction_order = 2` corresponds to a bi-additive sieve
  model. This argument is passed as the argument `interaction_order` to
  the function
  [`sieve_preprocess`](https://rdrr.io/pkg/Sieve/man/sieve_preprocess.html).

- `screen_basis_with_lasso`:

  EXPERIMENTAL. Whether to use a lasso-based DR-learner algorithm to
  screen sieve basis functions for EP-learner. There are no theoretical
  guarantees for EP-learner with `screen_basis_with_lasso = TRUE`.
  However, this argument may be useful in high dimensional settings. It
  also reduces the computational complexity of EP-learner, as the
  argument `sieve_num_basis` does not need to be externally
  cross-validated.

- `targeting_style`:

  One of `"dr"` or `"r"`. The default `"dr"` path keeps the current EP
  debiasing update. The `"r"` path uses an R-learner-style
  overlap-weighted targeting update and currently supports
  binary-treatment CATE tasks only. High-level wrappers such as
  [`fit_cate()`](https://larsvanderlaan.github.io/hte3/reference/fit_cate.md)
  compare both variants by expanding multiple `Lrnr_cate_EP` fits rather
  than by passing a vector here.

- `r_targeting_basis`:

  First-stage basis construction for `targeting_style = "r"`. The
  default `"v_plus_propensity"` builds the first-stage sieve on the
  modifier set `V` plus the treated-arm propensity score `e(W)`;
  `"full_w"` builds the first-stage EP-R sieve on the full confounder
  set `W`.

- `...`:

  Additional arguments to pass to the initialization function.

------------------------------------------------------------------------

### Method `get_pseudo_data()`

#### Usage

    Lrnr_cate_EP$get_pseudo_data(
      hte3_task,
      treatment_level = NULL,
      control_level = NULL,
      train = TRUE,
      ...
    )

#### Arguments

- `...`:

  Additional arguments to pass to the initialization function.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    Lrnr_cate_EP$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
