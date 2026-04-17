# Lrnr_crr_EP Class

This class defines an efficient plug-in (EP) meta-learner for the
conditional relative average treatment effect (crr).

## Format

An R6 class with public methods to initialize the learner, create a
regression task, and access the base learner.

## Super classes

[`sl3::Lrnr_base`](https://tlverse.org/sl3/reference/Lrnr_base.html) -\>
[`hte3::Lrnr_hte`](https://larsvanderlaan.github.io/hte3/reference/Lrnr_hte.md)
-\> `Lrnr_crr_EP`

## Methods

### Public methods

- [`Lrnr_crr_EP$new()`](#method-Lrnr_crr_EP-new)

- [`Lrnr_crr_EP$get_pseudo_data()`](#method-Lrnr_crr_EP-get_pseudo_data)

- [`Lrnr_crr_EP$clone()`](#method-Lrnr_crr_EP-clone)

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

    Lrnr_crr_EP$new(
      base_learner,
      sieve_num_basis = NULL,
      sieve_interaction_order = 3,
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

  The number of trigonometric basis functions used to construct the CRR
  EP sieve space through
  [`sieve_preprocess`](https://rdrr.io/pkg/Sieve/man/sieve_preprocess.html).
  Use `sieve_num_basis` for one fixed EP fit and `sieve_basis_grid` in
  [`fit_crr()`](https://larsvanderlaan.github.io/hte3/reference/fit_crr.md)
  for wrapper-level EP basis-size selection.

- `sieve_interaction_order`:

  The maximum interaction degree of tensor-product basis functions in
  the EP sieve basis.

- `treatment_level`:

  Optional treated level used for the contrast.

- `control_level`:

  Optional control level used for the contrast.

------------------------------------------------------------------------

### Method `get_pseudo_data()`

#### Usage

    Lrnr_crr_EP$get_pseudo_data(
      hte3_task,
      treatment_level = NULL,
      control_level = NULL,
      train = TRUE,
      ...
    )

#### Arguments

- `treatment_level`:

  Optional treated level used for the contrast.

- `control_level`:

  Optional control level used for the contrast.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    Lrnr_crr_EP$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
