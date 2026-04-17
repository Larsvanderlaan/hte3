# Lrnr_crr_IPW Class

This class defines a inverse probability weighted (IPW) meta-learner of
the conditional relative average treatment effect (crr).

## Format

An R6 class with public methods to initialize the learner, create a
regression task, and access the base learner.

## Super classes

[`sl3::Lrnr_base`](https://tlverse.org/sl3/reference/Lrnr_base.html) -\>
[`hte3::Lrnr_hte`](https://larsvanderlaan.github.io/hte3/reference/Lrnr_hte.md)
-\> `Lrnr_crr_IPW`

## Methods

### Public methods

- [`Lrnr_crr_IPW$new()`](#method-Lrnr_crr_IPW-new)

- [`Lrnr_crr_IPW$get_pseudo_data()`](#method-Lrnr_crr_IPW-get_pseudo_data)

- [`Lrnr_crr_IPW$clone()`](#method-Lrnr_crr_IPW-clone)

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

    Lrnr_crr_IPW$new(
      base_learner,
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

------------------------------------------------------------------------

### Method `get_pseudo_data()`

#### Usage

    Lrnr_crr_IPW$get_pseudo_data(
      hte3_task,
      treatment_level = NULL,
      control_level = NULL,
      train = TRUE,
      ...
    )

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    Lrnr_crr_IPW$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
