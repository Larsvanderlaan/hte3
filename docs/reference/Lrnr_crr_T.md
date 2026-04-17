# Lrnr_crr_T Class

This class constructs a T-learner for estimation of the conditional
relative average treatment effect (crr)

## Format

An R6 class with public methods to initialize the learner, create a
regression task, and access the base learner.

## Details

The learner first estimates arm-specific outcome regressions as
functions of the adjustment set `W`. By default, it then regresses the
resulting first-stage log-risk-ratio contrast onto the modifier set `V`.
When `stratify_by_treatment = FALSE`, the first-stage outcome model is
fit as a pooled S-learner-style regression and evaluated
counterfactually under each treatment level.

## Super classes

[`sl3::Lrnr_base`](https://tlverse.org/sl3/reference/Lrnr_base.html) -\>
[`hte3::Lrnr_hte`](https://larsvanderlaan.github.io/hte3/reference/Lrnr_hte.md)
-\> `Lrnr_crr_T`

## Methods

### Public methods

- [`Lrnr_crr_T$new()`](#method-Lrnr_crr_T-new)

- [`Lrnr_crr_T$get_pseudo_data()`](#method-Lrnr_crr_T-get_pseudo_data)

- [`Lrnr_crr_T$make_metalearner_task()`](#method-Lrnr_crr_T-make_metalearner_task)

- [`Lrnr_crr_T$clone()`](#method-Lrnr_crr_T-clone)

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

------------------------------------------------------------------------

### Method `new()`

#### Usage

    Lrnr_crr_T$new(
      base_learner,
      treatment_level = 1,
      control_level = 0,
      stratify_by_treatment = TRUE,
      second_stage_regression = TRUE,
      ...
    )

#### Arguments

- `base_learner`:

  A `sl3` learner object inheriting from
  [`Lrnr_base`](https://tlverse.org/sl3/reference/Lrnr_base.html) that
  specifies the base supervised learning algorithm used by the
  meta-learner.

- `stratify_by_treatment`:

  Logical indicating whether to estimate outcome regressions separately
  in each treatment arm (i.e., T-learner) or with a pooled
  S-learner-style outcome model.

- `second_stage_regression`:

  Logical indicating whether to regress the first-stage log-risk-ratio
  contrast onto the modifier set. This defaults to `TRUE`. Setting it to
  `FALSE` is only supported when the modifier set and confounder set are
  the same.

------------------------------------------------------------------------

### Method `get_pseudo_data()`

#### Usage

    Lrnr_crr_T$get_pseudo_data(hte3_task, ...)

------------------------------------------------------------------------

### Method `make_metalearner_task()`

#### Usage

    Lrnr_crr_T$make_metalearner_task(hte3_task, train = TRUE)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    Lrnr_crr_T$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
