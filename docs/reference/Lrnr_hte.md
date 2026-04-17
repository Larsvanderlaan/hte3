# Lrnr_hte Class

This class serves as a general template for constructing meta-learners
to estimate heterogeneous treatment effects (HTEs).

## Format

An R6 class with public methods to initialize the learner, create a
regression task, and access the base learner.

## Super class

[`sl3::Lrnr_base`](https://tlverse.org/sl3/reference/Lrnr_base.html) -\>
`Lrnr_hte`

## Methods

### Public methods

- [`Lrnr_hte$new()`](#method-Lrnr_hte-new)

- [`Lrnr_hte$get_pseudo_data()`](#method-Lrnr_hte-get_pseudo_data)

- [`Lrnr_hte$check_treatment_type()`](#method-Lrnr_hte-check_treatment_type)

- [`Lrnr_hte$get_modifiers()`](#method-Lrnr_hte-get_modifiers)

- [`Lrnr_hte$make_metalearner_task()`](#method-Lrnr_hte-make_metalearner_task)

- [`Lrnr_hte$clone()`](#method-Lrnr_hte-clone)

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

------------------------------------------------------------------------

### Method `new()`

#### Usage

    Lrnr_hte$new(
      params,
      base_learner,
      transform_function = NULL,
      pseudo_outcome_type = c("continuous", "binomial", "quasibinomial"),
      pseudo_family = gaussian(),
      ...
    )

#### Arguments

- `params`:

  A list of parameters for the meta-learning algorithm.

- `base_learner`:

  A `sl3` learner object inheriting from
  [`Lrnr_base`](https://tlverse.org/sl3/reference/Lrnr_base.html) that
  specifies the base supervised learning algorithm used by the
  meta-learner.

- `transform_function`:

  A function to transform the predictions of the base learner. Default
  is the identity transform.

- `pseudo_outcome_type`:

  The outcome type of the pseudo-outcome used by the meta-learner.
  Options are `c("continuous", "binomial", "quasibinomial")`. Default is
  `"continuous"`. For example, the DR-learner, EP-learner, and R-learner
  of the CATE involve (weighted) least-squares regression using a
  pseudo-outcome with `pseudo_outcome_type = "continuous"`. The CRATE
  EP-learner involves performing weighted logistic regression using a
  pseudo-outcome taking values in 0,1 with
  `pseudo_outcome_type = "quasibinomial"`.

- `pseudo_family`:

  A [`family`](https://rdrr.io/r/stats/family.html) object specifying
  the loss function (involving pseudo-weights and pseudo-outcomes) used
  to fit `base_learner` in the meta-learner algorithm. Default is
  [`gaussian()`](https://rdrr.io/r/stats/family.html).

------------------------------------------------------------------------

### Method `get_pseudo_data()`

#### Usage

    Lrnr_hte$get_pseudo_data(hte3_task, ...)

#### Arguments

- `hte3_task`:

  A `hte3_Task` object containing the data and necessary information for
  heterogeneous treatment effect estimation.

- `...`:

  Additional arguments in `params` needed to compute the pseudo-data.

#### Returns

A list containing attributes `pseudo_outcome` and `pseudo_weights`.
(Used internally) Check Compatibility with Treatment Type

This method checks whether this meta-learner is compatible with the
`treatment` variable type in `hte3_Task`.

------------------------------------------------------------------------

### Method `check_treatment_type()`

#### Usage

    Lrnr_hte$check_treatment_type(hte3_task)

#### Arguments

- `hte3_task`:

  A `hte3_Task` object containing the data and necessary information for
  heterogeneous treatment effect estimation.

------------------------------------------------------------------------

### Method `get_modifiers()`

#### Usage

    Lrnr_hte$get_modifiers(hte3_task, return_matrix = FALSE)

------------------------------------------------------------------------

### Method `make_metalearner_task()`

#### Usage

    Lrnr_hte$make_metalearner_task(hte3_task, train = TRUE)

#### Arguments

- `hte3_task`:

  A `hte3_Task` object containing the data and necessary information for
  heterogeneous treatment effect estimation.

- `train`:

  Logical indicating whether to create the task for training or
  prediction. Default is `TRUE`. If `FALSE` then a `hte3_Task` object
  for prediction is returned with covariates being the effect
  `modifiers`.

#### Returns

A `sl3_Task` object containing the pseudo-data of `get_pseudo_data` and
`outcome_Type`=`pseudo_outcome_type`.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    Lrnr_hte$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
