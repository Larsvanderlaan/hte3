# Internal use. Converts a single outcome learner into a multivariate outcome learner that predicts a matrix of predictions obtained by evaluating the single outcome learner at each possible value of variable_stratify.

Internal use. Converts a single outcome learner into a multivariate
outcome learner that predicts a matrix of predictions obtained by
evaluating the single outcome learner at each possible value of
variable_stratify.

Internal use. Converts a single outcome learner into a multivariate
outcome learner that predicts a matrix of predictions obtained by
evaluating the single outcome learner at each possible value of
variable_stratify.

## Super classes

[`sl3::Lrnr_base`](https://tlverse.org/sl3/reference/Lrnr_base.html) -\>
[`sl3::Lrnr_stratified`](https://tlverse.org/sl3/reference/Lrnr_stratified.html)
-\> `Lrnr_stratified_multivariate`

## Methods

### Public methods

- [`Lrnr_stratified_multivariate$new()`](#method-Lrnr_stratified_multivariate-new)

- [`Lrnr_stratified_multivariate$clone()`](#method-Lrnr_stratified_multivariate-clone)

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

    Lrnr_stratified_multivariate$new(learner, variable_stratify, ...)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    Lrnr_stratified_multivariate$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
