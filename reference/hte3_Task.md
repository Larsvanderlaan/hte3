# Task object for meta-learners in causal data structures.

Constructs a `hte3_Task` object for meta-learners in causal data
structures, containing relevant data and nuisance function estimators.

## Format

`R6Class` object.

## Value

An object of class `hte3_Task`.

## Details

Class for Storing Data, NPSEM, and Nuisance Function Estimators for hte3
Learners

This class inherits from
[`sl3_Task`](https://tlverse.org/sl3/reference/sl3_Task.html) and
[`tmle3_Task`](http://tlverse.org/tmle3/reference/tmle3_Task.md). In
addition to all the methods supported by
[`sl3_Task`](https://tlverse.org/sl3/reference/sl3_Task.html) and
[`tmle3_Task`](http://tlverse.org/tmle3/reference/tmle3_Task.md), it
supports the following functionalities specific to the heterogeneous
treatment effect estimation (hte3) framework.

## Super classes

[`sl3::sl3_Task`](https://tlverse.org/sl3/reference/sl3_Task.html) -\>
[`tmle3::tmle3_Task`](http://tlverse.org/tmle3/reference/tmle3_Task.md)
-\> `hte3_Task`

## Methods

### Public methods

- [`hte3_Task$new()`](#method-hte3_Task-new)

- [`hte3_Task$add_nuisance_estimator()`](#method-hte3_Task-add_nuisance_estimator)

- [`hte3_Task$get_nuisance_estimates()`](#method-hte3_Task-get_nuisance_estimates)

- [`hte3_Task$next_in_chain()`](#method-hte3_Task-next_in_chain)

- [`hte3_Task$subset_task()`](#method-hte3_Task-subset_task)

- [`hte3_Task$clone()`](#method-hte3_Task-clone)

Inherited methods

- [`sl3::sl3_Task$add_columns()`](https://tlverse.org/sl3/reference/sl3_Task.html#method-add_columns)
- [`sl3::sl3_Task$add_interactions()`](https://tlverse.org/sl3/reference/sl3_Task.html#method-add_interactions)
- [`sl3::sl3_Task$get_data()`](https://tlverse.org/sl3/reference/sl3_Task.html#method-get_data)
- [`sl3::sl3_Task$get_folds()`](https://tlverse.org/sl3/reference/sl3_Task.html#method-get_folds)
- [`sl3::sl3_Task$get_node()`](https://tlverse.org/sl3/reference/sl3_Task.html#method-get_node)
- [`sl3::sl3_Task$has_node()`](https://tlverse.org/sl3/reference/sl3_Task.html#method-has_node)
- [`sl3::sl3_Task$offset_transformed()`](https://tlverse.org/sl3/reference/sl3_Task.html#method-offset_transformed)
- [`sl3::sl3_Task$revere_fold_task()`](https://tlverse.org/sl3/reference/sl3_Task.html#method-revere_fold_task)
- [`tmle3::tmle3_Task$generate_counterfactual_task()`](http://tlverse.org/tmle3/reference/tmle3_Task.html#method-generate_counterfactual_task)
- [`tmle3::tmle3_Task$get_node_bounds()`](http://tlverse.org/tmle3/reference/tmle3_Task.html#method-get_node_bounds)
- [`tmle3::tmle3_Task$get_regression_task()`](http://tlverse.org/tmle3/reference/tmle3_Task.html#method-get_regression_task)
- [`tmle3::tmle3_Task$get_tmle_node()`](http://tlverse.org/tmle3/reference/tmle3_Task.html#method-get_tmle_node)
- [`tmle3::tmle3_Task$print()`](http://tlverse.org/tmle3/reference/tmle3_Task.html#method-print)
- [`tmle3::tmle3_Task$scale()`](http://tlverse.org/tmle3/reference/tmle3_Task.html#method-scale)
- [`tmle3::tmle3_Task$unscale()`](http://tlverse.org/tmle3/reference/tmle3_Task.html#method-unscale)

------------------------------------------------------------------------

### Method `new()`

#### Usage

    hte3_Task$new(data, npsem, likelihood = NULL, ...)

#### Arguments

- `data`:

  A named data frame or data.table containing treatment effect
  modifiers, potential confounders, treatment, outcome, and optionally
  weights and subject IDs. See the
  [`sl3_Task`](https://tlverse.org/sl3/reference/sl3_Task.html)
  documentation for further options.

- `npsem`:

  A list containing
  [`tmle3_Node`](http://tlverse.org/tmle3/reference/tmle3_Node.md)
  objects containing the non-parametric structural equation model
  (NPSEM) specifying the relationship between variables in the treatment
  effect estimation framework.

- `likelihood`:

  An [`Likelihood`](http://tlverse.org/tmle3/reference/Likelihood.md)
  object specifying the relevant likelihood/nuisance estimators for the
  meta-learner. used for estimating the parameters of the NPSEM.

- `...`:

  Additional arguments to pass to the initialization function.

------------------------------------------------------------------------

### Method `add_nuisance_estimator()`

#### Usage

    hte3_Task$add_nuisance_estimator(node, learner)

------------------------------------------------------------------------

### Method `get_nuisance_estimates()`

#### Usage

    hte3_Task$get_nuisance_estimates(
      nodes,
      hte3_task = NULL,
      fold_number = "validation"
    )

------------------------------------------------------------------------

### Method `next_in_chain()`

#### Usage

    hte3_Task$next_in_chain(
      covariates = NULL,
      outcome = NULL,
      id = NULL,
      weights = NULL,
      offset = NULL,
      time = NULL,
      folds = NULL,
      column_names = NULL,
      new_nodes = NULL,
      new_outcome_type = NULL,
      ...
    )

#### Arguments

- `...`:

  Additional arguments to pass to the initialization function.

------------------------------------------------------------------------

### Method `subset_task()`

#### Usage

    hte3_Task$subset_task(row_index, drop_folds = FALSE)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    hte3_Task$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
