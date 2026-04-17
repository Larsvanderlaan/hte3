# XGBoost Learner with Early Stopping

Package-local `sl3` learner for `xgboost` models with an internal
validation split. The learner uses early stopping to choose the number
of boosting rounds, then retrains the final model on the full training
data with the selected iteration count.

## Format

An R6 class inheriting from
[`sl3::Lrnr_base`](https://tlverse.org/sl3/reference/Lrnr_base.html).

## Arguments

- nrounds:

  Maximum number of boosting rounds considered before early stopping.

- validation_fraction:

  Fraction of the training sample reserved for the internal validation
  split.

- early_stopping_rounds:

  Number of validation rounds without improvement before stopping.

- seed:

  Random seed used for the internal validation split.

- nthread:

  Number of threads used by `xgboost`.

- ...:

  Additional arguments passed to
  [`xgb.train`](https://rdrr.io/pkg/xgboost/man/xgb.train.html).

## Value

An [`sl3::Lrnr_base`](https://tlverse.org/sl3/reference/Lrnr_base.html)
learner that selects the number of boosting rounds by early stopping and
then refits on the full training sample.

## Details

This learner is used by
[`get_autoML`](https://larsvanderlaan.github.io/hte3/reference/get_autoML.md)
to avoid manually fixing `nrounds` in the default boosted-tree
candidates.

## Methods

### Public methods

- `Lrnr_xgboost_early_stopping$new()`

- `Lrnr_xgboost_early_stopping$clone()`

### Method `new()`

#### Usage

    Lrnr_xgboost_early_stopping$new(
      nrounds = 500,
      validation_fraction = 0.2,
      early_stopping_rounds = 25,
      seed = 1,
      nthread = 1,
      ...
    )

#### Arguments

- `nrounds`:

  Maximum number of boosting rounds considered before early stopping.

- `validation_fraction`:

  Fraction of the training sample reserved for the internal validation
  split.

- `early_stopping_rounds`:

  Number of validation rounds without improvement before stopping.

- `seed`:

  Random seed used for the internal validation split.

- `nthread`:

  Number of threads used by `xgboost`.

- `...`:

  Additional arguments passed to
  [`xgb.train`](https://rdrr.io/pkg/xgboost/man/xgb.train.html).

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    Lrnr_xgboost_early_stopping$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
