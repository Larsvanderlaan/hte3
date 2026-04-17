# Create an hte3 Task with a Production-Oriented Interface

`hte_task()` is the recommended high-level constructor for new users. It
wraps
[`make_hte3_Task_tx()`](https://larsvanderlaan.github.io/hte3/reference/make_hte3_Task_tx.md)
with validated, production-oriented argument names.

## Usage

``` r
hte_task(
  data,
  modifiers,
  confounders = modifiers,
  treatment,
  outcome,
  id = NULL,
  weights = NULL,
  treatment_type = c("default", "binomial", "categorical", "continuous"),
  propensity = NULL,
  outcome_regression = NULL,
  outcome_mean = NULL,
  propensity_learner = get_autoML(),
  outcome_learner = get_autoML(),
  mean_learner = NULL,
  cross_fit = TRUE,
  folds = 10,
  for_prediction = FALSE,
  ...
)
```

## Arguments

- data:

  A data frame or `data.table` containing treatment, outcome, modifiers,
  and confounders.

- modifiers:

  Character vector of effect-modifier column names. These define the
  target modifier set `V`.

- confounders:

  Character vector of adjustment-variable column names used for nuisance
  estimation. These define the adjustment set `W`.

- treatment:

  Treatment column name.

- outcome:

  Outcome column name.

- id:

  Optional ID column name.

- weights:

  Optional weight column name.

- treatment_type:

  Treatment scale.

- propensity:

  Optional matrix of user-supplied propensity estimates.

- outcome_regression:

  Optional matrix of user-supplied outcome-regression estimates.

- outcome_mean:

  Optional vector of user-supplied marginal outcome estimates.

- propensity_learner:

  Learner used for propensity estimation when `propensity` is not
  supplied.

- outcome_learner:

  Learner used for outcome-regression estimation when
  `outcome_regression` is not supplied.

- mean_learner:

  Learner used for `m(W)` estimation when `outcome_mean` is not
  supplied.

- cross_fit:

  Whether to cross-fit nuisance learners.

- folds:

  Number of folds for nuisance estimation.

- for_prediction:

  Whether to construct a prediction-only task with no nuisance fitting.

- ...:

  Additional arguments passed through to
  [`make_hte3_Task_tx()`](https://larsvanderlaan.github.io/hte3/reference/make_hte3_Task_tx.md).

## Value

An `hte3_Task` object.

## Details

When `modifiers = V` and `confounders = W`, the natural CATE target is
`E[Y(1) - Y(0) | V] = E[tau(W) | V]`. In the supported
binary/categorical-treatment setting, DR-, EP-, and the default
two-stage T-learner align with that modifier-conditional target. The
current R-learner does not generally target that object when `V` is a
strict subset of `W`; it warns at fit time.

## Examples

``` r
if (FALSE) { # \dontrun{
library(sl3)

data <- hte3_example_data(n = 80, seed = 1)

task <- hte_task(
  data = data,
  modifiers = c("W1", "W2"),
  confounders = c("W1", "W2", "W3"),
  treatment = "A",
  outcome = "Y",
  propensity_learner = Lrnr_mean$new(),
  outcome_learner = Lrnr_mean$new(),
  mean_learner = Lrnr_mean$new(),
  cross_fit = FALSE
)

task
} # }
```
