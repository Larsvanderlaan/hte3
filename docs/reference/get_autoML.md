# Get Automated Machine Learning (AutoML) Learner

This function returns the default `sl3` learner stack used by the
high-level `hte3` wrappers. The stack always includes the core learners
`Lrnr_glmnet` and `Lrnr_gam`. It adds the following optional learners
only when their supporting runtime packages are available:

- `Lrnr_earth`: Multivariate Adaptive Regression Spline Model.

- `Lrnr_ranger`: Random Forest Model.

- `Lrnr_xgboost_early_stopping`: XGBoost models with early stopping
  across a small set of depth and shrinkage configurations.

If optional packages are unavailable, `get_autoML()` emits a warning
once and returns the available safe subset rather than failing at
runtime.

## Usage

``` r
get_autoML()
```

## Value

A [`sl3::Stack`](https://tlverse.org/sl3/reference/Stack.html) object
containing the available safe subset of the predefined default learner
library.
