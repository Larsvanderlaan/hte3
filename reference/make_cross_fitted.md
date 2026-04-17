# Create Cross-Fitted Learner

This function takes a learner and returns a cross-fitted version of it.
Cross-fitting involves fitting the learner to different subsets of the
data while using the complementary subsets for validation, to provide a
more robust estimate of model performance.

## Usage

``` r
make_cross_fitted(
  learner,
  calibrate = FALSE,
  cross_validate = inherits(learner, "Stack")
)
```

## Arguments

- learner:

  The learner to be cross-fitted. This can be a single learner or a
  stacked learner.

- calibrate:

  Currently not used. Logical indicating whether to calibrate the
  learner (default is `FALSE`).

- cross_validate:

  Logical indicating whether to perform cross-validation (default is
  `TRUE` if `learner` is a stacked learner).

## Value

A cross-fitted version of the input learner.
