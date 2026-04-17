# Predict Heterogeneous Treatment Effects using hte3 Learners

Generates treatment effect predictions using a trained hte3 learner.

## Usage

``` r
predict_hte3(hte3_learner, new_data)
```

## Arguments

- hte3_learner:

  A trained hte3 learner (`Lrnr_sl` object) used to make predictions.

- new_data:

  A data frame or data.table containing new data of effect `modifiers`
  for prediction.

## Value

A numeric vector of predicted heterogeneous treatment effects for the
provided new data.
