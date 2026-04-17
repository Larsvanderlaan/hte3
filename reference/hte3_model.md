# S3 Methods for Production hte3 Models

These methods provide the default prediction and summary interface for
the production `hte3_model` wrapper class.

## Usage

``` r
# S3 method for class 'hte3_model'
predict(object, new_data = NULL, ...)

# S3 method for class 'hte3_model'
summary(object, ...)
```

## Arguments

- object:

  An object created by
  [`fit_cate()`](https://larsvanderlaan.github.io/hte3/reference/fit_cate.md)
  or
  [`fit_crr()`](https://larsvanderlaan.github.io/hte3/reference/fit_crr.md).

- new_data:

  Optional prediction data. If omitted, predictions are returned on the
  training task.

- ...:

  Unused additional arguments.

## Value

[`predict()`](https://rdrr.io/r/stats/predict.html) returns model
predictions. [`summary()`](https://rdrr.io/r/base/summary.html) returns
a compact summary object with key training metadata, including
wrapper-level cross-validation selection metadata when available, such
as EP basis-size and EP targeting-style fields.
