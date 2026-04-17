# Truncate and Calibrate Propensity Scores

The truncation_method parameter controls how the propensity scores are
calibrated and truncated. - "isotonic": Performs isotonic calibration,
providing both data-adaptive truncation and calibration for the
propensity scores. - "adaptive": Adapts the truncation level using a
loss function for the inverse propensity score. - "deterministic":
Bounds the estimates away from 0 and 1 using the threshold
25/sqrt(n)/log(n). - "none": Bounds the estimates away from 0 and 1
using the threshold 1/n.

## Usage

``` r
truncate_propensity(
  pi.hat,
  A,
  treatment_level = max(A),
  truncation_method = c("isotonic", "adaptive", "deterministic", "none")
)
```

## Arguments

- pi.hat:

  A numeric vector containing estimates of the propensity score for
  `treatment_level`.

- A:

  A numeric vector of treatment values.

- treatment_level:

  A numeric value indicating the treatment level to calibrate.

- truncation_method:

  A character string indicating the truncation method: "isotonic",
  "adaptive", "deterministic", or "none".

## Value

A numeric vector of truncated and calibrated propensity scores.

## Details

This function truncates and calibrates propensity scores to ensure
bounded values and improve their reliability.
