# Advanced sl3 Integration

Read this article if you want the low-level API: explicit task
construction, manual learner choice, or direct calls to
[`cross_validate_cate()`](https://larsvanderlaan.github.io/hte3/reference/cross_validate_cate.md)
and
[`cross_validate_crr()`](https://larsvanderlaan.github.io/hte3/reference/cross_validate_crr.md).

## Build the Task Explicitly

Use
[`make_hte3_Task_tx()`](https://larsvanderlaan.github.io/hte3/reference/make_hte3_Task_tx.md)
when you want to supply the nuisance learners yourself or work directly
with the lower-level learner objects.

``` r
library(hte3)
library(sl3)

data <- hte3_example_data(n = 100, seed = 4)

task <- make_hte3_Task_tx(
  data = data,
  modifiers = c("W1", "W2"),
  confounders = c("W1", "W2", "W3"),
  treatment = "A",
  outcome = "Y",
  learner_pi = Lrnr_mean$new(),
  learner_mu = Lrnr_mean$new(),
  learner_m = Lrnr_mean$new(),
  cross_fit_and_cv = FALSE
)
```

## Train a Learner Directly

``` r
learner <- Lrnr_cate_DR$new(base_learner = Lrnr_mean$new())
trained <- learner$train(task)

head(predict_hte3(trained, data))
#> [1] 1.553732 1.553732 1.553732 1.553732 1.553732 1.553732
```

This is the right level when you want the exact learner object rather
than the wrapper returned by
[`fit_cate()`](https://larsvanderlaan.github.io/hte3/reference/fit_cate.md)
or
[`fit_crr()`](https://larsvanderlaan.github.io/hte3/reference/fit_crr.md).

## Low-Level Cross-Validation

You can cross-validate low-level learner portfolios directly.

``` r
cate_cv <- cross_validate_cate(
  list(
    Lrnr_cate_DR$new(base_learner = Lrnr_mean$new()),
    Lrnr_cate_R$new(base_learner = Lrnr_mean$new()),
    Lrnr_cate_EP$new(base_learner = Lrnr_mean$new(), sieve_num_basis = 4)
  ),
  task,
  cv_control = list(V = 2)
)

cate_cv$coefficients
#>                                   dr                                    r 
#>                                    0                                    0 
#> ep[style=dr, basis=4, interaction=3] 
#>                                    1
head(cate_cv$lrnr_sl$predict(task))
#> [1] 1.547718 1.547718 1.547718 1.547718 1.547718 1.547718
```

The analogous CRR workflow is:

``` r
crr_task <- make_hte3_Task_tx(
  data = data,
  modifiers = c("W1", "W2", "W3"),
  confounders = c("W1", "W2", "W3"),
  treatment = "A",
  outcome = "Y_binary",
  learner_pi = Lrnr_mean$new(),
  learner_mu = Lrnr_mean$new(),
  learner_m = Lrnr_mean$new(),
  cross_fit_and_cv = FALSE
)

crr_cv <- suppressWarnings(
  cross_validate_crr(
    list(
      Lrnr_crr_IPW$new(base_learner = Lrnr_mean$new()),
      Lrnr_crr_T$new(base_learner = Lrnr_mean$new()),
      Lrnr_crr_EP$new(base_learner = Lrnr_mean$new(), sieve_num_basis = 4)
    ),
    crr_task,
    cv_control = list(V = 2)
  )
)

crr_cv$coefficients
#>                        ipw                          t 
#>                          0                          1 
#> ep[basis=4, interaction=3] 
#>                          0
head(crr_cv$lrnr_sl$predict(crr_task))
#> [1] 0.1634243 0.1634243 0.1634243 0.1634243 0.1634243 0.1634243
```

## Important Constraints

- `cross_fit_and_cv` in
  [`make_hte3_Task_tx()`](https://larsvanderlaan.github.io/hte3/reference/make_hte3_Task_tx.md)
  controls nuisance estimation at the task stage.
- If `modifiers = V` and `confounders = W` with `V` a strict subset of
  `W`, the reduced-modifier target is
  `E[Y(1) - Y(0) | V] = E[tau(W) | V]`.
- In that reduced-modifier setting, the current low-level R-learner
  targets an overlap-weighted projection rather than the unweighted
  `V`-conditional CATE.
- For continuous-treatment CATE, the low-level `Lrnr_cate_R` path is the
  one currently implemented.

For EP-specific sieve construction and tuning, see [EP-Learner Sieve
Tuning](https://larsvanderlaan.github.io/hte3/articles/ep-learner.md).
If you are building a new learner, see [Extending hte3 with Custom
Learners](https://larsvanderlaan.github.io/hte3/articles/extending-hte3.md).
