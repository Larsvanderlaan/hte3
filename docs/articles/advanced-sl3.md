# Using sl3 directly

Construct an `hte3_Task`, train individual `sl3` learners, and compare
them with
[`cross_validate_cate()`](https://larsvanderlaan.github.io/hte3/reference/cross_validate_cate.md)
or
[`cross_validate_crr()`](https://larsvanderlaan.github.io/hte3/reference/cross_validate_crr.md).

## Create a task

Pass your nuisance learners to
[`make_hte3_Task_tx()`](https://larsvanderlaan.github.io/hte3/reference/make_hte3_Task_tx.md):

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

## Train a learner

``` r
learner <- Lrnr_cate_DR$new(base_learner = Lrnr_mean$new())
trained <- learner$train(task)

head(predict_hte3(trained, data))
#> [1] 1.553732 1.553732 1.553732 1.553732 1.553732 1.553732
```

`trained` is an `sl3` learner object. Use
[`predict_hte3()`](https://larsvanderlaan.github.io/hte3/reference/predict_hte3.md)
to predict effects on new data.

## Compare learners

Pass a list of learners to compare by cross-validation:

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

For conditional risk ratios:

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

## Targets and treatment types

- `cross_fit_and_cv` in
  [`make_hte3_Task_tx()`](https://larsvanderlaan.github.io/hte3/reference/make_hte3_Task_tx.md)
  controls nuisance estimation at the task stage.
- If `modifiers = V` and `confounders = W` with `V` a strict subset of
  `W`, the reduced-modifier target is
  `E[Y(1) - Y(0) | V] = E[tau(W) | V]`.
- In that reduced-modifier setting, the current low-level R-learner
  targets an overlap-weighted projection rather than the unweighted
  `V`-conditional CATE.
- For continuous-treatment CATE, use `Lrnr_cate_R`.

For EP-specific sieve construction and tuning, see
[EP-learner](https://larsvanderlaan.github.io/hte3/articles/ep-learner.md).
If you are building a new learner, see [Extending hte3 with Custom
Learners](https://larsvanderlaan.github.io/hte3/articles/extending-hte3.md).
