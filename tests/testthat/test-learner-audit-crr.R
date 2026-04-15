test_that("CRR learners fit, predict, and preserve a monotone log-risk-ratio signal", {
  skip_if_runtime_unavailable()
  data <- make_binary_crr_data()
  task <- make_crr_task_for_test(data)
  base_learner <- make_sl3_regression_learner(stats::binomial())
  modifier_data <- data[, c("W1", "W2"), with = FALSE]

  ipw_fit <- Lrnr_crr_IPW$new(base_learner = base_learner)$train(task)
  t_fit <- Lrnr_crr_T$new(base_learner = base_learner, stratify_by_treatment = TRUE)$train(task)
  t_pooled_fit <- Lrnr_crr_T$new(base_learner = base_learner, stratify_by_treatment = FALSE)$train(task)

  for (fit in list(ipw_fit, t_fit, t_pooled_fit)) {
    predictions <- predict_hte3(fit, modifier_data)
    expect_length(predictions, nrow(data))
    expect_true(all(is.finite(predictions)))
    expect_monotone_signal(predictions, data$W1)
  }
})

test_that("CRR EP learner and portfolio CV work end-to-end", {
  skip_if_runtime_unavailable(c("glmnet", "Sieve"))
  data <- make_binary_crr_data(seed = 30)
  task <- make_crr_task_for_test(data)
  base_learner <- make_sl3_regression_learner(stats::binomial())
  modifier_data <- data[, c("W1", "W2"), with = FALSE]

  ep_fit <- Lrnr_crr_EP$new(base_learner = base_learner, sieve_num_basis = 6)$train(task)
  ep_predictions <- predict_hte3(ep_fit, modifier_data)
  expect_length(ep_predictions, nrow(data))
  expect_true(all(is.finite(ep_predictions)))
  expect_monotone_signal(ep_predictions, data$W1)

  cv_fit <- cross_validate_crr(
    list(
      Lrnr_crr_IPW$new(base_learner = base_learner),
      Lrnr_crr_T$new(base_learner = base_learner)
    ),
    task,
    cv_control = list(V = 2)
  )
  expect_true(inherits(cv_fit$lrnr_sl, "Lrnr_sl"))
  expect_true(length(cv_fit$coefficients) >= 1L)

  model <- fit_crr(
    task,
    method = c("ipw", "t", "ep"),
    base_learner = base_learner,
    cross_validate = TRUE,
    cv_control = list(V = 2)
  )
  portfolio_predictions <- predict(model, modifier_data)
  expect_s3_class(model, "hte3_model")
  expect_match(model$method, "^portfolio\\(")
  expect_length(portfolio_predictions, nrow(data))
})

test_that("CRR learners handle categorical treatment contrasts and outcome validation", {
  skip_if_runtime_unavailable(c("glmnet", "Sieve"))
  data <- make_categorical_crr_data()
  task <- make_crr_task_for_test(data, treatment_type = "categorical")
  base_learner <- make_sl3_regression_learner(stats::binomial())
  modifier_data <- data[, c("W1", "W2"), with = FALSE]

  fits <- list(
    ipw = Lrnr_crr_IPW$new(base_learner = base_learner, treatment_level = "high", control_level = "control")$train(task),
    t = Lrnr_crr_T$new(base_learner = base_learner, treatment_level = "high", control_level = "control")$train(task),
    ep = Lrnr_crr_EP$new(base_learner = base_learner, sieve_num_basis = 6, treatment_level = "high", control_level = "control")$train(task)
  )

  for (fit in fits) {
    predictions <- predict_hte3(fit, modifier_data)
    expect_length(predictions, nrow(data))
    expect_true(all(is.finite(predictions)))
  }

  negative_outcome_data <- data.table::copy(make_binary_crr_data(seed = 31))
  negative_outcome_data$Y[1] <- -1
  negative_task <- hte_task(
    data = negative_outcome_data,
    modifiers = c("W1", "W2"),
    confounders = c("W1", "W2", "W3"),
    treatment = "A",
    outcome = "Y",
    propensity_learner = make_sl3_regression_learner(stats::binomial()),
    outcome_learner = make_sl3_regression_learner(stats::binomial()),
    mean_learner = make_sl3_regression_learner(stats::binomial()),
    cross_fit = FALSE
  )
  expect_error(fit_crr(negative_task, method = "ipw", base_learner = base_learner, cross_validate = FALSE), "non-negative")
})

test_that("CRR selector handles single-candidate portfolios and IPW guards zero positive weight", {
  skip_if_runtime_unavailable()
  data <- make_binary_crr_data(seed = 40)
  task <- make_crr_task_for_test(data)
  base_learner <- make_sl3_regression_learner(stats::binomial())

  selector_fit <- cross_validate_crr(
    list(Lrnr_crr_IPW$new(base_learner = base_learner)),
    task,
    cv_control = list(V = 2)
  )
  expect_true(inherits(selector_fit$lrnr_sl, "Lrnr_sl"))

  zero_outcome <- data.table::copy(data)
  zero_outcome$Y <- 0
  zero_task <- make_crr_task_for_test(zero_outcome)
  expect_error(
    Lrnr_crr_IPW$new(base_learner = base_learner)$train(zero_task),
    "positive CRR weight"
  )
})
