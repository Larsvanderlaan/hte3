test_that("CRR IPW and two-stage stratified T learners preserve a monotone log-risk-ratio signal", {
  skip_if_runtime_unavailable()
  data <- make_binary_crr_data(n = 1200)
  base_learner <- make_sl3_regression_learner(stats::binomial())
  modifier_data <- data[, c("W1", "W2"), with = FALSE]
  pi1_true <- stats::plogis(0.1 + 0.4 * data$W2 - 0.25 * data$W3)
  propensity_truth <- cbind("0" = 1 - pi1_true, "1" = pi1_true)
  learners <- list(
    Lrnr_crr_IPW$new(base_learner = base_learner),
    Lrnr_crr_T$new(base_learner = base_learner, stratify_by_treatment = TRUE)
  )
  tasks <- list(
    make_crr_task_for_test(data, propensity = propensity_truth),
    make_crr_task_for_test(data)
  )

  for (index in seq_along(learners)) {
    fit <- learners[[index]]$train(tasks[[index]])
    predictions <- predict_hte3(fit, modifier_data)
    expect_length(predictions, nrow(data))
    expect_true(all(is.finite(predictions)))
    expect_monotone_signal(predictions, data$W1)
  }
})

test_that("CRR T learner only skips the second stage when modifiers equal confounders", {
  skip_if_runtime_unavailable()
  data <- make_binary_crr_data(seed = 7)
  base_learner <- make_sl3_regression_learner(stats::binomial())
  reduced_task <- make_crr_task_for_test(data)
  full_task <- make_crr_task_for_test(
    data,
    modifiers = c("W1", "W2", "W3"),
    confounders = c("W1", "W2", "W3")
  )

  expect_error(
    Lrnr_crr_T$new(base_learner = base_learner, second_stage_regression = FALSE)$train(reduced_task),
    "only supported when `modifiers` and `confounders` are the same"
  )

  fit <- Lrnr_crr_T$new(base_learner = base_learner, second_stage_regression = FALSE)$train(full_task)
  predictions <- predict_hte3(fit, data[, c("W1", "W2", "W3"), with = FALSE])
  expect_length(predictions, nrow(data))
  expect_true(all(is.finite(predictions)))
})

test_that("CRR pooled T learner returns a valid positive average log-risk-ratio contrast", {
  skip_if_runtime_unavailable()
  data <- make_binary_crr_data(seed = 8)
  base_learner <- make_sl3_regression_learner(stats::binomial())
  fit <- Lrnr_crr_T$new(base_learner = base_learner, stratify_by_treatment = FALSE)$train(make_crr_task_for_test(data))
  predictions <- predict_hte3(fit, data[, c("W1", "W2"), with = FALSE])

  expect_length(predictions, nrow(data))
  expect_true(all(is.finite(predictions)))
  expect_positive_average_signal(predictions)
})

test_that("CRR EP learner and portfolio CV work end-to-end", {
  skip_if_runtime_unavailable(c("glmnet", "Sieve"))
  data <- make_binary_crr_data(seed = 30)
  base_learner <- make_sl3_regression_learner(stats::binomial())
  modifier_data <- data[, c("W1", "W2"), with = FALSE]

  ep_fit <- Lrnr_crr_EP$new(base_learner = base_learner, sieve_num_basis = 6)$train(make_crr_task_for_test(data))
  ep_predictions <- predict_hte3(ep_fit, modifier_data)
  expect_length(ep_predictions, nrow(data))
  expect_true(all(is.finite(ep_predictions)))
  expect_monotone_signal(ep_predictions, data$W1)

  cv_fit <- cross_validate_crr(
    list(
      Lrnr_crr_IPW$new(base_learner = base_learner),
      Lrnr_crr_T$new(base_learner = base_learner)
    ),
    make_crr_task_for_test(data),
    cv_control = list(V = 2)
  )
  expect_true(inherits(cv_fit$lrnr_sl, "Lrnr_sl"))
  expect_true(length(cv_fit$coefficients) >= 1L)

  model <- fit_crr(
    make_crr_task_for_test(data),
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
  base_learner <- make_sl3_regression_learner(stats::binomial())
  modifier_data <- data[, c("W1", "W2"), with = FALSE]
  learners <- list(
    Lrnr_crr_IPW$new(base_learner = base_learner, treatment_level = "high", control_level = "control"),
    Lrnr_crr_T$new(base_learner = base_learner, treatment_level = "high", control_level = "control"),
    Lrnr_crr_EP$new(base_learner = base_learner, sieve_num_basis = 6, treatment_level = "high", control_level = "control")
  )

  for (learner in learners) {
    fit <- learner$train(make_crr_task_for_test(data, treatment_type = "categorical"))
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
