test_that("CATE DR, R, and stratified T learners preserve a monotone signal on binary treatment data", {
  skip_if_runtime_unavailable()
  data <- make_binary_cate_data()
  base_learner <- make_sl3_regression_learner(stats::gaussian())
  modifier_data <- data[, c("W1", "W2"), with = FALSE]
  learners <- list(
    Lrnr_cate_DR$new(base_learner = base_learner),
    Lrnr_cate_R$new(base_learner = base_learner),
    Lrnr_cate_T$new(base_learner = base_learner, stratify_by_treatment = TRUE)
  )

  for (learner in learners) {
    fit <- learner$train(make_cate_task_for_test(data))
    predictions <- predict_hte3(fit, modifier_data)
    expect_length(predictions, nrow(data))
    expect_true(all(is.finite(predictions)))
    expect_monotone_signal(predictions, data$W1)
  }
})

test_that("CATE pooled T learner returns a valid average contrast with additive base learners", {
  skip_if_runtime_unavailable()
  data <- make_binary_cate_data(seed = 5)
  base_learner <- make_sl3_regression_learner(stats::gaussian())
  fit <- Lrnr_cate_T$new(base_learner = base_learner, stratify_by_treatment = FALSE)$train(make_cate_task_for_test(data))
  predictions <- predict_hte3(fit, data[, c("W1", "W2"), with = FALSE])

  expect_length(predictions, nrow(data))
  expect_true(all(is.finite(predictions)))
  expect_positive_average_signal(predictions)
})

test_that("CATE pooled T learner recovers heterogeneity with an interaction-capable base learner", {
  skip_if_runtime_unavailable("earth")
  data <- make_binary_cate_data(seed = 6)
  base_learner <- make_sl3_interaction_learner()
  fit <- Lrnr_cate_T$new(base_learner = base_learner, stratify_by_treatment = FALSE)$train(make_cate_task_for_test(data))
  predictions <- predict_hte3(fit, data[, c("W1", "W2"), with = FALSE])

  expect_length(predictions, nrow(data))
  expect_true(all(is.finite(predictions)))
  expect_monotone_signal(predictions, data$W1)
})

test_that("CATE EP learner and portfolio CV work end-to-end", {
  skip_if_runtime_unavailable(c("glmnet", "Sieve"))
  data <- make_binary_cate_data(seed = 10)
  base_learner <- make_sl3_regression_learner(stats::gaussian())
  modifier_data <- data[, c("W1", "W2"), with = FALSE]

  ep_fit <- Lrnr_cate_EP$new(base_learner = base_learner, sieve_num_basis = 6)$train(make_cate_task_for_test(data))
  ep_predictions <- predict_hte3(ep_fit, modifier_data)
  expect_length(ep_predictions, nrow(data))
  expect_true(all(is.finite(ep_predictions)))
  expect_monotone_signal(ep_predictions, data$W1)

  cv_fit <- cross_validate_cate(
    list(
      Lrnr_cate_DR$new(base_learner = base_learner),
      Lrnr_cate_R$new(base_learner = base_learner)
    ),
    make_cate_task_for_test(data),
    cv_control = list(V = 2)
  )
  expect_true(inherits(cv_fit$lrnr_sl, "Lrnr_sl"))
  expect_true(length(cv_fit$coefficients) >= 1L)

  model <- fit_cate(
    make_cate_task_for_test(data),
    method = c("dr", "r", "ep"),
    base_learner = base_learner,
    cross_validate = TRUE,
    cv_control = list(V = 2)
  )
  portfolio_predictions <- predict(model, modifier_data)
  expect_s3_class(model, "hte3_model")
  expect_match(model$method, "^portfolio\\(")
  expect_length(portfolio_predictions, nrow(data))
})

test_that("CATE learners handle categorical treatment contrasts consistently", {
  skip_if_runtime_unavailable(c("glmnet", "Sieve"))
  data <- make_categorical_cate_data()
  base_learner <- make_sl3_regression_learner(stats::gaussian())
  modifier_data <- data[, c("W1", "W2"), with = FALSE]
  learners <- list(
    Lrnr_cate_DR$new(base_learner = base_learner, treatment_level = "high", control_level = "control"),
    Lrnr_cate_T$new(base_learner = base_learner, treatment_level = "high", control_level = "control"),
    Lrnr_cate_EP$new(base_learner = base_learner, sieve_num_basis = 6, treatment_level = "high", control_level = "control")
  )

  for (learner in learners) {
    fit <- learner$train(make_cate_task_for_test(data, treatment_type = "categorical"))
    predictions <- predict_hte3(fit, modifier_data)
    expect_length(predictions, nrow(data))
    expect_true(all(is.finite(predictions)))
  }
})

test_that("CATE learner validation catches unsupported treatment codings and selector portfolios", {
  skip_if_runtime_unavailable()
  data <- make_binary_cate_data(seed = 20)
  data[, A_label := factor(ifelse(A == 1, "treated", "control"), levels = c("control", "treated"))]
  task <- make_cate_task_for_test(data, treatment = "A_label")
  base_learner <- make_sl3_regression_learner(stats::gaussian())

  expect_error(
    Lrnr_cate_R$new(base_learner = base_learner)$train(task),
    "numeric treatment coding"
  )

  selector_fit <- cross_validate_cate(
    list(Lrnr_cate_DR$new(base_learner = base_learner)),
    make_cate_task_for_test(make_binary_cate_data(seed = 21)),
    cv_control = list(V = 2)
  )
  expect_true(inherits(selector_fit$lrnr_sl, "Lrnr_sl"))
  expect_true(length(selector_fit$coefficients) >= 1L)
})
