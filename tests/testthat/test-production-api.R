skip_if_not_installed("sl3")
skip_if_not_installed("tmle3")

test_that("production CATE workflow fits and predicts", {
  skip_if_not("Lrnr_mean" %in% getNamespaceExports("sl3"), "sl3::Lrnr_mean is required for this test.")

  data <- hte3_example_data(n = 80, seed = 11)
  learner <- sl3::Lrnr_mean$new()

  task <- hte_task(
    data = data,
    modifiers = c("W1", "W2", "W3"),
    confounders = c("W1", "W2", "W3"),
    treatment = "A",
    outcome = "Y",
    propensity_learner = learner,
    outcome_learner = learner,
    mean_learner = learner,
    cross_fit = FALSE
  )

  model <- fit_cate(task, method = "dr", base_learner = learner, cross_validate = FALSE)
  preds <- predict(model, data)

  expect_s3_class(model, "hte3_model")
  expect_length(preds, nrow(data))
  expect_s3_class(summary(model), "summary.hte3_model")
})

test_that("production models predict on modifier-only new data", {
  skip_if_not("Lrnr_mean" %in% getNamespaceExports("sl3"), "sl3::Lrnr_mean is required for this test.")

  data <- hte3_example_data(n = 70, seed = 13)
  learner <- sl3::Lrnr_mean$new()

  task <- hte_task(
    data = data,
    modifiers = c("W1", "W2", "W3"),
    confounders = c("W1", "W2", "W3"),
    treatment = "A",
    outcome = "Y",
    propensity_learner = learner,
    outcome_learner = learner,
    mean_learner = learner,
    cross_fit = FALSE
  )

  model <- fit_cate(task, method = "t", base_learner = learner, cross_validate = FALSE)
  preds <- predict(model, data.table::data.table(W1 = data$W1, W2 = data$W2, W3 = data$W3))

  expect_length(preds, nrow(data))
  expect_true(all(is.finite(preds)))
})

test_that("production CRR workflow fits and predicts", {
  skip_if_not("Lrnr_mean" %in% getNamespaceExports("sl3"), "sl3::Lrnr_mean is required for this test.")

  data <- hte3_example_data(n = 80, seed = 12)
  learner <- sl3::Lrnr_mean$new()

  task <- hte_task(
    data = data,
    modifiers = c("W1", "W2", "W3"),
    confounders = c("W1", "W2", "W3"),
    treatment = "A",
    outcome = "Y_binary",
    propensity_learner = learner,
    outcome_learner = learner,
    mean_learner = learner,
    cross_fit = FALSE
  )

  model <- fit_crr(task, method = "ipw", base_learner = learner, cross_validate = FALSE)
  preds <- predict(model, data)

  expect_s3_class(model, "hte3_model")
  expect_length(preds, nrow(data))
})

test_that("hte_task validates weights, folds, and supplied nuisance estimates", {
  skip_if_not("Lrnr_mean" %in% getNamespaceExports("sl3"), "sl3::Lrnr_mean is required for this test.")

  data <- hte3_example_data(n = 50, seed = 14)
  learner <- sl3::Lrnr_mean$new()
  mu_hat <- cbind("0" = data$mu0, "1" = data$mu1)
  bad_pi <- cbind("0" = rep(0.7, nrow(data)), "1" = rep(0.4, nrow(data)))

  expect_error(
    hte_task(
      data = data,
      modifiers = c("W1", "W2"),
      confounders = c("W1", "W2", "W3"),
      treatment = "A",
      outcome = "Y",
      propensity = bad_pi,
      outcome_regression = mu_hat,
      cross_fit = FALSE
    ),
    "sum to 1"
  )

  weighted_data <- data.table::copy(data)
  weighted_data[, wt := -1]
  expect_error(
    hte_task(
      data = weighted_data,
      modifiers = c("W1", "W2"),
      confounders = c("W1", "W2", "W3"),
      treatment = "A",
      outcome = "Y",
      weights = "wt",
      propensity_learner = learner,
      outcome_learner = learner,
      mean_learner = learner,
      cross_fit = FALSE
    ),
    "weights"
  )

  expect_error(
    hte_task(
      data = data,
      modifiers = c("W1", "W2"),
      confounders = c("W1", "W2", "W3"),
      treatment = "A",
      outcome = "Y",
      propensity_learner = learner,
      outcome_learner = learner,
      mean_learner = learner,
      cross_fit = FALSE,
      folds = 1
    ),
    "greater than or equal to 2"
  )

  factor_weighted_data <- data.table::copy(data)
  factor_weighted_data[, wt_factor := factor(rep(c("10", "20"), length.out = .N))]
  factor_task <- hte_task(
    data = factor_weighted_data,
    modifiers = c("W1", "W2"),
    confounders = c("W1", "W2", "W3"),
    treatment = "A",
    outcome = "Y",
    weights = "wt_factor",
    propensity_learner = learner,
    outcome_learner = learner,
    mean_learner = learner,
    cross_fit = FALSE
  )

  expect_identical(
    factor_task$data[["wt_factor"]],
    rep(c(10, 20), length.out = nrow(factor_weighted_data))
  )
})

test_that("get_autoML resolves a safe subset of the default stack", {
  resolver <- getFromNamespace("resolve_automl_candidate_specs", "hte3")
  specs <- getFromNamespace("automl_spec_table", "hte3")()

  resolved <- resolver(
    specs = specs,
    sl3_exports = c("Lrnr_glmnet", "Lrnr_gam", "Lrnr_earth", "Lrnr_ranger"),
    package_available = function(pkg) !(pkg %in% c("ranger", "xgboost"))
  )

  expect_true(length(resolved) >= 2L)
  expect_true(all(vapply(resolved, function(spec) spec$learner_class, character(1)) %in% c("Lrnr_glmnet", "Lrnr_gam", "Lrnr_earth")))

  omitted <- attr(resolved, "omitted_specs")
  expect_true(length(omitted) >= 2L)
  expect_true(any(vapply(omitted, function(spec) identical(spec$required_package, "ranger"), logical(1))))
  expect_true(any(vapply(omitted, function(spec) identical(spec$required_package, "xgboost"), logical(1))))
})

test_that("fit_crr rejects continuous-treatment tasks", {
  skip_if_not("Lrnr_mean" %in% getNamespaceExports("sl3"), "sl3::Lrnr_mean is required for this test.")

  data <- hte3_example_data(n = 60, seed = 15)
  data[, A_cont := W1]
  learner <- sl3::Lrnr_mean$new()

  task <- hte_task(
    data = data,
    modifiers = c("W1", "W2"),
    confounders = c("W1", "W2", "W3"),
    treatment = "A_cont",
    outcome = "Y_binary",
    treatment_type = "continuous",
    propensity_learner = learner,
    outcome_learner = learner,
    mean_learner = learner,
    cross_fit = FALSE
  )

  expect_error(
    fit_crr(task, method = "ipw", base_learner = learner, cross_validate = FALSE),
    "binary or categorical treatments"
  )
})
