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
