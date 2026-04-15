skip_if_not_installed("sl3")
skip_if_not_installed("tmle3")

test_that("legacy task constructor still works on a simple example", {
  skip_if_not("Lrnr_mean" %in% getNamespaceExports("sl3"), "sl3::Lrnr_mean is required for this test.")

  data <- hte3_example_data(n = 60, seed = 21)
  learner <- sl3::Lrnr_mean$new()

  task <- make_hte3_Task_tx(
    data = data,
    modifiers = c("W1", "W2"),
    confounders = c("W1", "W2", "W3"),
    treatment = "A",
    outcome = "Y",
    learner_pi = learner,
    learner_mu = learner,
    learner_m = learner,
    cross_fit_and_cv = FALSE
  )

  fitted <- Lrnr_cate_DR$new(base_learner = learner)$train(task)

  expect_true(inherits(task, "hte3_Task"))
  expect_length(predict_hte3(fitted, data), nrow(data))
})

test_that("predict_hte3 checks for required modifier columns", {
  skip_if_not("Lrnr_mean" %in% getNamespaceExports("sl3"), "sl3::Lrnr_mean is required for this test.")

  data <- hte3_example_data(n = 40, seed = 22)
  learner <- sl3::Lrnr_mean$new()

  task <- make_hte3_Task_tx(
    data = data,
    modifiers = c("W1", "W2"),
    confounders = c("W1", "W2", "W3"),
    treatment = "A",
    outcome = "Y",
    learner_pi = learner,
    learner_mu = learner,
    learner_m = learner,
    cross_fit_and_cv = FALSE
  )

  fitted <- Lrnr_cate_DR$new(base_learner = learner)$train(task)

  expect_error(predict_hte3(fitted, data.table::data.table(W1 = data$W1)))
})
