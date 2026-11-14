test_that("nuisance predictions avoid sl3's default fold fallback", {
  skip_if_runtime_unavailable()
  data <- hte3_example_data(n = 80, seed = 730)
  task <- make_hte3_Task_tx(
    data, modifiers = c("W1", "W2", "W3"),
    confounders = c("W1", "W2", "W3"), treatment = "A", outcome = "Y",
    learner_pi = sl3::Lrnr_mean$new(), learner_mu = sl3::Lrnr_mean$new(),
    learner_m = sl3::Lrnr_mean$new(), cross_fit_and_cv = FALSE
  )

  for (node in c("pi", "mu", "m")) {
    factor <- task$likelihood$factor_list[[node]]
    expect_no_warning(actual <- factor$get_mean(task, "validation"))
    expect_equal(actual, factor$get_mean(task, "full"))
  }

  expect_no_warning(
    fit <- Lrnr_cate_DR$new(base_learner = sl3::Lrnr_mean$new())$train(task)
  )
  expect_true(all(is.finite(predict_hte3(fit, data))))
})

test_that("nuisance factors preserve cross-fitted and individual-fold predictions", {
  skip_if_runtime_unavailable()
  set.seed(731)
  regression_task <- sl3::sl3_Task$new(
    data = data.frame(W = 1:24, Y = (1:24)^2),
    covariates = "W", outcome = "Y", folds = 3
  )
  learner <- sl3::Lrnr_cv$new(sl3::Lrnr_mean$new(), full_fit = TRUE)$train(regression_task)
  factor <- LF_fit_hte3$new("m", learner, type = "mean")
  task <- list(
    get_regression_task = function(...) regression_task,
    unscale = function(preds, ...) preds
  )

  for (fold in list("validation", "full", 1L)) {
    expect_no_warning(actual <- factor$get_mean(task, fold))
    expect_equal(actual, learner$predict_fold(regression_task, fold))
  }
  expect_false(isTRUE(all.equal(
    factor$get_mean(task, "validation"), factor$get_mean(task, "full")
  )))
})

test_that("explicit continuous pseudo-outcomes do not trigger type-detection warnings", {
  skip_if_runtime_unavailable()
  data <- hte3_example_data(n = 80, seed = 732)
  task <- hte_task(
    data, modifiers = c("W1", "W2", "W3"),
    confounders = c("W1", "W2", "W3"), treatment = "A", outcome = "Y",
    propensity_learner = sl3::Lrnr_mean$new(),
    outcome_learner = sl3::Lrnr_mean$new(),
    mean_learner = sl3::Lrnr_mean$new(), cross_fit = FALSE
  )

  expect_no_warning(projection <- make_tlearner_projection_task(
    task, rep(c(0, 1), length.out = nrow(data))
  ))
  expect_identical(projection$outcome_type$type, "continuous")
})

test_that("genuine warnings from nuisance learners still reach the caller", {
  skip_if_runtime_unavailable()
  WarningLearner <- R6::R6Class(
    "Lrnr_warning_probe", inherit = sl3::Lrnr_mean,
    private = list(.predict = function(task) {
      warning("Nuisance model did not converge", call. = FALSE)
      super$.predict(task)
    })
  )
  regression_task <- sl3::sl3_Task$new(
    data = data.frame(W = 1:24, Y = (1:24)^2),
    covariates = "W", outcome = "Y"
  )
  learner <- WarningLearner$new()$train(regression_task)
  factor <- LF_fit_hte3$new("m", learner, type = "mean")
  task <- list(
    get_regression_task = function(...) regression_task,
    unscale = function(preds, ...) preds
  )

  expect_warning(factor$get_mean(task, "validation"), "Nuisance model did not converge")
})

test_that("forest and earth fits do not receive unused family arguments", {
  skip_if_runtime_unavailable()
  data <- hte3_example_data(n = 80, seed = 733)
  task <- hte_task(
    data, modifiers = c("W1", "W2", "W3"),
    confounders = c("W1", "W2", "W3"), treatment = "A", outcome = "Y",
    propensity_learner = sl3::Lrnr_mean$new(),
    outcome_learner = sl3::Lrnr_mean$new(),
    mean_learner = sl3::Lrnr_mean$new(), cross_fit = FALSE
  )
  learners <- list()
  if (requireNamespace("ranger", quietly = TRUE)) {
    learners$ranger <- sl3::Lrnr_ranger$new(num.trees = 20)
  }
  if (requireNamespace("earth", quietly = TRUE)) {
    learners$earth <- sl3::Lrnr_earth$new()
  }
  skip_if(length(learners) == 0L, "No optional forest or earth runtime available")

  for (learner in learners) {
    expect_no_message(fit <- fit_cate(
      task, method = "dr", base_learner = learner, cross_validate = FALSE
    ))
    expect_true(all(is.finite(predict(fit, data))))
  }
})
