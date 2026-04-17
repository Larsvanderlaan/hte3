skip_if_not_installed("sl3")
skip_if_not_installed("tmle3")

small_grf_params <- function() {
  list(num.trees = 40, num.threads = 1)
}

capture_warnings <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  list(value = value, warnings = warnings)
}

test_that("GRF helper factories return the expected learner types", {
  nuisance <- make_grf_nuisance_learners()
  expect_named(nuisance, c("propensity_learner", "outcome_learner", "mean_learner"))
  expect_true(inherits(nuisance$propensity_learner, "Lrnr_base"))

  cate_learners <- make_grf_cate_learners(method = c("ep", "r"), grf_params = small_grf_params(), tune = "none")
  expect_true(any(vapply(cate_learners, inherits, logical(1), "Lrnr_cate_EP")))
  expect_true(any(vapply(cate_learners, inherits, logical(1), "Lrnr_grf_causal_forest")))
  expect_identical(
    cate_learners[[which(vapply(cate_learners, inherits, logical(1), "Lrnr_cate_EP"))[1]]]$params$targeting_style,
    "r"
  )

  ep_style_learners <- make_grf_cate_learners(
    method = "ep",
    grf_params = small_grf_params(),
    tune = "none",
    ep_targeting_style = c("r", "dr", "r")
  )
  expect_length(ep_style_learners, 2L)
  expect_identical(vapply(ep_style_learners, function(learner) learner$params$targeting_style, character(1)), c("dr", "r"))

  crr_learners <- make_grf_crr_learners(method = c("ipw", "t", "ep"), grf_params = small_grf_params(), tune = "none")
  expect_true(any(vapply(crr_learners, inherits, logical(1), "Lrnr_crr_IPW")))
  expect_true(any(vapply(crr_learners, inherits, logical(1), "Lrnr_crr_T")))
  expect_true(any(vapply(crr_learners, inherits, logical(1), "Lrnr_crr_EP")))
})

test_that("GRF wrappers validate binary nuisance shapes before fitting", {
  data <- hte3_example_data(n = 50, seed = 901)

  expect_error(
    grf_cate(
      data = data,
      modifiers = c("W1", "W2"),
      confounders = c("W1", "W2", "W3"),
      treatment = "A",
      outcome = "Y",
      mu.hat = matrix(0, nrow(data), 3),
      cross_fit = FALSE
    ),
    "n x 2 matrix"
  )

  expect_error(
    grf_cate(
      data = data,
      modifiers = c("W1", "W2"),
      confounders = c("W1", "W2", "W3"),
      treatment = "A",
      outcome = "Y",
      pi.hat = rep(0.5, nrow(data) - 1L),
      cross_fit = FALSE
    ),
    "must have length"
  )
})

test_that("GRF quasibinomial forest predictions stay in (0, 1)", {
  skip_if_runtime_unavailable("grf")

  task <- sl3::sl3_Task$new(
    data = data.table::data.table(
      W1 = runif(60),
      W2 = runif(60),
      Y = pmin(0.99, pmax(0.01, stats::plogis(rnorm(60))))
    ),
    covariates = c("W1", "W2"),
    outcome = "Y",
    outcome_type = "quasibinomial"
  )

  learner <- Lrnr_grf_forest$new(tune = "none", num.trees = 40, num.threads = 1)
  fit <- learner$train(task)
  preds <- fit$predict(task)

  expect_true(all(preds > 0))
  expect_true(all(preds < 1))
})

test_that("public fitting wrappers suppress internal cv-fallback warnings", {
  data <- hte3_example_data(n = 80, seed = 900)
  task <- hte_task(
    data = data,
    modifiers = c("W1", "W2"),
    confounders = c("W1", "W2", "W3"),
    treatment = "A",
    outcome = "Y",
    propensity_learner = make_sl3_regression_learner(stats::binomial()),
    outcome_learner = make_sl3_regression_learner(stats::gaussian()),
    mean_learner = make_sl3_regression_learner(stats::gaussian()),
    cross_fit = FALSE
  )

  fit_result <- capture_warnings(
    fit_cate(
      task = task,
      method = "dr",
      base_learner = make_sl3_regression_learner(stats::gaussian()),
      cross_validate = FALSE
    )
  )

  expect_false(any(grepl("is not cv-aware: self\\$predict_fold reverts to self\\$predict", fit_result$warnings)))
})

test_that("grf_cate defaults to EP-vs-R selection and respects direct nuisances", {
  skip_if_runtime_unavailable("grf")

  data <- hte3_example_data(n = 80, seed = 902)
  mu.hat <- cbind(data$mu0, data$mu1)
  pi.hat <- data$pi1
  expected_m <- data$mu0 * (1 - pi.hat) + data$mu1 * pi.hat

  model <- grf_cate(
    data = data,
    modifiers = c("W1", "W2"),
    confounders = c("W1", "W2", "W3"),
    treatment = "A",
    outcome = "Y",
    mu.hat = mu.hat,
    pi.hat = pi.hat,
    cross_fit = FALSE,
    cv_control = list(V = 2),
    tune = "none",
    grf_params = small_grf_params()
  )

  preds <- predict(model, data.table::data.table(W1 = data$W1, W2 = data$W2))
  task <- model$training_task
  task_pi <- as.matrix(task$get_nuisance_estimates("pi"))
  task_m <- as.numeric(task$get_nuisance_estimates("m"))

  expect_identical(model$engine, "grf")
  expect_true(model$cross_validated)
  expect_true(model$selection_summary$selected_method %in% c("ep", "r"))
  expect_length(preds, nrow(data))
  expect_true(all(is.finite(preds)))
  expect_true(any(grepl("style=r", model$selection_summary$candidate_labels, fixed = TRUE)))
  expect_false(any(grepl("style=dr", model$selection_summary$candidate_labels, fixed = TRUE)))
  expect_equal(task_pi[, 2], pi.hat, tolerance = 1e-8)
  expect_equal(task_m, expected_m, tolerance = 1e-8)
})

test_that("grf public wrappers suppress internal cv-fallback warnings", {
  skip_if_runtime_unavailable("grf")

  data <- hte3_example_data(n = 80, seed = 9021)
  fit_result <- capture_warnings(
    grf_cate(
      data = data,
      modifiers = c("W1", "W2"),
      confounders = c("W1", "W2", "W3"),
      treatment = "A",
      outcome = "Y",
      method = "dr",
      cross_validate = FALSE,
      cross_fit = FALSE,
      tune = "none",
      grf_params = small_grf_params()
    )
  )

  expect_identical(fit_result$value$method, "dr")
  expect_false(any(grepl("is not cv-aware: self\\$predict_fold reverts to self\\$predict", fit_result$warnings)))
})

test_that("grf_cate supports explicit DR fits", {
  skip_if_runtime_unavailable("grf")

  data <- hte3_example_data(n = 80, seed = 903)

  model <- grf_cate(
    data = data,
    modifiers = c("W1", "W2"),
    confounders = c("W1", "W2", "W3"),
    treatment = "A",
    outcome = "Y",
    method = "dr",
    cross_validate = FALSE,
    cross_fit = FALSE,
    tune = "none",
    grf_params = small_grf_params()
  )

  preds <- predict(model, data.table::data.table(W1 = data$W1, W2 = data$W2))
  expect_identical(model$method, "dr")
  expect_true(all(is.finite(preds)))
})

test_that("grf_cate supports wrapper-level EP style selection", {
  skip_if_runtime_unavailable("grf")

  data <- hte3_example_data(n = 80, seed = 9031)

  model <- grf_cate(
    data = data,
    modifiers = c("W1", "W2"),
    confounders = c("W1", "W2", "W3"),
    treatment = "A",
    outcome = "Y",
    method = "ep",
    cross_fit = FALSE,
    cv_control = list(V = 2),
    ep_targeting_style = c("r", "dr"),
    tune = "none",
    grf_params = small_grf_params()
  )

  expect_true(model$cross_validated)
  expect_identical(model$selection_summary$selected_method, "ep")
  expect_true(model$selection_summary$selected_ep_targeting_style %in% c("dr", "r"))
  expect_true(any(grepl("style=dr", model$selection_summary$candidate_labels, fixed = TRUE)))
  expect_true(any(grepl("style=r", model$selection_summary$candidate_labels, fixed = TRUE)))
})

test_that("grf_crr supports IPW, T, and EP with bounded nuisance behavior", {
  skip_if_runtime_unavailable("grf")

  data <- hte3_example_data(n = 90, seed = 904)

  for (method in c("ipw", "t", "ep")) {
    model <- grf_crr(
      data = data,
      modifiers = c("W1", "W2"),
      confounders = c("W1", "W2", "W3"),
      treatment = "A",
      outcome = "Y_binary",
      method = method,
      cross_validate = FALSE,
      cross_fit = FALSE,
      tune = "none",
      grf_params = small_grf_params()
    )

    preds <- predict(model, data.table::data.table(W1 = data$W1, W2 = data$W2))
    expect_identical(model$method, method)
    expect_true(all(is.finite(preds)))
  }
})

test_that("GRF wrappers reject unsupported treatment types and CRR outcome violations", {
  skip_if_runtime_unavailable("grf")

  cate_data <- hte3_example_data(n = 70, seed = 905)
  cate_data[, A_cont := W1]

  expect_error(
    grf_cate(
      data = cate_data,
      modifiers = c("W1", "W2"),
      confounders = c("W1", "W2", "W3"),
      treatment = "A_cont",
      outcome = "Y",
      cross_fit = FALSE,
      tune = "none",
      grf_params = small_grf_params()
    ),
    "binary treatments"
  )

  expect_error(
    grf_cate(
      data = hte3_example_data(n = 70, seed = 9051),
      modifiers = c("W1", "W2"),
      confounders = c("W1", "W2", "W3"),
      treatment = "A",
      outcome = "Y",
      method = "ep",
      cross_validate = FALSE,
      ep_targeting_style = c("dr", "r"),
      cross_fit = FALSE,
      tune = "none",
      grf_params = small_grf_params()
    ),
    "must be TRUE when fitting a learner portfolio"
  )

  bad_crr <- hte3_example_data(n = 70, seed = 906)
  bad_crr[, Y_bad := Y_binary - 2]

  expect_error(
    grf_crr(
      data = bad_crr,
      modifiers = c("W1", "W2"),
      confounders = c("W1", "W2", "W3"),
      treatment = "A",
      outcome = "Y_bad",
      cross_fit = FALSE,
      tune = "none",
      grf_params = small_grf_params()
    ),
    "non-negative"
  )
})
