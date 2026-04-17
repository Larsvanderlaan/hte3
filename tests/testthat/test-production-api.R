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

test_that("EP wrapper CV exposes sieve grid metadata and preserves single-fit basis size", {
  skip_if_runtime_unavailable(c("glmnet", "Sieve"))
  data <- make_binary_cate_data(seed = 16)
  base_learner <- make_sl3_regression_learner(stats::gaussian())
  task <- make_cate_task_for_test(data)
  default_grid <- getFromNamespace("make_ep_basis_grid", "hte3")(task)
  explicit_grid <- c(2, 4, 6)

  default_model <- fit_cate(
    task,
    method = "ep",
    base_learner = base_learner,
    cross_validate = TRUE,
    cv_control = list(V = 2)
  )

  expect_true(default_model$cross_validated)
  expect_equal(default_model$selection_summary$ep_basis_grid, default_grid)
  expect_identical(default_model$selection_summary$selected_method, "ep")
  expect_identical(default_model$selection_summary$selected_ep_targeting_style, "dr")
  expect_null(default_model$selection_summary$selected_ep_r_targeting_basis)
  expect_true(default_model$selection_summary$selected_sieve_num_basis %in% default_grid)

  explicit_model <- fit_cate(
    task,
    method = "ep",
    base_learner = base_learner,
    cross_validate = TRUE,
    cv_control = list(V = 2),
    sieve_basis_grid = explicit_grid
  )

  expect_equal(explicit_model$selection_summary$ep_basis_grid, explicit_grid)
  expect_true(explicit_model$selection_summary$selected_sieve_num_basis %in% explicit_grid)

  summary_output <- capture.output(print(summary(explicit_model)))
  expect_true(any(grepl("Selected EP basis size", summary_output, fixed = TRUE)))
  expect_true(any(grepl("EP targeting style", summary_output, fixed = TRUE)))
  expect_true(any(grepl("EP basis grid", summary_output, fixed = TRUE)))

  single_model <- fit_cate(
    task,
    method = "ep",
    base_learner = base_learner,
    cross_validate = FALSE,
    sieve_num_basis = 9
  )

  expect_null(single_model$selection_summary)
  expect_identical(single_model$learner$params$sieve_num_basis, 9)
  single_summary_output <- capture.output(print(summary(single_model)))
  expect_true(any(grepl("EP targeting style: dr", single_summary_output, fixed = TRUE)))
})

test_that("EP-R wrapper CV exposes targeting metadata and basis defaults", {
  skip_if_runtime_unavailable(c("glmnet", "Sieve"))
  data <- make_binary_cate_data(seed = 27)
  base_learner <- make_sl3_regression_learner(stats::gaussian())
  task <- make_cate_task_for_test(
    data,
    modifiers = c("W1"),
    confounders = c("W1", "W2", "W3")
  )
  make_ep_basis_grid <- getFromNamespace("make_ep_basis_grid", "hte3")
  extract_ep_metadata <- getFromNamespace("extract_ep_learner_metadata", "hte3")

  default_dr_grid <- make_ep_basis_grid(task)
  full_w_grid <- make_ep_basis_grid(task, targeting_style = "r", r_targeting_basis = "full_w")
  propensity_grid <- make_ep_basis_grid(task, targeting_style = "r", r_targeting_basis = "v_plus_propensity")

  expect_equal(default_dr_grid, c(1L, 2L, 4L, 6L, 8L))
  expect_equal(full_w_grid, c(3L, 6L, 12L, 18L, 24L))
  expect_equal(propensity_grid, c(2L, 4L, 8L, 12L, 16L))

  default_ep_r_model <- suppressWarnings(fit_cate(
    task,
    method = "ep",
    base_learner = base_learner,
    cross_validate = TRUE,
    cv_control = list(V = 2),
    ep_targeting_style = "r"
  ))

  expect_equal(default_ep_r_model$selection_summary$ep_basis_grid, propensity_grid)
  expect_identical(default_ep_r_model$selection_summary$selected_ep_targeting_style, "r")
  expect_identical(default_ep_r_model$selection_summary$selected_ep_r_targeting_basis, "v_plus_propensity")

  full_w_model <- suppressWarnings(fit_cate(
    task,
    method = "ep",
    base_learner = base_learner,
    cross_validate = TRUE,
    cv_control = list(V = 2),
    ep_targeting_style = "r",
    ep_r_targeting_basis = "full_w"
  ))

  expect_equal(full_w_model$selection_summary$ep_basis_grid, full_w_grid)
  expect_identical(full_w_model$selection_summary$selected_ep_targeting_style, "r")
  expect_identical(full_w_model$selection_summary$selected_ep_r_targeting_basis, "full_w")
  expect_true(full_w_model$selection_summary$selected_sieve_num_basis %in% full_w_grid)

  propensity_model <- suppressWarnings(fit_cate(
    task,
    method = "ep",
    base_learner = base_learner,
    cross_validate = TRUE,
    cv_control = list(V = 2),
    sieve_interaction_order = 1,
    ep_targeting_style = "r",
    ep_r_targeting_basis = "v_plus_propensity"
  ))

  expect_equal(propensity_model$selection_summary$ep_basis_grid, propensity_grid)
  expect_identical(propensity_model$selection_summary$selected_ep_targeting_style, "r")
  expect_identical(propensity_model$selection_summary$selected_ep_r_targeting_basis, "v_plus_propensity")
  expect_true(propensity_model$selection_summary$selected_sieve_num_basis %in% propensity_grid)

  summary_output <- capture.output(print(summary(propensity_model)))
  expect_true(any(grepl("EP targeting style: r", summary_output, fixed = TRUE)))
  expect_true(any(grepl("EP-R first-stage basis: v_plus_propensity", summary_output, fixed = TRUE)))

  prop_metadata <- extract_ep_metadata(
    Lrnr_cate_EP$new(
      base_learner = base_learner,
      targeting_style = "r",
      r_targeting_basis = "v_plus_propensity",
      sieve_interaction_order = 1
    )
  )
  dr_metadata <- extract_ep_metadata(
    Lrnr_cate_EP$new(
      base_learner = base_learner,
      targeting_style = "dr",
      sieve_interaction_order = 1
    )
  )
  default_ep_r_metadata <- extract_ep_metadata(
    Lrnr_cate_EP$new(
      base_learner = base_learner,
      targeting_style = "r"
    )
  )
  expect_identical(prop_metadata$interaction_order, 2L)
  expect_identical(dr_metadata$interaction_order, 1L)
  expect_identical(default_ep_r_metadata$r_targeting_basis, "v_plus_propensity")
})

test_that("wrapper-level EP selection expands over targeting styles", {
  skip_if_runtime_unavailable(c("glmnet", "Sieve"))
  data <- make_binary_cate_data(seed = 271)
  base_learner <- make_sl3_regression_learner(stats::gaussian())
  task <- make_cate_task_for_test(
    data,
    modifiers = c("W1"),
    confounders = c("W1", "W2", "W3")
  )

  ep_only_model <- suppressWarnings(fit_cate(
    task,
    method = "ep",
    base_learner = base_learner,
    cross_validate = TRUE,
    cv_control = list(V = 2),
    sieve_basis_grid = 2,
    ep_targeting_style = c("r", "dr", "r")
  ))

  expect_identical(ep_only_model$selection_summary$selected_method, "ep")
  expect_equal(ep_only_model$selection_summary$ep_basis_grid, 2)
  expect_true(ep_only_model$selection_summary$selected_ep_targeting_style %in% c("dr", "r"))
  expect_true(any(grepl("style=dr", ep_only_model$selection_summary$candidate_labels, fixed = TRUE)))
  expect_true(any(grepl("style=r", ep_only_model$selection_summary$candidate_labels, fixed = TRUE)))

  mixed_model <- suppressWarnings(fit_cate(
    task,
    method = c("ep", "r"),
    base_learner = base_learner,
    cross_validate = TRUE,
    cv_control = list(V = 2),
    sieve_basis_grid = 2,
    ep_targeting_style = c("dr", "r")
  ))

  expect_true(mixed_model$selection_summary$selected_method %in% c("ep", "r"))
  expect_true(any(grepl("style=dr", mixed_model$selection_summary$candidate_labels, fixed = TRUE)))
  expect_true(any(grepl("style=r", mixed_model$selection_summary$candidate_labels, fixed = TRUE)))
  expect_true(any(mixed_model$selection_summary$candidate_labels == "r"))
})

test_that("EP wrapper CV validates sieve grid usage for CATE and CRR", {
  skip_if_runtime_unavailable(c("glmnet", "Sieve"))
  cate_data <- make_binary_cate_data(seed = 17)
  categorical_cate_data <- make_categorical_cate_data(seed = 171)
  crr_data <- make_binary_crr_data(seed = 18)
  cate_task <- make_cate_task_for_test(cate_data)
  categorical_task <- make_cate_task_for_test(
    categorical_cate_data,
    treatment_type = "categorical"
  )
  crr_task <- make_crr_task_for_test(crr_data)
  gaussian_learner <- make_sl3_regression_learner(stats::gaussian())
  binomial_learner <- make_sl3_regression_learner(stats::binomial())

  expect_error(
    fit_cate(
      cate_task,
      method = "dr",
      base_learner = gaussian_learner,
      cross_validate = TRUE,
      cv_control = list(V = 2),
      sieve_basis_grid = c(2, 4)
    ),
    "`sieve_basis_grid` can only be used"
  )

  expect_error(
    fit_cate(
      cate_task,
      method = "ep",
      base_learner = gaussian_learner,
      cross_validate = FALSE,
      sieve_basis_grid = c(2, 4)
    ),
    "`sieve_basis_grid` can only be used"
  )

  expect_error(
    fit_cate(
      cate_task,
      method = "ep",
      base_learner = gaussian_learner,
      cross_validate = TRUE,
      cv_control = list(V = 2),
      sieve_basis_grid = c(2, 2)
    ),
    "at least two distinct EP candidates"
  )

  expect_error(
    fit_cate(
      cate_task,
      method = "ep",
      base_learner = gaussian_learner,
      cross_validate = TRUE,
      cv_control = list(V = 2),
      sieve_basis_grid = 2,
      screen_basis_with_lasso = TRUE,
      ep_targeting_style = c("dr", "r")
    ),
    "does not include"
  )

  expect_error(
    fit_cate(
      categorical_task,
      method = "ep",
      base_learner = gaussian_learner,
      cross_validate = TRUE,
      cv_control = list(V = 2),
      sieve_basis_grid = 2,
      ep_targeting_style = c("dr", "r")
    ),
    "currently only supported for binary-treatment"
  )

  expect_error(
    fit_cate(
      cate_task,
      method = "ep",
      base_learner = gaussian_learner,
      cross_validate = TRUE,
      cv_control = list(V = 2),
      sieve_basis_grid = c(0, 2)
    ),
    "positive basis sizes"
  )

  expect_error(
    fit_crr(
      crr_task,
      method = "ep",
      base_learner = binomial_learner,
      cross_validate = TRUE,
      cv_control = list(V = 2),
      sieve_basis_grid = c(3, 3)
    ),
    "at least two distinct EP candidates"
  )
})

test_that("EP wrapper CV exposes sieve grid metadata for CRR and make_ep_stack accepts style portfolios", {
  skip_if_runtime_unavailable(c("glmnet", "Sieve"))
  crr_data <- make_binary_crr_data(seed = 19)
  crr_task <- make_crr_task_for_test(crr_data)
  base_learner <- make_sl3_regression_learner(stats::binomial())
  explicit_grid <- c(3, 6, 9)

  model <- fit_crr(
    crr_task,
    method = "ep",
    base_learner = base_learner,
    cross_validate = TRUE,
    cv_control = list(V = 2),
    sieve_basis_grid = explicit_grid
  )

  expect_equal(model$selection_summary$ep_basis_grid, explicit_grid)
  expect_true(model$selection_summary$selected_sieve_num_basis %in% explicit_grid)

  ep_stack <- make_ep_stack(
    base_learner = make_sl3_regression_learner(stats::gaussian()),
    hte3_task = make_cate_task_for_test(make_binary_cate_data(seed = 20)),
    sieve_basis_grid = c(2, 5),
    targeting_style = c("r", "dr"),
    r_targeting_basis = "v_plus_propensity"
  )
  expect_s3_class(ep_stack, "Stack")
  expect_length(ep_stack$params$learners, 4L)
  expect_identical(ep_stack$params$learners[[1]]$params$targeting_style, "dr")
  expect_identical(ep_stack$params$learners[[3]]$params$targeting_style, "r")
  expect_identical(ep_stack$params$learners[[1]]$params$r_targeting_basis, "v_plus_propensity")
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
