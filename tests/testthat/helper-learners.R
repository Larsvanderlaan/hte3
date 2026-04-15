skip_if_runtime_unavailable <- function(extra = character()) {
  packages <- c("sl3", "tmle3", extra)
  missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0L) {
    skip(sprintf("Missing runtime package(s): %s", paste(missing_packages, collapse = ", ")))
  }
}

make_sl3_regression_learner <- function(family = stats::gaussian()) {
  exports <- getNamespaceExports("sl3")

  if ("Lrnr_glm_fast" %in% exports) {
    learner <- tryCatch(
      sl3::Lrnr_glm_fast$new(family = family),
      error = function(...) NULL
    )
    if (!is.null(learner)) {
      return(learner)
    }
  }

  if ("Lrnr_glm" %in% exports) {
    return(sl3::Lrnr_glm$new(family = family))
  }

  skip("Need `sl3::Lrnr_glm_fast` or `sl3::Lrnr_glm` for learner audit tests.")
}

make_sl3_interaction_learner <- function() {
  exports <- getNamespaceExports("sl3")

  if (!("Lrnr_earth" %in% exports) || !requireNamespace("earth", quietly = TRUE)) {
    skip("Need `sl3::Lrnr_earth` and `earth` for pooled T-learner interaction tests.")
  }

  sl3::Lrnr_earth$new(degree = 2)
}

make_binary_cate_data <- function(n = 240, seed = 1) {
  set.seed(seed)
  W1 <- runif(n, -1, 1)
  W2 <- runif(n, -1, 1)
  W3 <- runif(n, -1, 1)
  pi1 <- stats::plogis(0.25 + 0.5 * W2 - 0.25 * W3)
  tau <- 0.4 + 0.9 * W1
  mu0 <- 0.25 + 0.5 * W2 - 0.3 * W3
  A <- stats::rbinom(n, 1, pi1)
  Y <- mu0 + A * tau + stats::rnorm(n, sd = 0.15)

  data.table::data.table(W1 = W1, W2 = W2, W3 = W3, A = A, Y = Y, tau = tau)
}

make_categorical_cate_data <- function(n = 260, seed = 2) {
  set.seed(seed)
  W1 <- runif(n, -1, 1)
  W2 <- runif(n, -1, 1)
  W3 <- runif(n, -1, 1)
  score_mid <- -0.2 + 0.3 * W2
  score_high <- 0.15 + 0.5 * W1 - 0.15 * W3
  exp_scores <- exp(cbind(control = 0, mid = score_mid, high = score_high))
  probs <- exp_scores / rowSums(exp_scores)
  treatment_levels <- c("control", "mid", "high")
  A <- factor(
    apply(probs, 1, function(p) sample(treatment_levels, size = 1, prob = p)),
    levels = treatment_levels
  )
  mu0 <- 0.2 + 0.4 * W2 - 0.2 * W3
  mu_mid <- mu0 + 0.2 + 0.2 * W1
  mu_high <- mu0 + 0.6 + 0.7 * W1
  Y_mean <- ifelse(A == "mid", mu_mid, ifelse(A == "high", mu_high, mu0))
  Y <- Y_mean + stats::rnorm(n, sd = 0.15)

  data.table::data.table(W1 = W1, W2 = W2, W3 = W3, A = A, Y = Y)
}

make_binary_crr_data <- function(n = 320, seed = 3) {
  set.seed(seed)
  W1 <- runif(n, -1, 1)
  W2 <- runif(n, -1, 1)
  W3 <- runif(n, -1, 1)
  pi1 <- stats::plogis(0.1 + 0.4 * W2 - 0.25 * W3)
  A <- stats::rbinom(n, 1, pi1)
  log_rr <- 0.12 + 0.35 * W1
  p0 <- stats::plogis(-2 + 0.35 * W2 - 0.2 * W3)
  p1 <- pmin(0.8, p0 * exp(log_rr))
  Y <- stats::rbinom(n, 1, ifelse(A == 1, p1, p0))

  data.table::data.table(W1 = W1, W2 = W2, W3 = W3, A = A, Y = Y, log_rr = log(p1 / p0))
}

make_categorical_crr_data <- function(n = 360, seed = 4) {
  set.seed(seed)
  W1 <- runif(n, -1, 1)
  W2 <- runif(n, -1, 1)
  W3 <- runif(n, -1, 1)
  score_mid <- -0.2 + 0.25 * W2
  score_high <- 0.1 + 0.45 * W1 - 0.2 * W3
  exp_scores <- exp(cbind(control = 0, mid = score_mid, high = score_high))
  probs <- exp_scores / rowSums(exp_scores)
  treatment_levels <- c("control", "mid", "high")
  A <- factor(
    apply(probs, 1, function(p) sample(treatment_levels, size = 1, prob = p)),
    levels = treatment_levels
  )
  p0 <- stats::plogis(-2 + 0.3 * W2 - 0.1 * W3)
  p_mid <- pmin(0.75, p0 * exp(0.08 + 0.15 * W1))
  p_high <- pmin(0.8, p0 * exp(0.18 + 0.3 * W1))
  probability <- ifelse(A == "mid", p_mid, ifelse(A == "high", p_high, p0))
  Y <- stats::rbinom(n, 1, probability)

  data.table::data.table(W1 = W1, W2 = W2, W3 = W3, A = A, Y = Y)
}

make_cate_task_for_test <- function(data,
                                    treatment = "A",
                                    outcome = "Y",
                                    treatment_type = "default",
                                    modifiers = c("W1", "W2"),
                                    confounders = c("W1", "W2", "W3"),
                                    propensity = NULL,
                                    outcome_regression = NULL,
                                    outcome_mean = NULL) {
  propensity_learner <- make_sl3_regression_learner(stats::binomial())
  outcome_learner <- make_sl3_regression_learner(stats::gaussian())
  mean_learner <- make_sl3_regression_learner(stats::gaussian())

  hte_task(
    data = data,
    modifiers = modifiers,
    confounders = confounders,
    treatment = treatment,
    outcome = outcome,
    treatment_type = treatment_type,
    propensity = propensity,
    outcome_regression = outcome_regression,
    outcome_mean = outcome_mean,
    propensity_learner = propensity_learner,
    outcome_learner = outcome_learner,
    mean_learner = mean_learner,
    cross_fit = FALSE
  )
}

make_crr_task_for_test <- function(data,
                                   treatment = "A",
                                   outcome = "Y",
                                   treatment_type = "default",
                                   modifiers = c("W1", "W2"),
                                   confounders = c("W1", "W2", "W3"),
                                   propensity = NULL,
                                   outcome_regression = NULL,
                                   outcome_mean = NULL) {
  propensity_learner <- make_sl3_regression_learner(stats::binomial())
  outcome_learner <- make_sl3_regression_learner(stats::binomial())
  mean_learner <- make_sl3_regression_learner(stats::binomial())

  hte_task(
    data = data,
    modifiers = modifiers,
    confounders = confounders,
    treatment = treatment,
    outcome = outcome,
    treatment_type = treatment_type,
    propensity = propensity,
    outcome_regression = outcome_regression,
    outcome_mean = outcome_mean,
    propensity_learner = propensity_learner,
    outcome_learner = outcome_learner,
    mean_learner = mean_learner,
    cross_fit = FALSE
  )
}

expect_monotone_signal <- function(predictions, signal, threshold = 0) {
  groups <- signal > stats::median(signal)
  expect_gt(mean(predictions[groups]), mean(predictions[!groups]) + threshold)
}

expect_positive_average_signal <- function(predictions, threshold = 0) {
  expect_gt(mean(predictions), threshold)
}
