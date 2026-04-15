hte3_example_data <- function(n = 200, seed = 123) {
  if (!is.numeric(n) || length(n) != 1L || n < 10) {
    stop("`n` must be a single integer greater than or equal to 10.", call. = FALSE)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  W1 <- runif(n, -1, 1)
  W2 <- runif(n, -1, 1)
  W3 <- runif(n, -1, 1)
  pi1 <- plogis(0.5 * W1 - 0.25 * W2 + 0.5 * W3)
  tau <- 0.75 + W1 - W2 + sin(2 * W3)
  mu0 <- 0.5 + W1 + 0.5 * W2^2 + cos(W3)
  mu1 <- mu0 + tau
  A <- rbinom(n, 1, pi1)
  Y <- rnorm(n, mean = mu0 + A * tau, sd = 0.25)
  Y_binary <- rbinom(n, 1, plogis(-0.5 + 0.5 * W1 + 0.25 * W2 + A * 0.4))

  data.table(
    W1 = W1,
    W2 = W2,
    W3 = W3,
    A = A,
    Y = Y,
    Y_binary = Y_binary,
    tau = tau,
    pi1 = pi1,
    mu0 = mu0,
    mu1 = mu1
  )
}

hte_task <- function(data,
                     modifiers,
                     confounders = modifiers,
                     treatment,
                     outcome,
                     id = NULL,
                     weights = NULL,
                     treatment_type = c("default", "binomial", "categorical", "continuous"),
                     propensity = NULL,
                     outcome_regression = NULL,
                     outcome_mean = NULL,
                     propensity_learner = get_autoML(),
                     outcome_learner = get_autoML(),
                     mean_learner = NULL,
                     cross_fit = TRUE,
                     folds = 10,
                     for_prediction = FALSE,
                     ...) {
  make_hte3_Task_tx(
    data = data,
    modifiers = modifiers,
    confounders = confounders,
    treatment = treatment,
    outcome = outcome,
    id = id,
    weights = weights,
    treatment_type = match.arg(treatment_type),
    pi.hat = propensity,
    mu.hat = outcome_regression,
    m.hat = outcome_mean,
    learner_pi = propensity_learner,
    learner_mu = outcome_learner,
    learner_m = mean_learner,
    cross_fit_and_cv = cross_fit,
    folds = folds,
    for_prediction = for_prediction,
    ...
  )
}

new_hte3_model <- function(kind, method, learner, training_task, cross_validated, call) {
  structure(
    list(
      kind = kind,
      method = method,
      learner = learner,
      training_task = training_task,
      cross_validated = cross_validated,
      call = call
    ),
    class = "hte3_model"
  )
}

validate_hte_task <- function(task) {
  if (!inherits(task, "hte3_Task")) {
    stop("`task` must be an `hte3_Task` created by `hte_task()` or `make_hte3_Task_tx()`.", call. = FALSE)
  }

  invisible(task)
}

normalize_method_spec <- function(method, choices) {
  method <- unique(as.character(method))
  if (length(method) == 0L) {
    stop("`method` must contain at least one supported method.", call. = FALSE)
  }
  unknown <- setdiff(method, choices)
  if (length(unknown) > 0L) {
    stop(
      sprintf("Unsupported method(s): %s. Choices are: %s.", paste(unknown, collapse = ", "), paste(choices, collapse = ", ")),
      call. = FALSE
    )
  }

  method
}

prepare_cv_flag <- function(base_learner, cross_validate, candidate_count) {
  explicit_flag <- !is.null(cross_validate)
  cross_validate <- default_cv_flag(base_learner, cross_validate)

  if (candidate_count > 1L) {
    if (explicit_flag && !isTRUE(cross_validate)) {
      stop("`cross_validate` must be TRUE when fitting a learner portfolio.", call. = FALSE)
    }
    cross_validate <- TRUE
  }

  cross_validate
}

format_method_label <- function(methods) {
  if (length(methods) == 1L) {
    return(methods)
  }

  sprintf("portfolio(%s)", paste(methods, collapse = ", "))
}

build_cate_candidates <- function(task,
                                  methods,
                                  base_learner,
                                  contrast,
                                  cross_validate,
                                  sieve_num_basis,
                                  sieve_interaction_order,
                                  screen_basis_with_lasso) {
  learners <- list()

  if ("dr" %in% methods) {
    learners <- c(learners, list(Lrnr_cate_DR$new(
      base_learner = base_learner,
      treatment_level = contrast$treatment_level,
      control_level = contrast$control_level
    )))
  }

  if ("r" %in% methods) {
    learners <- c(learners, list(Lrnr_cate_R$new(base_learner = base_learner)))
  }

  if ("t" %in% methods) {
    learners <- c(learners, list(Lrnr_cate_T$new(
      base_learner = base_learner,
      treatment_level = contrast$treatment_level,
      control_level = contrast$control_level
    )))
  }

  if ("ep" %in% methods) {
    ep_candidates <- if (cross_validate) {
      make_cate_ep_candidates(
        base_learner,
        task,
        contrast$treatment_level,
        contrast$control_level,
        sieve_interaction_order = sieve_interaction_order,
        screen_basis_with_lasso = screen_basis_with_lasso
      )
    } else {
      list(Lrnr_cate_EP$new(
        base_learner = base_learner,
        sieve_num_basis = sieve_num_basis,
        sieve_interaction_order = sieve_interaction_order,
        screen_basis_with_lasso = screen_basis_with_lasso,
        treatment_level = contrast$treatment_level,
        control_level = contrast$control_level
      ))
    }
    learners <- c(learners, ep_candidates)
  }

  learners
}

build_crr_candidates <- function(task,
                                 methods,
                                 base_learner,
                                 contrast,
                                 cross_validate,
                                 sieve_num_basis,
                                 sieve_interaction_order) {
  learners <- list()

  if ("ipw" %in% methods) {
    learners <- c(learners, list(Lrnr_crr_IPW$new(
      base_learner = base_learner,
      treatment_level = contrast$treatment_level,
      control_level = contrast$control_level
    )))
  }

  if ("t" %in% methods) {
    learners <- c(learners, list(Lrnr_crr_T$new(
      base_learner = base_learner,
      treatment_level = contrast$treatment_level,
      control_level = contrast$control_level
    )))
  }

  if ("ep" %in% methods) {
    ep_candidates <- if (cross_validate) {
      make_crr_ep_candidates(
        base_learner,
        task,
        contrast$treatment_level,
        contrast$control_level,
        sieve_interaction_order = sieve_interaction_order
      )
    } else {
      list(Lrnr_crr_EP$new(
        base_learner = base_learner,
        sieve_num_basis = sieve_num_basis,
        sieve_interaction_order = sieve_interaction_order,
        treatment_level = contrast$treatment_level,
        control_level = contrast$control_level
      ))
    }
    learners <- c(learners, ep_candidates)
  }

  learners
}

fit_cate <- function(task,
                     method = supported_cate_methods(),
                     base_learner = get_autoML(),
                     treatment_level = NULL,
                     control_level = NULL,
                     cross_validate = NULL,
                     cv_control = NULL,
                     sieve_num_basis = NULL,
                     sieve_interaction_order = 3,
                     screen_basis_with_lasso = FALSE) {
  validate_hte_task(task)
  methods <- if (missing(method)) "dr" else normalize_method_spec(method, supported_cate_methods())
  task_treatment_type <- get_task_treatment_type(task)

  if (task_treatment_type == "continuous" && any(methods %in% c("dr", "t", "ep"))) {
    stop("Only `method = \"r\"` is supported for continuous-treatment CATE tasks because the current high-level continuous-treatment path uses the partially linear R-learner effect model `A * tau(X)`.", call. = FALSE)
  }

  needs_contrast <- any(methods %in% c("dr", "t", "ep"))
  contrast <- if (needs_contrast) {
    resolve_treatment_levels(task, treatment_level, control_level)
  } else {
    NULL
  }
  cross_validate <- prepare_cv_flag(base_learner, cross_validate, candidate_count = length(methods))

  if (cross_validate && task_treatment_type == "continuous") {
    stop("Cross-validation for continuous-treatment CATE tasks is not supported because the current selector uses contrast-based DR loss.", call. = FALSE)
  }

  learners <- build_cate_candidates(
    task = task,
    methods = methods,
    base_learner = base_learner,
    contrast = contrast,
    cross_validate = cross_validate,
    sieve_num_basis = sieve_num_basis,
    sieve_interaction_order = sieve_interaction_order,
    screen_basis_with_lasso = screen_basis_with_lasso
  )

  trained <- if (cross_validate) {
    cross_validate_cate(
      learners,
      task,
      cv_control = cv_control,
      treatment_level = if (is.null(contrast)) NULL else contrast$treatment_level,
      control_level = if (is.null(contrast)) NULL else contrast$control_level
    )$lrnr_sl
  } else {
    learners[[1]]$train(task)
  }

  new_hte3_model("cate", format_method_label(methods), trained, task, cross_validate, match.call())
}

fit_crr <- function(task,
                    method = supported_crr_methods(),
                    base_learner = get_autoML(),
                    treatment_level = NULL,
                    control_level = NULL,
                    cross_validate = NULL,
                    cv_control = NULL,
                    sieve_num_basis = NULL,
                    sieve_interaction_order = 3) {
  validate_hte_task(task)
  methods <- if (missing(method)) "ep" else normalize_method_spec(method, supported_crr_methods())
  task_treatment_type <- get_task_treatment_type(task)
  if (identical(task_treatment_type, "continuous")) {
    stop("CRR models are only supported for binary or categorical treatments.", call. = FALSE)
  }
  contrast <- resolve_treatment_levels(task, treatment_level, control_level)
  cross_validate <- prepare_cv_flag(base_learner, cross_validate, candidate_count = length(methods))

  outcome <- task$get_tmle_node("outcome")
  if (any(outcome < 0, na.rm = TRUE)) {
    stop("CRR models require a non-negative outcome.", call. = FALSE)
  }

  learners <- build_crr_candidates(
    task = task,
    methods = methods,
    base_learner = base_learner,
    contrast = contrast,
    cross_validate = cross_validate,
    sieve_num_basis = sieve_num_basis,
    sieve_interaction_order = sieve_interaction_order
  )

  trained <- if (cross_validate) {
    cross_validate_crr(
      learners,
      task,
      cv_control = cv_control,
      treatment_level = contrast$treatment_level,
      control_level = contrast$control_level
    )$lrnr_sl
  } else {
    learners[[1]]$train(task)
  }

  new_hte3_model("crr", format_method_label(methods), trained, task, cross_validate, match.call())
}

predict.hte3_model <- function(object, new_data = NULL, ...) {
  if (is.null(new_data)) {
    return(object$learner$predict(object$training_task))
  }

  predict_hte3(object$learner, new_data)
}

print.hte3_model <- function(x, ...) {
  cat("<hte3_model>\n")
  cat(sprintf("  Target: %s\n", toupper(x$kind)))
  cat(sprintf("  Method: %s\n", x$method))
  cat(sprintf("  Cross-validated: %s\n", if (isTRUE(x$cross_validated)) "yes" else "no"))
  cat(sprintf("  Modifiers: %s\n", paste(x$training_task$npsem$modifiers$variables, collapse = ", ")))
  invisible(x)
}

summary.hte3_model <- function(object, ...) {
  structure(
    list(
      kind = object$kind,
      method = object$method,
      cross_validated = object$cross_validated,
      n = object$training_task$nrow,
      modifiers = object$training_task$npsem$modifiers$variables,
      treatment = object$training_task$npsem$treatment$variables
    ),
    class = "summary.hte3_model"
  )
}

print.summary.hte3_model <- function(x, ...) {
  cat("<summary.hte3_model>\n")
  cat(sprintf("  Target: %s\n", toupper(x$kind)))
  cat(sprintf("  Method: %s\n", x$method))
  cat(sprintf("  Rows: %d\n", x$n))
  cat(sprintf("  Modifiers: %s\n", paste(x$modifiers, collapse = ", ")))
  cat(sprintf("  Treatment variable: %s\n", paste(x$treatment, collapse = ", ")))
  invisible(x)
}
