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
  suppress_internal_cv_fallback_warnings(
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
  )
}

new_hte3_model <- function(kind, method, learner, training_task, cross_validated, call, selection_summary = NULL, engine = "sl3") {
  structure(
    list(
      kind = kind,
      method = method,
      learner = learner,
      training_task = training_task,
      cross_validated = cross_validated,
      call = call,
      selection_summary = selection_summary,
      engine = engine
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

validate_ep_cv_configuration <- function(methods, cross_validate, sieve_basis_grid, ep_targeting_style = "dr") {
  if (is.null(sieve_basis_grid)) {
    return(NULL)
  }

  normalized_grid <- normalize_sieve_basis_grid(sieve_basis_grid)

  if (!("ep" %in% methods) || !isTRUE(cross_validate)) {
    stop(
      "`sieve_basis_grid` can only be used when `method` includes \"ep\" and EP selection is being cross-validated.",
      call. = FALSE
    )
  }

  ep_variant_count <- if ("ep" %in% methods) length(normalize_ep_targeting_styles(ep_targeting_style)) else 1L

  if (length(methods) == 1L && identical(methods[[1L]], "ep") && (length(normalized_grid) * ep_variant_count) < 2L) {
    stop(
      paste(
        "`sieve_basis_grid` must yield at least two distinct EP candidates when cross-validating only EP learners.",
        "Provide multiple basis sizes or set `cross_validate = FALSE`."
      ),
      call. = FALSE
    )
  }

  normalized_grid
}

selection_summary_fields <- function(selection_summary, learner = NULL) {
  if (is.null(selection_summary)) {
    ep_metadata <- if (is.null(learner)) NULL else extract_ep_learner_metadata(learner)
    return(list(
      selected_candidate = NULL,
      selected_method = if (is.null(ep_metadata)) NULL else "ep",
      selected_sieve_num_basis = if (is.null(ep_metadata)) NULL else ep_metadata$sieve_num_basis,
      ep_basis_grid = NULL,
      selected_ep_targeting_style = if (is.null(ep_metadata)) NULL else ep_metadata$targeting_style,
      selected_ep_r_targeting_basis = if (is.null(ep_metadata)) NULL else ep_metadata$r_targeting_basis
    ))
  }

  list(
    selected_candidate = selection_summary$selected_candidate,
    selected_method = selection_summary$selected_method,
    selected_sieve_num_basis = selection_summary$selected_sieve_num_basis,
    ep_basis_grid = selection_summary$ep_basis_grid,
    selected_ep_targeting_style = selection_summary$selected_ep_targeting_style,
    selected_ep_r_targeting_basis = selection_summary$selected_ep_r_targeting_basis
  )
}

build_cate_candidates <- function(task,
                                  methods,
                                  base_learner,
                                  contrast,
                                  cross_validate,
                                  sieve_num_basis,
                                  sieve_basis_grid,
                                  sieve_interaction_order,
                                  screen_basis_with_lasso,
                                  ep_targeting_style,
                                  ep_r_targeting_basis) {
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
        sieve_basis_grid = sieve_basis_grid,
        sieve_interaction_order = sieve_interaction_order,
        screen_basis_with_lasso = screen_basis_with_lasso,
        targeting_style = ep_targeting_style,
        r_targeting_basis = ep_r_targeting_basis
      )
    } else {
      list(Lrnr_cate_EP$new(
        base_learner = base_learner,
        sieve_num_basis = sieve_num_basis,
        sieve_interaction_order = sieve_interaction_order,
        screen_basis_with_lasso = screen_basis_with_lasso,
        targeting_style = ep_targeting_style[[1L]],
        r_targeting_basis = ep_r_targeting_basis,
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
                                 sieve_basis_grid,
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
        sieve_basis_grid = sieve_basis_grid,
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
                     sieve_basis_grid = NULL,
                     sieve_interaction_order = 3,
                     screen_basis_with_lasso = FALSE,
                     ep_targeting_style = "dr",
                     ep_r_targeting_basis = "v_plus_propensity") {
  validate_hte_task(task)
  methods <- if (missing(method)) "dr" else normalize_method_spec(method, supported_cate_methods())
  task_treatment_type <- get_task_treatment_type(task)
  ep_targeting_styles <- if ("ep" %in% methods) normalize_ep_targeting_styles(ep_targeting_style) else ep_targeting_style

  if (task_treatment_type == "continuous" && any(methods %in% c("dr", "t", "ep"))) {
    stop("Only `method = \"r\"` is supported for continuous-treatment CATE tasks because the current high-level continuous-treatment path uses the partially linear R-learner effect model `A * tau(X)`.", call. = FALSE)
  }
  if ("ep" %in% methods && "r" %in% ep_targeting_styles && !identical(task_treatment_type, "binomial")) {
    stop("`ep_targeting_style = \"r\"` is currently only supported for binary-treatment CATE tasks.", call. = FALSE)
  }
  if ("ep" %in% methods && "r" %in% ep_targeting_styles && isTRUE(screen_basis_with_lasso)) {
    stop("`screen_basis_with_lasso = TRUE` is only supported when `ep_targeting_style` does not include \"r\".", call. = FALSE)
  }

  needs_contrast <- any(methods %in% c("dr", "t", "ep"))
  contrast <- if (needs_contrast) {
    resolve_treatment_levels(task, treatment_level, control_level)
  } else {
    NULL
  }
  cross_validate <- prepare_cv_flag(
    base_learner,
    cross_validate,
    candidate_count = count_cate_wrapper_candidates(methods, ep_targeting_styles)
  )
  sieve_basis_grid <- validate_ep_cv_configuration(
    methods,
    cross_validate,
    sieve_basis_grid,
    ep_targeting_style = ep_targeting_styles
  )

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
    sieve_basis_grid = sieve_basis_grid,
    sieve_interaction_order = sieve_interaction_order,
    screen_basis_with_lasso = screen_basis_with_lasso,
    ep_targeting_style = ep_targeting_styles,
    ep_r_targeting_basis = ep_r_targeting_basis
  )

  trained <- suppress_internal_cv_fallback_warnings(
    if (cross_validate) {
      cv_fit <- cross_validate_cate(
        learners,
        task,
        cv_control = cv_control,
        treatment_level = if (is.null(contrast)) NULL else contrast$treatment_level,
        control_level = if (is.null(contrast)) NULL else contrast$control_level
      )
      selection_summary <- cv_fit$selection_summary
      cv_fit$lrnr_sl
    } else {
      selection_summary <- NULL
      learners[[1]]$train(task)
    }
  )

  new_hte3_model("cate", format_method_label(methods), trained, task, cross_validate, match.call(), selection_summary = selection_summary)
}

fit_crr <- function(task,
                    method = supported_crr_methods(),
                    base_learner = get_autoML(),
                    treatment_level = NULL,
                    control_level = NULL,
                    cross_validate = NULL,
                    cv_control = NULL,
                    sieve_num_basis = NULL,
                    sieve_basis_grid = NULL,
                    sieve_interaction_order = 3) {
  validate_hte_task(task)
  methods <- if (missing(method)) "ep" else normalize_method_spec(method, supported_crr_methods())
  task_treatment_type <- get_task_treatment_type(task)
  if (identical(task_treatment_type, "continuous")) {
    stop("CRR models are only supported for binary or categorical treatments.", call. = FALSE)
  }
  contrast <- resolve_treatment_levels(task, treatment_level, control_level)
  cross_validate <- prepare_cv_flag(base_learner, cross_validate, candidate_count = length(methods))
  sieve_basis_grid <- validate_ep_cv_configuration(methods, cross_validate, sieve_basis_grid)

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
    sieve_basis_grid = sieve_basis_grid,
    sieve_interaction_order = sieve_interaction_order
  )

  trained <- suppress_internal_cv_fallback_warnings(
    if (cross_validate) {
      cv_fit <- cross_validate_crr(
        learners,
        task,
        cv_control = cv_control,
        treatment_level = contrast$treatment_level,
        control_level = contrast$control_level
      )
      selection_summary <- cv_fit$selection_summary
      cv_fit$lrnr_sl
    } else {
      selection_summary <- NULL
      learners[[1]]$train(task)
    }
  )

  new_hte3_model("crr", format_method_label(methods), trained, task, cross_validate, match.call(), selection_summary = selection_summary)
}

predict.hte3_model <- function(object, new_data = NULL, ...) {
  if (is.null(new_data)) {
    return(object$learner$predict(object$training_task))
  }

  predict_hte3(object$learner, new_data)
}

print.hte3_model <- function(x, ...) {
  selection_fields <- selection_summary_fields(x$selection_summary, learner = x$learner)
  cat("<hte3_model>\n")
  cat(sprintf("  Target: %s\n", toupper(x$kind)))
  if (!is.null(x$engine)) {
    cat(sprintf("  Engine: %s\n", x$engine))
  }
  cat(sprintf("  Method: %s\n", x$method))
  cat(sprintf("  Cross-validated: %s\n", if (isTRUE(x$cross_validated)) "yes" else "no"))
  if (!is.null(selection_fields$selected_candidate)) {
    cat(sprintf("  Selected candidate: %s\n", selection_fields$selected_candidate))
  }
  if (!is.null(selection_fields$selected_method)) {
    cat(sprintf("  Selected method: %s\n", selection_fields$selected_method))
  }
  if (!is.null(selection_fields$selected_sieve_num_basis)) {
    cat(sprintf("  Selected EP basis size: %s\n", selection_fields$selected_sieve_num_basis))
  }
  if (!is.null(selection_fields$selected_ep_targeting_style)) {
    cat(sprintf("  EP targeting style: %s\n", selection_fields$selected_ep_targeting_style))
  }
  if (!is.null(selection_fields$selected_ep_r_targeting_basis)) {
    cat(sprintf("  EP-R first-stage basis: %s\n", selection_fields$selected_ep_r_targeting_basis))
  }
  if (!is.null(selection_fields$ep_basis_grid)) {
    cat(sprintf("  EP basis grid: %s\n", paste(selection_fields$ep_basis_grid, collapse = ", ")))
  }
  cat(sprintf("  Modifiers: %s\n", paste(x$training_task$npsem$modifiers$variables, collapse = ", ")))
  invisible(x)
}

summary.hte3_model <- function(object, ...) {
  selection_fields <- selection_summary_fields(object$selection_summary, learner = object$learner)
  structure(
    list(
      kind = object$kind,
      engine = object$engine,
      method = object$method,
      cross_validated = object$cross_validated,
      n = object$training_task$nrow,
      modifiers = object$training_task$npsem$modifiers$variables,
      treatment = object$training_task$npsem$treatment$variables,
      selection_summary = object$selection_summary,
      selected_candidate = selection_fields$selected_candidate,
      selected_method = selection_fields$selected_method,
      selected_sieve_num_basis = selection_fields$selected_sieve_num_basis,
      ep_basis_grid = selection_fields$ep_basis_grid,
      selected_ep_targeting_style = selection_fields$selected_ep_targeting_style,
      selected_ep_r_targeting_basis = selection_fields$selected_ep_r_targeting_basis
    ),
    class = "summary.hte3_model"
  )
}

print.summary.hte3_model <- function(x, ...) {
  cat("<summary.hte3_model>\n")
  cat(sprintf("  Target: %s\n", toupper(x$kind)))
  if (!is.null(x$engine)) {
    cat(sprintf("  Engine: %s\n", x$engine))
  }
  cat(sprintf("  Method: %s\n", x$method))
  cat(sprintf("  Cross-validated: %s\n", if (isTRUE(x$cross_validated)) "yes" else "no"))
  cat(sprintf("  Rows: %d\n", x$n))
  if (!is.null(x$selected_candidate)) {
    cat(sprintf("  Selected candidate: %s\n", x$selected_candidate))
  }
  if (!is.null(x$selected_method)) {
    cat(sprintf("  Selected method: %s\n", x$selected_method))
  }
  if (!is.null(x$selected_sieve_num_basis)) {
    cat(sprintf("  Selected EP basis size: %s\n", x$selected_sieve_num_basis))
  }
  if (!is.null(x$selected_ep_targeting_style)) {
    cat(sprintf("  EP targeting style: %s\n", x$selected_ep_targeting_style))
  }
  if (!is.null(x$selected_ep_r_targeting_basis)) {
    cat(sprintf("  EP-R first-stage basis: %s\n", x$selected_ep_r_targeting_basis))
  }
  if (!is.null(x$ep_basis_grid)) {
    cat(sprintf("  EP basis grid: %s\n", paste(x$ep_basis_grid, collapse = ", ")))
  }
  cat(sprintf("  Modifiers: %s\n", paste(x$modifiers, collapse = ", ")))
  cat(sprintf("  Treatment variable: %s\n", paste(x$treatment, collapse = ", ")))
  invisible(x)
}
