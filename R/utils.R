call_with_args <- function(fun, args, silent = FALSE) {
  if (!is.function(fun)) {
    stop("`fun` must be a function.", call. = FALSE)
  }

  formal_names <- names(formals(fun))
  if (!is.null(formal_names) && !("..." %in% formal_names)) {
    args <- args[intersect(names(args), formal_names)]
  }

  runner <- function() {
    do.call(fun, args)
  }

  if (isTRUE(silent)) {
    tryCatch(runner(), error = function(err) stop(conditionMessage(err), call. = FALSE))
  } else {
    runner()
  }
}

coerce_prediction_matrix <- function(predictions) {
  if (inherits(predictions, "packed_predictions")) {
    return(as.matrix(predictions))
  }

  if (is.matrix(predictions)) {
    return(predictions)
  }

  if (is.data.frame(predictions)) {
    return(as.matrix(predictions))
  }

  if (is.list(predictions) && length(predictions) > 0) {
    if (length(predictions) == 1) {
      return(coerce_prediction_matrix(predictions[[1]]))
    }

    same_length <- vapply(predictions, length, integer(1))
    if (length(unique(same_length)) == 1L) {
      return(do.call(cbind, predictions))
    }
  }

  as.matrix(predictions)
}

subset_folds <- function(folds, row_index) {
  if (is.null(folds)) {
    return(NULL)
  }

  remap_indices <- function(indices) {
    mapped <- match(indices, row_index)
    mapped[!is.na(mapped)]
  }

  lapply(folds, function(fold) {
    if (!is.list(fold)) {
      return(fold)
    }

    for (slot in c("training_set", "validation_set", "assessment_set", "analysis_set")) {
      if (!is.null(fold[[slot]])) {
        fold[[slot]] <- remap_indices(fold[[slot]])
      }
    }

    fold
  })
}

assert_columns_present <- function(data, columns, label) {
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      sprintf("Missing %s column(s): %s", label, paste(missing_columns, collapse = ", ")),
      call. = FALSE
    )
  }
}

validate_vector_n <- function(x, n, label) {
  if (length(x) != n) {
    stop(sprintf("`%s` must have length %d.", label, n), call. = FALSE)
  }
}

validate_matrix_n <- function(x, n, label) {
  if (nrow(x) != n) {
    stop(sprintf("`%s` must have %d rows.", label, n), call. = FALSE)
  }
}

resolve_treatment_levels <- function(task, treatment_level = NULL, control_level = NULL) {
  treatment <- task$get_tmle_node("treatment")
  levels_found <- levels(factor(treatment))

  if (length(levels_found) < 2L) {
    stop("At least two treatment levels are required.", call. = FALSE)
  }

  if (is.null(control_level)) {
    control_level <- levels_found[[1]]
  }

  if (is.null(treatment_level)) {
    treatment_level <- levels_found[[2]]
  }

  if (!(as.character(treatment_level) %in% as.character(levels_found))) {
    stop("`treatment_level` was not found in the task treatment values.", call. = FALSE)
  }

  if (!(as.character(control_level) %in% as.character(levels_found))) {
    stop("`control_level` was not found in the task treatment values.", call. = FALSE)
  }

  if (identical(as.character(treatment_level), as.character(control_level))) {
    stop("`treatment_level` and `control_level` must be different.", call. = FALSE)
  }

  list(
    treatment_level = treatment_level,
    control_level = control_level,
    levels = levels_found
  )
}

canonicalize_treatment_type <- function(type) {
  if (is.null(type) || length(type) == 0L) {
    return(NA_character_)
  }

  type <- as.character(type[[1]])
  synonyms <- c(
    binary_treatment = "binomial",
    binomial_treatment = "binomial",
    binary = "binomial",
    categorical_treatment = "categorical",
    multinomial_treatment = "categorical",
    continuous_treatment = "continuous"
  )

  if (type %in% names(synonyms)) {
    return(unname(synonyms[[type]]))
  }

  type
}

infer_treatment_type <- function(treatment) {
  observed <- treatment[!is.na(treatment)]
  if (length(observed) == 0L) {
    return("continuous")
  }

  if (is.factor(observed)) {
    return(if (nlevels(observed) == 2L) "binomial" else "categorical")
  }

  unique_values <- unique(observed)
  if (length(unique_values) == 2L) {
    return("binomial")
  }

  if (is.numeric(observed)) {
    return("continuous")
  }

  "categorical"
}

get_task_treatment_type <- function(hte3_task) {
  declared <- canonicalize_treatment_type(hte3_task$npsem$treatment$variable_type$type)
  if (!is.na(declared) && declared %in% c("binomial", "categorical", "continuous")) {
    return(declared)
  }

  infer_treatment_type(hte3_task$get_tmle_node("treatment"))
}

validate_supported_treatment_type <- function(hte3_task, supported_types, learner_name) {
  supported_types <- vapply(supported_types, canonicalize_treatment_type, character(1))
  observed_type <- get_task_treatment_type(hte3_task)

  if (!(observed_type %in% supported_types)) {
    stop(
      sprintf(
        "%s is not compatible with `%s` treatments. Supported types: %s.",
        learner_name,
        observed_type,
        paste(sort(unique(supported_types)), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(observed_type)
}

validate_finite_vector <- function(x, label, lower = -Inf, upper = Inf, allow_na = FALSE) {
  if (!allow_na && anyNA(x)) {
    stop(sprintf("`%s` contains missing values.", label), call. = FALSE)
  }

  finite_values <- x[!is.na(x)]
  if (length(finite_values) > 0L && any(!is.finite(finite_values))) {
    stop(sprintf("`%s` contains non-finite values.", label), call. = FALSE)
  }

  if (length(finite_values) > 0L && any(finite_values < lower | finite_values > upper)) {
    stop(
      sprintf("`%s` must lie in [%s, %s].", label, format(lower), format(upper)),
      call. = FALSE
    )
  }

  invisible(x)
}

get_nuisance_matrix <- function(hte3_task, node, label = node) {
  estimates <- as.matrix(hte3_task$get_nuisance_estimates(node))
  validate_matrix_n(estimates, hte3_task$nrow, label)
  validate_finite_vector(as.vector(estimates), label)

  if (ncol(estimates) == 0L) {
    stop(sprintf("`%s` must have at least one column.", label), call. = FALSE)
  }

  estimates
}

get_nuisance_vector <- function(hte3_task, node, label = node) {
  estimates <- as.vector(hte3_task$get_nuisance_estimates(node))
  validate_vector_n(estimates, hte3_task$nrow, label)
  validate_finite_vector(estimates, label)
  estimates
}

extract_treatment_column <- function(estimates, treatment_level, label) {
  if (is.null(colnames(estimates))) {
    stop(sprintf("`%s` must have named columns for each treatment level.", label), call. = FALSE)
  }

  match_index <- match(as.character(treatment_level), as.character(colnames(estimates)))
  if (is.na(match_index)) {
    stop(
      sprintf("`%s` does not contain estimates for treatment level `%s`.", label, treatment_level),
      call. = FALSE
    )
  }

  as.vector(estimates[, match_index])
}

bound_probability <- function(x, eps = 1e-8) {
  pmin(pmax(x, eps), 1 - eps)
}

bound_positive <- function(x, eps = 1e-8) {
  pmax(x, eps)
}

bounded_qlogis <- function(x, eps = 1e-8) {
  stats::qlogis(bound_probability(x, eps = eps))
}

validate_nonnegative_outcome <- function(y, label = "outcome") {
  if (any(y < 0, na.rm = TRUE)) {
    stop(sprintf("%s must be non-negative.", label), call. = FALSE)
  }
}

resolve_contrast_nuisances <- function(hte3_task, treatment_level = NULL, control_level = NULL, truncation_method = "adaptive") {
  contrast <- resolve_treatment_levels(hte3_task, treatment_level, control_level)
  A <- hte3_task$get_tmle_node("treatment")
  Y <- hte3_task$get_tmle_node("outcome")
  pi.hat <- get_nuisance_matrix(hte3_task, "pi", "pi")
  mu.hat <- get_nuisance_matrix(hte3_task, "mu", "mu")

  pi.hat.1 <- truncate_propensity(
    extract_treatment_column(pi.hat, contrast$treatment_level, "pi"),
    A,
    treatment_level = contrast$treatment_level,
    truncation_method = truncation_method
  )
  pi.hat.0 <- truncate_propensity(
    extract_treatment_column(pi.hat, contrast$control_level, "pi"),
    A,
    treatment_level = contrast$control_level,
    truncation_method = truncation_method
  )

  list(
    A = A,
    Y = Y,
    treatment_level = contrast$treatment_level,
    control_level = contrast$control_level,
    pi.hat.1 = pi.hat.1,
    pi.hat.0 = pi.hat.0,
    mu.hat.1 = extract_treatment_column(mu.hat, contrast$treatment_level, "mu"),
    mu.hat.0 = extract_treatment_column(mu.hat, contrast$control_level, "mu")
  )
}

compute_cate_dr_pseudo_outcome <- function(hte3_task, treatment_level = NULL, control_level = NULL, truncation_method = "adaptive") {
  nuisance <- resolve_contrast_nuisances(
    hte3_task,
    treatment_level = treatment_level,
    control_level = control_level,
    truncation_method = truncation_method
  )

  with(nuisance, {
    mu.hat.1 - mu.hat.0 +
      (A == treatment_level) / pi.hat.1 * (Y - mu.hat.1) -
      (A == control_level) / pi.hat.0 * (Y - mu.hat.0)
  })
}

compute_crr_ratio_pseudo_data <- function(hte3_task, treatment_level = NULL, control_level = NULL, truncation_method = "adaptive") {
  nuisance <- resolve_contrast_nuisances(
    hte3_task,
    treatment_level = treatment_level,
    control_level = control_level,
    truncation_method = truncation_method
  )
  validate_nonnegative_outcome(nuisance$Y, label = "CRR outcomes")

  mu.tilde.1 <- nuisance$mu.hat.1 +
    (nuisance$A == nuisance$treatment_level) / nuisance$pi.hat.1 * (nuisance$Y - nuisance$mu.hat.1)
  mu.tilde.0 <- nuisance$mu.hat.0 +
    (nuisance$A == nuisance$control_level) / nuisance$pi.hat.0 * (nuisance$Y - nuisance$mu.hat.0)

  mu.tilde.1 <- bound_positive(mu.tilde.1)
  mu.tilde.0 <- bound_positive(mu.tilde.0)
  pseudo_weights <- bound_positive(mu.tilde.0 + mu.tilde.1)
  pseudo_outcome <- bound_probability(mu.tilde.1 / pseudo_weights)

  list(
    pseudo_outcome = pseudo_outcome,
    pseudo_weights = pseudo_weights,
    mu.tilde.1 = mu.tilde.1,
    mu.tilde.0 = mu.tilde.0
  )
}

make_tlearner_task <- function(hte3_task) {
  training_task <- hte3_task$next_in_chain(
    covariates = c(hte3_task$npsem$modifiers$variables, hte3_task$npsem$treatment$variables),
    outcome = hte3_task$npsem$outcome$variables
  )

  # Rebuild as an sl3 task so outcome typing does not force binomial models in
  # meta-learners that regress treatment-specific means on modifiers.
  sl3_Task$new(
    training_task$internal_data,
    column_names = training_task$column_names,
    row_index = training_task$row_index,
    nodes = training_task$nodes,
    folds = training_task$folds
  )
}

train_tlearner_models <- function(base_learner, learner_task, hte3_task, treatment_level, control_level, stratify_by_treatment = TRUE) {
  treatment <- hte3_task$get_tmle_node("treatment")
  modifiers <- hte3_task$npsem$modifiers$variables

  if (stratify_by_treatment) {
    learner_task <- learner_task$next_in_chain(covariates = modifiers)
    return(list(
      learner_trained_trt = base_learner$train(learner_task[treatment == treatment_level]),
      learner_trained_control = base_learner$train(learner_task[treatment == control_level])
    ))
  }

  list(learner_trained_pooled = base_learner$train(learner_task))
}

make_counterfactual_treatment_task <- function(hte3_task, treatment_level) {
  treatment_name <- hte3_task$npsem$treatment$variables
  treatment_values <- hte3_task$get_tmle_node("treatment")
  replacement <- rep(treatment_level, hte3_task$nrow)
  if (is.factor(treatment_values)) {
    replacement <- factor(as.character(replacement), levels = levels(treatment_values))
  } else if (is.numeric(treatment_values)) {
    replacement <- as.numeric(replacement)
  }
  counterfactual_data <- data.table(replacement)
  names(counterfactual_data) <- treatment_name
  hte3_task$generate_counterfactual_task(uuid::UUIDgenerate(), counterfactual_data)
}

predict_tlearner_means <- function(fit_object, hte3_task, learner_task, treatment_level, control_level, stratify_by_treatment = TRUE) {
  modifiers <- hte3_task$npsem$modifiers$variables

  if (stratify_by_treatment) {
    learner_task <- learner_task$next_in_chain(covariates = modifiers)
    return(list(
      mu1.hat = fit_object$learner_trained_trt$predict(learner_task),
      mu0.hat = fit_object$learner_trained_control$predict(learner_task)
    ))
  }

  cf0_hte3_task <- make_counterfactual_treatment_task(hte3_task, control_level)
  cf1_hte3_task <- make_counterfactual_treatment_task(hte3_task, treatment_level)

  list(
    mu0.hat = fit_object$learner_trained_pooled$predict(make_tlearner_task(cf0_hte3_task)),
    mu1.hat = fit_object$learner_trained_pooled$predict(make_tlearner_task(cf1_hte3_task))
  )
}

make_sieve_basis <- function(X, basisN = NULL, interaction_order = 3) {
  if (is.null(basisN)) {
    basisN <- ceiling((nrow(X))^(1 / 3) * ncol(X))
  }

  Sieve::sieve_preprocess(
    as.matrix(X),
    basisN = max(1L, as.integer(basisN)) + 1L,
    interaction_order = interaction_order,
    type = "cosine"
  )$Phi
}

make_ep_basis_grid <- function(hte3_task) {
  d <- max(1L, length(hte3_task$npsem$modifiers$variables))
  unique(as.integer(c(d, 2 * d, 4 * d, 6 * d, 8 * d)))
}

make_cate_ep_candidates <- function(base_learner,
                                    hte3_task,
                                    treatment_level,
                                    control_level,
                                    sieve_interaction_order = 3,
                                    screen_basis_with_lasso = FALSE) {
  grid_nbasis <- make_ep_basis_grid(hte3_task)
  lapply(grid_nbasis, function(basis_size) {
    Lrnr_cate_EP$new(
      base_learner = base_learner,
      sieve_num_basis = basis_size,
      sieve_interaction_order = sieve_interaction_order,
      screen_basis_with_lasso = screen_basis_with_lasso,
      treatment_level = treatment_level,
      control_level = control_level
    )
  })
}

make_crr_ep_candidates <- function(base_learner,
                                   hte3_task,
                                   treatment_level,
                                   control_level,
                                   sieve_interaction_order = 3) {
  grid_nbasis <- make_ep_basis_grid(hte3_task)
  lapply(grid_nbasis, function(basis_size) {
    Lrnr_crr_EP$new(
      base_learner = base_learner,
      sieve_num_basis = basis_size,
      sieve_interaction_order = sieve_interaction_order,
      treatment_level = treatment_level,
      control_level = control_level
    )
  })
}

default_cv_flag <- function(base_learner, cross_validate) {
  if (!is.null(cross_validate)) {
    return(isTRUE(cross_validate))
  }

  inherits(base_learner, "Stack")
}

make_crr_ep_stack <- function(base_learner, hte3_task, treatment_level, control_level) {
  do.call(Stack$new, make_crr_ep_candidates(base_learner, hte3_task, treatment_level, control_level))
}

get_autoML <- function() {
  Stack$new(
    Lrnr_glmnet$new(),
    Lrnr_gam$new(),
    Lrnr_earth$new(degree = 2),
    Lrnr_ranger$new(max.depth = 10),
    Lrnr_xgboost$new(min_child_weight = 15, max_depth = 3, nrounds = 40, eta = 0.15, subsample = 0.9),
    Lrnr_xgboost$new(min_child_weight = 15, max_depth = 4, nrounds = 40, eta = 0.15, subsample = 0.9),
    Lrnr_xgboost$new(min_child_weight = 15, max_depth = 5, nrounds = 40, eta = 0.15, subsample = 0.9)
  )
}

make_cross_fitted <- function(learner, calibrate = FALSE, cross_validate = inherits(learner, "Stack")) {
  if (is.null(learner)) {
    return(NULL)
  }

  learner <- Lrnr_cv$new(learner)
  if (cross_validate) {
    learner <- make_learner(Pipeline, learner, Lrnr_cv_selector$new(loss_squared_error))
  }

  learner
}

truncate_propensity <- function(pi.hat, A, treatment_level = max(A), truncation_method = c("isotonic", "adaptive", "deterministic", "none")) {
  truncation_method <- match.arg(truncation_method)
  n <- length(A)
  treatment_indicator <- as.numeric(A == treatment_level)
  pi.hat <- bound_probability(pi.hat)
  pi.hat.star <- pi.hat
  lower_bound <- max(1 / max(n, 2L), 1e-8)
  upper_bound <- min(0.5, 1 - lower_bound)

  if (truncation_method == "isotonic") {
    calibrator_pi <- as.stepfun(isoreg(pi.hat, treatment_indicator))
    pi.hat.star <- calibrator_pi(pi.hat)
  } else if (truncation_method == "adaptive" && lower_bound < upper_bound) {
    risk_function <- function(cutoff) {
      bounded <- pmax(pi.hat, cutoff)
      alpha.hat <- ifelse(treatment_indicator == 1, 1 / bounded, 0)
      alpha1.hat <- 1 / bounded
      mean(alpha.hat^2 - 2 * alpha1.hat)
    }
    cutoff <- optim(1 / sqrt(max(n, 2L)), fn = risk_function, method = "Brent", lower = lower_bound, upper = upper_bound)$par
    pi.hat.star <- pmax(pi.hat.star, cutoff)

    risk_function <- function(cutoff) {
      bounded <- pmin(pi.hat.star, 1 - cutoff)
      alpha.hat <- ifelse(treatment_indicator == 0, 1 / (1 - bounded), 0)
      alpha1.hat <- 1 / (1 - bounded)
      mean(alpha.hat^2 - 2 * alpha1.hat)
    }
    cutoff <- optim(1 / sqrt(max(n, 2L)), fn = risk_function, method = "Brent", lower = lower_bound, upper = upper_bound)$par
    pi.hat.star <- pmin(pi.hat.star, 1 - cutoff)
  } else if (truncation_method == "deterministic") {
    cutoff <- min(upper_bound, max(lower_bound, 25 / sqrt(max(n, 2L)) / log(max(n, 3L))))
    pi.hat.star <- pmax(pi.hat.star, cutoff)
    pi.hat.star <- pmin(pi.hat.star, 1 - cutoff)
  }

  pi.hat.star <- pmax(pi.hat.star, lower_bound)
  pi.hat.star <- pmin(pi.hat.star, 1 - lower_bound)
  pi.hat.star
}

Lrnr_stratified_multivariate <- R6Class(
  classname = "Lrnr_stratified_multivariate", inherit = Lrnr_stratified,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(learner, variable_stratify, ...) {
      super$initialize(learner, variable_stratify, ...)
    }
  ),
  active = list(
    strata_levels = function() {
      variable_stratify <- self$params$variable_stratify
      strata_ids <- self$training_task$data[, variable_stratify, with = FALSE][[1]]
      if (!is.factor(strata_ids)) {
        strata_ids <- factor(strata_ids)
      }
      levels(strata_ids)
    }
  ),
  private = list(
    .train = function(task) {
      learner <- self$params$learner
      learner_trained <- learner$train(task)
      list(learner_trained = learner_trained)
    },
    .predict = function(task = NULL) {
      learner_trained <- self$fit_object$learner_trained
      variable_stratify <- self$params$variable_stratify
      strata_levels <- self$strata_levels

      prediction_mat <- as.matrix(do.call(cbind, lapply(strata_levels, function(strata) {
        if (all(is.numeric(task$data[[variable_stratify]]))) {
          strata <- as.numeric(strata)
        } else {
          strata <- as.character(strata)
        }

        new_data <- data.table::copy(task$data)
        data.table::set(new_data, , variable_stratify, rep(strata, nrow(new_data)))
        new_task <- sl3_Task$new(new_data, nodes = task$nodes)
        learner_trained$predict(new_task)
      })))

      colnames(prediction_mat) <- as.character(strata_levels)
      prediction_mat
    },
    .required_packages = NULL
  )
)

calibrate <- function(predictor, outcome) {
  as.stepfun(isoreg(predictor, outcome))(predictor)
}

estimate_mu <- function(W, A, Y, treatment_level, learner = Lrnr_gam$new(family = "gaussian"), weights = NULL, cross_fit_and_cv = TRUE, stratified_by_trt = TRUE, return_learner = FALSE, folds = 10, ...) {
  W <- as.data.table(W)
  if (length(names(W)) == 0L) {
    names(W) <- paste0("W", seq_len(ncol(W)))
  }

  data <- cbind(W, data.table(A = A, Y = Y, weights = weights))

  if (cross_fit_and_cv) {
    learner <- make_learner(Pipeline, Lrnr_cv$new(learner), Lrnr_cv_selector$new(loss_squared_error()))
  }

  if (stratified_by_trt) {
    task <- sl3_Task$new(data, covariates = names(W), outcome = "Y", weights = weights, folds = folds, ...)
    learner_trained <- learner$train(task[A == treatment_level])
    mu.hat <- learner_trained$predict(task)
  } else {
    task <- sl3_Task$new(data, covariates = c(names(W), "A"), outcome = "Y", weights = weights, folds = folds, ...)
    learner_trained <- learner$train(task)
    new_data <- data.table::copy(data)
    new_data$A <- treatment_level
    mu.hat <- learner_trained$predict(sl3_Task$new(new_data, covariates = c(names(W), "A"), outcome = "Y", weights = weights, folds = folds, ...))
  }

  output <- list(mu.hat = mu.hat, folds = folds)
  if (return_learner) {
    output$learner <- learner_trained
  }

  output
}

estimate_m <- function(W, Y, learner = Lrnr_gam$new(family = "gaussian"), weights = NULL, cross_fit_and_cv = TRUE, return_learner = FALSE, folds = 10, ...) {
  W <- as.data.table(W)
  if (length(names(W)) == 0L) {
    names(W) <- paste0("W", seq_len(ncol(W)))
  }

  data <- cbind(W, data.table(Y = Y, weights = weights))

  if (cross_fit_and_cv) {
    learner <- make_learner(Pipeline, Lrnr_cv$new(learner), Lrnr_cv_selector$new(loss_squared_error()))
  }

  task <- sl3_Task$new(data, covariates = names(W), outcome = "Y", weights = weights, folds = folds, ...)
  learner_trained <- learner$train(task)
  m.hat <- learner_trained$predict(task)

  output <- list(m.hat = m.hat, folds = folds)
  if (return_learner) {
    output$learner <- learner_trained
  }

  output
}

estimate_pi <- function(W, A, binomial_learner = Lrnr_gam$new(family = "binomial"), weights = NULL, treatment_level = max(A), cross_fit_and_cv = TRUE, return_learner = FALSE, folds = 10, ...) {
  A <- factor(A, levels = sort(unique(A)))
  learner <- Lrnr_independent_binomial$new(binomial_learner)

  W <- as.data.table(W)
  if (length(names(W)) == 0L) {
    names(W) <- paste0("W", seq_len(ncol(W)))
  }

  data <- cbind(W, data.table(A = A, weights = weights))

  if (cross_fit_and_cv) {
    learner <- make_learner(Pipeline, Lrnr_cv$new(learner), Lrnr_cv_selector$new(loss_squared_error()))
  }

  task <- sl3_Task$new(data, covariates = names(W), outcome = "A", weights = weights, folds = folds, ...)
  learner_trained <- learner$train(task)
  pi.hat <- coerce_prediction_matrix(learner_trained$predict(task))
  colnames(pi.hat) <- levels(A)

  output <- list(pi.hat = pi.hat, levels = levels(A), folds = folds)
  if (return_learner) {
    output$learner <- learner_trained
  }

  output
}
