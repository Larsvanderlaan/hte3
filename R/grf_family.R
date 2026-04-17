#' GRF Forest Learner
#'
#' A family-aware `sl3` learner backed by the `grf` package.
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data
#' @return Learner object with methods for training and prediction.
Lrnr_grf_forest <- R6Class(
  classname = "Lrnr_grf_forest",
  inherit = Lrnr_base,
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(family = NULL,
                          tune = c("light", "none", "all"),
                          positive_level = "1",
                          ...) {
      params <- c(
        list(
          family = family,
          tune = match.arg(tune),
          positive_level = as.character(positive_level)
        ),
        list(...)
      )
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("continuous", "binomial", "quasibinomial", "weights"),
    .required_packages = c("grf"),
    .train = function(task) {
      validate_grf_available()

      outcome_type <- self$get_outcome_type(task)
      family <- resolve_grf_family(self$params$family, outcome_type$type)
      family_name <- grf_family_name(family)
      X <- grf_feature_matrix(task$X)
      y <- task$Y
      args <- prepare_grf_training_args(self$params)
      positive_level <- args$positive_level
      args$positive_level <- NULL

      if (task$has_node("weights")) {
        args$sample.weights <- task$weights
      }

      forest <- if (identical(family_name, "binomial")) {
        y_factor <- grf_binary_factor(y, label = "task outcome")
        call_with_args(
          grf::probability_forest,
          c(list(X = X, Y = y_factor), args),
          silent = TRUE
        )
      } else if (identical(family_name, "quasibinomial")) {
        y_link <- bounded_qlogis(coerce_numeric_vector_input(y, "task outcome"))
        call_with_args(
          grf::regression_forest,
          c(list(X = X, Y = y_link), args),
          silent = TRUE
        )
      } else {
        y_numeric <- coerce_numeric_vector_input(y, "task outcome")
        call_with_args(
          grf::regression_forest,
          c(list(X = X, Y = y_numeric), args),
          silent = TRUE
        )
      }

      list(
        forest = forest,
        x_names = colnames(X),
        family_name = family_name,
        positive_level = positive_level
      )
    },
    .predict = function(task) {
      fit_object <- private$.fit_object
      new_x <- align_grf_newx(task$X, fit_object$x_names)
      prediction <- stats::predict(fit_object$forest, newdata = new_x)
      values <- prediction$predictions

      if (identical(fit_object$family_name, "binomial")) {
        if (is.null(dim(values))) {
          return(as.numeric(values))
        }
        if (!is.null(colnames(values)) && fit_object$positive_level %in% colnames(values)) {
          return(as.numeric(values[, fit_object$positive_level]))
        }
        return(as.numeric(values[, ncol(values), drop = TRUE]))
      }

      values <- as.numeric(values)
      if (identical(fit_object$family_name, "quasibinomial")) {
        return(bound_probability(stats::plogis(values)))
      }

      values
    }
  )
)

#' GRF Causal Forest Learner
#'
#' A direct CATE learner backed by `grf::causal_forest`.
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data
#' @return Learner object with methods for training and prediction.
Lrnr_grf_causal_forest <- R6Class(
  classname = "Lrnr_grf_causal_forest",
  inherit = Lrnr_base,
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(treatment_level = NULL,
                          control_level = NULL,
                          modifiers = NULL,
                          tune = c("light", "none", "all"),
                          ...) {
      params <- c(
        list(
          treatment_level = treatment_level,
          control_level = control_level,
          modifiers = modifiers,
          tune = match.arg(tune)
        ),
        list(...)
      )
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("continuous", "binomial", "weights"),
    .required_packages = c("grf", "sl3"),
    .train = function(task) {
      validate_grf_available()
      validate_hte_task(task)
      validate_supported_treatment_type(task, "binomial", "Lrnr_grf_causal_forest")

      modifiers <- private$.resolve_modifiers(task)
      contrast <- resolve_treatment_levels(
        task,
        treatment_level = self$params$treatment_level,
        control_level = self$params$control_level
      )
      nuisance <- resolve_contrast_nuisances(
        task,
        treatment_level = contrast$treatment_level,
        control_level = contrast$control_level
      )
      m.hat <- get_nuisance_vector(task, "m", "m")
      treatment <- task$get_tmle_node("treatment")
      outcome <- coerce_numeric_vector_input(task$get_tmle_node("outcome"), "outcome")
      treatment_indicator <- as.numeric(as.character(treatment) == as.character(contrast$treatment_level))
      x_data <- task$data[, modifiers, with = FALSE]
      X <- grf_feature_matrix(x_data)
      args <- prepare_grf_training_args(self$params)
      args$treatment_level <- NULL
      args$control_level <- NULL
      args$modifiers <- NULL
      args$X <- X
      args$Y <- outcome
      args$W <- treatment_indicator
      args$Y.hat <- m.hat
      args$W.hat <- nuisance$pi.hat.1

      if (!is.null(task$weights)) {
        args$sample.weights <- task$weights
      }

      forest <- call_with_args(grf::causal_forest, args, silent = TRUE)
      list(
        forest = forest,
        x_names = colnames(X),
        prediction_template = make_prediction_template(task),
        contrast = contrast
      )
    },
    .predict = function(task) {
      validate_hte_task(task)
      fit_object <- private$.fit_object
      modifiers <- private$.resolve_modifiers(task)
      new_x <- align_grf_newx(task$data[, modifiers, with = FALSE], fit_object$x_names)
      prediction <- stats::predict(fit_object$forest, newdata = new_x)
      as.numeric(prediction$predictions)
    },
    .resolve_modifiers = function(task) {
      if (!is.null(self$params$modifiers)) {
        modifiers <- as.character(self$params$modifiers)
      } else {
        modifiers <- task$npsem$modifiers$variables
      }

      assert_columns_present(task$data, modifiers, "modifier")
      modifiers
    }
  )
)

validate_grf_available <- function() {
  if (!requireNamespace("grf", quietly = TRUE)) {
    stop(
      paste(
        "The GRF workflow requires the optional `grf` package.",
        "Install it with `install.packages(\"grf\")` and try again."
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

grf_family_name <- function(family) {
  if (is.null(family)) {
    return("gaussian")
  }

  if (is.character(family)) {
    family_name <- as.character(family[[1L]])
  } else if (inherits(family, "family")) {
    family_name <- family$family
  } else {
    stop("`family` must be NULL, a character string, or a family object.", call. = FALSE)
  }

  if (!(family_name %in% c("gaussian", "binomial", "quasibinomial"))) {
    stop("GRF learners support only gaussian, binomial, and quasibinomial families.", call. = FALSE)
  }

  family_name
}

resolve_grf_family <- function(family, outcome_type = NULL) {
  if (!is.null(family)) {
    family_name <- grf_family_name(family)
  } else if (identical(outcome_type, "binomial")) {
    family_name <- "binomial"
  } else if (identical(outcome_type, "quasibinomial")) {
    family_name <- "quasibinomial"
  } else {
    family_name <- "gaussian"
  }

  switch(
    family_name,
    gaussian = stats::gaussian(),
    binomial = stats::binomial(),
    quasibinomial = stats::quasibinomial()
  )
}

prepare_grf_training_args <- function(params) {
  args <- params
  tune <- args$tune
  args$tune <- NULL
  args$family <- NULL

  if (identical(tune, "light")) {
    if (is.null(args$tune.parameters)) {
      args$tune.parameters <- c("sample.fraction", "mtry", "min.node.size")
    }
    if (is.null(args$tune.num.trees)) {
      args$tune.num.trees <- 200
    }
    if (is.null(args$tune.num.reps)) {
      args$tune.num.reps <- 20
    }
    if (is.null(args$tune.num.draws)) {
      args$tune.num.draws <- 100
    }
  } else if (identical(tune, "all")) {
    if (is.null(args$tune.parameters)) {
      args$tune.parameters <- "all"
    }
  } else if (is.null(args$tune.parameters)) {
    args$tune.parameters <- "none"
  }

  args
}

grf_feature_matrix <- function(x) {
  matrix_x <- data.matrix(as.data.frame(x))
  if (is.null(colnames(matrix_x)) && ncol(matrix_x) > 0L) {
    colnames(matrix_x) <- names(x)
  }
  matrix_x
}

align_grf_newx <- function(x, train_names) {
  new_x <- as.data.table(x)

  missing_names <- setdiff(train_names, names(new_x))
  if (length(missing_names) > 0L) {
    stop(
      sprintf(
        "Prediction task covariates missing covariates in training task: %s",
        paste(missing_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  extra_names <- setdiff(names(new_x), train_names)
  if (length(extra_names) > 0L) {
    warning(
      sprintf(
        "Prediction task has new covariates not seen in training: %s; dropping these for prediction",
        paste(extra_names, collapse = ", ")
      ),
      call. = FALSE
    )
    new_x <- new_x[, train_names, with = FALSE]
  } else {
    new_x <- new_x[, train_names, with = FALSE]
  }

  grf_feature_matrix(new_x)
}

grf_binary_factor <- function(y, label = "outcome") {
  y <- coerce_numeric_vector_input(y, label)
  validate_finite_vector(y, label, lower = 0, upper = 1)

  y_binary <- round(y)
  if (any(abs(y - y_binary) > sqrt(.Machine$double.eps))) {
    stop(sprintf("`%s` must contain only 0/1 values for binomial GRF fits.", label), call. = FALSE)
  }

  factor(as.character(y_binary), levels = c("0", "1"))
}

infer_grf_outcome_family <- function(y) {
  observed <- y[!is.na(y)]
  if (length(observed) > 0L && all(observed %in% c(0, 1))) {
    return(stats::binomial())
  }

  stats::gaussian()
}

normalize_binary_mu_hat <- function(mu.hat, n, outcome_type, contrast) {
  mu.hat <- validate_outcome_matrix_input(mu.hat, n, "mu.hat", outcome_type = outcome_type)
  if (ncol(mu.hat) != 2L) {
    stop("`mu.hat` must be an n x 2 matrix ordered as (control, treatment).", call. = FALSE)
  }

  colnames(mu.hat) <- c(as.character(contrast$control_level), as.character(contrast$treatment_level))
  mu.hat
}

expand_binary_pi_hat <- function(pi.hat, n, contrast) {
  pi.hat <- coerce_numeric_vector_input(pi.hat, "pi.hat")
  validate_vector_n(pi.hat, n, "pi.hat")
  validate_finite_vector(pi.hat, "pi.hat", lower = 0, upper = 1)

  out <- cbind(
    1 - pi.hat,
    pi.hat
  )
  colnames(out) <- c(as.character(contrast$control_level), as.character(contrast$treatment_level))
  out
}

normalize_binary_m_hat <- function(m.hat, n, outcome_type) {
  validate_outcome_vector_input(m.hat, n, "m.hat", outcome_type = outcome_type)
}

compute_binary_m_hat <- function(mu.hat, pi.hat) {
  as.vector(rowSums(mu.hat * pi.hat))
}

resolve_binary_data_contrast <- function(data, treatment, treatment_level = 1, control_level = 0) {
  treatment_values <- data[[treatment]]
  treatment_levels <- levels(factor(treatment_values))

  if (length(treatment_levels) != 2L) {
    stop("GRF wrappers currently support only binary treatments.", call. = FALSE)
  }

  if (!(as.character(treatment_level) %in% as.character(treatment_levels)) &&
      identical(treatment_level, 1) &&
      length(treatment_levels) == 2L) {
    treatment_level <- treatment_levels[[2L]]
  }

  if (!(as.character(control_level) %in% as.character(treatment_levels)) &&
      identical(control_level, 0) &&
      length(treatment_levels) == 2L) {
    control_level <- treatment_levels[[1L]]
  }

  if (!(as.character(treatment_level) %in% as.character(treatment_levels))) {
    stop("`treatment_level` was not found in the treatment values.", call. = FALSE)
  }

  if (!(as.character(control_level) %in% as.character(treatment_levels))) {
    stop("`control_level` was not found in the treatment values.", call. = FALSE)
  }

  if (identical(as.character(treatment_level), as.character(control_level))) {
    stop("`treatment_level` and `control_level` must be different.", call. = FALSE)
  }

  list(
    treatment_level = treatment_level,
    control_level = control_level,
    levels = treatment_levels
  )
}

make_grf_base_learner <- function(tune = c("light", "none", "all"), grf_params = list(), family = NULL) {
  tune <- match.arg(tune)
  if (!is.list(grf_params)) {
    stop("`grf_params` must be a named list of GRF arguments.", call. = FALSE)
  }
  do.call(Lrnr_grf_forest$new, c(list(family = family, tune = tune), grf_params))
}

#' Build Default GRF Nuisance Learners
#'
#' @param outcome_family Optional family object for the outcome and mean models.
#' @param grf_params Optional named list of GRF arguments passed to the learners.
#'
#' @return A named list with propensity, outcome, and mean learners.
#' @export
make_grf_nuisance_learners <- function(outcome_family = NULL, grf_params = list()) {
  if (is.null(outcome_family)) {
    outcome_family <- stats::gaussian()
  }

  list(
    propensity_learner = make_grf_base_learner(tune = "none", grf_params = grf_params, family = stats::binomial()),
    outcome_learner = make_grf_base_learner(tune = "none", grf_params = grf_params, family = outcome_family),
    mean_learner = make_grf_base_learner(tune = "none", grf_params = grf_params, family = outcome_family)
  )
}

#' Build GRF-Backed CATE Learners
#'
#' @param method Character vector of GRF-supported CATE methods.
#' @param hte3_task Optional task used when expanding EP basis-size portfolios.
#' @param treatment_level Treated level for the target contrast.
#' @param control_level Control level for the target contrast.
#' @param tune Tuning mode for the final GRF learner(s).
#' @param grf_params Optional named list of GRF arguments passed to the learners.
#' @param sieve_num_basis Optional EP basis size for single-fit EP learners.
#' @param sieve_basis_grid Optional EP basis grid.
#' @param sieve_interaction_order EP interaction order.
#' @param screen_basis_with_lasso Whether EP basis screening is enabled.
#' @param ep_targeting_style EP targeting variant used when `method` includes `"ep"`.
#' @param ep_r_targeting_basis First-stage basis construction used for EP-R.
#'
#' @return A list of learner objects.
#' @export
make_grf_cate_learners <- function(method = c("ep", "r"),
                                   hte3_task = NULL,
                                   treatment_level = 1,
                                   control_level = 0,
                                   tune = c("light", "none", "all"),
                                   grf_params = list(),
                                   sieve_num_basis = NULL,
                                   sieve_basis_grid = NULL,
                                   sieve_interaction_order = 3,
                                   screen_basis_with_lasso = FALSE,
                                   ep_targeting_style = "r",
                                   ep_r_targeting_basis = "v_plus_propensity") {
  tune <- match.arg(tune)
  methods <- normalize_method_spec(method, supported_grf_cate_methods())
  ep_targeting_styles <- if ("ep" %in% methods) normalize_ep_targeting_styles(ep_targeting_style) else ep_targeting_style
  if ("ep" %in% methods && "r" %in% ep_targeting_styles && isTRUE(screen_basis_with_lasso)) {
    stop("`screen_basis_with_lasso = TRUE` is only supported when `ep_targeting_style` does not include \"r\".", call. = FALSE)
  }
  base_learner <- make_grf_base_learner(tune = tune, grf_params = grf_params)
  learners <- list()

  if ("dr" %in% methods) {
    learners <- c(learners, list(Lrnr_cate_DR$new(
      base_learner = base_learner,
      treatment_level = treatment_level,
      control_level = control_level
    )))
  }

  if ("r" %in% methods) {
    learners <- c(learners, list(do.call(
      Lrnr_grf_causal_forest$new,
      c(
        list(
          treatment_level = treatment_level,
          control_level = control_level,
          tune = tune
        ),
        grf_params
      )
    )))
  }

  if ("ep" %in% methods) {
    if (!is.null(sieve_basis_grid)) {
      if (is.null(hte3_task)) {
        stop("`hte3_task` must be supplied when `sieve_basis_grid` is used.", call. = FALSE)
      }
      learners <- c(
        learners,
        make_cate_ep_candidates(
          base_learner,
          hte3_task,
          treatment_level,
          control_level,
          sieve_basis_grid = sieve_basis_grid,
          sieve_interaction_order = sieve_interaction_order,
          screen_basis_with_lasso = screen_basis_with_lasso,
          targeting_style = ep_targeting_styles,
          r_targeting_basis = ep_r_targeting_basis
        )
      )
    } else {
      learners <- c(learners, lapply(ep_targeting_styles, function(style) {
        Lrnr_cate_EP$new(
          base_learner = base_learner,
          sieve_num_basis = sieve_num_basis,
          sieve_interaction_order = sieve_interaction_order,
          screen_basis_with_lasso = screen_basis_with_lasso,
          targeting_style = style,
          r_targeting_basis = ep_r_targeting_basis,
          treatment_level = treatment_level,
          control_level = control_level
        )
      }))
    }
  }

  learners
}

#' Build GRF-Backed CRR Learners
#'
#' @param method Character vector of GRF-supported CRR methods.
#' @param hte3_task Optional task used when expanding EP basis-size portfolios.
#' @param treatment_level Treated level for the target contrast.
#' @param control_level Control level for the target contrast.
#' @param tune Tuning mode for the final GRF learner(s).
#' @param grf_params Optional named list of GRF arguments passed to the learners.
#' @param sieve_num_basis Optional EP basis size for single-fit EP learners.
#' @param sieve_basis_grid Optional EP basis grid.
#' @param sieve_interaction_order EP interaction order.
#'
#' @return A list of learner objects.
#' @export
make_grf_crr_learners <- function(method = c("ep"),
                                  hte3_task = NULL,
                                  treatment_level = 1,
                                  control_level = 0,
                                  tune = c("light", "none", "all"),
                                  grf_params = list(),
                                  sieve_num_basis = NULL,
                                  sieve_basis_grid = NULL,
                                  sieve_interaction_order = 3) {
  tune <- match.arg(tune)
  methods <- normalize_method_spec(method, supported_grf_crr_methods())
  base_learner <- make_grf_base_learner(tune = tune, grf_params = grf_params)
  learners <- list()

  if ("ipw" %in% methods) {
    learners <- c(learners, list(Lrnr_crr_IPW$new(
      base_learner = base_learner,
      treatment_level = treatment_level,
      control_level = control_level
    )))
  }

  if ("t" %in% methods) {
    learners <- c(learners, list(Lrnr_crr_T$new(
      base_learner = base_learner,
      treatment_level = treatment_level,
      control_level = control_level
    )))
  }

  if ("ep" %in% methods) {
    if (!is.null(sieve_basis_grid)) {
      if (is.null(hte3_task)) {
        stop("`hte3_task` must be supplied when `sieve_basis_grid` is used.", call. = FALSE)
      }
      learners <- c(
        learners,
        make_crr_ep_candidates(
          base_learner,
          hte3_task,
          treatment_level,
          control_level,
          sieve_basis_grid = sieve_basis_grid,
          sieve_interaction_order = sieve_interaction_order
        )
      )
    } else {
      learners <- c(learners, list(Lrnr_crr_EP$new(
        base_learner = base_learner,
        sieve_num_basis = sieve_num_basis,
        sieve_interaction_order = sieve_interaction_order,
        treatment_level = treatment_level,
        control_level = control_level
      )))
    }
  }

  learners
}

resolve_grf_cate_defaults <- function(method = NULL, cross_validate = NULL) {
  if (is.null(method)) {
    if (isTRUE(cross_validate) || is.null(cross_validate)) {
      return(list(method = c("ep", "r"), cross_validate = TRUE))
    }
    return(list(method = "ep", cross_validate = FALSE))
  }

  list(method = method, cross_validate = cross_validate)
}

prepare_grf_binary_nuisances <- function(data,
                                         treatment,
                                         outcome,
                                         treatment_level = 1,
                                         control_level = 0,
                                         mu.hat = NULL,
                                         pi.hat = NULL,
                                         m.hat = NULL) {
  contrast <- resolve_binary_data_contrast(
    data,
    treatment = treatment,
    treatment_level = treatment_level,
    control_level = control_level
  )
  outcome_family <- infer_grf_outcome_family(data[[outcome]])
  outcome_type <- grf_family_name(outcome_family)
  normalized_mu <- NULL
  normalized_pi <- NULL
  normalized_m <- NULL

  if (!is.null(mu.hat)) {
    normalized_mu <- normalize_binary_mu_hat(mu.hat, nrow(data), outcome_type = outcome_type, contrast = contrast)
  }

  if (!is.null(pi.hat)) {
    normalized_pi <- expand_binary_pi_hat(pi.hat, nrow(data), contrast = contrast)
  }

  if (!is.null(m.hat)) {
    normalized_m <- normalize_binary_m_hat(m.hat, nrow(data), outcome_type = outcome_type)
  } else if (!is.null(normalized_mu) && !is.null(normalized_pi)) {
    normalized_m <- compute_binary_m_hat(normalized_mu, normalized_pi)
  }

  list(
    contrast = contrast,
    outcome_family = outcome_family,
    mu.hat = normalized_mu,
    pi.hat = normalized_pi,
    m.hat = normalized_m
  )
}

#' Fit GRF-Backed CATE Models
#'
#' @param data A data frame or data.table.
#' @param modifiers Effect modifiers.
#' @param confounders Adjustment covariates.
#' @param treatment Treatment column.
#' @param outcome Outcome column.
#' @param method Optional GRF CATE method specification.
#' @param cross_validate Optional logical flag for outer learner selection.
#' @param mu.hat Optional `n x 2` matrix of nuisance outcome regressions ordered as `(control, treatment)`.
#' @param pi.hat Optional length-`n` vector of treatment propensities `P(A = treatment_level | X)`.
#' @param m.hat Optional length-`n` vector of marginal outcome means.
#' @param treatment_level Treated level for the target contrast.
#' @param control_level Control level for the target contrast.
#' @param cross_fit Whether nuisance learners should be cross-fitted.
#' @param folds Number of folds used for nuisance cross-fitting.
#' @param cv_control Optional outer cross-validation control list.
#' @param ep_targeting_style EP targeting variant used when `method` includes `"ep"`.
#' @param ep_r_targeting_basis First-stage basis construction used for EP-R.
#' @param tune Tuning mode for the final GRF learner(s).
#' @param grf_params Optional named list of GRF arguments passed to the learners.
#' @param ... Additional arguments passed to `hte_task()`.
#'
#' @return An `hte3_model`.
#' @export
grf_cate <- function(data,
                     modifiers,
                     confounders = modifiers,
                     treatment,
                     outcome,
                     ...,
                     method = NULL,
                     cross_validate = NULL,
                     mu.hat = NULL,
                     pi.hat = NULL,
                     m.hat = NULL,
                     treatment_level = 1,
                     control_level = 0,
                     cross_fit = TRUE,
                     folds = 10,
                     cv_control = NULL,
                     ep_targeting_style = "r",
                     ep_r_targeting_basis = "v_plus_propensity",
                     tune = c("light", "none", "all"),
                     grf_params = list()) {
  tune <- match.arg(tune)
  data <- as.data.table(data)
  nuisance_spec <- prepare_grf_binary_nuisances(
    data = data,
    treatment = treatment,
    outcome = outcome,
    treatment_level = treatment_level,
    control_level = control_level,
    mu.hat = mu.hat,
    pi.hat = pi.hat,
    m.hat = m.hat
  )
  validate_grf_available()
  default_spec <- resolve_grf_cate_defaults(method = method, cross_validate = cross_validate)
  methods <- normalize_method_spec(default_spec$method, supported_grf_cate_methods())
  ep_targeting_styles <- if ("ep" %in% methods) normalize_ep_targeting_styles(ep_targeting_style) else ep_targeting_style
  nuisance_learners <- make_grf_nuisance_learners(
    outcome_family = nuisance_spec$outcome_family,
    grf_params = grf_params
  )

  task <- suppress_internal_cv_fallback_warnings(
    hte_task(
      data = data,
      modifiers = modifiers,
      confounders = confounders,
      treatment = treatment,
      outcome = outcome,
      treatment_type = "default",
      propensity = nuisance_spec$pi.hat,
      outcome_regression = nuisance_spec$mu.hat,
      outcome_mean = nuisance_spec$m.hat,
      propensity_learner = nuisance_learners$propensity_learner,
      outcome_learner = nuisance_learners$outcome_learner,
      mean_learner = nuisance_learners$mean_learner,
      cross_fit = cross_fit,
      folds = folds,
      ...
    )
  )

  if (!identical(get_task_treatment_type(task), "binomial")) {
    stop("GRF wrappers currently support only binary treatments.", call. = FALSE)
  }

  cross_validate_input <- default_spec$cross_validate
  candidate_count <- count_cate_wrapper_candidates(methods, ep_targeting_styles)
  if (candidate_count == 1L) {
    cross_validate_input <- FALSE
  }

  cross_validate <- prepare_cv_flag(
    make_grf_base_learner(tune = tune, grf_params = grf_params),
    cross_validate_input,
    candidate_count = candidate_count
  )

  learners <- make_grf_cate_learners(
    method = methods,
    treatment_level = nuisance_spec$contrast$treatment_level,
    control_level = nuisance_spec$contrast$control_level,
    tune = tune,
    grf_params = grf_params,
    ep_targeting_style = ep_targeting_styles,
    ep_r_targeting_basis = ep_r_targeting_basis
  )

  trained <- suppress_internal_cv_fallback_warnings(
    if (cross_validate) {
      cv_fit <- cross_validate_cate(
        learners,
        task,
        cv_control = cv_control,
        treatment_level = nuisance_spec$contrast$treatment_level,
        control_level = nuisance_spec$contrast$control_level
      )
      selection_summary <- cv_fit$selection_summary
      cv_fit$lrnr_sl
    } else {
      selection_summary <- NULL
      learners[[1L]]$train(task)
    }
  )

  new_hte3_model(
    "cate",
    format_method_label(methods),
    trained,
    task,
    cross_validate,
    match.call(),
    selection_summary = selection_summary,
    engine = "grf"
  )
}

#' Fit GRF-Backed CRR Models
#'
#' @inheritParams grf_cate
#' @return An `hte3_model`.
#' @export
grf_crr <- function(data,
                    modifiers,
                    confounders = modifiers,
                    treatment,
                    outcome,
                    ...,
                    method = NULL,
                    cross_validate = NULL,
                    mu.hat = NULL,
                    pi.hat = NULL,
                    m.hat = NULL,
                    treatment_level = 1,
                    control_level = 0,
                    cross_fit = TRUE,
                    folds = 10,
                    cv_control = NULL,
                    tune = c("light", "none", "all"),
                    grf_params = list()) {
  tune <- match.arg(tune)
  data <- as.data.table(data)
  nuisance_spec <- prepare_grf_binary_nuisances(
    data = data,
    treatment = treatment,
    outcome = outcome,
    treatment_level = treatment_level,
    control_level = control_level,
    mu.hat = mu.hat,
    pi.hat = pi.hat,
    m.hat = m.hat
  )
  validate_grf_available()
  methods <- if (is.null(method)) "ep" else normalize_method_spec(method, supported_grf_crr_methods())
  nuisance_learners <- make_grf_nuisance_learners(
    outcome_family = nuisance_spec$outcome_family,
    grf_params = grf_params
  )

  task <- suppress_internal_cv_fallback_warnings(
    hte_task(
      data = data,
      modifiers = modifiers,
      confounders = confounders,
      treatment = treatment,
      outcome = outcome,
      treatment_type = "default",
      propensity = nuisance_spec$pi.hat,
      outcome_regression = nuisance_spec$mu.hat,
      outcome_mean = nuisance_spec$m.hat,
      propensity_learner = nuisance_learners$propensity_learner,
      outcome_learner = nuisance_learners$outcome_learner,
      mean_learner = nuisance_learners$mean_learner,
      cross_fit = cross_fit,
      folds = folds,
      ...
    )
  )

  if (!identical(get_task_treatment_type(task), "binomial")) {
    stop("GRF wrappers currently support only binary treatments.", call. = FALSE)
  }

  validate_nonnegative_outcome(task$get_tmle_node("outcome"), label = "CRR outcomes")

  cross_validate_input <- cross_validate
  if (length(methods) == 1L) {
    cross_validate_input <- FALSE
  }

  cross_validate <- prepare_cv_flag(
    make_grf_base_learner(tune = tune, grf_params = grf_params),
    cross_validate_input,
    candidate_count = length(methods)
  )

  learners <- make_grf_crr_learners(
    method = methods,
    treatment_level = nuisance_spec$contrast$treatment_level,
    control_level = nuisance_spec$contrast$control_level,
    tune = tune,
    grf_params = grf_params
  )

  trained <- suppress_internal_cv_fallback_warnings(
    if (cross_validate) {
      cv_fit <- cross_validate_crr(
        learners,
        task,
        cv_control = cv_control,
        treatment_level = nuisance_spec$contrast$treatment_level,
        control_level = nuisance_spec$contrast$control_level
      )
      selection_summary <- cv_fit$selection_summary
      cv_fit$lrnr_sl
    } else {
      selection_summary <- NULL
      learners[[1L]]$train(task)
    }
  )

  new_hte3_model(
    "crr",
    format_method_label(methods),
    trained,
    task,
    cross_validate,
    match.call(),
    selection_summary = selection_summary,
    engine = "grf"
  )
}
