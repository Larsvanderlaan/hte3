#'  Task object for meta-learners in the point-treatment setting.
#'
#' Constructs a \code{hte3_Task} object for meta-learners in the point-treatment setting containing relevant data and nuisance function estimators.
#'
#' @param data A named data frame or data.table containing treatment effect modifiers, potential confounders, treatment, and outcome.
#'   Optionally, the dataset can contain weights and subject IDs. See the \code{data} argument of \code{\link{sl3_Task}} for further options.
#'
#' @param modifiers A character vector of variable names in \code{colnames(data)} for treatment effect moderators.
#'
#' @param confounders A character vector of variable names in \code{colnames(data)} for potential confounders \code{W} for which to adjust.
#'
#' @param treatment A character specifying the variable name in \code{colnames(data)} for the numeric treatment assignment \code{A}.
#'
#' @param outcome A character specifying the variable name in \code{colnames(data)} for the outcome variable \code{Y}.
#'
#' @param id An (optional) character specifying the variable name in \code{colnames(data)} for observation IDs.
#'
#' @param weights An (optional) character specifying the variable name in \code{colnames(data)} for observation weights.
#'
#' @param pi.hat An (optional) numeric matrix of dimension \code{n by nlevels(A)} containing estimates of the propensity score \code{a -> pi(a | W_i)}
#'   at each treatment level, where column \code{j} corresponds to treatment level \code{sort(unique(A))[j]}.
#'   This argument can be used by DR-type learners. Alternatively, the estimates can be learned internally by passing an \code{sl3_Learner} to the \code{learner_pi} argument.
#'
#' @param mu.hat An (optional) numeric matrix of dimension \code{n by nlevels(A)} containing estimates of the outcome regression \code{a -> mu(a, W_i)}
#'   at each treatment level, where column \code{j} corresponds to treatment level \code{sort(unique(A))[j]}.
#'   This argument can be used by DR-type learners. Alternatively, the estimates can be learned internally by passing an \code{sl3_Learner} to the \code{learner_m} argument.
#'
#' @param m.hat An (optional) numeric vector of size \code{n} containing estimates of the conditional mean outcome \code{E[Y | W_i]}.
#'   This argument can be used by R-type learners. Alternatively, the estimates can be learned internally by passing an \code{sl3_Learner} to the \code{learner_m} argument.
#'
#' @param learner_pi A binomial \code{sl3_Learner} or \code{Stack} object specifying the learning algorithm to estimate the propensity score \code{w -> pi(a_star | w)}
#'   at a given treatment level \code{a_star}. During training, this learner is passed an \code{sl3_Task} object \code{task} that contains a feature matrix \code{task$X}
#'   of \code{confounders} and an outcome vector \code{task$Y} corresponding to the binary treatment indicator \code{1(A = a_star)}.
#'   By default, \code{cross_fit_and_cv = TRUE}, and \code{learner_m} is fit with 10-fold cross-fitting using \code{\link{make_cross_fitted}}.
#'   If a \code{Stack} object, the best model is selected using cross-validation.
#'
#' @param learner_mu A \code{sl3_Learner} or \code{sl3_Task} object specifying the learning algorithm to estimate the outcome regression \code{w -> mu(a_star, w)}
#'   at a given treatment level \code{a_star}. The outcome regression is estimated via a treatment-stratified regression.
#'   During training, this learner is passed an \code{sl3_Task} object \code{task} that is subsetted to variables with treatment \code{A = a_star},
#'   and contains a feature matrix \code{task$X} of \code{confounders} and an outcome vector \code{task$Y} corresponding to \code{outcome}.
#'   By default, \code{cross_fit_and_cv = TRUE}, and \code{learner_m} is fit with 10-fold cross-fitting using \code{\link{make_cross_fitted}}.
#'   If a \code{Stack} object, the best model is selected using cross-validation.
#'
#' @param learner_m A \code{sl3_Learner} or \code{sl3_Task} object specifying the learning algorithm to estimate the conditional mean outcome \code{m}.
#'   During training, this learner is passed an \code{sl3_Task} object \code{task} that contains a feature matrix \code{task$X} of \code{confounders}
#'   and an outcome vector \code{task$Y} corresponding to \code{outcome}.
#'   By default, \code{cross_fit_and_cv = TRUE}, and \code{learner_m} is fit with 10-fold cross-fitting using \code{\link{make_cross_fitted}}.
#'   If a \code{Stack} object, the best model is selected using cross-validation.
#' @param multinomial_learner A \code{multinomial} \code{Lrnr_base} object to convert \code{learner_pi} from a \code{binomial} learner to a \code{multinomial} learner.
#' @param cross_fit_and_cv Whether to cross-fit the specified nuisance learners by applying \code{learner <- make_cross_fitted(learner)}.
#' @param for_prediction A \code{boolean} of whether to return an hte3_Task without the nuisance training and estimates. This can be useful when you wish to construct a task for predicting on a new dataset.
#' @param ... Additional arguments to pass to the \code{\link[sl3]{sl3_Task}} and \code{\link[tmle3]{tmle_Task}} constructors.
#' @return A \code{hte3_Task} object for the point-treatment data-structure
#'
#' @export
make_hte3_Task_tx <- function(data,
                              modifiers, confounders, treatment, outcome,
                              id = NULL, weights = NULL,
                              treatment_type = c("default", "binomial", "categorical", "continuous"),
                              pi.hat = NULL, mu.hat = NULL,
                              m.hat = NULL,
                              learner_pi = get_autoML(), learner_mu = get_autoML(),
                              learner_m = NULL,
                              multinomial_learner = Lrnr_independent_binomial,
                              cross_fit_and_cv = TRUE,
                              folds = 10,
                              warn = TRUE,
                              for_prediction = FALSE,
                              ...) {
  treatment_type <- match.arg(treatment_type)
  data <- as.data.table(data)
  assert_columns_present(data, modifiers, "modifier")
  assert_columns_present(data, confounders, "confounder")
  assert_columns_present(data, c(treatment, outcome), "required")

  if (!is.numeric(folds) || length(folds) != 1L || folds < 2 || folds != as.integer(folds)) {
    stop("`folds` must be a single integer greater than or equal to 2.", call. = FALSE)
  }

  if (!is.null(id)) {
    assert_columns_present(data, id, "id")
  }

  if (!is.null(weights)) {
    assert_columns_present(data, weights, "weight")
  }

  if (anyNA(data[[outcome]])) {
    stop("Missing outcomes are not currently supported in `make_hte3_Task_tx()`.", call. = FALSE)
  }

  if (anyNA(data[[treatment]])) {
    stop("Missing treatment values are not currently supported in `make_hte3_Task_tx()`.", call. = FALSE)
  }

  if (!is.null(weights)) {
    data[[weights]] <- coerce_numeric_vector_input(data[[weights]], "weights")
    validate_finite_vector(data[[weights]], "weights", lower = 0)
  }

  if (!is.null(id) && anyNA(data[[id]])) {
    stop("Missing IDs are not currently supported in `make_hte3_Task_tx()`.", call. = FALSE)
  }

  if (warn && !is.factor(data[[treatment]]) && treatment_type == "categorical") {
    warning("Treatment variable in `data` is specified as categorical but is not a factor. Forcefully converting treatment to a factor variable. The  control (or reference) level for the treatment is taken as the first element in levels(factor(data[[treatment]])).")
    data[[treatment]] <- factor(data[[treatment]])
  }

  treatment_levels <- levels(factor(data[[treatment]]))

  if (is.null(weights)) {
    weights <- paste0("weights_", uuid::UUIDgenerate())
    data[[weights]] <- 1
  }

  if (is.null(id)) {
    id <- paste0("id_", uuid::UUIDgenerate())
    data[[id]] <- 1:nrow(data)
  }

  node_list <- list(confounders = confounders, treatment = treatment, outcome = outcome, modifiers = modifiers)
  if (treatment_type == "default") {
    var_type_trt <- sl3::variable_type(type = NULL, x = data[[treatment]])
  } else {
    var_type_trt <- sl3::variable_type(type = treatment_type, levels = treatment_levels)
  }
  var_type_modifiers <- sl3::variable_type(x = data[, modifiers, with = FALSE])
  var_type_confounders <- sl3::variable_type(x = data[, confounders, with = FALSE])
  var_type_outcome <- sl3::variable_type(x = data[[outcome]])

  npsem <- list(
    define_node("modifiers", node_list$modifiers, variable_type = var_type_modifiers),
    define_node("confounders", node_list$confounders, variable_type = var_type_confounders),
    define_node("treatment", node_list$treatment, c("confounders"), variable_type = var_type_trt),
    define_node("outcome", node_list$outcome, c("treatment", "confounders"), variable_type = var_type_outcome),
    define_node("pi", node_list$treatment, c("confounders"), variable_type = sl3::variable_type(type = "categorical", levels = treatment_levels)),
    define_node("mu", node_list$outcome, c("treatment", "confounders"), variable_type = var_type_outcome),
    define_node("m", node_list$outcome, c("confounders"), variable_type = var_type_outcome)
  )

  if (for_prediction) {
    return(hte3_Task$new(data, npsem, likelihood = Likelihood$new(list()), id = id, weights = weights, folds = folds, ...))
  }

  if (!is.null(pi.hat)) {
    pi.hat <- validate_probability_matrix_input(pi.hat, nrow(data), "pi.hat")
    if (ncol(pi.hat) != length(treatment_levels)) {
      stop("pi.hat should be a matrix with ncol(pi.hat) = length(treatment_levels). The jth column should be estimated for the propensity score of the jth treatment level in treatment_levels.")
    }
    colnames(pi.hat) <- treatment_levels
  }
  if (!is.null(mu.hat)) {
    mu.hat <- validate_outcome_matrix_input(mu.hat, nrow(data), "mu.hat", var_type_outcome$type)
    if (ncol(mu.hat) != length(treatment_levels)) {
      stop("mu.hat should be a matrix with ncol(mu.hat) = length(treatment_levels). The jth column should be estimated for the outcome regression of the jth treatment level in treatment_levels.")
    }
    colnames(mu.hat) <- treatment_levels
  }

  if (!is.null(m.hat)) {
    m.hat <- validate_outcome_vector_input(m.hat, nrow(data), "m.hat", var_type_outcome$type)
  }

  if (cross_fit_and_cv) {
    learner_pi <- make_cross_fitted(learner_pi, cross_validate = TRUE)
    learner_mu <- make_cross_fitted(learner_mu, cross_validate = TRUE)
    if (!is.null(learner_m)) {
      learner_m <- make_cross_fitted(learner_m, cross_validate = TRUE)
    }
  }

  factor_list <- list()
  if (!is.null(pi.hat)) {
    mean_fun <- function(learner_task) {
      pi.hat[match(learner_task$id, data[[id]]), , drop = FALSE]
    }
    factor_list$pi <- LF_known$new("pi", mean_fun = mean_fun, type = "mean")
  } else if (!is.null(learner_pi)) {
    if (!is.null(multinomial_learner)) {
      learner_pi <- multinomial_learner$new(learner_pi)
    }
    factor_list$pi <- LF_fit_hte3$new("pi", learner_pi, type = "mean")
  }

  if (!is.null(mu.hat)) {
    mean_fun <- function(learner_task) {
      mu.hat[match(learner_task$id, data[[id]]), , drop = FALSE]
    }
    factor_list$mu <- LF_known$new("mu", mean_fun = mean_fun, type = "mean")
  } else if (!is.null(learner_mu)) {
    learner_mu <- Lrnr_stratified_multivariate$new(learner = learner_mu, variable_stratify = treatment)
    factor_list$mu <- LF_fit_hte3$new("mu", learner_mu, type = "mean")
  }

  if (!is.null(learner_m) && is.null(m.hat)) {
    factor_list$m <- LF_fit_hte3$new("m", learner_m, type = "mean")
  }

  likelihood <- Likelihood$new(factor_list = factor_list)
  task <- hte3_Task$new(data, npsem, likelihood, id = id, weights = weights, folds = folds, ...)

  if (is.null(learner_m) || !is.null(m.hat)) {
    pi.hat <- as.matrix(task$get_nuisance_estimates("pi"))
    mu.hat <- as.matrix(task$get_nuisance_estimates("mu"))
    if (is.null(m.hat)) {
      m.hat <- as.vector(rowSums(pi.hat * mu.hat))
    }

    mean_fun <- function(learner_task) {
      m.hat[match(learner_task$id, data[[id]])]
    }
    new_factor_list <- list(m = LF_known$new("m", mean_fun = mean_fun, type = "mean"))
    task$likelihood$add_factors(new_factor_list)
  }

  task
}
