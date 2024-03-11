
call_with_args <- sl3:::call_with_args


#' Get Automated Machine Learning (AutoML) Learner
#'
#' This function returns a predefined stack of off-the-shelf learning algorithms for Automated Machine Learning (AutoML).
# The stack includes a variety of learners with tuned parameters for quick experimentation and model selection.
# The emphasis is on speed and simplicity.
#
#' The following base learners are included in the stack:
#' \itemize{
#'   \item \code{Lrnr_glmnet}: Elastic Net Regularized Generalized Linear Model.
#'   \item \code{Lrnr_gam}: Generalized Additive Model.
#'   \item \code{Lrnr_earth}: Multivariate Adaptive Regression Spline Model.
#'   \item \code{Lrnr_ranger}: Random Forest Model.
#'   \item \code{Lrnr_xgboost}: XGBoost Model with different depth configurations.
#' }
#'
#' @return A `sl3::Stack` object containing a predefined set of off-the-shelf learners for AutoML.
#'
#' @import sl3
#' @export
get_autoML <- function() {
  # Simple off-the-shelf learning algorithm
  # Emphasis on speed.
  # xgboost tuning parameters are somewhat arbitrary but decent for off-the-shelf.
  learner <- Stack$new(Lrnr_glmnet$new(),
                       Lrnr_gam$new(),
                       Lrnr_earth$new(degree =2),
                       Lrnr_ranger$new(max.depth = 10),
                       Lrnr_xgboost$new(min_child_weight = 15, max_depth = 3, nrounds = 40, eta = 0.15, subsample = 0.9),
                       Lrnr_xgboost$new(min_child_weight = 15, max_depth = 4, nrounds = 40, eta = 0.15, subsample = 0.9) ,
                       Lrnr_xgboost$new(min_child_weight = 15, max_depth = 5, nrounds = 40, eta = 0.15, subsample = 0.9)
  )
  return(learner)


}




#' Create Cross-Fitted Learner
#'
#' This function takes a learner and returns a cross-fitted version of it. Cross-fitting involves
#' fitting the learner to different subsets of the data while using the complementary subsets for validation,
#' to provide a more robust estimate of model performance.
#'
#' @param learner The learner to be cross-fitted. This can be a single learner or a stacked learner.
#' @param calibrate Currently not used. Logical indicating whether to calibrate the learner (default is \code{FALSE}).
#' @param cross_validate Logical indicating whether to perform cross-validation (default is \code{TRUE} if \code{learner} is a stacked learner).
#'
#' @return A cross-fitted version of the input learner.
#'
#' @export
make_cross_fitted <- function(learner, calibrate = FALSE, cross_validate = inherits(learner, "Stack")) {
  if(is.null(learner)) {
    return(null)
  }
  learner <- Lrnr_cv$new(learner)
  if(cross_validate) {
    learner <- make_learner(Pipeline, learner, Lrnr_cv_selector$new(loss_squared_error))
  }
  return(learner)

}




#' Truncate and Calibrate Propensity Scores
#'
#' This function truncates and calibrates propensity scores to ensure bounded values and improve their reliability.
#'
#' @param pi.hat A numeric vector containing estimates of the propensity score \(\hat{\pi}(W)\) corresponding to \code{treatment_level}.
#' @param A A numeric vector of treatment values.
#' @param treatment_level A numeric value indicating the treatment level to calibrate (default is the maximum value of \(\mathbf{A}\)).
#' @param truncation_method A character string indicating the truncation method: "isotonic", "adaptive", "deterministic", or "none".
#'
#' @return A numeric vector of truncated and calibrated propensity scores.
#'
#'
#' @description
#' The truncation_method parameter controls how the propensity scores are calibrated and truncated.
#' - "isotonic": Performs isotonic calibration, providing both data-adaptive truncation and calibration for the propensity scores.
#' - "adaptive": Adapts the truncation level using a loss function for the inverse propensity score.
#' - "deterministic": Bounds the estimates away from 0 and 1 using the threshold 25/sqrt(n)/log(n).
#' - "none": Bounds the estimates away from 0 and 1 using the threshold 1/n.
#' @export
truncate_propensity <- function(pi.hat, A, treatment_level = max(A), truncation_method = c("isotonic", "adaptive" , "deterministic" , "none")) {
  n <- length(A)
  treatment_indicator <- as.numeric(A==treatment_level)
  pi.hat <- pmax(pi.hat, 1e-8)
  pi.hat <- pmin(pi.hat, 1 - 1e-8)

  if(truncation_method == "isotonic") {
    # Calibrates and automatically truncates using isotonic regression/calibration
    calibrator_pi <- as.stepfun(isoreg(pi.hat, treatment_indicator))
    pi.hat.star <- calibrator_pi(pi.hat)

  } else if(truncation_method == "adaptive") {
    # Choose lower truncation level by minimizing loss for inverse propensity score
    risk_function <- function(cutoff) {
      pi.hat <- pmax(pi.hat, cutoff)
      alpha.hat <- ifelse(treatment_indicator==1, 1/pi.hat, 0)
      alpha1.hat <- 1/pi.hat
      mean( alpha.hat^2 - 2*alpha1.hat)
    }
    cutoff <- optim(1/sqrt(n), fn = risk_function, method = "Brent", lower = 1/n, upper = 0.5)$par
    pi.hat.star <- pmax(pi.hat, cutoff)

    # Choose upper truncation level
    risk_function <- function(cutoff) {
      pi.hat <- pmin(pi.hat, 1-cutoff)
      alpha.hat <- ifelse(treatment_indicator == 0, 1/(1-pi.hat), 0)
      alpha1.hat <- 1/(1-pi.hat)
      mean( alpha.hat^2 - 2*alpha1.hat)
    }
    cutoff <- optim(1/sqrt(n), fn = risk_function, method = "Brent", lower = 1/n, upper = 0.5)$par
    pi.hat.star <- pmin(pi.hat.star, 1-cutoff)

  } else if(truncation_method == "deterministic") {
    # Rule of thumb truncation level based on sample size
    pi.hat.star <- pmax(pi.hat, 25/sqrt(n)/log(n))
    pi.hat.star <- pmin(pi.hat, 1 - 25/sqrt(n)/log(n))
  }

  pi.hat.star <- pmax(pi.hat, 1/n)
  pi.hat.star <- pmin(pi.hat, 1 - 1/n)
  return(pi.hat.star)

}


#' Internal use.
#' Converts a single outcome learner into a multivariate outcome learner
#' that predicts a matrix of predictions obtained by evaluating the
#' single outcome learner at each possible value of variable_stratify.
#
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
      if(!is.factor(strata_ids)) strata_ids <- factor(strata_ids)
      return(levels(strata_ids))
    }
  ),
  private = list(
    .train = function(task) {
      learner <- self$params$learner
      learner_trained <- learner$train(task)
      return(list(learner_trained = learner_trained))
    },
    .predict = function(task = NULL) {

      learner_trained <- self$fit_object$learner_trained
      variable_stratify <- self$params$variable_stratify

      strata_levels <- self$strata_levels

      prediction_mat <- as.matrix(do.call(cbind,lapply(strata_levels, function(strata) {
        if(all(is.numeric(task$data[[variable_stratify]]))) {
          strata <- as.numeric(strata)
        } else {
          strata <- as.character(strata)
        }
        # construct subtask
        new_data <- task$data
        set(new_data, , variable_stratify, rep(strata, nrow(new_data)))
        new_task <- sl3_Task$new(new_data, nodes = task$nodes)
        return(learner_trained$predict(new_task))
      })))
      colnames(prediction_mat) <- as.character(strata_levels)
      predictions <-  sl3:::pack_predictions(prediction_mat)
      return(predictions)
    },
    .required_packages = NULL
  )
)



#' Calibrate Predictor-Outcome Relationship
#'
#' This function calibrates the predictor-outcome relationship using isotonic regression.
#'
#' @param predictor A numeric vector representing the predictor variable.
#' @param outcome A numeric vector representing the outcome variable.
#'
#' @return A step function that represents the calibrated relationship between the predictor and the outcome.
#'
#' @export
calibrate <- function(predictor, outcome) {
  return(as.stepfun(isoreg(predictor, outcome))(predictor))
}





#' Helper function for estimation of outcome regression using \code{sl3}.
#' @import sl3 data.table
#' @export
estimate_mu <- function(W, A, Y, treatment_level, learner = Lrnr_gam$new(family = "gaussian"), weights = NULL, cross_fit_and_cv = TRUE, stratified_by_trt = TRUE, return_learner = FALSE, folds = 10,...) {

  W <- as.data.table(W)
  if(length(names(W))==0) {
    names(W) <- paste0("W", 1:ncol(W))
  }

  data <- cbind(W, data.table(A, Y, weights))

  if(stratified_by_trt) {
    learner <- Lrnr_stratified$new(learner = learner, variable_stratify = "A")
  }
  if(cross_fit_and_cv) {
    learner <- make_learner(Pipeline, Lrnr_cv$new(learner), Lrnr_cv_selector$new(loss_squared_error()))
  }

  task <- sl3_Task$new(data, covariates = c(names(W), "A"), outcome = "Y", weights = weights, folds = folds, ...)
  learner_trained <- learner$train(task[A == treatment_level])
  mu.hat <- learner_trained$predict(task)

  output <- list(mu.hat = mu,hat, folds = folds)
  if(return_learner) output$learner <- learner_trained
  return(output)
}

#' Helper function for estimation of treatment-averaged outcome regression using \code{sl3}.
#' @import sl3 data.table
#' @export
estimate_m <- function(W, Y, learner = Lrnr_gam$new(family = "gaussian"), weights = NULL, cross_fit_and_cv = TRUE,  return_learner = FALSE, folds = 10,...) {

  W <- as.data.table(W)
  if(length(names(W))==0) {
    names(W) <- paste0("W", 1:ncol(W))
  }

  data <- cbind(W, data.table( Y, weights))

  if(cross_fit_and_cv) {
    learner <- make_learner(Pipeline, Lrnr_cv$new(learner), Lrnr_cv_selector$new(loss_squared_error()))
  }

  task <- sl3_Task$new(data, covariates = names(W), outcome = "Y", weights = weights, folds = folds, ...)
  learner_trained <- learner$train(task)
  m.hat <- learner_trained$predict(task)

  output <- list(m.hat = m,hat, folds = folds)
  if(return_learner) output$learner <- learner_trained
  return(output)
}

#' Helper function for estimation of propensity score using \code{sl3}.
#' @import sl3 data.table
#' @export
estimate_pi <- function(W, A, binomial_learner = Lrnr_gam$new(family = "gaussian"), weights = NULL, treatment_level = max(A), cross_fit_and_cv = TRUE,  return_learner = FALSE, folds = 10,...) {

  A <- factor(A, levels = sort(unique(A)))
  learner <- Lrnr_independent_binomial$new(binomial_learner)


  W <- as.data.table(W)
  if(length(names(W))==0) {
    names(W) <- paste0("W", 1:ncol(W))
  }

  data <- cbind(W, data.table(A, weights))

  if(cross_fit_and_cv) {
    learner <- make_learner(Pipeline, Lrnr_cv$new(learner), Lrnr_cv_selector$new(loss_squared_error()))
  }

  task <- sl3_Task$new(data, covariates = names(W), outcome = "A", weights = weights, folds = folds, ...)
  learner_trained <- learner$train(task)
  pi.hat <- sl3:::unpack_predictions(learner_trained$predict(task))

  output <- list(pi.hat = pi.hat, levels = levels(A), folds = folds)
  if(return_learner) output$learner <- learner_trained
  return(output)
}


