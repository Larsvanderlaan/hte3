
#' @export
loss_crr_quasibinomial <- function (pred, observed)
{
  out <- log(1 + exp(pred)) - observed * pred
  attributes(out)$name <- "crr"
  return(out)
}


#' Lrnr_crr_DR_nonconvex Class
#'
#' This class defines an DR meta-learner for the conditional relative average treatment effect (crr).
#'
#' @format An R6 class with public methods to initialize the learner, create a regression task, and access the base learner.
#' @inheritParams Lrnr_cate_EP
#' @import Sieve
#' @export
Lrnr_crr_DR_selector <- R6Class(
  classname = "Lrnr_crr_DR_selector",
  inherit = Lrnr_hte,
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(
    treatment_level = NULL,
    control_level = NULL,
    ...
    ) {
      params <- sl3:::args_to_list()
      base_learner <- Lrnr_cv_selector$new(loss_crr_quasibinomial)
      super$initialize(params = params, base_learner = base_learner,
                       transform_function = NULL,
                       pseudo_outcome_type = "continuous",
                       pseudo_family = quasibinomial(),
                       ...)
    },
    get_pseudo_data = function(hte3_task, treatment_level = NULL, control_level = NULL, train = TRUE, ...) {
      data <- hte3_task$data
      A <- hte3_task$get_tmle_node("treatment")
      Y <- hte3_task$get_tmle_node("outcome")
      if(is.null(treatment_level)) treatment_level <- levels(factor(A))[2]
      if(is.null(control_level)) control_level <- levels(factor(A))[1]
      # should output a matrix where each column corresponds to a treatment level
      pi.hat <- as.matrix(hte3_task$get_nuisance_estimates("pi"))
      mu.hat <- as.matrix(hte3_task$get_nuisance_estimates("mu"))

      # get estimates for relevant treatment levels
      # assumes column names are treatment levels
      index.pi.1 <- match(as.character(treatment_level), as.character(colnames(pi.hat)))
      index.pi.0 <- match(as.character(control_level), as.character(colnames(pi.hat)))
      pi.hat.1 <- pi.hat[, index.pi.1]
      pi.hat.0 <- pi.hat[, index.pi.0]
      pi.hat.1 <- causalutils::truncate_propensity(pi.hat.1, A, treatment_level = treatment_level, truncation_method = "adaptive")
      pi.hat.0 <- causalutils::truncate_propensity(pi.hat.0, A, treatment_level = control_level, truncation_method = "adaptive")

      index.mu.1 <- match(as.character(treatment_level), as.character(colnames(mu.hat)))
      index.mu.0 <- match(as.character(control_level), as.character(colnames(mu.hat)))
      mu.hat.1 <- mu.hat[, index.mu.1]
      mu.hat.0 <- mu.hat[, index.mu.0]




      # get pseudo outcomes and weights for logistic regression.
      mu.tilde.1 <- mu.hat.1 + 1*(A == treatment_level) / pi.hat.1 * (Y - mu.hat.1)
      mu.tilde.0 <- mu.hat.0 + 1*(A == control_level) / pi.hat.0 * (Y - mu.hat.0)
      pseudo_outcome <- mu.tilde.1 /  (mu.tilde.0 + mu.tilde.1)
      pseudo_weights <- (mu.tilde.0 + mu.tilde.1) # note, observation weights are automatically added to the task and should not be added to pseudo_weights
      return(list(pseudo_outcome = pseudo_outcome, pseudo_weights = pseudo_weights))
    }
  ),
  active = list(
    coefficients = function() {
      base_learner <- self$fit_object$base_learner_trained
      return(base_learner$coefficients)
    },
    cv_risk = function() {
      base_learner <- self$fit_object$base_learner_trained
      return(base_learner$fit_object$cv_risk)
    }

  ),
  private = list(
    .treatment_type = c("binary_treatment", "categorical_treatment"),
    .properties = c(
      "crr", "EP"
    )
  )
)
