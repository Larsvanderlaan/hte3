
loss_crr_quasibinomial <- function (pred, observed)
{
  pred <- coerce_scalar_numeric(pred, "selector predictions")
  observed <- coerce_scalar_numeric(observed, "selector targets")

  if (length(pred) != length(observed)) {
    stop("CRR selector loss requires `pred` and `observed` to have matching lengths.", call. = FALSE)
  }

  out <- log1p(exp(pred)) - observed * pred
  attributes(out)$name <- "crr"
  return(out)
}


#' Lrnr_crr_DR_nonconvex Class
#'
#' This class defines an DR meta-learner for the conditional relative average treatment effect (crr).
#'
#' @format An R6 class with public methods to initialize the learner, create a regression task, and access the base learner.
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
      params <- list(treatment_level = treatment_level, control_level = control_level, ...)
      base_learner <- Lrnr_cv_selector$new(loss_crr_quasibinomial)
      super$initialize(params = params, base_learner = base_learner,
                       transform_function = NULL,
                       pseudo_outcome_type = "continuous",
                       pseudo_family = quasibinomial(),
                       ...)
    },
    get_pseudo_data = function(hte3_task, treatment_level = NULL, control_level = NULL, train = TRUE, ...) {
      compute_crr_ratio_pseudo_data(
        hte3_task,
        treatment_level = treatment_level,
        control_level = control_level
      )
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
    .treatment_type = c("binomial", "categorical"),
    .properties = c(
      "crr", "EP"
    )
  )
)
