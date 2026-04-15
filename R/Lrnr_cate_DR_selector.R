#' @export
Lrnr_cate_DR_selector <- R6Class(
  classname = "Lrnr_cate_DR_selector",
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
      base_learner <- Lrnr_cv_selector$new(loss_squared_error)
      super$initialize(params = params, base_learner = base_learner,
                       transform_function = NULL,
                       pseudo_outcome_type = "continuous",
                       pseudo_family = gaussian(),
                       ...)
    },
    get_pseudo_data = function(hte3_task, treatment_level = NULL, control_level = NULL, ...) {
      list(
        pseudo_outcome = compute_cate_dr_pseudo_outcome(
          hte3_task,
          treatment_level = treatment_level,
          control_level = control_level
        )
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
      "CATE", "selector"
    )
  )
)
