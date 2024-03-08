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
      params <- sl3:::args_to_list()
      base_learner <- Lrnr_cv_selector$new(loss_squared_error)
      super$initialize(params = params, base_learner = base_learner,
                       transform_function = NULL,
                       pseudo_outcome_type = "continuous",
                       pseudo_family = gaussian(),
                       ...)
    },
    get_pseudo_data = function(hte3_task, treatment_level = NULL, control_level = NULL, ...) {
      data <- hte3_task$data
      A <- hte3_task$get_tmle_node("treatment")
      Y <- hte3_task$get_tmle_node("outcome")
      print(treatment_level)
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
      pi.hat.1 <- truncate_propensity(pi.hat.1, A, treatment_level = treatment_level, truncation_method = "adaptive")
      pi.hat.0 <- truncate_propensity(pi.hat.0, A, treatment_level = control_level, truncation_method = "adaptive")

      index.mu.1 <- match(as.character(treatment_level), as.character(colnames(mu.hat)))
      index.mu.0 <- match(as.character(control_level), as.character(colnames(mu.hat)))
      mu.hat.1 <- mu.hat[, index.mu.1]
      mu.hat.0 <- mu.hat[, index.mu.0]
      # compute DR-learner pseudo-outcome based on EIF
      pseudo_outcome <- mu.hat.1 - mu.hat.0 + (A == treatment_level)/pi.hat.1 * (Y - mu.hat.1)  - (A == control_level)/pi.hat.0 * (Y - mu.hat.0)

      return(list(pseudo_outcome = pseudo_outcome))
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
      "CATE", "selector"
    )
  )
)
