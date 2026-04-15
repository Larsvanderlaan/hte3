
#' Lrnr_crr_T Class
#'
#' This class constructs a T-learner for estimation of the conditional relative average treatment effect (crr)
#'
#' @format An R6 class with public methods to initialize the learner, create a regression task, and access the base learner.
#' @param base_learner A \code{\link{sl3}} learner object inheriting from \code{\link[sl3]{Lrnr_base}} that specifies the base supervised learning algorithm used by the meta-learner.
#' @export
Lrnr_crr_T <- R6Class(
  classname = "Lrnr_crr_T", inherit = Lrnr_hte,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(base_learner, treatment_level = 1, control_level = 0, stratify_by_treatment = TRUE, ...) {
      params <- list(
        base_learner = base_learner,
        treatment_level = treatment_level,
        control_level = control_level,
        stratify_by_treatment = stratify_by_treatment,
        ...
      )
      super$initialize(params = params, base_learner = base_learner,
                       ...)
    },
    get_pseudo_data = function(hte3_task, ...) {
      stop("Not used by T-learner.")
    },
    make_metalearner_task = function(hte3_task, train = TRUE) {
      make_tlearner_task(hte3_task)
    }
  ),
  active = list(
    base_learner = function() {
      return(private$.base_learner)
    }
  ),
  private = list(
    .treatment_type = c("binomial", "categorical"),
    .properties = c(
      "continuous", "binomial", "categorical", "importance",
      "weights"),
    .train = function(hte3_task) {
      self$check_treatment_type(hte3_task)
      validate_nonnegative_outcome(hte3_task$get_tmle_node("outcome"), label = "CRR outcomes")
      learner_task <- self$make_metalearner_task(hte3_task)
      contrast <- resolve_treatment_levels(hte3_task, self$params$treatment_level, self$params$control_level)
      train_tlearner_models(
        self$base_learner,
        learner_task,
        hte3_task,
        treatment_level = contrast$treatment_level,
        control_level = contrast$control_level,
        stratify_by_treatment = self$params$stratify_by_treatment
      )
    },
    .predict = function(hte3_task) {
      fit_object <- self$fit_object
      learner_task <- self$make_metalearner_task(hte3_task)
      contrast <- resolve_treatment_levels(hte3_task, self$params$treatment_level, self$params$control_level)
      means <- predict_tlearner_means(
        fit_object,
        hte3_task,
        learner_task,
        treatment_level = contrast$treatment_level,
        control_level = contrast$control_level,
        stratify_by_treatment = self$params$stratify_by_treatment
      )
      log(bound_positive(means$mu1.hat)) - log(bound_positive(means$mu0.hat))
    },
    .required_packages = c("sl3"),
    .base_learner = NULL
  )
)
