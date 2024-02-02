
#' Lrnr_crr_T Class
#'
#' This class constructs a T-learner for estimation of the conditional relative average treatment effect (crr)
#'
#' @format An R6 class with public methods to initialize the learner, create a regression task, and access the base learner.
#' @inheritParams Lrnr_cate_T
#' @export
Lrnr_crr_T <- R6Class(
  classname = "Lrnr_crr_T", inherit = Lrnr_hte,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(base_learner, treatment_level = 1, control_level = 0, stratify_by_treatment = TRUE, ...) {
      params <- sl3:::args_to_list()
      super$initialize(params = params, base_learner = base_learner,
                       ...)
    },
    get_pseudo_data = function(hte3_task, ...) {
      stop("Not used by T-learner.")
    },
    make_metalearner_task = function(hte3_task, train = TRUE) {
      params <- self$params
      training_task <- hte3_task$next_in_chain(covariates = c(hte3_task$npsem$modifiers$variables, hte3_task$npsem$treatment$variables),
                                               outcome = hte3_task$npsem$outcome$variables)

      # remake task so outcome type is continuous
      training_task <- sl3_Task$new(training_task$internal_data,
                                    column_names =  training_task$column_names,
                                    row_index = training_task$row_index,
                                    nodes = training_task$nodes,
                                    folds = training_task$folds)
      return(training_task)
    }
  ),
  active = list(
    base_learner = function() {
      return(private$.base_learner)
    }
  ),
  private = list(
    .properties = c(
      "continuous", "binomial", "categorical", "importance",
      "weights"),
    .train = function(hte3_task) {
      args <- self$params
      learner_task <- self$make_metalearner_task(hte3_task)
      if(self$params$stratify_by_treatment) {
        trt <- hte3_task$get_tmle_node("treatment")
        learner_task <- learner_task$next_in_chain(covariates = hte3_task$npsem$modifiers$variables)
        learner_trained_trt <- self$base_learner$train(learner_task[trt == self$params$treatment_level])
        learner_trained_control <- self$base_learner$train(learner_task[trt == self$params$control_level])
        fit_object = list(learner_trained_trt = learner_trained_trt, learner_trained_control = learner_trained_control)
      } else{
        learner_trained_pooled <- self$base_learner$train(learner_task)
        fit_object = list(learner_trained_pooled = learner_trained_pooled)
      }

      return(fit_object)
    },
    .predict = function(hte3_task) {
      fit_object <- self$fit_object
      learner_task <- self$make_metalearner_task(hte3_task)
      if(self$params$stratify_by_treatment) {
        learner_trained_trt <- fit_object$learner_trained_trt
        learner_trained_control <- fit_object$learner_trained_control
        learner_task <- learner_task$next_in_chain(covariates = hte3_task$npsem$modifiers$variables)
        mu0.hat <- learner_trained_control$predict(learner_task)
        mu1.hat <- learner_trained_trt$predict(learner_task)
      } else{
        learner_trained_pooled <- fit_object$learner_trained_pooled
        # generate counterfactual tasks
        dat0 <- data.table(rep(self$params$control_level, hte3_task$nrow))
        dat1 <- data.table(rep(self$params$treatment_level, hte3_task$nrow))
        names(dat0) <- names(dat1) <-  hte3_task$npsem$treatment$variables
        cf0_hte3_task <- hte3_task$generate_counterfactual_task(uuid::UUIDgenerate(), dat0)
        cf1_hte3_task <- hte3_task$generate_counterfactual_task(uuid::UUIDgenerate(), dat1)
        mu0.hat <- learner_trained_pooled$predict(self$make_metalearner_task(cf0_hte3_task))
        mu1.hat <- learner_trained_pooled$predict(self$make_metalearner_task(cf1_hte3_task))
      }
      mu0.hat <- pmax(mu0.hat, 1e-10)
      mu1.hat <- pmax(mu1.hat, 1e-10)
      predictions <- log(mu1.hat) - log(mu0.hat)
      return(predictions)
    },
    .required_packages = c("sl3"),
    .base_learner = NULL
  )
)


#' Lrnr_cate_S Class
#'
#' This class constructs an S-learner for estimation of the conditional relative average treatment effect (crr).
#' This object is a wrapper for \code{Lrnr_crr_T} with \code{stratify_by_treatment} = \code{FALSE}.
#'
#' @format An R6 class with public methods to initialize the learner, create a regression task, and access the base learner.
#' @inheritParams Lrnr_crr_T
#' @export
Lrnr_crr_S <- R6Class(
  classname = "Lrnr_crr_S", inherit = Lrnr_crr_T,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(base_learner, treatment_level = 1, control_level = 0, ...) {
      super$initialize(base_learner = base_learner,
                       treatment_level = treatment_level,
                       control_level = control_level,
                       stratify_by_treatment = FALSE,
                       ...)
    }
  ))
