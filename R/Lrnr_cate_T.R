
#' Lrnr_cate_T Class
#'
#' This class constructs a T-learner (or S-learner) for estimation of the conditional average treatment effect (CATE)
#'
#' @format An R6 class with public methods to initialize the learner, create a regression task, and access the base learner.
#' @param stratify_by_treatment Logical indicating whether to estimate outcome regression separately in each treatment arm (i.e., T-learner) or pooled across treatment arms (i.e., S-learner).
#' @param second_stage_regression Logical indicating whether to regress the
#' first-stage contrast estimates onto the modifier set. This defaults to
#' \code{TRUE}. Setting it to \code{FALSE} is only supported when the modifier
#' set and confounder set are the same.
#' @param base_learner A \code{\link{sl3}} learner object inheriting from \code{\link[sl3]{Lrnr_base}} that specifies the base supervised learning algorithm used by the meta-learner.
#' @details The learner first estimates outcome regressions as functions of the
#' adjustment set \code{W}. By default, it then regresses the resulting
#' first-stage contrast \code{tau(W)} onto the modifier set \code{V}. When
#' \code{stratify_by_treatment = FALSE}, the first-stage outcome model is fit as
#' a pooled S-learner-style regression and evaluated counterfactually under each
#' treatment level.
#' @export
Lrnr_cate_T <- R6Class(
  classname = "Lrnr_cate_T", inherit = Lrnr_hte,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(base_learner, treatment_level = 1, control_level = 0,
                          stratify_by_treatment = TRUE,
                          second_stage_regression = TRUE,
                          ...) {
      params <- list(
        base_learner = base_learner,
        treatment_level = treatment_level,
        control_level = control_level,
        stratify_by_treatment = stratify_by_treatment,
        second_stage_regression = second_stage_regression,
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
      learner_task <- self$make_metalearner_task(hte3_task)
      contrast <- resolve_treatment_levels(hte3_task, self$params$treatment_level, self$params$control_level)
      fit_object <- train_tlearner_models(
        self$base_learner,
        learner_task,
        hte3_task,
        treatment_level = contrast$treatment_level,
        control_level = contrast$control_level,
        stratify_by_treatment = self$params$stratify_by_treatment
      )
      first_stage_means <- predict_tlearner_means(
        fit_object,
        hte3_task,
        learner_task,
        treatment_level = contrast$treatment_level,
        control_level = contrast$control_level,
        stratify_by_treatment = self$params$stratify_by_treatment
      )
      fit_object$use_second_stage <- isTRUE(self$params$second_stage_regression)
      if (fit_object$use_second_stage) {
        fit_object$effect_learner_trained <- fit_tlearner_projection(
          self$base_learner,
          hte3_task,
          outcome_values = first_stage_means$mu1.hat - first_stage_means$mu0.hat,
          family = stats::gaussian(),
          outcome_type = "continuous"
        )
      } else if (!modifiers_equal_confounders(hte3_task)) {
        stop(
          "T-learner without second-stage regression is only supported when `modifiers` and `confounders` are the same.",
          call. = FALSE
        )
      }
      fit_object$prediction_template <- make_prediction_template(hte3_task)
      fit_object$contrast_levels <- contrast
      fit_object
    },
    .predict = function(hte3_task) {
      fit_object <- self$fit_object
      if (isTRUE(fit_object$use_second_stage)) {
        return(predict_tlearner_projection(fit_object$effect_learner_trained, hte3_task))
      }
      contrast <- fit_object$contrast_levels
      learner_task <- self$make_metalearner_task(hte3_task)
      means <- predict_tlearner_means(
        fit_object,
        hte3_task,
        learner_task,
        treatment_level = contrast$treatment_level,
        control_level = contrast$control_level,
        stratify_by_treatment = self$params$stratify_by_treatment
      )
      means$mu1.hat - means$mu0.hat
    },
    .required_packages = c("sl3"),
    .base_learner = NULL
  )
)
