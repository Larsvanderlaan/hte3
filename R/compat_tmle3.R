# Internal compatibility with tmle3/sl3 training and prediction methods.
LF_fit_hte3 <- R6Class(
  classname = "LF_fit_hte3",
  inherit = LF_fit,
  portable = TRUE,
  class = TRUE,
  public = list(
    delayed_train = function(tmle_task) {
      if (self$learner$is_trained) {
        return(self$learner)
      }

      outcome_node <- self$name
      learner_task <- tmle_task$get_regression_task(
        outcome_node,
        scale = TRUE,
        drop_censored = TRUE,
        is_time_variant = self$is_time_variant
      )

      sl3::delayed_learner_train(self$learner, learner_task)
    },
    get_mean = function(tmle_task, fold_number) {
      # The inherited sl3 method ignores the fold and calls predict(). Request
      # that full-fit prediction directly; preserve overridden fold methods.
      if (identical(
        body(self$learner$predict_fold),
        body(sl3::Lrnr_base$public_methods$predict_fold)
      )) {
        fold_number <- "full"
      }
      super$get_mean(tmle_task, fold_number)
    }
  )
)
