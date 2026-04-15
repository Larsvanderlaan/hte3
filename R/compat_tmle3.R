# Internal compatibility shim for tmle3/sl3 versions where LF_fit does not
# namespace-qualify delayed_learner_train().
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
    }
  )
)
