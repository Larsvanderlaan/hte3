#' Lrnr_cate_R Class
#'
#' This class defines the R-learner of Xie and Wager (2021) for estimation of the conditional average treatment effect.
#'
#' @format An R6 class with public methods to initialize the learner, create a regression task, and access the base learner.
#' @inheritParams Lrnr_hte
#' @export
Lrnr_cate_R <- R6Class(
  classname = "Lrnr_cate_R", inherit = Lrnr_hte,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(base_learner, ...) {
      params <- sl3:::args_to_list()
      super$initialize(params = params, base_learner = base_learner,
                       ...)
    },
    get_pseudo_data = function(hte3_task, ...) {
      A <- hte3_task$get_tmle_node("treatment")
      Y <- hte3_task$get_tmle_node("outcome")
      A_levels <- as.numeric(levels(factor(A)))
      pi.hat <- as.matrix(hte3_task$get_nuisance_estimates("pi"))
      pi.hat <- do.call(cbind, lapply(seq_along(A_levels), function(index) {
        pi.hat.i <- causalutils::truncate_propensity(pi.hat[, index], A, treatment_level = A_levels[index], truncation_method = "adaptive")
      }))
      e.hat <- as.vector(pi.hat %*% A_levels)
      m.hat <- as.vector(hte3_task$get_nuisance_estimates("m"))

      # compute R-learner pseudo-outcome based on EIF
      pseudo_outcome <- (Y - m.hat) / (A - e.hat)
      pseudo_weights <- (A - e.hat)^2

      return(list(pseudo_outcome = pseudo_outcome, pseudo_weights = pseudo_weights))
    }
  ),
  active = list(
  ),
  private = list(
    .treatment_type = c("binary_treatment", "continuous_treatment"),
    .properties = c(
      "cate", "pcate", "Rlearner"
    )
  )
)
