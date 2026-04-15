#' Lrnr_cate_R Class
#'
#' This class defines the R-learner of Xie and Wager (2021) for estimation of the conditional average treatment effect.
#'
#' @format An R6 class with public methods to initialize the learner, create a regression task, and access the base learner.
#' @param base_learner A \code{\link{sl3}} learner object inheriting from \code{\link[sl3]{Lrnr_base}} that specifies the base supervised learning algorithm used by the meta-learner.
#' @export
Lrnr_cate_R <- R6Class(
  classname = "Lrnr_cate_R", inherit = Lrnr_hte,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(base_learner, ...) {
      params <- list(base_learner = base_learner, ...)
      super$initialize(params = params, base_learner = base_learner,
                       ...)
    },
    get_pseudo_data = function(hte3_task, ...) {
      A <- hte3_task$get_tmle_node("treatment")
      Y <- hte3_task$get_tmle_node("outcome")
      A_numeric <- suppressWarnings(as.numeric(as.character(A)))
      if (anyNA(A_numeric)) {
        stop("`Lrnr_cate_R` requires numeric treatment coding.", call. = FALSE)
      }

      A_levels <- suppressWarnings(as.numeric(as.character(levels(factor(A)))))
      if (anyNA(A_levels)) {
        stop("`Lrnr_cate_R` requires numeric treatment coding.", call. = FALSE)
      }
      pi.hat <- get_nuisance_matrix(hte3_task, "pi", "pi")
      pi.hat <- do.call(cbind, lapply(seq_along(A_levels), function(index) {
        truncate_propensity(pi.hat[, index], A_numeric, treatment_level = A_levels[index], truncation_method = "adaptive")
      }))
      e.hat <- as.vector(pi.hat %*% A_levels)
      m.hat <- get_nuisance_vector(hte3_task, "m", "m")
      residual_treatment <- A_numeric - e.hat
      safe_residual <- ifelse(abs(residual_treatment) < 1e-8, sign(residual_treatment + (residual_treatment == 0)) * 1e-8, residual_treatment)

      # compute R-learner pseudo-outcome based on EIF
      pseudo_outcome <- (Y - m.hat) / safe_residual
      pseudo_weights <- bound_positive(residual_treatment^2)

      return(list(pseudo_outcome = pseudo_outcome, pseudo_weights = pseudo_weights))
    }
  ),
  active = list(
  ),
  private = list(
    .treatment_type = c("binomial", "continuous"),
    .properties = c(
      "cate", "pcate", "Rlearner"
    )
  )
)
