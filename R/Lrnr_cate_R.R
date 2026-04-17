#' Lrnr_cate_R Class
#'
#' This class defines the R-learner of Xie and Wager (2021) for estimation of the conditional average treatment effect.
#'
#' For continuous treatment, this implementation follows the partially linear
#' effect-model view in which the outcome regression is decomposed using an
#' `A * tau(X)` term rather than a fully general treatment-response surface.
#'
#' When the chosen modifier set `V` is a strict subset of the adjustment set
#' `W`, the natural target is `E[Y(1) - Y(0) | V] = E[tau(W) | V]`. The current
#' implementation does not generally target that object in the reduced-modifier
#' setting. Instead, it learns the overlap-weighted projection
#' `f_R(V) = E[Var(A|W) tau(W) | V] / E[Var(A|W) | V]`, which simplifies to
#' `E[e(W)(1-e(W)) tau(W) | V] / E[e(W)(1-e(W)) | V]` for binary treatment.
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
      # R-learner residualization is fit on the modifier task, so reduced
      # modifier sets do not generally recover E[Y(1)-Y(0) | V].
      warn_rlearner_reduced_modifier_target(hte3_task)

      A <- hte3_task$get_tmle_node("treatment")
      Y <- hte3_task$get_tmle_node("outcome")
      A_numeric <- coerce_numeric_treatment_values(A, "`Lrnr_cate_R`")
      A_levels <- sort(unique(A_numeric))
      pi.hat <- get_nuisance_matrix(hte3_task, "pi", "pi")
      pi.hat <- do.call(cbind, lapply(seq_along(A_levels), function(index) {
        truncate_propensity(pi.hat[, index], A_numeric, treatment_level = A_levels[index], truncation_method = "adaptive")
      }))
      e.hat <- as.vector(pi.hat %*% A_levels)
      m.hat <- get_nuisance_vector(hte3_task, "m", "m")
      residual_treatment <- A_numeric - e.hat
      safe_residual <- stabilize_treatment_residual(residual_treatment)

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
