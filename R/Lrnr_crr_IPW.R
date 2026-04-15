#' Lrnr_crr_IPW Class
#'
#' This class defines a inverse probability weighted (IPW) meta-learner of the conditional relative average treatment effect (crr).
#'
#' @format An R6 class with public methods to initialize the learner, create a regression task, and access the base learner.
#' @param base_learner A \code{\link{sl3}} learner object inheriting from \code{\link[sl3]{Lrnr_base}} that specifies the base supervised learning algorithm used by the meta-learner.
#' @import Sieve
#' @export
Lrnr_crr_IPW <- R6Class(
  classname = "Lrnr_crr_IPW", inherit = Lrnr_hte,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(base_learner, treatment_level = NULL, control_level = NULL, ...) {
      params <- list(base_learner = base_learner, treatment_level = treatment_level, control_level = control_level, ...)
      super$initialize(params = params, base_learner = base_learner,
                       transform_function = bounded_qlogis,
                       pseudo_outcome_type = "binomial",
                       pseudo_family = binomial(),
                       ...)
    },
    get_pseudo_data = function(hte3_task, treatment_level = NULL, control_level = NULL, train = TRUE, ...) {
      A <- hte3_task$get_tmle_node("treatment")
      Y <- hte3_task$get_tmle_node("outcome")
      validate_nonnegative_outcome(Y, label = "CRR outcomes")

      contrast <- resolve_treatment_levels(hte3_task, treatment_level, control_level)
      pi.hat <- get_nuisance_matrix(hte3_task, "pi", "pi")
      pi.hat.1 <- extract_treatment_column(pi.hat, contrast$treatment_level, "pi")
      pi.hat.0 <- extract_treatment_column(pi.hat, contrast$control_level, "pi")
      n <- length(Y)
      truncation_method <- ifelse(n >= 500, "isotonic", "adaptive")
      pi.hat.1 <- truncate_propensity(pi.hat.1, A, treatment_level = contrast$treatment_level, truncation_method = truncation_method)
      pi.hat.0 <- truncate_propensity(pi.hat.0, A, treatment_level = contrast$control_level, truncation_method = truncation_method)

      pseudo_outcome <- as.numeric(A == contrast$treatment_level)
      invpi <- (A == contrast$treatment_level) / pi.hat.1 + (A == contrast$control_level) / pi.hat.0
      pseudo_weights <- Y * invpi
      if (!any(pseudo_weights > 0, na.rm = TRUE)) {
        stop("`Lrnr_crr_IPW` requires at least one observation with positive CRR weight.", call. = FALSE)
      }
      return(list(pseudo_outcome = pseudo_outcome, pseudo_weights = pseudo_weights))
    }
  ),
  active = list(

  ),
  private = list(
    .treatment_type = c("binomial", "categorical"),
    .properties = c(
      "crr", "IPW"
    )
  )
)
