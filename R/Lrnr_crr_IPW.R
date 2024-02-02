#' Lrnr_crr_IPW Class
#'
#' This class defines a inverse probability weighted (IPW) meta-learner of the conditional relative average treatment effect (crr).
#'
#' @format An R6 class with public methods to initialize the learner, create a regression task, and access the base learner.
#' @inheritParams Lrnr_cate_EP
#' @import Sieve
#' @export
Lrnr_crr_IPW <- R6Class(
  classname = "Lrnr_crr_IPW", inherit = Lrnr_hte,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(base_learner, treatment_level = NULL, control_level = NULL, ...) {
      params <- sl3:::args_to_list()
      super$initialize(params = params, base_learner = base_learner,
                       transform_function = stats::qlogis,
                       pseudo_outcome_type = "binomial",
                       pseudo_family = binomial(),
                       ...)
    },
    get_pseudo_data = function(hte3_task, treatment_level = NULL, control_level = NULL, train = TRUE, ...) {
      data <- hte3_task$data
      n <- nrow(data)
      A <- hte3_task$get_tmle_node("treatment")
      Y <- hte3_task$get_tmle_node("outcome")
      if(is.null(treatment_level)) treatment_level <- levels(factor(A))[2]
      if(is.null(control_level)) control_level <- levels(factor(A))[1]
      # should output a matrix where each column corresponds to a treatment level
      pi.hat <- as.matrix(hte3_task$get_nuisance_estimates("pi"))

      # get estimates for relevant treatment levels
      # assumes column names are treatment levels
      index.pi.1 <- match(as.character(treatment_level), as.character(colnames(pi.hat)))
      index.pi.0 <- match(as.character(control_level), as.character(colnames(pi.hat)))
      pi.hat.1 <- pi.hat[, index.pi.1]
      pi.hat.0 <- pi.hat[, index.pi.0]
      truncation_method <- ifelse(n >= 500, "isotonic", "adaptive")
      truncation_method <- "adaptive"
      pi.hat.1 <- causalutils::truncate_propensity(pi.hat.1, A, treatment_level = treatment_level, truncation_method = truncation_method)
      pi.hat.0 <- causalutils::truncate_propensity(pi.hat.0, A, treatment_level = control_level, truncation_method = truncation_method)
      pi.hat.1 <- pmax(pi.hat.1, 1e-10)
      pi.hat.0 <- pmax(pi.hat.0, 1e-10)



      pseudo_outcome <- as.numeric(A == treatment_level)
      invpi <- (A==treatment_level)/pi.hat.1 + (A==control_level)/pi.hat.0
      pseudo_weights <- Y * invpi
      return(list(pseudo_outcome = pseudo_outcome, pseudo_weights = pseudo_weights))
    }
  ),
  active = list(

  ),
  private = list(
    .treatment_type = c("binary_treatment", "categorical_treatment"),
    .properties = c(
      "crr", "IPW"
    )
  )
)
