#' Lrnr_cate_DR Class
#'
#' This class defines a doubly-robust (DR) meta-learner of the conditional average treatment effect.
#'
#' @format An R6 class with public methods to initialize the learner, create a regression task, and access the base learner.
#' @inheritParams Lrnr_hte
#' @param base_learner A \code{\link{sl3}} learner object inheriting from \code{\link[sl3]{Lrnr_base}} that specifies the base supervised learning algorithm used by the meta-learner.
#' @param treatment_level A treatment level encoding the treatment assignment of interest.
#' @param control_level A treatment level encoding the control (or reference) treatment assignment.
#' @export
Lrnr_cate_DR <- R6Class(
  classname = "Lrnr_cate_DR", inherit = Lrnr_hte,
  portable = TRUE, class = TRUE,
  public = list(
    #' @inheritParams Lrnr_hte
    initialize = function(base_learner, treatment_level = NULL, control_level = NULL, ...) {
      params <- sl3:::args_to_list()
      super$initialize(params = params, base_learner = base_learner,
                       ...)
    },
    #' @inheritParams Lrnr_hte
    get_pseudo_data = function(hte3_task, treatment_level = NULL, control_level = NULL, ...) {
      data <- hte3_task$data
      A <- hte3_task$get_tmle_node("treatment")
      Y <- hte3_task$get_tmle_node("outcome")
      if(is.null(treatment_level)) treatment_level <- levels(factor(A))[2]
      if(is.null(control_level)) control_level <- levels(factor(A))[1]



      # should output a matrix where each column corresponds to a treatment level
      pi.hat <- as.matrix(hte3_task$get_nuisance_estimates("pi"))
      mu.hat <- as.matrix(hte3_task$get_nuisance_estimates("mu"))

      # get estimates for relevant treatment levels
      # assumes column names are treatment levels
      index.pi.1 <- match(as.character(treatment_level), as.character(colnames(pi.hat)))
      index.pi.0 <- match(as.character(control_level), as.character(colnames(pi.hat)))
      pi.hat.1 <- pi.hat[, index.pi.1]
      pi.hat.0 <- pi.hat[, index.pi.0]
      pi.hat.1 <- truncate_propensity(pi.hat.1, A, treatment_level = treatment_level, truncation_method = "adaptive")
      pi.hat.0 <- truncate_propensity(pi.hat.0, A, treatment_level = control_level, truncation_method = "adaptive")

      index.mu.1 <- match(as.character(treatment_level), as.character(colnames(mu.hat)))
      index.mu.0 <- match(as.character(control_level), as.character(colnames(mu.hat)))
      mu.hat.1 <- mu.hat[, index.mu.1]
      mu.hat.0 <- mu.hat[, index.mu.0]
      # compute DR-learner pseudo-outcome based on EIF
      pseudo_outcome <- mu.hat.1 - mu.hat.0 + (A == treatment_level)/pi.hat.1 * (Y - mu.hat.1)  - (A == control_level)/pi.hat.0 * (Y - mu.hat.0)

      return(list(pseudo_outcome = pseudo_outcome))
    }
  ),
  active = list(

  ),
  private = list(
    .properties = c(
      "cate", "binary treatment", "categorical treatment", "DRlearner"
    )
  )
)

