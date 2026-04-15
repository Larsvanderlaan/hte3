#' Lrnr_cate_DR Class
#'
#' This class defines a doubly-robust (DR) meta-learner of the conditional average treatment effect.
#'
#' @format An R6 class with public methods to initialize the learner, create a regression task, and access the base learner.
#' @param base_learner A \code{\link{sl3}} learner object inheriting from \code{\link[sl3]{Lrnr_base}} that specifies the base supervised learning algorithm used by the meta-learner.
#' @param treatment_level A treatment level encoding the treatment assignment of interest.
#' @param control_level A treatment level encoding the control (or reference) treatment assignment.
#' @export
Lrnr_cate_DR <- R6Class(
  classname = "Lrnr_cate_DR", inherit = Lrnr_hte,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(base_learner, treatment_level = NULL, control_level = NULL, ...) {
      params <- list(base_learner = base_learner, treatment_level = treatment_level, control_level = control_level, ...)
      super$initialize(params = params, base_learner = base_learner,
                       ...)
    },
    get_pseudo_data = function(hte3_task, treatment_level = NULL, control_level = NULL, ...) {
      list(
        pseudo_outcome = compute_cate_dr_pseudo_outcome(
          hte3_task,
          treatment_level = treatment_level,
          control_level = control_level
        )
      )
    }
  ),
  active = list(

  ),
  private = list(
    .treatment_type = c("binomial", "categorical"),
    .properties = c(
      "cate", "binary treatment", "categorical treatment", "DRlearner"
    )
  )
)
