#' Lrnr_cate_EP Class
#'
#' This class defines the EP-learner of van der Laan et al. (2023) for estimation of the conditional average treatment effect.
#'
#' The EP-learner is a robust and doubly-robust meta-learner that inherits desirable properties of both T-learner and DR-learner.
#' In the supported binary/categorical-treatment setting, it targets the
#' conditional mean difference over the chosen modifier set `V`, namely
#' `E[Y(1) - Y(0) | V]`.
#' @format An R6 class with public methods to initialize the learner, create a regression task, and access the base learner.
#' @import Sieve
#' @param base_learner A \code{\link{sl3}} learner object inheriting from \code{\link[sl3]{Lrnr_base}} that specifies the base supervised learning algorithm used by the meta-learner.
#' @param sieve_num_basis The number of trigonometric basis functions to construct the EP-learner sieve space.
#' This argument is passed as the argument \code{basisN} to the function \code{\link[Sieve]{sieve_preprocess}}.
#' By default, \code{sieve_num_basis = ceiling((n^(1/3)*d)} where \code{n} denotes the sample size and \code{d} is the number of effect modifiers.
#' @param sieve_interaction_order The maximum interaction degree of tensor-product basis functions in the EP-learner sieve basis. Default is 3.
#' \code{sieve_interaction_order = 1} corresponds to an additive sieve model and \code{sieve_interaction_order = 2} corresponds to a bi-additive sieve model.
#' This argument is passed as the argument \code{interaction_order} to the function \code{\link[Sieve]{sieve_preprocess}}.
#' @param screen_basis_with_lasso EXPERIMENTAL. Whether to use a lasso-based DR-learner algorithm to screen sieve basis functions for EP-learner.
#' There are no theoretical guarantees for EP-learner with \code{screen_basis_with_lasso = TRUE}. However, this argument may be useful in high dimensional settings.
#' It also reduces the computational complexity of EP-learner, as the argument \code{sieve_num_basis} does not need to be externally cross-validated.
#' @param ... Additional arguments to pass to the initialization function.
#'
#' @export
Lrnr_cate_EP <- R6Class(
  classname = "Lrnr_cate_EP", inherit = Lrnr_hte,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(base_learner, sieve_num_basis = NULL,  sieve_interaction_order = 3, screen_basis_with_lasso = FALSE, treatment_level = NULL, control_level = NULL, ...) {
      params <- list(
        base_learner = base_learner,
        sieve_num_basis = sieve_num_basis,
        sieve_interaction_order = sieve_interaction_order,
        screen_basis_with_lasso = screen_basis_with_lasso,
        treatment_level = treatment_level,
        control_level = control_level,
        ...
      )
      super$initialize(params = params, base_learner = base_learner,
                       ...)
    },
    get_pseudo_data = function(hte3_task, treatment_level = NULL, control_level = NULL, train = TRUE, ...) {
      nuisance <- resolve_contrast_nuisances(
        hte3_task,
        treatment_level = treatment_level,
        control_level = control_level
      )

      # get treatment effect modifiers to regress on.
      X <- self$get_modifiers(hte3_task, return_matrix = TRUE)
      sieve_basis <- make_sieve_basis(
        X,
        basisN = self$params$sieve_num_basis,
        interaction_order = self$params$sieve_interaction_order
      )

      # if true, use DR-learner approach to screen basis functions of sieve.
      if(self$params$screen_basis_with_lasso) {
        # construct DR-learner pseudo-outcome
        pseudo_outcome_DR <- compute_cate_dr_pseudo_outcome(
          hte3_task,
          treatment_level = nuisance$treatment_level,
          control_level = nuisance$control_level
        )
        # use lasso DR-learner to find basis functions predictive of CATE
        glmnet_fit <- glmnet::cv.glmnet(sieve_basis, pseudo_outcome_DR,
                                        family = "gaussian",
                                        standardize = FALSE,
                                        intercept = FALSE,
                                        alpha = 1)
        # extract columns with nonzero coefficients
        keep <- union(c(1), which(coef(glmnet_fit, s = "lambda.min")[-1] != 0))
        sieve_basis <- sieve_basis[, keep, drop = FALSE]
      }

      # perform EP-learner sieve-based debiasing of outcome regression nuisance.
      sieve_basis_train <- sieve_basis * ((nuisance$A == nuisance$treatment_level) - (nuisance$A == nuisance$control_level))
      mu.hat <- (nuisance$A == nuisance$treatment_level) * nuisance$mu.hat.1 +
        (nuisance$A == nuisance$control_level) * nuisance$mu.hat.0

      glmnet_fit <- glmnet::glmnet(sieve_basis_train, nuisance$Y,
                                   offset = mu.hat,
                                   weights = (nuisance$A == nuisance$treatment_level) / nuisance$pi.hat.1 +
                                     (nuisance$A == nuisance$control_level) / nuisance$pi.hat.0,
                                   family = "gaussian",
                                   standardize = FALSE,
                                   lambda = 1e-8,
                                   intercept = FALSE,
                                   alpha = 0)
      beta <- as.matrix(glmnet_fit$beta)


      # get debiased outcome regression
      mu.hat.1.star <- nuisance$mu.hat.1 + sieve_basis %*% beta
      mu.hat.0.star <- nuisance$mu.hat.0 - sieve_basis %*% beta
      # construct EP-learner pseudo-outcome
      pseudo_outcome <- mu.hat.1.star  - mu.hat.0.star

      return(list(pseudo_outcome = pseudo_outcome))
    }
  ),
  active = list(

  ),
  private = list(
    .treatment_type = c("binomial", "categorical"),
    .properties = c(
      "cate", "EP"
    )
  )
)


#' Create an Ensemble of CATE EP-learners with Varying Sieve Basis Sizes
#'
#' This function creates an ensemble of CATE EP-learners of varying sieve dimensions for use with cross-validation.
#'
#' @param base_learner The base learner of \code{Lrnr_cate_EP}.
#' @param hte3_task A \code{hte3_Task} object containing the data and necessary information for the heterogeneous treatment effect estimation.
#' @inheritParams Lrnr_cate_DR
#' @return A stack of CATE EP meta-learners with varying sieve basis sizes.
#' @export
make_ep_stack <- function(base_learner, hte3_task, treatment_level = 1, control_level = 0) {
  do.call(Stack$new, make_cate_ep_candidates(base_learner, hte3_task, treatment_level, control_level))
}
