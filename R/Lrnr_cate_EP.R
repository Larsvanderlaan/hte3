#' Lrnr_cate_EP Class
#'
#' This class defines the EP-learner of van der Laan et al. (2023) for estimation of the conditional average treatment effect.
#'
#' The EP-learner is a robust and doubly-robust meta-learner that inherits desirable properties of both T-learner and DR-learner.
#' @format An R6 class with public methods to initialize the learner, create a regression task, and access the base learner.
#' @inheritParams Lrnr_hte
#' @inheritParams Lrnr_cate_DR
#' @import Sieve
#' @param sieve_num_basis The number of trignometric basis functions to construct the EP-learner sieve space.
#' This argument is passed as the argument \code{basisN} to the function \code{\link[Sieve]{sieve_preprocess}}.
#' By default, \code{sieve_num_basis = ceiling((n^(1/3)*d)} where \code{n} denotes the sample size and \code{d} is the number of effect modifiers.
#' @param sieve_interaction_order The maximum interaction degree of tensor-product basis functions in the EP-learner sieve basis. Default is 3.
#' \code{sieve_interaction_order = 1} corresponds to an additive sieve model and \code{sieve_interaction_order = 2} corresponds to a bi-additive sieve model.
#' This argument is passed as the argument \code{interaction_order} to the function \code{\link[Sieve]{sieve_preprocess}}.
#' @param screen_basis_with_lasso EXPERIMENTAL. Whether to use a lasso-based DR-learner algorithm to screen sieve basis functions for EP-learner.
#' There are no theoretical gaurantees for EP-learner with \code{screen_basis_with_lasso = TRUE}. However, this argument may be useful in high dimensional settings.
#' It also reduces the computational complexity of EP-learner, as the argument \code{sieve_num_basis} does not need to be externally cross-validated.
#' @param ... Additional arguments to pass to the initialization function.
#'
#' @export
Lrnr_cate_EP <- R6Class(
  classname = "Lrnr_cate_EP", inherit = Lrnr_hte,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(base_learner, sieve_num_basis = NULL,  sieve_interaction_order = 3, screen_basis_with_lasso = FALSE, treatment_level = NULL, control_level = NULL, ...) {
      params <- sl3:::args_to_list()
      super$initialize(params = params, base_learner = base_learner,
                       ...)
    },
    get_pseudo_data = function(hte3_task, treatment_level = NULL, control_level = NULL, train = TRUE, ...) {
      # get data, treatment, and outcome
      data <- hte3_task$data
      A <- hte3_task$get_tmle_node("treatment")
      Y <- hte3_task$get_tmle_node("outcome")
      # if control and treatment levels of `A` not provided, use defaults.
      if(is.null(treatment_level)) treatment_level <- levels(factor(A))[2]
      if(is.null(control_level)) control_level <- levels(factor(A))[1]
      # Get matrix of nuisance estimates for propensity score and outcome regression
      pi.hat <- as.matrix(hte3_task$get_nuisance_estimates("pi"))
      mu.hat <- as.matrix(hte3_task$get_nuisance_estimates("mu"))

      # get propensity estimates for relevant treatment levels
      index.pi.1 <- match(as.character(treatment_level), as.character(colnames(pi.hat)))
      index.pi.0 <- match(as.character(control_level), as.character(colnames(pi.hat)))
      pi.hat.1 <- pi.hat[, index.pi.1]
      pi.hat.0 <- pi.hat[, index.pi.0]
      # use data-adaptive propensity score truncation to bound away from 0 and 1.
      pi.hat.1 <- truncate_propensity(pi.hat.1, A, treatment_level = treatment_level, truncation_method = "adaptive")
      pi.hat.0 <- truncate_propensity(pi.hat.0, A, treatment_level = control_level, truncation_method = "adaptive")

      # get outcome regression estimates for relevant treatment levels
      index.mu.1 <- match(as.character(treatment_level), as.character(colnames(mu.hat)))
      index.mu.0 <- match(as.character(control_level), as.character(colnames(mu.hat)))
      mu.hat.1 <- mu.hat[, index.mu.1]
      mu.hat.0 <- mu.hat[, index.mu.0]

      # get treatment effect modifiers to regress on.
      X <- self$get_modifiers(hte3_task, return_matrix = TRUE)


      args <- self$params
      basisN <- self$params$sieve_num_basis
      if(is.null(basisN)) {
        basisN <- ceiling((nrow(X))^(1/3)*ncol(X))
      }
      basisN <- basisN + 1
      interaction_order <- self$params$sieve_interaction_order
      # construct sieve basis using the asymmetric cosine basis.
      sieve_basis <- Sieve::sieve_preprocess(as.matrix(X),
                                             basisN = basisN,
                                             interaction_order = interaction_order,
                                             type = "cosine")$Phi



      # if true, use DR-learner approach to screen basis functions of sieve.
      if(self$params$screen_basis_with_lasso) {
        # construct DR-learner pseudo-outcome
        pseudo_outcome_DR <- mu.hat.1 - mu.hat.0
        + (A == treatment_level)/pi.hat.1 * (Y - mu.hat.1)
        - (A == control_level)/pi.hat.0 * (Y - mu.hat.0)
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
      sieve_basis_train <- sieve_basis * ((A == treatment_level) - (A == control_level))
      mu.hat <- (A == treatment_level) * (mu.hat.1) + (A == control_level) * (mu.hat.0)

      glmnet_fit <- glmnet::glmnet(sieve_basis_train, Y,
                                   offset = mu.hat,
                                   weights = (A == treatment_level)/pi.hat.1 + (A == control_level)/pi.hat.0,
                                   family = "gaussian",
                                   standardize = FALSE,
                                   lambda = 1e-8,
                                   intercept = FALSE,
                                   alpha = 0)
      beta <- as.matrix(glmnet_fit$beta)


      # get debiased outcome regression
      mu.hat.1.star <- mu.hat.1 + sieve_basis %*% beta
      mu.hat.0.star <- mu.hat.0 - sieve_basis %*% beta
      mu.hat.star <- ifelse(A == treatment_level, mu.hat.1.star, mu.hat.0.star)
      #construct EP-learner pseudo-outcome
      pseudo_outcome <- mu.hat.1.star  - mu.hat.0.star

      return(list(pseudo_outcome = pseudo_outcome))
    }
  ),
  active = list(

  ),
  private = list(
    .treatment_type = c("binary_treatment", "categorical_treatment"),
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
  d <- length(hte3_task$npsem$modifiers$variables)
  n <- hte3_task$nrow
  # 3 per, 5 per,
  grid_nbasis <- c(d, 2*d, 4*d, 6*d, 8*d)

  lrnr_cate_epstack <- Stack$new(Lrnr_cate_EP$new(base_learner = base_learner, sieve_num_basis  = d,
                                                  treatment_level = treatment_level, control_level = control_level)
                                 , Lrnr_cate_EP$new(base_learner = base_learner, sieve_num_basis  = 5*d,
                                                    treatment_level = treatment_level, control_level = control_level),
                                 Lrnr_cate_EP$new(base_learner = base_learner, sieve_num_basis = 7*d,
                                                  treatment_level = treatment_level, control_level = control_level),
                                 Lrnr_cate_EP$new(base_learner = base_learner, sieve_num_basis = 9*d,
                                                  treatment_level = treatment_level, control_level = control_level),
                                 Lrnr_cate_EP$new(base_learner = base_learner, sieve_num_basis = 12*d,
                                                  treatment_level = treatment_level, control_level = control_level),
                                 Lrnr_cate_EP$new(base_learner = base_learner, sieve_num_basis = 15*d,
                                                  treatment_level = treatment_level, control_level = control_level),
                                 Lrnr_cate_EP$new(base_learner = base_learner, sieve_num_basis = 20*d,
                                                  treatment_level = treatment_level, control_level = control_level))
  return(lrnr_cate_epstack)
}

