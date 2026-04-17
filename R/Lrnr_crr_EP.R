#' Lrnr_crr_EP Class
#'
#' This class defines an efficient plug-in (EP) meta-learner for the conditional relative average treatment effect (crr).
#'
#' @format An R6 class with public methods to initialize the learner, create a regression task, and access the base learner.
#' @param base_learner A \code{\link{sl3}} learner object inheriting from \code{\link[sl3]{Lrnr_base}} that specifies the base supervised learning algorithm used by the meta-learner.
#' @param sieve_num_basis The number of trigonometric basis functions used to
#'   construct the CRR EP sieve space through
#'   \code{\link[Sieve]{sieve_preprocess}}. Use \code{sieve_num_basis} for one
#'   fixed EP fit and \code{sieve_basis_grid} in \code{fit_crr()} for
#'   wrapper-level EP basis-size selection.
#' @param sieve_interaction_order The maximum interaction degree of
#'   tensor-product basis functions in the EP sieve basis.
#' @param treatment_level Optional treated level used for the contrast.
#' @param control_level Optional control level used for the contrast.
#' @import Sieve
#' @export
Lrnr_crr_EP <- R6Class(
  classname = "Lrnr_crr_EP",
  inherit = Lrnr_hte,
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(
    base_learner,
    sieve_num_basis = NULL,
    sieve_interaction_order = 3,
    treatment_level = NULL,
    control_level = NULL,
    ...
    ) {
      params <- list(
        base_learner = base_learner,
        sieve_num_basis = sieve_num_basis,
        sieve_interaction_order = sieve_interaction_order,
        treatment_level = treatment_level,
        control_level = control_level,
        ...
      )
      super$initialize(params = params, base_learner = base_learner,
                       transform_function = bounded_qlogis,
                       pseudo_outcome_type = "quasibinomial",
                       pseudo_family = quasibinomial(),
                       ...)
    },
    get_pseudo_data = function(hte3_task, treatment_level = NULL, control_level = NULL, train = TRUE, ...) {
      nuisance <- resolve_contrast_nuisances(
        hte3_task,
        treatment_level = treatment_level,
        control_level = control_level
      )
      validate_nonnegative_outcome(nuisance$Y, label = "CRR outcomes")
      X <- self$get_modifiers(hte3_task, return_matrix = TRUE)
      sieve_basis <- make_sieve_basis(
        X,
        basisN = self$params$sieve_num_basis,
        interaction_order = self$params$sieve_interaction_order
      )

      weights_invpi <- (nuisance$A == nuisance$treatment_level) / nuisance$pi.hat.1 +
        (nuisance$A == nuisance$control_level) / nuisance$pi.hat.0

      # if(self$params$screen_basis_with_lasso) {
      #
      #   # use DR-learner to find basis functions predictive of CATE
      #   glmnet_fit <- glmnet::cv.glmnet(sieve_basis[Y > 0,], as.numeric(A == treatment_level)[Y > 0],
      #                                   weights = (weights_invpi * hte3_task$weights)[Y > 0],
      #                                   family = "binomial",
      #                                   standardize = FALSE,
      #                                   intercept = FALSE,
      #                                   alpha = 1)
      #
      #   keep <- union(c(1), which(coef(glmnet_fit, s = "lambda.min")[-1] != 0))
      #   sieve_basis <- sieve_basis[, keep, drop = FALSE]
      # }
      sieve_basis_train <- cbind(
        sieve_basis * (nuisance$A == nuisance$treatment_level),
        sieve_basis * (nuisance$A == nuisance$control_level)
      )
      mu.hat <- (nuisance$A == nuisance$treatment_level) * nuisance$mu.hat.1 +
        (nuisance$A == nuisance$control_level) * nuisance$mu.hat.0

      # Determine link function for debiasing step of outcome regression
      if(all(nuisance$Y %in% c(0,1))) {
        glmnet_family <- "binomial"
        family <- binomial()
      } else if(all(nuisance$Y >= 0 & nuisance$Y <= 1)) {
        glmnet_family <- binomial()
        family <- binomial()
      } else if(all(nuisance$Y >= 0)) {
        glmnet_family <- "poisson"
        family <- poisson()
      } else {
        stop("Outcome must be nonnegative for crr estimation.")
      }

      if (identical(family$family, "binomial")) {
        mu.hat <- bound_probability(mu.hat)
        nuisance$mu.hat.1 <- bound_probability(nuisance$mu.hat.1)
        nuisance$mu.hat.0 <- bound_probability(nuisance$mu.hat.0)
      } else {
        mu.hat <- bound_positive(mu.hat)
        nuisance$mu.hat.1 <- bound_positive(nuisance$mu.hat.1)
        nuisance$mu.hat.0 <- bound_positive(nuisance$mu.hat.0)
      }

      # perform EP-learner debiasing algorithm of outcome regression nuisance
      glmnet_fit <- glmnet::glmnet(sieve_basis_train, nuisance$Y,
                                   offset = family$linkfun(mu.hat),
                                   weights = weights_invpi,
                                   family = glmnet_family,
                                   standardize = FALSE,
                                   lambda = 1e-8,
                                   intercept = FALSE,
                                   alpha = 0)
      beta <- as.matrix(glmnet_fit$beta)

      # get sieve-debiased estimates.
      sieve_basis_1 <- cbind(sieve_basis * 1, sieve_basis * 0)
      sieve_basis_0 <- cbind(sieve_basis * 0, sieve_basis * 1)
      mu.hat.1.star <- family$linkinv(family$linkfun(nuisance$mu.hat.1) + sieve_basis_1 %*% beta)
      mu.hat.0.star <- family$linkinv(family$linkfun(nuisance$mu.hat.0) + sieve_basis_0 %*% beta)
      mu.hat.1.star <- if (identical(family$family, "binomial")) bound_probability(mu.hat.1.star) else bound_positive(mu.hat.1.star)
      mu.hat.0.star <- if (identical(family$family, "binomial")) bound_probability(mu.hat.0.star) else bound_positive(mu.hat.0.star)
      # get pseudo outcomes and weights for logistic regression.
      pseudo_weights <- bound_positive(mu.hat.1.star + mu.hat.0.star)
      pseudo_outcome <- bound_probability(mu.hat.1.star / pseudo_weights)
      return(list(pseudo_outcome = pseudo_outcome, pseudo_weights = pseudo_weights))
    }
  ),
  active = list(

  ),
  private = list(
    .treatment_type = c("binomial", "categorical"),
    .properties = c(
      "crr", "EP"
    )
  )
)
