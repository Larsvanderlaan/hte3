#' Lrnr_crr_EP Class
#'
#' This class defines an efficient plug-in (EP) meta-learner for the conditional relative average treatment effect (crr).
#'
#' @format An R6 class with public methods to initialize the learner, create a regression task, and access the base learner.
#' @param base_learner A \code{\link{sl3}} learner object inheriting from \code{\link[sl3]{Lrnr_base}} that specifies the base supervised learning algorithm used by the meta-learner.
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
      params <- sl3:::args_to_list()
      super$initialize(params = params, base_learner = base_learner,
                       transform_function = stats::qlogis,
                       pseudo_outcome_type = "quasibinomial",
                       pseudo_family = quasibinomial(),
                       ...)
    },
    get_pseudo_data = function(hte3_task, treatment_level = NULL, control_level = NULL, train = TRUE, ...) {
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
      pi.hat.1 <- causalutils::truncate_propensity(pi.hat.1, A, treatment_level = treatment_level, truncation_method = "adaptive")
      pi.hat.0 <- causalutils::truncate_propensity(pi.hat.0, A, treatment_level = control_level, truncation_method = "adaptive")

      index.mu.1 <- match(as.character(treatment_level), as.character(colnames(mu.hat)))
      index.mu.0 <- match(as.character(control_level), as.character(colnames(mu.hat)))
      mu.hat.1 <- mu.hat[, index.mu.1]
      mu.hat.0 <- mu.hat[, index.mu.0]

      #pseudo_outcome <- mu.hat.1 - mu.hat.0
      # + (A == treatment_level)/pi.hat.1 * (Y - mu.hat.1)
      #- (A == control_level)/pi.hat.0 * (Y - mu.hat.0)

      X <- self$get_modifiers(hte3_task, return_matrix = TRUE)
      args <- self$params
      args$X <- X
      basisN <- self$params$sieve_num_basis
      if(is.null(basisN)) {
        basisN <- ceiling((nrow(X))^(1/3)*ncol(X))
      }
      basisN <- basisN + 1
      interaction_order <- self$params$sieve_interaction_order
      sieve_basis <- Sieve::sieve_preprocess(as.matrix(X),
                                             basisN = basisN,
                                             interaction_order = interaction_order,
                                             type = "cosine")$Phi



      weights_invpi <- (A == treatment_level)/pi.hat.1 +  (A == control_level)/pi.hat.0

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
      print(dim(sieve_basis))

      sieve_basis_train <- cbind(sieve_basis * (A == treatment_level), sieve_basis * (A == control_level))
      mu.hat <- (A == treatment_level) * (mu.hat.1) + (A == control_level) * (mu.hat.0)

      # Determine link function for debiasing step of outcome regression
      if(all(Y %in% c(0,1))) {
        glmnet_family <- "binomial"
        family <- binomial()
      } else if(all(Y >= 0 & Y <= 1)) {
        glmnet_family <- binomial()
        family <- binomial()
      } else if(all(Y >= 0)) {
        glmnet_family <- "poisson"
        family <- poisson()
      } else {
        stop("Outcome must be nonnegative for crr estimation.")
      }
      # perform EP-learner debiasing algorithm of outcome regression nuisance
      glmnet_fit <- glmnet::glmnet(sieve_basis_train, Y,
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
      mu.hat.1.star <- family$linkinv(family$linkfun(mu.hat.1) + sieve_basis_1 %*% beta)
      mu.hat.0.star <- family$linkinv(family$linkfun(mu.hat.0) + sieve_basis_0 %*% beta)
      # get pseudo outcomes and weights for logistic regression.
      pseudo_outcome <- mu.hat.1.star  /  (mu.hat.1.star + mu.hat.0.star)
      pseudo_weights <- (mu.hat.1.star + mu.hat.0.star) # note, observation weights are automatically added to the task and should not be added to pseudo_weights
      return(list(pseudo_outcome = pseudo_outcome, pseudo_weights = pseudo_weights))
    }
  ),
  active = list(

  ),
  private = list(
    .treatment_type = c("binary_treatment", "categorical_treatment"),
    .properties = c(
      "crr", "EP"
    )
  )
)
