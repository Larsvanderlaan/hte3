compute_cate_ep_dr_pseudo_outcome <- function(hte3_task, nuisance, modifier_matrix, params) {
  sieve_basis <- make_sieve_basis(
    modifier_matrix,
    basisN = params$sieve_num_basis,
    interaction_order = params$sieve_interaction_order
  )

  if (isTRUE(params$screen_basis_with_lasso)) {
    pseudo_outcome_DR <- compute_cate_dr_pseudo_outcome(
      hte3_task,
      treatment_level = nuisance$treatment_level,
      control_level = nuisance$control_level
    )
    glmnet_fit <- glmnet::cv.glmnet(
      sieve_basis,
      pseudo_outcome_DR,
      family = "gaussian",
      standardize = FALSE,
      intercept = FALSE,
      alpha = 1
    )
    keep <- union(c(1), which(coef(glmnet_fit, s = "lambda.min")[-1] != 0))
    sieve_basis <- sieve_basis[, keep, drop = FALSE]
  }

  sieve_basis_train <- sieve_basis * ((nuisance$A == nuisance$treatment_level) - (nuisance$A == nuisance$control_level))
  mu.hat <- (nuisance$A == nuisance$treatment_level) * nuisance$mu.hat.1 +
    (nuisance$A == nuisance$control_level) * nuisance$mu.hat.0

  glmnet_fit <- glmnet::glmnet(
    sieve_basis_train,
    nuisance$Y,
    offset = mu.hat,
    weights = (nuisance$A == nuisance$treatment_level) / nuisance$pi.hat.1 +
      (nuisance$A == nuisance$control_level) / nuisance$pi.hat.0,
    family = "gaussian",
    standardize = FALSE,
    lambda = 1e-8,
    intercept = FALSE,
    alpha = 0
  )
  beta <- as.matrix(glmnet_fit$beta)

  mu.hat.1.star <- nuisance$mu.hat.1 + sieve_basis %*% beta
  mu.hat.0.star <- nuisance$mu.hat.0 - sieve_basis %*% beta
  list(pseudo_outcome = as.vector(mu.hat.1.star - mu.hat.0.star))
}

compute_cate_ep_r_pseudo_outcome <- function(hte3_task, nuisance, params) {
  if (!identical(get_task_treatment_type(hte3_task), "binomial")) {
    stop(
      "`targeting_style = \"r\"` is currently only supported for binary-treatment CATE tasks.",
      call. = FALSE
    )
  }

  warn_ep_r_reduced_modifier_target(hte3_task, r_targeting_basis = params$r_targeting_basis)

  D <- as.numeric(nuisance$A == nuisance$treatment_level)
  e.hat <- nuisance$pi.hat.1
  residual_treatment <- D - e.hat
  safe_residual <- stabilize_treatment_residual(residual_treatment)
  first_stage_covariates <- get_ep_r_first_stage_covariates(
    hte3_task,
    e.hat = e.hat,
    r_targeting_basis = params$r_targeting_basis
  )
  interaction_order <- effective_ep_interaction_order(
    targeting_style = params$targeting_style,
    r_targeting_basis = params$r_targeting_basis,
    interaction_order = params$sieve_interaction_order
  )
  sieve_basis <- make_sieve_basis(
    first_stage_covariates,
    basisN = params$sieve_num_basis,
    interaction_order = interaction_order
  )
  tau0.hat <- nuisance$mu.hat.1 - nuisance$mu.hat.0

  glmnet_fit <- glmnet::glmnet(
    sieve_basis,
    (nuisance$Y - nuisance$m.hat) / safe_residual,
    offset = tau0.hat,
    weights = bound_positive(residual_treatment^2),
    family = "gaussian",
    standardize = FALSE,
    lambda = 1e-8,
    intercept = FALSE,
    alpha = 0
  )
  beta <- as.matrix(glmnet_fit$beta)
  tau.hat <- tau0.hat + sieve_basis %*% beta
  list(pseudo_outcome = as.vector(tau.hat))
}

#' Lrnr_cate_EP Class
#'
#' This class defines the EP-learner of van der Laan et al. (2023) for estimation of the conditional average treatment effect.
#'
#' The EP-learner is a robust and doubly-robust meta-learner that inherits desirable properties of both T-learner and DR-learner.
#' The current implementation exposes two CATE EP variants through
#' \code{targeting_style}: \code{"dr"} is the plug-in analogue of the
#' DR-learner and \code{"r"} is the plug-in analogue of the R-learner.
#' The \code{"r"} variant uses an overlap-weighted targeting objective that can
#' be more stable under severe overlap violations, but it generally targets an
#' overlap-weighted projection in reduced-modifier settings.
#' @format An R6 class with public methods to initialize the learner, create a regression task, and access the base learner.
#' @import Sieve
#' @param base_learner A \code{\link{sl3}} learner object inheriting from \code{\link[sl3]{Lrnr_base}} that specifies the base supervised learning algorithm used by the meta-learner.
#' @param sieve_num_basis The number of trigonometric basis functions to construct the EP-learner sieve space.
#' This argument is passed as the argument \code{basisN} to the function \code{\link[Sieve]{sieve_preprocess}}.
#' By default, \code{sieve_num_basis = ceiling((n^(1/3)*d)} where \code{n}
#' denotes the sample size and \code{d} is the dimension of the first-stage EP
#' sieve covariates.
#' For wrapper-level EP cross-validation, use \code{sieve_basis_grid} in
#' \code{fit_cate()} or \code{fit_crr()} to compare multiple basis sizes.
#' @param sieve_interaction_order The maximum interaction degree of tensor-product basis functions in the EP-learner sieve basis. Default is 3.
#' \code{sieve_interaction_order = 1} corresponds to an additive sieve model and \code{sieve_interaction_order = 2} corresponds to a bi-additive sieve model.
#' This argument is passed as the argument \code{interaction_order} to the function \code{\link[Sieve]{sieve_preprocess}}.
#' @param screen_basis_with_lasso EXPERIMENTAL. Whether to use a lasso-based DR-learner algorithm to screen sieve basis functions for EP-learner.
#' There are no theoretical guarantees for EP-learner with \code{screen_basis_with_lasso = TRUE}. However, this argument may be useful in high dimensional settings.
#' It also reduces the computational complexity of EP-learner, as the argument \code{sieve_num_basis} does not need to be externally cross-validated.
#' @param targeting_style One of \code{"dr"} or \code{"r"}. The default
#'   \code{"dr"} path keeps the current EP debiasing update. The
#'   \code{"r"} path uses an R-learner-style overlap-weighted targeting update
#'   and currently supports binary-treatment CATE tasks only. High-level
#'   wrappers such as \code{fit_cate()} compare both variants by expanding
#'   multiple \code{Lrnr_cate_EP} fits rather than by passing a vector here.
#' @param r_targeting_basis First-stage basis construction for
#'   \code{targeting_style = "r"}. The default
#'   \code{"v_plus_propensity"} builds the first-stage sieve on the modifier
#'   set \code{V} plus the treated-arm propensity score \code{e(W)};
#'   \code{"full_w"} builds the first-stage EP-R sieve on the full
#'   confounder set \code{W}.
#' @param ... Additional arguments to pass to the initialization function.
#'
#' @export

Lrnr_cate_EP <- R6Class(
  classname = "Lrnr_cate_EP", inherit = Lrnr_hte,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(base_learner,
                          sieve_num_basis = NULL,
                          sieve_interaction_order = 3,
                          screen_basis_with_lasso = FALSE,
                          targeting_style = c("dr", "r"),
                          r_targeting_basis = c("v_plus_propensity", "full_w"),
                          treatment_level = NULL,
                          control_level = NULL,
                          ...) {
      targeting_style <- match.arg(targeting_style)
      r_targeting_basis <- match.arg(r_targeting_basis)
      if (identical(targeting_style, "r") && isTRUE(screen_basis_with_lasso)) {
        stop(
          "`screen_basis_with_lasso = TRUE` is only supported for `targeting_style = \"dr\"`.",
          call. = FALSE
        )
      }
      params <- list(
        base_learner = base_learner,
        sieve_num_basis = sieve_num_basis,
        sieve_interaction_order = sieve_interaction_order,
        screen_basis_with_lasso = screen_basis_with_lasso,
        targeting_style = targeting_style,
        r_targeting_basis = r_targeting_basis,
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

      if (identical(self$params$targeting_style, "r")) {
        return(compute_cate_ep_r_pseudo_outcome(
          hte3_task = hte3_task,
          nuisance = nuisance,
          params = self$params
        ))
      }

      modifier_matrix <- self$get_modifiers(hte3_task, return_matrix = TRUE)
      compute_cate_ep_dr_pseudo_outcome(
        hte3_task = hte3_task,
        nuisance = nuisance,
        modifier_matrix = modifier_matrix,
        params = self$params
      )
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
#' @param treatment_level Treatment level used for the treated arm in each EP learner in the stack.
#' @param control_level Reference treatment level used for the control arm in each EP learner in the stack.
#' @param sieve_basis_grid Optional vector of sieve basis sizes to include in
#'   the ensemble. If \code{NULL}, the default heuristic grid uses the
#'   first-stage EP sieve dimension for the chosen targeting style.
#' @param targeting_style EP targeting style or styles to include in the stack.
#'   Use \code{c("dr", "r")} to build both EP variants across the requested
#'   basis grid.
#' @param r_targeting_basis First-stage basis construction used for EP-R fits
#'   in the stack.
#' @return A stack of CATE EP meta-learners with varying sieve basis sizes.
#' @examples
#' \dontrun{
#' library(sl3)
#'
#' data <- hte3_example_data(n = 80, seed = 1)
#' task <- hte_task(
#'   data = data,
#'   modifiers = c("W1", "W2"),
#'   confounders = c("W1", "W2", "W3"),
#'   treatment = "A",
#'   outcome = "Y",
#'   propensity_learner = Lrnr_mean$new(),
#'   outcome_learner = Lrnr_mean$new(),
#'   mean_learner = Lrnr_mean$new(),
#'   cross_fit = FALSE
#' )
#'
#' ep_stack <- make_ep_stack(
#'   base_learner = Lrnr_mean$new(),
#'   hte3_task = task,
#'   sieve_basis_grid = c(2, 4, 6),
#'   targeting_style = c("dr", "r")
#' )
#'
#' ep_stack
#' }
#' @export
make_ep_stack <- function(base_learner,
                          hte3_task,
                          treatment_level = 1,
                          control_level = 0,
                          sieve_basis_grid = NULL,
                          targeting_style = "dr",
                          r_targeting_basis = "v_plus_propensity") {
  do.call(
    Stack$new,
    make_cate_ep_candidates(
      base_learner,
      hte3_task,
      treatment_level,
      control_level,
      sieve_basis_grid = sieve_basis_grid,
      targeting_style = targeting_style,
      r_targeting_basis = r_targeting_basis
    )
  )
}
