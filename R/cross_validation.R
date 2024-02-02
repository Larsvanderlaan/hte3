
#' Cross-Validate Heterogeneous Treatment Effect Models
#'
#' Cross-validates a collection of heterogeneous treatment effect (hte3) learners using the specified loss function.
#'
#' @param hte_learners A single \code{Lrnr_hte} meta-learner or a \code{Stack} of \code{Lrnr_hte} meta-learners to cross-validate.
#' @param hte3_task An \code{hte3_Task} object containing the data and necessary information for heterogeneous treatment effect estimation.
#' @param cv_metalearner An optional metalearner (\code{Lrnr_base} object) used to combine the cross-validated learners. Default is \code{Lrnr_cv_selector$new(loss_squared_error)}.
#' @param cv_loss_spec A function specifying the loss function to use for cross-validation. Default is \code{loss_spec_cate_DR_binary}.
#' @param cv_control A list of control parameters for cross-validation passed to \code{\link[sl3]{Lrnr_sl}}. Default is \code{NULL}.
#' @param ... Additional arguments to pass to the loss function and other functions.
#'
#' @return A trained \code{Lrnr_sl} object containing the cross-validated ensemble of heterogeneous treatment effect learners.
#'
#' @importFrom data.table data.table
#' @export
cross_validate <- function(hte_learners, hte3_task, cv_metalearner = Lrnr_cv_selector$new(loss_squared_error), cv_control = NULL, ...) {
  args <- list(hte3_task = hte3_task, ...)
  if(is.list(hte_learners)) {
    hte_learners <- Stack$new(hte_learners)
  }
  lrnr_sl <- Lrnr_sl$new(learners = hte_learners, metalearner = cv_metalearner, cv_control = cv_control, ...)$train(hte3_task)
  return(lrnr_sl)
}




#' Cross-Validate CATE Models
#'
#' Cross-validates a collection of CATE (hte3) learners using the DR-learner loss.
#'
#' @inheritParams cross_validate
#' @export
cross_validate_cate <- function(hte_learners, hte3_task, cv_control = NULL, treatment_level = NULL, control_level = NULL, ...) {
  cv_metalearner <- Lrnr_cate_DR_selector$new(treatment_level = treatment_level, control_level = control_level)
  lrnr_sl <- cross_validate(hte_learners, hte3_task, cv_metalearner = cv_metalearner, cv_control = cv_control)
  out <- list(lrnr_sl = lrnr_sl, cv_risk = lrnr_sl$fit_object$cv_meta_fit$cv_risk, coefficients = lrnr_sl$fit_object$cv_meta_fit$coefficients)
  return(out)
}



#' Cross-Validate CRR Models
#'
#' Cross-validates a collection of CRR (hte3) learners using the DR-learner loss.
#'
#' @inheritParams cross_validate
#' @export
cross_validate_crr <- function(hte_learners, hte3_task, cv_control = NULL, treatment_level = NULL, control_level = NULL, ...) {
  cv_metalearner <- Lrnr_crr_DR_selector$new(treatment_level = treatment_level, control_level = control_level)
  lrnr_sl <- cross_validate(hte_learners, hte3_task, cv_metalearner = cv_metalearner, cv_control = cv_control)
  out <- list(lrnr_sl = lrnr_sl, cv_risk = lrnr_sl$fit_object$cv_meta_fit$cv_risk, coefficients = lrnr_sl$fit_object$cv_meta_fit$coefficients)
  return(out)
}






# cross_validate <- function(hte_learners, hte3_task, cv_metalearner = Lrnr_cv_selector$new(loss_squared_error), cv_loss_spec = loss_spec_cate_DR_binary, cv_control = NULL, ...) {
#   args <- list(hte3_task = hte3_task, ...)
#   if(is.list(hte_learners)) {
#     hte_learners <- Stack$new(hte_learners)
#   }
#   pseudo_data <- sl3:::call_with_args(cv_loss_spec, args, silent = TRUE)
#   if(is.character(pseudo_data$family)) pseudo_data$family <- get(pseudo_data$family)
#   pseudo_outcome <- pseudo_data$pseudo_outcome
#   pseudo_weights <- pseudo_data$pseudo_weights
#   if(is.null(pseudo_weights)) {
#     pseudo_weights <- hte3_task$weights
#   } else {
#     pseudo_weights <- hte3_task$weights * pseudo_weights
#   }
#   new_data <- data.table(pseudo_outcome, pseudo_weights)
#   names(new_data) <- c("pseudo_outcome", "pseudo_weights")
#   column_names <- hte3_task$add_columns(new_data)
#   if(pseudo_data$family$family == "gaussian") {
#     new_outcome_type <- "continuous"
#   } else if(pseudo_data$family$family == "binomial") {
#     new_outcome_type <- "binomial"
#   }
#   new_hte3_task<- hte3_task$next_in_chain(covariates = c(),
#                                           outcome = "pseudo_outcome",
#                                           weights = "pseudo_weights",
#                                           column_names = column_names,
#                                           new_outcome_type = new_outcome_type)
#   lrnr_sl <- Lrnr_sl$new(learners = hte_learners, cv_metalearner = cv_metalearner, cv_control = cv_control, ...)$train(new_hte3_task)
#   return(lrnr_sl)
# }


