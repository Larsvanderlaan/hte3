describe_hte_candidate <- function(learner) {
  method <- NA_character_
  label <- class(learner)[1]
  sieve_num_basis <- NULL
  ep_targeting_style <- NULL
  ep_r_targeting_basis <- NULL

  if (inherits(learner, "Lrnr_cate_DR")) {
    method <- "dr"
    label <- "dr"
  } else if (inherits(learner, "Lrnr_cate_R")) {
    method <- "r"
    label <- "r"
  } else if (inherits(learner, "Lrnr_grf_causal_forest")) {
    method <- "r"
    label <- "r"
  } else if (inherits(learner, "Lrnr_cate_T") || inherits(learner, "Lrnr_crr_T")) {
    method <- "t"
    label <- "t"
  } else if (inherits(learner, "Lrnr_crr_IPW")) {
    method <- "ipw"
    label <- "ipw"
  } else if (inherits(learner, "Lrnr_cate_EP") || inherits(learner, "Lrnr_crr_EP")) {
    method <- "ep"
    ep_metadata <- extract_ep_learner_metadata(learner)
    sieve_num_basis <- ep_metadata$sieve_num_basis
    interaction_order <- ep_metadata$interaction_order
    basis_label <- if (is.null(sieve_num_basis)) "default" else sieve_num_basis
    if (inherits(learner, "Lrnr_cate_EP")) {
      ep_targeting_style <- ep_metadata$targeting_style
      ep_r_targeting_basis <- ep_metadata$r_targeting_basis
      label <- sprintf(
        "ep[style=%s, basis=%s, interaction=%s",
        ep_targeting_style,
        basis_label,
        interaction_order
      )
      if (!is.null(ep_r_targeting_basis)) {
        label <- sprintf("%s, r_basis=%s", label, ep_r_targeting_basis)
      }
      label <- paste0(label, "]")
    } else {
      label <- sprintf("ep[basis=%s, interaction=%s]", basis_label, interaction_order)
    }
    if (isTRUE(ep_metadata$screen_basis_with_lasso)) {
      label <- sprintf("%s, lasso-screen=TRUE", sub("\\]$", "", label))
      label <- paste0(label, "]")
    }
  }

  list(
    label = label,
    method = method,
    sieve_num_basis = sieve_num_basis,
    ep_targeting_style = ep_targeting_style,
    ep_r_targeting_basis = ep_r_targeting_basis
  )
}

prepare_cv_learner_stack <- function(hte_learners) {
  keep_original_stack <- inherits(hte_learners, "Stack") && length(hte_learners$learners) > 1L

  if (inherits(hte_learners, "Stack")) {
    learners <- hte_learners$learners
  } else if (is.list(hte_learners)) {
    learners <- hte_learners
  } else {
    learners <- list(hte_learners)
  }

  descriptors <- lapply(learners, describe_hte_candidate)

  if (length(learners) == 1L) {
    learners[[2L]] <- learners[[1L]]$clone(deep = TRUE)
    descriptors[[2L]] <- descriptors[[1L]]
  }

  stack <- if (keep_original_stack) {
    hte_learners
  } else {
    do.call(Stack$new, learners)
  }

  list(stack = stack, descriptors = descriptors)
}

#' Cross-Validate Heterogeneous Treatment Effect Models
#'
#' Cross-validates a collection of heterogeneous treatment effect learners.
#'
#' @param hte_learners A single \code{Lrnr_hte} learner, a list of
#'   \code{Lrnr_hte} learners, or a \code{Stack} of \code{Lrnr_hte} learners
#'   to cross-validate.
#' @param hte3_task An \code{hte3_Task} object containing the data and necessary information for heterogeneous treatment effect estimation.
#' @param cv_metalearner An optional metalearner (\code{Lrnr_base} object) used to combine the cross-validated learners. Default is \code{Lrnr_cv_selector$new(loss_squared_error)}.
#' @param cv_control A list of control parameters for cross-validation passed to \code{\link[sl3]{Lrnr_sl}}. Default is \code{NULL}.
#' @param ... Additional arguments to pass to the loss function and other functions.
#'
#' @return A trained \code{Lrnr_sl} object containing the cross-validated ensemble of heterogeneous treatment effect learners.
#'
#' @importFrom data.table data.table
#' @export
cross_validate <- function(hte_learners, hte3_task, cv_metalearner = Lrnr_cv_selector$new(loss_squared_error), cv_control = NULL, ...) {
  prepared <- prepare_cv_learner_stack(hte_learners)
  lrnr_sl <- Lrnr_sl$new(learners = prepared$stack, metalearner = cv_metalearner, cv_control = cv_control, ...)$train(hte3_task)
  attr(lrnr_sl, "hte3_candidate_descriptors") <- prepared$descriptors
  return(lrnr_sl)
}

build_cv_selection_summary <- function(lrnr_sl) {
  descriptors <- attr(lrnr_sl, "hte3_candidate_descriptors")
  coefficients <- lrnr_sl$fit_object$cv_meta_fit$coefficients
  coefficients <- as.numeric(coefficients)
  cv_risk <- lrnr_sl$fit_object$cv_meta_fit$cv_risk

  if (is.null(descriptors) || length(descriptors) != length(coefficients)) {
    candidate_labels <- paste0("candidate_", seq_along(coefficients))
    named_coefficients <- stats::setNames(coefficients, candidate_labels)
    selected_candidate <- if (length(named_coefficients) > 0L) names(named_coefficients)[[which.max(named_coefficients)]] else NULL
    return(list(
      candidate_labels = candidate_labels,
      coefficients = named_coefficients,
      cv_risk = cv_risk,
      selected_candidate = selected_candidate,
      selected_method = selected_candidate,
      selected_sieve_num_basis = NULL,
      ep_basis_grid = NULL,
      selected_ep_targeting_style = NULL,
      selected_ep_r_targeting_basis = NULL
    ))
  }

  raw_labels <- vapply(descriptors, `[[`, character(1), "label")
  raw_methods <- vapply(
    descriptors,
    function(descriptor) {
      if (is.null(descriptor$method)) {
        return(NA_character_)
      }
      descriptor$method
    },
    character(1)
  )
  raw_basis <- vapply(
    descriptors,
    function(descriptor) {
      if (is.null(descriptor$sieve_num_basis)) {
        return(NA_integer_)
      }
      as.integer(descriptor$sieve_num_basis)
    },
    integer(1)
  )
  raw_targeting_styles <- vapply(
    descriptors,
    function(descriptor) {
      if (is.null(descriptor$ep_targeting_style)) {
        return(NA_character_)
      }
      descriptor$ep_targeting_style
    },
    character(1)
  )
  raw_r_basis <- vapply(
    descriptors,
    function(descriptor) {
      if (is.null(descriptor$ep_r_targeting_basis)) {
        return(NA_character_)
      }
      descriptor$ep_r_targeting_basis
    },
    character(1)
  )

  candidate_labels <- unique(raw_labels)
  aggregated_coefficients <- numeric(length(candidate_labels))
  aggregated_methods <- rep(NA_character_, length(candidate_labels))
  aggregated_basis <- rep(NA_integer_, length(candidate_labels))
  aggregated_targeting_styles <- rep(NA_character_, length(candidate_labels))
  aggregated_r_basis <- rep(NA_character_, length(candidate_labels))

  for (index in seq_along(candidate_labels)) {
    label <- candidate_labels[[index]]
    matches <- raw_labels == label
    aggregated_coefficients[[index]] <- sum(coefficients[matches], na.rm = TRUE)

    method_values <- unique(raw_methods[matches])
    method_values <- method_values[!is.na(method_values)]
    if (length(method_values) > 0L) {
      aggregated_methods[[index]] <- method_values[[1L]]
    }

    basis_values <- unique(raw_basis[matches])
    basis_values <- basis_values[!is.na(basis_values)]
    if (length(basis_values) > 0L) {
      aggregated_basis[[index]] <- basis_values[[1L]]
    }

    targeting_style_values <- unique(raw_targeting_styles[matches])
    targeting_style_values <- targeting_style_values[!is.na(targeting_style_values)]
    if (length(targeting_style_values) > 0L) {
      aggregated_targeting_styles[[index]] <- targeting_style_values[[1L]]
    }

    r_basis_values <- unique(raw_r_basis[matches])
    r_basis_values <- r_basis_values[!is.na(r_basis_values)]
    if (length(r_basis_values) > 0L) {
      aggregated_r_basis[[index]] <- r_basis_values[[1L]]
    }
  }

  names(aggregated_coefficients) <- candidate_labels
  if (length(aggregated_coefficients) > 0L) {
    selected_index <- which.max(aggregated_coefficients)[[1L]]
    selected_candidate <- candidate_labels[[selected_index]]
    selected_method <- aggregated_methods[[selected_index]]
    if (is.na(selected_method)) {
      selected_method <- selected_candidate
    }
    selected_sieve_num_basis <- aggregated_basis[[selected_index]]
    if (is.na(selected_sieve_num_basis) || !identical(selected_method, "ep")) {
      selected_sieve_num_basis <- NULL
    }
    selected_ep_targeting_style <- aggregated_targeting_styles[[selected_index]]
    if (is.na(selected_ep_targeting_style) || !identical(selected_method, "ep")) {
      selected_ep_targeting_style <- NULL
    }
    selected_ep_r_targeting_basis <- aggregated_r_basis[[selected_index]]
    if (is.na(selected_ep_r_targeting_basis) || !identical(selected_ep_targeting_style, "r")) {
      selected_ep_r_targeting_basis <- NULL
    }
  } else {
    selected_candidate <- NULL
    selected_method <- NULL
    selected_sieve_num_basis <- NULL
    selected_ep_targeting_style <- NULL
    selected_ep_r_targeting_basis <- NULL
  }

  ep_basis_grid <- unique(aggregated_basis[aggregated_methods == "ep" & !is.na(aggregated_basis)])
  if (length(ep_basis_grid) == 0L) {
    ep_basis_grid <- NULL
  }

  list(
    candidate_labels = candidate_labels,
    coefficients = aggregated_coefficients,
    cv_risk = cv_risk,
    selected_candidate = selected_candidate,
    selected_method = selected_method,
    selected_sieve_num_basis = selected_sieve_num_basis,
    ep_basis_grid = ep_basis_grid,
    selected_ep_targeting_style = selected_ep_targeting_style,
    selected_ep_r_targeting_basis = selected_ep_r_targeting_basis
  )
}




#' Cross-Validate CATE Models
#'
#' Cross-validates a collection of CATE learners using the CATE selector loss.
#'
#' @param treatment_level The treatment level to be considered the treated
#'   group in the contrast used for selection.
#' @param control_level The treatment level to be considered the control
#'   group in the contrast used for selection.
#' @inheritParams cross_validate
#' @return A list containing the following elements:
#'   - `lrnr_sl`: A `Lrnr_sl` object that represents the cross-validated ensemble of CATE learners.
#'   - `cv_risk`: The cross-validation risk associated with the ensemble, which serves as a measure of the ensemble's predictive performance.
#'   - `coefficients`: The coefficients derived from the cross-validation process, providing insights into the relative importance of different learners within the ensemble.
#'   - `selection_summary`: A compact summary of candidate labels, selected candidate, selected method, and any EP basis-size metadata.
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
#' cv_fit <- cross_validate_cate(
#'   list(
#'     Lrnr_cate_DR$new(base_learner = Lrnr_mean$new()),
#'     Lrnr_cate_EP$new(base_learner = Lrnr_mean$new(), sieve_num_basis = 4)
#'   ),
#'   task,
#'   cv_control = list(V = 2)
#' )
#'
#' cv_fit$selection_summary
#' }
#' @export
cross_validate_cate <- function(hte_learners, hte3_task, cv_control = NULL, treatment_level = NULL, control_level = NULL, ...) {
  cv_metalearner <- Lrnr_cate_DR_selector$new(treatment_level = treatment_level, control_level = control_level)
  lrnr_sl <- cross_validate(hte_learners, hte3_task, cv_metalearner = cv_metalearner, cv_control = cv_control)
  selection_summary <- build_cv_selection_summary(lrnr_sl)
  attr(lrnr_sl, "hte3_selection_summary") <- selection_summary
  out <- list(
    lrnr_sl = lrnr_sl,
    cv_risk = lrnr_sl$fit_object$cv_meta_fit$cv_risk,
    coefficients = selection_summary$coefficients,
    selection_summary = selection_summary
  )
  return(out)
}



#' Cross-Validate CRR Models
#'
#' Cross-validates a collection of CRR learners using the CRR selector loss.
#'
#' @param treatment_level The treatment level to be considered the treated
#'   group in the contrast used for selection.
#' @param control_level The treatment level to be considered the control
#'   group in the contrast used for selection.
#' @inheritParams cross_validate
#' @return A list containing the following elements:
#'   - `lrnr_sl`: A `Lrnr_sl` object that represents the cross-validated ensemble of CRR learners.
#'   - `cv_risk`: The cross-validation risk associated with the ensemble, which serves as a measure of the ensemble's predictive performance.
#'   - `coefficients`: The coefficients derived from the cross-validation process, providing insights into the relative importance of different learners within the ensemble.
#'   - `selection_summary`: A compact summary of candidate labels, selected candidate, selected method, and any EP basis-size metadata.
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
#'   outcome = "Y_binary",
#'   propensity_learner = Lrnr_mean$new(),
#'   outcome_learner = Lrnr_mean$new(),
#'   mean_learner = Lrnr_mean$new(),
#'   cross_fit = FALSE
#' )
#'
#' cv_fit <- cross_validate_crr(
#'   list(
#'     Lrnr_crr_IPW$new(base_learner = Lrnr_mean$new()),
#'     Lrnr_crr_EP$new(base_learner = Lrnr_mean$new(), sieve_num_basis = 4)
#'   ),
#'   task,
#'   cv_control = list(V = 2)
#' )
#'
#' cv_fit$selection_summary
#' }
#' @export
cross_validate_crr <- function(hte_learners, hte3_task, cv_control = NULL, treatment_level = NULL, control_level = NULL, ...) {
  cv_metalearner <- Lrnr_crr_DR_selector$new(treatment_level = treatment_level, control_level = control_level)
  lrnr_sl <- cross_validate(hte_learners, hte3_task, cv_metalearner = cv_metalearner, cv_control = cv_control)
  selection_summary <- build_cv_selection_summary(lrnr_sl)
  attr(lrnr_sl, "hte3_selection_summary") <- selection_summary
  out <- list(
    lrnr_sl = lrnr_sl,
    cv_risk = lrnr_sl$fit_object$cv_meta_fit$cv_risk,
    coefficients = selection_summary$coefficients,
    selection_summary = selection_summary
  )
  return(out)
}






# cross_validate <- function(hte_learners, hte3_task, cv_metalearner = Lrnr_cv_selector$new(loss_squared_error), cv_loss_spec = loss_spec_cate_DR_binary, cv_control = NULL, ...) {
#   args <- list(hte3_task = hte3_task, ...)
#   if(is.list(hte_learners)) {
#     hte_learners <- Stack$new(hte_learners)
#   }
#   pseudo_data <- call_with_args(cv_loss_spec, args, silent = TRUE)
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
