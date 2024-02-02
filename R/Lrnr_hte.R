#' Lrnr_hte Class
#'
#' This class serves as a general template for constructing meta-learners to estimate heterogeneous treatment effects (HTEs).
#'
#' @format An R6 class with public methods to initialize the learner, create a regression task, and access the base learner.
#'
#' @param params A list of parameters for the meta-learning algorithm.
#' @param modifers A character vector of variable names for treatment effect moderators.
#' If \code{NULL}, then the \code{modifiers} variable of the \code{hte3_Task} object is used.
#' See \code{make_hte3_Task_tx} for more details.
#' @param base_learner A \code{\link{sl3}} learner object inheriting from \code{\link[sl3]{Lrnr_base}} that specifies the base supervised learning algorithm used by the meta-learner.
#' @param transform_function A function to transform the predictions of the base learner. Default is the identity transform.
#' @param pseudo_outcome_type The outcome type of the pseudo-outcome used by the meta-learner. Options are \code{c("continuous", "binomial", "quasibinomial")}. Default is \code{"continuous"}.
#'   For example, the DR-learner, EP-learner, and R-learner of the CATE involve (weighted) least-squares regression using a pseudo-outcome with \code{pseudo_outcome_type = "continuous"}.
#'   The CRATE EP-learner involves performing weighted logistic regression using a pseudo-outcome taking values in [0,1] with \code{pseudo_outcome_type = "quasibinomial"}.
#' @param pseudo_family A \code{\link[stats]{family}} object specifying the loss function (involving pseudo-weights and pseudo-outcomes) used to fit \code{base_learner} in the meta-learner algorithm. Default is \code{gaussian()}.
#'
#' @inheritParams Lrnr_base
#' @importFrom sl3 Lrnr_base
#' @export
Lrnr_hte <- R6Class(
  classname = "Lrnr_hte", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(

    initialize = function(params, base_learner,   transform_function = NULL, pseudo_outcome_type = c("continuous", "binomial", "quasibinomial"), pseudo_family = gaussian(),  ...) {
      if(!is.null(pseudo_outcome_type)) {
        pseudo_outcome_type <- match.arg(pseudo_outcome_type)
      }
      private$.base_learner <- base_learner

      super$initialize(params = params,  ...)
      private$.transform_function <- transform_function
      private$.pseudo_outcome_type <- pseudo_outcome_type
      private$.pseudo_family <- pseudo_family
    },
    #' (Used internally) Provides Pseudo-Data for Meta-Learner Training
    #'
    #' This method provides pseudo-data, such as pseudo-outcomes and pseudo-weights, which are used to fit \code{base_learner} with the loss function implied by \code{pseudo_family} in the meta-learner algorithm.
    #' NOTE: Observation weights provided by \code{hte3_task$weights} should not be included in the pseudo-weights of this function. The weights from \code{hte3_task$weights} are incorporated in the meta-learning fitting process downstream.
    #' @param hte3_task A \code{hte3_Task} object containing the data and necessary information for heterogeneous treatment effect estimation.
    #' @param ... Additional arguments in \code{params} needed to compute the pseudo-data.
    #'
    #' @return A list containing attributes \code{pseudo_outcome} and \code{pseudo_weights}.
    get_pseudo_data = function(hte3_task, ...) {
      msg <- stop("Lrnr_hte is a template class and does not implement get_pseudo_data.")
      return(list(pseudo_outcome = msg, pseudo_weights = msg))
    },
    #' (Used internally) Check Compatibility with Treatment Type
    #'
    #' This method checks whether this meta-learner is compatible with the \code{treatment} variable type in \code{hte3_Task}.
    #'
    #' @param hte3_task A \code{hte3_Task} object containing the data and necessary information for heterogeneous treatment effect estimation.
    #'
    #' @keywords internal
    check_treatment_type = function(hte3_task) {
      type <- hte3_task$npsem$treatment$variable_type$type
      valid_types <- paste0("c(", paste0(private$.treatment_type, collapse = ", "), ")")
      if(!(type %in% private$.treatment_type)) {
        stop(paste0(self$name, "is not compatible with ", type, " treatments. Only treatments of the following type are permitted for this learner: ", valid_types))
      }
    },
    get_modifiers = function(hte3_task, return_matrix = FALSE) {
      if("modifiers" %in% names(self$params)) {
        modifiers <- self$params$modifiers
      } else if ("covariates" %in% names(self$params)) {
        modifiers <- self$params$covariates
      } else {
        modifiers <- hte3_task$nodes$covariates
      }
      if(is.null(modifiers) || length(modifiers) == 0) {
        modifiers <- hte3_task$npsem$modifiers$variables
      }
      if(return_matrix) {
        X <- as.matrix(hte3_task$data[, modifiers, with = FALSE])
        return(X)
      } else {
        return(modifiers)
      }

    },
    #' (Used internally) Construct Meta-Learner Task
    #'
    #' This method constructs a \code{sl3_Task} object used to train \code{base_learner} in the meta-learner algorithm.
    #'
    #' @param hte3_task A \code{hte3_Task} object containing the data and necessary information for heterogeneous treatment effect estimation.
    #' @param train Logical indicating whether to create the task for training or prediction. Default is \code{TRUE}.
    #' If \code{FALSE} then a \code{hte3_Task} object for prediction is returned with covariates being the effect \code{modifiers}.
    #' @importFrom sl3 make_sl3_Task
    #'
    #' @return A \code{sl3_Task} object containing the pseudo-data of \code{get_pseudo_data} and \code{outcome_Type}=\code{pseudo_outcome_type}.
    #' @keywords internal
    make_metalearner_task = function(hte3_task, train = TRUE) {
      params <- self$params
      get_pseudo_data <- self$get_pseudo_data
      params$hte3_task <- hte3_task
      modifiers <- self$get_modifiers(hte3_task)

      if(!train) {
        # if for prediction, just send task.
        return(hte3_task$next_in_chain(covariates = modifiers))
      }
      # use specified meta learner to get pseudo outcome and pseudo weight data
      pseudo_data <- call_with_args(get_pseudo_data, params, silent = TRUE)
      pseudo_outcome <- pseudo_data$pseudo_outcome
      pseudo_weights <- pseudo_data$pseudo_weights
      outcome_type <- self$pseudo_outcome_info$outcome_type
      if(is.null(outcome_type)) outcome_type <- "continuous"
      # incorporate observational weights into pseudo_weights
      if(is.null(pseudo_weights)) {
        pseudo_weights <- hte3_task$weights
      } else {
        pseudo_weights <- hte3_task$weights * pseudo_weights
      }

      new_data <- data.table(pseudo_outcome, pseudo_weights)
      names(new_data) <- c("pseudo_outcome", "pseudo_weights")
      column_names <- hte3_task$add_columns(new_data)
      training_task <- hte3_task$next_in_chain(covariates = modifiers,
                                               outcome = "pseudo_outcome",
                                               weights = "pseudo_weights",
                                               column_names = column_names,
                                               new_outcome_type = outcome_type)

      return(training_task)
    }
  ),
  active = list(
    #' Access the Base Learner
    #'
    #' This method returns the (untrained) base learner specified during initialization.
    #'
    #' @return The base learner object used in the meta-learner.
    base_learner = function() {
      return(private$.base_learner)
    },
    #' Retrieve Pseudo-Outcome Information
    #'
    #' This method returns information about the pseudo-outcome used in the meta-learner.
    #'
    #' @return A list containing attributes \code{outcome_type} and \code{family}.
    pseudo_outcome_info = function() {
      return(list(outcome_type = private$.pseudo_outcome_type, family =private$.pseudo_family ))
    },
    #' Access the Transform Function
    #'
    #' This method returns the transform function used to transform predictions of \code{base_learner}.
    #'
    #' @return The transform function for predictions.
    transform_function = function() {
      fun <- private$.transform_function
      if(is.function(fun)) fun <- Vectorize(fun)
      return(fun)
    }
  ),
  private = list(
    .properties = c(
      "continuous", "binomial", "categorical", "importance",
      "weights"),
    .train = function(hte3_task) {
      args <- self$params
      learner_task <- self$make_metalearner_task(hte3_task)
      base_learner <- self$base_learner
      # make sure right family is being used.
      family <- self$pseudo_outcome_info$family
      base_learner <- base_learner$reparameterize(list(family = family))
      # some learners dont handle zero weights well for training
      learner_trained <- base_learner$train(learner_task)
      private$.base_learner <- learner_trained
      fit_object = list(base_learner_trained = learner_trained)
      return(fit_object)
    },
    .predict = function(hte3_task) {
      fit_object <- self$fit_object
      learner_trained <- fit_object$base_learner_trained
      learner_task <- self$make_metalearner_task(hte3_task, train = FALSE)
      predictions <- learner_trained$predict(learner_task)
      # transform function

      transform <- self$transform_function
      if(!is.null(transform)) {
        if(is.data.table(predictions)) {
          predictions <- as.data.table(apply(predictions, 2, transform))
        } else {
          predictions <- transform(predictions)
        }
      }
      # compute predictions

      return(predictions)
    },
    .chain = function(hte3_task) {
      # TODO make CV use chain
      predictions <- self$predict(hte3_task)
      predictions <- as.data.table(predictions)
      # Add predictions as new columns
      learner_task <- self$make_metalearner_task(hte3_task, train = TRUE)#$revere_fold_task("full")

      new_col_names <- learner_task$add_columns(predictions, self$fit_uuid)
      return(learner_task$next_in_chain(
        covariates = names(predictions),
        column_names = new_col_names
      ))
    },
    .treatment_type = c("binomial_treatment", "categorical_treatment", "continuous_treatment"),
    .required_packages = c("sl3"),
    .base_learner = NULL,
    .pseudo_outcome_type = "continuous",
    .pseudo_family = gaussian(),
    .transform_function = NULL
  )
)
