#' Task object for meta-learners in causal data structures.
#'
#' Constructs a \code{hte3_Task} object for meta-learners in causal data structures, containing relevant data and nuisance function estimators.
#'
#' Class for Storing Data, NPSEM, and Nuisance Function Estimators for hte3 Learners
#'
#' This class inherits from \code{\link[sl3]{sl3_Task}} and \code{\link[tmle3]{tmle3_Task}}.
#' In addition to all the methods supported by \code{\link[sl3]{sl3_Task}} and \code{\link[tmle3]{tmle3_Task}},
#' it supports the following functionalities specific to the heterogeneous treatment effect estimation (hte3) framework.
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#' @importFrom sl3 sl3_Task
#' @importFrom digest digest
#' @import data.table tmle3
#'
#' @keywords data
#'
#' @return An object of class \code{hte3_Task}.
#'
#' @format \code{\link{R6Class}} object.
#'
#' @param data A named data frame or data.table containing treatment effect modifiers,
#'   potential confounders, treatment, outcome, and optionally weights and subject IDs.
#'   See the \code{\link{sl3_Task}} documentation for further options.
#'
#' @param npsem A list containing \code{\link[tmle3]{tmle3_Node}} objects containing the non-parametric structural equation model (NPSEM)
#'   specifying the relationship between variables in the treatment effect estimation framework.
#'
#' @param likelihood An \code{\link[tmle3]{Likelihood}} object specifying the relevant likelihood/nuisance estimators for the meta-learner.
#'   used for estimating the parameters of the NPSEM.
#'
#' @param ... Additional arguments to pass to the initialization function.
#'
#' @export
hte3_Task <- R6Class(
  classname = "hte3_Task",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Task,
  public = list(
    initialize = function(data, npsem, likelihood = NULL,   ...) {
      super$initialize(data, npsem, ...)
      dot_args <- list(...)


      if(is.null(likelihood)) {
        likelihood <- Likelihood$new(factor_list = list())
      } else {
        likelihood$train(self)
      }

      private$.likelihood <- likelihood
    },
    add_nuisance_estimator = function(node, learner) {
      factor <- LF_fit$new(node, learner, type = "mean")
      self$likelihood$add_factors(list(factor))
    },
    get_nuisance_estimates = function(nodes, hte3_task = NULL, fold_number = "validation") {
      if(is.null(hte3_task)) hte3_task <- self
      estimates <- self$likelihood$get_likelihoods(hte3_task, nodes, fold_number = fold_number)

      if(class(estimates[[1]]) %in% c("packed_predictions")) {
        estimates <- sl3:::unpack_predictions(as.matrix(estimates))
      }
      return(estimates)
    },
    next_in_chain = function(covariates = NULL, outcome = NULL, id = NULL,
                             weights = NULL, offset = NULL, time = NULL,
                             folds = NULL, column_names = NULL,
                             new_nodes = NULL, new_outcome_type = NULL, ...) {
      if (is.null(new_nodes)) {
        new_nodes <- self$nodes

        if (!is.null(covariates)) {
          new_nodes$covariates <- covariates
        }

        if (!is.null(outcome)) {
          new_nodes$outcome <- outcome
        }

        if (!missing(id)) {
          new_nodes$id <- id
        }

        if (!missing(weights)) {
          new_nodes$weights <- weights
        }

        if (!missing(offset)) {
          new_nodes$offset <- offset
        }

        if (!missing(time)) {
          new_nodes$time <- time
        }
      }

      if (is.null(column_names)) {
        column_names <- private$.column_names
      }

      if (is.null(folds)) {
        folds <- private$.folds
      }

      all_nodes <- unlist(new_nodes)

      # verify nodes are contained in dataset
      missing_cols <- setdiff(all_nodes, names(column_names))

      assertthat::assert_that(
        length(missing_cols) == 0,
        msg = sprintf(
          "Couldn't find %s",
          paste(missing_cols, collapse = " ")
        )
      )
      new_task <- self$clone()

      if(is.null(new_outcome_type)) {
        if ((is.null(new_nodes$outcome) &&
             is.null(self$nodes$outcome)) ||
            all(new_nodes$outcome == self$nodes$outcome)) {
          # if we have the same outcome, transfer outcome properties
          new_outcome_type <- self$outcome_type
        } else {
          # otherwise, let the new task guess
          new_outcome_type <- NULL
        }
      }
      new_task$initialize(
        data = private$.shared_data,
        npsem = self$npsem,
        likelihood = self$likelihood,
        nodes = new_nodes,
        folds = folds,
        column_names = column_names,
        row_index = private$.row_index,
        outcome_type = new_outcome_type,
        ...
      )
      return(new_task)
    },
    subset_task = function(row_index, drop_folds = FALSE) {
      if (is.logical(row_index)) {
        row_index <- which(row_index)
      }
      old_row_index <- private$.row_index
      if (!is.null(old_row_index)) {
        # index into the logical rows of this task
        row_index <- old_row_index[row_index]
      }

      must_reindex <- any(duplicated(row_index))
      if (must_reindex) {
        new_shared_data <- private$.shared_data$clone()
        new_shared_data$reindex(row_index)
        row_index <- seq_along(row_index)
      } else {
        new_shared_data <- private$.shared_data
      }

      new_task <- self$clone()
      if (drop_folds) {
        new_folds <- NULL
      } else {
        if (must_reindex) {
          stop("subset indices have copies, this requires dropping folds.")
        }
        new_folds <- subset_folds(private$.folds, row_index)
      }

      new_task$initialize(
        data = new_shared_data,
        npsem = self$npsem,
        likelihood = self$likelihood,
        nodes = private$.nodes,
        folds = new_folds,
        column_names = private$.column_names,
        row_index = row_index,
        outcome_type = self$outcome_type
      )
      return(new_task)
    }

  ),
  active = list(
    likelihood = function() {
      return(private$.likelihood)
    },
    data = function() {
      all_variables <- unique(c(unlist(lapply(self$npsem, `[[`, "variables")), unlist(self$nodes )))
      self$get_data(columns = all_variables)
    }
  ),
  private = list(
    .likelihood = NULL
  )
)
