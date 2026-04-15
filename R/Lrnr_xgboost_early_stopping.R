#' XGBoost Learner with Early Stopping
#'
#' This learner fits gradient-boosted trees with an internal validation split
#' to select the number of boosting rounds by early stopping, then refits on
#' the full training sample at the selected iteration count.
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data
#' @return Learner object with methods for training and prediction. See
#'   \code{\link[sl3]{Lrnr_base}} for documentation on learners.
#' @format \code{\link{R6Class}} object.
#' @family Learners
#'
#' @section Parameters:
#' \describe{
#'   \item{\code{nrounds = 500}}{Maximum number of boosting rounds considered
#'   before early stopping.}
#'   \item{\code{validation_fraction = 0.2}}{Fraction of the training sample
#'   reserved for the internal validation split.}
#'   \item{\code{early_stopping_rounds = 25}}{Number of validation rounds
#'   without improvement before stopping.}
#'   \item{\code{seed = 1}}{Random seed used for the internal validation split.}
#'   \item{\code{...}}{Additional arguments passed to
#'   \code{\link[xgboost]{xgb.train}}.}
#' }
Lrnr_xgboost_early_stopping <- R6Class(
  classname = "Lrnr_xgboost_early_stopping", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(nrounds = 500,
                          validation_fraction = 0.2,
                          early_stopping_rounds = 25,
                          seed = 1,
                          nthread = 1,
                          ...) {
      params <- c(
        list(
          nrounds = nrounds,
          validation_fraction = validation_fraction,
          early_stopping_rounds = early_stopping_rounds,
          seed = seed,
          nthread = nthread
        ),
        list(...)
      )
      super$initialize(params = params, ...)
    },

    importance = function(...) {
      self$assert_trained()

      args <- c(list(model = self$fit_object), list(...))
      importance_result <- call_with_args(xgboost::xgb.importance, args)
      rownames(importance_result) <- importance_result[["Feature"]]
      importance_result
    }
  ),

  private = list(
    .properties = c(
      "continuous", "binomial", "categorical", "weights",
      "offset", "importance"
    ),

    .train = function(task) {
      verbose <- getOption("sl3.verbose")
      if (is.null(verbose)) {
        verbose <- 0L
      }
      args <- self$params
      outcome_type <- self$get_outcome_type(task)

      Xmat <- as.matrix(task$X)
      if (is.integer(Xmat)) {
        Xmat[, 1] <- as.numeric(Xmat[, 1])
      }

      Y <- outcome_type$format(task$Y)
      if (outcome_type$type == "categorical") {
        Y <- as.numeric(Y) - 1
      }

      weights <- if (task$has_node("weights")) task$weights else NULL

      training_offset <- FALSE
      base_margin <- NULL
      link_fun <- NULL
      if (task$has_node("offset")) {
        if (outcome_type$type == "categorical") {
          stop("offsets are not supported for outcome_type = 'categorical'")
        }

        family <- args$family
        if (is.null(family)) {
          family <- outcome_type$glm_family(return_object = TRUE)
        }
        link_fun <- family$linkfun
        base_margin <- task$offset_transformed(link_fun)
        training_offset <- TRUE
      }

      if (is.null(args$objective)) {
        if (outcome_type$type == "binomial") {
          args$objective <- "binary:logistic"
        } else if (outcome_type$type == "quasibinomial") {
          args$objective <- "reg:logistic"
        } else if (outcome_type$type == "categorical") {
          args$objective <- "multi:softprob"
          args$num_class <- length(outcome_type$levels)
        }
      }

      validation_fraction <- args$validation_fraction
      early_stopping_rounds <- args$early_stopping_rounds
      split_seed <- args$seed
      args$validation_fraction <- NULL
      args$seed <- NULL
      args$family <- NULL
      args$verbose <- as.integer(verbose)
      args$print_every_n <- 1000

      split <- private$.make_validation_split(
        y = task$Y,
        validation_fraction = validation_fraction,
        seed = split_seed,
        stratify = outcome_type$type %in% c("binomial", "categorical")
      )

      best_nrounds <- args$nrounds
      if (!is.null(split) && !is.null(early_stopping_rounds) && early_stopping_rounds > 0) {
        dtrain <- private$.make_dmatrix(
          x = Xmat[split$training, , drop = FALSE],
          y = Y[split$training],
          weights = if (is.null(weights)) NULL else weights[split$training],
          base_margin = if (is.null(base_margin)) NULL else base_margin[split$training]
        )
        dvalid <- private$.make_dmatrix(
          x = Xmat[split$validation, , drop = FALSE],
          y = Y[split$validation],
          weights = if (is.null(weights)) NULL else weights[split$validation],
          base_margin = if (is.null(base_margin)) NULL else base_margin[split$validation]
        )

        train_args <- args
        train_args$data <- dtrain
        train_args$watchlist <- list(train = dtrain, validation = dvalid)
        train_args$early_stopping_rounds <- early_stopping_rounds

        validation_fit <- call_with_args(xgboost::xgb.train, train_args)
        if (!is.null(validation_fit$best_iteration)) {
          best_nrounds <- validation_fit$best_iteration
        }
      }

      dfull <- private$.make_dmatrix(
        x = Xmat,
        y = Y,
        weights = weights,
        base_margin = base_margin
      )

      final_args <- args
      final_args$data <- dfull
      final_args$nrounds <- best_nrounds
      final_args$watchlist <- NULL
      final_args$early_stopping_rounds <- NULL

      fit_object <- call_with_args(xgboost::xgb.train, final_args)
      fit_object$best_iteration <- best_nrounds
      fit_object$best_ntreelimit <- best_nrounds
      fit_object$selected_nrounds <- best_nrounds
      fit_object$training_offset <- training_offset
      fit_object$link_fun <- link_fun

      fit_object
    },

    .predict = function(task = NULL) {
      outcome_type <- private$.training_outcome_type

      Xmat <- as.matrix(task$X)
      if (is.integer(Xmat)) {
        Xmat[, 1] <- as.numeric(Xmat[, 1])
      }

      xgb_data <- xgboost::xgb.DMatrix(Xmat)

      if (self$fit_object$training_offset) {
        offset <- task$offset_transformed(
          self$fit_object$link_fun,
          for_prediction = TRUE
        )
        xgboost::setinfo(xgb_data, "base_margin", offset)
      }

      fit_object <- private$.fit_object
      predictions <- rep.int(list(numeric()), 1)

      if (nrow(Xmat) > 0) {
        ntreelimit <- 0
        booster <- fit_object[["params"]][["booster"]]
        is_gblinear <- !is.null(booster) && identical(booster, "gblinear")
        if (!is.null(fit_object[["best_ntreelimit"]]) && !is_gblinear) {
          ntreelimit <- fit_object[["best_ntreelimit"]]
        }

        predictions <- stats::predict(
          fit_object,
          newdata = xgb_data,
          ntreelimit = ntreelimit,
          reshape = TRUE
        )
      }

      if (outcome_type$type == "categorical") {
        predictions <- sl3::pack_predictions(predictions)
      }

      predictions
    },

    .make_dmatrix = function(x, y, weights = NULL, base_margin = NULL) {
      dmatrix <- xgboost::xgb.DMatrix(x, label = y)

      if (!is.null(weights)) {
        xgboost::setinfo(dmatrix, "weight", weights)
      }

      if (!is.null(base_margin)) {
        xgboost::setinfo(dmatrix, "base_margin", base_margin)
      }

      dmatrix
    },

    .make_validation_split = function(y, validation_fraction, seed, stratify = FALSE) {
      n <- length(y)
      if (is.null(validation_fraction) ||
          validation_fraction <= 0 ||
          validation_fraction >= 1 ||
          n < 10) {
        return(NULL)
      }

      old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      if (old_seed_exists) {
        old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      }

      if (!is.null(seed)) {
        set.seed(seed)
        on.exit({
          if (old_seed_exists) {
            assign(".Random.seed", old_seed, envir = .GlobalEnv)
          } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv)
          }
        }, add = TRUE)
      }

      if (isTRUE(stratify)) {
        groups <- split(seq_len(n), as.character(y))
        validation_idx <- unlist(lapply(groups, function(group_idx) {
          if (length(group_idx) <= 1L) {
            return(integer())
          }

          n_group_validation <- floor(length(group_idx) * validation_fraction)
          if (n_group_validation == 0L) {
            n_group_validation <- 1L
          }

          sample(group_idx, size = min(n_group_validation, length(group_idx) - 1L))
        }), use.names = FALSE)
      } else {
        n_validation <- floor(n * validation_fraction)
        if (n_validation < 1L || n_validation >= n) {
          return(NULL)
        }
        validation_idx <- sample(seq_len(n), size = n_validation)
      }

      validation_idx <- sort(unique(validation_idx))
      if (length(validation_idx) < 1L || length(validation_idx) >= n) {
        return(NULL)
      }

      list(
        training = setdiff(seq_len(n), validation_idx),
        validation = validation_idx
      )
    },

    .required_packages = c("xgboost")
  )
)
