

#' Predict Heterogeneous Treatment Effects using hte3 Learners
#'
#' Generates treatment effect predictions using a trained hte3 learner.
#'
#' @param hte3_learner A trained hte3 learner (\code{Lrnr_sl} object) used to make predictions.
#' @param new_data A data frame or data.table containing new data of effect \code{modifiers} for prediction.
#'
#' @return A numeric vector of predicted heterogeneous treatment effects for the provided new data.
#'
#' @importFrom data.table as.data.table
#' @export
predict_hte3 <- function(hte3_learner, new_data) {
  if (inherits(new_data, "hte3_Task")) {
    return(hte3_learner$predict(new_data))
  }

  new_data <- as.data.table(new_data)
  if (!hte3_learner$is_trained) {
    stop("hte3_learner must be a trained hte3 learner object.")
  }
  task <- hte3_learner$training_task
  required_covariates <- task$npsem$modifiers$variables
  assert_columns_present(new_data, required_covariates, "modifier")

  id_var <- task$nodes$id
  weight_var <- task$nodes$weights
  npsem_columns <- unique(unlist(lapply(task$npsem, `[[`, "variables")))

  missing_npsem_cols <- setdiff(npsem_columns, names(new_data))
  for (column in missing_npsem_cols) {
    new_data[[column]] <- NA
  }

  if (!is.null(id_var) && !(id_var %in% names(new_data))) {
    new_data[[id_var]] <- seq_len(nrow(new_data))
  }

  if (!is.null(weight_var) && !(weight_var %in% names(new_data))) {
    new_data[[weight_var]] <- 1
  }

  new_task <- hte3_Task$new(
    new_data,
    npsem = task$npsem,
    likelihood = Likelihood$new(factor_list = list()),
    id = id_var,
    weights = weight_var
  )
  hte3_learner$predict(new_task)
}
