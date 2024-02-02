

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
  new_data <- as.data.table(new_data)
  if(!hte3_learner$is_trained) {
    stop("hte3_learner must be a trained hte3 learner object.")
  }
  task <- hte3_learner$training_task
  id_var <- task$nodes$id
  if(!is.null(id_var) && !(id_var %in% names(new_data))) {
    new_data[[id_var]] <- 1:nrow(new_data)
  }


  new_task <- hte3_Task$new(new_data, npsem = task$npsem, likelihood = task$likelihood, id = id_var)
  return(hte3_learner$predict(new_task))
}
