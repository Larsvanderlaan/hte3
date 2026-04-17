hte3_repository <- function() {
  "Larsvanderlaan/hte3"
}

hte3_install_command <- function() {
  sprintf('remotes::install_github("%s", dependencies = TRUE)', hte3_repository())
}

supported_cate_methods <- function() {
  c("dr", "r", "t", "ep")
}

supported_crr_methods <- function() {
  c("ep", "ipw", "t")
}

supported_grf_cate_methods <- function() {
  c("dr", "r", "ep")
}

supported_grf_crr_methods <- function() {
  c("ep", "ipw", "t")
}

default_cate_method <- function() {
  "dr"
}

default_crr_method <- function() {
  "ep"
}

automl_spec_table <- function() {
  list(
    list(
      id = "glmnet",
      learner_class = "Lrnr_glmnet",
      required_package = NULL,
      constructor_args = list()
    ),
    list(
      id = "gam",
      learner_class = "Lrnr_gam",
      required_package = NULL,
      constructor_args = list()
    ),
    list(
      id = "earth",
      learner_class = "Lrnr_earth",
      required_package = "earth",
      constructor_args = list(degree = 2)
    ),
    list(
      id = "ranger",
      learner_class = "Lrnr_ranger",
      required_package = "ranger",
      constructor_args = list(max.depth = 10)
    ),
    list(
      id = "xgboost_depth2",
      learner_class = "Lrnr_xgboost_early_stopping",
      required_package = "xgboost",
      constructor_args = list(
        min_child_weight = 15,
        max_depth = 2,
        eta = 0.20,
        subsample = 0.8,
        colsample_bytree = 0.8
      )
    ),
    list(
      id = "xgboost_depth3",
      learner_class = "Lrnr_xgboost_early_stopping",
      required_package = "xgboost",
      constructor_args = list(
        min_child_weight = 15,
        max_depth = 3,
        eta = 0.15,
        subsample = 0.9
      )
    ),
    list(
      id = "xgboost_depth4",
      learner_class = "Lrnr_xgboost_early_stopping",
      required_package = "xgboost",
      constructor_args = list(
        min_child_weight = 15,
        max_depth = 4,
        eta = 0.15,
        subsample = 0.9
      )
    ),
    list(
      id = "xgboost_depth5",
      learner_class = "Lrnr_xgboost_early_stopping",
      required_package = "xgboost",
      constructor_args = list(
        min_child_weight = 15,
        max_depth = 5,
        eta = 0.15,
        subsample = 0.9
      )
    ),
    list(
      id = "xgboost_depth4_slow",
      learner_class = "Lrnr_xgboost_early_stopping",
      required_package = "xgboost",
      constructor_args = list(
        min_child_weight = 15,
        max_depth = 4,
        eta = 0.08,
        subsample = 0.8,
        colsample_bytree = 0.8
      )
    )
  )
}
