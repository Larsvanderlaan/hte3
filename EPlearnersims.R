


# generate synthetic dataset
nsims <- 1000
key <- paste0("simsEP_n=",n, "_pos=",pos, "_hard=", hard, "_base_learner=" , base_learner, "_sim_type=", sim_type, "_nsims=", nsims)

n <- as.numeric(n)
hard <- as.logical(hard)
pos <- as.logical(pos)
sim_type <- as.character(sim_type)
base_learner <- as.character(base_learner)
library(data.table)
library(sl3)
library(hte3)
set.seed(12345)
# Generate data
if(sim_type == "CATEhigh") {
  source("~/EPlearner/FinalSimulationCode/simCATEHighDim.R")
  sim.fun <- sim.CATEHighDim
} else if(sim_type == "CATElow") {
  source("~/EPlearner/FinalSimulationCode/simCATE.R")
  sim.fun <- sim.CATE
} else if (sim_type == "CRATE") {
  sim.fun <- sim.LRR
}

make_ep_stack_sims <- function(base_learner, hte3_task, treatment_level = 1, control_level = 0, sieve_interaction_order = interaction_order) {
  d <- length(hte3_task$npsem$modifiers$variables)
  n <- hte3_task$nrow
  lrnr_cate_epstack <- Stack$new(Lrnr_cate_EP$new(base_learner = base_learner, sieve_num_basis  = d,
                                                  sieve_interaction_order = sieve_interaction_order,
                                                  treatment_level = treatment_level, control_level = control_level)
                                 , Lrnr_cate_EP$new(base_learner = base_learner, sieve_num_basis  = 2*d,
                                                    sieve_interaction_order = sieve_interaction_order,
                                                    treatment_level = treatment_level, control_level = control_level),
                                 Lrnr_cate_EP$new(base_learner = base_learner, sieve_num_basis = 3*d,
                                                  sieve_interaction_order = sieve_interaction_order,
                                                  treatment_level = treatment_level, control_level = control_level),
                                 Lrnr_cate_EP$new(base_learner = base_learner, sieve_num_basis = 5*d,
                                                  sieve_interaction_order = sieve_interaction_order,
                                                  treatment_level = treatment_level, control_level = control_level),
                                 Lrnr_cate_EP$new(base_learner = base_learner, sieve_num_basis = 8*d,
                                                  sieve_interaction_order = sieve_interaction_order,
                                                  treatment_level = treatment_level, control_level = control_level))
  return(lrnr_cate_epstack)
}


one_sim <- function(iter) {
  try({
    print(iter)
    print(key)
    data <- as.data.frame(sim.fun(n, hard, pos))
    data_val <- as.data.frame(sim.fun(n = 10000, hard, pos))
    # confounders for estimation of CATE
    confounders <- grep("[W].+", colnames(data), value = TRUE)
    # nuisance learning ensemble for outcome regression and propensity score.
    lrnr_nuisance <- Stack$new(Lrnr_xgboost$new(max_depth = 2, min_child_weight = 15),
                               Lrnr_xgboost$new(max_depth = 3, min_child_weight = 15),
                               Lrnr_xgboost$new(max_depth = 4, min_child_weight = 15),
                               Lrnr_xgboost$new(max_depth = 5, min_child_weight = 15),
                               Lrnr_earth$new(degree = 2))

    # make hte3 task (with cross-fitted nuisances)
    hte3_task <- make_hte3_Task_tx(data, modifiers = confounders, confounders = confounders, treatment = "A", outcome = "Y",
                                   learner_pi = lrnr_nuisance,
                                   learner_mu = lrnr_nuisance,
                                   learner_m = lrnr_nuisance)



    # base learning algorithms for HTE meta-learners
    lrnr_xg_stack <- Stack$new(Lrnr_xgboost$new(max_depth = 1, min_child_weight = 5),
                               Lrnr_xgboost$new(max_depth = 2, min_child_weight = 5),
                               Lrnr_xgboost$new(max_depth = 3, min_child_weight = 5),
                               Lrnr_xgboost$new(max_depth = 4, min_child_weight = 5),
                               Lrnr_xgboost$new(max_depth = 5, min_child_weight = 5))
    lrnr_rf_stack <- Stack$new(Lrnr_ranger$new(max.depth = 4),
                               Lrnr_ranger$new(max.depth = 6),
                               Lrnr_ranger$new(max.depth = 8))
    lrnr_gam <- Lrnr_gam$new()
    lrnr_earth <-  Lrnr_earth$new(degree = 1)

    if(sim_type == "CATEhigh") {
      interaction_order <- 1
    } else {
      interaction_order <- 3
    }
    if(base_learner %in% c("xg", "rf")) {
      # xg is a s tack, so cross-validate everything
      base_lrnr <- get(paste0("lrnr_", base_learner, "_stack"))
      ep_lrnr <- cross_validate(make_ep_stack_sims(base_lrnr, hte3_task, sieve_interaction_order = interaction_order), hte3_task, cv_control = list(V = 5))
      ep_lrnr_lasso <- cross_validate(Lrnr_cate_EP$new(base_learner = base_lrnr, sieve_num_basis = 300,
                                                       sieve_interaction_order = interaction_order,
                                                       treatment_level = 1, control_level = 0, screen_basis_with_lasso = TRUE), hte3_task, cv_control = list(V = 5))
      dr_lrnr <- cross_validate(Lrnr_cate_DR$new(base_lrnr), hte3_task, cv_control = list(V = 5))
      r_lrnr <- cross_validate(Lrnr_cate_R$new(base_lrnr), hte3_task, cv_control = list(V = 5))
      t_lrnr <- cross_validate(Lrnr_cate_T$new(base_lrnr), hte3_task, cv_control = list(V = 5))
    } else {
      # Only need to cross-validate the EP-learner sieve
      base_lrnr <- get(paste0("lrnr_", base_learner))
      ep_lrnr <- cross_validate(make_ep_stack_sims(base_lrnr, hte3_task, sieve_interaction_order = interaction_order), hte3_task, cv_control = list(V = 5))
      ep_lrnr_lasso <- Lrnr_cate_EP$new(base_learner = base_lrnr, sieve_num_basis = 300,
                                        sieve_interaction_order = interaction_order,
                                        treatment_level = 1, control_level = 0, screen_basis_with_lasso = TRUE)$train(hte3_task)
      dr_lrnr <- Lrnr_cate_DR$new(base_lrnr)$train(hte3_task)
      if(base_learner == "earth") {
        r_lrnr <- Lrnr_cate_R$new(Lrnr_mean$new())$train(hte3_task)
      } else {
        r_lrnr <- Lrnr_cate_R$new(base_lrnr)$train(hte3_task)
      }

      t_lrnr <- Lrnr_cate_T$new(base_lrnr)$train(hte3_task)
    }
    true_cate <- data_val$EY1W - data_val$EY0W
    ep_preds <- predict_hte3(ep_lrnr, new_data = data_val)
    ep_lasso_preds <- predict_hte3(ep_lrnr_lasso, new_data = data_val)
    dr_preds <- predict_hte3(dr_lrnr, new_data = data_val)
    r_preds <- predict_hte3(r_lrnr, new_data = data_val)
    t_preds <- predict_hte3(t_lrnr, new_data = data_val)

    # ep_lrnr$cv_risk(loss_squared_error)
    ep_lrnr_stack <- ep_lrnr$fit_object$cv_fit
    ep_preds_stack <- predict_hte3(ep_lrnr_stack, new_data = data_val)
    print(colMeans((as.matrix(ep_preds_stack) - true_cate)^2))
    # mean((ep_preds - true_cate)^2)
    # mean((ep_lasso_preds - true_cate)^2)
    # mean((true_cate - dr_preds)^2)

    ep_risk <- mean((true_cate - ep_preds)^2)
    ep_lasso_risk <- mean((true_cate - ep_lasso_preds)^2)
    dr_risk <- mean((true_cate - dr_preds)^2)
    r_risk <- mean((true_cate - r_preds)^2)
    t_risk <- mean((true_cate - t_preds)^2)

    if(base_learner %in% c("rf", "xg")) {
      X <- hte3_task$get_tmle_node("modifiers")
      W <- hte3_task$get_tmle_node("treatment")
      Y <- hte3_task$get_tmle_node("outcome")
      cf_fit <- grf::causal_forest(X = as.matrix(X),
                                   Y = as.vector(as.matrix(Y)),
                                   W = as.vector(as.matrix(W)),
                                   , tune.parameters = "all")
      preds_cf <- as.vector(predict(cf_fit, newdata = as.matrix(data_val)[, confounders])$predictions)
      causal_forest_risk <- mean((true_cate - preds_cf)^2)
    } else{
      causal_forest_risk <- NULL
    }


    out <- data.table(iter = iter, ep_lasso_risk = ep_lasso_risk, ep_risk = ep_risk, dr_risk = dr_risk, r_risk = r_risk, t_risk = t_risk, causal_forest_risk = causal_forest_risk, lrnr = base_learner)
    return(out)
  })
  return(data.table())
}

results <- rbindlist(lapply(1:nsims, one_sim))
fwrite(results, file = paste0("~/EPlearner/results/", key))
