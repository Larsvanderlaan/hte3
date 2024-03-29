


# generate synthetic dataset
nsims <- 1000
key <- paste0("simsEP_n=",n, "_pos=",pos, "_hard=", hard, "_base_learner=" , base_learner, "_sim_type=", "crr", "_nsims=", nsims)

n <- as.numeric(n)
hard <- as.logical(hard)
pos <- as.logical(pos)
base_learner <- as.character(base_learner)
library(data.table)
library(sl3)
library(hte3)
library(SuperLearner)
set.seed(12345)
# Generate data
source("~/EPlearner/FinalSimulationCode/simRR.R")
sim.fun <- sim.RR



SL.gam_2 <- function(Y, X, newX, family, obsWeights, cts.num = 4,...) {
  deg.gam <- 2
  SuperLearner::SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num,... )
}
SL.gam_3 <- function(Y, X, newX, family, obsWeights, cts.num = 4,...) {
  deg.gam <- 3
  SuperLearner::SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num,... )
}
SL.gam_5 <- function(Y, X, newX, family, obsWeights, cts.num = 4,...) {
  deg.gam <- 5
  SuperLearner::SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num,... )
}

lrnr_gam_stack <- Stack$new(
  Lrnr_pkg_SuperLearner$new("SL.gam_2", name = paste0("Lrnr_gam_s", 2, "_x")),
  Lrnr_pkg_SuperLearner$new("SL.gam_3", name = paste0("Lrnr_gam_s", 3, "_x")),
  Lrnr_pkg_SuperLearner$new("SL.gam_5", name = paste0("Lrnr_gam_s", 5, "_x"))
)


make_ep_stack_sims <- function(base_learner, hte3_task, treatment_level = 1, control_level = 0, sieve_interaction_order = interaction_order) {
  d <- length(hte3_task$npsem$modifiers$variables)
  n <- hte3_task$nrow
  lrnr_crr_epstack <- Stack$new(Lrnr_crr_EP$new(base_learner = base_learner, sieve_num_basis  = 3*d,
                                                    sieve_interaction_order = sieve_interaction_order,
                                                    treatment_level = treatment_level, control_level = control_level),
                                  Lrnr_crr_EP$new(base_learner = base_learner, sieve_num_basis  = 4*d,
                                                  sieve_interaction_order = sieve_interaction_order,
                                                  treatment_level = treatment_level, control_level = control_level)
                                 , Lrnr_crr_EP$new(base_learner = base_learner, sieve_num_basis  = 5*d,
                                                    sieve_interaction_order = sieve_interaction_order,
                                                    treatment_level = treatment_level, control_level = control_level),
                                 Lrnr_crr_EP$new(base_learner = base_learner, sieve_num_basis = 6*d,
                                                  sieve_interaction_order = sieve_interaction_order,
                                                  treatment_level = treatment_level, control_level = control_level),
                                 Lrnr_crr_EP$new(base_learner = base_learner, sieve_num_basis = 7*d,
                                                  sieve_interaction_order = sieve_interaction_order,
                                                  treatment_level = treatment_level, control_level = control_level),
                                 Lrnr_crr_EP$new(base_learner = base_learner, sieve_num_basis = 8*d,
                                                  sieve_interaction_order = sieve_interaction_order,
                                                  treatment_level = treatment_level, control_level = control_level))
  return(lrnr_crr_epstack)
}


one_sim <- function(iter) {
  try({
    print(iter)
    print(key)
    set.seed(98103 + iter)
    data <- as.data.frame(sim.fun(n, hard, pos))
    set.seed(12345 + iter)
    data_val <- as.data.frame(sim.fun(n = 10000, hard, pos))
    # confounders for estimation of crr
    confounders <- grep("[W].+", colnames(data), value = TRUE)
    # nuisance learning ensemble for outcome regression and propensity score.
    lrnr_nuisance <- Stack$new(Lrnr_xgboost$new(max_depth = 2, min_child_weight = 10),
                               Lrnr_xgboost$new(max_depth = 3, min_child_weight = 10),
                               Lrnr_xgboost$new(max_depth = 4, min_child_weight = 10),
                               Lrnr_xgboost$new(max_depth = 5, min_child_weight = 10),
                               Lrnr_earth$new(degree = 2))

    # make hte3 task (with cross-fitted nuisances)
    hte3_task <- make_hte3_Task_tx(data, modifiers = confounders, confounders = confounders, treatment = "A", outcome = "Y",
                                   learner_pi = lrnr_nuisance,
                                   learner_mu = lrnr_nuisance,
                                   learner_m = lrnr_nuisance)



    # base learning algorithms for HTE meta-learners
    lrnr_xg_stack <- Stack$new(Lrnr_xgboost$new(max_depth = 1, min_child_weight = 5, objective = "reg:logistic"),
                               Lrnr_xgboost$new(max_depth = 2, min_child_weight = 5, objective = "reg:logistic"),
                               Lrnr_xgboost$new(max_depth = 3, min_child_weight = 5, objective = "reg:logistic"),
                               Lrnr_xgboost$new(max_depth = 4, min_child_weight = 5, objective = "reg:logistic"),
                               Lrnr_xgboost$new(max_depth = 5, min_child_weight = 5, objective = "reg:logistic"))

    lrnr_rf_stack <- Stack$new(Lrnr_xgboost$new(
      max_depth = 3, nrounds = 1,
      subsample = 0.5,
      colsample_bytree = 0.5,
      eta=1,
      lambda=0,
      num_parallel_tree = 500,
      verbosity = 0, name = "Lrnr_rf_3_xg", objective = "reg:logistic"),
      Lrnr_xgboost$new(
        max_depth = 5, nrounds = 1,
        subsample = 0.5,
        eta=1,
        lambda=0,
        colsample_bytree = 0.5,
        num_parallel_tree = 500,
        verbosity = 0, name = "Lrnr_rf_5_xg", objective = "reg:logistic"),
      Lrnr_xgboost$new(
        max_depth = 9, nround = 1,
        subsample = 0.5,
        eta=1,
        lambda=0,
        colsample_bytree = 0.5,
        num_parallel_tree = 500,
        verbosity = 0, name = "Lrnr_rf_9_xg", objective = "reg:logistic"),
      Lrnr_xgboost$new(
        max_depth = 7, nrounds = 1,
        subsample = 0.5,
        eta=1,
        lambda=0,
        colsample_bytree = 0.5,
        num_parallel_tree = 500,
        verbosity = 0, name = "Lrnr_rf_7_xg", objective = "reg:logistic") )




    lrnr_glm <-  Lrnr_glm$new(family = binomial())

    interaction_order <- 3
    if(base_learner %in% c("xg", "rf", "gam")) {
      # xg is a s tack, so cross-validate everything
      base_lrnr <- get(paste0("lrnr_", base_learner, "_stack"))
      ep_lrnr <- cross_validate_crr(make_ep_stack_sims(base_lrnr, hte3_task, sieve_interaction_order = interaction_order), hte3_task, cv_control = list(V = 5))$lrnr_sl
      ipw_lrnr <- cross_validate_crr(Lrnr_crr_IPW$new(base_lrnr), hte3_task, cv_control = list(V = 5))$lrnr_sl
      t_lrnr <- cross_validate_crr(Lrnr_crr_T$new(base_lrnr), hte3_task, cv_control = list(V = 5))$lrnr_sl
    } else {
      # Only need to cross-validate the EP-learner sieve
      base_lrnr <- get(paste0("lrnr_", base_learner))
      ep_lrnr <- cross_validate_crr(make_ep_stack_sims(base_lrnr, hte3_task, sieve_interaction_order = interaction_order), hte3_task,  cv_control = list(V = 5))$lrnr_sl
      ipw_lrnr <- Lrnr_crr_IPW$new(base_lrnr)$train(hte3_task)
      base_lrnr_poisson <- base_lrnr$reparameterize(list(family = poisson()))
      t_lrnr <- Lrnr_crr_T$new(base_lrnr_poisson)$train(hte3_task)
    }
    true_crr <- log(data_val$EY1W) - log(data_val$EY0W)
    ep_preds <- predict_hte3(ep_lrnr, new_data = data_val)
    ipw_preds <- predict_hte3(ipw_lrnr, new_data = data_val)
    t_preds <- predict_hte3(t_lrnr, new_data = data_val)

    ep_risk <- mean((true_crr - ep_preds)^2)
    ipw_risk <- mean((true_crr - ipw_preds)^2)
    t_risk <- mean((true_crr - t_preds)^2)



    out <- data.table(iter = iter, ep_risk = ep_risk, ipw_risk = ipw_risk, t_risk = t_risk , lrnr = base_learner)
    return(out)
  })
  return(data.table())
}

results <- rbindlist(lapply(1:nsims, one_sim))
fwrite(results, file = paste0("~/EPlearner/results/", key))
