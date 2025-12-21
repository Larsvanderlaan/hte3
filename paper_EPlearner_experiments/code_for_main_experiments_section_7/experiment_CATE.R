# ----------------------------
# Reproducible paths via BASE_PATH
# (matches your new convention)
#   SIM_CODE_DIR <- file.path(BASE_PATH, "code_for_experiments_and_figures")
#   RESULTS_DIR  <- file.path(BASE_PATH, "experiment_results")
# ----------------------------

# generate synthetic dataset
nsims <- 1000

n <- as.numeric(n)
hard <- as.logical(hard)
pos <- as.logical(pos)
sim_type <- as.character(sim_type)
base_learner <- as.character(base_learner)

key <- paste0(
  "simsEP_n=", n,
  "_pos=", pos,
  "_hard=", hard,
  "_base_learner=", base_learner,
  "_sim_type=", sim_type,
  "_nsims=", nsims
)

library(data.table)
library(sl3)
library(hte3)
library(grf)
library(SuperLearner)

# ---- base path (single source of truth) ----
BASE_PATH <- Sys.getenv("BASE_PATH", unset = NA_character_)
if (is.na(BASE_PATH) || !nzchar(BASE_PATH)) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/"), error = function(e) NA_character_)
  if (is.na(this_file)) {
    stop("Set BASE_PATH env var, e.g. Sys.setenv(BASE_PATH='/path/to/project_root')")
  }
  BASE_PATH <- normalizePath(dirname(this_file), winslash = "/")
} else {
  BASE_PATH <- normalizePath(path.expand(BASE_PATH), winslash = "/")
}

SIM_CODE_DIR <- file.path(BASE_PATH, "code_for_experiments_and_figures")
RESULTS_DIR  <- file.path(BASE_PATH, "experiment_results")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Generate data: source from SIM_CODE_DIR (no ~ paths)
# ----------------------------
if (sim_type == "CATEhigh") {
  source(file.path(SIM_CODE_DIR, "simCATEHighDim.R"))
  nsims <- 500
  sim.fun <- sim.CATEHighDim
} else if (sim_type == "CATElow") {
  source(file.path(SIM_CODE_DIR, "simCATE.R"))
  sim.fun <- sim.CATE
} else if (sim_type == "CRATE") {
  # NOTE: assumes sim.LRR is defined in some sourced file.
  # If it's also in SIM_CODE_DIR, source it here similarly.
  sim.fun <- sim.LRR
} else {
  stop(sprintf("Unknown sim_type: %s", sim_type))
}

# (Optional) since nsims can change above, you may want to refresh key:
key <- paste0(
  "simsEP_n=", n,
  "_pos=", pos,
  "_hard=", hard,
  "_base_learner=", base_learner,
  "_sim_type=", sim_type,
  "_nsims=", nsims
)

# ----------------------------
# SuperLearner wrappers
# ----------------------------
SL.gam_2 <- function(Y, X, newX, family, obsWeights, cts.num = 4, ...) {
  deg.gam <- 2
  SuperLearner::SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num, ...)
}
SL.gam_3 <- function(Y, X, newX, family, obsWeights, cts.num = 4, ...) {
  deg.gam <- 3
  SuperLearner::SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num, ...)
}
SL.gam_5 <- function(Y, X, newX, family, obsWeights, cts.num = 4, ...) {
  deg.gam <- 5
  SuperLearner::SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num, ...)
}

lrnr_gam_stack <- Stack$new(
  Lrnr_pkg_SuperLearner$new("SL.gam_2", name = paste0("Lrnr_gam_s", 2, "_x")),
  Lrnr_pkg_SuperLearner$new("SL.gam_3", name = paste0("Lrnr_gam_s", 3, "_x")),
  Lrnr_pkg_SuperLearner$new("SL.gam_3", name = paste0("Lrnr_gam_s", 4, "_x")),
  Lrnr_pkg_SuperLearner$new("SL.gam_5", name = paste0("Lrnr_gam_s", 5, "_x"))
)

make_ep_stack_sims_low <- function(base_learner, hte3_task,
                                   treatment_level = 1, control_level = 0,
                                   sieve_interaction_order = interaction_order) {
  d <- length(hte3_task$npsem$modifiers$variables)
  n <- hte3_task$nrow
  Stack$new(
    Lrnr_cate_EP$new(base_learner = base_learner, sieve_num_basis = 3 * d,
                     sieve_interaction_order = sieve_interaction_order,
                     treatment_level = treatment_level, control_level = control_level),
    Lrnr_cate_EP$new(base_learner = base_learner, sieve_num_basis = 4 * d,
                     sieve_interaction_order = sieve_interaction_order,
                     treatment_level = treatment_level, control_level = control_level),
    Lrnr_cate_EP$new(base_learner = base_learner, sieve_num_basis = 5 * d,
                     sieve_interaction_order = sieve_interaction_order,
                     treatment_level = treatment_level, control_level = control_level),
    Lrnr_cate_EP$new(base_learner = base_learner, sieve_num_basis = 6 * d,
                     sieve_interaction_order = sieve_interaction_order,
                     treatment_level = treatment_level, control_level = control_level),
    Lrnr_cate_EP$new(base_learner = base_learner, sieve_num_basis = 7 * d,
                     sieve_interaction_order = sieve_interaction_order,
                     treatment_level = treatment_level, control_level = control_level),
    Lrnr_cate_EP$new(base_learner = base_learner, sieve_num_basis = 8 * d,
                     sieve_interaction_order = sieve_interaction_order,
                     treatment_level = treatment_level, control_level = control_level)
  )
}

make_ep_stack_sims_high <- function(base_learner, hte3_task,
                                    treatment_level = 1, control_level = 0,
                                    sieve_interaction_order = interaction_order) {
  d <- length(hte3_task$npsem$modifiers$variables)
  n <- hte3_task$nrow
  Stack$new(
    Lrnr_cate_EP$new(base_learner = base_learner, sieve_num_basis = 1 * d,
                     sieve_interaction_order = sieve_interaction_order,
                     treatment_level = treatment_level, control_level = control_level),
    Lrnr_cate_EP$new(base_learner = base_learner, sieve_num_basis = 2 * d,
                     sieve_interaction_order = sieve_interaction_order,
                     treatment_level = treatment_level, control_level = control_level),
    Lrnr_cate_EP$new(base_learner = base_learner, sieve_num_basis = 3 * d,
                     sieve_interaction_order = sieve_interaction_order,
                     treatment_level = treatment_level, control_level = control_level),
    Lrnr_cate_EP$new(base_learner = base_learner, sieve_num_basis = 4 * d,
                     sieve_interaction_order = sieve_interaction_order,
                     treatment_level = treatment_level, control_level = control_level),
    Lrnr_cate_EP$new(base_learner = base_learner, sieve_num_basis = 5 * d,
                     sieve_interaction_order = sieve_interaction_order,
                     treatment_level = treatment_level, control_level = control_level),
    Lrnr_cate_EP$new(base_learner = base_learner, sieve_num_basis = 6 * d,
                     sieve_interaction_order = sieve_interaction_order,
                     treatment_level = treatment_level, control_level = control_level)
  )
}

one_sim <- function(iter) {
  ({
    print(iter)
    print(key)

    set.seed(98103 + iter)
    data <- as.data.frame(sim.fun(n, hard, pos))

    set.seed(12345 + iter)
    data_val <- as.data.frame(sim.fun(n = 10000, hard, pos))

    confounders <- grep("[W].+", colnames(data), value = TRUE)

    lrnr_nuisance <- Stack$new(
      Lrnr_xgboost$new(max_depth = 2, min_child_weight = 10),
      Lrnr_xgboost$new(max_depth = 3, min_child_weight = 10),
      Lrnr_xgboost$new(max_depth = 4, min_child_weight = 10),
      Lrnr_xgboost$new(max_depth = 5, min_child_weight = 10),
      Lrnr_earth$new(degree = 2)
    )

    hte3_task <- make_hte3_Task_tx(
      data,
      modifiers = confounders, confounders = confounders,
      treatment = "A", outcome = "Y",
      learner_pi = lrnr_nuisance,
      learner_mu = lrnr_nuisance,
      learner_m  = lrnr_nuisance
    )

    lrnr_xg_stack <- Stack$new(
      Lrnr_xgboost$new(max_depth = 1, min_child_weight = 3),
      Lrnr_xgboost$new(max_depth = 2, min_child_weight = 3),
      Lrnr_xgboost$new(max_depth = 3, min_child_weight = 3),
      Lrnr_xgboost$new(max_depth = 4, min_child_weight = 3),
      Lrnr_xgboost$new(max_depth = 5, min_child_weight = 3)
    )
    lrnr_rf_stack <- Stack$new(
      Lrnr_ranger$new(max.depth = 6),
      Lrnr_ranger$new(max.depth = 8),
      Lrnr_ranger$new(max.depth = 10),
      Lrnr_ranger$new(max.depth = 12)
    )

    if (sim_type == "CATEhigh") {
      make_ep_stack_sims <- make_ep_stack_sims_high
      interaction_order <- 1
      lrnr_earth <- Lrnr_earth$new(degree = 1)
    } else {
      make_ep_stack_sims <- make_ep_stack_sims_low
      interaction_order <- 3
      lrnr_earth <- Lrnr_earth$new(degree = 2)
    }

    if (base_learner %in% c("xg", "rf", "gam")) {
      base_lrnr <- get(paste0("lrnr_", base_learner, "_stack"))
      ep_lrnr <- cross_validate_cate(
        make_ep_stack_sims(base_lrnr, hte3_task, sieve_interaction_order = interaction_order),
        hte3_task, cv_control = list(V = 5)
      )$lrnr_sl
      dr_lrnr <- cross_validate_cate(Lrnr_cate_DR$new(base_lrnr), hte3_task, cv_control = list(V = 5))$lrnr_sl
      r_lrnr  <- cross_validate_cate(Lrnr_cate_R$new(base_lrnr),  hte3_task, cv_control = list(V = 5))$lrnr_sl
      t_lrnr  <- cross_validate_cate(Lrnr_cate_T$new(base_lrnr),  hte3_task, cv_control = list(V = 5))$lrnr_sl
    } else {
      base_lrnr <- get(paste0("lrnr_", base_learner))
      dr_lrnr <- Lrnr_cate_DR$new(base_lrnr)$train(hte3_task)
      ep_fit <- cross_validate_cate(
        make_ep_stack_sims(base_lrnr, hte3_task, sieve_interaction_order = interaction_order),
        hte3_task, cv_control = list(V = 5)
      )
      print(ep_fit$coef)
      print(ep_fit$cv_risk)
      ep_lrnr <- ep_fit$lrnr_sl

      if (base_learner == "earth") {
        r_lrnr <- Lrnr_cate_R$new(Lrnr_mean$new())$train(hte3_task)
      } else {
        r_lrnr <- Lrnr_cate_R$new(base_lrnr)$train(hte3_task)
      }
      t_lrnr <- Lrnr_cate_T$new(base_lrnr)$train(hte3_task)
    }

    true_cate <- data_val$EY1W - data_val$EY0W

    ep_preds <- predict_hte3(ep_lrnr, new_data = data_val)
    dr_preds <- predict_hte3(dr_lrnr, new_data = data_val)
    r_preds  <- predict_hte3(r_lrnr,  new_data = data_val)
    t_preds  <- predict_hte3(t_lrnr,  new_data = data_val)

    ep_risk <- mean((true_cate - ep_preds)^2)
    dr_risk <- mean((true_cate - dr_preds)^2)
    r_risk  <- mean((true_cate - r_preds)^2)
    t_risk  <- mean((true_cate - t_preds)^2)

    if (base_learner %in% c("rf", "xg")) {
      X <- hte3_task$get_tmle_node("modifiers")
      W <- hte3_task$get_tmle_node("treatment")
      Y <- hte3_task$get_tmle_node("outcome")

      W.hat <- as.matrix(hte3_task$get_nuisance_estimates("pi"))[, 2]
      W.hat <- causalutils::truncate_propensity(W.hat, W, treatment_level = 1, truncation_method = "adaptive")
      Y.hat <- as.vector(hte3_task$get_nuisance_estimates("m"))

      cf_fit <- grf::causal_forest(
        X = as.matrix(X),
        Y = as.vector(as.matrix(Y)),
        W = as.vector(as.matrix(W)),
        W.hat = W.hat,
        Y.hat = Y.hat,
        tune.parameters = "all"
      )

      preds_cf <- as.vector(predict(cf_fit, newdata = as.matrix(data_val)[, confounders])$predictions)
      causal_forest_risk <- mean((true_cate - preds_cf)^2)
    } else {
      causal_forest_risk <- NA_real_
    }

    data.table(
      iter = iter,
      ep_risk = ep_risk,
      dr_risk = dr_risk,
      r_risk = r_risk,
      t_risk = t_risk,
      causal_forest_risk = causal_forest_risk,
      lrnr = base_learner
    )
  })
}

results <- rbindlist(lapply(seq_len(nsims), one_sim))

# Write using RESULTS_DIR (no hard-coded ~ paths)
fwrite(results, file = file.path(RESULTS_DIR, paste0(key, ".csv")))
