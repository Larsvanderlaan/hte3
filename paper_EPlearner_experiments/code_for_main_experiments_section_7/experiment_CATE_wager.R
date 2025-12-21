## ---- parse CLI args and set params ----
args <- commandArgs(trailingOnly = TRUE)

parse_kv <- function(x) {
  kv <- strsplit(x, "=", fixed = TRUE)[[1]]
  if (length(kv) == 2) setNames(list(kv[2]), kv[1]) else NULL
}
args_named <- do.call(c, Filter(Negate(is.null), lapply(args, parse_kv)))

# Required
n <- as.numeric(args_named[["n"]])
stopifnot(!is.na(n))

# Optional with defaults
d            <- if (!is.null(args_named[["d"]]))            as.numeric(args_named[["d"]])            else 5
sigma        <- if (!is.null(args_named[["sigma"]]))        as.numeric(args_named[["sigma"]])        else 1
base_learner <- if (!is.null(args_named[["base_learner"]])) as.character(args_named[["base_learner"]]) else "xg"
sim_type     <- if (!is.null(args_named[["sim_type"]]))     as.character(args_named[["sim_type"]])     else "A"
nsims        <- if (!is.null(args_named[["nsims"]]))        as.integer(args_named[["nsims"]])        else 250

message(sprintf(
  "Parsed args: n=%s, d=%s, sigma=%s, base_learner=%s, sim_type=%s, nsims=%s",
  n, d, sigma, base_learner, sim_type, nsims
))

# Parameter-specific key
key <- sprintf(
  "simsEP_n=%s_d=%s_sigma=%s_base_learner=%s_sim_type=%s_nsims=%s",
  n, d, sigma, base_learner, sim_type, nsims
)

# Libraries
library(data.table)
library(sl3)
library(hte3)
library(grf)
library(SuperLearner)

# ----------------------------
# Reproducible paths via BASE_PATH
# ----------------------------
BASE_PATH <- Sys.getenv("BASE_PATH", unset = NA_character_)
if (is.na(BASE_PATH) || !nzchar(BASE_PATH)) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/"), error = function(e) NA_character_)
  if (is.na(this_file)) {
    stop("Set BASE_PATH env var, e.g. Sys.setenv(BASE_PATH='/path/to/paper_EPlearner_experiments')")
  }
  BASE_PATH <- normalizePath(dirname(this_file), winslash = "/")
} else {
  BASE_PATH <- normalizePath(path.expand(BASE_PATH), winslash = "/")
}

# You said you changed these:
SIM_CODE_DIR <- file.path(BASE_PATH, "code_for_experiments_and_figures")
RESULTS_DIR  <- file.path(BASE_PATH, "experiment_results")

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# (Optional) if you source any helpers from SIM_CODE_DIR, do it here, e.g.:
# source(file.path(SIM_CODE_DIR, "some_helper.R"))

# ---------- Wager/Athey simulator (Setups A/B/C/D) ----------
simulate_wager <- function(n, d = 5, sigma = 1, setup = c("A","B","C","D")) {
  setup <- match.arg(setup)

  trimo <- function(x) pmin(pmax(x, 0.1), 0.9)

  if (setup == "A") {
    X <- matrix(runif(n * d, 0, 1), n, d)
    b_fun   <- function(x) sin(pi * x[,1] * x[,2]) + 2 * (x[,3] - 0.5)^2 + x[,4] + 0.5 * x[,5]
    e_fun   <- function(x) trimo(sin(pi * x[,1] * x[,2]))
    tau_fun <- function(x) (x[,1] + x[,2]) / 2
  } else if (setup == "B") {
    X <- matrix(rnorm(n * d), n, d)
    b_fun   <- function(x) pmax(x[,1] + x[,2], x[,3], 0) + pmax(x[,4] + x[,5], 0)
    e_fun   <- function(x) rep(0.5, n)
    tau_fun <- function(x) x[,1] + log1p(exp(x[,2]))
  } else if (setup == "C") {
    X <- matrix(rnorm(n * d), n, d)
    b_fun   <- function(x) 2 * log1p(exp(x[,1] + x[,2] + x[,3]))
    e_fun   <- function(x) 1 / (1 + exp(-(x[,2] + x[,3])))
    tau_fun <- function(x) rep(1, n)
  } else {
    X <- matrix(rnorm(n * d), n, d)
    b_fun   <- function(x) (pmax(x[,1] + x[,2] + x[,3], 0) + pmax(x[,4] + x[,5], 0)) / 2
    e_fun   <- function(x) 1 / (1 + exp(-x[,1] - x[,2]))
    tau_fun <- function(x) pmax(x[,1] + x[,2] + x[,3], 0) - pmax(x[,4] + x[,5], 0)
  }

  bX   <- b_fun(X)
  eX   <- e_fun(X)
  tauX <- tau_fun(X)

  A <- rbinom(n, 1, eX)
  Y <- bX + (A - 0.5) * tauX + sigma * rnorm(n)

  EY1W <- bX + 0.5 * tauX
  EY0W <- bX - 0.5 * tauX

  dat <- data.frame(
    Y = Y, A = A,
    setNames(as.data.frame(X), paste0("W", seq_len(ncol(X)))),
    EY1W = EY1W, EY0W = EY0W
  )
  return(dat)
}

# Generate data wrapper (kept signature-compatible)
sim.fun <- function(n, hard, pos, d = 5, sigma = 1) {
  simulate_wager(n = n, d = d, sigma = sigma, setup = sim_type)
}

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

lrnr_mean_stack <- Stack$new(Lrnr_mean$new(), Lrnr_mean$new())

make_ep_stack_sims <- function(base_learner, hte3_task, treatment_level = 1, control_level = 0,
                               sieve_interaction_order = 2) {
  d <- length(hte3_task$npsem$modifiers$variables)
  n <- hte3_task$nrow

  lrnr_cate_epstack <- Stack$new(
    Lrnr_cate_EP$new(base_learner = base_learner, sieve_num_basis = 1,
                     sieve_interaction_order = 1,
                     treatment_level = treatment_level, control_level = control_level),
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
                     treatment_level = treatment_level, control_level = control_level)
  )
  return(lrnr_cate_epstack)
}

# ---- parameterized one_sim ----
one_sim <- function(iter, sim_type, n, d = 5, sigma = 1, base_learner = "xg") {
  ({
    key_local <- paste0(
      "simsEP_n=", n,
      "_base_learner=", base_learner,
      "_sim_type=", sim_type,
      "_nsims=", nsims
    )

    print(paste("iter", iter, "key", key_local))

    set.seed(981037 + iter)
    sim.fun_local <- function(n) simulate_wager(n = n, d = d, sigma = sigma, setup = sim_type)

    data <- sim.fun_local(n)
    set.seed(123457 + iter)
    data_val <- sim.fun_local(10000)

    confounders <- grep("[W].+", colnames(data), value = TRUE)

    lrnr_nuisance <- Stack$new(
      Lrnr_xgboost$new(max_depth = 2, min_child_weight = 10),
      Lrnr_xgboost$new(max_depth = 3, min_child_weight = 10),
      Lrnr_xgboost$new(max_depth = 4, min_child_weight = 10),
      Lrnr_xgboost$new(max_depth = 5, min_child_weight = 10),
      Lrnr_earth$new(degree = 2)
    )

    hte3_task <- make_hte3_Task_tx(
      data, modifiers = confounders, confounders = confounders,
      treatment = "A", outcome = "Y",
      learner_pi = lrnr_nuisance,
      learner_mu = lrnr_nuisance,
      learner_m  = lrnr_nuisance
    )

    lrnr_xg_stack <- Stack$new(
      Lrnr_xgboost$new(max_depth = 1, min_child_weight = 5),
      Lrnr_xgboost$new(max_depth = 2, min_child_weight = 5),
      Lrnr_xgboost$new(max_depth = 3, min_child_weight = 5),
      Lrnr_xgboost$new(max_depth = 4, min_child_weight = 5),
      Lrnr_xgboost$new(max_depth = 5, min_child_weight = 5)
    )
    lrnr_rf_stack <- Stack$new(
      Lrnr_ranger$new(max.depth = 1),
      Lrnr_ranger$new(max.depth = 6),
      Lrnr_ranger$new(max.depth = 8),
      Lrnr_ranger$new(max.depth = 10)
    )

    interaction_order <- 2

    if (base_learner %in% c("xg", "rf", "gam", "mean")) {
      base_lrnr <- get(paste0("lrnr_", base_learner, "_stack"))
      ep_lrnr <- cross_validate_cate(
        make_ep_stack_sims(base_lrnr, hte3_task, sieve_interaction_order = interaction_order),
        hte3_task, cv_control = list(V = 10)
      )$lrnr_sl
      dr_lrnr <- cross_validate_cate(Lrnr_cate_DR$new(base_lrnr), hte3_task, cv_control = list(V = 10))$lrnr_sl
      r_lrnr  <- cross_validate_cate(Lrnr_cate_R$new(base_lrnr),  hte3_task, cv_control = list(V = 10))$lrnr_sl
      t_lrnr  <- cross_validate_cate(Lrnr_cate_T$new(base_lrnr),  hte3_task, cv_control = list(V = 10))$lrnr_sl
    } else {
      base_lrnr <- get(paste0("lrnr_", base_learner))
      dr_lrnr <- Lrnr_cate_DR$new(base_lrnr)$train(hte3_task)
      ep_lrnr <- cross_validate_cate(
        make_ep_stack_sims(base_lrnr, hte3_task, sieve_interaction_order = interaction_order),
        hte3_task, cv_control = list(V = 10)
      )$lrnr_sl
      r_lrnr <- Lrnr_cate_R$new(base_lrnr)$train(hte3_task)
      t_lrnr <- Lrnr_cate_T$new(base_lrnr)$train(hte3_task)
    }

    true_cate <- data_val$EY1W - data_val$EY0W
    ep_preds  <- predict_hte3(ep_lrnr, new_data = data_val)
    dr_preds  <- predict_hte3(dr_lrnr, new_data = data_val)
    r_preds   <- predict_hte3(r_lrnr,  new_data = data_val)
    t_preds   <- predict_hte3(t_lrnr,  new_data = data_val)

    if (identical(base_learner, "rf")) {
      X <- as.matrix(hte3_task$get_tmle_node("modifiers"))
      W <- as.vector(hte3_task$get_tmle_node("treatment"))
      Y <- as.vector(hte3_task$get_tmle_node("outcome"))

      nuis_pi <- hte3_task$get_nuisance_estimates("pi")
      W.hat <- if (is.matrix(nuis_pi)) nuis_pi[, 2] else as.numeric(nuis_pi)
      W.hat <- pmin(pmax(W.hat, 1e-3), 1 - 1e-3)

      Y.hat <- as.vector(hte3_task$get_nuisance_estimates("m"))

      cf_fit <- grf::causal_forest(
        X = X, Y = Y, W = W,
        W.hat = W.hat, Y.hat = Y.hat,
        tune.parameters = "all"
      )

      preds_cf <- as.vector(predict(cf_fit, newdata = as.matrix(data_val[, confounders]))$predictions)
      causal_forest_risk <- mean((true_cate - preds_cf)^2)
    } else {
      causal_forest_risk <- NA_real_
    }

    data.table(
      iter = iter,
      sim_type = sim_type,
      n = n,
      d = d,
      sigma = sigma,
      lrnr = base_learner,
      ep_risk = mean((true_cate - ep_preds)^2),
      dr_risk = mean((true_cate - dr_preds)^2),
      r_risk  = mean((true_cate - r_preds)^2),
      t_risk  = mean((true_cate - t_preds)^2),
      causal_forest_risk = causal_forest_risk
    )
  })
}

# ---- run nsims with parsed args ----
results <- rbindlist(lapply(seq_len(nsims), function(iter) {
  one_sim(iter, sim_type = sim_type, n = n, d = d, sigma = sigma, base_learner = base_learner)
}))

# ----------------------------
# WRITE OUTPUT using RESULTS_DIR (no ~, no hard-coded paths)
# ----------------------------
fwrite(results, file = file.path(RESULTS_DIR, paste0(key, ".csv")))

# ---- ad hoc local run (also uses RESULTS_DIR if you write) ----
nsims <- 20
results <- rbindlist(lapply(seq_len(nsims), function(iter) {
  one_sim(iter, sim_type = "C", n = 1000, d = 5, sigma = 1, base_learner = "rf")
}))

DT_avg <- results[, .(
  ep_MSE = mean(ep_risk, na.rm = TRUE),
  dr_MSE = mean(dr_risk, na.rm = TRUE),
  r_MSE  = mean(r_risk,  na.rm = TRUE),
  t_MSE  = mean(t_risk,  na.rm = TRUE),
  cf_MSE = mean(causal_forest_risk, na.rm = TRUE)
)]

print(DT_avg)

# If you want to save the averaged table too:
# fwrite(DT_avg, file = file.path(RESULTS_DIR, paste0(key, "_avg.csv")))
