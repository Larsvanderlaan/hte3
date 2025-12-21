library(SuperLearner)
library(sl3)
library(future)
library(npcausalML)
source("FinalSimulationCode/simRR.R")
nsims = 1000


library(earth)
SL.gam1 <- function(Y, X, newX, family, obsWeights, cts.num = 4,...) {
  deg.gam <- 1
  SuperLearner::SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num,... )
}
SL.gam2 <- function(Y, X, newX, family, obsWeights, cts.num = 4,...) {
  deg.gam <- 3
  SuperLearner::SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num,... )
}
SL.gam3 <- function(Y, X, newX, family, obsWeights, cts.num = 4,...) {
  deg.gam <- 5
  SuperLearner::SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num,... )
}
SL.gam4 <- function(Y, X, newX, family, obsWeights, cts.num = 4,...) {
  deg.gam <- 7
  SuperLearner::SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num,... )
}
SL.gam5 <- function(Y, X, newX, family, obsWeights, cts.num = 4,...) {
  deg.gam <- 9
  SuperLearner::SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num,... )
}


list_of_sieves_uni   <- list(
  "no_sieve" = NULL,
  fourier_basis$new(orders = c(1,0)),
  fourier_basis$new(orders = c(2,0)),
  fourier_basis$new(orders = c(3,0)),
  fourier_basis$new(orders = c(4,0))
)




onesim <- function(n) {

  sieve_list <- list_of_sieves_uni

  data <- as.data.frame(sim.RR(n, hard, pos))
  W <- data[,grep("^W", colnames(data))]
  A <- data$A
  Y <- data$Y
  EY1Wtrue <- data$EY1W
  EY0Wtrue <- data$EY0W
  pA1Wtrue <- data$pA1W
  EYWtrue <- ifelse(A==1, EY1Wtrue, EY0Wtrue)

  LRR <- log(EY1Wtrue/EY0Wtrue)



  # sieve method
  lrnr_Y <- make_learner(Pipeline, Lrnr_cv$new(Stack$new(
    Lrnr_xgboost$new(max_depth =2, verbosity = 0, nrounds = 10, objective = "binary:logistic"),
    Lrnr_xgboost$new(max_depth =3, verbosity = 0, nrounds = 10, objective = "binary:logistic"),
     Lrnr_xgboost$new(max_depth =5, verbosity = 0, nrounds = 10, objective = "binary:logistic"),
    Lrnr_xgboost$new(max_depth =6, verbosity = 0, nrounds = 10, objective = "binary:logistic")
  ))
  , Lrnr_cv_selector$new(loss_squared_error))


  lrnr_A <- make_learner(Pipeline, Lrnr_cv$new(
    Stack$new(


      Lrnr_xgboost$new(max_depth =2, verbosity = 0, nrounds = 10),
      Lrnr_xgboost$new(max_depth =3, verbosity = 0, nrounds = 10),
       Lrnr_xgboost$new(max_depth =5, verbosity = 0, nrounds = 10),
      Lrnr_xgboost$new(max_depth =6, verbosity = 0, nrounds = 10)
    )
  ), Lrnr_cv_selector$new(loss_squared_error))


  data_train <-  data #as.data.frame(sim.LRR(n, hard, pos))

  initial_likelihood <- npcausalML:::estimate_initial_likelihood(W=data_train[,c("W1", "W2", "W3"), drop = F], data_train$A, data_train$Y,  weights = rep(1,n), lrnr_A, lrnr_Y, folds = 10, outcome_type = "continuous")
  data1 <- data
  data0 <- data
  data1$A <- 1
  data0$A <- 0
  taskY <- sl3_Task$new(data, covariates = c("W1", "W2", "W3", "A"), outcome = "Y", folds = origami::folds_vfold(n), outcome_type = "continuous")
  folds <- taskY$folds
  taskY0 <- sl3_Task$new(data0, covariates = c("W1", "W2", "W3", "A"), outcome = "Y", folds = folds, outcome_type = "continuous")
  taskY1 <- sl3_Task$new(data1, covariates = c("W1", "W2", "W3", "A"), outcome = "Y", folds = folds, outcome_type = "continuous")
  taskA <- sl3_Task$new(data, covariates = c("W1", "W2", "W3"), outcome = "A", folds = folds)

  pA1W_est <- initial_likelihood$internal$sl3_Learner_pA1W_trained$predict(taskA)
  EY1W_est <- initial_likelihood$internal$sl3_Learner_EYAW_trained$predict(taskY1)
  EY0W_est <- initial_likelihood$internal$sl3_Learner_EYAW_trained$predict(taskY0)

  pA1W_est <- pmax(pA1W_est, 0.01)
  pA1W_est <- pmin(pA1W_est, 0.99)

  data$weightsIPW <- data$Y/ifelse(data$A==1,pA1W_est, 1 - pA1W_est)
  sl3_Task_IPW <- sl3_Task$new(data, covariates = c("W1", "W2", "W3"), outcome = "A", weights = "weightsIPW", folds = 5)


  lrnr_xg <- list(    Lrnr_xgboost$new(max_depth = 1, verbosity = 0, nrounds = 10, objective = "reg:logistic"),    Lrnr_xgboost$new(max_depth = 2, verbosity = 0, nrounds = 10, objective = "reg:logistic"), Lrnr_xgboost$new(max_depth = 3, verbosity = 0, nrounds = 10, objective = "reg:logistic"), Lrnr_xgboost$new(max_depth = 4, verbosity = 0, nrounds = 10, objective = "reg:logistic"),   Lrnr_xgboost$new(max_depth = 5, verbosity = 0, nrounds = 10, objective = "reg:logistic" ))
  lrnr_xg_sl <-  Lrnr_sl$new(Stack$new(lrnr_xg), metalearner = Lrnr_cv_selector$new(loss_loglik_binomial))


  lrnr_rf <- list(Lrnr_xgboost$new(
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


  lrnr_rf_sl <-  Lrnr_sl$new(Stack$new(lrnr_rf), metalearner = Lrnr_cv_selector$new(loss_loglik_binomial))




  LRR_library <- c(lrnr_rf,   lrnr_xg  )
  LRR_library_sl <-  c(lrnr_rf , lrnr_xg, list(lrnr_xg_sl , lrnr_rf_sl  ))
  LRR_library_subst <- c(lrnr_rf , lrnr_xg, list(lrnr_xg_sl , lrnr_rf_sl  ))



  IPW_learner <- Stack$new(LRR_library_sl) #TODO
  IPW_learner <- IPW_learner$train(sl3_Task_IPW)
  preds_IPW <- apply(IPW_learner$predict(sl3_Task_IPW), 2, qlogis)


  subst_compare <- Stack$new(LRR_library_subst)
   subst_EY1W_trained <-subst_compare$train(taskY1[A==1]$next_in_chain(covariates = c("W1", "W2", "W3")))
  subst_EY0W_trained <- subst_compare$train(taskY0[A==0]$next_in_chain(covariates = c("W1", "W2", "W3")))

  subst_EY1W <-subst_EY1W_trained$predict(taskY1$next_in_chain(covariates = c("W1", "W2", "W3")))
  subst_EY0W <- subst_EY0W_trained$predict(taskY0$next_in_chain(covariates = c("W1", "W2", "W3")))


  subst_EY1W <- pmax(subst_EY1W, 1e-5)
  subst_EY0W <- pmax(subst_EY0W, 1e-5)
  subst_LRR <- log(subst_EY1W/ subst_EY0W)

    print("EP")

  fit_npcausalML <- EP_learn(LRR_library,V =  W, A = A, Y = Y, EY1W = EY1W_est  , EY0W = EY0W_est  , pA1W = pA1W_est, sieve_basis_generator_list = sieve_list ,EP_learner_spec = EP_learner_spec_LRR, cross_validate = TRUE, nfolds = 2)
  preds <- fit_npcausalML$full_predictions
print("done")



  # Compute least-squares risk of predictions using oracle loss function.
  risks_oracle <- as.vector(apply(preds, 2, function(theta) {
    mean((theta -  LRR)^2)
  }))

  # Compute estimated cross-validated one-step risk of predictions
  cvrisksDR <- as.vector(apply(fit_npcausalML$cv_predictions, 2, function(theta) {
    loss <- efficient_loss_function_LRR(theta, A, Y, EY1W_est,EY0W_est, pA1W_est )
    mean(loss)
  }))#[-grep("IPW", colnames(fit_npcausalML$cv_predictions))])

  # Compute estimated cross-validated oracle one-step risk of predictions

  cvrisksDRoracle <- as.vector(apply(fit_npcausalML$cv_predictions, 2, function(theta) {
    loss <- efficient_loss_function_LRR(theta, A, Y, EY1Wtrue,EY0Wtrue, pA1Wtrue )
    mean(loss)
  }))#[-grep("IPW", colnames(fit_npcausalML$cv_predictions))])
  lrnrs_full <-  colnames(fit_npcausalML$cv_predictions)
  lrnrs <- gsub("[._]fourier.+", "", lrnrs_full)
  lrnrs <- gsub("_no_sieve", "", lrnrs)
  degree <- as.numeric(stringr::str_match(lrnrs_full, "fourier_basis_([0-9]+)")[,2])
  degree[grep("no_sieve", lrnrs_full)] <- 0



  tmp <- data.table(lrnrs_full, lrnrs , degree, risk = cvrisksDR, risks_oracle = risks_oracle, cvrisksDR = cvrisksDR, cvrisksDRoracle)

  rf_keep <- tmp[grep("rf", lrnrs_full),risks_oracle[which.min(risk)], by = degree]$V1
  xg_keep <- tmp[grep("xgboost", lrnrs_full),risks_oracle[which.min(risk)], by = degree]$V1
  names(xg_keep) <- paste0("Lrnr_xgboost_cv", "_fourier.basis_", 0:4, "_plugin")
  names(rf_keep) <- paste0("Lrnr_rf_cv", "_fourier.basis_", 0:4, "_plugin")
  risks_oracle <- c(risks_oracle, xg_keep, rf_keep)

  rf_keep <- tmp[grep("rf", lrnrs_full),cvrisksDR[which.min(risk)], by = degree]$V1
  xg_keep <- tmp[grep("xgboost", lrnrs_full),cvrisksDR[which.min(risk)], by = degree]$V1
  names(xg_keep) <- paste0("Lrnr_xgboost_cv", "_fourier.basis_", 0:4, "_plugin")
  names(rf_keep) <- paste0("Lrnr_rf_cv", "_fourier.basis_", 0:4, "_plugin")
  cvrisksDR <- c(cvrisksDR, xg_keep, rf_keep)

  rf_keep <- tmp[grep("rf", lrnrs_full),cvrisksDRoracle[which.min(risk)], by = degree]$V1
  xg_keep <- tmp[grep("xgboost", lrnrs_full),cvrisksDRoracle[which.min(risk)], by = degree]$V1
  names(xg_keep) <- paste0("Lrnr_xgboost_cv", "_fourier.basis_", 0:4, "_plugin")
  names(rf_keep) <- paste0("Lrnr_rf_cv", "_fourier.basis_", 0:4, "_plugin")
  cvrisksDRoracle <- c(cvrisksDRoracle, xg_keep, rf_keep)

  sieve_names <- c(colnames(fit_npcausalML$cv_predictions),   names(xg_keep),   names(rf_keep))



  risk_subst_cv <- mean((log(EY1W_est / EY0W_est) - LRR)^2)






  risk_subst<-  apply(subst_LRR, 2, function(pred) {
    mean((pred - LRR)^2)
  })

  risk_IPW <-  apply(preds_IPW, 2, function(pred) {
    mean((pred - LRR)^2)
  })


  list(risk_subst_cv = risk_subst_cv, risk_IPW = risk_IPW,  risk_subst = risk_subst,    sieve =data.frame(sieve_names, cvrisksDRoracle, cvrisksDR, risks_oracle))
}


hard <- hard == "TRUE"
pos <- pos == "TRUE"
n <- as.numeric(n)

simresults <- lapply(1:nsims, function(i){try({
  print(i)
  onesim(n)
})
})
save(simresults, file = paste0("mainSimResults/","simsLRR", hard,pos, "n", n, "_xgboost"))


