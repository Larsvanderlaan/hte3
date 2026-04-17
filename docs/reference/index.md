# Package index

## Production API

- [`hte_task()`](https://larsvanderlaan.github.io/hte3/reference/hte_task.md)
  : Create an hte3 Task with a Production-Oriented Interface
- [`fit_cate()`](https://larsvanderlaan.github.io/hte3/reference/fit_cate.md)
  : Fit a CATE Model
- [`fit_crr()`](https://larsvanderlaan.github.io/hte3/reference/fit_crr.md)
  : Fit a CRR Model
- [`grf_cate()`](https://larsvanderlaan.github.io/hte3/reference/grf_cate.md)
  : Fit GRF-Backed CATE Models
- [`grf_crr()`](https://larsvanderlaan.github.io/hte3/reference/grf_crr.md)
  : Fit GRF-Backed CRR Models
- [`predict(`*`<hte3_model>`*`)`](https://larsvanderlaan.github.io/hte3/reference/hte3_model.md)
  [`summary(`*`<hte3_model>`*`)`](https://larsvanderlaan.github.io/hte3/reference/hte3_model.md)
  : S3 Methods for Production hte3 Models
- [`hte3_example_data()`](https://larsvanderlaan.github.io/hte3/reference/hte3_example_data.md)
  : Simulate Tiny Example Data for hte3

## Advanced API

- [`hte3_Task`](https://larsvanderlaan.github.io/hte3/reference/hte3_Task.md)
  : Task object for meta-learners in causal data structures.

- [`make_hte3_Task_tx()`](https://larsvanderlaan.github.io/hte3/reference/make_hte3_Task_tx.md)
  : Task object for meta-learners in the point-treatment setting.

- [`make_cross_fitted()`](https://larsvanderlaan.github.io/hte3/reference/make_cross_fitted.md)
  : Create Cross-Fitted Learner

- [`calibrate()`](https://larsvanderlaan.github.io/hte3/reference/calibrate.md)
  : Calibrate Predictor-Outcome Relationship

- [`estimate_m()`](https://larsvanderlaan.github.io/hte3/reference/estimate_m.md)
  :

  Helper function for estimation of treatment-averaged outcome
  regression using `sl3`.

- [`estimate_mu()`](https://larsvanderlaan.github.io/hte3/reference/estimate_mu.md)
  :

  Helper function for estimation of outcome regression using `sl3`.

- [`estimate_pi()`](https://larsvanderlaan.github.io/hte3/reference/estimate_pi.md)
  :

  Helper function for estimation of propensity score using `sl3`.

- [`make_ep_stack()`](https://larsvanderlaan.github.io/hte3/reference/make_ep_stack.md)
  : Create an Ensemble of CATE EP-learners with Varying Sieve Basis
  Sizes

- [`make_grf_nuisance_learners()`](https://larsvanderlaan.github.io/hte3/reference/make_grf_nuisance_learners.md)
  : Build Default GRF Nuisance Learners

- [`make_grf_cate_learners()`](https://larsvanderlaan.github.io/hte3/reference/make_grf_cate_learners.md)
  : Build GRF-Backed CATE Learners

- [`make_grf_crr_learners()`](https://larsvanderlaan.github.io/hte3/reference/make_grf_crr_learners.md)
  : Build GRF-Backed CRR Learners

- [`Lrnr_cate_DR`](https://larsvanderlaan.github.io/hte3/reference/Lrnr_cate_DR.md)
  : Lrnr_cate_DR Class

- [`Lrnr_cate_DR_selector`](https://larsvanderlaan.github.io/hte3/reference/Lrnr_cate_DR_selector.md)
  : Lrnr_cate_DR_selector Class

- [`Lrnr_cate_EP`](https://larsvanderlaan.github.io/hte3/reference/Lrnr_cate_EP.md)
  : Lrnr_cate_EP Class

- [`Lrnr_cate_R`](https://larsvanderlaan.github.io/hte3/reference/Lrnr_cate_R.md)
  : Lrnr_cate_R Class

- [`Lrnr_cate_T`](https://larsvanderlaan.github.io/hte3/reference/Lrnr_cate_T.md)
  : Lrnr_cate_T Class

- [`Lrnr_crr_DR_selector`](https://larsvanderlaan.github.io/hte3/reference/Lrnr_crr_DR_selector.md)
  : Lrnr_crr_DR_nonconvex Class

- [`Lrnr_crr_EP`](https://larsvanderlaan.github.io/hte3/reference/Lrnr_crr_EP.md)
  : Lrnr_crr_EP Class

- [`Lrnr_crr_IPW`](https://larsvanderlaan.github.io/hte3/reference/Lrnr_crr_IPW.md)
  : Lrnr_crr_IPW Class

- [`Lrnr_crr_T`](https://larsvanderlaan.github.io/hte3/reference/Lrnr_crr_T.md)
  : Lrnr_crr_T Class

- [`Lrnr_grf_causal_forest`](https://larsvanderlaan.github.io/hte3/reference/Lrnr_grf_causal_forest.md)
  : GRF Causal Forest Learner

- [`Lrnr_grf_forest`](https://larsvanderlaan.github.io/hte3/reference/Lrnr_grf_forest.md)
  : GRF Forest Learner

- [`Lrnr_stratified_multivariate`](https://larsvanderlaan.github.io/hte3/reference/Lrnr_stratified_multivariate.md)
  : Internal use. Converts a single outcome learner into a multivariate
  outcome learner that predicts a matrix of predictions obtained by
  evaluating the single outcome learner at each possible value of
  variable_stratify.

- [`Lrnr_xgboost_early_stopping`](https://larsvanderlaan.github.io/hte3/reference/Lrnr_xgboost_early_stopping.md)
  : XGBoost Learner with Early Stopping

- [`cross_validate()`](https://larsvanderlaan.github.io/hte3/reference/cross_validate.md)
  : Cross-Validate Heterogeneous Treatment Effect Models

- [`cross_validate_cate()`](https://larsvanderlaan.github.io/hte3/reference/cross_validate_cate.md)
  : Cross-Validate CATE Models

- [`cross_validate_crr()`](https://larsvanderlaan.github.io/hte3/reference/cross_validate_crr.md)
  : Cross-Validate CRR Models

- [`predict_hte3()`](https://larsvanderlaan.github.io/hte3/reference/predict_hte3.md)
  : Predict Heterogeneous Treatment Effects using hte3 Learners

- [`get_autoML()`](https://larsvanderlaan.github.io/hte3/reference/get_autoML.md)
  : Get Automated Machine Learning (AutoML) Learner

- [`truncate_propensity()`](https://larsvanderlaan.github.io/hte3/reference/truncate_propensity.md)
  : Truncate and Calibrate Propensity Scores
