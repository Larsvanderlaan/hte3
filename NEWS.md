# hte3 0.1.1

- Removed internal fold-prediction and pseudo-outcome type warnings, and
  avoided unused `family` arguments for forest and earth learners. Cross-fitted
  predictions and warnings about model fitting or target changes are unchanged.
- Hardened the default `get_autoML()` stack so the high-level API uses the
  available safe subset of learners and warns once when optional runtime
  packages are missing.
- Finished the two-stage T-learner update for reduced-modifier settings so
  first-stage regressions use the confounder set `W` and the default second
  stage projects onto the modifier set `V`.
- Tightened task validation for folds, weights, treatment values, and
  user-supplied nuisance estimates.
- Added release checks for the custom website, vignette inventory, and
  package/docs sync, while excluding deployed site assets from package builds.
