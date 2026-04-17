# Changelog

## hte3 0.1.1

- Hardened the default
  [`get_autoML()`](https://larsvanderlaan.github.io/hte3/reference/get_autoML.md)
  stack so the high-level API uses the available safe subset of learners
  and warns once when optional runtime packages are missing.
- Finished the two-stage T-learner update for reduced-modifier settings
  so first-stage regressions use the confounder set `W` and the default
  second stage projects onto the modifier set `V`.
- Tightened task validation for folds, weights, treatment values, and
  user-supplied nuisance estimates.
- Added release checks for the custom website, vignette inventory, and
  package/docs sync, while excluding deployed site assets from package
  builds.
