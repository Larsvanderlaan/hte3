# Extending hte3 with Custom Learners

This guide is for researchers who want to prototype a new `hte3` learner
and for contributors who want to add a new learner to the package
itself.

The key design idea is that `hte3` separates:

- data and nuisance estimates, stored in an `hte3_Task`
- meta-learner logic, implemented as `Lrnr_*` classes
- cross-validation and selection, handled by `cross_validate_*()` or the
  high-level wrappers

## Framework Overview

Most custom learners in `hte3` inherit from the internal base class
`hte3:::Lrnr_hte`. That class already handles:

- training a base `sl3` learner on a pseudo-outcome task
- prediction on new modifier data
- passing through observation weights
- basic treatment-type validation

Your learner supplies the parts that are specific to the estimand and
the learning algorithm.

The minimal ingredients are:

1.  `initialize()`: declare the learner parameters and pass the base
    learner to `Lrnr_hte`.
2.  `get_pseudo_data()`: compute the pseudo-outcome and, optionally,
    pseudo-weights from an `hte3_Task`.
3.  Treatment compatibility: set the private `.treatment_type` field so
    the learner rejects unsupported task types.
4.  Pseudo-outcome type and family: choose the right regression target
    for the downstream base learner.

## The Data Flow

At fit time, the typical flow is:

1.  Build an `hte3_Task` with modifiers `V`, confounders `W`, treatment
    `A`, outcome `Y`, and nuisance estimates.
2.  Inside `get_pseudo_data()`, read the nuisance quantities you need
    from the task.
3.  Construct a pseudo-outcome and optional pseudo-weights.
4.  Let `Lrnr_hte` build the derived `sl3_Task` and train the supplied
    `base_learner`.
5.  Predict on the modifier set `V` through the inherited prediction
    path.

When you are adding a learner to the package, try to keep the
pseudo-data step as simple and testable as possible. Small helper
functions are easier to test than one large learner method.

## Minimal Contract

Use this checklist when you create a new learner:

- The learner inherits from `hte3:::Lrnr_hte`.
- `initialize()` stores every tuning knob in `params`.
- `get_pseudo_data()` returns a list with `pseudo_outcome` and,
  optionally, `pseudo_weights`.
- The returned pseudo-outcome matches the declared
  `pseudo_outcome_type`/`pseudo_family`.
- The learner sets an appropriate `.treatment_type`.
- Predictions only depend on the modifier variables that define the
  target.

## Research Prototype Skeleton

The example below is intentionally minimal and is shown with
`eval = FALSE`. It is meant as a template for a research prototype, not
as a recommended final estimator.

``` r
library(R6)
library(hte3)
library(sl3)

Lrnr_cate_prototype <- R6Class(
  classname = "Lrnr_cate_prototype",
  inherit = hte3:::Lrnr_hte,
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(base_learner, treatment_level = NULL, control_level = NULL, ...) {
      params <- list(
        base_learner = base_learner,
        treatment_level = treatment_level,
        control_level = control_level,
        ...
      )

      super$initialize(
        params = params,
        base_learner = base_learner,
        pseudo_outcome_type = "continuous",
        pseudo_family = stats::gaussian(),
        ...
      )
    },
    get_pseudo_data = function(hte3_task, treatment_level = NULL, control_level = NULL, ...) {
      pseudo_outcome <- hte3:::compute_cate_dr_pseudo_outcome(
        hte3_task,
        treatment_level = treatment_level,
        control_level = control_level
      )

      list(pseudo_outcome = pseudo_outcome)
    }
  ),
  private = list(
    .treatment_type = c("binomial", "categorical"),
    .properties = c("cate", "prototype")
  )
)
```

That is the minimum viable pattern:

- parameters live in `self$params`
- the nuisance-aware pseudo-outcome is computed inside
  `get_pseudo_data()`
- the base class takes care of task construction and prediction

## Choosing the Right Family

Use the regression family that matches the pseudo-outcome you are
fitting:

- Gaussian pseudo-outcomes: CATE-style least-squares meta-learners such
  as DR, R, and EP.
- Quasibinomial pseudo-outcomes with weights: bounded ratio-style
  targets such as the current CRR EP learner.
- Binomial pseudo-outcomes: only when the pseudo-outcome itself is a
  genuine Bernoulli-style target.

If you are unsure, inspect the existing `Lrnr_cate_DR`, `Lrnr_cate_EP`,
and `Lrnr_crr_EP` implementations and match the family to the
pseudo-outcome scale, not to the original observed outcome scale.

## When to Expose a New Learner in the Package

Before adding a research prototype to the package API, lock down these
decisions:

- What estimand does it target?
- Which treatment types does it support?
- Does it require contrast levels?
- What nuisance quantities does it assume are available?
- How should users tune it?
- Does it need its own selector loss for outer cross-validation?

If any of those are still unstable, keep the learner as a research-only
prototype first.

## Promotion Checklist

When a learner is ready to become part of `hte3`, do all of the
following:

- Add the learner class and roxygen docs under `R/`.
- Add tests for fit/predict behavior, treatment validation, and any
  special numerical guards.
- Add at least one end-to-end wrapper or cross-validation test if the
  learner participates in a public workflow.
- Document the learner in the relevant workflow vignette and on the
  website.
- Decide whether it belongs in the high-level wrapper API or only in the
  advanced interface.
- Run the release checks described in the contribution guide.

## Related Guides

- The `advanced-sl3` vignette shows how low-level learner portfolios
  plug into
  [`cross_validate_cate()`](https://larsvanderlaan.github.io/hte3/reference/cross_validate_cate.md)
  and
  [`cross_validate_crr()`](https://larsvanderlaan.github.io/hte3/reference/cross_validate_crr.md).
- The `ep-learner` vignette shows how a more complex learner uses the
  `Sieve`-based basis construction inside `hte3`.
- The `contributing` vignette describes the repo workflow, site sync
  commands, and release checks for package contributors.
