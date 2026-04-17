# Lrnr_cate_DR_selector Class

Cross-validated selector learner for CATE portfolios.

## Format

An R6 class that inherits from
[`Lrnr_hte`](https://larsvanderlaan.github.io/hte3/reference/Lrnr_hte.md).

## Details

This learner is primarily used internally by
[`cross_validate_cate`](https://larsvanderlaan.github.io/hte3/reference/cross_validate_cate.md)
to score and combine candidate CATE learners using the DR
pseudo-outcome. It remains exported for backward compatibility and
advanced workflows.
