% Manually maintained to document the exported selector class.
\name{Lrnr_cate_DR_selector}
\alias{Lrnr_cate_DR_selector}
\title{Lrnr_cate_DR_selector Class}
\format{
An R6 class that inherits from \code{\link{Lrnr_hte}}.
}
\description{
Cross-validated selector learner for CATE portfolios.
}
\details{
This learner is primarily used internally by \code{\link{cross_validate_cate}}
to score and combine candidate CATE learners using the DR pseudo-outcome. It
remains exported for backward compatibility and advanced workflows.
}
