\name{fit_crr}
\alias{fit_crr}
\title{Fit a CRR Model}
\usage{
fit_crr(
  task,
  method = supported_crr_methods(),
  base_learner = get_autoML(),
  treatment_level = NULL,
  control_level = NULL,
  cross_validate = NULL,
  cv_control = NULL,
  sieve_num_basis = NULL,
  sieve_basis_grid = NULL,
  sieve_interaction_order = 3
)
}
\arguments{
\item{task}{An \code{hte3_Task}.}

\item{method}{Either a single method or a character vector of methods drawn
from \code{"ep"}, \code{"ipw"}, and \code{"t"}. A vector creates a
learner portfolio that is selected by cross-validation.}

\item{base_learner}{Base supervised learner used by the meta-learner.}

\item{treatment_level}{Optional treated level used for the contrast.}

\item{control_level}{Optional control level used for the contrast.}

\item{cross_validate}{Whether to cross-validate across candidate learners.
Defaults to \code{TRUE} for stacks and is automatically enabled for
learner portfolios.}

\item{cv_control}{Optional \code{sl3} cross-validation control list.}

\item{sieve_num_basis}{Optional sieve basis size for EP-learner fits.}

\item{sieve_basis_grid}{Optional EP basis-size grid used when
\code{method} includes \code{"ep"} and wrapper-level cross-validation is
enabled. If \code{NULL}, the wrapper uses the heuristic default grid
\code{c(d, 2d, 4d, 6d, 8d)}, where \code{d} is the number of modifiers.}

\item{sieve_interaction_order}{Interaction order for EP-learner sieve construction.}
}
\value{
An object of class \code{hte3_model}.
}
\description{
High-level wrapper for fitting conditional risk-ratio models with
non-negative outcomes. The wrapper supports both single-method fits and
cross-validated portfolios across multiple learner families.
}
\details{
CRR models are supported for binary or categorical treatment tasks and
require a non-negative outcome.

For EP-learner fits, \code{sieve_num_basis} controls a single fixed sieve size
when \code{cross_validate = FALSE}. When \code{cross_validate = TRUE},
\code{sieve_basis_grid} controls which EP basis sizes are compared by the
wrapper. After fitting, \code{summary()} reports the selected candidate and
the selected EP basis size when applicable.
}
\examples{
\dontrun{
library(sl3)

data <- hte3_example_data(n = 80, seed = 1)
task <- hte_task(
  data = data,
  modifiers = c("W1", "W2"),
  confounders = c("W1", "W2", "W3"),
  treatment = "A",
  outcome = "Y_binary",
  propensity_learner = Lrnr_mean$new(),
  outcome_learner = Lrnr_mean$new(),
  mean_learner = Lrnr_mean$new(),
  cross_fit = FALSE
)

fit <- fit_crr(
  task,
  method = "ep",
  base_learner = Lrnr_mean$new(),
  cross_validate = TRUE,
  cv_control = list(V = 2),
  sieve_basis_grid = c(3, 6, 9)
)

summary(fit)
}
}
