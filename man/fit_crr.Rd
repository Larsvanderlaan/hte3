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
}
