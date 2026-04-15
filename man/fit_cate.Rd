\name{fit_cate}
\alias{fit_cate}
\title{Fit a CATE Model}
\usage{
fit_cate(
  task,
  method = c("dr", "r", "t", "ep"),
  base_learner = get_autoML(),
  treatment_level = NULL,
  control_level = NULL,
  cross_validate = NULL,
  cv_control = NULL,
  sieve_num_basis = NULL,
  sieve_interaction_order = 3,
  screen_basis_with_lasso = FALSE
)
}
\arguments{
\item{task}{An \code{hte3_Task}.}

\item{method}{Either a single method or a character vector of methods drawn
from \code{"dr"}, \code{"r"}, \code{"t"}, and \code{"ep"}. A
vector creates a learner portfolio that is selected by cross-validation.}

\item{base_learner}{Base supervised learner used by the meta-learner.}

\item{treatment_level}{Optional treated level used for binary or categorical contrasts.}

\item{control_level}{Optional control level used for binary or categorical contrasts.}

\item{cross_validate}{Whether to cross-validate across candidate learners.
Defaults to \code{TRUE} for stacks and is automatically enabled for
learner portfolios.}

\item{cv_control}{Optional \code{sl3} cross-validation control list.}

\item{sieve_num_basis}{Optional sieve basis size for EP-learner fits.}

\item{sieve_interaction_order}{Interaction order for EP-learner sieve construction.}

\item{screen_basis_with_lasso}{Whether to screen EP sieve terms with lasso.}
}
\value{
An object of class \code{hte3_model}.
}
\description{
High-level wrapper for fitting conditional average treatment effect models.
The wrapper supports both single-method fits and cross-validated portfolios
across multiple learner families.
}
\details{
For continuous-treatment CATE tasks, the supported high-level method is
\code{"r"} only. That path uses the partially linear R-learner effect model
\code{A * tau(X)}.
}
