\name{fit_cate}
\alias{fit_cate}
\title{Fit a CATE Model}
\usage{
fit_cate(
  task,
  method = supported_cate_methods(),
  base_learner = get_autoML(),
  treatment_level = NULL,
  control_level = NULL,
  cross_validate = NULL,
  cv_control = NULL,
  sieve_num_basis = NULL,
  sieve_basis_grid = NULL,
  sieve_interaction_order = 3,
  screen_basis_with_lasso = FALSE,
  ep_targeting_style = "dr",
  ep_r_targeting_basis = "v_plus_propensity"
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

\item{sieve_basis_grid}{Optional EP basis-size grid used when
\code{method} includes \code{"ep"} and wrapper-level cross-validation is
enabled. If \code{NULL}, the wrapper uses the heuristic default grid
\code{c(d, 2d, 4d, 6d, 8d)}, where \code{d} is the number of modifiers.}

\item{sieve_interaction_order}{Interaction order for EP-learner sieve construction.}

\item{screen_basis_with_lasso}{Whether to screen EP sieve terms with lasso.}

\item{ep_targeting_style}{EP targeting variant or variants used when
\code{method} includes \code{"ep"}. Use \code{"dr"} for EP-DR,
\code{"r"} for EP-R, or \code{c("dr", "r")} to let wrapper-level
cross-validation compare both EP variants.}

\item{ep_r_targeting_basis}{First-stage basis construction used for
\code{ep_targeting_style = "r"}. The default \code{"v_plus_propensity"}
uses the modifier set plus the treated-arm propensity score; \code{"full_w"}
uses the full confounder set.}
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

If \code{modifiers = V} and \code{confounders = W} with \code{V} a strict
subset of \code{W}, the natural target is \code{E[Y(1) - Y(0) | V]}. In the
supported binary/categorical-treatment setting, DR-, EP-, and the default
two-stage T-learner target that modifier-conditional CATE surface. The
current R-learner does not generally target that object in the
reduced-modifier setting and warns at fit time.

For EP-learner fits, \code{sieve_num_basis} controls a single fixed sieve size
when \code{cross_validate = FALSE}. When \code{cross_validate = TRUE},
\code{sieve_basis_grid} controls which EP basis sizes are compared by the
wrapper. After fitting, \code{summary()} reports the selected candidate, the
selected EP targeting style, and the selected EP basis size when applicable.
Wrapper-level EP tuning continues to use the DR selector loss from
\code{cross_validate_cate()}, including for the EP-R variant.

The CATE EP family now has two variants. EP-DR
(\code{ep_targeting_style = "dr"}) is the standard plug-in DR-style EP
variant. EP-R (\code{ep_targeting_style = "r"}) uses an overlap-weighted
targeting objective, can be more stable under severe overlap violations, and
currently supports binary-treatment CATE tasks only. For EP-R, the default
first-stage basis is \code{ep_r_targeting_basis = "v_plus_propensity"}, which
uses the modifier set plus the treated-arm propensity score as a
dimension-reduced summary of the adjustment set. To compare EP-DR and EP-R
directly at the wrapper level, set
\code{ep_targeting_style = c("dr", "r")}.
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
  outcome = "Y",
  propensity_learner = Lrnr_mean$new(),
  outcome_learner = Lrnr_mean$new(),
  mean_learner = Lrnr_mean$new(),
  cross_fit = FALSE
)

fit <- fit_cate(
  task,
  method = "ep",
  base_learner = Lrnr_mean$new(),
  cross_validate = TRUE,
  cv_control = list(V = 2),
  ep_targeting_style = c("dr", "r"),
  sieve_basis_grid = c(2, 4, 6, 8)
)

summary(fit)
}
}
