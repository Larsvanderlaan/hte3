\name{hte_task}
\alias{hte_task}
\title{Create an hte3 Task with a Production-Oriented Interface}
\usage{
hte_task(
  data,
  modifiers,
  confounders = modifiers,
  treatment,
  outcome,
  id = NULL,
  weights = NULL,
  treatment_type = c("default", "binomial", "categorical", "continuous"),
  propensity = NULL,
  outcome_regression = NULL,
  outcome_mean = NULL,
  propensity_learner = get_autoML(),
  outcome_learner = get_autoML(),
  mean_learner = NULL,
  cross_fit = TRUE,
  folds = 10,
  for_prediction = FALSE,
  ...
)
}
\arguments{
\item{data}{A data frame or \code{data.table} containing treatment, outcome, modifiers, and confounders.}

\item{modifiers}{Character vector of effect-modifier column names. These define the target modifier set \code{V}.}

\item{confounders}{Character vector of adjustment-variable column names used for nuisance estimation. These define the adjustment set \code{W}.}

\item{treatment}{Treatment column name.}

\item{outcome}{Outcome column name.}

\item{id}{Optional ID column name.}

\item{weights}{Optional weight column name.}

\item{treatment_type}{Treatment scale.}

\item{propensity}{Optional matrix of user-supplied propensity estimates.}

\item{outcome_regression}{Optional matrix of user-supplied outcome-regression estimates.}

\item{outcome_mean}{Optional vector of user-supplied marginal outcome estimates.}

\item{propensity_learner}{Learner used for propensity estimation when \code{propensity} is not supplied.}

\item{outcome_learner}{Learner used for outcome-regression estimation when \code{outcome_regression} is not supplied.}

\item{mean_learner}{Learner used for \code{m(W)} estimation when \code{outcome_mean} is not supplied.}

\item{cross_fit}{Whether to cross-fit nuisance learners.}

\item{folds}{Number of folds for nuisance estimation.}

\item{for_prediction}{Whether to construct a prediction-only task with no nuisance fitting.}

\item{...}{Additional arguments passed through to \code{make_hte3_Task_tx()}.}
}
\value{
An \code{hte3_Task} object.
}
\description{
\code{hte_task()} is the recommended high-level constructor for new users. It
wraps \code{make_hte3_Task_tx()} with validated, production-oriented argument
names.
}
\details{
When \code{modifiers = V} and \code{confounders = W}, the natural CATE target is
\code{E[Y(1) - Y(0) | V] = E[tau(W) | V]}. In the supported
binary/categorical-treatment setting, DR-, EP-, and the default two-stage
T-learner align with that modifier-conditional target. The current R-learner
does not generally target that object when \code{V} is a strict subset of
\code{W}; it warns at fit time.
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

task
}
}
