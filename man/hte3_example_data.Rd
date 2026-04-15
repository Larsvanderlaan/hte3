\name{hte3_example_data}
\alias{hte3_example_data}
\title{Simulate Tiny Example Data for hte3}
\usage{
hte3_example_data(n = 200, seed = 123)
}
\arguments{
\item{n}{Number of rows to simulate.}

\item{seed}{Optional random seed.}
}
\value{
A \code{data.table} with columns \code{W1}, \code{W2}, \code{W3}, \code{A},
\code{Y}, \code{Y_binary}, \code{tau}, \code{pi1}, \code{mu0}, and
\code{mu1}.
}
\description{
Convenience generator used in the package examples, tests, and vignettes.
}
