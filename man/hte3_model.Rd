\name{hte3_model}
\alias{predict.hte3_model}
\alias{summary.hte3_model}
\title{S3 Methods for Production hte3 Models}
\usage{
\method{predict}{hte3_model}(object, new_data = NULL, ...)

\method{summary}{hte3_model}(object, ...)
}
\arguments{
\item{object}{An object created by \code{fit_cate()} or \code{fit_crr()}.}

\item{new_data}{Optional prediction data. If omitted, predictions are returned on the training task.}

\item{...}{Unused additional arguments.}
}
\value{
\code{predict()} returns model predictions. \code{summary()} returns a compact
summary object with key training metadata.
}
\description{
These methods provide the default prediction and summary interface for the
production \code{hte3_model} wrapper class.
}
