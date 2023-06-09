% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_ssf.R
\name{predict_ssf}
\alias{predict_ssf}
\title{Predict SSF for new data}
\usage{
predict_ssf(mod, new_data)
}
\arguments{
\item{mod}{Fitted model returned by \code{\link{hmmSSF}}}

\item{new_data}{Data frame with covariate values used for prediction}
}
\value{
Matrix of predictions in each state
}
\description{
Predict SSF for new data
}
