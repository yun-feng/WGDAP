% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feval.R
\name{feval}
\alias{feval}
\title{Evaluate Loss}
\usage{
feval(data, par, var)
}
\arguments{
\item{data}{the CNV dataset (see \code{\link[graphics]{load_data}})}

\item{par}{The parameter used for fixed point method. (see \code{\link[graphics]{set_par}})}

\item{var}{The current value for all variables (see \code{\link[graphics]{init_var}})}
}
\value{
Value for the current loss
}
\description{
Evaluate current loss for fixed-point method based on current variable values.
}
\examples{

init_loss<-feval(wkdata,CNV_par,common_var)

}
