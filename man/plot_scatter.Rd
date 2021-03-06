% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_scatter.R
\name{plot_scatter}
\alias{plot_scatter}
\title{Plot Recurrent CNV Before v.s. After WGD for each loci}
\usage{
plot_scatter(data, var, loci = NULL)
}
\arguments{
\item{data}{the CNV dataset (see \code{\link[graphics]{load_data}})}

\item{var}{The optimal value for all variables (see \code{\link[graphics]{Run_lasso}})}

\item{loci}{The particular loci to be emphasized}
}
\value{
ggplot2 graph object
}
\description{
Plot the values for the common events happenning before v.s. after WGD for each loci, by scatterplot
ggplot2 is needed for this function
}
\examples{

library(ggplot2)
p_scatter<-plot_scatter(wkdata,Lasso_res$var,loci=c(3606,864,4858,5261,5570))

}
