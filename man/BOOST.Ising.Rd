% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BOOST.Ising.R
\name{BOOST.Ising}
\alias{BOOST.Ising}
\title{BOOST-Ising Model for Dichotomised Expression Levels}
\usage{
BOOST.Ising(
  bin.expr,
  neighbor.info,
  gene.name = NULL,
  mean.omega0 = 1,
  sigma.omega0 = 2.5,
  mean.theta = 0,
  sigma.theta = 1,
  n.iter = 10000,
  burn.prop = 0.5
)
}
\arguments{
\item{bin.expr}{A numeric vector \eqn{p} of length \eqn{n} that denotes the
dichotomised gene expression levels. Each entry is one if the gene is
highly expressed at spot \eqn{i} and zero otherwise.}

\item{neighbor.info}{An \eqn{n}-by-\eqn{K} numeric matrix \eqn{A} that
denotes the long format of the adjacency matrix. Each entry denotes the
neighbor for spot \eqn{i}.}

\item{gene.name}{A character string that specifies the name of the gene
passed. To be used when storing the results. The default value is \code{NULL}
to keep the gene expression levels unnamed.}

\item{mean.omega0}{A numeric value that denotes the prior mean of the
normally-distributed first-order intensity parameter. The default is a
mean of one.}

\item{sigma.omega0}{A numeric value that denotes the prior standard deviation of
the normally-distributed first-order intensity parameter. The default is
a standard deviation of 2.5.}

\item{mean.theta}{A numeric value that denotes the prior mean of the
normally-distributed interaction parameter. The default is a mean of 0.}

\item{sigma.theta}{A numeric value that denotes the prior standard deviation of
the normally-distributed interaction parameter. The default is a
standard deviation of 1.}

\item{n.iter}{An integer value to specify the number of iterations for the
DMH algorithm. The default is 10,000 iterations.}

\item{burn.prop}{A numeric value to specify the proportion of iterations to
use as warm-up. The default is 0.50 to use half of the iterations
for warm-up.}
}
\value{
\code{BOOST.Ising} returns an object of class "\code{BOOST.Ising}".
The function \code{\link[base:print]{base::print()}} i.e., \code{\link[=print.BOOST.Ising]{print.BOOST.Ising()}}, can be used to
print a summary of the results.

An object of class "\code{BOOST.Ising}" is a list containing the following components:

\item{call}{the function call in which all of the specified arguments are specified by their full names.}
\item{model}{the name of statistical model or technique.}
\item{gene.name}{the name of gene evaluated.}
\item{summary}{a summary table that contains a summary of the estimated parameters.}
\item{measures}{the estimated Bayes factors and corresponding \eqn{p}-values.}
\item{time}{the execution time of the function.}
}
\description{
Fit the BOOST-Ising model to detect whether the gene is spatially
variable (SV). The fit is done within a double Metropolis-Hastings (DMH)
algorithm. Only one gene must be present and the expression levels must
be dichotomised.
}
\details{
The primary interest lies in the identification of SV genes via making
inferences on the interaction parameter between the low and
high-expression states. See Jiang et al. (2021) for more information on
the model fitting and posterior inference procedures.
}
\examples{
\dontrun{
library(boost)

## load sample dataset
data(mob)

## extract a sample gene, dichotomise expression levels, and get spatial network
g <- binarize.st(mob, "Apoe", cluster.method = "GMC")
A <- get.neighbors(mob.spots, 4, method = "distance")

## fit the model
res <- BOOST.Ising(g, A, gene.name = "Apoe", n.iter = 500)
print(res)
}
}
\references{
Jiang, X., Li, Q., & Xiao, G. (2021). Bayesian Modeling of
Spatial Transcriptomics Data via a Modified Ising Model.
\emph{arXiv preprint arXiv:2104.13957}.
}
\seealso{
\code{\link[=normalize.st]{normalize.st()}} for normalizing sequence count data;
\code{\link[=binarize.st]{binarize.st()}} for dichotomising relative expression levels;
\link{print.BOOST.Ising} for printing a summary of results to console.
}
\keyword{BOOST-Ising}
\keyword{method}
