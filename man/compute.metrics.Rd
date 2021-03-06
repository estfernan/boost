% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metrics.R
\name{compute.metrics}
\alias{compute.metrics}
\title{Performance of Gene Identification}
\usage{
compute.metrics(
  predictor,
  truth,
  predictor.type = c("BF", "p-value"),
  threshold = NULL
)
}
\arguments{
\item{predictor}{A numeric vector of length \eqn{n} that denotes the
\eqn{p}-values or Bayes factors (BFs).}

\item{truth}{A logical vector of length \eqn{n} that represents the ground
truth corresponding to the predictor.}

\item{predictor.type}{A character string that specifies whether
\eqn{p}-values of Bayes factors (BFs) were provided.}

\item{threshold}{A numeric value that specifies the cutoff for
defining SV genes.}
}
\value{
A list object that contains six performance metrics
(Sensitivity, Specificity, F1_score, FDR, AUC, and MCC).
}
\description{
Calculate the performance of spatially variable (SV) gene identification
on simulated data.
}
\references{
Li, X., Wang, X., & Xiao, G. (2019). A comparative study of rank
aggregation methods for partial and top ranked lists in genomic
applications. \emph{Briefings in bioinformatics}, \emph{20}(1), 178–189.
\url{https://doi.org/10.1093/bib/bbx101}.

Robin, X., Turck, N., Hainard, A. et al. pROC: an open-source package for
R and S+ to analyze and compare ROC curves.
\emph{BMC Bioinformatics} \strong{12}, 77 (2011).
\url{https://doi.org/10.1186/1471-2105-12-77}.
}
\keyword{metrics}
