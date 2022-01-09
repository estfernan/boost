##
## R package boost by Esteban Fernández, Xi Jiang, Suhana Bedi, and Qiwei Li
## Copyright (C) 2021
##
## This file is part of the R package boost.
##
## The R package boost is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package boost is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##

##' Statistical Test to Identify Genes with Spatial Patterns
##'
##' SpatialDE is a method to identify and characterize spatially variable (SV)
##' genes via a statistical test.
##'
##' This method builds on Gaussian process regression, a class of models.
##' For a given gene, the expression variability is decomposed into spatial and
##' non-spatial components, using two random effect terms: a spatial variance
##' term that parametrizes gene expression covariance by pairwise distances of
##' samples, and a noise term that models non-spatial variability. The ratio
##' of the variance explained by these two components quantifies the
##' Fraction of Spatial Variance (FSV). Significant SV genes can be identified
##' by comparing this full model to a model without the
##' spatial variance component.
##'
##' See <https://bioconductor.org/packages/release/bioc/html/spatialDE.html>
##' for more information.
##'
##' @param norm.expr A numeric vector \eqn{p} of length \eqn{n} that denotes the
##'   relative gene expression levels. Each entry is a numeric value that
##'   denotes the normalized gene count at spot \eqn{i}.
##' @param spots An \eqn{n}-by-\eqn{2} numeric matrix \eqn{T} to represent the
##'   geospatial profile, where each row indicates the spot location in
##'   the grid.
##' @param gene.name A character string that specifies the name of the gene
##'   passed. To be used when storing the results. The default value is `NULL`
##'   to keep the gene expression levels unnamed.
##'
##' @return `SpatialDE` returns an object of class "`SpatialDE`".
##'   The function [base::print()] i.e., [print.SpatialDE()], can be used to
##'   print a summary of the results.
##'
##' An object of class "`SpatialDE`" is a list containing the following components:
##'
##' \item{call}{the function call in which all of the specified arguments are specified by their full names.}
##' \item{model}{the name of statistical model or technique.}
##' \item{gene.name}{the name of gene evaluated.}
##' \item{summary}{a summary table that contains a summary of the model.}
##' \item{measures}{the estimated \eqn{p}-value.}
##' \item{time}{the execution time of the function.}
##'
##' @references
##'
##' Corso D, Malfait M, Moses L (2021). *spatialDE: R wrapper for SpatialDE*.
##' doi: [10.18129/B9.bioc.spatialDE](https://doi.org/10.18129/B9.bioc.spatialDE),
##' R package version 1.0.0, <http://www.bioconductor.org/packages/spatialDE>.
##'
##' Svensson V, Teichmann SA, Stegle O (2018).
##' "SpatialDE: identification of spatially variable genes."
##' *Nature Methods*, **15**(5), 343–346. ISSN 1548-7105,
##' doi: [10.1038/nmeth.4636](https://doi.org/10.1038/nmeth.4636),
##' <https://www.nature.com/articles/nmeth.4636>.
##'
##' @example inst/examples/ex_SpatialDE.R
##'
##' @seealso
##' [normalize.st()] for normalizing the expression levels;
##' [print.SpatialDE()] for printing a summary of results to console.
##'
##' @export
##' @keywords method
##'
##' @importFrom spatialDE run
##'
SpatialDE <- function(norm.expr, spots, gene.name = NULL)
{
  if (!is.vector(norm.expr))
  {
    stop("value passed to 'norm.expr' must be a vector")
  }

  if (length(norm.expr) != dim(spots)[1])
  {
    stop("length of 'norm.expr' does not match sample size of 'spots'")
  }

  spots <- as.data.frame(spots)
  x     <- matrix(norm.expr, nrow = 1)

  colnames(spots) <- c("x", "y")
  rownames(x)     <- gene.name

  res <- spatialDE::run(x, spots, verbose = FALSE)

  rownames(res) <- ""

  structure(
    list(
      call      = match.call(),
      model     = "SpatialDE",
      gene.name = gene.name,
      summary   = res[c("g", "n", "FSV", "l", "BIC")],
      measures  = list(p.value = res$pval),
      time      = res$time
    ),
    class = "SpatialDE"
  )
}

##' Print SpatialDE Results
##'
##' @param x An object of class `SpatialDE`.
##' @param ... Additional arguments passed to [base::print()].
##'
##' @export
##'
print.SpatialDE <- function(x, ...)
{
  with(
    x,
    {
      cat("\nCall:\n")
      dput(call)

      cat("\nModel:", model, "\n")

      cat("\nSummary:\n")

      print(summary, digits = 2, ...)

      cat("\np-value in favor of a spatially-variable pattern: ")
      cat(pval.format(measures$p.value), "\n")
    }
  )

  cat("\n")
  invisible(x)
}
