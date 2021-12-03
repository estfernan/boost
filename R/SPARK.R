##
## R package boost by Esteban Fernández, Xi Jiang, Suhana Bedi, and Qiwei Li
## Copyright (C) 2021
##
## This file is part of the R package boost.
##
## The R package boost is free software: You can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by the
## Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package boost is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY without even the implied warranty of MERCHANTABILITY
## or FITNESS FOR A PARTICULAR PURPOSE.
##

##' Spatially Resolved Transcriptomic Analysis
##'
##' SPARK is an efficient method to identify spatially-variable genes.
##'
##' This method directly models count data generated from various spatial
##' resolved transcriptomic techniques through generalized spatial
##' linear models. It relies on penalized quasi-likelihood algorithm for
##' scalable computation and recently developed statistical formulas for
##' hypothesis testing, providing effective control of type I errors and
##' yielding high statistical power.
##'
##' See <https://xzhoulab.github.io/SPARK/> for more information.
##'
##' @param abs.expr A numeric vector \eqn{p} of length \eqn{n} that denotes the
##'   absolute gene expression levels. Each entry is an integer that denotes
##'   the gene count at spot \eqn{i}.
##' @param spots An \eqn{n}-by-\eqn{2} numeric matrix \eqn{T} to represent the
##'   geospatial profile, where each row indicates the spot location in
##'   the grid.
##' @param size.factor A numeric vector \eqn{s} of length \eqn{n} to compute
##'   the relative gene expression levels. Each entry denotes the size factor
##'   of sample \eqn{i} that captures all nuisance effects.
##' @param gene.name A character string that specifies the name of the gene
##'   passed. To be used when storing the results. The default value is `NULL`
##'   to keep the gene expression levels unnamed.
##'
##' @return `SPARK` returns an object of class "`SPARK`".
##'   The function [base::print()] i.e., [print.SPARK()], can be used to
##'   print a summary of the results.
##'
##' An object of class "`SPARK`" is a list containing the following components:
##'   \item{call}{the function call in which all of the specified arguments are specified by their full names.}
##'   \item{model}{the name of statistical model or technique.}
##'   \item{gene.name}{the name of gene evaluated.}
##'   \item{summary}{a summary table that contains the p-values for the different tests.}
##'   \item{measures}{the adjusted and combined \eqn{p}-values.}
##'   \item{time}{the execution time of the function.}
##'
##' @references Sun, S., Zhu, J. & Zhou, X. Statistical analysis of
##' spatial expression patterns for spatially resolved transcriptomic studies.
##' *Nat Methods* **17**, 193–200 (2020).
##' <https://doi.org/10.1038/s41592-019-0701-7>.
##'
##' @examples
##' ## Need to implement this example.
##'
##' @seealso
##' [get.size.factor()] for obtaining the size factors.
##'
##' @export
##' @keywords method
##'
##' @importFrom SPARK CreateSPARKObject spark.vc spark.test
##' @importFrom utils capture.output
##'
SPARK <- function(abs.expr, spots, size.factor, gene.name = NULL)
{
  if (!is.vector(abs.expr))
  {
    stop("value passed to 'abs.expr' must be a vector")
  }

  if (length(abs.expr) != dim(spots)[1])
  {
    stop("length of 'abs.expr' does not match sample size of 'spots'")
  }

  if (length(size.factor) != length(abs.expr))
  {
    stop("length of 'size.factor' does not match the length of 'abs.expr'")
  }

  loc_spark <- data.frame(spots)
  count_spark <- matrix(abs.expr, nrow = 1)
  count_spark <- as.data.frame(count_spark)

  colnames(count_spark) <- rownames(loc_spark)
  rownames(count_spark) <- gene.name
  stt <- Sys.time()

  spark <- CreateSPARKObject(counts = count_spark,
                             location = loc_spark,
                             percentage = 0,
                             min_total_counts = 0)

  ## total counts for each cell/spot
  spark@lib_size <- size.factor

  invisible(
    capture.output(
      {
        spark <- spark.vc(spark,
                          covariates = NULL,
                          lib_size = spark@lib_size,
                          num_core = 1,
                          verbose = FALSE)

        spark <- spark.test(spark,
                            check_positive = TRUE,
                            verbose = FALSE)
      }
    )
  )

  re <- spark@res_mtest
  endd <- Sys.time()

  rownames(re) <- NULL

  run_time <- as.numeric(difftime(endd, stt, units = "secs"))

  adjusted.pval <- re["adjusted_pvalue"]
  combined.pval <- re["combined_pvalue"]

  names(adjusted.pval) <- NULL
  names(combined.pval) <- NULL

  re[c("adjusted_pvalue", "combined_pvalue")] <- NULL

  structure(
    list(
      call      = match.call(),
      model     = "SPARK",
      gene.name = gene.name,
      summary   = re,
      measures  = c(
        adjusted.pval = adjusted.pval, combined.pval = combined.pval
      ),
      time      = run_time
    ),
    class = "SPARK"
  )
}

##' Print SPARK Results
##'
##' @param x An object of class `SPARK`.
##' @param ... Additional arguments passed to [base::print()].
##'
##' @export
##'
print.SPARK <- function(x, ...)
{
  print("this function has not been implemented yet")
}
