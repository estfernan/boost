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

### Incorporate into print(): "Key: 0=Low-expressed, 1=Highly-expressed"

##' Binary Spatial Extraction of Genes
##'
##' BinSpect is a fast computational method that identifies spatially variable
##' genes using the binarized expression levels.
##'
##' This method uses a graph-based approach to compute the probability of
##' encountering two physically neighboring cells being both
##' highly-expressed (=1). It is based on the provided neighborhood information,
##' which can be interpreted as a graph structure. Based on the dichotomised
##' gene expression levels, an edge in the "graph" connecting two cells is
##' classified as 1-0, 1-1, 0-1, or 0-0. A contingency table that counts
##' all the classified edges for the entire neighborhood network is created.
##' A \eqn{p}-value is then reported based on the hypergeometric
##' distribution test.
##'
##' See <http://spatialgiotto.rc.fas.harvard.edu/giotto_spatial_genes.html>
##' for more information.
##'
##' @param bin.expr A numeric vector \eqn{p} of length \eqn{n} that denotes the
##'   dichotomised gene expression levels. Each entry is one if the gene is
##'   highly expressed at spot \eqn{i} and zero otherwise.
##' @param neighbor.info An \eqn{n}-by-\eqn{K} numeric matrix \eqn{A} that
##'   denotes the long format of the adjacency matrix. Each entry denotes the
##'   neighbor for spot \eqn{i}.
##' @param do.fisher.test A logical value that indicates if the Fisher test
##'   should be used. The default value is `FALSE` to use a normal
##'   approximation test.
##' @param gene.name A character string that specifies the name of the gene
##'   passed. To be used when storing the results. The default value is `NULL`
##'   to keep the gene expression levels unnamed.
##'
##' @return `binSpect` returns an object of class "`binSpect`".
##'   The function [base::print()] i.e., [print.binSpect()], can be used to
##'   print a summary of the results.
##'
##' An object of class "`binSpect`" is a list containing the following components:
##'
##' \item{call}{the function call in which all of the specified arguments are specified by their full names.}
##' \item{model}{the name of statistical model or technique.}
##' \item{gene.name}{the name of gene evaluated.}
##' \item{summary}{a summary table that contains the contingency table of low-expressed and highly-expressed interactions.}
##' \item{measures}{the estimated odds ratio and corresponding \eqn{p}-value.}
##' \item{time}{the execution time of the function.}
##'
##' @references Dries, R., Zhu, Q., Dong, R. et al. Giotto: a toolbox for
##' integrative analysis and visualization of spatial expression data.
##' *Genome Biol* **22**, 78 (2021).
##' <https://doi.org/10.1186/s13059-021-02286-2>
##'
##' @example inst/examples/ex_BinSpect.R
##'
##' @seealso
##' [binarize.st()] for dichotomising the expression levels;
##' [get.neighbors()] for getting the neighborhood information.
##'
##' @export
##' @keywords method
##'
##' @importFrom stats fisher.test pnorm
##'
binSpect <- function(bin.expr, neighbor.info, do.fisher.test = FALSE, gene.name = NULL)
{
  ##
  ## NOTE: The code below has not been thoroughly checked and reformatted for
  ##         a better understanding
  ##

  if (!is.vector(bin.expr))
  {
    stop("value passed to 'bin.expr' must be a vector")
  }

  if (length(sort(unique(bin.expr))) != 2)
  {
    stop("vector passed to 'bin.expr' must be binary")
  }
  else if (sum(sort(unique(bin.expr)) == 0:1) < 2)
  {
    stop("vector passed to 'bin.expr' must only containe 0's and 1's")
  }

  if (length(bin.expr) != dim(neighbor.info)[1])
  {
    stop("length of 'bin.expr' does not match sample size of 'neighbor.info'")
  }

  st_time <- Sys.time()
  contingency_table <- matrix(0, ncol = 2, nrow = 2)

  sample_num <- nrow(neighbor.info)
  n_neighbor <- ncol(neighbor.info)

  for (i in 1:sample_num)
  {
    for (j in 1:n_neighbor)
    {
      if (neighbor.info[i, j] != 0)
      {
        nei_temp <- neighbor.info[i, j]

        row <- bin.expr[i] + 1
        col <- bin.expr[nei_temp] + 1

        contingency_table[row, col] <- contingency_table[row, col] + 1
      }
    }
  }

  colnames(contingency_table) <- c("0", "1")
  rownames(contingency_table) <- c("0", "1")

  if (do.fisher.test)
  {
    re <- fisher.test(contingency_table)
    p_value <- re$p.value
    estimate <- re$estimate
    names(estimate) <- NULL
  }
  else
  {
    a <- contingency_table[1, 1]
    b <- contingency_table[2, 1]
    c <- contingency_table[1, 2]
    d <- contingency_table[2, 2]
    estimate <- (a*d) / (b*c)

    p_value <- 2 * (
      1 - pnorm(
        log(ifelse(estimate > 1, estimate, 1/estimate)) / sqrt(1/a + 1/b + 1/c + 1/d)
      )
    )
  }

  end_time <- Sys.time()
  run_time <- as.numeric(difftime(end_time, st_time, units = "secs"))

  structure(
    list(
      call      = match.call(),
      model     = "BinSpect",
      gene.name = gene.name,
      summary   = as.data.frame(contingency_table),
      measures  = list(OR = estimate, p.val = p_value),
      time      = run_time
    ),
    class = "binSpect"
  )
}

##' Print BinSpect Results
##'
##' @param x An object of class `binSpect`.
##' @param ... Additional arguments passed to [base::print()].
##'
##' @export
##'
print.binSpect <- function(x, ...)
{
  with(
    x,
    {
      cat("\nCall:\n")
      dput(call)

      cat("\nModel:", model, "\n")

      cat("\nContingency Table for Classified Edges:\n")

      print(summary, digits = 2, ...)

      cat("\nOdds ratio in favor of a spatially-variable pattern: ")
      dput(round(measures$OR, 2))

      cat("p-value in favor of a spatially-variable pattern: ")
      cat(pval.format(measures$p.val), "\n")
    }
  )

  cat("\n")
  invisible(x)
}
