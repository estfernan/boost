##' Run SPARK-X
##'
##' Wrapper for [SPARK::sparkx()] (Zhu et al., 2021). Following STAr
##' (Supplementary Notes S1), the model is fitted to the raw counts with
##' `option = "mixture"`, and the adjusted combined p-values are used to
##' select SVGs.
##'
##' @param count An \eqn{n}-by-\eqn{p} count matrix with spots in rows and
##'   genes in columns.
##' @param spots An \eqn{n}-by-2 matrix of spatial coordinates.
##' @param X Optional \eqn{n}-by-\eqn{q} matrix of spot-level covariates.
##' @param numCores Number of CPU cores.
##' @param option Kernel option passed to SPARK-X (default `"mixture"`).
##' @param verbose Logical; print SPARK-X progress messages.
##'
##' @return A data frame with one row per gene and columns `gene`, `method`,
##'   `score`, `score_type`, `pvalue` (combined p-value), and `padj` (adjusted
##'   combined p-value reported by SPARK-X). The native SPARK-X output is
##'   attached as attribute `"native"`.
##'
##' @references Zhu, J., Sun, S., & Zhou, X. (2021). SPARK-X: non-parametric
##'   modeling enables scalable and robust detection of spatial expression
##'   patterns for large spatial transcriptomic studies. *Genome Biology*,
##'   22, 184.
##'
##' @export
##'
SPARK.X <- function(count, spots, X = NULL, numCores = 1,
                    option = c("mixture", "single"), verbose = FALSE)
{
  .check_svg_input(count, spots)
  option <- match.arg(option)

  spot.id <- .spot_names(count, spots)

  ## SPARK-X expects a genes x spots count matrix and a spots x d location matrix
  counts_gs <- Matrix::Matrix(t(as.matrix(count)), sparse = TRUE)
  dimnames(counts_gs) <- list(colnames(count), spot.id)
  loc <- as.matrix(spots)[, 1:2, drop = FALSE]
  rownames(loc) <- spot.id

  fit <- SPARK::sparkx(count_in = counts_gs, locus_in = loc, X_in = X,
                       numCores = numCores, option = option, verbose = verbose)

  res <- as.data.frame(fit$res_mtest)
  out <- .svg_table(rownames(res), "SPARK-X", pvalue = res$combinedPval,
                    padj = res$adjustedPval)
  attr(out, "native") <- fit
  out
}
