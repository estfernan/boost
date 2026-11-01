##' Run Trendsceek
##'
##' Wrapper for the R package trendsceek (Edsgard et al., 2018), following
##' STAr (Supplementary Notes S1): spatial positions are converted to a point
##' pattern with `trendsceek::pos2pp()`, the log10-transformed expression is
##' set as marks with `trendsceek::set_marks()`, and the permutation tests are
##' run with `trendsceek::trendsceek_test()`. The BH-adjusted p-values of the
##' mark-variogram statistic, which performed best in the STAr simulations,
##' are used by default.
##'
##' @param count An \eqn{n}-by-\eqn{p} count matrix with spots in rows and
##'   genes in columns.
##' @param spots An \eqn{n}-by-2 matrix of spatial coordinates.
##' @param statistic Summary statistic used for SVG selection: `"markvario"`
##'   (default), `"markcorr"`, `"Emark"`, or `"Vmark"`.
##' @param nrand Number of permutations.
##' @param ncores Number of CPU cores.
##' @param log.fcn Transformation applied to the marks (default `log10`).
##' @param ... Additional arguments passed to `trendsceek::trendsceek_test()`.
##'
##' @details The trendsceek package calls functions such as `spatstat::ppp()`
##'   and `spatstat::markvario()` from the legacy, single-package version of
##'   spatstat (< 2.0). With newer spatstat versions, which split these
##'   functions into sub-packages, trendsceek fails; install a legacy spatstat
##'   version (e.g., `remotes::install_version("spatstat", "1.64-1")`) in a
##'   separate library to run Trendsceek.
##'
##' @return A data frame with one row per gene and columns `gene`, `method`,
##'   `score`, `score_type`, `pvalue` (minimum permutation p-value), and
##'   `padj`. The native trendsceek output is attached as attribute `"native"`.
##'
##' @references Edsgard, D., Johnsson, P., & Sandberg, R. (2018).
##'   Identification of spatial expression trends in single-cell gene
##'   expression data. *Nature Methods*, 15, 339-342.
##'
##' @export
##'
Trendsceek <- function(count, spots,
                       statistic = c("markvario", "markcorr", "Emark", "Vmark"),
                       nrand = 1000, ncores = 1, log.fcn = log10, ...)
{
  .check_svg_input(count, spots)
  .require_svg_package("trendsceek", "remotes::install_github('edsgard/trendsceek')")
  statistic <- match.arg(statistic)
  if (!requireNamespace("spatstat", quietly = TRUE) ||
      !("ppp" %in% getNamespaceExports("spatstat")))
    stop("trendsceek requires the legacy spatstat package (< 2.0); ",
         "see ?Trendsceek", call. = FALSE)

  spot.id <- .spot_names(count, spots)
  counts_gs <- t(as.matrix(count))
  dimnames(counts_gs) <- list(colnames(count), spot.id)

  pp  <- trendsceek::pos2pp(as.matrix(spots)[, 1:2, drop = FALSE])
  pp  <- trendsceek::set_marks(pp, counts_gs, log.fcn = log.fcn)
  fit <- trendsceek::trendsceek_test(pp, nrand = nrand, ncores = ncores, ...)

  st  <- as.data.frame(fit$supstats[[statistic]])
  out <- .svg_table(st$gene, "Trendsceek", score_type = statistic,
                    pvalue = st$min.pval, padj = st$p.bh)
  attr(out, "native") <- fit
  out
}
