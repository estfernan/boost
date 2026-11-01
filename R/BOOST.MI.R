##' BOOST-MI Model for Dichotomised Expression Levels
##'
##' `BOOST.MI` is the current name used for the modified-Ising implementation
##' originally exposed by [BOOST.Ising()]. It is a backwards-compatible alias:
##' both functions fit exactly the same model and return a `BOOST.Ising` object.
##' The default prior standard deviation of the interaction parameter is 1/3,
##' as used in STAr (Supplementary Notes S1). For real datasets, STAr reports
##' the larger of the attraction and repulsion Bayes factors (`BF.max`).
##'
##' @inheritParams BOOST.Ising
##' @return An object of class `BOOST.Ising`; see [BOOST.Ising()] for details.
##' @seealso [BOOST.Ising()]
##' @export
BOOST.MI <- function(
  bin.expr, neighbor.info,
  gene.name = NULL,
  mean.omega0 = 1, sigma.omega0 = 2.5,
  mean.theta = 0, sigma.theta = 1/3,
  n.iter = 1e4, burn.prop = 0.50
)
{
  res <- BOOST.Ising(
    bin.expr = bin.expr,
    neighbor.info = neighbor.info,
    gene.name = gene.name,
    mean.omega0 = mean.omega0,
    sigma.omega0 = sigma.omega0,
    mean.theta = mean.theta,
    sigma.theta = sigma.theta,
    n.iter = n.iter,
    burn.prop = burn.prop
  )
  res$call  <- match.call()
  res$model <- "BOOST-MI"
  res
}
