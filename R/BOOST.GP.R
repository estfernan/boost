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

##' Fit BOOST-GP Model for Gene Expression Levels
##'
##' Fit the BOOST-GP model to detect whether the gene is spatially
##' variable (SV). The fit is done within a Metropolis (SSVS) search variable
##' selection algorithm. Only one gene must be present, but no normalization
##' is necessary.
##'
##' The primary interest lies in the identification of SV genes via a
##' selection indicator. See Li et al. (2020) for more information on the
##' model fitting and posterior inference procedures.
##'
##' @param abs.expr A numeric vector \eqn{p} of length \eqn{n} that denotes the
##'   absolute gene expression levels. Each entry is an integer that denotes
##'   the gene count at spot \eqn{i}.
##' @param spots An \eqn{n}-by-\eqn{2} numeric matrix \eqn{T} to represent the
##'   geospatial profile, where each row indicates the spot location in
##'   the grid.
##' @param size.factor A numeric vector \eqn{s} of length \eqn{n} to compute
##'   the relative gene expression levels. Each entry denotes the size factor
##'   of sample \eqn{i} that captures all nuisance effects. The default
##'   is `NULL` to evaluate the absolute expression levels.
##' @param gene.name A character string that specifies the name of the gene
##'   passed. To be used when storing the results. The default value is `NULL`
##'   to keep the gene expression levels unnamed.
##' @param n.iter An integer value to specify the number of iterations for the
##'   DMH algorithm. The default is 10,000 iterations.
##' @param burn.prop A numeric value to specify the proportion of iterations to
##'   use as warm-up. The default is 0.50 to use half of the iterations
##'   for warm-up.
##' @param update.prop A numeric value to specify the proportion of samples to
##'   update in each iteration. The default is 0.2 to update one-fifth of the
##'   total samples.
##' @param init.b.sigma A numeric value to specify the initial value of the
##'   scale parameter in the inverse-gamma prior on the variance of the
##'   multivariate normal distribution prior for the log-expression levels.
##'   The default is `NULL` to set the value as twice the sample variance of
##'   the log-expression levels.
##' @param init.h A numeric value to specify the scaling of the variance for
##'   the normal prior set on each coefficient. The default is one to not scale
##'   the variance.
##'
##' @return `BOOST.GP` returns an object of class "`BOOST.GP`".
##'   The function [base::print()] i.e., [print.BOOST.GP()], can be used to
##'   print a summary of the results.
##'
##' An object of class "`BOOST.GP`" is a list containing the following components:
##'
##' \item{call}{the function call in which all of the specified arguments are specified by their full names.}
##' \item{model}{the name of statistical model or technique.}
##' \item{gene.name}{the name of gene evaluated.}
##' \item{summary}{a summary table that contains a summary of the estimated parameters.}
##' \item{measures}{the estimated Bayes factor and corresponding \eqn{p}-value}
##' \item{time}{the execution time of the function.}
##'
##' @references Li, Q., Zhang M., Xie Y., & Xiao, G. (2020). Bayesian Modeling
##' of Spatial Molecular Profiling Data via Gaussian Process.
##' *arXiv preprint arXiv:2012.03326*.
##'
##' @example inst/examples/ex_GP.R
##'
##' @seealso
##' [get.size.factor()] for estimating the size factor;
##' [print.BOOST.GP()] for printing a summary of results to console.
##'
##' @export
##' @keywords method BOOST-GP
##'
##' @importFrom stats pchisq
##'
BOOST.GP <- function(
  abs.expr, spots,
  size.factor = NULL,
  gene.name = NULL,
  n.iter = 1e4, burn.prop = 0.50,
  update.prop = 0.2,
  init.b.sigma = NULL, init.h = 1
)
{
  ##
  ## NOTE: The code below has not been thoroughly checked and reformatted for a
  ##         better understanding
  ##

  invalid.obj     <- !is.vector(abs.expr)
  absolute.levels <- is.null(size.factor)
  default.b.sigma <- is.null(init.b.sigma)

  if (invalid.obj)
  {
    stop("value passed to 'count' is not a vector for the expression levels")
  }

  n <- length(abs.expr)

  if (absolute.levels)
  {
    size.factor <- rep(1, n)
  }

  if (default.b.sigma)
  {
    log.expr     <- log(abs.expr)[abs.expr > 0]
    init.b.sigma <- var(log.expr) * 2
    init.b.sigma <- round(init.b.sigma, 3)
  }

  Y    <- matrix(abs.expr, ncol = 1)
  D    <- dist.GP(spots, 2)
  burn <- ceiling(n.iter * burn.prop)

  MCMC <- .boost_gp(
    Y, D$dist, D$nei,
    size.factor,
    n.iter, burn,
    init.b.sigma, init.h,
    update.prop
  )

  iter <- seq(burn + 1, n.iter)
  N    <- length(iter)

  log.lr             <- MCMC$logBF[iter, 1]
  log.lr[log.lr < 0] <- 0

  n        <- (n.iter - burn) * update.prop
  log.lr.f <- sort(log.lr, decreasing = TRUE)[1:n]

  BF    <- 19 * MCMC$gamma_ppi / (1 - MCMC$gamma_ppi + 10^(-10))
  p.val <- pchisq(BF * 2, df = 1, lower.tail = FALSE)

  l <- MCMC$l[iter, 1]
  l <- sum(l * log.lr) / sum(log.lr)

  log.lambda <- round(MCMC$logLambda[, 1], 3)

  est <- with(
    MCMC,
    {
      phi     <- phi[iter]
      lambda  <- exp(log.lambda)[iter]
      l       <- l[iter]

      data.frame(
        "Estimate" = c(
          mean(phi, na.rm = TRUE),
          mean(lambda, na.rm = TRUE),
          mean(l, na.rm = TRUE)
        ),
        "Std. Error" = c(
          sd(phi, na.rm = TRUE),
          sd(lambda, na.rm = TRUE),
          sd(l, na.rm = TRUE)
        ),
        "Lower 95% CI" = c(
          quantile(phi, 0.025, na.rm = TRUE),
          quantile(lambda, 0.025, na.rm = TRUE),
          quantile(l, 0.025, na.rm = TRUE)
        ),
        "Upper 95% CI" = c(
          quantile(phi, 0.975, na.rm = TRUE),
          quantile(lambda, 0.975, na.rm = TRUE),
          quantile(l, 0.975, na.rm = TRUE)
        ),
        row.names = c("Dispersion", "Normalized gene expression levels", "Kernel"),
        check.names = FALSE
      )
    }
  )

  ## add option to output this to an RData file
  # list(
  #   expr.levels   = abs.expr,
  #   spots         = spots,
  #   distance.info = D$dist,
  #   neighbor.info = D$nei,
  #   size.factor   = size.factor,
  #   SSVS          = list(b.sigma.init = init.b.sigma, h.init = init.h),
  #   control       = list(n.iter = n.iter, burn.prop = burn.prop, burn = burn, update.prop = update.prop, verbose = verbose),
  #   nSpots        = n,
  #   spots.loc     = apply(spots, 1, paste, collapse = " x "),
  #   estimate      = list(phi  = MCMC$phi_s[, 1], gamma = MCMC$gamma_s[, 1], lambda = exp(log.lambda), l = MCMC$l_s[, 1]),
  #   PPI           = list(H = MCMC$H_ppi[, 1], gamma = MCMC$gamma_ppi),
  #   accept        = MCMC$accept
  # )

  structure(
    list(
      call          = match.call(),
      model         = "BOOST-GP",
      gene.name     = gene.name,
      summary       = est,
      measures      = list(
        BF    = BF,
        p.val = p.val
      ),
      time          = MCMC$time
    ),
    class = "BOOST.GP"
  )
}

##' Print BOOST-GP Model Fitting Results
##'
##' S3 class method printing MCMC results from BOOST-GP fitted model. The
##' console displays the function call, model type, parameter estimates, and
##' Bayes Factor (BF) measurement for testing wether the gene is spatially
##' variable.
##'
##' @param x An object of class `BOOST.GP`.
##' @param ... Additional arguments passed to [base::print()].
##'
##' @return No return value, called for side effects i.e. printing to console.
##'
##' @examples
##' ## See example given in function for BOOST-GP
##'
##' @seealso
##' [BOOST.GP()] for BOOST-GP model fitting.
##'
##' @export
##' @keywords diagnostics BOOST-GP
##'
print.BOOST.GP <- function(x, ...)
{
  with(
    x,
    {
      cat("\nCall:\n")
      dput(call)

      cat("\nModel:", model, "\n")

      cat("\nParameters:\n")

      print(summary, digits = 2, ...)

      cat("\nBayes Factor in favor of a spatially variable gene: ")
      dput(measures$BF)

      cat("p-value: ")
      dput(measures$p.val)
    }
  )

  cat("\n")
  invisible(x)
}
