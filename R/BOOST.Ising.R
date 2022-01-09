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

##' BOOST-Ising Model for Dichotomised Expression Levels
##'
##' Fit the BOOST-Ising model to detect whether the gene is spatially
##' variable (SV). The fit is done within a double Metropolis-Hastings (DMH)
##' algorithm. Only one gene must be present and the expression levels must
##' be dichotomised.
##'
##' The primary interest lies in the identification of SV genes via making
##' inferences on the interaction parameter between the low and
##' high-expression states. See Jiang et al. (2021) for more information on
##' the model fitting and posterior inference procedures.
##'
##' @param bin.expr A numeric vector \eqn{p} of length \eqn{n} that denotes the
##'   dichotomised gene expression levels. Each entry is one if the gene is
##'   highly expressed at spot \eqn{i} and zero otherwise.
##' @param neighbor.info An \eqn{n}-by-\eqn{K} numeric matrix \eqn{A} that
##'   denotes the long format of the adjacency matrix. Each entry denotes the
##'   neighbor for spot \eqn{i}.
##' @param gene.name A character string that specifies the name of the gene
##'   passed. To be used when storing the results. The default value is `NULL`
##'   to keep the gene expression levels unnamed.
##' @param mean.omega0 A numeric value that denotes the prior mean of the
##'   normally-distributed first-order intensity parameter. The default is a
##'   mean of one.
##' @param sigma.omega0 A numeric value that denotes the prior standard deviation of
##'   the normally-distributed first-order intensity parameter. The default is
##'   a standard deviation of 2.5.
##' @param mean.theta A numeric value that denotes the prior mean of the
##'   normally-distributed interaction parameter. The default is a mean of 0.
##' @param sigma.theta A numeric value that denotes the prior standard deviation of
##'   the normally-distributed interaction parameter. The default is a
##'   standard deviation of 1.
##' @param n.iter An integer value to specify the number of iterations for the
##'   DMH algorithm. The default is 10,000 iterations.
##' @param burn.prop A numeric value to specify the proportion of iterations to
##'   use as warm-up. The default is 0.50 to use half of the iterations
##'   for warm-up.
##'
##' @return `BOOST.Ising` returns an object of class "`BOOST.Ising`".
##'   The function [base::print()] i.e., [print.BOOST.Ising()], can be used to
##'   print a summary of the results.
##'
##' An object of class "`BOOST.Ising`" is a list containing the following components:
##'
##' \item{call}{the function call in which all of the specified arguments are specified by their full names.}
##' \item{model}{the name of statistical model or technique.}
##' \item{gene.name}{the name of gene evaluated.}
##' \item{summary}{a summary table that contains a summary of the estimated parameters.}
##' \item{measures}{the estimated Bayes factors and corresponding \eqn{p}-values.}
##' \item{time}{the execution time of the function.}
##'
##' @references Jiang, X., Li, Q., & Xiao, G. (2021). Bayesian Modeling of
##' Spatial Transcriptomics Data via a Modified Ising Model.
##' *arXiv preprint arXiv:2104.13957*.
##'
##' @example inst/examples/ex_Ising.R
##'
##' @seealso
##' [normalize.st()] for normalizing sequence count data;
##' [binarize.st()] for dichotomising relative expression levels;
##' [print.BOOST.Ising] for printing a summary of results to console.
##'
##' @export
##' @keywords method BOOST-Ising
##'
##' @importFrom stats rnorm
##'
BOOST.Ising <- function(
  bin.expr, neighbor.info,
  gene.name = NULL,
  mean.omega0 = 1, sigma.omega0 = 2.5,
  mean.theta = 0, sigma.theta = 1,
  n.iter = 1e4, burn.prop = 0.50
)
{
  invalid.obj  <- !is.vector(bin.expr)
  dichotomized <- max(bin.expr) <= 1

  if (invalid.obj)
  {
    stop("value passed to 'count' is not a vector for the expression levels")
  }

  if (!dichotomized)
  {
    stop("values passed to 'count' must be the dichotomized expression levels")
  }

  burn        <- ceiling(n.iter * burn.prop)
  M           <- 3
  tau         <- 0.25
  omega0.init <- rnorm(1, 1, sigma.omega0)
  theta.init  <- rnorm(1, 0, sigma.theta)

  MCMC <- .BOOST_Ising_MCMC_cpp(
    bin.expr, neighbor.info,
    mean.omega0, sigma.omega0,
    mean.theta, sigma.theta,
    n.iter, burn,
    M, tau,
    omega0.init, theta.init
  )

  iter <- seq(burn, n.iter + 1)
  N    <- length(iter)

  pval.neg <- sum(MCMC$theta_s[iter] > 0) / N
  pval.pos <- sum(MCMC$theta_s[iter] < 0) / N

  est <- with(
    MCMC,
    {
      omega0 <- omega0_s[iter]
      theta  <- theta_s[iter]

      data.frame(
        "Prior Mean" = c(
          mean.omega0, mean.theta
        ),
        "Prior S.D." = c(
          sigma.omega0, sigma.theta
        ),
        "Estimate" = c(
          mean(omega0), mean(theta)
        ),
        "Std. Error" = c(
          sd(omega0), sd(theta)
        ),
        "Lower 95% CI" = c(
          quantile(omega0, 0.025), quantile(theta, 0.025)
        ),
        "Upper 95% CI" = c(
          quantile(omega0, 0.975), quantile(theta, 0.975)
        ),
        row.names = c("First-Order Intensity", "Interaction"),
        check.names = FALSE
      )
    }
  )

  ## add option to output this to an RData file
  # list(
  #   expr.levels   = bin.expr,
  #   neighbor.info = neighbor.info,
  #   base.prior    = list(mean = mean.omega0, sd = sigma.omega0),
  #   coef.prior    = list(mean = mean.theta,  sd = sigma.theta),
  #   DMH           = list(M = M, tau = tau, omega0.init = omega0.init, theta.init = theta.init),
  #   control       = list(n.iter = n.iter, burn.prop = burn.prop, burn = burn, verbose = verbose),
  #   nSpots        = length(bin.expr),
  #   estimate      = list(omega0 = MCMC$omega0_s, theta  = MCMC$theta_s),
  #   accept        = MCMC$accept
  # )

  structure(
    list(
      call          = match.call(),
      model         = "BOOST-Ising",
      gene.name     = gene.name,
      summary       = est,
      measures      = list(
        BF.neg   = (1 - pval.neg) / (pval.neg + 1e-10),  # repulsion
        BF.pos   = (1 - pval.pos) / (pval.pos + 1e-10),  # attractions
        pval.neg = pval.neg, pval.pos = pval.pos
      ),
      time          = MCMC$time
    ),
    class = "BOOST.Ising"
  )
}

##' Print BOOST-Ising Model Fitting Results
##'
##' S3 class method printing MCMC results from BOOST-Ising fitted model. The
##' console displays the function call, model type, parameter estimates, and
##' Bayes Factor (BF) measurement for testing an attraction pattern in the gene.
##'
##' @param x An object of class `BOOST.Ising`.
##' @param ... Additional arguments passed to [base::print()].
##'
##' @return No return value, called for side effects i.e. printing to console.
##'
##' @examples
##' ## Need to implement the example for the procedure.
##'
##' @seealso
##' [BOOST.Ising()] for BOOST-Ising model fitting.
##'
##' @export
##' @keywords diagnostics BOOST-Ising
##'
print.BOOST.Ising <- function(x, ...)
{
  with(
    x,
    {
      cat("\nCall:\n")
      dput(call)

      cat("\nModel:", model, "\n")

      cat("\nParameters:\n")

      print(summary, digits = 2, ...)

      cat("\nBayes Factor in favor of an attraction pattern: ")
      dput(measures$BF.neg)

      cat("p-value: ")
      cat(pval.format(measures$pval.neg))
    }
  )

  cat("\n")
  invisible(x)
}
