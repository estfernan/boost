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

# //////////////////////////////////////////////////////////////////////////////
# Algorithms ----
# //////////////////////////////////////////////////////////////////////////////

##' Compute the Euclidean distance across two matrices
##'
##' @noRd
##' @keywords internal
##'
vectorized_pdist <- function(A,B)
{
  ##
  ## NOTE: The code below has not been thoroughly checked and reformatted for a
  ##         better understanding
  ##

  an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
  bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))
  m = nrow(A)
  n = nrow(B)
  tmp = matrix(rep(an, n), nrow = m)
  tmp = tmp +  matrix(rep(bn, m), nrow = m, byrow = TRUE)
  sqrt( tmp - 2 * tcrossprod(A,B) )
}

# //////////////////////////////////////////////////////////////////////////////
# Formatting functions ----
# //////////////////////////////////////////////////////////////////////////////

#'
#' Utility function for formatting p-values to use in the printing functions.
#'
#' @noRd
#' @keywords internal
#'
pval.format <- function(x, accuracy = 0.001, digits = 3)
{
  p.val <- round(x, digits = digits)

  if (p.val < accuracy)
  {
    p.val <- paste0("<", accuracy)
  }

  else if (p.val > 1 - accuracy)
  {
    p.val <- paste0(">", 1 - accuracy)
  }

  return(p.val)
}

# //////////////////////////////////////////////////////////////////////////////
# Graphics functions ----
# //////////////////////////////////////////////////////////////////////////////

##' Density Plot for MCMC Samples
##'
##' @noRd
##' @keywords internal plot
##'
##' @importFrom ggplot2 ggplot aes geom_histogram geom_density geom_vline labs
##'   coord_fixed theme_light theme element_blank element_rect element_text
##' @importFrom rlang .data
##' @importFrom stats quantile
##'
density.plot <- function(s, ind = NULL, title = "", subtitle = "")
{
  thin <- !is.null(ind)

  if (thin)
  {
    s <- s[ind]
  }

  ggplot(mapping = aes(x = s)) +
    geom_histogram(
      mapping = aes(y = .data$..density..),
      bins = 30,
      color = "#000000", fill = "#E69F00"
    ) +
    geom_density() +
    geom_vline(xintercept = mean(s), color = "#D55E00", linetype = 2) +
    geom_vline(xintercept = quantile(s, 0.025), color = "#D55E00") +
    geom_vline(xintercept = quantile(s, 0.975), color = "#D55E00") +
    labs(
      title    = title,
      subtitle = subtitle
    ) +
    coord_fixed() +
    theme_light() +
    theme(
      aspect.ratio = 1,
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
      plot.title   = element_text(face = 2)
    )
}

##' Trace Plot for MCMC Samples
##'
##' @noRd
##' @keywords internal plot
##'
##' @importFrom ggplot2 ggplot aes geom_line geom_hline labs scale_x_continuous
##'   coord_fixed theme_light theme element_blank element_rect element_text
##' @importFrom scales comma
##'
trace.plot <- function(s, iter, title = "", subtitle = "")
{
  ggplot(mapping = aes(x = iter, y = s)) +
    geom_line() +
    geom_hline(yintercept = mean(s), color = "#D55E00", linetype = 2) +
    labs(title = title, subtitle = subtitle) +
    scale_x_continuous(labels = comma) +
    coord_fixed() +
    theme_light() +
    theme(
      aspect.ratio = 1,
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
      plot.title   = element_text(face = 2)
    )
}

# //////////////////////////////////////////////////////////////////////////////
# Helper functions ----
# //////////////////////////////////////////////////////////////////////////////

#'
#' Utility function for BOOST-GP model to obtain the distance matrix and
#'   neighbor information for the geospatial profile
#'
#' @noRd
#' @keywords internal
#'
dist.GP <- function(loc, cutoff)
{
  ##
  ## NOTE: The code below has not been thoroughly checked and reformatted for a
  ##         better understanding
  ##

  dist <- as.matrix(dist(loc, method = "euclidean", diag = TRUE, upper = TRUE))
  rownames(dist) <- NULL
  colnames(dist) <- NULL

  n <- dim(dist)[1]
  neighbors <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    temp <- c(which(dist[i,] <= cutoff) - 1, -1)
    neighbors[i, 1:length(temp)] <- temp
  }
  neighbors <- neighbors[, which(colSums(neighbors) != 0)]
  return (list("dist" = dist, "nei" = neighbors))
}
