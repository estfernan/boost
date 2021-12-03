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
