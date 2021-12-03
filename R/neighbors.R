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

##' Get Neighborhood Information
##'
##' Obtain a long format of the adjacency matrix, specifically, for the
##' geospatial profile of the spatial transcriptomics (ST) data.
##'
##' Generally, for 10X and ST data, the number of neighbors are 6 and 4,
##' respectively. These values represent how the spots are arranged in a
##' lattice. For example, 6 and 4 indicate a hexagonal and square tiling,
##' respectively.
##'
##' The neighboring spots are obtained by two possible methods, where the
##' neighbors are obtained via a (1) distance threshold, which sets a cutoff
##' value, or the (2) k-nearest neighbors algorithm, which chooses the closest
##' neighbors in distance. Additionally, for method (1), a spot might not have
##' a full connectivity.
##'
##' See Jiang et al. (2021) for an introduction to the lattice grid system.
##'
##' @param spots An \eqn{n}-by-\eqn{2} numeric matrix \eqn{T} to represent the
##'   geospatial profile, where each row indicates the spot location in
##'   the grid.
##' @param n.neighbor An integer value that denotes the number of neighbors or
##'   connectivity. See "Details" for more information.
##' @param method An optional character string that specifies the method to
##'   obtain the neighboring spots. The default is "distance" to set
##'   a cutoff value.
##'
##' @return A numeric matrix that denotes the long format of the
##'   adjacency matrix. Each entry denotes the neighbor for spot \eqn{i}. If a
##'   spot does not have full connectivity, then the latest entry or entries
##'   are zero.
##'
##' @references Jiang, X., Li, Q., & Xiao, G. (2021). Bayesian Modeling of
##'   Spatial Transcriptomics Data via a Modified Ising Model.
##'   *arXiv preprint arXiv:2104.13957*.
##'
##' @examples
##' ## Need to implement the example.
##'
##' @seealso
##' [st.plot()] for plotting the expression levels and, hence, the lattice grid.
##'
##' @export
##' @keywords preprocessing
##'
##' @importFrom dbscan kNN
##'
get.neighbors <- function(spots, n.neighbor, method = c("distance", "KNN"))
{
  if (length(method) > 1)
  {
    method <- method[1]
  }

  n <- nrow(spots)

  P     <- matrix(0, nrow = n, ncol = n.neighbor)
  spots <- as.matrix(spots)

  if (method == "distance")
  {
    if (n.neighbor == 4)
    {
      spots <- round(spots)
      aa    <- sqrt(2)
    }
    else if (n.neighbor == 6)
    {
      aa <- sqrt(3)
    }
    else
    {
      aa <- 1.2
    }

    dist_matrix <- vectorized_pdist(spots, spots)
    min_dist    <- min(dist_matrix[dist_matrix > 0])

    dist_threshold <- min_dist*(aa - 1)*0.5 + min_dist

    for (i in 1:n)
    {
      k <- 1

      for (j in 1:n)
      {
        if (dist_matrix[i, j] > 0 & dist_matrix[i, j] < dist_threshold)
        {
          P[i, k] <- j
          k       <- k + 1
        }
      }
    }
  }
  else if (method == "KNN")
  {
    re_temp <- kNN(spots, k = n.neighbor)
    P       <- re_temp$id
  }
  else
  {
    stop("value passed to 'method' is not a valid option")
  }

  return(P)
}
