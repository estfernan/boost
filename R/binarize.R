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

##' Dichotomise Gene Expression Levels
##'
##' Cluster the gene expression levels into highly-expressed and
##' low-expressed groups.
##'
##' After some filtering steps, a clustering method is applied to de-noise the
##' expression levels of the given gene by dichotomising the spots into
##' low-expressed (=0) and highly-expressed (=1) groups. This preprocessing
##' step outputs the suitable data type required for some of the model fitting
##' procedures, as well as making it more robust to over-dispersion and
##' zero-inflation.
##'
##' For the gene at a given spot, the expression level was modified to a binary
##' format. Instead of a numeric value, a logical value is used to indicate the
##' group affiliation. This method provides three clustering methods: (1) rank,
##' a quantile-based approach which simply applies a cutoff at the specified
##' percentage, (2) Gaussian mixture clustering (GMC), which is a model-based approach
##' that fits a two-component Gaussian mixture model (GMM) with unequal
##' variances, and (3) k-means (k-means), a distance-based approach which is implicitly
##' based on the pairwise distances of the expression levels.
##'
##' See Jiang et al. (2021) for more information on the filtering steps and
##' last two clustering methods.
##'
##' @param count An \eqn{n}-by-\eqn{p} numeric matrix \eqn{Y} that denotes
##'   the relative gene expression count table. Each entry is the relative
##'   count for gene \eqn{j} collected at spot \eqn{i}.
##' @param gene.name A character string that specifies the gene in the
##'   expression count table to dichotomise.
##' @param cluster.method  An optional character string to specify the
##'   clustering technique. The default is "rank" for
##'   rank-based clustering.
##' @param percentage.rank An optional numeric value to specify the cutoff for
##'   the rank-based clustering method. The default is \eqn{0.30} to set the
##'   cell counts above the 70% quantile as highly-expressed.
##'
##' @return A binary and numeric vector to represent dichotomisation results.
##'   See "Details" for more information on how to interpret the entries.
##'
##' @references Jiang, X., Li, Q., & Xiao, G. (2021). Bayesian Modeling of
##'   Spatial Transcriptomics Data via a Modified Ising Model.
##'   *arXiv preprint arXiv:2104.13957*.
##'
##' @seealso
##' [st.plot()] for plotting the dichotomised expression levels.
##'
##' @export
##' @keywords preprocessing
##'
##' @importFrom mclust Mclust mclustBIC
##' @importFrom stats IQR kmeans quantile
##'
binarize.st <- function(
  count, gene.name,
  cluster.method = c("rank", "GMC", "k-means"),
  percentage.rank = 0.30
)
{
  ##
  ## NOTE: The code below has not been thoroughly checked and reformatted for
  ##         a better understanding
  ##

  if (!(gene.name %in% colnames(count)))
  {
    stop("value passed to 'gene.name' is not in the expression count table")
  }

  if (length(cluster.method) > 1)
  {
    cluster.method <- cluster.method[1]
  }

  count <- count[, gene.name]
  count <- as.matrix(count)

  count_nor <- count
  gene_num <- ncol(count)
  sample_num <- nrow(count)

  # clustering
  dout <- 3

  if (cluster.method == "GMC")
  {
    count_binary <- matrix(0, nrow = sample_num, ncol = gene_num)

    for (i in 1:gene_num)
    {
      g <- count_nor[,i]

      mIQR <- IQR(g)
      list_large <- which(g > quantile(g, 0.75) + dout * mIQR)
      list_small <- which(g  == 0)

      if (mIQR == 0)
      {
        list_large <- NULL
      }

      list_remain <- setdiff(1:sample_num, union(list_large, list_small))

      if (length(list_remain) > 2)
      {
        if (min(count_nor[list_remain,i]) == max(count_nor[list_remain,i]))
        {
          count_binary[list_remain,i] <- 1
        }
        else
        {
          k <- Mclust(count_nor[list_remain,i], G = 2, verbose = FALSE)
          if (k$parameters$mean[1] > k$parameters$mean[2])
          {
            count_binary[list_remain,i] <- 3 - k$classification
          }
          else
          {
            count_binary[list_remain,i] <- k$classification
          }
        }
      }
      else if (length(list_remain) == 1)
      {
        count_binary[list_remain,i] <- sample(c(1,2), 1)
      }
      else if (length(list_remain) == 2)
      {
        count_binary[list_remain,i] <- order(count_nor[list_remain,i])
      }

      count_binary[list_large,i] <- 2
      count_binary[list_small,i] <- 1

      # If cannot split into two groups, assume that there is no spatial pattern
      if (sum(count_binary[,i]) == 2*sample_num)
      {
        count_binary[,i] <- sample(c(1,2),sample_num, replace = TRUE)
        print(i)
      }
      if (sum(count_binary[,i]) == sample_num)
      {
        count_binary[,i] <- sample(c(1,2),sample_num, replace = TRUE)
        print(i)
      }
    }
  }
  else if (cluster.method == "k-means")
  {
    count_binary <- matrix(0, nrow = sample_num, ncol = gene_num)

    for (i in 1:gene_num)
    {
      g <- count_nor[,i]

      mIQR <- IQR(g)
      list_large <- which(g > quantile(g, 0.75) + dout * mIQR)
      list_small <- NULL

      if (mIQR == 0)
      {
        list_large <- NULL
      }

      list_remain <- setdiff(1:sample_num, list_large)

      if (length(list_remain) > 2)
      {
        if (min(count_nor[list_remain,i]) == max(count_nor[list_remain,i]))
        {
          count_binary[list_remain,i] <- 1
        }
        else
        {
          k <- kmeans(count_nor[list_remain,i], 2, nstart = 25)
          if (k$centers[1] > k$centers[2])
          {
            count_binary[list_remain,i] <- 3 - k$cluster
          }
          else {
            count_binary[list_remain,i] <- k$cluster
          }
        }
      }
      else if (length(list_remain) == 1)
      {
        count_binary[list_remain,i] <- sample(c(1,2), 1)
      }
      else if (length(list_remain) == 1)
      {
        count_binary[list_remain,i] <- order(count_nor[list_remain,i])
      }

      count_binary[list_large,i] <- 2
      count_binary[list_small,i] <- 1

      if (max(count_binary[,i]) == 1 )
      {
        count_binary[sample(1:sample_num, floor(sample_num/2)), i] <- 2
      }
      else if (min(count_binary[,i]) == 2 )
      {
        count_binary[sample(1:sample_num, floor(sample_num/2)), i] <- 1
      }
    }
  }
  else if (cluster.method == "rank")
  {
    if (percentage.rank < 1 & percentage.rank > 0)
    {
      count_binary <- matrix(0, nrow = sample_num, ncol = gene_num)

      threshold_temp <- apply(count_nor, 2,
                              quantile,
                              probs = 1 - percentage.rank)

      for (i in 1:gene_num)
      {
        list_large <- which(count_nor[,i] > threshold_temp[i])
        list_small <- which(count_nor[,i] <= threshold_temp[i])
        if (length(list_large) == 0)
        {
          list_large <- sample(which(count_nor[,i] == threshold_temp[i]),
                               floor(sample_num*percentage.rank),
                               replace = FALSE)
        }

        count_binary[list_small,i] <- 1
        count_binary[list_large,i] <- 2

        if (max(count_binary[,i]) == 1 )
        {
          count_binary[sample(1:sample_num, floor(sample_num/2)), i] <- 2
        }
        else if (min(count_binary[,i]) == 2 )
        {
          count_binary[sample(1:sample_num, floor(sample_num/2)), i] <- 1
        }
      }
    }
    else
    {
      stop("value passed to 'percentage.rank' must be between 0 to 1")
    }
  }
  else
  {
    stop("value passed to 'cluster.method' is not a valid option")
  }

  count_binary <- as.vector(count_binary) - 1

  return(count_binary)
}
