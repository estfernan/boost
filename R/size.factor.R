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

##' Get Size Factor Estimate
##'
##' Obtain the size factor estimates to capture all nuisance effects and
##' compute the relative gene expression levels.
##'
##' Normalization is critical to sequence count data analysis. As a result, to
##' counteract various artifacts and bias due to biological and technical
##' reasons, the read counts should be converted to their relative gene
##' expression levels i.e., divide the counts by the result of this function.
##' If the main interest is in the absolute gene expression level, then
##' the size factor is set to the same value; otherwise, the size factor is
##' computed directly from the gene expression count data.
##'
##' This method offers the following options for computing size factors:
##' total sum scaling (TSS), upper-quartiles (Q75), relative log
##' expression (RLE), and weighted trimmed mean by M-values (TMM).
##'
##' See Jiang et al. (2021) for more information and references on
##' these options.
##'
##' @param count An \eqn{n}-by-\eqn{p} numeric matrix \eqn{Y} that denotes the
##'   absolute gene expression count table. Each entry is the absolute read
##'   count for gene \eqn{j} collected at spot \eqn{i}.
##' @param estimation.method An optional character string to specify the
##'   size factor computation technique. The default is "TSS" for
##'   total sum scaling.
##'
##' @return A numeric vector, where each entry denotes the size factor for
##'   sample \eqn{i}.
##'
##' @references Jiang, X., Li, Q., & Xiao, G. (2021). Bayesian Modeling of
##'   Spatial Transcriptomics Data via a Modified Ising Model.
##'   *arXiv preprint arXiv:2104.13957*.
##'
##' @export
##' @keywords preprocessing
##'
##' @importFrom stats quantile
##'
get.size.factor <- function(
  count,
  estimation.method = c("TSS", "Q75", "RLE", "TMM")
)
{
  ##
  ## NOTE: The code below has not been thoroughly checked and reformatted for
  ##         a better understanding
  ##

  if (length(estimation.method) > 1)
  {
    estimation.method <- estimation.method[1]
  }

  sample_num <- nrow(count)
  count_rowsum <- rowSums(count)

  if (estimation.method == "TSS")
  {
    ##
    ## TSS (Total Sum Scaling)
    ##

    raw_s_factors <- count_rowsum
  }
  else if (estimation.method == "Q75")
  {
    ##
    ## Q75 (Upper Quantile normalization)
    ##

    ## function for calculating non-zero quantile
    non_zero_quantile <- function(x, probs)
    {
      quantile(x[x > 0], probs = probs)
    }

    count_q75 <- apply(count, 1, non_zero_quantile, probs = 0.75)

    count_N <- count_rowsum / sample_num
    raw_s_factors <- count_q75/count_N
  }
  else if (estimation.method == "RLE")
  {
    ##
    ## RLE (Relative Log Expression normalization)
    ##

    ## function for calculating the geometric mean
    geo_mean <- function(x)
    {
      exp(sum(log(x[x > 0])) / length(x))
    }

    ## function for calculating non-zero median
    non_zero_median <- function(x)
    {
      median(x[x > 0])
    }
    ref_sample <- apply(count, 2, geo_mean)
    norm_rle_1 <- sweep(count, 2, ref_sample, FUN = "/")
    raw_s_factors <- apply(as.matrix(norm_rle_1), 1, non_zero_median)
  }
  else if (estimation.method == "TMM")
  {
    ##
    ## TMM (Trimmed Mean Method)
    ##

    count_t <- t(count)
    # raw_s_factors <- calcNormFactors(count_t, method = "TMM")

    raw_s_factors <- calcNormFactors_copy(count_t, method = estimation.method,
                                          refColumn = NULL,
                                          logratioTrim = 0.3, sumTrim = 0.05,
                                          doWeighting = TRUE, Acutoff = -1e10)
  }
  else
  {
    warning(
      paste0("value passed to 'estimation.method' is not a valid ",
             "normalization technique,", "\n  ",
             "the size factor for relative gene expression levels",
             "will be returned")
    )

    raw_s_factors <- rep(1, sample_num)
  }

  return(raw_s_factors)
}
