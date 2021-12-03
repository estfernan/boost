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

##' Normalize Sequence Count Data
##'
##' Obtain the relative expression levels for the gene expression count table.
##'
##' Normalization is critical to the development of analysis techniques for any
##' sequence count data that suffer from various sequence artifacts and biases.
##' This method provides eight choices of normalization methods in two
##' categories that attempt to counteract those biases due to biological and
##' technical reasons. The first type is based on size factor estimation
##' (i.e., TSS, Q75, RLE, TMM) where a quantity is obtained, for each
##' sample, that captures all nuisance effects. The other type of normalization
##' method  is based on variance-stabilizing transformation (VST)
##' (i.e., A-VST, N-VST, log-VST), which aims to transform a random variable
##' with a negative binomial distribution into one with an approximately
##' normal distribution.
##'
##' See Jiang et al. (2021) for more information on the normalization methods.
##'
##' @param count An \eqn{n}-by-\eqn{p} integer matrix \eqn{Y} that denotes the
##'   absolute gene expression count table. Each entry is the read count for
##'   gene \eqn{j} collected at spot \eqn{i}.
##' @param scaling.method  An optional character string to specify the
##'   normalization technique. The default is "TSS" for the for total sum
##'   scaling.
##'
##' @return A numeric matrix that denotes the relative expression levels.
##'
##' @references Jiang, X., Li, Q., & Xiao, G. (2021). Bayesian Modeling of
##'   Spatial Transcriptomics Data via a Modified Ising Model.
##'   *arXiv preprint arXiv:2104.13957*.
##'
##' @examples
##' ## Need to implement the example.
##'
##' @seealso
##' [st.plot()] for plotting the relative expression levels.
##'
##' @export
##' @keywords preprocessing
##'
##' @importFrom stats quantile median var coef nls resid lm
##'
normalize.st <- function(
  count,
  scaling.method = c("TSS", "Q75", "RLE", "TMM", "A-VST", "N-VST", "log-VST")
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

  gene_num   <- ncol(count)
  sample_num <- nrow(count)

  N <- rowSums(count)

  if (scaling.method == "TSS")
  {
    ##
    ## TSS(Total Sum Scaling)
    ##

    ### scale-factors
    raw_s_factors <- N
    scale_coeff <- exp((-1/nrow(count)) * sum(log(raw_s_factors)))
    scaled_s_factors <- scale_coeff * raw_s_factors

    ### normalized count matrix
    db.norm <- sweep(count, 1, N, FUN = "/")
    count_nor <- db.norm
  }
  else if (scaling.method == "Q75")
  {
    ##
    ## Q75(Upper Quantile normalization)
    ##

    ### scale factors
    count_q75 <- apply(count, 1,function(x) { quantile(x[x > 0],0.75) })
    count_N <- N / nrow(count)
    raw_s_factors <- count_q75/count_N
    scale_coeff <- exp((-1/nrow(count)) * sum(log(raw_s_factors)))
    scaled_s_factors <- scale_coeff * raw_s_factors

    ### normalized count matrix
    db.norm <- sweep(count, 1, raw_s_factors, FUN = "/")
    count_nor <- db.norm
  }
  else if (scaling.method == "RLE")
  {
    ##
    ## RLE(Relative Log Expression normalization)
    ##

    ### scale_factors

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
    scale_coeff <- exp((-1/nrow(count)) * sum(log(raw_s_factors)))
    scaled_s_factors <- scale_coeff * raw_s_factors

    ### normalized count matrix
    db.norm <- sweep(count, 1, raw_s_factors, FUN = "/")
    count_nor <- db.norm
  }
  else if (scaling.method == "TMM")
  {
    ##
    ## TMM(Trimmed Mean Method)
    ##

    ### scale_factors
    count_t <- t(count)
    # raw_s_factors <- calcNormFactors(count_t,method = "TMM")
    raw_s_factors <- calcNormFactors_copy(count_t, method = "TMM", # scaling.method,
                                          refColumn = NULL,
                                          logratioTrim = 0.3, sumTrim = 0.05,
                                          doWeighting = TRUE, Acutoff = -1e10)
    scale_coeff <- exp((-1/nrow(count)) * sum(log(raw_s_factors)))
    scaled_s_factors <- scale_coeff * raw_s_factors

    # normalized count matrix
    db.norm <- sweep(count, 1, raw_s_factors, FUN = "/")
    count_nor <- db.norm
  }

  else if (scaling.method == "N-VST") {
    ##
    ## Naive Transformation (VST)
    ##

    varx = apply(count, 2, var)
    meanx = apply(count, 2, mean)
    phi = coef(nls(varx ~ meanx + phi * meanx^2, start = list(phi = 1)))
    y <- sqrt( phi * count)
    naive_norm <- log(y + sqrt(1 + y^2))
    total_n_counts <- apply(count, 1, sum)
    log_total_n_counts <- log(total_n_counts)
    db.norm <- apply(naive_norm, 2, function(x){resid(lm(x ~ log_total_n_counts))} )

    ## All the above normalized counts were negative so reversed their signs
    count_nor <- db.norm
  }
  else if (scaling.method == "A-VST")
  {
    ##
    ## Anscombe's Transformation (VST)
    ##

    varx = apply(count, 2, var)
    meanx = apply(count, 2, mean)
    phi = coef(nls(varx ~ meanx + phi * meanx^2, start = list(phi = 1)))

    ### Took the absolute value of counts as negative under the square root not on Real subspace
    y <- sqrt( abs( (count + (3/8))/((1/phi)-(3/4)) ) )
    anscombe_norm <- log(y + sqrt(1 + y^2))
    total_a_counts <- apply(count, 1, sum)
    log_total_a_counts <- log(total_a_counts)
    db.norm <- apply(anscombe_norm, 2, function(x){resid(lm(x ~ log_total_a_counts))} )

    ## All the above normalized counts were negative so reversed their signs
    count_nor <- db.norm
  }
  else if(scaling.method == "log-VST")
  {
    ##
    ## Log Transformation
    ##
    varx = apply(count, 2, var)
    meanx = apply(count, 2, mean)
    phi = coef(nls(varx ~ meanx + phi * meanx^2, start = list(phi = 1)))
    log_norm <- log2(count + (1/(2*phi)))
    total_l_counts <- apply(count, 1, sum)
    log_total_l_counts <- log(total_l_counts)
    db.norm <- apply(log_norm, 2, function(x){resid(lm(x ~ log_total_l_counts))} )

    ## All the above normalized counts were negative so reversed their signs
    count_nor <- db.norm
  }
  else
  {
    warning(
      paste0("value passed to 'estimation.method' is not a valid ",
             "normalization technique,", "\n  ",
             "the relative gene expression levels will be returned")
    )

    count_nor <- count
  }

  dimnames(count_nor) <- dimnames(count)

  return(count_nor)
}
