##
## R package boost by Esteban Fernández, Xi Jiang, Suhana Bedi, and Qiwei Li
## Copyright (C) 2021
##
## This file is part of the R package boost.
##
## The R package boost is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package boost is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##

##' boost: Bayesian Modeling of Spatial Transcriptomics Data
##'
##' This package provides functions to detect genes with spatial expression
##' pattern, also known as spatially variable (SV) genes, in spatial
##' transcriptomics (ST) data. Two major and novel Bayesian models are
##' implemented via a Gaussian process or Ising model. In addition, it also
##' provides other standard statistical tools such as SPARK, binSpect, etc.
##' Utilities are available to normalize count data, dichotomise expression
##' levels, get spatial neighbors, and view the results of the procedures.
##'
##' See Jiang et al. (2021) and Li et al. (2020) for details and illustrations
##' of how the Bayesian models are fitted and their applications.
###'
##' @references
##'
##' Jiang, X., Li, Q., & Xiao, G. (2021). Bayesian Modeling of
##' Spatial Transcriptomics Data via a Modified Ising Model.
##' *arXiv preprint arXiv:2104.13957*.
##'
##' Li, Q., Zhang M., Xie Y., & Xiao, G. (2020). Bayesian Modeling of
##' Spatial Molecular Profiling Data via Gaussian Process.
##' *arXiv preprint arXiv:2012.03326*.
##'
##' @useDynLib boost, .registration = TRUE
##' @importFrom Rcpp sourceCpp
##'
##' @docType package
##' @name boost
##'
NULL
