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
##' pattern, also known as spatially variable genes (SVGs), in spatial
##' transcriptomics data, and is the companion R package of the Spatial
##' Transcriptomics Arena (STAr, <https://lce.biohpc.swmed.edu/star/>). The
##' Bayesian models BOOST-GP and BOOST-MI are implemented natively together
##' with SPARK, SpatialDE, and BinSpect; nnSVG, SPARK-X, MERINGUE, Trendsceek,
##' Seurat Moran's I, SOMDE, and SpaGCN are available through wrappers that
##' call the original software. GEO and MC2 rank aggregation, data filtering,
##' normalization, binarization, neighbor construction, visualization, and
##' performance metrics are also provided. See [svg.methods()] for an
##' overview and <https://boost-r.readthedocs.io/en/latest/> for tutorials.
##'
##' See Jiang et al. (2022) and Li et al. (2021) for details of the Bayesian
##' models.
##'
##' @references
##'
##' Jiang, X., Xiao, G., & Li, Q. (2022). A Bayesian modified Ising
##' model for identifying spatially variable genes from spatial transcriptomics
##' data. *Statistics in Medicine*, 41(23), 4647-4665.
##'
##' Li, Q., Zhang, M., Xie, Y., & Xiao, G. (2021). Bayesian modeling of
##' spatial molecular profiling data via Gaussian process.
##' *Bioinformatics*, 37(22), 4129-4136.
##'
##' @useDynLib boost, .registration = TRUE
##' @importFrom Rcpp sourceCpp
##'
##' @name boost
##' @aliases boost-package
##'
"_PACKAGE"
