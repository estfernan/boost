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

##' Mouse Olfactory Bulb (MOB) Dataset
##'
##' The mouse olfactory bulb (MOB) study is a publicly available ST dataset
##' accessible on the Spatial Research Lab <http://www.spatialresearch.org>.
##' There are twelve replicates in the MOB study. Following the previous
##' studies (Li et al., 2020; Sun et al., 2020; Svensson et al., 2018), the
##' MOB replicate eleven, which contains 16,218 genes measured on 262 spots,
##' was saved.
##'
##' @docType data
##'
##' @name mob
##' @aliases mob mob.spots
##'
##' @format `mob` is a matrix with 262 spots and 16,218 genes that denotes
##' the expression count table. `mob.spots` is a matrix of the
##' geospatial profile.
##'
##' @references
##'
##' Li, Q., M. Zhang, Y. Xie, and G. Xiao (2020). Bayesian modeling of
##' spatial molecular profiling data via Gaussian process.
##' *arXiv preprint arXiv:2012.03326*.
##'
##' Sun, S., J. Zhu, and X. Zhou (2020). Statistical analysis of spatial
##' expression patterns for spatially resolved transcriptomic studies.
##' *Nature Methods 17* (2), 193–200.
##'
##' Svensson, V., S. A. Teichmann, and O. Stegle (2018).
##' SpatialDE: Identification of spatially variable genes.
##' *Nature Methods 15* (5), 343–346.
##'
##' @source Publicly available at <http://www.spatialresearch.org>.
##'
##' @examples
##' ## mouse olfactory bulb (MOB) dataset
##' data(mob)
##'
##' @keywords dataset
##'
NULL

##' Human Breast Cancer (BC) Dataset
##'
##' The human breast cancer (BC) study is a publicly available ST dataset
##' accessible on the Spatial Research Lab <http://www.spatialresearch.org>.
##' There are four layers in the BC study. Following the previous
##' studies (Li et al., 2020; Sun et al., 2020; Svensson et al., 2018), the
##' BC layer 2, which contains 14,789 genes measured on 250 spots, was saved.
##'
##' @docType data
##'
##' @name bc
##' @aliases bc bc.spots
##'
##' @format `bc` is a matrix with 250 spots and 14,789 genes that denotes
##' the expression count table. `bc.spots` is a matrix of the
##' geospatial profile.
##'
##' @references
##'
##' Li, Q., M. Zhang, Y. Xie, and G. Xiao (2020). Bayesian modeling of
##' spatial molecular profiling data via Gaussian process.
##' *arXiv preprint arXiv:2012.03326*.
##'
##' Sun, S., J. Zhu, and X. Zhou (2020). Statistical analysis of spatial
##' expression patterns for spatially resolved transcriptomic studies.
##' *Nature Methods 17* (2), 193–200.
##'
##' Svensson, V., S. A. Teichmann, and O. Stegle (2018).
##' SpatialDE: Identification of spatially variable genes.
##' *Nature Methods 15* (5), 343–346.
##'
##' @source Publicly available at <http://www.spatialresearch.org>.
##'
##' @examples
##' ## human breast cancer (BC) dataset
##' data(bc)
##'
##' @keywords dataset
##'
NULL
