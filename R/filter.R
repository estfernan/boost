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

##' Gene-Filtering Criteria Used in STAr
##'
##' Returns the spot- and gene-level quality-control criteria applied before
##' each SVG detection method in STAr (Section 3.4 and Supplementary Notes S1
##' of the STAr paper). When the original publication of a method specifies a
##' gene-filtering rule, that rule is used; otherwise, genes expressed in
##' fewer than 10% of spatial locations are excluded.
##'
##' @param method An optional character string naming a method (e.g.,
##'   `"BOOST-MI"`). If `NULL` (the default), the criteria for all methods are
##'   returned.
##'
##' @return A data frame with columns `method`, `min.total` (minimum total
##'   count per spot), `min.percentage` (minimum proportion of spots in which a
##'   gene must be expressed), `min.spots` (minimum number of spots in which a
##'   gene must be expressed), and `source`.
##'
##' @seealso [filter.st()] for applying the criteria.
##'
##' @examples
##' svg.filter.criteria()
##' svg.filter.criteria("BOOST-MI")
##'
##' @export
##'
svg.filter.criteria <- function(method = NULL)
{
  tab <- data.frame(
    method = c("BOOST-GP", "BOOST-MI", "SPARK", "SpatialDE", "BinSpect",
               "SpaGCN", "default"),
    min.total = rep(10, 7),
    min.percentage = c(0.10, 0.20, 0.10, 0.10, 0.10, 0, 0.10),
    min.spots = c(0, 0, 0, 0, 0, 3, 0),
    source = c("original publication", "original publication",
               "original publication", "STAr default", "STAr default",
               "original publication", "STAr default"),
    stringsAsFactors = FALSE
  )

  if (is.null(method)) return(tab)

  method <- .svg_method_name(method)
  if (!(method %in% tab$method)) method <- "default"
  tab[tab$method == method, , drop = FALSE]
}

##' Filter Spots and Genes
##'
##' Remove low-quality spots and rarely expressed genes from a spot-by-gene
##' count matrix before SVG detection. Spots are filtered first, and the
##' gene-level criteria are then evaluated on the retained spots.
##'
##' The thresholds can be given directly or taken from the criteria used in
##' STAr for a given method (see [svg.filter.criteria()]): spots with fewer
##' than ten total counts are removed; genes
##' expressed in fewer than 10% of spots are removed for BOOST-GP, SPARK,
##' SpatialDE and BinSpect; genes with zero counts in more than 80% of spots
##' are removed for BOOST-MI; and genes expressed in fewer than three spots are
##' removed for SpaGCN.
##'
##' @param count An \eqn{n}-by-\eqn{p} count matrix with spots in rows and
##'   genes in columns.
##' @param spots An \eqn{n}-by-2 matrix of spatial coordinates, in the same
##'   spot order as `count`.
##' @param method An optional character string naming an SVG detection method.
##'   When given, the STAr criteria for that method are used and the other
##'   threshold arguments are ignored.
##' @param min.total Minimum total count across genes for a spot to be kept.
##' @param min.percentage Minimum proportion of (retained) spots in which a
##'   gene must have a non-zero count to be kept.
##' @param min.spots Minimum number of (retained) spots in which a gene must
##'   have a non-zero count to be kept.
##'
##' @return A list with components `count` (the filtered count matrix),
##'   `spots` (the filtered coordinates), `kept.spots` and `kept.genes`
##'   (logical vectors over the input spots and genes), and `criteria`.
##'
##' @seealso [svg.filter.criteria()]
##'
##' @examples
##' data(mob)
##' res <- filter.st(mob, mob.spots, method = "BOOST-MI")
##' dim(res$count)
##'
##' @export
##' @keywords preprocessing
##'
filter.st <- function(
  count, spots,
  method = NULL,
  min.total = 10,
  min.percentage = 0.10,
  min.spots = 0
)
{
  count <- as.matrix(count)
  spots <- as.matrix(spots)

  if (nrow(count) != nrow(spots))
  {
    stop("'count' (spots x genes) and 'spots' must have the same number of rows")
  }

  if (!is.null(method))
  {
    crit           <- svg.filter.criteria(method)
    min.total      <- crit$min.total
    min.percentage <- crit$min.percentage
    min.spots      <- crit$min.spots
  }

  kept.spots <- rowSums(count) >= min.total
  count_f    <- count[kept.spots, , drop = FALSE]

  n.expr     <- colSums(count_f > 0)
  kept.genes <- (n.expr >= min.percentage * nrow(count_f)) & (n.expr >= min.spots)

  list(
    count      = count_f[, kept.genes, drop = FALSE],
    spots      = spots[kept.spots, , drop = FALSE],
    kept.spots = kept.spots,
    kept.genes = kept.genes,
    criteria   = list(min.total = min.total,
                      min.percentage = min.percentage,
                      min.spots = min.spots)
  )
}

.svg_method_name <- function(method)
{
  key <- toupper(gsub("[^A-Za-z0-9]", "", method))
  map <- c(BOOSTGP = "BOOST-GP", BOOSTMI = "BOOST-MI", BOOSTISING = "BOOST-MI",
           SPARK = "SPARK", SPATIALDE = "SpatialDE", BINSPECT = "BinSpect",
           SPAGCN = "SpaGCN", NNSVG = "nnSVG", SPARKX = "SPARK-X",
           MERINGUE = "MERINGUE", TRENDSCEEK = "Trendsceek", SOMDE = "SOMDE",
           SEURATMORANSI = "Seurat-Moran's I", MORANSI = "Seurat-Moran's I")
  if (key %in% names(map)) unname(map[key]) else method
}
