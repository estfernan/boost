##' Run nnSVG
##'
##' Wrapper for the Bioconductor implementation of nnSVG (Weber et al., 2023),
##' following the STAr settings (Supplementary Notes S1): low-expressed and
##' mitochondrial genes are removed with [nnSVG::filter_genes()] using its
##' default parameters, library size factors and log-normalized expression are
##' computed with `scuttle::computeLibraryFactors()` and
##' `scuttle::logNormCounts()`, and nnSVG is then run with its default
##' settings. Genes with BH-adjusted p-values below 0.05 are SVGs.
##'
##' @param count An \eqn{n}-by-\eqn{p} count matrix with spots in rows and
##'   genes in columns (column names are gene names).
##' @param spots An \eqn{n}-by-2 matrix of spatial coordinates.
##' @param n_neighbors Number of nearest neighbors used by nnSVG (default 10).
##' @param filter Logical; apply [nnSVG::filter_genes()] with its default
##'   settings before fitting (default `TRUE`).
##' @param n_threads Number of threads passed to [nnSVG::nnSVG()].
##' @param ... Additional arguments passed to [nnSVG::nnSVG()].
##'
##' @return A data frame with one row per gene and columns `gene`, `method`,
##'   `score` (likelihood-ratio statistic), `score_type`, `pvalue`, and `padj`.
##'   The native nnSVG `rowData` is attached as attribute `"native"`.
##'
##' @references Weber, L. M., et al. (2023). nnSVG for the scalable
##'   identification of spatially variable genes using nearest-neighbor
##'   Gaussian processes. *Nature Communications*, 14, 4059.
##'
##' @export
##'
nnSVG <- function(count, spots, n_neighbors = 10, filter = TRUE,
                  n_threads = 1, ...)
{
  .check_svg_input(count, spots)
  .require_svg_package("nnSVG", "BiocManager::install('nnSVG')")
  .require_svg_package("SpatialExperiment", "BiocManager::install('SpatialExperiment')")
  .require_svg_package("SummarizedExperiment", "BiocManager::install('SummarizedExperiment')")
  .require_svg_package("scuttle", "BiocManager::install('scuttle')")

  spot.id <- .spot_names(count, spots)
  genes   <- colnames(count)

  counts_gs <- t(as.matrix(count))
  dimnames(counts_gs) <- list(genes, spot.id)
  coords <- as.matrix(spots)[, 1:2, drop = FALSE]
  dimnames(coords) <- list(spot.id, c("x", "y"))

  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts_gs),
    rowData = data.frame(gene_id = genes, gene_name = genes, row.names = genes),
    colData = data.frame(sample_id = rep("sample01", length(spot.id)),
                         row.names = spot.id),
    spatialCoords = coords
  )

  if (filter) spe <- nnSVG::filter_genes(spe)
  spe <- scuttle::computeLibraryFactors(spe)
  spe <- scuttle::logNormCounts(spe)
  spe <- nnSVG::nnSVG(spe, n_neighbors = n_neighbors, n_threads = n_threads, ...)

  rd  <- as.data.frame(SummarizedExperiment::rowData(spe))
  out <- .svg_table(rownames(spe), "nnSVG", score = rd$LR_stat,
                    score_type = "LR_stat", pvalue = rd$pval, padj = rd$padj)
  attr(out, "native") <- rd
  out
}
