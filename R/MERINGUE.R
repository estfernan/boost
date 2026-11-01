##' Run MERINGUE
##'
##' Wrapper for the R package MERINGUE (Miller et al., 2021). By default
##' (`approach = "differential"`), the STAr implementation in Supplementary
##' Notes S1 is followed: spots and genes with fewer than 100 total counts are
##' removed, counts are CPM-normalized and log-transformed, a spatially unaware
##' analysis is performed by PCA and graph-based clustering
##' (\eqn{k = 30} nearest neighbors), and differentially expressed genes for
##' each cluster are identified with `MERINGUE::getDifferentialGenes()`
##' (Wilcoxon rank-sum test). The minimum adjusted p-value across clusters is
##' used to define SVGs. With `approach = "spatial"`, MERINGUE's spatial
##' autocorrelation test (`MERINGUE::getSpatialPatterns()`) is used instead.
##'
##' @param count An \eqn{n}-by-\eqn{p} count matrix with spots in rows and
##'   genes in columns.
##' @param spots An \eqn{n}-by-2 matrix of spatial coordinates.
##' @param approach Either `"differential"` (STAr default) or `"spatial"`.
##' @param min.reads Minimum total counts for a gene to be retained.
##' @param min.lib.size Minimum total counts for a spot to be retained.
##' @param n.pcs Number of principal components used for clustering.
##' @param k Number of nearest neighbors for graph-based clustering.
##' @param filterDist Optional maximum distance for spatial neighbors
##'   (`approach = "spatial"` only).
##' @param verbose Logical; print progress messages.
##'
##' @return A data frame with one row per gene and columns `gene`, `method`,
##'   `score`, `score_type`, `pvalue`, and `padj`. For the differential
##'   approach, `padj` is the minimum adjusted p-value across clusters. The
##'   native output is attached as attribute `"native"`.
##'
##' @references Miller, B. F., et al. (2021). Characterizing spatial gene
##'   expression heterogeneity in spatially resolved single-cell
##'   transcriptomic data with nonuniform cellular densities.
##'   *Genome Research*, 31, 1843-1855.
##'
##' @export
##'
MERINGUE <- function(count, spots,
                     approach = c("differential", "spatial"),
                     min.reads = 100, min.lib.size = 100,
                     n.pcs = 5, k = 30, filterDist = NA,
                     verbose = FALSE)
{
  .check_svg_input(count, spots)
  .require_svg_package("MERINGUE", "remotes::install_github('JEFworks-Lab/MERINGUE')")
  approach <- match.arg(approach)

  spot.id <- .spot_names(count, spots)
  cd <- t(as.matrix(count))
  dimnames(cd) <- list(colnames(count), spot.id)
  pos <- as.matrix(spots)[, 1:2, drop = FALSE]
  dimnames(pos) <- list(spot.id, c("x", "y"))

  cd  <- MERINGUE::cleanCounts(counts = cd, min.reads = min.reads,
                               min.lib.size = min.lib.size, plot = FALSE,
                               verbose = verbose)
  mat <- MERINGUE::normalizeCounts(counts = cd, log = TRUE, verbose = verbose)
  mat <- as.matrix(mat)

  if (approach == "differential")
  {
    pcs <- stats::prcomp(t(mat), center = TRUE)$x[, seq_len(n.pcs), drop = FALSE]
    com <- MERINGUE::getClusters(pcs, k = k, verbose = verbose)
    dg  <- MERINGUE::getDifferentialGenes(mat, com, verbose = verbose)
    padj <- do.call(cbind, lapply(dg, function(d) d[rownames(mat), "p.adj"]))
    out <- .svg_table(rownames(mat), "MERINGUE", score_type = "pvalue",
                      padj = apply(padj, 1, min))
    attr(out, "native") <- list(clusters = com, differential = dg)
  }
  else
  {
    w   <- MERINGUE::getSpatialNeighbors(pos[colnames(mat), , drop = FALSE],
                                         filterDist = filterDist)
    res <- MERINGUE::getSpatialPatterns(mat, w, verbose = verbose)
    out <- .svg_table(rownames(res), "MERINGUE", score = res$observed,
                      score_type = "MoransI", pvalue = res$p.value,
                      padj = res$p.adj)
    attr(out, "native") <- res
  }

  out
}
