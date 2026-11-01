##' Run Seurat Moran's I
##'
##' Computes Moran's I for each gene with `Seurat::RunMoransI()` (Hao et al.,
##' 2024), which uses inverse squared distance weights. Following STAr
##' (Supplementary Notes S1), p-values are adjusted with the Benjamini-Hochberg
##' procedure.
##'
##' @param expr An \eqn{n}-by-\eqn{p} expression matrix with spots in rows and
##'   genes in columns (e.g., normalized expression).
##' @param spots An \eqn{n}-by-2 matrix of spatial coordinates.
##' @param verbose Logical; display Seurat progress messages.
##'
##' @return A data frame with one row per gene and columns `gene`, `method`,
##'   `score` (observed Moran's I), `score_type`, `pvalue`, and `padj`.
##'
##' @details `Seurat::RunMoransI()` requires the 'Rfast2' or 'ape' package.
##'
##' @export
##'
Seurat.MoransI <- function(expr, spots, verbose = FALSE)
{
  .check_svg_input(expr, spots)
  .require_svg_package("Seurat", "install.packages('Seurat')")

  data_gs <- t(as.matrix(expr))
  rownames(data_gs) <- colnames(expr)

  ans <- Seurat::RunMoransI(
    data = data_gs,
    pos = as.matrix(spots)[, 1:2, drop = FALSE],
    verbose = verbose
  )

  .svg_table(rownames(ans), "Seurat-Moran's I", score = ans$observed,
             score_type = "MoransI", pvalue = ans$p.value,
             padj = stats::p.adjust(ans$p.value, method = "BH"))
}
