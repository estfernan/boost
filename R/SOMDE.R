##' Run SOMDE through reticulate
##'
##' SOMDE (Hao et al., 2021) is distributed as a Python package. Following
##' STAr (Supplementary Notes S1), a self-organizing map is trained on the
##' spatial coordinates with `somn = 20`, node-level expression is generated
##' with `mtx()`, normalized, and SOMDE is run with its default settings.
##' Q-values are used as the criterion for SVG detection.
##'
##' @param count An \eqn{n}-by-\eqn{p} count matrix with spots in rows and
##'   genes in columns.
##' @param spots An \eqn{n}-by-2 matrix of spatial coordinates.
##' @param k SOM size parameter (`somn`; default 20).
##' @param python Optional path to the Python executable used by reticulate.
##'
##' @return A data frame with one row per gene and columns `gene`, `method`,
##'   `score` (log-likelihood ratio), `score_type`, `pvalue`, and `padj`
##'   (q-value). The native SOMDE output is attached as attribute `"native"`.
##'
##' @details Requires the Python package 'somde' (`pip install somde`) in the
##'   reticulate environment.
##'
##' @references Hao, M., Hua, K., & Zhang, X. (2021). SOMDE: a scalable method
##'   for identifying spatially variable genes with self-organizing map.
##'   *Bioinformatics*, 37(23), 4392-4398.
##'
##' @export
##'
SOMDE <- function(count, spots, k = 20, python = NULL)
{
  .check_svg_input(count, spots)
  .require_python_module("somde", python)

  helper <- .star_python()
  res <- helper$run_somde(
    counts = as.matrix(count),
    genes  = colnames(count),
    spots  = .spot_names(count, spots),
    coords = as.matrix(spots)[, 1:2, drop = FALSE],
    k      = as.integer(k)
  )
  res <- as.data.frame(res)

  out <- .svg_table(res$g, "SOMDE", score = res$LLR, score_type = "LLR",
                    pvalue = res$pval, padj = res$qval)
  attr(out, "native") <- res
  out
}

.star_python <- function()
{
  path <- system.file("python", package = "boost")
  reticulate::import_from_path("star_svg", path = path, delay_load = FALSE)
}
