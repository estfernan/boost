##' Run SpaGCN through reticulate
##'
##' SpaGCN (Hu et al., 2021) is Python software that identifies spatial
##' domains and then detects SVGs by domain-guided differential expression.
##' By default, the STAr implementation in Supplementary Notes S1 is used:
##' genes expressed in fewer than three spots and special genes (ERCC spike-ins
##' and mitochondrial genes) are removed; spatial domains are identified with
##' two clusters (to detect global SVGs) and refined; and, for each domain,
##' neighboring domains are found with `find_neighbor_clusters()` and
##' differential expression is tested with `rank_genes_groups()`. The minimum
##' adjusted p-value across the two domains is used to select SVGs.
##' Histology images are used when supplied (e.g., for 10x Visium data).
##'
##' @param count An \eqn{n}-by-\eqn{p} count matrix with spots in rows and
##'   genes in columns.
##' @param spots An \eqn{n}-by-2 matrix of spot array coordinates.
##' @param pixels Optional \eqn{n}-by-2 matrix of pixel coordinates on the
##'   histology image.
##' @param image Optional path to the histology image. When given together
##'   with `pixels`, histology information is used in the adjacency matrix.
##' @param n.clusters Number of spatial domains (default 2).
##' @param shape Spot arrangement used for domain refinement: `"square"` for
##'   ST data or `"hexagon"` for 10x Visium data.
##' @param seed Random seed.
##' @param max.epochs Maximum number of training epochs.
##' @param python Optional path to the Python executable used by reticulate.
##' @param runner Optional R function `(spagcn, count, spots, ...)` that
##'   replaces the default workflow; it receives the imported SpaGCN module.
##' @param ... Additional arguments passed to `runner`.
##'
##' @return For the default workflow, a data frame with one row per gene and
##'   columns `gene`, `method`, `score`, `score_type`, `pvalue`, and `padj`
##'   (minimum adjusted p-value across domains); the per-domain results and
##'   the refined domain labels are attached as attribute `"native"`. If
##'   `runner` is given, its result is returned.
##'
##' @details Requires the Python packages 'SpaGCN', 'scanpy', 'anndata', and
##'   'torch' in the reticulate environment.
##'
##' @references Hu, J., et al. (2021). SpaGCN: Integrating gene expression,
##'   spatial location and histology to identify spatial domains and spatially
##'   variable genes by graph convolutional network. *Nature Methods*, 18,
##'   1342-1351.
##'
##' @export
##'
SpaGCN <- function(count, spots, pixels = NULL, image = NULL,
                   n.clusters = 2, shape = c("square", "hexagon"),
                   seed = 100, max.epochs = 200, python = NULL,
                   runner = NULL, ...)
{
  .check_svg_input(count, spots)
  .require_python_module("SpaGCN", python)
  shape <- match.arg(shape)

  if (!is.null(runner))
  {
    if (!is.function(runner)) stop("'runner' must be a function", call. = FALSE)
    mod <- reticulate::import("SpaGCN", delay_load = FALSE)
    return(runner(mod, count, spots, ...))
  }

  histology <- !is.null(image) && !is.null(pixels)
  xy <- as.matrix(spots)[, 1:2, drop = FALSE]

  helper <- .star_python()
  res <- helper$run_spagcn_svg(
    counts     = as.matrix(count),
    genes      = colnames(count),
    spots      = .spot_names(count, spots),
    x_array    = as.numeric(xy[, 1]),
    y_array    = as.numeric(xy[, 2]),
    x_pixel    = if (histology) as.integer(pixels[, 1]) else NULL,
    y_pixel    = if (histology) as.integer(pixels[, 2]) else NULL,
    image      = if (histology) image else NULL,
    histology  = histology,
    n_clusters = as.integer(n.clusters),
    shape      = shape,
    seed       = as.integer(seed),
    max_epochs = as.integer(max.epochs)
  )

  tab <- as.data.frame(res$result)
  out <- .svg_table(tab$gene, "SpaGCN", score_type = "pvalue", padj = tab$padj)
  attr(out, "native") <- list(result = tab, domains = unlist(res$domains))
  out
}
