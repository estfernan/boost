##' List SVG Detection and Rank-Aggregation Methods in boost
##'
##' Returns a catalogue of the 12 SVG detection methods benchmarked in STAr and
##' the two rank-aggregation procedures, together with the package function,
##' implementation type, expected input, native evidence measure, and the SVG
##' selection criterion used in STAr.
##'
##' @return A data frame with columns `method`, `function_name`, `category`,
##'   `implementation` (`"native"`, `"R wrapper"`, or `"Python wrapper"`),
##'   `input`, `statistic`, and `svg_criterion`.
##'
##' @examples
##' svg.methods()
##'
##' @export
svg.methods <- function()
{
  data.frame(
    method = c("BOOST-GP", "BOOST-MI", "SPARK", "SpatialDE", "BinSpect",
               "SpaGCN", "nnSVG", "SOMDE", "Trendsceek", "SPARK-X",
               "MERINGUE", "Seurat-Moran's I", "GEO", "MC2"),
    function_name = c("BOOST.GP", "BOOST.MI", "SPARK", "SpatialDE", "binSpect",
                      "SpaGCN", "nnSVG", "SOMDE", "Trendsceek", "SPARK.X",
                      "MERINGUE", "Seurat.MoransI", "rank.aggregation",
                      "rank.aggregation"),
    category = c(rep("SVG detection", 12), rep("Rank aggregation", 2)),
    implementation = c(rep("native", 5), "Python wrapper", "R wrapper",
                       "Python wrapper", "R wrapper", "R wrapper", "R wrapper",
                       "R wrapper", "native", "native"),
    input = c("counts (one gene) + size factors", "binarized expression (one gene) + neighbors",
              "counts (one gene) + size factors", "normalized expression (one gene)",
              "binarized expression (one gene) + neighbors", "counts (all genes)",
              "counts (all genes)", "counts (all genes)", "counts (all genes)",
              "counts (all genes)", "counts (all genes)", "normalized expression (all genes)",
              "ranked method results", "ranked method results"),
    statistic = c("Bayes factor", "Bayes factor", "p-value", "p-value", "p-value",
                  "p-value", "p-value", "q-value", "p-value", "p-value",
                  "p-value", "p-value", "aggregated rank", "aggregated rank"),
    svg_criterion = c("BF > 150", "BF > 150", rep("adjusted p-value < 0.05", 10),
                      "top-ranked genes", "top-ranked genes"),
    stringsAsFactors = FALSE
  )
}

##' Convert a boost Method Result to the STAr Standard Format
##'
##' Converts a single-gene result returned by BOOST-GP, BOOST-MI/BOOST-Ising,
##' SPARK, SpatialDE, or BinSpect into a one-row data frame with the same
##' columns as the multi-gene wrappers. The native evidence measure of each
##' method is retained: Bayes factors for BOOST-GP and BOOST-MI (stored in
##' `score`, with `pvalue` left missing) and p-values for the other methods.
##' The original model object is not modified.
##'
##' @param x A fitted object returned by one of the gene-level SVG methods.
##' @param bf For BOOST-MI, which Bayes factor to report: `"max"` (default),
##'   the larger of the attraction and repulsion Bayes factors, as reported
##'   for real datasets in STAr; `"attraction"`; or `"repulsion"`.
##'
##' @return A one-row data frame with columns `gene`, `method`, `score`,
##'   `score_type`, `pvalue`, and `padj` (`padj` is missing because
##'   multiple-testing adjustment requires all genes; see [svg.select()]).
##'
##' @export
as.svg.result <- function(x, bf = c("max", "attraction", "repulsion"))
{
  bf   <- match.arg(bf)
  gene <- if (is.null(x$gene.name)) NA_character_ else as.character(x$gene.name)

  if (inherits(x, "BOOST.GP")) {
    return(.svg_table(gene, "BOOST-GP", score = x$measures$BF, score_type = "BF"))
  }

  if (inherits(x, "BOOST.Ising")) {
    m <- x$measures
    score <- switch(bf,
                    max        = max(m$BF.neg, m$BF.pos),
                    attraction = m$BF.neg,
                    repulsion  = m$BF.pos)
    return(.svg_table(gene, "BOOST-MI", score = score, score_type = "BF"))
  }

  if (inherits(x, "SPARK")) {
    return(.svg_table(gene, "SPARK", pvalue = x$measures[["p.value"]]))
  }

  if (inherits(x, "SpatialDE")) {
    return(.svg_table(gene, "SpatialDE", pvalue = x$measures$p.value))
  }

  if (inherits(x, "binSpect")) {
    return(.svg_table(gene, "BinSpect", score = x$measures$OR,
                      score_type = "pvalue", pvalue = x$measures$p.val))
  }

  stop("unsupported object class: expected BOOST.GP, BOOST.Ising, SPARK, ",
       "SpatialDE, or binSpect", call. = FALSE)
}

##' Select SVGs Using Method-Appropriate Criteria
##'
##' Applies the SVG selection criteria used in STAr to a table of results in
##' the standard format (from [as.svg.result()] or the method wrappers):
##' Bayes factor greater than 150 for BOOST-GP and BOOST-MI, and
##' Benjamini-Hochberg adjusted p-value below 0.05 for the other methods.
##' Missing adjusted p-values are computed within each method from `pvalue`.
##'
##' @param results A data frame in the standard format, possibly containing
##'   several methods.
##' @param bf.threshold Bayes-factor threshold (default 150).
##' @param alpha Significance level for adjusted p-values (default 0.05).
##'
##' @return `results` with the `padj` column completed and a logical column
##'   `is.svg`.
##'
##' @export
svg.select <- function(results, bf.threshold = 150, alpha = 0.05)
{
  results <- as.data.frame(results, stringsAsFactors = FALSE)
  results$is.svg <- FALSE

  for (m in unique(results$method))
  {
    idx <- which(results$method == m)
    if (all(results$score_type[idx] == "BF"))
    {
      results$is.svg[idx] <- results$score[idx] > bf.threshold
    }
    else
    {
      miss <- is.na(results$padj[idx])
      if (any(miss))
      {
        results$padj[idx] <- stats::p.adjust(results$pvalue[idx], method = "BH")
      }
      results$is.svg[idx] <- results$padj[idx] < alpha
    }
  }

  results$is.svg[is.na(results$is.svg)] <- FALSE
  results
}

##' Rank Genes Within Each Method and Aggregate with GEO or MC2
##'
##' Converts results in the standard format into within-method ranks (larger
##' Bayes factors and smaller p-values rank higher) and aggregates them with
##' [rank.aggregation()]. In STAr, GEO and MC2 aggregate BOOST-GP, BOOST-MI,
##' BinSpect, SPARK, and SpatialDE using the union of the top 1,000 genes from
##' each method; SpaGCN is excluded because its domain-guided rankings are not
##' directly comparable.
##'
##' @param results A data frame in the standard format with results from two
##'   or more methods.
##' @param K Number of top genes considered from each method (1,000 in STAr).
##' @param method Either `"GEO"` or `"MC2"`.
##' @param ties.method Tie handling passed to [rank.aggregation()].
##'
##' @return A data frame with columns `gene`, `method`, and `rank`.
##'
##' @export
svg.aggregate <- function(results, K = 1000, method = c("GEO", "MC2"),
                          ties.method = c("average", "first", "last",
                                          "random", "max", "min"))
{
  method      <- match.arg(method)
  ties.method <- match.arg(ties.method)

  ranks <- svg.rank.table(results)
  K     <- min(K, nrow(ranks))
  out   <- rank.aggregation(data = ranks, K = K, method = method,
                            ties.method = ties.method)
  out$method <- method
  rownames(out) <- NULL
  out[, c("gene", "method", "rank")]
}

##' Build a Gene-by-Method Rank Table
##'
##' Converts results in the standard format into the wide rank table expected
##' by [rank.aggregation()]: one row per gene, a first column `gene`, and one
##' column of within-method ranks per method (1 = strongest evidence).
##'
##' @param results A data frame in the standard format.
##'
##' @return A data frame with column `gene` followed by one rank column per
##'   method.
##'
##' @export
svg.rank.table <- function(results)
{
  results <- as.data.frame(results, stringsAsFactors = FALSE)
  genes   <- sort(unique(results$gene))
  out     <- data.frame(gene = genes, stringsAsFactors = FALSE)

  for (m in unique(results$method))
  {
    r <- results[results$method == m, , drop = FALSE]
    r <- r[!duplicated(r$gene), , drop = FALSE]
    evidence <- if (all(r$score_type == "BF")) -r$score
                else if (!all(is.na(r$pvalue))) r$pvalue
                else r$padj
    rk <- rank(evidence, ties.method = "average", na.last = "keep")
    out[[make.names(m)]] <- rk[match(genes, r$gene)]
  }

  out
}
