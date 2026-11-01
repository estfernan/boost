## Internal helpers shared by the SVG method wrappers.
##
## All wrappers follow the package-wide convention: expression/count matrices
## have spots in rows and genes in columns, and `spots` has one row per spot.

.check_svg_input <- function(count, spots)
{
  if (is.null(dim(count)) || length(dim(count)) != 2L)
    stop("'count' must be a two-dimensional spots x genes matrix", call. = FALSE)
  if (is.null(dim(spots)) || ncol(as.matrix(spots)) < 2L)
    stop("'spots' must have one row per spot and at least two coordinate columns",
         call. = FALSE)
  if (nrow(count) != nrow(spots))
    stop("nrow(count) must equal nrow(spots): 'count' must be spots x genes ",
         "(spots in rows, genes in columns)", call. = FALSE)
  if (is.null(colnames(count)))
    stop("'count' must have gene names as column names", call. = FALSE)
  invisible(TRUE)
}

.spot_names <- function(count, spots)
{
  nm <- rownames(count)
  if (is.null(nm)) nm <- rownames(spots)
  if (is.null(nm)) nm <- paste0("spot", seq_len(nrow(count)))
  make.unique(as.character(nm))
}

.require_svg_package <- function(package, install)
{
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. Install with `%s`.", package, install),
         call. = FALSE)
  }
  invisible(TRUE)
}

.require_python_module <- function(module, python = NULL)
{
  .require_svg_package("reticulate", "install.packages('reticulate')")
  if (!is.null(python)) reticulate::use_python(python, required = TRUE)
  ## check that the module is installed without importing it (importing may
  ## require compatibility shims that the helper functions apply first)
  util  <- reticulate::import("importlib.util", convert = TRUE)
  found <- !is.null(util$find_spec(module))
  if (!found) {
    stop(sprintf(paste0("Python module '%s' is required in the reticulate ",
                        "environment (e.g., `pip install %s`)."), module, module),
         call. = FALSE)
  }
  invisible(TRUE)
}

## Standard multi-gene output: one row per gene with method, native score,
## score type, raw p-value, and adjusted p-value.
.svg_table <- function(gene, method, score = NA_real_, score_type = "pvalue",
                       pvalue = NA_real_, padj = NA_real_)
{
  data.frame(gene = as.character(gene), method = method,
             score = as.numeric(score), score_type = score_type,
             pvalue = as.numeric(pvalue), padj = as.numeric(padj),
             stringsAsFactors = FALSE, row.names = NULL)
}
