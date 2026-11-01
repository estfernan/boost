Additional SV Gene Detection Methods
====================================

Besides the five methods implemented natively (BOOST-GP, BOOST-MI, BinSpect, SPARK, and SpatialDE; see :doc:`sv_gene_detection`), ``boost`` provides wrappers for the other seven SVG detection methods benchmarked in STAr. The wrappers call the original software, which is installed as an optional dependency (see :doc:`installation`), and follow the settings described in Supplementary Notes S1 of the STAr paper.

* **nnSVG**: nearest-neighbor Gaussian processes (R/Bioconductor).
* **SPARK-X**: scalable non-parametric covariance test (R package SPARK).
* **MERINGUE**: clustering-based differential expression or spatial autocorrelation (R).
* **Trendsceek**: marked point-process permutation tests (R).
* **Seurat Moran's I**: Moran's I spatial autocorrelation through Seurat (R).
* **SOMDE**: self-organizing map followed by SpatialDE on the nodes (Python, through ``reticulate``).
* **SpaGCN**: graph convolutional network-based spatial domains followed by domain-guided differential expression (Python, through ``reticulate``).

Common interface
----------------

All wrappers analyze all genes at once and take the same inputs as the rest of the package: a **spot-by-gene** count (or expression) matrix with gene names as column names, and a spot-by-2 coordinate matrix. They return a data frame with one row per gene and the columns

``gene | method | score | score_type | pvalue | padj``

which is the same format returned by ``as.svg.result()`` for the native methods. The native output of the original software is attached as ``attr(result, "native")``. SVGs are selected with ``svg.select()`` (adjusted p-value below 0.05)::

    library(boost)
    data(mob)

    flt   <- filter.st(mob, mob.spots)   # STAr default criterion
    count <- flt$count
    spots <- flt$spots

    svg.methods()

nnSVG
-----

Low-expressed and mitochondrial genes are removed with ``nnSVG::filter_genes()`` (default settings), log-normalized expression is computed with ``scuttle``, and nnSVG is run with its default settings (10 nearest neighbors)::

    res_nnsvg <- nnSVG(count, spots)

SPARK-X
-------

SPARK-X is fitted to the raw counts with ``option = "mixture"``; ``padj`` is the adjusted combined p-value reported by SPARK-X::

    res_sparkx <- SPARK.X(count, spots)

MERINGUE
--------

By default, MERINGUE follows the STAr implementation: spots and genes with fewer than 100 total counts are removed, counts are CPM-normalized and log-transformed, spots are clustered by PCA and graph-based clustering (:math:`k = 30`), and genes differentially expressed between clusters are identified by Wilcoxon rank-sum tests; ``padj`` is the minimum adjusted p-value across clusters. MERINGUE's spatial autocorrelation test is available with ``approach = "spatial"``::

    res_meringue <- MERINGUE(count, spots)
    res_meringue_spatial <- MERINGUE(count, spots, approach = "spatial")

Trendsceek
----------

Positions are converted to a point pattern, log10-transformed expression is set as marks, and permutation tests are run. The BH-adjusted p-values of the mark-variogram statistic, which performed best in the STAr simulations, are used by default (``statistic = "markvario"``). Trendsceek is computationally intensive, so start with a small gene set::

    res_trendsceek <- Trendsceek(count[, 1:100], spots, nrand = 1000)

Seurat Moran's I
----------------

Moran's I is computed with ``Seurat::RunMoransI()`` on normalized expression, and p-values are adjusted with the Benjamini-Hochberg procedure. In STAr, Moran's I serves as a reference spatial-autocorrelation measure::

    norm <- log1p(normalize.st(count, scaling.method = "TSS") * 1e4)
    res_moran <- Seurat.MoransI(norm, spots)

SOMDE
-----

A self-organizing map is trained on the coordinates with ``somn = 20``, node-level expression is generated and normalized, and SpatialDE is run on the nodes. ``padj`` is the SOMDE q-value::

    res_somde <- SOMDE(count, spots, k = 20)

SpaGCN
------

Genes expressed in fewer than three spots and special genes (ERCC spike-ins and mitochondrial genes) are removed, two spatial domains are identified and refined, and differential expression is tested for each domain against its neighboring domain; ``padj`` is the minimum adjusted p-value across the two domains. For 10x Visium data, supply the pixel coordinates and the histology image, and use ``shape = "hexagon"``::

    res_spagcn <- SpaGCN(count, spots)

    ## 10x Visium with histology
    ## res_spagcn <- SpaGCN(count, array_coords, pixels = pixel_coords,
    ##                      image = "tissue_hires_image.png", shape = "hexagon")

A custom SpaGCN workflow can be supplied with ``runner``, a function that receives the imported SpaGCN module, the count matrix, and the coordinates.

Combining results
-----------------

Because all methods return the same columns, their results can be stacked and compared::

    all_results <- rbind(res_nnsvg, res_sparkx, res_meringue, res_moran)
    all_results <- svg.select(all_results)
    tapply(all_results$is.svg, all_results$method, sum)

In STAr, GEO and MC2 rank aggregation combine BOOST-GP, BOOST-MI, BinSpect, SPARK, and SpatialDE (see :doc:`tutorial`); ``svg.aggregate()`` can also be applied to any other combination of methods.
