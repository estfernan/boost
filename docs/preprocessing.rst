Data Pre-processing
=============================

Before SV gene detection, there are several steps to pre-process the data.

Data Filtration
------------------------------

Spots with fewer than ten total counts across all genes are excluded. Genes are then filtered on the retained spots. In STAr, the gene-filtering criterion specified in the original publication of each SVG detection method is followed; for methods without an explicit rule, genes expressed in fewer than 10% of spots are excluded (Section 3.4 and Supplementary Notes S1 of the STAr paper). ``svg.filter.criteria()`` lists these criteria::

    svg.filter.criteria()
    ##      method min.total min.percentage min.spots               source
    ## 1  BOOST-GP        10            0.1         0 original publication
    ## 2  BOOST-MI        10            0.2         0 original publication
    ## 3     SPARK        10            0.1         0 original publication
    ## 4 SpatialDE        10            0.1         0         STAr default
    ## 5  BinSpect        10            0.1         0         STAr default
    ## 6    SpaGCN        10            0.0         3 original publication
    ## 7   default        10            0.1         0         STAr default

``filter.st()`` applies the criteria of a given method, or user-specified thresholds. Here we use the default criterion (genes expressed in at least 10% of spots), which is also used for BOOST-GP, SPARK, SpatialDE, and BinSpect::

    filter_result <- filter.st(count, sample_info)
    count <- filter_result$count
    sample_info <- filter_result$spots

    dim(count)
    ## [1]   260 11274
    dim(sample_info)
    ## [1] 260   2

After filtration, the Mouse Olfactory Bulb data (replicate 11) has 260 spots and 11,274 genes. For BOOST-MI, genes with zero counts in more than 80% of spots are excluded instead::

    filter_mi <- filter.st(count_raw, sample_info_raw, method = "BOOST-MI")
    dim(filter_mi$count)
    ## [1]  260 9769

where ``count_raw`` and ``sample_info_raw`` are the unfiltered data from :doc:`loaddata`. For SpaGCN, which performs its own filtering, genes expressed in fewer than three spots are excluded.


Size Factor Estimation
-----------------------------

Size factors are inputs for some SV gene detection approaches (BOOST-GP and SPARK). ``get.size.factor()`` estimates them from the count matrix :math:`Y` with one of four estimation methods ("TSS", "Q75", "RLE" and "TMM"). The output is a vector of length :math:`n` with the estimated size factor of each spot. In STAr, total sum scaling (TSS) is used::

    size_factor <- get.size.factor(count, estimation.method = "TSS")
    size_factor <- size_factor / mean(size_factor)


Expression Counts Normalization
------------------------------------

Some methods need normalized gene expression levels as input. ``normalize.st()`` obtains normalized expression levels, including adjusting for the library size, stabilizing the variance, and log-transforming. Seven normalization methods are provided ("TSS", "Q75", "RLE", "TMM", "A-VST", "N-VST" and "log-VST"). The output has the same shape as the input count matrix :math:`Y`::

    normalized_count <- normalize.st(count, scaling.method = "TSS")

The "log-VST" option is the approximate Anscombe variance-stabilizing transformation followed by regression on the log total counts, described in Supplementary Notes S2 of the STAr paper; it is the input used for SpatialDE.


Method-specific preprocessing
-----------------------------

SV detection methods do not all use the same representation of the data. BOOST-GP and SPARK model raw counts with size factors; BOOST-MI uses TSS-normalized expression binarized by Gaussian mixture clustering (``binarize.st(..., cluster.method = "GMC")``); BinSpect uses binarized expression (k-means or rank-based binarization); and SpatialDE uses "log-VST"-normalized expression. The wrappers for the other methods perform their own transformations internally (see :doc:`new_methods`). For this reason, STAr preserves method-specific preprocessing rather than forcing every method to use one common transformed matrix. See :doc:`tutorial` for a complete example.
