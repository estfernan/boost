Installation
============

Install boost
-------------

The latest version of ``boost`` can be installed from GitHub::

    install.packages("remotes")
    remotes::install_github("estfernan/boost")

Then load the package::

    library(boost)

SPARK and SpatialDE are required dependencies. If they are not installed automatically, install them first::

    remotes::install_github("xzhoulab/SPARK")       # SPARK and SPARK-X

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("spatialDE")

Optional dependencies for the other benchmarked methods
--------------------------------------------------------

The wrappers for the other seven methods call the original software, which is not copied into ``boost``. Install only the packages for the methods you plan to run.

nnSVG::

    BiocManager::install(c("nnSVG", "SpatialExperiment", "scuttle"))

MERINGUE and Trendsceek::

    remotes::install_github("JEFworks-Lab/MERINGUE")
    remotes::install_github("edsgard/trendsceek")

trendsceek relies on the legacy, single-package version of spatstat (< 2.0); with newer spatstat versions it fails. To run Trendsceek, install a legacy spatstat version (e.g., ``remotes::install_version("spatstat", "1.64-1")``), preferably in a separate library.

Seurat Moran's I (``RunMoransI`` also needs ``Rfast2`` or ``ape``)::

    install.packages(c("Seurat", "ape"))

SOMDE and SpaGCN are Python packages, accessed from R through ``reticulate``::

    install.packages("reticulate")
    reticulate::py_install(c("somde", "SpaGCN", "scanpy"), pip = TRUE)

SOMDE was written for older SciPy and pandas versions; ``SOMDE()`` restores the few functions it needs (``scipy.misc.derivative``, NumPy functions in the ``scipy`` namespace, and ``pandas.Series.ravel``) when they are missing. SOMDE is designed for large datasets: its self-organizing map has about :math:`\sqrt{n/20}` nodes per side, so with a few hundred spots only a handful of nodes are used. The Louvain clustering used by SpaGCN requires ``setuptools`` < 80 (``pip install "setuptools<80"``) with recent Python versions.

To use a specific Python installation, pass ``python = "/path/to/python"`` to ``SOMDE()`` or ``SpaGCN()``, or call ``reticulate::use_python()`` before running them.

Check available methods
-----------------------

After installation, ``svg.methods()`` lists the SVG detection and rank-aggregation methods, their implementation, input, native statistic, and the SVG selection criterion used in STAr::

    svg.methods()

For reproducibility, record the R, Python, and package versions used in each analysis (e.g., ``sessionInfo()`` and ``reticulate::py_config()``).
