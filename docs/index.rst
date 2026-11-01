.. boost-r documentation master file, created by
   sphinx-quickstart on Fri Apr  8 14:45:30 2022.

boost R Package
==========================

``boost`` is the companion R package of the `Spatial Transcriptomics Arena (STAr) <https://lce.biohpc.swmed.edu/star/>`_, described in "Spatial Transcriptomics Arena: An Integrated Platform for Methodology Research with a Focus on Spatially Variable Gene Detection" by Xi Jiang, Danni Luo, Huimin Li, Esteban Fernández, Kevin Lutz, Suhana Bedi, Jie Yang, Yuanchun Zhan, Bo Yao, Guanghua Xiao, Yang Xie, Qiwei Li, and Xiaowei Zhan.

The package provides the 12 spatially variable gene (SVG) detection methods benchmarked in STAr. BOOST-GP, BOOST-MI, BinSpect, SPARK, and SpatialDE are implemented natively; nnSVG, SPARK-X, MERINGUE, Trendsceek, Seurat Moran's I, SOMDE, and SpaGCN are available through wrappers that call the original software, which is installed as an optional dependency. GEO and MC2 are provided for rank aggregation, together with utilities for filtering, normalization, binarization, spatial neighbor construction, visualization, and performance evaluation.

The :doc:`tutorial` demonstrates the complete workflow on a real STAr dataset, and :doc:`simulation` shows how the STAr performance metrics are computed on simulated data with known ground truth.

.. figure:: website.png
    :align: center
    :figwidth: 700px

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   loaddata
   preprocessing
   sv_gene_detection
   new_methods
   tutorial
   simulation
