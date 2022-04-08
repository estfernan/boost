Load Data
==========================

The current version of boost requires two input data:

* The gene expression count matrix :math:`Y: n \times p` (:math:`n` - number of spots; :math:`p` - number of genes)
* The location information matrix :math:`T: n \times 2`. It includes the :math:`x` and :math:`y` coordinate for each sample point.

Both data should be stored in R matrix format. For gene expression count matrix :math:`Y`, column names should be gene names.

In the following, we conduct the analysis on Mouse Olfactory Bulb data as an example to show the functions of boost package.

Files originally retrieved from the `SMP-Gym website <https://lce.biohpc.swmed.edu/smp_gym/explorer.php>`_ .

