Installation
==========================

The latest version of the boost package is under development at GitHub, and it can be installed from GitHub by::

    devtools::install_github("estfernan/boost")

SPARK is a dependency and is necessary for using the SPARK model. It can be installed from GitHub by::

    devtools::install_github('xzhoulab/SPARK')

SpatialDE is also necessary for using the SpatialDE model. It can be installed from Bioconductor by::

    BiocManager::install('spatialDE')

Import the package into the environment::

    library(boost)
