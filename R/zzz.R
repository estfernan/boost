##
## R package boost by Esteban Fernández, Xi Jiang, Suhana Bedi, and Qiwei Li
## Copyright (C) 2021
##
## This file is part of the R package boost.
##
## The R package boost is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package boost is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## This loading message should be temporary until SPARK can be installed
## through CRAN without any issues. Currently, SPARK can only be installed
## through GitHub. Although, SpatialDE needs to be installed through
## Bioconductor, so the message might need to stay.
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.onAttach <- function(libname, pkgname)
{
  packageStartupMessage(
    paste0(
      "boost: Bayesian modeling of spatial transcriptomics data\n",
      "\n",
      "To use the SPARK model, please install the source package by:\n",
      "\n",
      "  #> install.packages('devtools')\n",
      "  #> devtools::install_github('xzhoulab/SPARK')\n",
      "\n",
      "To use the SpatialDE model, please install the source package by:\n",
      "\n",
      "  #> install.packages('BiocManager')\n",
      "  #> BiocManager::install('spatialDE')\n",
      "\n",
      "More information about the package can be found in ",
      "<https://github.com/estfernan/boost>."
    )
  )
}
