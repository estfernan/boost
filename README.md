boost
=====

<!-- [![Build Status][gha-icon]][gha-url] -->

The R package **boost** (Bayesian Modeling of Spatial Transcriptomics Data)
provides functions to detect spatially variable (SV) genes in 
spatial transcriptomics (ST) data. This package provides two novel Bayesian 
methods, as well as other standard procedures, for facilitating studies in 
spatial molecular profiling (SMP). For the step-by-step tutorial of **boost**, please refer to https://boost-r.readthedocs.io/en/latest/

## Development

The latest version of the package is under development at [GitHub][github-url],
which can be installed from GitHub by:

```R
devtools::install_github("estfernan/boost")
library(boost)
```

`SPARK` is a dependency and is necessary for using the *SPARK* model. 
It can be installed from [GitHub][spark-url] by:

```R
devtools::install_github('xzhoulab/SPARK')
```

`SpatialDE` is also necessary for using the *SpatialDE* model. 
It can be installed from [Bioconductor][spatialde-url] by:

```R
BiocManager::install('spatialDE')
```

## References

- Jiang, X., Li, Q., & Xiao, G. (2021).
  Bayesian Modeling of Spatial Transcriptomics Data via a Modified Ising Model. 
  *arXiv preprint arXiv:2104.13957*.

- Li, Q., Zhang, M., Xie, Y., & Xiao, G. (2020). 
  Bayesian Modeling of Spatial Molecular Profiling Data via Gaussian Process. 
  *arXiv preprint arXiv:2012.03326*.

## License

[GNU General Public License][gpl] (≥ 3)

[gha-icon]: https://github.com/estfernan/boost/workflows/R-CMD-check/badge.svg
[gha-url]: https://github.com/estfernan/boost/actions
[github-url]: https://github.com/estfernan/boost
[spark-url]: https://github.com/xzhoulab/SPARK
[spatialde-url]: https://bioconductor.org/packages/release/bioc/html/spatialDE.html
[gpl]: https://www.gnu.org/licenses/
