boost
=====

The R package **boost** (Bayesian Modeling of Spatial Transcriptomics Data)
provides functions to detect spatially variable (SV) genes in 
spatial transcriptomics (ST) data. This package provides two novel Bayesian 
methods, as well as other standard procedures, for facilitating studies in 
spatial molecular profiling (SMP).

## Development

The latest version of the package is under development at [GitHub][github-url],
which can be installed from GitHub by:

```R
devtools::install_github("estfernan/boost")
library(boost)
```

`SPARK` is a dependency in this package and is necessary for using the 
*SPARK* model. It can be installed from [GitHub][spark-url] by:

```R
devtools::install_github('xzhoulab/SPARK')
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

[github-url]: https://github.com/estfernan/boost
[spark-url]: https://github.com/xzhoulab/SPARK
[gpl]: https://www.gnu.org/licenses/
