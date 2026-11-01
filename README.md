# boost

**boost** is the companion R package of the [Spatial Transcriptomics Arena (STAr)](https://lce.biohpc.swmed.edu/star/), an integrated platform for methodology research on spatially variable gene (SVG) detection. It provides the 12 SVG detection methods benchmarked in STAr, GEO and MC2 rank aggregation, and utilities for filtering, normalization, binarization, spatial neighbor construction, visualization, and performance evaluation.

Step-by-step tutorials are available at **<https://boost-r.readthedocs.io/en/latest/>**.

## Methods

| Method | Function | Implementation | Native statistic | SVG criterion in STAr |
|---|---|---|---|---|
| BOOST-GP | `BOOST.GP()` | native | Bayes factor | BF > 150 |
| BOOST-MI | `BOOST.MI()` | native | Bayes factor | BF > 150 |
| BinSpect | `binSpect()` | native | p-value | adjusted p < 0.05 |
| SPARK | `SPARK()` | native (calls SPARK) | p-value | adjusted p < 0.05 |
| SpatialDE | `SpatialDE()` | native (calls spatialDE) | p-value | adjusted p < 0.05 |
| nnSVG | `nnSVG()` | R wrapper | p-value | adjusted p < 0.05 |
| SPARK-X | `SPARK.X()` | R wrapper | p-value | adjusted p < 0.05 |
| MERINGUE | `MERINGUE()` | R wrapper | p-value | adjusted p < 0.05 |
| Trendsceek | `Trendsceek()` | R wrapper | p-value | adjusted p < 0.05 |
| Seurat Moran's I | `Seurat.MoransI()` | R wrapper | p-value | adjusted p < 0.05 |
| SOMDE | `SOMDE()` | Python wrapper (reticulate) | q-value | q < 0.05 |
| SpaGCN | `SpaGCN()` | Python wrapper (reticulate) | p-value | adjusted p < 0.05 |
| GEO | `svg.aggregate(..., method = "GEO")` | native | aggregated rank | top-ranked genes |
| MC2 | `svg.aggregate(..., method = "MC2")` | native | aggregated rank | top-ranked genes |

The wrappers call the original software, which is installed as an optional dependency and is not copied into boost; their settings follow Supplementary Notes S1 of the STAr paper. `svg.methods()` returns this table in R.

## Installation

```r
install.packages("remotes")
remotes::install_github("estfernan/boost")
library(boost)
```

SPARK and SpatialDE are required dependencies; if they are not installed automatically:

```r
remotes::install_github("xzhoulab/SPARK")
BiocManager::install("spatialDE")
```

The packages for the other methods are optional; see the [installation page](https://boost-r.readthedocs.io/en/latest/installation.html).

## Quick start

All functions use a spot-by-gene count matrix (spots in rows, genes in columns) and a spot-by-2 coordinate matrix.

```r
data(mob)                                   # mouse olfactory bulb (ST)

## method-specific filtering used in STAr
flt <- filter.st(mob, mob.spots, method = "SPARK")
sf  <- get.size.factor(flt$count, estimation.method = "TSS")

## single-gene methods, converted to a standard result table
res <- rbind(
  as.svg.result(SPARK(flt$count[, "Penk"], flt$spots, size.factor = sf, gene.name = "Penk")),
  as.svg.result(SPARK(flt$count[, "Nrf1"], flt$spots, size.factor = sf, gene.name = "Nrf1"))
)

## multi-gene wrappers return the same columns
res_sparkx <- SPARK.X(flt$count, flt$spots)

## SVG selection with the STAr criteria (BF > 150; BH-adjusted p < 0.05)
res_sparkx <- svg.select(res_sparkx)
sum(res_sparkx$is.svg)
```

The standardized output has the columns `gene | method | score | score_type | pvalue | padj`. Bayes factors and p-values are kept as they are and never converted into a common score. `svg.aggregate()` combines the rankings of several methods with GEO or MC2.

## References

- Jiang, X., Luo, D., Li, H., Fernández, E., Lutz, K., Bedi, S., Yang, J., Zhan, Y., Yao, B., Xiao, G., Xie, Y., Li, Q., & Zhan, X. Spatial Transcriptomics Arena: An Integrated Platform for Methodology Research with a Focus on Spatially Variable Gene Detection. Manuscript under revision.
- Jiang, X., Xiao, G., & Li, Q. (2022). A Bayesian modified Ising model for identifying spatially variable genes from spatial transcriptomics data. *Statistics in Medicine*, 41(23), 4647–4665.
- Li, Q., Zhang, M., Xie, Y., & Xiao, G. (2021). Bayesian modeling of spatial molecular profiling data via Gaussian process. *Bioinformatics*, 37(22), 4129–4136.

References for the other methods are given in the corresponding help pages.

## License

[GNU General Public License](https://www.gnu.org/licenses/) (≥ 3)
