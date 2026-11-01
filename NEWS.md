# boost 0.4.0

Aligned with the revised STAr paper.

## New features

- Wrappers for the seven SVG detection methods benchmarked in STAr that are not implemented natively: `nnSVG()`, `SPARK.X()`, `MERINGUE()`, `Trendsceek()`, `Seurat.MoransI()`, `SOMDE()`, and `SpaGCN()`. The wrappers call the original software (optional dependencies, not copied into boost) with the settings described in Supplementary Notes S1. `SOMDE()` and `SpaGCN()` use `reticulate` and the helper module in `inst/python/star_svg.py`.
- `filter.st()` and `svg.filter.criteria()` apply the method-specific spot and gene filtering used in STAr (genes expressed in at least 10% of spots by default, at least 20% for BOOST-MI, and in at least three spots for SpaGCN).
- `svg.select()` applies the STAr SVG criteria (Bayes factor > 150 for BOOST-GP and BOOST-MI; Benjamini-Hochberg adjusted p-value < 0.05 otherwise).
- `svg.rank.table()` and `svg.aggregate()` build within-method rankings and aggregate them with GEO or MC2 (K = 1,000 in STAr).
- `svg.methods()` now lists the implementation, input, native statistic, and SVG criterion of each method.

## Changes

- All functions, including the wrappers, use a spot-by-gene matrix (spots in rows, genes in columns).
- `BOOST.MI()` / `BOOST.Ising()`: the default prior standard deviation of the interaction parameter is now 1/3 (previously 1), as used in STAr. The results also include `BF.max`, the larger of the attraction and repulsion Bayes factors reported for real datasets in STAr, and `print()` shows both Bayes factors.
- `as.svg.result()` reports `BF.max` for BOOST-MI by default (`bf = "attraction"` for the attraction pattern only), keeps Bayes factors without a p-value for the Bayesian methods, and returns the same columns as the wrappers (including `padj`).
- `aggregate.svg()` was renamed to `svg.aggregate()`, which now takes results in the standard format.
- References updated to the published BOOST-GP (Bioinformatics, 2021) and BOOST-MI (Statistics in Medicine, 2022) papers; `citation("boost")` now includes the STAr paper.
- Documentation and tutorials rewritten for the STAr workflow (data loading, filtering, method-specific inputs, SVG detection, standardization, selection, rank aggregation, and visualization).

# boost 0.3.0

- Added `BOOST.MI()` as the current-name alias for `BOOST.Ising()` while preserving backward compatibility.
- Added `svg.methods()` and `as.svg.result()`.

# boost 0.2.0

- Incorporated the *SpatialDE* method for detecting spatially variable genes.

# boost 0.1.0

- The first version of `boost` providing functions for fitting the novel Bayesian models, other standard techniques, and preprocessing tools for spatial transcriptomics (ST) data.
