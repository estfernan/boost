"""Python helpers used by the boost wrappers for SOMDE and SpaGCN.

The functions follow the STAr implementation described in Supplementary
Notes S1 of the STAr paper. They are called from R through reticulate by
boost::SOMDE() and boost::SpaGCN(); the third-party packages themselves are
not copied into boost and must be installed in the Python environment.
"""

import random

import numpy as np
import pandas as pd


def _ensure_somde_compat():
    """Compatibility shims for running SOMDE with current SciPy/pandas.

    SOMDE imports scipy.misc.derivative (removed in SciPy 1.12), uses NumPy
    functions through the scipy namespace (e.g., scipy.arange, removed in
    SciPy 1.x), and calls pandas.Series.ravel (removed in pandas 3.0). If
    these are missing, the original behavior is restored so that SOMDE runs
    unchanged.
    """
    import scipy
    import scipy.misc

    for name in dir(np):
        if not name.startswith("_") and not hasattr(scipy, name):
            setattr(scipy, name, getattr(np, name))

    if not hasattr(pd.Series, "ravel"):
        pd.Series.ravel = lambda self, *args, **kwargs: \
            np.asarray(self).ravel(*args, **kwargs)

    if hasattr(scipy.misc, "derivative"):
        return

    def central_diff_weights(Np, ndiv=1):
        ho = Np >> 1
        x = np.arange(-ho, ho + 1.0)[:, np.newaxis]
        X = x ** 0.0
        for k in range(1, Np):
            X = np.hstack([X, x ** k])
        w = np.prod(np.arange(1, ndiv + 1), axis=0) * np.linalg.inv(X)[ndiv]
        return w

    def derivative(func, x0, dx=1.0, n=1, args=(), order=3):
        weights = central_diff_weights(order, n)
        val = 0.0
        ho = order >> 1
        for k in range(order):
            val += weights[k] * func(x0 + (k - ho) * dx, *args)
        return val / np.prod((dx,) * n, axis=0)

    scipy.misc.derivative = derivative


def run_somde(counts, genes, spots, coords, k=20, alpha=0.5):
    """Run SOMDE (Hao et al., 2021).

    counts : spots x genes array of raw counts
    coords : spots x 2 array of spatial coordinates
    k      : SOM size parameter ('somn'), 20 in STAr
    """
    _ensure_somde_compat()
    from somde import SomNode

    df = pd.DataFrame(np.asarray(counts, dtype=float).T,
                      index=list(genes), columns=list(spots))
    X = np.asarray(coords, dtype=np.float32)

    som = SomNode(X, int(k))
    som.mtx(df, alpha=alpha)
    som.norm()
    result, _ = som.run()
    return result.reset_index(drop=True)


def run_spagcn_svg(counts, genes, spots, x_array, y_array,
                   x_pixel=None, y_pixel=None, image=None,
                   histology=False, n_clusters=2, p=0.5,
                   shape="square", seed=100, max_epochs=200,
                   min_cells=3):
    """Detect SVGs with SpaGCN (Hu et al., 2021) as in STAr.

    Genes expressed in fewer than `min_cells` spots and special genes
    (ERCC spike-ins and mitochondrial genes) are removed, spatial domains are
    identified with `n_clusters` clusters (two for global SVGs) and refined,
    and domain-guided differential expression is performed for each domain
    against its neighboring domain(s). The minimum adjusted p-value across
    domains is returned for each gene.
    """
    import anndata as ad
    import scanpy as sc
    import torch
    import SpaGCN as spg

    adata = ad.AnnData(X=np.asarray(counts, dtype=np.float32))
    adata.var_names = [str(g) for g in genes]
    adata.obs_names = [str(s) for s in spots]
    x_array = [float(v) for v in x_array]
    y_array = [float(v) for v in y_array]
    adata.obs["x_array"] = x_array
    adata.obs["y_array"] = y_array

    if x_pixel is None or y_pixel is None:
        x_pixel, y_pixel = x_array, y_array
    else:
        x_pixel = [int(v) for v in x_pixel]
        y_pixel = [int(v) for v in y_pixel]

    img = None
    if histology:
        if image is None:
            raise ValueError("histology=True requires an image")
        import cv2
        img = cv2.imread(image) if isinstance(image, str) else np.asarray(image)

    adj = spg.calculate_adj_matrix(x=x_pixel, y=y_pixel, x_pixel=x_pixel,
                                   y_pixel=y_pixel, image=img, beta=49,
                                   alpha=1, histology=histology)

    spg.prefilter_genes(adata, min_cells=int(min_cells))
    spg.prefilter_specialgenes(adata)
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)

    l = spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)
    res = spg.search_res(adata, adj, l, int(n_clusters), start=0.7, step=0.1,
                         tol=5e-3, lr=0.05, max_epochs=20, r_seed=seed,
                         t_seed=seed, n_seed=seed)

    clf = spg.SpaGCN()
    clf.set_l(l)
    random.seed(seed)
    torch.manual_seed(seed)
    np.random.seed(seed)
    clf.train(adata, adj, init_spa=True, init="louvain", res=res, tol=5e-3,
              lr=0.05, max_epochs=int(max_epochs))
    y_pred, _ = clf.predict()

    adj_2d = spg.calculate_adj_matrix(x=x_array, y=y_array, histology=False)
    refined = spg.refine(sample_id=adata.obs.index.tolist(),
                         pred=[int(v) for v in y_pred], dis=adj_2d,
                         shape=shape)
    refined = [int(v) for v in refined]
    adata.obs["pred"] = refined

    domains = sorted(set(refined))
    out = pd.DataFrame({"gene": list(adata.var_names)})
    nz = adj_2d[adj_2d != 0]
    start, end = np.quantile(nz, q=0.001), np.quantile(nz, q=0.1)

    for d in domains:
        others = [o for o in domains if o != d]
        if not others:
            continue
        try:
            r = spg.search_radius(target_cluster=d,
                                  cell_id=adata.obs.index.tolist(),
                                  x=x_array, y=y_array, pred=refined,
                                  start=start, end=end, num_min=10,
                                  num_max=14, max_run=100)
            nbr = spg.find_neighbor_clusters(target_cluster=d,
                                             cell_id=adata.obs.index.tolist(),
                                             x=x_array, y=y_array,
                                             pred=refined, radius=r,
                                             ratio=1 / 2)
        except Exception:
            nbr = None
        if not nbr:
            nbr = others
        nbr = [int(v) for v in list(nbr)[0:3]]
        de = spg.rank_genes_groups(input_adata=adata, target_cluster=d,
                                   nbr_list=nbr, label_col="pred",
                                   adj_nbr=True, log=True)
        de = de[["genes", "pvals_adj"]].rename(
            columns={"genes": "gene", "pvals_adj": "padj_domain%d" % d})
        out = out.merge(de, on="gene", how="left")

    cols = [c for c in out.columns if c.startswith("padj_domain")]
    out["padj"] = out[cols].min(axis=1) if cols else np.nan
    return {"result": out, "domains": refined}
