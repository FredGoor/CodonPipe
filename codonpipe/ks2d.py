"""codonpipe.ks2d

2D Kolmogorov–Smirnov (Fasano–Franceschini) with permutation p-values.
"""

# codonpipe/ks2d.py
import math
import itertools
import numpy as np
import pandas as pd

from scipy.ndimage import gaussian_filter


def _four_quadrant_cdfs(counts):
    total = counts.sum()
    if total <= 0:
        z = np.zeros_like(counts, dtype=float)
        return z, z, z, z
    F_ll = counts.cumsum(axis=0).cumsum(axis=1) / total
    F_ul = np.flipud(counts).cumsum(axis=0)
    F_ul = np.flipud(F_ul.cumsum(axis=1)) / total
    F_lu = np.fliplr(counts).cumsum(axis=1)
    F_lu = np.fliplr(F_lu.cumsum(axis=0)) / total
    F_uu = np.flipud(np.fliplr(counts)).cumsum(axis=0).cumsum(axis=1)
    F_uu = np.flipud(np.fliplr(F_uu)) / total
    return F_ll, F_lu, F_ul, F_uu


def _D_from_counts(countsA, countsB):
    FA = _four_quadrant_cdfs(countsA)
    FB = _four_quadrant_cdfs(countsB)
    D = 0.0
    for a, b in zip(FA, FB):
        d = float(np.max(np.abs(a - b)))
        if d > D:
            D = d
    return D


def _binned_stat(A, B, bins=151):
    X = np.vstack([A, B])
    nA = len(A)
    N = len(X)
    bx = np.linspace(X[:, 0].min(), X[:, 0].max(), bins + 1)
    by = np.linspace(X[:, 1].min(), X[:, 1].max(), bins + 1)
    ix = np.clip(np.searchsorted(bx, X[:, 0], side='right') - 1, 0, bins - 1)
    iy = np.clip(np.searchsorted(by, X[:, 1], side='right') - 1, 0, bins - 1)
    lin = ix * bins + iy
    counts_all = np.bincount(lin, minlength=bins * bins).reshape(bins, bins)
    maskA = np.zeros(N, dtype=bool)
    maskA[:nA] = True
    countsA = np.bincount(lin[maskA], minlength=bins * bins).reshape(bins, bins)
    countsB = counts_all - countsA
    D = _D_from_counts(countsA, countsB)
    cache = dict(lin=lin, N=N, nA=nA, bins=bins, counts_all=counts_all)
    return D, cache


def _binned_perm(cache, n_perm=1000, seed=None):
    rng = np.random.default_rng(seed)
    lin = cache['lin']
    N = cache['N']
    nA = cache['nA']
    bins = cache['bins']
    counts_all = cache['counts_all']
    for _ in range(int(n_perm)):
        idxA = rng.choice(N, size=nA, replace=False)
        countsA = np.bincount(lin[idxA], minlength=bins * bins).reshape(bins, bins)
        countsB = counts_all - countsA
        yield _D_from_counts(countsA, countsB)


def ks2d_ff(A, B, n_perm=2000, method='binned', bins=151, seed=None):
    A = np.asarray(A, float)
    B = np.asarray(B, float)
    assert A.ndim == 2 and A.shape[1] == 2 and B.ndim == 2 and B.shape[1] == 2

    if method == 'binned':
        Dobs, cache = _binned_stat(A, B, bins=bins)
        ge = 0
        for Dp in _binned_perm(cache, n_perm=n_perm, seed=seed):
            if Dp >= Dobs - 1e-12:
                ge += 1
        p = (ge + 1) / (n_perm + 1)
        return Dobs, p, dict(method='binned', n_perm=n_perm, bins=bins, seed=seed)

    try:
        from astroML.density_estimation import fasano_franceschini
    except Exception as e:
        raise ImportError("astroML not available for exact FF; use method='binned'.") from e

    Dobs, _ = fasano_franceschini(A, B)
    rng = np.random.default_rng(seed)
    X = np.vstack([A, B])
    N = len(X)
    nA = len(A)
    ge = 0
    for _ in range(int(n_perm)):
        idx = rng.choice(N, size=nA, replace=False)
        mask = np.zeros(N, dtype=bool)
        mask[idx] = True
        Dp, _ = fasano_franceschini(X[mask], X[~mask])
        if Dp >= Dobs - 1e-12:
            ge += 1
    p = (ge + 1) / (n_perm + 1)
    return Dobs, p, dict(method='exact', n_perm=n_perm, bins=None, seed=seed)


def adjust_pvalues_bh(pvals):
    p = np.asarray(pvals, float)
    mask = np.isfinite(p)
    m = mask.sum()
    if m == 0:
        return np.full_like(p, np.nan, dtype=float)

    order = np.argsort(np.where(mask, p, np.inf))
    ranked = p[order]
    adj_rev = np.minimum.accumulate((ranked[::-1] * m) / np.arange(1, m + 1)[::-1])
    adj = np.clip(adj_rev[::-1], 0, 1)

    out = np.full_like(p, np.nan, dtype=float)
    out[order] = adj
    return out


def _clean_text(s):
    if not isinstance(s, str):
        return s
    return s.strip().strip('"').strip("'")


def compute_2d_ks_for_clusters(Y, row_names, cluster_txt_df, ks_settings):
    """
    Compute pairwise 2D KS comparisons between cluster distributions in a 2D embedding.

    Returns a long-form table with ALL unordered cluster pairs:
      cluster_A, cluster_B, n_A, n_B, D_2D_KS, p_perm, p_adj_BH, method, bins
    """
    if Y is None or getattr(Y, "shape", (0,))[0] == 0:
        print("[WARN] 2D KS: no DR coordinates (Y) available; skipping.")
        return None

    if cluster_txt_df is None or getattr(cluster_txt_df, "empty", True):
        print("[WARN] 2D KS: cluster table is empty; skipping.")
        return None

    Y2 = np.asarray(Y, float)
    if Y2.ndim != 2 or Y2.shape[1] < 2:
        raise ValueError(f"2D KS requires at least 2 columns in Y; got shape {Y2.shape}")
    Y2 = Y2[:, :2]

    idx_map = {str(tag): i for i, tag in enumerate(row_names)}

    cluster_coords = {}
    for col in cluster_txt_df.columns:
        tags = [_clean_text(v) for v in cluster_txt_df[col].tolist() if isinstance(v, str)]
        tags = [t for t in tags if t and t.lower() != "nan"]
        tags_unique = sorted(set(tags))

        indices, missing = [], 0
        for t in tags_unique:
            idx = idx_map.get(t)
            if idx is None:
                missing += 1
            else:
                indices.append(idx)

        if missing > 0:
            print(f"[INFO] 2D KS: cluster '{col}' has {missing} locus_tag(s) not found in DR coordinates (ignored).")
        if len(indices) == 0:
            print(f"[WARN] 2D KS: cluster '{col}' has no genes with DR coordinates; skipping.")
            continue

        cluster_coords[col] = np.asarray(Y2[np.array(indices), :], float)

    names = sorted(cluster_coords.keys())
    if len(names) < 2:
        print("[WARN] 2D KS: need at least two non-empty clusters; skipping.")
        return None

    method = ks_settings.get("method", "binned")
    bins = ks_settings.get("bins", 151)
    n_perm = ks_settings.get("n_perm", 2000)
    seed = ks_settings.get("random_seed", None)

    results = []
    for a, b in itertools.combinations(names, 2):
        A = cluster_coords[a]
        B = cluster_coords[b]
        if A.shape[0] == 0 or B.shape[0] == 0:
            continue

        D, p, info = ks2d_ff(A, B, n_perm=n_perm, method=method, bins=bins, seed=seed)

        results.append(dict(
            cluster_A=a, cluster_B=b,
            n_A=int(A.shape[0]), n_B=int(B.shape[0]),
            D_2D_KS=float(D),
            p_perm=float(p),
            method=info.get("method", method),
            bins=info.get("bins", bins),
        ))

    if not results:
        print("[WARN] 2D KS: no pairwise results; skipping.")
        return None

    ks_df = pd.DataFrame(results)
    ks_df["p_adj_BH"] = adjust_pvalues_bh(ks_df["p_perm"].values)

    alpha = ks_settings.get("alpha", 0.01)
    sig_mask = (ks_df["p_adj_BH"] < alpha) & np.isfinite(ks_df["p_adj_BH"])
    total_sig = int(sig_mask.sum())
    print(f"[INFO] 2D KS: computed {len(ks_df)} pairwise comparisons across {len(names)} clusters "
          f"({total_sig} with p_adj_BH < {alpha}).")

    ks_df = ks_df.sort_values(["p_adj_BH", "D_2D_KS"], ascending=[True, False]).reset_index(drop=True)
    cols = ["cluster_A", "cluster_B", "n_A", "n_B", "D_2D_KS", "p_perm", "p_adj_BH", "method", "bins"]
    return ks_df[cols]
