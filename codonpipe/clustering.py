"""codonpipe.clustering

Core computations for codon/AA usage embedding and clustering.
"""

# codonpipe/clustering.py
import os
import re
import math
import itertools
import textwrap
from datetime import datetime
from tkinter import Tk
from tkinter.filedialog import askopenfilename

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, optimal_leaf_ordering, leaves_list
from scipy.ndimage import gaussian_filter1d, gaussian_filter
from scipy.stats import rankdata, ttest_ind, mannwhitneyu

from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.cluster import SpectralClustering, DBSCAN, KMeans

from matplotlib.colors import ListedColormap

from .colormaps import get_cmap_any
from .excel_outputs import (
    compute_acu_from_counts_df,
    compute_rcu_from_counts_df,
    compute_rcu_devz,
    compute_relative_usage_from_counts_df,
    compute_devz,
    read_trna_decoding_table,
    compute_trna_counts_from_codon_counts_df,
    compute_rtu_from_counts_df,
    compute_ztu,
)

try:
    from sklearn_extra.cluster import KMedoids
except Exception as e:
    print("[WARN] sklearn-extra KMedoids unavailable; falling back to KMeans:", e)
    KMedoids = None

try:
    import umap
except ImportError:
    umap = None


def _show_figure_nonblocking(fig):
    """Show a Matplotlib figure when the backend supports it.

    Do not return early for Agg-like backends: in Spyder, depending on the
    selected graphics backend, ``plt.show()`` can still be needed to push the
    saved figure to the Plots pane. Calling it is harmless when the backend
    cannot display interactively.
    """
    try:
        if fig is not None and getattr(fig, "canvas", None) is not None:
            try:
                fig.canvas.draw_idle()
            except Exception:
                pass
        try:
            plt.show(block=False)
        except TypeError:
            plt.show()
    except Exception:
        pass


def wrap_label_no_break(text, width=12):
    s = str(text)
    if width is None or int(width) <= 0:
        return s
    lines = textwrap.wrap(s, width=int(width), break_long_words=False, break_on_hyphens=False)
    return "\n".join(lines) if lines else s



# -----------------------------
# Small file-selection helpers
# -----------------------------
def choose_excel(initialdir, title):
    """Open a file dialog to select an Excel workbook.

    Falls back to a console prompt if a GUI is not available.
    """
    try:
        root = Tk()
        root.withdraw()
        path = askopenfilename(
            title=title,
            initialdir=initialdir or None,
            filetypes=[("Excel files", "*.xlsx")]
        )
        root.destroy()
    except Exception:
        path = input(f"{title} (.xlsx): ").strip().strip('"').strip("'")
    if not path:
        raise SystemExit("No file selected.")
    return path


def choose_cluster_file(initialdir):
    """Open a file dialog to select a locus-tag cluster file.

    Falls back to a console prompt if a GUI is not available.
    """
    try:
        root = Tk()
        root.withdraw()
        path = askopenfilename(
            title="Select locus-tag cluster file (TXT/CSV/TSV/Excel)",
            initialdir=initialdir or None,
            filetypes=[
                ("Cluster files", "*.txt *.csv *.tsv *.xlsx *.xls *.xlsm"),
                ("Text files", "*.txt *.csv *.tsv"),
                ("Excel files", "*.xlsx *.xls *.xlsm"),
                ("All files", "*.*"),
            ]
        )
        root.destroy()
    except Exception:
        path = input("Cluster file path: ").strip().strip('"').strip("'")
    if not path:
        raise SystemExit("No cluster file selected.")
    return path


def _clean_text(s):
    if not isinstance(s, str):
        return s
    return s.strip().strip('"').strip("'")


def _normalize_colname(c):
    return str(c).strip().lower().replace(" ", "").replace("-", "").replace(".", "").replace("_", "")


def _is_nonempty_text(x):
    s = _clean_text(x) if isinstance(x, str) else str(x).strip() if x is not None else ""
    return bool(s) and s.lower() not in {"nan", "na", "none"}


def _make_unique_preferred_ids(preferred_ids, fallback_ids):
    out = []
    seen = set()
    for pref, fb in zip(list(preferred_ids), list(fallback_ids)):
        pref_s = _clean_text(pref) if isinstance(pref, str) else str(pref).strip() if pref is not None else ""
        fb_s = _clean_text(fb) if isinstance(fb, str) else str(fb).strip() if fb is not None else ""

        key = pref_s if _is_nonempty_text(pref_s) else fb_s
        if (not key) or (key in seen):
            base = fb_s if _is_nonempty_text(fb_s) else (pref_s if _is_nonempty_text(pref_s) else "gene")
            key = base
            i = 2
            while (not key) or (key in seen):
                key = f"{base}__dup{i}"
                i += 1

        seen.add(key)
        out.append(key)
    return out


def _prepare_export_geneids_df(geneids_df):
    if geneids_df is None:
        return None

    df = geneids_df.copy()
    if df.empty:
        return df

    if "RefSeq_LocusTag_RS" in df.columns:
        refseq = df["RefSeq_LocusTag_RS"].astype(str).map(_clean_text)
        if "LocusTag" in df.columns:
            locus_old = df["LocusTag"].astype(str).map(_clean_text)
        else:
            locus_old = pd.Series([""] * len(df), index=df.index, dtype=object)
        df["LocusTag"] = [r if _is_nonempty_text(r) else l for r, l in zip(refseq.tolist(), locus_old.tolist())]

    for col in ["Old_LocusTag", "PrimaryID"]:
        if col in df.columns:
            df = df.drop(columns=[col])

    preferred = [
        "LocusTag", "GeneSymbol", "EntrezGeneID", "ProteinDescription",
        "RefSeqProteinID", "UniProtID", "RefSeq_LocusTag_RS"
    ]
    ordered = [c for c in preferred if c in df.columns] + [c for c in df.columns if c not in preferred]
    return df.loc[:, ordered]


def _prepare_export_codon_df(codon_df, geneids_df):
    df = codon_df.copy()
    if df.empty:
        return df

    fallback_ids = [str(x).strip() for x in df.index.tolist()]
    if geneids_df is not None and (not geneids_df.empty) and ("RefSeq_LocusTag_RS" in geneids_df.columns):
        preferred_ids = geneids_df["RefSeq_LocusTag_RS"].astype(str).map(_clean_text).tolist()
        if len(preferred_ids) != len(fallback_ids):
            preferred_ids = fallback_ids
    else:
        preferred_ids = fallback_ids

    df.index = _make_unique_preferred_ids(preferred_ids, fallback_ids)
    return df


def infer_prefix_from_codon_basename(base_name):
    s = str(base_name)
    s = re.sub(r'(?i)codonusage', '', s)
    s = re.sub(r'(?i)codons?', '', s)
    s = re.sub(r'__+', '_', s)
    s = s.strip('_').strip()
    return s



def normalize_usage_basis(value):
    s = str(value or "").strip().upper()
    if s in ("RCU", "ACU", "AA"):
        return s
    raise ValueError(f"Unknown usage_basis: {value}")


def auto_codon_set_for_usage_basis(usage_basis):
    ub = normalize_usage_basis(usage_basis)
    return {"AA": "64_withSTOP", "ACU": "61_noSTOP", "RCU": "59_noSTOP_MW"}[ub]

def normalize_codon_set(value):
    s = str(value or "").strip().lower().replace(" ", "").replace("-", "_")
    s = s.replace("withstop", "with_stop").replace("nostop", "no_stop")
    aliases = {
        "59": "59_noSTOP_MW",
        "59_no_stop_mw": "59_noSTOP_MW",
        "59_no_stop": "59_noSTOP_MW",
        "59_nostop_mw": "59_noSTOP_MW",
        "59_nostop": "59_noSTOP_MW",
        "61": "61_noSTOP",
        "61_no_stop": "61_noSTOP",
        "61_nostop": "61_noSTOP",
        "64": "64_withSTOP",
        "64_all": "64_withSTOP",
        "64_with_stop": "64_withSTOP",
        "64_withstop": "64_withSTOP",
    }
    if s in aliases:
        return aliases[s]
    raise ValueError(f"Unknown codon_set: {value}")

def normalize_dimred_method(value):
    s = str(value or "").strip().lower().replace("-", "").replace("_", "")
    aliases = {
        "umap": "umap",
        "tsne": "tsne",
        "t-sne": "tsne",
        "pca": "pca",
        "none": "none",
    }
    if s in aliases:
        return aliases[s]
    raise ValueError(f"Unknown dimred_method: {value}")

def normalize_cluster_method(value):
    s = str(value or "").strip().lower().replace("-", "").replace("_", "")
    aliases = {
        "hierarchical": "hierarchical",
        "kmedoids": "kmedoids",
        "kmeans": "kmeans",
        "kmean": "kmeans",
        "dbscan": "dbscan",
        "spectral": "spectral",
    }
    if s in aliases:
        return aliases[s]
    raise ValueError(f"Unknown cluster_method: {value}")

def codon_set_kind(value):
    return normalize_codon_set(value).lower()

def auto_find_geneids_excel(codon_file):
    folder = os.path.dirname(codon_file)
    base = os.path.splitext(os.path.basename(codon_file))[0]

    for pat in (r'(?i)codons', r'(?i)codon'):
        new_base = re.sub(pat, 'geneIDs', base)
        if new_base != base:
            candidate = os.path.join(folder, new_base + ".xlsx")
            if os.path.exists(candidate):
                return candidate

    for suffix in ("_geneIDs", "_geneids", " geneIDs", " geneids"):
        candidate = os.path.join(folder, base + suffix + ".xlsx")
        if os.path.exists(candidate):
            return candidate

    hint = infer_prefix_from_codon_basename(base).lower().strip('_')
    candidates = []
    for fname in os.listdir(folder):
        if not fname.lower().endswith(".xlsx"):
            continue
        if os.path.abspath(os.path.join(folder, fname)) == os.path.abspath(codon_file):
            continue
        lname = fname.lower()
        if ("geneid" in lname) or ("gene_ids" in lname) or ("geneids" in lname):
            score = 0
            if hint and hint in lname:
                score -= 10
            score -= len(os.path.commonprefix([hint, lname]))
            candidates.append((score, fname))

    if candidates:
        candidates.sort()
        return os.path.join(folder, candidates[0][1])

    return None


def _infer_locus_tag_series(df, cols_l):
    cols = list(df.columns)
    locus_col = None
    for i, c in enumerate(cols_l):
        key = _normalize_colname(c)
        if key in ("locustag", "locustags", "locus"):
            locus_col = cols[i]
            break
    if locus_col is not None:
        return df[locus_col].astype(str).map(_clean_text)

    try:
        idx = pd.Series(df.index.astype(str))
        frac_like = (idx.str.contains(r"[A-Za-z]", regex=True)).mean()
        if frac_like > 0.5 and not idx.str.fullmatch(r"\d+").all():
            return idx.map(_clean_text)
    except Exception:
        pass

    if len(cols) > 0:
        return df[cols[0]].astype(str).map(_clean_text)

    return pd.Series([], dtype=str)


def load_geneids_maps(gene_file):
    if not gene_file or (not os.path.exists(gene_file)):
        print("[WARN] GeneIDs Excel not found; GeneSymbol/Description maps will be empty.")
        return {}, {}

    try:
        df = pd.read_excel(gene_file, dtype=str)
    except Exception:
        df = pd.read_excel(gene_file)

    if df is None or df.empty:
        print("[WARN] GeneIDs Excel is empty; GeneSymbol/Description maps will be empty.")
        return {}, {}

    cols = list(df.columns)
    cols_l = [str(c).strip().lower() for c in cols]
    locus_series = _infer_locus_tag_series(df, cols_l)

    gene_col = None
    for i, c in enumerate(cols_l):
        if _normalize_colname(c) == "genesymbol":
            gene_col = cols[i]
            break

    desc_col = None
    preferred = ["description", "product", "annotation", "function", "geneproduct", "genedescription"]
    for pref in preferred:
        for i, c in enumerate(cols_l):
            key = _normalize_colname(c)
            if key == pref:
                desc_col = cols[i]
                break
        if desc_col is not None:
            break

    if desc_col is None:
        for i, c in enumerate(cols_l):
            key = _normalize_colname(c)
            if ("desc" in key) or ("product" in key) or ("annot" in key) or ("function" in key):
                desc_col = cols[i]
                break

    gene_series = df[gene_col].astype(str).map(_clean_text) if gene_col is not None else pd.Series([""] * len(df), dtype=str)
    desc_series = df[desc_col].astype(str).map(_clean_text) if desc_col is not None else pd.Series([""] * len(df), dtype=str)

    gene_map, desc_map = {}, {}
    for lt, gs, ds in zip(locus_series.tolist(), gene_series.tolist(), desc_series.tolist()):
        if not isinstance(lt, str):
            continue
        lt = lt.strip()
        if lt == "" or lt.lower() == "nan":
            continue
        gs = "" if (not isinstance(gs, str) or gs.lower() == "nan") else gs.strip()
        ds = "" if (not isinstance(ds, str) or ds.lower() == "nan") else ds.strip()

        if lt not in gene_map or ((not gene_map[lt]) and gs):
            gene_map[lt] = gs
        if lt not in desc_map or ((not desc_map[lt]) and ds):
            desc_map[lt] = ds

    return gene_map, desc_map


# -----------------------------
# Codon usage tables
# -----------------------------
def compute_AA_ACU_RCU(count_df):
    """
    Build harmonized gene-level AA, ACU, and RCU tables from raw codon counts.

    count_df must contain per-gene codon counts. ACU and RCU are computed
    using the corrected logic used for the exported codon-usage workbooks.
    """
    df = count_df.copy()
    var_names = [str(v).replace('END', 'STOP') for v in df.columns.tolist()]
    df.columns = var_names
    row_names = df.index.to_numpy()

    aa_names = np.array([v.split('_')[0] for v in var_names], dtype=object)
    unique20 = [a for a in pd.unique(aa_names) if a != 'STOP']

    counts = df.to_numpy(dtype=float)
    AA = np.zeros((counts.shape[0], len(unique20)), dtype=float)
    for i, aa in enumerate(unique20):
        cols = (aa_names == aa)
        AA[:, i] = counts[:, cols].sum(axis=1)
    AA_df = pd.DataFrame(AA, index=row_names, columns=unique20)

    acu_df = compute_acu_from_counts_df(df)
    rcu_df = compute_rcu_from_counts_df(df)

    return acu_df, acu_df.shape[1], AA_df, rcu_df, aa_names, unique20, row_names, var_names


def subset_usage(SET, usage_basis, C_abs, C_rel_df, AA_df, aa_names, var_names):
    n_codons_total = C_abs.shape[1]
    mask64 = np.ones(n_codons_total, dtype=bool)
    mask61 = (aa_names != 'STOP')
    mask59 = mask61 & (aa_names != 'Met') & (aa_names != 'Trp')

    codon_set = codon_set_kind(SET['codon_set'])
    if codon_set == '64_withstop':
        feat_mask = mask64
    elif codon_set == '61_nostop':
        feat_mask = mask61
    elif codon_set == '59_nostop_mw':
        feat_mask = mask59
    else:
        raise ValueError(f"Unknown codon_set: {SET['codon_set']}")

    usage_basis = normalize_usage_basis(usage_basis)
    if usage_basis == 'AA':
        Usage = AA_df.to_numpy(dtype=float)
        feature_labels = AA_df.columns.to_numpy()
        isAA = True
    elif usage_basis == 'ACU':
        Usage = C_abs
        feature_labels = np.array(var_names, dtype=object)
        isAA = False
    elif usage_basis == 'RCU':
        Usage = C_rel_df.to_numpy(dtype=float)
        feature_labels = np.array(var_names, dtype=object)
        isAA = False
    else:
        raise ValueError(f"Unknown usage_basis: {usage_basis}")

    if not isAA:
        Usage = Usage[:, feat_mask]
        feature_labels = feature_labels[feat_mask]

    return Usage, feature_labels, isAA


def _impute_missing_usage_for_embedding(Usage):
    """
    Replace NaN/inf entries before dimensionality reduction.

    RCU-style tables can legitimately contain NaN when a multi-codon amino-acid
    family is absent from a gene. Those NaNs are correct for exported per-gene
    tables and cluster means, but UMAP/t-SNE/PCA cannot consume them.

    For embedding/clustering only, impute missing values column-wise with the
    whole-genome mean of the corresponding feature. This keeps undefined values
    neutral with respect to genome-centered analyses, while preserving the raw
    exported tables unchanged.
    """
    values = np.asarray(Usage, dtype=float).copy()
    bad = ~np.isfinite(values)
    if not bad.any():
        return values

    col_means = np.nanmean(np.where(np.isfinite(values), values, np.nan), axis=0)
    col_means = np.where(np.isfinite(col_means), col_means, 0.0)
    rows, cols = np.where(bad)
    values[rows, cols] = col_means[cols]
    return values


def normalize_values(SET, Usage):
    values = _impute_missing_usage_for_embedding(Usage)
    if SET['center_features']:
        values -= values.mean(axis=0, keepdims=True)
    if SET['scale_features']:
        sd = values.std(axis=0, ddof=0, keepdims=True)
        sd[sd == 0] = 1.0
        values /= sd
    return values


# -----------------------------
# DimRed
# -----------------------------
def run_dimred(SET, values, feature_labels):
    method = normalize_dimred_method(SET['dimred_method'])

    if method == 'umap':
        if umap is None:
            raise ImportError("umap-learn not installed but dimred_method='umap'.")
        X_umap = values.copy()
        clip = SET.get('umap_clip_abs', 0.0)
        if clip > 0:
            X_umap = np.clip(X_umap, -clip, clip)
        random_state = None if SET.get('umap_randomize', False) else 0
        umap_init = str(SET.get('umap_init', 'spectral')).strip().lower()
        reducer = umap.UMAP(
            n_neighbors=SET['umap_neighbors'],
            min_dist=SET['umap_min_dist'],
            n_components=SET['umap_components'],
            metric=SET['umap_metric'],
            random_state=random_state,
            init=umap_init,
        )
        return reducer.fit_transform(X_umap)

    if method == 'tsne':
        tsne = TSNE(
            n_components=SET['tsne_dims'],
            perplexity=SET['tsne_perplexity'],
            metric=SET['tsne_distance'],
            learning_rate=SET['tsne_learnrate'],
            early_exaggeration=SET['tsne_exaggeration'],
            init='random',
            random_state=0,
        )
        return tsne.fit_transform(values)

    if method == 'pca':
        n_features = values.shape[1]
        k = max(2, min(SET['pca_npcs'], n_features))
        Xpca = values.copy()
        if SET['pca_center'] and not SET['center_features']:
            Xpca -= Xpca.mean(axis=0, keepdims=True)
        if SET['pca_scale'] and not SET['scale_features']:
            s = Xpca.std(axis=0, ddof=0, keepdims=True)
            s[s == 0] = 1.0
            Xpca /= s
        pca = PCA(n_components=k, svd_solver='auto', random_state=0)
        return pca.fit_transform(Xpca)

    if method == 'none':
        return None

    raise ValueError(f"Unknown dimred_method: {SET['dimred_method']}")


# -----------------------------
# Clustering
# -----------------------------
def _order_members_fast_or_hier(Z, members, dist_metric, SET):
    members = np.asarray(members, dtype=int)
    m = len(members)
    if m <= 2:
        return members

    # For very large clusters, avoid O(n^2)-memory / O(n^3)-time exact ordering.
    fast_threshold = int(SET.get('cluster_fast_order_threshold', 800))
    olo_threshold = int(SET.get('cluster_optimal_leaf_max_size', 250))
    use_olo = bool(SET.get('cluster_use_optimal_leaf_ordering', True))

    if m > fast_threshold:
        return members[np.argsort(Z[members, 0], kind='stable')]

    Dloc = pdist(Z[members, :], metric=dist_metric)
    tree = linkage(Dloc, method='single')
    if use_olo and m <= olo_threshold:
        tree = optimal_leaf_ordering(tree, Dloc)
    return members[leaves_list(tree)]


def _order_by_cluster_then_hier(Z, labels, dist_metric, SET=None):
    if SET is None:
        SET = {}
    unique_labels = np.unique(labels)
    k = len(unique_labels)
    centroids = np.zeros((k, Z.shape[1]), dtype=float)
    for i, c in enumerate(unique_labels):
        centroids[i, :] = Z[labels == c, :].mean(axis=0)
    seq = np.argsort(centroids[:, 0])

    order = []
    for i in seq:
        c = unique_labels[i]
        members = np.where(labels == c)[0]
        sub_order = _order_members_fast_or_hier(Z, members, dist_metric, SET)
        order.extend(sub_order.tolist())
    return np.array(order, dtype=int), labels


def cluster_genes(SET, Z):
    n = Z.shape[0]
    method = normalize_cluster_method(SET['cluster_method'])

    if method == 'hierarchical':
        hier_fast_threshold = int(SET.get('hierarchical_fast_order_threshold', 2000))
        use_olo = bool(SET.get('cluster_use_optimal_leaf_ordering', True))
        olo_threshold = int(SET.get('hierarchical_optimal_leaf_max_size', 400))
        if n > hier_fast_threshold:
            order = np.argsort(Z[:, 0], kind='stable')
        else:
            Dg = pdist(Z, metric=SET['gene_dist_metric'])
            tree = linkage(Dg, method=SET['gene_linkage'])
            if use_olo and n <= olo_threshold:
                tree = optimal_leaf_ordering(tree, Dg)
            order = leaves_list(tree)
        labels = np.zeros(n, dtype=int)
        return order, labels

    if method == 'kmedoids':
        k = max(2, min(int(SET.get('kmedoids_k', 12)), n - 1))
        if KMedoids is not None:
            model = KMedoids(
                n_clusters=k,
                metric=SET['kmedoids_dist'],
                random_state=0,
                init='k-medoids++',
                max_iter=300
            )
            labels = model.fit_predict(Z)
        else:
            print("[WARN] KMedoids unavailable; falling back to KMeans for 'kmedoids'.")
            km = KMeans(n_clusters=k, random_state=0, n_init='auto')
            labels = km.fit_predict(Z)
        return _order_by_cluster_then_hier(Z, labels, SET['gene_dist_metric'], SET)

    if method == 'kmeans':
        k = max(2, min(int(SET.get('kmeans_k', SET.get('kmedoids_k', 12))), n - 1))
        km = KMeans(n_clusters=k, random_state=0, n_init='auto')
        labels = km.fit_predict(Z)
        return _order_by_cluster_then_hier(Z, labels, SET['gene_dist_metric'], SET)

    if method == 'spectral':
        k = max(2, min(SET['spectral_k'], n - 1))
        spec = SpectralClustering(n_clusters=k, affinity='nearest_neighbors', random_state=0)
        labels = spec.fit_predict(Z)
        return _order_by_cluster_then_hier(Z, labels, SET['gene_dist_metric'], SET)

    if method == 'dbscan':
        db = DBSCAN(eps=SET['dbscan_eps'], min_samples=SET['dbscan_minpts'], metric=SET['dbscan_dist'])
        labels = db.fit_predict(Z)

        unique_labels = sorted(l for l in np.unique(labels) if l != -1)
        order = []
        for c in unique_labels:
            members = np.where(labels == c)[0]
            if len(members) > 2:
                Dloc = pdist(Z[members, :], metric=SET['gene_dist_metric'])
                tree = linkage(Dloc, method='single')
                tree_opt = optimal_leaf_ordering(tree, Dloc)
                sub_order = members[leaves_list(tree_opt)]
            else:
                sub_order = members
            order.extend(sub_order.tolist())

        noise = np.where(labels == -1)[0]
        if len(noise) > 0:
            idx = np.argsort(Z[noise, 0])
            order.extend(noise[idx].tolist())
        return np.array(order, dtype=int), labels

    raise ValueError(f"Unknown cluster_method: {SET['cluster_method']}")


# -----------------------------
# Feature reordering / smoothing
# -----------------------------
def reorder_features(SET, values, feature_labels):
    metric = SET['feature_dist_metric'].lower()
    if metric == 'spearman':
        ranks = np.apply_along_axis(rankdata, 1, values.T)
        Df = pdist(ranks, metric='correlation')
    else:
        Df = pdist(values.T, metric=metric)

    # Short CDS windows (for example codons 1-20) can make some codon features
    # constant across all genes. Correlation/Spearman distances are undefined for
    # constant vectors, producing NaN distances and breaking feature ordering.
    # Treat those undefined distances as neutral/large finite distances so the
    # clustering run still completes while preserving all codon features.
    Df = np.asarray(Df, dtype=float)
    if Df.size and not np.isfinite(Df).all():
        finite = Df[np.isfinite(Df)]
        fill_value = float(np.nanmax(finite)) if finite.size else 0.0
        Df = np.where(np.isfinite(Df), Df, fill_value)

    tree = linkage(Df, method=SET['feature_linkage'])
    tree_opt = optimal_leaf_ordering(tree, Df)
    order = leaves_list(tree_opt)
    return feature_labels[order], order


def smooth_and_bin(SET, V):
    if SET['apply_smoothing']:
        window = SET['smooth_window_genes']
        sigma = max(window / 3.0, 1.0)
        V = gaussian_filter1d(V, sigma=sigma, axis=1, mode='nearest')

    bin_size = max(1, int(round(SET['bin_size_genes'])))
    if SET['apply_binning']:
        ncols = V.shape[1]
        nbins = int(np.ceil(ncols / bin_size))
        Vb = np.zeros((V.shape[0], nbins), dtype=float)
        for i in range(nbins):
            L = i * bin_size
            R = min((i + 1) * bin_size, ncols)
            Vb[:, i] = V[:, L:R].mean(axis=1)
        V = Vb

    return V, bin_size


# -----------------------------
# Custom colormaps
# -----------------------------
def load_colormap_sheet(excel_path, sheet_name):
    df = pd.read_excel(excel_path, sheet_name=sheet_name)
    if df.empty:
        raise ValueError(f"Sheet '{sheet_name}' in '{excel_path}' is empty.")

    cols_lower = [str(c).strip().lower() for c in df.columns]

    def find_channel(ch):
        for i, name in enumerate(cols_lower):
            if name.startswith(ch):
                return df.columns[i]
        raise ValueError(f"No column starting with '{ch}' in sheet '{sheet_name}'.")

    r_col = find_channel('r')
    g_col = find_channel('g')
    b_col = find_channel('b')

    df = df[[r_col, g_col, b_col]].dropna()
    colors = df.to_numpy(dtype=float)
    if colors.max() > 1.0:
        colors = colors / 255.0

    return ListedColormap(colors, name=str(sheet_name))


def load_custom_colormaps(SET):
    cmap_dict = {}
    if not SET.get('use_custom_colormaps', False):
        return cmap_dict

    excel_path = SET.get('custom_cmap_excel', '')
    if not excel_path:
        print("Custom colormaps enabled but 'custom_cmap_excel' is empty; skipping.")
        return cmap_dict
    if not os.path.exists(excel_path):
        print(f"Custom colormap Excel not found: {excel_path}; skipping.")
        return cmap_dict

    sheet_map = SET.get('custom_cmap_sheets', {})
    loaded = []
    for name, sheet in sheet_map.items():
        try:
            cmap_dict[name] = load_colormap_sheet(excel_path, sheet)
            loaded.append(name)
        except Exception as e:
            print(f"Could not load colormap '{name}' from sheet '{sheet}': {e}")

    if loaded:
        print("[INFO] Loaded custom colormaps:", ", ".join(loaded))
    return cmap_dict


# -----------------------------
# Density colors (scatter)
# -----------------------------
def _compute_density_colors(Y, nbins=150, sigma=4.0, use_log=True,
                           min_rel=0.0, max_rel=1.0, cmap=None,
                           metric="density", enrichment_eps=1e-12,
                           enrichment_use_log=True):
    if cmap is None:
        cmap = plt.get_cmap("plasma")
    if Y is None or Y.shape[0] == 0:
        return None

    metric = (metric or "density").strip().lower()
    x = Y[:, 0]
    y = Y[:, 1]

    x_range = x.max() - x.min()
    y_range = y.max() - y.min()
    x_pad = 0.05 * x_range if x_range > 0 else 1.0
    y_pad = 0.05 * y_range if y_range > 0 else 1.0

    x_min = x.min() - x_pad
    x_max = x.max() + x_pad
    y_min = y.min() - y_pad
    y_max = y.max() + y_pad

    xedges = np.linspace(x_min, x_max, nbins + 1)
    yedges = np.linspace(y_min, y_max, nbins + 1)
    H, _, _ = np.histogram2d(x, y, bins=[xedges, yedges])
    Hs_lin = gaussian_filter(H, sigma=sigma, mode="nearest")

    if not np.isfinite(Hs_lin).any() or Hs_lin.max() <= 0:
        return cmap(np.zeros_like(x))

    ix = np.searchsorted(xedges, x, side="right") - 1
    iy = np.searchsorted(yedges, y, side="right") - 1
    ix = np.clip(ix, 0, nbins - 1)
    iy = np.clip(iy, 0, nbins - 1)
    dens_lin = Hs_lin[ix, iy]

    if metric == "density":
        dens_score = np.log1p(dens_lin) if use_log else dens_lin
    elif metric == "enrichment":
        base = float(np.nanmean(dens_lin)) if np.isfinite(dens_lin).any() else 0.0
        base = max(base, 0.0)
        ratio = (dens_lin + enrichment_eps) / (base + enrichment_eps)
        dens_score = np.log1p(ratio) if enrichment_use_log else ratio
    else:
        raise ValueError("metric must be 'density' or 'enrichment'")

    if not np.isfinite(dens_score).any() or dens_score.max() <= 0:
        return cmap(np.zeros_like(x))

    dmax = float(np.nanmax(dens_score))
    min_rel = max(0.0, min(1.0, float(min_rel)))
    max_rel = max(min_rel, min(1.0, float(max_rel)))

    vmin = min_rel * dmax
    vmax = max(vmin + 1e-12, max_rel * dmax)

    dens_clipped = np.clip(dens_score, vmin, vmax)
    norm = (dens_clipped - vmin) / (vmax - vmin)
    norm = np.clip(norm, 0.0, 1.0)
    return cmap(norm)


# -----------------------------
# Plotting
# -----------------------------
def plot_heatmap(
    SET,
    V,
    features_reorder,
    n_genes_full,
    bin_size,
    title_txt,
    custom_cmaps,
    x_label='Genes ordered by codon usage similarity',
    y_label='Codons',
    colorbar_label=None,
    caxis_limits=None,
):
    fig, ax = plt.subplots(figsize=SET.get('heatmap_fig_size', (18, 4)), dpi=int(SET.get('figure_dpi', 300)))

    cmap_name = SET.get('heatmap_colormap_name', 'parula')
    cmap = get_cmap_any(cmap_name, custom_maps=custom_cmaps, fallback='plasma')

    im = ax.imshow(V, aspect='auto', origin='lower', cmap=cmap)
    use_caxis = SET.get('heatmap_caxis_limits', None) if caxis_limits is None else caxis_limits
    if use_caxis is not None:
        im.set_clim(*use_caxis)
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.tick_params(direction='out')
    cbar.set_label(colorbar_label or SET.get('colorbar_title_string', r'$\sigma_{\mathrm{codon}}$'), fontsize=int(SET.get('colorbar_title_size', 11)))

    ax.set_xlabel(str(x_label), fontname=SET.get('font_name', 'Arial'))
    ax.set_ylabel(str(y_label), fontname=SET.get('font_name', 'Arial'))
    ax.tick_params(direction='out')
    ax.tick_params(axis='x', labelsize=int(SET.get('font_size_xticks', 8)))
    ax.tick_params(axis='y', labelsize=int(SET.get('font_size_yticks', 7)))

    step = max(1, int(SET.get('xtick_every_genes', 500)))
    last_multiple = (n_genes_full - 1) // step * step
    ticks_genes = np.arange(0, last_multiple + 1, step)
    tick_pos = 1 + (ticks_genes // bin_size if bool(SET.get('apply_binning', False)) else ticks_genes)
    ax.set_xticks(tick_pos - 1)
    ax.set_xticklabels([str(x) for x in ticks_genes])

    ax.set_yticks(np.arange(len(features_reorder)))
    ax.set_yticklabels(features_reorder, fontsize=int(SET.get('font_size_yticks', 7)))

    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    try:
        if SET.get('heatmap_xmin', None) not in (None, "") or SET.get('heatmap_xmax', None) not in (None, ""):
            ax.set_xlim(left=(float(SET.get('heatmap_xmin')) if SET.get('heatmap_xmin', None) not in (None, "") else None),
                        right=(float(SET.get('heatmap_xmax')) if SET.get('heatmap_xmax', None) not in (None, "") else None))
    except Exception:
        pass
    try:
        if SET.get('heatmap_ymin', None) not in (None, "") or SET.get('heatmap_ymax', None) not in (None, ""):
            ax.set_ylim(bottom=(float(SET.get('heatmap_ymin')) if SET.get('heatmap_ymin', None) not in (None, "") else None),
                        top=(float(SET.get('heatmap_ymax')) if SET.get('heatmap_ymax', None) not in (None, "") else None))
    except Exception:
        pass
    ax.set_title(title_txt, fontsize=int(SET.get('font_size_titles', 10)), fontname=SET.get('font_name', 'Arial'))
    return fig, ax


def _plot_feature_by_gene_heatmap(
    SET,
    values_df: pd.DataFrame,
    ordered_genes,
    title_txt: str,
    y_label: str,
    out_path: str,
    custom_cmaps,
    show_fig: bool = True,
    colorbar_label: str | None = None,
    caxis_limits=None,
    feature_label_map: dict | None = None,
    reorder_feature_rows: bool = True,
):
    if values_df is None or values_df.empty:
        return ""

    ordered_gene_list = [str(g) for g in (list(ordered_genes) if ordered_genes is not None else [])]
    if not ordered_gene_list:
        return ""

    df_ord = values_df.reindex(index=ordered_gene_list)
    if df_ord.empty or df_ord.shape[1] == 0:
        return ""

    values_raw = df_ord.to_numpy(dtype=float)
    values_for_order = _impute_missing_usage_for_embedding(values_raw)

    feature_labels = np.asarray(df_ord.columns.astype(str).tolist(), dtype=object)
    if reorder_feature_rows and values_for_order.shape[1] > 1:
        features_reorder, feat_order = reorder_features(SET, values_for_order, feature_labels)
    else:
        features_reorder = feature_labels
        feat_order = np.arange(values_for_order.shape[1], dtype=int)

    V = values_for_order[:, feat_order].T
    V_smooth, bin_size = smooth_and_bin(SET, V)

    display_labels = [str(feature_label_map.get(str(f), str(f))) if feature_label_map else str(f) for f in list(features_reorder)]

    local_set = dict(SET)
    use_custom = bool(SET.get('trna_supp_heatmaps_customize', False))
    default_width = float((SET.get('heatmap_fig_size', (18, 4)) or (18, 4))[0])
    main_rows = max(1, _main_heatmap_row_count_from_settings(SET))
    main_height = float((SET.get('heatmap_fig_size', (18, 4)) or (18, 4))[1])
    mirrored_cell_h = main_height / float(main_rows)
    nrows_here = max(1, len(display_labels))

    fig_w = _optional_float_value(SET.get('trna_supp_heatmaps_fig_width')) if use_custom else None
    fig_h = _optional_float_value(SET.get('trna_supp_heatmaps_fig_height')) if use_custom else None
    cell_h = _optional_float_value(SET.get('trna_supp_heatmaps_cell_height')) if use_custom else None
    dpi = _optional_int_value(SET.get('trna_supp_heatmaps_dpi')) if use_custom else None
    xtick_every = _optional_int_value(SET.get('trna_supp_heatmaps_xtick_every_genes')) if use_custom else None
    ytick_fs = _optional_int_value(SET.get('trna_supp_heatmaps_ytick_fontsize')) if use_custom else None
    title_fs = _optional_int_value(SET.get('trna_supp_heatmaps_title_fontsize')) if use_custom else None
    x_min = _optional_float_value(SET.get('trna_supp_heatmaps_xmin')) if use_custom else None
    x_max = _optional_float_value(SET.get('trna_supp_heatmaps_xmax')) if use_custom else None
    y_min = _optional_float_value(SET.get('trna_supp_heatmaps_ymin')) if use_custom else None
    y_max = _optional_float_value(SET.get('trna_supp_heatmaps_ymax')) if use_custom else None

    effective_cell_h = float(cell_h if cell_h is not None else mirrored_cell_h)
    effective_fig_h = float(fig_h if fig_h is not None else max(2.0, effective_cell_h * nrows_here))
    local_set['heatmap_fig_size'] = (float(fig_w if fig_w is not None else default_width), effective_fig_h)
    if dpi is not None:
        local_set['figure_dpi'] = int(dpi)
    if xtick_every is not None:
        local_set['xtick_every_genes'] = int(max(1, xtick_every))
    if ytick_fs is not None:
        local_set['font_size_yticks'] = int(ytick_fs)
    if title_fs is not None:
        local_set['font_size_titles'] = int(title_fs)
    if x_min is not None:
        local_set['heatmap_xmin'] = x_min
    if x_max is not None:
        local_set['heatmap_xmax'] = x_max
    if y_min is not None:
        local_set['heatmap_ymin'] = y_min
    if y_max is not None:
        local_set['heatmap_ymax'] = y_max

    fig, _ = plot_heatmap(
        local_set,
        V_smooth,
        display_labels,
        len(ordered_gene_list),
        bin_size,
        '',
        custom_cmaps,
        x_label='Genes ordered by codon-usage clustering',
        y_label=y_label,
        colorbar_label=colorbar_label,
        caxis_limits=caxis_limits,
    )
    fig.savefig(out_path, dpi=local_set['figure_dpi'], bbox_inches='tight')
    if show_fig:
        _show_figure_nonblocking(fig)
    else:
        plt.close(fig)
    return out_path


_RNA_WC_PAIRS = {
    ('A', 'U'), ('U', 'A'), ('G', 'C'), ('C', 'G')
}
_RNA_WOBBLE_PAIRS = {
    ('G', 'U'), ('U', 'G'),
    ('I', 'A'), ('I', 'U'), ('I', 'C'),
}


def _split_anticodon_group(label: str):
    s = str(label or '').strip()
    if not s:
        return []
    return [x.strip().upper().replace('T', 'U') for x in s.split('/') if x.strip()]


def _pairing_mode_for_anticodon_codon(anticodon_5to3: str, codon_internal: str) -> str:
    anti = str(anticodon_5to3 or '').strip().upper().replace('T', 'U')
    cod = str(codon_internal or '').strip().upper()
    if '_' in cod:
        cod = cod.split('_', 1)[1]
    cod = cod.replace('T', 'U')
    if len(anti) != 3 or len(cod) != 3:
        return 'Other'
    pair = (anti[0], cod[2])
    if pair in _RNA_WC_PAIRS:
        return 'WC'
    if pair in _RNA_WOBBLE_PAIRS:
        return 'Wobble'
    return 'Other'


def _split_pair_label_family_fields(label: str):
    s = str(label or '')
    if '__' not in s:
        return s, ''
    return s.split('__', 1)[0], s.split('__', 1)[1]


def _pair_label_to_display(label: str) -> str:
    s = str(label or '')
    if '__' not in s:
        return s
    base, mode = s.split('__', 1)
    pretty = 'WC' if str(mode).strip().lower() == 'wc' else ('Wobble' if str(mode).strip().lower() == 'wobble' else str(mode))
    return f"{base} | {pretty}"


def _metric_caxis_from_name(SET, metric_name: str):
    metric = str(metric_name or '').strip().upper()
    if metric == 'ATU':
        return None
    if metric == 'RTU':
        return (0.0, 1.0)
    return SET.get('heatmap_caxis_limits', None)


def _optional_float_value(x):
    try:
        if x in (None, ''):
            return None
        return float(x)
    except Exception:
        return None


def _positive_optional_float_value(x):
    try:
        if x in (None, ''):
            return None
        val = float(x)
        return val if val > 0 else None
    except Exception:
        return None


def _optional_int_value(x):
    try:
        if x in (None, ''):
            return None
        return int(float(x))
    except Exception:
        return None


def _main_heatmap_row_count_from_settings(SET):
    usage = str(SET.get('usage_basis', 'ACU') or 'ACU').upper().strip()
    codon_set = str(SET.get('codon_set', '61') or '61').strip()
    if usage == 'AA':
        return 20
    if codon_set.startswith('59'):
        return 59
    if codon_set.startswith('64'):
        return 64
    return 61


def _select_requested_clusters(available_clusters, requested):
    available = [str(c) for c in list(available_clusters or []) if str(c).strip()]
    if not available:
        return []
    s = str(requested or '').strip()
    if s == '' or s.lower() == 'all':
        return available
    tokens = [t.strip() for t in re.split(r'[;,]+', s) if t.strip()]
    lower_map = {str(c).strip().lower(): str(c) for c in available}
    chosen = []
    for tok in tokens:
        if tok.isdigit():
            idx = int(tok) - 1
            if 0 <= idx < len(available):
                chosen.append(available[idx])
                continue
        mapped = lower_map.get(tok.lower())
        if mapped is not None:
            chosen.append(mapped)
    dedup = []
    seen = set()
    for c in chosen:
        if c not in seen:
            dedup.append(c)
            seen.add(c)
    return dedup or available


def _cluster_feature_mean_matrix(values_df: pd.DataFrame, cluster_df: pd.DataFrame, selected_clusters):
    if values_df is None or values_df.empty or cluster_df is None or getattr(cluster_df, 'empty', True):
        return pd.DataFrame()
    if selected_clusters is None:
        selected = [str(c) for c in cluster_df.columns]
    else:
        selected = [str(c) for c in list(selected_clusters or []) if str(c) in list(cluster_df.columns)]
    if not selected:
        return pd.DataFrame(index=values_df.columns)
    out = pd.DataFrame(index=values_df.columns, columns=selected, dtype=float)
    for cname in selected:
        genes = [str(v).strip() for v in cluster_df[cname].replace({np.nan: ''}).astype(str).tolist()]
        genes = [g for g in genes if g and g.lower() != 'nan' and g in values_df.index]
        if not genes:
            continue
        sub = values_df.loc[genes].apply(pd.to_numeric, errors='coerce')
        out[cname] = sub.mean(axis=0, skipna=True)
    return out


def _feature_group_bounds(labels, sep='_'):
    bounds = []
    start = 0
    current = None
    labels = [str(x) for x in list(labels or [])]
    for i, lab in enumerate(labels):
        grp = lab.split(sep, 1)[0] if sep in lab else lab
        if current is None:
            current = grp
            start = i
        elif grp != current:
            bounds.append((start, i - 1, current))
            current = grp
            start = i
    if labels:
        bounds.append((start, len(labels) - 1, current))
    return bounds


def _plot_cluster_shift_heatmap(
    SET,
    matrix_df: pd.DataFrame,
    out_path: str,
    title_txt: str,
    show_fig: bool = True,
    colorbar_label: str = 'Mean z-score',
    group_bounds=None,
    row_display_labels=None,
    x_label: str = '',
    y_label: str = 'Features',
    log2_colorbar: bool = False,
    bracket_type: str = 'square',
    bracket_x: float | None = None,
    label_x: float | None = None,
    bracket_lw: float = 1.4,
    red_row_mask=None,
):
    if matrix_df is None or matrix_df.empty:
        return ''

    use_custom = bool(SET.get('trna_shift_heatmaps_customize', False))
    dpi = _optional_int_value(SET.get('trna_shift_heatmaps_dpi')) if use_custom else None
    fig_w = _positive_optional_float_value(SET.get('trna_shift_heatmaps_fig_width')) if use_custom else None
    fig_h = _positive_optional_float_value(SET.get('trna_shift_heatmaps_fig_height')) if use_custom else None
    cell_w = _positive_optional_float_value(SET.get('trna_shift_heatmaps_cell_width')) if use_custom else None
    cell_h_scale = _positive_optional_float_value(SET.get('trna_shift_heatmaps_cell_height')) if use_custom else None
    xtick_fs = _optional_int_value(SET.get('trna_shift_heatmaps_xtick_fontsize')) if use_custom else None
    ytick_fs = _optional_int_value(SET.get('trna_shift_heatmaps_ytick_fontsize')) if use_custom else None
    title_fs = _optional_int_value(SET.get('trna_shift_heatmaps_title_fontsize')) if use_custom else None
    x_min = _optional_float_value(SET.get('trna_shift_heatmaps_xmin')) if use_custom else None
    x_max = _optional_float_value(SET.get('trna_shift_heatmaps_xmax')) if use_custom else None
    y_min = _optional_float_value(SET.get('trna_shift_heatmaps_ymin')) if use_custom else None
    y_max = _optional_float_value(SET.get('trna_shift_heatmaps_ymax')) if use_custom else None

    dpi = int(dpi if dpi is not None else SET.get('figure_dpi', 300))
    nrows, ncols = matrix_df.shape

    main_rows = max(1, _main_heatmap_row_count_from_settings(SET))
    main_fig_w, main_fig_h = (SET.get('heatmap_fig_size', (18, 4)) or (18, 4))
    ref_cell_h = float(main_fig_h) / float(main_rows)

    effective_cell_w = float(cell_w if cell_w is not None else 0.8)
    effective_cell_h = float(ref_cell_h * (cell_h_scale if cell_h_scale is not None else 1.5))

    fig_w = float(fig_w if fig_w is not None else max(5.4, 3.2 + ncols * effective_cell_w))
    fig_h = float(fig_h if fig_h is not None else max(4.0, 1.8 + nrows * effective_cell_h))
    xtick_fs = int(xtick_fs if xtick_fs is not None else max(7, SET.get('font_size_xticks', 8) - 1))
    ytick_fs = int(ytick_fs if ytick_fs is not None else max(8, SET.get('font_size_yticks', 8) + 3))
    title_fs = int(title_fs if title_fs is not None else max(10, SET.get('font_size_titles', 10) + 2))

    cmap = plt.get_cmap('RdBu_r')
    arr = np.asarray(matrix_df.values, dtype=float)
    if bool(log2_colorbar):
        arr = np.sign(arr) * np.log2(np.abs(arr) + 1.0)
    finite = arr[np.isfinite(arr)]
    vmax = float(np.nanmax(np.abs(finite))) if finite.size else 1.0
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = 1.0
    vmin = -vmax

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)
    im = ax.imshow(arr, aspect='auto', origin='upper', cmap=cmap, vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(im, ax=ax)
    cbar_text = str(colorbar_label)
    if bool(log2_colorbar):
        cbar_text = f'log2-transformed\n{cbar_text}'
    cbar.set_label(cbar_text, fontsize=max(9, title_fs - 1))
    cbar.ax.tick_params(labelsize=max(8, xtick_fs - 1))

    xlabels = [wrap_label_no_break(str(c), 12) for c in matrix_df.columns]
    ylabels = list(row_display_labels) if row_display_labels is not None else [str(i) for i in matrix_df.index]
    ax.set_xticks(np.arange(len(xlabels)))
    ax.set_xticklabels(xlabels, rotation=0, fontsize=xtick_fs)
    ax.set_yticks(np.arange(len(ylabels)))
    red_mask = list(red_row_mask) if red_row_mask is not None else [False] * len(ylabels)
    red_mask = red_mask + [False] * max(0, len(ylabels) - len(red_mask))
    if any(red_mask):
        ax.set_yticklabels([''] * len(ylabels))
        trans = ax.get_yaxis_transform()
        for i_lab, lab in enumerate(ylabels):
            ax.text(-0.02, i_lab, str(lab), transform=trans, ha='right', va='center',
                    fontsize=ytick_fs, color=('red' if red_mask[i_lab] else 'black'))
    else:
        ax.set_yticklabels(ylabels, fontsize=ytick_fs)

    # Titles are intentionally suppressed for Decoding strategies figures.
    ax.set_title('')
    ax.set_xlabel(str(x_label), fontsize=max(9, xtick_fs + 1))
    if str(y_label).strip():
        ax.set_ylabel(str(y_label), fontsize=max(9, ytick_fs + 1))
    else:
        ax.set_ylabel('')

    if group_bounds:
        for start, end, _lab in group_bounds:
            rect = Rectangle((-0.5, start - 0.5), len(xlabels), end - start + 1,
                             fill=False, edgecolor='black', linewidth=0.9)
            ax.add_patch(rect)
        eff_label_x = float(label_x if label_x is not None else -0.34)
        eff_bracket_x = float(bracket_x if bracket_x is not None else -0.18)
        cap_len = 0.03
        btype = str(bracket_type or 'square').strip().lower()
        trans = ax.get_yaxis_transform()
        for start, end, lab in group_bounds:
            yc = 0.5 * (start + end)
            y0 = start - 0.35
            y1 = end + 0.35
            ax.text(eff_label_x, yc, str(lab), ha='right', va='center', fontsize=max(10, ytick_fs + 1),
                    fontweight='bold', clip_on=False, transform=trans)
            if btype in {'brace', 'curly'}:
                ym = 0.5 * (y0 + y1)
                ax.plot([eff_bracket_x + cap_len, eff_bracket_x, eff_bracket_x + cap_len], [y0, y0, ym],
                        transform=trans, color='black', linewidth=bracket_lw, solid_capstyle='butt', clip_on=False, zorder=7)
                ax.plot([eff_bracket_x + cap_len, eff_bracket_x, eff_bracket_x + cap_len], [ym, y1, y1],
                        transform=trans, color='black', linewidth=bracket_lw, solid_capstyle='butt', clip_on=False, zorder=7)
            else:
                ax.plot([eff_bracket_x, eff_bracket_x], [y0, y1], transform=trans, color='black', linewidth=bracket_lw,
                        solid_capstyle='butt', clip_on=False, zorder=7)
                ax.plot([eff_bracket_x, eff_bracket_x + cap_len], [y0, y0], transform=trans, color='black', linewidth=bracket_lw,
                        solid_capstyle='butt', clip_on=False, zorder=7)
                ax.plot([eff_bracket_x, eff_bracket_x + cap_len], [y1, y1], transform=trans, color='black', linewidth=bracket_lw,
                        solid_capstyle='butt', clip_on=False, zorder=7)
        fig.subplots_adjust(left=0.33, right=0.93, bottom=0.16, top=0.92)
    else:
        fig.subplots_adjust(left=0.18, right=0.93, bottom=0.16, top=0.92)

    if x_min is not None or x_max is not None:
        ax.set_xlim(left=x_min, right=x_max)
    if y_min is not None or y_max is not None:
        cur_bottom, cur_top = ax.get_ylim()
        ax.set_ylim(bottom=(y_min if y_min is not None else cur_bottom),
                    top=(y_max if y_max is not None else cur_top))

    fig.savefig(out_path, dpi=dpi, bbox_inches='tight')
    if show_fig:
        _show_figure_nonblocking(fig)
    else:
        plt.close(fig)
    return out_path


def _build_trna_family_maps(trna_rules: dict):
    ordered_labels = [str(v) for v in list(trna_rules.get('ordered_trna_labels', []) or [])]
    aa_to_labels = {}
    aa_order = []
    for lab in ordered_labels:
        aa = str(lab).split('_', 1)[0].strip()
        if not aa:
            continue
        aa_to_labels.setdefault(aa, []).append(lab)
        if aa not in aa_order:
            aa_order.append(aa)

    aa_to_codons = {}
    table_df = trna_rules.get('table_df')
    if isinstance(table_df, pd.DataFrame) and (not table_df.empty):
        for _, row in table_df.iterrows():
            aa = str(row.get('AA', '')).strip()
            codon_internal = str(row.get('Codon_internal', '')).strip()
            if not aa or not codon_internal:
                continue
            aa_to_codons.setdefault(aa, [])
            if codon_internal not in aa_to_codons[aa]:
                aa_to_codons[aa].append(codon_internal)
            if aa not in aa_order:
                aa_order.append(aa)
    return aa_order, aa_to_labels, aa_to_codons



def _adaptive_cluster_colors(n):
    n = int(max(1, n))
    cmap_name = 'tab10' if n <= 10 else ('tab20' if n <= 20 else 'nipy_spectral')
    cmap = plt.get_cmap(cmap_name)
    if n == 1:
        return [cmap(0)]
    if cmap_name in {'tab10', 'tab20'}:
        return [cmap(i % cmap.N) for i in range(n)]
    return [cmap(i / max(1, n - 1)) for i in range(n)]


def _safe_total_sense_counts(count_df: pd.DataFrame):
    if count_df is None or count_df.empty:
        return pd.Series(dtype=float)
    cols = [c for c in count_df.columns if not str(c).upper().startswith('STOP')]
    return count_df.loc[:, cols].apply(pd.to_numeric, errors='coerce').sum(axis=1).replace(0, np.nan)


def _genome_codon_frequency_percent(count_df: pd.DataFrame, codon_labels):
    labels = [str(x) for x in list(codon_labels or [])]
    if count_df is None or count_df.empty:
        return pd.Series([np.nan] * len(labels), index=labels, dtype=float)
    cols = [c for c in count_df.columns if not str(c).upper().startswith('STOP')]
    total = pd.to_numeric(count_df.loc[:, cols].sum(axis=0), errors='coerce').sum()
    out = []
    for lab in labels:
        val = pd.to_numeric(count_df[lab], errors='coerce').sum() if lab in count_df.columns else np.nan
        out.append(float(val) / float(total) * 100.0 if pd.notna(val) and total and total > 0 else np.nan)
    return pd.Series(out, index=labels, dtype=float, name='Genomic codon frequency (%)')

def _wrap_title_max_words(text, max_words=6):
    words = str(text or '').split()
    if len(words) <= int(max_words):
        return str(text or '')
    mid = int(math.ceil(len(words) / 2.0))
    return ' '.join(words[:mid]) + '\n' + ' '.join(words[mid:])


def _rolling_nanmean(values, window=5):
    ser = pd.Series(np.asarray(values, dtype=float))
    window = max(1, int(window or 1))
    if window <= 1:
        return ser.to_numpy(dtype=float)
    return ser.rolling(window=window, min_periods=1, center=True).mean().to_numpy(dtype=float)


def _smooth_profile_values(values, method='running average', window=5):
    arr = np.asarray(values, dtype=float)
    method = str(method or '').strip().lower().replace('_', ' ').replace('-', ' ')
    window = max(1, int(window or 1))
    if method in {'', 'none', 'no', 'false', 'off', '0'} or window <= 1:
        return arr
    ser = pd.Series(arr, dtype=float)
    if 'median' in method:
        return ser.rolling(window=window, min_periods=1, center=True).median().to_numpy(dtype=float)
    if 'gauss' in method:
        # Gaussian smoothing while preserving NaN gaps as much as possible.
        try:
            from scipy.ndimage import gaussian_filter1d
            valid = np.isfinite(arr).astype(float)
            filled = np.where(np.isfinite(arr), arr, 0.0)
            sigma = max(0.5, float(window) / 2.0)
            num = gaussian_filter1d(filled, sigma=sigma, mode='nearest')
            den = gaussian_filter1d(valid, sigma=sigma, mode='nearest')
            out = num / np.where(den <= 1e-12, np.nan, den)
            return out.astype(float)
        except Exception:
            return ser.rolling(window=window, min_periods=1, center=True).mean().to_numpy(dtype=float)
    # Default: centered running average.
    return ser.rolling(window=window, min_periods=1, center=True).mean().to_numpy(dtype=float)


def _plot_ordered_profile_line_surface(
    SET,
    values,
    ordered_genes,
    out_path: str,
    title_txt: str,
    y_label: str,
    plot_kind: str = 'line',
    smooth: bool = True,
    smooth_window: int = 5,
    smooth_method: str = 'running average',
    show_fig: bool = True,
    caption_size: int | None = None,
):
    ordered = [str(g) for g in list(ordered_genes)]
    if isinstance(values, pd.DataFrame):
        data = values.reindex(index=ordered).apply(pd.to_numeric, errors='coerce')
    else:
        data = pd.DataFrame({'Profile': pd.Series(values, index=ordered, dtype=float).reindex(ordered)})
    if data.dropna(how='all').empty:
        print(f'[WARN] {title_txt}: no finite values available; figure was not saved.')
        return ''

    smooth_window = int(max(1, smooth_window or 1))
    smooth_method = str(smooth_method or 'running average').strip().lower()
    if bool(smooth):
        if smooth_method in {'', 'yes', 'true', 'on', '1'}:
            smooth_method = 'running average'
    else:
        smooth_method = 'none'

    if smooth_method in {'none', 'no', 'false', 'off', '0'} or smooth_window <= 1:
        plot_data = data.copy()
        smooth_note = ''
    else:
        plot_data = data.copy()
        for col in plot_data.columns:
            plot_data[col] = _smooth_profile_values(plot_data[col].to_numpy(dtype=float), method=smooth_method, window=smooth_window)
        smooth_note = f' ({smooth_method}, {smooth_window}-gene window)'

    x = np.arange(len(plot_data), dtype=float)
    dpi = int(SET.get('figure_dpi', 300))
    fig_w = float((SET.get('heatmap_fig_size', (18, 4)) or (18, 4))[0])
    fig_h = max(3.2, float((SET.get('heatmap_fig_size', (18, 4)) or (18, 4))[1]) * 0.9)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)
    kind = str(plot_kind or 'line').strip().lower().replace('surface', 'area')
    caption_fs = int(caption_size if caption_size is not None else SET.get('trna_profile_caption_size', 13))
    caption_fs = max(6, caption_fs)
    cmap = plt.get_cmap('tab10' if plot_data.shape[1] <= 10 else 'tab20')

    for i, col in enumerate(plot_data.columns):
        y_plot = plot_data[col].to_numpy(dtype=float)
        valid = np.isfinite(y_plot)
        if not valid.any():
            continue
        color = cmap(i % cmap.N)
        label = str(col)
        if kind in {'area', 'surface'}:
            baseline = 0.0
            finite = y_plot[valid]
            if finite.size and not (np.nanmin(finite) <= 0 <= np.nanmax(finite)):
                baseline = float(np.nanmin(finite))
            ax.fill_between(x, baseline, y_plot, where=valid, alpha=0.18, linewidth=0, color=color)
            ax.plot(x[valid], y_plot[valid], lw=1.45, color=color, label=label)
        else:
            ax.plot(x[valid], y_plot[valid], lw=1.65, color=color, label=label)

    # Titles are intentionally suppressed for Decoding strategies figures;
    # legends/axis labels carry the necessary information and avoid clutter above the plot.
    ax.set_title('')
    ax.set_xlabel('Genes ordered by codon-usage clustering', fontsize=caption_fs)
    ax.set_ylabel(str(y_label), fontsize=caption_fs)
    ax.grid(axis='y', alpha=0.25, lw=0.6)
    if plot_data.shape[1] > 1:
        ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1.0), frameon=True, fontsize=caption_fs)
    step = max(1, int(SET.get('xtick_every_genes', 500)))
    ticks = np.arange(0, max(1, len(plot_data)), step)
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(int(t)) for t in ticks], fontsize=caption_fs)
    ax.tick_params(axis='both', labelsize=caption_fs)
    if len(plot_data) > 0:
        ax.set_xlim(0, max(0, len(plot_data) - 1))
        try:
            ax.margins(x=0)
        except Exception:
            pass
    fig.tight_layout()
    os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches='tight')
    if show_fig:
        _show_figure_nonblocking(fig)
    else:
        plt.close(fig)
    return out_path


def _compute_wobble_percent_profile(count_df: pd.DataFrame, ordered_genes, requested_pairs):
    ordered = [str(g) for g in list(ordered_genes or [])]
    pair_values = []
    for aa, wc, wob in list(requested_pairs or []):
        wc_lab = f'{aa}_{str(wc).upper().replace("U", "T")}'
        wob_lab = f'{aa}_{str(wob).upper().replace("U", "T")}'
        if wc_lab not in count_df.columns or wob_lab not in count_df.columns:
            continue
        wc_vals = pd.to_numeric(count_df.reindex(index=ordered)[wc_lab], errors='coerce')
        wob_vals = pd.to_numeric(count_df.reindex(index=ordered)[wob_lab], errors='coerce')
        denom = wc_vals + wob_vals
        pct = wob_vals.div(denom.replace(0, np.nan)) * 100.0
        pair_values.append(pct)
    if not pair_values:
        return pd.Series([np.nan] * len(ordered), index=ordered, dtype=float)
    mat = pd.concat(pair_values, axis=1)
    return mat.mean(axis=1, skipna=True).reindex(ordered)


def _compute_weighted_trna_abundance_profile(count_df: pd.DataFrame, trna_rules: dict, ordered_genes):
    ordered = [str(g) for g in list(ordered_genes or [])]
    if count_df is None or count_df.empty:
        return pd.Series([np.nan] * len(ordered), index=ordered, dtype=float)
    abundance = trna_rules.get('abundance_series')
    if abundance is None or getattr(abundance, 'empty', True):
        return pd.Series([np.nan] * len(ordered), index=ordered, dtype=float)
    abundance = pd.to_numeric(pd.Series(abundance), errors='coerce')
    codon_to_decoders = dict(trna_rules.get('codon_to_decoders', {}) or {})
    total = _safe_total_sense_counts(count_df).reindex(ordered)
    out = pd.Series(0.0, index=ordered, dtype=float)
    any_contrib = False
    for codon in count_df.columns:
        codon_s = str(codon)
        if codon_s not in codon_to_decoders:
            continue
        counts = pd.to_numeric(count_df.reindex(index=ordered)[codon_s], errors='coerce').fillna(0.0)
        frac = counts.div(total)
        decs = list(codon_to_decoders.get(codon_s, []) or [])
        weighted_ab = 0.0
        weight_sum = 0.0
        for lab, weight in decs:
            val = abundance.get(str(lab), np.nan)
            if pd.notna(val):
                weighted_ab += float(val) * float(weight)
                weight_sum += float(weight)
        if weight_sum > 0:
            out = out.add(frac * (weighted_ab / weight_sum), fill_value=0.0)
            any_contrib = True
    out[~np.isfinite(total)] = np.nan
    if not any_contrib:
        out[:] = np.nan
    return out.reindex(ordered)



def _normalize_anti_key_for_profile(value):
    return str(value or '').strip().upper().replace('T', 'U').replace(' ', '')


def _count_df_column_for_aa_codon(count_df: pd.DataFrame, aa: str, codon: str):
    """Return the column name for an amino-acid/codon pair in a count table.

    Direct decoding replots usually use internal CodonPipe columns such as
    ``Ile_ATA``.  Some ad-hoc raw tables may instead contain bare codons such as
    ``ATA``.  The unified decoding parser stores ``Codon_internal`` as
    ``AA_DNAcodon`` (for example ``Ile_ATA``); this helper therefore accepts
    either a bare RNA/DNA codon or an already-internalized label.
    """
    raw = str(codon or '').strip()
    aa_s = str(aa or '').strip()
    col_lookup = {str(c).strip(): c for c in count_df.columns}

    candidates = []
    if raw:
        # Already-internalized label from the parser, e.g. Ile_ATA.
        raw_internal = raw.replace('U', 'T')
        candidates.append(raw_internal)
        # Bare codon extracted from AA_codon.
        if '_' in raw_internal:
            candidates.append(raw_internal.split('_', 1)[1])
            candidates.append(raw_internal.rsplit('_', 1)[-1])
        else:
            candidates.append(raw_internal)

    codons = []
    for c in candidates:
        c = str(c or '').strip().upper().replace('U', 'T')
        if not c:
            continue
        # If c is still AA_CODON, keep it as a direct candidate and also try the final token.
        if '_' in c:
            codons.append(c)
            codons.append(c.rsplit('_', 1)[-1])
        else:
            codons.append(c)

    final_candidates = []
    for c in codons:
        if not c:
            continue
        if '_' in c:
            final_candidates.append(c)
        else:
            final_candidates.append(f'{aa_s}_{c}')
            final_candidates.append(c)

    seen = set()
    for cand in final_candidates:
        if cand in seen:
            continue
        seen.add(cand)
        if cand in col_lookup:
            return col_lookup[cand]
    return None


def _compute_rare_trna_fraction_profiles(count_df: pd.DataFrame, trna_rules: dict, ordered_genes):
    """Compute rare/(rare+abundant) decoding fractions for selected AA families.

    For each gene and each selected amino acid, codons are first counted, then
    assigned to their decoding tRNA/anticodon group using the uploaded decoding
    table.  The plotted value is::

        rare tRNA codon counts / (rare + abundant tRNA codon counts) * 100

    Thus Ile, for example, is computed as CAU-decoded AUA counts divided by all
    Ile codons assigned to CAU or GAU, multiplied by 100.
    """
    ordered = [str(g) for g in list(ordered_genes or [])]
    if count_df is None or count_df.empty:
        return pd.DataFrame(index=ordered)
    table = trna_rules.get('table_df') if isinstance(trna_rules, dict) else None
    if table is None or getattr(table, 'empty', True):
        return pd.DataFrame(index=ordered)

    # User-defined contrasts. Rare/abundant is based on the intended biological
    # comparison, not automatically ranked by absolute abundance.
    # Gly includes UCC/CCC (actual table) and UGC/CCC as a tolerant alias.
    contrasts = {
        'Arg': dict(rare=['CCG', 'UCU/CCU'], abundant=['ACG']),
        'Gly': dict(rare=['UCC/CCC', 'UGC/CCC'], abundant=['GCC']),
        # Corrected: Ile CAU decodes AUA and is the rare decoding group;
        # GAU decodes AUU/AUC and is the abundant/common decoding group.
        'Ile': dict(rare=['CAU'], abundant=['GAU']),
        'Leu': dict(rare=['UAA/CAA'], abundant=['UAG/GAG/CAG']),
        'Ser': dict(rare=['UGA/GGA/CGA'], abundant=['GCU']),
    }

    df = table.copy()
    aa_col = 'AA' if 'AA' in df.columns else df.columns[0]
    cod_col = 'Codon_internal' if 'Codon_internal' in df.columns else None
    anti_col = 'Anticodon' if 'Anticodon' in df.columns else None
    if cod_col is None or anti_col is None:
        return pd.DataFrame(index=ordered)

    out = pd.DataFrame(index=ordered, dtype=float)
    for aa, spec in contrasts.items():
        rare_keys = {_normalize_anti_key_for_profile(x) for x in spec.get('rare', [])}
        abundant_keys = {_normalize_anti_key_for_profile(x) for x in spec.get('abundant', [])}
        sub = df[df[aa_col].astype(str).str.strip() == aa]
        rare_cols, abundant_cols = [], []
        for _, row in sub.iterrows():
            cod = str(row.get(cod_col, '')).strip().upper().replace('U', 'T')
            anti = _normalize_anti_key_for_profile(row.get(anti_col, ''))
            col = _count_df_column_for_aa_codon(count_df, aa, cod)
            if col is None:
                continue
            if anti in rare_keys:
                rare_cols.append(col)
            if anti in abundant_keys:
                abundant_cols.append(col)
        rare_cols = list(dict.fromkeys(rare_cols))
        abundant_cols = list(dict.fromkeys(abundant_cols))
        if not rare_cols or not abundant_cols:
            print(f'[WARN] Rare-tRNA fraction: {aa} skipped because rare or abundant codon columns were not found.')
            continue
        rare_counts = count_df.reindex(index=ordered).loc[:, rare_cols].apply(pd.to_numeric, errors='coerce').fillna(0.0).sum(axis=1)
        abundant_counts = count_df.reindex(index=ordered).loc[:, abundant_cols].apply(pd.to_numeric, errors='coerce').fillna(0.0).sum(axis=1)
        denom = rare_counts + abundant_counts
        out[aa] = rare_counts.div(denom.replace(0, np.nan)).astype(float) * 100.0
        print(f'[INFO] Rare-tRNA fraction {aa}: rare={list(rare_cols)}; abundant={list(abundant_cols)}.')
    return out.reindex(index=ordered)

def _read_rna_stability_from_decoding_workbook(path: str, ordered_genes, sheet_name: str = 'RNA stability'):
    ordered = [str(g) for g in list(ordered_genes or [])]
    try:
        df = pd.read_excel(path, sheet_name=sheet_name, dtype=object)
    except Exception as e:
        print(f"[WARN] mRNA stability sheet '{sheet_name}' could not be read from decoding workbook: {e}")
        return pd.Series([np.nan] * len(ordered), index=ordered, dtype=float)
    if df is None or df.empty or df.shape[1] < 2:
        print(f"[WARN] mRNA stability sheet '{sheet_name}' is empty or has fewer than two columns.")
        return pd.Series([np.nan] * len(ordered), index=ordered, dtype=float)
    locus = df.iloc[:, 0].astype(str).str.strip()
    vals = pd.to_numeric(df.iloc[:, 1], errors='coerce')
    ser = pd.Series(vals.to_numpy(dtype=float), index=locus, dtype=float)
    ser = ser[~ser.index.duplicated(keep='first')]
    aligned = ser.reindex(ordered)
    print(f"[INFO] mRNA stability overlap: {aligned.notna().sum()} / {len(ordered)} reordered genes have half-life values.")
    return aligned



def _normalize_optional_string_list(value):
    """Normalize a GUI/list setting to a de-duplicated list, or None for 'all'."""
    if value is None:
        return None
    if isinstance(value, str):
        text = value.strip()
        if text == "" or text.lower() in {"all", "any", "none/all", "default"}:
            return None
        raw = re.split(r"[,;\n\r\t]+", text)
    else:
        raw = list(value or [])
    out, seen = [], set()
    for item in raw:
        val = str(item or "").strip()
        if not val or val.lower() in {"nan", "none", "na", "n/a"}:
            continue
        key = val.lower()
        if key in seen:
            continue
        seen.add(key)
        out.append(val)
    return out



PLOT6_DEFAULT_EXCLUDED_MODIFICATION_KEYS = {"ac4c34"}
PLOT6_DEFAULT_AAS = ["Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Phe", "Pro", "Ser", "Thr", "Tyr", "Val"]
PLOT6_DEFAULT_EXCLUDED_AA_KEYS = {"sec", "selenocysteine", "trp", "tryptophan", "met", "methionine"}


def _plot6_feature_key_ascii(value):
    """Return a permissive ASCII key for Plot 6 feature matching."""
    s = str(value or "").strip()
    trans = str.maketrans({
        "⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4",
        "⁵": "5", "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9",
        "₀": "0", "₁": "1", "₂": "2", "₃": "3", "₄": "4",
        "₅": "5", "₆": "6", "₇": "7", "₈": "8", "₉": "9",
    })
    s = s.translate(trans)
    s = s.replace("′", "'").replace("’", "'").replace("`", "'")
    s = re.sub(r"[^A-Za-z0-9]+", "", s)
    return s.lower()


def _canonicalize_plot6_modification_feature_label(feature):
    """Return the canonical Plot 6 label for one modification feature."""
    s = str(feature or "").strip()
    if not s:
        return ""
    key = _plot6_feature_key_ascii(s)
    if key in {"cmo5u34", "mcmo5u34"}:
        return "(m)cmo5U34"
    if key in {"cm32", "um32", "cm32um32"}:
        return "Cm32/Um32"
    if key in {"q34", "gluq34"}:
        return "Q34"
    if key in {"cm34", "cmnm5um34", "cmmn5um34"}:
        return "Cm34"
    if key == "mnm5u34":
        return "mnm5U34"
    if key in {"mnm5s2u34", "cmnm5s2u34"}:
        return "mnm5s2U34"
    if key in {"m6t6a37", "t6a37", "ct6a37"}:
        return "ct6A37"
    if key == "m6a37":
        return "m6A37"
    return s


def _canonicalize_plot6_modification_feature_list(values):
    """Canonicalize a user-selected Plot 6 modification list, preserving order."""
    vals = _normalize_optional_string_list(values)
    if vals is None:
        return None
    out, seen = [], set()
    for v in vals:
        c = _canonicalize_plot6_modification_feature_label(v)
        if not c:
            continue
        key = _plot6_feature_key_ascii(c)
        if key in seen:
            continue
        seen.add(key)
        out.append(c)
    return out


def _is_default_excluded_plot6_modification(feature):
    return _plot6_feature_key_ascii(feature) in PLOT6_DEFAULT_EXCLUDED_MODIFICATION_KEYS


def _plot6_aa_key(value):
    return re.sub(r"[^A-Za-z0-9]+", "", str(value or "").strip()).lower()


def _is_default_excluded_plot6_aa(aa):
    return _plot6_aa_key(aa) in PLOT6_DEFAULT_EXCLUDED_AA_KEYS


def _apply_default_plot6_aa_exclusions(include_aas):
    """Return explicit Plot 6 amino-acid filtering.

    ``None`` means the user kept the default Plot 6 setting. The default
    manuscript-oriented set is Ala/Arg/Asn/Asp/Cys/Gly/His/Ile/Leu/Phe/Pro/Ser/Thr/Tyr.
    An explicit list provided by the GUI or script is respected verbatim, so
    choosing 'Select all' can include any available amino-acid family again.
    """
    aas = _normalize_optional_string_list(include_aas)
    if aas is not None:
        return aas
    return list(PLOT6_DEFAULT_AAS)


def _format_plot6_amino_acid_subtitle(SET) -> str:
    """Return a concise subtitle describing amino-acid families used by Plot 6."""
    aas = _normalize_optional_string_list((SET or {}).get('trna_modifications_include_aas', None))
    if aas is None:
        return 'Amino acids used: ' + ', '.join(PLOT6_DEFAULT_AAS)
    clean = []
    seen = set()
    for aa in aas:
        val = str(aa or '').strip()
        if not val:
            continue
        key = val.lower()
        if key in seen:
            continue
        seen.add(key)
        clean.append(val)
    if not clean:
        return 'Amino acids used: all amino acid families'
    return 'Amino acids used: ' + ', '.join(clean)


def _codon_internal_aa(label):
    s = str(label or "").strip()
    if "_" not in s:
        return ""
    return s.split("_", 1)[0].strip()


def _safe_total_counts_for_selected_aas(count_df: pd.DataFrame, include_aas=None):
    aas = _apply_default_plot6_aa_exclusions(include_aas)
    keep = {str(a).strip().lower() for a in aas if str(a).strip()}
    cols = [c for c in count_df.columns if _codon_internal_aa(c).lower() in keep]
    if not cols:
        return pd.Series(np.nan, index=count_df.index, dtype=float)
    total = count_df.loc[:, cols].apply(pd.to_numeric, errors='coerce').fillna(0.0).sum(axis=1)
    return total.replace(0.0, np.nan)



def _modification_position_key(feature_label: str) -> str:
    """Return '32', '34', '37' or '' for anticodon-loop modification labels."""
    s = str(feature_label or '').strip()
    if not s:
        return ''
    # Modification labels are usually compact names such as s2C32, mnm5s2U34,
    # cmo5U34, i6A37, ms2i6A37. Keep this deliberately permissive so
    # user-selected aliases and mathtext display labels still group correctly.
    for pos in ('32', '34', '37'):
        if re.search(rf'(?<!\d){pos}(?!\d)', s):
            return pos
    return ''


def _format_plot6_modification_label_mathtext(label: str) -> str:
    """Format Plot 6 modification labels with chemical superscripts and position subscripts.

    Examples:
      mnm5s2U34 -> mnm$^{5}$s$^{2}$U$_{34}$
      s²C32     -> s$^{2}$C$_{32}$
      Cm32/Um32 -> Cm$_{32}$/Um$_{32}$
    """
    raw = str(label or '').strip()
    if not raw:
        return raw

    trans = str.maketrans({
        '⁰': '0', '¹': '1', '²': '2', '³': '3', '⁴': '4',
        '⁵': '5', '⁶': '6', '⁷': '7', '⁸': '8', '⁹': '9',
        '₀': '0', '₁': '1', '₂': '2', '₃': '3', '₄': '4',
        '₅': '5', '₆': '6', '₇': '7', '₈': '8', '₉': '9',
    })

    def _format_part(part: str) -> str:
        p = str(part or '').strip().translate(trans)
        if not p:
            return p
        m = re.search(r'(32|34|37)\s*$', p)
        if m:
            pos = m.group(1)
            prefix = p[:m.start()].strip()
        else:
            pos = ''
            prefix = p
        # Remaining numbers denote chemical positions (m1, s2, ac4, i6, etc.).
        prefix = re.sub(r'(\d+)', r'$^{\1}$', prefix)
        return prefix + (f'$_{{{pos}}}$' if pos else '')

    # Keep simple composite labels readable without merging the components.
    pieces = re.split(r'(/)', raw)
    return ''.join(_format_part(x) if x != '/' else '/' for x in pieces)




def _format_plot6_modification_label_for_axis(label: str) -> str:
    """Return Plot 6 modification labels, splitting a few long names over two lines."""
    raw = str(label or '').strip()
    if not raw:
        return raw
    key = re.sub(r'[^A-Za-z0-9]+', '', raw).lower()
    manual = {
        'cm32um32': 'Cm$_{32}$\nUm$_{32}$',
        'mcmo5u34': '(m)cmo$^{5}$\nU$_{34}$',
        'mnm5s2u34': 'mnm$^{5}$s$^{2}$\nU$_{34}$',
        'mnm5u34': 'mnm$^{5}$\nU$_{34}$',
    }
    if key in manual:
        return manual[key]
    return _format_plot6_modification_label_mathtext(raw)


def _plot6_modification_label_line_count(label: str) -> int:
    s = str(label or '')
    return max(1, s.count('\n') + 1)

def _order_modification_features_by_position(features, feature_display=None):
    """Order Plot 6 modification features: position 32, then 34, then 37, then other."""
    feature_display = feature_display or {}

    def _display(f):
        if isinstance(feature_display, dict):
            return str(feature_display.get(f, f))
        return str(f)

    def _sort_key(f):
        lab = _display(f)
        pos = _modification_position_key(lab) or _modification_position_key(f)
        rank = {'32': 0, '34': 1, '37': 2}.get(pos, 3)
        return (rank, lab.lower(), str(f).lower())

    return sorted([str(f) for f in ([] if features is None else list(features))], key=_sort_key)


def _build_modification_position_metadata(features, feature_display=None):
    """Build x-axis label/group metadata for Plot 6 modification features."""
    feature_display = feature_display or {}
    meta = {}
    for f in ([] if features is None else list(features)):
        f = str(f)
        label = str(feature_display.get(f, f) if isinstance(feature_display, dict) else f)
        pos = _modification_position_key(label) or _modification_position_key(f)
        top_label = _format_plot6_modification_label_for_axis(label)
        n_lines = _plot6_modification_label_line_count(top_label)
        row = {
            'top': top_label,
            # Keep the *top* of one-line and two-line modification labels aligned.
            # The second line should extend downward, not shift the whole label block lower.
            'top_y': -0.074,
            'top_color': 'black',
            'n_top_lines': n_lines,
        }
        if pos in {'32', '34', '37'}:
            row['group'] = f'position_{pos}_modifications'
            row['group_label'] = f'Position {pos}'
        else:
            row['group'] = ''
            row['group_label'] = ''
        meta[f] = row
    return meta


def _compute_modification_usage_z(count_df: pd.DataFrame, trna_rules: dict, feature_mode: str = 'modifications', assignment_model: str = 'conservative', selected_features=None, include_aas=None, exclude_default_modifications: bool = True):
    """Return per-gene z-scored usage of codons annotated with each modification or tRME.

    feature_mode='modifications' uses modification names.
    feature_mode='enzymes' uses tRNA-modification enzyme names.

    Plot 6 is intentionally full-table-only: modification/tRME codon
    assignments are read exclusively from ``modification_assignment_models``,
    which is generated from the workbook sheet named ``Decoding table (full)``.
    The compact pooled decoding table may still define tRNA usage/ZTU, but it is
    never used as a fallback source for Plot 6 modification or tRME assignment.

    assignment_model controls how ambiguous full-table decoders are handled:
      - 'conservative': a codon is assigned to a feature only if all compatible
        decoders listed for that codon carry that feature.
      - 'permissive': a codon is assigned when at least one compatible decoder
        carries that feature.
      - 'fractional': implemented internally in the reader but disabled for
        Plot 6 export; only conservative and permissive plots are generated.
    """
    if count_df is None or count_df.empty or not trna_rules:
        return pd.DataFrame(), {}

    feature_mode = str(feature_mode or 'modifications').strip().lower()
    assignment_model = str(assignment_model or 'conservative').strip().lower()
    using_modifications = feature_mode not in {'enzyme', 'enzymes', 'trme', 'trmes', 'trme enzymes'}
    if using_modifications:
        selected_features_list = _canonicalize_plot6_modification_feature_list(selected_features)
    else:
        selected_features_list = _normalize_optional_string_list(selected_features)
    selected_feature_keys = None if selected_features_list is None else {_plot6_feature_key_ascii(x) for x in selected_features_list if str(x).strip()}
    include_aas_list = _apply_default_plot6_aa_exclusions(include_aas)
    include_aa_keys = {str(x).strip().lower() for x in include_aas_list if str(x).strip()}
    model_maps = trna_rules.get('modification_assignment_models') or {}
    if assignment_model not in {'conservative', 'permissive', 'fractional'}:
        print(f"[WARN] Plot 6 assignment model '{assignment_model}' is disabled because Plot 6 now uses only the full decoder-level table. Use conservative or permissive.")
        return pd.DataFrame(), {}
    model_payload = model_maps.get(assignment_model) or {}

    if feature_mode in {'enzyme', 'enzymes', 'trme', 'trmes', 'trme enzymes'}:
        by_codon = model_payload.get('trmes_by_codon') or {}
        ylabel = 'tRNA modification enzymes'
    else:
        by_codon = model_payload.get('modifications_by_codon') or {}
        ylabel = 'tRNA modifications'

    feature_to_codon_weights = {}
    for codon, items in by_codon.items():
        codon = str(codon)
        if codon not in count_df.columns:
            continue
        codon_aa = _codon_internal_aa(codon)
        if codon_aa.lower() not in include_aa_keys:
            continue
        if isinstance(items, dict):
            iterable = list(items.items())
        else:
            iterable = [(item, 1.0) for item in list(items or [])]
        for item, weight in iterable:
            feat = str(item).strip()
            if using_modifications:
                feat = _canonicalize_plot6_modification_feature_label(feat)
            if not feat:
                continue
            if selected_feature_keys is None and using_modifications and bool(exclude_default_modifications) and _is_default_excluded_plot6_modification(feat):
                continue
            if selected_feature_keys is not None and _plot6_feature_key_ascii(feat) not in selected_feature_keys:
                continue
            try:
                w = float(weight)
            except Exception:
                w = 1.0
            if not np.isfinite(w) or w <= 0:
                continue
            w = float(max(0.0, min(1.0, w)))
            feature_to_codon_weights.setdefault(feat, {})
            feature_to_codon_weights[feat][codon] = feature_to_codon_weights[feat].get(codon, 0.0) + w
    if not feature_to_codon_weights:
        return pd.DataFrame(), {}

    total = _safe_total_counts_for_selected_aas(count_df, include_aas=include_aas_list)
    raw = pd.DataFrame(index=count_df.index)
    for feat, codon_weights in feature_to_codon_weights.items():
        vals = pd.Series(0.0, index=count_df.index, dtype=float)
        for codon, weight in codon_weights.items():
            if codon not in count_df.columns:
                continue
            vals = vals.add(pd.to_numeric(count_df[codon], errors='coerce').fillna(0.0) * float(weight), fill_value=0.0)
        raw[feat] = vals.div(total) * 100.0
    if raw.empty:
        return raw, {}

    mu = raw.mean(axis=0, skipna=True)
    sd = raw.std(axis=0, skipna=True, ddof=0).replace(0.0, np.nan)
    z = (raw - mu) / sd
    for col in raw.columns:
        if not np.isfinite(sd.get(col, np.nan)):
            z.loc[raw[col].notna(), col] = 0.0
    label_map = {str(col): str(col) for col in raw.columns}
    return z, label_map






def _safe_xlsx_sheet_name(name, used=None):
    """Return a valid unique Excel worksheet name derived from a cluster name."""
    used = used if used is not None else set()
    used_lc = {str(u).lower() for u in used}
    base = str(name or 'Cluster').strip() or 'Cluster'
    base = re.sub(r"[\[\]\:\*\?\/\\]", "_", base)
    base = base.replace("'", "")
    base = re.sub(r"\s+", " ", base).strip() or 'Cluster'
    base = base[:31]
    candidate = base
    i = 2
    while str(candidate).lower() in used_lc:
        suffix = f"_{i}"
        candidate = (base[:31 - len(suffix)] + suffix)[:31]
        i += 1
    used.add(candidate)
    return candidate


def _metadata_text_for_plot6_aas(include_aas):
    aas = _apply_default_plot6_aa_exclusions(include_aas)
    clean = []
    seen = set()
    for aa in aas:
        val = str(aa or '').strip()
        if not val:
            continue
        key = val.lower()
        if key in seen:
            continue
        seen.add(key)
        clean.append(val)
    if _normalize_optional_string_list(include_aas) is None:
        prefix = 'default Plot 6 amino-acid set'
    else:
        prefix = 'explicit Plot 6 amino-acid selection'
    if clean:
        return prefix + ': ' + ', '.join(clean)
    return prefix + ': none selected'


def _build_gene_group_membership(cluster_df: pd.DataFrame, genes=None):
    """Map each gene to all cluster columns in which it appears."""
    gene_set = None if genes is None else {str(g).strip() for g in list(genes or []) if str(g).strip()}
    membership = {g: [] for g in (gene_set or [])}
    if cluster_df is None or getattr(cluster_df, 'empty', True):
        return membership
    for cname in list(cluster_df.columns):
        cname_s = str(cname)
        try:
            vals = cluster_df[cname].replace({np.nan: ''}).astype(str).tolist()
        except Exception:
            vals = []
        for v in vals:
            gene = str(v or '').strip()
            if not gene or gene.lower() == 'nan':
                continue
            if gene_set is not None and gene not in gene_set:
                continue
            membership.setdefault(gene, [])
            if cname_s not in membership[gene]:
                membership[gene].append(cname_s)
    return membership


def _cluster_gene_lists_for_selected(values_df: pd.DataFrame, cluster_df: pd.DataFrame, selected_clusters):
    if selected_clusters is None:
        selected = [str(c) for c in cluster_df.columns]
    else:
        selected = [str(c) for c in list(selected_clusters or []) if str(c) in list(cluster_df.columns)]
    if not selected:
        return {}
    out = {}
    for cname in selected:
        genes = [str(v).strip() for v in cluster_df[cname].replace({np.nan: ''}).astype(str).tolist()]
        genes = [g for g in genes if g and g.lower() != 'nan' and g in values_df.index]
        out[cname] = genes
    return out


def _trna_modifications_output_dir(output_dir: str) -> str:
    """Return the dedicated Plot 6 tRNA-modification output folder.

    Standard pipeline runs pass the project output folder here, whereas the
    GUI's "plot from existing workbook" route historically passed the
    project's Figures folder.  In both cases, place Plot 6 tRNA-modification
    plots and tables in a sibling/child folder named exactly "tRNA modifications"
    rather than inside the generic Figures folder.
    """
    base = str(output_dir or '').strip() or '.'
    base_norm = os.path.basename(os.path.normpath(base)).strip().lower()
    if base_norm == 'figures':
        base = os.path.dirname(os.path.normpath(base)) or base
    return os.path.join(base, 'tRNA modifications')


def _export_trna_modification_enrichment_table(
    output_dir: str,
    count_df: pd.DataFrame,
    cluster_df: pd.DataFrame,
    selected_clusters,
    trna_rules: dict,
    assignment_models,
    include_aas=None,
    use_log2: bool = True,
):
    """Export per-gene Plot 6 tRNA-modification enrichment values.

    The export intentionally ignores the Plot 6 selected-feature list so that it
    contains all modification features compatible with the selected amino-acid
    set.  It still uses the same full-table assignment models and amino-acid
    filtering as the plotted tRNA-modification analysis.
    """
    if output_dir is None or count_df is None or count_df.empty or cluster_df is None or getattr(cluster_df, 'empty', True):
        return ''
    models = [str(m).strip().lower() for m in list(assignment_models or []) if str(m).strip().lower() in {'conservative', 'permissive'}]
    if not models:
        models = ['conservative', 'permissive']

    model_results = {}
    label_maps = {}
    shared_features = []
    shared_label_map = {}
    # Use permissive first to define the broadest feature order, then add any
    # features detected only in stricter models.
    scan_order = [m for m in ['permissive', 'conservative'] if m in models] + [m for m in models if m not in {'permissive', 'conservative'}]
    for model in models:
        z_df, label_map = _compute_modification_usage_z(
            count_df,
            trna_rules,
            feature_mode='modifications',
            assignment_model=model,
            selected_features=None,
            include_aas=include_aas,
            exclude_default_modifications=False,
        )
        model_results[model] = z_df
        label_maps[model] = label_map or {}
    for model in scan_order:
        z_df = model_results.get(model)
        label_map = label_maps.get(model, {})
        if z_df is None or z_df.empty:
            continue
        for feat in list(z_df.columns):
            feat = str(feat)
            if feat not in shared_features:
                shared_features.append(feat)
            if feat not in shared_label_map:
                shared_label_map[feat] = str(label_map.get(feat, feat))
    if not shared_features:
        print('[WARN] Plot 6 export: no tRNA-modification features available for the selected amino-acid set.')
        return ''
    shared_features = _order_modification_features_by_position(shared_features, shared_label_map)
    shared_label_map = {feat: str(shared_label_map.get(feat, feat)) for feat in shared_features}

    cluster_genes = _cluster_gene_lists_for_selected(count_df, cluster_df, selected_clusters)
    cluster_genes = {str(c): genes for c, genes in cluster_genes.items() if genes}
    if not cluster_genes:
        print('[WARN] Plot 6 export: no selected cluster contains genes present in the codon-count table.')
        return ''

    all_genes = [str(g) for g in list(count_df.index)]
    gene_to_groups = _build_gene_group_membership(cluster_df, genes=all_genes)

    # Prepare model matrices once for the whole genome.  These are then reused
    # for the master genome sheet and the per-cluster sheets.
    export_model_tables = {}
    for model in models:
        z_df = model_results.get(model)
        if z_df is None or z_df.empty:
            model_table = pd.DataFrame(0.0, index=count_df.index, columns=shared_features)
        else:
            model_table = z_df.copy()
            for feat in shared_features:
                if feat not in model_table.columns:
                    model_table[feat] = 0.0
            model_table = model_table.reindex(index=count_df.index, columns=shared_features)
        if bool(use_log2):
            model_table = model_table.apply(lambda col: pd.Series(_signed_log2_values(col.values), index=col.index), axis=0)
        export_model_tables[model] = model_table

    def _build_export_data(genes):
        genes = [str(g) for g in list(genes or []) if str(g) in count_df.index]
        data = pd.DataFrame({
            'group': [', '.join(gene_to_groups.get(g, [])) for g in genes],
            'gene': genes,
        })
        for model in models:
            model_table = export_model_tables.get(model)
            if model_table is None or model_table.empty:
                continue
            for feat in shared_features:
                display = shared_label_map.get(feat, feat)
                col = str(display) if len(models) == 1 else f'{display} [{model}]'
                data[col] = pd.to_numeric(model_table.loc[genes, feat], errors='coerce').astype(float).values
        return data

    out_path = os.path.join(output_dir, 'tRNA modification enrichment table.xlsx')
    os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
    aa_text = _metadata_text_for_plot6_aas(include_aas)
    value_desc = 'signed log2-transformed z-score enrichment values: sign(z) × log2(|z| + 1)' if bool(use_log2) else 'genome-normalized z-score enrichment values'
    selected_clusters_text = ', '.join(str(c) for c in cluster_genes.keys())
    method_text = (
        'For each gene and modification, codons associated with that modification in the full decoder-level decoding table '
        'are summed within the selected amino-acid families, divided by the total selected-amino-acid codon count for that gene, multiplied by 100, '
        'and z-scored across all genes in the input dataset. Conservative = all compatible decoders carry the modification; permissive = at least one compatible decoder carries it. '
        f'Table values are {value_desc}; all tRNA modifications compatible with the selected amino-acid set are exported regardless of Plot 6 feature selection.'
    )
    try:
        with pd.ExcelWriter(out_path, engine='xlsxwriter') as writer:
            used_sheet_names = set()
            workbook = writer.book
            header_fmt = workbook.add_format({'bold': True, 'bg_color': '#D9EAF7', 'border': 1, 'text_wrap': True, 'valign': 'top'})
            number_fmt = workbook.add_format({'num_format': '0.000', 'border': 0})
            text_fmt = workbook.add_format({'border': 0, 'valign': 'top'})
            meta_key_fmt = workbook.add_format({'bold': True, 'bg_color': '#EAF2F8', 'border': 1, 'text_wrap': True, 'valign': 'top'})
            meta_val_fmt = workbook.add_format({'border': 1, 'text_wrap': True, 'valign': 'top'})

            # Metadata sheet contains all explanations, keeping data sheets clean.
            metadata_rows = [
                ('Export', 'Plot 6 per-gene tRNA modification enrichment table'),
                ('Amino acids used', aa_text),
                ('Selected clusters', selected_clusters_text),
                ('Assignment models exported', ', '.join(models)),
                ('Value scale', value_desc),
                ('Data sheets', 'decoding summary lists codons assigned to each modification in conservative/permissive regimes; genome contains all genes in the input codon-count table; one additional sheet is exported for each selected cluster.'),
                ('group column', 'Comma-separated list of all cluster columns in the input cluster table containing the gene; empty means the gene was not found in any cluster column.'),
                ('Modification columns', 'All modifications compatible with the selected amino-acid set are exported, including modifications not selected for plotting.'),
                ('Computation', method_text),
            ]
            meta_df = pd.DataFrame(metadata_rows, columns=['field', 'value'])
            meta_df.to_excel(writer, sheet_name='metadata', index=False)
            used_sheet_names.add('metadata')
            ws_meta = writer.sheets['metadata']
            ws_meta.set_column(0, 0, 26, meta_key_fmt)
            ws_meta.set_column(1, 1, 120, meta_val_fmt)
            ws_meta.set_row(0, 24)
            for col_idx, col_name in enumerate(meta_df.columns):
                ws_meta.write(0, col_idx, col_name, header_fmt)
            for r in range(1, len(meta_df) + 1):
                ws_meta.set_row(r, 42 if r in {2, 6, 7, 8, 9} else 28)
            ws_meta.freeze_panes(1, 0)

            def _write_data_sheet(sheet_name, data):
                data.to_excel(writer, sheet_name=sheet_name, index=False)
                ws = writer.sheets[sheet_name]
                ncols = max(1, len(data.columns))
                ws.set_row(0, 36)
                for col_idx, col_name in enumerate(data.columns):
                    ws.write(0, col_idx, col_name, header_fmt)
                    if col_idx == 0:
                        max_len = max([len(str(col_name))] + [len(str(x)) for x in data.iloc[:500, col_idx].fillna('').astype(str).tolist()]) if not data.empty else len(str(col_name))
                        ws.set_column(col_idx, col_idx, min(max(12, max_len + 2), 48), text_fmt)
                    elif col_idx == 1:
                        max_len = max([len(str(col_name))] + [len(str(x)) for x in data.iloc[:500, col_idx].fillna('').astype(str).tolist()]) if not data.empty else len(str(col_name))
                        ws.set_column(col_idx, col_idx, min(max(14, max_len + 2), 36), text_fmt)
                    else:
                        ws.set_column(col_idx, col_idx, 14, number_fmt)
                ws.freeze_panes(1, 2)
                if len(data) > 0 and ncols > 0:
                    ws.autofilter(0, 0, len(data), ncols - 1)

            def _codon_display_for_summary(codon_internal):
                s_codon = str(codon_internal or '').strip()
                if '_' in s_codon:
                    aa, cod = s_codon.split('_', 1)
                else:
                    aa, cod = '', s_codon
                cod = str(cod or '').upper().replace('T', 'U')
                return f'{aa}-{cod}' if aa and cod else s_codon

            def _build_decoding_summary_df(model):
                model_payload = (trna_rules.get('modification_assignment_models') or {}).get(model) or {}
                by_codon = model_payload.get('modifications_by_codon') or {}
                include_aas_list = _apply_default_plot6_aa_exclusions(include_aas)
                include_aa_keys = {str(x).strip().lower() for x in include_aas_list if str(x).strip()}
                feature_to_codons = {feat: [] for feat in shared_features}
                feature_key_to_feat = {_plot6_feature_key_ascii(feat): feat for feat in shared_features}
                for codon in list(count_df.columns):
                    codon_s = str(codon)
                    if codon_s not in by_codon:
                        continue
                    codon_aa = _codon_internal_aa(codon_s)
                    if include_aa_keys and codon_aa.lower() not in include_aa_keys:
                        continue
                    items = by_codon.get(codon_s) or []
                    iterable = items.keys() if isinstance(items, dict) else list(items or [])
                    for item in iterable:
                        feat = _canonicalize_plot6_modification_feature_label(str(item).strip())
                        mapped = feature_key_to_feat.get(_plot6_feature_key_ascii(feat))
                        if mapped is None:
                            continue
                        disp_codon = _codon_display_for_summary(codon_s)
                        if disp_codon not in feature_to_codons[mapped]:
                            feature_to_codons[mapped].append(disp_codon)
                display_cols = [shared_label_map.get(feat, feat) for feat in shared_features]
                max_len = max([len(v) for v in feature_to_codons.values()] + [0])
                rows = []
                for i in range(max_len):
                    row = []
                    for feat in shared_features:
                        vals = feature_to_codons.get(feat, [])
                        row.append(vals[i] if i < len(vals) else '')
                    rows.append(row)
                return pd.DataFrame(rows, columns=display_cols)

            def _write_decoding_summary_sheet():
                sheet_name = _safe_xlsx_sheet_name('decoding summary', used_sheet_names)
                ws = workbook.add_worksheet(sheet_name)
                used_sheet_names.add(sheet_name)
                writer.sheets[sheet_name] = ws
                title_fmt = workbook.add_format({'bold': True, 'bg_color': '#1F4E79', 'font_color': '#FFFFFF', 'border': 1})
                note_fmt = workbook.add_format({'italic': True, 'text_wrap': True, 'valign': 'top'})
                ws.write(0, 0, 'Decoding summary', title_fmt)
                ws.write(1, 0, 'Each block lists the codons assigned to each grouped mature tRNA modification under the indicated decoder-assignment regime. Only the explicit Modifications position 32/34/37 columns are used; blank cells are treated as no modification.', note_fmt)
                ws.set_row(1, 42)
                start_row = 3
                for model in [m for m in ['conservative', 'permissive'] if m in models]:
                    summary_df = _build_decoding_summary_df(model)
                    ws.write(start_row, 0, f'{model} regime', title_fmt)
                    start_row += 1
                    for col_idx, col_name in enumerate(summary_df.columns):
                        ws.write(start_row, col_idx, col_name, header_fmt)
                        ws.set_column(col_idx, col_idx, 18, text_fmt)
                    for r_idx, row_vals in enumerate(summary_df.values.tolist(), start=start_row + 1):
                        for c_idx, value in enumerate(row_vals):
                            ws.write(r_idx, c_idx, value, text_fmt)
                    if summary_df.empty:
                        ws.write(start_row + 1, 0, '(no codons assigned under this regime)', text_fmt)
                        start_row += 3
                    else:
                        start_row += len(summary_df) + 3
                ws.freeze_panes(4, 0)

            _write_decoding_summary_sheet()

            genome_sheet = _safe_xlsx_sheet_name('genome', used_sheet_names)
            _write_data_sheet(genome_sheet, _build_export_data(all_genes))

            for cluster_name, genes in cluster_genes.items():
                genes = [g for g in genes if g in count_df.index]
                if not genes:
                    continue
                sheet_name = _safe_xlsx_sheet_name(cluster_name, used_sheet_names)
                _write_data_sheet(sheet_name, _build_export_data(genes))

        print(f'[INFO] Exported Plot 6 per-gene tRNA modification enrichment table → {out_path}')
        return out_path
    except Exception as e:
        print(f'[WARN] Could not export Plot 6 tRNA modification enrichment table: {e}')
        return ''

def _setting_bool(SET, key, default=False):
    v = SET.get(key, default)
    if isinstance(v, str):
        return v.strip().lower() in {'1', 'true', 'yes', 'y', 'on'}
    return bool(v)


def _signed_log2_values(vals):
    arr = np.asarray(vals, dtype=float)
    return np.sign(arr) * np.log2(np.abs(arr) + 1.0)


def _filter_outliers_mean_sd(vals, n_sd=3.0):
    arr = np.asarray(vals, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size < 3:
        return arr
    mu = float(np.nanmean(arr))
    sd = float(np.nanstd(arr, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        return arr
    n_sd = max(0.1, float(n_sd))
    keep = (arr >= mu - n_sd * sd) & (arr <= mu + n_sd * sd)
    return arr[keep]


def _draw_mean_sd_box(ax, x, vals, width, color):
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return
    mean = float(np.nanmean(vals))
    sd = float(np.nanstd(vals, ddof=1)) if vals.size > 1 else 0.0
    y0 = mean - sd
    height = 2.0 * sd
    if not np.isfinite(height) or height <= 0:
        height = 1e-9
    rect = Rectangle((x - width / 2.0, y0), width, height,
                     facecolor=color, edgecolor='black', linewidth=0.8, alpha=0.68, zorder=2)
    ax.add_patch(rect)
    ax.plot([x - width / 2.0, x + width / 2.0], [mean, mean], color='black', lw=1.0, zorder=3)



def _p_to_stars(p):
    try:
        p = float(p)
    except Exception:
        return 'ns'
    if not np.isfinite(p):
        return 'ns'
    if p < 0.001:
        return '***'
    if p < 0.01:
        return '**'
    if p < 0.05:
        return '*'
    return 'ns'


def _compute_two_group_test(a, b, method='Student t-test'):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    method_in = str(method or 'none').strip()
    method_l = method_in.lower()
    if method_l in {'', 'none', 'no', 'false', 'off'}:
        return np.nan, np.nan, 'none'
    if a.size < 1 or b.size < 1:
        return np.nan, np.nan, method_in
    try:
        if method_l in {'mann-whitney u', 'mann whitney u', 'mann-whitney', 'mann whitney', 'nonparametric', 'non-parametric'}:
            if a.size < 1 or b.size < 1:
                return np.nan, np.nan, 'Mann-Whitney U'
            stat, pval = mannwhitneyu(a, b, alternative='two-sided')
            return float(stat), float(pval), 'Mann-Whitney U'
        equal_var = not (method_l in {'welch t-test', 'welch', 'welch ttest'})
        if a.size < 2 or b.size < 2:
            return np.nan, np.nan, 'Welch t-test' if not equal_var else 'Student t-test'
        with np.errstate(all='ignore'):
            stat, pval = ttest_ind(a, b, equal_var=equal_var, nan_policy='omit')
        return float(stat), float(pval), 'Welch t-test' if not equal_var else 'Student t-test'
    except Exception:
        return np.nan, np.nan, method_in


def _export_decoding_stats_excel(rows, out_path):
    if not rows:
        return ''
    try:
        df = pd.DataFrame(rows)
        os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
        with pd.ExcelWriter(out_path, engine='xlsxwriter') as writer:
            df.to_excel(writer, sheet_name='statistics', index=False)
        return out_path
    except Exception as e:
        print(f'[WARN] Could not export decoding statistics: {e}')
        return ''

def _plot_cluster_feature_boxplots(
    SET,
    values_df: pd.DataFrame,
    cluster_df: pd.DataFrame,
    selected_clusters,
    feature_order,
    feature_display,
    out_path: str,
    title_txt: str,
    y_label: str = 'Enrichment vs genome (z-score)',
    secondary_values: pd.Series | None = None,
    secondary_std: pd.Series | None = None,
    secondary_label: str = '',
    show_fig: bool = True,
    plot_key: str = 'shift',
    feature_metadata: dict | None = None,
):
    if values_df is None or values_df.empty or cluster_df is None or getattr(cluster_df, 'empty', True):
        return ''
    features = [str(f) for f in list(feature_order or []) if str(f) in values_df.columns]
    if not features:
        return ''
    display = [str(feature_display.get(f, f) if isinstance(feature_display, dict) else feature_display[i] if feature_display is not None and i < len(feature_display) else f) for i, f in enumerate(features)]
    feature_metadata = feature_metadata or {}

    cluster_genes = _cluster_gene_lists_for_selected(values_df, cluster_df, selected_clusters)
    clusters = [c for c, genes in cluster_genes.items() if genes]
    if not clusters:
        return ''

    plot_key = str(plot_key or 'shift').strip().lower()
    is_modification_plot = plot_key in {'modifications', 'modification'}
    if is_modification_plot:
        # Plot 6 figure-detail settings are used as the default behavior.
        # The GUI checkbox remains for backward compatibility, but the preset
        # values are always honored whenever they are provided.
        use_custom = True
        dpi = _optional_int_value(SET.get('trna_modification_plots_dpi'))
        fig_w = _positive_optional_float_value(SET.get('trna_modification_plots_fig_width'))
        fig_h = _positive_optional_float_value(SET.get('trna_modification_plots_fig_height'))
        xtick_fs = _optional_int_value(SET.get('trna_modification_plots_caption_size'))
        ytick_fs = xtick_fs
        title_fs = xtick_fs
    else:
        use_custom = bool(SET.get('trna_shift_heatmaps_customize', False))
        dpi = _optional_int_value(SET.get('trna_shift_heatmaps_dpi')) if use_custom else None
        fig_w = _positive_optional_float_value(SET.get('trna_shift_heatmaps_fig_width')) if use_custom else None
        fig_h = _positive_optional_float_value(SET.get('trna_shift_heatmaps_fig_height')) if use_custom else None
        xtick_fs = _optional_int_value(SET.get('trna_shift_heatmaps_xtick_fontsize')) if use_custom else None
        ytick_fs = _optional_int_value(SET.get('trna_shift_heatmaps_ytick_fontsize')) if use_custom else None
        title_fs = _optional_int_value(SET.get('trna_shift_heatmaps_title_fontsize')) if use_custom else None

    key_prefix = {
        'wobble': 'trna_wobble',
        'shift': 'trna_shift',
        'modifications': 'trna_modifications',
        'modification': 'trna_modifications',
    }.get(plot_key, 'trna_shift')
    style = str(SET.get(f'{key_prefix}_boxplot_style', SET.get(f'{key_prefix}_plot_kind', 'boxplot')) or 'boxplot').strip().lower()
    if style not in {'boxplot', 'violin'}:
        style = 'boxplot'
    use_log2 = _setting_bool(SET, f'{key_prefix}_boxplot_log2', True)
    exclude_outliers = _setting_bool(SET, f'{key_prefix}_exclude_outliers', False)
    outlier_sd = float(SET.get(f'{key_prefix}_outlier_sd', 3.0) or 3.0)
    if is_modification_plot:
        y_min = _optional_float_value(SET.get('trna_modification_plots_ymin'))
        y_max = _optional_float_value(SET.get('trna_modification_plots_ymax'))
        if y_min is None:
            y_min = _optional_float_value(SET.get(f'{key_prefix}_boxplot_ymin'))
        if y_max is None:
            y_max = _optional_float_value(SET.get(f'{key_prefix}_boxplot_ymax'))
    else:
        y_min = _optional_float_value(SET.get(f'{key_prefix}_boxplot_ymin'))
        y_max = _optional_float_value(SET.get(f'{key_prefix}_boxplot_ymax'))
    if is_modification_plot and SET.get('trna_modification_plots_caption_size') is not None:
        caption_size = int(SET.get('trna_modification_plots_caption_size') or 17)
    else:
        caption_size = int(SET.get(f'{key_prefix}_boxplot_caption_size', 17 if is_modification_plot else 13) or (17 if is_modification_plot else 13))
    caption_size = max(6, caption_size)
    stats_test = str(SET.get(f'{key_prefix}_stats_test', 'none') or 'none').strip()
    pair_stats_test = str(SET.get(f'{key_prefix}_pair_stats_test', 'none') or 'none').strip()

    # User-adjustable vertical clearance between unbracketed
    # cluster-vs-reference stars and bracketed within-feature-pair
    # comparisons.  Values are fractions of the current y-axis span.
    default_pair_gap = 0.30 if plot_key == 'shift' else 0.26
    pair_stats_gap_frac = _optional_float_value(SET.get(f'{key_prefix}_pair_stats_gap', default_pair_gap))
    if pair_stats_gap_frac is None or not np.isfinite(pair_stats_gap_frac):
        pair_stats_gap_frac = default_pair_gap
    pair_stats_gap_frac = float(np.clip(pair_stats_gap_frac, 0.02, 1.50))

    reference_cluster = str(SET.get('decoding_reference_cluster', '') or '').strip()

    dpi = int(dpi if dpi is not None else SET.get('figure_dpi', 300))
    n_features = len(features)
    n_clusters = len(clusters)

    group_values = [str((feature_metadata.get(f, {}) or {}).get('group', '')) for f in features]
    has_groups = any(group_values)
    x_positions = []
    pos = 0.0
    for i in range(n_features):
        if i == 0:
            pos = 0.0
        else:
            prev = group_values[i - 1] if has_groups else ''
            cur = group_values[i] if has_groups else ''
            pos += 1.0 + (1.05 if has_groups and cur and prev and cur != prev else (0.30 if plot_key in {'modifications', 'modification'} else 0.0))
        x_positions.append(pos)
    xbase = np.asarray(x_positions, dtype=float)

    n_groups = len([1 for i, g in enumerate(group_values) if g and (i == 0 or g != group_values[i-1])]) if has_groups else n_features
    if is_modification_plot:
        # Plot 6 needs extra room for multi-line modification labels and
        # position-group captions below the x-axis, plus the title/legend above.
        default_fig_w = max(8.6, 0.78 * (xbase[-1] + 1 if len(xbase) else n_features) + 3.8)
        default_fig_h = max(8.0, 5.9 + 0.12 * n_clusters)
    else:
        default_fig_w = max(8.2, 0.78 * (xbase[-1] + 1 if len(xbase) else n_features) + 3.8)
        default_fig_h = max(5.3, 4.1 + 0.08 * n_clusters)
    fig_w = float(fig_w if fig_w is not None else default_fig_w)
    fig_h = float(fig_h if fig_h is not None else default_fig_h)
    # A single per-plot caption size controls title, labels, ticks, legend, and
    # Plot 6 position-group captions.
    xtick_fs = caption_size
    ytick_fs = caption_size
    title_fs = caption_size
    if is_modification_plot:
        # Plot 6 gets a slightly wider default box width than the other tRNA
        # boxplots.  A dedicated Figure details control can override it.
        raw_box_width = SET.get('trna_modification_plots_box_width', None)
        box_width = float(raw_box_width if raw_box_width not in (None, '') else 0.18)
    else:
        box_width = float(SET.get('trna_boxplot_width', 0.12) or 0.12)
    box_width = max(0.04, min(0.8 / max(1, n_clusters), box_width))
    if is_modification_plot and n_clusters > 1:
        # Keep approximately the same empty space between neighboring cluster
        # boxes as before, while allowing the boxes themselves to be 25% wider.
        previous_width = 0.12
        previous_span = min(0.84, max(0.24, n_clusters * previous_width * 1.18))
        previous_spacing = previous_span / float(max(1, n_clusters - 1))
        previous_edge_gap = max(0.015, previous_spacing - previous_width)
        group_span = min(0.84, max(0.24, (n_clusters - 1) * (box_width + previous_edge_gap)))
    else:
        group_span = min(0.84, max(0.24, n_clusters * box_width * 1.18))
    offsets = np.linspace(-group_span / 2.0, group_span / 2.0, n_clusters) if n_clusters > 1 else np.array([0.0])
    colors = _adaptive_cluster_colors(n_clusters)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)
    ax2 = None
    secondary_style = str(SET.get('trna_secondary_axis_style', 'bars') or 'bars').strip().lower()
    if secondary_values is not None and secondary_style != 'none':
        sec_vals = pd.to_numeric(pd.Series(secondary_values), errors='coerce').reindex(features)
        if sec_vals.notna().any():
            ax2 = ax.twinx()
            ax2.set_zorder(ax.get_zorder() - 1)
            ax.set_zorder(ax2.get_zorder() + 1)
            ax.patch.set_visible(False)
            alpha = max(0.0, min(1.0, float(SET.get('trna_secondary_axis_alpha', 0.22) or 0.22)))
            bw = max(0.05, min(1.0, float(SET.get('trna_secondary_axis_bar_width', 0.72) or 0.72)))
            sec_y = sec_vals.to_numpy(dtype=float)
            if secondary_style == 'dots':
                ax2.scatter(xbase, sec_y, s=36, color='grey', alpha=alpha, zorder=0)
                if secondary_std is not None:
                    sec_std = pd.to_numeric(pd.Series(secondary_std), errors='coerce').reindex(features).to_numpy(dtype=float)
                    ax2.errorbar(xbase, sec_y, yerr=sec_std, fmt='none', ecolor='grey', alpha=min(0.8, alpha + 0.25), capsize=3, lw=0.8, zorder=0)
            elif secondary_style == 'line':
                ax2.plot(xbase, sec_y, color='grey', alpha=alpha, marker='o', lw=1.3, zorder=0)
                if secondary_std is not None:
                    sec_std = pd.to_numeric(pd.Series(secondary_std), errors='coerce').reindex(features).to_numpy(dtype=float)
                    ax2.errorbar(xbase, sec_y, yerr=sec_std, fmt='none', ecolor='grey', alpha=min(0.8, alpha + 0.25), capsize=3, lw=0.8, zorder=0)
            else:
                yerr = None
                if secondary_std is not None:
                    yerr = pd.to_numeric(pd.Series(secondary_std), errors='coerce').reindex(features).to_numpy(dtype=float)
                ax2.bar(xbase, sec_y, width=bw, color='grey', alpha=alpha, yerr=yerr,
                        error_kw={'ecolor': 'grey', 'lw': 0.8, 'capsize': 3, 'alpha': min(0.8, alpha + 0.25)}, zorder=0)
            ax2.set_ylabel(str(secondary_label or ''), fontsize=caption_size)
            ax2.tick_params(axis='y', labelsize=caption_size)
            ax2.grid(False)

    rng = np.random.default_rng(12345)
    per_cluster_n = {}
    all_plot_values = []
    plotted_values = {}
    for j, cname in enumerate(clusters):
        color = colors[j]
        genes = cluster_genes[cname]
        per_cluster_n[cname] = len(genes)
        data = []
        positions = xbase + offsets[j]
        for f in features:
            vals = pd.to_numeric(values_df.loc[genes, f], errors='coerce').dropna().to_numpy(dtype=float)
            if exclude_outliers:
                vals = _filter_outliers_mean_sd(vals, n_sd=outlier_sd)
            if use_log2:
                vals = _signed_log2_values(vals)
            plotted_values[(str(cname), str(f))] = vals.copy()
            data.append(vals)
            if vals.size:
                all_plot_values.extend(vals.tolist())

        if style == 'violin':
            for pos_i, vals in zip(positions, data):
                if vals.size == 0:
                    continue
                parts = ax.violinplot([vals], positions=[pos_i], widths=box_width * 1.75,
                                      showmeans=False, showmedians=False, showextrema=False)
                for body in parts.get('bodies', []):
                    body.set_facecolor(color)
                    body.set_edgecolor('black')
                    body.set_alpha(0.42)
                    body.set_linewidth(0.75)
                mean = float(np.nanmean(vals))
                sd = float(np.nanstd(vals, ddof=1)) if vals.size > 1 else 0.0
                ax.plot([pos_i - box_width / 2.0, pos_i + box_width / 2.0], [mean, mean], color='black', lw=1.0, zorder=4)
                if vals.size > 1 and np.isfinite(sd):
                    ax.plot([pos_i, pos_i], [mean - sd, mean + sd], color='black', lw=0.8, zorder=4)
        else:
            for pos_i, vals in zip(positions, data):
                _draw_mean_sd_box(ax, pos_i, vals, box_width, color)

        if bool(SET.get('trna_boxplot_show_points', True)):
            pt_alpha = max(0.0, min(1.0, float(SET.get('trna_boxplot_point_alpha', 0.35) or 0.35)))
            pt_size = max(1.0, float(SET.get('trna_boxplot_point_size', 10.5) or 10.5))
            for pos_i, vals in zip(positions, data):
                if vals.size == 0:
                    continue
                jitter = rng.normal(0, box_width * 0.13, size=vals.size)
                ax.scatter(np.full(vals.size, pos_i) + jitter, vals, s=pt_size, color=color,
                           alpha=pt_alpha, edgecolors='none', zorder=3)

    if not any(np.isfinite(np.asarray(all_plot_values, dtype=float))):
        print(f'[WARN] {title_txt}: no finite cluster enrichment values were found. The figure was not saved.')
        plt.close(fig)
        return ''

    ax.axhline(0, color='black', lw=0.8, alpha=0.45)
    ax.set_xticks(xbase)
    if is_modification_plot and n_features > 0:
        # Extra right-side data margin keeps the final position-group cap visible.
        ax.set_xlim(float(xbase[0]) - 0.55, float(xbase[-1]) + 0.85)

    if feature_metadata:
        ax.set_xticklabels([''] * n_features)
        trans = ax.get_xaxis_transform()
        top_y_values = []
        top_text_artists = []
        top_label_fs = max(6, caption_size - 2) if is_modification_plot else caption_size
        for i, f in enumerate(features):
            meta = feature_metadata.get(f, {}) or {}
            top = str(meta.get('top', display[i]) or display[i])
            top_color = str(meta.get('top_color', meta.get('bottom_color', 'black')) or 'black')
            top_y = _optional_float_value(meta.get('top_y', -0.070))
            if top_y is None or not np.isfinite(top_y):
                top_y = -0.070
            top_y_values.append(float(top_y))
            txt = ax.text(xbase[i], float(top_y), top, ha='center', va='top', transform=trans,
                          fontsize=top_label_fs, color=top_color, clip_on=False, linespacing=(0.95 if is_modification_plot else 1.0))
            top_text_artists.append(txt)

        # Use actual rendered x-label extents to draw amino-acid group bars. This
        # makes the bottom bars cover the full width of long labels such as Gly,
        # Leu and Ser anticodon groups instead of relying on a fixed half-width.
        label_x_bounds = []
        try:
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            inv_data = ax.transData.inverted()
            for txt, xb in zip(top_text_artists, xbase):
                bbox = txt.get_window_extent(renderer=renderer)
                x0 = inv_data.transform((bbox.x0, bbox.y0))[0]
                x1 = inv_data.transform((bbox.x1, bbox.y1))[0]
                if not (np.isfinite(x0) and np.isfinite(x1)):
                    raise ValueError('non-finite label extent')
                label_x_bounds.append((min(x0, x1), max(x0, x1)))
        except Exception:
            label_x_bounds = [(float(x) - 0.42, float(x) + 0.42) for x in xbase]

        min_top_y = min(top_y_values) if top_y_values else -0.070
        # Keep the group bar visibly below the lowest anticodon/tRNA label. The
        # larger offset is especially useful when Leu/Ser labels are staggered
        # on two rows.
        if is_modification_plot:
            custom_group_y = _optional_float_value(SET.get('trna_modification_plots_group_bar_y'))
            custom_group_gap = _optional_float_value(SET.get('trna_modification_plots_group_label_gap'))
            group_bracket_y = float(custom_group_y) if custom_group_y is not None else min(-0.255, float(min_top_y) - 0.135)
            group_label_y = group_bracket_y - float(custom_group_gap if custom_group_gap is not None else 0.070)
        else:
            group_bracket_y = min(-0.175, float(min_top_y) - 0.085)
            group_label_y = group_bracket_y - 0.070
        if has_groups:
            start = 0
            for i in range(1, n_features + 1):
                if i == n_features or group_values[i] != group_values[start]:
                    grp = group_values[start]
                    if grp:
                        # Span the full rendered label width and add a small
                        # proportional margin so the bar is centered under the
                        # whole amino-acid group.
                        x0 = min(label_x_bounds[k][0] for k in range(start, i))
                        x1 = max(label_x_bounds[k][1] for k in range(start, i))
                        group_width = max(0.1, float(x1 - x0))
                        pad = max(0.10, 0.06 * group_width)
                        x0 -= pad
                        x1 += pad
                        yline = group_bracket_y
                        cap = 0.026
                        ax.plot([x0, x1], [yline, yline], transform=trans, color='black', lw=1.20, clip_on=False)
                        ax.plot([x0, x0], [yline - cap, yline + cap], transform=trans, color='black', lw=1.20, clip_on=False)
                        ax.plot([x1, x1], [yline - cap, yline + cap], transform=trans, color='black', lw=1.20, clip_on=False)
                        group_label = str((feature_metadata.get(features[start], {}) or {}).get('group_label', grp))
                        group_label_weight = 'normal' if is_modification_plot else 'bold'
                        group_fs = max(6, caption_size - 1) if is_modification_plot else caption_size
                        ax.text(0.5 * (x0 + x1), group_label_y, group_label, ha='center', va='top', transform=trans,
                                fontsize=group_fs, fontweight=group_label_weight, clip_on=False)
                    start = i
        if is_modification_plot:
            bottom_margin = 0.39 if min_top_y < -0.120 else (0.36 if min_top_y < -0.090 else 0.33)
            fig.subplots_adjust(bottom=bottom_margin, right=0.91, top=0.78)
        else:
            bottom_margin = 0.41 if min_top_y < -0.120 else (0.38 if min_top_y < -0.090 else 0.33)
            fig.subplots_adjust(bottom=bottom_margin, right=0.94, top=0.82)
    else:
        ax.set_xticklabels([wrap_label_no_break(d, 13) for d in display], rotation=0, fontsize=caption_size)
        fig.subplots_adjust(bottom=0.20, right=0.94, top=0.82)

    ylabel = y_label
    if use_log2:
        if plot_key == 'shift':
            ylabel = 'log$_2$(ZTU)'
        elif plot_key == 'wobble':
            ylabel = 'log$_2$(ZCU)'
        elif plot_key in {'modifications', 'modification'}:
            ylabel = 'tRNA modification enrichment\n(signed log2 z-score)'
        else:
            ylabel = f'log$_2$({ylabel})'
    if is_modification_plot and use_log2:
        ylabel = 'tRNA modification enrichment\n(signed log2 z-score)'
    ax.set_ylabel(ylabel, fontsize=caption_size)
    # Most Decoding strategies figures suppress titles because the legend sits
    # above the graph. Plot 6 uses the title to show the assignment model
    # (conservative / permissive / estimated_fraction).
    if is_modification_plot:
        # Plot 6 title is placed as a figure-level title above the legend.
        ax.set_title('')
    else:
        ax.set_title('')
    ax.tick_params(axis='both', labelsize=caption_size)
    if y_min is not None or y_max is not None:
        ax.set_ylim(bottom=y_min if y_min is not None else ax.get_ylim()[0], top=y_max if y_max is not None else ax.get_ylim()[1])
    elif all_plot_values:
        finite_vals = np.asarray([v for v in all_plot_values if np.isfinite(v)], dtype=float)
        if finite_vals.size:
            lo, hi = float(np.nanpercentile(finite_vals, 1)), float(np.nanpercentile(finite_vals, 99))
            span = hi - lo
            if np.isfinite(span) and span > 0:
                ax.set_ylim(lo - 0.20 * span, hi + 0.25 * span)
    # Statistical comparisons.
    # 1) Cluster-vs-reference stars are drawn without horizontal bars.
    # 2) Within-group feature-pair stars are drawn above horizontal bars, so the
    #    two statistical layers remain visually distinguishable when both are enabled.
    stats_rows = []
    annotation_top_by_feature = {str(f): -np.inf for f in features}
    ref_stats_enabled = reference_cluster and str(stats_test).strip().lower() not in {'', 'none', 'no', 'false', 'off'}
    pair_stats_enabled = str(pair_stats_test).strip().lower() not in {'', 'none', 'no', 'false', 'off'} and plot_key in {'wobble', 'shift'}

    if ref_stats_enabled and reference_cluster in clusters:
        ylo, yhi = ax.get_ylim()
        span = yhi - ylo if np.isfinite(yhi - ylo) and (yhi - ylo) > 0 else 1.0
        fixed_plot6_star_rows = plot_key in {'modifications', 'modification'}
        star_y_by_cluster = {}
        if fixed_plot6_star_rows:
            # Add some headroom once, then place all cluster-vs-reference stars
            # on fixed horizontal rows just inside the upper frame. This keeps
            # Plot 6 significance annotations vertically aligned across all
            # modification features instead of following each local data maximum.
            if y_max is None:
                n_comp = max(1, len([c for c in clusters if str(c) != str(reference_cluster)]))
                ax.set_ylim(top=yhi + max(0.24, 0.065 * n_comp) * span)
                ylo, yhi = ax.get_ylim()
                span = yhi - ylo if np.isfinite(yhi - ylo) and (yhi - ylo) > 0 else span
            comp_clusters = [str(c) for c in clusters if str(c) != str(reference_cluster)]
            star_offset_frac = _optional_float_value(SET.get('trna_modification_plots_star_offset')) if use_custom else None
            if star_offset_frac is None or not np.isfinite(star_offset_frac):
                star_offset_frac = 0.070
            star_offset_frac = float(np.clip(star_offset_frac, 0.010, 0.250))
            row_gap = 0.050 * span
            first_y = yhi - star_offset_frac * span
            for comp_rank, cname in enumerate(comp_clusters):
                star_y_by_cluster[cname] = first_y - comp_rank * row_gap
        max_annotation_y = yhi
        for i_f, f in enumerate(features):
            ref_vals = plotted_values.get((str(reference_cluster), str(f)), np.asarray([], dtype=float))
            feature_vals = []
            for c in clusters:
                vals = plotted_values.get((str(c), str(f)), np.asarray([], dtype=float))
                if vals is not None and len(vals):
                    feature_vals.extend(np.asarray(vals, dtype=float)[np.isfinite(vals)].tolist())
            local_top = float(np.nanmax(feature_vals)) if feature_vals else yhi
            comp_rank = 0
            for j, c in enumerate(clusters):
                if str(c) == str(reference_cluster):
                    continue
                vals = plotted_values.get((str(c), str(f)), np.asarray([], dtype=float))
                stat, pval, test_name = _compute_two_group_test(ref_vals, vals, method=stats_test)
                stars = _p_to_stars(pval)
                comp_rank += 1
                if fixed_plot6_star_rows:
                    ypos = float(star_y_by_cluster.get(str(c), yhi - 0.035 * span))
                    star_va = 'top'
                    star_clip = True
                else:
                    ypos = local_top + (0.06 + 0.055 * (comp_rank - 1)) * span
                    star_va = 'bottom'
                    star_clip = False
                xpos = float(xbase[i_f] + offsets[j])
                ax.text(xpos, ypos, stars, ha='center', va=star_va, fontsize=caption_size, fontweight='bold', clip_on=star_clip, zorder=8)
                max_annotation_y = max(max_annotation_y, ypos)
                annotation_top_by_feature[str(f)] = max(annotation_top_by_feature.get(str(f), -np.inf), ypos)
                stats_rows.append({
                    'plot': str(title_txt),
                    'comparison_type': 'cluster_vs_reference',
                    'feature_group': str((feature_metadata.get(f, {}) or {}).get('group', '')),
                    'feature': str(f),
                    'feature_label': str((feature_display.get(f, f) if isinstance(feature_display, dict) else display[i_f] if i_f < len(display) else f)),
                    'reference_cluster': str(reference_cluster),
                    'comparison_cluster': str(c),
                    'feature_A': str(f),
                    'feature_B': '',
                    'n_reference': int(len(ref_vals) if ref_vals is not None else 0),
                    'n_comparison': int(len(vals) if vals is not None else 0),
                    'mean_reference': float(np.nanmean(ref_vals)) if ref_vals is not None and len(ref_vals) else np.nan,
                    'mean_comparison': float(np.nanmean(vals)) if vals is not None and len(vals) else np.nan,
                    'test': test_name,
                    'statistic': stat,
                    'p_value': pval,
                    'significance': stars,
                    'values_are_log2_transformed': bool(use_log2),
                    'outliers_excluded': bool(exclude_outliers),
                    'outlier_sd_multiplier': float(outlier_sd),
                })
        if (not fixed_plot6_star_rows) and max_annotation_y > yhi:
            ax.set_ylim(top=max_annotation_y + 0.08 * span)
    elif reference_cluster and reference_cluster not in clusters and str(stats_test).strip().lower() not in {'', 'none', 'no', 'false', 'off'}:
        print(f'[WARN] Reference cluster "{reference_cluster}" is not among the plotted clusters; statistical annotations were skipped.')

    if pair_stats_enabled:
        # Compare x-axis features within each amino-acid/feature group.
        # Earlier builds annotated only groups with exactly two displayed features.
        # This version annotates all pairwise comparisons for groups with >=2
        # features, which is required for cases such as Arg isoacceptors in Plot 5.
        group_to_indices = {}
        for i_f, grp in enumerate(group_values):
            grp_s = str(grp or '').strip()
            if not grp_s:
                continue
            group_to_indices.setdefault(grp_s, []).append(i_f)

        if plot_key == 'shift':
            # Plot 5 (tRNA/isoacceptor shift): keep the annotation readable by
            # comparing only simple two-isoacceptor groups. Multi-isoacceptor
            # families such as Arg are intentionally skipped rather than drawing
            # all pairwise brackets above the graph.
            pair_groups = {grp: idxs for grp, idxs in group_to_indices.items() if len(idxs) == 2}
            skipped_groups_few = [grp for grp, idxs in group_to_indices.items() if len(idxs) < 2]
            skipped_groups_many = [grp for grp, idxs in group_to_indices.items() if len(idxs) > 2]
            if skipped_groups_few:
                print('[INFO] Within-group tRNA-pair statistics skipped groups with fewer than two displayed features: ' + ', '.join(map(str, skipped_groups_few)))
            if skipped_groups_many:
                print('[INFO] Within-group tRNA-pair statistics skipped multi-isoacceptor groups (>2 displayed features; no all-pairwise brackets drawn): ' + ', '.join(map(str, skipped_groups_many)))
        else:
            # Plot 4 (wobble/codon-pair) and other feature-pair plots can still
            # annotate all pairwise comparisons when more than two features are shown.
            pair_groups = {grp: idxs for grp, idxs in group_to_indices.items() if len(idxs) >= 2}
            skipped_groups = [grp for grp, idxs in group_to_indices.items() if len(idxs) < 2]
            if skipped_groups:
                print('[INFO] Within-group pair statistics skipped groups with fewer than two displayed features: ' + ', '.join(map(str, skipped_groups)))

        if pair_groups:
            ylo, yhi = ax.get_ylim()
            span = yhi - ylo if np.isfinite(yhi - ylo) and (yhi - ylo) > 0 else 1.0
            max_annotation_y = yhi
            bar_lw = max(0.9, caption_size * 0.095)
            cap_h = 0.018 * span
            text_gap = 0.012 * span

            # The GUI "Bracket gap" is the base distance between
            # successive bracket rows.  The first clearance, between the
            # unbracketed cluster-vs-reference stars and the first bracketed
            # within-feature comparison, is intentionally 50% larger so the
            # two statistical layers remain visually distinct.
            base_bracket_gap = max(
                pair_stats_gap_frac * span,
                1.35 * caption_size * span / max(fig_h * 72.0, 1.0),
            )
            layer_gap = 1.5 * base_bracket_gap
            pair_stack_gap = base_bracket_gap
            cluster_stack_gap = base_bracket_gap

            for grp, idxs in pair_groups.items():
                pair_index_list = list(itertools.combinations(list(idxs), 2))
                for pair_rank, (i_a, i_b) in enumerate(pair_index_list):
                    f_a, f_b = features[i_a], features[i_b]
                    for j, c in enumerate(clusters):
                        vals_a = plotted_values.get((str(c), str(f_a)), np.asarray([], dtype=float))
                        vals_b = plotted_values.get((str(c), str(f_b)), np.asarray([], dtype=float))
                        stat, pval, test_name = _compute_two_group_test(vals_a, vals_b, method=pair_stats_test)
                        stars = _p_to_stars(pval)

                        finite_a = np.asarray(vals_a, dtype=float)
                        finite_a = finite_a[np.isfinite(finite_a)]
                        finite_b = np.asarray(vals_b, dtype=float)
                        finite_b = finite_b[np.isfinite(finite_b)]
                        local_vals = []
                        if finite_a.size:
                            local_vals.extend(finite_a.tolist())
                        if finite_b.size:
                            local_vals.extend(finite_b.tolist())
                        local_top = float(np.nanmax(local_vals)) if local_vals else yhi

                        # Place the bracket well above any existing cluster-vs-reference
                        # stars for either feature, then stack additional pairwise
                        # and cluster-specific brackets above that.
                        prev_top = max(
                            annotation_top_by_feature.get(str(f_a), -np.inf),
                            annotation_top_by_feature.get(str(f_b), -np.inf),
                        )
                        if np.isfinite(prev_top):
                            local_top = max(
                                local_top + 0.14 * span,
                                float(prev_top) + layer_gap,
                                yhi + 0.055 * span,
                            )
                        else:
                            local_top = max(local_top + 0.24 * span, yhi + 0.055 * span)

                        ypos = local_top + pair_rank * pair_stack_gap + j * cluster_stack_gap

                        x1 = float(xbase[i_a] + offsets[j])
                        x2 = float(xbase[i_b] + offsets[j])
                        if x2 < x1:
                            x1, x2 = x2, x1
                        ax.plot([x1, x1, x2, x2], [ypos - cap_h, ypos, ypos, ypos - cap_h],
                                color='black', lw=bar_lw, clip_on=False, zorder=9)
                        ax.text(0.5 * (x1 + x2), ypos + text_gap, stars, ha='center', va='bottom',
                                fontsize=caption_size, fontweight='bold', clip_on=False, zorder=10)
                        max_annotation_y = max(max_annotation_y, ypos + text_gap)
                        annotation_top_by_feature[str(f_a)] = max(annotation_top_by_feature.get(str(f_a), -np.inf), ypos + text_gap)
                        annotation_top_by_feature[str(f_b)] = max(annotation_top_by_feature.get(str(f_b), -np.inf), ypos + text_gap)

                        label_a = str((feature_display.get(f_a, f_a) if isinstance(feature_display, dict) else display[i_a] if i_a < len(display) else f_a))
                        label_b = str((feature_display.get(f_b, f_b) if isinstance(feature_display, dict) else display[i_b] if i_b < len(display) else f_b))
                        stats_rows.append({
                            'plot': str(title_txt),
                            'comparison_type': 'within_group_feature_pair',
                            'feature_group': str(grp),
                            'feature': '',
                            'feature_label': f'{label_a} vs {label_b}',
                            'reference_cluster': '',
                            'comparison_cluster': str(c),
                            'feature_A': str(f_a),
                            'feature_B': str(f_b),
                            'n_reference': int(len(vals_a) if vals_a is not None else 0),
                            'n_comparison': int(len(vals_b) if vals_b is not None else 0),
                            'mean_reference': float(np.nanmean(vals_a)) if vals_a is not None and len(vals_a) else np.nan,
                            'mean_comparison': float(np.nanmean(vals_b)) if vals_b is not None and len(vals_b) else np.nan,
                            'test': test_name,
                            'statistic': stat,
                            'p_value': pval,
                            'significance': stars,
                            'values_are_log2_transformed': bool(use_log2),
                            'outliers_excluded': bool(exclude_outliers),
                            'outlier_sd_multiplier': float(outlier_sd),
                        })
            if max_annotation_y > yhi:
                ax.set_ylim(top=max_annotation_y + 0.18 * span)
        else:
            print('[INFO] Within-group pair statistics were requested, but no displayed feature group had at least two features.')

    if stats_rows:
        stats_path = os.path.splitext(out_path)[0] + ' - statistics.xlsx'
        exported = _export_decoding_stats_excel(stats_rows, stats_path)
        if exported:
            print(f'[INFO] Decoding statistics exported: {exported}')

    ax.set_xlim(float(np.nanmin(xbase)) - 0.9, float(np.nanmax(xbase)) + 0.9)
    ax.grid(axis='y', alpha=0.22, lw=0.6)
    handles = [Patch(facecolor=colors[i], edgecolor='black', alpha=0.68,
                     label=f'{str(c)} (n={per_cluster_n.get(c, 0)})') for i, c in enumerate(clusters)]
    if ax2 is not None and secondary_label:
        # Keep the secondary-axis ylabel as requested by the plot, but make the
        # legend entry one-line so it does not consume extra vertical space above
        # the graph.
        legend_secondary_label = re.sub(r'\s+', ' ', str(secondary_label).replace('\n', ' ')).strip()
        handles.append(Patch(facecolor='grey', edgecolor='grey', alpha=max(0.18, float(SET.get('trna_secondary_axis_alpha', 0.22) or 0.22)), label=legend_secondary_label))
    # Legend is placed above the plot. Up to four entries are shown per
    # row by default, so cluster legends wrap cleanly when many clusters are
    # displayed. Plot 6 can override the number of columns in Figure details.
    custom_legend_ncol = None
    if is_modification_plot:
        custom_legend_ncol = _optional_int_value(SET.get('trna_modification_plots_legend_ncol'))
    if custom_legend_ncol is not None and custom_legend_ncol > 0:
        legend_ncol = max(1, min(int(custom_legend_ncol), max(1, len(handles))))
    else:
        legend_ncol = max(1, min(4, len(handles)))
    legend_rows = int(math.ceil(len(handles) / float(legend_ncol))) if handles else 1
    if is_modification_plot:
        # Keep the legend immediately above the axes. The figure-level title is
        # placed above the legend, preventing title/legend overlap.
        legend_y = 1.020
    else:
        legend_y = 1.045 + 0.050 * max(0, legend_rows - 1)
    ax.legend(
        handles=handles,
        loc='lower center',
        bbox_to_anchor=(0.5, legend_y),
        ncol=legend_ncol,
        borderaxespad=0.0,
        fontsize=caption_size,
        frameon=True,
        columnspacing=1.1,
        handletextpad=0.45,
    )
    try:
        # Make room for one or more legend rows above the axis while keeping
        # the secondary y-axis visible on the right.
        if feature_metadata and is_modification_plot:
            # Lower than ordinary x-labels, but far less compressed than v43.
            bottom_margin = 0.34
        elif feature_metadata:
            bottom_margin = 0.32
        else:
            bottom_margin = 0.20
        if is_modification_plot:
            top_margin = max(0.62, 0.765 - 0.045 * max(0, legend_rows - 1))
            title_y = min(0.985, top_margin + 0.185 + 0.030 * max(0, legend_rows - 1))
            if title_txt:
                fig.suptitle(str(title_txt), fontsize=caption_size, y=title_y)
            aa_subtitle = _format_plot6_amino_acid_subtitle(SET)
            # Put the amino-acid reference text in the whitespace between title
            # and legend, intentionally closer to the title than to the legend.
            subtitle_y = max(top_margin + 0.092, title_y - 0.060)
            fig.text(0.5, subtitle_y, aa_subtitle, ha='center', va='center',
                     fontsize=max(6, caption_size - 1))
        else:
            top_margin = max(0.66, 0.89 - 0.050 * max(0, legend_rows - 1))
        fig.subplots_adjust(bottom=bottom_margin, right=(0.91 if is_modification_plot else 0.94), top=top_margin)
    except Exception:
        pass
    try:
        os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
        fig.savefig(out_path, dpi=dpi, bbox_inches='tight')
    finally:
        if show_fig:
            _show_figure_nonblocking(fig)
        else:
            plt.close(fig)
    return out_path

def _compute_wc_wobble_pair_counts(count_df: pd.DataFrame, trna_rules: dict, allowed_aas=None):
    codon_to_decoders = dict(trna_rules.get('codon_to_decoders', {}) or {})
    aa_order, aa_to_labels, aa_to_codons = _build_trna_family_maps(trna_rules)
    allowed = {str(a) for a in list(allowed_aas or [])}

    ordered_pair_labels = []
    for aa in aa_order:
        if allowed and aa not in allowed:
            continue
        labels = aa_to_labels.get(aa, [])
        if len(labels) != 1:
            continue
        label = labels[0]
        seen_modes = set()
        for codon_internal in aa_to_codons.get(aa, []):
            for _lab, _weight in list(codon_to_decoders.get(codon_internal, []) or []):
                if str(_lab) != str(label):
                    continue
                component_modes = []
                anti_group = str(label).split('_', 1)[1] if '_' in str(label) else ''
                for anti in _split_anticodon_group(anti_group):
                    mode = _pairing_mode_for_anticodon_codon(anti, codon_internal)
                    if mode in {'WC', 'Wobble'}:
                        component_modes.append(mode)
                for mode in ['WC', 'Wobble']:
                    if mode in component_modes and mode not in seen_modes:
                        ordered_pair_labels.append(f'{label}__{mode}')
                        seen_modes.add(mode)

    out = pd.DataFrame(0.0, index=count_df.index, columns=ordered_pair_labels, dtype=float)
    if out.empty:
        return out

    for codon in count_df.columns:
        codon_s = str(codon)
        aa = codon_s.split('_', 1)[0] if '_' in codon_s else ''
        if allowed and aa not in allowed:
            continue
        decoders = list(codon_to_decoders.get(codon_s, []) or [])
        if not decoders:
            continue
        codon_counts = pd.to_numeric(count_df[codon_s], errors='coerce').fillna(0.0)
        for label, weight in decoders:
            label = str(label)
            anti_group = label.split('_', 1)[1] if '_' in label else ''
            modes = []
            for anti in _split_anticodon_group(anti_group):
                mode = _pairing_mode_for_anticodon_codon(anti, codon_s)
                if mode in {'WC', 'Wobble'}:
                    modes.append(mode)
            modes = sorted(set(modes), key=lambda x: 0 if x == 'WC' else 1)
            if not modes:
                continue
            frac = float(weight) / float(len(modes))
            for mode in modes:
                col = f'{label}__{mode}'
                if col in out.columns:
                    out[col] = pd.to_numeric(out[col], errors='coerce').fillna(0.0) + codon_counts * frac
    return out


_SUPP1_TRNA_ORDER = [
    "Gly_UCC/CCC",
    "Leu_UAA/CAA",
    "Ile_CAU",
    "Arg_CCG",
    "Arg_UCU/CCU",
    "Ser_GCU",
    "Ser_UGA/GGA/CGA",
    "Ile_GAU",
    "Arg_ACG",
    "Leu_UAG/GAG/CAG",
    "Gly_GCC",
]

def _anticodon_group_key(label):
    lab = str(label or '')
    if '_' not in lab:
        return lab.strip().lower()
    aa, anti = lab.split('_', 1)
    parts = [p.strip().upper().replace('T', 'U') for p in re.split(r'[/,+;]+', anti) if p.strip()]
    return aa.strip().lower() + '_' + '/'.join(sorted(parts))

def _pick_ordered_labels_with_group_aliases(preferred_order, available_columns):
    available = [str(c) for c in list(available_columns) if str(c).strip()] if available_columns is not None else []
    by_key = {_anticodon_group_key(c): c for c in available}
    selected = []
    seen = set()
    for preferred in list(preferred_order or []):
        lab = preferred if preferred in available else by_key.get(_anticodon_group_key(preferred), '')
        if lab and lab not in seen:
            selected.append(lab)
            seen.add(lab)
    return selected

_SUPP1_TRNA_LABEL_MAP = {lab: lab.replace("_", "-", 1) for lab in _SUPP1_TRNA_ORDER}

_SUPP2_SINGLE_TRNA_CODON_ORDER = [
    "Tyr_TAT", "Phe_TTT", "Asn_AAT", "His_CAT", "Asp_GAT", "Cys_TGT", "Glu_GAG", "Lys_AAA",
    "Lys_AAG", "Glu_GAA", "Cys_TGC", "Asp_GAC", "His_CAC", "Asn_AAC", "Phe_TTC", "Tyr_TAC",
]

def render_trna_gene_ordered_heatmaps(SET, count_df: pd.DataFrame, ordered_genes, output_dir: str, cluster_df: pd.DataFrame | None = None, selected_clusters=None):
    # This function can now be called directly from the GUI to replot decoding
    # figures without rerunning clustering.  Fill in the plotting defaults that
    # are normally provided by the full pipeline settings dictionary.
    _defaults = dict(
        export_trna_usage_enable=True,
        heatmap_fig_size=(18, 4),
        figure_dpi=300,
        heatmap_colormap_name='parula',
        heatmap_caxis_limits=(-2.5, 2.5),
        font_name='Arial',
        font_size_xticks=8,
        font_size_yticks=7,
        font_size_titles=10,
        colorbar_title_size=11,
        colorbar_title_string='Enrichment',
        xtick_every_genes=500,
        apply_binning=False,
        bin_size_genes=50,
        apply_smoothing=False,
        smooth_window_genes=6,
        feature_dist_metric='spearman',
        feature_linkage='single',
        cluster_use_optimal_leaf_ordering=True,
        cluster_optimal_leaf_max_size=250,
        cluster_fast_order_threshold=800,
        usage_basis='RCU',
        codon_set='59',
        trna_gene_wobble_plot_kind='heatmap',
        trna_gene_trna_plot_kind='heatmap',
        trna_mrna_stability_enable=False,
        trna_mrna_stability_plot_kind='line',
        trna_gene_wobble_smooth=True,
        trna_gene_wobble_smooth_method='running average',
        trna_gene_wobble_smooth_window=5,
        trna_gene_trna_smooth=True,
        trna_gene_trna_smooth_method='running average',
        trna_gene_trna_smooth_window=5,
        trna_mrna_stability_smooth=True,
        trna_mrna_stability_smooth_method='running average',
        trna_mrna_stability_smooth_window=5,
        trna_wobble_boxplot_log2=True,
        trna_shift_boxplot_log2=True,
        trna_modifications_boxplot_log2=True,
        trna_modifications_boxplot_caption_size=17,
        trna_modification_plots_caption_size=17,
        trna_modification_plots_box_width=0.18,
        trna_modifications_selected_features=None,
        trna_modifications_include_aas=None,  # None -> default Plot 6 amino-acid set: Ala/Arg/Asn/Asp/Cys/Gln/Glu/Gly/His/Ile/Leu/Lys/Phe/Pro/Ser/Thr/Tyr/Val
    )
    _tmp = dict(_defaults)
    _tmp.update(dict(SET or {}))
    SET = _tmp

    out = {}
    if not bool(SET.get('export_trna_usage_enable', False)):
        return out
    trna_path = str(SET.get('trna_decoding_table_path', '') or '').strip()
    if not trna_path:
        return out
    if count_df is None or getattr(count_df, 'empty', True):
        return out

    trna_rules = read_trna_decoding_table(trna_path, sheet_name=SET.get('trna_decoding_table_sheet', ''))
    ordered_gene_list = [str(g).strip() for g in (list(ordered_genes) if ordered_genes is not None else []) if str(g).strip() and str(g).strip().lower() != 'nan']
    if not ordered_gene_list:
        return out

    # Make direct decoding robust to old/new raw workbook formats and U/T codon labels.
    count_df = count_df.copy()
    count_df.index = pd.Index([str(i).strip() for i in count_df.index], dtype=object)
    norm_cols = []
    for c in count_df.columns:
        cs = str(c).strip()
        if '_' in cs:
            aa, cod = cs.split('_', 1)
            cod = cod.strip().upper().replace('U', 'T')
            cs = f'{aa.strip()}_{cod}'
        norm_cols.append(cs)
    count_df.columns = norm_cols
    count_df = count_df.apply(pd.to_numeric, errors='coerce').fillna(0.0)
    if count_df.columns.duplicated().any():
        count_df = count_df.T.groupby(level=0).sum().T
    count_df = count_df[~count_df.index.duplicated(keep='first')]

    overlap_order = [g for g in ordered_gene_list if g in count_df.index]
    print(f'[INFO] Decoding direct plot gene-ID overlap: {len(overlap_order)} / {len(ordered_gene_list)} reordered genes found in codon-count table.')
    if len(overlap_order) == 0:
        preview_counts = ', '.join(list(map(str, count_df.index[:5])))
        preview_order = ', '.join(list(map(str, ordered_gene_list[:5])))
        raise ValueError(
            'No reordered genes from the clustering workbook match the raw codon-count table.\n'
            'This usually means the raw workbook was read with the wrong header row, or the project mixes different identifier modes.\n'
            f'First reordered genes: {preview_order}\n'
            f'First count-table genes: {preview_counts}'
        )

    os.makedirs(output_dir or '.', exist_ok=True)
    custom_cmaps = load_custom_colormaps(SET)

    trna_counts_df = compute_trna_counts_from_codon_counts_df(count_df, trna_rules)
    atu_df = compute_acu_from_counts_df(trna_counts_df)
    rtu_df = compute_rtu_from_counts_df(trna_counts_df)
    ztu_df = compute_ztu(rtu_df, baseline_genes=overlap_order)

    # Shared codon-usage tables for gene-ordered and cluster-level decoding plots.
    rcu_df = compute_rcu_from_counts_df(count_df)
    zcu_df = compute_rcu_devz(rcu_df, baseline_genes=overlap_order)

    aa_full_name = {
        'Ala': 'Alanine', 'Arg': 'Arginine', 'Asn': 'Asparagine', 'Asp': 'Aspartate',
        'Cys': 'Cysteine', 'Gln': 'Glutamine', 'Glu': 'Glutamate', 'Gly': 'Glycine',
        'His': 'Histidine', 'Ile': 'Isoleucine', 'Leu': 'Leucine', 'Lys': 'Lysine',
        'Met': 'Methionine', 'Phe': 'Phenylalanine', 'Pro': 'Proline', 'Ser': 'Serine',
        'Thr': 'Threonine', 'Trp': 'Tryptophan', 'Tyr': 'Tyrosine', 'Val': 'Valine',
    }
    requested_pairs = [
        ('Tyr', 'TAC', 'TAT'),
        ('Phe', 'TTC', 'TTT'),
        ('Asn', 'AAC', 'AAT'),
        ('His', 'CAC', 'CAT'),
        ('Asp', 'GAC', 'GAT'),
        ('Cys', 'TGC', 'TGT'),
        ('Glu', 'GAA', 'GAG'),
        ('Lys', 'AAA', 'AAG'),
    ]

    # ---------------- Plot 1: wobble decoding along reordered genome ----------------
    if bool(SET.get('trna_single_box_codon_heatmap_enable', True)):
        plot_kind = str(SET.get('trna_gene_wobble_plot_kind', 'heatmap') or 'heatmap').strip().lower()
        if plot_kind in {'line', 'surface', 'area'}:
            profile = _compute_wobble_percent_profile(count_df, ordered_gene_list, requested_pairs)
            out_path = os.path.join(output_dir, f'Plot 1 - wobble decoding along reordered genome {plot_kind}.png')
            saved = _plot_ordered_profile_line_surface(
                SET,
                values=profile,
                ordered_genes=ordered_gene_list,
                out_path=out_path,
                title_txt='Plot 1: wobble decoding along reordered genome',
                y_label='Mean wobble codon usage (%)',
                plot_kind=plot_kind,
                smooth=bool(SET.get('trna_gene_wobble_smooth', True)),
                smooth_window=int(SET.get('trna_gene_wobble_smooth_window', 5) or 5),
                smooth_method=str(SET.get('trna_gene_wobble_smooth_method', 'running average') or 'running average'),
                show_fig=bool(SET.get('trna_single_box_codon_heatmap_show_fig', True)),
                caption_size=int(SET.get('trna_gene_wobble_caption_size', 13) or 13),
            )
            if saved:
                out['wobble_decoding_profile'] = saved
        else:
            wc_rows = []
            wobble_rows = []
            row_labels = {}
            for aa, wc_codon, wobble_codon in requested_pairs:
                wc_lab = f'{aa}_{wc_codon}'
                wob_lab = f'{aa}_{wobble_codon}'
                if wc_lab in zcu_df.columns:
                    wc_rows.append(wc_lab)
                    row_labels[wc_lab] = f"{aa_full_name.get(aa, aa)} {wc_codon}"
                if wob_lab in zcu_df.columns:
                    wobble_rows.append(wob_lab)
                    row_labels[wob_lab] = f"{aa_full_name.get(aa, aa)} {wobble_codon}"
            # Sort within WC and wobble blocks by decreasing mean ZCU over the reordered genome.
            means = zcu_df.reindex(index=ordered_gene_list).mean(axis=0, skipna=True)
            wc_rows = sorted(wc_rows, key=lambda lab: float(means.get(lab, -np.inf)), reverse=True)
            wobble_rows = sorted(wobble_rows, key=lambda lab: float(means.get(lab, -np.inf)), reverse=True)
            selected_single_codons = wc_rows + wobble_rows
            if selected_single_codons:
                zcu_df_gene = zcu_df.loc[:, selected_single_codons].copy()
                out_path = os.path.join(output_dir, 'Plot 1 - wobble decoding along reordered genome heatmap.png')
                saved = _plot_feature_by_gene_heatmap(
                    SET,
                    values_df=zcu_df_gene,
                    ordered_genes=ordered_gene_list,
                    title_txt='Plot 1 — wobble decoding along reordered genome',
                    y_label='Codons',
                    out_path=out_path,
                    custom_cmaps=custom_cmaps,
                    show_fig=bool(SET.get('trna_single_box_codon_heatmap_show_fig', True)),
                    colorbar_label='ZCU',
                    caxis_limits=SET.get('heatmap_caxis_limits', None),
                    feature_label_map=row_labels,
                    reorder_feature_rows=False,
                )
                if saved:
                    out['wobble_decoding_heatmap'] = saved

    # ---------------- Plot 2: tRNA usage shift along reordered genome ----------------
    heatmap_metric = str(SET.get('trna_gene_heatmap_metric', 'ZTU') or 'ZTU').upper()
    metric_df_map = {'ATU': atu_df, 'RTU': rtu_df, 'ZTU': ztu_df}
    selected_trna_labels = _pick_ordered_labels_with_group_aliases(_SUPP1_TRNA_ORDER, trna_counts_df.columns)
    if bool(SET.get('trna_gene_heatmap_enable', True)):
        plot_kind = str(SET.get('trna_gene_trna_plot_kind', 'heatmap') or 'heatmap').strip().lower()
        if plot_kind in {'line', 'surface', 'area'}:
            profile = _compute_rare_trna_fraction_profiles(count_df, trna_rules, ordered_gene_list)
            out_path = os.path.join(output_dir, f'Plot 2 - rare tRNA fraction along reordered genome {plot_kind}.png')
            saved = _plot_ordered_profile_line_surface(
                SET,
                values=profile,
                ordered_genes=ordered_gene_list,
                out_path=out_path,
                title_txt='Plot 2: rare tRNA fraction along reordered genome',
                y_label='Fraction of rare tRNA used for decoding (%)',
                plot_kind=plot_kind,
                smooth=bool(SET.get('trna_gene_trna_smooth', True)),
                smooth_window=int(SET.get('trna_gene_trna_smooth_window', 5) or 5),
                smooth_method=str(SET.get('trna_gene_trna_smooth_method', 'running average') or 'running average'),
                show_fig=bool(SET.get('trna_gene_heatmap_show_fig', True)),
                caption_size=int(SET.get('trna_gene_trna_caption_size', 13) or 13),
            )
            if saved:
                out['trna_rare_fraction_profile'] = saved
        elif selected_trna_labels:
            metric_df = metric_df_map.get(heatmap_metric, ztu_df)
            metric_df = metric_df.loc[:, selected_trna_labels].copy()
            if metric_df.shape[1] > 0:
                out_path = os.path.join(output_dir, f'Plot 2 - tRNA usage shift along reordered genome heatmap ({heatmap_metric}).png')
                saved = _plot_feature_by_gene_heatmap(
                    SET,
                    values_df=metric_df,
                    ordered_genes=ordered_gene_list,
                    title_txt=f'Plot 2 — tRNA usage shift along reordered genome ({heatmap_metric})',
                    y_label='tRNAs',
                    out_path=out_path,
                    custom_cmaps=custom_cmaps,
                    show_fig=bool(SET.get('trna_gene_heatmap_show_fig', True)),
                    colorbar_label=heatmap_metric,
                    caxis_limits=_metric_caxis_from_name(SET, heatmap_metric),
                    feature_label_map={lab: lab.replace("_", "-", 1) for lab in selected_trna_labels},
                    reorder_feature_rows=False,
                )
                if saved:
                    out['trna_usage_heatmap'] = saved

    # ---------------- Plot 3: mRNA stability along reordered genome ----------------
    if bool(SET.get('trna_mrna_stability_enable', False)):
        plot_kind = str(SET.get('trna_mrna_stability_plot_kind', 'line') or 'line').strip().lower()
        if plot_kind not in {'line', 'surface', 'area'}:
            plot_kind = 'line'
        stability = _read_rna_stability_from_decoding_workbook(trna_path, ordered_gene_list, sheet_name='RNA stability')
        out_path = os.path.join(output_dir, f'Plot 3 - mRNA stability along reordered genome {plot_kind}.png')
        saved = _plot_ordered_profile_line_surface(
            SET,
            values=stability,
            ordered_genes=ordered_gene_list,
            out_path=out_path,
            title_txt='Plot 3: mRNA stability along reordered genome',
            y_label='mRNA half-life',
            plot_kind=plot_kind,
            smooth=bool(SET.get('trna_mrna_stability_smooth', True)),
            smooth_window=int(SET.get('trna_mrna_stability_smooth_window', 5) or 5),
            smooth_method=str(SET.get('trna_mrna_stability_smooth_method', 'running average') or 'running average'),
            show_fig=True,
            caption_size=int(SET.get('trna_mrna_stability_caption_size', 13) or 13),
        )
        if saved:
            out['mrna_stability_profile'] = saved

    # ---------------- Cluster-level plots 4, 5 and 6 ----------------
    if cluster_df is None or getattr(cluster_df, 'empty', True):
        return out

    # Clean cluster table values and report whether selected cluster genes match the count table.
    cluster_df = cluster_df.copy().fillna('')
    for _c in cluster_df.columns:
        cluster_df[_c] = cluster_df[_c].astype(str).str.strip()

    if selected_clusters is None:
        available_clusters = [str(c) for c in cluster_df.columns]
    else:
        available_clusters = [str(c) for c in list(selected_clusters or []) if str(c) in list(cluster_df.columns)]
    if not available_clusters:
        print('[INFO] No decoding clusters selected; skipping cluster-level decoding plots 4–6.')
        return out
    cluster_gene_hits = 0
    cluster_gene_total = 0
    for _c in available_clusters:
        _genes = [str(v).strip() for v in cluster_df[_c].tolist() if str(v).strip() and str(v).strip().lower() != 'nan']
        cluster_gene_total += len(_genes)
        cluster_gene_hits += sum(1 for _g in _genes if _g in count_df.index)
    print(f'[INFO] Decoding selected-cluster overlap: {cluster_gene_hits} / {cluster_gene_total} cluster gene entries found in codon-count table.')
    if cluster_gene_total > 0 and cluster_gene_hits == 0:
        print('[WARN] None of the selected cluster gene IDs match the codon-count table; cluster boxplots/heatmaps will be empty.')

    aa_order, aa_to_labels, aa_to_codons = _build_trna_family_maps(trna_rules)

    aa_full_name = {
        'Ala': 'Alanine', 'Arg': 'Arginine', 'Asn': 'Asparagine', 'Asp': 'Aspartate',
        'Cys': 'Cysteine', 'Gln': 'Glutamine', 'Glu': 'Glutamate', 'Gly': 'Glycine',
        'His': 'Histidine', 'Ile': 'Isoleucine', 'Leu': 'Leucine', 'Lys': 'Lysine',
        'Met': 'Methionine', 'Phe': 'Phenylalanine', 'Pro': 'Proline', 'Ser': 'Serine',
        'Thr': 'Threonine', 'Trp': 'Tryptophan', 'Tyr': 'Tyrosine', 'Val': 'Valine',
    }

    # Plot 4: cluster-level z-score synonymous codon usage. Include both the
    # Watson-Crick codon and the wobble codon for each single-tRNA/two-codon box,
    # with the WC codon plotted first and the wobble codon second.
    if bool(SET.get('trna_wobble_heatmap_enable', True)):
        requested_pairs = [
            ('Tyr', 'TAC', 'TAT'),
            ('Phe', 'TTC', 'TTT'),
            ('Asn', 'AAC', 'AAT'),
            ('His', 'CAC', 'CAT'),
            ('Asp', 'GAC', 'GAT'),
            ('Cys', 'TGC', 'TGT'),
            ('Glu', 'GAA', 'GAG'),
            ('Lys', 'AAA', 'AAG'),
        ]

        codon_rows = []
        codon_display = []
        codon_meta = {}
        missing_requested = []
        for aa, wc_codon, wobble_codon in requested_pairs:
            for codon, mode in [(wc_codon, 'WC'), (wobble_codon, 'Wobble')]:
                internal_label = f'{aa}_{codon}'
                if internal_label in zcu_df.columns:
                    codon_rows.append(internal_label)
                    codon_display.append(codon)
                    codon_meta[internal_label] = {
                        'top': codon,
                        'top_color': 'red' if mode == 'Wobble' else 'black',
                        'bottom': aa_full_name.get(aa, aa),
                        'bottom_color': 'black',
                        'group': aa,
                        'group_label': aa_full_name.get(aa, aa),
                    }
                else:
                    missing_requested.append(internal_label)

        if missing_requested:
            print('[WARN] Plot 4: requested WC/wobble codon(s) absent from ZCU table and skipped: '
                  + ', '.join(missing_requested))

        if codon_rows:
            chosen = list(available_clusters)
            plot_kind = str(SET.get('trna_wobble_plot_kind', 'boxplot') or 'boxplot').strip().lower()
            zcu_manual_only = zcu_df.loc[:, codon_rows].copy()
            if plot_kind in {'boxplot', 'violin'}:
                secondary = _genome_codon_frequency_percent(count_df, codon_rows)
                out_path = os.path.join(output_dir, f'Plot 4 - enrichment in wobble decoding within clusters {plot_kind}s.png')
                saved = _plot_cluster_feature_boxplots(
                    SET,
                    values_df=zcu_manual_only,
                    cluster_df=cluster_df,
                    selected_clusters=chosen,
                    feature_order=codon_rows,
                    feature_display={k: v for k, v in zip(codon_rows, codon_display)},
                    out_path=out_path,
                    title_txt='z-score synonymous codon usage',
                    y_label='z-score synonymous codon usage',
                    secondary_values=secondary,
                    secondary_std=None,
                    secondary_label='Codon freq.\nin genome (%)',
                    show_fig=bool(SET.get('trna_wobble_heatmap_show_fig', True)),
                    plot_key='wobble',
                    feature_metadata=codon_meta,
                )
                if saved:
                    out['wobble_cluster_enrichment'] = saved
            else:
                mat = _cluster_feature_mean_matrix(zcu_manual_only, cluster_df, chosen)
                if not mat.empty:
                    mat = mat.reindex(codon_rows)
                    row_labels = [f"{aa_full_name.get(r.split('_', 1)[0], r.split('_', 1)[0])} {r.split('_', 1)[1]}" for r in codon_rows]
                    mat.index = row_labels
                    out_path = os.path.join(output_dir, 'Plot 4 - enrichment in wobble decoding within clusters heatmap.png')
                    wobble_red_mask = [bool((codon_meta.get(r, {}) or {}).get('top_color') == 'red') for r in codon_rows]
                    saved = _plot_cluster_shift_heatmap(
                        SET,
                        mat,
                        out_path=out_path,
                        title_txt='z-score synonymous codon usage',
                        show_fig=bool(SET.get('trna_wobble_heatmap_show_fig', True)),
                        colorbar_label='Codon usage enrichment\nvs genome',
                        group_bounds=None,
                        row_display_labels=row_labels,
                        x_label='Gene clusters',
                        y_label='',
                        log2_colorbar=bool(SET.get('trna_wobble_heatmap_log2_colorbar', True)),
                        red_row_mask=wobble_red_mask,
                    )
                    if saved:
                        out['wobble_cluster_enrichment'] = saved

    # Plot 5: tRNA usage shift within clusters (ZTU values), only multi-tRNA families.
    if bool(SET.get('trna_shift_heatmap_enable', True)):
        trna_rows = []
        trna_display = []
        trna_meta = {}
        for aa in aa_order:
            labels = [lab for lab in aa_to_labels.get(aa, []) if lab in ztu_df.columns]
            if len(labels) <= 1:
                continue
            for lab_idx, lab in enumerate(labels):
                anti = lab.split('_', 1)[1] if '_' in lab else lab
                trna_rows.append(lab)
                trna_display.append(anti)
                # Leu and Ser have multiple adjacent isoacceptor labels that can
                # overlap; stagger these labels onto two vertical rows.
                top_y = -0.070
                if aa in {'Leu', 'Ser'}:
                    # Stagger adjacent Leu/Ser anticodon labels with a slightly
                    # larger separation so the two text rows remain readable.
                    top_y = -0.070 if (lab_idx % 2 == 0) else -0.150
                trna_meta[lab] = {
                    'top': anti,
                    'top_color': 'black',
                    'top_y': top_y,
                    'bottom': aa_full_name.get(aa, aa),
                    'bottom_color': 'black',
                    'group': aa,
                    'group_label': aa_full_name.get(aa, aa),
                }
        if trna_rows:
            chosen = list(available_clusters)
            plot_kind = str(SET.get('trna_shift_plot_kind', 'boxplot') or 'boxplot').strip().lower()
            ztu_selected = ztu_df.loc[:, trna_rows].copy()
            if plot_kind in {'boxplot', 'violin'}:
                abundance = trna_rules.get('abundance_series')
                abundance_std = trna_rules.get('abundance_std_series')
                if abundance is not None:
                    abundance = pd.to_numeric(pd.Series(abundance), errors='coerce')
                if abundance_std is not None:
                    abundance_std = pd.to_numeric(pd.Series(abundance_std), errors='coerce')
                out_path = os.path.join(output_dir, f'Plot 5 - tRNA usage shift within clusters {plot_kind}s.png')
                saved = _plot_cluster_feature_boxplots(
                    SET,
                    values_df=ztu_selected,
                    cluster_df=cluster_df,
                    selected_clusters=chosen,
                    feature_order=trna_rows,
                    feature_display={k: v for k, v in zip(trna_rows, trna_display)},
                    out_path=out_path,
                    title_txt='tRNA usage shift vs genome',
                    y_label='tRNA usage shift vs genome (ZTU)',
                    secondary_values=abundance,
                    secondary_std=abundance_std,
                    secondary_label='tRNA abundance\n(molecules/cell)',
                    show_fig=bool(SET.get('trna_shift_heatmap_show_fig', True)),
                    plot_key='shift',
                    feature_metadata=trna_meta,
                )
                if saved:
                    out['trna_usage_cluster_shift'] = saved
            else:
                mat = _cluster_feature_mean_matrix(ztu_selected, cluster_df, chosen)
                if not mat.empty:
                    row_labels = [lab.replace('_', ' ', 1) for lab in trna_rows]
                    mat = mat.reindex(trna_rows)
                    mat.index = row_labels
                    out_path = os.path.join(output_dir, 'Plot 5 - tRNA usage shift within clusters heatmap.png')
                    saved = _plot_cluster_shift_heatmap(
                        SET,
                        mat,
                        out_path=out_path,
                        title_txt='tRNA usage shift vs genome',
                        show_fig=bool(SET.get('trna_shift_heatmap_show_fig', True)),
                        colorbar_label='tRNA enrichment\nvs genome',
                        group_bounds=None,
                        row_display_labels=row_labels,
                        x_label='Gene clusters',
                        y_label='',
                        log2_colorbar=bool(SET.get('trna_shift_heatmap_log2_colorbar', True)),
                    )
                    if saved:
                        out['trna_usage_cluster_shift'] = saved

    # Plot 6: shift in codons associated with tRNA modifications or tRNA
    # modification enzymes.  This is the only branch using the conservative /
    # permissive decoder-assignment models; all tRNA usage / ZTU plots keep the
    # original pooled-decoder behavior.
    if bool(SET.get('trna_modification_heatmap_enable', True)):
        feature_mode = str(SET.get('trna_modifications_feature_mode', 'modifications') or 'modifications').strip().lower()
        models_raw = SET.get('trna_modifications_assignment_models', 'conservative,permissive')
        if isinstance(models_raw, str):
            assignment_models = [m.strip().lower() for m in re.split(r'[,;/]+', models_raw) if m.strip()]
        else:
            assignment_models = [str(m).strip().lower() for m in list(models_raw or []) if str(m).strip()]
        disabled_models = [m for m in assignment_models if m == 'legacy']
        disabled_fractional = [m for m in assignment_models if m == 'fractional']
        assignment_models = [m for m in assignment_models if m in {'conservative', 'permissive'}]
        if disabled_models:
            print("[WARN] Plot 6 legacy row-level model was requested but is disabled: Plot 6 now uses only 'Decoding table (full)'.")
        if disabled_fractional:
            print("[INFO] Plot 6 fractional/estimated-fraction model is disabled and will not be exported.")
        if not assignment_models:
            assignment_models = ['conservative', 'permissive']

        chosen = list(available_clusters)
        selected_mod_features = _normalize_optional_string_list(SET.get('trna_modifications_selected_features', None))
        selected_mod_aas = _normalize_optional_string_list(SET.get('trna_modifications_include_aas', None))
        if selected_mod_features is not None:
            print(f"[INFO] Plot 6: restricting plotted features to {len(selected_mod_features)} selected modification/enzyme item(s).")
        if selected_mod_aas is not None:
            print(f"[INFO] Plot 6: restricting calculation to amino acid(s): {', '.join(selected_mod_aas)}.")
        else:
            print('[INFO] Plot 6: using default amino-acid set: ' + ', '.join(PLOT6_DEFAULT_AAS) + '.')
        plot_kind = str(SET.get('trna_modifications_plot_kind', 'boxplot') or 'boxplot').strip().lower()
        mode_label = 'tRNA modification enzymes' if feature_mode in {'enzyme', 'enzymes', 'trme', 'trmes'} else 'tRNA modifications'
        model_display = {
            'conservative': 'conservative decoder model',
            'permissive': 'permissive decoder model',
            
        }
        model_suffix = {
            'conservative': 'conservative',
            'permissive': 'permissive',
            
        }

        # Compute all requested assignment models first, then force them onto
        # one shared feature axis. Without this, the conservative plot can lose
        # modifications that exist only in the permissive model,
        # causing Plot 6 variants to have different x labels and ordering.
        model_results = {}
        label_maps = {}
        shared_feature_order = []
        shared_label_map = {}

        # Use the broadest models first to define a stable, comparable x-axis.
        # Preserve user-requested models, then add any remaining results.
        feature_order_scan = [m for m in ['permissive', 'conservative'] if m in assignment_models]
        feature_order_scan += [m for m in assignment_models if m not in feature_order_scan]

        for assignment_model in assignment_models:
            mod_z_df, mod_label_map = _compute_modification_usage_z(
                count_df,
                trna_rules,
                feature_mode=feature_mode,
                assignment_model=assignment_model,
                selected_features=selected_mod_features,
                include_aas=selected_mod_aas,
            )
            model_results[assignment_model] = mod_z_df
            label_maps[assignment_model] = mod_label_map or {}

        for assignment_model in feature_order_scan:
            mod_z_df = model_results.get(assignment_model)
            mod_label_map = label_maps.get(assignment_model, {})
            if mod_z_df is None or mod_z_df.empty:
                continue
            for feat in list(mod_z_df.columns):
                feat = str(feat)
                if feat not in shared_feature_order:
                    shared_feature_order.append(feat)
                if feat not in shared_label_map:
                    shared_label_map[feat] = str(mod_label_map.get(feat, feat))

        plot6_is_modifications = feature_mode not in {'enzyme', 'enzymes', 'trme', 'trmes', 'trme enzymes'}
        if plot6_is_modifications and shared_feature_order:
            shared_feature_order = _order_modification_features_by_position(shared_feature_order, shared_label_map)
            shared_label_map = {feat: str(shared_label_map.get(feat, feat)) for feat in shared_feature_order}
            plot6_feature_metadata = _build_modification_position_metadata(shared_feature_order, shared_label_map)
        else:
            plot6_feature_metadata = None

        plot6_output_dir = _trna_modifications_output_dir(output_dir) if plot6_is_modifications else output_dir
        if plot6_is_modifications:
            os.makedirs(plot6_output_dir, exist_ok=True)
            exported_mod_table = _export_trna_modification_enrichment_table(
                output_dir=plot6_output_dir,
                count_df=count_df,
                cluster_df=cluster_df,
                selected_clusters=chosen,
                trna_rules=trna_rules,
                assignment_models=assignment_models,
                include_aas=selected_mod_aas,
                use_log2=_setting_bool(SET, 'trna_modifications_boxplot_log2', True),
            )
            if exported_mod_table:
                out['trna_modification_enrichment_table'] = exported_mod_table

        if not shared_feature_order:
            print(f"[WARN] Plot 6: no {mode_label} features could be computed in any assignment model.")
        else:
            for assignment_model in assignment_models:
                mod_z_df = model_results.get(assignment_model)
                if mod_z_df is None or mod_z_df.empty:
                    # Keep the plot axis comparable even when a strict model has
                    # no feature-specific codons. All-zero values indicate no
                    # genome-normalized enrichment signal under that model.
                    mod_z_df = pd.DataFrame(0.0, index=count_df.index, columns=shared_feature_order)
                    print(f"[INFO] Plot 6 ({assignment_model}): no direct {mode_label} features; plotting shared labels with zero enrichment.")
                else:
                    mod_z_df = mod_z_df.copy()
                    for feat in shared_feature_order:
                        if feat not in mod_z_df.columns:
                            mod_z_df[feat] = 0.0
                    mod_z_df = mod_z_df.reindex(columns=shared_feature_order)

                mod_rows = list(shared_feature_order)
                mod_label_map = dict(shared_label_map)
                disp_model = model_display.get(assignment_model, assignment_model)
                suffix = model_suffix.get(assignment_model, assignment_model)
                if plot_kind in {'boxplot', 'violin'}:
                    out_path = os.path.join(plot6_output_dir, f'Plot 6 - {mode_label} shift within clusters - {suffix} {plot_kind}s.png')
                    saved = _plot_cluster_feature_boxplots(
                        SET,
                        values_df=mod_z_df,
                        cluster_df=cluster_df,
                        selected_clusters=chosen,
                        feature_order=mod_rows,
                        feature_display=mod_label_map,
                        out_path=out_path,
                        title_txt=(suffix if plot6_is_modifications else f'{mode_label} enrichment vs genome ({disp_model})'),
                        y_label='tRNA enrichment' if plot6_is_modifications else f'{mode_label} enrichment\nvs genome (z-score; {suffix})',
                        secondary_values=None,
                        secondary_std=None,
                        secondary_label='',
                        show_fig=bool(SET.get('trna_modification_heatmap_show_fig', True)),
                        plot_key='modifications',
                        feature_metadata=plot6_feature_metadata,
                    )
                    if saved:
                        out[f'trna_modification_cluster_shift_{suffix}'] = saved
                        # Keep the historical key pointing to the first generated
                        # modification plot for backward compatibility with older GUI logs.
                        out.setdefault('trna_modification_cluster_shift', saved)
                else:
                    mat = _cluster_feature_mean_matrix(mod_z_df, cluster_df, chosen)
                    if not mat.empty:
                        display = [mod_label_map.get(m, m) for m in mod_rows]
                        mat = mat.reindex(mod_rows)
                        mat.index = display
                        out_path = os.path.join(plot6_output_dir, f'Plot 6 - {mode_label} shift within clusters - {suffix} heatmap.png')
                        saved = _plot_cluster_shift_heatmap(
                            SET,
                            mat,
                            out_path=out_path,
                            title_txt=(f'tRNA modifications enrichment vs genome ({disp_model})' if 'modification' in mode_label else f'{mode_label} enrichment vs genome ({disp_model})'),
                            show_fig=bool(SET.get('trna_modification_heatmap_show_fig', True)),
                            colorbar_label=f'{mode_label}\nenrichment vs genome\n({suffix})',
                            group_bounds=None,
                            row_display_labels=display,
                            x_label='Gene clusters',
                            y_label=mode_label,
                            log2_colorbar=False,
                        )
                        if saved:
                            out[f'trna_modification_cluster_shift_{suffix}'] = saved
                            out.setdefault('trna_modification_cluster_shift', saved)

    return out


def plot_scatter(SET, Y, usage_basis, title_cset, custom_cmaps):
    if Y is None:
        return None, None

    fig, ax = plt.subplots(figsize=SET['scatter_fig_size'], dpi=SET['figure_dpi'])

    point_size = SET.get('scatter_point_size', 2.0)
    alpha = SET.get('scatter_point_alpha', 0.6)
    edge_width = SET.get('scatter_edge_width', 0.0)
    color_mode = SET.get('scatter_color_mode', 'plain').lower().strip()

    colors = None
    if color_mode in ("density", "enrichment"):
        density_cmap_name = SET.get('density_cmap_name', 'plasma')
        density_cmap = get_cmap_any(density_cmap_name, custom_maps=custom_cmaps, fallback='plasma')

        colors = _compute_density_colors(
            Y,
            nbins=SET.get('density_nbins', 150),
            sigma=SET.get('density_sigma', 4.0),
            use_log=SET.get('density_use_log', True),
            min_rel=SET.get('density_min_rel', 0.0),
            max_rel=SET.get('density_max_rel', 1.0),
            cmap=density_cmap,
            metric=color_mode,
            enrichment_eps=SET.get("density_enrichment_eps", 1e-12),
            enrichment_use_log=SET.get("density_enrichment_use_log", True),
        )

    ax.scatter(
        Y[:, 0], Y[:, 1],
        s=point_size,
        c=colors if colors is not None else None,
        alpha=alpha,
        linewidths=edge_width
    )

    method = normalize_dimred_method(SET['dimred_method'])
    if method == 'tsne':
        xlab, ylab = 'tSNE1', 'tSNE2'
    elif method == 'umap':
        xlab, ylab = 'UMAP1', 'UMAP2'
    elif method == 'pca':
        xlab, ylab = 'PC1', 'PC2'
    else:
        xlab, ylab = 'Dim 1', 'Dim 2'

    ax.set_xlabel(xlab, fontname=SET['font_name'], fontsize=SET['font_size_axes'])
    ax.set_ylabel(ylab, fontname=SET['font_name'], fontsize=SET['font_size_axes'])
    ax.tick_params(labelsize=max(1, SET['font_size_axes'] - 1))

    tag_usage = usage_basis.upper()
    tag_dim = SET['dimred_method'].upper()
    tag_clu = SET['cluster_method'].upper()
    title_txt = SET['scatter_title_template'].format(
        USAGE=tag_usage, CSET=title_cset, DIMRED=tag_dim, CLUSTER=tag_clu
    )

    if method == 'tsne':
        param = (f"tSNE params: perplexity={SET['tsne_perplexity']}, "
                 f"distance={SET['tsne_distance']}, "
                 f"exaggeration={SET['tsne_exaggeration']}, "
                 f"learnrate={SET['tsne_learnrate']}")
    elif method == 'umap':
        param = (f"UMAP params: n_neighbors={SET['umap_neighbors']}, "
                 f"min_dist={SET['umap_min_dist']:.3f}, "
                 f"metric={SET['umap_metric']}, "
                 f"randomize={int(SET['umap_randomize'])}")
    elif method == 'pca':
        param = (f"PCA params: nPCs={Y.shape[1]}, "
                 f"center={int(SET['pca_center'])}, "
                 f"scale={int(SET['pca_scale'])}")
    else:
        param = ""

    extra = ""
    cm = SET.get("scatter_color_mode", "plain").lower().strip()
    if cm == "enrichment":
        extra = " | color=enrichment(vs mean density)"
    elif cm == "density":
        extra = " | color=density"

    try:
        if SET.get('scatter_xmin', None) not in (None, "") or SET.get('scatter_xmax', None) not in (None, ""):
            ax.set_xlim(left=(float(SET.get('scatter_xmin')) if SET.get('scatter_xmin', None) not in (None, "") else None),
                        right=(float(SET.get('scatter_xmax')) if SET.get('scatter_xmax', None) not in (None, "") else None))
    except Exception:
        pass
    try:
        if SET.get('scatter_ymin', None) not in (None, "") or SET.get('scatter_ymax', None) not in (None, ""):
            ax.set_ylim(bottom=(float(SET.get('scatter_ymin')) if SET.get('scatter_ymin', None) not in (None, "") else None),
                        top=(float(SET.get('scatter_ymax')) if SET.get('scatter_ymax', None) not in (None, "") else None))
    except Exception:
        pass
    ax.set_title(f"{param}\n{title_txt}{extra}" if param else f"{title_txt}{extra}",
                 fontsize=SET['font_size_titles'], fontname=SET['font_name'])
    return fig, ax


def make_run_output_dir(base_folder, SET):
    method = SET['dimred_method'].lower()
    if method == 'tsne':
        per = SET.get('tsne_perplexity', 30)
        ex = SET.get('tsne_exaggeration', 12)
        lr = SET.get('tsne_learnrate', 200)
        name = f"tsne per{per} ex{ex} lr{lr}"
    elif method == 'umap':
        nn = SET.get('umap_neighbors', 15)
        md = SET.get('umap_min_dist', 0.1)
        md_str = f"{md:.3g}".replace('.', 'p')
        name = f"umap nn{nn} md{md_str}"
    elif method == 'pca':
        pcs = SET.get('pca_npcs', 2)
        name = f"pca npc{pcs}"
    else:
        name = "nodimred"

    requested_root = str(SET.get('default_root', '') or '').strip()
    if requested_root:
        try:
            out_dir = os.path.abspath(requested_root)
            os.makedirs(out_dir, exist_ok=True)
            return out_dir, "direct export folder"
        except Exception:
            pass

    out_dir = os.path.join(base_folder, name)
    os.makedirs(out_dir, exist_ok=True)
    return out_dir, name

# =========================================================
# BIG FUNCTION MOVED OUT: run_codon_clustering
# =========================================================

def run_codon_clustering(SET):
    SET['usage_basis'] = normalize_usage_basis(SET.get('usage_basis', 'RCU'))
    SET['codon_set'] = normalize_codon_set(SET.get('codon_set', '59_noSTOP_MW'))
    SET['dimred_method'] = normalize_dimred_method(SET.get('dimred_method', 'umap'))
    SET['cluster_method'] = normalize_cluster_method(SET.get('cluster_method', 'kmeans'))

    if SET['usage_basis'] == 'ACU' and codon_set_kind(SET['codon_set']) != '61_nostop':
        print("Warning: for ACU, recommended codon_set is '61_noSTOP'.")
    auto_cset = auto_codon_set_for_usage_basis(SET['usage_basis'])
    if normalize_codon_set(SET.get('codon_set', auto_cset)) != auto_cset:
        print(f"[INFO] Codon set auto-adjusted to {auto_cset} for usage basis {SET['usage_basis']}.")
    SET['codon_set'] = auto_cset

    print("---------------------- Open input data ----------------------")

    from .fasta_metrics import choose_fasta, compute_codon_usage_tables_from_cds_fasta

    fasta_path_for_input = str(SET.get("fasta_path", "") or SET.get("input_fasta", "")).strip()
    if not fasta_path_for_input:
        fasta_path_for_input = choose_fasta(SET.get('default_root', ''))
    if not os.path.isfile(fasta_path_for_input):
        raise FileNotFoundError(f"FASTA not found: {fasta_path_for_input}")

    base = os.path.splitext(os.path.basename(fasta_path_for_input))[0]
    prefix_hint = infer_prefix_from_codon_basename(base)
    strain_prefix = prefix_hint
    if strain_prefix and not strain_prefix.endswith("_"):
        strain_prefix += "_"

    include_stops_default = (codon_set_kind(SET.get("codon_set", "61_noSTOP")) == "64_withstop")

    fasta_codon_range = str(SET.get("fasta_codon_range", "all") or "all").strip() or "all"
    freq_df, geneids_df, abs_df, gc_pct = compute_codon_usage_tables_from_cds_fasta(
        fasta_path=fasta_path_for_input,
        row_id_mode=str(SET.get("fasta_row_id_mode", "primary")),
        trim_to_multiple_of_3=bool(SET.get("fasta_trim_to_multiple_of_3", True)),
        include_stops=bool(SET.get("fasta_include_stops", include_stops_default)),
        keep_first_duplicate=True,
        organism_mode='prokaryote',
        codon_range=fasta_codon_range,
    )
    if geneids_df is not None and (not geneids_df.empty) and "CodonRange" in geneids_df.columns:
        _range_series = geneids_df["CodonRange"].dropna().astype(str)
        codon_range_label = str(_range_series.iloc[0]) if not _range_series.empty else fasta_codon_range
    else:
        codon_range_label = fasta_codon_range
    n_empty_region = int((pd.to_numeric(geneids_df.get("SelectedCodonsCounted", pd.Series(dtype=float)), errors="coerce").fillna(0) == 0).sum()) if geneids_df is not None else 0
    print(f"[INFO] FASTA codon range used for codon counts: {codon_range_label}")
    if n_empty_region:
        print(f"[WARN] {n_empty_region} CDS record(s) had no codons in the requested range and were kept with zero counts.")

    count_df = abs_df.copy()
    if str(SET.get("fasta_codon_table_mode", "abs")).strip().lower() not in ("abs", "counts", "count"):
        print("[INFO] Harmonized codon-feature construction uses raw codon counts internally; FASTA codon table mode is ignored for clustering features.")

    gene_symbol_map = {}
    gene_desc_map = {}
    if geneids_df is not None and (not geneids_df.empty):
        for _, r in geneids_df.iterrows():
            lt = str(r.get("LocusTag", "")).strip()
            if not lt or lt.lower() == "nan":
                continue
            gs = str(r.get("GeneSymbol", "") or "").strip()
            ds = str(r.get("ProteinDescription", "") or "").strip()
            if gs and gs.lower() != "nan" and gs != "NA":
                gene_symbol_map[lt] = gs
            if ds and ds.lower() != "nan" and ds != "NA":
                gene_desc_map[lt] = ds

    codon_file = fasta_path_for_input
    gene_file = None
    geneids_xlsx = None
    gc_txt_path = None

    base_folder = os.path.dirname(fasta_path_for_input)
    output_dir, output_subfolder = make_run_output_dir(base_folder, SET)
    print(f"[INFO] Output subfolder: {output_subfolder}")

    if bool(SET.get("fasta_write_intermediate_excels", False)):
        try:
            geneids_xlsx = os.path.join(output_dir, f"{base}_geneIDs.xlsx")
            export_geneids_df = _prepare_export_geneids_df(geneids_df) if geneids_df is not None else None
            if export_geneids_df is not None:
                export_geneids_df.to_excel(geneids_xlsx, index=False)
            gc_txt_path = os.path.join(output_dir, "GC_content.txt")
            with open(gc_txt_path, "w", encoding="utf-8") as fh:
                fh.write(f"Codon range used for codon counts: {codon_range_label}\n")
                fh.write(f"GC content of concatenated selected CDS regions: {gc_pct:.2f}%\n")
            print(f"[INFO] Wrote intermediate tables in: {output_dir}")
        except Exception as e:
            print(f"[WARN] Could not write intermediate FASTA-derived tables: {e}")
            geneids_xlsx = None
            gc_txt_path = None

    SET['usage_basis'] = normalize_usage_basis(SET.get('usage_basis', 'ACU'))
    auto_cset = auto_codon_set_for_usage_basis(SET['usage_basis'])
    if normalize_codon_set(SET.get('codon_set', auto_cset)) != auto_cset:
        print(f"[INFO] Codon set auto-adjusted to {auto_cset} for usage basis {SET['usage_basis']}.")
    SET['codon_set'] = auto_cset

    C_abs_df, _, AA_df, C_rel_df, aa_names, _, RowNames, var_names = compute_AA_ACU_RCU(count_df)

    Usage, feature_labels, _ = subset_usage(
        SET, SET['usage_basis'], C_abs_df.to_numpy(dtype=float), C_rel_df, AA_df, aa_names, var_names
    )
    values = normalize_values(SET, Usage)
    nGenes_full = values.shape[0]

    Y = run_dimred(SET, values, feature_labels)
    if Y is None:
        print("[WARN] dimred_method=='none'; falling back to PCA(2) for coordinates.")
        SET_fallback = dict(SET)
        SET_fallback['dimred_method'] = 'pca'
        Y = run_dimred(SET_fallback, values, feature_labels)

    gene_order, labels = cluster_genes(SET, Y)
    ordered_genes = RowNames[gene_order]

    features_reorder, feat_order = reorder_features(SET, values, feature_labels)
    V = values[gene_order, :].T[feat_order, :]
    V_smooth, bin_size = smooth_and_bin(SET, V)

    custom_cmaps = load_custom_colormaps(SET)

    tag_usage = SET['usage_basis'].upper()
    tag_dim = SET['dimred_method'].upper()
    tag_clu = SET['cluster_method'].upper()
    title_cset = f"{AA_df.shape[1]} AAs" if tag_usage == 'AA' else f"{Usage.shape[1]} codons"

    heatmap_title = SET['heatmap_title_template'].format(
        USAGE=tag_usage, CSET=title_cset, DIMRED=tag_dim, CLUSTER=tag_clu
    )

    fig_h = None
    if bool(SET.get('plot_codon_gene_heatmap_enable', True)):
        fig_h, _ = plot_heatmap(SET, V_smooth, features_reorder, nGenes_full, bin_size, heatmap_title, custom_cmaps)
    else:
        print("[INFO] Main codon vs genes heatmap disabled by user settings.")

    generate_scatter_figure = not bool(SET.get('ordered_show_mode', True))
    if generate_scatter_figure:
        fig_s, _ = plot_scatter(SET, Y, SET['usage_basis'], title_cset, custom_cmaps)
    else:
        fig_s = None

    tag_cset = title_cset.replace(' ', '')
    tag_pipe = f"{tag_dim}_{tag_clu}"
    if SET['dimred_method'].lower() == 'pca' and Y is not None:
        tag_pipe = f"{tag_pipe}_PC{Y.shape[1]}"
    out_base = os.path.join(output_dir, f"{base}_{tag_usage}_{tag_cset}_{tag_pipe}")

    fmt = str(SET.get('figure_output_format', '') or '').strip().lstrip('.').lower()
    if fmt == 'jpg':
        fmt = 'jpeg'
    if fmt:
        if fig_h is not None:
            fig_h.savefig(out_base + f"_heatmap.{fmt}", dpi=SET['figure_dpi'], bbox_inches='tight')
        if fig_s is not None:
            fig_s.savefig(out_base + f"_scatter.{fmt}", dpi=SET['figure_dpi'], bbox_inches='tight')
    else:
        if SET.get('save_png', True):
            if fig_h is not None:
                fig_h.savefig(out_base + "_heatmap.png", dpi=SET['figure_dpi'], bbox_inches='tight')
            if fig_s is not None:
                fig_s.savefig(out_base + "_scatter.png", dpi=SET['figure_dpi'], bbox_inches='tight')
        if SET.get('save_pdf', False):
            if fig_h is not None:
                fig_h.savefig(out_base + "_heatmap.pdf", bbox_inches='tight')
            if fig_s is not None:
                fig_s.savefig(out_base + "_scatter.pdf", bbox_inches='tight')
        if SET.get('save_jpeg', False):
            if fig_h is not None:
                fig_h.savefig(out_base + "_heatmap.jpeg", dpi=SET['figure_dpi'], bbox_inches='tight')
            if fig_s is not None:
                fig_s.savefig(out_base + "_scatter.jpeg", dpi=SET['figure_dpi'], bbox_inches='tight')

    if fig_h is not None:
        if bool(SET.get('plot_codon_gene_heatmap_show_fig', SET.get('show_main_pipeline_figures', True))):
            _show_figure_nonblocking(fig_h)
        else:
            plt.close(fig_h)
    if fig_s is not None and bool(SET.get('show_main_pipeline_scatter', False)):
        _show_figure_nonblocking(fig_s)
    elif fig_s is not None:
        plt.close(fig_s)

    completed_outputs = []
    if fig_h is not None:
        completed_outputs.append('heatmap')
    if fig_s is not None:
        completed_outputs.append('scatter')
    suffix = ' + '.join(completed_outputs) if completed_outputs else 'tables only'
    print(f"\n[INFO] Codon-usage clustering completed ({suffix}).")

    return dict(
        input_mode='fasta',
        fasta_path=fasta_path_for_input,
        fasta_codon_range=fasta_codon_range,
        fasta_codon_range_label=codon_range_label,
        fasta_codon_range_empty_genes=n_empty_region,
        codon_file=codon_file,
        gene_file=gene_file,
        geneids_df=geneids_df,
        geneids_xlsx=geneids_xlsx,
        strain_prefix=strain_prefix,
        prefix_hint=prefix_hint,
        gene_symbol_map=gene_symbol_map,
        gene_desc_map=gene_desc_map,
        RowNames=RowNames,
        Usage=Usage,
        values=values,
        Y=Y,
        ordered_genes=ordered_genes,
        features_reorder=features_reorder,
        heatmap_matrix=V_smooth,
        heatmap_bin_size=bin_size,
        heatmap_n_genes_full=nGenes_full,
        heatmap_title=heatmap_title,
        heatmap_x_label='Genes ordered by codon usage similarity',
        heatmap_y_label='Amino acids' if tag_usage == 'AA' else 'Codons',
        AA_df=AA_df,
        C_abs_df=C_abs_df,
        C_rel_df=C_rel_df,
        codon_count_df=count_df,
        output_dir=output_dir,
        output_subfolder=output_subfolder,
        codon_xlsx='',
        gc_txt_path=gc_txt_path,
    )
