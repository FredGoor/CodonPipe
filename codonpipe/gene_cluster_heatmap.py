# codonpipe/gene_cluster_heatmap.py
"""Cluster landscape along the reordered genome axis (PLOT 2).

Replaces the earlier per-row min-max density heatmap with KS stars. Two things
changed, both because the old display could not be compared across rows:

1.  Colour is now log2(observed / expected) genes per bin on ONE shared scale,
    so a given colour means the same enrichment in every row. The old display
    rescaled each row to its own maximum, which made a 1.9-fold enrichment in
    one row look identical to a 6-fold enrichment in another.
2.  The axis is divided into equal-count bins of genes (not equal spans of
    locus tags), so bin area is directly comparable along the axis.

The per-row statistic is the overlap coefficient against a reference cluster:
sum of min(p, q) over bins, where p and q are the two binned densities each
normalised to sum 1. It is threshold-free, runs from 0 (disjoint) to 1
(identical), and is tested against random gene sets of the same size.

Pairwise statistics between clusters are reported separately (short menu of
named comparisons and the full n x n matrix), not drawn on the figure.

Expected workbook structure (default: sheet 'Locus Tags'):
- Column 1: reordered genome locus tags (defines the axis)
- Columns 2..N: each column is a gene cluster listing locus tags (header row)
"""

from __future__ import annotations

import os
import json
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy.ndimage import gaussian_filter1d

from .colormaps import get_cmap_any


DEFAULT_BINS = 120
DEFAULT_SMOOTHING = 1.5
DEFAULT_COLOR_CAP = 2.6
DEFAULT_PERMUTATIONS = 2000
DEFAULT_FRAME_COLORS = ["#762a83", "#e08214", "#b2182b", "#2166ac", "#1b7837"]
DEFAULT_CORE_SHARE = 0.5


# ----------------------------
# Colormap helpers
# ----------------------------

def _load_custom_colormaps(xlsx_path: str) -> Dict[str, LinearSegmentedColormap]:
    """Load custom colormaps from an Excel file.

    Excel format:
      - one sheet per colormap (sheet name = cmap name)
      - no header; columns: [index, R, G, B] with RGB in 0-255 or 0-1
    """
    maps: Dict[str, LinearSegmentedColormap] = {}
    if not xlsx_path:
        return maps
    if not os.path.isfile(xlsx_path):
        print(f"[WARN] Custom cmap Excel not found: {xlsx_path}")
        return maps

    try:
        xls = pd.ExcelFile(xlsx_path)
        for sheet in xls.sheet_names:
            df = xls.parse(sheet, header=None)
            if df.shape[1] < 4:
                print(f"[WARN] Sheet '{sheet}' has <4 columns; expected index + R,G,B. Skipped.")
                continue

            rgb_raw = df.iloc[:, 1:4].apply(pd.to_numeric, errors="coerce")
            rgb = rgb_raw.values
            rgb = rgb[~np.isnan(rgb).any(axis=1)]
            if rgb.shape[0] < 2:
                print(f"[WARN] Sheet '{sheet}' has <2 valid color rows. Skipped.")
                continue

            if np.nanmax(rgb) > 1.0:
                rgb = np.clip(rgb / 255.0, 0.0, 1.0)
            else:
                rgb = np.clip(rgb, 0.0, 1.0)

            maps[sheet] = LinearSegmentedColormap.from_list(sheet, [tuple(c) for c in rgb], N=256)

        if maps:
            print(f"[INFO] Loaded custom colormaps: {', '.join(maps.keys())}")
        else:
            print("[INFO] No valid custom colormaps found.")
    except Exception as e:
        print(f"[WARN] Failed to parse custom colormaps: {e}")
    return maps


def _choose_cmap(name: str, custom_maps: Dict[str, LinearSegmentedColormap]):
    """Pick a colormap by name, preferring custom maps when available.

    The landscape plot is a signed obs/exp display, so the sensible default is a
    diverging map centred on 1.0 (no enrichment). 'RdBu_r' is used when the
    requested name is not found.
    """
    return get_cmap_any(name, custom_maps=custom_maps, fallback="RdBu_r")


def get_significance_stars(p: float) -> str:
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def _clean_tag(x) -> str:
    if x is None:
        return ""
    s = str(x)
    if s.lower() == "nan":
        return ""
    return s.strip().strip('"').strip("'").strip()


# ----------------------------
# Workbook input
# ----------------------------

SPACER_PREFIX = "__spacer__"

#: Tokens a caller may put in ``include_clusters`` to request a blank row that
#: does not exist as a column in the workbook. The GUI's cluster-rearrangement
#: dialog uses these so a user can insert separators anywhere in the row order.
BLANK_REQUEST_TOKENS = {"__blank__", "<blank>", "(blank)", "blank", "spacer", "space", "separator", ""}


def is_spacer_column(col_name, values) -> bool:
    """True when a cluster column is meant as a visual separator, not a cluster.

    A column counts as a spacer when its header is blank (pandas names those
    'Unnamed: N'), or is literally 'blank' or 'spacer', and it lists no locus
    tags. The requirement that it be empty is deliberate: a real cluster whose
    header someone forgot to fill in should raise attention as an oddly named
    row, not silently disappear into a gap.
    """
    name = str(col_name).strip().lower()
    if values:
        return False
    return (name == "" or name.startswith("unnamed:")
            or name in ("blank", "spacer", "space", "separator"))


def load_locus_tags_sheet(
    workbook_path: str,
    sheet_name: str = "Locus Tags",
    include_clusters: Optional[List[str]] = None,
) -> Tuple[List[str], Dict[str, List[str]]]:
    """Load genome-ordered locus tags and cluster columns from the workbook sheet."""
    if not os.path.isfile(workbook_path):
        raise FileNotFoundError(f"Workbook not found: {workbook_path}")

    df = pd.read_excel(workbook_path, sheet_name=sheet_name, dtype=str)
    if df.shape[1] < 2:
        raise ValueError(
            f"Sheet '{sheet_name}' must have at least 2 columns (genome axis + >=1 cluster). Found {df.shape[1]}."
        )

    ref_col = df.columns[0]
    reference_order = [_clean_tag(v) for v in df[ref_col].tolist()]
    reference_order = [v for v in reference_order if v != ""]
    if not reference_order:
        raise ValueError(f"No locus tags found in first column of sheet '{sheet_name}'.")

    clusters: Dict[str, List[str]] = {}
    data_cols = list(df.columns[1:])
    if include_clusters is not None:
        lower_to_actual = {str(col).strip().lower(): str(col) for col in data_cols}
        ordered_cols = []
        seen_cols = set()
        for req in include_clusters:
            key = str(req).strip().lower()
            if key in BLANK_REQUEST_TOKENS:
                # Requested separator: keep a placeholder in the row order even
                # though no such column exists in the workbook.
                ordered_cols.append(None)
                continue
            mapped = lower_to_actual.get(key)
            if mapped is not None and mapped not in seen_cols:
                ordered_cols.append(mapped)
                seen_cols.add(mapped)
        if not [c for c in ordered_cols if c is not None]:
            raise ValueError("No requested clusters were found in the locus-tags sheet.")
    else:
        ordered_cols = [str(col) for col in data_cols]

    spacer_n = 0
    for col_name in ordered_cols:
        if col_name is None:
            spacer_n += 1
            clusters[f"{SPACER_PREFIX}{spacer_n}"] = []
            continue
        vals = [_clean_tag(v) for v in df[col_name].tolist()]
        vals = [v for v in vals if v != ""]
        if is_spacer_column(col_name, vals):
            spacer_n += 1
            clusters[f"{SPACER_PREFIX}{spacer_n}"] = []
            continue
        clusters[col_name] = vals

    return reference_order, clusters


# ----------------------------
# Binning and statistics
# ----------------------------

def _bin_index(positions: np.ndarray, n_genes: int, nbins: int) -> np.ndarray:
    """Map gene rank positions to equal-count bin indices."""
    if n_genes <= 0:
        return np.zeros(0, dtype=int)
    return np.clip((np.asarray(positions, float) / float(n_genes) * nbins).astype(int), 0, nbins - 1)


def _density(bins: np.ndarray, nbins: int, smoothing: float) -> np.ndarray:
    counts = np.bincount(bins, minlength=nbins).astype(float)
    if smoothing and smoothing > 0:
        counts = gaussian_filter1d(counts, float(smoothing), mode="nearest")
    return counts


def _enrichment(bins: np.ndarray, nbins: int, smoothing: float, pseudo: float = 0.5) -> np.ndarray:
    """log2(observed / expected) genes per bin, with a pseudo-count."""
    obs = _density(bins, nbins, smoothing)
    exp = len(bins) / float(nbins)
    return np.log2((obs + pseudo) / (exp + pseudo))


def _normalised_density(bins: np.ndarray, nbins: int, smoothing: float) -> np.ndarray:
    d = _density(bins, nbins, smoothing)
    total = d.sum()
    return d / total if total > 0 else d


def overlap_coefficient(bins_a: np.ndarray, bins_b: np.ndarray, nbins: int, smoothing: float) -> float:
    """Sum of min(p, q) over bins; 0 = disjoint, 1 = identical."""
    p = _normalised_density(bins_a, nbins, smoothing)
    q = _normalised_density(bins_b, nbins, smoothing)
    return float(np.minimum(p, q).sum())


def overlap_coefficient_p(
    bins_a: np.ndarray,
    bins_b: np.ndarray,
    n_genes: int,
    nbins: int,
    smoothing: float,
    permutations: int = DEFAULT_PERMUTATIONS,
    rng: Optional[np.random.Generator] = None,
) -> Tuple[float, float]:
    """Observed overlap and its one-sided P against random gene sets of the same size.

    The null draws a random set of the same size as cluster A from the genome
    axis and recomputes the overlap with cluster B, so the test controls for
    cluster size (small clusters overlap less by construction).
    """
    rng = rng or np.random.default_rng(0)
    obs = overlap_coefficient(bins_a, bins_b, nbins, smoothing)
    if permutations <= 0 or len(bins_a) == 0 or len(bins_b) == 0:
        return obs, float("nan")
    q = _normalised_density(bins_b, nbins, smoothing)
    n_a = len(bins_a)
    hits = 0
    for _ in range(int(permutations)):
        draw = rng.choice(n_genes, size=n_a, replace=False)
        p = _normalised_density(_bin_index(draw, n_genes, nbins), nbins, smoothing)
        if float(np.minimum(p, q).sum()) >= obs:
            hits += 1
    return obs, (hits + 1) / (float(permutations) + 1)


def separation_test(
    bins_a: np.ndarray,
    bins_b: np.ndarray,
    nbins: int,
    smoothing: float,
    permutations: int = DEFAULT_PERMUTATIONS,
    rng: Optional[np.random.Generator] = None,
) -> Tuple[float, float, float]:
    """Are two clusters more separated along the axis than chance?

    Null: pool the two clusters and re-split at the same two sizes. This is
    unbiased with respect to cluster size, unlike drawing both sets from the
    whole genome. Returns (observed distance, z, P).
    """
    rng = rng or np.random.default_rng(0)
    pa = _normalised_density(bins_a, nbins, smoothing)
    pb = _normalised_density(bins_b, nbins, smoothing)
    obs = float(np.abs(pa - pb).sum())
    if permutations <= 0 or len(bins_a) == 0 or len(bins_b) == 0:
        return obs, float("nan"), float("nan")
    pool = np.concatenate([bins_a, bins_b])
    n_a = len(bins_a)
    null = np.empty(int(permutations), dtype=float)
    for i in range(int(permutations)):
        perm = rng.permutation(pool)
        qa = _normalised_density(perm[:n_a], nbins, smoothing)
        qb = _normalised_density(perm[n_a:], nbins, smoothing)
        null[i] = float(np.abs(qa - qb).sum())
    sd = float(null.std(ddof=1))
    z = (obs - float(null.mean())) / sd if sd > 0 else float("nan")
    p = (int((null >= obs).sum()) + 1) / (float(permutations) + 1)
    return obs, z, p


def core_band(bins: np.ndarray, share: float = DEFAULT_CORE_SHARE) -> Tuple[int, int]:
    """Narrowest contiguous bin stretch holding `share` of a cluster's genes."""
    if len(bins) == 0:
        return (0, 0)
    b = np.sort(bins)
    k = max(1, int(np.ceil(float(share) * len(b))))
    if k > len(b):
        k = len(b)
    best = None
    for i in range(0, len(b) - k + 1):
        width = b[i + k - 1] - b[i]
        if best is None or width < best[0]:
            best = (width, int(b[i]), int(b[i + k - 1]))
    return (best[1], best[2])


# ----------------------------
# Parameter parsing helpers
# ----------------------------

def _parse_pairs(spec) -> List[Tuple[str, str]]:
    """Parse 'A vs B; C vs D' (or a list of 2-tuples) into pairs."""
    if not spec:
        return []
    if isinstance(spec, (list, tuple)):
        out = []
        for item in spec:
            if isinstance(item, (list, tuple)) and len(item) == 2:
                out.append((str(item[0]).strip(), str(item[1]).strip()))
        return out
    pairs = []
    for chunk in str(spec).split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        for sep in (" vs ", " VS ", " v. ", "|"):
            if sep in chunk:
                a, b = chunk.split(sep, 1)
                pairs.append((a.strip(), b.strip()))
                break
    return pairs


def _parse_groups(spec) -> List[Tuple[str, List[str]]]:
    """Parse 'Label: A, B | Label2: C, D' (or a dict) into row groups."""
    if not spec:
        return []
    if isinstance(spec, dict):
        return [(str(k), [str(v) for v in vals]) for k, vals in spec.items()]
    text = str(spec).strip()
    if text.startswith("{"):
        try:
            return _parse_groups(json.loads(text))
        except Exception:
            pass
    groups = []
    for chunk in text.split("|"):
        chunk = chunk.strip()
        if not chunk or ":" not in chunk:
            continue
        label, members = chunk.split(":", 1)
        names = [m.strip() for m in members.split(",") if m.strip()]
        if names:
            groups.append((label.strip(), names))
    return groups


def _parse_frames(spec) -> List[Tuple[str, str]]:
    """Parse 'Label: Cluster | Label2: Cluster2' into frame definitions."""
    return [(lab, names[0]) for lab, names in _parse_groups(spec) if names]


# ----------------------------
# Statistics tables
# ----------------------------

def compute_pair_statistics(
    cluster_bins: Dict[str, np.ndarray],
    cluster_members: Dict[str, Sequence[str]],
    n_genes: int,
    nbins: int,
    smoothing: float,
    key_pairs: Optional[Sequence[Tuple[str, str]]] = None,
    full_matrix: bool = True,
    permutations: int = DEFAULT_PERMUTATIONS,
    seed: int = 0,
) -> Dict[str, pd.DataFrame]:
    """Pairwise colocalisation statistics: named short menu and full n x n matrices."""
    names = list(cluster_bins.keys())
    tables: Dict[str, pd.DataFrame] = {}
    rng = np.random.default_rng(seed)

    if key_pairs:
        rows = []
        for a, b in key_pairs:
            if a not in cluster_bins or b not in cluster_bins:
                print(f"[WARN] Key comparison skipped, cluster not plotted: {a} vs {b}")
                continue
            ov = overlap_coefficient(cluster_bins[a], cluster_bins[b], nbins, smoothing)
            _, z, p_sep = separation_test(cluster_bins[a], cluster_bins[b], nbins, smoothing,
                                          permutations=permutations, rng=rng)
            shared = len(set(cluster_members[a]) & set(cluster_members[b]))
            rows.append({
                "cluster_A": a,
                "cluster_B": b,
                "n_A": len(cluster_bins[a]),
                "n_B": len(cluster_bins[b]),
                "shared_genes": shared,
                "overlap_coefficient": round(ov, 4),
                "separation_z": round(z, 2) if np.isfinite(z) else np.nan,
                "separation_p": p_sep,
                "separation_stars": get_significance_stars(p_sep) if np.isfinite(p_sep) else "na",
            })
        if rows:
            tables["Key comparisons"] = pd.DataFrame(rows)

    if full_matrix and len(names) > 1:
        ov = np.eye(len(names))
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                ov[i, j] = ov[j, i] = overlap_coefficient(
                    cluster_bins[names[i]], cluster_bins[names[j]], nbins, smoothing)
        tables["Overlap coefficient matrix"] = pd.DataFrame(np.round(ov, 4), index=names, columns=names)

        shared = np.zeros((len(names), len(names)), dtype=int)
        for i in range(len(names)):
            si = set(cluster_members[names[i]])
            shared[i, i] = len(si)
            for j in range(i + 1, len(names)):
                shared[i, j] = shared[j, i] = len(si & set(cluster_members[names[j]]))
        tables["Shared genes matrix"] = pd.DataFrame(shared, index=names, columns=names)

    return tables


def write_statistics(tables: Dict[str, pd.DataFrame], output_path: str) -> Optional[str]:
    """Write the statistics tables next to the figure (xlsx if possible, else csv)."""
    if not tables or not output_path:
        return None
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    try:
        with pd.ExcelWriter(output_path) as writer:
            for sheet, df in tables.items():
                keep_index = "matrix" in sheet.lower()
                df.to_excel(writer, sheet_name=sheet[:31], index=keep_index)
        print(f"[INFO] Cluster landscape statistics written: {output_path}")
        return output_path
    except Exception as exc:
        print(f"[WARN] Could not write statistics workbook ({exc}); falling back to CSV.")
        base = os.path.splitext(output_path)[0]
        written = []
        for sheet, df in tables.items():
            path = f"{base}_{sheet.replace(' ', '_').lower()}.csv"
            df.to_csv(path, index="matrix" in sheet.lower())
            written.append(path)
        return written[0] if written else None


# ----------------------------
# Core plotting logic
# ----------------------------

def plot_cluster_landscape(
    reference_order: List[str],
    clusters_dict: Dict[str, List[str]],
    output_file: str,
    cmap_obj=None,
    bins: int = DEFAULT_BINS,
    smoothing: float = DEFAULT_SMOOTHING,
    color_cap: float = DEFAULT_COLOR_CAP,
    reference_cluster: str = "",
    show_overlap_column: bool = True,
    row_groups=None,
    frame_clusters=None,
    core_share: float = DEFAULT_CORE_SHARE,
    permutations: int = DEFAULT_PERMUTATIONS,
    height_per_cluster: float = 0.30,
    label_fontsize: int = 10,
    dpi: int = 300,
    show: bool = True,
    title: str = "Cluster landscape along the reordered genome axis",
    seed: int = 0,
) -> Tuple[str, Dict[str, float], Dict[str, float]]:
    """Draw the landscape. Returns (output_file, overlap values, overlap P values)."""
    reference_index = {locus: i for i, locus in enumerate(reference_order)}
    n_genes = len(reference_order)
    nbins = max(4, int(bins))
    layout_names = list(clusters_dict.keys())
    names = [n for n in layout_names if not str(n).startswith(SPACER_PREFIX)]

    cluster_bins: Dict[str, np.ndarray] = {}
    cluster_members: Dict[str, List[str]] = {}
    for name in names:
        positions = [reference_index[t] for t in clusters_dict[name] if t in reference_index]
        cluster_members[name] = [t for t in clusters_dict[name] if t in reference_index]
        cluster_bins[name] = _bin_index(np.asarray(positions, float), n_genes, nbins)

    ref_name = str(reference_cluster or "").strip()
    if ref_name:
        matched = {n.lower(): n for n in names}.get(ref_name.lower())
        if matched is None:
            print(f"[WARN] Reference cluster '{ref_name}' is not plotted; overlap column disabled.")
            ref_name = ""
        else:
            ref_name = matched
    elif show_overlap_column and names:
        ref_name = names[-1]
        print(f"[INFO] No reference cluster given; using the last plotted cluster '{ref_name}'.")

    overlaps: Dict[str, float] = {}
    overlap_p: Dict[str, float] = {}
    if show_overlap_column and ref_name:
        rng = np.random.default_rng(seed)
        for name in names:
            if name == ref_name:
                overlaps[name] = 1.0
                overlap_p[name] = 0.0
                continue
            ov, p = overlap_coefficient_p(
                cluster_bins[name], cluster_bins[ref_name], n_genes, nbins, smoothing,
                permutations=permutations, rng=rng)
            overlaps[name] = ov
            overlap_p[name] = p

    # Row layout: optional group brackets insert a blank spacer row between groups.
    groups = _parse_groups(row_groups)
    if groups:
        ordered: List[Optional[str]] = []
        spans: List[Tuple[str, int, int]] = []
        lower_to_actual = {n.lower(): n for n in names}
        used = set()
        for label, members in groups:
            start = len(ordered)
            for member in members:
                actual = lower_to_actual.get(str(member).strip().lower())
                if actual is None:
                    print(f"[WARN] Row group '{label}' refers to a cluster that is not plotted: {member}")
                    continue
                ordered.append(actual)
                used.add(actual)
            if len(ordered) > start:
                spans.append((label, start, len(ordered) - 1))
                ordered.append(None)
        leftover = [n for n in names if n not in used]
        if leftover:
            if ordered and ordered[-1] is None:
                ordered = ordered[:-1]
                ordered.append(None)
            start = len(ordered)
            ordered.extend(leftover)
            spans.append(("Other", start, len(ordered) - 1))
        while ordered and ordered[-1] is None:
            ordered = ordered[:-1]
    else:
        # A spacer column from the workbook becomes a None, i.e. a group break.
        ordered = [None if str(n).startswith(SPACER_PREFIX) else n for n in layout_names]
        while ordered and ordered[0] is None:
            ordered = ordered[1:]
        while ordered and ordered[-1] is None:
            ordered = ordered[:-1]
        spans = []

    # ------------------------------------------------------------------
    # Layout. Each cluster is drawn as its own framed strip rather than as one
    # row of a single heatmap: strips separated by a hairline gap read as
    # independent tracks, which is what a reader needs when comparing rows that
    # contain very different numbers of genes. A None in `ordered` (a group
    # break, or a spacer column in the workbook) widens that gap instead of
    # inserting an empty row. All geometry is in inches and converted at the
    # end, so row height does not drift with the number of rows.
    # ------------------------------------------------------------------
    cap = max(0.2, float(color_cap))
    fig_w = 14.0
    h_row = max(0.12, float(height_per_cluster))
    gap_row = max(0.02, h_row * 0.12)
    gap_group = max(0.18, h_row * 0.85)
    top_in = 0.62 if title else 0.24
    bot_in = 1.26

    row_data, row_label = {}, {}
    for name in ordered:
        if name is None:
            continue
        row_data[name] = _enrichment(cluster_bins[name], nbins, smoothing)
        row_label[name] = f"{name}  ({len(cluster_members[name])})"

    plan, pending_break, first = [], False, True
    for name in ordered:
        if name is None:
            pending_break = True
            continue
        plan.append((0.0 if first else (gap_group if pending_break else gap_row), name))
        pending_break, first = False, False
    if not plan:
        raise ValueError("Nothing to plot: every cluster column resolved to a spacer.")
    n_draw = len(plan)
    fig_h = top_in + bot_in + n_draw * h_row + sum(g for g, _ in plan)

    char_w = 0.56 * float(label_fontsize) / 72.0
    label_block = char_w * max(len(v) for v in row_label.values()) + 0.10
    group_block = (char_w * max(len(g[0]) for g in spans) + 0.30) if spans else 0.0
    lab_r = 0.12 + group_block + label_block
    hm_l = lab_r + 0.14
    bar_on = bool(show_overlap_column and ref_name)
    bar_w = 1.15 if bar_on else 0.0
    # 1.75 in on the right holds the colorbar, its label, and the gap that
    # keeps a bar annotation reading "1.00 ***" clear of the colorbar.
    hm_r = fig_w - (bar_w + 0.34 if bar_on else 0.0) - 1.75
    if hm_r - hm_l < 2.0:
        hm_l = min(hm_l, fig_w * 0.45)
        hm_r = max(hm_l + 2.0, hm_r)

    fig = plt.figure(figsize=(fig_w, fig_h))
    norm = TwoSlopeNorm(vmin=-cap, vcenter=0.0, vmax=cap)
    fx, fy = lambda v: v / fig_w, lambda v: v / fig_h

    y = fig_h - top_in - h_row
    axes_by_name, y_by_name, im = {}, {}, None
    for k_i, (gap, name) in enumerate(plan):
        if k_i:
            y -= gap + h_row
        ax = fig.add_axes([fx(hm_l), fy(y), fx(hm_r - hm_l), fy(h_row)])
        im = ax.imshow(row_data[name][None, :], aspect="auto", cmap=cmap_obj, norm=norm,
                       extent=[0, 100, 0, 1], interpolation="nearest")
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_linewidth(0.8); spine.set_color("0.15")
        fig.text(fx(lab_r), fy(y + h_row / 2.0), row_label[name], ha="right", va="center",
                 fontsize=int(label_fontsize), color="black")
        axes_by_name[name] = ax
        y_by_name[name] = y

        if bar_on:
            axb = fig.add_axes([fx(hm_r + 0.34), fy(y), fx(bar_w), fy(h_row)])
            axb.set_xlim(0, 1.0); axb.set_ylim(0, 1)
            axb.set_xticks([]); axb.set_yticks([])
            for spine in axb.spines.values():
                spine.set_visible(False)
            value = overlaps.get(name, np.nan)
            if np.isfinite(value):
                axb.barh([0.5], [value], height=0.62, color="0.62", zorder=3)
                p = overlap_p.get(name, np.nan)
                stars = ("***" if p < 0.005 else ("*" if p < 0.05 else "")) if np.isfinite(p) else ""
                axb.text(min(value + 0.03, 1.0), 0.5, f"{value:.2f} {stars}".strip(),
                         va="center", ha="left", fontsize=max(5.0, int(label_fontsize) - 2.5),
                         color="0.35", clip_on=False)

    last_ax = axes_by_name[plan[-1][1]]
    last_ax.set_xticks([0, 25, 50, 75, 100])
    last_ax.set_xticklabels(["0", "25", "50", "75", "100"], fontsize=int(label_fontsize))
    last_ax.tick_params(axis="x", length=3.5, width=0.9, color="0.15", pad=3)
    last_ax.set_xlabel("position along the reordered genome axis (%)",
                       fontsize=int(label_fontsize), labelpad=5)

    if bar_on:
        y_top = y_by_name[plan[0][1]] + h_row
        fig.text(fx(hm_r + 0.34 + bar_w / 2.0), fy(y_top + 0.10),
                 f"Overlap with\n{ref_name}", ha="center", va="bottom",
                 fontsize=max(6.0, int(label_fontsize) - 1))

    if spans:
        import matplotlib.lines as _mlines
        bx = fx(0.12 + group_block)
        for label, start, end in spans:
            seg = [n for n in ordered[start:end + 1] if n is not None]
            if not seg:
                continue
            y_hi = y_by_name[seg[0]] + h_row
            y_lo = y_by_name[seg[-1]]
            fig.add_artist(_mlines.Line2D([bx, bx], [fy(y_lo), fy(y_hi)], lw=2.2, color="0.3"))
            fig.text(bx - 0.006, fy((y_lo + y_hi) / 2.0), label, ha="right", va="center",
                     fontsize=max(6.0, int(label_fontsize) - 1), fontweight="bold")

    frames = _parse_frames(frame_clusters)
    if frames:
        y_top = y_by_name[plan[0][1]] + h_row
        y_bot = y_by_name[plan[-1][1]]
        axf = fig.add_axes([fx(hm_l), fy(y_bot), fx(hm_r - hm_l), fy(y_top - y_bot)], zorder=5)
        axf.set_xlim(0, nbins); axf.set_ylim(0, 1); axf.set_axis_off()
        axf.patch.set_alpha(0.0)
        for i, (label, cluster_name) in enumerate(frames):
            actual = {n.lower(): n for n in names}.get(str(cluster_name).strip().lower())
            if actual is None:
                print(f"[WARN] Frame '{label}' refers to a cluster that is not plotted: {cluster_name}")
                continue
            lo, hi = core_band(cluster_bins[actual], core_share)
            color = DEFAULT_FRAME_COLORS[i % len(DEFAULT_FRAME_COLORS)]
            axf.add_patch(plt.Rectangle((lo, 0), max(hi - lo, 1), 1, fill=False,
                                        edgecolor=color, lw=2.0, clip_on=False))
            axf.text((lo + hi) / 2.0, 1.004, label, ha="center", va="bottom", color=color,
                     fontsize=int(label_fontsize), fontweight="bold")

    cb_x = hm_r + (0.34 + bar_w + 0.85 if bar_on else 0.40)
    cb_h = max(0.9, min(2.2, 0.35 * fig_h))
    cb_y = (fig_h - cb_h) / 2.0
    cbar = fig.colorbar(im, cax=fig.add_axes([fx(cb_x), fy(cb_y), fx(0.16), fy(cb_h)]),
                        ticks=[-cap, -1, 0, 1, cap])
    cbar.ax.set_yticklabels([f"1/{2 ** cap:.0f}", "1/2", "1", "2", f"{2 ** cap:.0f}"])
    cbar.ax.tick_params(labelsize=max(5.0, int(label_fontsize) - 2))
    cbar.outline.set_linewidth(0.6)
    cbar.set_label("genes per bin, obs/exp", fontsize=max(6.0, int(label_fontsize) - 1))

    if title:
        fig.text(0.012, 0.955, title, fontsize=int(label_fontsize) + 3, fontweight="bold", va="bottom")
    footer = (f"Resolution: {nbins} equal-count bins, Gaussian smoothing {float(smoothing):.1f} bins.")
    if show_overlap_column and ref_name:
        footer += (f"   Bar column: overlap coefficient with '{ref_name}' along this axis "
                   f"(0 = disjoint, 1 = identical); *** P < 0.005, * P < 0.05 against random "
                   f"gene sets of the same size, {int(permutations)} permutations.")
    fig.text(0.012, 0.025, footer, fontsize=int(label_fontsize) - 2, color="0.35", wrap=True)

    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
    fig.savefig(output_file, dpi=int(dpi))
    print(f"[INFO] Cluster landscape written: {output_file}")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return output_file, overlaps, overlap_p


def run_gene_cluster_localization_heatmap_from_workbook(
    workbook_path: str,
    sheet_name: str = "Locus Tags",
    output_file: Optional[str] = None,
    colormap: str = "RdBu_r",
    custom_cmaps_xlsx: str = "",
    bins: int = DEFAULT_BINS,
    smoothing: float = DEFAULT_SMOOTHING,
    color_cap: float = DEFAULT_COLOR_CAP,
    reference_cluster: str = "",
    show_overlap_column: bool = True,
    row_groups=None,
    frame_clusters=None,
    core_share: float = DEFAULT_CORE_SHARE,
    key_comparisons=None,
    full_matrix: bool = True,
    permutations: int = DEFAULT_PERMUTATIONS,
    statistics_file: Optional[str] = None,
    height_per_cluster: float = 0.30,
    label_fontsize: int = 10,
    dpi: int = 300,
    show: bool = True,
    include_clusters: Optional[List[str]] = None,
    title: str = "Cluster landscape along the reordered genome axis",
    seed: int = 0,
    # Accepted for backward compatibility with saved settings; no longer used.
    sigma: Optional[float] = None,
    spread_factor: Optional[float] = None,
    cmap_min_rel: Optional[float] = None,
    cmap_max_rel: Optional[float] = None,
    x_min=None,
    x_max=None,
    y_min=None,
    y_max=None,
) -> str:
    """Convenience wrapper used by the pipeline."""
    reference_order, clusters = load_locus_tags_sheet(
        workbook_path, sheet_name=sheet_name, include_clusters=include_clusters)

    if output_file is None or str(output_file).strip() == "":
        output_file = os.path.join(os.path.dirname(workbook_path), "cluster_landscape.png")

    for legacy_name, legacy_value in (("sigma", sigma), ("spread factor", spread_factor),
                                      ("min rel", cmap_min_rel), ("max rel", cmap_max_rel)):
        if legacy_value is not None:
            print(f"[INFO] Parameter '{legacy_name}' belonged to the previous per-row density display "
                  f"and is ignored by the cluster landscape plot.")
            break

    custom_maps = _load_custom_colormaps(custom_cmaps_xlsx) if custom_cmaps_xlsx else {}
    cmap_obj = _choose_cmap(colormap, custom_maps)

    path, _, _ = plot_cluster_landscape(
        reference_order=reference_order,
        clusters_dict=clusters,
        output_file=output_file,
        cmap_obj=cmap_obj,
        bins=bins,
        smoothing=smoothing,
        color_cap=color_cap,
        reference_cluster=reference_cluster,
        show_overlap_column=show_overlap_column,
        row_groups=row_groups,
        frame_clusters=frame_clusters,
        core_share=core_share,
        permutations=permutations,
        height_per_cluster=height_per_cluster,
        label_fontsize=label_fontsize,
        dpi=dpi,
        show=show,
        title=title,
        seed=seed,
    )

    key_pairs = _parse_pairs(key_comparisons)
    if key_pairs or full_matrix:
        reference_index = {locus: i for i, locus in enumerate(reference_order)}
        n_genes = len(reference_order)
        nbins = max(4, int(bins))
        cluster_bins = {}
        cluster_members = {}
        for name, loci in clusters.items():
            kept = [t for t in loci if t in reference_index]
            cluster_members[name] = kept
            cluster_bins[name] = _bin_index(
                np.asarray([reference_index[t] for t in kept], float), n_genes, nbins)
        tables = compute_pair_statistics(
            cluster_bins, cluster_members, n_genes, nbins, smoothing,
            key_pairs=key_pairs, full_matrix=full_matrix,
            permutations=permutations, seed=seed)
        if statistics_file is None or str(statistics_file).strip() == "":
            statistics_file = os.path.splitext(output_file)[0] + "_statistics.xlsx"
        write_statistics(tables, statistics_file)

    return path
