# codonpipe/gene_cluster_heatmap.py
"""Gene-cluster localization heatmap along a reordered genome axis (KS test).

Gene-cluster localization plotting along a reordered genome axis.

Expected workbook structure (default: sheet 'Locus Tags'):
- Column 1: reordered genome locus tags (defines genome axis)
- Columns 2..N: each column is a gene cluster listing locus tags (header in first row)

Output:
- PNG heatmap, one row per gene cluster, colored by smoothed (Gaussian) relative density.
- Row labels include (n=..., KS significance stars).

Notes:
- Density is normalized per cluster to [0,1] to highlight within-cluster peaks.
- KS test compares cluster positions to a uniform distribution along the genome axis.
"""

from __future__ import annotations

import os
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from .colormaps import get_cmap_any
from scipy.stats import ks_1samp, uniform


# ----------------------------
# Colormap helpers
# ----------------------------

def _load_custom_colormaps(xlsx_path: str) -> Dict[str, LinearSegmentedColormap]:
    """Load custom colormaps from an Excel file.

    Excel format:
      - one sheet per colormap (sheet name = cmap name)
      - no header; columns: [index, R, G, B] with RGB in 0–255 or 0–1
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

    Supports built-in 'parula' / 'parula_r' (see :mod:`codonpipe.colormaps`).
    """
    return get_cmap_any(name, custom_maps=custom_maps, fallback="plasma")


# ----------------------------
# Core plotting logic
# ----------------------------

def get_significance_stars(p: float) -> str:
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    else:
        return 'ns'


def _clean_tag(x) -> str:
    if x is None:
        return ""
    s = str(x)
    if s.lower() == "nan":
        return ""
    return s.strip().strip('"').strip("'").strip()


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

    # First column = genome axis
    ref_col = df.columns[0]
    reference_order = [_clean_tag(v) for v in df[ref_col].tolist()]
    reference_order = [v for v in reference_order if v != ""]
    if not reference_order:
        raise ValueError(f"No locus tags found in first column of sheet '{sheet_name}'.")

    # Remaining columns = clusters
    clusters: Dict[str, List[str]] = {}
    data_cols = list(df.columns[1:])
    if include_clusters is not None:
        lower_to_actual = {str(col).strip().lower(): str(col) for col in data_cols}
        ordered_cols = []
        seen_cols = set()
        for req in include_clusters:
            mapped = lower_to_actual.get(str(req).strip().lower())
            if mapped is not None and mapped not in seen_cols:
                ordered_cols.append(mapped)
                seen_cols.add(mapped)
        if not ordered_cols:
            raise ValueError("No requested clusters were found in the locus-tags sheet.")
    else:
        ordered_cols = [str(col) for col in data_cols]

    for col_name in ordered_cols:
        vals = [_clean_tag(v) for v in df[col_name].tolist()]
        vals = [v for v in vals if v != ""]
        clusters[col_name] = vals

    return reference_order, clusters


def plot_gene_density_heatmap(
    reference_order: List[str],
    clusters_dict: Dict[str, List[str]],
    sigma: float,
    spread_factor: float,
    cmap_obj,
    height_scale: float,
    label_fontsize: int,
    output_file: str,
    cmap_min_rel: float,
    cmap_max_rel: float,
    dpi: int = 300,
    show: bool = True,
    title: str = "Gene Cluster Localization Heatmap (with KS Test)",
    x_min=None,
    x_max=None,
    y_min=None,
    y_max=None,
) -> str:
    """Create and save the heatmap. Returns the output_file path."""
    reference_index = {locus: i for i, locus in enumerate(reference_order)}
    genome_length = len(reference_order)
    cluster_names = list(clusters_dict.keys())

    density_matrix = np.zeros((len(cluster_names), genome_length), dtype=float)
    cluster_labels: List[str] = []

    adjusted_sigma = float(sigma) * float(spread_factor)
    window = np.arange(genome_length)

    for i, cluster_name in enumerate(cluster_names):
        loci = clusters_dict[cluster_name]
        indices = [reference_index[l] for l in loci if l in reference_index]

        density = np.zeros(genome_length, dtype=float)
        for idx in indices:
            gaussian = np.exp(-0.5 * ((window - idx) / adjusted_sigma) ** 2)
            density += gaussian

        max_val = float(np.max(density)) if density.size else 0.0
        if max_val > 0:
            density_matrix[i, :] = density / max_val
        else:
            density_matrix[i, :] = 0.0

        # KS test vs uniform distribution
        norm_positions = (np.array(indices, dtype=float) / float(genome_length)) if len(indices) > 0 else np.array([])
        if norm_positions.size > 0:
            _, p_value = ks_1samp(norm_positions, uniform.cdf)
        else:
            p_value = 1.0

        label = f"{cluster_name} (n={len(indices)}, {get_significance_stars(p_value)})"
        cluster_labels.append(label)

    # Clamp colormap limits
    vmin = max(0.0, min(1.0, float(cmap_min_rel)))
    vmax = max(vmin + 1e-6, min(1.0, float(cmap_max_rel)))

    fig_height = float(height_scale) * len(cluster_names) + 1.5
    plt.figure(figsize=(14, fig_height))
    im = plt.imshow(
        density_matrix,
        aspect='auto',
        cmap=cmap_obj,
        interpolation='nearest',
        vmin=vmin,
        vmax=vmax
    )
    plt.colorbar(im, label='Relative Gene Density')
    plt.yticks(np.arange(len(cluster_labels)), cluster_labels, fontsize=label_fontsize)
    plt.xticks([], [])
    try:
        if x_min not in (None, "") or x_max not in (None, ""):
            plt.xlim(left=(float(x_min) if x_min not in (None, "") else None), right=(float(x_max) if x_max not in (None, "") else None))
    except Exception:
        pass
    try:
        if y_min not in (None, "") or y_max not in (None, ""):
            cur0, cur1 = plt.ylim()
            logical_min = min(cur0, cur1)
            logical_max = max(cur0, cur1)
            ylo = float(y_min) if y_min not in (None, "") else logical_min
            yhi = float(y_max) if y_max not in (None, "") else logical_max
            if cur0 > cur1:
                plt.ylim(yhi, ylo)
            else:
                plt.ylim(ylo, yhi)
    except Exception:
        pass
    plt.xlabel("Genome Position (Ordered Locus Tags)")
    plt.title(title)
    plt.tight_layout()

    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
    plt.savefig(output_file, dpi=int(dpi))

    if show:
        plt.show()
    else:
        plt.close()

    return output_file


def run_gene_cluster_localization_heatmap_from_workbook(
    workbook_path: str,
    sheet_name: str = "Locus Tags",
    output_file: Optional[str] = None,
    colormap: str = "plasma",
    custom_cmaps_xlsx: str = "",
    sigma: float = 10,
    spread_factor: float = 5,
    height_per_cluster: float = 0.3,
    label_fontsize: int = 10,
    dpi: int = 300,
    cmap_min_rel: float = 0.2,
    cmap_max_rel: float = 1.0,
    show: bool = True,
    include_clusters: Optional[List[str]] = None,
    title: str = "Gene Cluster Localization Heatmap (with KS Test)",
    x_min=None,
    x_max=None,
    y_min=None,
    y_max=None,
) -> str:
    """Convenience wrapper used by the pipeline."""
    reference_order, clusters = load_locus_tags_sheet(workbook_path, sheet_name=sheet_name, include_clusters=include_clusters)

    if output_file is None or str(output_file).strip() == "":
        output_file = os.path.join(os.path.dirname(workbook_path), "gene_cluster_heatmap_KS.png")

    custom_maps = _load_custom_colormaps(custom_cmaps_xlsx) if custom_cmaps_xlsx else {}
    cmap_obj = _choose_cmap(colormap, custom_maps)

    return plot_gene_density_heatmap(
        reference_order=reference_order,
        clusters_dict=clusters,
        sigma=sigma,
        spread_factor=spread_factor,
        cmap_obj=cmap_obj,
        height_scale=height_per_cluster,
        label_fontsize=label_fontsize,
        output_file=output_file,
        cmap_min_rel=cmap_min_rel,
        cmap_max_rel=cmap_max_rel,
        dpi=dpi,
        show=show,
        title=title,
        x_min=x_min,
        x_max=x_max,
        y_min=y_min,
        y_max=y_max,
    )
