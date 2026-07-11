"""Cluster overlay plotting for 2D embeddings (UMAP / t-SNE / PCA).

Inputs:
  - Excel workbook produced by Clustering_Pipeline.py
  - Coordinate sheet (default: "Genes reordered")
  - Binary membership sheet (default: "Binary")

Outputs:
  - A multi-panel figure showing one panel per cluster,
    with optional density or enrichment coloring.
  - Optionally, a genome-axis gene-cluster localization heatmap.
"""


# ======================== USER SETTINGS ===========================

EXCEL_PATH   = r""
UMAP_SHEET   = "Genes reordered"
DATA_SHEET   = "Binary"

USE_FILE_DIALOG = True

INCLUDE_COLUMNS  = None
HIGHLIGHT_COLUMNS = None

PLOT_MODE = "grid"

MAX_NROWS    = 2
PANEL_W_IN   = 5.0
PANEL_H_IN   = 5.0
DENSITY_XMIN = None
DENSITY_XMAX = None
DENSITY_YMIN = None
DENSITY_YMAX = None

X_LABEL = "t-SNE 1"
Y_LABEL = "t-SNE 2"

SHOW_BACKGROUND_SMALL   = True
BACKGROUND_COLOR        = "grey"
BACKGROUND_ALPHA        = 0.18
BACKGROUND_SIZE         = 8
BACKGROUND_LW           = 0.0

HIGHLIGHT_MARKER         = "o"
HIGHLIGHT_EDGE_COLOR     = "black"
HIGHLIGHT_EDGE_WIDTH     = 0.6
HIGHLIGHT_ALPHA          = 0.95

ALL_SCATTER_POINT_SIZE      = 3.5
ALL_SCATTER_ALPHA           = 0.7
ALL_SCATTER_EDGE_WIDTH      = 0.15
ALL_SCATTER_DENSITY_NBINS   = 150
ALL_SCATTER_DENSITY_SIGMA   = 4.0
ALL_SCATTER_DENSITY_USE_LOG = True
ALL_SCATTER_DENSITY_MIN_REL = 0.00
ALL_SCATTER_DENSITY_MAX_REL = 1.00

CLUSTER_SCATTER_PRESETS = [
    {"MIN_N": 0,    "MAX_N": 100,  "POINT_SIZE": 40.0, "ALPHA": 0.95, "EDGE_WIDTH": 0.6,  "DENSITY_NBINS": 50,  "DENSITY_SIGMA": 2.0, "DENSITY_USE_LOG": True, "DENSITY_MIN_REL": 0.00, "DENSITY_MAX_REL": 1.00},
    {"MIN_N": 100,  "MAX_N": 250,  "POINT_SIZE": 35.0, "ALPHA": 0.95, "EDGE_WIDTH": 0.6,  "DENSITY_NBINS": 80,  "DENSITY_SIGMA": 2.5, "DENSITY_USE_LOG": True, "DENSITY_MIN_REL": 0.00, "DENSITY_MAX_REL": 1.00},
    {"MIN_N": 250,  "MAX_N": 500,  "POINT_SIZE": 20.0, "ALPHA": 0.95, "EDGE_WIDTH": 0.3,  "DENSITY_NBINS": 110, "DENSITY_SIGMA": 3.0, "DENSITY_USE_LOG": True, "DENSITY_MIN_REL": 0.00, "DENSITY_MAX_REL": 1.00},
    {"MIN_N": 500,  "MAX_N": 1500, "POINT_SIZE": 14,   "ALPHA": 0.9,  "EDGE_WIDTH": 0.2,  "DENSITY_NBINS": 150, "DENSITY_SIGMA": 4.0, "DENSITY_USE_LOG": True, "DENSITY_MIN_REL": 0.00, "DENSITY_MAX_REL": 1.00},
    {"MIN_N": 1500, "MAX_N": 2500, "POINT_SIZE": 6.0,  "ALPHA": 0.85, "EDGE_WIDTH": 0.10, "DENSITY_NBINS": 180, "DENSITY_SIGMA": 4.5, "DENSITY_USE_LOG": True, "DENSITY_MIN_REL": 0.00, "DENSITY_MAX_REL": 1.00},
    {"MIN_N": 2500, "MAX_N": None, "POINT_SIZE": 3.5,  "ALPHA": 0.8,  "EDGE_WIDTH": 0.08, "DENSITY_NBINS": 150, "DENSITY_SIGMA": 4.0, "DENSITY_USE_LOG": True, "DENSITY_MIN_REL": 0.00, "DENSITY_MAX_REL": 1.00},
]

DRAW_GRID        = True
GRID_MAJOR_STEP  = None
GRID_MAJOR_ALPHA = 0.55
GRID_MAJOR_LW    = 1.0
GRID_MAJOR_COLOR = "black"

CAPTION_SIZE         = 12
CAPTION_WEIGHT       = "bold"
CAPTION_PAD          = 8
CAPTION_WRAP_ENABLED     = True
CAPTION_WRAP_CHAR_WIDTH  = 35
CAPTION_MAX_NAME_LINES   = 2

SUBPLOT_WSPACE = 0.20
SUBPLOT_HSPACE = 0.30

SHOW_COLORBAR = True

CUSTOM_CMAPS_XLSX = ""
DENSITY_CMAP_NAME = "plasma_r"

# ===================== Gene-cluster localization heatmap =====================
# Uses the workbook sheet produced by Clustering_Pipeline.py (default: 'Locus Tags').
# The sheet must be:
#   - Column 1: reordered genome locus tags (genome axis)
#   - Other columns: each column lists locus tags belonging to that cluster
RUN_GENE_CLUSTER_LOCALIZATION_HEATMAP = True
LOCUS_TAGS_SHEET = "Locus Tags"

GCHM_COLORMAP = "plasma"          # can be 'plasma', 'viridis', or a custom colormap sheet name
GCHM_CUSTOM_CMAPS_XLSX = CUSTOM_CMAPS_XLSX  # "" disables custom colormaps

GCHM_SIGMA = 10
GCHM_SPREAD_FACTOR = 5
GCHM_HEIGHT_PER_CLUSTER = 0.3
GCHM_LABEL_FONTSIZE = 10
GCHM_DPI = 300

# Colormap range in relative density units (0..1)
GCHM_CMAP_MIN_REL = 0.2
GCHM_CMAP_MAX_REL = 1.0

# Output file (saved next to the workbook)
GCHM_OUTPUT_FILENAME = "gene_cluster_heatmap_KS.png"
GCHM_OUTPUT_BASENAME = "gene_cluster_heatmap_KS"

# If running from Clustering_Pipeline, you may prefer show=False to avoid popping a second window
GCHM_SHOW_FIG = True
GCHM_INCLUDE_COLUMNS = None
GCHM_XMIN = None
GCHM_XMAX = None
GCHM_YMIN = None
GCHM_YMAX = None

# ---- Density vs enrichment coloring ----
COLOR_MODE = "enrichment"   # "density" or "enrichment"

# ---- Enrichment coloring options (used only when COLOR_MODE = "enrichment") ----
ENRICHMENT_SCALE = "log2"  # "ratio" or "log2"
ENRICHMENT_STYLE = "unilateral_diverging"
ENRICHMENT_CMAP_NAME = "plasma_r"
ENRICHMENT_VMAX = 3
ENRICHMENT_PERCENTILE = 99.0
ENRICHMENT_SYMMETRIC = True
ENRICHMENT_EPS = 1e-12

# CHANGED: force the top-left "All genes" panel to use DENSITY coloring,
# even when COLOR_MODE="enrichment".
ALL_GENES_PANEL_MODE = "density"  # "density" or "background"

ENRICHMENT_FORCE_NEUTRAL_BELOW_CENTER = False
ENRICHMENT_NEUTRAL_RGBA = (0.85, 0.85, 0.85, 1.0)

# ---- Output ----
SAVE_FIG = True
PNG_DPI = 220
FIGURE_FORMAT = "png"
RUN_2D_DENSITY_PLOTS = True
INCLUDE_GENOMIC_DENSITY_PANEL = True

# ---- Debugging ----
DEBUG_PRINT_CONFIG = False  # Set False to silence debug prints

# ---- Optional codon-usage multi-panel bar plot ----
CODON_USAGE_PLOT_MODE = "NONE"   # "ACU", "RCU", "Z", or "NONE"
CODON_USAGE_WORKBOOK = ""
CODON_USAGE_ORDER_SHEET = "Codons reordered"
CODON_USAGE_PNG_DPI = 220
CODON_USAGE_PANEL_W_IN = PANEL_W_IN
CODON_USAGE_PANEL_H_IN = PANEL_H_IN
CODON_USAGE_SHOW_FIG = True
CODON_USAGE_OUTPUT_FORMAT = FIGURE_FORMAT
CODON_USAGE_OUTPUT_BASENAME = "Average codon usage per cluster"
CODON_USAGE_XMIN = None
CODON_USAGE_XMAX = None
CODON_USAGE_YMIN = None
CODON_USAGE_YMAX = None

# ================================================================

import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Optional, Tuple

# Optional: allow Clustering_Pipeline.py to enforce consistent colormaps across the whole pipeline
# (it sets CODONPIPE_DENSITY_CMAP / CODONPIPE_ENRICHMENT_CMAP before importing this script).
_env_density = os.environ.get("CODONPIPE_DENSITY_CMAP", "").strip()
_env_enrich  = os.environ.get("CODONPIPE_ENRICHMENT_CMAP", "").strip()

if _env_density:
    DENSITY_CMAP_NAME = _env_density
if _env_enrich:
    ENRICHMENT_CMAP_NAME = _env_enrich
elif _env_density:
    # If only one cmap is specified, use it for both density + enrichment
    ENRICHMENT_CMAP_NAME = _env_density


# Optional: allow Clustering_Pipeline.py to drive gene-cluster heatmap parameters via env vars
# (keeps the bridge minimal; Plotting_Pipeline remains runnable standalone).
def _env_bool(key: str, default: bool) -> bool:
    v = os.environ.get(key, "").strip().lower()
    if v == "":
        return default
    return v in ("1", "true", "yes", "y", "on")

def _env_float(key: str, default: float) -> float:
    v = os.environ.get(key, "").strip()
    if v == "":
        return default
    try:
        return float(v)
    except Exception:
        return default

def _env_int(key: str, default: int) -> int:
    v = os.environ.get(key, "").strip()
    if v == "":
        return default
    try:
        return int(float(v))
    except Exception:
        return default

RUN_GENE_CLUSTER_LOCALIZATION_HEATMAP = _env_bool("CODONPIPE_GCHM_ENABLE", RUN_GENE_CLUSTER_LOCALIZATION_HEATMAP)
LOCUS_TAGS_SHEET = os.environ.get("CODONPIPE_GCHM_SHEET", LOCUS_TAGS_SHEET).strip() or LOCUS_TAGS_SHEET
GCHM_COLORMAP = os.environ.get("CODONPIPE_GCHM_COLORMAP", GCHM_COLORMAP).strip() or GCHM_COLORMAP
GCHM_CUSTOM_CMAPS_XLSX = os.environ.get("CODONPIPE_GCHM_CUSTOM_CMAPS_XLSX", GCHM_CUSTOM_CMAPS_XLSX).strip()
GCHM_SIGMA = _env_float("CODONPIPE_GCHM_SIGMA", GCHM_SIGMA)
GCHM_SPREAD_FACTOR = _env_float("CODONPIPE_GCHM_SPREAD_FACTOR", GCHM_SPREAD_FACTOR)
GCHM_HEIGHT_PER_CLUSTER = _env_float("CODONPIPE_GCHM_HEIGHT_PER_CLUSTER", GCHM_HEIGHT_PER_CLUSTER)
GCHM_LABEL_FONTSIZE = _env_int("CODONPIPE_GCHM_LABEL_FONTSIZE", GCHM_LABEL_FONTSIZE)
GCHM_DPI = _env_int("CODONPIPE_GCHM_DPI", GCHM_DPI)
GCHM_CMAP_MIN_REL = _env_float("CODONPIPE_GCHM_CMAP_MIN_REL", GCHM_CMAP_MIN_REL)
GCHM_CMAP_MAX_REL = _env_float("CODONPIPE_GCHM_CMAP_MAX_REL", GCHM_CMAP_MAX_REL)
GCHM_OUTPUT_FILENAME = os.environ.get("CODONPIPE_GCHM_OUTPUT_FILENAME", GCHM_OUTPUT_FILENAME).strip() or GCHM_OUTPUT_FILENAME
GCHM_SHOW_FIG = _env_bool("CODONPIPE_GCHM_SHOW", GCHM_SHOW_FIG)
GCHM_XMIN = os.environ.get("CODONPIPE_GCHM_XMIN", "").strip() or GCHM_XMIN
GCHM_XMAX = os.environ.get("CODONPIPE_GCHM_XMAX", "").strip() or GCHM_XMAX
GCHM_YMIN = os.environ.get("CODONPIPE_GCHM_YMIN", "").strip() or GCHM_YMIN
GCHM_YMAX = os.environ.get("CODONPIPE_GCHM_YMAX", "").strip() or GCHM_YMAX
CODON_USAGE_PLOT_MODE = os.environ.get("CODONPIPE_CODON_USAGE_PLOT_MODE", CODON_USAGE_PLOT_MODE).strip().upper() or CODON_USAGE_PLOT_MODE
CODON_USAGE_WORKBOOK = os.environ.get("CODONPIPE_CODON_USAGE_WORKBOOK", CODON_USAGE_WORKBOOK).strip() or CODON_USAGE_WORKBOOK
FIGURE_FORMAT = os.environ.get("CODONPIPE_FIGURE_FORMAT", FIGURE_FORMAT).strip().lstrip(".").lower() or FIGURE_FORMAT
CODON_USAGE_OUTPUT_FORMAT = os.environ.get("CODONPIPE_CODON_USAGE_OUTPUT_FORMAT", CODON_USAGE_OUTPUT_FORMAT).strip().lstrip(".").lower() or CODON_USAGE_OUTPUT_FORMAT
GCHM_OUTPUT_BASENAME = os.environ.get("CODONPIPE_GCHM_OUTPUT_BASENAME", GCHM_OUTPUT_BASENAME).strip() or GCHM_OUTPUT_BASENAME
CODON_USAGE_OUTPUT_BASENAME = os.environ.get("CODONPIPE_CODON_USAGE_OUTPUT_BASENAME", CODON_USAGE_OUTPUT_BASENAME).strip() or CODON_USAGE_OUTPUT_BASENAME
RUN_2D_DENSITY_PLOTS = _env_bool("CODONPIPE_RUN_2D_DENSITY_PLOTS", RUN_2D_DENSITY_PLOTS)
DENSITY_XMIN = os.environ.get("CODONPIPE_DENSITY_XMIN", "").strip() or DENSITY_XMIN
DENSITY_XMAX = os.environ.get("CODONPIPE_DENSITY_XMAX", "").strip() or DENSITY_XMAX
DENSITY_YMIN = os.environ.get("CODONPIPE_DENSITY_YMIN", "").strip() or DENSITY_YMIN
DENSITY_YMAX = os.environ.get("CODONPIPE_DENSITY_YMAX", "").strip() or DENSITY_YMAX
CODON_USAGE_XMIN = os.environ.get("CODONPIPE_CODON_USAGE_XMIN", "").strip() or CODON_USAGE_XMIN
CODON_USAGE_XMAX = os.environ.get("CODONPIPE_CODON_USAGE_XMAX", "").strip() or CODON_USAGE_XMAX
CODON_USAGE_YMIN = os.environ.get("CODONPIPE_CODON_USAGE_YMIN", "").strip() or CODON_USAGE_YMIN
CODON_USAGE_YMAX = os.environ.get("CODONPIPE_CODON_USAGE_YMAX", "").strip() or CODON_USAGE_YMAX

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = HERE
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from codonpipe.density_plot_core import run_density_dotplot
from codonpipe.gene_cluster_heatmap import run_gene_cluster_localization_heatmap_from_workbook



def _codon_order_from_workbooks(main_workbook_path: str, codon_usage_workbook: str, order_sheet: str = "Codons reordered") -> list:
    order = []
    try:
        xls_main = pd.ExcelFile(main_workbook_path)
        if order_sheet in xls_main.sheet_names:
            df = pd.read_excel(main_workbook_path, sheet_name=order_sheet)
            if not df.empty:
                vals = df.iloc[:, 0].astype(str).tolist()
                for v in vals:
                    s = str(v).strip()
                    if not s or s.lower() == "nan":
                        continue
                    if "_" not in s:
                        continue
                    if s == "Met_ATG":
                        s = "Met_AUG"
                    order.append(s)
    except Exception:
        pass

    # Append any codons present in the codon-usage workbook but absent from the reordered sheet.
    try:
        xls_cu = pd.ExcelFile(codon_usage_workbook)
        for sheet in ("ACU per cluster", "RCU per cluster", "z-scores per cluster", "ACU", "RCU", "ZCU", "Z"):
            if sheet in xls_cu.sheet_names:
                df = pd.read_excel(codon_usage_workbook, sheet_name=sheet)
                if not df.empty:
                    vals = [str(v).strip() for v in df.iloc[:, 0].tolist()]
                    vals = [v for v in vals if v and v.lower() != "nan"]
                    for v in vals:
                        if v not in order:
                            order.append(v)
                break
    except Exception:
        pass

    return order


def _read_described_table_sheet(workbook_path: str, sheet_name: str) -> pd.DataFrame:
    """Read a summary sheet that may start with a one-line description row."""
    for header_row in (0, 1):
        df = pd.read_excel(workbook_path, sheet_name=sheet_name, header=header_row)
        if df is None or df.empty:
            continue
        cols = [str(c).strip().lower() for c in df.columns]
        if "codon" in cols or "trna" in cols:
            return df
    # Fallback for unexpected but still valid sheets.
    return pd.read_excel(workbook_path, sheet_name=sheet_name)


def _read_codon_usage_summary_sheet(codon_usage_workbook: str, plot_mode: str) -> Tuple[str, pd.DataFrame]:
    mode = (plot_mode or "NONE").strip().upper()
    sheet_map = {
        "ACU": ["ACU per cluster", "ACU"],
        "RCU": ["RCU per cluster", "RCU"],
        "Z": ["z-scores per cluster", "ZCU", "Z"],
    }
    if mode not in sheet_map:
        raise ValueError("plot_mode must be one of ACU, RCU, Z.")
    xls = pd.ExcelFile(codon_usage_workbook)
    for sheet_name in sheet_map[mode]:
        if sheet_name in xls.sheet_names:
            df = _read_described_table_sheet(codon_usage_workbook, sheet_name)
            if df.empty:
                raise ValueError(f"Sheet {sheet_name!r} is empty in {codon_usage_workbook}.")
            return sheet_name, df
    raise ValueError(
        f"Could not find a summary sheet for mode {mode!r} in {codon_usage_workbook}. "
        f"Tried: {sheet_map[mode]}"
    )


def _plot_cluster_codon_usage_panels(main_workbook_path: str,
                                     codon_usage_workbook: str,
                                     plot_mode: str,
                                     include_columns,
                                     max_nrows: int,
                                     panel_w_in: float,
                                     panel_h_in: float,
                                     dpi: int = 220,
                                     show: bool = False,
                                     output_format: str = "png",
                                     output_basename: str = "",
                                     x_min=None,
                                     x_max=None,
                                     y_min=None,
                                     y_max=None) -> Optional[str]:
    raw_mode = (plot_mode or "NONE").strip()
    mode_map = {
        "ACU": "ACU",
        "ABS": "ACU",
        "ABSOLUTE": "ACU",
        "ABSOLUTE CODON USAGE": "ACU",
        "RCU": "RCU",
        "REL": "RCU",
        "RELATIVE": "RCU",
        "RELATIVE CODON USAGE": "RCU",
        "Z": "Z",
        "ZSCORE": "Z",
        "Z-SCORE": "Z",
        "Z SCORES": "Z",
        "Z-SCORES": "Z",
        "RELATIVE CODON USAGE Z-SCORES": "Z",
        "NONE": "NONE",
        "NO": "NONE",
        "NO PLOT": "NONE",
    }
    mode = mode_map.get(raw_mode.upper(), raw_mode.upper())
    if mode == "NONE":
        return None
    if not codon_usage_workbook or not os.path.exists(codon_usage_workbook):
        print(f"[WARN] Codon-usage workbook not found; skipping codon-usage panel plot:\n  {codon_usage_workbook}")
        return None

    sheet_name, df = _read_codon_usage_summary_sheet(codon_usage_workbook, mode)
    if "Codon" not in df.columns:
        raise ValueError(f"Sheet {sheet_name!r} must contain a first column named 'Codon'.")

    codon_order = _codon_order_from_workbooks(main_workbook_path, codon_usage_workbook, order_sheet=CODON_USAGE_ORDER_SHEET)
    if not codon_order:
        codon_order = [str(v).strip() for v in df["Codon"].tolist()]

    df["Codon"] = df["Codon"].astype(str).str.strip()
    df = df.drop_duplicates(subset=["Codon"], keep="first").set_index("Codon", drop=False)

    genome_names = {"whole genome", "genome", "all genes", "whole genome "}
    wanted = [str(c) for c in (include_columns or [])]
    panel_cols = []
    for cname in wanted:
        cname_s = str(cname).strip()
        if cname_s in df.columns and cname_s.lower() not in genome_names and cname_s not in panel_cols:
            panel_cols.append(cname_s)

    if not panel_cols:
        panel_cols = [
            c for c in df.columns
            if c != "Codon" and str(c).strip().lower() not in genome_names
        ]

    plot_df = df.reindex(codon_order)
    n_panels = len(panel_cols)
    if n_panels == 0:
        print("[WARN] No codon-usage columns selected; skipping codon-usage panel plot.")
        return None

    nrows = min(max(1, int(max_nrows)), int(np.ceil(np.sqrt(n_panels))))
    ncols = int(np.ceil(n_panels / nrows))

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(float(panel_w_in) * ncols, float(panel_h_in) * 1.5 * nrows),
        squeeze=False
    )
    plt.subplots_adjust(
        top=0.945,
        wspace=max(0.02, float(SUBPLOT_WSPACE) * 0.5),
        hspace=max(0.18, float(SUBPLOT_HSPACE) * 0.65),
    )

    y = np.arange(len(codon_order))
    xlab_map = {
        "ACU": "Average absolute codon usage frequency",
        "RCU": "Average relative codon usage frequency",
        "Z": "Average RCU z-score",
    }
    title_map = {
        "ACU": "Average absolute codon usage per cluster",
        "RCU": "Average relative codon usage per cluster",
        "Z": "Average RCU z-scores per cluster",
    }
    output_format = str(output_format or "png").strip().lstrip(".").lower() or "png"
    file_map = {
        "ACU": f"Average absolute codon usage per cluster.{output_format}",
        "RCU": f"Average relative codon usage per cluster.{output_format}",
        "Z": f"Average RCU z-scores per cluster.{output_format}",
    }

    finite_vals = []
    for col in panel_cols:
        vals = pd.to_numeric(plot_df[col], errors="coerce").to_numpy(dtype=float)
        finite_vals.extend(vals[np.isfinite(vals)].tolist())
    finite_vals = np.asarray(finite_vals, dtype=float)
    if finite_vals.size == 0:
        xmin, xmax = -1.0, 1.0
    else:
        vmin = float(np.nanmin(finite_vals))
        vmax = float(np.nanmax(finite_vals))
        if vmin < 0 and vmax > 0:
            lim = max(abs(vmin), abs(vmax))
            pad = 0.05 * lim if lim > 0 else 0.1
            xmin, xmax = -(lim + pad), (lim + pad)
        else:
            pad = 0.05 * max(abs(vmin), abs(vmax), 1e-12)
            xmin, xmax = min(0.0, vmin - pad), max(0.0, vmax + pad)

    ytick_fs = 5 if len(codon_order) >= 60 else (6 if len(codon_order) >= 40 else 7)
    panel_title_fs = max(int(CAPTION_SIZE) * 2, 22)
    super_title_fs = max(panel_title_fs + 4, 28)

    for i, col in enumerate(panel_cols):
        r, c = divmod(i, ncols)
        ax = axes[r, c]
        vals = pd.to_numeric(plot_df[col], errors="coerce").to_numpy(dtype=float)
        colors = ["red" if (np.isfinite(v) and v < 0) else "blue" for v in vals]
        ax.barh(y, vals, color=colors)
        ax.set_yticks(y)
        if c == 0:
            ax.set_yticklabels(codon_order, fontsize=ytick_fs)
            ax.set_ylabel("Codon")
            ax.tick_params(axis="y", length=2)
        else:
            ax.set_yticklabels([])
            ax.set_ylabel("")
            ax.tick_params(axis="y", length=0)
        ax.invert_yaxis()
        xlo = xmin
        xhi = xmax
        try:
            if x_min not in (None, ""):
                xlo = float(x_min)
        except Exception:
            pass
        try:
            if x_max not in (None, ""):
                xhi = float(x_max)
        except Exception:
            pass
        ax.set_xlim(xlo, xhi)
        try:
            cur0, cur1 = ax.get_ylim()
            logical_min = min(cur0, cur1)
            logical_max = max(cur0, cur1)
            ylo = float(y_min) if y_min not in (None, "") else logical_min
            yhi = float(y_max) if y_max not in (None, "") else logical_max
            if cur0 > cur1:
                ax.set_ylim(yhi, ylo)
            else:
                ax.set_ylim(ylo, yhi)
        except Exception:
            pass
        ax.set_xlabel(xlab_map.get(mode, "Value"))
        ax.axvline(0.0, linewidth=1.0, color="black", alpha=0.7)
        ax.grid(True, axis="x", alpha=0.25)
        ax.set_title(str(col), fontsize=panel_title_fs, fontweight=CAPTION_WEIGHT, pad=max(2, int(CAPTION_PAD * 0.6)))

    for j in range(n_panels, nrows * ncols):
        r, c = divmod(j, ncols)
        axes[r, c].axis("off")

    fig.suptitle(
        title_map.get(mode, "Average codon usage per cluster"),
        fontsize=super_title_fs,
        fontweight="bold",
        y=0.985,
    )
    chosen_name = output_basename.strip() if str(output_basename or "").strip() else os.path.splitext(file_map.get(mode, f"Average codon usage per cluster.{output_format}"))[0]
    figures_dir = os.path.join(os.path.dirname(codon_usage_workbook), "Figures")
    os.makedirs(figures_dir, exist_ok=True)
    out_png = os.path.join(figures_dir, f"{chosen_name}.{output_format}")
    fig.savefig(out_png, dpi=int(dpi), bbox_inches="tight")
    if show:
        try:
            plt.show()
        except Exception:
            try:
                plt.show(block=False)
            except Exception:
                pass
    else:
        plt.close(fig)

    print(f"[INFO] Codon-usage multi-panel plot saved to:\n  {out_png}")
    return out_png




def main():
    # ------------------------------------------------------------------
    # Auto-sync axis labels with the dimensionality-reduction method used
    # upstream (Clustering_Pipeline sets CODONPIPE_DIMRED_METHOD).
    # ------------------------------------------------------------------
    global X_LABEL, Y_LABEL

    dimred = os.environ.get("CODONPIPE_DIMRED_METHOD", "").strip().lower()
    x_from_env = os.environ.get("CODONPIPE_X_LABEL", "").strip()
    y_from_env = os.environ.get("CODONPIPE_Y_LABEL", "").strip()

    # Allow explicit axis labels via env vars (highest priority)
    if x_from_env:
        X_LABEL = x_from_env
    if y_from_env:
        Y_LABEL = y_from_env

    # Otherwise, auto-label based on dimred method, but only if the labels
    # are still the default t-SNE ones (so we don't clobber custom labels).
    if dimred and (not x_from_env) and (not y_from_env):
        default_x = str(X_LABEL).strip().lower() in {"t-sne 1", "tsne 1", "tsne1", "tsne_1"}
        default_y = str(Y_LABEL).strip().lower() in {"t-sne 2", "tsne 2", "tsne2", "tsne_2"}
        if default_x and default_y:
            if dimred == "tsne":
                X_LABEL, Y_LABEL = "tSNE1", "tSNE2"
            elif dimred == "umap":
                X_LABEL, Y_LABEL = "UMAP1", "UMAP2"
            elif dimred == "pca":
                X_LABEL, Y_LABEL = "PC1", "PC2"

    # Collect all UPPERCASE globals as config for the core plotter
    cfg = {k: v for k, v in globals().items() if k.isupper()}

    if cfg.get("DEBUG_PRINT_CONFIG", False):
        print("\n========== [PLOT DEBUG] Plotting_Pipeline.py ==========")
        print(f"[PLOT DEBUG] Script path: {os.path.abspath(__file__)}")
        print(f"[PLOT DEBUG] EXCEL_PATH={cfg.get('EXCEL_PATH')}")
        print(f"[PLOT DEBUG] UMAP_SHEET={cfg.get('UMAP_SHEET')} | DATA_SHEET={cfg.get('DATA_SHEET')}")
        print(f"[PLOT DEBUG] COLOR_MODE={cfg.get('COLOR_MODE')}")
        print(f"[PLOT DEBUG] DENSITY_CMAP_NAME={cfg.get('DENSITY_CMAP_NAME')}")
        print(
            "[PLOT DEBUG] ENRICH settings: "
            f"SCALE={cfg.get('ENRICHMENT_SCALE')}, "
            f"STYLE={cfg.get('ENRICHMENT_STYLE')}, "
            f"CMAP={cfg.get('ENRICHMENT_CMAP_NAME')}, "
            f"VMAX={cfg.get('ENRICHMENT_VMAX')}, "
            f"PCTL={cfg.get('ENRICHMENT_PERCENTILE')}, "
            f"SYM={cfg.get('ENRICHMENT_SYMMETRIC')}, "
            f"EPS={cfg.get('ENRICHMENT_EPS')}"
        )
        print(f"[PLOT DEBUG] ALL_GENES_PANEL_MODE={cfg.get('ALL_GENES_PANEL_MODE')}")
        print("=======================================================\n")

    # Run plotting
    if RUN_2D_DENSITY_PLOTS:
        run_density_dotplot(cfg)
    else:
        print("[INFO] 2D embedded density plots disabled; skipping.")

    # Optional: gene-cluster localization heatmap along genome axis
    if RUN_GENE_CLUSTER_LOCALIZATION_HEATMAP and cfg.get("EXCEL_PATH"):
        try:
            wb = cfg.get("EXCEL_PATH")
            figures_dir = os.path.join(os.path.dirname(wb), "Figures")
            os.makedirs(figures_dir, exist_ok=True)
            if str(GCHM_OUTPUT_FILENAME or "").strip():
                out_name = str(GCHM_OUTPUT_FILENAME).strip()
                out_png = out_name if os.path.isabs(out_name) else os.path.join(figures_dir, out_name)
            else:
                out_png = os.path.join(figures_dir, f"{GCHM_OUTPUT_BASENAME}.{FIGURE_FORMAT}")
            print("[INFO] Generating gene-cluster localization heatmap (KS)...")
            run_gene_cluster_localization_heatmap_from_workbook(
                workbook_path=wb,
                sheet_name=LOCUS_TAGS_SHEET,
                output_file=out_png,
                colormap=GCHM_COLORMAP,
                custom_cmaps_xlsx=GCHM_CUSTOM_CMAPS_XLSX,
                sigma=GCHM_SIGMA,
                spread_factor=GCHM_SPREAD_FACTOR,
                height_per_cluster=GCHM_HEIGHT_PER_CLUSTER,
                label_fontsize=GCHM_LABEL_FONTSIZE,
                dpi=GCHM_DPI,
                cmap_min_rel=GCHM_CMAP_MIN_REL,
                cmap_max_rel=GCHM_CMAP_MAX_REL,
                show=GCHM_SHOW_FIG,
                include_clusters=GCHM_INCLUDE_COLUMNS,
                x_min=GCHM_XMIN,
                x_max=GCHM_XMAX,
                y_min=GCHM_YMIN,
                y_max=GCHM_YMAX,
            )
            print(f"[INFO] Gene-cluster heatmap saved to:\n  {out_png}")
        except Exception as e:
            print(f"[WARN] Gene-cluster heatmap failed (skipping): {e}")


    # Optional: codon-usage multi-panel bar plots
    try:
        _plot_cluster_codon_usage_panels(
            main_workbook_path=cfg.get("EXCEL_PATH"),
            codon_usage_workbook=CODON_USAGE_WORKBOOK,
            plot_mode=CODON_USAGE_PLOT_MODE,
            include_columns=cfg.get("HIGHLIGHT_COLUMNS") or cfg.get("INCLUDE_COLUMNS"),
            max_nrows=int(cfg.get("MAX_NROWS", 2)),
            panel_w_in=float(CODON_USAGE_PANEL_W_IN),
            panel_h_in=float(CODON_USAGE_PANEL_H_IN),
            dpi=int(CODON_USAGE_PNG_DPI),
            show=bool(CODON_USAGE_SHOW_FIG),
            output_format=str(CODON_USAGE_OUTPUT_FORMAT or FIGURE_FORMAT),
            output_basename=str(CODON_USAGE_OUTPUT_BASENAME or ""),
            x_min=CODON_USAGE_XMIN,
            x_max=CODON_USAGE_XMAX,
            y_min=CODON_USAGE_YMIN,
            y_max=CODON_USAGE_YMAX,
        )
    except Exception as e:
        print(f"[WARN] Codon-usage multi-panel plot failed (skipping): {e}")

    # Methods/TXT sidecar exports are intentionally disabled in the GUI workflow.


if __name__ == "__main__":
    main()
