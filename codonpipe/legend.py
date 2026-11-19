#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: codonpipe/legend.py
"""
codonpipe.legend
================

Generate three manuscript-supporting TXT files:

1) Methods details (reproducibility-first)
   - exhaustive parameter listing + I/O paths + software versions

2) Legends (reader-facing)
   - describes each output figure/table so readers understand what they are seeing

3) Methods (paper-ready narrative)
   - written in the style of a Materials & Methods section, includes key parameters
     but avoids dumping every single knob unless relevant

Designed to be used by:
  * Clustering_Pipeline.py: writes all 3 files (overwrite each run)
  * Plotting_Pipeline.py: appends plotting-related sections to Legends + Methods details + Methods
"""

from __future__ import annotations

import os
import sys
import platform
from datetime import datetime
from typing import Any, Dict, List, Optional

try:
    from importlib import metadata as _metadata  # py>=3.8
except Exception:  # pragma: no cover
    _metadata = None


_RULE = "=" * 92
_RULE2 = "-" * 92


def _now_iso() -> str:
    try:
        return datetime.now().isoformat(timespec="seconds")
    except Exception:
        return str(datetime.now())


def _s(x: Any) -> str:
    try:
        return str(x)
    except Exception:
        return repr(x)


def _sec(title: str) -> str:
    return f"{_RULE}\n{title}\n{_RULE2}\n"


def _kv_lines(d: Dict[str, Any], keys: List[str], indent: str = "  ", redact_paths: bool = False) -> str:
    lines = []
    path_like_keys = {"default_root", "custom_cmap_excel", "plotting_pipeline_script_path", "trna_decoding_table_path"}
    for k in keys:
        if k in d:
            value = d.get(k)
            if redact_paths and k in path_like_keys:
                value = _format_path(value, True)
            lines.append(f"{indent}{k} = {_s(value)}")
    if not lines:
        return indent + "(none)"
    return "\n".join(lines)


def _write(path: str, text: str, mode: str = "w") -> None:
    folder = os.path.dirname(path)
    if folder:
        os.makedirs(folder, exist_ok=True)
    with open(path, mode, encoding="utf-8") as fh:
        fh.write(text)



def _format_path(p: str, redact: bool) -> str:
    """Format a file path for human-readable outputs.

    If redact=True, only the basename is shown to avoid embedding user-specific
    directory structures in text outputs meant for sharing.
    """
    if not p:
        return ""
    try:
        s = str(p)
    except Exception:
        return ""
    if not redact:
        return s
    try:
        return os.path.basename(s)
    except Exception:
        return s


def _format_config_value(key: str, value: Any, redact_paths: bool) -> str:
    """Format configuration values, redacting path-like keys when requested."""
    if not redact_paths:
        return _s(value)
    key_l = str(key or "").lower()
    path_markers = ("path", "file", "folder", "dir", "root", "xlsx", "excel", "fasta")
    if any(marker in key_l for marker in path_markers):
        return _format_path(value, True)
    return _s(value)

def _pkg_version(name: str) -> str:
    if _metadata is None:
        return "unknown"
    try:
        return _metadata.version(name)
    except Exception:
        return "unknown"


def _dimred_summary(SET: Dict[str, Any]) -> str:
    dim = str(SET.get("dimred_method", "umap")).lower().strip()
    if dim == "umap":
        return (
            "UMAP (2D) using "
            f"n_neighbors={SET.get('umap_neighbors')}, "
            f"min_dist={SET.get('umap_min_dist')}, "
            f"metric={SET.get('umap_metric')}, "
            f"components={SET.get('umap_components', 2)}, "
            f"randomize={SET.get('umap_randomize')}, "
            f"seed={SET.get('random_seed', 'NA')}."
        )
    if dim == "tsne":
        return (
            "t-SNE (2D) using "
            f"perplexity={SET.get('tsne_perplexity')}, "
            f"distance={SET.get('tsne_distance')}, "
            f"exaggeration={SET.get('tsne_exaggeration')}, "
            f"learning_rate={SET.get('tsne_learnrate')}, "
            f"dims={SET.get('tsne_dims', 2)}, "
            f"seed={SET.get('random_seed', 'NA')}."
        )
    return (
        "PCA using "
        f"n_components={SET.get('pca_npcs', 2)}, "
        f"center={SET.get('pca_center')}, "
        f"scale={SET.get('pca_scale')}, "
        f"seed={SET.get('random_seed', 'NA')}."
    )


def _clustering_summary(SET: Dict[str, Any]) -> str:
    cm = str(SET.get("cluster_method", "kmedoids")).lower().strip()
    if cm == "kmedoids":
        return f"K-medoids clustering (k={SET.get('kmedoids_k')}, distance={SET.get('kmedoids_dist')})."
    if cm == "hierarchical":
        return f"Hierarchical clustering (gene distance={SET.get('gene_dist_metric')}, linkage={SET.get('gene_linkage')})."
    if cm == "spectral":
        return f"Spectral clustering (k={SET.get('spectral_k')}, distance={SET.get('spectral_dist')})."
    if cm == "dbscan":
        return (
            f"DBSCAN clustering (eps={SET.get('dbscan_eps')}, minpts={SET.get('dbscan_minpts')}, distance={SET.get('dbscan_dist')})."
        )
    return f"Clustering method={SET.get('cluster_method')}."


def _usage_summary(SET: Dict[str, Any]) -> str:
    codon_range = str(SET.get('fasta_codon_range', 'all') or 'all').strip() or 'all'
    return f"Usage basis={SET.get('usage_basis')} with codon set={SET.get('codon_set')}; CDS codon range={codon_range}."


def write_all_texts(
    out_base: str,
    SET: Dict[str, Any],
    PIPELINE: Dict[str, Any],
    KS_SETTINGS: Dict[str, Any],
    clustering_results: Dict[str, Any],
    codon_file: str,
    cluster_file_path: str,
    fasta_path: str,
    out_xlsx: str,
    out_cluster_xlsx: str,
    out_codon_usage_xlsx: str = "",
    redact_paths: bool = True,
    sheet_name_quantitative: str = "Quantitative",
    sheet_name_quant_locustags: str = "Quantitative Locus tags",
    sheet_ks_raw: str = "2D KS - raw data",
    sheet_ks_cmp: str = "2D KS - comparison",
) -> Dict[str, str]:
    """
    Write three files:
      - <out_base>_Methods details.txt
      - <out_base>_Legends.txt
      - <out_base>_Methods.txt

    Returns a dict of paths.
    """
    paths = dict(
        methods_details=f"{out_base}_Methods details.txt",
        legends=f"{out_base}_Legends.txt",
        methods=f"{out_base}_Methods.txt",
    )

    _write(paths["methods_details"], _build_methods_details(
        SET, PIPELINE, KS_SETTINGS, clustering_results,
        codon_file=codon_file,
        cluster_file_path=cluster_file_path,
        fasta_path=fasta_path,
        out_xlsx=out_xlsx,
        out_cluster_xlsx=out_cluster_xlsx,
        out_codon_usage_xlsx=out_codon_usage_xlsx,
        redact_paths=redact_paths,
        sheet_name_quantitative=sheet_name_quantitative,
        sheet_name_quant_locustags=sheet_name_quant_locustags,
        sheet_ks_raw=sheet_ks_raw,
        sheet_ks_cmp=sheet_ks_cmp,
    ), mode="w")

    _write(paths["legends"], _build_legends(
        SET, PIPELINE, clustering_results,
        out_xlsx=out_xlsx,
        out_cluster_xlsx=out_cluster_xlsx,
        out_codon_usage_xlsx=out_codon_usage_xlsx,
        redact_paths=redact_paths,
        sheet_name_quantitative=sheet_name_quantitative,
        sheet_name_quant_locustags=sheet_name_quant_locustags,
        sheet_ks_raw=sheet_ks_raw,
        sheet_ks_cmp=sheet_ks_cmp,
    ), mode="w")

    _write(paths["methods"], _build_methods(
        SET, PIPELINE, KS_SETTINGS, clustering_results,
        codon_file=codon_file,
        cluster_file_path=cluster_file_path,
        fasta_path=fasta_path,
        redact_paths=redact_paths,
        sheet_name_quantitative=sheet_name_quantitative,
        sheet_name_quant_locustags=sheet_name_quant_locustags,
        sheet_ks_raw=sheet_ks_raw,
        sheet_ks_cmp=sheet_ks_cmp,
    ), mode="w")

    return paths


def _build_methods_details(
    SET: Dict[str, Any],
    PIPELINE: Dict[str, Any],
    KS_SETTINGS: Dict[str, Any],
    clustering_results: Dict[str, Any],
    codon_file: str,
    cluster_file_path: str,
    fasta_path: str,
    out_xlsx: str,
    out_cluster_xlsx: str,
    out_codon_usage_xlsx: str,
    redact_paths: bool,
    sheet_name_quantitative: str,
    sheet_name_quant_locustags: str,
    sheet_ks_raw: str,
    sheet_ks_cmp: str,
) -> str:
    gene_ids_file = clustering_results.get("gene_file", "") or clustering_results.get("gene_ids_file", "")
    output_subfolder = clustering_results.get("output_subfolder", "")

    txt = []
    txt.append(_sec("METHODS DETAILS (REPRODUCIBILITY-FIRST) — RUN METADATA"))
    txt.append(f"Generated: {_now_iso()}\n")
    txt.append(f"Python: {sys.version.split()[0]}\n")
    txt.append(f"Platform: {platform.platform()}\n")
    txt.append(f"Working directory: {_format_path(os.getcwd(), redact_paths) if redact_paths else os.getcwd()}\n\n")

    txt.append(_sec("INPUTS (FILE PATHS)"))
    txt.append(f"Codon usage workbook (counts): {_format_path(codon_file, redact_paths)}\n")
    if gene_ids_file:
        txt.append(f"GeneIDs workbook (annotations): {_format_path(gene_ids_file, redact_paths)}\n")
    txt.append(f"Cluster file (locus-tag lists): {_format_path(cluster_file_path, redact_paths)}\n")
    txt.append(f"CDS FASTA (canonicalization + metrics): {_format_path(fasta_path, redact_paths)}\n\n")

    txt.append(_sec("OUTPUTS (FILE PATHS)"))
    txt.append(f"Main workbook: {_format_path(out_xlsx, redact_paths)}\n")
    txt.append(f"Per-cluster workbook: {_format_path(out_cluster_xlsx, redact_paths)}\n")
    if out_codon_usage_xlsx:
        txt.append(f"Codon-usage workbook: {_format_path(out_codon_usage_xlsx, redact_paths)}\n")
    txt.append(f"Output subfolder (figures): {_format_path(output_subfolder, redact_paths)}\n\n")

    txt.append(_sec("CORE PIPELINE PARAMETERS — SET"))
    # Dump SET in a stable-ish order: show important blocks first, then everything.
    txt.append("Key parameters:\n")
    txt.append(_kv_lines(SET, [
        "default_root",
        "usage_basis", "codon_set", "fasta_codon_range",
        "dimred_method", "cluster_method", "do_2d_ks",
        "center_features", "scale_features",
        "umap_neighbors", "umap_min_dist", "umap_metric", "umap_components", "umap_randomize", "umap_clip_abs",
        "tsne_perplexity", "tsne_distance", "tsne_dims", "tsne_exaggeration", "tsne_learnrate",
        "pca_npcs", "pca_center", "pca_scale",
        "kmedoids_k", "kmedoids_dist",
        "spectral_k", "spectral_dist",
        "dbscan_eps", "dbscan_minpts", "dbscan_dist",
        "gene_dist_metric", "gene_linkage",
        "feature_dist_metric", "feature_linkage",
        "apply_smoothing", "smooth_window_genes",
        "apply_binning", "bin_size_genes",
        "heatmap_colormap_name", "heatmap_caxis_limits", "heatmap_fig_size", "xtick_every_genes",
        "scatter_fig_size", "scatter_point_size", "scatter_point_alpha", "scatter_edge_width", "scatter_color_mode",
        "density_nbins", "density_sigma", "density_use_log", "density_min_rel", "density_max_rel", "density_cmap_name",
        "density_enrichment_eps", "density_enrichment_use_log",
        "font_name", "font_size_axes", "font_size_xticks", "font_size_yticks", "font_size_titles", "colorbar_title_size", "figure_dpi",
        "heatmap_title_template", "scatter_title_template", "colorbar_title_string",
        "save_png", "save_pdf", "save_jpeg",
        "use_custom_colormaps", "custom_cmap_excel", "custom_cmap_sheets",
        "random_seed",
    ], redact_paths=redact_paths))
    txt.append("\n\n")

    txt.append(_sec("PIPELINE SETTINGS — PIPELINE"))
    for k, v in PIPELINE.items():
        txt.append(f"  {k} = {_s(v)}\n")
    txt.append("\n")

    txt.append(_sec("2D KS SETTINGS — KS_SETTINGS"))
    for k, v in KS_SETTINGS.items():
        txt.append(f"  {k} = {_s(v)}\n")
    txt.append("\n")

    txt.append(_sec("EXCEL OUTPUT STRUCTURE (SHEETS)"))
    txt.append(
        "Main workbook sheets (high level):\n"
        f"  - {PIPELINE.get('sheet_coordinates', 'Genes reordered')} (headered reordered-gene coordinate table)\n"
        f"  - {PIPELINE.get('sheet_locus_tags', 'Locus Tags')}\n"
        f"  - {PIPELINE.get('sheet_binary', 'Binary')}\n"
        f"  - {sheet_name_quantitative}\n"
        f"  - {sheet_name_quant_locustags}\n"
        "  - Meta\n"
        f"  - {sheet_ks_raw} / {sheet_ks_cmp} (if enabled)\n"
        "  - Missing (if applicable)\n\n"
    )
    if out_codon_usage_xlsx:
        txt.append(
            "Separate codon-usage workbook sheets (high level):\n"
            "  - ACU per cluster\n"
            "  - RCU\n"
            "  - ZCU\n\n"
            "A companion subfolder named 'Raw codon usage tables' is created next to this workbook. "
            "It contains one workbook per cluster (plus Whole genome), each with the sheets: counts, ACU, RCU, and z-scores.\n\n"
        )

    txt.append(_sec("SOFTWARE VERSIONS"))
    txt.append(
        "Versions captured at runtime (if available):\n"
        f"  numpy = {_pkg_version('numpy')}\n"
        f"  pandas = {_pkg_version('pandas')}\n"
        f"  matplotlib = {_pkg_version('matplotlib')}\n"
        f"  scipy = {_pkg_version('scipy')}\n"
        f"  scikit-learn = {_pkg_version('scikit-learn')}\n"
        f"  umap-learn = {_pkg_version('umap-learn')}\n"
        f"  xlsxwriter = {_pkg_version('XlsxWriter')}\n\n"
        "Note: optional dependencies (e.g., sklearn-extra for KMedoids) may be absent; the pipeline may fall back to alternatives.\n"
    )

    return "".join(txt)


def _build_legends(
    SET: Dict[str, Any],
    PIPELINE: Dict[str, Any],
    clustering_results: Dict[str, Any],
    out_xlsx: str,
    out_cluster_xlsx: str,
    out_codon_usage_xlsx: str,
    redact_paths: bool,
    sheet_name_quantitative: str,
    sheet_name_quant_locustags: str,
    sheet_ks_raw: str,
    sheet_ks_cmp: str,
) -> str:
    """
    Reader-facing legends: explain what each output contains and how to interpret it.
    Avoid full parameter dumps; keep it understandable to a reader.
    """
    dimred_line = _dimred_summary(SET)
    cluster_line = _clustering_summary(SET)
    usage_line = _usage_summary(SET)

    txt = []
    txt.append(_sec("LEGENDS (READER-FACING) — OVERVIEW"))
    txt.append(
        "This file provides reader-facing descriptions of each output figure and table produced by the codon-usage "
        "clustering pipeline. These descriptions are intended to be used as figure legends and table captions.\n\n"
    )

    txt.append(_sec("OUTPUT — MAIN EXCEL WORKBOOK"))
    txt.append(f"File: {_format_path(out_xlsx, redact_paths)}\n\n")
    txt.append(
        "This workbook consolidates the reordered-gene coordinate table, cluster membership tables, per-gene FASTA-derived "
        "metrics, DAVID-derived summaries when requested, and (optionally) 2D KS comparisons between clusters.\n\n"
    )

    if out_codon_usage_xlsx:
        txt.append(_sec("OUTPUT — CODON-USAGE WORKBOOK"))
        txt.append(f"File: {_format_path(out_codon_usage_xlsx, redact_paths)}\n\n")
        txt.append(
            "This workbook centralizes cluster-level codon-usage summaries. It contains the sheets 'ACU', 'RCU', and 'ZCU', "
            "and, when decoding strategy analyses are enabled, also 'ATU', 'RTU', and 'ZTU'. "
            "A companion subfolder named 'Raw codon usage tables' is created next to this workbook; it contains "
            "one workbook per cluster (plus Whole genome), each with the sheets 'Codon counts', 'ACU', 'RCU', 'ZCU', and optionally 'tRNA counts', 'ATU', 'RTU', and 'ZTU'. "
            "The decoding table can use either the original multi-column format with anticodon weights or a simplified 3-column format (AA, codon, anticodon). Pooled tRNA labels such as UCU/CCU are accepted. "
            "If the unified decoding table contains tRNA molecules/cell mean values in column 4, an additional workbook named 'tRNA abundance correlations.xlsx' can also be exported, together with cluster-by-amino-acid heatmaps for ATU, RTU, and ZTU abundance correlations.\n\n"
            f"Context for this run: {usage_line}\n\n"
        )

    txt.append(_sec(f"TABLE — '{PIPELINE.get('sheet_coordinates', 'Genes reordered')}' (reordered genes + 2D coordinates)"))
    txt.append(
        "This sheet provides the reordered genes together with gene names, 2D embedding coordinates, and protein descriptions. "
        "These coordinates are used by the plotting pipeline to generate cluster overlay figures.\n\n"
        f"Embedding used: {dimred_line}\n\n"
    )

    txt.append(_sec("TABLE — 'Codons reordered'"))
    txt.append(
        "This sheet lists the feature ordering (codons or amino acids) used for visualization in the heatmap. It provides a direct "
        "mapping between heatmap column positions and feature identifiers.\n\n"
    )

    txt.append(_sec(f"TABLE — '{PIPELINE.get('sheet_locus_tags', 'Locus Tags')}'"))
    txt.append(
        "This sheet merges (i) the genome-ordered list of locus tags produced by the clustering pipeline with (ii) cluster membership "
        "information from the user-provided cluster file. Each row corresponds to a gene in genome order and records which curated cluster(s) "
        "it belongs to.\n\n"
        "Interpretation: genes sharing a curated cluster label are expected to represent related functional groups, genomic islands, or "
        "user-defined categories.\n\n"
    )

    txt.append(_sec(f"TABLE — '{PIPELINE.get('sheet_binary', 'Binary')}' (binary cluster membership table)"))
    txt.append(
        "This sheet is a machine-friendly representation of curated cluster membership. Rows correspond to genome-ordered genes and columns "
        "correspond to curated clusters (0/1 entries). This table is the primary input for downstream plotting of per-cluster overlays.\n\n"
    )

    txt.append(_sec(f"TABLE — '{sheet_name_quantitative}' (per-gene FASTA metrics)"))
    txt.append(
        "This sheet contains per-gene sequence-derived metrics computed from the CDS FASTA file (e.g., sequence length, GC fraction, strand, "
        "and additional threshold-based flags). These metrics enable stratifying or interpreting clusters in relation to sequence composition "
        "and other derived features.\n\n"
    )

    txt.append(_sec(f"TABLE — '{sheet_name_quant_locustags}' (locus-tag lists derived from binary metrics)"))
    txt.append(
        "This sheet provides a convenience view for binary (0/1) metrics present in the 'Quantitative' sheet: for each binary metric column, "
        "it lists the locus tags of genes for which the value equals 1. This is useful for quickly extracting gene sets meeting specific criteria "
        "(e.g., low GC flags or membership in specific codon-bias categories).\n\n"
    )

    txt.append(_sec("TABLE — 'Meta'"))
    txt.append(
        "This sheet records run metadata, input file paths, and key configuration values. It provides provenance for the outputs "
        "and supports traceability.\n\n"
    )

    txt.append(_sec(f"TABLES — '{sheet_ks_raw}' and '{sheet_ks_cmp}' (2D KS comparisons)"))
    if bool(SET.get("do_2d_ks", True)):
        txt.append(
            "These sheets summarize pairwise differences between clusters in the 2D embedding using a 2D Kolmogorov–Smirnov (KS) statistic.\n\n"
            f"- '{sheet_ks_raw}' contains a long-form table of all cluster-pair comparisons, including BH-adjusted p-values (p_adj_BH).\n"
            f"- '{sheet_ks_cmp}' contains two matrices in the original cluster order:\n"
            "    (i) the p_adj_BH matrix (diagonal = 1), and\n"
            "    (ii) the -log10(p_adj_BH) matrix (diagonal = 0), which provides a more readable scale for small p-values.\n\n"
            "Interpretation: smaller p_adj_BH values (or larger -log10 values) indicate stronger evidence that two clusters occupy different "
            "regions in the 2D embedding.\n\n"
        )
    else:
        txt.append("2D KS analysis was disabled for this run.\n\n")

    txt.append(_sec("OUTPUT — PER-CLUSTER GENE LISTS WORKBOOK"))
    txt.append(f"File: {_format_path(out_cluster_xlsx, redact_paths)}\n\n")
    txt.append(
        "This workbook provides one sheet per curated cluster. Each sheet lists member genes along with their 2D embedding coordinates and "
        "available annotations (e.g., gene symbols and descriptions when present). This output supports manual review of cluster composition "
        "and targeted follow-up analyses.\n\n"
    )

    txt.append(_sec("FIGURE — CLUSTER OVERLAY PLOTS (Plotting pipeline)"))
    txt.append(
        "The plotting pipeline generates multi-panel cluster overlay figures in the 2D embedding. Each panel shows one curated cluster, "
        "where genes are plotted at their embedding coordinates and colored by either local density or enrichment relative to all genes "
        "(depending on the chosen coloring mode).\n\n"
        "Interpretation:\n"
        "- Density mode highlights where cluster genes concentrate in the embedding.\n"
        "- Enrichment mode highlights embedding regions where the cluster is over-represented relative to the global gene distribution.\n\n"
    )

    txt.append(_sec("RUN CONTEXT (ONE-LINE SUMMARY)"))
    txt.append(f"{usage_line} {dimred_line} {cluster_line}\n")

    return "".join(txt)


def _build_methods(
    SET: Dict[str, Any],
    PIPELINE: Dict[str, Any],
    KS_SETTINGS: Dict[str, Any],
    clustering_results: Dict[str, Any],
    codon_file: str,
    cluster_file_path: str,
    fasta_path: str,
    redact_paths: bool,
    sheet_name_quantitative: str,
    sheet_name_quant_locustags: str,
    sheet_ks_raw: str,
    sheet_ks_cmp: str,
) -> str:
    """
    Manuscript-style Methods: narrative, paper-ready.
    Includes key parameters for reproducibility but not a full dump.
    """
    dimred_line = _dimred_summary(SET)
    cluster_line = _clustering_summary(SET)
    usage_line = _usage_summary(SET)

    txt = []
    txt.append(_sec("METHODS (MANUSCRIPT-STYLE)"))

    txt.append(
        "Codon/AA usage preprocessing\n"
        "----------------------------\n"
        "Codon usage matrices were obtained from a codon-count Excel workbook and processed according to the selected usage basis and codon set. "
        "Features were optionally centered and scaled across genes prior to downstream analyses.\n\n"
        f"Input codon usage workbook: {codon_file}\n"
        f"Usage configuration: {usage_line}\n"
        f"Feature normalization: center_features={SET.get('center_features')}, scale_features={SET.get('scale_features')}.\n\n"
    )

    txt.append(
        "Dimensionality reduction\n"
        "------------------------\n"
        "To visualize genome-wide codon usage variation, genes were embedded into a two-dimensional space using the specified dimensionality "
        "reduction method. The resulting 2D coordinates were exported to the main output workbook and used for downstream visualization.\n\n"
        f"{dimred_line}\n\n"
    )

    txt.append(
        "Clustering and genome ordering\n"
        "------------------------------\n"
        "Genes were clustered based on their codon/AA usage profiles and/or their 2D embedding. The clustering output was used to derive a gene "
        "ordering for heatmap visualization, enabling detection of regional genomic structure.\n\n"
        f"{cluster_line}\n\n"
    )

    txt.append(
        "Heatmap visualization\n"
        "---------------------\n"
        "A heatmap of codon/AA usage was generated with genes ordered by the clustering-derived ordering and features reordered for visualization. "
        "Optional smoothing and/or binning along the gene axis were applied to emphasize broader regional patterns.\n\n"
        f"apply_smoothing={SET.get('apply_smoothing')}, smooth_window_genes={SET.get('smooth_window_genes')}; "
        f"apply_binning={SET.get('apply_binning')}, bin_size_genes={SET.get('bin_size_genes')}.\n"
        f"Colormap={SET.get('heatmap_colormap_name')}, caxis_limits={SET.get('heatmap_caxis_limits')}, figure_size={SET.get('heatmap_fig_size')}.\n\n"
    )

    txt.append(
        "Integration of curated gene sets and FASTA-derived metrics\n"
        "---------------------------------------------------------\n"
        "Curated gene sets (clusters) provided as locus-tag lists were integrated with the genome-ordered genes. A binary membership table "
        "was generated (0/1 membership per cluster) for downstream plotting. In parallel, a CDS FASTA file was used to canonicalize identifiers "
        "and compute per-gene sequence metrics (e.g., length, GC fraction, and derived flags), exported to the 'Quantitative' sheet. "
        "Binary metric columns were additionally summarized as locus-tag lists in the 'Quantitative Locus tags' sheet.\n\n"
        f"Cluster file: {cluster_file_path}\n"
        f"FASTA file: {fasta_path}\n"
        f"Sheets: '{PIPELINE.get('sheet_binary', 'Binary')}', '{sheet_name_quantitative}', '{sheet_name_quant_locustags}'.\n\n"
    )

    if bool(SET.get("do_2d_ks", True)):
        txt.append(
            "2D KS cluster comparisons\n"
            "-------------------------\n"
            "To quantify whether curated clusters occupy distinct regions in the 2D embedding, pairwise 2D Kolmogorov–Smirnov (KS) comparisons were "
            "performed between clusters using their 2D coordinate distributions. Statistical significance was assessed using permutation testing and "
            "corrected for multiple testing across all cluster pairs using the Benjamini–Hochberg procedure. Results were exported as a long-form table "
            f"('{sheet_ks_raw}') and as matrices of BH-adjusted p-values and -log10(p_adj_BH) values ('{sheet_ks_cmp}').\n\n"
            f"KS parameters: method={KS_SETTINGS.get('method')}, bins={KS_SETTINGS.get('bins')}, n_perm={KS_SETTINGS.get('n_perm')}, "
            f"alpha={KS_SETTINGS.get('alpha')}, random_seed={KS_SETTINGS.get('random_seed')}.\n\n"
        )
    else:
        txt.append(
            "2D KS cluster comparisons\n"
            "-------------------------\n"
            "Pairwise 2D KS comparisons were disabled for this run.\n\n"
        )

    txt.append(
        "Cluster overlay plotting\n"
        "------------------------\n"
        "Cluster overlay figures were generated from the main workbook by reading the 2D coordinates sheet and the binary cluster membership sheet. "
        "Each panel displays one curated cluster, with genes positioned by their 2D coordinates and colored by either local density or enrichment "
        "relative to all genes. Exact plotting parameters are recorded in the 'Methods details' file for full reproducibility.\n\n"
    )

    return "".join(txt)


def append_plotting_sections(
    out_base: str,
    cfg: Dict[str, Any],
    workbook_path: str,
    redact_paths: bool = True,
) -> None:
    """
    Append plotting-related info to:
      - <out_base>_Methods details.txt (parameter dump)
      - <out_base>_Legends.txt (reader-facing figure legend addendum)
      - <out_base>_Methods.txt (brief methods addendum)

    This is called from Plotting_Pipeline.py after plots are generated.
    """
    md_path = f"{out_base}_Methods details.txt"
    leg_path = f"{out_base}_Legends.txt"
    m_path = f"{out_base}_Methods.txt"

    # Derive expected output filename tag (best-effort)
    base = os.path.splitext(os.path.basename(workbook_path))[0]
    color_mode = str(cfg.get("COLOR_MODE", "density") or "density").lower().strip()
    enrich_scale = str(cfg.get("ENRICHMENT_SCALE", "ratio") or "ratio").lower().strip()
    enrich_style = str(cfg.get("ENRICHMENT_STYLE", "unilateral_diverging") or "unilateral_diverging").lower().strip()
    mode_tag = f"_ENRICH_{enrich_scale}_{enrich_style}" if color_mode == "enrichment" else "_DENSITY"
    out_png = os.path.join(os.path.dirname(workbook_path), f"{base}_UMAP_local{mode_tag}.png")

    # ---- Append to Methods details (full config dump)
    txt = []
    txt.append(_sec("PLOTTING PARAMETERS (APPENDED) — Plotting_Pipeline.py"))
    txt.append(f"Workbook used: {_format_path(workbook_path, redact_paths)}\n")
    txt.append(f"Expected output figure (best-effort): {_format_path(out_png, redact_paths)}\n\n")
    txt.append("Selected key plotting parameters:\n")
    txt.append(_kv_lines(cfg, [
        "UMAP_SHEET", "DATA_SHEET",
        "PLOT_MODE", "MAX_NROWS", "PANEL_W_IN", "PANEL_H_IN",
        "X_LABEL", "Y_LABEL",
        "COLOR_MODE",
        "DENSITY_CMAP_NAME",
        "ALL_SCATTER_DENSITY_NBINS", "ALL_SCATTER_DENSITY_SIGMA", "ALL_SCATTER_DENSITY_USE_LOG",
        "ENRICHMENT_SCALE", "ENRICHMENT_STYLE", "ENRICHMENT_CMAP_NAME",
        "ENRICHMENT_VMAX", "ENRICHMENT_PERCENTILE", "ENRICHMENT_SYMMETRIC", "ENRICHMENT_EPS",
        "SAVE_FIG", "PNG_DPI",
    ]))
    txt.append("\n\n")
    txt.append("Full configuration dump (all UPPERCASE globals captured):\n")
    for k in sorted([k for k in cfg.keys() if str(k).isupper()]):
        txt.append(f"  {k} = {_format_config_value(k, cfg.get(k), redact_paths)}\n")
    txt.append("\n")
    _write(md_path, "".join(txt), mode="a")

    # ---- Append to Legends (reader-facing)
    txt2 = []
    txt2.append(_sec("FIGURE LEGEND ADDENDUM — CLUSTER OVERLAY PLOTS (APPENDED)"))
    txt2.append(
        "Cluster overlay plots show genes in the 2D embedding, with one panel per curated cluster. Within each panel, genes belonging to the "
        "cluster are displayed at their embedding coordinates. Colors represent either local density of the cluster in embedding space or "
        "enrichment relative to all genes (depending on COLOR_MODE).\n\n"
    )
    txt2.append(f"Figure output (best-effort filename): {_format_path(out_png, redact_paths)}\n")
    txt2.append(
        f"Coloring mode: {color_mode}. "
        "In enrichment mode, enrichment is computed by comparing a smoothed local density field of cluster genes to the corresponding global density field.\n\n"
    )
    _write(leg_path, "".join(txt2), mode="a")

    # ---- Append to Methods (brief narrative)
    txt3 = []
    txt3.append(_sec("METHODS ADDENDUM — CLUSTER OVERLAY PLOTTING (APPENDED)"))
    txt3.append(
        "Cluster overlay figures were generated by plotting genes at their 2D embedding coordinates and highlighting curated cluster members. "
        "Local density fields were estimated on a grid and smoothed with a Gaussian kernel; enrichment was computed as the ratio (or log2-ratio) "
        "between cluster density and global density (where applicable). Key plotting settings are reported in the Methods details file.\n\n"
    )
    txt3.append(f"Output figure (best-effort filename): {_format_path(out_png, redact_paths)}\n\n")
    _write(m_path, "".join(txt3), mode="a")
