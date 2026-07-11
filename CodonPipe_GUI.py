#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CodonPipe GUI launcher.

Desktop tkinter interface for Clustering_Pipeline.py.
The GUI collects user settings and suppresses the console prompts from the
existing pipeline by monkey-patching a few interactive helper functions.

Assumptions for this GUI version:
- CodonPipe_GUI.py, Clustering_Pipeline.py and Plotting_Pipeline.py are in the
  same folder.
- An example folder may be present next to this script:
      Example (Salmonella SL1344)/
          SL1344 cds from genomic.fna
          SL1344 clusters.xlsx
"""

import copy
import io
import itertools
import json
import os
import queue
import shutil
import re
import sys
import threading
import traceback
import warnings
from contextlib import redirect_stdout, redirect_stderr

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import Clustering_Pipeline as CP
from codonpipe.fasta_metrics import (
    get_fasta_metric_group_defaults,
    build_fasta_metric_cluster_df,
    append_fasta_metric_clusters,
    build_locus_index,
)

try:
    from codonpipe.excel_outputs import read_trna_decoding_table as _read_trna_decoding_table_for_gui
except Exception:
    _read_trna_decoding_table_for_gui = None

CODONPIPE_GUI_BUILD = "direct-decoding-v46-prefilled-figure-details-2026-06-26"
GCHM_DEFAULT_SIGMA = 7.5
PLOT6_DEFAULT_DPI = 300
PLOT6_DEFAULT_FIG_WIDTH = 18.0
PLOT6_DEFAULT_FIG_HEIGHT = 8.0
PLOT6_DEFAULT_YMIN = -2.1
PLOT6_DEFAULT_YMAX = 3.5
PLOT6_DEFAULT_GROUP_BAR_Y = -0.26
PLOT6_DEFAULT_GROUP_LABEL_GAP = 0.07
PLOT6_DEFAULT_STAR_OFFSET = 0.07
# Defaults shown in the Figure details tab. These mirror the automatic values
# used by the plotting code whenever users leave the corresponding custom
# fields empty, so enabling customization starts from a known baseline.
TRNA_SUPP_DETAILS_DEFAULT_DPI = 300
TRNA_SUPP_DETAILS_DEFAULT_FIG_WIDTH = 18.0
TRNA_SUPP_DETAILS_DEFAULT_FIG_HEIGHT = 4.0
TRNA_SUPP_DETAILS_DEFAULT_CELL_HEIGHT = 0.07
TRNA_SUPP_DETAILS_DEFAULT_XTICK_EVERY = 500
TRNA_SUPP_DETAILS_DEFAULT_YTICK_FONTSIZE = 7
TRNA_SUPP_DETAILS_DEFAULT_TITLE_FONTSIZE = 10
TRNA_SHIFT_DETAILS_DEFAULT_DPI = 300
TRNA_SHIFT_DETAILS_DEFAULT_FIG_WIDTH = 8.2
TRNA_SHIFT_DETAILS_DEFAULT_FIG_HEIGHT = 5.3
TRNA_SHIFT_DETAILS_DEFAULT_CELL_WIDTH = 0.8
TRNA_SHIFT_DETAILS_DEFAULT_CELL_HEIGHT = 1.5
TRNA_SHIFT_DETAILS_DEFAULT_XTICK_FONTSIZE = 7
TRNA_SHIFT_DETAILS_DEFAULT_YTICK_FONTSIZE = 11
TRNA_SHIFT_DETAILS_DEFAULT_TITLE_FONTSIZE = 12
GCHM_DEFAULT_SPREAD_FACTOR = 4.0
CP.SET["gchm_sigma"] = GCHM_DEFAULT_SIGMA
CP.SET["gchm_spread_factor"] = GCHM_DEFAULT_SPREAD_FACTOR

PLOTTING_SCRIPT = os.path.join(HERE, "Plotting_Pipeline.py")
EXAMPLE_DIR = os.path.join(HERE, "Example (Salmonella SL1344)")
EXAMPLE_FASTA = os.path.join(EXAMPLE_DIR, "SL1344 cds from genomic.fna")
EXAMPLE_CLUSTERS = os.path.join(EXAMPLE_DIR, "SL1344 clusters.xlsx")

DEFAULT_BROWSER_ROOT = ""

USER_CLUSTER_CHOICES = [
    "Inferred from FASTA file",
    "Inferred from DAVID gene2terms txt file",
    "Provided by user, xlsx file with 1 column per cluster",
]
USER_CLUSTER_TO_INTERNAL = {
    "Inferred from FASTA file": "basic",
    "Inferred from DAVID gene2terms txt file": "basic",
    "Provided by user, xlsx file with 1 column per cluster": "refined",
    # Backward-compatible labels from older GUI builds
    "Not provided - naïve": "basic",
    "Not provided - naive": "basic",
    "Provided by user": "refined",
}
USER_CLUSTER_TO_TERMS_SOURCE = {
    "Inferred from FASTA file": "geneids",
    "Inferred from DAVID gene2terms txt file": "david_gene2terms",
    "Provided by user, xlsx file with 1 column per cluster": "geneids",
    "Provided by user": "geneids",
}
INTERNAL_CLUSTER_TO_USER = {
    "basic": "Inferred from FASTA file",
    "refined": "Provided by user, xlsx file with 1 column per cluster",
}

USAGE_METRIC_CHOICES = [
    "Relative codon usage",
    "Absolute codon usage",
    "Amino acid identity",
]
USAGE_DISPLAY_TO_INTERNAL = {
    "Relative codon usage": "RCU",
    "Absolute codon usage": "ACU",
    "Amino acid identity": "AA",
    "RCU": "RCU",
    "ACU": "ACU",
    "AA": "AA",
}
USAGE_INTERNAL_TO_DISPLAY = {
    "RCU": "Relative codon usage",
    "ACU": "Absolute codon usage",
    "AA": "Amino acid identity",
}

DIMRED_CHOICES = ["UMAP", "tSNE", "PCA"]
DIMRED_DISPLAY_TO_INTERNAL = {
    "UMAP": "umap",
    "umap": "umap",
    "tSNE": "tsne",
    "t-SNE": "tsne",
    "tsne": "tsne",
    "PCA": "pca",
    "pca": "pca",
}
DIMRED_INTERNAL_TO_DISPLAY = {"umap": "UMAP", "tsne": "tSNE", "pca": "PCA"}

CLUSTERING_METHOD_CHOICES = ["kmeans", "kmedoids", "hierarchical"]
STATISTICAL_TEST_CHOICES = ["disabled", "2D Kolmogorov-Smirnov"]
FUNCTIONAL_SCAN_CHOICES = ["disabled", "DAVID sliding-window scan"]

CODON_USAGE_MODE_CHOICES = ["ACU", "RCU", "ZCU"]
USER_CODON_MODE_TO_INTERNAL = {
    "ACU": "ACU",
    "RCU": "RCU",
    "ZCU": "Z",
    "RCU z-scores": "Z",
    "Z": "Z",
}
INTERNAL_CODON_MODE_TO_USER = {
    "ACU": "ACU",
    "RCU": "RCU",
    "Z": "ZCU",
    "NONE": "ZCU",
}

CODON_COMPARE_METRIC_CHOICES = [
    "Relative codon usage",
    "Absolute codon frequency",
    "Amino acid identity",
]

CODON_COMPARE_RCU_DISPLAY_CHOICES = [
    "Relative codon usage with genome",
    "Relative codon usage no genome",
    "z-scores",
]
CODON_COMPARE_STAT_CHOICES = [
    "None",
    "Student t test",
    "Welch t test",
    "Mann-Whitney U",
]
CODON_COMPARE_PLOT_STYLE_CHOICES = [
    "Mean ± SD",
    "Line plot",
    "Boxplot",
    "Violin plot",
]

FASTA_EXTENSIONS = (".fasta", ".fa", ".fna", ".ffn", ".fas")
PRELOADED_GENOME_DIRS = [
    os.path.join(HERE, "Preloaded genomes"),
    os.path.join(HERE, "preloaded_genomes"),
    os.path.join(HERE, "genomes"),
    EXAMPLE_DIR,
]
# Add future bundled genomes here. Paths can be absolute or relative to CodonPipe_GUI.py.
# The GUI also auto-discovers FASTA files placed in PRELOADED_GENOME_DIRS.
PRELOADED_GENOME_REGISTRY = [
    {
        "label": "Salmonella enterica serovar Typhimurium SL1344",
        "organism": "Salmonella SL1344",
        "fasta": EXAMPLE_FASTA,
    },
    {
        "label": "Shigella flexneri",
        "organism": "Shigella flexneri",
        "fasta": os.path.join(HERE, "Preloaded genomes", "Shigella flexneri cds from genomic.fna"),
    },
    {
        "label": "Escherichia coli K-12 MG1655",
        "organism": "E. coli MG1655",
        "fasta": os.path.join(HERE, "Preloaded genomes", "E coli MG1655 cds from genomic.fna"),
    },
    {
        "label": "Pseudomonas aeruginosa PAO1",
        "organism": "P. aeruginosa PAO1",
        "fasta": os.path.join(HERE, "Preloaded genomes", "Pseudomonas aeruginosa PAO1 cds from genomic.fna"),
    },
]


def _normalize_preloaded_name_for_match(value):
    """Normalize organism/file names for robust preloaded companion-file matching."""
    s = str(value or "").strip().lower()
    if not s:
        return ""
    s = os.path.splitext(os.path.basename(s))[0] if os.path.sep in s or "." in os.path.basename(s) else s
    s = re.sub(r"(?i)\b(cds|from|genomic|genome|sequence|sequences|fasta|ffn|fna)\b", " ", s)
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return " ".join(s.split())


def _preloaded_name_variants_from_record(rec):
    """Return likely species/strain name variants for companion files."""
    vals = []
    if rec:
        vals.extend([rec.get("organism", ""), rec.get("label", "")])
        fasta = rec.get("fasta", "")
        if fasta:
            vals.append(_clean_genome_label_from_filename(fasta))
            vals.append(os.path.splitext(os.path.basename(str(fasta)))[0])
    out, seen = [], set()
    for v in vals:
        raw = " ".join(str(v or "").replace("_", " ").replace("-", " ").split())
        norm = _normalize_preloaded_name_for_match(raw)
        for candidate in (raw, norm):
            candidate = str(candidate or "").strip()
            key = candidate.lower()
            if candidate and key not in seen:
                seen.add(key)
                out.append(candidate)
    return out


def _clean_genome_label_from_filename(path):
    base = os.path.splitext(os.path.basename(str(path)))[0]
    base = re.sub(r"(?i)\b(cds|from|genomic|genome|sequence|sequences|fasta|ffn|fna)\b", " ", base)
    base = re.sub(r"[_\-]+", " ", base)
    return " ".join(base.split()) or os.path.splitext(os.path.basename(str(path)))[0]


def _preferred_decoding_sheet_name(path):
    """Prefer the compact pooled decoding sheet when a workbook has multiple sheets."""
    wanted = {
        re.sub(r"[^a-z0-9]+", "", x.lower())
        for x in [
            "Decoding table (compact)", "Decoding table compact",
            "CA table", "Codon-anticodon table", "Codon anticodon table"
        ]
    }
    try:
        xls = pd.ExcelFile(path)
    except Exception:
        return ""
    for sheet in xls.sheet_names:
        key = re.sub(r"[^a-z0-9]+", "", str(sheet or "").strip().lower())
        if key in wanted:
            return sheet
    return xls.sheet_names[0] if xls.sheet_names else ""


def _discover_preloaded_genomes():
    """Return valid bundled genomes visible to the GUI.

    A genome is considered available only when its FASTA file exists. To add more
    genomes without editing code, place CDS FASTA files in a folder named
    'Preloaded genomes', 'preloaded_genomes', or 'genomes' next to this GUI.
    """
    records = []
    seen_paths = set()

    def _add(label, organism, fasta):
        fasta = os.path.abspath(os.path.expanduser(str(fasta or "")))
        if not fasta or not os.path.isfile(fasta):
            return
        key = os.path.normcase(fasta)
        if key in seen_paths:
            return
        seen_paths.add(key)
        records.append({
            "label": str(label or _clean_genome_label_from_filename(fasta)).strip(),
            "organism": str(organism or label or _clean_genome_label_from_filename(fasta)).strip(),
            "fasta": fasta,
        })

    for rec in PRELOADED_GENOME_REGISTRY:
        _add(rec.get("label"), rec.get("organism"), rec.get("fasta"))

    # Keep startup fast and predictable: only scan the top level of the
    # preloaded-genome folders. Recursive os.walk() can become very slow when
    # users accidentally place large folders in/near "Preloaded genomes".
    for folder in PRELOADED_GENOME_DIRS:
        folder = os.path.abspath(os.path.expanduser(str(folder)))
        if not os.path.isdir(folder):
            continue
        try:
            files = os.listdir(folder)
        except Exception:
            continue
        for fname in files:
            path = os.path.join(folder, fname)
            if os.path.isfile(path) and fname.lower().endswith(FASTA_EXTENSIONS):
                label = _clean_genome_label_from_filename(path)
                _add(label, label, path)

    records.sort(key=lambda r: r["label"].lower())
    return records


def _default_preloaded_genome_label(records):
    """Prefer Salmonella/SL1344 as the startup genome when it is available."""
    records = list(records or [])
    if not records:
        return "No preloaded genome found"

    def _score(rec):
        text = " ".join([
            str(rec.get("label", "")),
            str(rec.get("organism", "")),
            os.path.basename(str(rec.get("fasta", ""))),
        ]).lower()
        has_salmonella = "salmonella" in text
        has_sl1344 = "sl1344" in text
        if has_salmonella and has_sl1344:
            return 0
        if has_sl1344:
            return 1
        if has_salmonella:
            return 2
        return 10

    best = sorted(records, key=lambda rec: (_score(rec), str(rec.get("label", "")).lower()))[0]
    return str(best.get("label") or records[0].get("label") or "No preloaded genome found")


# -----------------------------------------------------------------------------
# Codon-usage comparison helpers
# -----------------------------------------------------------------------------

_CODON_TO_AA3 = {
    "TTT":"Phe", "TTC":"Phe", "TTA":"Leu", "TTG":"Leu",
    "CTT":"Leu", "CTC":"Leu", "CTA":"Leu", "CTG":"Leu",
    "ATT":"Ile", "ATC":"Ile", "ATA":"Ile", "ATG":"Met",
    "GTT":"Val", "GTC":"Val", "GTA":"Val", "GTG":"Val",
    "TCT":"Ser", "TCC":"Ser", "TCA":"Ser", "TCG":"Ser",
    "CCT":"Pro", "CCC":"Pro", "CCA":"Pro", "CCG":"Pro",
    "ACT":"Thr", "ACC":"Thr", "ACA":"Thr", "ACG":"Thr",
    "GCT":"Ala", "GCC":"Ala", "GCA":"Ala", "GCG":"Ala",
    "TAT":"Tyr", "TAC":"Tyr",
    "CAT":"His", "CAC":"His", "CAA":"Gln", "CAG":"Gln",
    "AAT":"Asn", "AAC":"Asn", "AAA":"Lys", "AAG":"Lys",
    "GAT":"Asp", "GAC":"Asp", "GAA":"Glu", "GAG":"Glu",
    "TGT":"Cys", "TGC":"Cys", "TGG":"Trp",
    "CGT":"Arg", "CGC":"Arg", "CGA":"Arg", "CGG":"Arg",
    "AGT":"Ser", "AGC":"Ser", "AGA":"Arg", "AGG":"Arg",
    "GGT":"Gly", "GGC":"Gly", "GGA":"Gly", "GGG":"Gly",
}

_CODON_ORDER_BY_AA = []
for _aa in sorted(set(_CODON_TO_AA3.values())):
    for _codon in sorted([c for c, a in _CODON_TO_AA3.items() if a == _aa]):
        _CODON_ORDER_BY_AA.append((_aa, _codon))
_CODON_ORDER = [c for _aa, c in _CODON_ORDER_BY_AA]
_AA_TO_CODONS = {}
for _aa, _codon in _CODON_ORDER_BY_AA:
    _AA_TO_CODONS.setdefault(_aa, []).append(_codon)
_SINGLE_CODON_AA3 = {aa for aa, codons in _AA_TO_CODONS.items() if len(codons) == 1}
_AA_ORDER = list(_AA_TO_CODONS.keys())

# Codon-display groups used by the "Codons to display" selector in the
# codon-usage comparison panel.  The grouping follows synonymous codon boxes:
# four-codon boxes stay together, two-codon boxes stay together, and six-codon
# amino acids are split into their 4-codon and 2-codon boxes (e.g. Leu CTN vs
# Leu TTR, Ser TCN vs Ser AGY, Arg CGN vs Arg AGR).
_CODON_DISPLAY_GROUPS = [
    ("Ala — GCN", ["GCA", "GCC", "GCG", "GCT"]),
    ("Arg — CGN", ["CGA", "CGC", "CGG", "CGT"]),
    ("Arg — AGR", ["AGA", "AGG"]),
    ("Asn — AAY", ["AAC", "AAT"]),
    ("Asp — GAY", ["GAC", "GAT"]),
    ("Cys — TGY", ["TGC", "TGT"]),
    ("Gln — CAR", ["CAA", "CAG"]),
    ("Glu — GAR", ["GAA", "GAG"]),
    ("Gly — GGN", ["GGA", "GGC", "GGG", "GGT"]),
    ("His — CAY", ["CAC", "CAT"]),
    ("Ile — ATH", ["ATA", "ATC", "ATT"]),
    ("Leu — CTN", ["CTA", "CTC", "CTG", "CTT"]),
    ("Leu — TTR", ["TTA", "TTG"]),
    ("Lys — AAR", ["AAA", "AAG"]),
    ("Met — ATG", ["ATG"]),
    ("Phe — TTY", ["TTC", "TTT"]),
    ("Pro — CCN", ["CCA", "CCC", "CCG", "CCT"]),
    ("Ser — TCN", ["TCA", "TCC", "TCG", "TCT"]),
    ("Ser — AGY", ["AGC", "AGT"]),
    ("Thr — ACN", ["ACA", "ACC", "ACG", "ACT"]),
    ("Trp — TGG", ["TGG"]),
    ("Tyr — TAY", ["TAC", "TAT"]),
    ("Val — GTN", ["GTA", "GTC", "GTG", "GTT"]),
]
_CODON_DISPLAY_GROUPS = [(name, [c for c in codons if c in _CODON_ORDER]) for name, codons in _CODON_DISPLAY_GROUPS]
_CODON_DISPLAY_GROUPS = [(name, codons) for name, codons in _CODON_DISPLAY_GROUPS if codons]


def _codon_display_group_summary(selected_codons):
    selected = {str(c).upper().replace("U", "T") for c in list(selected_codons or []) if str(c).strip()}
    full = []
    partial = []
    for name, codons in _CODON_DISPLAY_GROUPS:
        cset = set(codons)
        n = len(cset & selected)
        if n == len(cset):
            full.append(name)
        elif n > 0:
            partial.append(name)
    return full, partial


def _iter_fasta_records_for_compare(path):
    header = None
    chunks = []
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            yield header, "".join(chunks)


def _clean_dna_for_compare(seq):
    s = str(seq or "").upper().replace("U", "T")
    s = re.sub(r"[^ACGT]", "N", s)
    return s[:len(s) - (len(s) % 3)]


def _extract_fasta_record_id(header):
    """Extract the canonical CDS identifier used by CodonPipe from a FASTA header."""
    h = str(header or "").strip()
    m = re.search(r"\[locus_tag=([^\]]+)\]", h)
    if m:
        return m.group(1).strip()
    m = re.search(r"\[protein_id=([^\]]+)\]", h)
    if m:
        return m.group(1).strip()
    return h.split()[0].strip() if h else ""


def _extract_gene_symbol_from_header(header):
    h = str(header or "").strip()
    for key in ("gene", "gene_synonym"):
        m = re.search(r"\[" + key + r"=([^\]]+)\]", h)
        if m:
            return m.group(1).strip()
    return ""


def _extract_product_from_header(header):
    h = str(header or "").strip()
    for key in ("product", "protein"):
        m = re.search(r"\[" + key + r"=([^\]]+)\]", h)
        if m:
            return m.group(1).strip()
    # As a fallback, keep the descriptive text after the primary ID.
    parts = h.split(None, 1)
    return parts[1].strip() if len(parts) > 1 else ""


def _safe_output_name(text, fallback="output", max_len=120):
    s = str(text or fallback).strip()
    s = re.sub(r"[^A-Za-z0-9._ -]+", "_", s)
    s = re.sub(r"\s+", " ", s).strip(" ._-)")
    s = s[:int(max_len)].strip(" ._-")
    return s or fallback


def _suggest_group_name_from_paths(paths, fallback="FASTA group", max_len=42):
    """Suggest a readable comparison-group name from one or more FASTA paths."""
    paths = _split_custom_cds_paths_gui(paths)
    stems = []
    for path in paths:
        stem = os.path.splitext(os.path.basename(str(path)))[0].strip()
        stem = re.sub(r"(?i)(cds[_ -]*from[_ -]*genomic|cds|genomic|sequence|sequences|fasta|fna|ffn)", "", stem)
        stem = re.sub(r"[_-]+", " ", stem)
        stem = re.sub(r"\s+", " ", stem).strip(" ._-")
        if stem:
            stems.append(stem)
    if not stems:
        return fallback
    if len(stems) == 1:
        return _safe_output_name(stems[0], fallback=fallback, max_len=max_len)

    common = os.path.commonprefix(stems).strip(" ._-")
    common = re.sub(r"\s+", " ", common).strip(" ._-")
    if len(common) >= 5:
        label = f"{common} ({len(stems)} files)"
    elif len(stems) <= 3:
        label = " + ".join(stems)
    else:
        label = f"{len(stems)} FASTA files"
    return _safe_output_name(label, fallback=fallback, max_len=max_len)


def _make_unique_label(label, used_labels, fallback="Group"):
    base = str(label or fallback).strip() or fallback
    candidate = base
    suffix = 2
    used = {str(v).strip().lower() for v in list(used_labels or [])}
    while candidate.strip().lower() in used:
        candidate = f"{base} ({suffix})"
        suffix += 1
    return candidate


def _wrap_fasta(seq, width=80):
    s = str(seq or "")
    return "\n".join(s[i:i + width] for i in range(0, len(s), int(width)))


def _codon_count_rows_from_fastas(paths):
    """Return a DataFrame of per-CDS sense-codon counts for one or more FASTA files."""
    paths = _split_custom_cds_paths_gui(paths)
    rows = []
    names = []
    for path in paths:
        if not path:
            continue
        if not os.path.isfile(path):
            raise FileNotFoundError(path)
        for idx, (header, seq) in enumerate(_iter_fasta_records_for_compare(path), start=1):
            clean = _clean_dna_for_compare(seq)
            counts = dict.fromkeys(_CODON_ORDER, 0)
            valid = 0
            for i in range(0, len(clean), 3):
                codon = clean[i:i+3]
                if codon in counts:
                    counts[codon] += 1
                    valid += 1
            if valid <= 0:
                continue
            label = str(header or "").split()[0] or f"CDS_{idx}"
            names.append(label)
            rows.append(counts)
    if not rows:
        return pd.DataFrame(columns=_CODON_ORDER, dtype=float)
    return pd.DataFrame(rows, index=names, columns=_CODON_ORDER, dtype=float)


def _metric_values_from_codon_counts(counts_df, metric):
    """Convert per-CDS codon counts to ACU or RCU values."""
    metric = str(metric or "Relative codon usage").strip().lower()
    df = counts_df.copy().astype(float)
    if df.empty:
        return df
    if metric.startswith("absolute"):
        denom = df.sum(axis=1).replace(0, np.nan)
        return df.div(denom, axis=0)

    # Relative codon usage within each synonymous family. For Met and Trp, keep
    # the pipeline convention: value = 1 for all CDS.
    out = pd.DataFrame(index=df.index, columns=_CODON_ORDER, dtype=float)
    for aa, codons in _AA_TO_CODONS.items():
        if aa in _SINGLE_CODON_AA3:
            out[codons[0]] = 1.0
        else:
            denom = df[codons].sum(axis=1).replace(0, np.nan)
            out[codons] = df[codons].div(denom, axis=0)
    return out


def _aa_frequency_values_from_codon_counts(counts_df):
    """Return per-CDS absolute amino-acid frequencies for the 20 amino acids."""
    df = counts_df.copy().astype(float)
    if df.empty:
        return pd.DataFrame(columns=_AA_ORDER, dtype=float)
    aa_counts = pd.DataFrame(index=df.index, columns=_AA_ORDER, dtype=float)
    for aa in _AA_ORDER:
        codons = [c for c in _AA_TO_CODONS.get(aa, []) if c in df.columns]
        aa_counts[aa] = df[codons].sum(axis=1) if codons else 0.0
    denom = aa_counts.sum(axis=1).replace(0, np.nan)
    return aa_counts.div(denom, axis=0).loc[:, _AA_ORDER]


def _zcu_values_from_rcu(rcu_df, genome_rcu_df):
    """Compute RCU z-scores for custom CDS against the selected genome baseline."""
    if rcu_df is None or rcu_df.empty:
        return pd.DataFrame(columns=_CODON_ORDER, dtype=float)
    if genome_rcu_df is None or genome_rcu_df.empty:
        return pd.DataFrame(index=rcu_df.index, columns=_CODON_ORDER, dtype=float)
    mu = genome_rcu_df.loc[:, _CODON_ORDER].mean(axis=0, skipna=True)
    sd = genome_rcu_df.loc[:, _CODON_ORDER].std(axis=0, skipna=True, ddof=0)
    sd2 = sd.replace(0.0, np.nan)
    z = (rcu_df.loc[:, _CODON_ORDER] - mu) / sd2
    # Match the pipeline convention: zero-variance but defined features become 0.
    zero_cols = sd.index[sd == 0.0].tolist()
    for col in zero_cols:
        if col in z.columns:
            mask = rcu_df[col].notna()
            z.loc[mask, col] = 0.0
    return z.loc[:, _CODON_ORDER]


def _write_codon_usage_raw_tables_xlsx(out_path, counts_df, genome_rcu_df=None):
    """Write per-CDS codon count, ACU, RCU and ZCU tables for a comparison group."""
    counts = counts_df.copy().loc[:, _CODON_ORDER]
    acu = _metric_values_from_codon_counts(counts, "Absolute codon frequency")
    rcu = _metric_values_from_codon_counts(counts, "Relative codon usage")
    zcu = _zcu_values_from_rcu(rcu, genome_rcu_df)
    aa_freq = _aa_frequency_values_from_codon_counts(counts)

    def _prep(df, round_values=True):
        out = df.copy()
        out.index.name = "CDS"
        if round_values:
            out = out.round(6)
        return out.reset_index()

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with pd.ExcelWriter(out_path, engine="xlsxwriter") as writer:
        _prep(counts, round_values=False).to_excel(writer, sheet_name="Codon counts", index=False)
        _prep(acu).to_excel(writer, sheet_name="ACU", index=False)
        _prep(rcu).to_excel(writer, sheet_name="RCU", index=False)
        _prep(zcu).to_excel(writer, sheet_name="ZCU", index=False)
        _prep(aa_freq).to_excel(writer, sheet_name="AA frequency", index=False)
    return out_path


def _mean_sd_for_compare(values_df, feature_order=None):
    features = list(feature_order or _CODON_ORDER)
    if values_df is None or values_df.empty:
        return (
            np.full(len(features), np.nan, dtype=float),
            np.full(len(features), np.nan, dtype=float),
            0,
        )

    aligned = pd.DataFrame(index=values_df.index, columns=features, dtype=float)
    for feature in features:
        if feature in values_df.columns:
            aligned[feature] = pd.to_numeric(values_df[feature], errors="coerce")
    arr = aligned.to_numpy(dtype=float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        mean = np.nanmean(arr, axis=0)
        if arr.shape[0] > 1:
            sd = np.nanstd(arr, axis=0, ddof=1)
        else:
            sd = np.full(arr.shape[1], np.nan, dtype=float)
    return mean, sd, int(arr.shape[0])


def _show_matplotlib_figure_nonblocking(fig):
    try:
        backend = (plt.get_backend() or "").lower()
    except Exception:
        backend = ""
    noninteractive_agg = ("agg" in backend) and ("inline" not in backend) and ("nbagg" not in backend) and ("ipympl" not in backend)
    if noninteractive_agg:
        return
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


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def _split_keywords_text(value):
    text = str(value or "").replace("\n", ";")
    raw = re.split(r"[;,]+", text)
    out = []
    seen = set()
    for item in raw:
        kw = " ".join(str(item or "").strip().split())
        if not kw:
            continue
        key = kw.lower()
        if key in seen:
            continue
        seen.add(key)
        out.append(kw)
    return out


def _split_selection_text(value):
    """Parse a GUI/list setting into a de-duplicated list. None means all/default."""
    if value is None:
        return None
    if isinstance(value, str):
        text = value.strip()
        if text == "" or text.lower() in {"all", "any", "default"}:
            return None
        raw = re.split(r"[,;\n\r\t]+", text)
    else:
        raw = list(value or [])
    out = []
    seen = set()
    for item in raw:
        val = " ".join(str(item or "").strip().split())
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
PLOT6_DEFAULT_AA_KEYS = {re.sub(r"[^A-Za-z0-9]+", "", aa).lower() for aa in PLOT6_DEFAULT_AAS}
PLOT6_DEFAULT_EXCLUDED_AA_KEYS = {"sec", "selenocysteine", "trp", "tryptophan", "met", "methionine"}


def _plot6_feature_key_ascii_gui(value):
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


def _canonicalize_gui_plot6_modification_feature(feature):
    s = str(feature or "").strip()
    if not s:
        return ""
    key = _plot6_feature_key_ascii_gui(s)
    if key in {"cmo5u34", "mcmo5u34"}:
        return "(m)cmo5U34"
    if key in {"cm32", "um32", "cm32um32"}:
        return "Cm32/Um32"
    if key in {"q34", "gluq34"}:
        return "Q34"
    if key in {"cm34", "cmnm5um34", "cmmn5um34"}:
        # cmnm5Um34 is expanded to both Cm34 and mnm5U34 in the backend.
        # A single GUI text selection cannot expand to two boxes, so this helper
        # maps it to the methylation bin for checklist display/deduplication.
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




def _gui_plot6_modification_position_key(feature):
    s = str(feature or '').strip()
    if not s:
        return ''
    for pos in ('32', '34', '37'):
        if re.search(rf'(?<!\d){pos}(?!\d)', s):
            return pos
    return ''


def _gui_plot6_modification_sort_key(feature):
    s = str(feature or '').strip()
    pos = _gui_plot6_modification_position_key(s)
    return ({'32': 0, '34': 1, '37': 2}.get(pos, 3), s.lower())


def _canonicalize_gui_plot6_modification_selection(values):
    vals = _split_selection_text(values)
    if vals is None:
        return None
    out, seen = [], set()
    for v in vals:
        c = _canonicalize_gui_plot6_modification_feature(v)
        if not c:
            continue
        key = _plot6_feature_key_ascii_gui(c)
        if key in seen:
            continue
        seen.add(key)
        out.append(c)
    return out


def _is_default_excluded_gui_plot6_modification(feature):
    return _plot6_feature_key_ascii_gui(feature) in PLOT6_DEFAULT_EXCLUDED_MODIFICATION_KEYS


def _plot6_aa_key_gui(value):
    return re.sub(r"[^A-Za-z0-9]+", "", str(value or "").strip()).lower()


def _is_default_excluded_gui_plot6_aa(aa):
    return _plot6_aa_key_gui(aa) in PLOT6_DEFAULT_EXCLUDED_AA_KEYS


def _default_gui_plot6_included_aas(available):
    available_clean = [str(aa).strip() for aa in list(available or []) if str(aa).strip()]
    by_key = {_plot6_aa_key_gui(aa): aa for aa in available_clean}
    selected = [by_key[key] for key in [_plot6_aa_key_gui(aa) for aa in PLOT6_DEFAULT_AAS] if key in by_key]
    if selected:
        return selected
    return [aa for aa in available_clean if _plot6_aa_key_gui(aa) in PLOT6_DEFAULT_AA_KEYS]


def _cluster_mode_internal(mode):
    return USER_CLUSTER_TO_INTERNAL.get(str(mode or "").strip(), "basic")


def _cluster_mode_terms_source(mode):
    return USER_CLUSTER_TO_TERMS_SOURCE.get(str(mode or "").strip(), "geneids")


def _safe_int(value, default):
    try:
        return int(float(str(value).strip()))
    except Exception:
        return int(default)


def _safe_float(value, default):
    try:
        return float(str(value).strip())
    except Exception:
        return float(default)


def _optional_limit_value(value):
    s = str(value).strip()
    if s == "":
        return None
    try:
        return float(s)
    except Exception:
        return None


def _validate_codon_range_text(value):
    """Validate the GUI codon-range field.

    Returns the normalized text to pass to the pipeline. Coordinates are 1-based
    and inclusive; e.g. "1-20" means the first 20 codons.
    """
    s = str(value or "all").strip() or "all"
    sl = s.lower().replace("–", "-").replace("—", "-")
    sl = re.sub(r"\s+", " ", sl)
    if sl in {"", "all", "full", "whole", "entire", "entire gene", "entire cds"}:
        return "all"
    if re.match(r"^first\s+\d+$", sl):
        return s
    if re.match(r"^\d+$", sl):
        return s
    m = re.match(r"^(\d+)\s*(?:-|:|\.\.)\s*(\d+|end|all)?$", sl)
    if not m:
        raise ValueError("Codon range must be 'all', '1-20', '20-200', '20-end', or 'first 20'.")
    start = int(m.group(1))
    raw_end = m.group(2)
    if start < 1:
        raise ValueError("Codon range start must be >= 1.")
    if raw_end not in (None, "", "end", "all"):
        end = int(raw_end)
        if end < start:
            raise ValueError("Codon range end must be greater than or equal to the start.")
    return s


def _safe_bool(value):
    if isinstance(value, bool):
        return value
    s = str(value).strip().lower()
    return s in {"1", "true", "yes", "y", "on"}


def _split_custom_cds_paths_gui(value):
    """Parse one or more FASTA paths from GUI text or a Python sequence."""
    if isinstance(value, (list, tuple, set)):
        raw_items = []
        for item in value:
            if isinstance(item, (list, tuple, set)):
                raw_items.extend(list(item))
            else:
                raw_items.append(item)
    else:
        text = str(value or "").strip()
        if not text:
            return []
        raw_items = []
        for line in text.splitlines():
            raw_items.extend(line.split(";"))

    out, seen = [], set()
    for item in raw_items:
        path = str(item or "").strip().strip('"').strip("'")
        if not path:
            continue
        path = os.path.abspath(os.path.expanduser(path))
        key = os.path.normcase(path)
        if key in seen:
            continue
        seen.add(key)
        out.append(path)
    return out


def _usage_internal(value):
    return USAGE_DISPLAY_TO_INTERNAL.get(str(value).strip(), str(value).strip().upper())


def _dimred_internal(value):
    return DIMRED_DISPLAY_TO_INTERNAL.get(str(value).strip(), str(value).strip().lower())


def _bool_to_choice(value, enabled_label):
    return enabled_label if bool(value) else "disabled"

def _smooth_method_display(value, default="running average"):
    s = str(value or "").strip().lower().replace("_", " ").replace("-", " ")
    if s in {"", "true", "yes", "on", "1"}:
        return default
    if s in {"false", "no", "off", "0", "none", "no smoothing"}:
        return "none"
    if "median" in s:
        return "running median"
    if "gauss" in s:
        return "gaussian"
    return "running average"


def _smooth_method_to_bool(value):
    return str(value or "").strip().lower() not in {"none", "no", "false", "0", "off", "no smoothing"}


def _read_cluster_file(path, sheet_name=""):
    ext = os.path.splitext(path)[1].lower()
    if ext in (".xlsx", ".xls", ".xlsm"):
        if sheet_name:
            df = pd.read_excel(path, sheet_name=sheet_name, dtype=str)
        else:
            df = pd.read_excel(path, dtype=str)
    else:
        df = pd.read_csv(path, sep=None, engine="python", dtype=str)
    return df.fillna("")


def _cluster_sizes_from_df(cluster_df):
    sizes = {}
    for c in cluster_df.columns:
        vals = cluster_df[c].dropna().astype(str).str.strip()
        vals = vals[vals != ""]
        sizes[str(c)] = int(vals.shape[0])
    return sizes


class TextRedirector(io.TextIOBase):
    def __init__(self, q):
        self.q = q

    def write(self, s):
        if s:
            self.q.put(str(s))
        return len(s)

    def flush(self):
        return



class TeeTextRedirector(io.TextIOBase):
    def __init__(self, write_gui, *streams):
        self.write_gui = write_gui
        self.streams = [s for s in streams if s is not None]

    def write(self, s):
        if not s:
            return 0
        try:
            self.write_gui(str(s))
        except Exception:
            pass
        for st in self.streams:
            try:
                st.write(str(s))
                st.flush()
            except Exception:
                pass
        return len(s)

    def flush(self):
        for st in self.streams:
            try:
                st.flush()
            except Exception:
                pass
        return

class ToolTip:

    def __init__(self, widget, text, enabled_var, delay=450):
        self.widget = widget
        self.text = text
        self.enabled_var = enabled_var
        self.delay = delay
        self.tipwindow = None
        self._after_id = None
        widget.bind("<Enter>", self._schedule, add="+")
        widget.bind("<Leave>", self._hide, add="+")
        widget.bind("<ButtonPress>", self._hide, add="+")

    def _schedule(self, _event=None):
        self._hide()
        if not self.enabled_var.get() or not self.text:
            return
        self._after_id = self.widget.after(self.delay, self._show)

    def _show(self):
        if self.tipwindow is not None or not self.enabled_var.get() or not self.text:
            return
        try:
            x = self.widget.winfo_rootx() + 18
            y = self.widget.winfo_rooty() + self.widget.winfo_height() + 8
        except Exception:
            return
        self.tipwindow = tw = tk.Toplevel(self.widget)
        tw.wm_overrideredirect(True)
        tw.wm_geometry(f"+{x}+{y}")
        label = tk.Label(
            tw,
            text=self.text,
            justify="left",
            relief="solid",
            borderwidth=1,
            padx=8,
            pady=6,
            wraplength=420,
            background="#fff8dc",
            foreground="#202020",
            font=("Arial", 9),
        )
        label.pack()

    def _hide(self, _event=None):
        if self._after_id is not None:
            try:
                self.widget.after_cancel(self._after_id)
            except Exception:
                pass
            self._after_id = None
        if self.tipwindow is not None:
            self.tipwindow.destroy()
            self.tipwindow = None


class ScrollableFrame(ttk.Frame):
    """A vertical scrollable frame with safe mouse-wheel routing.

    Earlier GUI versions used ``canvas.bind_all(<MouseWheel>, self._on_mousewheel)``
    in every scrollable frame. That works until a temporary Toplevel window
    containing a ScrollableFrame is destroyed: Tk can keep the old global callback
    around briefly, and the next wheel event may try to access a widget path that
    no longer exists (``TclError: bad window path name``).

    This implementation installs one global mouse-wheel dispatcher per Tk root and
    routes wheel events only to the currently active ScrollableFrame. Destroyed
    frames are ignored safely.
    """

    _active_frame = None
    _wheel_bound_roots = set()

    def __init__(self, parent):
        super().__init__(parent)
        self.canvas = tk.Canvas(self, highlightthickness=0, borderwidth=0)
        self.vsb = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=self.vsb.set)
        self.vsb.pack(side="right", fill="y")
        self.canvas.pack(side="left", fill="both", expand=True)

        self.inner = ttk.Frame(self.canvas)
        self._window_id = self.canvas.create_window((0, 0), window=self.inner, anchor="nw")

        self.inner.bind("<Configure>", self._on_inner_configure)
        self.canvas.bind("<Configure>", self._on_canvas_configure)

        # Activate this scroller when the pointer is over it.  The actual wheel
        # callback is class-level and bound to the root, not to this widget, so
        # destroyed popups cannot leave broken callbacks behind.
        for widget in (self, self.canvas, self.inner):
            widget.bind("<Enter>", self._activate_mousewheel, add="+")
            widget.bind("<Leave>", self._maybe_deactivate_mousewheel, add="+")
        self.bind("<Destroy>", self._on_destroy, add="+")
        self._ensure_global_mousewheel_binding()

    def _ensure_global_mousewheel_binding(self):
        try:
            root = self._root()
        except Exception:
            root = self.winfo_toplevel()

        key = str(root)
        if key in ScrollableFrame._wheel_bound_roots:
            return

        ScrollableFrame._wheel_bound_roots.add(key)
        try:
            root.bind_all("<MouseWheel>", ScrollableFrame._dispatch_mousewheel, add="+")
            # Linux/X11 sends Button-4/Button-5 instead of MouseWheel.
            root.bind_all("<Button-4>", ScrollableFrame._dispatch_mousewheel, add="+")
            root.bind_all("<Button-5>", ScrollableFrame._dispatch_mousewheel, add="+")
        except Exception:
            pass

    def _activate_mousewheel(self, _event=None):
        ScrollableFrame._active_frame = self

    def _maybe_deactivate_mousewheel(self, _event=None):
        # Do not deactivate when moving between child widgets inside the same
        # scrollable frame.
        try:
            x, y = self.winfo_pointerxy()
            widget = self.winfo_containing(x, y)
            while widget is not None:
                if widget in (self, self.canvas, self.inner):
                    return
                widget = getattr(widget, "master", None)
        except Exception:
            pass

        if ScrollableFrame._active_frame is self:
            ScrollableFrame._active_frame = None

    def _on_destroy(self, _event=None):
        if ScrollableFrame._active_frame is self:
            ScrollableFrame._active_frame = None

    def _on_inner_configure(self, _event=None):
        try:
            if not self.winfo_exists():
                return
            self.canvas.configure(scrollregion=self.canvas.bbox("all"))
        except tk.TclError:
            return

    def _on_canvas_configure(self, event=None):
        try:
            if event is not None and self.winfo_exists():
                self.canvas.itemconfigure(self._window_id, width=event.width)
        except tk.TclError:
            return

    @classmethod
    def _dispatch_mousewheel(cls, event):
        sf = cls._active_frame
        if sf is None:
            return
        try:
            if not sf.winfo_exists() or not sf.winfo_ismapped():
                cls._active_frame = None
                return
        except tk.TclError:
            cls._active_frame = None
            return

        try:
            if getattr(event, "num", None) == 4:
                units = -1
            elif getattr(event, "num", None) == 5:
                units = 1
            else:
                delta = getattr(event, "delta", 0)
                units = int(-1 * (delta / 120)) if delta else 0
            if units:
                sf.canvas.yview_scroll(units, "units")
        except tk.TclError:
            cls._active_frame = None
        except Exception:
            pass


# -----------------------------------------------------------------------------
# GUI app
# -----------------------------------------------------------------------------

class CodonPipeGUI:
    def __init__(self, root):
        self.root = root
        self.root.title(f"CodonPipe — Python Interface [{CODONPIPE_GUI_BUILD}]")
        self.root.geometry("1580x980")
        self.root.minsize(1250, 820)

        self.log_queue = queue.Queue()
        self.worker = None
        self.is_running = False
        self.tooltip_refs = []
        self._size_pair_guard = False
        self.run_buttons = []
        self.last_clustering_workbook = ""
        self.last_clustering_output_dir = ""
        self._auto_loaded_cluster_file = ""
        self._preloaded_applied_key = ""
        self._preloaded_companion_cache = {}
        # Keep startup fast: never read Excel workbooks just to build labels.
        # Preloaded FASTA + companion xlsx paths are applied immediately, but
        # companion detection is path-based and top-level only.
        self._startup_building = True
        self._defer_preloaded_companion_autoload = False

        self._build_variables()
        self._build_ui()
        self._bind_size_lock_pairs()
        self._sync_codon_set_from_usage_basis()
        self._refresh_genome_source_widgets()
        self._refresh_cluster_source_widgets()
        self._refresh_dimred_params()
        self._refresh_cluster_params()
        self._refresh_optional_sections()
        self._update_active_clusters_status()
        self._startup_building = False
        self._poll_log_queue()

    # ------------------------- variables -------------------------
    def _build_variables(self):
        S = CP.SET
        P = CP.PIPELINE

        # Tooltips are always enabled; no user-facing toggle is displayed.
        self.tooltips_enabled = tk.BooleanVar(value=True)
        self.simplified_interface = tk.BooleanVar(value=True)

        self.default_root = tk.StringVar(value=str(S.get("default_root", "") or DEFAULT_BROWSER_ROOT))
        self.preloaded_genomes = _discover_preloaded_genomes()
        default_genome = _default_preloaded_genome_label(self.preloaded_genomes)
        self.preloaded_genome_choice = tk.StringVar(value=default_genome)
        self.genome_not_available = tk.BooleanVar(value=not bool(self.preloaded_genomes))
        self.fasta_path = tk.StringVar(value=str(S.get("fasta_path", "")))
        self.organism_name = tk.StringVar(value="Organism")
        self.add_custom_cds_enable = tk.BooleanVar(value=bool(S.get("custom_cds_enable", False)))
        custom_paths_default = S.get("custom_cds_paths", [])
        if isinstance(custom_paths_default, (list, tuple, set)):
            custom_paths_default = "; ".join([str(x) for x in custom_paths_default])
        self.custom_cds_paths = tk.StringVar(value=str(custom_paths_default or ""))
        self.include_custom_cds_cluster = tk.BooleanVar(value=True)
        self.custom_cds_cluster_name = tk.StringVar(value=str(S.get("custom_cds_cluster_name", "custom") or "custom"))

        self.codon_compare_metric = tk.StringVar(value="Relative codon usage")
        self.codon_compare_rcu_display = tk.StringVar(value="Relative codon usage with genome")
        self.codon_compare_statistics = tk.StringVar(value="None")
        self.codon_compare_plot_style = tk.StringVar(value="Mean ± SD")
        self.codon_compare_selected_codons = list(_CODON_ORDER)
        self.codon_compare_selected_codons_status = tk.StringVar(value="All codons selected")
        self.codon_compare_highlighted_codons = []
        self.codon_compare_highlighted_codons_status = tk.StringVar(value="No codon highlighted")
        self.codon_compare_log2_y_scale = tk.BooleanVar(value=False)
        self.codon_compare_groups = []

        # Figure-details controls for the lightweight plots generated from the
        # "Codon usage analyses" tab. These are intentionally separate from the
        # per-cluster codon-usage profile controls used by the main Figures tab.
        self.codon_compare_raw_custom_axes = tk.BooleanVar(value=False)
        self.codon_compare_raw_xmin = tk.StringVar(value="")
        self.codon_compare_raw_xmax = tk.StringVar(value="")
        self.codon_compare_raw_ymin = tk.StringVar(value="")
        self.codon_compare_raw_ymax = tk.StringVar(value="")
        self.codon_compare_corr_custom_axes = tk.BooleanVar(value=False)
        self.codon_compare_corr_xmin = tk.StringVar(value="")
        self.codon_compare_corr_xmax = tk.StringVar(value="")
        self.codon_compare_corr_ymin = tk.StringVar(value="")
        self.codon_compare_corr_ymax = tk.StringVar(value="")
        # Shared styling for lightweight codon/amino-acid analysis plots.
        # Caption size is intentionally large by default for publication/readability.
        self.codon_compare_caption_size = tk.StringVar(value="14")
        self.codon_compare_marker_size = tk.StringVar(value="8.0")
        self.codon_compare_line_width = tk.StringVar(value="2.0")
        self.codon_compare_aa_spacing = tk.StringVar(value="1.75")
        self.codon_compare_codon_aa_gap = tk.StringVar(value="0.65")
        self.codon_compare_codon_gap = tk.StringVar(value="0.33")
        self.codon_compare_highlight_color = tk.StringVar(value="#0057D9")
        self.codon_compare_legend_ncol = tk.StringVar(value="3")

        self.fasta_metric_default_defs = get_fasta_metric_group_defaults()
        self.fasta_metric_group_rows = []
        self.fasta_metric_cluster_df = pd.DataFrame()
        self.fasta_metric_cluster_scores_df = pd.DataFrame()
        self.fasta_metric_cluster_path = ""
        self.fasta_metric_cluster_status = tk.StringVar(value="No FASTA-derived metric group selected")
        self.extract_clusters_selection = None
        self.extract_clusters_status = tk.StringVar(value="All available clusters selected")

        # Default to FASTA-inferred clusters for casual users. The DAVID-gene2terms
        # mode remains available from the dropdown.
        default_cluster_mode = "Inferred from FASTA file"
        self.user_cluster_mode = tk.StringVar(value=default_cluster_mode)
        self.david_gene2terms_path = tk.StringVar(value=str(S.get("david_gene2terms_path", "")))
        self.refined_cluster_file = tk.StringVar(value="")
        self.refined_cluster_sheet = tk.StringVar(value="")
        self.active_clusters_selection = None  # None means all currently available clusters
        self.active_clusters_status = tk.StringVar(value="All clusters selected")
        self.figure_clusters_selection = None  # None means all clusters in the existing workbook are plotted from the Figures tab
        self.figure_clusters_status = tk.StringVar(value="All figure clusters selected")
        self.decoding_clusters_selection = []  # Default: auto-preselect Ribosomal proteins + Virulence when both exist
        self.decoding_clusters_selection_user_set = False
        self.decoding_clusters_status = tk.StringVar(value="No decoding clusters selected")
        self.decoding_reference_cluster = tk.StringVar(value=str(CP.SET.get("decoding_reference_cluster", "") or ""))
        self.decoding_reference_cluster_status = tk.StringVar(value="No reference cluster selected")
        self.keyword_group_vars = {}
        self._init_keyword_group_vars()

        self.usage_basis = tk.StringVar(value="Relative codon usage")
        self.fasta_codon_range = tk.StringVar(value=str(S.get("fasta_codon_range", "all") or "all"))
        self.codon_set = tk.StringVar(value=str(CP.RUNTIME_DEFAULTS.get("codon_set", "61")))
        self.dimred_method = tk.StringVar(value=DIMRED_INTERNAL_TO_DISPLAY.get(str(CP.RUNTIME_DEFAULTS.get("dimred_method", "umap")).lower(), "UMAP"))
        self.cluster_method = tk.StringVar(value=str(CP.RUNTIME_DEFAULTS.get("cluster_method", "kmeans")))
        self.statistical_test_method = tk.StringVar(value=_bool_to_choice(CP.RUNTIME_DEFAULTS.get("do_2d_ks", False), "2D Kolmogorov-Smirnov"))
        self.functional_scan_method = tk.StringVar(value=_bool_to_choice(CP.RUNTIME_DEFAULTS.get("run_david_scan", False), "DAVID sliding-window scan"))

        self.center_features = tk.BooleanVar(value=bool(S.get("center_features", True)))
        self.scale_features = tk.BooleanVar(value=bool(S.get("scale_features", True)))

        self.dimred_param_vars = {
            "umap": {
                "umap_neighbors": tk.StringVar(value=str(S.get("umap_neighbors", 10))),
                "umap_min_dist": tk.StringVar(value=str(S.get("umap_min_dist", 0.01))),
                "umap_metric": tk.StringVar(value=str(S.get("umap_metric", "cosine"))),
                "umap_components": tk.StringVar(value=str(S.get("umap_components", 2))),
                "umap_randomize": tk.BooleanVar(value=bool(S.get("umap_randomize", False))),
                "umap_clip_abs": tk.StringVar(value=str(S.get("umap_clip_abs", 0.0))),
                "umap_init": tk.StringVar(value=str(S.get("umap_init", "spectral"))),
            },
            "tsne": {
                "tsne_perplexity": tk.StringVar(value=str(S.get("tsne_perplexity", 10))),
                "tsne_distance": tk.StringVar(value=str(S.get("tsne_distance", "cosine"))),
                "tsne_dims": tk.StringVar(value=str(S.get("tsne_dims", 2))),
                "tsne_exaggeration": tk.StringVar(value=str(S.get("tsne_exaggeration", 10))),
                "tsne_learnrate": tk.StringVar(value=str(S.get("tsne_learnrate", 100))),
            },
            "pca": {
                "pca_npcs": tk.StringVar(value=str(S.get("pca_npcs", 3))),
                "pca_center": tk.BooleanVar(value=bool(S.get("pca_center", True))),
                "pca_scale": tk.BooleanVar(value=bool(S.get("pca_scale", True))),
            },
        }
        self.dimred_specs = {
            "umap": [
                ("umap_neighbors", "Neighbors", "entry", None, "Number of nearest neighbours used by UMAP."),
                ("umap_min_dist", "Min dist", "entry", None, "How tightly UMAP packs nearby points."),
                ("umap_metric", "Metric", "combo", ["cosine", "euclidean", "manhattan", "correlation"], "Distance metric used by UMAP."),
                ("umap_components", "Components", "entry", None, "Usually 2 for plotting."),
                ("umap_randomize", "Randomize seed", "check", None, "Allow some stochastic variation between runs."),
                ("umap_clip_abs", "Clip abs", "entry", None, "Optional absolute clipping before UMAP."),
                ("umap_init", "Init", "combo", ["spectral", "random"], "UMAP initialization strategy."),
            ],
            "tsne": [
                ("tsne_perplexity", "Perplexity", "entry", None, "Effective neighborhood size for t-SNE."),
                ("tsne_distance", "Distance", "combo", ["cosine", "euclidean", "manhattan", "correlation"], "Distance metric used by t-SNE."),
                ("tsne_dims", "Dimensions", "entry", None, "Usually 2 for plotting."),
                ("tsne_exaggeration", "Exaggeration", "entry", None, "Early exaggeration factor for t-SNE."),
                ("tsne_learnrate", "Learning rate", "entry", None, "t-SNE learning rate."),
            ],
            "pca": [
                ("pca_npcs", "n PCs", "entry", None, "Number of principal components to compute."),
                ("pca_center", "Center", "check", None, "Center variables before PCA."),
                ("pca_scale", "Scale", "check", None, "Scale variables before PCA."),
            ],
        }

        self.cluster_param_vars = {
            "kmeans": {
                "kmeans_k": tk.StringVar(value=str(S.get("kmeans_k", 12))),
                "gene_dist_metric": tk.StringVar(value=str(S.get("gene_dist_metric", "euclidean"))),
            },
            "kmedoids": {
                "kmedoids_k": tk.StringVar(value=str(S.get("kmedoids_k", 12))),
                "kmedoids_dist": tk.StringVar(value=str(S.get("kmedoids_dist", "euclidean"))),
                "gene_dist_metric": tk.StringVar(value=str(S.get("gene_dist_metric", "euclidean"))),
            },
            "hierarchical": {
                "gene_dist_metric": tk.StringVar(value=str(S.get("gene_dist_metric", "euclidean"))),
                "gene_linkage": tk.StringVar(value=str(S.get("gene_linkage", "single"))),
                "cluster_use_optimal_leaf_ordering": tk.BooleanVar(value=bool(S.get("cluster_use_optimal_leaf_ordering", True))),
                "hierarchical_optimal_leaf_max_size": tk.StringVar(value=str(S.get("hierarchical_optimal_leaf_max_size", 400))),
                "hierarchical_fast_order_threshold": tk.StringVar(value=str(S.get("hierarchical_fast_order_threshold", 2000))),
            },
            "spectral": {
                "spectral_k": tk.StringVar(value=str(S.get("spectral_k", 12))),
                "gene_dist_metric": tk.StringVar(value=str(S.get("gene_dist_metric", "euclidean"))),
            },
            "dbscan": {
                "dbscan_eps": tk.StringVar(value=str(S.get("dbscan_eps", 1.5))),
                "dbscan_minpts": tk.StringVar(value=str(S.get("dbscan_minpts", 10))),
                "dbscan_dist": tk.StringVar(value=str(S.get("dbscan_dist", "euclidean"))),
                "gene_dist_metric": tk.StringVar(value=str(S.get("gene_dist_metric", "euclidean"))),
            },
        }
        self.cluster_specs = {
            "kmeans": [
                ("kmeans_k", "K", "entry", None, "Number of clusters for k-means."),
                ("gene_dist_metric", "Ordering metric", "combo", ["euclidean", "cosine", "correlation", "cityblock"], "Metric used to order genes within clusters."),
            ],
            "kmedoids": [
                ("kmedoids_k", "K", "entry", None, "Number of clusters for k-medoids."),
                ("kmedoids_dist", "Medoid metric", "combo", ["euclidean", "cosine", "manhattan", "correlation"], "Distance metric used by k-medoids."),
                ("gene_dist_metric", "Ordering metric", "combo", ["euclidean", "cosine", "correlation", "cityblock"], "Metric used to order genes within clusters."),
            ],
            "hierarchical": [
                ("gene_dist_metric", "Gene distance", "combo", ["euclidean", "cosine", "correlation", "cityblock"], "Distance metric used for hierarchical clustering."),
                ("gene_linkage", "Linkage", "combo", ["single", "average", "complete", "ward"], "Linkage method."),
                ("cluster_use_optimal_leaf_ordering", "Use optimal leaf ordering", "check", None, "Improves leaf ordering for smaller datasets."),
                ("hierarchical_optimal_leaf_max_size", "OLO max size", "entry", None, "Maximum cluster size for optimal leaf ordering."),
                ("hierarchical_fast_order_threshold", "Fast threshold", "entry", None, "Switch threshold for faster ordering."),
            ],
            "spectral": [
                ("spectral_k", "K", "entry", None, "Number of clusters for spectral clustering."),
                ("gene_dist_metric", "Ordering metric", "combo", ["euclidean", "cosine", "correlation", "cityblock"], "Metric used to order genes within clusters."),
            ],
            "dbscan": [
                ("dbscan_eps", "eps", "entry", None, "Neighbourhood radius for DBSCAN."),
                ("dbscan_minpts", "min points", "entry", None, "Minimum local points for DBSCAN."),
                ("dbscan_dist", "DBSCAN metric", "combo", ["euclidean", "cosine", "manhattan", "correlation"], "Distance metric used by DBSCAN."),
                ("gene_dist_metric", "Ordering metric", "combo", ["euclidean", "cosine", "correlation", "cityblock"], "Metric used to order genes within clusters."),
            ],
        }

        self.enable_2d_ks = tk.BooleanVar(value=bool(CP.RUNTIME_DEFAULTS.get("do_2d_ks", False)))
        self.ks_alpha = tk.StringVar(value=str(CP.KS_SETTINGS.get("alpha", 0.01)))
        self.ks_method = tk.StringVar(value=str(CP.KS_SETTINGS.get("method", "binned")))
        self.ks_bins = tk.StringVar(value=str(CP.KS_SETTINGS.get("bins", 151)))
        self.ks_n_perm = tk.StringVar(value=str(CP.KS_SETTINGS.get("n_perm", 2000)))
        self.ks_seed = tk.StringVar(value=str(CP.KS_SETTINGS.get("random_seed", 42)))

        self.enable_david = tk.BooleanVar(value=bool(CP.RUNTIME_DEFAULTS.get("run_david_scan", False)))
        self.david_email = tk.StringVar(value=str(S.get("david_user_email", "")))
        self.david_window_size = tk.StringVar(value=str(S.get("david_window_size", 100)))
        self.david_step_size = tk.StringVar(value=str(S.get("david_step_size", 50)))
        self.david_wait_time = tk.StringVar(value=str(S.get("david_wait_time", 0.0)))
        self.david_max_clusters = tk.StringVar(value=str(S.get("david_max_clusters", 3)))
        self.david_min_valid_ids = tk.StringVar(value=str(S.get("david_min_valid_ids_per_window", 3)))
        self.david_top_n_hits = tk.StringVar(value=str(S.get("david_top_n_hits", 10)))

        self.figure_format = tk.StringVar(value="png")
        self.enable_main_heatmap = tk.BooleanVar(value=bool(S.get("plot_codon_gene_heatmap_enable", True)))
        self.main_heatmap_show_fig = tk.BooleanVar(value=bool(S.get("plot_codon_gene_heatmap_show_fig", S.get("show_main_pipeline_figures", True))))
        self.main_heatmap_custom_aesthetics = tk.BooleanVar(value=False)
        self.main_heatmap_custom_axes = tk.BooleanVar(value=False)
        self.main_heatmap_dpi = tk.StringVar(value=str(S.get("figure_dpi", 300)))
        self.main_heatmap_colormap = tk.StringVar(value=str(S.get("heatmap_colormap_name", "parula")))
        self.main_heatmap_fig_width = tk.StringVar(value=str((S.get("heatmap_fig_size", (18, 4)) or (18, 4))[0]))
        self.main_heatmap_fig_height = tk.StringVar(value=str((S.get("heatmap_fig_size", (18, 4)) or (18, 4))[1]))
        caxis = S.get("heatmap_caxis_limits", (-0.5, 2.5)) or (-0.5, 2.5)
        self.main_heatmap_caxis_min = tk.StringVar(value=str(caxis[0]))
        self.main_heatmap_caxis_max = tk.StringVar(value=str(caxis[1]))
        self.main_heatmap_xtick_every = tk.StringVar(value=str(S.get("xtick_every_genes", 500)))
        self.main_heatmap_xmin = tk.StringVar(value="")
        self.main_heatmap_xmax = tk.StringVar(value="")
        self.main_heatmap_ymin = tk.StringVar(value="")
        self.main_heatmap_ymax = tk.StringVar(value="")

        self.enable_2d_density_plots = tk.BooleanVar(value=True)
        self.density_figure_dpi = tk.StringVar(value=str(S.get("figure_dpi", 300)))
        self.density_panel_w_in = tk.StringVar(value="5.0")
        self.density_panel_h_in = tk.StringVar(value="5.0")
        self.plot_rows = tk.StringVar(value=str(P.get("plot_max_nrows", 4) or 4))
        self.show_colorbar = tk.BooleanVar(value=True)
        self.show_2d_fig = tk.BooleanVar(value=True)
        self.include_genomic_density_map = tk.BooleanVar(value=True)
        self.color_mode = tk.StringVar(value="enrichment")
        self.density_subplot_wspace = tk.StringVar(value="0.20")
        self.density_subplot_hspace = tk.StringVar(value="0.30")
        self.figure_suptitle = tk.StringVar(value="")
        self.density_cmap = tk.StringVar(value="plasma_r")
        self.enrichment_cmap = tk.StringVar(value="plasma_r")
        self.plot_cluster_min_genes = tk.StringVar(value="2")
        self.density_custom_axes = tk.BooleanVar(value=False)
        self.density_xmin = tk.StringVar(value="")
        self.density_xmax = tk.StringVar(value="")
        self.density_ymin = tk.StringVar(value="")
        self.density_ymax = tk.StringVar(value="")
        self.density_custom_aesthetics = tk.BooleanVar(value=False)

        self.gchm_enable = tk.BooleanVar(value=bool(S.get("gchm_enable", True)))
        self.heatmap_dpi = tk.StringVar(value=str(S.get("gchm_dpi", 300)))
        self.gchm_show_fig = tk.BooleanVar(value=bool(S.get("gchm_show_fig", True)))
        self.gchm_colormap = tk.StringVar(value=str(S.get("gchm_colormap", "plasma")))
        self.gchm_sigma = tk.StringVar(value=str(S.get("gchm_sigma", GCHM_DEFAULT_SIGMA)))
        self.gchm_spread_factor = tk.StringVar(value=str(S.get("gchm_spread_factor", GCHM_DEFAULT_SPREAD_FACTOR)))
        self.gchm_height_per_cluster = tk.StringVar(value=str(S.get("gchm_height_per_cluster", 0.3)))
        self.gchm_label_fontsize = tk.StringVar(value=str(S.get("gchm_label_fontsize", 10)))
        self.gchm_cmap_min_rel = tk.StringVar(value=str(S.get("gchm_cmap_min_rel", 0.2)))
        self.gchm_cmap_max_rel = tk.StringVar(value=str(S.get("gchm_cmap_max_rel", 1.0)))
        self.gchm_custom_axes = tk.BooleanVar(value=False)
        self.gchm_xmin = tk.StringVar(value="")
        self.gchm_xmax = tk.StringVar(value="")
        self.gchm_ymin = tk.StringVar(value="")
        self.gchm_ymax = tk.StringVar(value="")
        self.gchm_custom_aesthetics = tk.BooleanVar(value=False)

        self.apply_smoothing = tk.BooleanVar(value=bool(S.get("apply_smoothing", True)))
        self.smooth_window_genes = tk.StringVar(value=str(S.get("smooth_window_genes", 6)))
        self.apply_binning = tk.BooleanVar(value=bool(S.get("apply_binning", False)))
        self.bin_size_genes = tk.StringVar(value=str(S.get("bin_size_genes", 50)))

        self.enable_codon_usage_plot = tk.BooleanVar(value=True)
        self.codon_usage_plot_mode = tk.StringVar(value=INTERNAL_CODON_MODE_TO_USER.get(str(P.get("codon_usage_plot_mode", "Z")), "ZCU"))
        self.codon_usage_show_fig = tk.BooleanVar(value=True)
        self.codon_usage_dpi = tk.StringVar(value=str(S.get("figure_dpi", 300)))
        self.codon_panel_w_in = tk.StringVar(value="5.0")
        self.codon_panel_h_in = tk.StringVar(value="5.0")
        self.codon_custom_axes = tk.BooleanVar(value=False)
        self.codon_xmin = tk.StringVar(value="")
        self.codon_xmax = tk.StringVar(value="")
        self.codon_ymin = tk.StringVar(value="")
        self.codon_ymax = tk.StringVar(value="")
        self.codon_custom_aesthetics = tk.BooleanVar(value=False)
        self.enable_trna_usage = tk.BooleanVar(value=bool(S.get("export_trna_usage_enable", False) or str(S.get("trna_decoding_table_path", "")).strip()))
        self.enable_trna_abundance_corr = tk.BooleanVar(value=False)
        self.trna_decoding_table_path = tk.StringVar(value=str(S.get("trna_decoding_table_path", "")))
        self.trna_decoding_table_sheet = tk.StringVar(value=str(S.get("trna_decoding_table_sheet", "")))
        self.trna_abundance_sheet = tk.StringVar(value="")
        try:
            self.trna_decoding_table_path.trace_add("write", lambda *_args: self._on_decoding_table_path_changed())
        except Exception:
            pass

        # Gene-ordered supplementary heatmaps
        self.enable_trna_gene_heatmap = tk.BooleanVar(value=bool(S.get("trna_gene_heatmap_enable", True)))
        self.trna_gene_heatmap_metric = tk.StringVar(value=str(S.get("trna_gene_heatmap_metric", "ZTU")))
        self.trna_gene_heatmap_show_fig = tk.BooleanVar(value=bool(S.get("trna_gene_heatmap_show_fig", True)))
        self.enable_trna_single_box_codon_heatmap = tk.BooleanVar(value=bool(S.get("trna_single_box_codon_heatmap_enable", True)))
        self.trna_single_box_codon_heatmap_show_fig = tk.BooleanVar(value=bool(S.get("trna_single_box_codon_heatmap_show_fig", True)))

        # Cluster-level enrichment heatmaps (decodingpipe-like style)
        self.enable_trna_shift_heatmap = tk.BooleanVar(value=bool(S.get("trna_shift_heatmap_enable", True)))
        self.trna_shift_heatmap_show_fig = tk.BooleanVar(value=bool(S.get("trna_shift_heatmap_show_fig", True)))
        self.trna_shift_heatmap_clusters = tk.StringVar(value=str(S.get("trna_shift_heatmap_clusters", "all")))
        self.enable_trna_wobble_heatmap = tk.BooleanVar(value=bool(S.get("trna_wobble_heatmap_enable", True)))
        self.trna_wobble_heatmap_show_fig = tk.BooleanVar(value=bool(S.get("trna_wobble_heatmap_show_fig", True)))
        self.trna_wobble_heatmap_clusters = tk.StringVar(value=str(S.get("trna_wobble_heatmap_clusters", "all")))
        self.trna_shift_heatmap_log2_colorbar = tk.BooleanVar(value=bool(S.get("trna_shift_heatmap_log2_colorbar", True)))
        self.trna_wobble_heatmap_log2_colorbar = tk.BooleanVar(value=bool(S.get("trna_wobble_heatmap_log2_colorbar", True)))
        self.trna_gene_wobble_plot_kind = tk.StringVar(value=str(S.get("trna_gene_wobble_plot_kind", "line") or "line").replace("surface", "area"))
        self.trna_gene_trna_plot_kind = tk.StringVar(value=str(S.get("trna_gene_trna_plot_kind", "heatmap") or "heatmap").replace("surface", "area"))
        self.trna_mrna_stability_plot_kind = tk.StringVar(value=str(S.get("trna_mrna_stability_plot_kind", "line") or "line").replace("surface", "area"))
        self.trna_wobble_plot_kind = tk.StringVar(value=str(S.get("trna_wobble_plot_kind", "boxplot") or "boxplot"))
        self.trna_shift_plot_kind = tk.StringVar(value=str(S.get("trna_shift_plot_kind", "boxplot") or "boxplot"))
        self.trna_modifications_plot_kind = tk.StringVar(value=str(S.get("trna_modifications_plot_kind", "boxplot") or "boxplot"))
        self.trna_gene_wobble_smooth = tk.StringVar(value=_smooth_method_display(S.get("trna_gene_wobble_smooth_method", S.get("trna_gene_wobble_smooth", "running average"))))
        self.trna_gene_wobble_smooth_window = tk.StringVar(value=str(S.get("trna_gene_wobble_smooth_window", 40)))
        self.trna_gene_trna_smooth = tk.StringVar(value=_smooth_method_display(S.get("trna_gene_trna_smooth_method", S.get("trna_gene_trna_smooth", "running average"))))
        self.trna_gene_trna_smooth_window = tk.StringVar(value=str(S.get("trna_gene_trna_smooth_window", 40)))
        self.trna_mrna_stability_smooth = tk.StringVar(value=_smooth_method_display(S.get("trna_mrna_stability_smooth_method", S.get("trna_mrna_stability_smooth", "running average"))))
        self.trna_mrna_stability_smooth_window = tk.StringVar(value=str(S.get("trna_mrna_stability_smooth_window", 100)))
        self.enable_trna_mrna_stability = tk.BooleanVar(value=bool(S.get("trna_mrna_stability_enable", True)))
        self.trna_gene_wobble_caption_size = tk.StringVar(value=str(S.get("trna_gene_wobble_caption_size", 13)))
        self.trna_gene_trna_caption_size = tk.StringVar(value=str(S.get("trna_gene_trna_caption_size", 13)))
        self.trna_mrna_stability_caption_size = tk.StringVar(value=str(S.get("trna_mrna_stability_caption_size", 13)))
        self.trna_wobble_boxplot_caption_size = tk.StringVar(value=str(S.get("trna_wobble_boxplot_caption_size", 13)))
        self.trna_shift_boxplot_caption_size = tk.StringVar(value=str(S.get("trna_shift_boxplot_caption_size", 13)))
        self.trna_modifications_boxplot_caption_size = tk.StringVar(value=str(S.get("trna_modifications_boxplot_caption_size", 17)))
        stat_default = str(S.get("trna_wobble_stats_test", "Student t-test") or "Student t-test")
        self.trna_wobble_stats_test = tk.StringVar(value=stat_default)
        self.trna_shift_stats_test = tk.StringVar(value=str(S.get("trna_shift_stats_test", stat_default) or stat_default))
        self.trna_modifications_stats_test = tk.StringVar(value=str(S.get("trna_modifications_stats_test", stat_default) or stat_default))
        self.trna_wobble_pair_stats_test = tk.StringVar(value=str(S.get("trna_wobble_pair_stats_test", "Student t-test") or "Student t-test"))
        self.trna_shift_pair_stats_test = tk.StringVar(value=str(S.get("trna_shift_pair_stats_test", "Student t-test") or "Student t-test"))
        # Vertical spacing, as a fraction of the y-axis data span, between
        # unbracketed cluster-vs-reference stars and the first bracketed
        # within-codon / within-tRNA comparison.  Lower values move brackets
        # closer to the cluster-comparison stars.
        self.trna_wobble_pair_stats_gap = tk.StringVar(value=str(S.get("trna_wobble_pair_stats_gap", 0.05)))
        self.trna_shift_pair_stats_gap = tk.StringVar(value=str(S.get("trna_shift_pair_stats_gap", 0.05)))
        self.enable_trna_modification_heatmap = tk.BooleanVar(value=bool(S.get("trna_modification_heatmap_enable", True)))
        self.trna_modification_heatmap_show_fig = tk.BooleanVar(value=bool(S.get("trna_modification_heatmap_show_fig", True)))
        self.trna_secondary_axis_style = tk.StringVar(value=str(S.get("trna_secondary_axis_style", "bars") or "bars"))
        self.trna_secondary_axis_alpha = tk.StringVar(value=str(S.get("trna_secondary_axis_alpha", 0.22)))
        self.trna_secondary_axis_bar_width = tk.StringVar(value=str(S.get("trna_secondary_axis_bar_width", 0.72)))
        self.trna_boxplot_width = tk.StringVar(value=str(S.get("trna_boxplot_width", 0.12)))
        self.trna_boxplot_show_points = tk.BooleanVar(value=bool(S.get("trna_boxplot_show_points", True)))
        self.trna_boxplot_point_alpha = tk.StringVar(value=str(S.get("trna_boxplot_point_alpha", 0.35)))
        self.trna_boxplot_point_size = tk.StringVar(value=str(S.get("trna_boxplot_point_size", 10.5)))
        self.trna_wobble_boxplot_style = tk.StringVar(value=str(S.get("trna_wobble_boxplot_style", "boxplot") or "boxplot"))
        self.trna_shift_boxplot_style = tk.StringVar(value=str(S.get("trna_shift_boxplot_style", "boxplot") or "boxplot"))
        self.trna_modifications_boxplot_style = tk.StringVar(value=str(S.get("trna_modifications_boxplot_style", "boxplot") or "boxplot"))
        self.trna_wobble_boxplot_log2 = tk.StringVar(value="yes" if bool(S.get("trna_wobble_boxplot_log2", True)) else "no")
        self.trna_shift_boxplot_log2 = tk.StringVar(value="yes" if bool(S.get("trna_shift_boxplot_log2", True)) else "no")
        self.trna_modifications_boxplot_log2 = tk.StringVar(value="yes" if bool(S.get("trna_modifications_boxplot_log2", True)) else "no")
        self.trna_wobble_boxplot_ymin = tk.StringVar(value=str(S.get("trna_wobble_boxplot_ymin", "") or ""))
        self.trna_wobble_boxplot_ymax = tk.StringVar(value=str(S.get("trna_wobble_boxplot_ymax", "") or ""))
        self.trna_shift_boxplot_ymin = tk.StringVar(value=str(S.get("trna_shift_boxplot_ymin", "") or ""))
        self.trna_shift_boxplot_ymax = tk.StringVar(value=str(S.get("trna_shift_boxplot_ymax", "") or ""))
        self.trna_modifications_boxplot_ymin = tk.StringVar(value=str(S.get("trna_modifications_boxplot_ymin", "") or ""))
        self.trna_modifications_boxplot_ymax = tk.StringVar(value=str(S.get("trna_modifications_boxplot_ymax", "") or ""))
        self.trna_wobble_exclude_outliers = tk.StringVar(value="yes" if bool(S.get("trna_wobble_exclude_outliers", False)) else "no")
        self.trna_shift_exclude_outliers = tk.StringVar(value="yes" if bool(S.get("trna_shift_exclude_outliers", False)) else "no")
        self.trna_modifications_exclude_outliers = tk.StringVar(value="yes" if bool(S.get("trna_modifications_exclude_outliers", False)) else "no")
        self.trna_wobble_outlier_sd = tk.StringVar(value=str(S.get("trna_wobble_outlier_sd", 3.0)))
        self.trna_shift_outlier_sd = tk.StringVar(value=str(S.get("trna_shift_outlier_sd", 3.0)))
        self.trna_modifications_outlier_sd = tk.StringVar(value=str(S.get("trna_modifications_outlier_sd", 3.0)))
        self.trna_modifications_feature_mode = tk.StringVar(value=str(S.get("trna_modifications_feature_mode", "modifications") or "modifications"))
        self.trna_modifications_selection = _canonicalize_gui_plot6_modification_selection(S.get("trna_modifications_selected_features", None))  # None -> default feature set
        self.trna_modifications_selection_status = tk.StringVar(value="All available modifications/enzymes selected")
        self.trna_modification_aas_selection = _split_selection_text(S.get("trna_modifications_include_aas", None))  # None -> default manuscript amino-acid set
        self.trna_modification_aas_status = tk.StringVar(value="Default: Ala, Arg, Asn, Asp, Cys, Gln, Glu, Gly, His, Ile, Leu, Lys, Phe, Pro, Ser, Thr, Tyr, Val")
        self.trna_shift_heatmap_bracket_type = tk.StringVar(value=str(S.get("trna_shift_heatmap_bracket_type", "brace")))
        self.trna_wobble_heatmap_bracket_type = tk.StringVar(value=str(S.get("trna_wobble_heatmap_bracket_type", "brace")))
        self.trna_shift_heatmap_bracket_x = tk.StringVar(value=str(S.get("trna_shift_heatmap_bracket_x", -0.24)))
        self.trna_shift_heatmap_label_x = tk.StringVar(value=str(S.get("trna_shift_heatmap_label_x", -0.35)))
        self.trna_wobble_heatmap_bracket_x = tk.StringVar(value=str(S.get("trna_wobble_heatmap_bracket_x", -0.17)))
        self.trna_wobble_heatmap_label_x = tk.StringVar(value=str(S.get("trna_wobble_heatmap_label_x", -0.27)))

        # Customization for gene-ordered supplementary heatmaps 1 and 2
        self.trna_supp_heatmaps_customize = tk.BooleanVar(value=bool(S.get("trna_supp_heatmaps_customize", False)))
        self.trna_supp_heatmaps_dpi = tk.StringVar(value=str(S.get("trna_supp_heatmaps_dpi", TRNA_SUPP_DETAILS_DEFAULT_DPI) or TRNA_SUPP_DETAILS_DEFAULT_DPI))
        self.trna_supp_heatmaps_fig_width = tk.StringVar(value=str(S.get("trna_supp_heatmaps_fig_width", TRNA_SUPP_DETAILS_DEFAULT_FIG_WIDTH) or TRNA_SUPP_DETAILS_DEFAULT_FIG_WIDTH))
        self.trna_supp_heatmaps_fig_height = tk.StringVar(value=str(S.get("trna_supp_heatmaps_fig_height", TRNA_SUPP_DETAILS_DEFAULT_FIG_HEIGHT) or TRNA_SUPP_DETAILS_DEFAULT_FIG_HEIGHT))
        self.trna_supp_heatmaps_cell_height = tk.StringVar(value=str(S.get("trna_supp_heatmaps_cell_height", TRNA_SUPP_DETAILS_DEFAULT_CELL_HEIGHT) or TRNA_SUPP_DETAILS_DEFAULT_CELL_HEIGHT))
        self.trna_supp_heatmaps_xtick_every_genes = tk.StringVar(value=str(S.get("trna_supp_heatmaps_xtick_every_genes", TRNA_SUPP_DETAILS_DEFAULT_XTICK_EVERY) or TRNA_SUPP_DETAILS_DEFAULT_XTICK_EVERY))
        self.trna_supp_heatmaps_ytick_fontsize = tk.StringVar(value=str(S.get("trna_supp_heatmaps_ytick_fontsize", TRNA_SUPP_DETAILS_DEFAULT_YTICK_FONTSIZE) or TRNA_SUPP_DETAILS_DEFAULT_YTICK_FONTSIZE))
        self.trna_supp_heatmaps_title_fontsize = tk.StringVar(value=str(S.get("trna_supp_heatmaps_title_fontsize", TRNA_SUPP_DETAILS_DEFAULT_TITLE_FONTSIZE) or TRNA_SUPP_DETAILS_DEFAULT_TITLE_FONTSIZE))
        self.trna_supp_heatmaps_xmin = tk.StringVar(value=str(S.get("trna_supp_heatmaps_xmin", "") or ""))
        self.trna_supp_heatmaps_xmax = tk.StringVar(value=str(S.get("trna_supp_heatmaps_xmax", "") or ""))
        self.trna_supp_heatmaps_ymin = tk.StringVar(value=str(S.get("trna_supp_heatmaps_ymin", "") or ""))
        self.trna_supp_heatmaps_ymax = tk.StringVar(value=str(S.get("trna_supp_heatmaps_ymax", "") or ""))

        # Customization for cluster-level enrichment heatmaps 3 and 4
        self.trna_shift_heatmaps_customize = tk.BooleanVar(value=bool(S.get("trna_shift_heatmaps_customize", False)))
        self.trna_shift_heatmaps_dpi = tk.StringVar(value=str(S.get("trna_shift_heatmaps_dpi", TRNA_SHIFT_DETAILS_DEFAULT_DPI) or TRNA_SHIFT_DETAILS_DEFAULT_DPI))
        self.trna_shift_heatmaps_fig_width = tk.StringVar(value=str(S.get("trna_shift_heatmaps_fig_width", TRNA_SHIFT_DETAILS_DEFAULT_FIG_WIDTH) or TRNA_SHIFT_DETAILS_DEFAULT_FIG_WIDTH))
        self.trna_shift_heatmaps_fig_height = tk.StringVar(value=str(S.get("trna_shift_heatmaps_fig_height", TRNA_SHIFT_DETAILS_DEFAULT_FIG_HEIGHT) or TRNA_SHIFT_DETAILS_DEFAULT_FIG_HEIGHT))
        self.trna_shift_heatmaps_cell_width = tk.StringVar(value=str(S.get("trna_shift_heatmaps_cell_width", TRNA_SHIFT_DETAILS_DEFAULT_CELL_WIDTH) or TRNA_SHIFT_DETAILS_DEFAULT_CELL_WIDTH))
        self.trna_shift_heatmaps_cell_height = tk.StringVar(value=str(S.get("trna_shift_heatmaps_cell_height", TRNA_SHIFT_DETAILS_DEFAULT_CELL_HEIGHT) or TRNA_SHIFT_DETAILS_DEFAULT_CELL_HEIGHT))
        self.trna_shift_heatmaps_xtick_fontsize = tk.StringVar(value=str(S.get("trna_shift_heatmaps_xtick_fontsize", TRNA_SHIFT_DETAILS_DEFAULT_XTICK_FONTSIZE) or TRNA_SHIFT_DETAILS_DEFAULT_XTICK_FONTSIZE))
        self.trna_shift_heatmaps_ytick_fontsize = tk.StringVar(value=str(S.get("trna_shift_heatmaps_ytick_fontsize", TRNA_SHIFT_DETAILS_DEFAULT_YTICK_FONTSIZE) or TRNA_SHIFT_DETAILS_DEFAULT_YTICK_FONTSIZE))
        self.trna_shift_heatmaps_title_fontsize = tk.StringVar(value=str(S.get("trna_shift_heatmaps_title_fontsize", TRNA_SHIFT_DETAILS_DEFAULT_TITLE_FONTSIZE) or TRNA_SHIFT_DETAILS_DEFAULT_TITLE_FONTSIZE))
        self.trna_shift_heatmaps_xmin = tk.StringVar(value=str(S.get("trna_shift_heatmaps_xmin", "") or ""))
        self.trna_shift_heatmaps_xmax = tk.StringVar(value=str(S.get("trna_shift_heatmaps_xmax", "") or ""))
        self.trna_shift_heatmaps_ymin = tk.StringVar(value=str(S.get("trna_shift_heatmaps_ymin", "") or ""))
        self.trna_shift_heatmaps_ymax = tk.StringVar(value=str(S.get("trna_shift_heatmaps_ymax", "") or ""))

        # Dedicated customization for Plot 6 / tRNA modification boxplots.
        self.trna_modification_plots_customize = tk.BooleanVar(value=bool(S.get("trna_modification_plots_customize", True)))
        self.trna_modification_plots_dpi = tk.StringVar(value=str(S.get("trna_modification_plots_dpi", PLOT6_DEFAULT_DPI) or PLOT6_DEFAULT_DPI))
        self.trna_modification_plots_fig_width = tk.StringVar(value=str(S.get("trna_modification_plots_fig_width", PLOT6_DEFAULT_FIG_WIDTH) or PLOT6_DEFAULT_FIG_WIDTH))
        self.trna_modification_plots_fig_height = tk.StringVar(value=str(S.get("trna_modification_plots_fig_height", PLOT6_DEFAULT_FIG_HEIGHT) or PLOT6_DEFAULT_FIG_HEIGHT))
        self.trna_modification_plots_caption_size = tk.StringVar(value=str(S.get("trna_modification_plots_caption_size", S.get("trna_modifications_boxplot_caption_size", 17)) or 17))
        self.trna_modification_plots_ymin = tk.StringVar(value=str(S.get("trna_modification_plots_ymin", PLOT6_DEFAULT_YMIN) if S.get("trna_modification_plots_ymin", PLOT6_DEFAULT_YMIN) is not None else PLOT6_DEFAULT_YMIN))
        self.trna_modification_plots_ymax = tk.StringVar(value=str(S.get("trna_modification_plots_ymax", PLOT6_DEFAULT_YMAX) if S.get("trna_modification_plots_ymax", PLOT6_DEFAULT_YMAX) is not None else PLOT6_DEFAULT_YMAX))
        self.trna_modification_plots_group_bar_y = tk.StringVar(value=str(S.get("trna_modification_plots_group_bar_y", PLOT6_DEFAULT_GROUP_BAR_Y) or PLOT6_DEFAULT_GROUP_BAR_Y))
        self.trna_modification_plots_group_label_gap = tk.StringVar(value=str(S.get("trna_modification_plots_group_label_gap", PLOT6_DEFAULT_GROUP_LABEL_GAP) or PLOT6_DEFAULT_GROUP_LABEL_GAP))
        self.trna_modification_plots_star_offset = tk.StringVar(value=str(S.get("trna_modification_plots_star_offset", PLOT6_DEFAULT_STAR_OFFSET) or PLOT6_DEFAULT_STAR_OFFSET))
        self.trna_modification_plots_legend_ncol = tk.StringVar(value=str(S.get("trna_modification_plots_legend_ncol", "") or ""))
        self.trna_modification_plots_box_width = tk.StringVar(value=str(S.get("trna_modification_plots_box_width", 0.18) or 0.18))


    # ------------------------- gene-cluster inference state -------------------------
    def _init_keyword_group_vars(self):
        """Initialize editable keyword groups from Clustering_Pipeline.py."""
        self.keyword_group_vars = {}
        source = getattr(CP, "BASIC_CLUSTER_KEYWORD_GROUPS", {}) or {}
        for i, (name, keywords) in enumerate(source.items()):
            key = f"group_{i:03d}"
            kws = "; ".join([str(k) for k in list(keywords or [])])
            self.keyword_group_vars[key] = {
                "original_name": str(name),
                "enabled": tk.BooleanVar(value=True),
                "name": tk.StringVar(value=str(name)),
                "keywords": tk.StringVar(value=kws),
                "frame": None,
                "keywords_frame": None,
            }

    def _collect_keyword_groups(self):
        """Return the currently enabled GUI keyword groups as a plain dict."""
        groups = {}
        used_names = set()
        for payload in self.keyword_group_vars.values():
            if not bool(payload["enabled"].get()):
                continue
            name = " ".join(str(payload["name"].get() or "").strip().split())
            if not name:
                continue
            keywords = _split_keywords_text(payload["keywords"].get())
            if not keywords:
                continue
            base = name
            suffix = 2
            while name.lower() in used_names:
                name = f"{base} ({suffix})"
                suffix += 1
            used_names.add(name.lower())
            groups[name] = keywords
        return groups

    def _current_inferred_cluster_names(self):
        names = list(self._collect_keyword_groups().keys())
        if bool(self.add_custom_cds_enable.get()):
            custom_name = (self.custom_cds_cluster_name.get().strip() or "custom")
            if custom_name.lower() not in {n.lower() for n in names}:
                names.append(custom_name)
        return names

    def _mark_active_clusters_dirty(self):
        # Keep the user selection if possible, but update the status text to reflect
        # that available cluster names may have changed.
        self._update_active_clusters_status()
        try:
            self._update_extract_clusters_status()
        except Exception:
            pass

    def _update_active_clusters_status(self):
        """Update cluster-selection labels without reading Excel files.

        Older builds opened the cluster workbook and the latest clustering workbook
        repeatedly during startup just to count columns for these status labels.
        On large Dropbox folders or cloud-synced Excel files this could make the
        GUI take minutes to appear. Cluster names are now loaded only when the
        user opens a cluster picker or launches a plot.
        """
        if self.active_clusters_selection is None:
            self.active_clusters_status.set("All clusters selected")
        else:
            self.active_clusters_status.set(f"{len(self.active_clusters_selection)} active cluster(s) selected")
        self._update_decoding_clusters_status()
        self._update_figure_clusters_status()

    def _update_figure_clusters_status(self):
        status_var = getattr(self, "figure_clusters_status", None)
        if status_var is None:
            return
        if self.figure_clusters_selection is None:
            status_var.set("All figure clusters selected")
        else:
            status_var.set(f"{len(self.figure_clusters_selection)} figure cluster(s) selected")

    def _update_decoding_clusters_status(self):
        status_var = getattr(self, "decoding_clusters_status", None)
        if status_var is None:
            return
        if self.decoding_clusters_selection is None:
            status_var.set("All clusters selected for decoding")
        elif len(self.decoding_clusters_selection) == 0:
            status_var.set("No decoding clusters selected")
        else:
            status_var.set(f"{len(self.decoding_clusters_selection)} decoding cluster(s) selected")
        ref_var = getattr(self, "decoding_reference_cluster", None)
        ref_status = getattr(self, "decoding_reference_cluster_status", None)
        if ref_var is not None and ref_status is not None:
            ref = str(ref_var.get() or "").strip()
            if ref:
                ref_status.set(f"Reference cluster: {ref}")
            else:
                ref_status.set("No reference cluster selected")

        self._update_trna_modification_selection_status()
        self._update_trna_modification_aas_status()

    def _update_trna_modification_selection_status(self):
        status = getattr(self, "trna_modifications_selection_status", None)
        if status is None:
            return
        selected = getattr(self, "trna_modifications_selection", None)
        mode = str(getattr(self, "trna_modifications_feature_mode", tk.StringVar(value="modifications")).get() or "modifications").strip().lower()
        label = "enzymes" if mode in {"enzyme", "enzymes", "trme", "trmes"} else "modifications"
        if selected is None:
            if label == "modifications":
                status.set("Default: all available modifications except ac4C34")
            else:
                status.set(f"All available {label} selected")
        elif len(selected) == 0:
            status.set(f"No {label} selected")
        else:
            status.set(f"{len(selected)} {label} selected")

    def _update_trna_modification_aas_status(self):
        status = getattr(self, "trna_modification_aas_status", None)
        if status is None:
            return
        selected = getattr(self, "trna_modification_aas_selection", None)
        if selected is None:
            status.set("Default: " + ", ".join(PLOT6_DEFAULT_AAS))
        elif len(selected) == 0:
            status.set("No amino acid selected")
        else:
            status.set(f"{len(selected)} amino acid(s) considered: {', '.join(selected[:8])}" + ("..." if len(selected) > 8 else ""))

    def _get_available_cluster_names(self, show_errors=True):
        """Return cluster names visible to the Active clusters selector."""
        mode = self.user_cluster_mode.get().strip()
        internal = _cluster_mode_internal(mode)
        names = []
        if internal == "refined":
            path = self.refined_cluster_file.get().strip()
            if path and os.path.isfile(path):
                try:
                    df = _read_cluster_file(path, sheet_name=self.refined_cluster_sheet.get().strip())
                    names = [str(c).strip() for c in df.columns if str(c).strip()]
                except Exception as e:
                    if show_errors:
                        messagebox.showwarning("Active clusters", f"Could not read the cluster file:\n{e}")
            elif show_errors:
                messagebox.showwarning("Active clusters", "Please select the user cluster xlsx file first.")
        else:
            names = self._current_inferred_cluster_names()

        if bool(self.add_custom_cds_enable.get()):
            custom_name = (self.custom_cds_cluster_name.get().strip() or "custom")
            if custom_name.lower() not in {str(n).strip().lower(): n for n in names}:
                names.append(custom_name)

        for cfg in self._collect_fasta_metric_cluster_configs():
            label = str(cfg.get("label") or "").strip()
            if label and label.lower() not in {str(n).strip().lower(): n for n in names}:
                names.append(label)

        # De-duplicate while preserving order.
        out, seen = [], set()
        for n in names:
            key = str(n).strip().lower()
            if not key or key in seen:
                continue
            seen.add(key)
            out.append(str(n).strip())
        return out

    def _get_available_figure_cluster_names(self, show_errors=True):
        """Return clusters available for the Figures-tab picker.

        The Figures tab primarily replots from an existing clustering workbook,
        but GUI-computed FASTA-derived metric clusters may be generated after
        that workbook was created. Therefore, merge workbook columns with the
        current GUI cluster source instead of returning the workbook columns
        immediately. This keeps the Figures picker synchronized with the
        Input/Output and Decoding Strategies pickers.
        """
        names = []
        workbook = ""
        try:
            workbook = self._locate_clustering_workbook()
        except Exception:
            workbook = ""
        if workbook:
            try:
                df = pd.read_excel(workbook, sheet_name=CP.PIPELINE.get("sheet_locus_tags", "Locus Tags"), nrows=1, dtype=str)
                names.extend([str(c).strip() for c in list(df.columns)[1:] if str(c).strip()])
            except Exception as e:
                if show_errors:
                    messagebox.showwarning("Figure cluster picker", f"Could not read cluster names from the existing clustering workbook:\n{e}")

        # Add clusters that exist in the current GUI state but are not yet
        # present in the last workbook, especially newly computed FASTA-derived
        # metric groups. The plotting path below can create a temporary augmented
        # workbook so these clusters are usable immediately.
        try:
            names.extend(self._get_available_cluster_names(show_errors=False))
        except Exception:
            pass

        out, seen = [], set()
        for n in names:
            key = str(n).strip().lower()
            if not key or key in seen:
                continue
            seen.add(key)
            out.append(str(n).strip())
        return out

    def _active_cluster_selection_for_available(self, available):
        available = list(available or [])
        if self.active_clusters_selection is None:
            return available
        selected_lower = {str(c).strip().lower() for c in self.active_clusters_selection}
        return [c for c in available if str(c).strip().lower() in selected_lower]

    def _default_decoding_cluster_selection_for_available(self, available):
        """Default Plot 4-6 cluster selection for the Decoding strategies tab.

        When the standard manuscript clusters are present, start from only
        Ribosomal proteins and Virulence. If either cluster is absent, keep the
        historical empty default rather than guessing another cluster.
        """
        available = [str(c).strip() for c in list(available or []) if str(c).strip()]
        lower_to_actual = {str(c).strip().lower(): str(c) for c in available}
        rib = lower_to_actual.get("ribosomal proteins")
        vir = lower_to_actual.get("virulence")
        if rib and vir:
            return [rib, vir]
        return []

    def _ensure_default_decoding_preselection_for_available(self, available):
        available = [str(c).strip() for c in list(available or []) if str(c).strip()]
        if not available:
            return
        if not bool(getattr(self, "decoding_clusters_selection_user_set", False)):
            default = self._default_decoding_cluster_selection_for_available(available)
            if default:
                self.decoding_clusters_selection = default
                if not str(self.decoding_reference_cluster.get() or "").strip():
                    self.decoding_reference_cluster.set(default[0])

    def _decoding_cluster_selection_for_available(self, available):
        available = list(available or [])
        self._ensure_default_decoding_preselection_for_available(available)
        if self.decoding_clusters_selection is None:
            return available
        selected_lower = {str(c).strip().lower() for c in self.decoding_clusters_selection}
        return [c for c in available if str(c).strip().lower() in selected_lower]

    def _figure_cluster_selection_for_available(self, available):
        available = list(available or [])
        if self.figure_clusters_selection is None:
            return available
        selected_lower = {str(c).strip().lower() for c in self.figure_clusters_selection}
        return [c for c in available if str(c).strip().lower() in selected_lower]


    # ------------------------- UI -------------------------
    def _build_ui(self):
        style = ttk.Style()
        try:
            style.theme_use("clam")
        except Exception:
            pass

        outer = ttk.Frame(self.root, padding=10)
        outer.pack(fill="both", expand=True)

        title = ttk.Label(outer, text="CodonPipe — Desktop Interface", font=("Arial", 17, "bold"))
        title.pack(anchor="w", pady=(0, 8))
        self._tip(title, "Desktop launcher for the CodonPipe clustering and plotting workflow.")

        main_pane = ttk.Panedwindow(outer, orient="horizontal")
        main_pane.pack(fill="both", expand=True)

        left = ttk.Frame(main_pane)
        right = ttk.Frame(main_pane)
        main_pane.add(left, weight=4)
        main_pane.add(right, weight=2)

        simple_row = ttk.Frame(left)
        simple_row.pack(fill="x", pady=(0, 6))
        simple_chk = ttk.Checkbutton(
            simple_row,
            text="Simplified interface",
            variable=self.simplified_interface,
            command=self._refresh_simplified_interface,
        )
        simple_chk.pack(side="left")
        self._tip(simple_chk, "When enabled, only the most commonly used tabs are shown. Advanced decoding, keyword, clustering-parameter, and figure-detail tabs are hidden unless a decoding table is loaded.")

        notebook = ttk.Notebook(left)
        notebook.pack(fill="both", expand=True)
        self.notebook = notebook
        self._advanced_tabs = []

        self.tab_input = ScrollableFrame(notebook)
        self.tab_core = ScrollableFrame(notebook)
        self.tab_figures = ScrollableFrame(notebook)
        self.tab_codon_usage_analyses = ScrollableFrame(notebook)
        self.tab_trna = ScrollableFrame(notebook)
        self.tab_gene_clusters_inference = ScrollableFrame(notebook)
        self.tab_analysis_settings = ScrollableFrame(notebook)
        self.tab_plot_details = ScrollableFrame(notebook)
        notebook.add(self.tab_input, text="Input/Output")
        notebook.add(self.tab_core, text="Codon usage clustering")
        notebook.add(self.tab_figures, text="Figures")
        notebook.add(self.tab_codon_usage_analyses, text="Codon usage analyses")
        notebook.add(self.tab_trna, text="Decoding strategies")
        notebook.add(self.tab_gene_clusters_inference, text="Gene clusters keywords")
        notebook.add(self.tab_analysis_settings, text="Clustering analysis parameters")
        notebook.add(self.tab_plot_details, text="Figure details")
        self._advanced_tabs = [
            (self.tab_trna, "Decoding strategies"),
            (self.tab_gene_clusters_inference, "Gene clusters keywords"),
            (self.tab_analysis_settings, "Clustering analysis parameters"),
            (self.tab_plot_details, "Figure details"),
        ]

        self._build_input_tab(self.tab_input.inner)
        self._build_core_tab(self.tab_core.inner)
        self._build_figures_tab(self.tab_figures.inner)
        self._build_codon_usage_analyses_tab(self.tab_codon_usage_analyses.inner)
        self._build_trna_tab(self.tab_trna.inner)
        self._build_gene_clusters_inference_tab(self.tab_gene_clusters_inference.inner)
        self._build_analysis_settings_tab(self.tab_analysis_settings.inner)
        self._build_plots_details_tab(self.tab_plot_details.inner)
        self._build_right_panel(right)
        self._refresh_simplified_interface()
        self._update_decoding_status()

    def _tip(self, widget, text):
        if widget is None or not text:
            return
        try:
            tip = ToolTip(widget, text, self.tooltips_enabled)
            self.tooltip_refs.append(tip)
        except Exception:
            pass

    # ------------------------- tabs -------------------------
    def _build_input_tab(self, parent):
        files_box = ttk.LabelFrame(parent, text="Files and dataset context", padding=10)
        files_box.pack(fill="x", pady=5)
        files_box.columnconfigure(1, weight=1)

        self._add_file_row(
            files_box, 0, "Default export root", self.default_root, ask_dir=True,
            tip="Main output folder where CodonPipe will write its results."
        )

        genome_box = ttk.LabelFrame(parent, text="Upload genome", padding=10)
        genome_box.pack(fill="x", pady=6)
        genome_box.columnconfigure(1, weight=1)
        self.genome_box = genome_box

        self.preloaded_genome_label = ttk.Label(genome_box, text="Preloaded genome")
        self.preloaded_genome_label.grid(row=0, column=0, sticky="w", padx=(0, 8), pady=4)
        self._tip(self.preloaded_genome_label, "Choose one of the CDS FASTA genomes bundled with CodonPipe.")
        self.preloaded_genome_combo = ttk.Combobox(
            genome_box,
            textvariable=self.preloaded_genome_choice,
            state="readonly",
            values=[g["label"] for g in self.preloaded_genomes] or ["No preloaded genome found"],
            width=48,
        )
        self.preloaded_genome_combo.grid(row=0, column=1, sticky="ew", pady=4)
        self.preloaded_genome_combo.bind("<<ComboboxSelected>>", lambda _e: self._apply_preloaded_genome_choice())
        self._tip(self.preloaded_genome_combo, "The organism name and FASTA path are inferred from this selection.")

        self.genome_not_available_check = ttk.Checkbutton(
            genome_box,
            text="Genome not available / use my own FASTA",
            variable=self.genome_not_available,
            command=self._refresh_genome_source_widgets,
        )
        self.genome_not_available_check.grid(row=1, column=0, columnspan=3, sticky="w", pady=(6, 2))
        self._tip(self.genome_not_available_check, "Check this when the desired genome is not bundled with CodonPipe; manual FASTA and organism fields will appear.")

        self.manual_genome_frame = ttk.Frame(genome_box)
        self.manual_genome_frame.grid(row=2, column=0, columnspan=3, sticky="ew", pady=(4, 0))
        self.manual_genome_frame.columnconfigure(1, weight=1)
        self._add_file_row(
            self.manual_genome_frame, 0, "Input FASTA", self.fasta_path,
            filetypes=[("FASTA files", "*.fasta *.fa *.fna *.ffn *.fas"), ("All files", "*.*")],
            tip="Coding-sequence FASTA file used by the pipeline."
        )
        lbl = ttk.Label(self.manual_genome_frame, text="Organism name")
        lbl.grid(row=1, column=0, sticky="w", padx=(0, 8), pady=4)
        self._tip(lbl, "Short organism or strain name used in outputs.")
        ent = ttk.Entry(self.manual_genome_frame, textvariable=self.organism_name)
        ent.grid(row=1, column=1, sticky="ew", pady=4)
        self._tip(ent, "Short organism or strain name used in outputs.")

        note = ttk.Label(
            genome_box,
            text="To add more bundled genomes, place CDS FASTA files in a folder named 'Preloaded genomes' next to CodonPipe_GUI.py.",
            foreground="#444444",
            wraplength=760,
        )
        note.grid(row=3, column=0, columnspan=3, sticky="w", pady=(8, 0))

        cluster_box = ttk.LabelFrame(parent, text="Gene clusters", padding=10)
        cluster_box.pack(fill="x", pady=6)
        cluster_box.columnconfigure(1, weight=1)
        self.cluster_box = cluster_box

        lbl = ttk.Label(cluster_box, text="Gene clusters")
        lbl.grid(row=0, column=0, sticky="w", padx=(0, 8), pady=4)
        self._tip(lbl, "Choose how gene clusters are obtained for downstream visualizations and tests.")
        combo = ttk.Combobox(cluster_box, textvariable=self.user_cluster_mode, state="readonly", values=USER_CLUSTER_CHOICES, width=48)
        combo.grid(row=0, column=1, sticky="w", pady=4)
        combo.bind("<<ComboboxSelected>>", lambda _e: self._on_cluster_mode_changed())
        self._tip(combo, "Clusters can be inferred from FASTA annotations, inferred from a DAVID gene2terms TXT file, or loaded from a user xlsx file with one column per cluster.")

        self.david_terms_label = ttk.Label(cluster_box, text="DAVID gene2terms txt")
        self.david_terms_entry = ttk.Entry(cluster_box, textvariable=self.david_gene2terms_path)
        self.david_terms_button = ttk.Button(cluster_box, text="Browse", command=self._browse_david_gene2terms_file)
        self._tip(self.david_terms_entry, "Optional explicit DAVID gene2terms TXT file. If left blank, CodonPipe will look for 'DAVID gene2terms.txt' in the output/genome folders or use a newly generated one if the DAVID scan is enabled.")
        self._tip(self.david_terms_button, "Browse for a DAVID gene2terms TXT/TSV/CSV file.")

        self.refined_label = ttk.Label(cluster_box, text="Cluster xlsx file")
        self.refined_entry = ttk.Entry(cluster_box, textvariable=self.refined_cluster_file)
        self.refined_button = ttk.Button(cluster_box, text="Browse", command=self._browse_refined_cluster)
        self.refined_sheet_label = ttk.Label(cluster_box, text="Excel sheet (optional)")
        self.refined_sheet_entry = ttk.Entry(cluster_box, textvariable=self.refined_cluster_sheet)
        self._tip(self.refined_button, "Browse for a user-provided xlsx cluster file with one column per cluster.")
        self._tip(self.refined_sheet_entry, "Optional sheet name if the cluster workbook contains multiple sheets.")

        note_clusters = ttk.Label(
            cluster_box,
            text=(
                "If a cluster file is added to the 'Preloaded genomes' subfolder, it will be loaded by default. "
                "Name should be 'species name clusters.xlsx'."
            ),
            foreground="#444444",
            wraplength=760,
        )
        note_clusters.grid(row=3, column=0, columnspan=3, sticky="w", pady=(8, 0))

        self.active_clusters_button = ttk.Button(cluster_box, text="Active clusters", command=self._open_active_clusters_dialog)
        self.active_clusters_button.grid(row=5, column=0, sticky="w", padx=(0, 8), pady=(10, 4))
        self._tip(self.active_clusters_button, "Open the current cluster list and choose which clusters are active for plots and 2D KS analyses. By default, all available clusters are active.")
        self.active_clusters_label = ttk.Label(cluster_box, textvariable=self.active_clusters_status, foreground="#444444")
        self.active_clusters_label.grid(row=5, column=1, columnspan=2, sticky="w", pady=(10, 4))


        decoding_box = ttk.LabelFrame(parent, text="Decoding strategies", padding=10)
        decoding_box.pack(fill="x", pady=6)
        decoding_box.columnconfigure(1, weight=1)
        self.decoding_box = decoding_box

        info = ttk.Label(
            decoding_box,
            text=(
                "To use this analysis, either load a decoding table or place the Excel file in the 'Preloaded genomes' subfolder. "
                "Name should be 'species name decoding table.xlsx'."
            ),
            foreground="#444444",
            wraplength=760,
        )
        info.grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, 8))

        lbl = ttk.Label(decoding_box, text="Decoding table")
        lbl.grid(row=1, column=0, sticky="w", padx=(0, 8), pady=4)
        self._tip(lbl, "Excel workbook containing the unified decoding-strategy table: AA, codon, anticodon, optional tRNA abundance mean/STD, optional modifications and tRMEs.")
        ent = ttk.Entry(decoding_box, textvariable=self.trna_decoding_table_path)
        ent.grid(row=1, column=1, sticky="ew", pady=4)
        self._tip(ent, "Select the Excel workbook used for decoding strategy analyses.")
        btn = ttk.Button(decoding_box, text="Browse", command=self._browse_trna_decoding_table)
        btn.grid(row=1, column=2, sticky="w", padx=(8, 0), pady=4)
        self._tip(btn, "Choose the Excel workbook used for decoding strategy analyses.")

        lbl_sheet = ttk.Label(decoding_box, text="Sheet (optional)")
        lbl_sheet.grid(row=2, column=0, sticky="w", padx=(0, 8), pady=4)
        self._tip(lbl_sheet, "Optional sheet name. Leave blank to use the first workbook sheet.")
        ent_sheet = ttk.Entry(decoding_box, textvariable=self.trna_decoding_table_sheet)
        ent_sheet.grid(row=2, column=1, sticky="ew", pady=4)
        self._tip(ent_sheet, "Optional sheet name for the unified decoding table.")
        self.decoding_status_label = ttk.Label(decoding_box, text="", foreground="#444444")
        self.decoding_status_label.grid(row=3, column=0, columnspan=3, sticky="w", pady=(4, 0))

        custom_cds_box = ttk.LabelFrame(parent, text="Custom CDS", padding=10)
        custom_cds_box.pack(fill="x", pady=6)
        custom_cds_box.columnconfigure(1, weight=1)
        self.custom_cds_box = custom_cds_box

        self.add_custom_cds_check = ttk.Checkbutton(
            custom_cds_box,
            text="Add custom CDS to analysis",
            variable=self.add_custom_cds_enable,
            command=lambda: (self._refresh_optional_sections(), self._mark_active_clusters_dirty()),
        )
        self.add_custom_cds_check.grid(row=0, column=0, sticky="w", pady=(0, 4))
        self._tip(self.add_custom_cds_check, "Append one or more CDS FASTA files to the selected genome before codon-usage analysis. Each selected file may contain one or several CDS records.")

        cluster_name_frame = ttk.Frame(custom_cds_box)
        cluster_name_frame.grid(row=0, column=1, columnspan=2, sticky="ew", padx=(14, 0), pady=(0, 4))
        cluster_name_frame.columnconfigure(1, weight=1)
        lbl_cname = ttk.Label(cluster_name_frame, text="Gene cluster name")
        lbl_cname.grid(row=0, column=0, sticky="w", padx=(0, 6))
        self._tip(lbl_cname, "Name used if these custom CDS are included as a gene cluster.")
        ent_cname = ttk.Entry(cluster_name_frame, textvariable=self.custom_cds_cluster_name, width=28)
        ent_cname.grid(row=0, column=1, sticky="ew")
        ent_cname.bind("<KeyRelease>", lambda _e: self._mark_active_clusters_dirty())
        self._tip(ent_cname, "Default: custom. You can rename this to Reporter, Plasmid genes, Recoded CDS, etc.")

        self.custom_cds_content = ttk.Frame(custom_cds_box)
        self.custom_cds_content.grid(row=1, column=0, columnspan=3, sticky="ew")
        self.custom_cds_content.columnconfigure(1, weight=1)

        lbl = ttk.Label(self.custom_cds_content, text="Additional FASTA/.fna file(s)")
        lbl.grid(row=0, column=0, sticky="w", padx=(0, 8), pady=4)
        self._tip(lbl, "Select one or more FASTA/.fna files containing CDS records to append at the end of the selected genome CDS FASTA.")
        ent = ttk.Entry(self.custom_cds_content, textvariable=self.custom_cds_paths)
        ent.grid(row=0, column=1, sticky="ew", pady=4)
        self._tip(ent, "Multiple files are stored as a semicolon-separated list. Use the Browse button to select several files at once.")
        btn = ttk.Button(self.custom_cds_content, text="Browse", command=self._browse_custom_cds_files)
        btn.grid(row=0, column=2, sticky="w", padx=(8, 0), pady=4)
        self._tip(btn, "Select one or more FASTA/.fna files containing custom CDS records.")

        note_custom = ttk.Label(
            self.custom_cds_content,
            text="When custom CDS are added, they are automatically available as a gene cluster. Use Active clusters to include or exclude this group from plots.",
            foreground="#444444",
            wraplength=760,
        )
        note_custom.grid(row=1, column=0, columnspan=3, sticky="w", pady=(4, 0))


        state_box = ttk.LabelFrame(parent, text="Project state", padding=10)
        state_box.pack(fill="x", pady=6)
        state_box.columnconfigure(2, weight=1)
        save_btn = ttk.Button(state_box, text="Save project", command=self._save_project_state_dialog)
        save_btn.grid(row=0, column=0, sticky="w", padx=(0, 8), pady=2)
        self._tip(save_btn, "Save the current GUI choices to a reusable CodonPipe_session.json file.")
        load_btn = ttk.Button(state_box, text="Load project", command=self._load_project_state_dialog)
        load_btn.grid(row=0, column=1, sticky="w", padx=(0, 8), pady=2)
        self._tip(load_btn, "Reload a previously saved CodonPipe session JSON file.")
        ttk.Label(
            state_box,
            text="Session files store GUI choices only; they do not copy genome FASTA files or output folders.",
            foreground="#444444",
            wraplength=760,
        ).grid(row=0, column=2, sticky="w", pady=2)

    def _build_codon_usage_analyses_tab(self, parent):
        parent.columnconfigure(0, weight=1)

        ttk.Label(
            parent,
            text=f"Codon usage analyses — {CODONPIPE_GUI_BUILD}",
            foreground="#555555",
        ).pack(fill="x", pady=(0, 4))

        compare_box = ttk.LabelFrame(parent, text="Per cluster codon usage comparisons", padding=10)
        compare_box.pack(fill="x", pady=6)
        compare_box.columnconfigure(0, weight=1)
        self.codon_compare_box = compare_box

        info = ttk.Label(
            compare_box,
            text=(
                "Load one or more CDS FASTA/.fna files, or load a currently available gene cluster, as a comparison group. "
                "After each group is loaded, a new empty row appears so another group can be added. "
                "Suggested names are inferred automatically but can be edited before plotting. "
                "Tick one reference group for statistics; when a statistical test is selected, all other groups are compared against that reference. "
                "If no reference is selected, statistics compare each group against the genome."
            ),
            foreground="#444444",
            wraplength=980,
        )
        info.grid(row=0, column=0, columnspan=5, sticky="w", pady=(0, 6))

        controls = ttk.Frame(compare_box)
        controls.grid(row=1, column=0, columnspan=5, sticky="ew", pady=(0, 6))
        ttk.Label(controls, text="Metric").pack(side="left", padx=(0, 6))
        metric_combo = ttk.Combobox(
            controls,
            textvariable=self.codon_compare_metric,
            state="readonly",
            values=CODON_COMPARE_METRIC_CHOICES,
            width=28,
        )
        metric_combo.pack(side="left")
        self._tip(metric_combo, "Relative codon usage compares codons within each synonymous amino-acid family. Absolute codon frequency uses codon count divided by total sense codons in each CDS. Amino acid identity collapses synonymous codons and plots the 20 amino-acid frequencies.")
        self.codon_compare_rcu_display_combo = ttk.Combobox(
            controls,
            textvariable=self.codon_compare_rcu_display,
            state="readonly",
            values=CODON_COMPARE_RCU_DISPLAY_CHOICES,
            width=34,
        )
        self.codon_compare_rcu_display_combo.pack(side="left", padx=(8, 0))
        self._tip(self.codon_compare_rcu_display_combo, "Choose whether the raw codon-usage plot shows the genome mean, hides it, or uses genome-normalized RCU z-scores.")
        ttk.Label(controls, text="Statistics").pack(side="left", padx=(16, 6))
        stat_combo = ttk.Combobox(
            controls,
            textvariable=self.codon_compare_statistics,
            state="readonly",
            values=CODON_COMPARE_STAT_CHOICES,
            width=18,
        )
        stat_combo.pack(side="left")
        self._tip(stat_combo, "Statistical comparison for each displayed codon or amino acid. If one row is checked as Ref for statistics, all other groups are compared against that reference; otherwise each group is compared against the genome.")
        ttk.Label(controls, text="Display as").pack(side="left", padx=(16, 6))
        plot_style_combo = ttk.Combobox(
            controls,
            textvariable=self.codon_compare_plot_style,
            state="readonly",
            values=CODON_COMPARE_PLOT_STYLE_CHOICES,
            width=14,
        )
        plot_style_combo.pack(side="left")
        self._tip(plot_style_combo, "Choose whether the raw codon/amino-acid comparison is shown as mean±SD markers, line plot, boxplots, or violin plots.")
        self.codon_compare_metric.trace_add("write", lambda *_: self._refresh_codon_compare_metric_widgets())

        self.codon_compare_rows_frame = ttk.Frame(compare_box)
        self.codon_compare_rows_frame.grid(row=2, column=0, columnspan=5, sticky="ew")
        self.codon_compare_rows_frame.columnconfigure(2, weight=1)
        ttk.Label(self.codon_compare_rows_frame, text="Ref for\nstatistics", justify="center").grid(row=0, column=0, sticky="w", padx=(0, 8), pady=(0, 2))
        ttk.Label(self.codon_compare_rows_frame, text="Name").grid(row=0, column=1, sticky="w", padx=(0, 6), pady=(0, 2))
        ttk.Label(self.codon_compare_rows_frame, text="Gene(s)").grid(row=0, column=2, sticky="w", padx=(0, 6), pady=(0, 2))
        self._add_codon_compare_group_row()

        plot_row = ttk.Frame(compare_box)
        plot_row.grid(row=3, column=0, columnspan=5, sticky="w", pady=(10, 0))
        codon_pick_btn = ttk.Button(plot_row, text="Codons to display", command=self._open_codon_compare_codon_selector)
        self.codon_compare_codon_pick_btn = codon_pick_btn
        codon_pick_btn.pack(side="left")
        self._tip(codon_pick_btn, "Choose which codons are displayed in the raw usage and correlation plots. For Amino acid identity, all 20 amino acids are displayed.")
        highlight_btn = ttk.Button(plot_row, text="Codons highlight", command=self._open_codon_compare_highlight_selector)
        self.codon_compare_highlight_btn = highlight_btn
        highlight_btn.pack(side="left", padx=(8, 0))
        self._tip(highlight_btn, "Choose codons whose x-axis labels should be highlighted in the raw codon-usage plot. Highlighted codon labels are bold and blue; this does not change which codons are displayed.")
        ttk.Label(plot_row, textvariable=self.codon_compare_highlighted_codons_status, foreground="#0057D9").pack(side="left", padx=(6, 14))
        ttk.Label(plot_row, textvariable=self.codon_compare_selected_codons_status, foreground="#444444").pack(side="left", padx=(0, 18))
        plot_btn = ttk.Button(plot_row, text="Plot raw usage", command=self._plot_codon_usage_comparison)
        plot_btn.pack(side="left")
        self._tip(plot_btn, "Generate the genome-vs-group usage plot now. For Amino acid identity, this shows absolute amino-acid frequencies. This does not run the full clustering pipeline.")
        log2_chk = ttk.Checkbutton(plot_row, text="log2 y scale", variable=self.codon_compare_log2_y_scale)
        log2_chk.pack(side="left", padx=(8, 0))
        self._tip(log2_chk, "Use a base-2 logarithmic axis for Plot raw usage. Correlation plots use base-2 logarithmic X and Y axes so the y=x reference line remains meaningful. Non-positive values cannot be shown on log axes and are ignored by Matplotlib.")
        corr_btn = ttk.Button(plot_row, text="Plot correlations", command=self._plot_codon_usage_correlations)
        corr_btn.pack(side="left", padx=(8, 0))
        self._tip(corr_btn, "Generate all pairwise codon-usage or amino-acid frequency correlation plots between the loaded FASTA/cluster groups.")
        self._refresh_codon_compare_metric_widgets()
        self._update_codon_compare_selected_codons_status()

        extract_box = ttk.LabelFrame(parent, text="Extract fasta from clusters", padding=10)
        extract_box.pack(fill="x", pady=8)
        extract_box.columnconfigure(1, weight=1)
        ttk.Label(
            extract_box,
            text=(
                "Select one or more currently available clusters and export their CDS as FASTA files. "
                "This is useful for immediately re-loading those clusters in the comparison plot above."
            ),
            foreground="#444444",
            wraplength=980,
        ).grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, 8))

        pick_btn = ttk.Button(extract_box, text="Cluster picker", command=self._open_extract_clusters_dialog)
        pick_btn.grid(row=1, column=0, sticky="w", padx=(0, 8), pady=4)
        self._tip(pick_btn, "Choose which cluster(s) should be exported as FASTA. By default all available clusters are selected.")
        ttk.Label(extract_box, textvariable=self.extract_clusters_status, foreground="#444444").grid(row=1, column=1, sticky="w", pady=4)

        extract_btn = ttk.Button(extract_box, text="Extract fasta from cluster(s)", command=self._extract_fasta_from_clusters)
        extract_btn.grid(row=2, column=0, sticky="w", pady=(8, 0))
        self._tip(extract_btn, "Create a subfolder containing one FASTA per CDS and one combined FASTA for the selected cluster(s).")

        metric_box = ttk.LabelFrame(parent, text="FASTA-derived metric groups", padding=10)
        metric_box.pack(fill="x", pady=8)
        metric_box.columnconfigure(1, weight=1)
        ttk.Label(
            metric_box,
            text=(
                "Generate extra gene clusters directly from CDS sequence metrics. "
                "Codon-defined percentages are normalized only over the synonymous amino-acid families represented by the selected codons."
            ),
            foreground="#444444",
            wraplength=980,
        ).grid(row=0, column=0, columnspan=4, sticky="w", pady=(0, 8))
        metric_pick_btn = ttk.Button(metric_box, text="Choose metrics", command=self._open_fasta_metric_group_selector)
        metric_pick_btn.grid(row=1, column=0, sticky="w", padx=(0, 8), pady=4)
        self._tip(metric_pick_btn, "Select sequence-derived metric groups such as insensitive genes, tRNA-modification-dependent genes, AT-high genes, or GC-low genes.")
        ttk.Label(metric_box, textvariable=self.fasta_metric_cluster_status, foreground="#444444").grid(row=1, column=1, sticky="w", pady=4)

        self.fasta_metric_rows_frame = ttk.Frame(metric_box)
        self.fasta_metric_rows_frame.grid(row=2, column=0, columnspan=4, sticky="ew", pady=(6, 2))
        self.fasta_metric_rows_frame.columnconfigure(0, weight=1)
        metric_compute_btn = ttk.Button(metric_box, text="Compute clusters", command=self._compute_fasta_metric_clusters_from_gui)
        metric_compute_btn.grid(row=3, column=0, sticky="w", pady=(8, 0))
        self._tip(metric_compute_btn, "Compute FASTA-derived clusters with the current thresholds/top-hit settings and add them to the available cluster list.")
        self._refresh_fasta_metric_rows()

    # ------------------------- codon comparison controls -------------------------
    def _is_codon_compare_aa_metric(self, metric=None):
        metric = self.codon_compare_metric.get() if metric is None else metric
        return "amino" in str(metric or "").strip().lower()

    def _codon_compare_feature_order(self, metric=None):
        return list(_AA_ORDER) if self._is_codon_compare_aa_metric(metric) else list(_CODON_ORDER)

    def _refresh_codon_compare_metric_widgets(self):
        combo = getattr(self, "codon_compare_rcu_display_combo", None)
        metric = str(self.codon_compare_metric.get() or "").strip().lower()
        if combo is not None:
            if metric.startswith("relative"):
                try:
                    combo.pack(side="left", padx=(8, 0))
                except Exception:
                    pass
            else:
                try:
                    combo.pack_forget()
                except Exception:
                    pass

        disabled_for_aa = self._is_codon_compare_aa_metric(metric)
        pick_btn = getattr(self, "codon_compare_codon_pick_btn", None)
        if pick_btn is not None:
            try:
                pick_btn.configure(state=("disabled" if disabled_for_aa else "normal"))
            except Exception:
                pass
        highlight_btn = getattr(self, "codon_compare_highlight_btn", None)
        if highlight_btn is not None:
            try:
                highlight_btn.configure(state=("disabled" if disabled_for_aa else "normal"))
            except Exception:
                pass
        self._update_codon_compare_selected_codons_status()
        self._update_codon_compare_highlight_status()

    def _update_codon_compare_selected_codons_status(self):
        if self._is_codon_compare_aa_metric():
            self.codon_compare_selected_codons_status.set("Amino acid identity: all 20 amino acids displayed")
            return
        selected = [c for c in getattr(self, "codon_compare_selected_codons", []) if c in _CODON_ORDER]
        if len(selected) == len(_CODON_ORDER):
            self.codon_compare_selected_codons_status.set("All codon groups selected")
        elif len(selected) == 0:
            self.codon_compare_selected_codons_status.set("No codon group selected")
        else:
            full_groups, partial_groups = _codon_display_group_summary(selected)
            if partial_groups:
                self.codon_compare_selected_codons_status.set(
                    f"{len(full_groups)} full + {len(partial_groups)} partial groups ({len(selected)}/{len(_CODON_ORDER)} codons)"
                )
            else:
                self.codon_compare_selected_codons_status.set(
                    f"{len(full_groups)}/{len(_CODON_DISPLAY_GROUPS)} codon groups selected ({len(selected)} codons)"
                )

    def _update_codon_compare_highlight_status(self):
        """Update the status text for codon-label highlighting in raw-usage plots."""
        if self._is_codon_compare_aa_metric():
            self.codon_compare_highlighted_codons_status.set("Highlight disabled for amino-acid mode")
            return
        highlighted = [c for c in getattr(self, "codon_compare_highlighted_codons", []) if c in _CODON_ORDER]
        if not highlighted:
            self.codon_compare_highlighted_codons_status.set("No codon highlighted")
        elif len(highlighted) == len(_CODON_ORDER):
            self.codon_compare_highlighted_codons_status.set("All codons highlighted")
        else:
            self.codon_compare_highlighted_codons_status.set(f"{len(highlighted)} codon(s) highlighted")

    def _open_codon_compare_highlight_selector(self):
        """Open a codon selector for highlighting x-axis labels in raw-usage plots."""
        if self._is_codon_compare_aa_metric():
            messagebox.showinfo("Codons highlight", "Codon highlighting applies to codon-usage plots. Amino acid identity plots display amino-acid labels instead.")
            return
        chosen = self._open_codon_selector_dialog(
            "Codons highlight",
            getattr(self, "codon_compare_highlighted_codons", []),
            allow_empty=True,
            intro_text=(
                "Select codons whose x-axis labels should be highlighted in Plot raw usage. "
                "Highlighted codons are shown in bold blue. This does not change which codons are displayed."
            ),
        )
        if chosen is None:
            return
        self.codon_compare_highlighted_codons = [c for c in chosen if c in _CODON_ORDER]
        self._update_codon_compare_highlight_status()

    def _open_codon_selector_dialog(self, title, selected_codons, allow_empty=False, intro_text=None):
        selected = {str(c).upper().replace("U", "T") for c in list(selected_codons or [])}
        win = tk.Toplevel(self.root)
        win.title(title)
        win.geometry("520x720")
        win.transient(self.root)
        win.grab_set()

        outer = ttk.Frame(win, padding=10)
        outer.pack(fill="both", expand=True)
        ttk.Label(
            outer,
            text=(intro_text or "Select codons to include. Codons are grouped by amino acid."),
            wraplength=480,
            foreground="#444444",
        ).pack(anchor="w", pady=(0, 8))

        tools = ttk.Frame(outer)
        tools.pack(fill="x", pady=(0, 6))
        vars_by_codon = {c: tk.BooleanVar(value=(c in selected)) for c in _CODON_ORDER}

        def set_all(value):
            for v in vars_by_codon.values():
                v.set(bool(value))

        ttk.Button(tools, text="Select all", command=lambda: set_all(True)).pack(side="left")
        ttk.Button(tools, text="Select none", command=lambda: set_all(False)).pack(side="left", padx=(6, 0))

        scroller = ScrollableFrame(outer)
        scroller.pack(fill="both", expand=True)
        row_i = 0
        for aa, codons in _AA_TO_CODONS.items():
            lab = ttk.Label(scroller.inner, text=aa, font=("Arial", 9, "bold"))
            lab.grid(row=row_i, column=0, sticky="w", pady=(8 if row_i else 0, 2))
            cframe = ttk.Frame(scroller.inner)
            cframe.grid(row=row_i, column=1, sticky="w", pady=(8 if row_i else 0, 2))
            for j, codon in enumerate(codons):
                ttk.Checkbutton(cframe, text=codon, variable=vars_by_codon[codon]).grid(row=0, column=j, sticky="w", padx=(0, 8))
            row_i += 1

        result = {"selected": None}
        bottom = ttk.Frame(outer)
        bottom.pack(fill="x", pady=(8, 0))

        def apply_selection():
            chosen = [c for c in _CODON_ORDER if vars_by_codon[c].get()]
            if not chosen and not allow_empty:
                messagebox.showwarning(title, "Please select at least one codon.")
                return
            result["selected"] = chosen
            win.destroy()

        ttk.Button(bottom, text="Apply", command=apply_selection).pack(side="left")
        ttk.Button(bottom, text="Cancel", command=win.destroy).pack(side="left", padx=(8, 0))
        self.root.wait_window(win)
        return result["selected"]

    def _open_codon_compare_codon_selector(self):
        """Open a group-level codon selector for codon-usage comparison plots."""
        if self._is_codon_compare_aa_metric():
            messagebox.showinfo("Amino acid identity", "Amino acid identity plots always display the 20 amino acids; the codon selector is not used for this metric.")
            return
        selected = {str(c).upper().replace("U", "T") for c in getattr(self, "codon_compare_selected_codons", list(_CODON_ORDER))}
        win = tk.Toplevel(self.root)
        win.title("Codons to display")
        win.geometry("560x680")
        win.transient(self.root)
        win.grab_set()

        outer = ttk.Frame(win, padding=10)
        outer.pack(fill="both", expand=True)
        ttk.Label(
            outer,
            text=(
                "Select codon groups to display in the raw codon-usage and correlation plots. "
                "Groups are synonymous codon boxes; six-codon amino acids are split into their 4-codon and 2-codon boxes."
            ),
            wraplength=520,
            foreground="#444444",
        ).pack(anchor="w", pady=(0, 8))

        tools = ttk.Frame(outer)
        tools.pack(fill="x", pady=(0, 6))
        vars_by_group = {}
        for name, codons in _CODON_DISPLAY_GROUPS:
            cset = set(codons)
            vars_by_group[name] = tk.BooleanVar(value=bool(cset and cset.issubset(selected)))

        def set_all(value):
            for v in vars_by_group.values():
                v.set(bool(value))

        ttk.Button(tools, text="Select all groups", command=lambda: set_all(True)).pack(side="left")
        ttk.Button(tools, text="Select none", command=lambda: set_all(False)).pack(side="left", padx=(6, 0))

        scroller = ScrollableFrame(outer)
        scroller.pack(fill="both", expand=True)
        for row_i, (name, codons) in enumerate(_CODON_DISPLAY_GROUPS):
            label = f"{name}: " + ", ".join(codons)
            ttk.Checkbutton(scroller.inner, text=label, variable=vars_by_group[name]).grid(
                row=row_i, column=0, sticky="w", padx=(0, 4), pady=3
            )

        result = {"selected": None}
        bottom = ttk.Frame(outer)
        bottom.pack(fill="x", pady=(8, 0))

        def apply_selection():
            chosen = []
            for name, codons in _CODON_DISPLAY_GROUPS:
                if vars_by_group[name].get():
                    chosen.extend(codons)
            chosen = [c for c in _CODON_ORDER if c in set(chosen)]
            if not chosen:
                messagebox.showwarning("Codons to display", "Please select at least one codon group.")
                return
            result["selected"] = chosen
            win.destroy()

        ttk.Button(bottom, text="Apply", command=apply_selection).pack(side="left")
        ttk.Button(bottom, text="Cancel", command=win.destroy).pack(side="left", padx=(8, 0))
        self.root.wait_window(win)
        if result["selected"] is None:
            return
        self.codon_compare_selected_codons = [c for c in result["selected"] if c in _CODON_ORDER]
        self._update_codon_compare_selected_codons_status()

    def _selected_codon_indices_for_compare(self):
        selected = [c for c in getattr(self, "codon_compare_selected_codons", list(_CODON_ORDER)) if c in _CODON_ORDER]
        if not selected:
            selected = list(_CODON_ORDER)
        idx = [i for i, c in enumerate(_CODON_ORDER) if c in set(selected)]
        return selected, np.asarray(idx, dtype=int)

    # ------------------------- FASTA-derived metric groups -------------------------
    def _metric_mode_display(self, mode):
        return "top hits" if str(mode or "").strip().lower().replace("_", " ") == "top hits" else "bins"

    def _metric_mode_internal(self, mode):
        return "top_hits" if str(mode or "").strip().lower().replace("_", " ") == "top hits" else "bins"

    def _set_fasta_metric_group_rows_from_configs(self, configs):
        self.fasta_metric_group_rows = []
        for cfg in list(configs or []):
            key = str(cfg.get("key") or "").strip()
            default = dict(getattr(self, "fasta_metric_default_defs", {}).get(key, {}))
            merged = default
            merged.update(dict(cfg))
            merged["key"] = key
            merged["label"] = str(merged.get("label") or default.get("label") or key).strip()
            merged["metric_type"] = str(merged.get("metric_type") or default.get("metric_type") or "codon_group")
            merged["codons"] = [str(c).upper().replace("U", "T") for c in list(merged.get("codons", default.get("codons", [])) or []) if str(c).strip()]
            row = {
                "key": key,
                "label": merged["label"],
                "metric_type": merged["metric_type"],
                "direction": str(merged.get("direction", default.get("direction", "high")) or "high"),
                "codons": merged["codons"],
                "mode": tk.StringVar(value=self._metric_mode_display(merged.get("mode", default.get("mode", "bins")))),
                "cutoff": tk.StringVar(value=str(merged.get("cutoff", default.get("cutoff", 0.0)))),
                "top_n": tk.StringVar(value=str(merged.get("top_n", default.get("top_n", 100)))),
            }
            self.fasta_metric_group_rows.append(row)
        self._refresh_fasta_metric_rows()
        self._update_fasta_metric_cluster_status()
        self._mark_active_clusters_dirty()

    def _selected_fasta_metric_keys(self):
        return [str(r.get("key") or "") for r in getattr(self, "fasta_metric_group_rows", []) if str(r.get("key") or "")]

    def _open_fasta_metric_group_selector(self):
        defaults = getattr(self, "fasta_metric_default_defs", {}) or {}
        available_keys = list(defaults.keys())
        selected = set(self._selected_fasta_metric_keys())
        win = tk.Toplevel(self.root)
        win.title("Choose FASTA-derived metrics")
        win.geometry("520x520")
        win.transient(self.root)
        win.grab_set()

        outer = ttk.Frame(win, padding=10)
        outer.pack(fill="both", expand=True)
        ttk.Label(
            outer,
            text="Choose which FASTA-derived metric groups should be generated as cluster columns.",
            wraplength=480,
            foreground="#444444",
        ).pack(anchor="w", pady=(0, 8))
        vars_by_key = {k: tk.BooleanVar(value=(k in selected)) for k in available_keys}
        for k in available_keys:
            label = defaults[k].get("label", k)
            ttk.Checkbutton(outer, text=label, variable=vars_by_key[k]).pack(anchor="w", pady=3)

        bottom = ttk.Frame(outer)
        bottom.pack(fill="x", pady=(10, 0))
        result = {"apply": False}

        def apply_selection():
            result["apply"] = True
            win.destroy()

        ttk.Button(bottom, text="Apply", command=apply_selection).pack(side="left")
        ttk.Button(bottom, text="Cancel", command=win.destroy).pack(side="left", padx=(8, 0))
        self.root.wait_window(win)
        if not result["apply"]:
            return

        existing = {str(r.get("key")): r for r in getattr(self, "fasta_metric_group_rows", [])}
        configs = []
        for k in available_keys:
            if not vars_by_key[k].get():
                continue
            if k in existing:
                configs.append(self._fasta_metric_row_to_config(existing[k]))
            else:
                cfg = dict(defaults[k])
                cfg["key"] = k
                configs.append(cfg)
        self._set_fasta_metric_group_rows_from_configs(configs)

    def _fasta_metric_row_to_config(self, row):
        mode = self._metric_mode_internal(row.get("mode", tk.StringVar(value="bins")).get())
        try:
            cutoff = float(str(row.get("cutoff", tk.StringVar(value="0")).get()).strip())
        except Exception:
            cutoff = 0.0
        try:
            top_n = max(1, int(float(str(row.get("top_n", tk.StringVar(value="100")).get()).strip())))
        except Exception:
            top_n = 100
        return {
            "key": str(row.get("key") or ""),
            "label": str(row.get("label") or ""),
            "metric_type": str(row.get("metric_type") or "codon_group"),
            "direction": str(row.get("direction") or "high"),
            "codons": [c for c in list(row.get("codons", [])) if c in _CODON_ORDER],
            "mode": mode,
            "cutoff": cutoff,
            "top_n": top_n,
        }

    def _collect_fasta_metric_cluster_configs(self):
        return [self._fasta_metric_row_to_config(r) for r in getattr(self, "fasta_metric_group_rows", [])]

    def _update_fasta_metric_cluster_status(self):
        configs = self._collect_fasta_metric_cluster_configs() if hasattr(self, "fasta_metric_group_rows") else []
        if not configs:
            self.fasta_metric_cluster_status.set("No FASTA-derived metric group selected")
        elif getattr(self, "fasta_metric_cluster_path", ""):
            self.fasta_metric_cluster_status.set(f"{len(configs)} metric group(s) selected; clusters computed")
        else:
            self.fasta_metric_cluster_status.set(f"{len(configs)} metric group(s) selected; not computed yet")

    def _refresh_fasta_metric_row_entry_state(self, row):
        mode = self._metric_mode_internal(row.get("mode", tk.StringVar(value="bins")).get())
        cutoff_entry = row.get("cutoff_entry")
        top_entry = row.get("top_n_entry")
        try:
            if cutoff_entry is not None:
                cutoff_entry.configure(state="normal" if mode == "bins" else "disabled")
            if top_entry is not None:
                top_entry.configure(state="normal" if mode == "top_hits" else "disabled")
        except Exception:
            pass

    def _refresh_fasta_metric_rows(self):
        frame = getattr(self, "fasta_metric_rows_frame", None)
        if frame is None:
            return
        for child in frame.winfo_children():
            child.destroy()
        rows = getattr(self, "fasta_metric_group_rows", [])
        if not rows:
            ttk.Label(frame, text="No metric selected.", foreground="#666666").grid(row=0, column=0, sticky="w", pady=2)
            self._update_fasta_metric_cluster_status()
            return

        headers = ["Metric", "Codons", "Mode", "Cutoff (%)", "Top hits"]
        for j, txt in enumerate(headers):
            ttk.Label(frame, text=txt, font=("Arial", 9, "bold")).grid(row=0, column=j, sticky="w", padx=(0, 10), pady=(0, 4))

        for i, row in enumerate(rows, start=1):
            ttk.Label(frame, text=str(row.get("label") or row.get("key") or "Metric"), width=28).grid(row=i, column=0, sticky="w", padx=(0, 10), pady=3)
            if str(row.get("metric_type")) == "codon_group":
                codon_status = tk.StringVar(value=f"{len(row.get('codons', []))} codons")
                row["codon_status"] = codon_status
                cbtn = ttk.Button(frame, text="Select codons", command=lambda r=row: self._edit_fasta_metric_codons(r))
                cbtn.grid(row=i, column=1, sticky="w", padx=(0, 6), pady=3)
                ttk.Label(frame, textvariable=codon_status, foreground="#444444").grid(row=i, column=1, sticky="w", padx=(112, 10), pady=3)
                self._tip(cbtn, "Edit which codons define this sequence metric. The denominator will be the corresponding synonymous amino-acid families.")
            else:
                ttk.Label(frame, text="—", foreground="#666666").grid(row=i, column=1, sticky="w", padx=(0, 10), pady=3)

            mode_combo = ttk.Combobox(frame, textvariable=row["mode"], state="readonly", values=["bins", "top hits"], width=10)
            mode_combo.grid(row=i, column=2, sticky="w", padx=(0, 10), pady=3)
            def _on_mode_change(_e=None, r=row):
                self._refresh_fasta_metric_row_entry_state(r)
                self._update_fasta_metric_cluster_status()
            mode_combo.bind("<<ComboboxSelected>>", _on_mode_change)
            cutoff_entry = ttk.Entry(frame, textvariable=row["cutoff"], width=10)
            cutoff_entry.grid(row=i, column=3, sticky="w", padx=(0, 10), pady=3)
            top_entry = ttk.Entry(frame, textvariable=row["top_n"], width=10)
            top_entry.grid(row=i, column=4, sticky="w", padx=(0, 10), pady=3)
            row["cutoff_entry"] = cutoff_entry
            row["top_n_entry"] = top_entry
            self._refresh_fasta_metric_row_entry_state(row)
        self._update_fasta_metric_cluster_status()

    def _edit_fasta_metric_codons(self, row):
        chosen = self._open_codon_selector_dialog(
            f"Codons for {row.get('label', 'metric')}",
            row.get("codons", []),
            allow_empty=False,
        )
        if chosen is None:
            return
        row["codons"] = [c for c in chosen if c in _CODON_ORDER]
        status = row.get("codon_status")
        if status is not None:
            status.set(f"{len(row.get('codons', []))} codons")
        self.fasta_metric_cluster_path = ""
        self._update_fasta_metric_cluster_status()
        self._mark_active_clusters_dirty()

    def _compute_fasta_metric_clusters_from_gui(self):
        try:
            if not self.genome_not_available.get():
                self._apply_preloaded_genome_choice()
            fasta = self.fasta_path.get().strip()
            if not fasta or not os.path.isfile(fasta):
                raise ValueError("Please select a valid genome CDS FASTA first.")
            configs = self._collect_fasta_metric_cluster_configs()
            if not configs:
                raise ValueError("Please choose at least one FASTA-derived metric first.")
            codon_range = _validate_codon_range_text(self.fasta_codon_range.get())
            metric_cluster_df, scores_df = build_fasta_metric_cluster_df(
                fasta_path=fasta,
                metric_configs=configs,
                row_id_mode="locus",
                trim_to_multiple_of_3=bool(CP.SET.get("fasta_trim_to_multiple_of_3", True)),
                organism_mode="prokaryote",
                codon_range=codon_range,
            )
            if metric_cluster_df is None or metric_cluster_df.empty:
                raise ValueError("No FASTA-derived clusters were generated with the current settings.")

            # Export the preview/generated clusters in the same canonical locus-tag
            # namespace that the main pipeline uses for Gene lists per cluster.xlsx.
            locus_index, alias_map, _id_map_df, _missing, _dups = build_locus_index(
                fasta, organism_mode="prokaryote", codon_range=codon_range
            )
            metric_cluster_df, scores_df = CP._canonicalize_fasta_metric_outputs(
                metric_cluster_df, scores_df, alias_map=alias_map, ordered_genes=list(locus_index.keys())
            )

            out_root = self._codon_usage_analysis_output_dir()
            os.makedirs(out_root, exist_ok=True)
            cluster_path = os.path.join(out_root, "FASTA-derived metric groups.xlsx")
            scores_path = os.path.join(out_root, "FASTA-derived metric scores.xlsx")
            with pd.ExcelWriter(cluster_path, engine="xlsxwriter") as writer:
                metric_cluster_df.to_excel(writer, sheet_name="Clusters", index=False)
            if scores_df is not None and not scores_df.empty:
                with pd.ExcelWriter(scores_path, engine="xlsxwriter") as writer:
                    scores_df.to_excel(writer, sheet_name="Metric scores", index=False)
                    pd.DataFrame(configs).to_excel(writer, sheet_name="Metric settings", index=False)
            self.fasta_metric_cluster_df = metric_cluster_df
            self.fasta_metric_cluster_scores_df = scores_df
            self.fasta_metric_cluster_path = cluster_path
            self._update_fasta_metric_cluster_status()
            self.active_clusters_selection = None
            self.figure_clusters_selection = None
            self._update_active_clusters_status()
            self._append_log(
                f"[INFO] Computed FASTA-derived metric clusters: {', '.join(map(str, metric_cluster_df.columns))}\n"
                f"[INFO] Saved clusters: {cluster_path}\n"
                f"[INFO] Saved metric scores: {scores_path}\n"
            )
            messagebox.showinfo("FASTA-derived metric groups", "Clusters computed and added to the available cluster list.")
        except Exception as e:
            messagebox.showerror("FASTA-derived metric groups", str(e))

    def _build_fasta_metric_cluster_df_for_tools(self):
        configs = self._collect_fasta_metric_cluster_configs()
        if not configs:
            return pd.DataFrame()
        fasta = self.fasta_path.get().strip()
        if not fasta or not os.path.isfile(fasta):
            raise FileNotFoundError("Please select a valid genome CDS FASTA first.")
        codon_range = _validate_codon_range_text(self.fasta_codon_range.get())
        metric_cluster_df, scores_df = build_fasta_metric_cluster_df(
            fasta_path=fasta,
            metric_configs=configs,
            row_id_mode="locus",
            trim_to_multiple_of_3=bool(CP.SET.get("fasta_trim_to_multiple_of_3", True)),
            organism_mode="prokaryote",
            codon_range=codon_range,
        )
        try:
            locus_index, alias_map, _id_map_df, _missing, _dups = build_locus_index(
                fasta, organism_mode="prokaryote", codon_range=codon_range
            )
            metric_cluster_df, _scores_df = CP._canonicalize_fasta_metric_outputs(
                metric_cluster_df, scores_df, alias_map=alias_map, ordered_genes=list(locus_index.keys())
            )
        except Exception:
            pass
        return metric_cluster_df.fillna("") if metric_cluster_df is not None else pd.DataFrame()

    def _build_gene_clusters_inference_tab(self, parent):
        info = ttk.Label(
            parent,
            text=(
                "Edit the keyword groups used when gene clusters are inferred from FASTA annotations "
                "or from a DAVID gene2terms TXT file. Unchecked groups are ignored. "
                "Keywords are matched case-insensitively as text fragments. Separate keywords with semicolons or commas."
            ),
            wraplength=1150,
            foreground="#444444",
        )
        info.pack(fill="x", pady=(4, 8))

        controls = ttk.Frame(parent)
        controls.pack(fill="x", pady=(0, 6))
        ttk.Button(controls, text="Select all", command=lambda: self._set_all_keyword_groups(True)).pack(side="left", padx=(0, 6))
        ttk.Button(controls, text="Select none", command=lambda: self._set_all_keyword_groups(False)).pack(side="left")

        self.keyword_groups_box = ttk.LabelFrame(parent, text="Keyword groups for inferred gene clusters", padding=10)
        self.keyword_groups_box.pack(fill="both", expand=True, pady=5)
        self.keyword_groups_box.columnconfigure(0, weight=1)

        for i, (key, payload) in enumerate(self.keyword_group_vars.items()):
            row = ttk.Frame(self.keyword_groups_box)
            row.grid(row=i, column=0, sticky="ew", pady=(0, 6))
            row.columnconfigure(1, weight=1)
            payload["frame"] = row

            chk = ttk.Checkbutton(row, variable=payload["enabled"], command=lambda k=key: self._on_keyword_group_toggle(k))
            chk.grid(row=0, column=0, sticky="nw", padx=(0, 6), pady=2)
            self._tip(chk, "Enable or disable this inferred cluster group.")

            name_entry = ttk.Entry(row, textvariable=payload["name"], width=42)
            name_entry.grid(row=0, column=1, sticky="ew", pady=2)
            self._tip(name_entry, "Cluster name. This will become the column name in the inferred cluster table.")
            name_entry.bind("<KeyRelease>", lambda _e: self._mark_active_clusters_dirty())
            name_entry.bind("<FocusOut>", lambda _e: self._mark_active_clusters_dirty())

            kw_frame = ttk.Frame(row)
            kw_frame.grid(row=1, column=1, sticky="ew", pady=(2, 0))
            kw_frame.columnconfigure(1, weight=1)
            payload["keywords_frame"] = kw_frame
            ttk.Label(kw_frame, text="Keywords").grid(row=0, column=0, sticky="w", padx=(0, 8))
            kw_entry = ttk.Entry(kw_frame, textvariable=payload["keywords"])
            kw_entry.grid(row=0, column=1, sticky="ew")
            self._tip(kw_entry, "Semicolon- or comma-separated keywords used to infer this cluster.")
            kw_entry.bind("<KeyRelease>", lambda _e: self._mark_active_clusters_dirty())
            kw_entry.bind("<FocusOut>", lambda _e: self._mark_active_clusters_dirty())

        self._refresh_keyword_group_visibility()

    def _set_all_keyword_groups(self, value):
        for payload in self.keyword_group_vars.values():
            payload["enabled"].set(bool(value))
        self._refresh_keyword_group_visibility()
        self.active_clusters_selection = None
        self.figure_clusters_selection = None
        self.decoding_clusters_selection = []
        self._update_active_clusters_status()

    def _on_keyword_group_toggle(self, _key=None):
        self._refresh_keyword_group_visibility()
        self.active_clusters_selection = None
        self.figure_clusters_selection = None
        self.decoding_clusters_selection = []
        self._update_active_clusters_status()

    def _refresh_keyword_group_visibility(self):
        for payload in self.keyword_group_vars.values():
            frame = payload.get("keywords_frame")
            if frame is None:
                continue
            if bool(payload["enabled"].get()):
                frame.grid()
            else:
                frame.grid_remove()

    def _build_core_tab(self, parent):
        top = ttk.LabelFrame(parent, text="Codon usage clustering", padding=14)
        top.pack(fill="x", pady=5)
        top.columnconfigure(1, weight=1)

        self._add_labeled_widget(
            top, 0, 0, "Codon usage metric",
            ttk.Combobox(top, textvariable=self.usage_basis, state="readonly", values=USAGE_METRIC_CHOICES, width=28),
            "Choose the feature space used to compare genes before dimensional reduction."
        )
        top.grid_slaves(row=0, column=1)[0].bind("<<ComboboxSelected>>", lambda _e: self._sync_codon_set_from_usage_basis())

        self._add_labeled_widget(
            top, 1, 0, "CDS codon range",
            ttk.Entry(top, textvariable=self.fasta_codon_range, width=28),
            "Region of each CDS used to count codons before clustering. Use 'all' for the full CDS, '1-20' for the first 20 codons, '20-200' for codons 20 to 200, or '20-end'. Coordinates are 1-based and inclusive; ends beyond short genes are clipped to the gene end."
        )

        self._add_labeled_widget(
            top, 2, 0, "Dimensional reduction method",
            ttk.Combobox(top, textvariable=self.dimred_method, state="readonly", values=DIMRED_CHOICES, width=28),
            "Choose the method used to project genes into two dimensions."
        )
        top.grid_slaves(row=2, column=1)[0].bind("<<ComboboxSelected>>", lambda _e: self._refresh_dimred_params())

        self._add_labeled_widget(
            top, 3, 0, "Clustering method",
            ttk.Combobox(top, textvariable=self.cluster_method, state="readonly", values=CLUSTERING_METHOD_CHOICES, width=28),
            "Choose the method used to group genes in the embedded space."
        )
        top.grid_slaves(row=3, column=1)[0].bind("<<ComboboxSelected>>", lambda _e: self._refresh_cluster_params())

        self._add_labeled_widget(
            top, 4, 0, "Statistical significance of clusters",
            ttk.Combobox(top, textvariable=self.statistical_test_method, state="readonly", values=STATISTICAL_TEST_CHOICES, width=28),
            "Enable or disable 2D Kolmogorov-Smirnov comparisons between selected clusters."
        )
        top.grid_slaves(row=4, column=1)[0].bind("<<ComboboxSelected>>", lambda _e: self._sync_analysis_choice_booleans())

        ttk.Label(
            top,
            text="Detailed parameters for the selected methods are shown in the 'Clustering analysis parameters' tab.",
            foreground="#444444",
            wraplength=760,
        ).grid(row=5, column=0, columnspan=3, sticky="w", pady=(10, 0))

        run_row = ttk.Frame(top)
        run_row.grid(row=6, column=0, columnspan=3, sticky="w", pady=(14, 0))
        b1 = ttk.Button(run_row, text="Run clustering", command=lambda: self.run_pipeline(run_mode="clustering_only"))
        b1.pack(side="left")
        self.run_buttons.append(b1)
        self._tip(b1, "Run codon/AA usage computation, dimensional reduction, clustering, statistics, and workbook export, without plotting figures or decoding-strategy analyses.")
        b2 = ttk.Button(run_row, text="Run clustering followed by plotting", command=lambda: self.run_pipeline(run_mode="clustering_plots"))
        b2.pack(side="left", padx=(8, 0))
        self.run_buttons.append(b2)
        self._tip(b2, "Run the standard codon-usage clustering pipeline and then generate the figures selected in the Figures tab. Decoding-strategy analyses are not run from this button.")

        david_frame = ttk.LabelFrame(parent, text="Unbiased functional enrichment analysis using DAVID", padding=14)
        david_frame.pack(fill="x", pady=(12, 5))
        for c in range(6):
            david_frame.columnconfigure(c, weight=1)
        ttk.Label(
            david_frame,
            text=(
                "Please first run clustering analysis. This analysis extract genes using a sliding window from the reordered genome "
                "after clustering analysis and enquires for functional enrichment through DAVID bioinformatics."
            ),
            foreground="#444444",
            wraplength=900,
        ).grid(row=0, column=0, columnspan=6, sticky="w", pady=(0, 8))
        self._add_labeled_widget(david_frame, 1, 0, "DAVID email", ttk.Entry(david_frame, textvariable=self.david_email, width=28), "Registered DAVID email used for the web service.")
        self._add_labeled_widget(david_frame, 1, 2, "Window size", ttk.Entry(david_frame, textvariable=self.david_window_size, width=12), "Number of genes per sliding window.")
        self._add_labeled_widget(david_frame, 1, 4, "Step size", ttk.Entry(david_frame, textvariable=self.david_step_size, width=12), "Step between successive genome windows.")
        david_btn = ttk.Button(david_frame, text="Run DAVID sliding-window scan", command=self.run_david_sliding_window_scan)
        david_btn.grid(row=2, column=0, sticky="w", pady=(12, 0))
        self.run_buttons.append(david_btn)
        self._tip(david_btn, "Run DAVID functional enrichment on sliding windows from the last reordered genome. This does not rerun codon-usage clustering.")

    def _build_analysis_settings_tab(self, parent):
        parent.columnconfigure(0, weight=1)

        self.feature_settings_box = ttk.LabelFrame(parent, text="Codon-feature settings", padding=10)
        self.feature_settings_box.grid(row=0, column=0, sticky="ew", pady=5)
        self.feature_settings_box.columnconfigure(1, weight=1)
        self._add_labeled_widget(
            self.feature_settings_box, 0, 0, "Codon set (automatic)",
            ttk.Combobox(self.feature_settings_box, textvariable=self.codon_set, state="disabled", values=["64", "61", "59"], width=18),
            "Automatically selected from the codon usage metric: amino acid identity→64, absolute codon usage→61, relative codon usage→59."
        )
        chk = ttk.Checkbutton(self.feature_settings_box, text="Centered (by genomic average)", variable=self.center_features)
        chk.grid(row=1, column=0, columnspan=2, sticky="w", pady=3)
        self._tip(chk, "Subtract the whole-genome mean value for each feature before dimensional reduction and clustering.")
        chk = ttk.Checkbutton(self.feature_settings_box, text="Scaled to genomic deviation", variable=self.scale_features)
        chk.grid(row=2, column=0, columnspan=2, sticky="w", pady=3)
        self._tip(chk, "Divide each feature by its whole-genome standard deviation.")

        self.dimred_frame = ttk.LabelFrame(parent, text="Dimensional reduction settings", padding=10)
        self.dimred_frame.grid(row=1, column=0, sticky="ew", pady=6)
        self.dimred_frame.columnconfigure(0, weight=1)
        self.dimred_params_box = ttk.Frame(self.dimred_frame)
        self.dimred_params_box.grid(row=0, column=0, sticky="ew")

        self.cluster_frame = ttk.LabelFrame(parent, text="Clustering settings", padding=10)
        self.cluster_frame.grid(row=2, column=0, sticky="ew", pady=6)
        self.cluster_frame.columnconfigure(0, weight=1)
        self.cluster_params_box = ttk.Frame(self.cluster_frame)
        self.cluster_params_box.grid(row=0, column=0, sticky="ew")

        self.genome_axis_settings_box = ttk.LabelFrame(parent, text="Genome-axis preparation", padding=10)
        self.genome_axis_settings_box.grid(row=3, column=0, sticky="ew", pady=6)
        for c in range(4):
            self.genome_axis_settings_box.columnconfigure(c, weight=1)
        chk = ttk.Checkbutton(self.genome_axis_settings_box, text="Apply smoothing", variable=self.apply_smoothing)
        chk.grid(row=0, column=0, sticky="w", pady=4)
        self._tip(chk, "Smooth the ordered genome signal before downstream summaries.")
        self._add_labeled_widget(self.genome_axis_settings_box, 0, 1, "Smooth window genes", ttk.Entry(self.genome_axis_settings_box, textvariable=self.smooth_window_genes, width=12), "Window size used by the smoothing step.")
        chk = ttk.Checkbutton(self.genome_axis_settings_box, text="Apply binning", variable=self.apply_binning)
        chk.grid(row=1, column=0, sticky="w", pady=4)
        self._tip(chk, "Bin genes along the ordered genome axis for downstream summaries.")
        self._add_labeled_widget(self.genome_axis_settings_box, 1, 1, "Bin size genes", ttk.Entry(self.genome_axis_settings_box, textvariable=self.bin_size_genes, width=12), "Number of genes per bin.")

        self.ks_box = ttk.LabelFrame(parent, text="2D Kolmogorov-Smirnov settings", padding=10)
        self.ks_box.grid(row=4, column=0, sticky="ew", pady=6)
        for c in range(6):
            self.ks_box.columnconfigure(c, weight=1)
        self.ks_content = ttk.Frame(self.ks_box)
        self.ks_content.grid(row=0, column=0, columnspan=6, sticky="ew")
        for c in range(6):
            self.ks_content.columnconfigure(c, weight=1)
        self._add_labeled_widget(self.ks_content, 0, 0, "Alpha", ttk.Entry(self.ks_content, textvariable=self.ks_alpha, width=12), "Significance threshold.")
        self._add_labeled_widget(self.ks_content, 0, 2, "Method", ttk.Combobox(self.ks_content, textvariable=self.ks_method, state="readonly", values=["binned", "exact"], width=12), "KS implementation mode.")
        self._add_labeled_widget(self.ks_content, 0, 4, "Bins", ttk.Entry(self.ks_content, textvariable=self.ks_bins, width=12), "Number of bins for the binned implementation.")
        self._add_labeled_widget(self.ks_content, 1, 0, "Permutations", ttk.Entry(self.ks_content, textvariable=self.ks_n_perm, width=12), "Number of random permutations used for p-values.")
        self._add_labeled_widget(self.ks_content, 1, 2, "Random seed", ttk.Entry(self.ks_content, textvariable=self.ks_seed, width=12), "Random seed for reproducibility.")

        self.david_box = ttk.LabelFrame(parent, text="Advanced DAVID sliding-window scan settings", padding=10)
        self.david_box.grid(row=5, column=0, sticky="ew", pady=6)
        for c in range(6):
            self.david_box.columnconfigure(c, weight=1)
        self.david_content = ttk.Frame(self.david_box)
        self.david_content.grid(row=0, column=0, columnspan=6, sticky="ew")
        for c in range(6):
            self.david_content.columnconfigure(c, weight=1)
        self._add_labeled_widget(self.david_content, 0, 0, "Wait time", ttk.Entry(self.david_content, textvariable=self.david_wait_time, width=12), "Delay between DAVID requests, in seconds.")
        self._add_labeled_widget(self.david_content, 0, 2, "Max clusters", ttk.Entry(self.david_content, textvariable=self.david_max_clusters, width=12), "Number of top enrichment-derived cluster groups to keep.")
        self._add_labeled_widget(self.david_content, 0, 4, "Min valid IDs", ttk.Entry(self.david_content, textvariable=self.david_min_valid_ids, width=12), "Minimum number of valid IDs required per window.")
        self._add_labeled_widget(self.david_content, 1, 0, "Top N hits", ttk.Entry(self.david_content, textvariable=self.david_top_n_hits, width=12), "Top enrichment windows retained for output and term reconstruction.")

    def _build_figures_tab(self, parent):
        general = ttk.LabelFrame(parent, text="Figures", padding=14)
        general.pack(fill="x", pady=5)
        general.columnconfigure(1, weight=1)

        self._add_labeled_widget(
            general, 0, 0, "Image format",
            ttk.Combobox(general, textvariable=self.figure_format, state="readonly", values=["png", "tif", "tiff", "jpg", "jpeg", "pdf"], width=12),
            "Output image format used for saved figures."
        )

        picker_box = ttk.LabelFrame(parent, text="Figure cluster picker", padding=10)
        picker_box.pack(fill="x", pady=6)
        picker_box.columnconfigure(1, weight=1)
        ttk.Label(
            picker_box,
            text=(
                "Choose which clusters will be plotted from the Figures tab. "
                "This affects cluster-based figures only: clusters along the genome axis, "
                "2D cluster maps, and per-cluster codon-usage profiles. "
                "The codon-vs-genes heatmap remains a whole-genome heatmap."
            ),
            foreground="#606060",
            wraplength=820,
        ).grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, 4))
        b = ttk.Button(picker_box, text="Cluster picker", command=self._open_figure_clusters_dialog)
        b.grid(row=1, column=0, sticky="w", padx=(0, 8), pady=4)
        self.run_buttons.append(b)
        ttk.Label(picker_box, textvariable=self.figure_clusters_status, foreground="#444444").grid(row=1, column=1, columnspan=2, sticky="w", pady=4)

        list_box = ttk.LabelFrame(parent, text="Plot one figure", padding=14)
        list_box.pack(fill="x", pady=5)
        list_box.columnconfigure(0, weight=1)

        def add_plot_row(row, title, tooltip, command):
            lbl = ttk.Label(list_box, text=title, font=("Arial", 10, "bold"))
            lbl.grid(row=row, column=0, sticky="w", pady=(8 if row == 0 else 4, 4))
            self._tip(lbl, tooltip)
            btn = ttk.Button(list_box, text="PLOT", command=command)
            btn.grid(row=row, column=1, sticky="e", padx=(10, 0), pady=(8 if row == 0 else 4, 4))
            self.run_buttons.append(btn)

        plot_common_opts = ttk.Frame(list_box)
        plot_common_opts.grid(row=0, column=0, columnspan=2, sticky="w", pady=(0, 8))
        self._add_labeled_widget(
            plot_common_opts, 0, 0, "Max rows",
            ttk.Entry(plot_common_opts, textvariable=self.plot_rows, width=8),
            "Maximum number of subplot rows in the multi-panel 2D map. Default: 4."
        )
        self._add_labeled_widget(
            plot_common_opts, 0, 2, "Min genes threshold",
            ttk.Entry(plot_common_opts, textvariable=self.plot_cluster_min_genes, width=8),
            "Minimum number of genes required for a cluster to be included in cluster-based figures. Applies to Plot 2 and Plot 3, and to the equivalent plots generated after clustering. Default: 2."
        )

        add_plot_row(
            1,
            "PLOT 1. Codon vs genes heatmap",
            "Generate the main codon/feature × genes heatmap from the existing clustering workbook.",
            lambda: self.plot_existing_figure("figure_main_heatmap"),
        )
        add_plot_row(
            2,
            "PLOT 2. Clusters along genome axis",
            "Generate the gene-cluster localization heatmap along the reordered genome axis using the Figure cluster picker selection and the min-gene threshold above.",
            lambda: self.plot_existing_figure("figure_cluster_axis"),
        )
        add_plot_row(
            3,
            "PLOT 3. 2D dimensional reduction map",
            "Create the multi-panel UMAP/tSNE/PCA map using the Figure cluster picker selection, max rows, and min-gene threshold above.",
            lambda: self.plot_existing_figure("figure_2d_map"),
        )
        add_plot_row(
            4,
            "PLOT 4. Per-cluster codon usage profiles",
            "Create the multi-panel codon-usage summary figure using the Figure cluster picker selection.",
            lambda: self.plot_existing_figure("figure_codon_profiles"),
        )

        ttk.Label(
            parent,
            text="Detailed plot parameters are shown in the 'Figure details' tab.",
            foreground="#444444",
            wraplength=820,
        ).pack(anchor="w", pady=(8, 0))

    def _build_plots_details_tab(self, parent):
        parent.columnconfigure(0, weight=1)

        self.main_heatmap_details_box = ttk.LabelFrame(parent, text="Codon vs genes heatmap details", padding=10)
        self.main_heatmap_details_box.grid(row=0, column=0, sticky="ew", pady=6)
        for c in range(6):
            self.main_heatmap_details_box.columnconfigure(c, weight=1)
        chk = ttk.Checkbutton(self.main_heatmap_details_box, text="Custom figure aesthetics", variable=self.main_heatmap_custom_aesthetics, command=self._refresh_optional_sections)
        chk.grid(row=1, column=0, columnspan=3, sticky="w", pady=(6, 4))
        self._tip(chk, "Show optional controls for the main codon × genes heatmap.")
        self.main_heatmap_aesthetics_content = ttk.Frame(self.main_heatmap_details_box)
        self.main_heatmap_aesthetics_content.grid(row=2, column=0, columnspan=6, sticky="ew")
        for c in range(6):
            self.main_heatmap_aesthetics_content.columnconfigure(c, weight=1)
        self._add_labeled_widget(self.main_heatmap_aesthetics_content, 0, 0, "DPI", ttk.Entry(self.main_heatmap_aesthetics_content, textvariable=self.main_heatmap_dpi, width=12), "Saved resolution for the main heatmap.")
        self._add_labeled_widget(self.main_heatmap_aesthetics_content, 0, 2, "Colormap", ttk.Entry(self.main_heatmap_aesthetics_content, textvariable=self.main_heatmap_colormap, width=14), "Matplotlib or CodonPipe colormap name, e.g. parula, plasma, viridis.")
        self._add_labeled_widget(self.main_heatmap_aesthetics_content, 0, 4, "X tick every", ttk.Entry(self.main_heatmap_aesthetics_content, textvariable=self.main_heatmap_xtick_every, width=12), "Spacing between x-axis ticks along the gene axis.")
        self._add_labeled_widget(self.main_heatmap_aesthetics_content, 1, 0, "Figure width", ttk.Entry(self.main_heatmap_aesthetics_content, textvariable=self.main_heatmap_fig_width, width=12), "Main heatmap figure width in inches.")
        self._add_labeled_widget(self.main_heatmap_aesthetics_content, 1, 2, "Figure height", ttk.Entry(self.main_heatmap_aesthetics_content, textvariable=self.main_heatmap_fig_height, width=12), "Main heatmap figure height in inches.")
        self._add_labeled_widget(self.main_heatmap_aesthetics_content, 1, 4, "Color min", ttk.Entry(self.main_heatmap_aesthetics_content, textvariable=self.main_heatmap_caxis_min, width=12), "Lower color-axis limit for the main heatmap.")
        self._add_labeled_widget(self.main_heatmap_aesthetics_content, 2, 0, "Color max", ttk.Entry(self.main_heatmap_aesthetics_content, textvariable=self.main_heatmap_caxis_max, width=12), "Upper color-axis limit for the main heatmap.")
        chk = ttk.Checkbutton(self.main_heatmap_details_box, text="Custom axes limits", variable=self.main_heatmap_custom_axes, command=self._refresh_optional_sections)
        chk.grid(row=3, column=0, columnspan=3, sticky="w", pady=(6, 4))
        self.main_heatmap_axes_content = ttk.Frame(self.main_heatmap_details_box)
        self.main_heatmap_axes_content.grid(row=4, column=0, columnspan=6, sticky="ew")
        for c in range(6):
            self.main_heatmap_axes_content.columnconfigure(c, weight=1)
        self._add_labeled_widget(self.main_heatmap_axes_content, 0, 0, "X min", ttk.Entry(self.main_heatmap_axes_content, textvariable=self.main_heatmap_xmin, width=12), "Optional minimum x-axis value for the main heatmap.")
        self._add_labeled_widget(self.main_heatmap_axes_content, 0, 2, "X max", ttk.Entry(self.main_heatmap_axes_content, textvariable=self.main_heatmap_xmax, width=12), "Optional maximum x-axis value for the main heatmap.")
        self._add_labeled_widget(self.main_heatmap_axes_content, 0, 4, "Y min", ttk.Entry(self.main_heatmap_axes_content, textvariable=self.main_heatmap_ymin, width=12), "Optional minimum y-axis value for the main heatmap.")
        self._add_labeled_widget(self.main_heatmap_axes_content, 1, 0, "Y max", ttk.Entry(self.main_heatmap_axes_content, textvariable=self.main_heatmap_ymax, width=12), "Optional maximum y-axis value for the main heatmap.")

        self.heatmap_section = ttk.LabelFrame(parent, text="Clusters along genome axis details", padding=10)
        self.heatmap_section.grid(row=1, column=0, sticky="ew", pady=6)
        for c in range(6):
            self.heatmap_section.columnconfigure(c, weight=1)
        chk = ttk.Checkbutton(self.heatmap_section, text="Custom figure aesthetics", variable=self.gchm_custom_aesthetics, command=self._refresh_optional_sections)
        chk.grid(row=1, column=0, columnspan=3, sticky="w", pady=(6, 4))
        self.heatmap_aesthetics_content = ttk.Frame(self.heatmap_section)
        self.heatmap_aesthetics_content.grid(row=2, column=0, columnspan=6, sticky="ew")
        for c in range(6):
            self.heatmap_aesthetics_content.columnconfigure(c, weight=1)
        self._add_labeled_widget(self.heatmap_aesthetics_content, 0, 0, "DPI", ttk.Entry(self.heatmap_aesthetics_content, textvariable=self.heatmap_dpi, width=12), "Saved resolution for the genome-axis heatmap.")
        self._add_labeled_widget(self.heatmap_aesthetics_content, 0, 2, "Colormap", ttk.Combobox(self.heatmap_aesthetics_content, textvariable=self.gchm_colormap, state="readonly", values=["plasma", "viridis", "magma", "inferno", "cividis", "parula"], width=14), "Genome-axis heatmap colormap.")
        self._add_labeled_widget(self.heatmap_aesthetics_content, 0, 4, "Sigma", ttk.Entry(self.heatmap_aesthetics_content, textvariable=self.gchm_sigma, width=12), "Gaussian smoothing sigma along the genome axis.")
        self._add_labeled_widget(self.heatmap_aesthetics_content, 1, 0, "Spread factor", ttk.Entry(self.heatmap_aesthetics_content, textvariable=self.gchm_spread_factor, width=12), "Controls how broad each cluster signal appears.")
        self._add_labeled_widget(self.heatmap_aesthetics_content, 1, 2, "Height / cluster", ttk.Entry(self.heatmap_aesthetics_content, textvariable=self.gchm_height_per_cluster, width=12), "Vertical space per cluster row.")
        self._add_labeled_widget(self.heatmap_aesthetics_content, 1, 4, "Label fontsize", ttk.Entry(self.heatmap_aesthetics_content, textvariable=self.gchm_label_fontsize, width=12), "Cluster-label font size.")
        self._add_labeled_widget(self.heatmap_aesthetics_content, 2, 0, "Min rel", ttk.Entry(self.heatmap_aesthetics_content, textvariable=self.gchm_cmap_min_rel, width=12), "Minimum relative density used for the color range.")
        self._add_labeled_widget(self.heatmap_aesthetics_content, 2, 2, "Max rel", ttk.Entry(self.heatmap_aesthetics_content, textvariable=self.gchm_cmap_max_rel, width=12), "Maximum relative density used for the color range.")
        chk = ttk.Checkbutton(self.heatmap_section, text="Custom axes limits", variable=self.gchm_custom_axes, command=self._refresh_optional_sections)
        chk.grid(row=3, column=0, columnspan=3, sticky="w", pady=(6, 4))
        self.heatmap_axes_content = ttk.Frame(self.heatmap_section)
        self.heatmap_axes_content.grid(row=4, column=0, columnspan=6, sticky="ew")
        for c in range(6):
            self.heatmap_axes_content.columnconfigure(c, weight=1)
        self._add_labeled_widget(self.heatmap_axes_content, 0, 0, "X min", ttk.Entry(self.heatmap_axes_content, textvariable=self.gchm_xmin, width=12), "Optional minimum x-axis value for the genome-axis heatmap.")
        self._add_labeled_widget(self.heatmap_axes_content, 0, 2, "X max", ttk.Entry(self.heatmap_axes_content, textvariable=self.gchm_xmax, width=12), "Optional maximum x-axis value for the genome-axis heatmap.")
        self._add_labeled_widget(self.heatmap_axes_content, 0, 4, "Y min", ttk.Entry(self.heatmap_axes_content, textvariable=self.gchm_ymin, width=12), "Optional minimum y-axis value for the genome-axis heatmap.")
        self._add_labeled_widget(self.heatmap_axes_content, 1, 0, "Y max", ttk.Entry(self.heatmap_axes_content, textvariable=self.gchm_ymax, width=12), "Optional maximum y-axis value for the genome-axis heatmap.")

        self.density_section = ttk.LabelFrame(parent, text="2D dimensional reduction map details", padding=10)
        self.density_section.grid(row=2, column=0, sticky="ew", pady=6)
        for c in range(6):
            self.density_section.columnconfigure(c, weight=1)
        self.density_content = ttk.Frame(self.density_section)
        self.density_content.grid(row=0, column=0, columnspan=6, sticky="ew")
        for c in range(6):
            self.density_content.columnconfigure(c, weight=1)
        self._add_labeled_widget(self.density_content, 0, 0, "Color mode", ttk.Combobox(self.density_content, textvariable=self.color_mode, state="readonly", values=["enrichment", "density"], width=14), "Choose enrichment or raw density coloring for the 2D plots.")
        chk = ttk.Checkbutton(self.density_content, text="Show colorbar", variable=self.show_colorbar)
        chk.grid(row=0, column=2, sticky="w", pady=4)
        chk = ttk.Checkbutton(self.density_content, text="Custom figure aesthetics", variable=self.density_custom_aesthetics, command=self._refresh_optional_sections)
        chk.grid(row=2, column=0, columnspan=3, sticky="w", pady=(6, 4))
        self.density_aesthetics_content = ttk.Frame(self.density_content)
        self.density_aesthetics_content.grid(row=3, column=0, columnspan=6, sticky="ew")
        for c in range(6):
            self.density_aesthetics_content.columnconfigure(c, weight=1)
        self._add_labeled_widget(self.density_aesthetics_content, 0, 0, "DPI", ttk.Entry(self.density_aesthetics_content, textvariable=self.density_figure_dpi, width=12), "Saved resolution for the 2D map.")
        self._add_labeled_widget(self.density_aesthetics_content, 0, 2, "Optional title", ttk.Entry(self.density_aesthetics_content, textvariable=self.figure_suptitle, width=24), "Optional overall title for the 2D map.")
        self._add_labeled_widget(self.density_aesthetics_content, 1, 0, "Panel width (in)", ttk.Entry(self.density_aesthetics_content, textvariable=self.density_panel_w_in, width=12), "Width of each 2D map panel.")
        self._add_labeled_widget(self.density_aesthetics_content, 1, 2, "Panel height (in)", ttk.Entry(self.density_aesthetics_content, textvariable=self.density_panel_h_in, width=12), "Height of each 2D map panel.")
        self._add_labeled_widget(self.density_aesthetics_content, 1, 4, "Horizontal spacing", ttk.Entry(self.density_aesthetics_content, textvariable=self.density_subplot_wspace, width=12), "Spacing between panels horizontally.")
        self._add_labeled_widget(self.density_aesthetics_content, 2, 0, "Vertical spacing", ttk.Entry(self.density_aesthetics_content, textvariable=self.density_subplot_hspace, width=12), "Spacing between panels vertically.")
        self._add_labeled_widget(self.density_aesthetics_content, 2, 2, "Density cmap", ttk.Entry(self.density_aesthetics_content, textvariable=self.density_cmap, width=18), "Colormap used for density mode.")
        self._add_labeled_widget(self.density_aesthetics_content, 2, 4, "Enrichment cmap", ttk.Entry(self.density_aesthetics_content, textvariable=self.enrichment_cmap, width=18), "Colormap used for enrichment mode.")
        chk = ttk.Checkbutton(self.density_content, text="Custom axes limits", variable=self.density_custom_axes, command=self._refresh_optional_sections)
        chk.grid(row=4, column=0, columnspan=3, sticky="w", pady=(6, 4))
        self.density_axes_content = ttk.Frame(self.density_content)
        self.density_axes_content.grid(row=5, column=0, columnspan=6, sticky="ew")
        for c in range(6):
            self.density_axes_content.columnconfigure(c, weight=1)
        self._add_labeled_widget(self.density_axes_content, 0, 0, "X min", ttk.Entry(self.density_axes_content, textvariable=self.density_xmin, width=12), "Optional minimum x-axis value for the 2D map.")
        self._add_labeled_widget(self.density_axes_content, 0, 2, "X max", ttk.Entry(self.density_axes_content, textvariable=self.density_xmax, width=12), "Optional maximum x-axis value for the 2D map.")
        self._add_labeled_widget(self.density_axes_content, 0, 4, "Y min", ttk.Entry(self.density_axes_content, textvariable=self.density_ymin, width=12), "Optional minimum y-axis value for the 2D map.")
        self._add_labeled_widget(self.density_axes_content, 1, 0, "Y max", ttk.Entry(self.density_axes_content, textvariable=self.density_ymax, width=12), "Optional maximum y-axis value for the 2D map.")

        self.codon_section = ttk.LabelFrame(parent, text="Per cluster codon usage profile details", padding=10)
        self.codon_section.grid(row=3, column=0, sticky="ew", pady=6)
        for c in range(6):
            self.codon_section.columnconfigure(c, weight=1)
        self.codon_content = ttk.Frame(self.codon_section)
        self.codon_content.grid(row=0, column=0, columnspan=6, sticky="ew")
        for c in range(6):
            self.codon_content.columnconfigure(c, weight=1)
        self._add_labeled_widget(self.codon_content, 0, 0, "Plot mode", ttk.Combobox(self.codon_content, textvariable=self.codon_usage_plot_mode, state="readonly", values=CODON_USAGE_MODE_CHOICES, width=18), "Choose absolute codon usage, relative codon usage, or RCU z-scores.")
        chk = ttk.Checkbutton(self.codon_content, text="Custom figure aesthetics", variable=self.codon_custom_aesthetics, command=self._refresh_optional_sections)
        chk.grid(row=1, column=0, columnspan=3, sticky="w", pady=(6, 4))
        self.codon_aesthetics_content = ttk.Frame(self.codon_content)
        self.codon_aesthetics_content.grid(row=2, column=0, columnspan=6, sticky="ew")
        for c in range(6):
            self.codon_aesthetics_content.columnconfigure(c, weight=1)
        self._add_labeled_widget(self.codon_aesthetics_content, 0, 0, "DPI", ttk.Entry(self.codon_aesthetics_content, textvariable=self.codon_usage_dpi, width=12), "Saved resolution for the codon-usage plot.")
        self._add_labeled_widget(self.codon_aesthetics_content, 0, 2, "Panel width (in)", ttk.Entry(self.codon_aesthetics_content, textvariable=self.codon_panel_w_in, width=12), "Width of each codon-usage panel.")
        self._add_labeled_widget(self.codon_aesthetics_content, 0, 4, "Panel height (in)", ttk.Entry(self.codon_aesthetics_content, textvariable=self.codon_panel_h_in, width=12), "Height of each codon-usage panel.")
        chk = ttk.Checkbutton(self.codon_content, text="Custom axes limits", variable=self.codon_custom_axes, command=self._refresh_optional_sections)
        chk.grid(row=3, column=0, columnspan=3, sticky="w", pady=(6, 4))
        self.codon_axes_content = ttk.Frame(self.codon_content)
        self.codon_axes_content.grid(row=4, column=0, columnspan=6, sticky="ew")
        for c in range(6):
            self.codon_axes_content.columnconfigure(c, weight=1)
        self._add_labeled_widget(self.codon_axes_content, 0, 0, "X min", ttk.Entry(self.codon_axes_content, textvariable=self.codon_xmin, width=12), "Optional minimum x-axis value for the codon-usage plot.")
        self._add_labeled_widget(self.codon_axes_content, 0, 2, "X max", ttk.Entry(self.codon_axes_content, textvariable=self.codon_xmax, width=12), "Optional maximum x-axis value for the codon-usage plot.")
        self._add_labeled_widget(self.codon_axes_content, 0, 4, "Y min", ttk.Entry(self.codon_axes_content, textvariable=self.codon_ymin, width=12), "Optional minimum y-axis value for the codon-usage plot.")
        self._add_labeled_widget(self.codon_axes_content, 1, 0, "Y max", ttk.Entry(self.codon_axes_content, textvariable=self.codon_ymax, width=12), "Optional maximum y-axis value for the codon-usage plot.")

        self.codon_compare_details_box = ttk.LabelFrame(parent, text="Codon usage analyses plot details", padding=10)
        self.codon_compare_details_box.grid(row=4, column=0, sticky="ew", pady=6)
        for c in range(6):
            self.codon_compare_details_box.columnconfigure(c, weight=1)
        ttk.Label(
            self.codon_compare_details_box,
            text=(
                "Advanced settings for the two buttons in the Codon usage analyses tab: "
                "Plot raw codon usage and Plot correlations. These controls do not affect "
                "the main per-cluster codon-usage profiles generated from the Figures tab."
            ),
            foreground="#444444",
            wraplength=900,
        ).grid(row=0, column=0, columnspan=6, sticky="w", pady=(0, 6))

        chk = ttk.Checkbutton(
            self.codon_compare_details_box,
            text="Custom axes for Plot raw codon usage",
            variable=self.codon_compare_raw_custom_axes,
            command=self._refresh_optional_sections,
        )
        chk.grid(row=1, column=0, columnspan=4, sticky="w", pady=(6, 4))
        self._tip(chk, "Enable manual x/y limits for the raw codon-usage comparison plot generated from the Codon usage analyses tab.")
        self.codon_compare_raw_axes_content = ttk.Frame(self.codon_compare_details_box)
        self.codon_compare_raw_axes_content.grid(row=2, column=0, columnspan=6, sticky="ew")
        for c in range(6):
            self.codon_compare_raw_axes_content.columnconfigure(c, weight=1)
        self._add_labeled_widget(self.codon_compare_raw_axes_content, 0, 0, "X min", ttk.Entry(self.codon_compare_raw_axes_content, textvariable=self.codon_compare_raw_xmin, width=12), "Optional minimum x-axis value for the raw codon/amino-acid comparison plot. With codon AA gaps enabled, X is the plotted coordinate, not simply the codon index.")
        self._add_labeled_widget(self.codon_compare_raw_axes_content, 0, 2, "X max", ttk.Entry(self.codon_compare_raw_axes_content, textvariable=self.codon_compare_raw_xmax, width=12), "Optional maximum x-axis value for the raw codon/amino-acid comparison plot. With codon AA gaps enabled, X is the plotted coordinate, not simply the codon index.")
        self._add_labeled_widget(self.codon_compare_raw_axes_content, 0, 4, "Y min", ttk.Entry(self.codon_compare_raw_axes_content, textvariable=self.codon_compare_raw_ymin, width=12), "Optional minimum y-axis value for the raw codon-usage plot.")
        self._add_labeled_widget(self.codon_compare_raw_axes_content, 1, 0, "Y max", ttk.Entry(self.codon_compare_raw_axes_content, textvariable=self.codon_compare_raw_ymax, width=12), "Optional maximum y-axis value for the raw codon-usage plot.")

        self.codon_compare_style_content = ttk.Frame(self.codon_compare_details_box)
        self.codon_compare_style_content.grid(row=3, column=0, columnspan=6, sticky="ew", pady=(10, 2))
        for c in range(6):
            self.codon_compare_style_content.columnconfigure(c, weight=1)
        self._add_labeled_widget(self.codon_compare_style_content, 0, 0, "Caption size", ttk.Entry(self.codon_compare_style_content, textvariable=self.codon_compare_caption_size, width=12), "Font size for x/y tick labels, x/y labels, plot title, statistical stars and legend in the Codon usage analyses plots.")
        self._add_labeled_widget(self.codon_compare_style_content, 0, 2, "Dot size", ttk.Entry(self.codon_compare_style_content, textvariable=self.codon_compare_marker_size, width=12), "Marker size for mean/line/correlation points. Larger values make amino-acid comparison plots easier to read.")
        self._add_labeled_widget(self.codon_compare_style_content, 0, 4, "Line width", ttk.Entry(self.codon_compare_style_content, textvariable=self.codon_compare_line_width, width=12), "Line width for line plots, error bars, boxplot outlines and the y=x reference line.")
        self._add_labeled_widget(self.codon_compare_style_content, 1, 0, "AA spacing", ttk.Entry(self.codon_compare_style_content, textvariable=self.codon_compare_aa_spacing, width=12), "Horizontal spacing between amino-acid groups. Used only for Amino acid identity plots; larger values add more white space between amino acids.")
        self._add_labeled_widget(self.codon_compare_style_content, 1, 2, "Legend entries/row", ttk.Entry(self.codon_compare_style_content, textvariable=self.codon_compare_legend_ncol, width=12), "Maximum legend entries per row. Default is 3, so extra entries wrap onto a new row above the plot.")
        self._add_labeled_widget(self.codon_compare_style_content, 1, 4, "Codon AA gap", ttk.Entry(self.codon_compare_style_content, textvariable=self.codon_compare_codon_aa_gap, width=12), "Extra horizontal gap inserted between amino-acid families in raw codon-usage plots. Codons within the same synonymous family remain close together.")
        self._add_labeled_widget(self.codon_compare_style_content, 2, 0, "Codon spacing", ttk.Entry(self.codon_compare_style_content, textvariable=self.codon_compare_codon_gap, width=12), "Extra horizontal gap inserted between individual codons in raw codon-usage plots. Default 0.33 is approximately half of the default amino-acid-family gap.")
        self._add_labeled_widget(self.codon_compare_style_content, 2, 2, "Highlight color", ttk.Entry(self.codon_compare_style_content, textvariable=self.codon_compare_highlight_color, width=12), "Color used for highlighted codon x-axis labels in Plot raw usage. Default #0057D9 is a strong blue that contrasts with black labels.")

        chk = ttk.Checkbutton(
            self.codon_compare_details_box,
            text="Custom axes for Plot correlations",
            variable=self.codon_compare_corr_custom_axes,
            command=self._refresh_optional_sections,
        )
        chk.grid(row=4, column=0, columnspan=4, sticky="w", pady=(10, 4))
        self._tip(chk, "Enable manual x/y limits for pairwise codon-usage correlation plots. The dotted reference remains y = x.")
        self.codon_compare_corr_axes_content = ttk.Frame(self.codon_compare_details_box)
        self.codon_compare_corr_axes_content.grid(row=5, column=0, columnspan=6, sticky="ew")
        for c in range(6):
            self.codon_compare_corr_axes_content.columnconfigure(c, weight=1)
        self._add_labeled_widget(self.codon_compare_corr_axes_content, 0, 0, "X min", ttk.Entry(self.codon_compare_corr_axes_content, textvariable=self.codon_compare_corr_xmin, width=12), "Optional minimum x-axis value for correlation plots.")
        self._add_labeled_widget(self.codon_compare_corr_axes_content, 0, 2, "X max", ttk.Entry(self.codon_compare_corr_axes_content, textvariable=self.codon_compare_corr_xmax, width=12), "Optional maximum x-axis value for correlation plots.")
        self._add_labeled_widget(self.codon_compare_corr_axes_content, 0, 4, "Y min", ttk.Entry(self.codon_compare_corr_axes_content, textvariable=self.codon_compare_corr_ymin, width=12), "Optional minimum y-axis value for correlation plots.")
        self._add_labeled_widget(self.codon_compare_corr_axes_content, 1, 0, "Y max", ttk.Entry(self.codon_compare_corr_axes_content, textvariable=self.codon_compare_corr_ymax, width=12), "Optional maximum y-axis value for correlation plots.")

        self.decoding_details_box = ttk.LabelFrame(parent, text="Decoding strategies figure details", padding=10)
        self.decoding_details_box.grid(row=6, column=0, sticky="ew", pady=6)
        for c in range(6):
            self.decoding_details_box.columnconfigure(c, weight=1)
        ttk.Label(
            self.decoding_details_box,
            text="Advanced settings for decoding-strategy heatmaps and boxplots. These controls are available here rather than in the Decoding strategies tab.",
            foreground="#444444",
            wraplength=900,
        ).grid(row=0, column=0, columnspan=6, sticky="w", pady=(0, 6))

        chk = ttk.Checkbutton(self.decoding_details_box, text="Customize gene-ordered heatmaps 1 and 2", variable=self.trna_supp_heatmaps_customize, command=self._refresh_optional_sections)
        chk.grid(row=1, column=0, columnspan=4, sticky="w", pady=(6, 4))
        self.trna_gene_heatmap_custom_content = ttk.Frame(self.decoding_details_box)
        self.trna_gene_heatmap_custom_content.grid(row=2, column=0, columnspan=6, sticky="ew")
        for c in range(6):
            self.trna_gene_heatmap_custom_content.columnconfigure(c, weight=1)
        self._add_labeled_widget(self.trna_gene_heatmap_custom_content, 0, 0, "DPI", ttk.Entry(self.trna_gene_heatmap_custom_content, textvariable=self.trna_supp_heatmaps_dpi, width=12), "Saved resolution for gene-ordered heatmaps 1 and 2. Prefilled with the automatic default used by the pipeline.")
        self._add_labeled_widget(self.trna_gene_heatmap_custom_content, 0, 2, "Figure width (in)", ttk.Entry(self.trna_gene_heatmap_custom_content, textvariable=self.trna_supp_heatmaps_fig_width, width=12), "Figure width for gene-ordered heatmaps 1 and 2. Prefilled with the current main codon-heatmap width.")
        self.trna_supp_fig_height_entry = ttk.Entry(self.trna_gene_heatmap_custom_content, textvariable=self.trna_supp_heatmaps_fig_height, width=12)
        self._add_labeled_widget(self.trna_gene_heatmap_custom_content, 0, 4, "Figure height (in)", self.trna_supp_fig_height_entry, "Figure height for gene-ordered heatmaps 1 and 2. Prefilled with a useful baseline; adjust upward for many rows.")
        self.trna_supp_cell_height_entry = ttk.Entry(self.trna_gene_heatmap_custom_content, textvariable=self.trna_supp_heatmaps_cell_height, width=12)
        self._add_labeled_widget(self.trna_gene_heatmap_custom_content, 1, 0, "Cell height (in)", self.trna_supp_cell_height_entry, "Per-row height used when automatic sizing is derived from row count. Prefilled with the main heatmap row-height baseline.")
        self._add_labeled_widget(self.trna_gene_heatmap_custom_content, 1, 2, "X tick every genes", ttk.Entry(self.trna_gene_heatmap_custom_content, textvariable=self.trna_supp_heatmaps_xtick_every_genes, width=12), "Optional x-axis tick spacing, in genes, for gene-ordered heatmaps 1 and 2.")
        self._add_labeled_widget(self.trna_gene_heatmap_custom_content, 1, 4, "Y tick fontsize", ttk.Entry(self.trna_gene_heatmap_custom_content, textvariable=self.trna_supp_heatmaps_ytick_fontsize, width=12), "Optional y-axis label font size for gene-ordered heatmaps 1 and 2.")
        self._add_labeled_widget(self.trna_gene_heatmap_custom_content, 2, 0, "Title fontsize", ttk.Entry(self.trna_gene_heatmap_custom_content, textvariable=self.trna_supp_heatmaps_title_fontsize, width=12), "Optional title font size for gene-ordered heatmaps 1 and 2.")
        self._add_labeled_widget(self.trna_gene_heatmap_custom_content, 2, 2, "X min", ttk.Entry(self.trna_gene_heatmap_custom_content, textvariable=self.trna_supp_heatmaps_xmin, width=12), "Optional minimum x-axis limit for gene-ordered heatmaps 1 and 2.")
        self._add_labeled_widget(self.trna_gene_heatmap_custom_content, 2, 4, "X max", ttk.Entry(self.trna_gene_heatmap_custom_content, textvariable=self.trna_supp_heatmaps_xmax, width=12), "Optional maximum x-axis limit for gene-ordered heatmaps 1 and 2.")
        self._add_labeled_widget(self.trna_gene_heatmap_custom_content, 3, 0, "Y min", ttk.Entry(self.trna_gene_heatmap_custom_content, textvariable=self.trna_supp_heatmaps_ymin, width=12), "Optional minimum y-axis limit for gene-ordered heatmaps 1 and 2.")
        self._add_labeled_widget(self.trna_gene_heatmap_custom_content, 3, 2, "Y max", ttk.Entry(self.trna_gene_heatmap_custom_content, textvariable=self.trna_supp_heatmaps_ymax, width=12), "Optional maximum y-axis limit for gene-ordered heatmaps 1 and 2.")

        chk = ttk.Checkbutton(self.decoding_details_box, text="Customize cluster-level heatmaps/boxplots 3–5", variable=self.trna_shift_heatmaps_customize, command=self._refresh_optional_sections)
        chk.grid(row=3, column=0, columnspan=4, sticky="w", pady=(10, 4))
        self.trna_shift_heatmaps_custom_content = ttk.Frame(self.decoding_details_box)
        self.trna_shift_heatmaps_custom_content.grid(row=4, column=0, columnspan=6, sticky="ew")
        for c in range(6):
            self.trna_shift_heatmaps_custom_content.columnconfigure(c, weight=1)
        self._add_labeled_widget(self.trna_shift_heatmaps_custom_content, 0, 0, "DPI", ttk.Entry(self.trna_shift_heatmaps_custom_content, textvariable=self.trna_shift_heatmaps_dpi, width=12), "Optional saved resolution for cluster-level decoding plots 3–5.")
        self.trna_shift_fig_width_entry = ttk.Entry(self.trna_shift_heatmaps_custom_content, textvariable=self.trna_shift_heatmaps_fig_width, width=12)
        self._add_labeled_widget(self.trna_shift_heatmaps_custom_content, 0, 2, "Figure width (in)", self.trna_shift_fig_width_entry, "Optional figure width for cluster-level decoding plots 3–5.")
        self.trna_shift_fig_height_entry = ttk.Entry(self.trna_shift_heatmaps_custom_content, textvariable=self.trna_shift_heatmaps_fig_height, width=12)
        self._add_labeled_widget(self.trna_shift_heatmaps_custom_content, 0, 4, "Figure height (in)", self.trna_shift_fig_height_entry, "Optional figure height for cluster-level decoding plots 3–5.")
        self.trna_shift_cell_width_entry = ttk.Entry(self.trna_shift_heatmaps_custom_content, textvariable=self.trna_shift_heatmaps_cell_width, width=12)
        self._add_labeled_widget(self.trna_shift_heatmaps_custom_content, 1, 0, "Cell width (in / col)", self.trna_shift_cell_width_entry, "Optional cell width for cluster-level heatmaps. Leave blank to let figure width control the plot.")
        self.trna_shift_cell_height_entry = ttk.Entry(self.trna_shift_heatmaps_custom_content, textvariable=self.trna_shift_heatmaps_cell_height, width=12)
        self._add_labeled_widget(self.trna_shift_heatmaps_custom_content, 1, 2, "Cell height (relative)", self.trna_shift_cell_height_entry, "Relative cell height for cluster-level heatmaps 3–5.")
        self._add_labeled_widget(self.trna_shift_heatmaps_custom_content, 1, 4, "X tick fontsize", ttk.Entry(self.trna_shift_heatmaps_custom_content, textvariable=self.trna_shift_heatmaps_xtick_fontsize, width=12), "Optional x-axis label font size.")
        self._add_labeled_widget(self.trna_shift_heatmaps_custom_content, 2, 0, "Y tick fontsize", ttk.Entry(self.trna_shift_heatmaps_custom_content, textvariable=self.trna_shift_heatmaps_ytick_fontsize, width=12), "Optional y-axis label font size.")
        self._add_labeled_widget(self.trna_shift_heatmaps_custom_content, 2, 2, "Title fontsize", ttk.Entry(self.trna_shift_heatmaps_custom_content, textvariable=self.trna_shift_heatmaps_title_fontsize, width=12), "Optional title font size.")
        self._add_labeled_widget(self.trna_shift_heatmaps_custom_content, 2, 4, "X min", ttk.Entry(self.trna_shift_heatmaps_custom_content, textvariable=self.trna_shift_heatmaps_xmin, width=12), "Optional minimum x-axis limit.")
        self._add_labeled_widget(self.trna_shift_heatmaps_custom_content, 3, 0, "X max", ttk.Entry(self.trna_shift_heatmaps_custom_content, textvariable=self.trna_shift_heatmaps_xmax, width=12), "Optional maximum x-axis limit.")
        self._add_labeled_widget(self.trna_shift_heatmaps_custom_content, 3, 2, "Y min", ttk.Entry(self.trna_shift_heatmaps_custom_content, textvariable=self.trna_shift_heatmaps_ymin, width=12), "Optional minimum y-axis limit.")
        self._add_labeled_widget(self.trna_shift_heatmaps_custom_content, 3, 4, "Y max", ttk.Entry(self.trna_shift_heatmaps_custom_content, textvariable=self.trna_shift_heatmaps_ymax, width=12), "Optional maximum y-axis limit.")
        ttk.Checkbutton(self.trna_shift_heatmaps_custom_content, text="Log2 colorbar (tRNA usage heatmap)", variable=self.trna_shift_heatmap_log2_colorbar).grid(row=4, column=0, columnspan=2, sticky="w", pady=(4, 2))
        ttk.Checkbutton(self.trna_shift_heatmaps_custom_content, text="Log2 colorbar (wobble heatmap)", variable=self.trna_wobble_heatmap_log2_colorbar).grid(row=4, column=2, columnspan=2, sticky="w", pady=(4, 2))
        self._add_labeled_widget(self.trna_shift_heatmaps_custom_content, 5, 0, "Bracket type (tRNA)", ttk.Combobox(self.trna_shift_heatmaps_custom_content, textvariable=self.trna_shift_heatmap_bracket_type, state="readonly", values=["square", "brace"], width=12), "Bracket style for tRNA usage heatmap amino-acid group annotations.")
        self._add_labeled_widget(self.trna_shift_heatmaps_custom_content, 5, 2, "Bracket type (wobble)", ttk.Combobox(self.trna_shift_heatmaps_custom_content, textvariable=self.trna_wobble_heatmap_bracket_type, state="readonly", values=["square", "brace"], width=12), "Bracket style for wobble heatmap amino-acid group annotations.")
        self._add_labeled_widget(self.trna_shift_heatmaps_custom_content, 6, 0, "Bracket X (tRNA)", ttk.Entry(self.trna_shift_heatmaps_custom_content, textvariable=self.trna_shift_heatmap_bracket_x, width=12), "Horizontal bracket offset for tRNA usage heatmap.")
        self._add_labeled_widget(self.trna_shift_heatmaps_custom_content, 6, 2, "AA label X (tRNA)", ttk.Entry(self.trna_shift_heatmaps_custom_content, textvariable=self.trna_shift_heatmap_label_x, width=12), "Horizontal amino-acid label offset for tRNA usage heatmap.")
        self._add_labeled_widget(self.trna_shift_heatmaps_custom_content, 7, 0, "Bracket X (wobble)", ttk.Entry(self.trna_shift_heatmaps_custom_content, textvariable=self.trna_wobble_heatmap_bracket_x, width=12), "Horizontal bracket offset for wobble heatmap.")
        self._add_labeled_widget(self.trna_shift_heatmaps_custom_content, 7, 2, "AA label X (wobble)", ttk.Entry(self.trna_shift_heatmaps_custom_content, textvariable=self.trna_wobble_heatmap_label_x, width=12), "Horizontal amino-acid label offset for wobble heatmap.")

        self.trna_modification_plots_details_content = ttk.LabelFrame(self.decoding_details_box, text="tRNA modification plots", padding=8)
        self.trna_modification_plots_details_content.grid(row=5, column=0, columnspan=6, sticky="ew", pady=(10, 2))
        for c in range(6):
            self.trna_modification_plots_details_content.columnconfigure(c, weight=1)
        chk = ttk.Checkbutton(self.trna_modification_plots_details_content, text="Customize tRNA modification plots (Plot 6)", variable=self.trna_modification_plots_customize, command=self._refresh_optional_sections)
        chk.grid(row=0, column=0, columnspan=4, sticky="w", pady=(0, 4))
        self.trna_modification_plots_custom_content = ttk.Frame(self.trna_modification_plots_details_content)
        self.trna_modification_plots_custom_content.grid(row=1, column=0, columnspan=6, sticky="ew")
        for c in range(6):
            self.trna_modification_plots_custom_content.columnconfigure(c, weight=1)
        self._add_labeled_widget(self.trna_modification_plots_custom_content, 0, 0, "DPI", ttk.Entry(self.trna_modification_plots_custom_content, textvariable=self.trna_modification_plots_dpi, width=12), "Optional saved resolution for Plot 6 tRNA modification/enzyme boxplots.")
        self._add_labeled_widget(self.trna_modification_plots_custom_content, 0, 2, "Figure width (in)", ttk.Entry(self.trna_modification_plots_custom_content, textvariable=self.trna_modification_plots_fig_width, width=12), "Figure width used for Plot 6 by default. Prefilled with the default value; adjust if needed.")
        self._add_labeled_widget(self.trna_modification_plots_custom_content, 0, 4, "Figure height (in)", ttk.Entry(self.trna_modification_plots_custom_content, textvariable=self.trna_modification_plots_fig_height, width=12), "Figure height used for Plot 6 by default. Prefilled with the default value; increase this if the plot area still looks compressed.")
        self._add_labeled_widget(self.trna_modification_plots_custom_content, 1, 0, "Caption size", ttk.Entry(self.trna_modification_plots_custom_content, textvariable=self.trna_modification_plots_caption_size, width=12), "Font size for Plot 6 title, legend, y-axis title, ticks, modification labels, and position labels.")
        self._add_labeled_widget(self.trna_modification_plots_custom_content, 1, 2, "Y min", ttk.Entry(self.trna_modification_plots_custom_content, textvariable=self.trna_modification_plots_ymin, width=12), "Default y-axis minimum for Plot 6 boxplot/violin mode. Prefilled with the default value (-2.1).")
        self._add_labeled_widget(self.trna_modification_plots_custom_content, 1, 4, "Y max", ttk.Entry(self.trna_modification_plots_custom_content, textvariable=self.trna_modification_plots_ymax, width=12), "Default y-axis maximum for Plot 6 boxplot/violin mode. Prefilled with the default value (3.5).")
        self._add_labeled_widget(self.trna_modification_plots_custom_content, 2, 0, "Group line Y", ttk.Entry(self.trna_modification_plots_custom_content, textvariable=self.trna_modification_plots_group_bar_y, width=12), "Vertical position of the position-32/34/37 horizontal grouping lines in axis coordinates. More negative moves them lower; default is -0.22.")
        self._add_labeled_widget(self.trna_modification_plots_custom_content, 2, 2, "Group label gap", ttk.Entry(self.trna_modification_plots_custom_content, textvariable=self.trna_modification_plots_group_label_gap, width=12), "Vertical gap between the grouping line and its text label. Default is 0.06.")
        self._add_labeled_widget(self.trna_modification_plots_custom_content, 2, 4, "Star offset", ttk.Entry(self.trna_modification_plots_custom_content, textvariable=self.trna_modification_plots_star_offset, width=12), "Distance of significance-star rows below the upper y-axis frame, as a fraction of the y-axis span. Default is 0.07.")
        self._add_labeled_widget(self.trna_modification_plots_custom_content, 3, 0, "Legend columns", ttk.Entry(self.trna_modification_plots_custom_content, textvariable=self.trna_modification_plots_legend_ncol, width=12), "Optional number of legend columns for Plot 6. Leave blank to auto-wrap.")
        self._add_labeled_widget(self.trna_modification_plots_custom_content, 3, 2, "Box width", ttk.Entry(self.trna_modification_plots_custom_content, textvariable=self.trna_modification_plots_box_width, width=12), "Width of each Plot 6 boxplot. Default is 0.18.")

        self.trna_secondary_axis_details_content = ttk.LabelFrame(self.decoding_details_box, text="Secondary y-axis and boxplot overlay", padding=8)
        self.trna_secondary_axis_details_content.grid(row=6, column=0, columnspan=6, sticky="ew", pady=(10, 2))
        for c in range(6):
            self.trna_secondary_axis_details_content.columnconfigure(c, weight=1)
        self._add_labeled_widget(self.trna_secondary_axis_details_content, 0, 0, "Secondary axis style", ttk.Combobox(self.trna_secondary_axis_details_content, textvariable=self.trna_secondary_axis_style, state="readonly", values=["bars", "dots", "line", "none"], width=12), "How to display genomic codon frequency or tRNA abundance on boxplot secondary axes.")
        self._add_labeled_widget(self.trna_secondary_axis_details_content, 0, 2, "Grey overlay alpha", ttk.Entry(self.trna_secondary_axis_details_content, textvariable=self.trna_secondary_axis_alpha, width=12), "Transparency for grey secondary-axis bars/dots/line.")
        self._add_labeled_widget(self.trna_secondary_axis_details_content, 0, 4, "Grey bar width", ttk.Entry(self.trna_secondary_axis_details_content, textvariable=self.trna_secondary_axis_bar_width, width=12), "Width of the grey secondary-axis bars.")
        self._add_labeled_widget(self.trna_secondary_axis_details_content, 1, 0, "Box width", ttk.Entry(self.trna_secondary_axis_details_content, textvariable=self.trna_boxplot_width, width=12), "Width of each grouped boxplot.")
        chk = ttk.Checkbutton(self.trna_secondary_axis_details_content, text="Show jittered points", variable=self.trna_boxplot_show_points)
        chk.grid(row=1, column=2, sticky="w", pady=4)
        self._add_labeled_widget(self.trna_secondary_axis_details_content, 1, 4, "Point alpha", ttk.Entry(self.trna_secondary_axis_details_content, textvariable=self.trna_boxplot_point_alpha, width=12), "Transparency for jittered per-gene points over boxplots.")

    def _build_trna_tab(self, parent):
        box = ttk.LabelFrame(parent, text="Decoding strategies", padding=10)
        box.pack(fill="x", pady=6)
        for c in range(6):
            box.columnconfigure(c, weight=1)

        self.trna_content = ttk.Frame(box)
        self.trna_content.grid(row=0, column=0, columnspan=6, sticky="ew")
        for c in range(6):
            self.trna_content.columnconfigure(c, weight=1)

        note = ttk.Label(
            self.trna_content,
            text=(
                f"Build {CODONPIPE_GUI_BUILD}. "
                "Plots in this tab use the decoding table selected in Input/Output and reuse the existing clustering workbook. "
                "PLOT buttons do not rerun UMAP/tSNE/PCA or gene clustering. "
                "Use the cluster picker below only for cluster-level plots 4–6. Smoothing controls for plots 1–3 are independent and can be edited directly in each plot frame."
            ),
            foreground="#606060",
            wraplength=880,
        )
        note.grid(row=0, column=0, columnspan=6, sticky="w", pady=(0, 10))

        picker_box = ttk.LabelFrame(self.trna_content, text="Cluster picker for cluster-level decoding plots", padding=8)
        picker_box.grid(row=1, column=0, columnspan=6, sticky="ew", pady=(0, 8))
        picker_box.columnconfigure(1, weight=1)
        ttk.Label(
            picker_box,
            text="Choose which gene clusters will be used for plots 4–6 only. Plots 1–3 use the full reordered genome.",
            foreground="#606060",
            wraplength=760,
        ).grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, 4))
        b = ttk.Button(picker_box, text="Cluster picker", command=self._open_decoding_clusters_dialog)
        b.grid(row=1, column=0, sticky="w", pady=4)
        self.run_buttons.append(b)
        ttk.Label(picker_box, textvariable=self.decoding_clusters_status, foreground="#444444").grid(row=1, column=1, sticky="w", padx=(10, 0), pady=4)
        bref = ttk.Button(picker_box, text="Reference cluster for statistical analyses", command=self._open_decoding_reference_cluster_dialog)
        bref.grid(row=2, column=0, sticky="w", pady=(4, 2))
        self.run_buttons.append(bref)
        ttk.Label(picker_box, textvariable=self.decoding_reference_cluster_status, foreground="#444444").grid(row=2, column=1, sticky="w", padx=(10, 0), pady=(4, 2))

        # ------------------------------------------------------------------
        # Plot 1 — Wobble decoding along reordered genome
        # ------------------------------------------------------------------
        self.trna_wobble_gene_box = ttk.LabelFrame(self.trna_content, text="Plot 1 — wobble decoding along reordered genome", padding=8)
        self.trna_wobble_gene_box.grid(row=2, column=0, columnspan=6, sticky="ew", pady=(6, 6))
        for c in range(6):
            self.trna_wobble_gene_box.columnconfigure(c, weight=1)
        ttk.Label(
            self.trna_wobble_gene_box,
            text=("Displays wobble-associated synonymous codon usage along the reordered genome. "
                  "Heatmap mode shows codon-level ZCU values. Line/area mode shows the average wobble percentage across WC/wobble codon pairs."),
            foreground="#444444",
            wraplength=820,
        ).grid(row=0, column=0, columnspan=5, sticky="w", pady=(0, 4))
        self._add_labeled_widget(
            self.trna_wobble_gene_box, 1, 0, "Display as",
            ttk.Combobox(self.trna_wobble_gene_box, textvariable=self.trna_gene_wobble_plot_kind, state="readonly", values=["heatmap", "line", "area"], width=12),
            "Heatmap = codon x gene ZCU; line/area = average wobble percentage across selected synonymous codon pairs."
        )
        self.trna_wobble_gene_smoothing_frame = ttk.Frame(self.trna_wobble_gene_box)
        self.trna_wobble_gene_smoothing_frame.grid(row=2, column=0, columnspan=5, sticky="w", pady=4)
        self._add_labeled_widget(self.trna_wobble_gene_smoothing_frame, 0, 0, "Smoothing", ttk.Combobox(self.trna_wobble_gene_smoothing_frame, textvariable=self.trna_gene_wobble_smooth, state="readonly", values=["none", "running average", "running median", "gaussian"], width=16), "Smoothing method for line/area values along reordered genes.")
        self._add_labeled_widget(self.trna_wobble_gene_smoothing_frame, 0, 2, "Smoothing window size (in genes)", ttk.Entry(self.trna_wobble_gene_smoothing_frame, textvariable=self.trna_gene_wobble_smooth_window, width=8), "Editable number of genes used for smoothing; default 40. Use 1 for no effective smoothing.")
        self._add_labeled_widget(self.trna_wobble_gene_smoothing_frame, 0, 4, "Caption size", ttk.Entry(self.trna_wobble_gene_smoothing_frame, textvariable=self.trna_gene_wobble_caption_size, width=8), "Font size for title, axes, ticks, and legend in line/area mode.")
        b = ttk.Button(self.trna_wobble_gene_box, text="PLOT", command=lambda: self.plot_existing_decoding_figure("trna_figure_single_box"))
        b.grid(row=2, column=5, sticky="w", pady=4)
        self.run_buttons.append(b)

        # ------------------------------------------------------------------
        # Plot 2 — tRNA usage shift along reordered genome
        # ------------------------------------------------------------------
        self.trna_gene_heatmap_box = ttk.LabelFrame(self.trna_content, text="Plot 2 — tRNA usage shift along reordered genome", padding=8)
        self.trna_gene_heatmap_box.grid(row=3, column=0, columnspan=6, sticky="ew", pady=(6, 6))
        for c in range(6):
            self.trna_gene_heatmap_box.columnconfigure(c, weight=1)
        ttk.Label(
            self.trna_gene_heatmap_box,
            text=("Displays tRNA decoding-group usage along the reordered genome. "
                  "Heatmap mode shows selected tRNA ZTU values. Line/area mode shows the fraction of rare tRNA decoding within selected amino-acid families."),
            foreground="#444444",
            wraplength=820,
        ).grid(row=0, column=0, columnspan=5, sticky="w", pady=(0, 4))
        self._add_labeled_widget(
            self.trna_gene_heatmap_box, 1, 0, "Display as",
            ttk.Combobox(self.trna_gene_heatmap_box, textvariable=self.trna_gene_trna_plot_kind, state="readonly", values=["heatmap", "line", "area"], width=12),
            "Heatmap = tRNA x gene ZTU; line/area = rare-tRNA fraction profiles for selected amino-acid families."
        )
        self.trna_gene_heatmap_smoothing_frame = ttk.Frame(self.trna_gene_heatmap_box)
        self.trna_gene_heatmap_smoothing_frame.grid(row=2, column=0, columnspan=5, sticky="w", pady=4)
        self._add_labeled_widget(self.trna_gene_heatmap_smoothing_frame, 0, 0, "Smoothing", ttk.Combobox(self.trna_gene_heatmap_smoothing_frame, textvariable=self.trna_gene_trna_smooth, state="readonly", values=["none", "running average", "running median", "gaussian"], width=16), "Smoothing method for line/area values along reordered genes.")
        self._add_labeled_widget(self.trna_gene_heatmap_smoothing_frame, 0, 2, "Smoothing window size (in genes)", ttk.Entry(self.trna_gene_heatmap_smoothing_frame, textvariable=self.trna_gene_trna_smooth_window, width=8), "Editable number of genes used for smoothing; default 40. Use 1 for no effective smoothing.")
        self._add_labeled_widget(self.trna_gene_heatmap_smoothing_frame, 0, 4, "Caption size", ttk.Entry(self.trna_gene_heatmap_smoothing_frame, textvariable=self.trna_gene_trna_caption_size, width=8), "Font size for title, axes, ticks, and legend in line/area mode.")
        b = ttk.Button(self.trna_gene_heatmap_box, text="PLOT", command=lambda: self.plot_existing_decoding_figure("trna_figure_gene_heatmap"))
        b.grid(row=2, column=5, sticky="w", pady=4)
        self.run_buttons.append(b)

        # ------------------------------------------------------------------
        # Plot 3 — mRNA stability along reordered genome
        # ------------------------------------------------------------------
        self.trna_mrna_stability_box = ttk.LabelFrame(self.trna_content, text="Plot 3 — mRNA stability along reordered genome", padding=8)
        self.trna_mrna_stability_box.grid(row=4, column=0, columnspan=6, sticky="ew", pady=(6, 6))
        for c in range(6):
            self.trna_mrna_stability_box.columnconfigure(c, weight=1)
        ttk.Label(
            self.trna_mrna_stability_box,
            text=("Plots mRNA half-life values from the decoding workbook sheet named 'RNA stability'. "
                  "The sheet is matched by locus tag and then reordered according to the clustering workbook; missing genes are left blank."),
            foreground="#444444",
            wraplength=820,
        ).grid(row=0, column=0, columnspan=5, sticky="w", pady=(0, 4))
        self._add_labeled_widget(
            self.trna_mrna_stability_box, 1, 0, "Display as",
            ttk.Combobox(self.trna_mrna_stability_box, textvariable=self.trna_mrna_stability_plot_kind, state="readonly", values=["line", "area"], width=12),
            "Line or shaded-area plot of mRNA half-life along reordered genes."
        )
        self.trna_mrna_stability_smoothing_frame = ttk.Frame(self.trna_mrna_stability_box)
        self.trna_mrna_stability_smoothing_frame.grid(row=2, column=0, columnspan=5, sticky="w", pady=4)
        self._add_labeled_widget(self.trna_mrna_stability_smoothing_frame, 0, 0, "Smoothing", ttk.Combobox(self.trna_mrna_stability_smoothing_frame, textvariable=self.trna_mrna_stability_smooth, state="readonly", values=["none", "running average", "running median", "gaussian"], width=16), "Smoothing method for mRNA half-life values along reordered genes.")
        self._add_labeled_widget(self.trna_mrna_stability_smoothing_frame, 0, 2, "Smoothing window size (in genes)", ttk.Entry(self.trna_mrna_stability_smoothing_frame, textvariable=self.trna_mrna_stability_smooth_window, width=8), "Editable number of genes used for smoothing; default 100. Use 1 for no effective smoothing.")
        self._add_labeled_widget(self.trna_mrna_stability_smoothing_frame, 0, 4, "Caption size", ttk.Entry(self.trna_mrna_stability_smoothing_frame, textvariable=self.trna_mrna_stability_caption_size, width=8), "Font size for title, axes, ticks, and legend in line/area mode.")
        b = ttk.Button(self.trna_mrna_stability_box, text="PLOT", command=lambda: self.plot_existing_decoding_figure("trna_figure_rna_stability"))
        b.grid(row=2, column=5, sticky="w", pady=4)
        self.run_buttons.append(b)

        # ------------------------------------------------------------------
        # Plot 4 — Enrichment in wobble decoding within clusters
        # ------------------------------------------------------------------
        self.trna_wobble_shift_box = ttk.LabelFrame(self.trna_content, text="Plot 4 — enrichment in wobble decoding within clusters", padding=8)
        self.trna_wobble_shift_box.grid(row=5, column=0, columnspan=6, sticky="ew", pady=(6, 6))
        for c in range(6):
            self.trna_wobble_shift_box.columnconfigure(c, weight=1)
        ttk.Label(
            self.trna_wobble_shift_box,
            text="Cluster-level z-score synonymous codon usage for WC/wobble codon pairs. Grey secondary-axis bars show genomic codon frequency (%).",
            foreground="#444444",
            wraplength=820,
        ).grid(row=0, column=0, columnspan=5, sticky="w", pady=(0, 4))
        self._add_labeled_widget(
            self.trna_wobble_shift_box, 1, 0, "Display as",
            ttk.Combobox(self.trna_wobble_shift_box, textvariable=self.trna_wobble_plot_kind, state="readonly", values=["heatmap", "boxplot", "violin"], width=12),
            "Choose heatmap of cluster means, mean±SD boxplots, or violin plots of per-gene values."
        )
        self._add_labeled_widget(self.trna_wobble_shift_box, 2, 0, "Y min", ttk.Entry(self.trna_wobble_shift_box, textvariable=self.trna_wobble_boxplot_ymin, width=8), "Optional boxplot/violin y-axis minimum for plot 4.")
        self._add_labeled_widget(self.trna_wobble_shift_box, 2, 2, "Y max", ttk.Entry(self.trna_wobble_shift_box, textvariable=self.trna_wobble_boxplot_ymax, width=8), "Optional boxplot/violin y-axis maximum for plot 4.")
        self._add_labeled_widget(self.trna_wobble_shift_box, 2, 4, "Log2", ttk.Combobox(self.trna_wobble_shift_box, textvariable=self.trna_wobble_boxplot_log2, state="readonly", values=["yes", "no"], width=8), "Apply signed log2 transform: sign(x)*log2(abs(x)+1).")
        self._add_labeled_widget(self.trna_wobble_shift_box, 3, 0, "Exclude outliers", ttk.Combobox(self.trna_wobble_shift_box, textvariable=self.trna_wobble_exclude_outliers, state="readonly", values=["no", "yes"], width=8), "Exclude per-feature values outside mean ± SD multiplier before plotting.")
        self._add_labeled_widget(self.trna_wobble_shift_box, 3, 2, "Outlier SD x", ttk.Entry(self.trna_wobble_shift_box, textvariable=self.trna_wobble_outlier_sd, width=8), "Multiplier used for mean ± x*SD outlier filtering.")
        self._add_labeled_widget(self.trna_wobble_shift_box, 4, 0, "Caption size", ttk.Entry(self.trna_wobble_shift_box, textvariable=self.trna_wobble_boxplot_caption_size, width=8), "Font size for title, axes, ticks, and legend in boxplot/violin mode.")
        self._add_labeled_widget(self.trna_wobble_shift_box, 4, 2, "Stats test", ttk.Combobox(self.trna_wobble_shift_box, textvariable=self.trna_wobble_stats_test, state="readonly", values=["none", "Student t-test", "Welch t-test", "Mann-Whitney U"], width=16), "Statistical comparison of each selected cluster against the reference cluster for each x-axis feature. These stars are drawn without a horizontal bar.")
        self._add_labeled_widget(self.trna_wobble_shift_box, 5, 0, "Codon-pair stats", ttk.Combobox(self.trna_wobble_shift_box, textvariable=self.trna_wobble_pair_stats_test, state="readonly", values=["none", "Student t-test", "Welch t-test", "Mann-Whitney U"], width=16), "Within each plotted cluster, statistically compare all displayed synonymous-codon pairs within each amino-acid group. These stars are drawn above horizontal bars.")
        self._add_labeled_widget(self.trna_wobble_shift_box, 5, 2, "Bracket gap", ttk.Entry(self.trna_wobble_shift_box, textvariable=self.trna_wobble_pair_stats_gap, width=8), "Base vertical distance between bracket rows, expressed as a fraction of the y-axis data span. The first gap above cluster-vs-reference stars is 1.5× this value. Try 0.03–0.15.")
        b = ttk.Button(self.trna_wobble_shift_box, text="PLOT", command=lambda: self.plot_existing_decoding_figure("trna_figure_wobble"))
        b.grid(row=1, column=5, sticky="w", pady=4)
        self.run_buttons.append(b)

        # ------------------------------------------------------------------
        # Plot 5 — tRNA usage shift within clusters
        # ------------------------------------------------------------------
        self.trna_shift_heatmap_box = ttk.LabelFrame(self.trna_content, text="Plot 5 — tRNA usage shift within clusters", padding=8)
        self.trna_shift_heatmap_box.grid(row=6, column=0, columnspan=6, sticky="ew", pady=(6, 6))
        for c in range(6):
            self.trna_shift_heatmap_box.columnconfigure(c, weight=1)
        ttk.Label(
            self.trna_shift_heatmap_box,
            text="Cluster-level tRNA usage shift versus the genome. Grey secondary-axis values show absolute genomic tRNA abundance from the decoding table.",
            foreground="#444444",
            wraplength=820,
        ).grid(row=0, column=0, columnspan=5, sticky="w", pady=(0, 4))
        self._add_labeled_widget(
            self.trna_shift_heatmap_box, 1, 0, "Display as",
            ttk.Combobox(self.trna_shift_heatmap_box, textvariable=self.trna_shift_plot_kind, state="readonly", values=["heatmap", "boxplot", "violin"], width=12),
            "Choose heatmap of cluster means, mean±SD boxplots, or violin plots of per-gene values."
        )
        self._add_labeled_widget(self.trna_shift_heatmap_box, 2, 0, "Y min", ttk.Entry(self.trna_shift_heatmap_box, textvariable=self.trna_shift_boxplot_ymin, width=8), "Optional boxplot/violin y-axis minimum for plot 5.")
        self._add_labeled_widget(self.trna_shift_heatmap_box, 2, 2, "Y max", ttk.Entry(self.trna_shift_heatmap_box, textvariable=self.trna_shift_boxplot_ymax, width=8), "Optional boxplot/violin y-axis maximum for plot 5.")
        self._add_labeled_widget(self.trna_shift_heatmap_box, 2, 4, "Log2", ttk.Combobox(self.trna_shift_heatmap_box, textvariable=self.trna_shift_boxplot_log2, state="readonly", values=["yes", "no"], width=8), "Apply signed log2 transform: sign(x)*log2(abs(x)+1).")
        self._add_labeled_widget(self.trna_shift_heatmap_box, 3, 0, "Exclude outliers", ttk.Combobox(self.trna_shift_heatmap_box, textvariable=self.trna_shift_exclude_outliers, state="readonly", values=["no", "yes"], width=8), "Exclude per-feature values outside mean ± SD multiplier before plotting.")
        self._add_labeled_widget(self.trna_shift_heatmap_box, 3, 2, "Outlier SD x", ttk.Entry(self.trna_shift_heatmap_box, textvariable=self.trna_shift_outlier_sd, width=8), "Multiplier used for mean ± x*SD outlier filtering.")
        self._add_labeled_widget(self.trna_shift_heatmap_box, 4, 0, "Caption size", ttk.Entry(self.trna_shift_heatmap_box, textvariable=self.trna_shift_boxplot_caption_size, width=8), "Font size for title, axes, ticks, and legend in boxplot/violin mode.")
        self._add_labeled_widget(self.trna_shift_heatmap_box, 4, 2, "Stats test", ttk.Combobox(self.trna_shift_heatmap_box, textvariable=self.trna_shift_stats_test, state="readonly", values=["none", "Student t-test", "Welch t-test", "Mann-Whitney U"], width=16), "Statistical comparison of each selected cluster against the reference cluster for each x-axis feature. These stars are drawn without a horizontal bar.")
        self._add_labeled_widget(self.trna_shift_heatmap_box, 5, 0, "Isoacceptor-pair stats", ttk.Combobox(self.trna_shift_heatmap_box, textvariable=self.trna_shift_pair_stats_test, state="readonly", values=["none", "Student t-test", "Welch t-test", "Mann-Whitney U"], width=16), "Within each plotted cluster, statistically compare two-feature isoacceptor groups. Multi-isoacceptor groups such as Arg are skipped to keep the plot readable.")
        self._add_labeled_widget(self.trna_shift_heatmap_box, 5, 2, "Bracket gap", ttk.Entry(self.trna_shift_heatmap_box, textvariable=self.trna_shift_pair_stats_gap, width=8), "Base vertical distance between bracket rows, expressed as a fraction of the y-axis data span. The first gap above cluster-vs-reference stars is 1.5× this value. Try 0.03–0.15.")
        b = ttk.Button(self.trna_shift_heatmap_box, text="PLOT", command=lambda: self.plot_existing_decoding_figure("trna_figure_shift"))
        b.grid(row=1, column=5, sticky="w", pady=4)
        self.run_buttons.append(b)

        # ------------------------------------------------------------------
        # Plot 6 — tRNA modifications shift within clusters
        # ------------------------------------------------------------------
        self.trna_modification_box = ttk.LabelFrame(self.trna_content, text="Plot 6 — tRNA modifications shift within clusters", padding=8)
        self.trna_modification_box.grid(row=7, column=0, columnspan=6, sticky="ew", pady=(6, 6))
        for c in range(6):
            self.trna_modification_box.columnconfigure(c, weight=1)
        ttk.Label(
            self.trna_modification_box,
            text="Cluster-level enrichment of codons associated with tRNA modifications or enzymes. Plot 6 exports conservative and permissive decoder-assignment versions using only the full decoder-level table; the compact/pooled table is ignored for this plot.",
            foreground="#444444",
            wraplength=820,
        ).grid(row=0, column=0, columnspan=5, sticky="w", pady=(0, 4))
        self._add_labeled_widget(
            self.trna_modification_box, 1, 0, "Display as",
            ttk.Combobox(self.trna_modification_box, textvariable=self.trna_modifications_plot_kind, state="readonly", values=["heatmap", "boxplot", "violin"], width=12),
            "Choose heatmap of cluster means, mean±SD boxplots, or violin plots of per-gene values."
        )
        self._add_labeled_widget(self.trna_modification_box, 2, 0, "Show", ttk.Combobox(self.trna_modification_box, textvariable=self.trna_modifications_feature_mode, state="readonly", values=["modifications", "enzymes"], width=14), "Use column F modification names or column G tRNA-modification enzyme names for plot 6.")
        self._add_labeled_widget(self.trna_modification_box, 3, 0, "Log2", ttk.Combobox(self.trna_modification_box, textvariable=self.trna_modifications_boxplot_log2, state="readonly", values=["yes", "no"], width=8), "Apply signed log2 transform: sign(x)*log2(abs(x)+1).")
        self._add_labeled_widget(self.trna_modification_box, 3, 2, "Exclude outliers", ttk.Combobox(self.trna_modification_box, textvariable=self.trna_modifications_exclude_outliers, state="readonly", values=["no", "yes"], width=8), "Exclude per-feature values outside mean ± SD multiplier before plotting.")
        self._add_labeled_widget(self.trna_modification_box, 3, 4, "Outlier SD x", ttk.Entry(self.trna_modification_box, textvariable=self.trna_modifications_outlier_sd, width=8), "Multiplier used for mean ± x*SD outlier filtering.")
        self._add_labeled_widget(self.trna_modification_box, 4, 0, "Caption size", ttk.Entry(self.trna_modification_box, textvariable=self.trna_modifications_boxplot_caption_size, width=8), "Font size for title, axes, ticks, and legend in boxplot/violin mode.")
        self._add_labeled_widget(self.trna_modification_box, 4, 2, "Stats test", ttk.Combobox(self.trna_modification_box, textvariable=self.trna_modifications_stats_test, state="readonly", values=["none", "Student t-test", "Welch t-test", "Mann-Whitney U"], width=16), "Statistical comparison of each selected cluster against the reference cluster for each x-axis feature.")
        fbtn = ttk.Button(self.trna_modification_box, text="Select modifications/enzymes", command=self._open_trna_modification_feature_dialog)
        fbtn.grid(row=5, column=0, sticky="w", pady=(6, 2))
        self.run_buttons.append(fbtn)
        ttk.Label(self.trna_modification_box, textvariable=self.trna_modifications_selection_status, foreground="#444444").grid(row=5, column=1, columnspan=2, sticky="w", padx=(8, 0), pady=(6, 2))
        aabtn = ttk.Button(self.trna_modification_box, text="Select amino acids", command=self._open_trna_modification_aa_dialog)
        aabtn.grid(row=5, column=3, sticky="w", pady=(6, 2))
        self.run_buttons.append(aabtn)
        ttk.Label(self.trna_modification_box, textvariable=self.trna_modification_aas_status, foreground="#444444").grid(row=5, column=4, columnspan=2, sticky="w", padx=(8, 0), pady=(6, 2))
        b = ttk.Button(self.trna_modification_box, text="PLOT", command=lambda: self.plot_existing_decoding_figure("trna_figure_modifications"))
        b.grid(row=1, column=5, sticky="w", pady=4)
        self.run_buttons.append(b)

        trna_run_row = ttk.Frame(self.trna_content)
        trna_run_row.grid(row=8, column=0, columnspan=6, sticky="w", pady=(14, 0))
        b = ttk.Button(trna_run_row, text="Plot all figures", command=lambda: self.plot_existing_decoding_figure("trna_all"))
        b.pack(side="left")
        self.run_buttons.append(b)
        self._tip(b, "Run all decoding-strategy figures using the decoding table from Input/Output and the selected clusters from this tab.")

        # Show smoothing controls only for profile plots (line/area), not heatmaps.
        try:
            for _var in (self.trna_gene_wobble_plot_kind, self.trna_gene_trna_plot_kind, self.trna_mrna_stability_plot_kind):
                _var.trace_add("write", lambda *_args: self._refresh_decoding_profile_smoothing_controls())
            self.trna_modifications_feature_mode.trace_add("write", lambda *_args: self._update_trna_modification_selection_status())
        except Exception:
            pass
        self._refresh_decoding_profile_smoothing_controls()
        self._update_trna_modification_selection_status()
        self._update_trna_modification_aas_status()

    def _refresh_decoding_profile_smoothing_controls(self):
        """Show smoothing controls only for line/area profile displays in Decoding strategies plots 1–3."""
        configs = [
            (getattr(self, "trna_wobble_gene_smoothing_frame", None), getattr(self, "trna_gene_wobble_plot_kind", None)),
            (getattr(self, "trna_gene_heatmap_smoothing_frame", None), getattr(self, "trna_gene_trna_plot_kind", None)),
            (getattr(self, "trna_mrna_stability_smoothing_frame", None), getattr(self, "trna_mrna_stability_plot_kind", None)),
        ]
        for frame, var in configs:
            if frame is None or var is None:
                continue
            try:
                kind = str(var.get() or "").strip().lower().replace("surface", "area")
            except Exception:
                kind = "line"
            try:
                if kind in {"line", "area"}:
                    frame.grid()
                else:
                    frame.grid_remove()
            except Exception:
                pass

    def _sync_codon_set_from_usage_basis(self):
        mapping = {'AA': '64', 'ACU': '61', 'RCU': '59'}
        self.codon_set.set(mapping.get(_usage_internal(self.usage_basis.get()), '61'))

    def _sync_analysis_choice_booleans(self):
        self.enable_2d_ks.set(self.statistical_test_method.get().strip() == "2D Kolmogorov-Smirnov")
        self._refresh_optional_sections()

    def _selected_preloaded_genome(self):
        label = self.preloaded_genome_choice.get().strip()
        for rec in self.preloaded_genomes:
            if rec.get("label") == label:
                return rec
        return self.preloaded_genomes[0] if self.preloaded_genomes else None

    def _apply_preloaded_genome_choice(self, force=False):
        if self.genome_not_available.get():
            return
        rec = self._selected_preloaded_genome()
        if rec:
            key = str(rec.get("label") or rec.get("fasta") or "")
            if (not force) and key and key == getattr(self, "_preloaded_applied_key", ""):
                return
            self._preloaded_applied_key = key
            self.fasta_path.set(rec.get("fasta", ""))
            self.organism_name.set(rec.get("organism") or rec.get("label") or "Organism")
            # Keep old behavior: preloaded FASTA, cluster workbook and decoding
            # workbook are visible as soon as the interface opens. The expensive
            # part from recent builds was not the path assignment itself, but
            # recursive folder scans and Excel sheet parsing. Companion lookup is
            # now fast path-based matching only.
            self._auto_load_companion_files_for_current_genome()

    def _auto_load_companion_files_for_current_genome(self):
        """Load preloaded cluster/decoding companion paths without slow sheet parsing."""
        self._auto_load_cluster_file_for_current_genome()
        if not self.trna_decoding_table_path.get().strip():
            auto_decoding = self._auto_find_decoding_table_for_current_genome()
            if auto_decoding:
                self.trna_decoding_table_path.set(auto_decoding)
        self._update_decoding_status()

    def _deferred_startup_companion_autoload(self):
        """Run lightweight companion-file autoloading after the GUI is visible."""
        try:
            self._defer_preloaded_companion_autoload = False
            if not self.genome_not_available.get():
                self._apply_preloaded_genome_choice(force=True)
            self._refresh_cluster_source_widgets()
            self._update_active_clusters_status()
        except Exception as e:
            try:
                self._append_log(f"[WARN] Deferred preloaded companion-file detection failed: {e}\n")
            except Exception:
                pass

    def _refresh_genome_source_widgets(self):
        if not self.preloaded_genomes:
            self.genome_not_available.set(True)
        use_manual = bool(self.genome_not_available.get())
        if hasattr(self, "preloaded_genome_label"):
            if use_manual:
                self.preloaded_genome_label.grid_remove()
                self.preloaded_genome_combo.grid_remove()
            else:
                self.preloaded_genome_label.grid()
                self.preloaded_genome_combo.grid()
                self._apply_preloaded_genome_choice()
        if hasattr(self, "manual_genome_frame"):
            if use_manual:
                self.manual_genome_frame.grid()
            else:
                self.manual_genome_frame.grid_remove()

    def _preloaded_search_dirs_for_current_genome(self):
        """Folders searched for preloaded companion files for the selected genome."""
        folders = []
        try:
            fasta = self.fasta_path.get().strip()
            if fasta:
                folders.append(os.path.dirname(os.path.abspath(fasta)))
        except Exception:
            pass
        for d in PRELOADED_GENOME_DIRS:
            try:
                if d:
                    folders.append(os.path.abspath(os.path.expanduser(str(d))))
            except Exception:
                pass
        out, seen = [], set()
        for folder in folders:
            folder = os.path.abspath(os.path.expanduser(str(folder or "")))
            key = os.path.normcase(folder)
            if folder and key not in seen and os.path.isdir(folder):
                seen.add(key)
                out.append(folder)
        return out

    def _current_preloaded_name_variants(self):
        rec = None if self.genome_not_available.get() else self._selected_preloaded_genome()
        variants = _preloaded_name_variants_from_record(rec)
        if self.organism_name.get().strip():
            variants.append(self.organism_name.get().strip())
        if self.fasta_path.get().strip():
            variants.append(_clean_genome_label_from_filename(self.fasta_path.get().strip()))
        out, seen = [], set()
        for v in variants:
            norm = _normalize_preloaded_name_for_match(v)
            if norm and norm not in seen:
                seen.add(norm)
                out.append(norm)
        return out
    def _current_preloaded_raw_name_variants(self):
        """Readable species/strain name variants used for exact companion-file paths."""
        rec = None if self.genome_not_available.get() else self._selected_preloaded_genome()
        variants = _preloaded_name_variants_from_record(rec)
        if self.organism_name.get().strip():
            variants.append(self.organism_name.get().strip())
        if self.fasta_path.get().strip():
            variants.append(_clean_genome_label_from_filename(self.fasta_path.get().strip()))
        out, seen = [], set()
        for v in variants:
            raw = " ".join(str(v or "").replace("_", " ").split()).strip()
            key = raw.lower()
            if raw and key not in seen:
                seen.add(key)
                out.append(raw)
        return out

    def _fast_find_preloaded_companion(self, kind):
        """Fast companion workbook lookup for startup.

        First tries exact expected names such as
        'Salmonella enterica SL1344 clusters.xlsx' or
        'Salmonella enterica SL1344 decoding table.xlsx'. If no exact match is
        found, it scans only the top level of the preloaded folders with
        os.scandir(). It never opens Excel files and never scans recursively.
        """
        cache_key = (
            str(kind),
            str(self.preloaded_genome_choice.get() if hasattr(self, "preloaded_genome_choice") else ""),
            str(self.fasta_path.get() if hasattr(self, "fasta_path") else ""),
            str(self.organism_name.get() if hasattr(self, "organism_name") else ""),
        )
        cache = getattr(self, "_preloaded_companion_cache", {})
        if cache_key in cache:
            cached = cache[cache_key]
            if cached == "" or os.path.isfile(cached):
                return cached

        if kind == "clusters":
            suffixes = ["clusters", "cluster"]
        elif kind == "decoding":
            suffixes = ["decoding table", "decoding", "codon anticodon table", "trna decoding table"]
        else:
            suffixes = [str(kind)]

        folders = self._preloaded_search_dirs_for_current_genome()
        raw_names = self._current_preloaded_raw_name_variants()
        extensions = [".xlsx", ".xlsm", ".xls"]

        # Exact expected names first: this keeps normal startup essentially O(1).
        for folder in folders:
            for raw in raw_names:
                for suffix in suffixes:
                    for ext in extensions:
                        candidate = os.path.join(folder, f"{raw} {suffix}{ext}")
                        if os.path.isfile(candidate):
                            cache[cache_key] = candidate
                            self._preloaded_companion_cache = cache
                            return candidate

        # Fallback: top-level scan only, scored like previous versions.
        hits = []
        fasta_dir = ""
        try:
            fasta_dir = os.path.dirname(os.path.abspath(self.fasta_path.get().strip()))
        except Exception:
            pass
        for folder in folders:
            try:
                with os.scandir(folder) as it:
                    for entry in it:
                        try:
                            if not entry.is_file():
                                continue
                            path = entry.path
                        except Exception:
                            continue
                        score = self._score_preloaded_companion_file(path, kind)
                        if score is None:
                            continue
                        if fasta_dir and os.path.abspath(folder) == os.path.abspath(fasta_dir):
                            score += 10
                        hits.append((-score, path))
            except Exception:
                continue
        if not hits:
            cache[cache_key] = ""
            self._preloaded_companion_cache = cache
            return ""
        hits.sort(key=lambda x: (x[0], x[1].lower()))
        cache[cache_key] = hits[0][1]
        self._preloaded_companion_cache = cache
        return hits[0][1]


    def _score_preloaded_companion_file(self, path, kind):
        """Score a likely preloaded cluster/decoding workbook for the current genome."""
        fname = os.path.basename(str(path or ""))
        low = fname.lower()
        if low.startswith("~$") or not low.endswith((".xlsx", ".xls", ".xlsm")):
            return None
        stem_norm = _normalize_preloaded_name_for_match(os.path.splitext(fname)[0])
        variants = self._current_preloaded_name_variants()
        if kind == "clusters":
            if "cluster" not in stem_norm:
                return None
            if any(x in stem_norm for x in ("decoding", "anticodon", "trna")):
                return None
            exact_suffixes = [f"{v} clusters" for v in variants] + [f"{v} cluster" for v in variants]
            if stem_norm in exact_suffixes:
                return 100
            score = 10
            if stem_norm == "clusters" or stem_norm == "cluster":
                score += 20
            for v in variants:
                if v and v in stem_norm:
                    score += 40
                    break
            return score
        if kind == "decoding":
            if not any(x in stem_norm for x in ("decoding", "decode", "anticodon", "trna", "codon")):
                return None
            exact_suffixes = [f"{v} decoding table" for v in variants]
            if stem_norm in exact_suffixes:
                return 100
            score = 0
            for token, pts in [
                ("decoding", 30), ("decode", 24), ("trna", 18),
                ("anticodon", 16), ("codon", 12), ("table", 4),
            ]:
                if token in stem_norm:
                    score += pts
            for v in variants:
                if v and v in stem_norm:
                    score += 30
                    break
            return score if score >= 20 else None
        return None

    def _auto_find_cluster_file_for_current_genome(self):
        """Find a preloaded '<species name> clusters.xlsx' workbook for the selected genome."""
        return self._fast_find_preloaded_companion("clusters")

    def _set_refined_cluster_file(self, path, read_sheet_names=True):
        path = str(path or "").strip()
        if not path:
            return
        self.user_cluster_mode.set("Provided by user, xlsx file with 1 column per cluster")
        self.refined_cluster_file.set(path)
        # Do not open Excel workbooks during automatic startup autoloading. On
        # Dropbox/OneDrive this can be surprisingly slow. An empty sheet name is
        # valid: pandas will read the first sheet when the pipeline actually runs.
        if read_sheet_names and path.lower().endswith((".xlsx", ".xls", ".xlsm")):
            try:
                xls = pd.ExcelFile(path)
                if xls.sheet_names:
                    self.refined_cluster_sheet.set(xls.sheet_names[0])
            except Exception:
                pass
        elif not read_sheet_names:
            self.refined_cluster_sheet.set("")
        self.active_clusters_selection = None
        self._refresh_cluster_source_widgets()
        self._update_active_clusters_status()

    def _auto_load_cluster_file_for_current_genome(self):
        """Automatically switch to a preloaded user-provided cluster workbook, if found."""
        auto_cluster = self._auto_find_cluster_file_for_current_genome()
        if auto_cluster:
            self._auto_loaded_cluster_file = auto_cluster
            self._set_refined_cluster_file(auto_cluster, read_sheet_names=False)
            return auto_cluster
        old_auto = getattr(self, "_auto_loaded_cluster_file", "")
        if old_auto and os.path.normcase(self.refined_cluster_file.get().strip()) == os.path.normcase(old_auto):
            self._auto_loaded_cluster_file = ""
            self.refined_cluster_file.set("")
            self.refined_cluster_sheet.set("")
            self.user_cluster_mode.set("Inferred from FASTA file")
            self.active_clusters_selection = None
            self._refresh_cluster_source_widgets()
            self._update_active_clusters_status()
        return ""

    def _auto_find_decoding_table_for_current_genome(self):
        """Find a preloaded '<species name> decoding table.xlsx' workbook for the selected genome."""
        return self._fast_find_preloaded_companion("decoding")

    def _on_decoding_table_path_changed(self):
        path = self.trna_decoding_table_path.get().strip()
        self.enable_trna_usage.set(bool(path))
        self._update_decoding_status()
        self._refresh_simplified_interface()

    def _update_decoding_status(self):
        lbl = getattr(self, "decoding_status_label", None)
        if lbl is None:
            return
        path = self.trna_decoding_table_path.get().strip()
        if path and os.path.isfile(path):
            lbl.configure(text="Decoding strategies enabled. Cluster-level plots will use the decoding-tab cluster picker.")
        elif path:
            lbl.configure(text="Decoding table path is set but the file was not found.")
        else:
            lbl.configure(text="No decoding table loaded. The decoding tab stays hidden in simplified mode.")

    def _browse_trna_decoding_table(self):
        path = filedialog.askopenfilename(
            title="Select decoding-strategy Excel table",
            initialdir=self.default_root.get().strip() or HERE,
            filetypes=[("Excel files", "*.xlsx *.xls *.xlsm"), ("All files", "*.*")],
        )
        if path:
            self.trna_decoding_table_path.set(path)
            try:
                preferred_sheet = _preferred_decoding_sheet_name(path)
                if preferred_sheet:
                    self.trna_decoding_table_sheet.set(preferred_sheet)
                    self.trna_abundance_sheet.set("")
            except Exception:
                pass


    # ------------------------- project state -------------------------
    def _refresh_simplified_interface(self):
        """Show only common tabs, but reveal decoding analyses when a table is loaded."""
        nb = getattr(self, "notebook", None)
        if nb is None:
            return
        advanced = list(getattr(self, "_advanced_tabs", []) or [])
        if not advanced:
            return

        # Companion decoding files are applied by _apply_preloaded_genome_choice().
        # Do not scan folders from this frequently called display-refresh method.
        decoding_loaded = bool(self.trna_decoding_table_path.get().strip())
        if bool(self.simplified_interface.get()):
            for tab, title in advanced:
                is_decoding_tab = (tab is getattr(self, "tab_trna", None))
                try:
                    if is_decoding_tab and decoding_loaded:
                        nb.add(tab, text="Decoding strategies")
                    else:
                        nb.hide(tab)
                except tk.TclError:
                    try:
                        if is_decoding_tab and decoding_loaded:
                            nb.insert("end", tab, text="Decoding strategies")
                    except Exception:
                        pass
                except Exception:
                    pass
        else:
            for tab, title in advanced:
                try:
                    nb.add(tab, text=title)
                except tk.TclError:
                    try:
                        nb.insert("end", tab, text=title)
                    except Exception:
                        pass
                except Exception:
                    pass
        self._update_decoding_status()

    def _iter_top_level_tk_variables(self):
        for name, value in vars(self).items():
            try:
                is_var = isinstance(value, tk.Variable)
            except Exception:
                is_var = False
            if not is_var:
                continue
            yield name, value

    def _collect_project_state(self):
        variables = {}
        for name, var in self._iter_top_level_tk_variables():
            try:
                variables[name] = var.get()
            except Exception:
                pass

        dimred_params = {}
        for method, payload in getattr(self, "dimred_param_vars", {}).items():
            dimred_params[method] = {}
            for key, var in payload.items():
                try:
                    dimred_params[method][key] = var.get()
                except Exception:
                    pass

        cluster_params = {}
        for method, payload in getattr(self, "cluster_param_vars", {}).items():
            cluster_params[method] = {}
            for key, var in payload.items():
                try:
                    cluster_params[method][key] = var.get()
                except Exception:
                    pass

        keyword_groups = {}
        for key, payload in getattr(self, "keyword_group_vars", {}).items():
            keyword_groups[key] = {
                "enabled": bool(payload.get("enabled").get()) if payload.get("enabled") is not None else True,
                "name": payload.get("name").get() if payload.get("name") is not None else str(key),
                "keywords": payload.get("keywords").get() if payload.get("keywords") is not None else "",
            }

        codon_compare_groups = []
        for payload in getattr(self, "codon_compare_groups", []):
            try:
                name = payload["name"].get()
                paths = payload["paths"].get()
                ref_for_statistics = bool(payload.get("ref", tk.BooleanVar(value=False)).get())
            except Exception:
                continue
            # Keep both filled and empty rows so the visual state is preserved.
            codon_compare_groups.append({
                "name": name,
                "paths": paths,
                "ref_for_statistics": ref_for_statistics,
            })

        return {
            "file_type": "CodonPipe GUI project state",
            "version": 1,
            "variables": variables,
            "dimred_param_vars": dimred_params,
            "cluster_param_vars": cluster_params,
            "keyword_group_vars": keyword_groups,
            "active_clusters_selection": self.active_clusters_selection,
            "figure_clusters_selection": self.figure_clusters_selection,
            "decoding_clusters_selection": self.decoding_clusters_selection,
            "decoding_clusters_selection_user_set": bool(getattr(self, "decoding_clusters_selection_user_set", False)),
            "trna_modifications_selection": self.trna_modifications_selection,
            "trna_modification_aas_selection": self.trna_modification_aas_selection,
            "extract_clusters_selection": self.extract_clusters_selection,
            "codon_compare_groups": codon_compare_groups,
            "codon_compare_selected_codons": list(getattr(self, "codon_compare_selected_codons", list(_CODON_ORDER))),
            "codon_compare_highlighted_codons": list(getattr(self, "codon_compare_highlighted_codons", [])),
            "fasta_metric_cluster_configs": self._collect_fasta_metric_cluster_configs(),
            "last_clustering_workbook": getattr(self, "last_clustering_workbook", ""),
            "last_clustering_output_dir": getattr(self, "last_clustering_output_dir", ""),
        }

    def _safe_set_tk_var(self, var, value):
        try:
            if isinstance(var, tk.BooleanVar):
                if isinstance(value, str):
                    var.set(value.strip().lower() in {"1", "true", "yes", "y", "on"})
                else:
                    var.set(bool(value))
            else:
                var.set("" if value is None else value)
        except Exception:
            pass

    def _apply_project_state(self, state):
        if not isinstance(state, dict):
            raise ValueError("The selected file is not a valid CodonPipe project state JSON file.")

        variables = state.get("variables", {}) or {}
        for name, value in variables.items():
            var = getattr(self, name, None)
            try:
                is_var = isinstance(var, tk.Variable)
            except Exception:
                is_var = False
            if is_var:
                self._safe_set_tk_var(var, value)

        for method, payload in (state.get("dimred_param_vars", {}) or {}).items():
            if method not in getattr(self, "dimred_param_vars", {}):
                continue
            for key, value in (payload or {}).items():
                var = self.dimred_param_vars[method].get(key)
                if var is not None:
                    self._safe_set_tk_var(var, value)

        for method, payload in (state.get("cluster_param_vars", {}) or {}).items():
            if method not in getattr(self, "cluster_param_vars", {}):
                continue
            for key, value in (payload or {}).items():
                var = self.cluster_param_vars[method].get(key)
                if var is not None:
                    self._safe_set_tk_var(var, value)

        for key, payload in (state.get("keyword_group_vars", {}) or {}).items():
            current = getattr(self, "keyword_group_vars", {}).get(key)
            if not current:
                continue
            if "enabled" in payload:
                self._safe_set_tk_var(current.get("enabled"), payload.get("enabled"))
            if "name" in payload:
                self._safe_set_tk_var(current.get("name"), payload.get("name"))
            if "keywords" in payload:
                self._safe_set_tk_var(current.get("keywords"), payload.get("keywords"))

        self.active_clusters_selection = state.get("active_clusters_selection", None)
        self.figure_clusters_selection = state.get("figure_clusters_selection", None)
        self.decoding_clusters_selection = state.get("decoding_clusters_selection", [])
        self.decoding_clusters_selection_user_set = bool(state.get("decoding_clusters_selection_user_set", False))
        self.trna_modifications_selection = _canonicalize_gui_plot6_modification_selection(state.get("trna_modifications_selection", None))
        self.trna_modification_aas_selection = _split_selection_text(state.get("trna_modification_aas_selection", None))
        self.extract_clusters_selection = state.get("extract_clusters_selection", None)
        self.last_clustering_workbook = str(state.get("last_clustering_workbook", "") or "")
        self.last_clustering_output_dir = str(state.get("last_clustering_output_dir", "") or "")
        self._set_codon_compare_groups_from_state(state.get("codon_compare_groups", []))
        saved_codons = [str(c).upper().replace("U", "T") for c in list(state.get("codon_compare_selected_codons", []) or [])]
        self.codon_compare_selected_codons = [c for c in saved_codons if c in _CODON_ORDER] or list(_CODON_ORDER)
        saved_highlights = [str(c).upper().replace("U", "T") for c in list(state.get("codon_compare_highlighted_codons", []) or [])]
        self.codon_compare_highlighted_codons = [c for c in saved_highlights if c in _CODON_ORDER]
        self._set_fasta_metric_group_rows_from_configs(state.get("fasta_metric_cluster_configs", []))

        # Keep dependent widgets and internal settings synchronized after loading.
        if not bool(self.genome_not_available.get()):
            self._apply_preloaded_genome_choice()
        self._sync_codon_set_from_usage_basis()
        self._sync_analysis_choice_booleans()
        self._refresh_genome_source_widgets()
        self._refresh_cluster_source_widgets()
        self._refresh_dimred_params()
        self._refresh_cluster_params()
        self._refresh_optional_sections()
        self._refresh_keyword_group_visibility()
        self._refresh_codon_compare_metric_widgets()
        self._update_codon_compare_selected_codons_status()
        self._update_codon_compare_highlight_status()
        self._update_fasta_metric_cluster_status()
        self._update_active_clusters_status()
        self._update_trna_modification_selection_status()
        self._update_trna_modification_aas_status()
        self._update_extract_clusters_status()
        self._refresh_simplified_interface()

    def _set_codon_compare_groups_from_state(self, groups):
        if not hasattr(self, "codon_compare_rows_frame"):
            return
        for payload in list(getattr(self, "codon_compare_groups", [])):
            try:
                payload.get("frame").destroy()
            except Exception:
                pass
        self.codon_compare_groups = []
        groups = [g for g in (groups or []) if isinstance(g, dict)]
        if not groups:
            self._add_codon_compare_group_row()
            return
        for g in groups:
            self._add_codon_compare_group_row(
                name=g.get("name", ""),
                paths=g.get("paths", ""),
                ref_for_statistics=bool(g.get("ref_for_statistics", g.get("ref", False))),
            )
        # Keep old sessions robust: if multiple references were saved, keep only the first.
        seen_ref = False
        for payload in self.codon_compare_groups:
            try:
                ref_var = payload.get("ref")
                if bool(ref_var.get()):
                    if seen_ref:
                        ref_var.set(False)
                    else:
                        seen_ref = True
            except Exception:
                pass
        # Ensure one blank row remains available for adding a new comparison group.
        if self.codon_compare_groups:
            last_paths = self.codon_compare_groups[-1]["paths"].get().strip()
            if last_paths:
                self._add_codon_compare_group_row()

    def _save_project_state_dialog(self):
        initial_dir = self.default_root.get().strip() or HERE
        path = filedialog.asksaveasfilename(
            initialdir=initial_dir,
            title="Save CodonPipe project state",
            defaultextension=".json",
            initialfile="CodonPipe_session.json",
            filetypes=[("CodonPipe session JSON", "*.json"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            state = self._collect_project_state()
            with open(path, "w", encoding="utf-8") as fh:
                json.dump(state, fh, indent=2, ensure_ascii=False)
            self._append_log(f"[INFO] Saved CodonPipe project state: {path}\n")
            messagebox.showinfo("Project state saved", f"Saved project state:\n{path}")
        except Exception as e:
            messagebox.showerror("Project state", f"Could not save project state:\n{e}")

    def _load_project_state_dialog(self):
        initial_dir = self.default_root.get().strip() or HERE
        path = filedialog.askopenfilename(
            initialdir=initial_dir,
            title="Load CodonPipe project state",
            filetypes=[("CodonPipe session JSON", "*.json"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            with open(path, "r", encoding="utf-8") as fh:
                state = json.load(fh)
            self._apply_project_state(state)
            self._append_log(f"[INFO] Loaded CodonPipe project state: {path}\n")
            messagebox.showinfo("Project state loaded", f"Loaded project state:\n{path}")
        except Exception as e:
            messagebox.showerror("Project state", f"Could not load project state:\n{e}")

    def _build_right_panel(self, parent):
        log_box = ttk.LabelFrame(parent, text="Console log", padding=8)
        log_box.pack(fill="both", expand=True)

        self.log_text = tk.Text(log_box, wrap="word", font=("Consolas", 10), background="#111111", foreground="#f0f0f0")
        self.log_text.pack(side="left", fill="both", expand=True)
        scroll = ttk.Scrollbar(log_box, orient="vertical", command=self.log_text.yview)
        scroll.pack(side="right", fill="y")
        self.log_text.configure(yscrollcommand=scroll.set)

    # ------------------------- helpers -------------------------
    def _tip(self, widget, text):
        if widget is None or not text:
            return
        try:
            tip = ToolTip(widget, text, self.tooltips_enabled)
            self.tooltip_refs.append(tip)
        except Exception:
            pass

    def _add_file_row(self, parent, row, label, var, filetypes=None, ask_dir=False, tip=""):
        parent.columnconfigure(1, weight=1)
        lbl = ttk.Label(parent, text=label)
        lbl.grid(row=row, column=0, sticky="w", padx=(0, 8), pady=4)
        self._tip(lbl, tip)
        ent = ttk.Entry(parent, textvariable=var)
        ent.grid(row=row, column=1, sticky="ew", pady=4)
        self._tip(ent, tip)
        if ask_dir:
            cmd = lambda v=var: self._browse_directory(v)
        else:
            cmd = lambda v=var, ft=filetypes: self._browse_file(v, ft)
        btn = ttk.Button(parent, text="Browse", command=cmd)
        btn.grid(row=row, column=2, sticky="w", padx=(8, 0), pady=4)
        self._tip(btn, tip)

    def _add_labeled_widget(self, parent, row, col, label_text, widget, tip=""):
        lbl = ttk.Label(parent, text=label_text)
        lbl.grid(row=row, column=col, sticky="w", padx=(0, 8), pady=4)
        widget.grid(row=row, column=col + 1, sticky="w", pady=4)
        self._tip(lbl, tip)
        self._tip(widget, tip)

    def _add_param_widget(self, parent, row, spec, var):
        key, label, widget_type, options, tip = spec
        lbl = ttk.Label(parent, text=label)
        lbl.grid(row=row, column=0, sticky="w", padx=(0, 8), pady=3)
        self._tip(lbl, tip)
        if widget_type == "entry":
            w = ttk.Entry(parent, textvariable=var, width=18)
        elif widget_type == "combo":
            w = ttk.Combobox(parent, textvariable=var, state="readonly", values=options or [], width=16)
        elif widget_type == "check":
            w = ttk.Checkbutton(parent, variable=var)
        else:
            w = ttk.Entry(parent, textvariable=var, width=18)
        w.grid(row=row, column=1, sticky="w", pady=3)
        self._tip(w, tip)
        return w

    def _on_cluster_mode_changed(self):
        self._auto_loaded_cluster_file = ""
        self.active_clusters_selection = None
        self._refresh_cluster_source_widgets()
        self._update_active_clusters_status()

    def _refresh_cluster_source_widgets(self):
        mode = self.user_cluster_mode.get().strip()
        is_refined = _cluster_mode_internal(mode) == "refined"
        is_david_terms = _cluster_mode_terms_source(mode) == "david_gene2terms"

        for widget in (
            self.david_terms_label, self.david_terms_entry, self.david_terms_button,
            self.refined_label, self.refined_entry, self.refined_button,
            self.refined_sheet_label, self.refined_sheet_entry,
        ):
            try:
                widget.grid_remove()
            except Exception:
                pass

        if is_david_terms:
            self.david_terms_label.grid(row=1, column=0, sticky="w", padx=(0, 8), pady=4)
            self.david_terms_entry.grid(row=1, column=1, sticky="ew", pady=4)
            self.david_terms_button.grid(row=1, column=2, sticky="w", padx=(8, 0), pady=4)

        if is_refined:
            self.refined_label.grid(row=1, column=0, sticky="w", padx=(0, 8), pady=4)
            self.refined_entry.grid(row=1, column=1, sticky="ew", pady=4)
            self.refined_button.grid(row=1, column=2, sticky="w", padx=(8, 0), pady=4)
            self.refined_sheet_label.grid(row=2, column=0, sticky="w", padx=(0, 8), pady=4)
            self.refined_sheet_entry.grid(row=2, column=1, sticky="ew", pady=4)

    def _refresh_dimred_params(self):
        for w in self.dimred_params_box.winfo_children():
            w.destroy()
        method = _dimred_internal(self.dimred_method.get())
        specs = self.dimred_specs.get(method, [])
        vars_dict = self.dimred_param_vars.get(method, {})
        for i, spec in enumerate(specs):
            self._add_param_widget(self.dimred_params_box, i, spec, vars_dict[spec[0]])

    def _refresh_cluster_params(self):
        for w in self.cluster_params_box.winfo_children():
            w.destroy()
        method = self.cluster_method.get().strip().lower()
        specs = self.cluster_specs.get(method, [])
        vars_dict = self.cluster_param_vars.get(method, {})
        for i, spec in enumerate(specs):
            self._add_param_widget(self.cluster_params_box, i, spec, vars_dict[spec[0]])

    def _set_section_enabled(self, frame, enabled):
        if enabled:
            frame.grid()
        else:
            frame.grid_remove()

    def _refresh_optional_sections(self):
        if hasattr(self, "custom_cds_content"):
            self._set_section_enabled(self.custom_cds_content, bool(self.add_custom_cds_enable.get()))
        if hasattr(self, "ks_box"):
            self._set_section_enabled(self.ks_box, bool(self.enable_2d_ks.get()))
        if hasattr(self, "david_box"):
            self._set_section_enabled(self.david_box, True)
        if hasattr(self, "main_heatmap_details_box"):
            self._set_section_enabled(self.main_heatmap_details_box, bool(self.enable_main_heatmap.get()))
        if hasattr(self, "main_heatmap_aesthetics_content"):
            self._set_section_enabled(self.main_heatmap_aesthetics_content, bool(self.enable_main_heatmap.get() and self.main_heatmap_custom_aesthetics.get()))
        if hasattr(self, "main_heatmap_axes_content"):
            self._set_section_enabled(self.main_heatmap_axes_content, bool(self.enable_main_heatmap.get() and self.main_heatmap_custom_axes.get()))
        if hasattr(self, "heatmap_section"):
            self._set_section_enabled(self.heatmap_section, bool(self.gchm_enable.get()))
        if hasattr(self, "heatmap_aesthetics_content"):
            self._set_section_enabled(self.heatmap_aesthetics_content, bool(self.gchm_enable.get() and self.gchm_custom_aesthetics.get()))
        if hasattr(self, "heatmap_axes_content"):
            self._set_section_enabled(self.heatmap_axes_content, bool(self.gchm_enable.get() and self.gchm_custom_axes.get()))
        if hasattr(self, "density_section"):
            self._set_section_enabled(self.density_section, bool(self.enable_2d_density_plots.get()))
        if hasattr(self, "density_aesthetics_content"):
            self._set_section_enabled(self.density_aesthetics_content, bool(self.enable_2d_density_plots.get() and self.density_custom_aesthetics.get()))
        if hasattr(self, "density_axes_content"):
            self._set_section_enabled(self.density_axes_content, bool(self.enable_2d_density_plots.get() and self.density_custom_axes.get()))
        if hasattr(self, "codon_section"):
            self._set_section_enabled(self.codon_section, bool(self.enable_codon_usage_plot.get()))
        if hasattr(self, "codon_aesthetics_content"):
            self._set_section_enabled(self.codon_aesthetics_content, bool(self.enable_codon_usage_plot.get() and self.codon_custom_aesthetics.get()))
        if hasattr(self, "codon_axes_content"):
            self._set_section_enabled(self.codon_axes_content, bool(self.enable_codon_usage_plot.get() and self.codon_custom_axes.get()))
        if hasattr(self, "codon_compare_raw_axes_content"):
            self._set_section_enabled(self.codon_compare_raw_axes_content, bool(self.codon_compare_raw_custom_axes.get()))
        if hasattr(self, "codon_compare_corr_axes_content"):
            self._set_section_enabled(self.codon_compare_corr_axes_content, bool(self.codon_compare_corr_custom_axes.get()))
        if hasattr(self, "trna_content"):
            self._set_section_enabled(self.trna_content, bool(self.enable_trna_usage.get()))
        if hasattr(self, "trna_input_box"):
            self._set_section_enabled(self.trna_input_box, bool(self.enable_trna_usage.get()))
        if hasattr(self, "trna_gene_heatmap_box"):
            self._set_section_enabled(self.trna_gene_heatmap_box, bool(self.enable_trna_usage.get()))
        if hasattr(self, "trna_wobble_gene_box"):
            self._set_section_enabled(self.trna_wobble_gene_box, bool(self.enable_trna_usage.get()))
        if hasattr(self, "trna_shift_heatmap_box"):
            self._set_section_enabled(self.trna_shift_heatmap_box, bool(self.enable_trna_usage.get()))
        if hasattr(self, "trna_wobble_shift_box"):
            self._set_section_enabled(self.trna_wobble_shift_box, bool(self.enable_trna_usage.get()))
        if hasattr(self, "trna_modification_box"):
            self._set_section_enabled(self.trna_modification_box, bool(self.enable_trna_usage.get()))
        if hasattr(self, "decoding_details_box"):
            self._set_section_enabled(self.decoding_details_box, bool(self.enable_trna_usage.get()))
        if hasattr(self, "trna_gene_heatmap_custom_content"):
            self._set_section_enabled(self.trna_gene_heatmap_custom_content, bool(self.enable_trna_usage.get() and self.trna_supp_heatmaps_customize.get()))
        if hasattr(self, "trna_shift_heatmaps_custom_content"):
            self._set_section_enabled(self.trna_shift_heatmaps_custom_content, bool(self.enable_trna_usage.get() and self.trna_shift_heatmaps_customize.get()))
        if hasattr(self, "trna_modification_plots_details_content"):
            self._set_section_enabled(self.trna_modification_plots_details_content, bool(self.enable_trna_usage.get()))
        if hasattr(self, "trna_modification_plots_custom_content"):
            self._set_section_enabled(self.trna_modification_plots_custom_content, bool(self.enable_trna_usage.get()))
        if hasattr(self, "trna_secondary_axis_details_content"):
            self._set_section_enabled(self.trna_secondary_axis_details_content, bool(self.enable_trna_usage.get()))
        self._refresh_size_lock_widgets()

    def _set_widget_enabled_state(self, widget, enabled):
        if widget is None:
            return
        try:
            widget.configure(state=("normal" if enabled else "disabled"))
        except Exception:
            pass

    def _make_exclusive_pair_callback(self, source_name):
        def _callback(*_args):
            if getattr(self, "_size_pair_guard", False):
                return
            self._size_pair_guard = True
            try:
                if source_name == "trna_supp_fig_h":
                    if self.trna_supp_heatmaps_fig_height.get().strip() and self.trna_supp_heatmaps_cell_height.get().strip():
                        self.trna_supp_heatmaps_cell_height.set("")
                elif source_name == "trna_supp_cell_h":
                    if self.trna_supp_heatmaps_cell_height.get().strip() and self.trna_supp_heatmaps_fig_height.get().strip():
                        self.trna_supp_heatmaps_fig_height.set("")
                elif source_name == "trna_shift_fig_w":
                    if self.trna_shift_heatmaps_fig_width.get().strip() and self.trna_shift_heatmaps_cell_width.get().strip():
                        self.trna_shift_heatmaps_cell_width.set("")
                elif source_name == "trna_shift_cell_w":
                    if self.trna_shift_heatmaps_cell_width.get().strip() and self.trna_shift_heatmaps_fig_width.get().strip():
                        self.trna_shift_heatmaps_fig_width.set("")
                elif source_name == "trna_shift_fig_h":
                    if self.trna_shift_heatmaps_fig_height.get().strip() and self.trna_shift_heatmaps_cell_height.get().strip():
                        self.trna_shift_heatmaps_cell_height.set("")
                elif source_name == "trna_shift_cell_h":
                    if self.trna_shift_heatmaps_cell_height.get().strip() and self.trna_shift_heatmaps_fig_height.get().strip():
                        self.trna_shift_heatmaps_fig_height.set("")
            finally:
                self._size_pair_guard = False
            self._refresh_size_lock_widgets()
        return _callback

    def _bind_size_lock_pairs(self):
        self.trna_supp_heatmaps_fig_height.trace_add("write", self._make_exclusive_pair_callback("trna_supp_fig_h"))
        self.trna_supp_heatmaps_cell_height.trace_add("write", self._make_exclusive_pair_callback("trna_supp_cell_h"))
        self.trna_shift_heatmaps_fig_width.trace_add("write", self._make_exclusive_pair_callback("trna_shift_fig_w"))
        self.trna_shift_heatmaps_cell_width.trace_add("write", self._make_exclusive_pair_callback("trna_shift_cell_w"))
        self.trna_shift_heatmaps_fig_height.trace_add("write", self._make_exclusive_pair_callback("trna_shift_fig_h"))
        self.trna_shift_heatmaps_cell_height.trace_add("write", self._make_exclusive_pair_callback("trna_shift_cell_h"))
        self._refresh_size_lock_widgets()

    def _refresh_size_lock_widgets(self):
        supp_active = bool(self.enable_trna_usage.get() and self.trna_supp_heatmaps_customize.get())
        shift_active = bool(self.enable_trna_usage.get() and self.trna_shift_heatmaps_customize.get())

        supp_fig_h_filled = bool(self.trna_supp_heatmaps_fig_height.get().strip())
        supp_cell_h_filled = bool(self.trna_supp_heatmaps_cell_height.get().strip())
        self._set_widget_enabled_state(getattr(self, "trna_supp_fig_height_entry", None), supp_active and not supp_cell_h_filled)
        self._set_widget_enabled_state(getattr(self, "trna_supp_cell_height_entry", None), supp_active and not supp_fig_h_filled)

        shift_fig_w_filled = bool(self.trna_shift_heatmaps_fig_width.get().strip())
        shift_cell_w_filled = bool(self.trna_shift_heatmaps_cell_width.get().strip())
        shift_fig_h_filled = bool(self.trna_shift_heatmaps_fig_height.get().strip())
        shift_cell_h_filled = bool(self.trna_shift_heatmaps_cell_height.get().strip())
        self._set_widget_enabled_state(getattr(self, "trna_shift_fig_width_entry", None), shift_active and not shift_cell_w_filled)
        self._set_widget_enabled_state(getattr(self, "trna_shift_cell_width_entry", None), shift_active and not shift_fig_w_filled)
        self._set_widget_enabled_state(getattr(self, "trna_shift_fig_height_entry", None), shift_active and not shift_cell_h_filled)
        self._set_widget_enabled_state(getattr(self, "trna_shift_cell_height_entry", None), shift_active and not shift_fig_h_filled)

    def _export_readme_setup(self):
        target_dir = filedialog.askdirectory(
            initialdir=self.default_root.get() or HERE,
            title="Choose folder for readme/setup export",
        )
        if not target_dir:
            return
        try:
            paths = []
            readme_path = os.path.join(target_dir, "CodonPipe_Readme.txt")
            with open(readme_path, "w", encoding="utf-8") as fh:
                fh.write(self._readme_text())
            paths.append(readme_path)

            deps_path = os.path.join(target_dir, "CodonPipe_setup_and_dependencies.txt")
            with open(deps_path, "w", encoding="utf-8") as fh:
                fh.write(self._setup_dependencies_text())
            paths.append(deps_path)

            quick_path = os.path.join(target_dir, "CodonPipe_Quick_start.txt")
            with open(quick_path, "w", encoding="utf-8") as fh:
                fh.write(self._quick_start_text())
            paths.append(quick_path)

            messagebox.showinfo(
                "CodonPipe",
                "Exported helper files to:\n\n" + "\n".join(paths),
            )
        except Exception as e:
            messagebox.showerror(
                "CodonPipe",
                f"Could not export helper files:\n{e}",
            )

    def _readme_text(self):
        return """CodonPipe desktop interface
============================

This folder contains the graphical launcher for the CodonPipe clustering pipeline.

Expected local layout
---------------------
Put these files in the same folder:
- CodonPipe_GUI.py
- Clustering_Pipeline.py
- Plotting_Pipeline.py

And keep the codonpipe package folder next to them, containing modules such as:
- clustering.py
- fasta_metrics.py
- excel_outputs.py
- density_bridge.py
- density_plot_core.py
- gene_cluster_heatmap.py
- david_window_scan.py
- legend.py
- ks2d.py

Important notes
---------------
- Run the GUI in the same conda environment where the pipeline already works.
- The GUI suppresses most console prompts by sending your choices directly into the pipeline.
- Optional sections in the GUI stay hidden until enabled, to keep the interface lighter for new users.
"""

    def _setup_dependencies_text(self):
        return """CodonPipe setup and dependencies
===============================

Recommended environment
-----------------------
Use the same conda environment that already runs your CodonPipe pipeline.

Typical Python modules used by the GUI and pipeline
---------------------------------------------------
python
numpy
pandas
matplotlib
scipy
scikit-learn
umap-learn
openpyxl
xlsxwriter
threadpoolctl
joblib
tkinter (usually included with Python)
suds-community   [needed only if DAVID web-service access is used]

Example conda / pip commands
---------------------------
conda install numpy pandas matplotlib scipy scikit-learn openpyxl xlsxwriter joblib threadpoolctl
conda install -c conda-forge umap-learn
pip install suds-community

If you do not plan to use DAVID, suds-community is not mandatory.
"""

    def _quick_start_text(self):
        return """CodonPipe quick start
=====================

1. Place CodonPipe_GUI.py, Clustering_Pipeline.py and Plotting_Pipeline.py in the same folder.
2. Make sure the codonpipe package folder is next to them.
3. Open the correct conda environment.
4. Run CodonPipe_GUI.py from Spyder or with: python CodonPipe_GUI.py
5. In the Input/Output tab, choose a preloaded genome or your own CDS FASTA and an export root.
6. Optionally add custom CDS to the analysis or compare separate FASTA groups against the selected genome.
7. In Gene clusters, choose inferred clusters or provide an xlsx file with one column per cluster.
8. In Codon usage clustering, choose the codon metric, dimensional reduction, clustering, statistics, and functional scan options.
9. In Figures, enable only the plots you want to generate.
10. Use the buttons at the bottom of the Codon usage clustering, Figures, decoding strategy analyses, or Codon usage analyses tabs depending on the task.

Fastest first test
------------------
Use a preloaded genome or your own CDS FASTA file to verify that the GUI, pipeline, and plotting scripts are all wired correctly on your computer.
"""

    def _browse_file(self, var, filetypes=None):
        path = filedialog.askopenfilename(
            initialdir=self.default_root.get() or HERE,
            filetypes=filetypes or [("All files", "*.*")],
        )
        if path:
            var.set(path)

    def _browse_directory(self, var):
        path = filedialog.askdirectory(initialdir=var.get() or HERE)
        if path:
            var.set(path)

    def _browse_refined_cluster(self):
        path = filedialog.askopenfilename(
            initialdir=self.default_root.get() or HERE,
            filetypes=[("Cluster files", "*.txt *.csv *.tsv *.xlsx *.xls *.xlsm"), ("All files", "*.*")],
        )
        if path:
            self._auto_loaded_cluster_file = ""
            self._set_refined_cluster_file(path)


    def _browse_david_gene2terms_file(self):
        path = filedialog.askopenfilename(
            initialdir=self.default_root.get() or HERE,
            title="Select DAVID gene2terms TXT file",
            filetypes=[
                ("DAVID gene2terms files", "*.txt *.tsv *.csv"),
                ("Text files", "*.txt *.tsv *.csv"),
                ("All files", "*.*"),
            ],
        )
        if path:
            self.david_gene2terms_path.set(path)
            self.user_cluster_mode.set("Inferred from DAVID gene2terms txt file")
            self.active_clusters_selection = None
            self._refresh_cluster_source_widgets()
            self._update_active_clusters_status()

    def _browse_custom_cds_files(self):
        paths = filedialog.askopenfilenames(
            initialdir=self.default_root.get() or HERE,
            title="Select custom CDS FASTA/.fna file(s)",
            filetypes=[
                ("FASTA/.fna CDS files", "*.fna *.ffn *.fa *.fas *.fasta"),
                ("FNA files", "*.fna"),
                ("All files", "*.*"),
            ],
        )
        if paths:
            self.custom_cds_paths.set("; ".join(paths))
            self.add_custom_cds_enable.set(True)
            if not self.custom_cds_cluster_name.get().strip():
                self.custom_cds_cluster_name.set("custom")
            self._refresh_optional_sections()
            self._mark_active_clusters_dirty()

    def _add_codon_compare_group_row(self, name="", paths="", ref_for_statistics=False):
        if not hasattr(self, "codon_compare_rows_frame"):
            return
        idx = len(self.codon_compare_groups)
        row = ttk.Frame(self.codon_compare_rows_frame)
        # Row 0 is the compact header. Data rows start at row 1.
        row.grid(row=idx + 1, column=0, columnspan=8, sticky="ew", pady=3)
        row.columnconfigure(2, weight=1)

        default_name = str(name or (self._suggest_unique_compare_group_name_from_paths(paths) if paths else "")).strip()
        ref_var = tk.BooleanVar(value=bool(ref_for_statistics))
        name_var = tk.StringVar(value=default_name)
        paths_var = tk.StringVar(value=str(paths or ""))
        payload = {"frame": row, "ref": ref_var, "name": name_var, "paths": paths_var}
        self.codon_compare_groups.append(payload)

        ref_chk = ttk.Checkbutton(row, text="", variable=ref_var, command=lambda p=payload: self._on_codon_compare_reference_changed(p))
        ref_chk.grid(row=0, column=0, sticky="w", padx=(10, 12))
        name_entry = ttk.Entry(row, textvariable=name_var, width=22)
        name_entry.grid(row=0, column=1, sticky="w", padx=(0, 10))
        path_entry = ttk.Entry(row, textvariable=paths_var)
        path_entry.grid(row=0, column=2, sticky="ew", padx=(0, 8))
        btn = ttk.Button(row, text="Load FASTA", command=lambda p=payload: self._browse_codon_compare_files(p))
        btn.grid(row=0, column=3, sticky="w", padx=(0, 6))
        cluster_btn = ttk.Button(row, text="Load cluster", command=lambda p=payload: self._load_cluster_as_codon_compare_group(p))
        cluster_btn.grid(row=0, column=4, sticky="w", padx=(0, 6))
        kill_btn = ttk.Button(row, text="×", width=3, command=lambda p=payload: self._remove_codon_compare_group_row(p))
        kill_btn.grid(row=0, column=5, sticky="w")
        payload["remove_button"] = kill_btn
        payload["cluster_button"] = cluster_btn
        payload["ref_checkbox"] = ref_chk
        self._tip(ref_chk, "Tick one group as the reference for statistical tests. Only one reference can be active; selecting this row automatically clears any other reference checkbox. Correlation plots also keep this reference on the X axis.")
        self._tip(name_entry, "Custom label shown in the plot legend and in correlation axis titles.")
        self._tip(path_entry, "One or more FASTA/.fna files, or a cluster converted to FASTA. Files in the same row are analyzed together as one gene group.")
        self._tip(btn, "Select one or more CDS FASTA/.fna files for this comparison group.")
        self._tip(cluster_btn, "Select one or more currently available gene clusters and load their CDS as this comparison group.")
        self._tip(kill_btn, "Delete this comparison group.")

    def _on_codon_compare_reference_changed(self, selected_payload):
        """Keep at most one comparison group checked as statistics reference."""
        try:
            selected_ref = selected_payload.get("ref")
            if selected_ref is None or not bool(selected_ref.get()):
                return
        except Exception:
            return
        for payload in list(getattr(self, "codon_compare_groups", [])):
            if payload is selected_payload:
                continue
            try:
                payload.get("ref").set(False)
            except Exception:
                pass

    def _refresh_codon_compare_group_rows(self):
        """Re-grid comparison rows and guarantee one empty loader row remains."""
        if not hasattr(self, "codon_compare_rows_frame"):
            return
        for idx, payload in enumerate(list(getattr(self, "codon_compare_groups", []))):
            frame = payload.get("frame")
            try:
                frame.grid(row=idx + 1, column=0, columnspan=8, sticky="ew", pady=3)
            except Exception:
                pass
        if not getattr(self, "codon_compare_groups", []):
            self._add_codon_compare_group_row()
            return
        try:
            has_blank = any(not p["paths"].get().strip() for p in self.codon_compare_groups)
        except Exception:
            has_blank = False
        if not has_blank:
            self._add_codon_compare_group_row()

    def _remove_codon_compare_group_row(self, group_payload):
        try:
            if group_payload in self.codon_compare_groups:
                self.codon_compare_groups.remove(group_payload)
        except Exception:
            pass
        try:
            group_payload.get("frame").destroy()
        except Exception:
            pass
        self._refresh_codon_compare_group_rows()

    def _current_compare_group_names(self, exclude_payload=None):
        names = []
        for p in getattr(self, "codon_compare_groups", []):
            if exclude_payload is not None and p is exclude_payload:
                continue
            try:
                name = p["name"].get().strip()
            except Exception:
                name = ""
            if name:
                names.append(name)
        return names

    def _suggest_unique_compare_group_name_from_paths(self, paths, exclude_payload=None):
        label = _suggest_group_name_from_paths(paths, fallback="FASTA group", max_len=42)
        return _make_unique_label(label, self._current_compare_group_names(exclude_payload=exclude_payload), fallback="FASTA group")

    def _browse_codon_compare_files(self, group_payload):
        paths = filedialog.askopenfilenames(
            initialdir=self.default_root.get() or HERE,
            title="Select comparison CDS FASTA/.fna file(s)",
            filetypes=[
                ("FASTA/.fna CDS files", "*.fna *.ffn *.fa *.fas *.fasta"),
                ("FNA files", "*.fna"),
                ("All files", "*.*"),
            ],
        )
        if not paths:
            return
        joined = "; ".join(paths)
        group_payload["paths"].set(joined)
        suggested = self._suggest_unique_compare_group_name_from_paths(paths, exclude_payload=group_payload)
        # Always refresh the suggested name after file loading. Users can still edit it afterwards.
        group_payload["name"].set(suggested)
        # Always keep one blank row at the bottom for the next comparison group.
        if self.codon_compare_groups and self.codon_compare_groups[-1] is group_payload:
            self._add_codon_compare_group_row()

    def _choose_clusters_for_compare_dialog(self):
        available = self._get_available_cluster_names(show_errors=True)
        if not available:
            messagebox.showinfo("Load cluster", "No clusters are currently available. Check the selected cluster source and keyword groups.")
            return None

        win = tk.Toplevel(self.root)
        win.title("Load cluster as codon-usage group")
        win.geometry("560x640")
        win.transient(self.root)
        win.grab_set()

        outer = ttk.Frame(win, padding=10)
        outer.pack(fill="both", expand=True)
        ttk.Label(
            outer,
            text="Select one or more clusters to load as a single codon-usage comparison group.",
            wraplength=520,
            foreground="#444444",
        ).pack(anchor="w", pady=(0, 8))

        tools = ttk.Frame(outer)
        tools.pack(fill="x", pady=(0, 6))
        vars_by_name = {name: tk.BooleanVar(value=False) for name in available}

        def set_all(value):
            for v in vars_by_name.values():
                v.set(bool(value))

        ttk.Button(tools, text="Select all", command=lambda: set_all(True)).pack(side="left")
        ttk.Button(tools, text="Select none", command=lambda: set_all(False)).pack(side="left", padx=(6, 0))

        scroller = ScrollableFrame(outer)
        scroller.pack(fill="both", expand=True)
        for name in available:
            ttk.Checkbutton(scroller.inner, text=name, variable=vars_by_name[name]).pack(anchor="w", pady=2)

        result = {"selected": None}
        bottom = ttk.Frame(outer)
        bottom.pack(fill="x", pady=(8, 0))

        def apply_selection():
            chosen = [name for name in available if vars_by_name[name].get()]
            if not chosen:
                messagebox.showwarning("Load cluster", "Please select at least one cluster.")
                return
            result["selected"] = chosen
            win.destroy()

        ttk.Button(bottom, text="Load selected cluster(s)", command=apply_selection).pack(side="left")
        ttk.Button(bottom, text="Cancel", command=win.destroy).pack(side="left", padx=(8, 0))
        self.root.wait_window(win)
        return result["selected"]

    def _load_cluster_as_codon_compare_group(self, group_payload):
        try:
            selected_clusters = self._choose_clusters_for_compare_dialog()
            if not selected_clusters:
                return
            fasta_path, label, n_found, n_missing = self._write_clusters_to_comparison_fasta(selected_clusters)
            group_payload["paths"].set(fasta_path)
            group_payload["name"].set(_make_unique_label(label, self._current_compare_group_names(exclude_payload=group_payload), fallback="Cluster group"))
            if self.codon_compare_groups and self.codon_compare_groups[-1] is group_payload:
                self._add_codon_compare_group_row()
            msg = f"Loaded {n_found} CDS from {len(selected_clusters)} cluster(s) as comparison group: {group_payload['name'].get()}"
            if n_missing:
                msg += f"\nWarning: {n_missing} cluster ID(s) were not found in the current FASTA and were skipped."
            self._append_log("\n[GUI] " + msg.replace("\n", "\n[GUI] ") + "\n")
        except Exception as e:
            messagebox.showerror("Load cluster failed", str(e))

    def _write_clusters_to_comparison_fasta(self, selected_clusters):
        cluster_df = self._build_current_cluster_df_for_tools()
        available_cols = {str(c).strip().lower(): str(c) for c in cluster_df.columns if str(c).strip()}
        resolved_clusters = []
        for cname in list(selected_clusters or []):
            mapped = available_cols.get(str(cname).strip().lower())
            if mapped is not None and mapped not in resolved_clusters:
                resolved_clusters.append(mapped)
        if not resolved_clusters:
            raise ValueError("No selected cluster was found in the current cluster table.")

        seq_records = self._current_sequence_records_for_export()
        if not seq_records:
            raise ValueError("No CDS sequences were found in the selected genome FASTA.")

        members = []
        for cname in resolved_clusters:
            vals = cluster_df[cname].fillna("").astype(str).str.strip().tolist()
            vals = [v for v in vals if v and v.lower() != "nan"]
            members.extend(vals)
        members = list(dict.fromkeys(members))
        found = [m for m in members if m in seq_records]
        missing = [m for m in members if m not in seq_records]
        if not found:
            raise ValueError("The selected cluster(s) contained no IDs matching the current genome FASTA records.")

        out_root = self._codon_usage_analysis_output_dir()
        out_dir = os.path.join(out_root, "Raw codon usage tables", "Cluster comparison FASTA")
        os.makedirs(out_dir, exist_ok=True)
        label = " + ".join(resolved_clusters)
        if len(resolved_clusters) > 1:
            file_base = "Multiple clusters"
        else:
            file_base = resolved_clusters[0]
        safe_base = _safe_output_name(file_base, fallback="Cluster group", max_len=80)
        out_path = os.path.join(out_dir, safe_base + ".fna")
        suffix = 2
        while os.path.exists(out_path):
            out_path = os.path.join(out_dir, f"{safe_base}_{suffix}.fna")
            suffix += 1
        with open(out_path, "w", encoding="utf-8") as fh:
            for rec_id in found:
                header, seq = seq_records[rec_id]
                fh.write(f">{header}\n{_wrap_fasta(seq)}\n")
        return out_path, _safe_output_name(label, fallback="Cluster group", max_len=42), len(found), len(missing)

    def _collect_codon_compare_group_records(self):
        """Return non-empty comparison groups with names, paths and statistics-reference flags."""
        records = []
        used = set()
        for i, payload in enumerate(getattr(self, "codon_compare_groups", []), start=1):
            paths = _split_custom_cds_paths_gui(payload["paths"].get())
            if not paths:
                continue
            proposed_name = payload["name"].get().strip() or _suggest_group_name_from_paths(paths, fallback=f"Group {i}", max_len=42)
            name = proposed_name
            suffix = 2
            while name.lower() in used:
                name = f"{proposed_name} ({suffix})"
                suffix += 1
            used.add(name.lower())
            try:
                ref_for_statistics = bool(payload.get("ref", tk.BooleanVar(value=False)).get())
            except Exception:
                ref_for_statistics = False
            records.append({
                "name": name,
                "paths": paths,
                "ref_for_statistics": ref_for_statistics,
            })
        return records

    def _collect_codon_compare_groups(self):
        return [(r["name"], r["paths"]) for r in self._collect_codon_compare_group_records()]

    def _codon_compare_values_from_counts(self, counts_df, metric, genome_rcu_df=None):
        """Return values used by the lightweight Codon usage analyses plots.

        Fraction-style metrics are displayed as percentages in the figures:
          - RCU = % of synonymous codons for the same amino acid
          - ACU = % of all sense codons in the CDS
          - amino-acid identity = % of all amino acids in the CDS

        The raw exported workbooks are intentionally left in fraction units so
        downstream analyses remain backward-compatible. RCU z-scores are not
        percentage-scaled.
        """
        metric_l = str(metric or "").strip().lower()
        if "amino" in metric_l:
            return _aa_frequency_values_from_codon_counts(counts_df) * 100.0, "Amino acid frequency (%)", True
        if metric_l.startswith("relative"):
            rcu = _metric_values_from_codon_counts(counts_df, "Relative codon usage")
            submode = str(self.codon_compare_rcu_display.get() or "Relative codon usage with genome").strip().lower()
            if "z" in submode:
                if genome_rcu_df is None or genome_rcu_df.empty:
                    raise ValueError("Genome RCU values are required to compute RCU z-scores.")
                return _zcu_values_from_rcu(rcu, genome_rcu_df), "RCU z-score", False
            return rcu * 100.0, "Relative codon usage (%)", ("with genome" in submode)
        return _metric_values_from_codon_counts(counts_df, "Absolute codon frequency") * 100.0, "Absolute codon frequency (%)", True

    def _pvalue_to_stars(self, pval):
        try:
            pval = float(pval)
        except Exception:
            return ""
        if not np.isfinite(pval):
            return ""
        if pval < 0.001:
            return "***"
        if pval < 0.01:
            return "**"
        if pval < 0.05:
            return "*"
        return "ns"

    def _apply_custom_axis_limits(self, ax, enabled_attr, xmin_attr, xmax_attr, ymin_attr, ymax_attr, context_label="plot"):
        """Apply optional GUI-defined axis limits to an existing Matplotlib axis."""
        enabled_var = getattr(self, enabled_attr, None)
        try:
            enabled = bool(enabled_var.get()) if enabled_var is not None else False
        except Exception:
            enabled = False
        if not enabled:
            return False

        def _value(attr):
            var = getattr(self, attr, None)
            try:
                raw = var.get() if var is not None else ""
            except Exception:
                raw = ""
            return _optional_limit_value(raw)

        xmin = _value(xmin_attr)
        xmax = _value(xmax_attr)
        ymin = _value(ymin_attr)
        ymax = _value(ymax_attr)

        if xmin is not None and xmax is not None and xmax <= xmin:
            raise ValueError(f"For {context_label}, X max must be greater than X min.")
        if ymin is not None and ymax is not None and ymax <= ymin:
            raise ValueError(f"For {context_label}, Y max must be greater than Y min.")

        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        if xmin is not None or xmax is not None:
            ax.set_xlim(left=(xmin if xmin is not None else x0), right=(xmax if xmax is not None else x1))
        if ymin is not None or ymax is not None:
            ax.set_ylim(bottom=(ymin if ymin is not None else y0), top=(ymax if ymax is not None else y1))
        return True

    def _apply_codon_compare_raw_axis_limits(self, ax):
        """Apply Figure-details axis limits to Codon usage analyses comparison plots.

        The dedicated controls live in:
            Figure details > Codon usage analyses plot details >
            Custom axes for Plot raw codon usage

        As a convenience/backward-compatible fallback, if those dedicated controls
        are not enabled, the older "Per cluster codon usage profile details" axis
        controls are also honored. This makes amino-acid identity plots respond
        to the same y-axis fields users naturally try first in the Figure details tab.
        """
        applied = self._apply_custom_axis_limits(
            ax,
            "codon_compare_raw_custom_axes",
            "codon_compare_raw_xmin",
            "codon_compare_raw_xmax",
            "codon_compare_raw_ymin",
            "codon_compare_raw_ymax",
            context_label="Plot raw codon usage",
        )
        if applied:
            return True

        # Fallback to the pre-existing codon-profile controls if the user enabled
        # those rather than the dedicated Codon usage analyses controls.
        return self._apply_custom_axis_limits(
            ax,
            "codon_custom_axes",
            "codon_xmin",
            "codon_xmax",
            "codon_ymin",
            "codon_ymax",
            context_label="Plot raw codon usage",
        )

    def _codon_compare_log2_enabled(self):
        """Return True when the lightweight codon/amino-acid plots should use log2 axes."""
        try:
            return bool(self.codon_compare_log2_y_scale.get())
        except Exception:
            return False

    def _positive_axis_floor(self, values, fallback=1e-9):
        """Small positive lower bound for log axes, based on the plotted data."""
        try:
            arr = np.asarray(list(values or []), dtype=float)
            arr = arr[np.isfinite(arr) & (arr > 0)]
            if arr.size:
                return max(float(np.nanmin(arr)) * 0.5, float(fallback))
        except Exception:
            pass
        return float(fallback)

    def _apply_log2_y_axis_for_codon_compare(self, ax, positive_values=None, context_label="Plot raw usage"):
        """Apply a base-2 logarithmic Y axis with safe positive limits."""
        if not self._codon_compare_log2_enabled():
            return False
        vals = [] if positive_values is None else list(positive_values)
        positives = [float(v) for v in vals if np.isfinite(v) and float(v) > 0]
        if not positives:
            try:
                self._append_log(f"[WARN] {context_label}: log2 scale requested but no positive values were available; keeping linear Y axis.\n")
            except Exception:
                pass
            return False
        floor = self._positive_axis_floor(positives)
        y0, y1 = ax.get_ylim()
        if not np.isfinite(y1) or y1 <= 0:
            y1 = max(positives) * 1.25
        if not np.isfinite(y0) or y0 <= 0:
            y0 = floor
        if y1 <= y0:
            y1 = max(max(positives) * 1.25, y0 * 2.0)
        ax.set_ylim(bottom=y0, top=y1)
        try:
            ax.set_yscale("log", base=2)
        except TypeError:
            ax.set_yscale("log", basey=2)
        return True

    def _apply_log2_xy_axes_for_codon_compare_corr(self, ax, x_values, y_values):
        """Apply base-2 logarithmic X and Y axes to correlation plots.

        The checkbox is labeled as a Y-scale control because the raw-usage plot
        only has one quantitative axis. For correlations, both axes are switched
        to log2 so the dotted y=x reference line remains a true identity line.
        """
        if not self._codon_compare_log2_enabled():
            return False
        x_pos = [float(v) for v in np.asarray(x_values, dtype=float) if np.isfinite(v) and float(v) > 0]
        y_pos = [float(v) for v in np.asarray(y_values, dtype=float) if np.isfinite(v) and float(v) > 0]
        if not x_pos or not y_pos:
            try:
                self._append_log("[WARN] Plot correlations: log2 scale requested but one axis has no positive values; keeping linear axes for that panel.\n")
            except Exception:
                pass
            return False
        xmin_floor = self._positive_axis_floor(x_pos)
        ymin_floor = self._positive_axis_floor(y_pos)
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        if not np.isfinite(x0) or x0 <= 0:
            x0 = xmin_floor
        if not np.isfinite(y0) or y0 <= 0:
            y0 = ymin_floor
        if not np.isfinite(x1) or x1 <= x0:
            x1 = max(max(x_pos) * 1.25, x0 * 2.0)
        if not np.isfinite(y1) or y1 <= y0:
            y1 = max(max(y_pos) * 1.25, y0 * 2.0)
        ax.set_xlim(left=x0, right=x1)
        ax.set_ylim(bottom=y0, top=y1)
        try:
            ax.set_xscale("log", base=2)
            ax.set_yscale("log", base=2)
        except TypeError:
            ax.set_xscale("log", basex=2)
            ax.set_yscale("log", basey=2)
        return True

    def _compute_group_vs_genome_codons_stats(self, genome_values_df, group_values_df, codons, method):
        method_l = str(method or "None").strip().lower()
        if method_l in {"", "none", "disabled"}:
            return {}
        try:
            from scipy.stats import ttest_ind, mannwhitneyu
        except Exception as e:
            raise RuntimeError("Statistics require scipy. Please install scipy or choose Statistics = None.") from e

        out = {}
        for codon in list(codons or []):
            if codon not in genome_values_df.columns or codon not in group_values_df.columns:
                continue
            x = pd.to_numeric(genome_values_df[codon], errors="coerce").dropna().to_numpy(dtype=float)
            y = pd.to_numeric(group_values_df[codon], errors="coerce").dropna().to_numpy(dtype=float)
            if x.size < 2 or y.size < 2:
                out[codon] = np.nan
                continue
            try:
                if method_l.startswith("student"):
                    _stat, pval = ttest_ind(y, x, equal_var=True, nan_policy="omit")
                elif method_l.startswith("welch"):
                    _stat, pval = ttest_ind(y, x, equal_var=False, nan_policy="omit")
                elif "mann" in method_l:
                    _stat, pval = mannwhitneyu(y, x, alternative="two-sided")
                else:
                    pval = np.nan
            except Exception:
                pval = np.nan
            out[codon] = pval
        return out

    def _plot_codon_usage_comparison(self):
        try:
            if not self.genome_not_available.get():
                self._apply_preloaded_genome_choice()
            genome_fasta = self.fasta_path.get().strip()
            if not genome_fasta or not os.path.isfile(genome_fasta):
                raise ValueError("Please select a valid genome CDS FASTA first.")

            records = self._collect_codon_compare_group_records()
            if not records:
                raise ValueError("Please load at least one comparison FASTA/.fna file or group of files.")
            groups = [(rec["name"], rec["paths"]) for rec in records]
            for _name, paths in groups:
                missing = [p for p in paths if not os.path.isfile(p)]
                if missing:
                    raise FileNotFoundError("Missing comparison FASTA file(s):\n" + "\n".join(missing))

            metric = self.codon_compare_metric.get().strip()
            feature_order = self._codon_compare_feature_order(metric)
            if self._is_codon_compare_aa_metric(metric):
                selected_features = list(feature_order)
            else:
                selected_features, _idx = self._selected_codon_indices_for_compare()
            genome_counts = _codon_count_rows_from_fastas([genome_fasta])
            if genome_counts.empty:
                raise ValueError("No valid CDS codons were found in the selected genome FASTA.")
            genome_rcu_values = _metric_values_from_codon_counts(genome_counts, "Relative codon usage")
            genome_values, metric_label, show_genome = self._codon_compare_values_from_counts(genome_counts, metric, genome_rcu_df=genome_rcu_values)
            genome_mean, genome_sd, genome_n = _mean_sd_for_compare(genome_values, feature_order=feature_order)

            group_stats = []
            raw_exports = []
            stat_method = self.codon_compare_statistics.get().strip() or "None"
            ref_indices = [i for i, rec in enumerate(records) if rec.get("ref_for_statistics")]
            # The GUI enforces a single checked reference, but keep this robust for
            # sessions loaded from older JSON files.
            if len(ref_indices) > 1:
                ref_indices = ref_indices[:1]
            stat_reference_values = genome_values
            stat_reference_label = "genome"

            for i, (name, paths) in enumerate(groups):
                counts = _codon_count_rows_from_fastas(paths)
                if counts.empty:
                    raise ValueError(f"No valid CDS codons were found in comparison group: {name}")
                values, _metric_label2, _show_genome2 = self._codon_compare_values_from_counts(counts, metric, genome_rcu_df=genome_rcu_values)
                mean, sd, n = _mean_sd_for_compare(values, feature_order=feature_order)
                group_stats.append(dict(name=name, mean=mean, sd=sd, n=n, values_df=values, pvals={}, is_reference=False))
                raw_exports.append((name, counts))

            if ref_indices:
                ref_idx = ref_indices[0]
                if 0 <= ref_idx < len(group_stats):
                    stat_reference_values = group_stats[ref_idx]["values_df"]
                    stat_reference_label = group_stats[ref_idx]["name"]
                    group_stats[ref_idx]["is_reference"] = True

            for stat in group_stats:
                if stat.get("is_reference"):
                    stat["pvals"] = {}
                else:
                    stat["pvals"] = self._compute_group_vs_genome_codons_stats(
                        stat_reference_values,
                        stat.get("values_df"),
                        selected_features,
                        stat_method,
                    )

            self._export_codon_usage_comparison_raw_tables(raw_exports, genome_rcu_values)
            self._render_codon_usage_comparison_plot(metric_label, genome_mean, genome_sd, genome_n, group_stats, selected_codons=selected_features, show_genome=show_genome, stat_method=stat_method, stat_reference_label=stat_reference_label)
        except Exception as e:
            messagebox.showerror("Codon usage comparison", str(e))

    def _codon_usage_analysis_output_dir(self):
        """Return the project output root used by lightweight codon-usage utilities."""
        candidates = [
            getattr(self, "last_clustering_output_dir", ""),
            self.default_root.get().strip() if hasattr(self, "default_root") else "",
            os.getcwd(),
        ]
        for d in candidates:
            d = str(d or "").strip()
            if not d:
                continue
            try:
                os.makedirs(d, exist_ok=True)
                return d
            except Exception:
                continue
        return os.getcwd()

    def _export_codon_usage_comparison_raw_tables(self, group_counts, genome_rcu_values):
        """Export one raw codon-usage workbook per FASTA comparison group."""
        if not group_counts:
            return []
        out_root = self._codon_usage_analysis_output_dir()
        raw_dir = os.path.join(out_root, "Raw codon usage tables")
        os.makedirs(raw_dir, exist_ok=True)
        written = []
        used = set()
        for name, counts in group_counts:
            base = _safe_output_name(f"Codon usage comparison - {name}", fallback="Codon usage comparison", max_len=90)
            fname = f"{base}.xlsx"
            suffix = 2
            while fname.lower() in used:
                fname = f"{base} ({suffix}).xlsx"
                suffix += 1
            used.add(fname.lower())
            out_path = os.path.join(raw_dir, fname)
            try:
                _write_codon_usage_raw_tables_xlsx(out_path, counts, genome_rcu_df=genome_rcu_values)
                written.append(out_path)
                self._append_log(f"[INFO] Saved raw codon-usage tables for '{name}': {out_path}\n")
            except Exception as e:
                self._append_log(f"[WARN] Could not save raw codon-usage tables for '{name}': {e}\n")
        return written

    def _codon_compare_style_value(self, attr, default, min_value=None, max_value=None, cast=float):
        """Read a numeric Codon usage analyses style control safely."""
        try:
            var = getattr(self, attr, None)
            raw = var.get() if var is not None else default
            val = cast(float(str(raw).strip())) if cast is int else cast(str(raw).strip())
        except Exception:
            val = default
        try:
            if min_value is not None:
                val = max(min_value, val)
            if max_value is not None:
                val = min(max_value, val)
        except Exception:
            val = default
        return val

    def _codon_compare_style(self):
        """Return shared plotting style for the lightweight codon/amino-acid plots."""
        caption = float(self._codon_compare_style_value("codon_compare_caption_size", 14.0, min_value=6.0, max_value=40.0, cast=float))
        marker = float(self._codon_compare_style_value("codon_compare_marker_size", 8.0, min_value=1.0, max_value=40.0, cast=float))
        line = float(self._codon_compare_style_value("codon_compare_line_width", 2.0, min_value=0.2, max_value=10.0, cast=float))
        aa_spacing = float(self._codon_compare_style_value("codon_compare_aa_spacing", 1.75, min_value=1.0, max_value=4.0, cast=float))
        codon_aa_gap = float(self._codon_compare_style_value("codon_compare_codon_aa_gap", 0.65, min_value=0.0, max_value=3.0, cast=float))
        codon_gap = float(self._codon_compare_style_value("codon_compare_codon_gap", 0.33, min_value=0.0, max_value=2.0, cast=float))
        legend_ncol = int(self._codon_compare_style_value("codon_compare_legend_ncol", 3, min_value=1, max_value=8, cast=int))
        try:
            highlight_color = str(self.codon_compare_highlight_color.get()).strip() or "#0057D9"
        except Exception:
            highlight_color = "#0057D9"
        return dict(
            caption_size=caption,
            tick_size=max(6.0, caption),
            label_size=max(7.0, caption),
            axis_title_size=max(7.0, caption * 1.40),
            title_size=max(8.0, caption + 2.0),
            legend_size=max(6.0, caption * 1.40),
            aa_group_label_size=max(8.5, caption * 1.40),
            marker_size=marker,
            line_width=line,
            aa_spacing=aa_spacing,
            codon_aa_gap=codon_aa_gap,
            codon_gap=codon_gap,
            legend_ncol=legend_ncol,
            highlight_color=highlight_color,
        )

    def _place_codon_compare_legend(self, ax, handles=None, labels=None, style=None):
        """Place legends above Codon usage analyses plots, wrapping after N entries."""
        style = style or self._codon_compare_style()
        if handles is None or labels is None:
            handles, labels = ax.get_legend_handles_labels()
        if not handles:
            return None
        ncol = max(1, min(int(style.get("legend_ncol", 3)), len(handles)))
        return ax.legend(
            handles=handles,
            labels=labels,
            loc="lower center",
            bbox_to_anchor=(0.5, 1.045),
            ncol=ncol,
            frameon=False,
            fontsize=style.get("legend_size", 14),
            borderaxespad=0.0,
            handlelength=1.8,
            columnspacing=1.2,
            handletextpad=0.5,
        )

    def _codon_compare_offsets(self, n_series, is_aa_metric=False, feature_spacing=1.0):
        """Series offsets within a codon/amino-acid group, scaled for AA plots."""
        n_series = max(1, int(n_series))
        if is_aa_metric:
            span = min(1.05, max(0.55, float(feature_spacing) * 0.58))
            return np.linspace(-0.5 * span, 0.5 * span, n_series) if n_series > 1 else np.array([0.0])
        return np.linspace(-0.36, 0.36, n_series) if n_series > 1 else np.array([0.0])

    def _codon_compare_feature_positions(self, selected_features, is_aa_metric=False, feature_spacing=1.0, codon_aa_gap=0.65, codon_gap=0.33):
        """Return plotted x positions and amino-acid group spans for raw-usage plots.

        Amino-acid identity plots use a uniform spacing controlled by AA spacing.
        Codon-usage plots keep synonymous codons grouped, insert a small
        gap between individual codons, and add a larger extra gap when the
        amino-acid family changes. This makes synonymous codon boxes easier
        to read while preserving clear separation between neighboring amino acids.
        """
        features = list(selected_features or [])
        if not features:
            return np.asarray([], dtype=float), {}

        if is_aa_metric:
            x = np.arange(len(features), dtype=float) * float(feature_spacing)
            return x, {str(aa): [float(pos)] for aa, pos in zip(features, x)}

        positions = []
        group_positions = {}
        xpos = 0.0
        prev_aa = None
        for feature in features:
            codon = str(feature).upper().replace("U", "T")
            aa = _CODON_TO_AA3.get(codon, "")
            if prev_aa is not None and aa != prev_aa:
                xpos += float(codon_aa_gap)
            positions.append(float(xpos))
            if aa:
                group_positions.setdefault(aa, []).append(float(xpos))
            xpos += 1.0 + float(codon_gap)
            prev_aa = aa
        return np.asarray(positions, dtype=float), group_positions

    def _render_codon_usage_comparison_plot(self, metric_label, genome_mean, genome_sd, genome_n, group_stats, selected_codons=None, show_genome=True, stat_method="None", stat_reference_label="genome"):
        # Robustly detect amino-acid plots from either the GUI label ("Amino acid identity")
        # or the computed metric label ("Amino acid frequency" / "Amino-acid frequency").
        is_aa_metric = self._is_codon_compare_aa_metric(metric_label)
        feature_order = list(_AA_ORDER) if is_aa_metric else list(_CODON_ORDER)
        selected_codons = [c for c in list(selected_codons or feature_order) if c in feature_order]
        if not selected_codons:
            selected_codons = list(feature_order)
        idx = np.asarray([feature_order.index(c) for c in selected_codons], dtype=int)
        # In amino-acid identity mode, give each amino acid a little more
        # horizontal breathing room. This keeps several clusters visually grouped
        # around the same amino acid while leaving a small gap before the next
        # amino-acid group. Codon plots use a smaller configurable within-codon gap.
        style = self._codon_compare_style()
        feature_spacing = float(style["aa_spacing"]) if is_aa_metric else 1.0
        x, aa_group_positions = self._codon_compare_feature_positions(
            selected_codons,
            is_aa_metric=is_aa_metric,
            feature_spacing=feature_spacing,
            codon_aa_gap=float(style.get("codon_aa_gap", 0.65)),
            codon_gap=float(style.get("codon_gap", 0.33)),
        )
        n_groups = len(group_stats)
        plot_style_var = getattr(self, "codon_compare_plot_style", None)
        try:
            plot_style = str(plot_style_var.get() if plot_style_var is not None else "Mean ± SD").strip().lower()
        except Exception:
            plot_style = "mean ± sd"
        is_line = "line" in plot_style
        is_box = "box" in plot_style
        is_violin = "violin" in plot_style
        is_distribution = is_box or is_violin

        width_per_feature = 0.58 if is_aa_metric else 0.36
        x_extent = float(np.nanmax(x) - np.nanmin(x) + 1.0) if len(x) else float(len(selected_codons))
        fig_w = max(10.0, min(34.0, width_per_feature * x_extent + 3.2 + 0.20 * max(0, n_groups - 1)))
        fig_h = 8.4
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))

        gm = np.asarray(genome_mean, dtype=float)[idx]
        gs = np.asarray(genome_sd, dtype=float)[idx]
        yerr_genome = np.where(np.isfinite(gs), gs, 0.0)

        colors = ["red", "#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e", "#8c564b", "#17becf", "#bcbd22"]
        markers = ["o", "^", "D", "v", "P", "X", "h", "*"]
        plotted_values = []

        def _finite_array(values):
            arr = np.asarray(values, dtype=float)
            return arr[np.isfinite(arr)]

        def _selected_values_by_codon(values_df):
            vals = []
            for codon in selected_codons:
                if values_df is None or codon not in values_df.columns:
                    vals.append(np.asarray([], dtype=float))
                else:
                    vals.append(_finite_array(pd.to_numeric(values_df[codon], errors="coerce").to_numpy(dtype=float)))
            return vals

        # Keep a per-group/per-codon y position for significance annotations.
        stat_y = {}

        if is_distribution:
            # Genome + each selected group are drawn side-by-side for every codon.
            series = []
            if show_genome:
                # Reconstruct the genome distribution from the already computed genome values when available.
                # The full genome values are not stored in group_stats, so the genome distribution is shown as
                # mean±SD markers; cluster groups are shown as true per-gene box/violin distributions.
                series.append({"name": f"Genome mean ± SD (n={genome_n})", "kind": "genome", "color": "white"})
            for i, stat in enumerate(group_stats):
                series.append({"name": f"{stat['name']} (n={stat['n']})", "kind": "group", "stat": stat, "color": colors[i % len(colors)]})
            n_series = max(1, len(series))
            offsets = self._codon_compare_offsets(n_series, is_aa_metric=is_aa_metric, feature_spacing=feature_spacing)
            span = float(np.nanmax(offsets) - np.nanmin(offsets)) if len(offsets) > 1 else (0.72 if not is_aa_metric else min(1.05, feature_spacing * 0.58))
            width = min(0.24 if is_aa_metric else 0.16, max(0.04, span * 0.72 / max(n_series, 1)))

            legend_handles = []
            try:
                from matplotlib.lines import Line2D
                from matplotlib.patches import Patch
            except Exception:
                Line2D = None
                Patch = None

            for si, ser in enumerate(series):
                off = float(offsets[si])
                if ser["kind"] == "genome":
                    ax.errorbar(
                        x + off, gm, yerr=yerr_genome,
                        fmt="s", ms=style["marker_size"], mfc="white", mec="black", mew=max(0.8, style["line_width"] * 0.55),
                        ecolor="black", elinewidth=style["line_width"], capsize=3.5, linestyle="none",
                        zorder=5, label=ser["name"],
                    )
                    plotted_values.extend((gm + yerr_genome)[np.isfinite(gm + yerr_genome)].tolist())
                    plotted_values.extend((gm - yerr_genome)[np.isfinite(gm - yerr_genome)].tolist())
                    if Line2D is not None:
                        legend_handles.append(Line2D([0], [0], marker="s", color="black", markerfacecolor="white", linestyle="none", markersize=style["marker_size"], label=ser["name"]))
                    continue

                stat = ser["stat"]
                data = _selected_values_by_codon(stat.get("values_df"))
                safe_data = [arr if arr.size else np.asarray([np.nan]) for arr in data]
                positions = x + off
                if is_box:
                    bp = ax.boxplot(
                        safe_data, positions=positions, widths=width, patch_artist=True,
                        manage_ticks=False, showfliers=False, medianprops={"color": "black", "linewidth": style["line_width"]},
                        whiskerprops={"linewidth": max(0.8, style["line_width"] * 0.8)}, capprops={"linewidth": max(0.8, style["line_width"] * 0.8)},
                    )
                    for patch in bp.get("boxes", []):
                        patch.set_facecolor(ser["color"])
                        patch.set_alpha(0.45)
                        patch.set_edgecolor("black")
                        patch.set_linewidth(max(0.8, style["line_width"] * 0.8))
                else:
                    for arr, pos in zip(safe_data, positions):
                        if not np.isfinite(arr).any():
                            continue
                        vp = ax.violinplot([arr], positions=[pos], widths=width * 1.35, showmeans=False, showmedians=True, showextrema=False)
                        for body in vp.get("bodies", []):
                            body.set_facecolor(ser["color"])
                            body.set_edgecolor("black")
                            body.set_alpha(0.45)
                        for part_name in ("cmedians",):
                            part = vp.get(part_name)
                            if part is not None:
                                part.set_color("black")
                                part.set_linewidth(style["line_width"])
                for j, arr in enumerate(data):
                    if arr.size:
                        plotted_values.extend(arr.tolist())
                        stat_y[(stat["name"], selected_codons[j])] = float(np.nanpercentile(arr, 95)) if np.isfinite(arr).any() else np.nan
                if Patch is not None:
                    legend_handles.append(Patch(facecolor=ser["color"], edgecolor="black", alpha=0.45, label=ser["name"]))
            if legend_handles:
                self._place_codon_compare_legend(ax, handles=legend_handles, labels=[h.get_label() for h in legend_handles], style=style)

        elif is_line:
            if show_genome:
                ax.plot(x, gm, marker="s", ms=style["marker_size"], mfc="white", mec="black", color="black", linewidth=style["line_width"],
                        label=f"Genome mean (n={genome_n})", zorder=4)
                plotted_values.extend(gm[np.isfinite(gm)].tolist())
            for i, stat in enumerate(group_stats):
                mean = np.asarray(stat["mean"], dtype=float)[idx]
                color = colors[i % len(colors)]
                ax.plot(x, mean, marker=markers[i % len(markers)], ms=style["marker_size"], color=color, linewidth=style["line_width"],
                        markeredgecolor="black", markeredgewidth=max(0.4, style["line_width"] * 0.25), label=f"{stat['name']} (n={stat['n']})", zorder=5 + i)
                plotted_values.extend(mean[np.isfinite(mean)].tolist())
                for j, codon in enumerate(selected_codons):
                    stat_y[(stat["name"], codon)] = mean[j]
            self._place_codon_compare_legend(ax, style=style)

        else:
            if show_genome:
                all_offsets = self._codon_compare_offsets(n_groups + 1, is_aa_metric=is_aa_metric, feature_spacing=feature_spacing)
                genome_offset = float(all_offsets[0]) if len(all_offsets) else (-0.18)
                offsets = all_offsets[1:] if len(all_offsets) > 1 else np.array([0.0])
            else:
                genome_offset = 0.0
                offsets = self._codon_compare_offsets(n_groups, is_aa_metric=is_aa_metric, feature_spacing=feature_spacing)
            if show_genome:
                ax.errorbar(
                    x + genome_offset, gm, yerr=yerr_genome,
                    fmt="s", ms=style["marker_size"], mfc="white", mec="black", mew=max(0.8, style["line_width"] * 0.55),
                    ecolor="black", elinewidth=style["line_width"], capsize=3.5, linestyle="none",
                    label=f"Genome mean ± SD (n={genome_n})", zorder=4,
                )
                plotted_values.extend((gm + yerr_genome)[np.isfinite(gm + yerr_genome)].tolist())
                plotted_values.extend((gm - yerr_genome)[np.isfinite(gm - yerr_genome)].tolist())
            for i, stat in enumerate(group_stats):
                mean = np.asarray(stat["mean"], dtype=float)[idx]
                sd = np.asarray(stat["sd"], dtype=float)[idx]
                n = stat["n"]
                yerr = np.where(np.isfinite(sd), sd, 0.0) if n > 1 else None
                ax.errorbar(
                    x + float(offsets[i]), mean, yerr=yerr,
                    fmt=markers[i % len(markers)], ms=style["marker_size"],
                    mfc=colors[i % len(colors)], mec="black", mew=max(0.4, style["line_width"] * 0.25),
                    ecolor=colors[i % len(colors)], elinewidth=style["line_width"], capsize=3.5,
                    linestyle="none", label=f"{stat['name']} (n={n})", zorder=5 + i,
                )
                err = np.where(np.isfinite(sd), sd, 0.0)
                plotted_values.extend((mean + err)[np.isfinite(mean + err)].tolist())
                plotted_values.extend((mean - err)[np.isfinite(mean - err)].tolist())
                for j, codon in enumerate(selected_codons):
                    stat_y[(stat["name"], codon)] = mean[j] + err[j] if np.isfinite(mean[j]) else np.nan
            self._place_codon_compare_legend(ax, style=style)

        axis_pad = max(0.75, 0.55 * feature_spacing) if is_aa_metric else max(0.75, 0.50 + 0.50 * (float(style.get("codon_aa_gap", 0.65)) + float(style.get("codon_gap", 0.33))))
        if len(x) > 0:
            ax.set_xlim(float(x[0]) - axis_pad, float(x[-1]) + axis_pad)
        else:
            ax.set_xlim(-axis_pad, axis_pad)
        metric_lower = str(metric_label).lower()
        finite_plotted = [v for v in plotted_values if np.isfinite(v)]
        if "z-score" in metric_lower:
            ax.set_ylabel("RCU z-score vs genome", fontsize=style["axis_title_size"])
            lo = min(finite_plotted) if finite_plotted else -1.0
            hi = max(finite_plotted) if finite_plotted else 1.0
            pad = max(0.25, 0.15 * (hi - lo if hi > lo else 1.0))
            ax.axhline(0, color="0.5", lw=style["line_width"], ls=(0, (3, 3)), zorder=1)
            ax.set_ylim(lo - pad, hi + pad)
        elif "relative" in metric_lower:
            ax.set_ylabel("Relative codon usage (%)", fontsize=style["axis_title_size"])
            ymax = 105.0
            if finite_plotted:
                ymax = max(ymax, float(np.nanmax(finite_plotted)) * 1.12)
            ax.set_ylim(bottom=0, top=ymax)
        elif "amino acid" in metric_lower:
            ax.set_ylabel("Amino-acid frequency (%)", fontsize=style["axis_title_size"])
            upper = max(finite_plotted) if finite_plotted else 0.0
            ax.set_ylim(bottom=0, top=max(upper * 1.22, 5.0))
        else:
            ax.set_ylabel("Absolute codon frequency (%)", fontsize=style["axis_title_size"])
            upper = max(finite_plotted) if finite_plotted else 0.0
            ax.set_ylim(bottom=0, top=max(upper * 1.22, 2.0))

        stat_method_l = str(stat_method or "None").strip().lower()
        if stat_method_l not in {"", "none", "disabled"}:
            y0, y1 = ax.get_ylim()
            pad = 0.025 * (y1 - y0 if y1 > y0 else 1.0)
            if is_distribution:
                n_series = n_groups + (1 if show_genome else 0)
                series_offsets = self._codon_compare_offsets(n_series, is_aa_metric=is_aa_metric, feature_spacing=feature_spacing)
                group_offsets = series_offsets[(1 if show_genome else 0):]
            elif is_line:
                group_offsets = np.zeros(max(n_groups, 1), dtype=float)
            else:
                group_offsets = offsets if 'offsets' in locals() else self._codon_compare_offsets(n_groups, is_aa_metric=is_aa_metric, feature_spacing=feature_spacing)
            for i, stat in enumerate(group_stats):
                for j, codon in enumerate(selected_codons):
                    pval = (stat.get("pvals") or {}).get(codon, np.nan)
                    stars = self._pvalue_to_stars(pval)
                    y_base = stat_y.get((stat["name"], codon), np.nan)
                    if not stars or not np.isfinite(y_base):
                        continue
                    y_text = y_base + pad * (1.0 + 0.65 * i)
                    off = float(group_offsets[i]) if i < len(group_offsets) else 0.0
                    ax.text(
                        x[j] + off, y_text, stars,
                        ha="center", va="bottom", fontsize=max(8.0, style["caption_size"] * 0.85), rotation=90 if stars == "ns" else 0,
                        color="0.35" if stars == "ns" else "black", clip_on=False, zorder=20,
                    )
            y0, y1 = ax.get_ylim()
            ax.set_ylim(y0, y1 + 0.10 * (y1 - y0 if y1 > y0 else 1.0))

        # Apply Figure-details custom axes after automatic scaling and after
        # optional statistical annotations, so user y-limits also work for
        # amino-acid identity plots.
        self._apply_codon_compare_raw_axis_limits(ax)

        title_org = self.organism_name.get().strip() or "selected genome"
        style_var = getattr(self, "codon_compare_plot_style", None)
        try:
            style_label = str(style_var.get() if style_var is not None else "Mean ± SD") or "Mean ± SD"
        except Exception:
            style_label = "Mean ± SD"
        stat_method_l_for_title = str(stat_method or "None").strip().lower()
        ref_txt = ""
        if stat_method_l_for_title not in {"", "none", "disabled"}:
            ref_txt = f" | stats ref: {stat_reference_label}"
        # Keep the raw-usage plot clean for dense codon labels: no title and no
        # x-axis title. The legend and output filename still document the plotted metric.
        ax.set_title("")
        ax.set_xlabel("")
        ax.set_xticks(x)
        ax.set_xticklabels(selected_codons, rotation=90, fontsize=style["tick_size"])
        if not is_aa_metric:
            highlighted_codons = {str(c).upper().replace("U", "T") for c in getattr(self, "codon_compare_highlighted_codons", [])}
            if highlighted_codons:
                for tick_label in ax.get_xticklabels():
                    txt = str(tick_label.get_text()).upper().replace("U", "T")
                    if txt in highlighted_codons:
                        tick_label.set_fontweight("bold")
                        tick_label.set_color(style.get("highlight_color", "#0057D9"))
        ax.tick_params(axis="x", length=0, pad=3, labelsize=style["tick_size"])
        ax.tick_params(axis="y", labelsize=style["tick_size"])
        ax.grid(axis="y", alpha=0.25, linewidth=max(0.7, style["line_width"] * 0.4))
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        if not is_aa_metric:
            # Synonymous-family bars and amino-acid labels below codon tick labels.
            trans = ax.get_xaxis_transform()
            # aa_group_positions uses the true x coordinates, including the
            # extra gap inserted between amino-acid families.
            for aa, present in aa_group_positions.items():
                present = [float(v) for v in present if np.isfinite(v)]
                if not present:
                    continue
                start = min(present)
                end = max(present)
                center = (start + end) / 2
                ax.plot([start - 0.38, end + 0.38], [-0.18, -0.18], transform=trans, color="black", lw=style["line_width"], clip_on=False)
                ax.text(center, -0.27, aa, transform=trans, ha="center", va="top", fontsize=style["aa_group_label_size"], clip_on=False)

        if not ax.get_legend_handles_labels()[0]:
            self._place_codon_compare_legend(ax, style=style)
        legend_rows = int(np.ceil((len(ax.get_legend_handles_labels()[0]) or 1) / max(1, int(style["legend_ncol"]))))
        # With titles removed, keep the legend close to the plotting area instead
        # of preserving the older title-sized vertical gap.
        fig.subplots_adjust(left=0.09, right=0.98, bottom=(0.20 if is_aa_metric else 0.31),
                            top=max(0.64, 0.89 - 0.050 * max(0, legend_rows - 1)))

        # Re-apply once at the end to protect the user-defined limits from any
        # future plot-formatting code added above this point.
        self._apply_codon_compare_raw_axis_limits(ax)
        self._apply_log2_y_axis_for_codon_compare(ax, positive_values=finite_plotted, context_label="Plot raw usage")

        try:
            out_root = self._codon_usage_analysis_output_dir()
            out_dir = os.path.join(out_root, "Figures")
            os.makedirs(out_dir, exist_ok=True)
            safe_metric = re.sub(r"[^A-Za-z0-9]+", "_", str(metric_label)).strip("_")
            safe_style = re.sub(r"[^A-Za-z0-9]+", "_", str(style_label)).strip("_")
            base_name = "Amino acid identity comparison" if is_aa_metric else "Codon usage comparison"
            out_path = os.path.join(out_dir, f"{base_name} - {safe_metric} - {safe_style}.png")
            fig.savefig(out_path, dpi=300, bbox_inches="tight")
            self._append_log(f"[INFO] Saved codon usage comparison plot: {out_path}\n")
        except Exception as e:
            self._append_log(f"[WARN] Could not save codon usage comparison plot: {e}\n")

        _show_matplotlib_figure_nonblocking(fig)

    def _plot_codon_usage_correlations(self):
        """Plot codon-usage correlations between selected reference groups and other groups.

        If one group is ticked as "Ref for statistics", correlation plots keep that group on the X axis.
        If no reference is ticked, all pairwise comparisons are plotted. Each point represents
        one sense codon. A dotted y=x reference line and the Pearson correlation R²
        are shown in each panel.
        """
        try:
            records = self._collect_codon_compare_group_records()
            if len(records) < 2:
                raise ValueError(
                    "Please load at least two FASTA/cluster comparison groups before plotting correlations."
                )
            for rec in records:
                missing = [p for p in rec["paths"] if not os.path.isfile(p)]
                if missing:
                    raise FileNotFoundError("Missing comparison FASTA file(s):\n" + "\n".join(missing))

            metric = self.codon_compare_metric.get().strip()
            feature_order = self._codon_compare_feature_order(metric)
            genome_rcu_values = None
            if str(metric or "").strip().lower().startswith("relative"):
                submode = str(self.codon_compare_rcu_display.get() or "").strip().lower()
                if "z" in submode:
                    if not self.genome_not_available.get():
                        self._apply_preloaded_genome_choice()
                    genome_fasta = self.fasta_path.get().strip()
                    if not genome_fasta or not os.path.isfile(genome_fasta):
                        raise ValueError("Relative codon usage z-score correlations require a valid genome CDS FASTA baseline.")
                    genome_counts = _codon_count_rows_from_fastas([genome_fasta])
                    genome_rcu_values = _metric_values_from_codon_counts(genome_counts, "Relative codon usage")
            group_means = []
            group_ns = []
            for rec in records:
                name = rec["name"]
                paths = rec["paths"]
                counts = _codon_count_rows_from_fastas(paths)
                if counts.empty:
                    raise ValueError(f"No valid CDS codons were found in comparison group: {name}")
                values, _metric_label, _show_genome = self._codon_compare_values_from_counts(counts, metric, genome_rcu_df=genome_rcu_values)
                mean, _sd, n = _mean_sd_for_compare(values, feature_order=feature_order)
                group_means.append((name, mean))
                group_ns.append(n)

            ref_indices = [i for i, rec in enumerate(records) if rec.get("ref_for_statistics")]
            if ref_indices:
                pair_indices = []
                seen_unordered = set()
                for ref_idx in ref_indices:
                    for other_idx in range(len(records)):
                        if other_idx == ref_idx:
                            continue
                        key = tuple(sorted((ref_idx, other_idx)))
                        if key in seen_unordered:
                            continue
                        seen_unordered.add(key)
                        # The checked reference group is kept on the X axis.
                        pair_indices.append((ref_idx, other_idx))
                if not pair_indices:
                    raise ValueError("Please select at least two non-empty comparison groups.")
            else:
                # Backward-compatible fallback: no reference selected means all pairwise comparisons.
                pair_indices = list(itertools.combinations(range(len(records)), 2))

            if self._is_codon_compare_aa_metric(metric):
                selected_features = list(feature_order)
            else:
                selected_features, _idx = self._selected_codon_indices_for_compare()
            self._render_codon_usage_correlation_plot(metric, group_means, group_ns, pair_indices=pair_indices, selected_codons=selected_features)
        except Exception as e:
            messagebox.showerror("Codon usage correlations", str(e))

    def _render_codon_usage_correlation_plot(self, metric, group_means, group_ns, pair_indices=None, selected_codons=None):
        if len(group_means) < 2:
            raise ValueError("At least two comparison groups are required.")

        if pair_indices is None:
            pair_indices = list(itertools.combinations(range(len(group_means)), 2))
        pair_indices = list(pair_indices)
        if not pair_indices:
            raise ValueError("No codon-usage correlation comparisons were selected.")
        n_panels = len(pair_indices)
        n_cols = min(3, n_panels)
        n_rows = int(np.ceil(n_panels / n_cols))
        style = self._codon_compare_style()
        fig_w = max(6.8 * n_cols, 7.2)
        fig_h = max(6.4 * n_rows, 6.2)
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_w, fig_h), squeeze=False)
        axes_flat = axes.ravel()

        is_aa_metric = self._is_codon_compare_aa_metric(metric)
        feature_order = list(_AA_ORDER) if is_aa_metric else list(_CODON_ORDER)
        selected_codons = [c for c in list(selected_codons or feature_order) if c in feature_order]
        if not selected_codons:
            selected_codons = list(feature_order)
        selected_idx = np.asarray([feature_order.index(c) for c in selected_codons], dtype=int)

        # Use a stable, pleasant amino-acid coloring so synonymous codons are easy to inspect.
        cmap = plt.get_cmap("tab20")
        aa_names = list(_AA_TO_CODONS.keys())
        aa_color = {aa: cmap(i % cmap.N) for i, aa in enumerate(aa_names)}
        if is_aa_metric:
            codon_color = [aa_color.get(aa, cmap(i % cmap.N)) for i, aa in enumerate(selected_codons)]
        else:
            codon_color = [aa_color[_CODON_TO_AA3[c]] for c in selected_codons]

        for panel_idx, (idx_x, idx_y) in enumerate(pair_indices):
            ax = axes_flat[panel_idx]
            x_name, x_mean = group_means[idx_x]
            y_name, y_mean = group_means[idx_y]
            x_n = group_ns[idx_x]
            y_n = group_ns[idx_y]
            x_all = np.asarray(x_mean, dtype=float)[selected_idx]
            y_all = np.asarray(y_mean, dtype=float)[selected_idx]
            mask = np.isfinite(x_all) & np.isfinite(y_all)
            x = x_all[mask]
            y = y_all[mask]
            codons = np.asarray(selected_codons, dtype=object)[mask]
            colors = np.asarray(codon_color, dtype=object)[mask]
            if x.size < 2:
                ax.text(0.5, 0.5, "Not enough valid codons", ha="center", va="center", transform=ax.transAxes)
                ax.set_axis_off()
                continue

            ax.scatter(x, y, s=float(style["marker_size"]) ** 2, c=list(colors), edgecolors="black", linewidths=max(0.35, style["line_width"] * 0.25), alpha=0.92, zorder=4)

            # Codon labels make the scatter directly interpretable.
            for xi, yi, codon in zip(x, y, codons):
                ax.annotate(str(codon), (xi, yi), xytext=(4, 3), textcoords="offset points", fontsize=max(7.0, style["caption_size"] * 0.65), alpha=0.85)

            # The displayed R² remains a correlation statistic: Pearson r squared
            # across the plotted codons. The dotted line is a fixed y=x reference,
            # not a fitted trendline.
            if x.size >= 2 and np.nanstd(x) > 0 and np.nanstd(y) > 0:
                r = float(np.corrcoef(x, y)[0, 1])
                r2 = r ** 2 if np.isfinite(r) else np.nan
            else:
                r2 = np.nan

            lim_max = float(np.nanmax([np.nanmax(x), np.nanmax(y), 1e-12])) * 1.10
            metric_l = str(metric).lower()
            if metric_l.startswith("relative"):
                lim_max = max(lim_max, 102.0)
            elif "amino" in metric_l:
                lim_max = max(lim_max, 5.0)
            elif "z" in metric_l:
                lim_max = max(lim_max, 1.0)
            else:
                lim_max = max(lim_max, 2.0)
            ax.set_xlim(0, lim_max)
            ax.set_ylim(0, lim_max)
            self._apply_custom_axis_limits(
                ax,
                "codon_compare_corr_custom_axes",
                "codon_compare_corr_xmin",
                "codon_compare_corr_xmax",
                "codon_compare_corr_ymin",
                "codon_compare_corr_ymax",
                context_label="Plot correlations",
            )
            self._apply_log2_xy_axes_for_codon_compare_corr(ax, x, y)

            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            line_min = min(float(xlim[0]), float(ylim[0]))
            line_max = max(float(xlim[1]), float(ylim[1]))
            line_x = np.linspace(line_min, line_max, 100)
            ax.plot(line_x, line_x, color="black", linestyle=(0, (2, 2)), linewidth=style["line_width"], zorder=3)
            ax.set_aspect("equal", adjustable="box")
            ax.grid(True, alpha=0.25, linewidth=max(0.7, style["line_width"] * 0.4))
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            if is_aa_metric:
                axis_metric = "amino-acid frequency (%)"
            elif str(metric).strip().lower().startswith("relative") and "z" not in str(metric).strip().lower():
                axis_metric = "relative codon usage (%)"
            elif "z" in str(metric).strip().lower():
                axis_metric = "RCU z-score"
            else:
                axis_metric = "absolute codon frequency (%)"
            ax.set_xlabel(f"{x_name} mean {axis_metric} (n={x_n})", fontsize=style["axis_title_size"])
            ax.set_ylabel(f"{y_name} mean {axis_metric} (n={y_n})", fontsize=style["axis_title_size"])
            # Plot titles are intentionally suppressed in Codon usage analysis figures;
            # the legend and axis labels already identify the comparison.
            ax.set_title("")
            ax.tick_params(axis="both", labelsize=style["tick_size"])

            r2_txt = "NA" if not np.isfinite(r2) else f"{r2:.3f}"
            ax.text(
                0.04, 0.96,
                f"y = x reference\nR² = {r2_txt}",
                ha="left", va="top", transform=ax.transAxes,
                fontsize=max(8.0, style["caption_size"] * 0.85),
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="0.75", alpha=0.9),
            )

        for j in range(n_panels, len(axes_flat)):
            axes_flat[j].set_axis_off()

        metric_label = str(metric).strip() or "Codon usage"
        # No figure title: keep the space between the graph and the legend clean.

        # Compact amino-acid legend using colored dots. It follows the same
        # "entries per row" control as the raw comparison plot.
        handles = []
        for aa in aa_names:
            handles.append(plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=aa_color[aa],
                                      markeredgecolor="black", markeredgewidth=max(0.35, style["line_width"] * 0.20), markersize=style["marker_size"], label=aa))
        legend_ncol = max(1, min(int(style["legend_ncol"]), len(handles)))
        fig.legend(handles=handles, title=("Amino acid" if not is_aa_metric else "Amino acids"), loc="upper center", bbox_to_anchor=(0.5, 0.945),
                   frameon=False, fontsize=style["legend_size"], title_fontsize=style["legend_size"], ncol=legend_ncol,
                   columnspacing=1.1, handletextpad=0.4)
        legend_rows = int(np.ceil(len(handles) / max(1, legend_ncol)))
        # Titles are removed, so the panels can move upward and sit closer to
        # the legend without creating the old title-sized gap.
        fig.subplots_adjust(left=0.08, right=0.98, bottom=0.09,
                            top=max(0.58, 0.90 - 0.040 * max(0, legend_rows - 1)),
                            wspace=0.34, hspace=0.34)

        try:
            out_root = self._codon_usage_analysis_output_dir()
            out_dir = os.path.join(out_root, "Figures")
            os.makedirs(out_dir, exist_ok=True)
            safe_metric = re.sub(r"[^A-Za-z0-9]+", "_", str(metric)).strip("_")
            base_name = "Amino acid identity correlations" if is_aa_metric else "Codon usage correlations"
            out_path = os.path.join(out_dir, f"{base_name} - {safe_metric}.png")
            fig.savefig(out_path, dpi=300, bbox_inches="tight")
            self._append_log(f"[INFO] Saved codon usage correlation plot: {out_path}\n")
        except Exception as e:
            self._append_log(f"[WARN] Could not save codon usage correlation plot: {e}\n")

        _show_matplotlib_figure_nonblocking(fig)

    def _close_interface(self):
        try:
            self.root.after(0, self.root.destroy)
        except Exception:
            try:
                self.root.destroy()
            except Exception:
                pass

    def _clear_log(self):
        self.log_text.delete("1.0", "end")

    def _append_log(self, text):
        self.log_text.insert("end", text)
        self.log_text.see("end")
        self.root.update_idletasks()

    def _poll_log_queue(self):
        try:
            while True:
                msg = self.log_queue.get_nowait()
                self._append_log(msg)
        except queue.Empty:
            pass
        self.root.after(100, self._poll_log_queue)

    def _apply_default_settings(self):
        self.default_root.set(DEFAULT_BROWSER_ROOT)
        if self.preloaded_genomes:
            self.genome_not_available.set(False)
            self.preloaded_genome_choice.set(_default_preloaded_genome_label(self.preloaded_genomes))
            self._apply_preloaded_genome_choice()
        else:
            self.genome_not_available.set(True)
            self.organism_name.set("Organism")
        self.user_cluster_mode.set("Inferred from FASTA file")
        self.refined_cluster_file.set("")
        self.refined_cluster_sheet.set("")
        self._auto_loaded_cluster_file = ""
        self.add_custom_cds_enable.set(False)
        self.custom_cds_paths.set("")
        self.include_custom_cds_cluster.set(True)
        self.custom_cds_cluster_name.set("custom")
        if self.preloaded_genomes and not self.genome_not_available.get():
            self._apply_preloaded_genome_choice()
        self.usage_basis.set("Relative codon usage")
        self.fasta_codon_range.set(str(CP.SET.get("fasta_codon_range", "all") or "all"))
        self.codon_compare_metric.set("Relative codon usage")
        self.codon_compare_rcu_display.set("Relative codon usage with genome")
        self.codon_compare_statistics.set("None")
        self.codon_compare_plot_style.set("Mean ± SD")
        self.codon_compare_raw_custom_axes.set(False)
        self.codon_compare_raw_xmin.set("")
        self.codon_compare_raw_xmax.set("")
        self.codon_compare_raw_ymin.set("")
        self.codon_compare_raw_ymax.set("")
        self.codon_compare_corr_custom_axes.set(False)
        self.codon_compare_corr_xmin.set("")
        self.codon_compare_corr_xmax.set("")
        self.codon_compare_corr_ymin.set("")
        self.codon_compare_corr_ymax.set("")
        self.codon_compare_caption_size.set("14")
        self.codon_compare_marker_size.set("8.0")
        self.codon_compare_line_width.set("2.0")
        self.codon_compare_aa_spacing.set("1.75")
        self.codon_compare_codon_aa_gap.set("0.65")
        self.codon_compare_codon_gap.set("0.33")
        self.codon_compare_highlight_color.set("#0057D9")
        self.codon_compare_legend_ncol.set("3")
        self.codon_compare_log2_y_scale.set(False)
        self.codon_compare_selected_codons = list(_CODON_ORDER)
        self.codon_compare_highlighted_codons = []
        self._update_codon_compare_selected_codons_status()
        self._update_codon_compare_highlight_status()
        self.fasta_metric_cluster_df = pd.DataFrame()
        self.fasta_metric_cluster_scores_df = pd.DataFrame()
        self.fasta_metric_cluster_path = ""
        self._set_fasta_metric_group_rows_from_configs([])
        self._sync_codon_set_from_usage_basis()
        self.dimred_method.set("UMAP")
        self.cluster_method.set("kmeans")

        self.statistical_test_method.set(_bool_to_choice(CP.RUNTIME_DEFAULTS.get("do_2d_ks", False), "2D Kolmogorov-Smirnov"))
        self.enable_2d_ks.set(bool(CP.RUNTIME_DEFAULTS.get("do_2d_ks", False)))
        self.ks_alpha.set(str(CP.KS_SETTINGS.get("alpha", 0.01)))
        self.ks_method.set(str(CP.KS_SETTINGS.get("method", "binned")))
        self.ks_bins.set(str(CP.KS_SETTINGS.get("bins", 151)))
        self.ks_n_perm.set(str(CP.KS_SETTINGS.get("n_perm", 2000)))
        self.ks_seed.set(str(CP.KS_SETTINGS.get("random_seed", 42)))

        self.functional_scan_method.set("disabled")
        self.enable_david.set(False)
        self.david_email.set(str(CP.SET.get("david_user_email", "")))
        self.david_window_size.set(str(CP.SET.get("david_window_size", 100)))
        self.david_step_size.set(str(CP.SET.get("david_step_size", 50)))
        self.david_wait_time.set(str(CP.SET.get("david_wait_time", 0.0)))
        self.david_max_clusters.set(str(CP.SET.get("david_max_clusters", 3)))
        self.david_min_valid_ids.set(str(CP.SET.get("david_min_valid_ids_per_window", 3)))
        self.david_top_n_hits.set(str(CP.SET.get("david_top_n_hits", 10)))

        self.figure_format.set("png")
        self.enable_main_heatmap.set(True)
        self.main_heatmap_show_fig.set(True)
        self.main_heatmap_custom_aesthetics.set(False)
        self.main_heatmap_custom_axes.set(False)
        self.main_heatmap_dpi.set(str(CP.SET.get("figure_dpi", 300)))
        self.main_heatmap_colormap.set(str(CP.SET.get("heatmap_colormap_name", "parula")))
        self.main_heatmap_fig_width.set(str((CP.SET.get("heatmap_fig_size", (18, 4)) or (18, 4))[0]))
        self.main_heatmap_fig_height.set(str((CP.SET.get("heatmap_fig_size", (18, 4)) or (18, 4))[1]))
        _caxis = CP.SET.get("heatmap_caxis_limits", (-0.5, 2.5)) or (-0.5, 2.5)
        self.main_heatmap_caxis_min.set(str(_caxis[0]))
        self.main_heatmap_caxis_max.set(str(_caxis[1]))
        self.main_heatmap_xtick_every.set(str(CP.SET.get("xtick_every_genes", 500)))
        self.main_heatmap_xmin.set("")
        self.main_heatmap_xmax.set("")
        self.main_heatmap_ymin.set("")
        self.main_heatmap_ymax.set("")
        self.density_figure_dpi.set(str(CP.SET.get("figure_dpi", 300)))
        self.density_panel_w_in.set("5.0")
        self.density_panel_h_in.set("5.0")
        self.plot_rows.set("4")
        self.show_colorbar.set(True)
        self.show_2d_fig.set(True)
        self.color_mode.set("enrichment")
        self.density_subplot_wspace.set("0.20")
        self.density_subplot_hspace.set("0.30")
        self.figure_suptitle.set("")
        self.density_cmap.set("plasma_r")
        self.enrichment_cmap.set("plasma_r")
        self.enable_2d_density_plots.set(True)
        self.plot_cluster_min_genes.set("2")
        self.density_custom_aesthetics.set(False)
        self.density_custom_axes.set(False)
        self.density_xmin.set("")
        self.density_xmax.set("")
        self.density_ymin.set("")
        self.density_ymax.set("")

        self.gchm_enable.set(bool(CP.SET.get("gchm_enable", True)))
        self.heatmap_dpi.set(str(CP.SET.get("gchm_dpi", 300)))
        self.gchm_show_fig.set(True)
        self.gchm_colormap.set(str(CP.SET.get("gchm_colormap", "plasma")))
        self.gchm_sigma.set(str(CP.SET.get("gchm_sigma", GCHM_DEFAULT_SIGMA)))
        self.gchm_spread_factor.set(str(CP.SET.get("gchm_spread_factor", GCHM_DEFAULT_SPREAD_FACTOR)))
        self.gchm_height_per_cluster.set(str(CP.SET.get("gchm_height_per_cluster", 0.3)))
        self.gchm_label_fontsize.set(str(CP.SET.get("gchm_label_fontsize", 10)))
        self.gchm_cmap_min_rel.set(str(CP.SET.get("gchm_cmap_min_rel", 0.2)))
        self.gchm_cmap_max_rel.set(str(CP.SET.get("gchm_cmap_max_rel", 1.0)))
        self.gchm_custom_aesthetics.set(False)
        self.gchm_custom_axes.set(False)
        self.gchm_xmin.set("")
        self.gchm_xmax.set("")
        self.gchm_ymin.set("")
        self.gchm_ymax.set("")

        self.apply_smoothing.set(bool(CP.SET.get("apply_smoothing", True)))
        self.smooth_window_genes.set(str(CP.SET.get("smooth_window_genes", 6)))
        self.apply_binning.set(bool(CP.SET.get("apply_binning", False)))
        self.bin_size_genes.set(str(CP.SET.get("bin_size_genes", 50)))
        self.center_features.set(bool(CP.SET.get("center_features", True)))
        self.scale_features.set(bool(CP.SET.get("scale_features", True)))

        self.enable_codon_usage_plot.set(True)
        self.codon_usage_plot_mode.set("RCU z-scores")
        self.codon_usage_show_fig.set(True)
        self.codon_usage_dpi.set(str(CP.SET.get("figure_dpi", 300)))
        self.codon_panel_w_in.set("5.0")
        self.codon_panel_h_in.set("5.0")
        self.codon_custom_aesthetics.set(False)
        self.codon_custom_axes.set(False)
        self.codon_xmin.set("")
        self.codon_xmax.set("")
        self.codon_ymin.set("")
        self.codon_ymax.set("")

        self.trna_decoding_table_path.set(str(CP.SET.get("trna_decoding_table_path", "")))
        self.trna_decoding_table_sheet.set(str(CP.SET.get("trna_decoding_table_sheet", "")))
        _has_trna_path = bool(self.trna_decoding_table_path.get().strip())
        self.enable_trna_usage.set(bool(CP.SET.get("export_trna_usage_enable", _has_trna_path)))

        self.enable_trna_gene_heatmap.set(bool(CP.SET.get("trna_gene_heatmap_enable", True)))
        self.trna_gene_heatmap_metric.set(str(CP.SET.get("trna_gene_heatmap_metric", "ZTU")))
        self.trna_gene_heatmap_show_fig.set(True)
        self.enable_trna_single_box_codon_heatmap.set(bool(CP.SET.get("trna_single_box_codon_heatmap_enable", True)))
        self.trna_single_box_codon_heatmap_show_fig.set(True)

        self.enable_trna_shift_heatmap.set(bool(CP.SET.get("trna_shift_heatmap_enable", True)))
        self.trna_shift_heatmap_show_fig.set(True)
        self.trna_shift_heatmap_clusters.set(str(CP.SET.get("trna_shift_heatmap_clusters", "all")))
        self.enable_trna_wobble_heatmap.set(bool(CP.SET.get("trna_wobble_heatmap_enable", True)))
        self.trna_wobble_heatmap_show_fig.set(True)
        self.trna_wobble_heatmap_clusters.set(str(CP.SET.get("trna_wobble_heatmap_clusters", "all")))
        self.trna_gene_wobble_plot_kind.set(str(CP.SET.get("trna_gene_wobble_plot_kind", "heatmap") or "heatmap"))
        self.trna_gene_trna_plot_kind.set(str(CP.SET.get("trna_gene_trna_plot_kind", "heatmap") or "heatmap"))
        self.trna_mrna_stability_plot_kind.set(str(CP.SET.get("trna_mrna_stability_plot_kind", "line") or "line"))
        self.trna_gene_wobble_smooth.set(_smooth_method_display(CP.SET.get("trna_gene_wobble_smooth_method", CP.SET.get("trna_gene_wobble_smooth", "running average"))))
        self.trna_gene_wobble_smooth_window.set(str(CP.SET.get("trna_gene_wobble_smooth_window", 40)))
        self.trna_gene_trna_smooth.set(_smooth_method_display(CP.SET.get("trna_gene_trna_smooth_method", CP.SET.get("trna_gene_trna_smooth", "running average"))))
        self.trna_gene_trna_smooth_window.set(str(CP.SET.get("trna_gene_trna_smooth_window", 40)))
        self.trna_mrna_stability_smooth.set(_smooth_method_display(CP.SET.get("trna_mrna_stability_smooth_method", CP.SET.get("trna_mrna_stability_smooth", "running average"))))
        self.trna_mrna_stability_smooth_window.set(str(CP.SET.get("trna_mrna_stability_smooth_window", 40)))
        self.trna_wobble_plot_kind.set(str(CP.SET.get("trna_wobble_plot_kind", "boxplot") or "boxplot"))
        self.trna_shift_plot_kind.set(str(CP.SET.get("trna_shift_plot_kind", "boxplot") or "boxplot"))
        self.trna_modifications_plot_kind.set(str(CP.SET.get("trna_modifications_plot_kind", "boxplot") or "boxplot"))
        self.trna_wobble_pair_stats_test.set(str(CP.SET.get("trna_wobble_pair_stats_test", "Student t-test") or "Student t-test"))
        self.trna_shift_pair_stats_test.set(str(CP.SET.get("trna_shift_pair_stats_test", "Student t-test") or "Student t-test"))
        self.trna_wobble_pair_stats_gap.set(str(CP.SET.get("trna_wobble_pair_stats_gap", 0.05)))
        self.trna_shift_pair_stats_gap.set(str(CP.SET.get("trna_shift_pair_stats_gap", 0.05)))
        self.trna_wobble_boxplot_log2.set("yes" if bool(CP.SET.get("trna_wobble_boxplot_log2", True)) else "no")
        self.trna_shift_boxplot_log2.set("yes" if bool(CP.SET.get("trna_shift_boxplot_log2", True)) else "no")
        self.trna_modifications_boxplot_log2.set("yes" if bool(CP.SET.get("trna_modifications_boxplot_log2", True)) else "no")
        self.enable_trna_mrna_stability.set(bool(CP.SET.get("trna_mrna_stability_enable", True)))
        self.enable_trna_modification_heatmap.set(bool(CP.SET.get("trna_modification_heatmap_enable", True)))
        self.trna_modifications_selection = _canonicalize_gui_plot6_modification_selection(CP.SET.get("trna_modifications_selected_features", None))
        self.trna_modification_aas_selection = _split_selection_text(CP.SET.get("trna_modifications_include_aas", None))  # None -> default manuscript amino-acid set
        self._update_trna_modification_selection_status()
        self._update_trna_modification_aas_status()
        self.trna_modification_heatmap_show_fig.set(True)
        self.trna_secondary_axis_style.set(str(CP.SET.get("trna_secondary_axis_style", "bars") or "bars"))
        self.trna_secondary_axis_alpha.set(str(CP.SET.get("trna_secondary_axis_alpha", 0.22)))
        self.trna_secondary_axis_bar_width.set(str(CP.SET.get("trna_secondary_axis_bar_width", 0.72)))
        self.trna_boxplot_width.set(str(CP.SET.get("trna_boxplot_width", 0.12)))
        self.trna_boxplot_show_points.set(bool(CP.SET.get("trna_boxplot_show_points", True)))
        self.trna_boxplot_point_alpha.set(str(CP.SET.get("trna_boxplot_point_alpha", 0.35)))

        self.trna_supp_heatmaps_customize.set(bool(CP.SET.get("trna_supp_heatmaps_customize", False)))
        self.trna_supp_heatmaps_dpi.set(str(CP.SET.get("trna_supp_heatmaps_dpi", TRNA_SUPP_DETAILS_DEFAULT_DPI) or TRNA_SUPP_DETAILS_DEFAULT_DPI))
        self.trna_supp_heatmaps_fig_width.set(str(CP.SET.get("trna_supp_heatmaps_fig_width", TRNA_SUPP_DETAILS_DEFAULT_FIG_WIDTH) or TRNA_SUPP_DETAILS_DEFAULT_FIG_WIDTH))
        self.trna_supp_heatmaps_fig_height.set(str(CP.SET.get("trna_supp_heatmaps_fig_height", TRNA_SUPP_DETAILS_DEFAULT_FIG_HEIGHT) or TRNA_SUPP_DETAILS_DEFAULT_FIG_HEIGHT))
        self.trna_supp_heatmaps_cell_height.set(str(CP.SET.get("trna_supp_heatmaps_cell_height", TRNA_SUPP_DETAILS_DEFAULT_CELL_HEIGHT) or TRNA_SUPP_DETAILS_DEFAULT_CELL_HEIGHT))
        self.trna_supp_heatmaps_xtick_every_genes.set(str(CP.SET.get("trna_supp_heatmaps_xtick_every_genes", TRNA_SUPP_DETAILS_DEFAULT_XTICK_EVERY) or TRNA_SUPP_DETAILS_DEFAULT_XTICK_EVERY))
        self.trna_supp_heatmaps_ytick_fontsize.set(str(CP.SET.get("trna_supp_heatmaps_ytick_fontsize", TRNA_SUPP_DETAILS_DEFAULT_YTICK_FONTSIZE) or TRNA_SUPP_DETAILS_DEFAULT_YTICK_FONTSIZE))
        self.trna_supp_heatmaps_title_fontsize.set(str(CP.SET.get("trna_supp_heatmaps_title_fontsize", TRNA_SUPP_DETAILS_DEFAULT_TITLE_FONTSIZE) or TRNA_SUPP_DETAILS_DEFAULT_TITLE_FONTSIZE))
        self.trna_supp_heatmaps_xmin.set(str(CP.SET.get("trna_supp_heatmaps_xmin", "") or ""))
        self.trna_supp_heatmaps_xmax.set(str(CP.SET.get("trna_supp_heatmaps_xmax", "") or ""))
        self.trna_supp_heatmaps_ymin.set(str(CP.SET.get("trna_supp_heatmaps_ymin", "") or ""))
        self.trna_supp_heatmaps_ymax.set(str(CP.SET.get("trna_supp_heatmaps_ymax", "") or ""))

        self.trna_shift_heatmaps_customize.set(bool(CP.SET.get("trna_shift_heatmaps_customize", False)))
        self.trna_shift_heatmaps_dpi.set(str(CP.SET.get("trna_shift_heatmaps_dpi", TRNA_SHIFT_DETAILS_DEFAULT_DPI) or TRNA_SHIFT_DETAILS_DEFAULT_DPI))
        self.trna_shift_heatmaps_fig_width.set(str(CP.SET.get("trna_shift_heatmaps_fig_width", TRNA_SHIFT_DETAILS_DEFAULT_FIG_WIDTH) or TRNA_SHIFT_DETAILS_DEFAULT_FIG_WIDTH))
        self.trna_shift_heatmaps_fig_height.set(str(CP.SET.get("trna_shift_heatmaps_fig_height", TRNA_SHIFT_DETAILS_DEFAULT_FIG_HEIGHT) or TRNA_SHIFT_DETAILS_DEFAULT_FIG_HEIGHT))
        self.trna_shift_heatmaps_cell_width.set(str(CP.SET.get("trna_shift_heatmaps_cell_width", TRNA_SHIFT_DETAILS_DEFAULT_CELL_WIDTH) or TRNA_SHIFT_DETAILS_DEFAULT_CELL_WIDTH))
        self.trna_shift_heatmaps_cell_height.set(str(CP.SET.get("trna_shift_heatmaps_cell_height", TRNA_SHIFT_DETAILS_DEFAULT_CELL_HEIGHT) or TRNA_SHIFT_DETAILS_DEFAULT_CELL_HEIGHT))
        self.trna_shift_heatmaps_xtick_fontsize.set(str(CP.SET.get("trna_shift_heatmaps_xtick_fontsize", TRNA_SHIFT_DETAILS_DEFAULT_XTICK_FONTSIZE) or TRNA_SHIFT_DETAILS_DEFAULT_XTICK_FONTSIZE))
        self.trna_shift_heatmaps_ytick_fontsize.set(str(CP.SET.get("trna_shift_heatmaps_ytick_fontsize", TRNA_SHIFT_DETAILS_DEFAULT_YTICK_FONTSIZE) or TRNA_SHIFT_DETAILS_DEFAULT_YTICK_FONTSIZE))
        self.trna_shift_heatmaps_title_fontsize.set(str(CP.SET.get("trna_shift_heatmaps_title_fontsize", TRNA_SHIFT_DETAILS_DEFAULT_TITLE_FONTSIZE) or TRNA_SHIFT_DETAILS_DEFAULT_TITLE_FONTSIZE))
        self.trna_shift_heatmaps_xmin.set(str(CP.SET.get("trna_shift_heatmaps_xmin", "") or ""))
        self.trna_shift_heatmaps_xmax.set(str(CP.SET.get("trna_shift_heatmaps_xmax", "") or ""))
        self.trna_shift_heatmaps_ymin.set(str(CP.SET.get("trna_shift_heatmaps_ymin", "") or ""))
        self.trna_shift_heatmaps_ymax.set(str(CP.SET.get("trna_shift_heatmaps_ymax", "") or ""))

        self.trna_modification_plots_customize.set(bool(CP.SET.get("trna_modification_plots_customize", True)))
        self.trna_modification_plots_dpi.set(str(CP.SET.get("trna_modification_plots_dpi", PLOT6_DEFAULT_DPI) or PLOT6_DEFAULT_DPI))
        self.trna_modification_plots_fig_width.set(str(CP.SET.get("trna_modification_plots_fig_width", PLOT6_DEFAULT_FIG_WIDTH) or PLOT6_DEFAULT_FIG_WIDTH))
        self.trna_modification_plots_fig_height.set(str(CP.SET.get("trna_modification_plots_fig_height", PLOT6_DEFAULT_FIG_HEIGHT) or PLOT6_DEFAULT_FIG_HEIGHT))
        self.trna_modification_plots_caption_size.set(str(CP.SET.get("trna_modification_plots_caption_size", CP.SET.get("trna_modifications_boxplot_caption_size", 17)) or 17))
        self.trna_modification_plots_ymin.set(str(CP.SET.get("trna_modification_plots_ymin", PLOT6_DEFAULT_YMIN) if CP.SET.get("trna_modification_plots_ymin", PLOT6_DEFAULT_YMIN) is not None else PLOT6_DEFAULT_YMIN))
        self.trna_modification_plots_ymax.set(str(CP.SET.get("trna_modification_plots_ymax", PLOT6_DEFAULT_YMAX) if CP.SET.get("trna_modification_plots_ymax", PLOT6_DEFAULT_YMAX) is not None else PLOT6_DEFAULT_YMAX))
        self.trna_modification_plots_group_bar_y.set(str(CP.SET.get("trna_modification_plots_group_bar_y", PLOT6_DEFAULT_GROUP_BAR_Y) or PLOT6_DEFAULT_GROUP_BAR_Y))
        self.trna_modification_plots_group_label_gap.set(str(CP.SET.get("trna_modification_plots_group_label_gap", PLOT6_DEFAULT_GROUP_LABEL_GAP) or PLOT6_DEFAULT_GROUP_LABEL_GAP))
        self.trna_modification_plots_star_offset.set(str(CP.SET.get("trna_modification_plots_star_offset", PLOT6_DEFAULT_STAR_OFFSET) or PLOT6_DEFAULT_STAR_OFFSET))
        self.trna_modification_plots_legend_ncol.set(str(CP.SET.get("trna_modification_plots_legend_ncol", "") or ""))
        self.trna_modification_plots_box_width.set(str(CP.SET.get("trna_modification_plots_box_width", 0.18) or 0.18))

        self._refresh_genome_source_widgets()
        self._refresh_cluster_source_widgets()
        self._refresh_dimred_params()
        self._refresh_cluster_params()
        self._refresh_optional_sections()

    # ------------------------- configuration extraction -------------------------
    def _validate_inputs(self):
        # Do not call _apply_preloaded_genome_choice() here: locating an existing
        # workbook must not trigger companion-file scans or Excel reads.
        fasta = self.fasta_path.get().strip()
        if not fasta:
            raise ValueError("Please select a FASTA file.")
        if not os.path.isfile(fasta):
            raise ValueError(f"FASTA file not found:\n{fasta}")

        _validate_codon_range_text(self.fasta_codon_range.get())

        if self.add_custom_cds_enable.get():
            custom_paths = _split_custom_cds_paths_gui(self.custom_cds_paths.get())
            if not custom_paths:
                raise ValueError("'Add custom CDS to analysis' is enabled, but no additional FASTA/.fna file was selected.")
            if self.add_custom_cds_enable.get() and not self.custom_cds_cluster_name.get().strip():
                raise ValueError("Please enter a gene cluster name for the custom CDS cluster.")
            missing = [p for p in custom_paths if not os.path.isfile(p)]
            if missing:
                raise ValueError("One or more custom CDS FASTA files were not found:\n" + "\n".join(missing))
            bad_ext = [p for p in custom_paths if not p.lower().endswith(FASTA_EXTENSIONS)]
            if bad_ext:
                raise ValueError("Custom CDS files should be FASTA/.fna files:\n" + "\n".join(bad_ext))

        mode = self.user_cluster_mode.get().strip()
        is_refined = _cluster_mode_internal(mode) == "refined"
        is_david_terms = _cluster_mode_terms_source(mode) == "david_gene2terms"
        if is_refined:
            cf = self.refined_cluster_file.get().strip()
            if not cf:
                raise ValueError("Gene clusters is set to 'Provided by user', but no cluster file was selected.")
            if not os.path.isfile(cf):
                raise ValueError(f"Cluster file not found:\n{cf}")

        if is_david_terms and self.david_gene2terms_path.get().strip():
            david_txt = self.david_gene2terms_path.get().strip()
            if not os.path.isfile(david_txt):
                raise ValueError(f"DAVID gene2terms TXT file not found:\n{david_txt}")

        if not self._collect_keyword_groups() and not is_refined and not self._collect_fasta_metric_cluster_configs():
            raise ValueError("No active keyword group or FASTA-derived metric group is available. Enable at least one keyword group or choose a FASTA-derived metric group.")


        trna_required = bool(self.enable_trna_usage.get()) or str(getattr(self, "_current_run_mode", "")).startswith("trna_")
        if trna_required:
            trna = self.trna_decoding_table_path.get().strip() or self._auto_find_decoding_table_for_current_genome()
            if trna:
                self.trna_decoding_table_path.set(trna)
            if not trna:
                raise ValueError("Decoding strategies are enabled, but no decoding Excel workbook was selected or found in the 'Preloaded genomes' subfolder.")
            if not os.path.isfile(trna):
                raise ValueError(f"Decoding strategy Excel workbook not found:\n{trna}")
            try:
                xls = pd.ExcelFile(trna)
                available_sheets = list(xls.sheet_names)
            except Exception as e:
                raise ValueError(f"Could not read the tRNA Excel workbook:\n{trna}\n\n{e}") from e

            decoding_sheet = self.trna_decoding_table_sheet.get().strip()
            if decoding_sheet and decoding_sheet not in available_sheets:
                raise ValueError(
                    "The decoding sheet name was not found in the selected decoding workbook.\n\n"
                    f"Requested: {decoding_sheet}\n"
                    f"Available: {', ' .join(available_sheets)}"
                )
        if not os.path.isfile(PLOTTING_SCRIPT):
            raise ValueError(f"Plotting_Pipeline.py was not found next to the GUI:\n{PLOTTING_SCRIPT}")

    def _collect_runtime_choices(self):
        return dict(
            usage_basis=_usage_internal(self.usage_basis.get()),
            codon_set=self.codon_set.get().strip(),
            dimred_method=_dimred_internal(self.dimred_method.get()),
            cluster_method=self.cluster_method.get().strip().lower(),
            do_2d_ks=(self.statistical_test_method.get().strip() == "2D Kolmogorov-Smirnov"),
            run_david_scan=False,
            compute_trna_usage=bool(self.enable_trna_usage.get()),
            compute_trna_abundance_correlations=False,
            run_2d_density_plots=bool(self.enable_2d_density_plots.get()),
            run_gchm=bool(self.gchm_enable.get()),
            cluster_source=_cluster_mode_internal(self.user_cluster_mode.get()),
            basic_cluster_input_source=_cluster_mode_terms_source(self.user_cluster_mode.get()),
            organism_name=self.organism_name.get().strip() or "Organism",
        )

    def _collect_set_overrides(self):
        fmt = str(self.figure_format.get() or "png").strip().lstrip(".").lower() or "png"
        S = {}
        S["default_root"] = self.default_root.get().strip()
        if not self.genome_not_available.get():
            self._apply_preloaded_genome_choice()
        S["fasta_path"] = self.fasta_path.get().strip()
        S["fasta_codon_range"] = _validate_codon_range_text(self.fasta_codon_range.get())
        metric_cfgs = self._collect_fasta_metric_cluster_configs()
        S["fasta_metric_clusters_enable"] = bool(metric_cfgs)
        S["fasta_metric_cluster_configs"] = metric_cfgs
        S["fasta_metric_cluster_file_path"] = getattr(self, "fasta_metric_cluster_path", "")
        S["custom_cds_enable"] = bool(self.add_custom_cds_enable.get())
        S["custom_cds_paths"] = _split_custom_cds_paths_gui(self.custom_cds_paths.get())
        S["custom_cds_include_as_cluster"] = bool(self.add_custom_cds_enable.get())
        S["custom_cds_cluster_name"] = self.custom_cds_cluster_name.get().strip() or "custom"
        S["basic_cluster_keyword_groups"] = self._collect_keyword_groups()
        S["david_gene2terms_path"] = self.david_gene2terms_path.get().strip()
        self._sync_codon_set_from_usage_basis()
        S["usage_basis"] = _usage_internal(self.usage_basis.get())
        S["codon_set"] = self.codon_set.get().strip()
        S["dimred_method"] = _dimred_internal(self.dimred_method.get())
        S["cluster_method"] = self.cluster_method.get().strip().lower()
        S["organism_name"] = self.organism_name.get().strip() or "Organism"
        S["center_features"] = bool(self.center_features.get())
        S["scale_features"] = bool(self.scale_features.get())
        main_heatmap_custom = bool(self.main_heatmap_custom_aesthetics.get())
        main_heatmap_axes = bool(self.main_heatmap_custom_axes.get())
        S["plot_codon_gene_heatmap_enable"] = bool(self.enable_main_heatmap.get())
        S["plot_codon_gene_heatmap_show_fig"] = True
        S["show_main_pipeline_figures"] = bool(self.enable_main_heatmap.get())
        S["figure_output_format"] = fmt
        S["figure_dpi"] = _safe_int(self.main_heatmap_dpi.get(), CP.SET.get("figure_dpi", 300)) if main_heatmap_custom else int(CP.SET.get("figure_dpi", 300))
        S["heatmap_colormap_name"] = (self.main_heatmap_colormap.get().strip() or "parula") if main_heatmap_custom else str(CP.SET.get("heatmap_colormap_name", "parula"))
        S["heatmap_fig_size"] = (
            _safe_float(self.main_heatmap_fig_width.get(), 18.0),
            _safe_float(self.main_heatmap_fig_height.get(), 4.0),
        ) if main_heatmap_custom else tuple(CP.SET.get("heatmap_fig_size", (18, 4)))
        S["heatmap_caxis_limits"] = (
            _safe_float(self.main_heatmap_caxis_min.get(), -0.5),
            _safe_float(self.main_heatmap_caxis_max.get(), 2.5),
        ) if main_heatmap_custom else tuple(CP.SET.get("heatmap_caxis_limits", (-0.5, 2.5)))
        S["xtick_every_genes"] = _safe_int(self.main_heatmap_xtick_every.get(), CP.SET.get("xtick_every_genes", 500)) if main_heatmap_custom else int(CP.SET.get("xtick_every_genes", 500))
        S["heatmap_xmin"] = _optional_limit_value(self.main_heatmap_xmin.get()) if main_heatmap_axes else None
        S["heatmap_xmax"] = _optional_limit_value(self.main_heatmap_xmax.get()) if main_heatmap_axes else None
        S["heatmap_ymin"] = _optional_limit_value(self.main_heatmap_ymin.get()) if main_heatmap_axes else None
        S["heatmap_ymax"] = _optional_limit_value(self.main_heatmap_ymax.get()) if main_heatmap_axes else None
        S["apply_smoothing"] = bool(self.apply_smoothing.get())
        S["smooth_window_genes"] = _safe_int(self.smooth_window_genes.get(), CP.SET.get("smooth_window_genes", 6))
        S["apply_binning"] = bool(self.apply_binning.get())
        S["bin_size_genes"] = _safe_int(self.bin_size_genes.get(), CP.SET.get("bin_size_genes", 50))
        S["plot_cluster_min_genes"] = _safe_int(self.plot_cluster_min_genes.get(), 2)
        density_custom_axes = bool(self.density_custom_axes.get())
        S["scatter_xmin"] = _optional_limit_value(self.density_xmin.get()) if density_custom_axes else None
        S["scatter_xmax"] = _optional_limit_value(self.density_xmax.get()) if density_custom_axes else None
        S["scatter_ymin"] = _optional_limit_value(self.density_ymin.get()) if density_custom_axes else None
        S["scatter_ymax"] = _optional_limit_value(self.density_ymax.get()) if density_custom_axes else None
        S["do_2d_ks"] = (self.statistical_test_method.get().strip() == "2D Kolmogorov-Smirnov")
        S["export_trna_usage_enable"] = bool(self.enable_trna_usage.get())
        S["export_trna_abundance_correlation_enable"] = False
        S["trna_decoding_table_path"] = self.trna_decoding_table_path.get().strip()
        S["trna_decoding_table_sheet"] = self.trna_decoding_table_sheet.get().strip()
        S["trna_abundance_sheet"] = ""
        S["trna_abundance_corr_show_fig"] = False
        S["trna_abundance_corr_dpi"] = int(CP.SET.get("figure_dpi", 300))
        S["trna_abundance_heatmap_metric"] = "ZTU"
        S["trna_abundance_scatter_metric"] = "ZTU"
        S["trna_abundance_heatmap_clusters"] = "all"
        S["trna_abundance_scatter_clusters"] = "all"
        S["trna_abundance_scatter_yscale"] = "linear"
        S["trna_abundance_scatter_show_fig"] = False
        S["trna_gene_heatmap_enable"] = bool(self.enable_trna_usage.get() and self.enable_trna_gene_heatmap.get())
        S["trna_gene_heatmap_metric"] = str(self.trna_gene_heatmap_metric.get() or "ZTU").upper()
        S["trna_gene_heatmap_show_fig"] = True
        S["trna_pairing_heatmap_enable"] = False
        S["trna_pairing_heatmap_show_fig"] = False
        S["trna_single_box_codon_heatmap_enable"] = bool(self.enable_trna_usage.get() and self.enable_trna_single_box_codon_heatmap.get())
        S["trna_single_box_codon_heatmap_show_fig"] = True
        S["trna_shift_heatmap_enable"] = bool(self.enable_trna_usage.get() and self.enable_trna_shift_heatmap.get())
        S["trna_shift_heatmap_show_fig"] = True
        S["trna_shift_heatmap_clusters"] = "all"
        S["trna_wobble_heatmap_enable"] = bool(self.enable_trna_usage.get() and self.enable_trna_wobble_heatmap.get())
        S["trna_wobble_heatmap_show_fig"] = True
        S["trna_wobble_heatmap_clusters"] = "all"
        S["trna_modification_heatmap_enable"] = bool(self.enable_trna_usage.get() and self.enable_trna_modification_heatmap.get())
        S["trna_modification_heatmap_show_fig"] = True
        S["trna_gene_wobble_plot_kind"] = str(self.trna_gene_wobble_plot_kind.get() or "heatmap").strip().lower()
        S["trna_gene_trna_plot_kind"] = str(self.trna_gene_trna_plot_kind.get() or "heatmap").strip().lower()
        S["trna_mrna_stability_enable"] = bool(self.enable_trna_usage.get() and self.enable_trna_mrna_stability.get())
        S["trna_mrna_stability_plot_kind"] = str(self.trna_mrna_stability_plot_kind.get() or "line").strip().lower()
        S["trna_gene_wobble_smooth_method"] = _smooth_method_display(self.trna_gene_wobble_smooth.get())
        S["trna_gene_wobble_smooth"] = _smooth_method_to_bool(self.trna_gene_wobble_smooth.get())
        S["trna_gene_wobble_smooth_window"] = _safe_int(self.trna_gene_wobble_smooth_window.get(), 40)
        S["trna_gene_trna_smooth_method"] = _smooth_method_display(self.trna_gene_trna_smooth.get())
        S["trna_gene_trna_smooth"] = _smooth_method_to_bool(self.trna_gene_trna_smooth.get())
        S["trna_gene_trna_smooth_window"] = _safe_int(self.trna_gene_trna_smooth_window.get(), 40)
        S["trna_mrna_stability_smooth_method"] = _smooth_method_display(self.trna_mrna_stability_smooth.get())
        S["trna_mrna_stability_smooth"] = _smooth_method_to_bool(self.trna_mrna_stability_smooth.get())
        S["trna_mrna_stability_smooth_window"] = _safe_int(self.trna_mrna_stability_smooth_window.get(), 100)
        S["trna_gene_wobble_caption_size"] = _safe_int(self.trna_gene_wobble_caption_size.get(), 13)
        S["trna_gene_trna_caption_size"] = _safe_int(self.trna_gene_trna_caption_size.get(), 13)
        S["trna_mrna_stability_caption_size"] = _safe_int(self.trna_mrna_stability_caption_size.get(), 13)
        S["trna_wobble_boxplot_caption_size"] = _safe_int(self.trna_wobble_boxplot_caption_size.get(), 13)
        S["trna_shift_boxplot_caption_size"] = _safe_int(self.trna_shift_boxplot_caption_size.get(), 13)
        S["trna_modifications_boxplot_caption_size"] = _safe_int(self.trna_modifications_boxplot_caption_size.get(), 17)
        S["decoding_reference_cluster"] = str(self.decoding_reference_cluster.get() or "").strip()
        S["trna_wobble_stats_test"] = str(self.trna_wobble_stats_test.get() or "none").strip()
        S["trna_shift_stats_test"] = str(self.trna_shift_stats_test.get() or "none").strip()
        S["trna_modifications_stats_test"] = str(self.trna_modifications_stats_test.get() or "none").strip()
        S["trna_wobble_pair_stats_test"] = str(self.trna_wobble_pair_stats_test.get() or "Student t-test").strip()
        S["trna_shift_pair_stats_test"] = str(self.trna_shift_pair_stats_test.get() or "Student t-test").strip()
        S["trna_wobble_pair_stats_gap"] = _safe_float(self.trna_wobble_pair_stats_gap.get(), 0.05)
        S["trna_shift_pair_stats_gap"] = _safe_float(self.trna_shift_pair_stats_gap.get(), 0.05)
        S["trna_wobble_plot_kind"] = str(self.trna_wobble_plot_kind.get() or "boxplot").strip().lower()
        S["trna_shift_plot_kind"] = str(self.trna_shift_plot_kind.get() or "boxplot").strip().lower()
        S["trna_modifications_plot_kind"] = str(self.trna_modifications_plot_kind.get() or "boxplot").strip().lower()
        S["trna_secondary_axis_style"] = str(self.trna_secondary_axis_style.get() or "bars").strip().lower()
        S["trna_secondary_axis_alpha"] = _safe_float(self.trna_secondary_axis_alpha.get(), 0.22)
        S["trna_secondary_axis_bar_width"] = _safe_float(self.trna_secondary_axis_bar_width.get(), 0.72)
        S["trna_boxplot_width"] = _safe_float(self.trna_boxplot_width.get(), 0.12)
        S["trna_boxplot_show_points"] = bool(self.trna_boxplot_show_points.get())
        S["trna_boxplot_point_alpha"] = _safe_float(self.trna_boxplot_point_alpha.get(), 0.35)
        S["trna_boxplot_point_size"] = _safe_float(self.trna_boxplot_point_size.get(), 10.5)
        S["trna_wobble_boxplot_style"] = str(self.trna_wobble_plot_kind.get() or "boxplot").strip().lower()
        S["trna_shift_boxplot_style"] = str(self.trna_shift_plot_kind.get() or "boxplot").strip().lower()
        S["trna_modifications_boxplot_style"] = str(self.trna_modifications_plot_kind.get() or "boxplot").strip().lower()
        S["trna_wobble_boxplot_log2"] = str(self.trna_wobble_boxplot_log2.get()).strip().lower() in {"yes", "true", "1", "on"}
        S["trna_shift_boxplot_log2"] = str(self.trna_shift_boxplot_log2.get()).strip().lower() in {"yes", "true", "1", "on"}
        S["trna_modifications_boxplot_log2"] = str(self.trna_modifications_boxplot_log2.get()).strip().lower() in {"yes", "true", "1", "on"}
        S["trna_wobble_boxplot_ymin"] = _optional_limit_value(self.trna_wobble_boxplot_ymin.get())
        S["trna_wobble_boxplot_ymax"] = _optional_limit_value(self.trna_wobble_boxplot_ymax.get())
        S["trna_shift_boxplot_ymin"] = _optional_limit_value(self.trna_shift_boxplot_ymin.get())
        S["trna_shift_boxplot_ymax"] = _optional_limit_value(self.trna_shift_boxplot_ymax.get())
        S["trna_modifications_boxplot_ymin"] = _optional_limit_value(self.trna_modifications_boxplot_ymin.get())
        S["trna_modifications_boxplot_ymax"] = _optional_limit_value(self.trna_modifications_boxplot_ymax.get())
        S["trna_wobble_exclude_outliers"] = str(self.trna_wobble_exclude_outliers.get()).strip().lower() in {"yes", "true", "1", "on"}
        S["trna_shift_exclude_outliers"] = str(self.trna_shift_exclude_outliers.get()).strip().lower() in {"yes", "true", "1", "on"}
        S["trna_modifications_exclude_outliers"] = str(self.trna_modifications_exclude_outliers.get()).strip().lower() in {"yes", "true", "1", "on"}
        S["trna_wobble_outlier_sd"] = _safe_float(self.trna_wobble_outlier_sd.get(), 3.0)
        S["trna_shift_outlier_sd"] = _safe_float(self.trna_shift_outlier_sd.get(), 3.0)
        S["trna_modifications_outlier_sd"] = _safe_float(self.trna_modifications_outlier_sd.get(), 3.0)
        S["trna_modifications_feature_mode"] = str(self.trna_modifications_feature_mode.get() or "modifications").strip().lower()
        S["trna_modifications_selected_features"] = None if self.trna_modifications_selection is None else list(self.trna_modifications_selection)
        S["trna_modifications_include_aas"] = None if self.trna_modification_aas_selection is None else list(self.trna_modification_aas_selection)
        S["trna_modifications_assignment_models"] = "conservative,permissive"
        S["trna_shift_heatmap_log2_colorbar"] = bool(self.trna_shift_heatmap_log2_colorbar.get())
        S["trna_wobble_heatmap_log2_colorbar"] = bool(self.trna_wobble_heatmap_log2_colorbar.get())
        S["trna_shift_heatmap_bracket_type"] = self.trna_shift_heatmap_bracket_type.get().strip() or "square"
        S["trna_wobble_heatmap_bracket_type"] = self.trna_wobble_heatmap_bracket_type.get().strip() or "square"
        S["trna_shift_heatmap_bracket_x"] = _safe_float(self.trna_shift_heatmap_bracket_x.get(), -0.20)
        S["trna_shift_heatmap_label_x"] = _safe_float(self.trna_shift_heatmap_label_x.get(), -0.31)
        S["trna_wobble_heatmap_bracket_x"] = _safe_float(self.trna_wobble_heatmap_bracket_x.get(), -0.17)
        S["trna_wobble_heatmap_label_x"] = _safe_float(self.trna_wobble_heatmap_label_x.get(), -0.31)

        trna_supp_custom = bool(self.trna_supp_heatmaps_customize.get())
        S["trna_supp_heatmaps_customize"] = trna_supp_custom
        S["trna_supp_heatmaps_dpi"] = _safe_int(self.trna_supp_heatmaps_dpi.get(), CP.SET.get("figure_dpi", 300)) if trna_supp_custom and self.trna_supp_heatmaps_dpi.get().strip() else None
        S["trna_supp_heatmaps_fig_width"] = _safe_float(self.trna_supp_heatmaps_fig_width.get(), CP.SET.get("heatmap_fig_size", (18, 4))[0]) if trna_supp_custom and self.trna_supp_heatmaps_fig_width.get().strip() else None
        S["trna_supp_heatmaps_fig_height"] = _safe_float(self.trna_supp_heatmaps_fig_height.get(), CP.SET.get("heatmap_fig_size", (18, 4))[1]) if trna_supp_custom and self.trna_supp_heatmaps_fig_height.get().strip() else None
        if S["trna_supp_heatmaps_fig_height"] is not None and float(S["trna_supp_heatmaps_fig_height"]) <= 0:
            S["trna_supp_heatmaps_fig_height"] = None
        S["trna_supp_heatmaps_cell_height"] = _safe_float(self.trna_supp_heatmaps_cell_height.get(), 0.0) if trna_supp_custom and self.trna_supp_heatmaps_cell_height.get().strip() else None
        if S["trna_supp_heatmaps_cell_height"] is not None and float(S["trna_supp_heatmaps_cell_height"]) <= 0:
            S["trna_supp_heatmaps_cell_height"] = None
        if S["trna_supp_heatmaps_fig_height"] is not None:
            S["trna_supp_heatmaps_cell_height"] = None
        S["trna_supp_heatmaps_xtick_every_genes"] = _safe_int(self.trna_supp_heatmaps_xtick_every_genes.get(), CP.SET.get("xtick_every_genes", 500)) if trna_supp_custom and self.trna_supp_heatmaps_xtick_every_genes.get().strip() else None
        S["trna_supp_heatmaps_ytick_fontsize"] = _safe_int(self.trna_supp_heatmaps_ytick_fontsize.get(), CP.SET.get("font_size_yticks", 3)) if trna_supp_custom and self.trna_supp_heatmaps_ytick_fontsize.get().strip() else None
        S["trna_supp_heatmaps_title_fontsize"] = _safe_int(self.trna_supp_heatmaps_title_fontsize.get(), CP.SET.get("font_size_titles", 10)) if trna_supp_custom and self.trna_supp_heatmaps_title_fontsize.get().strip() else None
        S["trna_supp_heatmaps_xmin"] = _optional_limit_value(self.trna_supp_heatmaps_xmin.get()) if trna_supp_custom else None
        S["trna_supp_heatmaps_xmax"] = _optional_limit_value(self.trna_supp_heatmaps_xmax.get()) if trna_supp_custom else None
        S["trna_supp_heatmaps_ymin"] = _optional_limit_value(self.trna_supp_heatmaps_ymin.get()) if trna_supp_custom else None
        S["trna_supp_heatmaps_ymax"] = _optional_limit_value(self.trna_supp_heatmaps_ymax.get()) if trna_supp_custom else None

        trna_shift_custom = bool(self.trna_shift_heatmaps_customize.get())
        S["trna_shift_heatmaps_customize"] = trna_shift_custom
        S["trna_shift_heatmaps_dpi"] = _safe_int(self.trna_shift_heatmaps_dpi.get(), CP.SET.get("figure_dpi", 300)) if trna_shift_custom and self.trna_shift_heatmaps_dpi.get().strip() else None
        S["trna_shift_heatmaps_fig_width"] = _safe_float(self.trna_shift_heatmaps_fig_width.get(), TRNA_SHIFT_DETAILS_DEFAULT_FIG_WIDTH) if trna_shift_custom and self.trna_shift_heatmaps_fig_width.get().strip() else None
        if S["trna_shift_heatmaps_fig_width"] is not None and float(S["trna_shift_heatmaps_fig_width"]) <= 0:
            S["trna_shift_heatmaps_fig_width"] = None
        S["trna_shift_heatmaps_fig_height"] = _safe_float(self.trna_shift_heatmaps_fig_height.get(), TRNA_SHIFT_DETAILS_DEFAULT_FIG_HEIGHT) if trna_shift_custom and self.trna_shift_heatmaps_fig_height.get().strip() else None
        if S["trna_shift_heatmaps_fig_height"] is not None and float(S["trna_shift_heatmaps_fig_height"]) <= 0:
            S["trna_shift_heatmaps_fig_height"] = None
        S["trna_shift_heatmaps_cell_width"] = _safe_float(self.trna_shift_heatmaps_cell_width.get(), TRNA_SHIFT_DETAILS_DEFAULT_CELL_WIDTH) if trna_shift_custom and self.trna_shift_heatmaps_cell_width.get().strip() else None
        if S["trna_shift_heatmaps_cell_width"] is not None and float(S["trna_shift_heatmaps_cell_width"]) <= 0:
            S["trna_shift_heatmaps_cell_width"] = None
        if S["trna_shift_heatmaps_fig_width"] is not None:
            S["trna_shift_heatmaps_cell_width"] = None
        S["trna_shift_heatmaps_cell_height"] = _safe_float(self.trna_shift_heatmaps_cell_height.get(), TRNA_SHIFT_DETAILS_DEFAULT_CELL_HEIGHT) if trna_shift_custom and self.trna_shift_heatmaps_cell_height.get().strip() else None
        if S["trna_shift_heatmaps_cell_height"] is not None and float(S["trna_shift_heatmaps_cell_height"]) <= 0:
            S["trna_shift_heatmaps_cell_height"] = None
        if S["trna_shift_heatmaps_fig_height"] is not None:
            S["trna_shift_heatmaps_cell_height"] = None
        S["trna_shift_heatmaps_xtick_fontsize"] = _safe_int(self.trna_shift_heatmaps_xtick_fontsize.get(), TRNA_SHIFT_DETAILS_DEFAULT_XTICK_FONTSIZE) if trna_shift_custom and self.trna_shift_heatmaps_xtick_fontsize.get().strip() else None
        S["trna_shift_heatmaps_ytick_fontsize"] = _safe_int(self.trna_shift_heatmaps_ytick_fontsize.get(), TRNA_SHIFT_DETAILS_DEFAULT_YTICK_FONTSIZE) if trna_shift_custom and self.trna_shift_heatmaps_ytick_fontsize.get().strip() else None
        S["trna_shift_heatmaps_title_fontsize"] = _safe_int(self.trna_shift_heatmaps_title_fontsize.get(), TRNA_SHIFT_DETAILS_DEFAULT_TITLE_FONTSIZE) if trna_shift_custom and self.trna_shift_heatmaps_title_fontsize.get().strip() else None
        S["trna_shift_heatmaps_xmin"] = _optional_limit_value(self.trna_shift_heatmaps_xmin.get()) if trna_shift_custom else None
        S["trna_shift_heatmaps_xmax"] = _optional_limit_value(self.trna_shift_heatmaps_xmax.get()) if trna_shift_custom else None
        S["trna_shift_heatmaps_ymin"] = _optional_limit_value(self.trna_shift_heatmaps_ymin.get()) if trna_shift_custom else None
        S["trna_shift_heatmaps_ymax"] = _optional_limit_value(self.trna_shift_heatmaps_ymax.get()) if trna_shift_custom else None

        trna_mod_custom = bool(self.trna_modification_plots_customize.get())
        S["trna_modification_plots_customize"] = trna_mod_custom
        # Plot 6 figure-detail presets are always saved and used as the default
        # behavior, even when the checkbox is left unchecked.
        S["trna_modification_plots_dpi"] = _safe_int(self.trna_modification_plots_dpi.get(), CP.SET.get("figure_dpi", 300)) if self.trna_modification_plots_dpi.get().strip() else PLOT6_DEFAULT_DPI
        S["trna_modification_plots_fig_width"] = _safe_float(self.trna_modification_plots_fig_width.get(), PLOT6_DEFAULT_FIG_WIDTH) if self.trna_modification_plots_fig_width.get().strip() else PLOT6_DEFAULT_FIG_WIDTH
        if S["trna_modification_plots_fig_width"] is not None and float(S["trna_modification_plots_fig_width"]) <= 0:
            S["trna_modification_plots_fig_width"] = PLOT6_DEFAULT_FIG_WIDTH
        S["trna_modification_plots_fig_height"] = _safe_float(self.trna_modification_plots_fig_height.get(), PLOT6_DEFAULT_FIG_HEIGHT) if self.trna_modification_plots_fig_height.get().strip() else PLOT6_DEFAULT_FIG_HEIGHT
        if S["trna_modification_plots_fig_height"] is not None and float(S["trna_modification_plots_fig_height"]) <= 0:
            S["trna_modification_plots_fig_height"] = PLOT6_DEFAULT_FIG_HEIGHT
        S["trna_modification_plots_caption_size"] = _safe_int(self.trna_modification_plots_caption_size.get(), 17) if self.trna_modification_plots_caption_size.get().strip() else 17
        S["trna_modification_plots_ymin"] = _optional_limit_value(self.trna_modification_plots_ymin.get())
        S["trna_modification_plots_ymax"] = _optional_limit_value(self.trna_modification_plots_ymax.get())
        S["trna_modification_plots_group_bar_y"] = _optional_limit_value(self.trna_modification_plots_group_bar_y.get())
        S["trna_modification_plots_group_label_gap"] = _optional_limit_value(self.trna_modification_plots_group_label_gap.get())
        S["trna_modification_plots_star_offset"] = _optional_limit_value(self.trna_modification_plots_star_offset.get())
        S["trna_modification_plots_legend_ncol"] = _safe_int(self.trna_modification_plots_legend_ncol.get(), 4) if self.trna_modification_plots_legend_ncol.get().strip() else None
        if S["trna_modification_plots_legend_ncol"] is not None and int(S["trna_modification_plots_legend_ncol"]) <= 0:
            S["trna_modification_plots_legend_ncol"] = None
        S["trna_modification_plots_box_width"] = _safe_float(self.trna_modification_plots_box_width.get(), 0.18) if self.trna_modification_plots_box_width.get().strip() else 0.18
        if S["trna_modification_plots_box_width"] is not None and float(S["trna_modification_plots_box_width"]) <= 0:
            S["trna_modification_plots_box_width"] = 0.18

        S["david_user_email"] = self.david_email.get().strip()
        S["david_window_size"] = _safe_int(self.david_window_size.get(), CP.SET.get("david_window_size", 100))
        S["david_step_size"] = _safe_int(self.david_step_size.get(), CP.SET.get("david_step_size", 50))
        S["david_wait_time"] = _safe_float(self.david_wait_time.get(), CP.SET.get("david_wait_time", 0.0))
        S["david_max_clusters"] = _safe_int(self.david_max_clusters.get(), CP.SET.get("david_max_clusters", 3))
        S["david_min_valid_ids_per_window"] = _safe_int(self.david_min_valid_ids.get(), CP.SET.get("david_min_valid_ids_per_window", 3))
        S["david_top_n_hits"] = _safe_int(self.david_top_n_hits.get(), CP.SET.get("david_top_n_hits", 10))
        S["david_plot_format"] = fmt

        gchm_custom_aesthetics = bool(self.gchm_custom_aesthetics.get())
        S["gchm_enable"] = bool(self.gchm_enable.get())
        S["gchm_show_fig"] = True
        S["gchm_colormap"] = (self.gchm_colormap.get().strip() or "plasma") if gchm_custom_aesthetics else str(CP.SET.get("gchm_colormap", "plasma"))
        S["gchm_sigma"] = _safe_float(self.gchm_sigma.get(), CP.SET.get("gchm_sigma", GCHM_DEFAULT_SIGMA)) if gchm_custom_aesthetics else float(CP.SET.get("gchm_sigma", GCHM_DEFAULT_SIGMA))
        S["gchm_spread_factor"] = _safe_float(self.gchm_spread_factor.get(), CP.SET.get("gchm_spread_factor", GCHM_DEFAULT_SPREAD_FACTOR)) if gchm_custom_aesthetics else float(CP.SET.get("gchm_spread_factor", GCHM_DEFAULT_SPREAD_FACTOR))
        S["gchm_height_per_cluster"] = _safe_float(self.gchm_height_per_cluster.get(), CP.SET.get("gchm_height_per_cluster", 0.3)) if gchm_custom_aesthetics else float(CP.SET.get("gchm_height_per_cluster", 0.3))
        S["gchm_label_fontsize"] = _safe_int(self.gchm_label_fontsize.get(), CP.SET.get("gchm_label_fontsize", 10)) if gchm_custom_aesthetics else int(CP.SET.get("gchm_label_fontsize", 10))
        S["gchm_dpi"] = _safe_int(self.heatmap_dpi.get(), CP.SET.get("gchm_dpi", 300)) if gchm_custom_aesthetics else int(CP.SET.get("gchm_dpi", 300))
        S["gchm_cmap_min_rel"] = _safe_float(self.gchm_cmap_min_rel.get(), CP.SET.get("gchm_cmap_min_rel", 0.2)) if gchm_custom_aesthetics else float(CP.SET.get("gchm_cmap_min_rel", 0.2))
        S["gchm_cmap_max_rel"] = _safe_float(self.gchm_cmap_max_rel.get(), CP.SET.get("gchm_cmap_max_rel", 1.0)) if gchm_custom_aesthetics else float(CP.SET.get("gchm_cmap_max_rel", 1.0))
        S["gchm_output_filename"] = f"gene_cluster_heatmap_KS.{fmt}"

        dimred_method = _dimred_internal(self.dimred_method.get())
        for key, var in self.dimred_param_vars.get(dimred_method, {}).items():
            if isinstance(var, tk.BooleanVar):
                S[key] = bool(var.get())
            else:
                default = CP.SET.get(key)
                if isinstance(default, bool):
                    S[key] = _safe_bool(var.get())
                elif isinstance(default, int):
                    S[key] = _safe_int(var.get(), default)
                elif isinstance(default, float):
                    S[key] = _safe_float(var.get(), default)
                else:
                    S[key] = str(var.get()).strip()

        cluster_method = self.cluster_method.get().strip().lower()
        for key, var in self.cluster_param_vars.get(cluster_method, {}).items():
            if isinstance(var, tk.BooleanVar):
                S[key] = bool(var.get())
            else:
                default = CP.SET.get(key)
                if isinstance(default, bool):
                    S[key] = _safe_bool(var.get())
                elif isinstance(default, int):
                    S[key] = _safe_int(var.get(), default)
                elif isinstance(default, float):
                    S[key] = _safe_float(var.get(), default)
                else:
                    S[key] = str(var.get()).strip()
        return S

    def _collect_ks_overrides(self):
        return dict(
            alpha=_safe_float(self.ks_alpha.get(), CP.KS_SETTINGS.get("alpha", 0.01)),
            method=self.ks_method.get().strip() or CP.KS_SETTINGS.get("method", "binned"),
            bins=_safe_int(self.ks_bins.get(), CP.KS_SETTINGS.get("bins", 151)),
            n_perm=_safe_int(self.ks_n_perm.get(), CP.KS_SETTINGS.get("n_perm", 2000)),
            random_seed=_safe_int(self.ks_seed.get(), CP.KS_SETTINGS.get("random_seed", 42)),
        )

    def _collect_pipeline_overrides(self):
        fmt = str(self.figure_format.get() or "png").strip().lstrip(".").lower() or "png"
        any_plotting = bool(self.enable_2d_density_plots.get() or self.gchm_enable.get() or self.enable_codon_usage_plot.get())
        codon_mode = USER_CODON_MODE_TO_INTERNAL.get(self.codon_usage_plot_mode.get().strip(), "Z") if self.enable_codon_usage_plot.get() else "NONE"
        density_custom_aesthetics = bool(self.density_custom_aesthetics.get())
        density_custom_axes = bool(self.density_custom_axes.get())
        gchm_custom_aesthetics = bool(self.gchm_custom_aesthetics.get())
        gchm_custom_axes = bool(self.gchm_custom_axes.get())
        codon_custom_aesthetics = bool(self.codon_custom_aesthetics.get())
        codon_custom_axes = bool(self.codon_custom_axes.get())
        return dict(
            auto_run_plotting_pipeline=any_plotting,
            plotting_pipeline_script_path=PLOTTING_SCRIPT,
            plot_max_nrows=_safe_int(self.plot_rows.get(), 4),
            MAX_NROWS=_safe_int(self.plot_rows.get(), 4),
            codon_usage_plot_mode=codon_mode,
            RUN_2D_DENSITY_PLOTS=bool(self.enable_2d_density_plots.get()),
            INCLUDE_GENOMIC_DENSITY_PANEL=bool(self.include_genomic_density_map.get()),
            PNG_DPI=_safe_int(self.density_figure_dpi.get(), CP.SET.get("figure_dpi", 300)) if density_custom_aesthetics else int(CP.SET.get("figure_dpi", 300)),
            FIGURE_FORMAT=fmt,
            OUTPUT_EXT=fmt,
            PANEL_W_IN=_safe_float(self.density_panel_w_in.get(), 5.0) if density_custom_aesthetics else 5.0,
            PANEL_H_IN=_safe_float(self.density_panel_h_in.get(), 5.0) if density_custom_aesthetics else 5.0,
            SHOW_FIG=True,
            SHOW_COLORBAR=bool(self.show_colorbar.get()),
            COLOR_MODE=self.color_mode.get().strip().lower() or "enrichment",
            SUBPLOT_WSPACE=_safe_float(self.density_subplot_wspace.get(), 0.20) if density_custom_aesthetics else 0.20,
            SUBPLOT_HSPACE=_safe_float(self.density_subplot_hspace.get(), 0.30) if density_custom_aesthetics else 0.30,
            FIGURE_SUPTITLE=self.figure_suptitle.get().strip() if density_custom_aesthetics else "",
            DENSITY_CMAP_NAME=(self.density_cmap.get().strip() or "plasma_r") if density_custom_aesthetics else "plasma_r",
            ENRICHMENT_CMAP_NAME=(self.enrichment_cmap.get().strip() or "plasma_r") if density_custom_aesthetics else "plasma_r",
            DENSITY_XMIN=_optional_limit_value(self.density_xmin.get()) if density_custom_axes else None,
            DENSITY_XMAX=_optional_limit_value(self.density_xmax.get()) if density_custom_axes else None,
            DENSITY_YMIN=_optional_limit_value(self.density_ymin.get()) if density_custom_axes else None,
            DENSITY_YMAX=_optional_limit_value(self.density_ymax.get()) if density_custom_axes else None,
            RUN_GENE_CLUSTER_LOCALIZATION_HEATMAP=bool(self.gchm_enable.get()),
            GCHM_COLORMAP=(self.gchm_colormap.get().strip() or "plasma") if gchm_custom_aesthetics else str(CP.SET.get("gchm_colormap", "plasma")),
            GCHM_SIGMA=_safe_float(self.gchm_sigma.get(), CP.SET.get("gchm_sigma", GCHM_DEFAULT_SIGMA)) if gchm_custom_aesthetics else float(CP.SET.get("gchm_sigma", GCHM_DEFAULT_SIGMA)),
            GCHM_SPREAD_FACTOR=_safe_float(self.gchm_spread_factor.get(), CP.SET.get("gchm_spread_factor", GCHM_DEFAULT_SPREAD_FACTOR)) if gchm_custom_aesthetics else float(CP.SET.get("gchm_spread_factor", GCHM_DEFAULT_SPREAD_FACTOR)),
            GCHM_HEIGHT_PER_CLUSTER=_safe_float(self.gchm_height_per_cluster.get(), CP.SET.get("gchm_height_per_cluster", 0.3)) if gchm_custom_aesthetics else float(CP.SET.get("gchm_height_per_cluster", 0.3)),
            GCHM_LABEL_FONTSIZE=_safe_int(self.gchm_label_fontsize.get(), CP.SET.get("gchm_label_fontsize", 10)) if gchm_custom_aesthetics else int(CP.SET.get("gchm_label_fontsize", 10)),
            GCHM_DPI=_safe_int(self.heatmap_dpi.get(), CP.SET.get("gchm_dpi", 300)) if gchm_custom_aesthetics else int(CP.SET.get("gchm_dpi", 300)),
            GCHM_CMAP_MIN_REL=_safe_float(self.gchm_cmap_min_rel.get(), CP.SET.get("gchm_cmap_min_rel", 0.2)) if gchm_custom_aesthetics else float(CP.SET.get("gchm_cmap_min_rel", 0.2)),
            GCHM_CMAP_MAX_REL=_safe_float(self.gchm_cmap_max_rel.get(), CP.SET.get("gchm_cmap_max_rel", 1.0)) if gchm_custom_aesthetics else float(CP.SET.get("gchm_cmap_max_rel", 1.0)),
            GCHM_OUTPUT_FILENAME=f"gene_cluster_heatmap_KS.{fmt}",
            GCHM_OUTPUT_BASENAME="gene_cluster_heatmap_KS",
            GCHM_SHOW_FIG=True,
            GCHM_XMIN=_optional_limit_value(self.gchm_xmin.get()) if gchm_custom_axes else None,
            GCHM_XMAX=_optional_limit_value(self.gchm_xmax.get()) if gchm_custom_axes else None,
            GCHM_YMIN=_optional_limit_value(self.gchm_ymin.get()) if gchm_custom_axes else None,
            GCHM_YMAX=_optional_limit_value(self.gchm_ymax.get()) if gchm_custom_axes else None,
            CODON_USAGE_PLOT_MODE=codon_mode,
            CODON_USAGE_PNG_DPI=_safe_int(self.codon_usage_dpi.get(), CP.SET.get("figure_dpi", 300)) if codon_custom_aesthetics else int(CP.SET.get("figure_dpi", 300)),
            CODON_USAGE_PANEL_W_IN=_safe_float(self.codon_panel_w_in.get(), 5.0) if codon_custom_aesthetics else 5.0,
            CODON_USAGE_PANEL_H_IN=_safe_float(self.codon_panel_h_in.get(), 5.0) if codon_custom_aesthetics else 5.0,
            CODON_USAGE_SHOW_FIG=True,
            CODON_USAGE_OUTPUT_FORMAT=fmt,
            CODON_USAGE_OUTPUT_BASENAME="Average codon usage per cluster",
            CODON_USAGE_XMIN=_optional_limit_value(self.codon_xmin.get()) if codon_custom_axes else None,
            CODON_USAGE_XMAX=_optional_limit_value(self.codon_xmax.get()) if codon_custom_axes else None,
            CODON_USAGE_YMIN=_optional_limit_value(self.codon_ymin.get()) if codon_custom_axes else None,
            CODON_USAGE_YMAX=_optional_limit_value(self.codon_ymax.get()) if codon_custom_axes else None,
        )

    # ------------------------- cluster FASTA extraction -------------------------
    def _open_extract_clusters_dialog(self):
        available = self._get_available_cluster_names(show_errors=True)
        if not available:
            messagebox.showinfo("Extract FASTA", "No clusters are currently available. Check the selected cluster source and keyword groups.")
            return

        win = tk.Toplevel(self.root)
        win.title("Clusters to extract as FASTA")
        win.geometry("560x640")
        win.transient(self.root)
        win.grab_set()

        outer = ttk.Frame(win, padding=10)
        outer.pack(fill="both", expand=True)
        ttk.Label(
            outer,
            text="Select the cluster(s) to export as FASTA. All clusters are selected by default.",
            wraplength=520,
            foreground="#444444",
        ).pack(anchor="w", pady=(0, 8))

        tools = ttk.Frame(outer)
        tools.pack(fill="x", pady=(0, 6))
        selected_initial = set(available if self.extract_clusters_selection is None else [c for c in self.extract_clusters_selection if c in available])
        vars_by_name = {name: tk.BooleanVar(value=(name in selected_initial)) for name in available}

        def set_all(value):
            for v in vars_by_name.values():
                v.set(bool(value))

        ttk.Button(tools, text="Select all", command=lambda: set_all(True)).pack(side="left")
        ttk.Button(tools, text="Select none", command=lambda: set_all(False)).pack(side="left", padx=(6, 0))

        scroller = ScrollableFrame(outer)
        scroller.pack(fill="both", expand=True)
        for name in available:
            ttk.Checkbutton(scroller.inner, text=name, variable=vars_by_name[name]).pack(anchor="w", pady=2)

        bottom = ttk.Frame(outer)
        bottom.pack(fill="x", pady=(8, 0))

        def apply_selection():
            chosen = [name for name in available if vars_by_name[name].get()]
            self.extract_clusters_selection = chosen
            self._update_extract_clusters_status(available)
            win.destroy()

        ttk.Button(bottom, text="Use selected clusters", command=apply_selection).pack(side="left")
        ttk.Button(bottom, text="Cancel", command=win.destroy).pack(side="left", padx=(8, 0))

    def _update_extract_clusters_status(self, available=None):
        if available is None:
            try:
                available = self._get_available_cluster_names(show_errors=False)
            except Exception:
                available = []
        if self.extract_clusters_selection is None:
            n = len(available or [])
            self.extract_clusters_status.set("All available clusters selected" if n == 0 else f"All available clusters selected ({n})")
        else:
            n = len([c for c in self.extract_clusters_selection if (not available or c in available)])
            self.extract_clusters_status.set(f"{n} cluster(s) selected for FASTA extraction")

    def _current_genome_fasta_path(self):
        if not self.genome_not_available.get():
            self._apply_preloaded_genome_choice()
        path = self.fasta_path.get().strip()
        if not path or not os.path.isfile(path):
            raise FileNotFoundError(f"Selected genome FASTA not found:\n{path}")
        return path

    def _current_sequence_records_for_export(self):
        """Return {id: (header_without_gt, clean_sequence)} for selected genome plus optional custom CDS."""
        base_fasta = self._current_genome_fasta_path()
        records = {}
        for header, seq in _iter_fasta_records_for_compare(base_fasta):
            rec_id = _extract_fasta_record_id(header)
            if not rec_id:
                continue
            clean = _clean_dna_for_compare(seq)
            if not clean:
                continue
            records.setdefault(rec_id, (header, clean))

        if bool(self.add_custom_cds_enable.get()):
            counter = 0
            for path in _split_custom_cds_paths_gui(self.custom_cds_paths.get()):
                if not os.path.isfile(path):
                    continue
                for header, seq in _iter_fasta_records_for_compare(path):
                    clean = _clean_dna_for_compare(seq)
                    if not clean:
                        continue
                    counter += 1
                    custom_id = f"custom_CDS_{counter:03d}"
                    product = f"Custom CDS from {os.path.basename(path)}; original header: {header}".replace("[", "(").replace("]", ")")
                    custom_header = f"{custom_id} [locus_tag={custom_id}] [gene={custom_id}] [product={product}]"
                    records[custom_id] = (custom_header, clean)
        return records

    def _current_fasta_annotation_df(self):
        rows = []
        for rec_id, (header, _seq) in self._current_sequence_records_for_export().items():
            rows.append({
                "PrimaryID": str(header).split()[0] if header else rec_id,
                "LocusTag": rec_id,
                "RefSeq_LocusTag_RS": rec_id,
                "GeneSymbol": _extract_gene_symbol_from_header(header),
                "ProteinDescription": _extract_product_from_header(header),
                "Header": header,
            })
        return pd.DataFrame(rows)

    def _build_current_cluster_df_for_tools(self):
        """Build the currently available cluster table without running the full pipeline."""
        mode = self.user_cluster_mode.get().strip()
        internal = _cluster_mode_internal(mode)
        if internal == "refined":
            path = self.refined_cluster_file.get().strip()
            if not path or not os.path.isfile(path):
                raise FileNotFoundError("Please select the user cluster xlsx file first.")
            cluster_df = _read_cluster_file(path, sheet_name=self.refined_cluster_sheet.get().strip())
        else:
            keyword_groups = self._collect_keyword_groups()
            if not keyword_groups:
                cluster_df = pd.DataFrame()
                source = "geneids"
            else:
                source = _cluster_mode_terms_source(mode)
            if keyword_groups and source == "david_gene2terms":
                david_path = self.david_gene2terms_path.get().strip()
                if not david_path or not os.path.isfile(david_path):
                    raise FileNotFoundError("Please select a DAVID gene2terms TXT file first.")
                annotation_df = CP._read_david_gene2terms_txt(david_path)
                cluster_df, _long_df, _summary_df = CP.build_basic_clusters_from_annotation_df(
                    annotation_df=annotation_df,
                    keyword_groups=keyword_groups,
                    output_id_column="DisplayLocusTag",
                    output_id_fallbacks=["GeneID", "EntrezGeneID", "DisplayLocusTag"],
                    search_columns_preferred=["DAVID_Terms", "CleanTerm", "Term"],
                )
            elif keyword_groups:
                annotation_df = self._current_fasta_annotation_df()
                cluster_df, _long_df, _summary_df = CP.build_basic_clusters_from_annotation_df(
                    annotation_df=annotation_df,
                    keyword_groups=keyword_groups,
                    output_id_column="RefSeq_LocusTag_RS",
                    output_id_fallbacks=["LocusTag", "PrimaryID", "GeneSymbol"],
                    search_columns_preferred=["GeneSymbol", "ProteinDescription", "Header"],
                )

        if bool(self.add_custom_cds_enable.get()):
            custom_name = self.custom_cds_cluster_name.get().strip() or "custom"
            custom_ids = [k for k in self._current_sequence_records_for_export().keys() if k.startswith("custom_CDS_")]
            if custom_ids:
                df = cluster_df.copy().reset_index(drop=True)
                col = custom_name
                existing_cols_lower = {str(c).lower(): c for c in df.columns}
                if col.lower() in existing_cols_lower:
                    col = existing_cols_lower[col.lower()]
                max_len = max(len(df), len(custom_ids))
                df = df.reindex(range(max_len)).fillna("")
                existing = [str(v).strip() for v in df[col].tolist()] if col in df.columns else []
                values = list(dict.fromkeys([v for v in existing if v] + custom_ids))
                df[col] = pd.Series(values + [""] * (max_len - len(values)), index=range(max_len))
                cluster_df = df

        metric_cluster_df = self._build_fasta_metric_cluster_df_for_tools()
        if metric_cluster_df is not None and not metric_cluster_df.empty:
            cluster_df = append_fasta_metric_clusters(cluster_df, metric_cluster_df)
        return cluster_df.fillna("")

    def _selected_extract_cluster_names(self, cluster_df):
        available = [str(c).strip() for c in cluster_df.columns if str(c).strip()]
        if self.extract_clusters_selection is None:
            return available
        wanted = {str(c).strip().lower() for c in self.extract_clusters_selection}
        return [c for c in available if c.lower() in wanted]

    def _extract_fasta_from_clusters(self):
        try:
            cluster_df = self._build_current_cluster_df_for_tools()
            selected_clusters = self._selected_extract_cluster_names(cluster_df)
            if not selected_clusters:
                raise ValueError("No cluster selected for FASTA extraction.")
            seq_records = self._current_sequence_records_for_export()
            if not seq_records:
                raise ValueError("No CDS sequences were found in the selected genome FASTA.")

            members = []
            for cname in selected_clusters:
                vals = cluster_df[cname].fillna("").astype(str).str.strip().tolist()
                vals = [v for v in vals if v and v.lower() != "nan"]
                members.extend(vals)
            members = list(dict.fromkeys(members))
            found = [m for m in members if m in seq_records]
            missing = [m for m in members if m not in seq_records]
            if not found:
                raise ValueError("The selected cluster(s) contained no IDs matching the current genome FASTA records.")

            root = self.default_root.get().strip() or os.path.dirname(self._current_genome_fasta_path())
            root = os.path.abspath(os.path.expanduser(root))
            short_label = " + ".join(selected_clusters)
            folder_name = _safe_output_name(short_label, fallback="clusters", max_len=30)
            if len(_safe_output_name(short_label, fallback="clusters", max_len=300)) > 30:
                folder_name = _safe_output_name(short_label, fallback="clusters", max_len=26).rstrip("._- ") + " etc"
            out_dir = os.path.join(root, folder_name)
            os.makedirs(out_dir, exist_ok=True)

            # One FASTA per gene/CDS.
            used_files = set()
            for rec_id in found:
                header, seq = seq_records[rec_id]
                base = _safe_output_name(rec_id, fallback="CDS", max_len=90)
                fname = base + ".fna"
                i = 2
                while fname.lower() in used_files or os.path.exists(os.path.join(out_dir, fname)):
                    fname = f"{base}_{i}.fna"
                    i += 1
                used_files.add(fname.lower())
                with open(os.path.join(out_dir, fname), "w", encoding="utf-8") as fh:
                    fh.write(f">{header}\n{_wrap_fasta(seq)}\n")

            combined_name = "Multiple clusters" if len(selected_clusters) > 1 else selected_clusters[0]
            combined_name = _safe_output_name(combined_name, fallback="cluster", max_len=90) + ".fna"
            combined_path = os.path.join(out_dir, combined_name)
            with open(combined_path, "w", encoding="utf-8") as fh:
                for rec_id in found:
                    header, seq = seq_records[rec_id]
                    fh.write(f">{header}\n{_wrap_fasta(seq)}\n")

            msg = f"Exported {len(found)} CDS from {len(selected_clusters)} cluster(s) to:\n{out_dir}"
            if missing:
                msg += f"\n\nWarning: {len(missing)} cluster ID(s) were not found in the current FASTA and were skipped."
            self._append_log("\n[GUI] " + msg.replace("\n", "\n[GUI] ") + "\n")
            messagebox.showinfo("Extract FASTA", msg)
        except Exception as e:
            messagebox.showerror("Extract FASTA failed", str(e))

    def _flush_log_queue_now(self):
        try:
            while True:
                msg = self.log_queue.get_nowait()
                self._append_log(msg)
        except queue.Empty:
            pass

    # ------------------------- active-cluster selector -------------------------
    def _open_active_clusters_dialog(self):
        available = self._get_available_cluster_names(show_errors=True)
        if not available:
            messagebox.showinfo("Active clusters", "No clusters are currently available. Check the selected cluster source and keyword groups.")
            return

        win = tk.Toplevel(self.root)
        win.title("Active clusters")
        win.geometry("560x640")
        win.transient(self.root)
        win.grab_set()

        outer = ttk.Frame(win, padding=10)
        outer.pack(fill="both", expand=True)
        ttk.Label(
            outer,
            text="Select the clusters to use for cluster plots and 2D KS analyses. All clusters are selected by default.",
            wraplength=520,
            foreground="#444444",
        ).pack(anchor="w", pady=(0, 8))

        tools = ttk.Frame(outer)
        tools.pack(fill="x", pady=(0, 6))

        selected_initial = set(self._active_cluster_selection_for_available(available))
        vars_by_name = {name: tk.BooleanVar(value=(name in selected_initial)) for name in available}

        def set_all(value):
            for v in vars_by_name.values():
                v.set(bool(value))

        ttk.Button(tools, text="Select all", command=lambda: set_all(True)).pack(side="left", padx=(0, 6))
        ttk.Button(tools, text="Select none", command=lambda: set_all(False)).pack(side="left")

        sf = ScrollableFrame(outer)
        sf.pack(fill="both", expand=True, pady=(4, 8))
        for i, name in enumerate(available):
            chk = ttk.Checkbutton(sf.inner, text=name, variable=vars_by_name[name])
            chk.grid(row=i, column=0, sticky="w", padx=4, pady=2)

        buttons = ttk.Frame(outer)
        buttons.pack(fill="x")

        def apply_and_close():
            chosen = [name for name in available if vars_by_name[name].get()]
            if not chosen:
                if not messagebox.askyesno("Active clusters", "No cluster is selected. Continue with no active cluster-specific plots/tests?"):
                    return
            # Store an explicit selection. If all are checked, reset to None so newly added
            # clusters are automatically included later.
            self.active_clusters_selection = None if len(chosen) == len(available) else chosen
            self._update_active_clusters_status()
            win.destroy()

        def cancel():
            win.destroy()

        ttk.Button(buttons, text="Apply", command=apply_and_close).pack(side="right", padx=(6, 0))
        ttk.Button(buttons, text="Cancel", command=cancel).pack(side="right")

    def _open_figure_clusters_dialog(self):
        available = self._get_available_figure_cluster_names(show_errors=True)
        if not available:
            messagebox.showinfo("Figure cluster picker", "No clusters are currently available. Check the selected cluster source and keyword groups.")
            return

        win = tk.Toplevel(self.root)
        win.title("Figure cluster picker")
        win.geometry("560x640")
        win.transient(self.root)
        win.grab_set()

        outer = ttk.Frame(win, padding=10)
        outer.pack(fill="both", expand=True)
        ttk.Label(
            outer,
            text=(
                "Select the clusters to plot from the Figures tab. This affects cluster-based figures only: "
                "clusters along the genome axis, 2D cluster maps, and per-cluster codon-usage profiles. "
                "It does not change the Active clusters used in Input/Output or the Decoding strategies cluster picker."
            ),
            wraplength=520,
            foreground="#444444",
        ).pack(anchor="w", pady=(0, 8))

        tools = ttk.Frame(outer)
        tools.pack(fill="x", pady=(0, 6))

        selected_initial = set(self._figure_cluster_selection_for_available(available))
        vars_by_name = {name: tk.BooleanVar(value=(name in selected_initial)) for name in available}

        def set_all(value):
            for v in vars_by_name.values():
                v.set(bool(value))

        ttk.Button(tools, text="Select all", command=lambda: set_all(True)).pack(side="left", padx=(0, 6))
        ttk.Button(tools, text="Select none", command=lambda: set_all(False)).pack(side="left")

        sf = ScrollableFrame(outer)
        sf.pack(fill="both", expand=True, pady=(4, 8))
        for i, name in enumerate(available):
            chk = ttk.Checkbutton(sf.inner, text=name, variable=vars_by_name[name])
            chk.grid(row=i, column=0, sticky="w", padx=4, pady=2)

        buttons = ttk.Frame(outer)
        buttons.pack(fill="x")

        def apply_and_close():
            chosen = [name for name in available if vars_by_name[name].get()]
            if not chosen:
                if not messagebox.askyesno("Figure cluster picker", "No cluster is selected. Continue with no figure cluster-specific plots?"):
                    return
            self.figure_clusters_selection = None if len(chosen) == len(available) else chosen
            self._update_figure_clusters_status()
            win.destroy()

        def cancel():
            win.destroy()

        ttk.Button(buttons, text="Apply", command=apply_and_close).pack(side="right", padx=(6, 0))
        ttk.Button(buttons, text="Cancel", command=cancel).pack(side="right")

    def _open_decoding_cluster_picker(self):
        """Backward-compatible alias used by the Decoding strategies tab button."""
        return self._open_decoding_clusters_dialog()

    def _open_decoding_clusters_dialog(self):
        available = self._get_available_cluster_names(show_errors=True)
        if not available:
            messagebox.showinfo("Decoding cluster picker", "No clusters are currently available. Check the selected cluster source and keyword groups.")
            return

        win = tk.Toplevel(self.root)
        win.title("Decoding cluster picker")
        win.geometry("560x640")
        win.transient(self.root)
        win.grab_set()

        outer = ttk.Frame(win, padding=10)
        outer.pack(fill="both", expand=True)
        ttk.Label(
            outer,
            text="Select the clusters to use for decoding-strategy analyses only. This does not change the Active clusters used by the main clustering figures.",
            wraplength=520,
            foreground="#444444",
        ).pack(anchor="w", pady=(0, 8))

        tools = ttk.Frame(outer)
        tools.pack(fill="x", pady=(0, 6))

        selected_initial = set(self._decoding_cluster_selection_for_available(available))
        vars_by_name = {name: tk.BooleanVar(value=(name in selected_initial)) for name in available}

        def set_all(value):
            for v in vars_by_name.values():
                v.set(bool(value))

        ttk.Button(tools, text="Select all", command=lambda: set_all(True)).pack(side="left", padx=(0, 6))
        ttk.Button(tools, text="Select none", command=lambda: set_all(False)).pack(side="left")

        sf = ScrollableFrame(outer)
        sf.pack(fill="both", expand=True, pady=(4, 8))
        for i, name in enumerate(available):
            chk = ttk.Checkbutton(sf.inner, text=name, variable=vars_by_name[name])
            chk.grid(row=i, column=0, sticky="w", padx=4, pady=2)

        buttons = ttk.Frame(outer)
        buttons.pack(fill="x")

        def apply_and_close():
            chosen = [name for name in available if vars_by_name[name].get()]
            if not chosen:
                if not messagebox.askyesno("Decoding cluster picker", "No cluster is selected. Continue with no decoding cluster-specific plots?"):
                    return
            self.decoding_clusters_selection = None if len(chosen) == len(available) else chosen
            self.decoding_clusters_selection_user_set = True
            if (not str(self.decoding_reference_cluster.get() or "").strip()) and any(str(c).strip().lower() == "ribosomal proteins" for c in chosen):
                self.decoding_reference_cluster.set(next(str(c) for c in chosen if str(c).strip().lower() == "ribosomal proteins"))
            self._update_decoding_clusters_status()
            win.destroy()

        def cancel():
            win.destroy()

        ttk.Button(buttons, text="Apply", command=apply_and_close).pack(side="right", padx=(6, 0))
        ttk.Button(buttons, text="Cancel", command=cancel).pack(side="right")


    def _open_decoding_reference_cluster_dialog(self):
        available = self._get_available_cluster_names(show_errors=True)
        if not available:
            messagebox.showinfo("Reference cluster", "No clusters are currently available. Check the selected cluster source and keyword groups.")
            return
        self._ensure_default_decoding_preselection_for_available(available)
        selected = self._decoding_cluster_selection_for_available(available)
        choices = selected if selected else available
        win = tk.Toplevel(self.root)
        win.title("Reference cluster for statistical analyses")
        win.geometry("520x220")
        win.transient(self.root)
        win.grab_set()
        outer = ttk.Frame(win, padding=12)
        outer.pack(fill="both", expand=True)
        ttk.Label(
            outer,
            text="Choose the reference cluster used for statistical comparisons in decoding boxplots/violin plots. Comparisons are made against this cluster independently for each x-axis feature.",
            wraplength=480,
            foreground="#444444",
        ).pack(anchor="w", pady=(0, 10))
        current = str(self.decoding_reference_cluster.get() or "").strip()
        if current not in choices and choices:
            current = choices[0]
        ref_var = tk.StringVar(value=current)
        row = ttk.Frame(outer)
        row.pack(fill="x", pady=4)
        ttk.Label(row, text="Reference cluster").pack(side="left", padx=(0, 8))
        combo = ttk.Combobox(row, textvariable=ref_var, state="readonly", values=choices, width=42)
        combo.pack(side="left", fill="x", expand=True)
        buttons = ttk.Frame(outer)
        buttons.pack(fill="x", pady=(14, 0))
        def apply_and_close():
            self.decoding_reference_cluster.set(str(ref_var.get() or "").strip())
            self._update_decoding_clusters_status()
            win.destroy()
        def clear_and_close():
            self.decoding_reference_cluster.set("")
            self._update_decoding_clusters_status()
            win.destroy()
        ttk.Button(buttons, text="Apply", command=apply_and_close).pack(side="right", padx=(6, 0))
        ttk.Button(buttons, text="Clear", command=clear_and_close).pack(side="right", padx=(6, 0))
        ttk.Button(buttons, text="Cancel", command=win.destroy).pack(side="right")

    def _current_decoding_table_path_for_selection(self):
        trna = self.trna_decoding_table_path.get().strip()
        if not trna:
            try:
                trna = self._auto_find_decoding_table_for_current_genome()
            except Exception:
                trna = ""
            if trna:
                self.trna_decoding_table_path.set(trna)
        return trna

    def _read_decoding_rules_for_selection(self):
        trna = self._current_decoding_table_path_for_selection()
        if not trna or not os.path.isfile(trna):
            messagebox.showinfo(
                "Decoding table required",
                "Please select a decoding table in the Input/Output tab first. The selector uses this workbook to list available modifications and amino acids."
            )
            return None
        try:
            reader = _read_trna_decoding_table_for_gui or getattr(CP, "read_trna_decoding_table", None)
            if reader is None:
                raise AttributeError(
                    "read_trna_decoding_table is not available. "
                    "Please check that codonpipe/excel_outputs.py is present next to the GUI."
                )
            return reader(trna, sheet_name=self.trna_decoding_table_sheet.get().strip())
        except Exception as e:
            messagebox.showwarning("Decoding selector", f"Could not read the decoding table:\n{e}")
            return None

    def _available_trna_modification_features(self):
        """Return selectable Plot 6 modification/enzyme features.

        The selector follows the same rule as Plot 6 itself: feature names are
        taken from the full decoder-level table when available. The compact
        pooled table may still define tRNA usage/ZTU, but it should not define
        the Plot 6 modification/enzyme checklist.
        """
        rules = self._read_decoding_rules_for_selection()
        if not rules:
            return []
        mode = str(self.trna_modifications_feature_mode.get() or "modifications").strip().lower()
        use_enzymes = mode in {"enzyme", "enzymes", "trme", "trmes"}
        map_key = "trmes_by_codon" if use_enzymes else "modifications_by_codon"
        candidate_col = "tRMEs" if use_enzymes else "tRNA_modifications"
        features = []
        seen = set()

        def add_feature(x):
            val = " ".join(str(x or "").strip().split())
            if not val or val.lower() in {"nan", "none", "na", "n/a"}:
                return
            if not use_enzymes:
                val = _canonicalize_gui_plot6_modification_feature(val)
                if not val:
                    return
            key = _plot6_feature_key_ascii_gui(val) if not use_enzymes else val.lower()
            if key not in seen:
                seen.add(key)
                features.append(val)

        def add_feature_text(value):
            text = str(value or "")
            for item in re.split(r"[,;\n\r]+", text):
                add_feature(item)

        # Primary source: candidate table generated from "Decoding table (full)".
        candidate_df = rules.get("modification_candidate_table")
        if candidate_df is not None and not getattr(candidate_df, "empty", True) and candidate_col in candidate_df.columns:
            for value in candidate_df[candidate_col].astype(str).tolist():
                add_feature_text(value)

        # Secondary source: the actual Plot 6 model maps, also generated from
        # the full decoder-level table in the patched backend.
        def collect_from_map(payload):
            for _codon, items in (payload or {}).items():
                iterable = items.keys() if isinstance(items, dict) else list(items or [])
                for item in iterable:
                    add_feature(item)

        models = rules.get("modification_assignment_models") or {}
        for model_name in ["permissive", "conservative"]:
            collect_from_map((models.get(model_name) or {}).get(map_key) or {})

        if use_enzymes:
            return sorted(features, key=lambda x: x.lower())
        return sorted(features, key=_gui_plot6_modification_sort_key)

    def _available_trna_modification_aas(self):
        """Return selectable amino-acid families for Plot 6.

        Prefer the full-table Plot 6 candidate table over the compact decoding
        table, so this button remains synchronized with the modification/enzyme
        plot even when the pooled compact table is selected for tRNA usage.
        """
        rules = self._read_decoding_rules_for_selection()
        if not rules:
            return []
        seen = set()
        present = []

        def add_aa(value):
            aa = str(value or "").strip()
            if not aa or aa.lower() in {"nan", "none", "na", "n/a"}:
                return
            key = aa.lower()
            if key not in seen:
                seen.add(key)
                present.append(aa)

        candidate_df = rules.get("modification_candidate_table")
        if candidate_df is not None and not getattr(candidate_df, "empty", True) and "AA" in candidate_df.columns:
            for aa in candidate_df["AA"].astype(str).tolist():
                add_aa(aa)

        # Fallbacks for legacy workbooks without a full Plot 6 candidate table.
        if not present:
            try:
                meta = rules.get("table_df")
                if meta is not None and not getattr(meta, "empty", True) and "AA" in meta.columns:
                    for aa in meta["AA"].astype(str).tolist():
                        add_aa(aa)
            except Exception:
                pass
        if not present:
            try:
                for codon in (rules.get("codon_to_decoders") or {}).keys():
                    aa = str(codon).split("_", 1)[0].strip() if "_" in str(codon) else ""
                    add_aa(aa)
            except Exception:
                pass

        canonical = list(_AA_TO_CODONS.keys())
        ordered = [aa for aa in canonical if aa.lower() in seen]
        ordered.extend([aa for aa in present if aa not in ordered])
        return ordered

    def _open_generic_checklist_dialog(self, title, available, selected, allow_empty=True, intro_text=""):
        available = [str(x).strip() for x in list(available or []) if str(x).strip()]
        if not available:
            messagebox.showinfo(title, "No selectable item was found in the current decoding table.")
            return None
        selected_set = {str(x).strip().lower() for x in list(selected or [])} if selected is not None else {x.lower() for x in available}
        win = tk.Toplevel(self.root)
        win.title(title)
        win.geometry("560x680")
        win.transient(self.root)
        win.grab_set()
        outer = ttk.Frame(win, padding=10)
        outer.pack(fill="both", expand=True)
        ttk.Label(outer, text=intro_text or "Select items to include.", wraplength=520, foreground="#444444").pack(anchor="w", pady=(0, 8))
        tools = ttk.Frame(outer)
        tools.pack(fill="x", pady=(0, 6))
        vars_by_item = {item: tk.BooleanVar(value=(item.lower() in selected_set)) for item in available}
        def set_all(value):
            for v in vars_by_item.values():
                v.set(bool(value))
        ttk.Button(tools, text="Select all", command=lambda: set_all(True)).pack(side="left", padx=(0, 6))
        ttk.Button(tools, text="Select none", command=lambda: set_all(False)).pack(side="left")
        sf = ScrollableFrame(outer)
        sf.pack(fill="both", expand=True, pady=(4, 8))
        for i, item in enumerate(available):
            ttk.Checkbutton(sf.inner, text=item, variable=vars_by_item[item]).grid(row=i, column=0, sticky="w", padx=4, pady=2)
        result = {"selected": None, "applied": False}
        buttons = ttk.Frame(outer)
        buttons.pack(fill="x")
        def apply_and_close():
            chosen = [item for item in available if vars_by_item[item].get()]
            if not chosen and not allow_empty:
                messagebox.showwarning(title, "Please select at least one item.")
                return
            result["selected"] = None if len(chosen) == len(available) else chosen
            result["applied"] = True
            win.destroy()
        ttk.Button(buttons, text="Apply", command=apply_and_close).pack(side="right", padx=(6, 0))
        ttk.Button(buttons, text="Cancel", command=win.destroy).pack(side="right")
        self.root.wait_window(win)
        return result["selected"] if result.get("applied") else "__CANCELLED__"

    def _open_trna_modification_feature_dialog(self):
        mode = str(self.trna_modifications_feature_mode.get() or "modifications").strip().lower()
        label = "enzymes" if mode in {"enzyme", "enzymes", "trme", "trmes"} else "modifications"
        available = self._available_trna_modification_features()
        selected_for_dialog = getattr(self, "trna_modifications_selection", None)
        if label == "modifications" and selected_for_dialog is None:
            # Default state: keep minor/specific modifications visible but unchecked.
            selected_for_dialog = [f for f in available if not _is_default_excluded_gui_plot6_modification(f)]
        result = self._open_generic_checklist_dialog(
            f"Select tRNA {label} for Plot 6",
            available,
            selected_for_dialog,
            allow_empty=True,
            intro_text=(
                f"Select which tRNA {label} should be plotted in Plot 6. The list is read from the selected decoding Excel workbook. "
                "By default, ac4C34 is unchecked; m6A37 is selected. Select all includes ac4C34 again. Select none disables Plot 6 feature output until you select at least one item again."
            ),
        )
        if result == "__CANCELLED__":
            return
        if label == "modifications":
            # In the generic dialog, None means every available item was selected.
            # Store this explicitly so a user click on 'Select all' overrides the default exclusions.
            if result is None:
                result = list(available)
            else:
                result = _canonicalize_gui_plot6_modification_selection(result)
        self.trna_modifications_selection = result
        self._update_trna_modification_selection_status()

    def _open_trna_modification_aa_dialog(self):
        available = self._available_trna_modification_aas()
        selected_for_dialog = getattr(self, "trna_modification_aas_selection", None)
        if selected_for_dialog is None:
            # Default state: select the manuscript amino-acid set and leave other amino acids visible but unchecked.
            selected_for_dialog = _default_gui_plot6_included_aas(available)
        result = self._open_generic_checklist_dialog(
            "Amino acids considered for Plot 6",
            available,
            selected_for_dialog,
            allow_empty=False,
            intro_text=(
                "Select which amino-acid families are considered when calculating Plot 6 modification/enzyme enrichment. "
                "Codons from unselected amino acids are ignored and the percentage denominator is restricted to the selected amino-acid families. "
                "By default, Ala, Arg, Asn, Asp, Cys, Gly, His, Ile, Leu, Phe, Pro, Ser, Thr and Tyr are selected. Select all includes every available amino acid."
            ),
        )
        if result == "__CANCELLED__":
            return
        # In the generic dialog, None means every available amino acid was selected.
        # Store this explicitly so a user click on 'Select all' overrides the default exclusions.
        self.trna_modification_aas_selection = list(available) if result is None else result
        self._update_trna_modification_aas_status()


    # ------------------------- monkey patching -------------------------
    def _cluster_sizes_from_cluster_df(self, cluster_df):
        """Return unique, non-empty locus-tag counts for each cluster column."""
        sizes = {}
        try:
            for col in list(cluster_df.columns):
                vals = []
                for v in cluster_df[col].replace({np.nan: ""}).astype(str).tolist():
                    vv = str(v).strip()
                    if vv and vv.lower() != "nan":
                        vals.append(vv)
                sizes[str(col)] = len(dict.fromkeys(vals))
        except Exception:
            sizes = {}
        return sizes

    def _filter_clusters_by_min_genes(self, cluster_names, cluster_sizes, min_genes, fallback=None):
        """Apply the GUI min-gene threshold to cluster names, with a safe fallback."""
        names = [str(c).strip() for c in list(cluster_names or []) if str(c).strip()]
        try:
            min_n = max(0, int(float(str(min_genes).strip())))
        except Exception:
            min_n = 0
        filtered = [c for c in names if int(cluster_sizes.get(c, 0) or 0) >= min_n]
        if filtered:
            return list(dict.fromkeys(filtered))
        fb = [str(c).strip() for c in list(fallback or names) if str(c).strip()]
        return list(dict.fromkeys(fb))

    def _make_choose_clusters_func(self, active_selection, min_genes):
        def _patched_choose_clusters(cluster_df, _min_genes_default):
            ordered = [str(col) for col in list(cluster_df.columns)]
            sizes = self._cluster_sizes_from_cluster_df(cluster_df)
            lower_map = {str(col).strip().lower(): str(col) for col in ordered}

            # None means the default GUI state: all currently available clusters
            # are active, but the min-gene threshold must still be applied.
            if active_selection is None:
                return self._filter_clusters_by_min_genes(ordered, sizes, min_genes, fallback=ordered)

            requested = [str(x).strip() for x in list(active_selection or []) if str(x).strip()]
            chosen = []
            for name in requested:
                mapped = lower_map.get(name.lower())
                if mapped is not None:
                    chosen.append(mapped)
            chosen = list(dict.fromkeys(chosen))
            return self._filter_clusters_by_min_genes(chosen, sizes, min_genes, fallback=chosen or ordered)
        return _patched_choose_clusters

    def _make_load_refined_func(self, cluster_path, sheet_name):
        def _patched_load_refined_cluster_df(_base_folder):
            cluster_df = _read_cluster_file(cluster_path, sheet_name=sheet_name)
            return cluster_df, cluster_path
        return _patched_load_refined_cluster_df

    def _apply_run_mode_overrides(self, run_mode, runtime_choices, set_overrides, pipeline_overrides):
        """Map GUI buttons to the subset of analyses/figures that should run."""
        mode = str(run_mode or "clustering_plots")
        runtime_choices["run_david_scan"] = False

        def disable_all_standard_figures():
            set_overrides["plot_codon_gene_heatmap_enable"] = False
            set_overrides["show_main_pipeline_figures"] = False
            runtime_choices["run_2d_density_plots"] = False
            runtime_choices["run_gchm"] = False
            pipeline_overrides["auto_run_plotting_pipeline"] = False
            pipeline_overrides["RUN_2D_DENSITY_PLOTS"] = False
            pipeline_overrides["RUN_GENE_CLUSTER_LOCALIZATION_HEATMAP"] = False
            pipeline_overrides["CODON_USAGE_PLOT_MODE"] = "NONE"
            pipeline_overrides["codon_usage_plot_mode"] = "NONE"

        # Core analysis buttons should never run decoding strategy analyses.
        if mode in {"clustering_only", "clustering_plots", "figure_main_heatmap", "figure_cluster_axis", "figure_2d_map", "figure_codon_profiles"}:
            runtime_choices["compute_trna_usage"] = False
            runtime_choices["compute_trna_abundance_correlations"] = False
            set_overrides["export_trna_usage_enable"] = False
            set_overrides["export_trna_abundance_correlation_enable"] = False

        if mode == "clustering_only":
            disable_all_standard_figures()
            return

        if mode == "clustering_plots":
            # "Run clustering followed by plotting" runs the full clustering
            # workflow and then plots the standard figure set. Cluster
            # subsetting for this run is controlled by the Input/Output tab's
            # "Active clusters" selector. The Figures-tab "Cluster picker" is
            # only used when replotting individual figures from an existing
            # workbook.
            set_overrides["plot_codon_gene_heatmap_enable"] = True
            set_overrides["show_main_pipeline_figures"] = True
            runtime_choices["run_2d_density_plots"] = True
            runtime_choices["run_gchm"] = True
            pipeline_overrides["auto_run_plotting_pipeline"] = True
            pipeline_overrides["RUN_2D_DENSITY_PLOTS"] = True
            pipeline_overrides["RUN_GENE_CLUSTER_LOCALIZATION_HEATMAP"] = True
            mode_key = USER_CODON_MODE_TO_INTERNAL.get(self.codon_usage_plot_mode.get().strip(), "Z")
            pipeline_overrides["CODON_USAGE_PLOT_MODE"] = mode_key
            pipeline_overrides["codon_usage_plot_mode"] = mode_key
            return

        if mode == "figure_main_heatmap":
            disable_all_standard_figures()
            set_overrides["plot_codon_gene_heatmap_enable"] = True
            set_overrides["show_main_pipeline_figures"] = True
            return

        if mode == "figure_cluster_axis":
            disable_all_standard_figures()
            runtime_choices["run_gchm"] = True
            pipeline_overrides["auto_run_plotting_pipeline"] = True
            pipeline_overrides["RUN_GENE_CLUSTER_LOCALIZATION_HEATMAP"] = True
            return

        if mode == "figure_2d_map":
            disable_all_standard_figures()
            runtime_choices["run_2d_density_plots"] = True
            pipeline_overrides["auto_run_plotting_pipeline"] = True
            pipeline_overrides["RUN_2D_DENSITY_PLOTS"] = True
            return

        if mode == "figure_codon_profiles":
            disable_all_standard_figures()
            pipeline_overrides["auto_run_plotting_pipeline"] = True
            mode_key = USER_CODON_MODE_TO_INTERNAL.get(self.codon_usage_plot_mode.get().strip(), "Z")
            pipeline_overrides["CODON_USAGE_PLOT_MODE"] = mode_key
            pipeline_overrides["codon_usage_plot_mode"] = mode_key
            return

        if mode.startswith("trna_"):
            disable_all_standard_figures()
            runtime_choices["compute_trna_usage"] = True
            runtime_choices["compute_trna_abundance_correlations"] = False
            set_overrides["export_trna_usage_enable"] = True
            set_overrides["export_trna_abundance_correlation_enable"] = False
            # Disable all tRNA figure toggles first, then re-enable the requested one(s).
            for key in [
                "trna_gene_heatmap_enable",
                "trna_single_box_codon_heatmap_enable",
                "trna_shift_heatmap_enable",
                "trna_wobble_heatmap_enable",
                "trna_modification_heatmap_enable",
            ]:
                set_overrides[key] = False
            if mode == "trna_all":
                set_overrides["trna_gene_heatmap_enable"] = True
                set_overrides["trna_single_box_codon_heatmap_enable"] = True
                set_overrides["trna_shift_heatmap_enable"] = True
                set_overrides["trna_wobble_heatmap_enable"] = True
                set_overrides["trna_modification_heatmap_enable"] = True
            elif mode == "trna_figure_gene_heatmap":
                set_overrides["trna_gene_heatmap_enable"] = True
            elif mode == "trna_figure_single_box":
                set_overrides["trna_single_box_codon_heatmap_enable"] = True
            elif mode == "trna_figure_shift":
                set_overrides["trna_shift_heatmap_enable"] = True
            elif mode == "trna_figure_wobble":
                set_overrides["trna_wobble_heatmap_enable"] = True
            elif mode == "trna_figure_modifications":
                set_overrides["trna_modification_heatmap_enable"] = True
            return

        # Default: clustering followed by selected plots, no tRNA.
        if mode == "clustering_plots":
            return


    # ------------------------- replot from existing clustering outputs -------------------------
    def _method_output_subfolder_name(self):
        method = _dimred_internal(self.dimred_method.get()).lower()
        if method == "tsne":
            per = _safe_float(self.dimred_param_vars.get("tsne", {}).get("tsne_perplexity", tk.StringVar(value="30")).get(), 30)
            ex = _safe_float(self.dimred_param_vars.get("tsne", {}).get("tsne_exaggeration", tk.StringVar(value="12")).get(), 12)
            lr = _safe_float(self.dimred_param_vars.get("tsne", {}).get("tsne_learnrate", tk.StringVar(value="200")).get(), 200)
            return f"tsne per{per:g} ex{ex:g} lr{lr:g}"
        if method == "umap":
            nn = _safe_int(self.dimred_param_vars.get("umap", {}).get("umap_neighbors", tk.StringVar(value="15")).get(), 15)
            md = _safe_float(self.dimred_param_vars.get("umap", {}).get("umap_min_dist", tk.StringVar(value="0.1")).get(), 0.1)
            md_str = f"{md:.3g}".replace(".", "p")
            return f"umap nn{nn} md{md_str}"
        if method == "pca":
            pcs = _safe_int(self.dimred_param_vars.get("pca", {}).get("pca_npcs", tk.StringVar(value="2")).get(), 2)
            return f"pca npc{pcs}"
        return "nodimred"

    def _candidate_output_dirs(self):
        dirs = []
        def add(d):
            d = str(d or "").strip()
            if not d:
                return
            d = os.path.abspath(os.path.expanduser(d))
            if d not in dirs:
                dirs.append(d)
        add(getattr(self, "last_clustering_output_dir", ""))
        add(self.default_root.get().strip())
        # Do not call _apply_preloaded_genome_choice() here: locating output
        # workbooks should not trigger genome/companion-file autoloading.
        fasta = self.fasta_path.get().strip()
        if fasta:
            base_dir = os.path.dirname(os.path.abspath(os.path.expanduser(fasta)))
            if self.default_root.get().strip():
                add(self.default_root.get().strip())
            else:
                add(os.path.join(base_dir, self._method_output_subfolder_name()))
                add(os.path.join(base_dir, "CodonPipe custom CDS merged FASTA", self._method_output_subfolder_name()))
        return dirs

    def _find_latest_file(self, candidates):
        existing = []
        for c in candidates:
            try:
                if c and os.path.isfile(c):
                    existing.append(os.path.abspath(c))
            except Exception:
                pass
        if not existing:
            return ""
        existing.sort(key=lambda p: os.path.getmtime(p), reverse=True)
        return existing[0]

    def _locate_clustering_workbook(self):
        candidates = []
        if getattr(self, "last_clustering_workbook", ""):
            candidates.append(self.last_clustering_workbook)
        for d in self._candidate_output_dirs():
            candidates.extend([
                os.path.join(d, "Clustering analysis results.xlsx"),
                os.path.join(d, f"{self.organism_name.get().strip() or 'Organism'}_ClusteringAnalysis.xlsx"),
            ])
            try:
                for fname in os.listdir(d):
                    if fname.lower().endswith("_clusteringanalysis.xlsx"):
                        candidates.append(os.path.join(d, fname))
            except Exception:
                pass
        return self._find_latest_file(candidates)

    def _set_last_clustering_outputs_from_disk(self):
        workbook = self._locate_clustering_workbook()
        if workbook:
            self.last_clustering_workbook = workbook
            self.last_clustering_output_dir = os.path.dirname(workbook)
        return workbook

    def _locate_codon_usage_workbook(self, output_dir):
        candidates = [
            os.path.join(output_dir, "Codon usage tables per cluster.xlsx"),
        ]
        try:
            for fname in os.listdir(output_dir):
                low = fname.lower()
                if low.endswith(".xlsx") and ("codon_usage" in low or "codon usage" in low) and "cluster" in low:
                    candidates.append(os.path.join(output_dir, fname))
        except Exception:
            pass
        return self._find_latest_file(candidates)

    def _active_clusters_for_workbook(self, workbook_path):
        try:
            df = pd.read_excel(workbook_path, sheet_name=CP.PIPELINE.get("sheet_locus_tags", "Locus Tags"), nrows=1, dtype=str)
            available = [str(c) for c in list(df.columns)[1:]]
        except Exception:
            available = []
        if not available:
            return []
        if self.active_clusters_selection is None:
            return available
        requested = [str(x).strip() for x in list(self.active_clusters_selection or []) if str(x).strip()]
        lower_map = {a.strip().lower(): a for a in available}
        chosen = [lower_map[r.lower()] for r in requested if r.lower() in lower_map]
        return list(dict.fromkeys(chosen)) or available

    def _figure_clusters_for_workbook(self, workbook_path):
        """Return figure-tab clusters after applying the min-gene threshold.

        This uses the actual cluster columns present in the already-generated
        clustering workbook, so PLOT buttons can replot a subset without rerunning
        clustering. The GUI min-gene threshold is applied to both Plot 2 and
        Plot 3 by filtering the cluster columns sent to the plotting pipeline.
        """
        sheet = CP.PIPELINE.get("sheet_locus_tags", "Locus Tags")
        try:
            df = pd.read_excel(workbook_path, sheet_name=sheet, dtype=str).fillna("")
            available = [str(c).strip() for c in list(df.columns)[1:] if str(c).strip()]
        except Exception:
            df = pd.DataFrame()
            available = []
        if not available:
            return []

        sizes = {}
        if not df.empty:
            for col in available:
                vals = [str(v).strip() for v in df[col].astype(str).tolist() if str(v).strip() and str(v).strip().lower() != "nan"]
                sizes[col] = len(dict.fromkeys(vals))

        if self.figure_clusters_selection is None:
            chosen = available
        else:
            requested = [str(x).strip() for x in list(self.figure_clusters_selection or []) if str(x).strip()]
            lower_map = {a.strip().lower(): a for a in available}
            chosen = [lower_map[r.lower()] for r in requested if r.lower() in lower_map]
            chosen = list(dict.fromkeys(chosen))

        min_genes = _safe_int(self.plot_cluster_min_genes.get(), 2)
        return self._filter_clusters_by_min_genes(chosen, sizes, min_genes, fallback=chosen or available)

    def _current_fasta_metric_cluster_df(self):
        """Return the current computed FASTA-derived metric cluster table, if any."""
        df = getattr(self, "fasta_metric_cluster_df", None)
        if isinstance(df, pd.DataFrame) and not df.empty:
            return df.fillna("").copy()
        path = str(getattr(self, "fasta_metric_cluster_path", "") or "").strip()
        if path and os.path.isfile(path):
            try:
                df = pd.read_excel(path, sheet_name="Clusters", dtype=str).fillna("")
                if not df.empty:
                    return df
            except Exception:
                try:
                    df = pd.read_excel(path, dtype=str).fillna("")
                    if not df.empty:
                        return df
                except Exception:
                    pass
        return pd.DataFrame()

    def _current_fasta_metric_cluster_members(self, ordered_tags=None):
        """Return {cluster_name: ordered locus tags} for computed metric clusters."""
        df = self._current_fasta_metric_cluster_df()
        if df is None or df.empty:
            return {}
        order = [str(x).strip() for x in list(ordered_tags or []) if str(x).strip()]
        order_index = {tag: i for i, tag in enumerate(order)}
        out = {}
        for col in df.columns:
            cname = str(col).strip()
            if not cname:
                continue
            vals = []
            seen = set()
            for v in df[col].astype(str).tolist():
                tag = str(v).strip()
                if not tag or tag.lower() == "nan" or tag in seen:
                    continue
                seen.add(tag)
                vals.append(tag)
            if order_index:
                vals = [v for v in vals if v in order_index]
                vals.sort(key=lambda t: order_index.get(t, 10**12))
            if vals:
                out[cname] = vals
        return out

    def _current_custom_cds_ids_for_gui(self):
        """Return stable custom_CDS_* IDs generated from the currently selected custom FASTA files."""
        if not bool(getattr(self, "add_custom_cds_enable", tk.BooleanVar(value=False)).get()):
            return []
        try:
            records = self._current_sequence_records_for_export()
        except Exception:
            return []
        ids, seen = [], set()
        for rec_id in records.keys():
            rec_id = str(rec_id).strip()
            if rec_id.startswith("custom_CDS_") and rec_id not in seen:
                seen.add(rec_id)
                ids.append(rec_id)
        return ids

    def _current_custom_cds_cluster_members(self, ordered_tags=None, warn_if_absent=False):
        """Return {cluster_name: custom_CDS IDs} for custom CDS already present in a workbook.

        This helper intentionally returns only custom CDS that are already present
        in the workbook order. For post-analysis custom CDS, the Figures-tab
        augmentation path uses _project_custom_cds_onto_existing_workbook() to add
        projected coordinates before plotting.
        """
        ids = self._current_custom_cds_ids_for_gui()
        if not ids:
            return {}
        cluster_name = str(getattr(self, "custom_cds_cluster_name", tk.StringVar(value="custom")).get() or "custom").strip() or "custom"
        order = [str(x).strip() for x in list(ordered_tags or []) if str(x).strip()]
        if order:
            order_index = {tag: i for i, tag in enumerate(order)}
            present = [x for x in ids if x in order_index]
            present.sort(key=lambda t: order_index.get(t, 10**12))
            if ids and not present and warn_if_absent:
                self._append_log(
                    "[INFO] The custom CDS cluster is not present in the current workbook; "
                    "CodonPipe will try to project the custom CDS onto the existing 2D map for Figures-tab plotting.\n"
                )
            ids = present
        return {cluster_name: ids} if ids else {}

    def _current_gui_extra_cluster_members(self, ordered_tags=None, warn_if_custom_absent=False):
        """Return GUI-generated cluster memberships that can augment an existing workbook.

        This includes FASTA-derived metric clusters and any custom CDS already
        present in the workbook. Post-analysis custom CDS are added separately by
        _augment_workbook_with_current_metric_clusters(), because they first need
        projected coordinates in the existing 2D map.
        """
        combined = {}

        def _merge(members_dict):
            for cname, vals in (members_dict or {}).items():
                cname = str(cname).strip()
                if not cname:
                    continue
                existing = combined.setdefault(cname, [])
                seen = set(existing)
                for v in vals:
                    tag = str(v).strip()
                    if tag and tag.lower() != "nan" and tag not in seen:
                        seen.add(tag)
                        existing.append(tag)

        _merge(self._current_custom_cds_cluster_members(ordered_tags=ordered_tags, warn_if_absent=warn_if_custom_absent))
        _merge(self._current_fasta_metric_cluster_members(ordered_tags=ordered_tags))

        order = [str(x).strip() for x in list(ordered_tags or []) if str(x).strip()]
        if order:
            order_index = {tag: i for i, tag in enumerate(order)}
            for cname in list(combined.keys()):
                vals = [v for v in combined[cname] if v in order_index]
                vals.sort(key=lambda t: order_index.get(t, 10**12))
                if vals:
                    combined[cname] = vals
                else:
                    combined.pop(cname, None)
        return combined

    def _custom_cds_count_df_for_projection(self):
        """Return per-custom-CDS codon counts with plain DNA codon columns."""
        records = self._current_sequence_records_for_export()
        rows = []
        ids = []
        for rec_id, (_header, seq) in records.items():
            rec_id = str(rec_id).strip()
            if not rec_id.startswith("custom_CDS_"):
                continue
            counts = dict.fromkeys(_CODON_ORDER, 0)
            clean = _clean_dna_for_compare(seq)
            for i in range(0, len(clean), 3):
                codon = clean[i:i + 3]
                if codon in counts:
                    counts[codon] += 1
            if sum(counts.values()) > 0:
                ids.append(rec_id)
                rows.append(counts)
        if not rows:
            return pd.DataFrame(columns=_CODON_ORDER, dtype=float)
        return pd.DataFrame(rows, index=ids, columns=_CODON_ORDER, dtype=float)

    def _codon_column_to_plain_dna(self, col):
        """Convert workbook codon labels such as Ala_GCA or Met_AUG to GCA/ATG."""
        s = str(col or "").strip()
        if "_" in s:
            s = s.split("_", 1)[1].strip()
        s = s.upper().replace("U", "T")
        if len(s) == 3 and all(ch in "ACGT" for ch in s):
            return s
        return ""

    def _load_genome_count_df_for_projection(self, workbook_path):
        """Load original genome codon counts for post-hoc custom-CDS projection."""
        output_dir = os.path.dirname(os.path.abspath(workbook_path)) if workbook_path else ""
        raw_wb = self._locate_whole_genome_raw_codon_usage_workbook(output_dir) if output_dir else ""
        errors = []
        if raw_wb and os.path.isfile(raw_wb):
            for header_row in (1, 0):
                try:
                    df = pd.read_excel(raw_wb, sheet_name="Codon counts", dtype=object, header=header_row)
                except Exception as e:
                    errors.append(str(e))
                    continue
                if df is None or df.empty or df.shape[1] < 10:
                    continue
                locus_col = df.columns[0]
                renamed = {}
                for c in list(df.columns)[1:]:
                    cod = self._codon_column_to_plain_dna(c)
                    if cod:
                        renamed[c] = cod
                if len(set(renamed.values())) < 20:
                    continue
                out = df[[locus_col] + list(renamed.keys())].copy()
                out[locus_col] = out[locus_col].astype(str).str.strip()
                out = out[out[locus_col].astype(str).str.lower() != "nan"]
                out = out[out[locus_col].astype(str).str.strip() != ""]
                out = out.rename(columns=renamed).set_index(locus_col)
                out = out.apply(pd.to_numeric, errors="coerce").fillna(0.0)
                out = out.T.groupby(level=0).sum().T
                for cod in _CODON_ORDER:
                    if cod not in out.columns:
                        out[cod] = 0.0
                out = out.loc[:, _CODON_ORDER]
                out.index = out.index.astype(str).str.strip()
                out = out[~out.index.duplicated(keep="first")]
                out = out[[not str(idx).startswith("custom_CDS_") for idx in out.index]]
                if not out.empty:
                    return out

        # Fallback: parse the selected genome FASTA. This is slower but keeps the
        # projection available for older workbooks that lack raw codon-count exports.
        try:
            base_fasta = self._current_genome_fasta_path()
            df = _codon_count_rows_from_fastas([base_fasta])
            if df is not None and not df.empty:
                df = df.loc[:, _CODON_ORDER]
                df.index = df.index.astype(str).str.strip()
                return df
        except Exception as e:
            errors.append(str(e))
        raise RuntimeError("Could not load genome codon counts for custom-CDS projection. " + ("; ".join(errors[-3:]) if errors else ""))

    def _workbook_meta_value(self, workbook_path, key, default=None):
        try:
            meta = pd.read_excel(workbook_path, sheet_name="Meta", dtype=object)
            if meta is None or meta.empty or not {"Key", "Value"}.issubset(set(meta.columns)):
                return default
            key_l = str(key).strip().lower()
            for _, row in meta.iterrows():
                if str(row.get("Key", "")).strip().lower() == key_l:
                    val = row.get("Value", default)
                    if pd.isna(val):
                        return default
                    return val
        except Exception:
            pass
        return default

    def _bool_from_meta(self, value, default):
        if value is None:
            return bool(default)
        s = str(value).strip().lower()
        if s in {"1", "true", "yes", "y", "on"}:
            return True
        if s in {"0", "false", "no", "n", "off"}:
            return False
        return bool(default)

    def _plain_counts_to_feature_df_for_projection(self, counts_df, workbook_path):
        """Convert plain codon counts to the feature matrix used for nearest-neighbour projection."""
        usage_basis = str(self._workbook_meta_value(workbook_path, "usage_basis", _usage_internal(self.usage_basis.get())) or "RCU").upper()
        codon_set = str(self._workbook_meta_value(workbook_path, "codon_set", self.codon_set.get().strip() or "59") or "59")
        df = counts_df.copy().astype(float)
        for cod in _CODON_ORDER:
            if cod not in df.columns:
                df[cod] = 0.0
        df = df.loc[:, _CODON_ORDER]

        if usage_basis == "AA":
            aa_cols = sorted({aa for aa in _CODON_TO_AA3.values() if aa})
            feat = pd.DataFrame(index=df.index, columns=aa_cols, dtype=float)
            for aa in aa_cols:
                codons = [c for c, a in _CODON_TO_AA3.items() if a == aa and c in df.columns]
                feat[aa] = df[codons].sum(axis=1) if codons else 0.0
            return feat

        if usage_basis == "ACU":
            feat = _metric_values_from_codon_counts(df, "Absolute codon frequency")
        else:
            feat = _metric_values_from_codon_counts(df, "Relative codon usage")

        cset = str(codon_set).lower()
        if "59" in cset:
            keep = [c for c in _CODON_ORDER if _CODON_TO_AA3.get(c) not in {"Met", "Trp"}]
        else:
            keep = list(_CODON_ORDER)
        return feat.loc[:, keep]

    def _normalize_projection_features(self, genome_feat, custom_feat, workbook_path):
        """Apply CodonPipe-like imputation/centering/scaling using genome statistics."""
        cols = [c for c in genome_feat.columns if c in custom_feat.columns]
        genome = genome_feat.loc[:, cols].copy().astype(float)
        custom = custom_feat.loc[:, cols].copy().astype(float)
        mu_for_impute = genome.replace([np.inf, -np.inf], np.nan).mean(axis=0, skipna=True).fillna(0.0)
        genome = genome.replace([np.inf, -np.inf], np.nan).fillna(mu_for_impute)
        custom = custom.replace([np.inf, -np.inf], np.nan).fillna(mu_for_impute)

        center = self._bool_from_meta(self._workbook_meta_value(workbook_path, "center_features", self.center_features.get()), self.center_features.get())
        scale = self._bool_from_meta(self._workbook_meta_value(workbook_path, "scale_features", self.scale_features.get()), self.scale_features.get())
        if center:
            mu = genome.mean(axis=0)
            genome = genome - mu
            custom = custom - mu
        if scale:
            sd = genome.std(axis=0, ddof=0).replace(0.0, 1.0).fillna(1.0)
            genome = genome / sd
            custom = custom / sd
        return genome.to_numpy(dtype=float), custom.to_numpy(dtype=float), list(genome.index), cols

    def _distance_vector_for_projection(self, genome_arr, custom_vec, metric):
        metric = str(metric or "euclidean").strip().lower()
        if metric in {"cosine", "correlation"}:
            A = np.asarray(genome_arr, dtype=float)
            b = np.asarray(custom_vec, dtype=float)
            if metric == "correlation":
                A = A - A.mean(axis=1, keepdims=True)
                b = b - b.mean()
            An = np.linalg.norm(A, axis=1)
            bn = float(np.linalg.norm(b))
            denom = An * bn
            sim = np.zeros(A.shape[0], dtype=float)
            ok = denom > 0
            sim[ok] = (A[ok] @ b) / denom[ok]
            return 1.0 - np.clip(sim, -1.0, 1.0)
        # cityblock/manhattan are treated as Manhattan; otherwise Euclidean.
        if metric in {"cityblock", "manhattan", "l1"}:
            return np.sum(np.abs(np.asarray(genome_arr, dtype=float) - np.asarray(custom_vec, dtype=float)), axis=1)
        return np.sqrt(np.sum((np.asarray(genome_arr, dtype=float) - np.asarray(custom_vec, dtype=float)) ** 2, axis=1))

    def _project_custom_cds_onto_existing_workbook(self, workbook_path, ordered_tags=None):
        """Project post-analysis custom CDS onto the existing embedding.

        UMAP/t-SNE coordinates are not available for custom CDS added after a run.
        To make them plottable without perturbing the original map, place each
        custom CDS at a weighted average of its nearest genome genes in the same
        codon-usage feature space used for the run. This is an overlay/projection,
        not a rerun of UMAP or clustering.
        """
        custom_counts = self._custom_cds_count_df_for_projection()
        if custom_counts is None or custom_counts.empty:
            return {}, pd.DataFrame(), []

        sheet_coords = CP.PIPELINE.get("sheet_coordinates", "Genes reordered")
        try:
            coord_df = pd.read_excel(workbook_path, sheet_name=sheet_coords, dtype=object).fillna("")
        except Exception as e:
            self._append_log(f"[WARN] Could not read coordinate sheet for custom-CDS projection: {e}\n")
            return {}, pd.DataFrame(), []
        if coord_df.empty or coord_df.shape[1] < 3:
            return {}, pd.DataFrame(), []

        cols = list(coord_df.columns)
        norm = {str(c).strip().lower().replace(" ", "").replace("_", ""): c for c in cols}
        locus_col = None
        for k in ("locustags", "locustag", "locus", "genomelocustags"):
            if k in norm:
                locus_col = norm[k]
                break
        if locus_col is None:
            locus_col = cols[0]
        x_col = None
        y_col = None
        for k in ("coordinates1", "x", "umapx", "tsnex", "pc1", "dim1"):
            if k in norm:
                x_col = norm[k]
                break
        for k in ("coordinates2", "y", "umapy", "tsney", "pc2", "dim2"):
            if k in norm:
                y_col = norm[k]
                break
        if x_col is None or y_col is None:
            # Legacy fallback: first three columns are locus/X/Y.
            locus_col, x_col, y_col = cols[0], cols[1], cols[2]

        coord_df = coord_df.copy()
        coord_df[locus_col] = coord_df[locus_col].astype(str).str.strip()
        existing_coord_tags = {str(v).strip() for v in coord_df[locus_col].astype(str).tolist() if str(v).strip() and str(v).strip().lower() != "nan"}
        coord_map = {}
        for _, row in coord_df.iterrows():
            tag = str(row.get(locus_col, "")).strip()
            if not tag or tag.lower() == "nan" or tag.startswith("custom_CDS_"):
                continue
            try:
                coord_map[tag] = (float(row.get(x_col)), float(row.get(y_col)))
            except Exception:
                continue
        if not coord_map:
            return {}, pd.DataFrame(), []

        try:
            genome_counts = self._load_genome_count_df_for_projection(workbook_path)
            common = [g for g in genome_counts.index.astype(str).tolist() if g in coord_map]
            genome_counts = genome_counts.loc[common]
            if genome_counts.empty:
                raise RuntimeError("No overlap between genome codon counts and workbook coordinates.")
            genome_feat = self._plain_counts_to_feature_df_for_projection(genome_counts, workbook_path)
            custom_feat = self._plain_counts_to_feature_df_for_projection(custom_counts, workbook_path)
            genome_arr, custom_arr, genome_ids, _feature_cols = self._normalize_projection_features(genome_feat, custom_feat, workbook_path)
        except Exception as e:
            self._append_log(f"[WARN] Could not prepare codon-usage features for custom-CDS projection: {e}\n")
            return {}, pd.DataFrame(), []

        dimred = str(self._workbook_meta_value(workbook_path, "dimred_method", _dimred_internal(self.dimred_method.get())) or "umap").lower()
        if dimred == "umap":
            metric = str(self._workbook_meta_value(workbook_path, "umap_metric", self.dimred_param_vars.get("umap", {}).get("umap_metric", tk.StringVar(value="cosine")).get()) or "cosine")
        elif dimred == "tsne":
            metric = str(self._workbook_meta_value(workbook_path, "tsne_distance", self.dimred_param_vars.get("tsne", {}).get("tsne_distance", tk.StringVar(value="cosine")).get()) or "cosine")
        else:
            metric = str(self._workbook_meta_value(workbook_path, "gene_dist_metric", self.cluster_param_vars.get("kmeans", {}).get("gene_dist_metric", tk.StringVar(value="euclidean")).get()) or "euclidean")

        rows = []
        members = []
        order = [str(x).strip() for x in list(ordered_tags or []) if str(x).strip()]
        order_index = {tag: i for i, tag in enumerate(order)}
        for custom_id, vec in zip(custom_feat.index.astype(str).tolist(), custom_arr):
            if custom_id in existing_coord_tags:
                # The custom CDS was already included when the workbook was generated;
                # keep its original embedding coordinate instead of projecting it.
                continue
            d = self._distance_vector_for_projection(genome_arr, vec, metric)
            if d.size == 0 or not np.isfinite(d).any():
                continue
            d = np.where(np.isfinite(d), d, np.inf)
            k = int(min(12, max(1, len(d))))
            nn = np.argsort(d)[:k]
            nearest_idx = int(nn[0])
            nearest_tag = genome_ids[nearest_idx]
            if d[nearest_idx] <= 1e-12:
                x, y = coord_map[nearest_tag]
            else:
                w = 1.0 / np.maximum(d[nn], 1e-9) ** 2
                coords = np.array([coord_map[genome_ids[i]] for i in nn], dtype=float)
                xy = np.sum(coords * w[:, None], axis=0) / np.sum(w)
                x, y = float(xy[0]), float(xy[1])
            insert_after = order_index.get(nearest_tag, 10**12)
            rows.append({
                "locus_tag": custom_id,
                "X": x,
                "Y": y,
                "nearest_locus_tag": nearest_tag,
                "insert_after_index": insert_after,
            })
            members.append(custom_id)

        if rows:
            self._append_log(f"[INFO] Projected {len(rows)} custom CDS record(s) onto the existing 2D map without rerunning clustering.\n")
        proj_df = pd.DataFrame(rows)
        cluster_name = str(getattr(self, "custom_cds_cluster_name", tk.StringVar(value="custom")).get() or "custom").strip() or "custom"
        return ({cluster_name: members} if members else {}), proj_df, members

    def _pad_dataframe_rows(self, df, n_rows):
        df = df.copy() if isinstance(df, pd.DataFrame) else pd.DataFrame()
        n_rows = int(max(n_rows, len(df)))
        if len(df) < n_rows:
            extra = pd.DataFrame("", index=range(n_rows - len(df)), columns=df.columns)
            df = pd.concat([df, extra], ignore_index=True)
        return df

    def _unique_gui_augmented_path(self, base_path, suffix):
        folder = os.path.dirname(os.path.abspath(base_path)) or "."
        stem, ext = os.path.splitext(os.path.basename(base_path))
        ext = ext or ".xlsx"
        candidate = os.path.join(folder, f"{stem}{suffix}{ext}")
        if os.path.abspath(candidate) == os.path.abspath(base_path):
            candidate = os.path.join(folder, f"{stem}{suffix}_copy{ext}")
        return candidate

    def _augment_workbook_with_current_metric_clusters(self, workbook_path):
        """Create a temporary workbook copy containing GUI-generated clusters.

        The original workbook is never modified. For FASTA-derived metric
        clusters, member genes already exist in the workbook. For custom CDS
        added after clustering, CodonPipe now projects them onto the existing 2D
        embedding using nearest-neighbour interpolation in codon-usage feature
        space, then adds them to the coordinate, locus-tag, and binary sheets so
        Figures-tab Plot 2 and Plot 3 can display the custom cluster.
        """
        if not workbook_path or not os.path.isfile(workbook_path):
            return workbook_path

        sheet_locus = CP.PIPELINE.get("sheet_locus_tags", "Locus Tags")
        sheet_binary = CP.PIPELINE.get("sheet_binary", "Binary")
        sheet_coords = CP.PIPELINE.get("sheet_coordinates", "Genes reordered")
        try:
            locus_df = pd.read_excel(workbook_path, sheet_name=sheet_locus, dtype=str).fillna("")
        except Exception:
            return workbook_path
        if locus_df.empty or locus_df.shape[1] < 1:
            return workbook_path

        genome_col = locus_df.columns[0]
        ordered_tags = [str(v).strip() for v in locus_df[genome_col].astype(str).tolist() if str(v).strip() and str(v).strip().lower() != "nan"]

        # Existing GUI-extra clusters whose genes already exist in the workbook
        # plus post-hoc custom CDS that can be projected into the existing map.
        gui_members = self._current_gui_extra_cluster_members(ordered_tags=ordered_tags, warn_if_custom_absent=True)
        projected_custom_members, projected_custom_df, _custom_ids = self._project_custom_cds_onto_existing_workbook(workbook_path, ordered_tags=ordered_tags)
        for cname, vals in (projected_custom_members or {}).items():
            existing = gui_members.setdefault(cname, [])
            seen = set(existing)
            for v in vals:
                if v not in seen:
                    existing.append(v)
                    seen.add(v)

        if not gui_members:
            return workbook_path

        existing_lower = {str(c).strip().lower() for c in list(locus_df.columns)[1:]}
        missing = {c: vals for c, vals in gui_members.items() if str(c).strip().lower() not in existing_lower}

        # Even if the custom-cluster column already exists, a post-hoc custom CDS
        # projection still requires adding coordinate/binary rows.
        has_projected_custom = isinstance(projected_custom_df, pd.DataFrame) and not projected_custom_df.empty
        if not missing and not has_projected_custom:
            return workbook_path

        # ----- Build augmented genome/cluster-order axis for the Locus Tags sheet.
        new_order = list(ordered_tags)
        if has_projected_custom:
            inserts = []
            for _, row in projected_custom_df.iterrows():
                tag = str(row.get("locus_tag", "")).strip()
                if not tag or tag in new_order:
                    continue
                pos = row.get("insert_after_index", 10**12)
                try:
                    pos = int(pos)
                except Exception:
                    pos = 10**12
                inserts.append((pos, tag))
            inserts.sort(key=lambda x: (x[0], x[1]))
            offset = 0
            for pos, tag in inserts:
                if pos >= len(new_order):
                    new_order.append(tag)
                else:
                    new_order.insert(pos + 1 + offset, tag)
                    offset += 1

        locus_cols = list(locus_df.columns)
        locus_aug_data = {genome_col: new_order}
        n_rows = len(new_order)
        for col in locus_cols[1:]:
            vals = [str(v).strip() for v in locus_df[col].astype(str).tolist() if str(v).strip() and str(v).strip().lower() != "nan"]
            locus_aug_data[col] = vals + [""] * max(0, n_rows - len(vals))
        for cname, members in missing.items():
            vals = list(members)
            locus_aug_data[cname] = vals + [""] * max(0, n_rows - len(vals))
        # If a custom cluster column already exists, extend it with the projected IDs.
        for cname, members in (projected_custom_members or {}).items():
            actual_col = None
            for col in list(locus_aug_data.keys())[1:]:
                if str(col).strip().lower() == str(cname).strip().lower():
                    actual_col = col
                    break
            if actual_col is not None:
                vals = [str(v).strip() for v in locus_aug_data[actual_col] if str(v).strip()]
                vals = list(dict.fromkeys(vals + list(members)))
                locus_aug_data[actual_col] = vals + [""] * max(0, n_rows - len(vals))
        # Normalize all column lengths.
        for col in list(locus_aug_data.keys()):
            vals = list(locus_aug_data[col])
            locus_aug_data[col] = vals[:n_rows] + [""] * max(0, n_rows - len(vals))
        locus_aug = pd.DataFrame(locus_aug_data)

        # ----- Build augmented Binary sheet, row-aligned to the new order.
        try:
            binary_df = pd.read_excel(workbook_path, sheet_name=sheet_binary, dtype=str).fillna("")
        except Exception:
            binary_df = pd.DataFrame({"Genome locus tags": ordered_tags, "Gene name": [""] * len(ordered_tags)})
        if binary_df.empty:
            binary_df = pd.DataFrame({"Genome locus tags": ordered_tags, "Gene name": [""] * len(ordered_tags)})
        binary_locus_col = binary_df.columns[0]
        binary_gene_col = binary_df.columns[1] if binary_df.shape[1] > 1 else None
        binary_df = binary_df.copy()
        binary_df[binary_locus_col] = binary_df[binary_locus_col].astype(str).str.strip()
        old_rows = {str(row.get(binary_locus_col, "")).strip(): row for _, row in binary_df.iterrows() if str(row.get(binary_locus_col, "")).strip()}
        binary_cols = list(binary_df.columns)
        for cname in missing.keys():
            if not any(str(c).strip().lower() == str(cname).strip().lower() for c in binary_cols):
                binary_cols.append(cname)
        for cname in projected_custom_members.keys():
            if not any(str(c).strip().lower() == str(cname).strip().lower() for c in binary_cols):
                binary_cols.append(cname)
        cluster_member_sets = {str(c): {str(v).strip() for v in vals} for c, vals in gui_members.items()}
        lower_cluster_cols = {str(c).strip().lower(): c for c in binary_cols}
        binary_rows = []
        for tag in new_order:
            row_out = {c: "" for c in binary_cols}
            row_out[binary_locus_col] = tag
            if tag in old_rows:
                old = old_rows[tag]
                for c in binary_cols:
                    if c in old.index:
                        row_out[c] = old[c]
            elif binary_gene_col is not None:
                row_out[binary_gene_col] = tag if tag.startswith("custom_CDS_") else ""
            for cname, member_set in cluster_member_sets.items():
                actual = lower_cluster_cols.get(str(cname).strip().lower(), cname)
                row_out[actual] = 1 if tag in member_set else row_out.get(actual, 0)
            binary_rows.append(row_out)
        binary_aug = pd.DataFrame(binary_rows, columns=binary_cols)

        # ----- Build augmented coordinate sheet for projected custom CDS.
        coord_aug = None
        if has_projected_custom:
            try:
                coord_df = pd.read_excel(workbook_path, sheet_name=sheet_coords, dtype=object).fillna("")
                cols = list(coord_df.columns)
                norm = {str(c).strip().lower().replace(" ", "").replace("_", ""): c for c in cols}
                locus_col = None
                for k in ("locustags", "locustag", "locus", "genomelocustags"):
                    if k in norm:
                        locus_col = norm[k]
                        break
                if locus_col is None:
                    locus_col = cols[0]
                x_col = None
                y_col = None
                for k in ("coordinates1", "x", "umapx", "tsnex", "pc1", "dim1"):
                    if k in norm:
                        x_col = norm[k]
                        break
                for k in ("coordinates2", "y", "umapy", "tsney", "pc2", "dim2"):
                    if k in norm:
                        y_col = norm[k]
                        break
                if x_col is None or y_col is None:
                    x_col, y_col = cols[1], cols[2]
                row_by_tag = {str(r.get(locus_col, "")).strip(): r for _, r in coord_df.iterrows() if str(r.get(locus_col, "")).strip()}
                proj_by_tag = {str(r.get("locus_tag", "")).strip(): r for _, r in projected_custom_df.iterrows()}
                rows = []
                for tag in new_order:
                    if tag in row_by_tag:
                        rows.append(dict(row_by_tag[tag]))
                    elif tag in proj_by_tag:
                        pr = proj_by_tag[tag]
                        d = {c: "" for c in cols}
                        d[locus_col] = tag
                        d[x_col] = float(pr.get("X"))
                        d[y_col] = float(pr.get("Y"))
                        # Fill optional human-readable fields when present.
                        for c in cols:
                            key = str(c).strip().lower().replace(" ", "").replace("_", "")
                            if key in {"genename", "gene"}:
                                d[c] = tag
                            elif key in {"proteindescription", "description", "product"}:
                                d[c] = f"Custom CDS projected near {pr.get('nearest_locus_tag', '')}"
                        rows.append(d)
                coord_aug = pd.DataFrame(rows, columns=cols)
            except Exception as e:
                self._append_log(f"[WARN] Could not augment coordinate sheet with projected custom CDS: {e}\n")
                coord_aug = None

        out_path = self._unique_gui_augmented_path(workbook_path, "__with_GUI_extra_clusters")
        try:
            if os.path.exists(out_path):
                os.remove(out_path)
            shutil.copy2(workbook_path, out_path)
            with pd.ExcelWriter(out_path, engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
                locus_aug.to_excel(writer, sheet_name=sheet_locus, index=False)
                binary_aug.to_excel(writer, sheet_name=sheet_binary, index=False)
                if coord_aug is not None:
                    coord_aug.to_excel(writer, sheet_name=sheet_coords, index=False)
                max_len = max((len(v) for v in gui_members.values()), default=0)
                if max_len > 0:
                    extra_df = pd.DataFrame({c: list(vals) + [""] * (max_len - len(vals)) for c, vals in gui_members.items()})
                    extra_df.to_excel(writer, sheet_name="GUI extra clusters", index=False)
                if has_projected_custom:
                    projected_custom_df.to_excel(writer, sheet_name="Projected custom CDS", index=False)
                metric_df = self._current_fasta_metric_cluster_df()
                if isinstance(metric_df, pd.DataFrame) and not metric_df.empty:
                    metric_df.to_excel(writer, sheet_name="FASTA metric clusters", index=False)
                scores_df = getattr(self, "fasta_metric_cluster_scores_df", None)
                if isinstance(scores_df, pd.DataFrame) and not scores_df.empty:
                    scores_df.to_excel(writer, sheet_name="FASTA metric scores", index=False)
            return out_path
        except Exception as e:
            self._append_log(f"[WARN] Could not create augmented GUI-cluster workbook; using original workbook. {e}\n")
            return workbook_path

    def _read_described_genelevel_sheet(self, workbook_path, sheet_name):
        for header_row in (0, 1):
            try:
                df = pd.read_excel(workbook_path, sheet_name=sheet_name, header=header_row)
            except Exception:
                continue
            if df is None or df.empty:
                continue
            cols = [str(c).strip().lower() for c in df.columns]
            if "locustag" in cols or "locus tag" in cols or "locus tags" in cols:
                return df
        return pd.DataFrame()

    def _read_described_summary_sheet(self, workbook_path, sheet_names):
        for sheet_name in sheet_names:
            try:
                for header_row in (0, 1):
                    df = pd.read_excel(workbook_path, sheet_name=sheet_name, header=header_row)
                    if df is None or df.empty:
                        continue
                    cols = [str(c).strip().lower() for c in df.columns]
                    if "codon" in cols:
                        return sheet_name, df
            except Exception:
                continue
        return "", pd.DataFrame()

    def _locate_whole_genome_raw_codon_usage_workbook(self, output_dir):
        raw_dir = os.path.join(output_dir, str(CP.SET.get("export_cluster_codon_usage_raw_subdir", "Raw codon usage tables")))
        candidates = []
        if os.path.isdir(raw_dir):
            for fname in os.listdir(raw_dir):
                low = fname.lower()
                if low.endswith(('.xlsx', '.xlsm', '.xls')) and ("whole genome" in low or "whole_genome" in low or "genome" in low):
                    candidates.append(os.path.join(raw_dir, fname))
        return self._find_latest_file(candidates)

    def _augment_codon_usage_workbook_with_current_metric_clusters(self, codon_usage_workbook, output_dir):
        """Create a temporary codon-usage summary workbook with GUI-extra cluster columns."""
        if not codon_usage_workbook or not os.path.isfile(codon_usage_workbook):
            return codon_usage_workbook
        gui_members = self._current_gui_extra_cluster_members()
        if not gui_members:
            return codon_usage_workbook
        raw_wb = self._locate_whole_genome_raw_codon_usage_workbook(output_dir)
        if not raw_wb:
            return codon_usage_workbook

        mode_specs = {
            "ACU": (["ACU per cluster", "ACU"], "ACU"),
            "RCU": (["RCU per cluster", "RCU"], "RCU"),
            "ZCU": (["z-scores per cluster", "ZCU", "Z"], "ZCU"),
        }
        replacements = {}
        for _mode, (summary_sheets, raw_sheet) in mode_specs.items():
            sheet_name, summary_df = self._read_described_summary_sheet(codon_usage_workbook, summary_sheets)
            raw_df = self._read_described_genelevel_sheet(raw_wb, raw_sheet)
            if not sheet_name or summary_df.empty or raw_df.empty:
                continue
            locus_col = None
            for col in raw_df.columns:
                if str(col).strip().lower().replace(" ", "") in {"locustag", "locustags", "locus"}:
                    locus_col = col
                    break
            if locus_col is None:
                locus_col = raw_df.columns[0]
            raw_df = raw_df.copy()
            raw_df[locus_col] = raw_df[locus_col].astype(str).str.strip()
            raw_df = raw_df.set_index(locus_col, drop=True)
            codon_col = None
            for col in summary_df.columns:
                if str(col).strip().lower() == "codon":
                    codon_col = col
                    break
            if codon_col is None:
                codon_col = summary_df.columns[0]
            out_df = summary_df.copy()
            codons = [str(v).strip() for v in out_df[codon_col].astype(str).tolist()]
            for cname, members in gui_members.items():
                if cname in out_df.columns:
                    continue
                genes_present = [g for g in members if g in raw_df.index]
                if not genes_present:
                    continue
                vals = []
                for codon in codons:
                    if codon not in raw_df.columns:
                        vals.append(np.nan)
                    else:
                        vals.append(pd.to_numeric(raw_df.loc[genes_present, codon], errors="coerce").mean(skipna=True))
                out_df[cname] = vals
            replacements[sheet_name] = out_df

        if not replacements:
            return codon_usage_workbook
        out_path = self._unique_gui_augmented_path(codon_usage_workbook, "__with_GUI_extra_clusters")
        try:
            if os.path.exists(out_path):
                os.remove(out_path)
            shutil.copy2(codon_usage_workbook, out_path)
            with pd.ExcelWriter(out_path, engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
                for sheet_name, df in replacements.items():
                    df.to_excel(writer, sheet_name=sheet_name, index=False)
            return out_path
        except Exception as e:
            self._append_log(f"[WARN] Could not create augmented codon-usage workbook with GUI-extra clusters; using original workbook. {e}\n")
            return codon_usage_workbook

    def _show_no_clustering_popup(self):
        messagebox.showinfo(
            "Run clustering first",
            "No completed CodonPipe clustering workbook was found for the current project.\n\n"
            "Please run 'Run clustering' in the Codon usage clustering tab first. "
            "After that, figure buttons will replot from the existing workbook without rerunning clustering or rewriting Excel files."
        )

    def _make_single_figure_pipeline_cfg(self, figure_mode, workbook_path):
        cfg = self._collect_pipeline_overrides()
        figure_clusters = self._figure_clusters_for_workbook(workbook_path)
        cfg["plotting_pipeline_script_path"] = PLOTTING_SCRIPT
        cfg["auto_run_plotting_pipeline"] = True
        cfg["plot_include_columns"] = list(figure_clusters)
        cfg["plot_highlight_columns"] = list(figure_clusters)
        cfg["gchm_include_columns"] = list(figure_clusters)
        cfg["RUN_2D_DENSITY_PLOTS"] = False
        cfg["RUN_GENE_CLUSTER_LOCALIZATION_HEATMAP"] = False
        cfg["CODON_USAGE_PLOT_MODE"] = "NONE"
        cfg["codon_usage_plot_mode"] = "NONE"
        cfg["SHOW_FIG"] = True
        cfg["GCHM_SHOW_FIG"] = True
        cfg["CODON_USAGE_SHOW_FIG"] = True
        cfg["MAX_NROWS"] = _safe_int(self.plot_rows.get(), 4)
        cfg["plot_max_nrows"] = _safe_int(self.plot_rows.get(), 4)
        if figure_mode == "figure_cluster_axis":
            cfg["RUN_GENE_CLUSTER_LOCALIZATION_HEATMAP"] = True
        elif figure_mode == "figure_2d_map":
            cfg["RUN_2D_DENSITY_PLOTS"] = True
        elif figure_mode == "figure_codon_profiles":
            mode_key = USER_CODON_MODE_TO_INTERNAL.get(self.codon_usage_plot_mode.get().strip(), "Z")
            cfg["CODON_USAGE_PLOT_MODE"] = mode_key
            cfg["codon_usage_plot_mode"] = mode_key
            output_dir = os.path.dirname(workbook_path)
            codon_wb = self._locate_codon_usage_workbook(output_dir)
            if codon_wb:
                codon_wb = self._augment_codon_usage_workbook_with_current_metric_clusters(codon_wb, output_dir)
                cfg["CODON_USAGE_WORKBOOK"] = codon_wb
        return cfg

    def _plot_main_heatmap_from_workbook(self, workbook_path):
        try:
            matrix_df = pd.read_excel(workbook_path, sheet_name="Heatmap matrix", dtype=object)
        except Exception as e:
            raise RuntimeError(
                "The selected clustering workbook does not contain the 'Heatmap matrix' sheet needed for replotting.\n\n"
                "Please rerun 'Run clustering' once with this updated CodonPipe version, then press PLOT again."
            ) from e
        if matrix_df.empty or matrix_df.shape[1] < 2:
            raise RuntimeError("The 'Heatmap matrix' sheet is empty or malformed.")
        feature_col = matrix_df.columns[0]
        features = matrix_df[feature_col].astype(str).tolist()
        data = matrix_df.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        if data.size == 0:
            raise RuntimeError("The 'Heatmap matrix' sheet does not contain numeric heatmap values.")
        n_genes_full = data.shape[1]
        bin_size = 1
        title = "Codons vs genes heatmap"
        try:
            meta = pd.read_excel(workbook_path, sheet_name="Heatmap metadata", dtype=object)
            if not meta.empty and {"Key", "Value"}.issubset(set(meta.columns)):
                meta_map = {str(r["Key"]): r["Value"] for _, r in meta.iterrows()}
                title = str(meta_map.get("title", title) or title)
                bin_size = _safe_int(meta_map.get("bin_size", 1), 1)
                n_genes_full = _safe_int(meta_map.get("n_genes_full", n_genes_full), n_genes_full)
        except Exception:
            pass

        set_overrides = self._collect_set_overrides()
        local_set = copy.deepcopy(CP.SET)
        local_set.update(set_overrides)
        custom_cmaps = {}
        try:
            from codonpipe.clustering import load_custom_colormaps, plot_heatmap
            custom_cmaps = load_custom_colormaps(local_set)
            fig, _ = plot_heatmap(
                local_set,
                data,
                features,
                n_genes_full,
                bin_size,
                title,
                custom_cmaps,
            )
        except Exception:
            # Lightweight fallback if the plotting helper is unavailable.
            fig, ax = plt.subplots(figsize=local_set.get("heatmap_fig_size", (18, 4)), dpi=local_set.get("figure_dpi", 300))
            im = ax.imshow(data, aspect="auto", origin="lower", cmap=local_set.get("heatmap_colormap_name", "viridis"))
            caxis = local_set.get("heatmap_caxis_limits", None)
            if caxis is not None:
                try:
                    im.set_clim(*caxis)
                except Exception:
                    pass
            fig.colorbar(im, ax=ax)
            ax.set_yticks(np.arange(len(features)))
            ax.set_yticklabels(features, fontsize=local_set.get("font_size_yticks", 3))
            ax.set_title(title)
            ax.xaxis.set_ticks_position('top')
            ax.xaxis.set_label_position('top')

        out_dir = os.path.dirname(workbook_path)
        fig_dir = os.path.join(out_dir, "Figures")
        os.makedirs(fig_dir, exist_ok=True)
        fmt = str(self.figure_format.get() or "png").strip().lstrip(".").lower() or "png"
        if fmt == "jpg":
            fmt = "jpeg"
        out_path = os.path.join(fig_dir, f"Codons vs genes 2D heatmap.{fmt}")
        fig.savefig(out_path, dpi=local_set.get("figure_dpi", 300), bbox_inches="tight")
        try:
            plt.show(block=False)
        except Exception:
            pass
        return out_path

    def _ordered_genes_from_clustering_workbook(self, workbook_path):
        """Read reordered locus tags from the clustering workbook."""
        sheet = CP.PIPELINE.get("sheet_coordinates", "Genes reordered")
        try:
            df = pd.read_excel(workbook_path, sheet_name=sheet, dtype=object)
        except Exception as e:
            raise RuntimeError(
                "Could not read the reordered genes from the clustering workbook.\n\n"
                f"Workbook: {workbook_path}\nSheet: {sheet}\n\n{e}"
            ) from e
        if df is None or df.empty:
            raise RuntimeError(f"The '{sheet}' sheet is empty in the clustering workbook.")
        cols = list(df.columns)
        locus_col = None
        for col in cols:
            key = str(col).strip().lower().replace(" ", "").replace("_", "")
            if key in {"locustags", "locustag", "locus", "genomelocustags"}:
                locus_col = col
                break
        if locus_col is None:
            locus_col = cols[0]
        ordered = [str(x).strip() for x in df[locus_col].tolist() if str(x).strip() and str(x).strip().lower() != "nan"]
        if not ordered:
            raise RuntimeError(f"No locus tags could be read from the '{sheet}' sheet.")
        return ordered

    def _locate_geneids_workbook_for_david(self, output_dir):
        candidates = [os.path.join(output_dir, "Gene IDs.xlsx")]
        try:
            for fname in os.listdir(output_dir):
                low = fname.lower()
                if low.endswith((".xlsx", ".xls", ".xlsm")) and ("gene id" in low or "geneids" in low or "gene_ids" in low):
                    candidates.append(os.path.join(output_dir, fname))
        except Exception:
            pass
        return self._find_latest_file(candidates)

    def _default_david_output_prefix(self, workbook_path):
        base = os.path.splitext(os.path.basename(str(workbook_path)))[0]
        if base.lower() == "clustering analysis results":
            org = self.organism_name.get().strip() or "Organism"
            return f"{org}_ClusteringAnalysis"
        return base

    def _looks_like_david_ssl_certificate_error(self, exc):
        txt = (str(exc or "") + "\n" + repr(exc or "")).lower()
        markers = [
            "certificate_verify_failed",
            "unable to get local issuer certificate",
            "ssl",
            "certificate verify failed",
        ]
        if "ssl" in txt and any(m in txt for m in markers[1:]):
            return True
        if any(m in txt for m in markers[:2]):
            return True
        cause = getattr(exc, "__cause__", None) or getattr(exc, "reason", None)
        if cause is not None and cause is not exc:
            return self._looks_like_david_ssl_certificate_error(cause)
        return False


    def _looks_like_david_identifier_error(self, exc):
        txt = (str(exc or "") + "\n" + repr(exc or "")).lower()
        markers = [
            "david did not recognize any identifier set",
            "david returned unrecognized",
            "unrecognized / non-annotated identifiers",
            "server raised fault: 'index: 0, size: 0'",
            "index: 0, size: 0",
        ]
        if any(m in txt for m in markers):
            return True
        cause = getattr(exc, "__cause__", None) or getattr(exc, "reason", None)
        if cause is not None and cause is not exc:
            return self._looks_like_david_identifier_error(cause)
        return False

    def run_david_sliding_window_scan(self):
        """Run DAVID sliding-window enrichment from the last clustering workbook only."""
        if getattr(self, "is_running", False):
            messagebox.showinfo("CodonPipe", "A run is already in progress.")
            return
        workbook = self._set_last_clustering_outputs_from_disk()
        if not workbook:
            self._show_no_clustering_popup()
            return
        email = self.david_email.get().strip()
        if not email:
            messagebox.showerror("DAVID email required", "Please enter a DAVID-registered email address before running the DAVID sliding-window scan.")
            return

        self.is_running = True
        self._set_run_buttons_state("disabled")
        self._append_log("\n" + "=" * 84 + "\nStarting DAVID sliding-window scan from existing clustering results...\n")
        self.root.update_idletasks()

        orig_stdout = sys.stdout
        orig_stderr = sys.stderr
        live_stdout = TeeTextRedirector(self._append_log, orig_stdout)
        live_stderr = TeeTextRedirector(self._append_log, orig_stderr)
        old_set = copy.deepcopy(CP.SET)
        try:
            output_dir = os.path.dirname(os.path.abspath(workbook))
            geneids_xlsx = self._locate_geneids_workbook_for_david(output_dir)
            if not geneids_xlsx:
                raise RuntimeError(
                    "The FASTA-derived Gene IDs workbook was not found in the clustering output folder.\n\n"
                    "Please rerun 'Run clustering' once, then run the DAVID scan again."
                )
            ordered_genes = self._ordered_genes_from_clustering_workbook(workbook)
            geneids_df = pd.read_excel(geneids_xlsx, dtype=str).fillna("")
            if geneids_df.empty:
                raise RuntimeError(f"The Gene IDs workbook is empty:\n{geneids_xlsx}")

            fmt = str(self.figure_format.get() or "png").strip().lstrip(".").lower() or "png"
            prefix = self._default_david_output_prefix(workbook)
            min_genes = _safe_int(self.plot_cluster_min_genes.get(), 2)

            CP.SET.update({
                "david_user_email": email,
                "david_window_size": _safe_int(self.david_window_size.get(), CP.SET.get("david_window_size", 100)),
                "david_step_size": _safe_int(self.david_step_size.get(), CP.SET.get("david_step_size", 50)),
                "david_wait_time": _safe_float(self.david_wait_time.get(), CP.SET.get("david_wait_time", 0.0)),
                "david_max_clusters": _safe_int(self.david_max_clusters.get(), CP.SET.get("david_max_clusters", 3)),
                "david_min_valid_ids_per_window": _safe_int(self.david_min_valid_ids.get(), CP.SET.get("david_min_valid_ids_per_window", 3)),
                "david_top_n_hits": _safe_int(self.david_top_n_hits.get(), CP.SET.get("david_top_n_hits", 10)),
                "david_plot_format": fmt,
            })

            with warnings.catch_warnings():
                with redirect_stdout(live_stdout), redirect_stderr(live_stderr):
                    print(f"[INFO] Using clustering workbook:\n  {workbook}")
                    print(f"[INFO] Using Gene IDs workbook:\n  {geneids_xlsx}")
                    print(f"[INFO] Ordered genes scanned by DAVID: {len(ordered_genes)}")
                    david_results = CP.run_david_window_scan_from_ordered_genes(
                        ordered_genes=list(ordered_genes),
                        geneids_df=geneids_df,
                        output_folder=output_dir,
                        output_prefix=prefix,
                        user_email=email,
                        alias_map={},
                        geneids_xlsx_path=geneids_xlsx,
                        window_size=int(CP.SET.get("david_window_size", 100)),
                        step_size=int(CP.SET.get("david_step_size", 50)),
                        wait_time=float(CP.SET.get("david_wait_time", 0.0)),
                        top_n_hits=int(CP.SET.get("david_top_n_hits", 10)),
                        max_clusters=int(CP.SET.get("david_max_clusters", 3)),
                        report_subdir_name=str(CP.SET.get("david_report_subdir_name", "DAVID window reports")),
                        min_valid_ids_per_window=int(CP.SET.get("david_min_valid_ids_per_window", 3)),
                        plot_format=fmt,
                        chart_threshold=float(CP.SET.get("david_chart_threshold", 1.0)),
                        chart_count=int(CP.SET.get("david_chart_count", 1)),
                        manual_term_queries=list(CP.SET.get("david_manual_term_queries", [])),
                        term_match_mode=str(CP.SET.get("david_term_match_mode", "contains")),
                        append_to_geneids_excel=bool(CP.SET.get("david_append_terms_to_geneids_excel", True)),
                    )
                    species_wb = CP._species_clusters_workbook_path(output_dir, self.organism_name.get().strip() or "Organism")
                    species_wb, david_sheet = CP._append_david_derived_clusters_sheet(
                        workbook_path=species_wb,
                        output_dir=output_dir,
                        organism_name=self.organism_name.get().strip() or "Organism",
                        david_results=david_results,
                        min_locus_tags=min_genes,
                    )
                    CP._organize_pipeline_outputs(
                        output_dir,
                        {"david_results": david_results},
                        move_root_files=False,
                        move_text_files=False,
                        move_figure_files=False,
                        move_david_outputs=True,
                    )
                    gene2terms = (david_results or {}).get("gene2terms_txt_path", "")
                    if gene2terms and os.path.isfile(gene2terms):
                        self.david_gene2terms_path.set(gene2terms)
                    print("[INFO] DAVID sliding-window scan completed.")
                    if species_wb:
                        msg = f"[INFO] DAVID-derived cluster table available in:\n  {species_wb}"
                        if david_sheet:
                            msg += f"\n  sheet: {david_sheet}"
                        print(msg)

            self._flush_log_queue_now()
            self._append_log("\n[GUI] DAVID sliding-window scan completed successfully.\n")
            self.root.update_idletasks()
            messagebox.showinfo("DAVID sliding-window scan", "DAVID sliding-window scan completed successfully.")
        except Exception as e:
            self._flush_log_queue_now()
            tb = traceback.format_exc()
            self._append_log("\n[GUI] DAVID sliding-window scan failed:\n" + tb + "\n")
            self.root.update_idletasks()
            if self._looks_like_david_ssl_certificate_error(e):
                messagebox.showerror(
                    "DAVID SSL certificate error",
                    "CodonPipe could not establish a verified HTTPS connection to the DAVID web service.\n\n"
                    "The patched DAVID connector first tries Python's normal certificate store, then the certifi CA bundle, "
                    "and finally a compatibility fallback with certificate verification disabled. If this message still appears, "
                    "please update your conda certificates with:\n\n"
                    "conda install -c conda-forge certifi ca-certificates openssl\n\n"
                    "Then restart Spyder and run the DAVID scan again. The full traceback is kept in the log."
                )
            elif self._looks_like_david_identifier_error(e):
                messagebox.showerror(
                    "DAVID identifiers not recognized",
                    "DAVID is reachable, but it did not recognize/annotate the identifiers in the current Gene IDs workbook.\n\n"
                    "For bacterial genomes, locus tags such as SL1344_* are often rejected by DAVID. The most reliable fix is to rerun CodonPipe with an NCBI RefSeq CDS FASTA that contains db_xref=GeneID entries, or to provide a Gene IDs workbook containing EntrezGeneID values.\n\n"
                    "The detailed DAVID attempts and identifier counts are kept in the log panel."
                )
            else:
                messagebox.showerror("DAVID sliding-window scan", str(e) if str(e) else "DAVID scan failed. See the log panel for details.")
        finally:
            CP.SET.clear()
            CP.SET.update(old_set)
            try:
                sys.stdout = orig_stdout
                sys.stderr = orig_stderr
            except Exception:
                pass
            self.is_running = False
            self._set_run_buttons_state("normal")

    def plot_existing_figure(self, figure_mode):
        if getattr(self, "is_running", False):
            messagebox.showinfo("CodonPipe", "A run is already in progress.")
            return
        workbook = self._set_last_clustering_outputs_from_disk()
        if not workbook:
            self._show_no_clustering_popup()
            return
        self.is_running = True
        self._set_run_buttons_state("disabled")
        self._append_log("\n" + "=" * 84 + f"\nReplotting from existing clustering workbook [{figure_mode}]...\n")
        self.root.update_idletasks()
        orig_stdout = sys.stdout
        orig_stderr = sys.stderr
        live_stdout = TeeTextRedirector(self._append_log, orig_stdout)
        live_stderr = TeeTextRedirector(self._append_log, orig_stderr)
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message=r"n_jobs value 1 overridden to 1 by setting random_state.*")
                with redirect_stdout(live_stdout), redirect_stderr(live_stderr):
                    if figure_mode == "figure_main_heatmap":
                        out_path = self._plot_main_heatmap_from_workbook(workbook)
                        print(f"[INFO] Codons vs genes heatmap replotted from existing workbook:\n  {out_path}")
                    else:
                        from codonpipe.density_bridge import run_density_plot_script
                        workbook_for_plot = self._augment_workbook_with_current_metric_clusters(workbook)
                        if workbook_for_plot != workbook:
                            print(f"[INFO] Using temporary workbook augmented with GUI-generated extra clusters:\n  {workbook_for_plot}")
                        cfg = self._make_single_figure_pipeline_cfg(figure_mode, workbook_for_plot)
                        run_density_plot_script(workbook_for_plot, cfg)
            self._flush_log_queue_now()
            self._append_log("\n[GUI] Figure plotting completed successfully.\n")
            self.root.update_idletasks()
            messagebox.showinfo("CodonPipe", "Figure plotting completed successfully.")
        except Exception as e:
            self._flush_log_queue_now()
            tb = traceback.format_exc()
            self._append_log("\n[GUI] Figure plotting failed:\n" + tb + "\n")
            self.root.update_idletasks()
            messagebox.showerror("CodonPipe", str(e) if str(e) else "Figure plotting failed. See the log panel for details.")
        finally:
            try:
                sys.stdout = orig_stdout
                sys.stderr = orig_stderr
            except Exception:
                pass
            self.is_running = False
            self._set_run_buttons_state("normal")


    # ------------------------- direct decoding replotting -------------------------
    def _display_codon_label_to_internal_for_decoding(self, label):
        """Convert workbook display codon labels to internal CodonPipe labels."""
        s = str(label or "").strip()
        if not s or "_" not in s:
            return s
        aa, cod = s.split("_", 1)
        cod = cod.strip().upper().replace("U", "T")
        return f"{aa.strip()}_{cod}"

    def _locate_raw_whole_genome_codon_counts_workbook(self, output_dir):
        candidates = []
        raw_dir = os.path.join(str(output_dir or ""), "Raw codon usage tables")
        for fname in ["Whole genome.xlsx", "Whole genome 2.xlsx", "whole genome.xlsx"]:
            candidates.append(os.path.join(raw_dir, fname))
        try:
            for fname in os.listdir(raw_dir):
                low = fname.lower()
                if low.endswith((".xlsx", ".xls", ".xlsm")) and "whole" in low and "genome" in low:
                    candidates.append(os.path.join(raw_dir, fname))
        except Exception:
            pass
        return self._find_latest_file(candidates)

    def _load_count_df_for_decoding_from_existing_outputs(self, workbook_path):
        """Load per-gene codon counts without rerunning clustering.

        Preferred path: ``Raw codon usage tables/Whole genome.xlsx`` created by a
        previous clustering run. CodonPipe raw workbooks contain a description row
        above the real header, so this loader first reads ``header=1`` and only
        falls back to ``header=0`` for older files. This is critical: reading the
        description row as the header produces a table whose gene IDs and codon
        columns cannot match the clustering workbook, which then yields empty
        decoding plots.

        Fallback: parse the selected FASTA only to count codons. This does not run
        dimensional reduction, clustering, 2D density plots, or any full pipeline
        step.
        """
        def _looks_like_codon_count_table(df):
            if df is None or df.empty or df.shape[1] < 2:
                return False
            cols = [self._display_codon_label_to_internal_for_decoding(c) for c in list(df.columns)[1:]]
            codon_like = 0
            for c in cols:
                ss = str(c)
                if "_" in ss:
                    aa, cod = ss.split("_", 1)
                    cod = cod.strip().upper().replace("U", "T")
                    if len(cod) == 3 and all(ch in "ACGT" for ch in cod):
                        codon_like += 1
            return codon_like >= 20

        def _read_raw_codon_counts_workbook(raw_wb):
            errors = []
            # New CodonPipe raw workbooks: row 0 = description, row 1 = header.
            # Older/ad-hoc workbooks: row 0 = header.
            for header_row in (1, 0):
                try:
                    df = pd.read_excel(raw_wb, sheet_name="Codon counts", dtype=object, header=header_row)
                except Exception as e:
                    errors.append(str(e))
                    continue
                if not _looks_like_codon_count_table(df):
                    continue
                locus_col = df.columns[0]
                out = df.copy()
                out[locus_col] = out[locus_col].astype(str).str.strip()
                out = out[out[locus_col].astype(str).str.lower() != "nan"]
                out = out[out[locus_col].astype(str).str.strip() != ""]
                out = out.set_index(locus_col)
                out.columns = [self._display_codon_label_to_internal_for_decoding(c) for c in out.columns]
                # Drop duplicated codon columns after U/T normalization by summing them.
                out = out.apply(pd.to_numeric, errors="coerce").fillna(0.0)
                out = out.T.groupby(level=0).sum().T
                out.index = out.index.astype(str).str.strip()
                out = out[~out.index.duplicated(keep="first")]
                return out, header_row
            raise RuntimeError(
                "Could not recognize the 'Codon counts' sheet in the raw workbook. "
                "Expected a first column of locus tags followed by codon columns such as Ala_GCA or Phe_TTT."
            )

        output_dir = os.path.dirname(os.path.abspath(workbook_path))
        raw_wb = self._locate_raw_whole_genome_codon_counts_workbook(output_dir)
        if raw_wb:
            out, header_row = _read_raw_codon_counts_workbook(raw_wb)
            print(f"[INFO] Loaded codon counts from existing raw workbook (header row={header_row}):\n  {raw_wb}")
            print(f"[INFO] Raw count table: {out.shape[0]} genes x {out.shape[1]} codon columns.")
            return out

        # Fallback: this is not a clustering rerun. It only counts codons from FASTA.
        fasta = self.fasta_path.get().strip()
        if not fasta or not os.path.isfile(fasta):
            raise RuntimeError(
                "Could not find 'Raw codon usage tables/Whole genome.xlsx' and the FASTA fallback is unavailable.\n\n"
                "Please rerun CodonPipe once with codon-usage table export enabled, or select the FASTA used for this project."
            )
        try:
            from codonpipe.fasta_metrics import compute_codon_usage_tables_from_cds_fasta
            _freq_df, _geneids_df, abs_df, _gc = compute_codon_usage_tables_from_cds_fasta(
                fasta_path=fasta,
                row_id_mode="locus",
                trim_to_multiple_of_3=bool(CP.SET.get("fasta_trim_to_multiple_of_3", True)),
                include_stops=False,
                keep_first_duplicate=True,
                organism_mode="prokaryote",
                codon_range=_validate_codon_range_text(self.fasta_codon_range.get()),
            )
            abs_df = abs_df.copy()
            abs_df.columns = [self._display_codon_label_to_internal_for_decoding(c) for c in abs_df.columns]
            abs_df = abs_df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
            print("[INFO] Raw codon-count workbook was not found; counted codons directly from FASTA without rerunning clustering.")
            return abs_df
        except Exception as e:
            raise RuntimeError(
                "Could not reconstruct codon counts for decoding plots.\n\n"
                "Preferred fix: run CodonPipe once so the 'Raw codon usage tables/Whole genome.xlsx' file is generated."
            ) from e

    def _cluster_df_from_clustering_workbook(self, workbook_path):
        sheet = CP.PIPELINE.get("sheet_locus_tags", "Locus Tags")
        try:
            df = pd.read_excel(workbook_path, sheet_name=sheet, dtype=object)
        except Exception as e:
            raise RuntimeError(
                f"Could not read cluster columns from the existing clustering workbook.\nWorkbook: {workbook_path}\nSheet: {sheet}"
            ) from e
        if df is None or df.empty or df.shape[1] < 2:
            raise RuntimeError(f"The '{sheet}' sheet does not contain cluster columns.")
        cluster_df = df.iloc[:, 1:].copy().fillna("")
        cluster_df.columns = [str(c).strip() for c in cluster_df.columns]
        return cluster_df

    def _decoding_selected_clusters_for_workbook(self, cluster_df):
        available = [str(c) for c in list(cluster_df.columns)]
        self._ensure_default_decoding_preselection_for_available(available)
        if self.decoding_clusters_selection is None:
            return available
        lower = {str(c).strip().lower(): str(c) for c in available}
        chosen = []
        for c in list(self.decoding_clusters_selection or []):
            mapped = lower.get(str(c).strip().lower())
            if mapped and mapped not in chosen:
                chosen.append(mapped)
        return chosen

    def _validate_decoding_existing_inputs(self):
        workbook = self._set_last_clustering_outputs_from_disk()
        if not workbook:
            self._show_no_clustering_popup()
            return "", ""
        trna = self.trna_decoding_table_path.get().strip() or self._auto_find_decoding_table_for_current_genome()
        if trna:
            self.trna_decoding_table_path.set(trna)
        if not trna or not os.path.isfile(trna):
            messagebox.showerror(
                "Decoding table required",
                "Please select a decoding table in the Input/Output tab, or place a 'species name decoding table.xlsx' file in the Preloaded genomes subfolder."
            )
            return "", ""
        return workbook, trna

    def _apply_decoding_run_mode_to_set(self, mode, S):
        for key in [
            "trna_gene_heatmap_enable",
            "trna_single_box_codon_heatmap_enable",
            "trna_mrna_stability_enable",
            "trna_shift_heatmap_enable",
            "trna_wobble_heatmap_enable",
            "trna_modification_heatmap_enable",
        ]:
            S[key] = False
        if mode == "trna_all":
            S["trna_single_box_codon_heatmap_enable"] = True   # Plot 1
            S["trna_gene_heatmap_enable"] = True               # Plot 2
            S["trna_mrna_stability_enable"] = True             # Plot 3
            S["trna_wobble_heatmap_enable"] = True             # Plot 4
            S["trna_shift_heatmap_enable"] = True              # Plot 5
            S["trna_modification_heatmap_enable"] = True       # Plot 6
        elif mode == "trna_figure_single_box":
            S["trna_single_box_codon_heatmap_enable"] = True
        elif mode == "trna_figure_gene_heatmap":
            S["trna_gene_heatmap_enable"] = True
        elif mode == "trna_figure_rna_stability":
            S["trna_mrna_stability_enable"] = True
        elif mode == "trna_figure_wobble":
            S["trna_wobble_heatmap_enable"] = True
        elif mode == "trna_figure_shift":
            S["trna_shift_heatmap_enable"] = True
        elif mode == "trna_figure_modifications":
            S["trna_modification_heatmap_enable"] = True

    def plot_existing_decoding_figure(self, run_mode="trna_all"):
        """Render decoding-strategy plots from existing clustering outputs only."""
        if getattr(self, "is_running", False):
            messagebox.showinfo("CodonPipe", "A run is already in progress.")
            return
        workbook, trna = self._validate_decoding_existing_inputs()
        if not workbook:
            return
        self.is_running = True
        self._set_run_buttons_state("disabled")
        self._append_log("\n" + "=" * 84 + f"\nPlotting decoding strategies from existing clustering workbook [{run_mode}]...\n")
        self._append_log(f"[GUI build] {CODONPIPE_GUI_BUILD} — direct decoding plotting; no clustering rerun.\n")
        self.root.update_idletasks()
        orig_stdout = sys.stdout
        orig_stderr = sys.stderr
        live_stdout = TeeTextRedirector(self._append_log, orig_stdout)
        live_stderr = TeeTextRedirector(self._append_log, orig_stderr)
        try:
            with warnings.catch_warnings():
                with redirect_stdout(live_stdout), redirect_stderr(live_stderr):
                    output_dir = os.path.dirname(os.path.abspath(workbook))
                    figures_dir = os.path.join(output_dir, "Figures")
                    os.makedirs(figures_dir, exist_ok=True)
                    print(f"[INFO] Reusing clustering workbook:\n  {workbook}")
                    print(f"[INFO] Reusing decoding table:\n  {trna}")
                    ordered_genes = self._ordered_genes_from_clustering_workbook(workbook)
                    cluster_df = self._cluster_df_from_clustering_workbook(workbook)
                    selected_clusters = self._decoding_selected_clusters_for_workbook(cluster_df)
                    count_df = self._load_count_df_for_decoding_from_existing_outputs(workbook)
                    S = self._collect_set_overrides()
                    S["export_trna_usage_enable"] = True
                    S["trna_decoding_table_path"] = trna
                    S["trna_decoding_table_sheet"] = self.trna_decoding_table_sheet.get().strip()
                    S["figure_output_format"] = str(self.figure_format.get() or "png").strip().lstrip(".").lower() or "png"
                    self._apply_decoding_run_mode_to_set(str(run_mode), S)
                    outputs = CP.render_trna_gene_ordered_heatmaps(
                        SET=S,
                        count_df=count_df,
                        ordered_genes=list(ordered_genes),
                        output_dir=figures_dir,
                        cluster_df=cluster_df,
                        selected_clusters=list(selected_clusters),
                    )
                    if outputs:
                        for key, path in outputs.items():
                            print(f"[INFO] {key}: {path}")
                    else:
                        print("[WARN] No decoding-strategy figure was generated. Check the decoding table, selected clusters and selected plot options.")
            self._flush_log_queue_now()
            self._append_log("\n[GUI] Decoding-strategy figure plotting completed successfully.\n")
            self.root.update_idletasks()
            messagebox.showinfo("CodonPipe", "Decoding-strategy figure plotting completed successfully.")
        except Exception as e:
            self._flush_log_queue_now()
            tb = traceback.format_exc()
            self._append_log("\n[GUI] Decoding plotting failed:\n" + tb + "\n")
            self.root.update_idletasks()
            messagebox.showerror("CodonPipe", str(e) if str(e) else "Decoding plotting failed. See the log panel for details.")
        finally:
            try:
                sys.stdout = orig_stdout
                sys.stderr = orig_stderr
            except Exception:
                pass
            self.is_running = False
            self._set_run_buttons_state("normal")

    # ------------------------- run pipeline -------------------------
    def _set_run_buttons_state(self, state):
        for btn in getattr(self, "run_buttons", []):
            try:
                btn.configure(state=state)
            except Exception:
                pass

    def run_pipeline(self, run_mode="clustering_plots"):
        # Decoding-strategy PLOT buttons must reuse existing clustering outputs.
        # Keep this guard so older saved GUI states or external callbacks cannot
        # accidentally restart the full codon-usage clustering pipeline.
        if str(run_mode or "").startswith("trna_"):
            self.plot_existing_decoding_figure(run_mode)
            return
        if getattr(self, "is_running", False):
            messagebox.showinfo("CodonPipe", "A run is already in progress.")
            return

        self._current_run_mode = str(run_mode or "clustering_plots")
        try:
            self._validate_inputs()
        except Exception as e:
            messagebox.showerror("Invalid settings", str(e))
            return

        self.is_running = True
        self._set_run_buttons_state("disabled")
        self._append_log("\n" + "=" * 84 + f"\nStarting CodonPipe run from GUI [{self._current_run_mode}]...\n")
        self.root.update_idletasks()
        self._run_pipeline_worker()

    def _looks_like_excel_file_lock_error(self, exc, traceback_text=""):
        """Return True for common Windows/Excel locked-workbook save failures."""
        msg = f"{exc}\n{traceback_text}".lower()
        lock_markers = (
            "filecreateerror",
            "permissionerror",
            "winerror 32",
            "being used by another process",
            "cannot access the file",
            "xlsxwriter",
        )
        return any(m in msg for m in lock_markers) and (
            "xlsx" in msg or "excelwriter" in msg or "tmp" in msg or "temp" in msg
        )

    def _show_excel_file_lock_error(self):
        messagebox.showerror(
            "Excel file is open or locked",
            "CodonPipe could not save one of the Excel output files.\n\n"
            "Most often, this happens because the previous CodonPipe Excel "
            "workbook is still open in Excel, or because Windows/Dropbox/OneDrive "
            "is temporarily locking the file.\n\n"
            "Please close any CodonPipe Excel workbooks in the output folder "
            "and run the analysis again. If the problem persists, try using a "
            "fresh empty output folder or pause cloud synchronization briefly."
        )

    def _run_pipeline_worker(self):
        runtime_choices = self._collect_runtime_choices()
        set_overrides = self._collect_set_overrides()
        pipeline_overrides = self._collect_pipeline_overrides()
        ks_overrides = self._collect_ks_overrides()
        run_mode = getattr(self, "_current_run_mode", "clustering_plots")
        self._apply_run_mode_overrides(run_mode, runtime_choices, set_overrides, pipeline_overrides)
        if str(run_mode).startswith("trna_"):
            # Decoding-strategy plots have their own cluster picker.
            active_clusters = self.decoding_clusters_selection
        elif str(run_mode) == "clustering_plots":
            # Full clustering runs must obey the Input/Output tab's
            # "Active clusters" selection. The Figures-tab "Cluster picker"
            # is intentionally reserved for replotting existing workbooks from
            # the Figures tab, so the two selectors do not interfere.
            active_clusters = self.active_clusters_selection
        else:
            active_clusters = self.active_clusters_selection
        min_genes = _safe_int(self.plot_cluster_min_genes.get(), 2)
        refined_cluster_path = self.refined_cluster_file.get().strip()
        refined_cluster_sheet = self.refined_cluster_sheet.get().strip()

        old_set = copy.deepcopy(CP.SET)
        old_pipeline = copy.deepcopy(CP.PIPELINE)
        old_ks = copy.deepcopy(CP.KS_SETTINGS)
        old_prompt_runtime_choices = CP._prompt_runtime_choices
        old_prompt_plot_grid_rows = CP._prompt_plot_grid_rows
        old_prompt_codon_usage_plot_mode = CP._prompt_codon_usage_cluster_plot_mode
        old_choose_clusters_for_display = CP._choose_clusters_for_display
        old_load_refined_cluster_df = CP._load_refined_cluster_df

        orig_stdout = sys.stdout
        orig_stderr = sys.stderr
        live_stdout = TeeTextRedirector(self._append_log, orig_stdout)
        live_stderr = TeeTextRedirector(self._append_log, orig_stderr)

        try:
            CP.SET.update(set_overrides)
            CP.PIPELINE.update(pipeline_overrides)
            CP.KS_SETTINGS.clear()
            CP.KS_SETTINGS.update(ks_overrides)

            CP._prompt_runtime_choices = lambda: runtime_choices
            CP._prompt_plot_grid_rows = lambda default_rows=2: pipeline_overrides["plot_max_nrows"]
            CP._prompt_codon_usage_cluster_plot_mode = lambda default_key="Z": pipeline_overrides["codon_usage_plot_mode"]
            CP._choose_clusters_for_display = self._make_choose_clusters_func(active_clusters, min_genes)
            if runtime_choices.get("cluster_source") == "refined":
                CP._load_refined_cluster_df = self._make_load_refined_func(refined_cluster_path, refined_cluster_sheet)

            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message=r"n_jobs value 1 overridden to 1 by setting random_state.*")
                warnings.filterwarnings("ignore", message=r".*Found Intel OpenMP.*LLVM OpenMP.*", category=RuntimeWarning)
                with redirect_stdout(live_stdout), redirect_stderr(live_stderr):
                    CP.main()

            self._flush_log_queue_now()
            self._append_log("\n[GUI] CodonPipe run completed successfully.\n")
            self.root.update_idletasks()
            messagebox.showinfo("CodonPipe", "Run completed successfully.")

        except Exception as e:
            self._flush_log_queue_now()
            tb = traceback.format_exc()
            self._append_log("\n[GUI] Run failed:\n" + tb + "\n")
            self.root.update_idletasks()
            if self._looks_like_excel_file_lock_error(e, tb):
                self._show_excel_file_lock_error()
            else:
                messagebox.showerror("CodonPipe", "Run failed. See the log panel for details.")
        finally:
            CP.SET.clear()
            CP.SET.update(old_set)
            CP.PIPELINE.clear()
            CP.PIPELINE.update(old_pipeline)
            CP.KS_SETTINGS.clear()
            CP.KS_SETTINGS.update(old_ks)
            CP._prompt_runtime_choices = old_prompt_runtime_choices
            CP._prompt_plot_grid_rows = old_prompt_plot_grid_rows
            CP._prompt_codon_usage_cluster_plot_mode = old_prompt_codon_usage_plot_mode
            CP._choose_clusters_for_display = old_choose_clusters_for_display
            CP._load_refined_cluster_df = old_load_refined_cluster_df
            try:
                sys.stdout = orig_stdout
                sys.stderr = orig_stderr
            except Exception:
                pass
            self.is_running = False
            self._set_run_buttons_state("normal")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def _bring_main_window_to_front(root):
    """Best-effort helper so Spyder/Windows users can see the Tk window immediately."""
    try:
        root.deiconify()
        root.update_idletasks()
        root.lift()
        root.focus_force()
        # On Windows, a short topmost pulse helps when Spyder keeps focus.
        try:
            root.attributes("-topmost", True)
            root.after(700, lambda: root.attributes("-topmost", False))
        except Exception:
            pass
    except Exception:
        pass


def main():
    print(f"[GUI build] {CODONPIPE_GUI_BUILD} — launching CodonPipe GUI.", flush=True)
    try:
        root = tk.Tk()
        print("[GUI] Tk root created; building interface...", flush=True)
        app = CodonPipeGUI(root)
        root.after(150, lambda: _bring_main_window_to_front(root))
        print("[GUI] Tkinter window created. If it is not visible, check the taskbar or behind Spyder.", flush=True)
        root.mainloop()
    except Exception:
        tb = traceback.format_exc()
        print("[GUI] Startup failed:\n" + tb, flush=True)
        try:
            messagebox.showerror("CodonPipe startup failed", tb)
        except Exception:
            pass
        raise


if __name__ == "__main__":
    main()
