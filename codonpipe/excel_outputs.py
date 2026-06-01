"""codonpipe.excel_outputs

Excel I/O utilities for the codon-usage clustering pipeline.
"""

# codonpipe/excel_outputs.py
import os
import re
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm


def read_locus_clusters(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Cluster file not found: {path}")
    ext = os.path.splitext(path)[1].lower()
    if ext in (".xlsx", ".xls", ".xlsm"):
        df = pd.read_excel(path, dtype=str)
    else:
        df = pd.read_csv(path, sep=None, engine="python", dtype=str)

    for col in df.columns:
        df[col] = df[col].replace({np.nan: ""}).astype(str).map(lambda x: x.strip().strip('"').strip("'") if x != "nan" else "")
    return df


def clusters_from_df(df_clusters):
    clusters = {}
    for col in df_clusters.columns:
        vals = [v for v in df_clusters[col].tolist() if isinstance(v, str) and v != ""]
        clusters[col] = set(vals)
    return clusters


def build_locus_tags_sheet(ordered_tags, df_clusters, genome_colname="Genome locus tags"):
    ordered_tags = list(ordered_tags)
    n_genome = len(ordered_tags)
    n_rows_clusters = len(df_clusters)
    n_rows = max(n_genome, n_rows_clusters)

    data = {genome_colname: [""] * n_rows}
    for i, tag in enumerate(ordered_tags):
        data[genome_colname][i] = tag

    for col in df_clusters.columns:
        col_vals = [v if isinstance(v, str) else "" for v in df_clusters[col].tolist()]
        if len(col_vals) < n_rows:
            col_vals += [""] * (n_rows - len(col_vals))
        else:
            col_vals = col_vals[:n_rows]
        data[col] = col_vals

    return pd.DataFrame(data)


def build_binary_sheet(ordered_tags, df_clusters, gene_symbol_map, locus_index):
    ordered_tags = list(ordered_tags)
    clusters = clusters_from_df(df_clusters)

    gene_names = []
    missing_gene = 0
    for lt in ordered_tags:
        g = (gene_symbol_map.get(lt, "") or "") if gene_symbol_map else ""
        if not g:
            rec = locus_index.get(lt)
            if rec is not None:
                g = rec.get("gene_name", "") or ""
        if not g:
            missing_gene += 1
        gene_names.append(g)

    data = {"Genome locus tags": ordered_tags, "Gene name": gene_names}
    for cname, members in clusters.items():
        data[cname] = [1 if lt in members else 0 for lt in ordered_tags]

    if missing_gene:
        print(f"[INFO] Gene names missing for {missing_gene} of {len(ordered_tags)} locus tags (left blank in 'Binary').")

    return pd.DataFrame(data)


def build_coordinates_df(row_names, Y, ordered_genes=None, gene_symbol_map=None, gene_desc_map=None, locus_index=None):
    """Build the coordinate sheet exported to the main workbook.

    New default behavior (used by Clustering_Pipeline.py): return a headered
    ``Genes reordered`` table with columns:
      - locus tags
      - gene name
      - coordinates 1
      - coordinates 2
      - protein description

    The row order follows ``ordered_genes`` when provided. Legacy behavior is
    preserved when ``ordered_genes`` is omitted.
    """
    if Y is None or Y.shape[0] == 0:
        raise ValueError("No DR coordinates available (Y is None or empty).")

    Y2 = np.column_stack([Y[:, 0], np.zeros_like(Y[:, 0])]) if Y.shape[1] < 2 else Y[:, :2]

    if ordered_genes is None:
        return pd.DataFrame(np.column_stack([np.asarray(row_names, dtype=object), Y2]))

    coords = {str(tag): (float(x), float(y)) for tag, (x, y) in zip(list(row_names), Y2.tolist())}
    rows = []
    for lt in list(ordered_genes):
        lt_s = str(lt)
        x, y = coords.get(lt_s, (np.nan, np.nan))
        gene_name = ""
        protein_desc = ""
        if gene_symbol_map:
            gene_name = (gene_symbol_map.get(lt_s, "") or "").strip()
        if gene_desc_map:
            protein_desc = (gene_desc_map.get(lt_s, "") or "").strip()
        rec = locus_index.get(lt_s) if locus_index else None
        if rec is not None:
            if not gene_name:
                gene_name = str(rec.get("gene_name", "") or "").strip()
            if not protein_desc:
                protein_desc = str(rec.get("protein_description", "") or rec.get("product", "") or "").strip()
        rows.append({
            "locus tags": lt_s,
            "gene name": gene_name,
            "coordinates 1": x,
            "coordinates 2": y,
            "protein description": protein_desc,
        })
    return pd.DataFrame(rows, columns=["locus tags", "gene name", "coordinates 1", "coordinates 2", "protein description"])


def build_meta_table(SET, PIPELINE, clustering_results,
                     fasta_path, locustags_path,
                     summary_df, output_subfolder):
    # Meta sheet writer (short implementation; behavior unchanged).
    # By default, paths are redacted to basenames so exported workbooks remain
    # suitable for sharing. Set SET["redact_output_paths"] = False to retain
    # full local paths for private reproducibility records.
    rows = []
    add = rows.append

    redact_paths = bool(SET.get("redact_output_paths", True)) if isinstance(SET, dict) else True

    def _path_value(value):
        if not value:
            return ""
        if not redact_paths:
            return str(value)
        try:
            return os.path.basename(str(value))
        except Exception:
            return str(value)

    add({"Category": "Run", "Key": "Timestamp", "Value": datetime.now().isoformat(timespec="seconds")})
    add({"Category": "Run", "Key": "Output subfolder", "Value": _path_value(output_subfolder)})

    add({"Category": "Files", "Key": "Codon source", "Value": _path_value(clustering_results.get("codon_file", ""))})
    add({"Category": "Files", "Key": "GeneIDs source", "Value": _path_value(clustering_results.get("gene_file", ""))})
    add({"Category": "Files", "Key": "CDS FASTA", "Value": _path_value(fasta_path)})
    add({"Category": "Files", "Key": "Cluster file", "Value": _path_value(locustags_path)})

    add({"Category": "Dataset", "Key": "Strain prefix", "Value": clustering_results["strain_prefix"]})
    add({"Category": "Dataset", "Key": "N genes (codon table)", "Value": len(clustering_results["RowNames"])})
    add({"Category": "Dataset", "Key": "N genes (ordered)", "Value": len(clustering_results["ordered_genes"])})

    for k in ("usage_basis", "codon_set", "fasta_codon_range"):
        add({"Category": "Usage", "Key": k, "Value": SET.get(k)})
    if clustering_results.get("fasta_codon_range_label"):
        add({"Category": "Usage", "Key": "fasta_codon_range_label", "Value": clustering_results.get("fasta_codon_range_label")})
    if clustering_results.get("fasta_codon_range_empty_genes") is not None:
        add({"Category": "Usage", "Key": "fasta_codon_range_empty_genes", "Value": clustering_results.get("fasta_codon_range_empty_genes")})

    dimred_keys = ["dimred_method", "tsne_perplexity", "tsne_exaggeration", "tsne_learnrate",
                   "umap_neighbors", "umap_min_dist", "umap_metric", "pca_npcs"]
    for k in dimred_keys:
        add({"Category": "DimRed", "Key": k, "Value": SET.get(k)})

    clu_keys = ["cluster_method", "kmeans_k", "kmedoids_k", "spectral_k", "dbscan_eps", "dbscan_minpts"]
    for k in clu_keys:
        add({"Category": "Clustering", "Key": k, "Value": SET.get(k)})

    for k in ("apply_smoothing", "smooth_window_genes"):
        add({"Category": "Smoothing", "Key": k, "Value": SET.get(k)})
    for k in ("apply_binning", "bin_size_genes"):
        add({"Category": "Binning", "Key": k, "Value": SET.get(k)})

    if summary_df is not None and not summary_df.empty:
        srow = summary_df.iloc[0]
        for col in summary_df.columns:
            add({"Category": "Quantitative summary", "Key": col, "Value": srow[col]})

    return pd.DataFrame(rows, columns=["Category", "Key", "Value"])


# =========================================================
# BIG FUNCTION MOVED OUT: write_per_cluster_workbook
# =========================================================
_INVALID_SHEET_CHARS = re.compile(r"[:\\/?*\[\]]")


def _safe_sheet_name(name: str, existing: set) -> str:
    base = str(name) if name is not None else "cluster"
    base = base.strip() or "cluster"
    base = _INVALID_SHEET_CHARS.sub("_", base)[:31]

    candidate = base
    if candidate not in existing:
        existing.add(candidate)
        return candidate

    i = 2
    while True:
        suffix = f"_{i}"
        max_len = 31 - len(suffix)
        candidate = (base[:max_len] if max_len > 0 else "") + suffix
        if candidate not in existing:
            existing.add(candidate)
            return candidate
        i += 1


def _coord_colnames(dimred_method: str):
    m = (dimred_method or "").lower().strip()
    if m == "umap":
        return "UMAP-X", "UMAP-Y"
    if m == "tsne":
        return "tSNE-X", "tSNE-Y"
    if m == "pca":
        return "PC1", "PC2"
    return "Dim1", "Dim2"


def _clean_text(s):
    if not isinstance(s, str):
        return s
    return s.strip().strip('"').strip("'")


def write_per_cluster_workbook(out_path,
                               cluster_df: pd.DataFrame,
                               ordered_genes: np.ndarray,
                               gene_symbol_map: dict,
                               gene_desc_map: dict,
                               row_names: np.ndarray,
                               Y: np.ndarray,
                               dimred_method: str):
    xcol, ycol = _coord_colnames(dimred_method)

    coords = {}
    if Y is not None and len(row_names) == Y.shape[0]:
        for tag, xy in zip(row_names.tolist(), Y[:, :2].tolist()):
            coords[str(tag)] = (float(xy[0]), float(xy[1]))

    order_index = {str(tag): i for i, tag in enumerate(ordered_genes.tolist())}

    existing = set()
    with pd.ExcelWriter(out_path, engine="xlsxwriter") as writer:
        # All Genes
        all_rows = []
        for lt in row_names.tolist():
            gs = (gene_symbol_map.get(lt, "") or "").strip() if gene_symbol_map else ""
            desc = (gene_desc_map.get(lt, "") or "").strip() if gene_desc_map else ""
            x, y = coords.get(str(lt), (np.nan, np.nan))
            all_rows.append({"Locus tags": lt, "GeneSymbol": gs, xcol: x, ycol: y, "Description": desc})
        pd.DataFrame(all_rows, columns=["Locus tags", "GeneSymbol", xcol, ycol, "Description"]).to_excel(
            writer, sheet_name=_safe_sheet_name("All Genes", existing), index=False
        )

        # Clusters
        for col in cluster_df.columns:
            raw = cluster_df[col].replace({np.nan: ""}).astype(str).tolist()
            members, seen = [], set()
            for v in raw:
                v = _clean_text(v)
                if not v or v.lower() == "nan" or v in seen:
                    continue
                seen.add(v)
                members.append(v)

            members.sort(key=lambda t: order_index.get(str(t), 10**12))
            rows = []
            for lt in members:
                gs = (gene_symbol_map.get(lt, "") or "").strip() if gene_symbol_map else ""
                desc = (gene_desc_map.get(lt, "") or "").strip() if gene_desc_map else ""
                x, y = coords.get(str(lt), (np.nan, np.nan))
                rows.append({"Locus tags": lt, "GeneSymbol": gs, xcol: x, ycol: y, "Description": desc})

            pd.DataFrame(rows, columns=["Locus tags", "GeneSymbol", xcol, ycol, "Description"]).to_excel(
                writer, sheet_name=_safe_sheet_name(col, existing), index=False
            )


# =========================================================
# Codon-usage / tRNA-usage workbook utilities
# =========================================================

UNIFIED_CODON_DISPLAY_ORDER = [
    "Ala_GCA", "Ala_GCC", "Ala_GCG", "Ala_GCT",
    "Arg_AGA", "Arg_AGG", "Arg_CGA", "Arg_CGC", "Arg_CGG", "Arg_CGT",
    "Asn_AAC", "Asn_AAT",
    "Asp_GAC", "Asp_GAT",
    "Cys_TGC", "Cys_TGT",
    "Gln_CAA", "Gln_CAG",
    "Glu_GAA", "Glu_GAG",
    "Gly_GGA", "Gly_GGC", "Gly_GGG", "Gly_GGT",
    "His_CAC", "His_CAT",
    "Ile_ATA", "Ile_ATC", "Ile_ATT",
    "Leu_CTA", "Leu_CTC", "Leu_CTG", "Leu_CTT", "Leu_TTA", "Leu_TTG",
    "Lys_AAA", "Lys_AAG",
    "Met_AUG",
    "Phe_TTC", "Phe_TTT",
    "Pro_CCA", "Pro_CCC", "Pro_CCG", "Pro_CCT",
    "Ser_AGC", "Ser_AGT", "Ser_TCA", "Ser_TCC", "Ser_TCG", "Ser_TCT",
    "Thr_ACA", "Thr_ACC", "Thr_ACG", "Thr_ACT",
    "Trp_TGG",
    "Tyr_TAC", "Tyr_TAT",
    "Val_GTA", "Val_GTC", "Val_GTG", "Val_GTT",
]

_SINGLE_CODON_AAS = {"Met", "Trp"}
_INVALID_FILENAME_CHARS = re.compile(r'[<>:"/\\|?*\x00-\x1f]+')


def _split_label_to_aa_codon(label: str):
    if label is None:
        return "", ""
    s = str(label).strip()
    if not s or "_" not in s:
        return "", ""
    aa, cod = s.split("_", 1)
    return aa.strip(), cod.strip().upper()


def _split_label_to_aa_trna(label: str):
    return _split_label_to_aa_codon(label)


def _dedup_keep_order(values):
    seen = set()
    out = []
    for v in values:
        v = _clean_text(v) if isinstance(v, str) else v
        if not v or str(v).lower() == "nan":
            continue
        if v in seen:
            continue
        seen.add(v)
        out.append(v)
    return out


def _cluster_to_gene_list(cluster_df: pd.DataFrame) -> dict:
    out = {}
    for col in cluster_df.columns:
        raw = cluster_df[col].replace({np.nan: ""}).astype(str).tolist()
        out[str(col)] = _dedup_keep_order(raw)
    return out


def _display_label_to_internal(label: str) -> str:
    s = str(label).strip()
    if s == "Met_AUG":
        return "Met_ATG"
    aa, cod = _split_label_to_aa_codon(s)
    if not aa or not cod:
        return s
    return f"{aa}_{cod.replace('U', 'T')}"


def _internal_label_to_display(label: str) -> str:
    aa, cod = _split_label_to_aa_codon(label)
    if not aa or not cod:
        return str(label)
    if aa == "Met" and cod == "ATG":
        return "Met_AUG"
    return f"{aa}_{cod}"


def _ordered_internal_labels(df: pd.DataFrame, display_order=None) -> list:
    display_order = list(display_order or UNIFIED_CODON_DISPLAY_ORDER)
    wanted = [_display_label_to_internal(v) for v in display_order]
    present = [v for v in wanted if v in df.columns]
    extras = [str(v) for v in df.columns if str(v) not in present]
    return present + extras


def _genelevel_codon_sheet(df: pd.DataFrame, display_order=None, round_decimals: int = 6) -> pd.DataFrame:
    ordered_internal = _ordered_internal_labels(df, display_order=display_order)
    out = pd.DataFrame(index=df.index)
    for internal in ordered_internal:
        disp = _internal_label_to_display(internal)
        out[disp] = pd.to_numeric(df[internal], errors="coerce")
    if round_decimals is not None:
        out = out.round(int(round_decimals))
    out.index.name = "LocusTag"
    return out


def _genelevel_generic_sheet(df: pd.DataFrame, ordered_labels=None, round_decimals: int = 6, index_name: str = "LocusTag") -> pd.DataFrame:
    ordered = [str(v) for v in list(ordered_labels or []) if str(v) in df.columns]
    extras = [str(v) for v in df.columns if str(v) not in ordered]
    final_cols = ordered + extras
    out = pd.DataFrame(index=df.index)
    for lab in final_cols:
        out[str(lab)] = pd.to_numeric(df[lab], errors="coerce")
    if round_decimals is not None:
        out = out.round(int(round_decimals))
    out.index.name = index_name
    return out


def _build_aa_to_cols(labels, splitter=_split_label_to_aa_codon) -> dict:
    aa_to_cols = {}
    for lab in list(labels):
        aa, _ = splitter(lab)
        if aa:
            aa_to_cols.setdefault(aa, []).append(str(lab))
    return aa_to_cols


def compute_acu_from_counts_df(count_df: pd.DataFrame) -> pd.DataFrame:
    df = count_df.copy()
    denom = pd.to_numeric(df.sum(axis=1), errors="coerce").replace(0, np.nan)
    return df.div(denom, axis=0)


def compute_relative_usage_from_counts_df(
    count_df: pd.DataFrame,
    splitter,
    single_family_value_if_absent: float = 1.0,
) -> pd.DataFrame:
    """
    Generic relative-usage computation from count tables.

    For multi-member synonymous families:
      value(feature) = count(feature) / sum(counts for the family)

    For single-member families:
      value(feature) = single_family_value_if_absent for every gene.

    If a multi-member family is absent from a gene, all values for that family
    remain NaN so they can be excluded from downstream means.
    """
    df = count_df.copy()
    labels = [str(v) for v in df.columns.tolist()]
    fam_to_cols = _build_aa_to_cols(labels, splitter=splitter)

    out = pd.DataFrame(index=df.index, columns=labels, dtype=float)
    for fam, cols in fam_to_cols.items():
        if len(cols) == 1:
            out[cols[0]] = float(single_family_value_if_absent)
            continue
        denom = pd.to_numeric(df[cols].sum(axis=1), errors="coerce").replace(0, np.nan)
        out[cols] = df[cols].div(denom, axis=0)
    return out


def compute_rcu_from_counts_df(
    count_df: pd.DataFrame,
    single_codon_value_if_absent: float = 1.0,
) -> pd.DataFrame:
    """
    Compute gene-level RCU from gene-level codon counts.

    Multi-codon synonymous families:
      - RCU = count(codon) / sum(counts for the synonymous family)
      - if the amino acid is absent from the gene, values remain NaN so they are
        excluded from downstream cluster averages.

    Single-codon families (Met/Trp):
      - RCU is fixed to 1.0 for every gene, including genes lacking that amino
        acid, exactly as requested.
    """
    return compute_relative_usage_from_counts_df(
        count_df=count_df,
        splitter=_split_label_to_aa_codon,
        single_family_value_if_absent=single_codon_value_if_absent,
    )


def compute_devz(relative_df: pd.DataFrame, baseline_genes: list, metric_name: str = "relative-usage") -> pd.DataFrame:
    """
    Compute gene-level z-scores relative to the whole-genome baseline.

    z = (value_gene - genome_mean_value) / genome_sd_value

    NaN input values remain NaN and are excluded from both baseline statistics
    and downstream cluster averages.
    """
    base = [g for g in baseline_genes if g in relative_df.index]
    if not base:
        raise ValueError(f"No baseline genes found in the provided {metric_name} table; cannot compute z-scores.")

    baseline = relative_df.loc[base]
    mu = baseline.mean(axis=0, skipna=True)
    sd = baseline.std(axis=0, skipna=True, ddof=0)

    sd2 = sd.replace(0.0, np.nan)
    z = (relative_df - mu) / sd2

    zero_var_cols = sd.index[sd == 0.0].tolist()
    if zero_var_cols:
        for c in zero_var_cols:
            mask = relative_df[c].notna()
            z.loc[mask, c] = 0.0
    return z


def compute_rcu_devz(rcu_df: pd.DataFrame, baseline_genes: list) -> pd.DataFrame:
    """
    Compute gene-level ZCU values versus the whole genome.

    ZCU = (RCU_gene - genome_mean_RCU) / genome_sd_RCU
    """
    return compute_devz(rcu_df, baseline_genes=baseline_genes, metric_name="RCU")


def _cluster_average_table(values_df: pd.DataFrame, cluster_to_genes: dict, feature_label: str = "Codon", round_decimals: int = 6) -> pd.DataFrame:
    rows = []
    for feature in values_df.columns.tolist():
        row = {feature_label: str(feature)}
        for cname, genes in cluster_to_genes.items():
            genes_present = [g for g in genes if g in values_df.index]
            if not genes_present:
                val = np.nan
            else:
                val = pd.to_numeric(values_df.loc[genes_present, feature], errors="coerce").mean(skipna=True)
            if round_decimals is not None and pd.notna(val):
                val = round(float(val), int(round_decimals))
            row[str(cname)] = val
        rows.append(row)
    out = pd.DataFrame(rows)
    return out.loc[:, [feature_label] + [c for c in out.columns if c != feature_label]]


def _prepend_whole_genome_cluster(cluster_to_genes: dict, ordered_genes: list, whole_genome_name: str = "Whole genome") -> dict:
    out = {str(whole_genome_name): _dedup_keep_order(list(ordered_genes))}
    for cname, genes in cluster_to_genes.items():
        out[str(cname)] = _dedup_keep_order(list(genes))
    return out


def _safe_filename_component(name: str, default: str = "cluster", max_len: int = 120) -> str:
    s = str(name or "").strip()
    s = _INVALID_FILENAME_CHARS.sub("_", s)
    s = re.sub(r"\s+", " ", s).strip().strip(".")
    if not s:
        s = default
    if len(s) > max_len:
        s = s[:max_len].rstrip()
    return s or default


def _unique_filepath(path: str) -> str:
    if not os.path.exists(path):
        return path
    root, ext = os.path.splitext(path)
    i = 2
    while True:
        candidate = f"{root}_{i}{ext}"
        if not os.path.exists(candidate):
            return candidate
        i += 1


def _subset_rows_preserve_order(df: pd.DataFrame, genes: list) -> pd.DataFrame:
    genes_present = [g for g in list(genes) if g in df.index]
    return df.loc[genes_present].copy()



def _normalize_triplet_rna(s):
    if pd.isna(s):
        return ""
    s = str(s or "").strip().upper().replace("T", "U")
    s = re.sub(r"\s+", "", s)
    if s in {"", "NAN", "NONE", "NA"}:
        return ""
    return s if len(s) == 3 and re.fullmatch(r"[ACGUN]+", s) else ""


def _split_pooled_anticodon_group(s):
    """
    Return one or more normalized RNA anticodons from a pooled label.

    Accepted examples:
      UCU
      UCU/CCU
      UCU+CCU
      UCU and CCU
      UCU, CCU
      UCU ; CCU
    """
    if pd.isna(s):
        return []
    raw = str(s or "").strip().upper().replace("T", "U")
    if raw in {"", "NAN", "NONE", "NA"}:
        return []
    raw = raw.replace(" AND ", "/").replace("&", "/").replace("+", "/").replace(";", "/").replace(",", "/")
    raw = re.sub(r"\s*/\s*", "/", raw)
    raw = re.sub(r"\s+", "", raw)
    parts = [p for p in raw.split("/") if p]
    out = []
    seen = set()
    for p in parts:
        trip = _normalize_triplet_rna(p)
        if trip and trip not in seen:
            seen.add(trip)
            out.append(trip)
    return out


def _normalize_anticodon_group_label(s):
    parts = _split_pooled_anticodon_group(s)
    if not parts:
        return ""
    return "/".join(parts)


def _normalize_aa_name(s):
    s = str(s or "").strip()
    if not s:
        return ""
    s = s.replace("_", " ").replace("-", " ")
    tokens = [tok for tok in s.split() if tok]
    if not tokens:
        return ""
    aa = tokens[0]
    if len(aa) <= 3:
        return aa[0].upper() + aa[1:].lower()
    return aa[:3][0].upper() + aa[:3][1:].lower()



def _split_optional_comma_list(value):
    if pd.isna(value):
        return []
    raw = str(value or "").strip()
    if raw.lower() in {"", "nan", "none", "na"}:
        return []
    out = []
    seen = set()
    for part in re.split(r"[,;]+", raw):
        item = " ".join(str(part or "").strip().split())
        if not item:
            continue
        key = item.lower()
        if key in seen:
            continue
        seen.add(key)
        out.append(item)
    return out


def _find_decoding_table_column(df, aliases, positional_index=None):
    """Find a decoding-table column by normalized header aliases, with positional fallback."""
    if df is None or df.empty:
        return None
    alias_norm = {re.sub(r"[^a-z0-9]+", "", str(a).lower()) for a in aliases}
    for col in df.columns:
        key = re.sub(r"[^a-z0-9]+", "", str(col).lower())
        if key in alias_norm:
            return col
    if positional_index is not None and positional_index < len(df.columns):
        return df.columns[positional_index]
    return None


def read_trna_decoding_table(path: str, sheet_name=None):
    """
    Read the unified CodonPipe decoding-strategy table.

    Required columns, in order or by matching header:
      1. AA
      2. Codon (5'-3')
      3. Anticodon (5'-3')

    Optional columns:
      4. tRNA molecules/cell (mean)
      5. tRNA molecules/cell (STD)
      6. tRNA modifications involved in decoding
      7. tRMEs involved in decoding

    The anticodon column may contain pooled decoding groups such as UGC/GGC,
    UAG/GAG/CAG, or UGG/CGG/GGG. Each row is treated as one decoding group with
    weight 1.0 for the corresponding codon. tRNA abundances are read from
    columns 4-5 of this same table; sheet 2 is no longer used.
    """
    if not path or not os.path.exists(path):
        raise FileNotFoundError(f"tRNA decoding table not found: {path}")

    if sheet_name is None or str(sheet_name).strip() == "":
        df = pd.read_excel(path, sheet_name=0)
    else:
        df = pd.read_excel(path, sheet_name=sheet_name)

    if df is None or df.empty:
        raise ValueError(f"tRNA decoding table is empty: {path}")
    if df.shape[1] < 3:
        raise ValueError(
            "The decoding-strategy table must contain the mandatory fields: "
            "AA, Codon (5'-3'), and Anticodon (5'-3'). Optional fields are "
            "tRNA molecules/cell mean, tRNA molecules/cell STD, tRNA modifications, and tRMEs."
        )

    aa_col = _find_decoding_table_column(df, ["AA", "Amino acid", "Amino-acid"], 0)
    codon_col = _find_decoding_table_column(df, ["Codon (5'-3')", "Codon", "Codon 5-3"], 1)
    anticodon_col = _find_decoding_table_column(df, ["Anticodon (5'-3')", "Anticodon", "Anticodon 5-3"], 2)
    mean_col = _find_decoding_table_column(df, ["tRNA molecules/cell (mean)", "tRNA molecules per cell mean", "tRNA abundance mean", "mean"], 3 if df.shape[1] >= 4 else None)
    std_col = _find_decoding_table_column(df, ["tRNA molecules/cell (STD)", "tRNA molecules per cell STD", "tRNA abundance STD", "STD", "SD"], 4 if df.shape[1] >= 5 else None)
    mod_col = _find_decoding_table_column(df, ["tRNA modifications involved in decoding", "tRNA modifications", "modifications"], 5 if df.shape[1] >= 6 else None)
    trme_col = _find_decoding_table_column(df, ["tRMEs involved in decoding", "tRMEs", "tRNA modification enzymes", "enzymes"], 6 if df.shape[1] >= 7 else None)

    ordered_trna_labels = []
    codon_to_decoders = {}
    decoding_rows = []
    label_components = {}
    abundance_records = []
    abundance_std_records = []
    modifications_by_label = {}
    trmes_by_label = {}
    modifications_by_codon = {}
    trmes_by_codon = {}

    for _, row in df.iterrows():
        aa = _normalize_aa_name(row.get(aa_col, "") if aa_col is not None else "")
        codon_rna = _normalize_triplet_rna(row.get(codon_col, "") if codon_col is not None else "")
        anticodon_group = _normalize_anticodon_group_label(row.get(anticodon_col, "") if anticodon_col is not None else "")

        if not aa or not codon_rna or not anticodon_group:
            continue

        label = f"{aa}_{anticodon_group}"
        if label not in ordered_trna_labels:
            ordered_trna_labels.append(label)
        if label not in label_components:
            aa_lab, anti_lab = _split_label_to_aa_trna(label)
            label_components[label] = [f"{aa_lab}_{x}" for x in _split_pooled_anticodon_group(anti_lab)] or [label]

        codon_internal = f"{aa}_{codon_rna.replace('U', 'T')}"
        codon_to_decoders[codon_internal] = [(label, 1.0)]

        mean_val = np.nan
        if mean_col is not None:
            mean_raw = row.get(mean_col, np.nan)
            if not pd.isna(mean_raw) and str(mean_raw).strip() != "":
                try:
                    mean_val = float(mean_raw)
                    abundance_records.append((label, mean_val))
                except Exception:
                    mean_val = np.nan

        std_val = np.nan
        if std_col is not None:
            std_raw = row.get(std_col, np.nan)
            if not pd.isna(std_raw) and str(std_raw).strip() != "":
                try:
                    std_val = float(std_raw)
                    abundance_std_records.append((label, std_val))
                except Exception:
                    std_val = np.nan

        mods = _split_optional_comma_list(row.get(mod_col, "") if mod_col is not None else "")
        trmes = _split_optional_comma_list(row.get(trme_col, "") if trme_col is not None else "")
        if mods:
            modifications_by_label.setdefault(label, [])
            modifications_by_codon[codon_internal] = mods
            for m in mods:
                if m not in modifications_by_label[label]:
                    modifications_by_label[label].append(m)
        if trmes:
            trmes_by_label.setdefault(label, [])
            trmes_by_codon[codon_internal] = trmes
            for e in trmes:
                if e not in trmes_by_label[label]:
                    trmes_by_label[label].append(e)

        decoding_rows.append({
            "AA": aa,
            "Codon": codon_rna,
            "Codon_internal": codon_internal,
            "Anticodon": anticodon_group,
            "tRNA_molecules_per_cell_mean": mean_val,
            "tRNA_molecules_per_cell_STD": std_val,
            "tRNA_modifications_involved_in_decoding": ", ".join(mods),
            "tRMEs_involved_in_decoding": ", ".join(trmes),
        })

    if not codon_to_decoders:
        raise ValueError(f"No valid codon-to-tRNA decoding rules were found in: {path}")

    meta_df = pd.DataFrame(decoding_rows)

    abundance_series = None
    if abundance_records:
        abundance_series = pd.Series(
            [v for _, v in abundance_records],
            index=[k for k, _ in abundance_records],
            dtype=float,
            name="Abundance",
        ).groupby(level=0).mean()

    abundance_std_series = None
    if abundance_std_records:
        abundance_std_series = pd.Series(
            [v for _, v in abundance_std_records],
            index=[k for k, _ in abundance_std_records],
            dtype=float,
            name="Abundance_STD",
        ).groupby(level=0).mean()

    return dict(
        table_df=meta_df,
        ordered_trna_labels=ordered_trna_labels,
        codon_to_decoders=codon_to_decoders,
        label_components=label_components,
        abundance_series=abundance_series,
        abundance_std_series=abundance_std_series,
        modifications_by_label=modifications_by_label,
        trmes_by_label=trmes_by_label,
        modifications_by_codon=modifications_by_codon,
        trmes_by_codon=trmes_by_codon,
    )


def read_trna_abundance_table(path: str, sheet_name=None):
    """
    Read tRNA abundance values from columns 4-5 of the unified decoding table.

    The same sheet used by :func:`read_trna_decoding_table` is used; the legacy
    second-sheet abundance format is intentionally no longer used. The returned
    series is indexed by the same pooled labels used in ATU/RTU/ZTU tables, e.g.
    ``Ala_UGC/GGC`` or ``Leu_UAG/GAG/CAG``.
    """
    try:
        rules = read_trna_decoding_table(path, sheet_name=sheet_name)
    except Exception:
        return None
    s = rules.get("abundance_series")
    if s is None or getattr(s, "empty", True):
        return None
    s = pd.to_numeric(s, errors="coerce").dropna()
    if s.empty:
        return None
    s.name = "Abundance"
    return s


def align_trna_abundance_to_usage_labels(abundance_series: pd.Series, usage_labels, trna_rules: dict | None = None) -> pd.Series:
    if abundance_series is None:
        return pd.Series(dtype=float)
    abundance_series = pd.to_numeric(abundance_series, errors="coerce")
    out = {}
    label_components = dict((trna_rules or {}).get("label_components", {}) or {})
    for lab in list(usage_labels or []):
        lab = str(lab)
        exact = abundance_series.get(lab, np.nan)
        if pd.notna(exact):
            out[lab] = float(exact)
            continue
        candidates = list(label_components.get(lab, []))
        if not candidates:
            aa, anti = _split_label_to_aa_trna(lab)
            candidates = [f"{aa}_{x}" for x in _split_pooled_anticodon_group(anti)] or [lab]
        vals = [float(abundance_series.get(c, np.nan)) for c in candidates if pd.notna(abundance_series.get(c, np.nan))]
        out[lab] = float(np.sum(vals)) if vals else np.nan
    return pd.Series(out, dtype=float, name="Aligned abundance")


def _read_summary_metric_sheet(path: str, sheet_name: str):
    try:
        df = pd.read_excel(path, sheet_name=sheet_name, skiprows=1)
    except Exception:
        return None
    if df is None or df.empty:
        return None
    return df


def _compute_trna_abundance_correlation_table(usage_cluster_df: pd.DataFrame, abundance_aligned: pd.Series, method: str = "pearson", include_aas=None, metric: str | None = None):
    if usage_cluster_df is None or usage_cluster_df.empty or abundance_aligned is None or abundance_aligned.empty:
        return None
    metric = str(metric or '').upper()
    feature_col = usage_cluster_df.columns[0]
    feature_labels = usage_cluster_df[feature_col].astype(str).tolist()
    cluster_cols = [c for c in usage_cluster_df.columns if c != feature_col]

    aa_order = []
    aa_to_labs = {}
    for lab in feature_labels:
        aa, _ = _split_label_to_aa_trna(lab)
        if not aa:
            continue
        if aa not in aa_to_labs:
            aa_order.append(aa)
            aa_to_labs[aa] = []
        aa_to_labs[aa].append(lab)

    if include_aas is not None:
        keep = {str(a) for a in include_aas}
        aa_order = [aa for aa in aa_order if aa in keep]

    rows = []
    for cluster in cluster_cols:
        if metric == 'ZTU' and str(cluster).strip().lower() == 'whole genome':
            continue
        row = {"Cluster": str(cluster)}
        series = pd.to_numeric(usage_cluster_df[cluster], errors="coerce")
        lab_to_usage = dict(zip(feature_labels, series))
        for aa in aa_order:
            labs = aa_to_labs.get(aa, [])
            x = []
            y = []
            for lab in labs:
                a = abundance_aligned.get(lab, np.nan)
                u = lab_to_usage.get(lab, np.nan)
                if pd.notna(a) and pd.notna(u):
                    x.append(float(a))
                    y.append(float(u))
            if len(x) < 2 or np.nanstd(x) == 0 or np.nanstd(y) == 0:
                row[aa] = np.nan
            else:
                row[aa] = float(np.corrcoef(np.asarray(x, dtype=float), np.asarray(y, dtype=float))[0, 1])
        rows.append(row)
    return pd.DataFrame(rows)




def _metric_relevant_trna_aas(usage_cluster_df: pd.DataFrame, metric: str):
    if usage_cluster_df is None or usage_cluster_df.empty:
        return []
    feature_col = usage_cluster_df.columns[0]
    feature_labels = usage_cluster_df[feature_col].astype(str).tolist()
    aa_to_labs = {}
    order = []
    for lab in feature_labels:
        aa, _anti = _split_label_to_aa_trna(lab)
        if not aa:
            continue
        aa_to_labs.setdefault(aa, []).append(lab)
        if aa not in order:
            order.append(aa)
    metric = str(metric or '').upper()
    if metric in {'RTU', 'ZTU'}:
        order = [aa for aa in order if len(aa_to_labs.get(aa, [])) > 1]
    return order


def _filter_trna_corr_df_for_plot(corr_df: pd.DataFrame, metric: str):
    if corr_df is None or corr_df.empty:
        return corr_df
    df = corr_df.copy()
    metric = str(metric or '').upper()
    if metric == 'ZTU' and 'Cluster' in df.columns:
        cluster_series = df['Cluster'].astype(str).str.strip().str.lower()
        df = df.loc[cluster_series != 'whole genome'].copy()
    return df

TRNA_CORR_SHEET_DESCRIPTIONS = {
    "ATU abundance correlation": (
        "Pearson correlation between aligned tRNA abundances and cluster-level mean ATU values, computed separately for each amino-acid family. "
        "For cluster k and amino acid A, r(k,A) = corr( abundance_A(t), mean_ATU_k(t) ) across tRNA groups t decoding A. "
        "If pooled tRNA labels such as UCU/CCU are used in the decoding table, abundances from sheet 2 are summed across the corresponding anticodons before correlation."
    ),
    "RTU abundance correlation": (
        "Pearson correlation between aligned tRNA abundances and cluster-level mean RTU values, computed separately for each amino-acid family. "
        "For cluster k and amino acid A, r(k,A) = corr( abundance_A(t), mean_RTU_k(t) ) across tRNA groups t decoding A. "
        "If pooled tRNA labels such as UCU/CCU are used in the decoding table, abundances from sheet 2 are summed across the corresponding anticodons before correlation."
    ),
    "ZTU abundance correlation": (
        "Pearson correlation between aligned tRNA abundances and cluster-level mean ZTU values, computed separately for each amino-acid family. "
        "For cluster k and amino acid A, r(k,A) = corr( abundance_A(t), mean_ZTU_k(t) ) across tRNA groups t decoding A. "
        "If pooled tRNA labels such as UCU/CCU are used in the decoding table, abundances from sheet 2 are summed across the corresponding anticodons before correlation."
    ),
}


def _write_trna_correlation_workbook(out_path: str, aligned_abundance: pd.Series, corr_tables: dict):
    with pd.ExcelWriter(out_path, engine="xlsxwriter") as writer:
        aligned_df = pd.DataFrame({"tRNA": aligned_abundance.index.astype(str), "Aligned abundance": pd.to_numeric(aligned_abundance, errors="coerce").values})
        _write_sheet_with_description(
            writer,
            "Aligned abundances",
            "Aligned tRNA abundances used for the correlation analyses. When pooled tRNA labels are used in the decoding table, abundances of the corresponding individual anticodons from sheet 2 are summed.",
            aligned_df,
            index=False,
        )
        for sheet_name, df in corr_tables.items():
            desc = TRNA_CORR_SHEET_DESCRIPTIONS.get(sheet_name, "Correlation table.")
            _write_sheet_with_description(writer, sheet_name, desc, df, index=False)
    return out_path


def _save_trna_corr_heatmap(corr_df: pd.DataFrame, out_path: str, title: str, show_fig: bool = False, dpi: int = 300,
                            xmin=None, xmax=None, ymin=None, ymax=None):
    if corr_df is None or corr_df.empty:
        return ""
    df = corr_df.copy()
    if "Cluster" in df.columns:
        df = df.set_index("Cluster")
    if df.empty:
        return ""
    vals = df.astype(float).to_numpy()
    masked = np.ma.masked_invalid(vals)

    fig_w = max(6.0, 0.9 * df.shape[1] + 2.5)
    fig_h = max(4.5, 0.45 * df.shape[0] + 1.8)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    cmap = cm.get_cmap("coolwarm").copy()
    try:
        cmap.set_bad(color="#d9d9d9")
    except Exception:
        pass
    im = ax.imshow(masked, aspect="auto", vmin=-1, vmax=1, cmap=cmap)
    ax.set_xticks(np.arange(df.shape[1]))
    ax.set_xticklabels(df.columns.tolist(), rotation=45, ha="right")
    ax.set_yticks(np.arange(df.shape[0]))
    ax.set_yticklabels(df.index.tolist())
    ax.set_xlabel("Amino acid")
    ax.set_ylabel("Cluster")
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Pearson r")
    _apply_optional_axis_limits(ax, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    fig.tight_layout()
    fig.savefig(out_path, dpi=int(dpi), bbox_inches="tight")
    if show_fig:
        try:
            plt.show()
        except Exception:
            pass
    else:
        plt.close(fig)
    return out_path


def _cluster_safe_name(name: str) -> str:
    s = str(name or '').strip() or 'cluster'
    s = re.sub(r'[\/:*?"<>|]+', '_', s)
    s = re.sub(r'\s+', ' ', s).strip().replace(' ', '_')
    return s[:120]


def _apply_optional_axis_limits(ax, xmin=None, xmax=None, ymin=None, ymax=None):
    try:
        if xmin not in (None, "") or xmax not in (None, ""):
            ax.set_xlim(left=(float(xmin) if xmin not in (None, "") else None), right=(float(xmax) if xmax not in (None, "") else None))
    except Exception:
        pass
    try:
        if ymin not in (None, "") or ymax not in (None, ""):
            ax.set_ylim(bottom=(float(ymin) if ymin not in (None, "") else None), top=(float(ymax) if ymax not in (None, "") else None))
    except Exception:
        pass


def _select_requested_clusters(cluster_cols, selection):
    cols = [str(c) for c in (cluster_cols or [])]
    sel = str(selection or '').strip()
    if not sel or sel.lower() == 'all':
        return cols
    tokens = [t.strip() for t in re.split(r"[;,]+", sel) if t.strip()]
    lower_map = {c.strip().lower(): c for c in cols}
    chosen = []
    for tok in tokens:
        if tok.isdigit():
            idx = int(tok) - 1
            if 0 <= idx < len(cols):
                chosen.append(cols[idx])
                continue
        mapped = lower_map.get(tok.lower())
        if mapped is not None:
            chosen.append(mapped)
    chosen = list(dict.fromkeys(chosen))
    return chosen if chosen else cols


def _save_trna_usage_abundance_scatter_panels(
    usage_cluster_df: pd.DataFrame,
    abundance_aligned: pd.Series,
    out_dir: str,
    metric: str = 'RTU',
    yscale: str = 'linear',
    show_fig: bool = False,
    dpi: int = 300,
    include_aas=None,
    xmin=None,
    xmax=None,
    ymin=None,
    ymax=None,
):
    if usage_cluster_df is None or usage_cluster_df.empty or abundance_aligned is None or abundance_aligned.empty:
        return []
    feature_col = usage_cluster_df.columns[0]
    feature_labels = usage_cluster_df[feature_col].astype(str).tolist()
    cluster_cols = [c for c in usage_cluster_df.columns if c != feature_col]

    aa_order = []
    aa_to_labs = {}
    for lab in feature_labels:
        aa, anti = _split_label_to_aa_trna(lab)
        if not aa:
            continue
        aa_to_labs.setdefault(aa, []).append(lab)
        if aa not in aa_order:
            aa_order.append(aa)
    if include_aas is not None:
        keep = {str(a) for a in include_aas}
        aa_order = [aa for aa in aa_order if aa in keep]

    metric = str(metric or 'RTU').upper()
    if metric in {'RTU', 'ZTU'}:
        aa_order = [aa for aa in aa_order if len(aa_to_labs.get(aa, [])) > 1]
    if not aa_order or not cluster_cols:
        return []

    scatter_dir = os.path.join(out_dir, f'tRNA abundance vs {metric} scatter panels')
    os.makedirs(scatter_dir, exist_ok=True)
    saved_paths = []
    n_panels = len(aa_order)
    ncols = min(4, max(1, int(np.ceil(np.sqrt(n_panels)))))
    nrows = int(np.ceil(n_panels / ncols))

    for cluster in cluster_cols:
        series = pd.to_numeric(usage_cluster_df[cluster], errors='coerce')
        lab_to_usage = dict(zip(feature_labels, series))
        fig, axes = plt.subplots(nrows, ncols, figsize=(4.0 * ncols, 3.4 * nrows), squeeze=False)
        axes_flat = axes.ravel()
        any_points = False
        for i, aa in enumerate(aa_order):
            ax = axes_flat[i]
            labs = aa_to_labs.get(aa, [])
            x = []
            y = []
            labels = []
            for lab in labs:
                a = abundance_aligned.get(lab, np.nan)
                u = lab_to_usage.get(lab, np.nan)
                if pd.notna(a) and pd.notna(u):
                    x.append(float(a))
                    y.append(float(u))
                    _, anti = _split_label_to_aa_trna(lab)
                    labels.append(anti)
            if x:
                any_points = True
                ax.scatter(x, y)
                for xi, yi, labtxt in zip(x, y, labels):
                    ax.annotate(str(labtxt), (xi, yi), xytext=(3, 3), textcoords='offset points', fontsize=7)
                try:
                    if str(yscale).lower() == 'log':
                        if np.all(np.asarray(y, dtype=float) > 0):
                            ax.set_yscale('log')
                        else:
                            ax.set_yscale('symlog', linthresh=1.0)
                    else:
                        ax.set_yscale('linear')
                except Exception:
                    ax.set_yscale('linear')
            ax.set_title(str(aa))
            ax.set_xlabel('tRNA abundance')
            ax.set_ylabel(metric)
            _apply_optional_axis_limits(ax, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
        for j in range(len(aa_order), len(axes_flat)):
            axes_flat[j].axis('off')
        fig.suptitle(f'{metric} vs tRNA abundance — {cluster}', y=0.995)
        fig.tight_layout(rect=[0, 0, 1, 0.97])
        if any_points:
            out_path = os.path.join(scatter_dir, f'{metric} vs abundance - {_cluster_safe_name(cluster)}.png')
            fig.savefig(out_path, dpi=int(dpi), bbox_inches='tight')
            saved_paths.append(out_path)
            if show_fig:
                try:
                    plt.show()
                except Exception:
                    pass
            else:
                plt.close(fig)
        else:
            plt.close(fig)
    return saved_paths


def write_trna_abundance_correlation_outputs(
    summary_workbook_path: str,
    trna_decoding_table_path: str,
    out_dir: str | None = None,
    trna_decoding_table_sheet=None,
    trna_abundance_sheet=None,
    write_workbook: bool = True,
    write_figures: bool = True,
    show_fig: bool = False,
    dpi: int = 300,
    heatmap_metric: str = "RTU",
    scatter_metric: str = "RTU",
    scatter_yscale: str = "linear",
    scatter_show_fig: bool = False,
    heatmap_clusters='all',
    scatter_clusters='all',
    heatmap_xmin=None,
    heatmap_xmax=None,
    heatmap_ymin=None,
    heatmap_ymax=None,
    scatter_xmin=None,
    scatter_xmax=None,
    scatter_ymin=None,
    scatter_ymax=None,
):
    """
    Create cluster-by-amino-acid correlation tables and optional heatmaps
    relating tRNA abundances (sheet 2, optional) to cluster-level ATU/RTU/ZTU.
    """
    if not summary_workbook_path or not os.path.exists(summary_workbook_path):
        raise FileNotFoundError(f"Summary workbook not found: {summary_workbook_path}")
    if not trna_decoding_table_path or not os.path.exists(trna_decoding_table_path):
        raise FileNotFoundError(f"tRNA workbook not found: {trna_decoding_table_path}")

    out_dir = out_dir or (os.path.dirname(summary_workbook_path) or ".")
    os.makedirs(out_dir, exist_ok=True)

    trna_rules = read_trna_decoding_table(trna_decoding_table_path, sheet_name=trna_decoding_table_sheet)
    abundance_series = read_trna_abundance_table(trna_decoding_table_path, sheet_name=(trna_decoding_table_sheet or trna_abundance_sheet))
    if abundance_series is None or abundance_series.empty:
        return {}

    corr_tables = {}
    usage_tables = {}
    aligned_abundance = None
    metric_to_sheet = {"ATU": "ATU abundance correlation", "RTU": "RTU abundance correlation", "ZTU": "ZTU abundance correlation"}
    for metric in ("ATU", "RTU", "ZTU"):
        df = _read_summary_metric_sheet(summary_workbook_path, metric)
        if df is None or df.empty:
            continue
        usage_tables[metric] = df
        feature_col = df.columns[0]
        usage_labels = df[feature_col].astype(str).tolist()
        aligned_abundance = align_trna_abundance_to_usage_labels(abundance_series, usage_labels, trna_rules=trna_rules)
        corr_df = _compute_trna_abundance_correlation_table(df, aligned_abundance, metric=metric)
        if corr_df is not None:
            corr_tables[metric_to_sheet[metric]] = corr_df

    if not corr_tables:
        return {}

    out = {}
    if write_workbook:
        wb_path = os.path.join(out_dir, "tRNA abundance correlations.xlsx")
        _write_trna_correlation_workbook(wb_path, aligned_abundance, corr_tables)
        out["workbook"] = wb_path

    if write_figures:
        heatmap_metric = str(heatmap_metric or 'RTU').upper()
        wanted_heatmaps = [heatmap_metric] if heatmap_metric in {'ATU', 'RTU', 'ZTU'} else ['ATU', 'RTU', 'ZTU']
        for metric in wanted_heatmaps:
            sheet_name = metric_to_sheet.get(metric)
            corr_df = corr_tables.get(sheet_name)
            usage_df = usage_tables.get(metric)
            if corr_df is None or usage_df is None:
                continue
            relevant_aas = _metric_relevant_trna_aas(usage_df, metric)
            plot_df = corr_df.copy()
            if relevant_aas:
                keep_cols = ['Cluster'] + [aa for aa in relevant_aas if aa in plot_df.columns]
                plot_df = plot_df.loc[:, [c for c in keep_cols if c in plot_df.columns]].copy()
            plot_df = _filter_trna_corr_df_for_plot(plot_df, metric)
            chosen_clusters = _select_requested_clusters(plot_df["Cluster"].tolist() if "Cluster" in plot_df.columns else [], heatmap_clusters)
            if "Cluster" in plot_df.columns and chosen_clusters:
                cluster_series = plot_df["Cluster"].astype(str)
                plot_df = plot_df.loc[cluster_series.isin(chosen_clusters)].copy()
                order_map = {str(c): i for i, c in enumerate(chosen_clusters)}
                plot_df["__order"] = plot_df["Cluster"].astype(str).map(lambda x: order_map.get(str(x), len(order_map)))
                plot_df = plot_df.sort_values("__order", kind="stable").drop(columns=["__order"])
            fig_path = os.path.join(out_dir, f"tRNA abundance correlations {metric}.png")
            saved = _save_trna_corr_heatmap(
                plot_df,
                fig_path,
                title=f"{metric} vs tRNA abundance correlations",
                show_fig=show_fig,
                dpi=dpi,
                xmin=heatmap_xmin,
                xmax=heatmap_xmax,
                ymin=heatmap_ymin,
                ymax=heatmap_ymax,
            )
            if saved:
                out[f"figure_{metric}"] = saved
        try:
            scatter_metric = str(scatter_metric or 'RTU').upper()
            usage_df = usage_tables.get(scatter_metric)
            if usage_df is not None and not usage_df.empty:
                relevant_aas = _metric_relevant_trna_aas(usage_df, scatter_metric)
                feature_col = usage_df.columns[0]
                selected_clusters = _select_requested_clusters([c for c in usage_df.columns if c != feature_col], scatter_clusters)
                usage_df_plot = usage_df.loc[:, [feature_col] + [c for c in selected_clusters if c in usage_df.columns]].copy()
                saved_panels = _save_trna_usage_abundance_scatter_panels(
                    usage_cluster_df=usage_df_plot,
                    abundance_aligned=aligned_abundance,
                    out_dir=out_dir,
                    metric=scatter_metric,
                    yscale=str(scatter_yscale or 'linear'),
                    show_fig=bool(scatter_show_fig),
                    dpi=int(dpi),
                    include_aas=relevant_aas,
                    xmin=scatter_xmin,
                    xmax=scatter_xmax,
                    ymin=scatter_ymin,
                    ymax=scatter_ymax,
                )
                for i, p in enumerate(saved_panels, start=1):
                    out[f"scatter_{scatter_metric}_{i}"] = p
        except Exception:
            pass
    return out


def compute_trna_counts_from_codon_counts_df(count_df: pd.DataFrame, trna_rules: dict) -> pd.DataFrame:
    ordered_labels = list(trna_rules.get("ordered_trna_labels", []))
    codon_to_decoders = dict(trna_rules.get("codon_to_decoders", {}))

    out = pd.DataFrame(0.0, index=count_df.index, columns=ordered_labels, dtype=float)
    for codon in count_df.columns:
        codon_s = str(codon)
        decoders = codon_to_decoders.get(codon_s, [])
        if not decoders:
            continue
        codon_counts = pd.to_numeric(count_df[codon_s], errors="coerce").fillna(0.0)
        for label, weight in decoders:
            if label not in out.columns:
                out[label] = 0.0
            out[label] = pd.to_numeric(out[label], errors="coerce").fillna(0.0) + codon_counts * float(weight)
    return out


def compute_rtu_from_counts_df(trna_count_df: pd.DataFrame, single_family_value_if_absent: float = 1.0) -> pd.DataFrame:
    """
    Compute gene-level RTU from gene-level tRNA counts.

    Multi-tRNA synonymous families (same amino acid):
      - RTU = count(tRNA) / sum(counts for the synonymous tRNA family)
      - if that amino-acid family is absent from the gene, values remain NaN.

    Single-tRNA families:
      - RTU is fixed to 1.0 for every gene, including genes lacking that amino acid.
    """
    return compute_relative_usage_from_counts_df(
        count_df=trna_count_df,
        splitter=_split_label_to_aa_trna,
        single_family_value_if_absent=single_family_value_if_absent,
    )


def compute_ztu(rtu_df: pd.DataFrame, baseline_genes: list) -> pd.DataFrame:
    return compute_devz(rtu_df, baseline_genes=baseline_genes, metric_name="RTU")


CODON_SHEET_DESCRIPTIONS = {
    "Codon counts": (
        "Codon counts: per-gene raw codon counts. "
        "For gene g and codon c, n(g,c) = number of occurrences of codon c in gene g."
    ),
    "ACU": (
        "ACU: per-gene absolute codon usage. "
        "ACU(g,c) = n(g,c) / sum_c' n(g,c')."
    ),
    "RCU": (
        "RCU: per-gene relative codon usage. "
        "RCU(g,c) = ACU(g,c) / sum_{k in Syn(c)} ACU(g,k) = n(g,c) / sum_{k in Syn(c)} n(g,k). "
        "For single-codon families (Met and Trp), RCU is fixed to 1. "
        "For absent multi-codon amino-acid families, RCU is left blank (not counted in averages)."
    ),
    "ZCU": (
        "ZCU: per-gene z-score of relative codon usage versus the whole-genome baseline. "
        "ZCU(g,c) = (RCU(g,c) - mu_genome(c)) / sigma_genome(c), where mu_genome(c) and sigma_genome(c) "
        "are the mean and standard deviation of gene-level RCU for codon c across the whole genome. "
        "If sigma_genome(c)=0, defined entries are set to 0."
    ),
}

TRNA_SHEET_DESCRIPTIONS = {
    "tRNA counts": (
        "tRNA counts: per-gene weighted tRNA decoding counts. "
        "For gene g and tRNA t, m(g,t) = sum_c n(g,c) * w(c,t), where n(g,c) is the codon count and w(c,t) "
        "is the decoding fraction assigned to tRNA t for codon c from the imported codon-anticodon table. "
        "Pooled tRNA labels such as UCU/CCU are allowed and are treated as one combined decoding group."
    ),
    "ATU": (
        "ATU: per-gene absolute tRNA usage. "
        "ATU(g,t) = m(g,t) / sum_u m(g,u)."
    ),
    "RTU": (
        "RTU: per-gene relative tRNA usage within synonymous tRNA families decoding the same amino acid. "
        "RTU(g,t) = ATU(g,t) / sum_{u in Syn_tRNA(t)} ATU(g,u) = m(g,t) / sum_{u in Syn_tRNA(t)} m(g,u). "
        "For single-tRNA families, RTU is fixed to 1. "
        "For absent multi-tRNA amino-acid families, RTU is left blank (not counted in averages)."
    ),
    "ZTU": (
        "ZTU: per-gene z-score of relative tRNA usage versus the whole-genome baseline. "
        "ZTU(g,t) = (RTU(g,t) - mu_genome(t)) / sigma_genome(t), where mu_genome(t) and sigma_genome(t) "
        "are the mean and standard deviation of gene-level RTU for tRNA t across the whole genome. "
        "If sigma_genome(t)=0, defined entries are set to 0."
    ),
}


def _write_sheet_with_description(writer, sheet_name: str, description: str, df: pd.DataFrame, index: bool = True):
    pd.DataFrame([[description]]).to_excel(
        writer,
        sheet_name=sheet_name,
        index=False,
        header=False,
        startrow=0,
        startcol=0,
    )
    df.to_excel(writer, sheet_name=sheet_name, index=index, startrow=1)


def _write_raw_cluster_workbook(
    out_path: str,
    counts_display: pd.DataFrame,
    acu_display: pd.DataFrame,
    rcu_display: pd.DataFrame,
    zcu_display: pd.DataFrame,
    trna_counts_display: pd.DataFrame | None = None,
    atu_display: pd.DataFrame | None = None,
    rtu_display: pd.DataFrame | None = None,
    ztu_display: pd.DataFrame | None = None,
) -> str:
    with pd.ExcelWriter(out_path, engine="xlsxwriter") as writer:
        _write_sheet_with_description(writer, "Codon counts", CODON_SHEET_DESCRIPTIONS["Codon counts"], counts_display, index=True)
        _write_sheet_with_description(writer, "ACU", CODON_SHEET_DESCRIPTIONS["ACU"], acu_display, index=True)
        _write_sheet_with_description(writer, "RCU", CODON_SHEET_DESCRIPTIONS["RCU"], rcu_display, index=True)
        _write_sheet_with_description(writer, "ZCU", CODON_SHEET_DESCRIPTIONS["ZCU"], zcu_display, index=True)
        if trna_counts_display is not None:
            _write_sheet_with_description(writer, "tRNA counts", TRNA_SHEET_DESCRIPTIONS["tRNA counts"], trna_counts_display, index=True)
            _write_sheet_with_description(writer, "ATU", TRNA_SHEET_DESCRIPTIONS["ATU"], atu_display, index=True)
            _write_sheet_with_description(writer, "RTU", TRNA_SHEET_DESCRIPTIONS["RTU"], rtu_display, index=True)
            _write_sheet_with_description(writer, "ZTU", TRNA_SHEET_DESCRIPTIONS["ZTU"], ztu_display, index=True)
    return out_path


def _write_summary_workbook(
    out_path: str,
    acu_cluster: pd.DataFrame,
    rcu_cluster: pd.DataFrame,
    zcu_cluster: pd.DataFrame,
    atu_cluster: pd.DataFrame | None = None,
    rtu_cluster: pd.DataFrame | None = None,
    ztu_cluster: pd.DataFrame | None = None,
) -> str:
    with pd.ExcelWriter(out_path, engine="xlsxwriter") as writer:
        _write_sheet_with_description(writer, "ACU", CODON_SHEET_DESCRIPTIONS["ACU"], acu_cluster, index=False)
        _write_sheet_with_description(writer, "RCU", CODON_SHEET_DESCRIPTIONS["RCU"], rcu_cluster, index=False)
        _write_sheet_with_description(writer, "ZCU", CODON_SHEET_DESCRIPTIONS["ZCU"], zcu_cluster, index=False)
        if atu_cluster is not None:
            _write_sheet_with_description(writer, "ATU", TRNA_SHEET_DESCRIPTIONS["ATU"], atu_cluster, index=False)
            _write_sheet_with_description(writer, "RTU", TRNA_SHEET_DESCRIPTIONS["RTU"], rtu_cluster, index=False)
            _write_sheet_with_description(writer, "ZTU", TRNA_SHEET_DESCRIPTIONS["ZTU"], ztu_cluster, index=False)
    return out_path


def write_codon_usage_by_cluster_workbook(
    out_path: str,
    ordered_genes: list,
    cluster_df: pd.DataFrame,
    count_df: pd.DataFrame,
    aa_df: pd.DataFrame | None = None,
    c_abs_df: pd.DataFrame | None = None,
    c_rel_df: pd.DataFrame | None = None,
    z_rcu_df: pd.DataFrame | None = None,
    round_decimals: int = 6,
    whole_genome_name: str = "Whole genome",
    raw_subdir_name: str = "Raw codon usage tables",
    write_raw_workbooks: bool = True,
    compute_trna_usage: bool = False,
    trna_decoding_table_path: str = "",
    trna_decoding_table_sheet=None,
) -> str:
    """
    Write codon-usage exports using the requested logic.

    Raw per-cluster workbooks (plus a Whole genome workbook) contain:
      1. Codon counts
      2. ACU
      3. RCU
      4. ZCU
      5. tRNA counts   (optional)
      6. ATU           (optional)
      7. RTU           (optional)
      8. ZTU           (optional)

    The summary workbook contains:
      - ACU
      - RCU
      - ZCU
      - ATU (optional)
      - RTU (optional)
      - ZTU (optional)
    """
    if count_df is None or getattr(count_df, "empty", True):
        raise ValueError("count_df is required and must contain gene-level codon counts.")

    cluster_to_genes = _prepend_whole_genome_cluster(
        _cluster_to_gene_list(cluster_df),
        ordered_genes=list(ordered_genes),
        whole_genome_name=whole_genome_name,
    )

    count_df = count_df.copy()

    acu_df = compute_acu_from_counts_df(count_df)
    rcu_df = compute_rcu_from_counts_df(count_df)
    zcu_df = compute_rcu_devz(rcu_df, baseline_genes=list(ordered_genes))

    counts_display = _genelevel_codon_sheet(count_df, round_decimals=None)
    acu_display = _genelevel_codon_sheet(acu_df, round_decimals=round_decimals)
    rcu_display = _genelevel_codon_sheet(rcu_df, round_decimals=round_decimals)
    zcu_display = _genelevel_codon_sheet(zcu_df, round_decimals=round_decimals)

    trna_counts_df = atu_df = rtu_df = ztu_df = None
    trna_counts_display = atu_display = rtu_display = ztu_display = None
    ordered_trna_labels = None
    if bool(compute_trna_usage):
        if not str(trna_decoding_table_path or "").strip():
            raise ValueError("compute_trna_usage=True but no tRNA decoding table path was provided.")
        trna_rules = read_trna_decoding_table(trna_decoding_table_path, sheet_name=trna_decoding_table_sheet)
        ordered_trna_labels = list(trna_rules.get("ordered_trna_labels", []))
        trna_counts_df = compute_trna_counts_from_codon_counts_df(count_df, trna_rules)
        atu_df = compute_acu_from_counts_df(trna_counts_df)
        rtu_df = compute_rtu_from_counts_df(trna_counts_df)
        ztu_df = compute_ztu(rtu_df, baseline_genes=list(ordered_genes))

        trna_counts_display = _genelevel_generic_sheet(trna_counts_df, ordered_labels=ordered_trna_labels, round_decimals=round_decimals)
        atu_display = _genelevel_generic_sheet(atu_df, ordered_labels=ordered_trna_labels, round_decimals=round_decimals)
        rtu_display = _genelevel_generic_sheet(rtu_df, ordered_labels=ordered_trna_labels, round_decimals=round_decimals)
        ztu_display = _genelevel_generic_sheet(ztu_df, ordered_labels=ordered_trna_labels, round_decimals=round_decimals)

    raw_dir = os.path.join(os.path.dirname(out_path) or ".", str(raw_subdir_name))
    if write_raw_workbooks:
        os.makedirs(raw_dir, exist_ok=True)
        used_paths = set()
        for cname, genes in cluster_to_genes.items():
            safe_name = _safe_filename_component(cname)
            raw_path = _unique_filepath(os.path.join(raw_dir, f"{safe_name}.xlsx"))
            while raw_path in used_paths:
                raw_path = _unique_filepath(raw_path)
            used_paths.add(raw_path)

            _write_raw_cluster_workbook(
                out_path=raw_path,
                counts_display=_subset_rows_preserve_order(counts_display, genes),
                acu_display=_subset_rows_preserve_order(acu_display, genes),
                rcu_display=_subset_rows_preserve_order(rcu_display, genes),
                zcu_display=_subset_rows_preserve_order(zcu_display, genes),
                trna_counts_display=_subset_rows_preserve_order(trna_counts_display, genes) if trna_counts_display is not None else None,
                atu_display=_subset_rows_preserve_order(atu_display, genes) if atu_display is not None else None,
                rtu_display=_subset_rows_preserve_order(rtu_display, genes) if rtu_display is not None else None,
                ztu_display=_subset_rows_preserve_order(ztu_display, genes) if ztu_display is not None else None,
            )

    acu_cluster = _cluster_average_table(acu_display, cluster_to_genes, feature_label="Codon", round_decimals=round_decimals)
    rcu_cluster = _cluster_average_table(rcu_display, cluster_to_genes, feature_label="Codon", round_decimals=round_decimals)
    zcu_cluster = _cluster_average_table(zcu_display, cluster_to_genes, feature_label="Codon", round_decimals=round_decimals)

    atu_cluster = rtu_cluster = ztu_cluster = None
    if trna_counts_display is not None:
        atu_cluster = _cluster_average_table(atu_display, cluster_to_genes, feature_label="tRNA", round_decimals=round_decimals)
        rtu_cluster = _cluster_average_table(rtu_display, cluster_to_genes, feature_label="tRNA", round_decimals=round_decimals)
        ztu_cluster = _cluster_average_table(ztu_display, cluster_to_genes, feature_label="tRNA", round_decimals=round_decimals)

    _write_summary_workbook(
        out_path=out_path,
        acu_cluster=acu_cluster,
        rcu_cluster=rcu_cluster,
        zcu_cluster=zcu_cluster,
        atu_cluster=atu_cluster,
        rtu_cluster=rtu_cluster,
        ztu_cluster=ztu_cluster,
    )
    return out_path
