"""codonpipe.fasta_metrics

FASTA parsing, identifier canonicalization, and per-gene sequence metrics.
"""

# codonpipe/fasta_metrics.py
import os
import re
import math
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

from tkinter import Tk
from tkinter.filedialog import askopenfilename


# -----------------------------
# FASTA parsing helpers
# -----------------------------
_LOCUS_PATTERNS = [
    re.compile(r"\[locus_tag=([^\]]+)\]"),
    re.compile(r"\blocus_tag=([^\s\]]+)"),
    re.compile(r"^>+\s*([\w.\-:]+)"),
]
_GENE_PATTERNS = [
    re.compile(r"\[gene=([^\]]+)\]"),
    re.compile(r"\bgene=([^\s\]]+)"),
]
_PRODUCT_PATTERNS = [
    re.compile(r"\[(?:product|protein)=([^\]]+)\]"),
    re.compile(r"\b(?:product|protein)=([^\[]+?)(?:\s*\[|$)"),
]

_PRIMARY_ID_RE = re.compile(r"^>\s*([^\s]+)")
_LOCUS_TAG_RE = re.compile(r"\[locus_tag=([^\]]+)\]")
_PROTEIN_ID_RE = re.compile(r"\[protein_id=([^\]]+)\]")


def _clean_text(s):
    if not isinstance(s, str):
        return s
    return s.strip().strip('"').strip("'")


def parse_fasta_header(header):
    if not header or not header.startswith(">"):
        return None, ""

    locus = None
    for pat in _LOCUS_PATTERNS:
        m = pat.search(header)
        if m:
            locus = _clean_text(m.group(1))
            break

    gene = ""
    for pat in _GENE_PATTERNS:
        m = pat.search(header)
        if m:
            gene = _clean_text(m.group(1))
            break

    if not gene:
        for pat in _PRODUCT_PATTERNS:
            m = pat.search(header)
            if m:
                gene = _clean_text(m.group(1))
                break

    return locus, gene


def fasta_records(fasta_path):
    header = None
    seq_chunks = []
    with open(fasta_path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line.rstrip("\n")
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, "".join(seq_chunks)


def extract_strand_from_header(header):
    loc_m = re.search(r"\[location=([^\]]+)\]", header)
    if loc_m:
        loc = loc_m.group(1)
        return (0, "-") if "complement(" in loc else (1, "+")
    if re.search(r"strand\s*=\s*-\b", header) or re.search(r"\(-\)", header):
        return 0, "-"
    if re.search(r"strand\s*=\s*\+\b", header) or re.search(r"\(\+\)", header):
        return 1, "+"
    return 1, "+"


def extract_primary_id_from_header(header: str) -> str:
    if not header:
        return ""
    m = _PRIMARY_ID_RE.search(header)
    if m:
        return _clean_text(m.group(1))
    return _clean_text(header.lstrip(">").strip().split()[0]) if header.startswith(">") else _clean_text(str(header).split()[0])


def extract_protein_id_from_header(header: str) -> str:
    if not header:
        return ""
    m = _PROTEIN_ID_RE.search(header)
    return _clean_text(m.group(1)) if m else ""


# -----------------------------
# Canonicalization utilities
# -----------------------------
def canonicalize_id(x, alias_map: dict):
    if x is None:
        return ""
    s = str(x).strip()
    if not s or s.lower() == "nan":
        return ""

    m = _LOCUS_TAG_RE.search(s)
    if m:
        return _clean_text(m.group(1))

    if s.startswith(">"):
        s = s[1:].strip()
    if " " in s:
        s = s.split()[0].strip()

    return alias_map.get(s, s)


def canonicalize_generic_map(value_map: dict, alias_map: dict) -> dict:
    if not value_map:
        return {}
    out = {}
    for k, v in value_map.items():
        kk = canonicalize_id(k, alias_map)
        if not kk:
            continue
        vv = "" if v is None else str(v).strip()
        if kk not in out or (not out[kk] and vv):
            out[kk] = vv
    return out


def canonicalize_cluster_df(cluster_df: pd.DataFrame, alias_map: dict) -> pd.DataFrame:
    df = cluster_df.copy()
    for col in df.columns:
        df[col] = df[col].replace({np.nan: ""}).astype(str).map(
            lambda v: canonicalize_id(v, alias_map) if (v and str(v).lower() != "nan") else ""
        )
    return df


def enforce_unique_after_canon(orig_ids, canon_ids):
    seen = set()
    out = []
    for o, c in zip(list(orig_ids), list(canon_ids)):
        key = c if (c and c not in seen) else str(o)
        if key in seen or not key:
            base = key if key else str(o)
            i = 2
            key2 = base
            while key2 in seen or not key2:
                key2 = f"{base}__dup{i}"
                i += 1
            key = key2
        out.append(key)
        seen.add(key)
    return np.array(out, dtype=object)




# -----------------------------
# Codon table generation from CDS FASTA
# -----------------------------

# Standard genetic code (DNA alphabet)
_CODON2AA_1 = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"STOP","TAG":"STOP",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"STOP","TGG":"W",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

_AA1_TO_AA3 = {
    "A":"Ala","R":"Arg","N":"Asn","D":"Asp","C":"Cys",
    "Q":"Gln","E":"Glu","G":"Gly","H":"His","I":"Ile",
    "L":"Leu","K":"Lys","M":"Met","F":"Phe","P":"Pro",
    "S":"Ser","T":"Thr","W":"Trp","Y":"Tyr","V":"Val",
    "STOP":"STOP",
}

_STOP_CODONS = {"TAA", "TAG", "TGA"}

# Lexicographic codon order: AAA, AAC, ..., TTT
_BASE_ORDER = "ACGT"
_ALL_CODONS_64 = [a + b + c for a in _BASE_ORDER for b in _BASE_ORDER for c in _BASE_ORDER]

_RE_KV = re.compile(r"\[([A-Za-z0-9_\/\-\s]+)=([^\]]+)\]")
_RE_GENEID = re.compile(r"GeneID:(\d+)")
_RE_UNIPROT = re.compile(r"UniProtKB\/Swiss-Prot:([A-Z0-9]+)")


def _fasta_iter_simple(path: str):
    header = None
    seq_chunks = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, "".join(seq_chunks)


def _parse_header_kv(header: str) -> dict:
    out = dict(
        primary_id="",
        locus_tag="NA",
        old_locus_tag="NA",
        gene="NA",
        product="NA",
        protein_id="NA",
        entrez="NA",
        uniprot="NA",
    )
    if not header:
        return out

    out["primary_id"] = header.split()[0].strip()

    kv = {}
    for m in _RE_KV.finditer(header):
        k = m.group(1).strip()
        v = m.group(2).strip()
        kv[k] = v

    if "locus_tag" in kv:
        out["locus_tag"] = kv["locus_tag"]
    if "old_locus_tag" in kv:
        out["old_locus_tag"] = kv["old_locus_tag"]
    if "gene" in kv:
        out["gene"] = kv["gene"]

    if "protein" in kv:
        out["product"] = kv["protein"]
    elif "product" in kv:
        out["product"] = kv["product"]

    if "protein_id" in kv:
        out["protein_id"] = kv["protein_id"]

    m = _RE_GENEID.search(header)
    if m:
        out["entrez"] = m.group(1)

    m = _RE_UNIPROT.search(header)
    if m:
        out["uniprot"] = m.group(1)

    # If old_locus_tag absent but primary looks like SL1344_0001, keep it
    if out["old_locus_tag"] == "NA" and re.match(r".+_\d+$", out["primary_id"]):
        out["old_locus_tag"] = out["primary_id"]

    return out


def _normalize_organism_mode(mode: str) -> str:
    """The current pipeline is restricted to prokaryotic CDS FASTA inputs."""
    return "prokaryote"


def _valid_meta_value(v) -> bool:
    s = str(v or "").strip()
    return bool(s) and s not in {"NA", "None", "nan"}


def _choose_row_id(meta: dict, row_id_mode: str, organism_mode: str = "prokaryote") -> str:
    mode = str(row_id_mode or "primary").lower().strip()

    primary_id = meta.get("primary_id", "")
    protein_id = meta.get("protein_id", "")
    old_locus_tag = meta.get("old_locus_tag", "NA")
    locus_tag = meta.get("locus_tag", "NA")

    if mode == "primary":
        return primary_id
    if mode == "protein":
        return protein_id if _valid_meta_value(protein_id) else primary_id
    if mode == "old":
        return old_locus_tag if _valid_meta_value(old_locus_tag) else primary_id
    if mode == "locus":
        return locus_tag if _valid_meta_value(locus_tag) else primary_id
    if mode == "auto":
        for v in (locus_tag, old_locus_tag, protein_id, primary_id):
            if _valid_meta_value(v):
                return v
        return primary_id
    raise ValueError("row_id_mode must be one of: 'primary', 'protein', 'old', 'locus', 'auto'")


def _clean_seq_dna(seq: str) -> str:
    s = (seq or "").upper().replace("U", "T")
    s = re.sub(r"[^ACGT]", "N", s)
    return s


def parse_codon_range(codon_range):
    """Parse a user-facing codon range into 1-based inclusive coordinates.

    Accepted inputs:
      - "all", "", None: use the full CDS
      - "1-20" or "1:20": codons 1 through 20, inclusive
      - "20-200": codons 20 through 200, inclusive
      - "20-" / "20:end": codon 20 through the end
      - "20" or "first 20": first 20 codons

    Returns
    -------
    (start, end, label)
        start and end are 1-based inclusive integers, or None for an open side.
    """
    if codon_range is None:
        return None, None, "all"

    if isinstance(codon_range, (list, tuple)) and len(codon_range) >= 2:
        start = None if codon_range[0] in (None, "", "all") else int(codon_range[0])
        end = None if codon_range[1] in (None, "", "end", "all") else int(codon_range[1])
    else:
        s = str(codon_range).strip().lower()
        s = s.replace("–", "-").replace("—", "-")
        s = re.sub(r"\s+", " ", s)
        if s in {"", "all", "full", "whole", "entire", "entire gene", "entire cds"}:
            return None, None, "all"

        m = re.match(r"^first\s+(\d+)$", s)
        if m:
            start, end = 1, int(m.group(1))
        else:
            m = re.match(r"^(\d+)\s*(?:-|:|\.\.)\s*(\d+|end|all)?$", s)
            if m:
                start = int(m.group(1))
                raw_end = m.group(2)
                end = None if raw_end in (None, "", "end", "all") else int(raw_end)
            else:
                m = re.match(r"^(\d+)$", s)
                if m:
                    # Convenience shortcut: "20" means the first 20 codons.
                    start, end = 1, int(m.group(1))
                else:
                    raise ValueError(
                        "Invalid codon range. Use 'all', '1-20', '20-200', '20-end', or 'first 20'."
                    )

    if start is not None and start < 1:
        raise ValueError("Codon range start must be >= 1.")
    if end is not None and end < 1:
        raise ValueError("Codon range end must be >= 1.")
    if start is not None and end is not None and end < start:
        raise ValueError("Codon range end must be greater than or equal to the start.")

    if start is None and end is None:
        label = "all"
    elif start is None:
        label = f"codons 1-{end}"
    elif end is None:
        label = f"codons {start}-end"
    else:
        label = f"codons {start}-{end}"
    return start, end, label


def _slice_seq_by_codon_range(seq: str, codon_range="all", trim_to_multiple_of_3: bool = True):
    """Return the CDS subsequence corresponding to a 1-based inclusive codon range."""
    s = _clean_seq_dna(seq)
    if trim_to_multiple_of_3:
        s = s[:len(s) - (len(s) % 3)]

    full_codons = len(s) // 3
    start, end, label = parse_codon_range(codon_range)
    if start is None and end is None:
        return s, full_codons, full_codons, label

    start_i = 0 if start is None else max(0, int(start) - 1)
    # End is inclusive for users but exclusive for slicing in codon units.
    end_i = full_codons if end is None else min(int(end), full_codons)

    if start_i >= full_codons or end_i <= start_i:
        return "", full_codons, 0, label

    sliced = s[start_i * 3:end_i * 3]
    return sliced, full_codons, len(sliced) // 3, label


def _codon_counts(seq: str, codon_list, trim_to_multiple_of_3: bool = True, codon_range="all"):
    s, full_codons, selected_codons, range_label = _slice_seq_by_codon_range(
        seq, codon_range=codon_range, trim_to_multiple_of_3=trim_to_multiple_of_3
    )
    counts = {c: 0 for c in codon_list}
    total_valid = 0
    for i in range(0, len(s), 3):
        c = s[i:i+3]
        if c in counts:
            counts[c] += 1
            total_valid += 1
    return counts, total_valid, full_codons, selected_codons, range_label, s


def compute_codon_usage_tables_from_cds_fasta(
    fasta_path: str,
    row_id_mode: str = "primary",
    trim_to_multiple_of_3: bool = True,
    include_stops: bool = True,
    keep_first_duplicate: bool = True,
    organism_mode: str = "prokaryote",
    codon_range="all",
):
    """Build per-gene codon usage tables from a CDS FASTA.

    Returns
    -------
    codon_freq_df : pandas.DataFrame
        Per-gene codon frequencies (rows=genes, cols=AA3_CODON) in AAA..TTT order.
    geneids_df : pandas.DataFrame
        Parsed identifiers/annotations from headers.
    codon_abs_df : pandas.DataFrame
        Per-gene absolute codon counts (same columns/order as codon_freq_df).
    gc_percent : float
        GC% of concatenated unique selected CDS regions.

    Notes
    -----
    codon_range is interpreted as 1-based inclusive codon coordinates. For example,
    "1-20" uses the first 20 codons and "20-200" uses codons 20 through 200.
    If the requested end extends beyond a CDS, the end is clipped to that CDS length.
    """
    if not fasta_path or (not os.path.isfile(fasta_path)):
        raise FileNotFoundError(f"FASTA not found: {fasta_path}")

    _range_start, _range_end, codon_range_label = parse_codon_range(codon_range)

    codon_list = list(_ALL_CODONS_64)
    if not include_stops:
        codon_list = [c for c in codon_list if c not in _STOP_CODONS]

    freq_rows = {}
    abs_rows = {}
    meta_rows = []
    seq_by_id = {}

    for header, seq in _fasta_iter_simple(fasta_path):
        meta = _parse_header_kv(header)
        row_id = _choose_row_id(meta, row_id_mode=row_id_mode, organism_mode=organism_mode)
        if not row_id:
            continue

        if keep_first_duplicate and row_id in freq_rows:
            continue

        counts, total_valid, full_codons, selected_codons, range_label, selected_seq = _codon_counts(
            seq, codon_list=codon_list, trim_to_multiple_of_3=trim_to_multiple_of_3, codon_range=codon_range
        )
        if total_valid == 0:
            freqs = {c: 0.0 for c in codon_list}
        else:
            freqs = {c: counts[c] / float(total_valid) for c in codon_list}

        freq_rows[row_id] = freqs
        abs_rows[row_id] = counts
        seq_by_id[row_id] = selected_seq

        meta_rows.append(dict(
            LocusTag=row_id,
            GeneSymbol=meta.get("gene", "NA"),
            EntrezGeneID=meta.get("entrez", "NA"),
            ProteinDescription=meta.get("product", "NA"),
            RefSeqProteinID=meta.get("protein_id", "NA"),
            UniProtID=meta.get("uniprot", "NA"),
            RefSeq_LocusTag_RS=meta.get("locus_tag", "NA"),
            Old_LocusTag=meta.get("old_locus_tag", "NA"),
            PrimaryID=meta.get("primary_id", ""),
            CodonRange=codon_range_label,
            FullCDSCodons=full_codons,
            SelectedCodonsCounted=selected_codons,
        ))

    if not freq_rows:
        raise ValueError("No CDS records could be parsed from the FASTA (empty file or invalid format).")

    freq_df = pd.DataFrame.from_dict(freq_rows, orient="index")[codon_list]
    abs_df  = pd.DataFrame.from_dict(abs_rows,  orient="index")[codon_list]

    new_cols = []
    for codon in codon_list:
        aa1 = _CODON2AA_1.get(codon, "X")
        aa3 = _AA1_TO_AA3.get(aa1, "Xxx")
        new_cols.append(f"{aa3}_{codon}")
    freq_df.columns = new_cols
    abs_df.columns = new_cols

    geneids_df = pd.DataFrame(meta_rows)

    concat = "".join(seq_by_id.values()).upper()
    gc = concat.count("G") + concat.count("C")
    gc_pct = 100.0 * gc / max(1, len(concat))

    return freq_df, geneids_df, abs_df, gc_pct


# -----------------------------
# Codon group metrics
# -----------------------------
CODON_GROUPS = {
    "sensitive": {
        "ACC", "ATC", "CCG", "ATT", "CTC", "GGC", "CTA", "CCA",
        "CGC", "CGA", "CGT", "GCC", "GTC", "CTT", "TCC", "CAA", "GTG", "GGT"
    },
    "insensitive": {
        "AGA", "AGG", "CGG", "GGA", "GGG", "ATA", "TTA", "TTG", "AGC", "AGT", "TCG", "ACG"
    },
    "MiaA": {
        "TGC", "TGT", "TTA", "TTG", "TTC", "TTT", "TCC", "TCT", "TCA", "TCG",
        "TGG", "TAC", "TAT"
    },
    "MnmEG": {
        "AGG", "AGA", "CAA", "CAG", "GAA", "GAG", "GGA", "GGG", "TTA", "TTG",
        "AAA", "AAG"
    },
    "Tgt": {
        "AAC", "AAT", "GAC", "GAT", "CAC", "CAT", "TAC", "TAT"
    },
    "regulatory": {
        "AGA", "AGG", "TTA", "TTG"
    },
    "regbis": {
        "AGA", "AGG", "TTA", "TTG", "GGA", "GGG"
    },
}

THRESHOLD_FLAGS = [
    ("perc_insensitive", 17, "bin_insensitive_gt17"),
    ("perc_MiaA",       20.0, "bin_MiaA_gt20"),
    ("perc_MnmEG",      25.0, "bin_MnmEG_gt25"),
    ("perc_Tgt",        17.0, "bin_Tgt_gt17"),
    ("perc_regulatory",  5.0, "bin_regulatory_gt5"),
    ("perc_regbis",      8.0, "bin_regbis_gt8"),
]


def gc_fraction(seq):
    if not seq:
        return math.nan
    s = seq.upper().replace("U", "T")
    a = s.count("A")
    t = s.count("T")
    g = s.count("G")
    c = s.count("C")
    denom = a + t + g + c
    if denom == 0:
        return math.nan
    return (g + c) / denom


def gc_flags(gc_pct):
    if gc_pct is None or (isinstance(gc_pct, float) and math.isnan(gc_pct)):
        return {"gc_lt_48": 0, "gc_lt_45": 0, "gc_lt_40": 0, "gc_lt_35": 0}
    return {
        "gc_lt_48": 1 if gc_pct < 48.0 else 0,
        "gc_lt_45": 1 if gc_pct < 45.0 else 0,
        "gc_lt_40": 1 if gc_pct < 40.0 else 0,
        "gc_lt_35": 1 if gc_pct < 35.0 else 0,
    }


def _valid_codon(c):
    if len(c) != 3:
        return False
    for ch in c:
        if ch not in ("A", "C", "G", "T"):
            return False
    return True


def _aa_family_for_codon(codon):
    """Return the amino-acid family used for denominator normalization."""
    c = str(codon or "").upper().replace("U", "T")
    return _CODON2AA_1.get(c, "")


def _normalize_codon_set(codons):
    out = []
    seen = set()
    for c in list(codons or []):
        cc = str(c or "").strip().upper().replace("U", "T")
        if _valid_codon(cc) and cc not in seen and _aa_family_for_codon(cc) not in {"", "STOP"}:
            seen.add(cc)
            out.append(cc)
    return out


def _codon_group_denominator_codons(group_codons):
    """Codons used as denominator for a codon-defined metric.

    The denominator is all valid synonymous codons belonging to the amino-acid
    families represented in ``group_codons``. This avoids diluting a metric by
    unrelated amino acids for which the selected codon group is irrelevant.
    """
    group_codons = _normalize_codon_set(group_codons)
    aa_families = {aa for aa in (_aa_family_for_codon(c) for c in group_codons) if aa and aa != "STOP"}
    return {c for c, aa in _CODON2AA_1.items() if aa in aa_families and aa != "STOP"}


def codon_group_percentage(seq, group_codons):
    """Percentage of selected codons normalized by their synonymous families.

    Example: if the group contains Leu-TTA/TGG-like codons, the numerator is the
    count of the selected codons, while the denominator is the count of all Leu
    codons in that CDS/window. Codons from unrelated amino acids are ignored.
    """
    selected = set(_normalize_codon_set(group_codons))
    denominator_codons = _codon_group_denominator_codons(selected)
    if not selected or not denominator_codons or not seq:
        return math.nan, 0, 0

    s = str(seq or "").upper().replace("U", "T")
    L = len(s) - (len(s) % 3)
    numerator = 0
    denominator = 0
    for i in range(0, L, 3):
        c = s[i:i+3]
        if not _valid_codon(c):
            continue
        if c in denominator_codons:
            denominator += 1
            if c in selected:
                numerator += 1
    pct = math.nan if denominator <= 0 else 100.0 * numerator / float(denominator)
    return pct, numerator, denominator


def compute_codon_percentages(seq, groups=None):
    """Compute codon-group percentages for one CDS/window.

    Percentages are normalized by the synonymous codon families represented in
    each group, not by the total number of codons in the gene. This is important
    for groups such as "insensitive", MiaA, MnmEG or Tgt codons: amino acids that
    cannot contribute to the metric should not affect the denominator.
    """
    if groups is None:
        groups = CODON_GROUPS

    if not seq:
        return 0, {f"perc_{name}": math.nan for name in groups.keys()}

    total_valid = 0
    s = str(seq or "").upper().replace("U", "T")
    L = len(s) - (len(s) % 3)
    for i in range(0, L, 3):
        if _valid_codon(s[i:i+3]):
            total_valid += 1

    pcts = {}
    for name, codon_set in groups.items():
        pct, _num, _den = codon_group_percentage(s, codon_set)
        pcts[f"perc_{name}"] = pct

    return total_valid, pcts


def perc_threshold_flags(pcts):
    flags = {}
    for key, cutoff, outname in THRESHOLD_FLAGS:
        val = pcts.get(key, math.nan)
        if val is None or (isinstance(val, float) and math.isnan(val)):
            flags[outname] = 0
        else:
            flags[outname] = 1 if (val > cutoff) else 0
    return flags


FASTA_METRIC_GROUP_DEFAULTS = {
    "insensitive": {
        "label": "Insensitive genes",
        "metric_type": "codon_group",
        "codons": sorted(CODON_GROUPS["insensitive"]),
        "mode": "bins",
        "cutoff": 17.0,
        "top_n": 100,
        "direction": "high",
    },
    "MiaA": {
        "label": "MiaA dependent genes",
        "metric_type": "codon_group",
        "codons": sorted(CODON_GROUPS["MiaA"]),
        "mode": "bins",
        "cutoff": 20.0,
        "top_n": 100,
        "direction": "high",
    },
    "MnmEG": {
        "label": "MnmEG dependent genes",
        "metric_type": "codon_group",
        "codons": sorted(CODON_GROUPS["MnmEG"]),
        "mode": "bins",
        "cutoff": 25.0,
        "top_n": 100,
        "direction": "high",
    },
    "Tgt": {
        "label": "Tgt dependent genes",
        "metric_type": "codon_group",
        "codons": sorted(CODON_GROUPS["Tgt"]),
        "mode": "bins",
        "cutoff": 17.0,
        "top_n": 100,
        "direction": "high",
    },
    "AT_high": {
        "label": "AT high genes",
        "metric_type": "AT_percent",
        "codons": [],
        "mode": "bins",
        "cutoff": 55.0,
        "top_n": 100,
        "direction": "high",
    },
    "GC_low": {
        "label": "GC low genes",
        "metric_type": "GC_percent",
        "codons": [],
        "mode": "bins",
        "cutoff": 45.0,
        "top_n": 100,
        "direction": "low",
    },
}


def get_fasta_metric_group_defaults():
    """Return copy-safe default definitions for GUI/pipeline metric clusters."""
    out = {}
    for key, cfg in FASTA_METRIC_GROUP_DEFAULTS.items():
        item = dict(cfg)
        item["codons"] = list(cfg.get("codons", []))
        out[key] = item
    return out


def _metric_config_with_defaults(config):
    cfg = dict(config or {})
    key = str(cfg.get("key") or cfg.get("metric_key") or cfg.get("name") or "").strip()
    default = dict(FASTA_METRIC_GROUP_DEFAULTS.get(key, {}))
    default.setdefault("key", key)
    merged = default
    merged.update(cfg)
    merged["key"] = key or str(merged.get("label", "metric")).strip()
    merged["label"] = str(merged.get("label") or merged.get("key") or "Metric group").strip()
    merged["metric_type"] = str(merged.get("metric_type") or "codon_group").strip()
    merged["mode"] = str(merged.get("mode") or "bins").strip().lower().replace(" ", "_")
    if merged["mode"] not in {"bins", "top_hits"}:
        merged["mode"] = "bins"
    merged["direction"] = str(merged.get("direction") or "high").strip().lower()
    try:
        merged["cutoff"] = float(merged.get("cutoff", 0.0))
    except Exception:
        merged["cutoff"] = 0.0
    try:
        merged["top_n"] = max(1, int(float(merged.get("top_n", 100))))
    except Exception:
        merged["top_n"] = 100
    merged["codons"] = _normalize_codon_set(merged.get("codons", []))
    return merged


def _at_percent(seq):
    s = str(seq or "").upper().replace("U", "T")
    a = s.count("A")
    t = s.count("T")
    g = s.count("G")
    c = s.count("C")
    denom = a + t + g + c
    return math.nan if denom <= 0 else 100.0 * (a + t) / float(denom)


def _gc_percent_from_seq(seq):
    frac = gc_fraction(seq)
    return math.nan if frac is None or (isinstance(frac, float) and math.isnan(frac)) else 100.0 * frac


def score_sequence_for_metric(seq, config):
    cfg = _metric_config_with_defaults(config)
    typ = str(cfg.get("metric_type") or "").strip()
    if typ == "codon_group":
        pct, _num, _den = codon_group_percentage(seq, cfg.get("codons", []))
        return pct
    if typ == "AT_percent":
        return _at_percent(seq)
    if typ == "GC_percent":
        return _gc_percent_from_seq(seq)
    return math.nan


def build_fasta_metric_cluster_df(
    fasta_path: str,
    metric_configs,
    row_id_mode: str = "primary",
    trim_to_multiple_of_3: bool = True,
    organism_mode: str = "prokaryote",
    codon_range="all",
):
    """Build one cluster column per selected FASTA-derived metric.

    ``mode='bins'`` uses the metric cutoff. For high-direction metrics, genes with
    score >= cutoff are retained; for low-direction metrics, score <= cutoff is
    retained. ``mode='top_hits'`` keeps the top N genes by score (or lowest N for
    low-direction metrics).
    """
    configs = [_metric_config_with_defaults(c) for c in list(metric_configs or [])]
    configs = [c for c in configs if c.get("label")]
    if not configs:
        return pd.DataFrame(), pd.DataFrame()
    if not fasta_path or (not os.path.isfile(fasta_path)):
        raise FileNotFoundError(f"FASTA not found: {fasta_path}")

    rows = []
    for header, seq in _fasta_iter_simple(fasta_path):
        meta = _parse_header_kv(header)
        row_id = _choose_row_id(meta, row_id_mode=row_id_mode, organism_mode=organism_mode)
        if not row_id:
            continue
        selected_seq, full_codons, selected_codons, range_label = _slice_seq_by_codon_range(
            seq, codon_range=codon_range, trim_to_multiple_of_3=trim_to_multiple_of_3
        )
        rec = {
            "LocusTag": row_id,
            "GeneSymbol": meta.get("gene", "NA"),
            "ProteinDescription": meta.get("product", "NA"),
            "CodonRange": range_label,
            "FullCDSCodons": full_codons,
            "SelectedCodonsCounted": selected_codons,
        }
        for cfg in configs:
            rec[str(cfg["label"])] = score_sequence_for_metric(selected_seq, cfg)
        rows.append(rec)

    scores_df = pd.DataFrame(rows)
    if scores_df.empty:
        return pd.DataFrame(), scores_df

    clusters = {}
    for cfg in configs:
        label = str(cfg["label"]).strip()
        series = pd.to_numeric(scores_df[label], errors="coerce")
        valid = scores_df.loc[series.notna(), ["LocusTag", label]].copy()
        valid[label] = pd.to_numeric(valid[label], errors="coerce")
        direction = str(cfg.get("direction", "high")).lower()
        if cfg.get("mode") == "top_hits":
            ascending = direction == "low"
            chosen = valid.sort_values(label, ascending=ascending).head(int(cfg.get("top_n", 100)))["LocusTag"].tolist()
        else:
            cutoff = float(cfg.get("cutoff", 0.0))
            if direction == "low":
                chosen = valid.loc[valid[label] <= cutoff, "LocusTag"].tolist()
            else:
                chosen = valid.loc[valid[label] >= cutoff, "LocusTag"].tolist()
        clusters[label] = pd.Series(list(dict.fromkeys([str(x).strip() for x in chosen if str(x).strip()])))

    cluster_df = pd.DataFrame(clusters) if clusters else pd.DataFrame()
    return cluster_df.fillna(""), scores_df


def append_fasta_metric_clusters(cluster_df: pd.DataFrame, metric_cluster_df: pd.DataFrame) -> pd.DataFrame:
    """Append FASTA-derived metric clusters to an existing cluster table."""
    if metric_cluster_df is None or getattr(metric_cluster_df, "empty", True):
        return cluster_df.copy().fillna("") if cluster_df is not None else pd.DataFrame()
    base = cluster_df.copy().fillna("") if cluster_df is not None else pd.DataFrame()
    metric_cluster_df = metric_cluster_df.copy().fillna("")
    max_len = max(len(base), len(metric_cluster_df))
    base = base.reindex(range(max_len)).fillna("")
    metric_cluster_df = metric_cluster_df.reindex(range(max_len)).fillna("")
    existing_lower = {str(c).strip().lower(): str(c) for c in base.columns}
    for col in metric_cluster_df.columns:
        out_col = str(col).strip() or "FASTA metric group"
        if out_col.lower() in existing_lower:
            base_col = out_col
            i = 2
            while f"{base_col} ({i})".lower() in existing_lower:
                i += 1
            out_col = f"{base_col} ({i})"
        base[out_col] = metric_cluster_df[col].tolist()
        existing_lower[out_col.lower()] = out_col
    return base.fillna("")


# =========================================================
# BIG FUNCTION MOVED OUT: build_locus_index
# =========================================================
def build_locus_index(fasta_path, organism_mode="prokaryote", codon_range="all"):
    index = {}
    dup_counts = {}
    missing_locus_headers = 0

    lt_re = re.compile(r"\[locus_tag=([^\]]+)\]")
    gb_re = re.compile(r"\[gbkey=([^\]]+)\]")

    counts = Counter()
    gb_by_lt = defaultdict(Counter)

    alias_map = {}
    idmap_rows = []

    print("[INFO] Parsing FASTA and indexing by locus_tag…")
    for header, seq in fasta_records(fasta_path):
        m = lt_re.search(header)
        if m:
            lt = m.group(1)
            counts[lt] += 1
            g = gb_re.search(header)
            gb = g.group(1) if g else "NA"
            gb_by_lt[lt][gb] += 1

        primary_id = extract_primary_id_from_header(header)
        protein_id = extract_protein_id_from_header(header)
        locus, gene_name = parse_fasta_header(header)

        canonical = locus or protein_id or primary_id
        if not canonical:
            missing_locus_headers += 1
            continue

        if primary_id:
            alias_map[primary_id] = canonical
        if protein_id:
            alias_map[protein_id] = canonical
        if locus:
            alias_map[locus] = locus

        idmap_rows.append({
            "canonical_id": canonical,
            "locus_tag": locus,
            "primary_id": primary_id,
            "protein_id": protein_id,
        })

        selected_seq, _full_codons, selected_codons, _range_label = _slice_seq_by_codon_range(seq, codon_range=codon_range, trim_to_multiple_of_3=True)
        gc = gc_fraction(selected_seq)
        gc_pct = None if (gc is None or (isinstance(gc, float) and math.isnan(gc))) else 100.0 * gc
        strand_bin, strand_sym = extract_strand_from_header(header)
        total_codons_valid, pcts = compute_codon_percentages(selected_seq, CODON_GROUPS)

        rec = {
            "gene_name":            gene_name,
            "gc_fraction":          gc,
            "gc_percent":           gc_pct,
            "seq_length":           len(seq) if seq else 0,
            "selected_seq_length":  len(selected_seq) if selected_seq else 0,
            "selected_codons_counted": selected_codons,
            "total_codons_valid":   total_codons_valid,
            "strand_binary":        strand_bin,
            "strand_symbol":        strand_sym,
        }
        rec.update(gc_flags(gc_pct))
        rec.update(pcts)
        rec.update(perc_threshold_flags(pcts))

        if canonical in index:
            dup_counts[canonical] = dup_counts.get(canonical, 1) + 1
            continue
        index[canonical] = rec

    id_map_df = pd.DataFrame(idmap_rows).drop_duplicates()

    print(f"[INFO] Indexed {len(index)} unique entries from FASTA (canonical IDs).")
    if missing_locus_headers:
        print(f"[WARN] {missing_locus_headers} FASTA headers without detectable locus_tag/ID (ignored).")
    if dup_counts:
        print(f"[WARN] Detected {len(dup_counts)} canonical ID(s) with duplicate entries in FASTA (first kept).")

    print("Unique locus tags:", len(counts))
    print("How many locus tags appear twice:", sum(1 for v in counts.values() if v == 2))
    for lt, n in counts.most_common(10):
        if n > 1:
            print(lt, n, dict(gb_by_lt[lt]))

    return index, alias_map, id_map_df, missing_locus_headers, dup_counts


# -----------------------------
# FASTA selection helpers (kept)
# -----------------------------
def auto_find_fasta(folder, prefix_hint=""):
    exts = (".fa", ".fna", ".fasta", ".ffn", ".faa")
    candidates = []
    ph = (prefix_hint or "").lower().strip("_")

    for fname in os.listdir(folder):
        lower = fname.lower()
        if not lower.endswith(exts):
            continue

        score = 0
        if ph and ph in lower:
            score -= 10
        if "cds_from_genomic" in lower:
            score -= 6
        if "cds" in lower:
            score -= 3
        if lower.endswith(".ffn") or lower.endswith(".fna"):
            score -= 1
        if lower.endswith(".faa") or "protein" in lower:
            score += 5

        candidates.append((score, fname))

    if not candidates:
        return None

    candidates.sort()
    return os.path.join(folder, candidates[0][1])


def choose_fasta(initialdir):
    try:
        root = Tk()
        root.withdraw()
        path = askopenfilename(
            title="Select CDS FASTA file (e.g. *.fna / *.ffn)…",
            initialdir=initialdir,
            filetypes=[("FASTA files", "*.fa *.fna *.fasta *.ffn *.faa"),
                       ("All files", "*.*")]
        )
        root.destroy()
    except Exception:
        print("[NO GUI] Please provide path to CDS FASTA file:")
        path = input("> ").strip().strip('"').strip("'")

    if not path:
        raise SystemExit("No FASTA file selected; aborting.")
    return path


# =========================================================
# BIG FUNCTIONS MOVED OUT: metrics tables
# =========================================================
def build_metrics_table(ordered_tags, locus_index, gene_symbol_map):
    rows = []
    missing = []
    zero_flags = {name: 0 for _, _, name in THRESHOLD_FLAGS}

    for lt in ordered_tags:
        rec = locus_index.get(lt)
        gene_sym = (gene_symbol_map.get(lt, "") or "") if gene_symbol_map else ""

        if rec is None:
            missing.append(lt)
            row = {
                "locus_tag": lt,
                "gene_name": gene_sym,
                "gc_fraction": math.nan,
                "gc_percent": math.nan,
                "seq_length": math.nan,
                "total_codons_valid": 0,
                "strand_binary": math.nan,
                "strand_symbol": "",
                "gc_lt_48": 0,
                "gc_lt_45": 0,
                "gc_lt_40": 0,
                "gc_lt_35": 0,
                "perc_sensitive": math.nan,
                "perc_insensitive": math.nan,
                "perc_MiaA": math.nan,
                "perc_MnmEG": math.nan,
                "perc_Tgt": math.nan,
                "perc_regulatory": math.nan,
                "perc_regbis": math.nan,
            }
            row.update(zero_flags)
        else:
            row = {"locus_tag": lt, "gene_name": gene_sym if gene_sym else rec.get("gene_name", "")}
            row.update(rec)
            if gene_sym:
                row["gene_name"] = gene_sym

        rows.append(row)

    df = pd.DataFrame(rows)

    base_cols = [
        "locus_tag", "gene_name", "gc_fraction", "gc_percent", "seq_length",
        "total_codons_valid", "strand_binary", "strand_symbol",
        "gc_lt_48", "gc_lt_45", "gc_lt_40", "gc_lt_35",
        "perc_sensitive", "perc_insensitive", "perc_MiaA", "perc_MnmEG",
        "perc_Tgt", "perc_regulatory", "perc_regbis",
    ]
    flag_cols = [flag_name for _, _, flag_name in THRESHOLD_FLAGS]

    ordered_cols = [c for c in base_cols + flag_cols if c in df.columns]
    extra_cols = [c for c in df.columns if c not in ordered_cols]
    df = df[ordered_cols + extra_cols]

    return df, missing


def build_summary(metrics_df, total_requested, missing_count):
    matched_df = metrics_df[metrics_df["gc_percent"].notna()] if "gc_percent" in metrics_df.columns else metrics_df.iloc[0:0]

    summary = {
        "total_requested": [total_requested],
        "found_in_fasta": [total_requested - missing_count],
        "missing_in_fasta": [missing_count],
        "count_gc_lt_48": [int(matched_df.get("gc_lt_48", pd.Series([0])).sum())] if not matched_df.empty else [0],
        "count_gc_lt_45": [int(matched_df.get("gc_lt_45", pd.Series([0])).sum())] if not matched_df.empty else [0],
        "count_gc_lt_40": [int(matched_df.get("gc_lt_40", pd.Series([0])).sum())] if not matched_df.empty else [0],
        "count_gc_lt_35": [int(matched_df.get("gc_lt_35", pd.Series([0])).sum())] if not matched_df.empty else [0],
    }

    for _, _, flag_name in THRESHOLD_FLAGS:
        summary[f"count_{flag_name}"] = [int(metrics_df[flag_name].sum())] if flag_name in metrics_df.columns else [0]

    return pd.DataFrame(summary)


def build_quantitative_bis(metrics_df):
    if metrics_df is None or metrics_df.empty:
        return pd.DataFrame()
    if "locus_tag" not in metrics_df.columns:
        raise ValueError("metrics_df must contain a 'locus_tag' column.")

    locus_tags = metrics_df["locus_tag"].astype(str).tolist()
    out_cols = {}

    for col in metrics_df.columns:
        if col in ("locus_tag", "gene_name"):
            continue

        s_num = pd.to_numeric(metrics_df[col], errors="coerce")
        uniq = pd.unique(s_num.dropna())
        if len(uniq) == 0:
            continue
        if not np.all(np.isin(uniq, [0, 1])):
            continue

        hits = [lt for lt, v in zip(locus_tags, s_num.fillna(0).astype(int).tolist()) if v == 1]
        out_cols[str(col)] = hits

    if not out_cols:
        return pd.DataFrame()

    max_len = max(len(v) for v in out_cols.values())
    data = {}
    for k, v in out_cols.items():
        vv = list(v)
        if len(vv) < max_len:
            vv += [""] * (max_len - len(vv))
        data[k] = vv

    return pd.DataFrame(data)
