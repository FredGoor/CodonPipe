#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DAVID sliding-window enrichment scan on a reordered genome
+ full-genome DAVID term table
+ whole-genome term-derived cluster reconstruction seeded ONLY from the top-N
  sliding-window enrichment hits.

Behavior:
1. DAVID is tried first with locus tags using the dedicated LOCUS_TAG id type.
2. If needed, legacy locus-tag-like symbol modes are also attempted.
3. If DAVID still does not recognize / annotate those locus tags, the scan
   automatically falls back to Entrez Gene IDs when available.
4. If neither locus tags nor Entrez Gene IDs are usable, the caller receives
   a DAVID-recognition error so the existing popup behavior can be triggered.
"""

from __future__ import annotations

import logging
import os
import re
import time
import ssl
import urllib.error
import urllib.request
from collections import OrderedDict, defaultdict
from typing import Callable, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_DAVID_WSDL_URL = 'https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService?wsdl'
DEFAULT_DAVID_ENDPOINT = 'https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/'


def _suppress_suds_logging():
    noisy_names = [
        "suds",
        "suds.client",
        "suds.transport",
        "suds.transport.http",
        "suds.metrics",
        "suds.xsd",
        "suds.xsd.query",
        "suds.wsdl",
        "suds.umx",
        "suds.mx",
    ]
    for name in noisy_names:
        lg = logging.getLogger(name)
        lg.handlers = []
        lg.propagate = False
        lg.setLevel(logging.CRITICAL + 1)


def _clean_string(x) -> str:
    if pd.isna(x):
        return ""
    s = str(x).strip()
    return "" if s.lower() == "nan" else s


def _clean_entrez(x) -> str:
    s = _clean_string(x)
    if s in {"", "NA", "None"}:
        return ""
    if s.endswith(".0"):
        s = s[:-2]
    return s


def _clean_generic_david_id(x) -> str:
    s = _clean_string(x)
    if s in {"", "NA", "None"}:
        return ""
    if s.endswith(".0") and s[:-2].isdigit():
        s = s[:-2]
    return s


def _canonicalize_token(tok: str, alias_map: Optional[Dict[str, str]] = None) -> str:
    s = _clean_string(tok)
    if not s:
        return ""
    if alias_map:
        return str(alias_map.get(s, s))
    return s


def _first_nonempty(values: List[str]) -> str:
    for v in values:
        vv = _clean_string(v)
        if vv and vv not in {"NA", "None"}:
            return vv
    return ""


def _clean_term_label(term: str) -> str:
    term = _clean_string(term)
    if not term:
        return ""
    cleaned = re.split(r'[~:]', term, maxsplit=1)[-1].strip()
    return cleaned if cleaned else term


def _split_generic_ids_field(gene_ids_field: str) -> List[str]:
    vals = []
    for tok in re.split(r"[,\s;]+", _clean_string(gene_ids_field)):
        t = _clean_generic_david_id(tok)
        if t:
            vals.append(t)
    return list(OrderedDict.fromkeys(vals))


def _make_unique_headers(headers: List[str], max_len: int = 60) -> List[str]:
    counts = OrderedDict()
    out = []
    for h in headers:
        key = _clean_string(h) or "Unnamed"
        key = re.sub(r"\s+", " ", key).strip()
        if len(key) > max_len:
            key = key[:max_len - 3].rstrip() + "..."
        n = counts.get(key, 0) + 1
        counts[key] = n
        out.append(key if n == 1 else f"{key} ({n})")
    return out


def _collect_unique_ids(values, cleaner: Callable[[object], str]) -> List[str]:
    out = []
    seen = OrderedDict()
    for v in list(values):
        vv = cleaner(v)
        if vv and vv not in seen:
            seen[vv] = True
            out.append(vv)
    return out


def _looks_like_david_unrecognized_ids_error(exc) -> bool:
    msg = str(exc or "").strip().lower()
    if msg == "":
        return False
    return (
        ("server raised fault" in msg and "index: 0, size: 0" in msg)
        or ("index: 0, size: 0" in msg)
        or ("no functional annotation" in msg)
        or ("unrecognized" in msg and "david" in msg)
        or ("no david annotation rows" in msg)
    )




def _looks_like_ssl_certificate_error(exc) -> bool:
    """Return True for the common Windows/conda DAVID SSL CA-chain failure."""
    msg = str(exc or "").lower()
    if "certificate_verify_failed" in msg or "certificateverifyfailed" in msg:
        return True
    if "unable to get local issuer certificate" in msg:
        return True
    if "ssl" in msg and "certificate" in msg and "verify" in msg:
        return True
    cause = getattr(exc, "__cause__", None) or getattr(exc, "reason", None)
    if cause is not None and cause is not exc:
        return _looks_like_ssl_certificate_error(cause)
    return False


class _SSLContextHttpTransportBase:
    """Mixin factory for suds HttpTransport with a custom urllib SSL context.

    This keeps the normal suds-py3 transport behavior but installs an
    HTTPSHandler using the provided SSL context.  It is intentionally built at
    runtime because suds-py3 is an optional dependency of CodonPipe.
    """

    pass


def _make_suds_transport_with_ssl_context(context):
    from suds.transport.http import HttpTransport

    class SSLContextHttpTransport(HttpTransport):
        def __init__(self, ssl_context, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self._codonpipe_ssl_context = ssl_context

        def u2handlers(self):
            handlers = []
            try:
                handlers.extend(super().u2handlers())
            except Exception:
                pass
            handlers.append(urllib.request.HTTPSHandler(context=self._codonpipe_ssl_context))
            return handlers

    return SSLContextHttpTransport(context)


def _create_david_suds_client(Client, wsdl_url: str, allow_insecure_ssl_fallback: bool = True):
    """Create a suds client, with robust Windows/conda SSL fallbacks.

    The normal verified HTTPS path is always tried first.  If Python cannot find
    a local issuer certificate, CodonPipe tries the certifi CA bundle.  As a last
    resort, it can create a client with certificate verification disabled; this
    is limited to connecting to DAVID's public WSDL/service endpoint and is
    reported clearly in the console.
    """
    try:
        return Client(wsdl_url)
    except Exception as first_exc:
        if not _looks_like_ssl_certificate_error(first_exc):
            raise

        print("[WARN] DAVID HTTPS certificate verification failed with the default Python/conda CA store.")
        print("[INFO] Retrying DAVID connection with the certifi CA bundle...")

        # First secure fallback: certifi bundle, if installed.
        try:
            import certifi
            context = ssl.create_default_context(cafile=certifi.where())
            transport = _make_suds_transport_with_ssl_context(context)
            return Client(wsdl_url, transport=transport)
        except Exception as certifi_exc:
            if not _looks_like_ssl_certificate_error(certifi_exc):
                # If certifi fails for another reason, keep going only if the
                # final fallback is allowed; otherwise surface the real issue.
                if not allow_insecure_ssl_fallback:
                    raise
            print("[WARN] certifi-based DAVID HTTPS connection also failed.")

        # Last-resort compatibility fallback for Windows/conda installations
        # whose CA store is broken.  DAVID does not transmit private data beyond
        # the user's submitted gene IDs, but the user should still be informed.
        if allow_insecure_ssl_fallback:
            print("[WARN] Retrying DAVID connection with SSL certificate verification disabled.")
            print("[WARN] This fallback is only used because the local Python/conda CA store could not validate DAVID's certificate.")
            context = ssl._create_unverified_context()
            transport = _make_suds_transport_with_ssl_context(context)
            return Client(wsdl_url, transport=transport)

        raise RuntimeError(
            "DAVID HTTPS certificate verification failed. Install/update certifi "
            "or update the conda CA certificates package, then retry."
        ) from first_exc

class DAVIDEnrichmentRunner:
    def __init__(self, user_email: str,
                 wsdl_url: str = DEFAULT_DAVID_WSDL_URL,
                 endpoint: str = DEFAULT_DAVID_ENDPOINT,
                 allow_insecure_ssl_fallback: bool = True):
        if not user_email or not str(user_email).strip():
            raise ValueError("A DAVID-registered email address is required.")

        try:
            from suds.client import Client
        except Exception as e:
            raise ImportError(
                "The DAVID scan requires suds-py3. Install it with: pip install suds-py3"
            ) from e

        _suppress_suds_logging()
        print("[INFO] Connecting to DAVID web service...")
        self.client = _create_david_suds_client(Client, wsdl_url, allow_insecure_ssl_fallback=allow_insecure_ssl_fallback)
        self.client.wsdl.services[0].setlocation(endpoint)
        self.client.service.authenticate(str(user_email).strip())

    def analyze_gene_list(self, gene_list: List[str], output_prefix: str,
                          id_type: str = 'ENTREZ_GENE_ID',
                          overlap: int = 3,
                          initial_seed: int = 3,
                          final_seed: int = 3,
                          linkage: float = 0.5,
                          kappa: int = 50):
        input_ids = ",".join([str(g) for g in gene_list if _clean_string(g) != ""])
        if input_ids == "":
            return None
        list_type = 0
        list_name = output_prefix
        self.client.service.addList(input_ids, id_type, list_name, list_type)
        return self.client.service.getTermClusterReport(overlap, initial_seed, final_seed, linkage, kappa)

    def get_chart_report(self, gene_list: List[str], output_prefix: str,
                         id_type: str = 'ENTREZ_GENE_ID',
                         threshold: float = 1.0,
                         count: int = 1):
        input_ids = ",".join([str(g) for g in gene_list if _clean_string(g) != ""])
        if input_ids == "":
            return []
        list_type = 0
        list_name = output_prefix
        self.client.service.addList(input_ids, id_type, list_name, list_type)
        return self.client.service.getChartReport(float(threshold), int(count))


def save_cluster_report(term_clustering_report, filename: str):
    def _to_list(obj):
        if obj is None:
            return []
        if isinstance(obj, (list, tuple)):
            return list(obj)
        for attr in ("item", "items", "annotationClusters", "annotationCluster"):
            try:
                val = getattr(obj, attr)
            except Exception:
                val = None
            if val is not None:
                if isinstance(val, (list, tuple)):
                    return list(val)
                return [val]
        return [obj]

    def _safe_records(cluster_obj):
        try:
            recs = getattr(cluster_obj, "simpleChartRecords", None)
        except Exception:
            recs = None
        if recs is None:
            return []
        if isinstance(recs, (list, tuple)):
            return list(recs)
        try:
            item = getattr(recs, "item", None)
        except Exception:
            item = None
        if item is None:
            return [recs]
        if isinstance(item, (list, tuple)):
            return list(item)
        return [item]

    with open(filename, 'w', encoding='utf-8') as fout:
        clusters = _to_list(term_clustering_report)
        if not clusters:
            fout.write("No DAVID clustering result returned.\n")
            return
        wrote_any = False
        for i, cluster in enumerate(clusters, start=1):
            try:
                score = getattr(cluster, "score", "")
            except Exception:
                score = ""
            fout.write(f"Annotation Cluster {i}\tEnrichmentScore:{score}\n")
            fout.write("Category\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\tFDR\n")
            for record in _safe_records(cluster):
                row = [
                    str(getattr(record, "categoryName", "")),
                    str(getattr(record, "termName", "")),
                    str(getattr(record, "listHits", "")),
                    str(getattr(record, "percent", "")),
                    str(getattr(record, "ease", "")),
                    str(getattr(record, "geneIds", "")),
                    str(getattr(record, "listTotals", "")),
                    str(getattr(record, "popHits", "")),
                    str(getattr(record, "popTotals", "")),
                    str(getattr(record, "foldEnrichment", "")),
                    str(getattr(record, "bonferroni", "")),
                    str(getattr(record, "benjamini", "")),
                    str(getattr(record, "afdr", "")),
                ]
                fout.write('\t'.join(row) + '\n')
                wrote_any = True
        if not wrote_any:
            fout.write("No annotation clusters with chart records were returned.\n")


def extract_clusters(filename: str, max_clusters: int = 3):
    scores, pvals, sizes = [], [], []
    terms_list, current_terms = [], []
    inside_data = False
    current_first_pval = None
    current_cluster_size = None

    with open(filename, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('Annotation Cluster'):
                if current_terms:
                    terms_list.append(current_terms[:3])
                    sizes.append(current_cluster_size)
                    pvals.append(current_first_pval)
                    current_terms = []
                    current_first_pval = None
                    current_cluster_size = None
                parts = line.strip().split('\t')
                for p in parts:
                    if p.startswith('EnrichmentScore:'):
                        try:
                            scores.append(float(p.split(':', 1)[1]))
                        except Exception:
                            scores.append(None)
                inside_data = False
            elif line.startswith('Category\tTerm'):
                inside_data = True
            elif inside_data and line.strip():
                fields = line.strip().split('\t')
                if len(fields) >= 13:
                    term = fields[1]
                    try:
                        pval = float(fields[4])
                    except Exception:
                        pval = None
                    if current_first_pval is None:
                        current_first_pval = pval
                        try:
                            current_cluster_size = int(float(fields[2]))
                        except Exception:
                            current_cluster_size = None
                    current_terms.append(term)

    if current_terms:
        terms_list.append(current_terms[:3])
        sizes.append(current_cluster_size)
        pvals.append(current_first_pval)

    scores += [None] * max(0, (max_clusters - len(scores)))
    sizes += [None] * max(0, (max_clusters - len(sizes)))
    terms_list += [''] * max(0, (max_clusters - len(terms_list)))
    pvals += [None] * max(0, (max_clusters - len(pvals)))

    cleaned_terms = []
    for t in terms_list[:max_clusters]:
        if isinstance(t, list):
            cleaned = [_clean_term_label(term) for term in t]
            cleaned_terms.append('; '.join(cleaned))
        else:
            cleaned_terms.append('')

    return scores[:max_clusters], pvals[:max_clusters], cleaned_terms, sizes[:max_clusters]


# Backward-compatible access through the class, as used elsewhere in the pipeline.
DAVIDEnrichmentRunner.save_cluster_report = staticmethod(save_cluster_report)
DAVIDEnrichmentRunner.extract_clusters = staticmethod(extract_clusters)


def extract_chart_report_to_df(chart_report) -> pd.DataFrame:
    def _to_list(obj):
        if obj is None:
            return []
        if isinstance(obj, (list, tuple)):
            return list(obj)
        for attr in ("item", "items", "simpleChartRecords", "simpleChartRecord"):
            try:
                val = getattr(obj, attr)
            except Exception:
                val = None
            if val is not None:
                if isinstance(val, (list, tuple)):
                    return list(val)
                return [val]
        return [obj]

    rows = []
    records = _to_list(chart_report)
    if not records:
        return pd.DataFrame(columns=[
            'Category', 'Term', 'CleanTerm', 'ListHits', 'Percent', 'Pvalue', 'Genes',
            'ListTotal', 'PopHits', 'PopTotal', 'FoldEnrichment', 'Bonferroni', 'Benjamini', 'FDR',
            'GeneListIDs', 'GeneListEntrez', 'NGeneIDsParsed', 'CleanTermLower'
        ])

    for rec in records:
        try:
            term = _clean_string(getattr(rec, 'termName', ''))
            genes = _clean_string(getattr(rec, 'geneIds', ''))
            rows.append({
                'Category': _clean_string(getattr(rec, 'categoryName', '')),
                'Term': term,
                'CleanTerm': _clean_term_label(term),
                'ListHits': getattr(rec, 'listHits', None),
                'Percent': getattr(rec, 'percent', None),
                'Pvalue': getattr(rec, 'ease', None),
                'Genes': genes,
                'ListTotal': getattr(rec, 'listTotals', None),
                'PopHits': getattr(rec, 'popHits', None),
                'PopTotal': getattr(rec, 'popTotals', None),
                'FoldEnrichment': getattr(rec, 'foldEnrichment', None),
                'Bonferroni': getattr(rec, 'bonferroni', None),
                'Benjamini': getattr(rec, 'benjamini', None),
                'FDR': getattr(rec, 'afdr', None),
            })
        except Exception:
            continue

    df = pd.DataFrame(rows)
    if df.empty:
        return pd.DataFrame(columns=[
            'Category', 'Term', 'CleanTerm', 'ListHits', 'Percent', 'Pvalue', 'Genes',
            'ListTotal', 'PopHits', 'PopTotal', 'FoldEnrichment', 'Bonferroni', 'Benjamini', 'FDR',
            'GeneListIDs', 'GeneListEntrez', 'NGeneIDsParsed', 'CleanTermLower'
        ])

    df['GeneListIDs'] = df['Genes'].map(_split_generic_ids_field)
    # Kept for backward compatibility with earlier exports.
    df['GeneListEntrez'] = df['GeneListIDs']
    df['NGeneIDsParsed'] = df['GeneListIDs'].map(len)
    df['CleanTermLower'] = df['CleanTerm'].astype(str).str.lower()
    return df


def build_reordered_gene_mapping(ordered_genes: List[str], geneids_df: pd.DataFrame,
                                 alias_map: Optional[Dict[str, str]] = None) -> pd.DataFrame:
    if geneids_df is None or geneids_df.empty:
        raise ValueError("No FASTA-derived GeneIDs table is available for DAVID mapping.")

    df = geneids_df.copy().fillna("")
    rows = []

    id_lookup: Dict[str, Dict[str, str]] = {}
    for _, r in df.iterrows():
        entrez = _clean_entrez(r.get("EntrezGeneID", ""))
        display_locus = _first_nonempty([
            r.get("RefSeq_LocusTag_RS", ""),
            r.get("Old_LocusTag", ""),
            r.get("LocusTag", ""),
            r.get("PrimaryID", ""),
            r.get("GeneSymbol", ""),
        ])
        refseq_protein = _first_nonempty([
            r.get("RefSeqProteinID", ""),
            r.get("ProteinID", ""),
            r.get("protein_id", ""),
        ])
        uniprot = _first_nonempty([
            r.get("UniProtID", ""),
            r.get("UniprotID", ""),
            r.get("UniProt", ""),
            r.get("uniprot", ""),
        ])
        payload = {
            "EntrezGeneID": entrez,
            "DisplayLocusTag": display_locus,
            "GeneSymbol": _clean_string(r.get("GeneSymbol", "")),
            "ProteinDescription": _clean_string(r.get("ProteinDescription", "")),
            "PrimaryID": _clean_string(r.get("PrimaryID", "")),
            "RefSeqProteinID": _clean_generic_david_id(refseq_protein),
            "UniProtID": _clean_generic_david_id(uniprot),
            "LocusTag": _clean_generic_david_id(r.get("LocusTag", "")),
            "Old_LocusTag": _clean_generic_david_id(r.get("Old_LocusTag", "")),
            "RefSeq_LocusTag_RS": _clean_generic_david_id(r.get("RefSeq_LocusTag_RS", "")),
        }

        candidates = [
            r.get("LocusTag", ""),
            r.get("RefSeq_LocusTag_RS", ""),
            r.get("Old_LocusTag", ""),
            r.get("PrimaryID", ""),
            r.get("GeneSymbol", ""),
            refseq_protein,
            uniprot,
        ]
        for c in candidates:
            cc = _clean_string(c)
            if not cc:
                continue
            for key in {cc, _canonicalize_token(cc, alias_map)}:
                if key and key not in id_lookup:
                    id_lookup[key] = payload

    for i, gene in enumerate(list(ordered_genes), start=1):
        g = _clean_string(gene)
        g_can = _canonicalize_token(g, alias_map)
        payload = id_lookup.get(g) or id_lookup.get(g_can) or {}
        entrez = _clean_entrez(payload.get("EntrezGeneID", ""))
        display_locus = _first_nonempty([payload.get("DisplayLocusTag", ""), g_can, g])

        refseq_protein = _clean_generic_david_id(payload.get("RefSeqProteinID", ""))
        uniprot = _clean_generic_david_id(payload.get("UniProtID", ""))

        rows.append({
            "OrderedIndex": i,
            "OrderedGene": g,
            "OrderedGeneCanonical": g_can,
            "DisplayLocusTag": display_locus,
            "EntrezGeneID": entrez,
            "RefSeqProteinID": refseq_protein,
            "UniProtID": uniprot,
            "GeneSymbol": _clean_string(payload.get("GeneSymbol", "")),
            "ProteinDescription": _clean_string(payload.get("ProteinDescription", "")),
            "MappedToEntrez": bool(entrez),
            "MappedToRefSeqProtein": bool(refseq_protein),
            "MappedToUniProt": bool(uniprot),
            "MappedToLocusTag": bool(_clean_generic_david_id(display_locus)),
        })

    return pd.DataFrame(rows)


def build_gene_term_table(mapping_df: pd.DataFrame, genome_chart_df: pd.DataFrame,
                          david_id_col: str = "EntrezGeneID",
                          david_id_type: str = "ENTREZ_GENE_ID") -> pd.DataFrame:
    term_map = defaultdict(list)
    cleaner = _clean_entrez if str(david_id_col) == "EntrezGeneID" else _clean_generic_david_id

    if genome_chart_df is not None and not genome_chart_df.empty:
        for _, row in genome_chart_df.iterrows():
            term = _clean_string(row.get("CleanTerm", ""))
            gids = row.get("GeneListIDs", [])
            if not isinstance(gids, list):
                gids = _split_generic_ids_field(gids)
            for gid in gids:
                gid_clean = cleaner(gid)
                if gid_clean and term and term not in term_map[gid_clean]:
                    term_map[gid_clean].append(term)

    out_rows = []
    for _, row in mapping_df.iterrows():
        query_id = cleaner(row.get(david_id_col, ""))
        terms = term_map.get(query_id, [])
        out_rows.append({
            "OrderedIndex": row.get("OrderedIndex", ""),
            "DisplayLocusTag": row.get("DisplayLocusTag", ""),
            "EntrezGeneID": _clean_entrez(row.get("EntrezGeneID", "")),
            "RefSeqProteinID": _clean_generic_david_id(row.get("RefSeqProteinID", "")),
            "UniProtID": _clean_generic_david_id(row.get("UniProtID", "")),
            "GeneSymbol": row.get("GeneSymbol", ""),
            "ProteinDescription": row.get("ProteinDescription", ""),
            "DAVID_QueryID": query_id,
            "DAVID_ID_Type": david_id_type,
            "N_DAVID_terms": len(terms),
            "DAVID_Terms": "; ".join(terms),
        })

    return pd.DataFrame(out_rows)


def _build_gene_term_table_fallback(mapping_df: pd.DataFrame,
                                    david_id_col: str = "EntrezGeneID",
                                    david_id_type: str = "ENTREZ_GENE_ID") -> pd.DataFrame:
    fallback = mapping_df.copy()
    for col in ["OrderedIndex", "DisplayLocusTag", "EntrezGeneID", "RefSeqProteinID", "UniProtID", "GeneSymbol", "ProteinDescription"]:
        if col not in fallback.columns:
            fallback[col] = ""
    cleaner = _clean_entrez if str(david_id_col) == "EntrezGeneID" else _clean_generic_david_id
    fallback["DAVID_QueryID"] = fallback.get(david_id_col, "").map(cleaner)
    fallback["DAVID_ID_Type"] = david_id_type
    fallback["N_DAVID_terms"] = 0
    fallback["DAVID_Terms"] = ""
    cols = [
        "OrderedIndex", "DisplayLocusTag", "EntrezGeneID", "RefSeqProteinID", "UniProtID", "GeneSymbol",
        "ProteinDescription", "DAVID_QueryID", "DAVID_ID_Type", "N_DAVID_terms", "DAVID_Terms"
    ]
    return fallback.loc[:, cols].fillna("")


def select_top_windows(summary_df: pd.DataFrame, top_n_hits: int) -> pd.DataFrame:
    work = summary_df.copy()
    work["Enrich1_num"] = pd.to_numeric(work.get("Enrich1", np.nan), errors="coerce")
    work["Pval1_num"] = pd.to_numeric(work.get("Pval1", np.nan), errors="coerce")
    work = work.dropna(subset=["Enrich1_num"]).sort_values(
        by=["Enrich1_num", "Pval1_num"],
        ascending=[False, True],
        na_position="last"
    )
    return work.head(int(top_n_hits)).copy()


def collect_auto_queries_from_top_windows(top_windows_df: pd.DataFrame,
                                          max_clusters: int = 3) -> pd.DataFrame:
    rows = []
    seen = OrderedDict()

    if top_windows_df is None or top_windows_df.empty:
        return pd.DataFrame(columns=[
            "Query", "Source", "TopWindowRank", "Window", "SourceCluster",
            "WindowEnrich1", "WindowPval1"
        ])

    for rank, (_, row) in enumerate(top_windows_df.iterrows(), start=1):
        window = _clean_string(row.get("Window", ""))
        enrich1 = row.get("Enrich1", np.nan)
        pval1 = row.get("Pval1", np.nan)

        for i in range(1, max_clusters + 1):
            txt = _clean_string(row.get(f"Cluster {i} Terms", ""))
            if not txt:
                continue
            for term in [x.strip() for x in txt.split(";") if x.strip()]:
                q = _clean_term_label(term)
                ql = q.lower()
                if q and ql not in seen:
                    seen[ql] = True
                    rows.append({
                        "Query": q,
                        "Source": "auto_from_top_windows",
                        "TopWindowRank": rank,
                        "Window": window,
                        "SourceCluster": i,
                        "WindowEnrich1": enrich1,
                        "WindowPval1": pval1,
                    })

    return pd.DataFrame(rows)


def build_term_cluster_outputs(genome_chart_df: pd.DataFrame,
                               mapping_df: pd.DataFrame,
                               term_queries: List[str],
                               match_mode: str = "contains",
                               david_id_col: str = "EntrezGeneID",
                               david_id_type: str = "ENTREZ_GENE_ID"):
    if genome_chart_df is None or genome_chart_df.empty:
        empty = pd.DataFrame()
        return empty, empty

    cleaner = _clean_entrez if str(david_id_col) == "EntrezGeneID" else _clean_generic_david_id

    mapped = mapping_df.copy()
    mapped["DAVID_QueryID"] = mapped.get(david_id_col, "").map(cleaner)
    mapped = mapped[mapped["DAVID_QueryID"] != ""].copy()
    mapped = mapped.drop_duplicates(subset=["DAVID_QueryID"], keep="first")

    gene_order = {gid: i for i, gid in enumerate(mapped["DAVID_QueryID"].tolist())}
    gid_to_locus = dict(zip(mapped["DAVID_QueryID"], mapped["DisplayLocusTag"]))

    queries_clean = []
    seen = set()
    for q in term_queries:
        qq = _clean_term_label(q)
        if qq and qq.lower() not in seen:
            queries_clean.append(qq)
            seen.add(qq.lower())

    cluster_cols = OrderedDict()
    match_rows = []
    raw_headers = []

    for query in queries_clean:
        ql = query.lower()
        if match_mode == "exact":
            hit_df = genome_chart_df[genome_chart_df["CleanTermLower"] == ql].copy()
        else:
            hit_df = genome_chart_df[genome_chart_df["CleanTermLower"].str.contains(re.escape(ql), na=False)].copy()

        all_gene_ids = []
        for _, row in hit_df.iterrows():
            gids = row.get("GeneListIDs", [])
            if not isinstance(gids, list):
                gids = _split_generic_ids_field(gids)
            gids = [cleaner(g) for g in gids if cleaner(g)]
            all_gene_ids.extend(gids)

            match_rows.append({
                "Query": query,
                "DAVID_ID_Type": david_id_type,
                "MatchedCategory": row.get("Category", ""),
                "MatchedTerm": row.get("CleanTerm", ""),
                "MatchedFullTerm": row.get("Term", ""),
                "MatchedGeneCount": len(gids),
                "MatchedGenes_DAVID_IDs": ", ".join(gids),
                "MatchedGenes_LocusTags": ", ".join([gid_to_locus.get(g, "") for g in gids if gid_to_locus.get(g, "")]),
                "Pvalue": row.get("Pvalue", ""),
                "FoldEnrichment": row.get("FoldEnrichment", ""),
            })

        all_gene_ids = list(OrderedDict.fromkeys([g for g in all_gene_ids if g in gid_to_locus]))
        all_gene_ids = sorted(all_gene_ids, key=lambda x: gene_order.get(x, 10**12))
        locus_tags = [gid_to_locus[g] for g in all_gene_ids if gid_to_locus.get(g, "")]

        raw_headers.append(query)
        cluster_cols[query] = pd.Series(locus_tags, dtype=object)

    unique_headers = _make_unique_headers(raw_headers)
    cluster_df = pd.DataFrame({
        unique_headers[i]: cluster_cols[list(cluster_cols.keys())[i]]
        for i in range(len(unique_headers))
    }) if unique_headers else pd.DataFrame()

    match_detail_df = pd.DataFrame(match_rows)
    return cluster_df, match_detail_df


def _prepare_gene_term_df_for_export(gene_term_df: pd.DataFrame) -> pd.DataFrame:
    if gene_term_df is None:
        return pd.DataFrame(columns=[
            "OrderedIndex", "DisplayLocusTag", "EntrezGeneID", "GeneSymbol",
            "ProteinDescription", "DAVID_QueryID", "DAVID_ID_Type", "N_DAVID_terms", "DAVID_Terms"
        ])

    df = gene_term_df.copy()
    if "DAVID_Categories" in df.columns:
        df = df.drop(columns=["DAVID_Categories"])

    desired = [
        "OrderedIndex", "DisplayLocusTag", "EntrezGeneID", "RefSeqProteinID", "UniProtID", "GeneSymbol",
        "ProteinDescription", "DAVID_QueryID", "DAVID_ID_Type", "N_DAVID_terms", "DAVID_Terms"
    ]
    present = [c for c in desired if c in df.columns]
    remaining = [c for c in df.columns if c not in present]
    df = df.loc[:, present + remaining]
    return df.fillna("")


def write_gene_term_txt(output_path: str, gene_term_df: pd.DataFrame):
    export_df = _prepare_gene_term_df_for_export(gene_term_df)
    export_df.to_csv(output_path, sep="\t", index=False, encoding="utf-8")


def append_gene_term_sheet_to_geneids_excel(geneids_xlsx_path: str, gene_term_df: pd.DataFrame,
                                            sheet_name: str = "DAVID Gene Terms"):
    if not geneids_xlsx_path or not os.path.isfile(geneids_xlsx_path):
        return
    export_df = _prepare_gene_term_df_for_export(gene_term_df)
    try:
        with pd.ExcelWriter(
            geneids_xlsx_path,
            engine="openpyxl",
            mode="a",
            if_sheet_exists="replace"
        ) as writer:
            export_df.to_excel(writer, sheet_name=sheet_name, index=False)
        print(f"[INFO] Appended DAVID gene-term sheet to:\n  {geneids_xlsx_path}")
    except Exception as e:
        print(f"[WARN] Could not append DAVID gene-term sheet to GeneIDs workbook: {e}")


def write_standalone_david_excel(output_path: str,
                                 summary_df: pd.DataFrame,
                                 filtered_df: pd.DataFrame,
                                 top_hits_df: pd.DataFrame,
                                 mapping_df: pd.DataFrame,
                                 genome_chart_df: pd.DataFrame,
                                 gene_term_df: pd.DataFrame,
                                 auto_queries_df: pd.DataFrame,
                                 term_cluster_df: pd.DataFrame,
                                 term_match_detail_df: pd.DataFrame):
    export_gene_term_df = _prepare_gene_term_df_for_export(gene_term_df)
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        summary_df.to_excel(writer, sheet_name="Window scan all", index=False)
        filtered_df.to_excel(writer, sheet_name="Window scan filtered", index=False)
        top_hits_df.to_excel(writer, sheet_name="Top enrichment hits", index=False)
        export_gene_term_df.to_excel(writer, sheet_name="Gene to DAVID terms", index=False)
        term_cluster_df.to_excel(writer, sheet_name="Term derived clusters", index=False)
        term_match_detail_df.to_excel(writer, sheet_name="Term match detail", index=False)
        auto_queries_df.to_excel(writer, sheet_name="Auto term queries", index=False)
        mapping_df.to_excel(writer, sheet_name="Gene mapping", index=False)
        genome_chart_df.to_excel(writer, sheet_name="Genome chart raw", index=False)


def _preflight_david_strategy(runner: DAVIDEnrichmentRunner,
                              genome_ids: List[str],
                              output_prefix: str,
                              id_type: str,
                              chart_threshold: float,
                              chart_count: int) -> Tuple[pd.DataFrame, object]:
    report = runner.get_chart_report(
        genome_ids,
        output_prefix=output_prefix,
        id_type=id_type,
        threshold=float(chart_threshold),
        count=int(chart_count),
    )
    chart_df = extract_chart_report_to_df(report)
    if chart_df.empty:
        raise ValueError("DAVID returned no functional annotation rows for this identifier scheme.")
    return chart_df, report


def _run_scan_for_identifier_scheme(
    mapping_df: pd.DataFrame,
    runner: DAVIDEnrichmentRunner,
    id_col: str,
    id_type: str,
    id_source_label: str,
    id_cleaner: Callable[[object], str],
    output_folder: str,
    output_prefix: str,
    geneids_xlsx_path: Optional[str],
    window_size: int,
    step_size: int,
    wait_time: float,
    top_n_hits: int,
    max_clusters: int,
    report_subdir_name: str,
    min_valid_ids_per_window: int,
    plot_format: str,
    chart_threshold: float,
    chart_count: int,
    manual_term_queries: Optional[List[str]],
    term_match_mode: str,
    append_to_geneids_excel: bool,
    write_standalone_excel: bool,
    precomputed_genome_chart_df: pd.DataFrame,
):
    os.makedirs(output_folder, exist_ok=True)
    report_dir = os.path.join(output_folder, report_subdir_name)
    os.makedirs(report_dir, exist_ok=True)

    summary_records = []
    total_genes = len(mapping_df)

    for start0 in range(0, total_genes - int(window_size) + 1, int(step_size)):
        end0 = start0 + int(window_size)
        subset = mapping_df.iloc[start0:end0].copy()
        start1, end1 = start0 + 1, end0
        window_label = f"{start1}-{end1}"
        output_prefix_window = f"{output_prefix}_{start1}to{end1}"
        print(f"[INFO] DAVID window scan ({id_source_label}, {id_type}): {window_label}")

        subset_ids = _collect_unique_ids(subset[id_col].astype(str).tolist(), id_cleaner)

        enrich_scores = [None] * max_clusters
        min_pvals = [None] * max_clusters
        cluster_terms = [""] * max_clusters
        cluster_sizes = [None] * max_clusters
        report_file = os.path.join(report_dir, f"{output_prefix_window}_fullReport.txt")
        status = "OK"

        if len(subset_ids) < int(min_valid_ids_per_window):
            status = f"Skipped: < {int(min_valid_ids_per_window)} mapped DAVID IDs"
            with open(report_file, 'w', encoding='utf-8') as fout:
                fout.write(status + "\n")
        else:
            try:
                clustering = runner.analyze_gene_list(subset_ids, output_prefix_window, id_type=id_type)
                DAVIDEnrichmentRunner.save_cluster_report(clustering, report_file)
                enrich_scores, min_pvals, cluster_terms, cluster_sizes = DAVIDEnrichmentRunner.extract_clusters(
                    report_file, max_clusters=max_clusters
                )
            except Exception as e:
                status = f"DAVID error: {e}"
                with open(report_file, 'w', encoding='utf-8') as fout:
                    fout.write(status + "\n")

        record = {
            "Window": window_label,
            "Start": start1,
            "End": end1,
            "NWindowGenes": len(subset),
            "NMappedDAVIDIDs": len(subset_ids),
            "DAVID_ID_Type": id_type,
            "DAVID_ID_Source": id_source_label,
            "ReportFile": report_file,
            "Status": status,
        }
        for i in range(1, max_clusters + 1):
            record[f"Enrich{i}"] = enrich_scores[i - 1]
            record[f"Pval{i}"] = min_pvals[i - 1]
            record[f"Cluster {i} Terms"] = cluster_terms[i - 1]
            record[f"Size{i}"] = cluster_sizes[i - 1]
        summary_records.append(record)

        if float(wait_time) > 0:
            time.sleep(float(wait_time))

    summary_df = pd.DataFrame(summary_records)
    if not summary_df.empty and 'Status' in summary_df.columns:
        status_series = summary_df['Status'].fillna('').astype(str)
        n_ok = int((status_series == 'OK').sum())
        n_skipped = int(status_series.str.startswith('Skipped').sum())
        n_failed = int(status_series.str.startswith('DAVID error').sum())
        print(f"[INFO] DAVID window-scan summary: OK={n_ok}, skipped={n_skipped}, failed={n_failed}")
    if summary_df.empty:
        raise ValueError("No sliding windows were evaluated. Check the chosen window size and genome length.")

    for idx in summary_df.index:
        for i in [2, 3]:
            p = summary_df.at[idx, f'Pval{i}'] if f'Pval{i}' in summary_df.columns else None
            try:
                pnum = float(p)
            except Exception:
                pnum = None
            if (pnum is None) or (pnum > 0.01):
                if f'Enrich{i}' in summary_df.columns:
                    summary_df.at[idx, f'Enrich{i}'] = np.nan
                if f'Pval{i}' in summary_df.columns:
                    summary_df.at[idx, f'Pval{i}'] = np.nan
                if f'Cluster {i} Terms' in summary_df.columns:
                    summary_df.at[idx, f'Cluster {i} Terms'] = ''
                if f'Size{i}' in summary_df.columns:
                    summary_df.at[idx, f'Size{i}'] = np.nan

    filtered_df = summary_df.dropna(subset=['Pval1', 'Pval2', 'Pval3'], how='all').copy()
    filtered_final_df = filtered_df[
        ['Window', 'Enrich1', 'Enrich2', 'Enrich3',
         'Cluster 1 Terms', 'Cluster 2 Terms', 'Cluster 3 Terms',
         'Size1', 'Size2', 'Size3']
    ].sort_values(by='Enrich1', ascending=False, na_position='last') if not filtered_df.empty else pd.DataFrame(
        columns=['Window', 'Enrich1', 'Enrich2', 'Enrich3',
                 'Cluster 1 Terms', 'Cluster 2 Terms', 'Cluster 3 Terms',
                 'Size1', 'Size2', 'Size3']
    )

    top_windows_df = select_top_windows(summary_df, top_n_hits=top_n_hits)

    raw_headers = []
    raw_values = []
    for _, hit in top_windows_df.iterrows():
        start1 = int(hit['Start'])
        end1 = int(hit['End'])
        subset = mapping_df.iloc[start1 - 1:end1].copy()
        locus_tags = [x for x in subset['DisplayLocusTag'].astype(str).tolist() if _clean_string(x) != ""]
        locus_tags = list(OrderedDict.fromkeys(locus_tags))
        header = _clean_string(hit.get('Cluster 1 Terms', '')) or _clean_string(hit.get('Window', ''))
        raw_headers.append(header)
        raw_values.append(locus_tags)

    top_headers = _make_unique_headers(raw_headers)
    top_hits_df = pd.DataFrame({
        top_headers[i]: pd.Series(raw_values[i], dtype=object)
        for i in range(len(top_headers))
    })

    genome_chart_df = precomputed_genome_chart_df.copy() if precomputed_genome_chart_df is not None else pd.DataFrame()
    gene_term_df = build_gene_term_table(
        mapping_df,
        genome_chart_df,
        david_id_col=id_col,
        david_id_type=id_type,
    ) if not genome_chart_df.empty else _build_gene_term_table_fallback(
        mapping_df,
        david_id_col=id_col,
        david_id_type=id_type,
    )

    manual_terms_clean = [x for x in (manual_term_queries or []) if _clean_string(x) != ""]
    if manual_terms_clean:
        auto_queries_df = pd.DataFrame({
            "Query": manual_terms_clean,
            "Source": ["manual"] * len(manual_terms_clean),
            "TopWindowRank": [np.nan] * len(manual_terms_clean),
            "Window": [""] * len(manual_terms_clean),
            "SourceCluster": [np.nan] * len(manual_terms_clean),
            "WindowEnrich1": [np.nan] * len(manual_terms_clean),
            "WindowPval1": [np.nan] * len(manual_terms_clean),
        })
    else:
        auto_queries_df = collect_auto_queries_from_top_windows(
            top_windows_df=top_windows_df,
            max_clusters=max_clusters,
        )

    query_terms = auto_queries_df["Query"].astype(str).tolist() if not auto_queries_df.empty else []

    term_cluster_df, term_match_detail_df = build_term_cluster_outputs(
        genome_chart_df=genome_chart_df,
        mapping_df=mapping_df,
        term_queries=query_terms,
        match_mode=term_match_mode,
        david_id_col=id_col,
        david_id_type=id_type,
    )

    summary_df['NegLog10_Pval1'] = pd.to_numeric(summary_df['Pval1'], errors='coerce')
    summary_df['NegLog10_Pval1'] = -np.log10(summary_df['NegLog10_Pval1'])
    plot_format = str(plot_format or 'png').lstrip('.')

    enrichment_plot_path = os.path.join(output_folder, f'{output_prefix}_david_enrichment_scores.{plot_format}')
    pvalue_plot_path = os.path.join(output_folder, f'{output_prefix}_david_neglog10_pval.{plot_format}')

    plt.figure(figsize=(10, 4))
    plt.bar(summary_df['Window'], pd.to_numeric(summary_df['Enrich1'], errors='coerce'), width=0.8)
    plt.xlabel('Window')
    plt.ylabel('Enrichment Score (Cluster 1)')
    plt.title('DAVID sliding-window enrichment scores along reordered genome')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(enrichment_plot_path, dpi=300)
    plt.close()

    plt.figure(figsize=(10, 4))
    plt.bar(summary_df['Window'], summary_df['NegLog10_Pval1'], width=0.8)
    plt.xlabel('Window')
    plt.ylabel('-log10(Pval) Cluster 1')
    plt.title('DAVID sliding-window -log10(Pvalue) along reordered genome')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(pvalue_plot_path, dpi=300)
    plt.close()

    if bool(append_to_geneids_excel):
        append_gene_term_sheet_to_geneids_excel(
            geneids_xlsx_path=geneids_xlsx_path,
            gene_term_df=gene_term_df,
            sheet_name="DAVID Gene Terms",
        )

    gene2terms_txt_path = os.path.join(output_folder, "DAVID gene2terms.txt")
    write_gene_term_txt(gene2terms_txt_path, gene_term_df)
    print(f"[INFO] DAVID gene-to-terms TXT saved:\n  {gene2terms_txt_path}")

    david_excel_path = os.path.join(output_folder, f"{output_prefix}_DAVID.xlsx")
    if bool(write_standalone_excel):
        write_standalone_david_excel(
            output_path=david_excel_path,
            summary_df=summary_df,
            filtered_df=filtered_final_df,
            top_hits_df=top_hits_df,
            mapping_df=mapping_df,
            genome_chart_df=genome_chart_df,
            gene_term_df=gene_term_df,
            auto_queries_df=auto_queries_df,
            term_cluster_df=term_cluster_df,
            term_match_detail_df=term_match_detail_df,
        )
        print(f"[INFO] DAVID standalone Excel saved:\n  {david_excel_path}")

    return {
        'summary_df': summary_df,
        'filtered_df': filtered_final_df,
        'top_hits_df': top_hits_df,
        'mapping_df': mapping_df,
        'gene_term_df': gene_term_df,
        'genome_chart_df': genome_chart_df,
        'term_cluster_df': term_cluster_df,
        'term_match_detail_df': term_match_detail_df,
        'auto_queries_df': auto_queries_df,
        'top_windows_df': top_windows_df,
        'report_dir': report_dir,
        'david_excel_path': david_excel_path,
        'gene2terms_txt_path': gene2terms_txt_path,
        'enrichment_plot_path': enrichment_plot_path,
        'pvalue_plot_path': pvalue_plot_path,
        'n_mapped_genes': len(_collect_unique_ids(mapping_df.get(id_col, []).tolist(), id_cleaner)),
        'n_total_genes': total_genes,
        'david_id_type': id_type,
        'david_id_source': id_source_label,
    }


def run_david_window_scan_from_ordered_genes(
    ordered_genes: List[str],
    geneids_df: pd.DataFrame,
    output_folder: str,
    output_prefix: str,
    user_email: str,
    alias_map: Optional[Dict[str, str]] = None,
    geneids_xlsx_path: Optional[str] = None,
    window_size: int = 100,
    step_size: int = 50,
    wait_time: float = 0.0,
    top_n_hits: int = 10,
    max_clusters: int = 3,
    report_subdir_name: str = "DAVID window reports",
    wsdl_url: str = DEFAULT_DAVID_WSDL_URL,
    endpoint: str = DEFAULT_DAVID_ENDPOINT,
    min_valid_ids_per_window: int = 3,
    plot_format: str = "png",
    chart_threshold: float = 1.0,
    chart_count: int = 1,
    manual_term_queries: Optional[List[str]] = None,
    term_match_mode: str = "contains",
    append_to_geneids_excel: bool = True,
    write_standalone_excel: bool = True,
    allow_insecure_ssl_fallback: bool = True,
):
    if int(window_size) <= 0 or int(step_size) <= 0:
        raise ValueError("window_size and step_size must be positive integers.")

    mapping_df = build_reordered_gene_mapping(ordered_genes, geneids_df, alias_map=alias_map)
    if mapping_df.empty:
        raise ValueError("The reordered genome could not be mapped to any DAVID-compatible identifiers.")

    def _ids_from_mapping(column_name: str, cleaner) -> List[str]:
        if column_name not in mapping_df.columns:
            return []
        return _collect_unique_ids(mapping_df[column_name].astype(str).tolist(), cleaner)

    locus_ids = _ids_from_mapping("DisplayLocusTag", _clean_generic_david_id)
    entrez_ids = _ids_from_mapping("EntrezGeneID", _clean_entrez)
    refseq_protein_ids = _ids_from_mapping("RefSeqProteinID", _clean_generic_david_id)
    uniprot_ids = _ids_from_mapping("UniProtID", _clean_generic_david_id)
    gene_symbol_ids = _ids_from_mapping("GeneSymbol", _clean_generic_david_id)

    id_count_msg = (
        f"locus tags={len(locus_ids)}, Entrez={len(entrez_ids)}, "
        f"RefSeq proteins={len(refseq_protein_ids)}, UniProt={len(uniprot_ids)}, "
        f"gene symbols={len(gene_symbol_ids)}"
    )
    print(f"[INFO] DAVID identifier availability in Gene IDs table: {id_count_msg}")

    if not any([locus_ids, entrez_ids, refseq_protein_ids, uniprot_ids, gene_symbol_ids]):
        raise ValueError(
            "DAVID could not be run because no potentially DAVID-compatible identifiers were available in the Gene IDs workbook.\n\n"
            "The online DAVID scan requires identifiers recognized by DAVID, ideally Entrez Gene IDs. "
            "Use an NCBI RefSeq CDS FASTA containing db_xref=GeneID entries, or provide a Gene IDs workbook with EntrezGeneID values."
        )

    os.makedirs(output_folder, exist_ok=True)
    runner = DAVIDEnrichmentRunner(user_email=user_email, wsdl_url=wsdl_url, endpoint=endpoint, allow_insecure_ssl_fallback=allow_insecure_ssl_fallback)

    attempt_errors = []
    strategies = []

    # Prefer stable database identifiers when available. Salmonella/other bacterial
    # locus tags are often not recognized by DAVID, whereas Entrez/RefSeq/UniProt
    # identifiers have a better chance of being accepted.
    if entrez_ids:
        strategies.append(("EntrezGeneID", "Entrez Gene IDs", "ENTREZ_GENE_ID", _clean_entrez, entrez_ids))
    if refseq_protein_ids:
        strategies.extend([
            ("RefSeqProteinID", "RefSeq protein IDs", "REFSEQ_PROTEIN", _clean_generic_david_id, refseq_protein_ids),
            ("RefSeqProteinID", "RefSeq protein IDs", "REFSEQ_PROTEIN_ID", _clean_generic_david_id, refseq_protein_ids),
        ])
    if uniprot_ids:
        strategies.extend([
            ("UniProtID", "UniProt IDs", "UNIPROT_ACCESSION", _clean_generic_david_id, uniprot_ids),
            ("UniProtID", "UniProt IDs", "UNIPROT_ID", _clean_generic_david_id, uniprot_ids),
        ])
    if gene_symbol_ids:
        strategies.extend([
            ("GeneSymbol", "gene symbols", "OFFICIAL_GENE_SYMBOL", _clean_generic_david_id, gene_symbol_ids),
            ("GeneSymbol", "gene symbols", "GENE_SYMBOL", _clean_generic_david_id, gene_symbol_ids),
        ])
    if locus_ids:
        strategies.extend([
            ("DisplayLocusTag", "locus tags", "LOCUS_TAG", _clean_generic_david_id, locus_ids),
            ("DisplayLocusTag", "locus tags", "GENE_SYMBOL", _clean_generic_david_id, locus_ids),
            ("DisplayLocusTag", "locus tags", "OFFICIAL_GENE_SYMBOL", _clean_generic_david_id, locus_ids),
        ])

    for id_col, id_source_label, id_type, cleaner, genome_ids in strategies:
        print(f"[INFO] Trying DAVID with {id_source_label} using id_type={id_type} ...")
        try:
            genome_chart_df, _ = _preflight_david_strategy(
                runner=runner,
                genome_ids=genome_ids,
                output_prefix=f"{output_prefix}_GENOME_PRECHECK_{id_type}",
                id_type=id_type,
                chart_threshold=chart_threshold,
                chart_count=chart_count,
            )
            print(f"[INFO] DAVID accepted {id_source_label} using id_type={id_type}.")
            return _run_scan_for_identifier_scheme(
                mapping_df=mapping_df,
                runner=runner,
                id_col=id_col,
                id_type=id_type,
                id_source_label=id_source_label,
                id_cleaner=cleaner,
                output_folder=output_folder,
                output_prefix=output_prefix,
                geneids_xlsx_path=geneids_xlsx_path,
                window_size=int(window_size),
                step_size=int(step_size),
                wait_time=float(wait_time),
                top_n_hits=int(top_n_hits),
                max_clusters=int(max_clusters),
                report_subdir_name=report_subdir_name,
                min_valid_ids_per_window=int(min_valid_ids_per_window),
                plot_format=plot_format,
                chart_threshold=float(chart_threshold),
                chart_count=int(chart_count),
                manual_term_queries=manual_term_queries,
                term_match_mode=term_match_mode,
                append_to_geneids_excel=bool(append_to_geneids_excel),
                write_standalone_excel=bool(write_standalone_excel),
                precomputed_genome_chart_df=genome_chart_df,
            )
        except Exception as e:
            attempt_errors.append(f"{id_source_label} via {id_type}: {e}")
            msg = str(e)
            if _looks_like_david_unrecognized_ids_error(e) or ("no functional annotation rows" in msg.lower()):
                print(f"[WARN] DAVID did not accept {id_source_label} with id_type={id_type}: {e}")
                continue
            raise

    details = " | ".join(attempt_errors) if attempt_errors else "no DAVID attempts succeeded"
    raise ValueError(
        "DAVID did not recognize any identifier set from the current Gene IDs workbook.\n\n"
        f"Identifier counts tested: {id_count_msg}.\n\n"
        "This is not a CodonPipe clustering failure: DAVID is reachable, but the IDs submitted for this organism are not accepted/annotated by DAVID. "
        "For bacterial genomes, local locus tags such as SL1344_* are often rejected. The most reliable solution is to rerun CodonPipe with an NCBI RefSeq CDS FASTA that contains db_xref=GeneID entries, so that the Gene IDs workbook contains EntrezGeneID values. "
        "Alternatively, provide a DAVID gene2terms TXT file and use it for cluster inference instead of the online sliding-window scan.\n\n"
        "DAVID attempts: " + details
    )
