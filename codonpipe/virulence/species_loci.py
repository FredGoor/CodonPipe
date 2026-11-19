#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Evidence layer E6 -- species-specific acquired virulence loci.

Why this layer exists
---------------------
The pan-bacterial evidence layers (E1 curated virulence-factor databases, E2
UniProtKB) are built from individually characterised proteins of reference
proteomes.  That construction covers some species far better than others, and
it covers *Klebsiella pneumoniae* badly: the defining determinants of
hypervirulent *K. pneumoniae* are acquired loci carried on an integrative
conjugative element (ICE*Kp*) or on a pLVPK-like virulence plasmid, not
proteins of a curated reference proteome.  Two measured consequences for
*K. pneumoniae* HS11286 (assembly GCF_000240185.1):

* of the 137 genes in the authors' manual curation, only 42 have any
  VFDB / VICTORS / PATRIC_VF record at all;
* the RefSeq annotation assigns **no gene symbols whatsoever** -- all 5,779
  ``gene`` fields simply repeat the locus tag -- so every symbol-based matching
  route in E1/E2/E4 is inert and the product name is the only signal left.

This layer closes that gap by applying published, species-scoped **locus
definitions** to the RefSeq product names, through the module's existing
identifier crosswalk (:meth:`~codonpipe.virulence.resolve.GeneTable.lookup`)
and its existing operon/island machinery
(:func:`codonpipe.virulence.context.extend_by_context`).  The definitions live
in the editable resource table ``resources/species_virulence_loci.tsv``, one
row per locus, scoped by taxon so that the layer is a strict no-op for any
genome outside the configured scopes.

What this is, and what it is not
--------------------------------
The *Klebsiella* rows reproduce the five acquired virulence loci and the
0-5 virulence-score rule published for **Kleborate** (Lam *et al.* 2021).
They are applied here **by annotation matching only**.  This module is *not*
Kleborate and is not a substitute for it:

* it performs **no sequence-identity or coverage scoring** against the
  Kleborate allele databases, and therefore assigns **no lineage types**
  (no ybST, cbST, AbST, SmST);
* it performs **no MLST**;
* it performs **no Kaptive K-locus or O-locus capsule typing**;
* it requires **no external binaries at all** -- no minimap2, no mash, no
  BLAST, no Kaptive;
* because a locus is called from the presence of *annotated* genes, a genome
  with sparse annotation, or one carrying a pseudogenised, frameshifted or
  partially deleted locus, can be called differently from how Kleborate --
  which reads the nucleotide sequence -- would call it.

A score produced by this module must be reported as a **"Kleborate-style
virulence score"**, with the locus definitions cited to Lam *et al.* (2021),
and never as Kleborate output.

How detection works
-------------------
1. **Scope.**  Each row of the resource table declares a scope
   (``tax_id`` / ``species_tax_id`` / ``genus_tax_id`` / ``lineage_name``)
   matched against the resolved :class:`~codonpipe.virulence.resolve.TargetInfo`.
   Rows that do not match are never evaluated; if no row matches, the layer
   returns an empty result and changes nothing downstream.
2. **Marker matching.**  A gene is a *marker* of a locus if its ``gene`` field
   is one of the locus's constituent symbols (resolved through the crosswalk),
   or if its RefSeq product name matches the locus's ``product_patterns``.
   Patterns are compiled **case-sensitively** by default; see the warning
   below.
3. **Co-location.**  Markers are grouped into clusters on one replicon,
   splitting wherever more than ``max_gene_gap`` genes intervene.  A cluster is
   accepted as the locus only if it holds at least ``min_markers`` distinct
   marker genes.  This is what prevents a lone annotation hit from fabricating
   a locus.
4. **Recruitment.**  For an accepted cluster, :func:`extend_by_context` is run
   with the cluster's markers as seeds to recruit the remaining locus members
   whose product names are too generic to match a pattern -- HS11286 annotates
   *ybtX* as "MFS superfamily transporter signal transducer" and *ybtA* as
   "AraC-type transcriptional regulator".  Recruitment is clamped to the span
   of the accepted cluster, so it cannot propagate across the replicon.
5. **Scoring.**  Present loci contribute their ``score_token`` to the scheme
   named in ``score_scheme``; see :data:`SCORE_SCHEMES`.

.. warning::
   **Case sensitivity in ``product_patterns`` is load-bearing.**  The
   salmochelin receptor *iroN*, written case-insensitively as
   ``\\biro[BCDEN]\\b``, also matches the ordinary English word "iron", which
   appears in dozens of unrelated iron-transport product names -- in HS11286
   that single mistake matches 39 genes and would fabricate a salmochelin
   locus.  Written case-sensitively as ``\\b[Ii]ro[BCDEN]\\b`` it matches
   IroE / iroB / iroN and not "iron".  Patterns are therefore compiled without
   ``re.IGNORECASE``, and alternatives that genuinely should ignore case are
   wrapped in Python's scoped inline flag ``(?i:...)`` in the table itself.

A note on "completeness"
------------------------
The reported completeness is ``n_marker_genes / expected_gene_count``: the
fraction of the published locus that the **annotation names recognisably**.
It is not sequence coverage and not an alignment statistic.  A locus can be
fully present yet report low completeness simply because its annotation is
terse, which is why the recruited block size is reported alongside it.

References
----------
Lam MMC, Wick RR, Watts SC, Cerdeira LT, Wyres KL, Holt KE (2021) A genomic
surveillance framework and genotyping tool for *Klebsiella pneumoniae* and its
related species complex. *Nat Commun* 12:4188.
Lam MMC, Wick RR, Wyres KL, Gorrie CL, Judd LM, Jenney AWJ, Brisse S, Holt KE
(2018) Genetic diversity, mobilisation and spread of the
yersiniabactin-encoding mobile element ICE*Kp* in *Klebsiella pneumoniae*
populations. *Microb Genom* 4:e000196.
Lam MMC, Wyres KL, Judd LM, Wick RR, Jenney A, Brisse S, Holt KE (2018)
Tracking key virulence loci encoding aerobactin and salmochelin siderophore
synthesis in *Klebsiella pneumoniae*. *Genome Med* 10:77.
Walker KA, Miner TA, Palacios M, Trzilova D, Frederick DR, Broberg CA,
Sepulveda VE, Quinn JD, Miller VL (2019) A *Klebsiella pneumoniae* regulatory
mutant has reduced capsule expression but retains hypermucoviscosity.
*mBio* 10:e00089-19.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Set, Tuple

import pandas as pd

from .context import extend_by_context
from .evidence import load_resource_table
from .resolve import GeneTable, TargetInfo

__all__ = [
    "SPECIES_LOCI_TABLE",
    "LocusDefinition",
    "LocusCall",
    "SpeciesLociResult",
    "SCORE_SCHEMES",
    "KLEBORATE_SCORE_RULE",
    "load_species_loci",
    "definitions_in_scope",
    "kleborate_virulence_score",
    "detect_species_loci",
]

#: Name of the shipped resource table.
SPECIES_LOCI_TABLE = "species_virulence_loci.tsv"

_TRUE = {"1", "true", "yes", "y", "t"}


# ---------------------------------------------------------------------------
# the published Kleborate virulence score
# ---------------------------------------------------------------------------
#: Verbatim statement of the rule implemented by
#: :func:`kleborate_virulence_score`, reproduced into the provenance report so
#: that a published result states its own decision rule.
KLEBORATE_SCORE_RULE: Tuple[str, ...] = (
    "Kleborate-style virulence score, 0-5, after Lam MMC et al. (2021) "
    "Nat Commun 12:4188:",
    "  0 = none of yersiniabactin (ybt), colibactin (clb), aerobactin (iuc)",
    "  1 = yersiniabactin only",
    "  2 = colibactin without aerobactin (regardless of yersiniabactin)",
    "  3 = aerobactin only",
    "  4 = aerobactin and yersiniabactin, without colibactin",
    "  5 = yersiniabactin, colibactin and aerobactin all present",
    "Salmochelin (iro) and the rmp capsule up-regulators are detected and "
    "reported but do NOT enter the score: Lam et al. exclude iro because its "
    "presence is very strongly associated with aerobactin, and prioritise iuc "
    "because aerobactin is the stronger predictor of the hypervirulence "
    "phenotype.",
    "Edge case: the published prose defines 3 as 'aerobactin only'. The rule "
    "as implemented returns 3 for the rare aerobactin-plus-colibactin, "
    "yersiniabactin-negative combination, following the order in which the "
    "published conditions are stated (2 is reachable only when aerobactin is "
    "absent).",
)


def kleborate_virulence_score(tokens: Set[str]) -> int:
    """The published 0-5 *K. pneumoniae* virulence score.

    Parameters
    ----------
    tokens:
        Score tokens of the loci called present.  Only ``ybt``, ``clb`` and
        ``iuc`` affect the result; ``iro`` and ``rmp`` tokens are ignored by
        design (see :data:`KLEBORATE_SCORE_RULE`).

    Returns
    -------
    int
        A value in ``0..5``.
    """
    ybt = "ybt" in tokens
    clb = "clb" in tokens
    iuc = "iuc" in tokens
    if iuc and ybt and clb:
        return 5
    if iuc and ybt:
        return 4
    if iuc:
        return 3
    if clb:
        return 2
    if ybt:
        return 1
    return 0


#: Registry of scoring schemes, keyed by the ``score_scheme`` column of the
#: resource table.  Add an entry here to score a new species' loci; nothing
#: else in the module is species-specific.
SCORE_SCHEMES: Dict[str, Callable[[Set[str]], int]] = {
    "kleborate_ybt_clb_iuc": kleborate_virulence_score,
}

#: Human-readable rule text per scheme, for the provenance report.
SCORE_SCHEME_RULES: Dict[str, Tuple[str, ...]] = {
    "kleborate_ybt_clb_iuc": KLEBORATE_SCORE_RULE,
}


# ---------------------------------------------------------------------------
# resource table
# ---------------------------------------------------------------------------
@dataclass
class LocusDefinition:
    """One row of ``species_virulence_loci.tsv``, parsed and compiled."""

    scope_type: str
    scope_value: str
    locus: str
    locus_label: str
    score_scheme: str
    score_token: str
    counts_toward_score: bool
    gene_symbols: Tuple[str, ...]
    product_pattern: str
    regex: Optional["re.Pattern"]
    case_sensitive: bool
    expected_gene_count: int
    min_markers: int
    max_gene_gap: int
    extend_locus: bool
    machinery: str
    axis_hint: str
    weight: float
    reference: str
    note: str


def _as_int(value, default: int) -> int:
    try:
        return int(float(str(value).strip()))
    except Exception:
        return default


def _as_float(value, default: float) -> float:
    try:
        return float(str(value).strip())
    except Exception:
        return default


def _as_bool(value, default: bool) -> bool:
    text = str(value).strip().lower()
    if not text:
        return default
    return text in _TRUE


def load_species_loci(directory: Optional[str] = None) -> List[LocusDefinition]:
    """Read and compile the species-locus definitions.

    A missing file is **not** an error here, unlike the other resource tables:
    this layer is an optional extension, and a deployment that has not shipped
    the table should simply run without it rather than fail.  An *invalid*
    regular expression in a table that does exist is an error, because
    silently skipping a locus would change results without saying so.
    """
    folder = directory or os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources")
    if not os.path.isfile(os.path.join(folder, SPECIES_LOCI_TABLE)):
        return []
    frame = load_resource_table(SPECIES_LOCI_TABLE, folder)

    out: List[LocusDefinition] = []
    for row in frame.itertuples(index=False):
        pattern = str(getattr(row, "product_patterns", "") or "").strip()
        case_sensitive = _as_bool(getattr(row, "case_sensitive", "1"), True)
        regex = None
        if pattern:
            try:
                regex = re.compile(pattern, 0 if case_sensitive else re.IGNORECASE)
            except re.error as exc:
                raise ValueError(
                    f"Invalid product_patterns regular expression for locus "
                    f"{getattr(row, 'locus', '?')!r} in {SPECIES_LOCI_TABLE}: "
                    f"{pattern!r} ({exc})"
                )
        symbols = tuple(
            s.strip()
            for s in str(getattr(row, "gene_symbols", "") or "").replace(";", ",").split(",")
            if s.strip()
        )
        out.append(
            LocusDefinition(
                scope_type=str(getattr(row, "scope_type", "") or "").strip().lower(),
                scope_value=str(getattr(row, "scope_value", "") or "").strip(),
                locus=str(getattr(row, "locus", "") or "").strip(),
                locus_label=str(getattr(row, "locus_label", "") or "").strip(),
                score_scheme=str(getattr(row, "score_scheme", "") or "").strip(),
                score_token=str(getattr(row, "score_token", "") or "").strip(),
                counts_toward_score=_as_bool(getattr(row, "counts_toward_score", "0"), False),
                gene_symbols=symbols,
                product_pattern=pattern,
                regex=regex,
                case_sensitive=case_sensitive,
                expected_gene_count=_as_int(getattr(row, "expected_gene_count", 0), 0),
                min_markers=max(1, _as_int(getattr(row, "min_markers", 1), 1)),
                max_gene_gap=max(0, _as_int(getattr(row, "max_gene_gap", 6), 6)),
                extend_locus=_as_bool(getattr(row, "extend_locus", "1"), True),
                machinery=str(getattr(row, "machinery", "") or "unclassified").strip().lower(),
                axis_hint=str(getattr(row, "axis_hint", "") or "").strip(),
                weight=_as_float(getattr(row, "weight", 2.0), 2.0),
                reference=str(getattr(row, "reference", "") or "").strip(),
                note=str(getattr(row, "note", "") or "").strip(),
            )
        )
    return out


def definitions_in_scope(
    definitions: Sequence[LocusDefinition], target: TargetInfo
) -> List[LocusDefinition]:
    """Select the definitions whose scope matches the resolved target.

    Scope is matched against the taxonomy that
    :func:`codonpipe.virulence.resolve.resolve_target` already established, so
    a genome outside every configured scope yields an empty list and the layer
    becomes a no-op without any special-casing elsewhere.
    """
    lineage = {str(n).strip().lower() for n in (target.lineage_names or []) if str(n).strip()}
    if target.scientific_name:
        lineage.add(str(target.scientific_name).strip().lower())
    if target.species_name:
        lineage.add(str(target.species_name).strip().lower())
    if target.genus_name:
        lineage.add(str(target.genus_name).strip().lower())

    def matches(defn: LocusDefinition) -> bool:
        value = defn.scope_value.strip()
        if not value:
            return False
        if defn.scope_type == "tax_id":
            return str(target.tax_id or "") == value
        if defn.scope_type == "species_tax_id":
            return str(target.species_tax_id or "") == value
        if defn.scope_type == "genus_tax_id":
            return str(target.genus_tax_id or "") == value
        if defn.scope_type == "lineage_name":
            return value.lower() in lineage
        return False

    return [d for d in definitions if matches(d)]


# ---------------------------------------------------------------------------
# detection result
# ---------------------------------------------------------------------------
@dataclass
class LocusCall:
    """One detected (or sub-threshold) locus cluster."""

    locus: str
    locus_label: str
    present: bool
    seqid: str
    n_markers: int
    min_markers: int
    expected_gene_count: int
    completeness: float
    marker_genes: List[str] = field(default_factory=list)
    marker_products: List[str] = field(default_factory=list)
    recruited_genes: List[str] = field(default_factory=list)
    score_token: str = ""
    counts_toward_score: bool = False
    machinery: str = ""
    reference: str = ""
    reason: str = ""

    @property
    def all_genes(self) -> List[str]:
        return list(self.marker_genes) + list(self.recruited_genes)


@dataclass
class SpeciesLociResult:
    """Outcome of the species-locus layer for one genome."""

    in_scope: bool = False
    scope_description: str = ""
    calls: List[LocusCall] = field(default_factory=list)
    hits: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=_HIT_COLUMNS))
    score: Optional[int] = None
    score_scheme: str = ""
    score_tokens: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    @property
    def present_calls(self) -> List[LocusCall]:
        return [c for c in self.calls if c.present]

    @property
    def loci_present(self) -> List[str]:
        """Distinct loci called present, in table order.

        A locus can be called more than once in one genome -- NTUH-K2044
        carries both a chromosomal and a plasmid-borne copy of *iro* and of
        *rmpADC* -- so the per-cluster calls are de-duplicated here while
        :attr:`present_calls` retains every cluster.
        """
        return list(dict.fromkeys(c.locus for c in self.present_calls))

    @property
    def n_clusters_present(self) -> int:
        """Number of accepted locus clusters, counting repeated copies."""
        return len(self.present_calls)

    def as_frame(self) -> pd.DataFrame:
        """One row per locus cluster considered, present or not."""
        rows = []
        for c in self.calls:
            rows.append(
                {
                    "locus": c.locus,
                    "locus_label": c.locus_label,
                    "present": c.present,
                    "replicon": c.seqid,
                    "n_marker_genes": c.n_markers,
                    "min_markers_required": c.min_markers,
                    "expected_gene_count": c.expected_gene_count,
                    "annotation_completeness": round(c.completeness, 3),
                    "n_recruited_genes": len(c.recruited_genes),
                    "marker_genes": "; ".join(c.marker_genes),
                    "marker_products": " | ".join(c.marker_products),
                    "recruited_genes": "; ".join(c.recruited_genes),
                    "score_token": c.score_token,
                    "counts_toward_score": c.counts_toward_score,
                    "machinery": c.machinery,
                    "reference": c.reference,
                    "reason": c.reason,
                }
            )
        return pd.DataFrame(
            rows,
            columns=[
                "locus", "locus_label", "present", "replicon", "n_marker_genes",
                "min_markers_required", "expected_gene_count", "annotation_completeness",
                "n_recruited_genes", "marker_genes", "marker_products",
                "recruited_genes", "score_token", "counts_toward_score",
                "machinery", "reference", "reason",
            ],
        )

    def summary_dict(self) -> dict:
        return {
            "species_loci_in_scope": bool(self.in_scope),
            "species_loci_scope": self.scope_description,
            "species_loci_present": list(self.loci_present),
            "species_loci_n_clusters_present": int(self.n_clusters_present),
            "species_loci_score_scheme": self.score_scheme,
            "species_loci_score": self.score,
            "species_loci_score_tokens": list(self.score_tokens),
            "species_loci_n_genes": int(
                0 if self.hits.empty else self.hits["locus_tag"].nunique()
            ),
            "species_loci_n_marker_genes": int(
                0 if self.hits.empty else (self.hits["role"] == "marker").sum()
            ),
        }


_HIT_COLUMNS = [
    "locus_tag", "locus", "locus_label", "role", "matched_field", "matched_text",
    "machinery", "axis_hint", "weight", "score_token", "counts_toward_score",
    "reference",
]


# ---------------------------------------------------------------------------
# detection
# ---------------------------------------------------------------------------
def _gene_order(gene_table: GeneTable) -> Dict[str, List[str]]:
    """Locus tags in coordinate order, per replicon."""
    frame = gene_table.frame
    usable = frame.dropna(subset=["start"]).copy()
    if usable.empty:
        return {}
    usable["seqid"] = usable["seqid"].fillna("unknown").astype(str)
    out: Dict[str, List[str]] = {}
    for seqid, group in usable.sort_values(["seqid", "start"]).groupby("seqid"):
        out[str(seqid)] = [str(t) for t in group["locus_tag"]]
    return out


def _find_markers(
    gene_table: GeneTable, defn: LocusDefinition
) -> Dict[str, Tuple[str, str]]:
    """``locus_tag -> (matched_field, matched_text)`` for one locus definition."""
    found: Dict[str, Tuple[str, str]] = {}

    # (a) constituent gene symbols, through the identifier crosswalk. Inert for
    # genomes whose RefSeq annotation carries no gene symbols, which is exactly
    # the HS11286 situation, but it is the primary route for genomes that do.
    for symbol in defn.gene_symbols:
        tag = gene_table.lookup(symbol)
        if tag and tag not in found:
            found[str(tag)] = ("gene_symbol", symbol)

    # (b) RefSeq product names, matched with the compiled pattern.
    if defn.regex is not None:
        frame = gene_table.frame
        for tag, product, gene in zip(
            frame["locus_tag"].astype(str),
            frame["product"].fillna("").astype(str),
            frame["gene"].fillna("").astype(str),
        ):
            if tag in found:
                continue
            match = defn.regex.search(product)
            if match:
                found[tag] = ("product", product[:200])
                continue
            # The 'gene' field is also worth a pattern match, because some
            # annotations put the symbol there without a product name. Guard
            # against the degenerate case where 'gene' merely repeats the
            # locus tag, which would let a pattern match a tag prefix.
            if gene and gene != tag:
                if defn.regex.search(gene):
                    found[tag] = ("gene", gene[:200])
    return found


def _cluster(indices: Sequence[int], max_gap: int) -> List[List[int]]:
    """Split sorted gene indices wherever more than ``max_gap`` genes intervene."""
    groups: List[List[int]] = []
    current: List[int] = []
    for idx in sorted(indices):
        if current and (idx - current[-1] - 1) > int(max_gap):
            groups.append(current)
            current = [idx]
        else:
            current.append(idx)
    if current:
        groups.append(current)
    return groups


def detect_species_loci(
    gene_table: GeneTable,
    target: TargetInfo,
    resource_dir: Optional[str] = None,
    definitions: Optional[Sequence[LocusDefinition]] = None,
) -> SpeciesLociResult:
    """Detect the species-specific acquired virulence loci of one genome.

    Parameters
    ----------
    gene_table:
        The gene universe, from
        :func:`codonpipe.virulence.resolve.build_crosswalk`.
    target:
        The resolved target; its taxonomy decides which definitions apply.
    resource_dir:
        Override the shipped resource folder.
    definitions:
        Pre-loaded definitions, to avoid re-reading the table.

    Returns
    -------
    SpeciesLociResult
        ``in_scope`` is ``False`` and every other field empty when no
        definition matches the target, which is the intended behaviour for
        genomes outside the configured scopes.
    """
    result = SpeciesLociResult()
    defs = list(definitions) if definitions is not None else load_species_loci(resource_dir)
    if not defs:
        result.notes.append(
            f"No species-locus definitions are available ({SPECIES_LOCI_TABLE} absent "
            "or empty), so the species-locus layer was skipped."
        )
        return result

    scoped = definitions_in_scope(defs, target)
    if not scoped:
        result.notes.append(
            f"None of the {len(defs)} species-locus definition(s) is scoped to "
            f"{target.scientific_name or target.query} (taxon {target.tax_id}; species "
            f"{target.species_tax_id}; genus {target.genus_tax_id}), so the "
            "species-locus layer is a no-op for this genome and contributes no "
            "evidence."
        )
        return result

    result.in_scope = True
    scopes = sorted({f"{d.scope_type}={d.scope_value}" for d in scoped})
    result.scope_description = ", ".join(scopes)
    result.notes.append(
        f"{len(scoped)} species-locus definition(s) are in scope for "
        f"{target.scientific_name or target.query} ({result.scope_description}): "
        + ", ".join(d.locus for d in scoped)
        + "."
    )

    order = _gene_order(gene_table)
    if not order:
        result.notes.append(
            "Genomic coordinates are unavailable for this genome, so locus "
            "co-location could not be assessed and no locus was called."
        )
        return result
    position = {
        tag: (seqid, i) for seqid, tags in order.items() for i, tag in enumerate(tags)
    }

    hit_rows: List[dict] = []
    schemes: Set[str] = set()

    for defn in scoped:
        markers = _find_markers(gene_table, defn)
        if not markers:
            result.calls.append(
                LocusCall(
                    locus=defn.locus, locus_label=defn.locus_label, present=False,
                    seqid="", n_markers=0, min_markers=defn.min_markers,
                    expected_gene_count=defn.expected_gene_count, completeness=0.0,
                    score_token=defn.score_token,
                    counts_toward_score=defn.counts_toward_score,
                    machinery=defn.machinery, reference=defn.reference,
                    reason="no gene of this genome matches the locus gene symbols or "
                           "product-name patterns",
                )
            )
            continue

        # group the markers by replicon, then into co-located clusters
        by_seq: Dict[str, List[int]] = {}
        for tag in markers:
            if tag in position:
                seqid, idx = position[tag]
                by_seq.setdefault(seqid, []).append(idx)

        clusters: List[Tuple[str, List[int]]] = []
        for seqid, idxs in by_seq.items():
            for group in _cluster(idxs, defn.max_gene_gap):
                clusters.append((seqid, group))
        clusters.sort(key=lambda sg: -len(sg[1]))

        if not clusters:
            result.calls.append(
                LocusCall(
                    locus=defn.locus, locus_label=defn.locus_label, present=False,
                    seqid="", n_markers=len(markers), min_markers=defn.min_markers,
                    expected_gene_count=defn.expected_gene_count, completeness=0.0,
                    score_token=defn.score_token,
                    counts_toward_score=defn.counts_toward_score,
                    machinery=defn.machinery, reference=defn.reference,
                    reason="marker genes have no genomic coordinates, so co-location "
                           "could not be assessed",
                )
            )
            continue

        for seqid, group in clusters:
            tags = [order[seqid][i] for i in group]
            n = len(tags)
            present = n >= defn.min_markers
            completeness = (
                n / float(defn.expected_gene_count) if defn.expected_gene_count else 0.0
            )
            recruited: List[str] = []

            if present and defn.extend_locus:
                # Reuse the module's own operon/island machinery, seeded only
                # with this cluster's markers, then clamp the result to the
                # cluster span so recruitment cannot propagate along the
                # replicon.
                ctx = extend_by_context(
                    gene_table,
                    tags,
                    operon_extension=True,
                    max_intergenic_gap=50,
                    require_same_strand=True,
                    island_bridging=True,
                    island_bridge_genes=defn.max_gene_gap,
                    min_island_seeds=defn.min_markers,
                )
                lo, hi = group[0] - defn.max_gene_gap, group[-1] + defn.max_gene_gap
                for tag in ctx.added:
                    where = position.get(str(tag))
                    if where and where[0] == seqid and lo <= where[1] <= hi:
                        recruited.append(str(tag))
                recruited = sorted(set(recruited), key=lambda t: position[t][1])

            call = LocusCall(
                locus=defn.locus,
                locus_label=defn.locus_label,
                present=present,
                seqid=seqid,
                n_markers=n,
                min_markers=defn.min_markers,
                expected_gene_count=defn.expected_gene_count,
                completeness=completeness,
                marker_genes=tags,
                marker_products=[markers[t][1] for t in tags],
                recruited_genes=recruited,
                score_token=defn.score_token,
                counts_toward_score=defn.counts_toward_score,
                machinery=defn.machinery,
                reference=defn.reference,
                reason=(
                    f"{n} co-located marker gene(s) on {seqid}, at or above the "
                    f"minimum of {defn.min_markers}"
                    if present
                    else f"only {n} co-located marker gene(s) on {seqid}, below the "
                         f"minimum of {defn.min_markers} required to call this locus; "
                         "reported for transparency and contributing no evidence"
                ),
            )
            result.calls.append(call)

            if not present:
                continue
            if defn.score_scheme:
                schemes.add(defn.score_scheme)
            for tag in tags:
                hit_rows.append(
                    {
                        "locus_tag": tag, "locus": defn.locus,
                        "locus_label": defn.locus_label, "role": "marker",
                        "matched_field": markers[tag][0],
                        "matched_text": markers[tag][1],
                        "machinery": defn.machinery, "axis_hint": defn.axis_hint,
                        "weight": defn.weight, "score_token": defn.score_token,
                        "counts_toward_score": defn.counts_toward_score,
                        "reference": defn.reference,
                    }
                )
            for tag in recruited:
                hit_rows.append(
                    {
                        "locus_tag": tag, "locus": defn.locus,
                        "locus_label": defn.locus_label, "role": "locus_context",
                        "matched_field": "genomic_context",
                        "matched_text": "recruited within the detected locus block",
                        "machinery": defn.machinery, "axis_hint": defn.axis_hint,
                        "weight": defn.weight, "score_token": defn.score_token,
                        "counts_toward_score": defn.counts_toward_score,
                        "reference": defn.reference,
                    }
                )

    frame = pd.DataFrame(hit_rows, columns=_HIT_COLUMNS)
    if not frame.empty:
        # One gene can sit in two overlapping locus definitions; keep the
        # marker row in preference to the recruited one.
        frame["_rank"] = (frame["role"] == "marker").map({True: 0, False: 1})
        frame = (
            frame.sort_values(["locus_tag", "_rank"])
            .drop_duplicates(["locus_tag", "locus"], keep="first")
            .drop(columns="_rank")
            .reset_index(drop=True)
        )
    result.hits = frame

    tokens = {
        c.score_token
        for c in result.present_calls
        if c.counts_toward_score and c.score_token
    }
    result.score_tokens = sorted(tokens)
    if schemes:
        scheme = sorted(schemes)[0]
        if len(schemes) > 1:
            result.notes.append(
                f"Several scoring schemes are in scope ({', '.join(sorted(schemes))}); "
                f"{scheme} was applied."
            )
        result.score_scheme = scheme
        fn = SCORE_SCHEMES.get(scheme)
        if fn is None:
            result.notes.append(
                f"Scoring scheme {scheme!r} is named in {SPECIES_LOCI_TABLE} but is not "
                "registered in SCORE_SCHEMES, so no score was computed."
            )
        else:
            result.score = int(fn(tokens))
    else:
        # In scope, but nothing present: the score of the scheme declared by
        # the in-scope definitions still applies, and it is zero.
        declared = sorted({d.score_scheme for d in scoped if d.score_scheme})
        if declared:
            result.score_scheme = declared[0]
            fn = SCORE_SCHEMES.get(declared[0])
            result.score = int(fn(set())) if fn else None

    present = result.present_calls
    if present:
        result.notes.append(
            "Loci called present: "
            + "; ".join(
                f"{c.locus} ({c.n_markers}/{c.expected_gene_count} genes named by the "
                f"annotation, {len(c.recruited_genes)} recruited by context, on {c.seqid})"
                for c in present
            )
            + "."
        )
    else:
        result.notes.append("No species-specific acquired virulence locus was called present.")
    if result.score is not None:
        result.notes.append(
            f"{result.score_scheme} score = {result.score} "
            f"(scoring loci present: {', '.join(result.score_tokens) or 'none'}). "
            "Computed from published locus definitions applied to the RefSeq "
            "annotation; this is NOT output of the Kleborate program and involves no "
            "MLST, no Kaptive K/O typing and no sequence-identity scoring."
        )
    return result
