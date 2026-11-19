#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Evidence harvesting for automated virulence-gene discovery.

Five independent evidence layers are collected, each as a flat table of
``EvidenceItem``-shaped rows so that every call in the final output can be
traced back to a database record:

``E1`` curated virulence-factor calls
    Retrieved from the BV-BRC specialty-gene service, which aggregates
    **VFDB**, **VICTORS** and **PATRIC_VF** and records, per gene, the source
    database, the source accession, the virulence-factor category, the percent
    identity, the query coverage, the e-value and whether the assignment rests
    on literature curation or on sequence homology.  This is the automated
    equivalent of screening a genome against those databases, with the
    advantage that the thresholds and provenance are fixed and citable rather
    than chosen per run.

``E2`` UniProtKB curation
    Keyword, Gene Ontology and subcellular-location annotation for the
    proteome, with the review status and the evidence codes retained.  Two
    scopes are used: the strain itself (mapped by RefSeq protein accession,
    exact) and, optionally, reviewed entries of the parent species (mapped by
    gene symbol, flagged as a cross-strain transfer) -- UniProt curation is
    concentrated in a handful of reference strains, so the second scope adds
    most of the experimentally supported annotation.

``E3`` Gene Ontology closures
    Descendant sets of the host-interaction and core-cellular-process roots,
    fetched from the EBI QuickGO service and validated against the live
    ontology.  Validation is not cosmetic: the widely used term ``GO:0009405``
    ("pathogenesis") is obsolete in current GO releases and silently returns
    nothing.

``E4`` controlled vocabulary
    Low-weight corroborating matches of RefSeq product names and gene symbols
    against a shipped, editable lexicon, which also supplies the functional
    subcategories used as CodonPipe cluster columns.

``E5`` genomic context
    Operon-level propagation of virulence evidence, implemented in
    :mod:`codonpipe.virulence.context`.

References
----------
Chen et al. (2005) VFDB: a reference database for bacterial virulence
factors. *Nucleic Acids Res* 33:D325-D328.
Sayers et al. (2019) Victors: a web-based knowledge base of virulence factors
in human and animal pathogens. *Nucleic Acids Res* 47:D693-D700.
Olson et al. (2023) Introducing the Bacterial and Viral Bioinformatics
Resource Center (BV-BRC). *Nucleic Acids Res* 51:D678-D689.
The UniProt Consortium (2025) UniProt: the Universal Protein Knowledgebase.
Binns et al. (2009) QuickGO: a web-based tool for Gene Ontology searching.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import pandas as pd

from .resolve import GeneTable, TargetInfo
from .webcache import WebCache, OfflineCacheMiss

__all__ = [
    "RESOURCE_DIR",
    "load_resource_table",
    "load_go_axes",
    "load_keyword_axes",
    "load_vf_category_axes",
    "load_vocabulary",
    "load_core_vocabulary",
    "fetch_go_closures",
    "bvbrc_genome_ids",
    "fetch_bvbrc_virulence",
    "fetch_uniprot_table",
    "map_bvbrc_to_genes",
    "map_uniprot_to_genes",
    "apply_vocabulary",
]

RESOURCE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources")

BVBRC = "https://www.bv-brc.org/api"
UNIPROT = "https://rest.uniprot.org"
QUICKGO = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms"

#: BV-BRC spells the virulence property two ways in its own data.
VF_PROPERTIES = ("Virulence Factor", "Virulance factor")

#: Fields pulled from UniProtKB.  Matching is done on stable identifiers
#: (accessions and keyword ids), never on display names.
UNIPROT_FIELDS = (
    "accession,reviewed,id,protein_name,gene_names,gene_primary,gene_oln,"
    "xref_refseq,keywordid,go_id,cc_subcellular_location,protein_existence,organism_id"
)


# ---------------------------------------------------------------------------
# resource files
# ---------------------------------------------------------------------------
def load_resource_table(name: str, directory: Optional[str] = None) -> pd.DataFrame:
    """Read one of the shipped, user-editable TSV resources.

    Lines beginning with ``#`` are comments.  A missing file is an error: the
    classification rules are defined by these tables and silently running
    without them would produce results that cannot be interpreted.
    """
    folder = directory or RESOURCE_DIR
    path = os.path.join(folder, name)
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"Required virulence-discovery resource {name!r} is missing from {folder}."
        )
    return pd.read_csv(path, sep="\t", comment="#", dtype=str).fillna("")


def _as_float(value, default: float = 0.0) -> float:
    try:
        return float(value)
    except Exception:
        return default


def load_go_axes(directory: Optional[str] = None) -> pd.DataFrame:
    frame = load_resource_table("go_axes.tsv", directory)
    frame["weight"] = frame["weight"].map(lambda v: _as_float(v, 1.0))
    frame["enabled"] = frame["enabled"].map(lambda v: str(v).strip() not in ("0", "", "false", "False"))
    return frame


def load_keyword_axes(directory: Optional[str] = None) -> pd.DataFrame:
    frame = load_resource_table("uniprot_keyword_axes.tsv", directory)
    frame["weight"] = frame["weight"].map(lambda v: _as_float(v, 1.0))
    return frame


def load_vf_category_axes(directory: Optional[str] = None) -> Dict[str, str]:
    frame = load_resource_table("vf_category_axes.tsv", directory)
    return {str(r.term).strip().lower(): str(r.axis).strip() for r in frame.itertuples()}


def _compile_patterns(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    compiled = []
    for pattern in frame["pattern"]:
        try:
            compiled.append(re.compile(str(pattern), re.IGNORECASE))
        except re.error as exc:
            raise ValueError(f"Invalid regular expression in resource file: {pattern!r} ({exc})")
    frame["regex"] = compiled
    return frame


def restrict_vocabulary_to_taxon(
    frame: pd.DataFrame, taxon_ids: Iterable[Optional[int]]
) -> Tuple[pd.DataFrame, List[str]]:
    """Drop vocabulary rows that name a taxon the target does not belong to.

    Some rules are gene-symbol patterns for a named pathogenicity island and are
    meaningful in one genus only: `sseA`-`sseL` is SPI-2 in *Salmonella*, but the
    same symbols occur in other enterobacteria without implying that island.
    Applying them everywhere mislabels clusters in the other species. A row may
    therefore declare `restrict_taxon`, a semicolon-separated list of NCBI taxon
    ids; it is kept only when one of them appears in the target's own lineage.
    Rows with an empty `restrict_taxon` are universal and always kept.
    """
    if "restrict_taxon" not in frame.columns:
        return frame, []
    keep_ids = {int(t) for t in taxon_ids if t}
    notes: List[str] = []
    mask = []
    dropped: List[str] = []
    for row in frame.itertuples():
        raw = str(getattr(row, "restrict_taxon", "") or "").strip()
        if not raw or raw.lower() == "nan":
            mask.append(True)
            continue
        wanted = {int(t) for t in re.split(r"[;,\s]+", raw) if t.strip().isdigit()}
        ok = bool(wanted & keep_ids)
        mask.append(ok)
        if not ok:
            dropped.append(str(getattr(row, "subcategory", "")))
    if dropped:
        notes.append(
            "Vocabulary rules restricted to another taxon were not applied to this "
            "genome: " + ", ".join(sorted(set(dropped)))
            + ". They are gene-symbol patterns tied to one lineage's "
            "nomenclature and would mislabel clusters here: the same symbols "
            "name unrelated systems elsewhere, as esp does for the ESX-1 "
            "secretion system in mycobacteria."
        )
    return frame[pd.Series(mask, index=frame.index)].copy(), notes


def load_vocabulary(directory: Optional[str] = None) -> pd.DataFrame:
    frame = load_resource_table("virulence_vocabulary.tsv", directory)
    frame["weight"] = frame["weight"].map(lambda v: _as_float(v, 0.5))
    if "as_cluster" not in frame.columns:
        frame["as_cluster"] = "1"
    frame["as_cluster"] = frame["as_cluster"].map(
        lambda v: str(v).strip() not in ("0", "", "false", "False", "no")
    )
    if "machinery" not in frame.columns:
        frame["machinery"] = "unclassified"
    frame["machinery"] = frame["machinery"].map(
        lambda v: str(v).strip().lower() or "unclassified"
    )
    return _compile_patterns(frame)


def load_core_vocabulary(directory: Optional[str] = None) -> pd.DataFrame:
    return _compile_patterns(load_resource_table("core_process_vocabulary.tsv", directory))


# ---------------------------------------------------------------------------
# E3 -- Gene Ontology closures
# ---------------------------------------------------------------------------
def fetch_go_closures(
    go_axes: pd.DataFrame, web: WebCache
) -> Tuple[Dict[str, Set[str]], List[str]]:
    """Expand each enabled GO root to its descendant closure.

    Returns ``(closures, warnings)`` where ``closures`` maps an axis name to
    the set of GO identifiers belonging to it (roots included) and
    ``warnings`` lists roots that are obsolete or unknown upstream.
    """
    closures: Dict[str, Set[str]] = {}
    warnings: List[str] = []
    roots = go_axes[go_axes["enabled"]]
    if roots.empty:
        return closures, ["No Gene Ontology roots are enabled in go_axes.tsv."]

    ids = sorted({str(x).strip() for x in roots["go_id"] if str(x).strip()})
    # Validate the roots first so obsolete terms are reported rather than
    # silently contributing an empty descendant set.
    try:
        meta = web.get_json(f"{QUICKGO}/{','.join(ids)}")
        for res in (meta or {}).get("results", []):
            if res.get("isObsolete"):
                warnings.append(
                    f"{res.get('id')} ('{res.get('name')}') is obsolete in the current "
                    "Gene Ontology release and contributes no terms."
                )
        seen = {r.get("id") for r in (meta or {}).get("results", [])}
        for missing in set(ids) - seen:
            warnings.append(f"{missing} is not present in the current Gene Ontology release.")
    except OfflineCacheMiss:
        warnings.append("Offline mode: Gene Ontology roots were not revalidated.")
    except Exception as exc:
        warnings.append(f"Gene Ontology root validation failed ({exc}); roots used as given.")

    for axis, group in roots.groupby("axis"):
        terms: Set[str] = set()
        group_ids = sorted({str(x).strip() for x in group["go_id"] if str(x).strip()})
        for start in range(0, len(group_ids), 20):
            chunk = group_ids[start : start + 20]
            try:
                payload = web.get_json(
                    f"{QUICKGO}/{','.join(chunk)}/descendants",
                    {"relations": "is_a,part_of"},
                )
            except OfflineCacheMiss:
                warnings.append(
                    f"Offline mode: descendants of {', '.join(chunk)} are not cached; "
                    f"the {axis} axis is incomplete."
                )
                terms.update(chunk)
                continue
            except Exception as exc:
                warnings.append(f"Could not expand {', '.join(chunk)} ({exc}).")
                terms.update(chunk)
                continue
            for res in (payload or {}).get("results", []):
                if res.get("id"):
                    terms.add(res["id"])
                for child in res.get("descendants") or []:
                    terms.add(child)
        closures[str(axis)] = terms
    return closures, warnings


# ---------------------------------------------------------------------------
# E1 -- BV-BRC curated virulence factors
# ---------------------------------------------------------------------------
def _bvbrc_query(core: str, query: str, web: WebCache) -> list:
    """Run one BV-BRC Data API query.

    BV-BRC expects RQL in the raw query string (``eq(...)``, ``in(...)``,
    ``select(...)``), which must reach the server unescaped, so the query is
    appended to the URL rather than passed as encoded parameters.
    """
    full = f"{BVBRC}/{core}/?{query}&http_accept=application/json"
    data = web.get_json(full)
    return data if isinstance(data, list) else []


def bvbrc_genome_ids(target: TargetInfo, web: WebCache) -> Tuple[List[str], List[str]]:
    """Find the BV-BRC genome records for the resolved assembly.

    BV-BRC indexes both GenBank (``GCA_``) and RefSeq (``GCF_``) accessions,
    and frequently splits one assembly into several records (chromosome plus
    plasmids).  All matching records are returned so that plasmid-encoded
    virulence factors are not missed.

    Both accessions are queried and their results pooled.  Returning on the
    first accession that matched anything is not safe: an assembly can carry a
    small, contig-level record under one accession and its full annotation
    under the other, in which case the first-match rule silently scores the
    genome on a fragment.  *K. pneumoniae* NTUH-K2044 is the worked example --
    its GenBank accession matches a 293-CDS plasmid-only record while the
    5,386-CDS record carries the RefSeq accession.
    """
    notes: List[str] = []
    accessions = [a for a in (target.genbank_accession, target.refseq_accession) if a]
    pooled: Dict[str, dict] = {}
    for acc in accessions:
        try:
            rows = _bvbrc_query(
                "genome",
                f"eq(assembly_accession,{acc})"
                "&select(genome_id,genome_name,genome_status,patric_cds,assembly_accession)"
                "&limit(200)",
                web,
            )
        except Exception as exc:
            notes.append(f"BV-BRC genome lookup for {acc} failed ({exc}).")
            continue
        for row in rows:
            gid = str(row.get("genome_id") or "")
            if gid:
                pooled.setdefault(gid, row)
    if pooled:
        records = sorted(
            pooled.values(),
            key=lambda r: (
                str(r.get("genome_status", "")).lower() != "complete",
                -int(r.get("patric_cds") or 0),
            ),
        )
        ids = [str(r["genome_id"]) for r in records]
        richest = records[0]
        detail = ", ".join(
            f"{r['genome_id']} ({int(r.get('patric_cds') or 0)} CDS,"
            f" {str(r.get('genome_status') or 'status unknown').lower()})"
            for r in records
        )
        notes.append(
            f"BV-BRC genome record(s) matched on assembly accession(s) "
            f"{', '.join(accessions)}: {detail}."
        )
        # A pool whose richest member holds very few genes means the assembly
        # is represented in BV-BRC only by a fragment, and the curated
        # evidence for this genome will be correspondingly partial. Say so
        # rather than letting a near-empty virulence table look like a
        # biological result.
        if int(richest.get("patric_cds") or 0) < 500:
            notes.append(
                f"The richest BV-BRC record for this assembly holds only "
                f"{int(richest.get('patric_cds') or 0)} CDS, so curated "
                f"virulence evidence will cover a fraction of the genome. "
                f"Treat the resulting counts as a lower bound and do not "
                f"compare them with a genome that has a complete record."
            )
        return ids, notes

    if target.tax_id:
        try:
            rows = _bvbrc_query(
                "genome",
                f"eq(taxon_id,{target.tax_id})"
                "&select(genome_id,genome_name,genome_status,patric_cds)"
                "&limit(200)",
                web,
            )
        except Exception as exc:
            notes.append(f"BV-BRC taxon lookup failed ({exc}).")
            rows = []
        complete = [r for r in rows if str(r.get("genome_status", "")).lower() == "complete"]
        pool = complete or rows
        if pool:
            pool.sort(key=lambda r: -int(r.get("patric_cds") or 0))
            best = str(pool[0]["genome_id"])
            notes.append(
                f"No BV-BRC record carried assembly accession {accessions}; fell back to "
                f"genome {best} ({pool[0].get('genome_name')}) matched on taxon "
                f"{target.tax_id}. Identifier mapping may be less complete."
            )
            return [best], notes

    notes.append(
        "No BV-BRC genome record could be matched, so VFDB/VICTORS/PATRIC_VF "
        "evidence is unavailable for this target."
    )
    return [], notes


def fetch_bvbrc_virulence(genome_ids: Sequence[str], web: WebCache) -> pd.DataFrame:
    """Retrieve all virulence-factor specialty-gene records for ``genome_ids``."""
    props = ",".join(f'"{p}"' for p in VF_PROPERTIES)
    select = (
        "select(genome_id,patric_id,refseq_locus_tag,alt_locus_tag,gene,product,property,"
        "source,source_id,evidence,classification,organism,identity,query_coverage,"
        "subject_coverage,e_value)"
    )
    rows: List[dict] = []
    for gid in genome_ids:
        rows.extend(
            _bvbrc_query(
                "sp_gene",
                f"in(property,({props}))&eq(genome_id,{gid})&{select}&limit(25000)",
                web,
            )
        )
    if not rows:
        return pd.DataFrame(
            columns=[
                "genome_id", "patric_id", "refseq_locus_tag", "alt_locus_tag", "gene",
                "product", "property", "source", "source_id", "evidence",
                "classification", "organism", "identity", "query_coverage",
                "subject_coverage", "e_value",
            ]
        )
    frame = pd.DataFrame(rows)
    if "classification" in frame.columns:
        frame["classification"] = frame["classification"].map(
            lambda v: "; ".join(str(x) for x in v) if isinstance(v, (list, tuple)) else (v or "")
        )
    for col in ("identity", "query_coverage", "subject_coverage"):
        if col in frame.columns:
            frame[col] = pd.to_numeric(frame[col], errors="coerce")
    return frame


def _bvbrc_feature_map(
    genome_ids: Sequence[str], web: WebCache, notes: Optional[List[str]] = None
) -> Dict[str, dict]:
    """``patric_id`` -> feature record, used to rescue unmapped VF calls.

    A failure here used to be swallowed, returning an empty map. The caller
    then reported the records it could not place as an "annotation mismatch
    between the database genome and the RefSeq assembly" -- a biological
    explanation for what may have been a network error. Failures are now
    reported through ``notes`` so the two cannot be confused.
    """
    out: Dict[str, dict] = {}
    failed: List[str] = []
    for gid in genome_ids:
        try:
            rows = _bvbrc_query(
                "genome_feature",
                f"eq(genome_id,{gid})&eq(feature_type,CDS)"
                "&select(patric_id,refseq_locus_tag,alt_locus_tag,gene,product,protein_id,"
                "start,end,strand,accession)&limit(50000)",
                web,
            )
        except Exception as exc:
            failed.append(f"{gid} ({type(exc).__name__})")
            continue
        for row in rows:
            if row.get("patric_id"):
                out[str(row["patric_id"])] = row
    if failed and notes is not None:
        notes.append(
            "The BV-BRC feature lookup used to place unmapped virulence-factor "
            f"records failed for genome(s) {', '.join(failed)}. Any records "
            "reported below as unmapped may be unmapped because of that failure "
            "rather than because of an annotation mismatch."
        )
    return out


def map_bvbrc_to_genes(
    vf_frame: pd.DataFrame,
    gene_table: GeneTable,
    genome_ids: Sequence[str],
    web: WebCache,
    min_identity: float = 80.0,
    min_coverage: float = 60.0,
) -> Tuple[pd.DataFrame, List[str]]:
    """Attach a canonical ``locus_tag`` to every virulence-factor record.

    Thresholds follow the convention used for ABRicate-style screens
    (>=80% nucleotide/protein identity, >=60% coverage).  Literature-curated
    assignments carry no alignment statistics and are never filtered out by
    them.
    """
    notes: List[str] = []
    if vf_frame.empty:
        return vf_frame.assign(locus_tag=pd.Series(dtype=str), mapped_via=pd.Series(dtype=str)), notes

    frame = vf_frame.copy()
    locus: List[Optional[str]] = []
    via: List[str] = []
    for row in frame.itertuples(index=False):
        hit = None
        source = ""
        for cand, label in (
            (getattr(row, "refseq_locus_tag", None), "refseq_locus_tag"),
            (getattr(row, "alt_locus_tag", None), "alt_locus_tag"),
            (getattr(row, "gene", None), "gene_symbol"),
        ):
            hit = gene_table.lookup(cand)
            if hit:
                source = label
                break
        locus.append(hit)
        via.append(source)
    frame["locus_tag"] = locus
    frame["mapped_via"] = via

    unmapped = frame["locus_tag"].isna()
    if unmapped.any():
        features = _bvbrc_feature_map(genome_ids, web, notes)
        if features:
            for idx in frame.index[unmapped]:
                feat = features.get(str(frame.at[idx, "patric_id"]))
                if not feat:
                    continue
                for cand, label in (
                    (feat.get("protein_id"), "protein_id_via_feature"),
                    (feat.get("refseq_locus_tag"), "refseq_locus_tag_via_feature"),
                    (feat.get("gene"), "gene_symbol_via_feature"),
                ):
                    hit = gene_table.lookup(cand)
                    if hit:
                        frame.at[idx, "locus_tag"] = hit
                        frame.at[idx, "mapped_via"] = label
                        break
        still = int(frame["locus_tag"].isna().sum())
        if still:
            notes.append(
                f"{still} of {len(frame)} virulence-factor records could not be mapped onto "
                "a gene of the analysed genome and were discarded (annotation mismatch "
                "between the database genome and the RefSeq assembly)."
            )
    frame = frame.dropna(subset=["locus_tag"]).copy()

    # Alignment-quality gate, applied only where statistics exist.
    is_literature = frame["evidence"].astype(str).str.lower().str.contains("literature")
    identity = pd.to_numeric(frame.get("identity"), errors="coerce")
    coverage = pd.to_numeric(frame.get("query_coverage"), errors="coerce")
    quality_ok = (
        is_literature
        | ((identity >= float(min_identity)) & (coverage >= float(min_coverage)))
        | (identity.isna() & coverage.isna())
    )
    dropped = int((~quality_ok).sum())
    if dropped:
        notes.append(
            f"{dropped} homology-based virulence-factor records fell below the "
            f"{min_identity:.0f}% identity / {min_coverage:.0f}% coverage threshold and were discarded."
        )
    frame = frame[quality_ok].copy()
    frame["is_literature"] = is_literature[quality_ok].values
    return frame, notes


# ---------------------------------------------------------------------------
# E2 -- UniProtKB
# ---------------------------------------------------------------------------
def fetch_uniprot_table(query: str, web: WebCache) -> pd.DataFrame:
    """Stream a UniProtKB query as a table.

    The ``stream`` endpoint returns the complete result set in one response,
    which keeps the cache a single file per query and avoids cursor paging.
    """
    import io as _io

    text = web.get_text(
        f"{UNIPROT}/uniprotkb/stream",
        {"query": query, "format": "tsv", "fields": UNIPROT_FIELDS},
        ext=".tsv",
    )
    if not text.strip():
        return pd.DataFrame()
    frame = pd.read_csv(_io.StringIO(text), sep="\t", dtype=str).fillna("")
    frame.columns = [c.strip() for c in frame.columns]
    return frame


_UNIPROT_COLS = {
    "Entry": "uniprot",
    "Reviewed": "reviewed",
    "Entry Name": "entry_name",
    "Protein names": "protein_name",
    "Gene Names": "gene_names",
    "Gene Names (primary)": "gene_primary",
    "Gene Names (ordered locus)": "gene_oln",
    "RefSeq": "refseq",
    "Keyword ID": "keyword_ids",
    "Gene Ontology IDs": "go_ids",
    "Subcellular location [CC]": "subcellular",
    "Protein existence": "protein_existence",
    "Organism (ID)": "organism_id",
}


def _normalise_uniprot(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(columns=list(_UNIPROT_COLS.values()))
    out = frame.rename(columns=_UNIPROT_COLS)
    for col in _UNIPROT_COLS.values():
        if col not in out.columns:
            out[col] = ""
    return out[list(_UNIPROT_COLS.values())].copy()


def map_uniprot_to_genes(
    target: TargetInfo,
    gene_table: GeneTable,
    web: WebCache,
    include_species_transfer: bool = True,
) -> Tuple[pd.DataFrame, List[str]]:
    """Collect UniProtKB annotation and attach canonical locus tags.

    Strain-scope entries are mapped by RefSeq protein accession, old locus tag
    or gene symbol.  Species-scope entries -- used only when the strain itself
    is sparsely curated -- are reviewed (Swiss-Prot) records mapped by gene
    symbol only, and are labelled ``species_reviewed`` so that the transfer is
    visible in the audit table and can be switched off.
    """
    notes: List[str] = []
    collected: List[pd.DataFrame] = []

    if target.tax_id:
        try:
            strain = _normalise_uniprot(
                fetch_uniprot_table(f"organism_id:{target.tax_id}", web)
            )
        except Exception as exc:
            strain = pd.DataFrame()
            notes.append(f"UniProtKB strain query failed ({exc}).")
        if not strain.empty:
            strain["scope"] = "strain"
            collected.append(strain)
            n_reviewed = int((strain["reviewed"].str.lower() == "reviewed").sum())
            notes.append(
                f"UniProtKB returned {len(strain)} entries for taxon {target.tax_id}, "
                f"{n_reviewed} of them reviewed."
            )

    if include_species_transfer and target.species_tax_id and target.species_tax_id != target.tax_id:
        try:
            species = _normalise_uniprot(
                fetch_uniprot_table(
                    f"taxonomy_id:{target.species_tax_id} AND reviewed:true", web
                )
            )
        except Exception as exc:
            species = pd.DataFrame()
            notes.append(f"UniProtKB species query failed ({exc}).")
        if not species.empty:
            species["scope"] = "species_reviewed"
            collected.append(species)
            notes.append(
                f"{len(species)} reviewed UniProtKB entries for {target.species_name} "
                f"(taxon {target.species_tax_id}) were available for cross-strain "
                "annotation transfer by gene symbol."
            )

    if not collected:
        return pd.DataFrame(
            columns=list(_UNIPROT_COLS.values()) + ["scope", "locus_tag", "mapped_via"]
        ), notes

    frame = pd.concat(collected, ignore_index=True)
    locus: List[Optional[str]] = []
    via: List[str] = []
    for row in frame.itertuples(index=False):
        hit = None
        label = ""
        if row.scope == "strain":
            for accession in str(row.refseq).replace(";", " ").split():
                hit = gene_table.lookup(accession.strip())
                if hit:
                    label = "refseq_protein"
                    break
            if not hit:
                for cand, lab in (
                    (row.gene_oln, "old_locus_tag"),
                    (row.gene_primary, "gene_symbol"),
                ):
                    hit = gene_table.lookup(cand)
                    if hit:
                        label = lab
                        break
        else:
            for cand in [row.gene_primary] + str(row.gene_names).split():
                hit = gene_table.lookup(cand)
                if hit:
                    label = "gene_symbol_species_transfer"
                    break
        locus.append(hit)
        via.append(label)
    frame["locus_tag"] = locus
    frame["mapped_via"] = via

    mapped = frame.dropna(subset=["locus_tag"]).copy()
    # A strain-scope mapping always wins over a species-scope transfer.
    mapped["_rank"] = (mapped["scope"] == "strain").map({True: 0, False: 1})
    mapped = (
        mapped.sort_values(["locus_tag", "_rank"])
        .drop_duplicates(["locus_tag", "uniprot"], keep="first")
        .drop(columns="_rank")
    )
    notes.append(
        f"{mapped['locus_tag'].nunique()} genes carry UniProtKB annotation "
        f"({int((mapped['scope'] == 'strain').sum())} strain-scope and "
        f"{int((mapped['scope'] == 'species_reviewed').sum())} transferred entries)."
    )
    return mapped, notes


# ---------------------------------------------------------------------------
# E4 -- controlled vocabulary
# ---------------------------------------------------------------------------
def apply_vocabulary(
    gene_table: GeneTable, vocabulary: pd.DataFrame
) -> pd.DataFrame:
    """Match product names and gene symbols against the shipped lexicon.

    Returns one row per (gene, vocabulary rule) hit with the subcategory, the
    suggested axis and the rule weight.
    """
    rows: List[dict] = []
    frame = gene_table.frame
    products = frame["product"].fillna("").astype(str).tolist()
    genes = frame["gene"].fillna("").astype(str).tolist()
    tags = frame["locus_tag"].astype(str).tolist()
    for rule in vocabulary.itertuples(index=False):
        field_name = str(rule.field).strip().lower()
        regex = rule.regex
        if field_name == "gene":
            values = genes
        elif field_name == "product":
            values = products
        else:
            values = [f"{g} {p}" for g, p in zip(genes, products)]
        for tag, value in zip(tags, values):
            if value and regex.search(value):
                rows.append(
                    {
                        "locus_tag": tag,
                        "subcategory": rule.subcategory,
                        "field": field_name,
                        "pattern": rule.pattern,
                        "axis_hint": rule.axis_hint,
                        "weight": float(rule.weight),
                        "as_cluster": bool(getattr(rule, "as_cluster", True)),
                        "machinery": str(getattr(rule, "machinery", "unclassified")),
                        "matched_text": value[:200],
                    }
                )
    return pd.DataFrame(
        rows,
        columns=[
            "locus_tag", "subcategory", "field", "pattern",
            "axis_hint", "weight", "as_cluster", "machinery", "matched_text",
        ],
    )


def apply_core_vocabulary(
    gene_table: GeneTable, core_vocabulary: pd.DataFrame
) -> pd.DataFrame:
    """Flag genes whose annotation indicates a core cellular function.

    This is the automated counterpart of the manual removal of housekeeping
    and general stress-response proteins that curated virulence screens
    require, and it is what makes Falkow's second postulate interpretable for
    the remaining genes.
    """
    rows: List[dict] = []
    frame = gene_table.frame
    products = frame["product"].fillna("").astype(str).tolist()
    genes = frame["gene"].fillna("").astype(str).tolist()
    tags = frame["locus_tag"].astype(str).tolist()
    for rule in core_vocabulary.itertuples(index=False):
        field_name = str(rule.field).strip().lower()
        if field_name == "gene":
            values = genes
        elif field_name == "product":
            values = products
        else:
            values = [f"{g} {p}" for g, p in zip(genes, products)]
        for tag, value in zip(tags, values):
            if value and rule.regex.search(value):
                rows.append(
                    {
                        "locus_tag": tag,
                        "field": field_name,
                        "pattern": rule.pattern,
                        "note": getattr(rule, "note", ""),
                        "matched_text": value[:200],
                    }
                )
    return pd.DataFrame(
        rows, columns=["locus_tag", "field", "pattern", "note", "matched_text"]
    )
