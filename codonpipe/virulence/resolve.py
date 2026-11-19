#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Target resolution and gene-identifier crosswalks.

The discovery pipeline accepts a free-form target specification -- a species
or strain name, an NCBI taxonomy identifier, or a genome assembly accession --
and must turn it into (i) an unambiguous NCBI taxon, (ii) a specific RefSeq
assembly, and (iii) a gene table whose identifiers match the coding-sequence
FASTA that CodonPipe analyses.

Identifier handling is the quiet source of most errors in this kind of
pipeline.  RefSeq CDS FASTA headers carry current ``locus_tag`` values
(e.g. ``SL1344_RS01405``) and RefSeq protein accessions (``WP_...``), whereas
the virulence databases were largely built against the original submitter
annotation and report the *old* locus tags (``SL1344_0276``) or GenBank
protein accessions.  :func:`build_crosswalk` reads both namespaces out of the
RefSeq GFF and exposes lookups in every direction, so evidence can be mapped
onto the exact genes CodonPipe will score.

All network access is mediated by :class:`~codonpipe.virulence.webcache.WebCache`.
"""

from __future__ import annotations

import io
import os
import re
import zipfile
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Tuple, Set

import pandas as pd

from .webcache import WebCache, NetworkUnavailable

__all__ = [
    "TargetInfo",
    "GeneTable",
    "resolve_target",
    "fetch_refseq_gff",
    "fetch_cds_fasta",
    "assembly_from_cds_fasta",
    "nucleotide_accessions_from_fasta",
    "parse_gff_genes",
    "parse_cds_fasta",
    "build_crosswalk",
    "ASSEMBLY_RE",
]

DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

ASSEMBLY_RE = re.compile(r"^GC[AF]_\d{9}(?:\.\d+)?$", re.IGNORECASE)
TAXID_RE = re.compile(r"^(?:txid[:=]?|taxid[:=]?|tax_id[:=]?)?(\d+)$", re.IGNORECASE)

#: Preference order when several assemblies exist for a taxon.
_LEVEL_RANK = {
    "complete genome": 0,
    "chromosome": 1,
    "scaffold": 2,
    "contig": 3,
}


# ---------------------------------------------------------------------------
# data containers
# ---------------------------------------------------------------------------
@dataclass
class TargetInfo:
    """Everything the pipeline needs to know about the resolved target."""

    query: str
    tax_id: Optional[int] = None
    scientific_name: Optional[str] = None
    rank: Optional[str] = None
    lineage_names: List[str] = field(default_factory=list)
    #: Every ancestor taxon id, root to strain. Vocabulary rules may be
    #: restricted to a taxon at any rank (class, order, family, genus), so the
    #: whole chain is needed, not just species and genus.
    lineage_tax_ids: List[int] = field(default_factory=list)
    species_tax_id: Optional[int] = None
    species_name: Optional[str] = None
    genus_tax_id: Optional[int] = None
    genus_name: Optional[str] = None
    assembly_accession: Optional[str] = None          # GCF_... preferred
    paired_accession: Optional[str] = None            # matching GCA_...
    assembly_name: Optional[str] = None
    assembly_level: Optional[str] = None
    refseq_category: Optional[str] = None
    organism_label: Optional[str] = None              # organism + strain
    strain: Optional[str] = None
    resolution_notes: List[str] = field(default_factory=list)

    @property
    def genbank_accession(self) -> Optional[str]:
        """The GCA_ accession, whichever slot it landed in."""
        for acc in (self.assembly_accession, self.paired_accession):
            if acc and acc.upper().startswith("GCA_"):
                return acc
        return None

    @property
    def refseq_accession(self) -> Optional[str]:
        for acc in (self.assembly_accession, self.paired_accession):
            if acc and acc.upper().startswith("GCF_"):
                return acc
        return None

    def to_dict(self) -> dict:
        d = asdict(self)
        d["genbank_accession"] = self.genbank_accession
        d["refseq_accession"] = self.refseq_accession
        return d


@dataclass
class GeneTable:
    """The gene universe for one target, plus identifier crosswalks.

    ``frame`` has one row per analysed coding sequence with the columns
    ``locus_tag, old_locus_tag, protein_id, gene, product, seqid, start,
    end, strand, order`` -- ``order`` being the position in the CDS FASTA so
    that output can preserve the input ordering.
    """

    frame: pd.DataFrame
    source: str                       # "cds_fasta+gff", "cds_fasta", "gff"
    maps: Dict[str, Dict[str, str]] = field(default_factory=dict)
    notes: List[str] = field(default_factory=list)

    @property
    def locus_tags(self) -> List[str]:
        return [t for t in self.frame["locus_tag"].tolist() if t]

    def lookup(self, value, namespaces=("locus_tag", "old_locus_tag",
                                        "protein_id", "protein_id_nover",
                                        "gene_upper")) -> Optional[str]:
        """Resolve any supported identifier to a canonical ``locus_tag``."""
        if value is None:
            return None
        raw = str(value).strip()
        if not raw:
            return None
        candidates = [raw, raw.upper()]
        # protein accessions are frequently quoted with or without version
        if "." in raw:
            candidates.append(raw.split(".")[0])
        for ns in namespaces:
            table = self.maps.get(ns) or {}
            for cand in candidates:
                hit = table.get(cand)
                if hit:
                    return hit
        return None


# ---------------------------------------------------------------------------
# taxonomy / assembly resolution
# ---------------------------------------------------------------------------
def _taxon_report(web: WebCache, token: str) -> Optional[dict]:
    from urllib.parse import quote

    try:
        payload = web.get_json(f"{DATASETS}/taxonomy/taxon/{quote(str(token))}/dataset_report")
    except Exception:
        return None
    reports = (payload or {}).get("reports") or []
    for rep in reports:
        tax = rep.get("taxonomy") or {}
        if tax.get("tax_id"):
            return tax
    return None


def _taxon_suggest(web: WebCache, text: str) -> Optional[int]:
    from urllib.parse import quote

    try:
        payload = web.get_json(f"{DATASETS}/taxonomy/taxon_suggest/{quote(str(text))}")
    except Exception:
        return None
    for hit in (payload or {}).get("sci_name_and_ids") or []:
        try:
            return int(hit["tax_id"])
        except Exception:
            continue
    return None


def _esearch_taxid(web: WebCache, text: str, email: str = "", api_key: str = "") -> Optional[int]:
    params = {"db": "taxonomy", "term": str(text), "retmode": "json", "retmax": "5"}
    if email:
        params["email"] = email
        params["tool"] = "CodonPipe"
    if api_key:
        params["api_key"] = api_key
    try:
        payload = web.get_json(f"{EUTILS}/esearch.fcgi", params)
    except Exception:
        return None
    ids = ((payload or {}).get("esearchresult") or {}).get("idlist") or []
    try:
        return int(ids[0])
    except Exception:
        return None


def _apply_lineage(info: "TargetInfo", tax: dict) -> None:
    """Record the lineage plus the species and genus identifiers.

    The species-level taxon matters because UniProtKB curation is concentrated
    in a few reference strains: for most pathogens the strain actually being
    analysed carries almost no reviewed entries, while its species does.
    """
    classification = tax.get("classification") or {}
    # NCBI returns the full ancestor chain in 'parents' and a rank-keyed
    # 'classification' map. Pool both: 'parents' is ordered root-to-strain but
    # is occasionally absent, and 'classification' carries the named ranks.
    ids: List[int] = []
    for value in tax.get("parents") or []:
        try:
            ids.append(int(value))
        except Exception:
            continue
    for node in classification.values():
        if isinstance(node, dict) and node.get("id"):
            try:
                ids.append(int(node["id"]))
            except Exception:
                continue
    if info.tax_id:
        ids.append(int(info.tax_id))
    seen: Set[int] = set()
    info.lineage_tax_ids = [i for i in ids if not (i in seen or seen.add(i))]
    info.lineage_names = [
        str(v.get("name"))
        for v in classification.values()
        if isinstance(v, dict) and v.get("name")
    ]
    for rank, slot_id, slot_name in (
        ("species", "species_tax_id", "species_name"),
        ("genus", "genus_tax_id", "genus_name"),
    ):
        node = classification.get(rank) or {}
        if isinstance(node, dict) and node.get("id"):
            try:
                setattr(info, slot_id, int(node["id"]))
            except Exception:
                pass
            setattr(info, slot_name, node.get("name"))
    if info.species_tax_id is None and str(tax.get("rank") or "").upper() == "SPECIES":
        info.species_tax_id = info.tax_id
        info.species_name = info.scientific_name


def _assembly_report(web: WebCache, accession: str) -> Optional[dict]:
    try:
        payload = web.get_json(
            f"{DATASETS}/genome/accession/{accession}/dataset_report",
            {"page_size": "1"},
        )
    except Exception:
        return None
    reports = (payload or {}).get("reports") or []
    return reports[0] if reports else None


def _genome_reports_for_taxon(web: WebCache, tax_id: int, reference_only: bool) -> List[dict]:
    params = {
        "filters.assembly_source": "refseq",
        "filters.exclude_atypical": "true",
        "page_size": "100",
    }
    if reference_only:
        params["filters.reference_only"] = "true"
    try:
        payload = web.get_json(f"{DATASETS}/genome/taxon/{tax_id}/dataset_report", params)
    except Exception:
        return []
    return (payload or {}).get("reports") or []


def _score_assembly(rep: dict) -> Tuple[int, int, str]:
    info = rep.get("assembly_info") or {}
    level = str(info.get("assembly_level") or "").lower()
    category = str(info.get("refseq_category") or "").lower()
    cat_rank = 0 if "reference" in category else (1 if "representative" in category else 2)
    return (cat_rank, _LEVEL_RANK.get(level, 9), str(rep.get("accession") or ""))


def _fill_from_report(info: TargetInfo, rep: dict) -> None:
    ainfo = rep.get("assembly_info") or {}
    org = rep.get("organism") or {}
    info.assembly_accession = rep.get("accession")
    info.paired_accession = rep.get("paired_accession")
    info.assembly_name = ainfo.get("assembly_name")
    info.assembly_level = ainfo.get("assembly_level")
    info.refseq_category = ainfo.get("refseq_category")
    strain = (ainfo.get("biosample") or {}).get("strain") or (
        (org.get("infraspecific_names") or {}).get("strain")
    )
    info.strain = strain
    name = org.get("organism_name") or info.scientific_name
    if strain and name and strain.lower() not in str(name).lower():
        info.organism_label = f"{name} {strain}"
    else:
        info.organism_label = name
    if info.tax_id is None:
        try:
            info.tax_id = int(org.get("tax_id"))
        except Exception:
            pass
    if not info.scientific_name:
        info.scientific_name = org.get("organism_name")


def resolve_target(
    query: str,
    web: WebCache,
    assembly_override: str = "",
    ncbi_email: str = "",
    ncbi_api_key: str = "",
) -> TargetInfo:
    """Resolve a free-form target into a taxon and a RefSeq assembly.

    Accepted ``query`` forms
    -----------------------
    * assembly accession, ``GCF_000210855.2`` or ``GCA_000210855.2``
    * NCBI taxonomy id, ``216597`` or ``txid:216597``
    * organism text, ``Salmonella Typhimurium SL1344``,
      ``Klebsiella pneumoniae``, ``Staphylococcus aureus NCTC 8325``

    ``assembly_override`` pins the assembly explicitly while still recording
    the taxonomy resolved from ``query``.
    """
    raw = str(query or "").strip()
    info = TargetInfo(query=raw)
    if not raw and not assembly_override:
        raise ValueError("No target specified: give a species name, NCBI taxid or assembly accession.")

    # --- 1. an explicit assembly accession short-circuits everything --------
    acc_candidates = [a for a in (assembly_override, raw) if a and ASSEMBLY_RE.match(str(a).strip())]
    if acc_candidates:
        acc = acc_candidates[0].strip()
        rep = _assembly_report(web, acc)
        if rep is None:
            raise ValueError(f"NCBI returned no assembly report for accession {acc!r}.")
        _fill_from_report(info, rep)
        info.resolution_notes.append(f"Assembly taken directly from accession {acc}.")
        if info.tax_id:
            tax = _taxon_report(web, str(info.tax_id))
            if tax:
                info.scientific_name = (tax.get("current_scientific_name") or {}).get("name") or info.scientific_name
                info.rank = tax.get("rank")
                _apply_lineage(info, tax)
        return info

    # --- 2. taxonomy ---------------------------------------------------------
    tax = None
    m = TAXID_RE.match(raw)
    if m:
        tax = _taxon_report(web, m.group(1))
        if tax:
            info.resolution_notes.append(f"Taxon resolved directly from identifier {m.group(1)}.")
    if tax is None:
        tax = _taxon_report(web, raw)
        if tax:
            info.resolution_notes.append("Taxon resolved by exact NCBI Taxonomy name match.")
    if tax is None:
        tid = _taxon_suggest(web, raw)
        if tid:
            tax = _taxon_report(web, str(tid))
            if tax:
                info.resolution_notes.append(
                    "Taxon resolved through the NCBI Taxonomy name-suggest service."
                )
    if tax is None:
        tid = _esearch_taxid(web, raw, ncbi_email, ncbi_api_key)
        if tid:
            tax = _taxon_report(web, str(tid))
            if tax:
                info.resolution_notes.append(
                    "Taxon resolved through an Entrez taxonomy search (synonym match)."
                )
    if tax is None:
        if getattr(web, "offline", False):
            raise ValueError(
                f"Offline mode is enabled and {raw!r} has not been resolved before, "
                "so there is nothing cached for it. Run this target once with "
                "network access, or point the cache folder at the one used by the "
                "original run."
            )
        raise ValueError(
            f"Could not resolve {raw!r} to an NCBI taxon. Try the exact scientific "
            "name, an NCBI taxonomy id, or a GCF_/GCA_ assembly accession."
        )

    info.tax_id = int(tax.get("tax_id"))
    info.scientific_name = (tax.get("current_scientific_name") or {}).get("name")
    info.rank = tax.get("rank")
    _apply_lineage(info, tax)

    # --- 3. assembly for that taxon -----------------------------------------
    reports = _genome_reports_for_taxon(web, info.tax_id, reference_only=True)
    if reports:
        info.resolution_notes.append("Assembly is the NCBI-designated reference genome for this taxon.")
    else:
        reports = _genome_reports_for_taxon(web, info.tax_id, reference_only=False)
        if reports:
            info.resolution_notes.append(
                "No NCBI reference genome for this taxon; the most complete "
                "RefSeq assembly was selected."
            )
    if not reports:
        raise ValueError(
            f"No RefSeq genome assembly is available for taxon {info.tax_id} "
            f"({info.scientific_name}). Supply an assembly accession instead."
        )
    reports.sort(key=_score_assembly)
    _fill_from_report(info, reports[0])
    if len(reports) > 1:
        info.resolution_notes.append(
            f"{len(reports)} RefSeq assemblies were available; "
            f"{info.assembly_accession} was chosen by reference status then assembly level."
        )
    return info


# ---------------------------------------------------------------------------
# annotation retrieval and parsing
# ---------------------------------------------------------------------------
_NUC_ACC_RE = re.compile(r"\b((?:N[CGTWZ]|A[CP]|C[PY]|[A-Z]{4})_?\d{5,8}(?:\.\d+)?)\b")


def nucleotide_accessions_from_fasta(path: str, max_records: int = 200) -> List[str]:
    """Collect the replicon accessions named in a CDS FASTA header block.

    RefSeq coding-sequence headers start ``>lcl|NC_016810.1_cds_...``, so the
    accession of the chromosome (and of each plasmid) can be read directly off
    the file the user already selected in CodonPipe.
    """
    found: List[str] = []
    seen = set()
    count = 0
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            count += 1
            if count > max_records:
                break
            head = line[1:].split()[0] if len(line) > 1 else ""
            head = head.split("|")[1] if "|" in head else head
            match = _NUC_ACC_RE.search(head)
            if match and match.group(1) not in seen:
                seen.add(match.group(1))
                found.append(match.group(1))
    return found


def assembly_from_cds_fasta(
    path: str, web: WebCache, ncbi_email: str = "", ncbi_api_key: str = ""
) -> Tuple[Optional[str], List[str]]:
    """Identify the assembly a coding-sequence FASTA was derived from.

    Returns ``(assembly_accession, notes)``.  The replicon accession is linked
    to its assembly through Entrez, which lets the desktop interface offer a
    one-click 'use the genome I already selected' action without the user
    having to look anything up.
    """
    notes: List[str] = []
    accessions = nucleotide_accessions_from_fasta(path)
    if not accessions:
        notes.append(
            "No INSDC/RefSeq replicon accession could be read from the FASTA headers; "
            "the genome has to be named explicitly."
        )
        return None, notes
    base = {"db": "assembly", "retmode": "json"}
    if ncbi_email:
        base["email"] = ncbi_email
        base["tool"] = "CodonPipe"
    if ncbi_api_key:
        base["api_key"] = ncbi_api_key
    for accession in accessions[:3]:
        try:
            link = web.get_json(
                f"{EUTILS}/elink.fcgi",
                {**base, "dbfrom": "nuccore", "id": accession},
            )
        except Exception as exc:
            notes.append(f"Entrez link lookup failed for {accession} ({exc}).")
            continue
        uids: List[str] = []
        for linkset in (link or {}).get("linksets", []):
            for db in linkset.get("linksetdbs", []) or []:
                uids.extend(str(u) for u in (db.get("links") or []))
        if not uids:
            continue
        try:
            payload = web.get_json(f"{EUTILS}/esummary.fcgi", {**base, "id": uids[0]})
        except Exception as exc:
            notes.append(f"Entrez summary failed for assembly uid {uids[0]} ({exc}).")
            continue
        record = ((payload or {}).get("result") or {}).get(str(uids[0])) or {}
        synonym = record.get("synonym") or {}
        accession_out = synonym.get("refseq") or record.get("assemblyaccession")
        if accession_out:
            notes.append(
                f"Assembly {accession_out} identified from replicon {accession} "
                "in the selected FASTA."
            )
            return str(accession_out), notes
    notes.append(
        "The replicon accessions in the FASTA could not be linked to an assembly; "
        "the genome has to be named explicitly."
    )
    return None, notes


def fetch_refseq_gff(accession: str, web: WebCache) -> str:
    """Download and cache the RefSeq GFF3 annotation for ``accession``."""
    data = web.get_bytes(
        f"{DATASETS}/genome/accession/{accession}/download",
        {"include_annotation_type": "GENOME_GFF", "filename": f"{accession}_gff.zip"},
        ext=".zip",
    )
    try:
        archive = zipfile.ZipFile(io.BytesIO(data))
    except Exception as exc:
        raise NetworkUnavailable(
            f"The annotation download for {accession} was not a readable zip archive: {exc}"
        )
    names = [n for n in archive.namelist() if n.endswith(".gff") or n.endswith(".gff3")]
    if not names:
        raise NetworkUnavailable(
            f"No GFF annotation inside the NCBI download for {accession}. "
            "The assembly may be unannotated."
        )
    return archive.read(names[0]).decode("utf-8", errors="replace")


def fetch_cds_fasta(accession: str, web: WebCache, dest_path: str) -> str:
    """Download the RefSeq coding-sequence FASTA for ``accession``.

    Used when the button is driven purely from a species name, so that the
    gene universe and the downstream CodonPipe analysis can share exactly the
    same coding sequences.  The file is written to ``dest_path`` in the
    standard NCBI ``cds_from_genomic`` form.
    """
    data = web.get_bytes(
        f"{DATASETS}/genome/accession/{accession}/download",
        {
            "include_annotation_type": "CDS_FASTA",
            "filename": f"{accession}_cds.zip",
        },
        ext=".zip",
    )
    try:
        archive = zipfile.ZipFile(io.BytesIO(data))
    except Exception as exc:
        raise NetworkUnavailable(
            f"The coding-sequence download for {accession} was not a readable zip: {exc}"
        )
    names = [n for n in archive.namelist() if n.endswith(".fna") or n.endswith(".fasta")]
    if not names:
        raise NetworkUnavailable(
            f"No coding-sequence FASTA inside the NCBI download for {accession}."
        )
    os.makedirs(os.path.dirname(os.path.abspath(dest_path)) or ".", exist_ok=True)
    payload = archive.read(names[0])
    tmp = dest_path + ".part"
    with open(tmp, "wb") as fh:
        fh.write(payload)
    os.replace(tmp, dest_path)
    return dest_path


def _gff_attrs(field_text: str) -> Dict[str, str]:
    from urllib.parse import unquote

    out: Dict[str, str] = {}
    for chunk in str(field_text or "").split(";"):
        if "=" not in chunk:
            continue
        key, _, value = chunk.partition("=")
        out[key.strip()] = unquote(value.strip())
    return out


def parse_gff_genes(gff_text: str) -> pd.DataFrame:
    """Extract protein-coding gene records from a RefSeq GFF3.

    Returns one row per ``locus_tag`` combining the ``gene`` feature (which
    carries ``old_locus_tag``) with its ``CDS`` child (which carries
    ``protein_id`` and the functional ``product``).
    """
    genes: Dict[str, dict] = {}
    cds: Dict[str, dict] = {}
    for line in gff_text.splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        ftype = parts[2]
        if ftype not in ("gene", "CDS"):
            continue
        attrs = _gff_attrs(parts[8])
        locus = attrs.get("locus_tag")
        if not locus:
            continue
        rec = {
            "locus_tag": locus,
            "seqid": parts[0],
            "start": int(parts[3]) if parts[3].isdigit() else None,
            "end": int(parts[4]) if parts[4].isdigit() else None,
            "strand": parts[6],
            "gene": attrs.get("gene") or attrs.get("Name"),
            "product": attrs.get("product"),
            "old_locus_tag": attrs.get("old_locus_tag"),
            "protein_id": attrs.get("protein_id"),
            "gene_biotype": attrs.get("gene_biotype"),
            "dbxref": attrs.get("Dbxref"),
        }
        target = genes if ftype == "gene" else cds
        if locus not in target:
            target[locus] = rec

    rows = []
    for locus, g in genes.items():
        c = cds.get(locus) or {}
        if g.get("gene_biotype") not in (None, "protein_coding") and not c:
            continue
        if not c and g.get("gene_biotype") != "protein_coding":
            continue
        rows.append(
            {
                "locus_tag": locus,
                "old_locus_tag": g.get("old_locus_tag") or c.get("old_locus_tag"),
                "protein_id": c.get("protein_id"),
                "gene": g.get("gene") or c.get("gene"),
                "product": c.get("product") or g.get("product"),
                "seqid": g.get("seqid") or c.get("seqid"),
                "start": g.get("start") if g.get("start") is not None else c.get("start"),
                "end": g.get("end") if g.get("end") is not None else c.get("end"),
                "strand": g.get("strand") or c.get("strand"),
                "dbxref": g.get("dbxref") or c.get("dbxref"),
            }
        )
    # CDS without a parent gene feature (rare, but do not silently drop genes)
    for locus, c in cds.items():
        if locus in genes:
            continue
        rows.append(
            {
                "locus_tag": locus,
                "old_locus_tag": c.get("old_locus_tag"),
                "protein_id": c.get("protein_id"),
                "gene": c.get("gene"),
                "product": c.get("product"),
                "seqid": c.get("seqid"),
                "start": c.get("start"),
                "end": c.get("end"),
                "strand": c.get("strand"),
                "dbxref": c.get("dbxref"),
            }
        )
    frame = pd.DataFrame(rows)
    if not frame.empty:
        frame = frame.sort_values(["seqid", "start"], na_position="last").reset_index(drop=True)
    return frame


_HEADER_FIELD_RE = re.compile(r"\[([A-Za-z_]+)=([^\]]*)\]")


def parse_cds_fasta(path: str) -> pd.DataFrame:
    """Read identifiers out of a RefSeq-style coding-sequence FASTA.

    This defines the gene universe actually analysed by CodonPipe, so the
    discovery output is guaranteed to use the same identifiers in the same
    order.
    """
    rows = []
    order = 0
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            fields = {k: v for k, v in _HEADER_FIELD_RE.findall(header)}
            record_id = header.split()[0] if header else ""
            protein = fields.get("protein_id")
            if not protein:
                # ``lcl|NC_016810.1_cds_WP_001575544.1_1`` -> WP_001575544.1
                m = re.search(r"_cds_([A-Za-z]{2,3}_?\d+\.\d+)_\d+$", record_id)
                if m:
                    protein = m.group(1)
            rows.append(
                {
                    "order": order,
                    "fasta_id": record_id,
                    "locus_tag": fields.get("locus_tag"),
                    "old_locus_tag": fields.get("old_locus_tag"),
                    "protein_id": protein,
                    "gene": fields.get("gene"),
                    "product": fields.get("protein"),
                    "location": fields.get("location"),
                    "pseudo": fields.get("pseudo"),
                }
            )
            order += 1
    return pd.DataFrame(rows)


def _add_map(
    maps: Dict[str, Dict[str, str]],
    namespace: str,
    key,
    locus: str,
    collisions: Optional[Dict[str, Set[str]]] = None,
) -> None:
    """Record one identifier -> locus_tag edge, remembering genuine collisions.

    First writer wins, which is fine for identifiers that are unique by
    construction (locus tags, protein accessions). Gene SYMBOLS are not unique:
    a genome can carry several copies of the same symbol, and silently keeping
    whichever came first in FASTA order makes the lookup order-dependent. Those
    keys are recorded here and dropped by :func:`_drop_ambiguous_symbols`, so an
    ambiguous symbol resolves to nothing rather than to an arbitrary one of its
    copies.
    """
    if key is None or locus is None:
        return
    key = str(key).strip()
    if not key or key.lower() in ("nan", "none"):
        return
    table = maps.setdefault(namespace, {})
    existing = table.get(key)
    if existing is None:
        table[key] = locus
    elif existing != locus and collisions is not None:
        collisions.setdefault(namespace, set()).add(key)


def _drop_ambiguous_symbols(
    maps: Dict[str, Dict[str, str]], collisions: Dict[str, Set[str]]
) -> List[str]:
    """Remove gene-symbol keys that matched more than one locus."""
    notes: List[str] = []
    ambiguous = collisions.get("gene_upper") or set()
    if not ambiguous:
        return notes
    table = maps.get("gene_upper") or {}
    for key in ambiguous:
        table.pop(key, None)
    notes.append(
        f"{len(ambiguous)} gene symbol(s) occur on more than one locus in this "
        "genome and were removed from the symbol lookup, so they resolve to "
        "nothing rather than to an arbitrary copy. Locus tags and protein "
        "accessions are unaffected."
    )
    return notes


def build_crosswalk(
    cds_frame: Optional[pd.DataFrame],
    gff_frame: Optional[pd.DataFrame],
) -> GeneTable:
    """Merge the FASTA gene universe with the RefSeq annotation.

    The FASTA is authoritative for *which* genes exist and in what order; the
    GFF supplies coordinates, old locus tags and curated product names.  When
    only one of the two is available the table is built from that alone.
    """
    notes: List[str] = []
    have_cds = cds_frame is not None and not cds_frame.empty
    have_gff = gff_frame is not None and not gff_frame.empty

    if have_cds and have_gff:
        source = "cds_fasta+gff"
        gff = gff_frame.copy()
        # Join on locus_tag first, then fall back to protein accession.
        base = cds_frame.copy()
        merged = base.merge(
            gff.add_suffix("_gff"),
            left_on="locus_tag",
            right_on="locus_tag_gff",
            how="left",
        )
        missing = merged["locus_tag_gff"].isna()
        if missing.any():
            prot_index = gff.dropna(subset=["protein_id"]).drop_duplicates("protein_id")
            prot_index = prot_index.set_index("protein_id")
            for idx in merged.index[missing]:
                prot = merged.at[idx, "protein_id"]
                if prot in prot_index.index:
                    row = prot_index.loc[prot]
                    for col in gff.columns:
                        merged.at[idx, f"{col}_gff"] = row[col]
            still = merged["locus_tag_gff"].isna().sum()
            if still:
                notes.append(
                    f"{int(still)} of {len(merged)} FASTA records could not be matched "
                    "to the RefSeq annotation; coordinate-based context evidence is "
                    "unavailable for those genes."
                )
        frame = pd.DataFrame(
            {
                "order": merged["order"],
                "locus_tag": merged["locus_tag"].where(
                    merged["locus_tag"].notna(), merged.get("locus_tag_gff")
                ),
                "old_locus_tag": merged["old_locus_tag"].where(
                    merged["old_locus_tag"].notna(), merged.get("old_locus_tag_gff")
                ),
                "protein_id": merged["protein_id"].where(
                    merged["protein_id"].notna(), merged.get("protein_id_gff")
                ),
                "gene": merged["gene"].where(merged["gene"].notna(), merged.get("gene_gff")),
                "product": merged["product"].where(
                    merged["product"].notna(), merged.get("product_gff")
                ),
                "seqid": merged.get("seqid_gff"),
                "start": merged.get("start_gff"),
                "end": merged.get("end_gff"),
                "strand": merged.get("strand_gff"),
            }
        )
    elif have_cds:
        source = "cds_fasta"
        frame = cds_frame.reindex(
            columns=["order", "locus_tag", "old_locus_tag", "protein_id", "gene", "product"]
        ).copy()
        for col in ("seqid", "start", "end", "strand"):
            frame[col] = None
        notes.append(
            "No RefSeq annotation was available, so genomic-context evidence "
            "(operon and island extension) is disabled for this run."
        )
    elif have_gff:
        source = "gff"
        frame = gff_frame.reindex(
            columns=["locus_tag", "old_locus_tag", "protein_id", "gene", "product",
                     "seqid", "start", "end", "strand"]
        ).copy()
        frame.insert(0, "order", range(len(frame)))
        notes.append("Gene universe taken from the RefSeq annotation (no CDS FASTA supplied).")
    else:
        raise ValueError("A CDS FASTA or a RefSeq annotation is required to build the gene table.")

    frame = frame.dropna(subset=["locus_tag"]).copy()
    frame["locus_tag"] = frame["locus_tag"].astype(str)
    frame = frame.drop_duplicates("locus_tag", keep="first").reset_index(drop=True)
    for col in ("start", "end"):
        frame[col] = pd.to_numeric(frame[col], errors="coerce")

    maps: Dict[str, Dict[str, str]] = {}
    collisions: Dict[str, Set[str]] = {}
    for row in frame.itertuples(index=False):
        locus = row.locus_tag
        _add_map(maps, "locus_tag", locus, locus, collisions)
        _add_map(maps, "locus_tag", str(locus).upper(), locus, collisions)
        _add_map(maps, "old_locus_tag", row.old_locus_tag, locus, collisions)
        if row.old_locus_tag:
            _add_map(maps, "old_locus_tag", str(row.old_locus_tag).upper(), locus, collisions)
        _add_map(maps, "protein_id", row.protein_id, locus, collisions)
        if row.protein_id:
            _add_map(maps, "protein_id_nover", str(row.protein_id).split(".")[0], locus, collisions)
        if row.gene:
            _add_map(maps, "gene_upper", str(row.gene).upper(), locus, collisions)
    notes.extend(_drop_ambiguous_symbols(maps, collisions))

    return GeneTable(frame=frame, source=source, maps=maps, notes=notes)
