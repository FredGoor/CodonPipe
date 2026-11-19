#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Output writers for automated virulence-gene discovery.

Four artefacts are produced per run:

``<organism> virulence clusters.xlsx``
    One column per gene cluster, each column a list of locus tags.  This is
    the format CodonPipe already reads from the *Gene clusters* selector, so
    the result can be used directly as the cluster file without any
    reformatting.

``<organism> virulence evidence.xlsx``
    The audit workbook: every gene of the genome with its evidence score, the
    individual score components, the per-postulate flags, the database records
    behind the call and the reason it was included or excluded, plus the raw
    evidence tables on separate sheets.

``<organism> virulence discovery report.txt``
    A human-readable provenance record: target resolution, parameters, source
    databases with access dates, per-layer counts, warnings and citations.
    Written to be pasted into a methods section or a reviewer response.

``<organism> virulence run.json``
    The machine-readable manifest: the complete parameter set, resolved
    accessions, resource-file checksums and the network access log, so the run
    can be repeated exactly.
"""

from __future__ import annotations

import hashlib
import json
import os
import time
from typing import Dict, List, Optional, Sequence

import pandas as pd

from .classify import (
    CATEGORY_GENERAL,
    CATEGORY_HOST,
    IN_VIVO_SPECIFIC_LABEL,
    LADDER,
    LADDER_LABELS,
    TIERS,
    TIER_LABELS,
    TIER_DEFINITIONS,
)

__all__ = [
    "CITATIONS",
    "safe_name",
    "build_cluster_columns",
    "write_cluster_workbook",
    "write_evidence_workbook",
    "write_report",
    "write_manifest",
]

CITATIONS: List[str] = [
    "Falkow S (1988) Molecular Koch's postulates applied to microbial pathogenicity. "
    "Rev Infect Dis 10(Suppl 2):S274-S276.",
    "Falkow S (2004) Molecular Koch's postulates applied to bacterial pathogenicity - "
    "a personal recollection 15 years later. Nat Rev Microbiol 2:67-72.",
    "Chen L, Yang J, Yu J, Yao Z, Sun L, Shen Y, Jin Q (2005) VFDB: a reference database "
    "for bacterial virulence factors. Nucleic Acids Res 33:D325-D328.",
    "Sayers S, Li L, Ong E, et al. (2019) Victors: a web-based knowledge base of virulence "
    "factors in human and animal pathogens. Nucleic Acids Res 47:D693-D700.",
    "Olson RD, Assaf R, Brettin T, et al. (2023) Introducing the Bacterial and Viral "
    "Bioinformatics Resource Center (BV-BRC): a resource combining PATRIC, IRD and ViPR. "
    "Nucleic Acids Res 51:D678-D689.",
    "The UniProt Consortium (2025) UniProt: the Universal Protein Knowledgebase in 2025. "
    "Nucleic Acids Res 53:D609-D617.",
    "Binns D, Dimmer E, Huntley R, Barrell D, O'Donovan C, Apweiler R (2009) QuickGO: a "
    "web-based tool for Gene Ontology searching. Bioinformatics 25:3045-3046.",
    "Lam MMC, Wick RR, Watts SC, Cerdeira LT, Wyres KL, Holt KE (2021) A genomic "
    "surveillance framework and genotyping tool for Klebsiella pneumoniae and its related "
    "species complex. Nat Commun 12:4188.",
    "Lam MMC, Wick RR, Wyres KL, Gorrie CL, Judd LM, Jenney AWJ, Brisse S, Holt KE (2018) "
    "Genetic diversity, mobilisation and spread of the yersiniabactin-encoding mobile "
    "element ICEKp in Klebsiella pneumoniae populations. Microb Genom 4:e000196.",
    "Lam MMC, Wyres KL, Judd LM, Wick RR, Jenney A, Brisse S, Holt KE (2018) Tracking key "
    "virulence loci encoding aerobactin and salmochelin siderophore synthesis in "
    "Klebsiella pneumoniae. Genome Med 10:77.",
    "Walker KA, Miner TA, Palacios M, Trzilova D, Frederick DR, Broberg CA, Sepulveda VE, "
    "Quinn JD, Miller VL (2019) A Klebsiella pneumoniae regulatory mutant has reduced "
    "capsule expression but retains hypermucoviscosity. mBio 10:e00089-19.",
    "Sayers EW, Beck J, Bolton EE, et al. (2024) Database resources of the National Center "
    "for Biotechnology Information. Nucleic Acids Res 52:D33-D43.",
    "O'Leary NA, Wright MW, Brister JR, et al. (2016) Reference sequence (RefSeq) database "
    "at NCBI: current status, taxonomic expansion, and functional annotation. "
    "Nucleic Acids Res 44:D733-D745.",
]


def safe_name(text: str, fallback: str = "genome", max_len: int = 90) -> str:
    """Make a string safe for use in a Windows file name."""
    cleaned = "".join(
        c if (c.isalnum() or c in " -_.()") else "_" for c in str(text or "")
    ).strip()
    cleaned = " ".join(cleaned.split())
    return (cleaned[:max_len] or fallback)


# ---------------------------------------------------------------------------
# cluster workbook
# ---------------------------------------------------------------------------
def build_cluster_columns(
    genes: pd.DataFrame,
    tiers: Sequence[str] = ("permissive", "standard", "strict"),
    include_subcategories: bool = True,
    min_subcategory_size: int = 2,
    subcategory_tier: str = "standard",
) -> "dict[str, List[str]]":
    """Assemble the ``{cluster name: [locus tags]}`` mapping.

    One column per requested tier, in ladder order, followed by the functional
    subcategories.  Genes keep the order of the input coding-sequence FASTA, so
    the workbook stays aligned with the rest of a CodonPipe project.
    """
    ordered = genes.sort_values("order", na_position="last")
    columns: "dict[str, List[str]]" = {}

    # The five-rung ladder, in order. A rung the caller did not ask for is
    # skipped; a rung that is empty for this genome (rungs 4 and 5 without
    # mutant-screen data) is dropped at the end with the other empties.
    legacy = {"permissive": "permissive", "standard": "standard", "strict": "host_interaction"}
    wanted = set()
    for name in tiers:
        key = str(name).strip().lower()
        wanted.add(legacy.get(key, key))
    if "all" in wanted or not wanted:
        wanted = set(LADDER)
    # Asking for the legacy "strict" means the rung the run actually used.
    if "strict" in {str(t).strip().lower() for t in tiers} and "tier_host_interaction" in ordered.columns:
        wanted.add("host_interaction")
    for rung in LADDER:
        if rung not in wanted:
            continue
        flag = f"tier_{rung}"
        if flag not in ordered.columns:
            continue
        columns[LADDER_LABELS[rung]] = ordered.loc[
            ordered[flag].fillna(False).astype(bool), "locus_tag"
        ].tolist()

    # The experimental set, where measurements exist for this organism. It is
    # written after the ladder and before the subcategories, and only when it
    # is non-empty, so a genome without screen data is unaffected.
    if "in_vivo_specific" in ordered.columns:
        members = ordered.loc[
            ordered["in_vivo_specific"].fillna(False).astype(bool), "locus_tag"
        ].tolist()
        if members:
            columns[IN_VIVO_SPECIFIC_LABEL] = members

    if include_subcategories:
        flag = f"tier_{subcategory_tier}"
        if flag not in ordered.columns:
            flag = "tier_standard"
        pool = ordered[ordered[flag].fillna(False).astype(bool)]
        counts: "dict[str, List[str]]" = {}
        for row in pool.itertuples(index=False):
            source = getattr(row, "cluster_subcategories", None)
            if source is None:
                source = getattr(row, "subcategories", "")
            for name in str(source or "").split(";"):
                name = name.strip()
                if name:
                    counts.setdefault(name, []).append(row.locus_tag)
        for name in sorted(counts):
            members = counts[name]
            if len(members) >= int(min_subcategory_size):
                columns[name] = members
    return {k: v for k, v in columns.items() if v}


def _write_columns_sheet(writer, sheet_name: str, columns: "dict[str, List[str]]") -> None:
    if not columns:
        pd.DataFrame({"(no genes)": []}).to_excel(writer, sheet_name=sheet_name, index=False)
        return
    height = max(len(v) for v in columns.values())
    padded = {k: list(v) + [None] * (height - len(v)) for k, v in columns.items()}
    pd.DataFrame(padded).to_excel(writer, sheet_name=sheet_name, index=False)


def _cluster_guide(columns: "dict[str, List[str]]") -> pd.DataFrame:
    """One row per column of the cluster sheet: name, size, and what it means.

    The sheet is built to be pasted alongside a user's own gene clusters, where
    the column header is all the context a reader gets. Without this, a column
    headed "6 In vivo specific" in someone else's workbook is a list of locus
    tags with no recoverable definition, and the rung's two clauses -- which are
    the whole point of it -- are invisible.
    """
    from .classify import LADDER, LADDER_LABELS, LADDER_DEFINITIONS

    label_to_rung = {LADDER_LABELS[r]: r for r in LADDER if r in LADDER_LABELS}
    rows = []
    for name, members in columns.items():
        rung = label_to_rung.get(name, "")
        if rung:
            kind = "ladder rung"
            meaning = LADDER_DEFINITIONS.get(rung, "")
        elif name == IN_VIVO_SPECIFIC_LABEL:
            kind = "not a rung"
            meaning = (
                "Every gene attenuated in an animal infection screen that is not "
                "required for growth in vitro. Defined by measurement alone, with "
                "no annotation input, and it stands OUTSIDE the ladder. It is NOT "
                "a host-interaction set: 'required in vivo and dispensable in "
                "vitro' is weaker than 'acts on the host' -- an auxotroph "
                "qualifies, because the host is nutrient-poor and rich medium is "
                "not."
            )
        else:
            kind = "functional category"
            meaning = (
                "Genes of the curation-anchored rung carrying this category in the "
                "vocabulary table. Categories are descriptive labels, not rungs: "
                "they overlap, and membership implies nothing about evidence "
                "strength."
            )
        rows.append({"column": name, "kind": kind, "genes": len(members),
                     "definition": meaning})
    return pd.DataFrame(rows, columns=["column", "kind", "genes", "definition"])


def write_cluster_workbook(
    path: str, columns: "dict[str, List[str]]", sheet_name: str = "Virulence clusters"
) -> str:
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with pd.ExcelWriter(path, engine="xlsxwriter") as writer:
        _write_columns_sheet(writer, sheet_name[:31], columns)
        guide = _cluster_guide(columns)
        guide.to_excel(writer, sheet_name="Column guide", index=False)
        sheet = writer.sheets["Column guide"]
        wrap = writer.book.add_format({"text_wrap": True, "valign": "top"})
        sheet.set_column(0, 0, 44)
        sheet.set_column(1, 1, 18)
        sheet.set_column(2, 2, 8)
        sheet.set_column(3, 3, 96, wrap)
        sheet.freeze_panes(1, 0)
    return path


# ---------------------------------------------------------------------------
# evidence workbook
# ---------------------------------------------------------------------------
_AUDIT_COLUMNS = [
    "locus_tag", "old_locus_tag", "gene", "product", "protein_id",
    "seqid", "start", "end", "strand",
    "tier", "tier_permissive", "tier_standard",
    "tier_host_interaction", "tier_in_vivo_contributing", "tier_in_vivo_required",
    "tier_in_vitro_neutral", "in_vivo_outside_rung2", "in_vivo_discordant",
    "in_vivo_effect_log2", "in_vivo_effect_rank", "in_vivo_significance",
    "in_vivo_specific",
    "virulence_general", "host_interaction", "evidence_score",
    "host_interaction_axis", "core_process_axis",
    "koch_postulate_1", "koch_postulate_2", "koch_postulate_2_basis",
    "koch_postulate_2_satisfied", "koch_postulate_2_satisfied_basis",
    "koch_postulate_3", "koch_postulate_3_basis", "in_vivo_fitness", "in_vivo_model",
    "in_vitro_requirement_source",
    "machinery", "species_locus", "species_locus_role",
    "vf_axes", "subcategories", "cluster_subcategories",
    "context_layer", "context_detail",
    "database_evidence", "host_evidence", "core_evidence",
    "score_components", "classification_reason",
]


def write_evidence_workbook(
    path: str,
    genes: pd.DataFrame,
    evidence_tables: Optional[Dict[str, pd.DataFrame]] = None,
    summary_rows: Optional[pd.DataFrame] = None,
) -> str:
    """Write the per-gene audit table plus the raw evidence layers."""
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    audit = genes.reindex(columns=[c for c in _AUDIT_COLUMNS if c in genes.columns])
    flag = "tier_permissive" if "tier_permissive" in audit.columns else "virulence_general"
    called = audit[audit[flag].fillna(False).astype(bool)]
    with pd.ExcelWriter(path, engine="xlsxwriter") as writer:
        if summary_rows is not None and not summary_rows.empty:
            summary_rows.to_excel(writer, sheet_name="Run summary", index=False)
        called.sort_values("evidence_score", ascending=False).to_excel(
            writer, sheet_name="Virulence calls", index=False
        )
        audit.to_excel(writer, sheet_name="All genes", index=False)
        for name, frame in (evidence_tables or {}).items():
            if frame is None or getattr(frame, "empty", True):
                continue
            frame.to_excel(writer, sheet_name=str(name)[:31], index=False)
    return path


# ---------------------------------------------------------------------------
# provenance report
# ---------------------------------------------------------------------------
def _fmt_counts(columns: "dict[str, List[str]]") -> List[str]:
    return [f"    {name:<44s} {len(members):>6d} genes" for name, members in columns.items()]


def write_report(
    path: str,
    target,
    params,
    summary: dict,
    columns: "dict[str, List[str]]",
    gene_table_source: str,
    layer_counts: Dict[str, int],
    warnings: Sequence[str],
    cache_summary: Dict[str, int],
    access_log: Sequence[dict],
    extra_sections: Optional[Dict[str, Sequence[str]]] = None,
) -> str:
    """Write the human-readable provenance report."""
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    now = time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())
    hosts: Dict[str, str] = {}
    for entry in access_log:
        from urllib.parse import urlparse

        host = urlparse(str(entry.get("url"))).netloc
        if host and host not in hosts:
            hosts[host] = str(entry.get("accessed_utc") or "")

    lines: List[str] = []
    add = lines.append
    add("=" * 78)
    add("CodonPipe - automated virulence-gene discovery")
    add("=" * 78)
    add(f"Report generated      : {now}")
    add("")
    add("1. TARGET")
    add("-" * 78)
    add(f"  Query               : {target.query}")
    add(f"  NCBI taxon          : {target.tax_id}  {target.scientific_name} (rank {target.rank})")
    if target.species_tax_id:
        add(f"  Species taxon       : {target.species_tax_id}  {target.species_name}")
    add(f"  Assembly            : {target.assembly_accession}  ({target.assembly_name}, {target.assembly_level})")
    add(f"  Paired accession    : {target.paired_accession}")
    if target.refseq_category:
        add(f"  RefSeq category     : {target.refseq_category}")
    add(f"  Gene universe       : {summary.get('n_genes')} coding sequences, from {gene_table_source}")
    for note in target.resolution_notes:
        add(f"  - {note}")
    add("")
    add("2. TIER DEFINITIONS (molecular Koch's postulates)")
    add("-" * 78)
    add("  The three tiers form a strictly nested ladder, computed in one scoring")
    add("  pass. Each adds a predicate to the previous one, so the nesting holds by")
    add("  construction and is asserted at run time.")
    add("")
    add("  permissive  -- postulate 1, read inclusively")
    add("    Any virulence evidence at all: a curated virulence-factor record, a")
    add("    UniProtKB virulence annotation, operon context around one of those, or")
    add("    a supporting annotation match on its own. Deliberately over-inclusive,")
    add("    so that a result can be tested against the most generous definition.")
    add("")
    add("  standard  -- postulate 1, anchored in curation")
    add("    Evidence anchored in database curation, or in operon context around a")
    add("    curated gene, above the evidence threshold. Annotation-only candidates")
    add("    are excluded.")
    add("")
    add("  strict  -- postulate 2 interpretable (the host-interaction set)")
    add("    The standard set restricted to genes that carry positive evidence of")
    add("    acting on the host and show no evidence of a core cellular function, so")
    add("    that a loss of virulence upon specific inactivation cannot be explained")
    add("    by a generic growth defect.")
    add("")
    add("  Postulate 3 asks for a measured loss of virulence when the gene is")
    add("  inactivated. No annotation can supply that, but a published mutant")
    add("  screen can, and this run reports one of four states per gene:")
    add("    demonstrated_in_vivo          the mutant is under-represented in an")
    add("                                  animal infection screen")
    add("    refuted_in_vivo               the mutant was assayed in such a screen")
    add("                                  and showed no fitness defect")
    add("    reported_in_source_literature asserted by a literature-curated record,")
    add("                                  without a measurement in this genome")
    add("    not_assessed                  no screen covers this gene")
    add("")
    add("  'refuted_in_vivo' and 'not_assessed' are kept apart deliberately.")
    add("  Screens cover part of a genome and effector families are redundant, so")
    add("  a gene with no attenuation phenotype has not thereby been shown to be")
    add("  irrelevant. Rung 3 is the annotation definition; experimental evidence")
    add("  enters at rungs 4 to 6, so rung 3 stays comparable between genomes that")
    add("  differ in how well studied they are.")

    add("")
    exp = str(summary.get("experimental_reference", "none") or "none")
    if exp != "none":
        add(f"  Experimental gene-requirement table applied: {exp}")
        add(f"    coverage of the gene universe  {summary.get('experimental_table_coverage', 0):.1%}")
        add(f"    postulate 2 demonstrated       {summary.get('n_postulate2_demonstrated', 0)}")
        add(f"    postulate 2 refuted            {summary.get('n_postulate2_refuted', 0)}")
        add(f"    gate: in vivo contributing     {summary.get('n_tier_in_vivo_contributing', 0)}")
        add(f"    tier 3 in vivo required        {summary.get('n_tier_in_vivo_required', 0)}")
        add(f"    attenuated but outside the curation gate "
            f"{summary.get('n_in_vivo_outside_rung2', 0)}")
        add("")
        add(f"  Also written: '{IN_VIVO_SPECIFIC_LABEL}', "
            f"{summary.get('n_in_vivo_specific', 0)} genes.")
        add("    Every gene attenuated in an animal infection screen that is not")
        add("    required for growth in vitro. Defined by measurement alone, with no")
        add("    annotation input, and it stands OUTSIDE the tier ladder.")
        add("    It is not a host-interaction set: 'required in vivo and dispensable")
        add("    in vitro' is weaker than 'acts on the host'. An auxotroph is")
        add("    attenuated because the host is nutrient-poor, and in rich medium it")
        add("    carries no in vitro fitness cost, so it is retained here. The")
        add("    host-interaction tier is a subset of this set.")
        add(f"    excluded as required in vitro  "
            f"{summary.get('n_in_vivo_specific_excluded_as_essential', 0)}")
        add("  A measured value overrides the annotation-based inference in both")
        add("  directions: an essential gene cannot give an interpretable")
        add("  loss-of-virulence phenotype however host-directed it looks, and a")
        add("  gene whose viable mutant has been built and assayed is dispensable")
        add("  in vitro however much core-process annotation it carries.")
        add("")
    add("3. PARAMETERS")
    add("-" * 78)
    add(f"  Stringency          : {params.stringency}")
    for key, value in sorted(params.to_dict().items()):
        if key == "stringency":
            continue
        add(f"    {key:<32s} {value}")
    add("")
    add("4. EVIDENCE LAYERS")
    add("-" * 78)
    for name, count in layer_counts.items():
        add(f"    {name:<52s} {count:>6d}")
    add("")
    add("5. RESULT")
    add("-" * 78)
    add(f"    {'Genes in the genome':<44s} {summary.get('n_genes', 0):>6d} genes")
    lines.extend(_fmt_counts(columns))
    add("")
    add(f"    Tier ladder nested (asserted)               : {summary.get('tier_ladder_nested', False)}")
    add(f"    Annotation-only candidates (permissive only) : {summary.get('n_annotation_only_candidates', 0)}")
    add(f"    Genes with any database evidence            : {summary.get('n_with_database_evidence', 0)}")
    add(f"    Literature-curated records                  : {summary.get('n_literature_curated', 0)}")
    add(f"    UniProtKB Virulence keyword                 : {summary.get('n_uniprot_virulence_keyword', 0)}")
    add(f"    Added by genomic context                    : {summary.get('n_context_added', 0)}")
    add(f"    Virulence-associated but not host-directed  : {summary.get('n_general_only', 0)}")
    add("")
    for note in summary.get("notes", []) or []:
        add(f"  - {note}")
    add("")
    if extra_sections:
        number = 6
        for title, body in extra_sections.items():
            add(f"{number}. {title.upper()}")
            add("-" * 78)
            for row in body:
                add(f"  {row}")
            add("")
            number += 1
    add("DATA SOURCES ACCESSED")
    add("-" * 78)
    for host, when in hosts.items():
        add(f"    {host:<44s} first access {when}")
    add(f"    HTTP requests: {cache_summary.get('requests', 0)} "
        f"({cache_summary.get('from_network', 0)} from the network, "
        f"{cache_summary.get('from_cache', 0)} replayed from cache)")
    add("")
    if warnings:
        add("WARNINGS")
        add("-" * 78)
        for item in warnings:
            add(f"    - {item}")
        add("")
    add("CITATIONS")
    add("-" * 78)
    for ref in CITATIONS:
        add(f"    {ref}")
    add("")

    text = "\n".join(lines)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    return path


# ---------------------------------------------------------------------------
# manifest
# ---------------------------------------------------------------------------
def _file_digest(path: str) -> Optional[str]:
    try:
        digest = hashlib.sha256()
        with open(path, "rb") as fh:
            for block in iter(lambda: fh.read(1 << 20), b""):
                digest.update(block)
        return digest.hexdigest()[:32]
    except Exception:
        return None


def write_manifest(
    path: str,
    target,
    params,
    summary: dict,
    columns: "dict[str, List[str]]",
    resource_dir: str,
    access_log: Sequence[dict],
    outputs: Dict[str, str],
    warnings: Sequence[str],
    module_version: str,
    extra: Optional[Dict[str, object]] = None,
) -> str:
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    resources = {}
    if os.path.isdir(resource_dir):
        for name in sorted(os.listdir(resource_dir)):
            full = os.path.join(resource_dir, name)
            if os.path.isfile(full):
                resources[name] = _file_digest(full)
    payload = {
        "module": "codonpipe.virulence",
        "module_version": module_version,
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "target": target.to_dict(),
        "parameters": params.to_dict(),
        "summary": {k: v for k, v in summary.items() if k != "notes"},
        "notes": summary.get("notes", []),
        "warnings": list(warnings),
        "cluster_sizes": {k: len(v) for k, v in columns.items()},
        "cluster_members": columns,
        "resource_sha256_prefix": resources,
        "outputs": outputs,
        "access_log": list(access_log),
    }
    for key, value in (extra or {}).items():
        payload[key] = value
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1)
    return path


# ---------------------------------------------------------------------------
# merging with a user's own cluster workbook
# ---------------------------------------------------------------------------
def merge_cluster_workbooks(
    base_path: str,
    addition_path: str,
    out_path: str = "",
    base_sheet=0,
    addition_sheet=0,
) -> Dict[str, object]:
    """Combine two one-column-per-cluster workbooks into a third.

    CodonPipe's clustering reads a workbook with one column per gene cluster.
    The discovery run writes a workbook containing *only* its own virulence
    columns, so handing that file to the clustering step in place of an
    existing cluster workbook would silently discard every other cluster the
    user had defined. This builds the union instead.

    Neither input is modified, and an existing ``out_path`` is never
    overwritten -- a numeric suffix is added until the name is free. Column
    names that collide are disambiguated rather than merged, because two
    columns with the same heading in different workbooks are not necessarily
    the same gene set and silently unioning them would be a guess.
    """

    base = pd.read_excel(base_path, sheet_name=base_sheet)
    addition = pd.read_excel(addition_path, sheet_name=addition_sheet)

    taken = {str(c) for c in base.columns}
    renamed: Dict[str, str] = {}
    columns: "dict[str, List[str]]" = {
        str(c): [str(v) for v in base[c].dropna().astype(str) if str(v).strip()]
        for c in base.columns
    }
    for col in addition.columns:
        name = str(col)
        if name in taken:
            candidate = f"{name} (virulence discovery)"
            n = 2
            while candidate in taken:
                candidate = f"{name} (virulence discovery {n})"
                n += 1
            renamed[name] = candidate
            name = candidate
        taken.add(name)
        columns[name] = [
            str(v) for v in addition[col].dropna().astype(str) if str(v).strip()
        ]

    if not out_path:
        stem = os.path.splitext(base_path)[0]
        out_path = f"{stem} + virulence clusters.xlsx"
    root, ext = os.path.splitext(out_path)
    n = 2
    while os.path.exists(out_path):
        out_path = f"{root} ({n}){ext}"
        n += 1

    write_cluster_workbook(out_path, columns, sheet_name="Gene clusters")
    return {
        "path": out_path,
        "n_columns_base": int(len(base.columns)),
        "n_columns_added": int(len(addition.columns)),
        "n_columns_total": int(len(columns)),
        "renamed": renamed,
        "base_path": base_path,
        "addition_path": addition_path,
    }
