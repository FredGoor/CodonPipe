#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Orchestration of the automated virulence-gene discovery pipeline.

One call, :func:`run_discovery`, takes a target specification -- a species or
strain name, an NCBI taxonomy identifier or a genome assembly accession -- and
produces the two virulence-gene categories together with a complete audit
trail.  The module is deliberately free of any user-interface dependency: it
is driven identically from the CodonPipe desktop interface, from a script or
from the command line, which is what makes the result reproducible.

Usage from a script::

    from codonpipe.virulence import DiscoveryConfig, run_discovery

    result = run_discovery(DiscoveryConfig(
        target="Salmonella Typhimurium SL1344",
        stringency="standard",
        output_dir="virulence_output",
    ))
    print(result.outputs["clusters_xlsx"])

From the command line::

    python -m codonpipe.virulence "Salmonella Typhimurium SL1344" --stringency standard

The pipeline never writes into the directory holding the input genome and
never modifies any CodonPipe file; all artefacts go to ``output_dir``.
"""

from __future__ import annotations

import os
import time
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional

import pandas as pd

from . import classify as _classify
from . import evidence as _evidence
from . import outputs as _outputs
from . import reference as _reference
from . import resolve as _resolve
from . import species_loci as _species_loci
from .webcache import WebCache, OfflineCacheMiss, NetworkUnavailable

__all__ = ["MODULE_VERSION", "DiscoveryConfig", "DiscoveryResult", "run_discovery"]

MODULE_VERSION = "1.2.0"


@dataclass
class DiscoveryConfig:
    """Everything that defines one discovery run."""

    #: Species or strain name, NCBI taxid (``216597`` / ``txid:216597``) or
    #: assembly accession (``GCF_000210855.2``).
    target: str = ""
    #: Pin the assembly explicitly while still resolving taxonomy from ``target``.
    assembly: str = ""
    #: Coding-sequence FASTA defining the gene universe.  When empty the
    #: RefSeq CDS FASTA of the resolved assembly is downloaded.
    cds_fasta: str = ""
    #: Which tier the run treats as its headline set. All three tiers are
    #: always computed; this only affects which one is highlighted and which
    #: tier supplies the functional subcategory columns.
    stringency: str = "standard"
    #: Which tiers to write as cluster columns: ``all`` (default) or a
    #: comma-separated subset such as ``"standard,strict"``. The legacy values
    #: ``both``/``general``/``host`` are still accepted.
    tiers: str = "all"
    #: Deprecated alias of ``tiers``, retained for saved sessions.
    categories: str = ""
    #: Where the four output artefacts are written.
    output_dir: str = ""
    #: Cache folder; defaults to ``<output_dir>/virulence_cache``.
    cache_dir: str = ""
    offline: bool = False
    refresh: bool = False
    include_species_transfer: bool = True
    include_subcategories: bool = True
    min_subcategory_size: int = 2
    #: Optional experimental in vitro requirement table (xlsx/csv/tsv) with
    #: columns ``locus_tag`` and ``in_vitro_requirement``.
    essentiality_table: str = ""
    #: Load the bundled experimental gene-requirement table when one exists
    #: for the resolved organism. A table named in ``essentiality_table``
    #: always takes precedence over the bundled one.
    use_bundled_reference: bool = True
    #: What the host-interaction tier demands of Falkow's third postulate.
    #: Retired. Rung 3 has one definition; kept so an older script that passes
    #: it keeps working, and ignored.
    strict_requires: str = ""
    #: Optional contact details forwarded to NCBI Entrez, which asks for them.
    ncbi_email: str = ""
    ncbi_api_key: str = ""
    #: Apply the species-specific acquired-virulence-locus layer (E6, see
    #: :mod:`codonpipe.virulence.species_loci`). Harmless to leave on: the
    #: layer is a strict no-op for any genome outside the taxonomic scopes
    #: configured in ``resources/species_virulence_loci.tsv``.
    include_species_loci: bool = True
    #: Override the shipped resource folder to use edited rule tables.
    resource_dir: str = ""
    #: Explicit scoring overrides, applied on top of the stringency preset.
    scoring_overrides: Dict[str, object] = field(default_factory=dict)
    output_prefix: str = ""

    def resolved_output_dir(self) -> str:
        return os.path.abspath(
            os.path.expanduser(self.output_dir or os.path.join(os.getcwd(), "CodonPipe virulence"))
        )

    def resolved_cache_dir(self) -> str:
        if self.cache_dir:
            return os.path.abspath(os.path.expanduser(self.cache_dir))
        return os.path.join(self.resolved_output_dir(), "virulence_cache")


@dataclass
class DiscoveryResult:
    """Outcome of a run, including the data needed to inspect it."""

    config: DiscoveryConfig
    target: "_resolve.TargetInfo"
    gene_table: "_resolve.GeneTable"
    genes: pd.DataFrame
    summary: Dict[str, object]
    params: "_classify.ScoringParams"
    columns: Dict[str, List[str]]
    evidence_tables: Dict[str, pd.DataFrame]
    outputs: Dict[str, str] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)
    layer_counts: Dict[str, int] = field(default_factory=dict)
    #: Result of the species-specific acquired-virulence-locus layer (E6).
    #: Always present; ``in_scope`` is ``False`` for genomes outside the
    #: configured scopes, in which case the layer changed nothing.
    species_loci: Optional["_species_loci.SpeciesLociResult"] = None

    @property
    def virulence_score(self) -> Optional[int]:
        """The species-specific virulence score, when one is defined.

        For *Klebsiella* this is the Kleborate-style 0-5 score computed from
        the published locus definitions by annotation matching. It is **not**
        output of the Kleborate program: no MLST, no Kaptive K/O typing and no
        sequence-identity scoring are performed. ``None`` when the genome is
        outside every configured scope.
        """
        return None if self.species_loci is None else self.species_loci.score

    def tier(self, name: str) -> List[str]:
        """Locus tags in one tier of the ladder, in input order."""
        flag = f"tier_{str(name).strip().lower()}"
        if flag not in self.genes.columns:
            raise ValueError(f"Unknown tier {name!r}; choose from {', '.join(_classify.TIERS)}.")
        sub = self.genes[self.genes[flag].fillna(False).astype(bool)]
        return sub.sort_values("order", na_position="last")["locus_tag"].tolist()

    @property
    def permissive(self) -> List[str]:
        return self.tier("permissive")

    @property
    def standard(self) -> List[str]:
        return self.tier("standard")

    @property
    def strict(self) -> List[str]:
        return self.tier("strict")

    # Retained for continuity with CodonPipe 1.2.0 scripts.
    @property
    def virulence_general(self) -> List[str]:
        return self.tier("standard")

    @property
    def host_interaction(self) -> List[str]:
        return self.tier("strict")


_LEGACY_CATEGORY_TO_TIERS = {
    "both": ("permissive", "standard", "strict"),
    "general": ("standard",),
    "host": ("strict",),
}


#: Spellings accepted for a ladder rung, including the legacy three-tier names
#: and the two names the third rung goes by.
_TIER_SPELLINGS = {
    "permissive": "permissive",
    "standard": "standard",
    "strict": "host_interaction",
    "host-interaction": "host_interaction",
    "host_interaction": "host_interaction",
    "in-vivo-contributing": "in_vivo_contributing",
    "in_vivo_contributing": "in_vivo_contributing",
    "in-vivo-required": "in_vivo_required",
    "in_vivo_required": "in_vivo_required",
    "in-vitro-neutral": "in_vitro_neutral",
    "in_vitro_neutral": "in_vitro_neutral",
    "tier_1": "permissive",
    "tier_2": "host_interaction",
    "tier_3": "in_vivo_required",
    "tier_4": "in_vitro_neutral",
}


def _requested_tiers(config: "DiscoveryConfig") -> tuple:
    """Resolve which ladder rungs a run should report.

    Accepts the five rung names, the legacy three-tier names, and the legacy
    category strings. Rungs that are empty for the genome in hand (4 and 5
    without mutant-screen data) are dropped when the workbook is written, so
    asking for all five is always safe.
    """
    raw = str(getattr(config, "tiers", "") or "").strip().lower()
    if raw in {"published", "4-tier", "four", "tiers"}:
        return _classify.PUBLISHED_TIERS
    if not raw:
        # Default: the four published tiers. The two gates are computed either
        # way and can be written by asking for "all"; they are left out here so
        # a default run produces a four-tier workbook rather than six columns a
        # reader has to know how to discount.
        legacy = str(getattr(config, "categories", "") or "").strip().lower()
        if legacy in _LEGACY_CATEGORY_TO_TIERS:
            return tuple(
                _TIER_SPELLINGS.get(t, t) for t in _LEGACY_CATEGORY_TO_TIERS[legacy]
            )
        return _classify.PUBLISHED_TIERS
    if raw == "all":
        legacy = str(getattr(config, "categories", "") or "").strip().lower()
        if legacy in _LEGACY_CATEGORY_TO_TIERS:
            return tuple(
                _TIER_SPELLINGS.get(t, t) for t in _LEGACY_CATEGORY_TO_TIERS[legacy]
            )
        return _classify.LADDER
    wanted = [t.strip() for t in raw.replace(";", ",").split(",") if t.strip()]
    unknown = [t for t in wanted if t not in _TIER_SPELLINGS]
    if unknown:
        raise ValueError(
            f"Unknown tier(s) {sorted(unknown)}; choose from "
            f"{', '.join(_classify.LADDER)} (or the legacy names permissive, "
            f"standard, strict)."
        )
    resolved = {_TIER_SPELLINGS[t] for t in wanted}
    return tuple(r for r in _classify.LADDER if r in resolved)


def _read_table(path: str) -> pd.DataFrame:
    lower = str(path).lower()
    if lower.endswith((".xlsx", ".xlsm", ".xls")):
        return pd.read_excel(path)
    if lower.endswith((".tsv", ".tab", ".txt")):
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)



def _is_template_workbook(path: str) -> bool:
    """True when `path` is a filled-in user template rather than a ready table.

    The marker is a "Screen definitions" sheet. Non-Excel files and workbooks
    that cannot be opened return False and are read as canonical tables, so a
    malformed file still raises in the place that reports it properly.
    """
    if os.path.splitext(path)[1].lower() not in (".xlsx", ".xlsm", ".xls"):
        return False
    try:
        return "Screen definitions" in pd.ExcelFile(path).sheet_names
    except Exception:
        return False


def run_discovery(
    config: DiscoveryConfig,
    log: Optional[Callable[[str], None]] = None,
    write_outputs: bool = True,
) -> DiscoveryResult:
    """Run the full discovery pipeline for one target.

    Parameters
    ----------
    config:
        See :class:`DiscoveryConfig`.
    log:
        Optional ``callable(str)`` receiving progress lines.  The desktop
        interface passes its console logger here.
    write_outputs:
        Set to ``False`` to compute the result without writing files, which is
        what the benchmark harness does.
    """
    emit = log if callable(log) else (lambda _m: None)
    started = time.time()
    warnings: List[str] = []
    out_dir = config.resolved_output_dir()
    cache_dir = config.resolved_cache_dir()

    web = WebCache(
        cache_dir,
        offline=config.offline,
        refresh=config.refresh,
        log=emit,
        user_agent=f"CodonPipe-virulence-discovery/{MODULE_VERSION}",
    )
    resource_dir = config.resource_dir or _evidence.RESOURCE_DIR

    # ------------------------------------------------------------- resources
    emit("[1/7] Loading classification rule tables ...")
    go_axes = _evidence.load_go_axes(resource_dir)
    keyword_axes = _evidence.load_keyword_axes(resource_dir)
    vf_category_axes = _evidence.load_vf_category_axes(resource_dir)
    vocabulary = _evidence.load_vocabulary(resource_dir)
    core_vocabulary = _evidence.load_core_vocabulary(resource_dir)
    locus_definitions = (
        _species_loci.load_species_loci(resource_dir) if config.include_species_loci else []
    )
    emit(
        f"        {len(go_axes)} Gene Ontology roots, {len(keyword_axes)} UniProt keywords, "
        f"{len(vf_category_axes)} virulence-factor categories, "
        f"{len(vocabulary)} vocabulary rules, {len(core_vocabulary)} core-process rules, "
        f"{len(locus_definitions)} species-specific virulence-locus definitions."
    )

    # ---------------------------------------------------------------- target
    emit("[2/7] Resolving target against NCBI Taxonomy and Assembly ...")
    target = _resolve.resolve_target(
        config.target,
        web,
        assembly_override=config.assembly,
        ncbi_email=config.ncbi_email,
        ncbi_api_key=config.ncbi_api_key,
    )
    emit(
        f"        taxon {target.tax_id} - {target.scientific_name}; "
        f"assembly {target.assembly_accession} ({target.assembly_level})."
    )
    for note in target.resolution_notes:
        emit(f"        {note}")

    # Vocabulary rules that name a taxon apply only inside that lineage. The
    # SPI-1 and SPI-2 gene-symbol rules are the case this exists for: sseA-sseL
    # and sopB are Salmonella island symbols, and applying them to Sakai, UTI89
    # or KPPR1 labelled unrelated genes as SPI-2.
    vocabulary, _vocab_notes = _evidence.restrict_vocabulary_to_taxon(
        vocabulary,
        target.lineage_tax_ids
        or (target.tax_id, target.species_tax_id, target.genus_tax_id),
    )
    for note in _vocab_notes:
        emit(f"        {note}")
        warnings.append(note)

    # ------------------------------------------------------------ annotation
    emit("[3/7] Retrieving the RefSeq annotation and building the gene table ...")
    gff_frame = None
    accession = target.refseq_accession or target.assembly_accession
    try:
        gff_frame = _resolve.parse_gff_genes(_resolve.fetch_refseq_gff(accession, web))
        emit(f"        {len(gff_frame)} annotated genes in {accession}.")
    except (NetworkUnavailable, OfflineCacheMiss) as exc:
        warnings.append(f"RefSeq annotation unavailable for {accession}: {exc}")
        emit(f"        annotation unavailable: {exc}")

    cds_path = os.path.abspath(os.path.expanduser(config.cds_fasta)) if config.cds_fasta else ""
    cds_frame = None
    if cds_path and os.path.isfile(cds_path):
        cds_frame = _resolve.parse_cds_fasta(cds_path)
        emit(f"        gene universe from the supplied FASTA: {len(cds_frame)} coding sequences.")
    else:
        if cds_path:
            warnings.append(f"The supplied CDS FASTA was not found: {cds_path}")
        dest = os.path.join(cache_dir, "genomes", f"{accession}_cds_from_genomic.fna")
        try:
            cds_path = _resolve.fetch_cds_fasta(accession, web, dest)
            cds_frame = _resolve.parse_cds_fasta(cds_path)
            emit(f"        downloaded RefSeq coding sequences: {len(cds_frame)} records.")
        except (NetworkUnavailable, OfflineCacheMiss) as exc:
            warnings.append(f"Could not obtain coding sequences for {accession}: {exc}")
            emit(f"        coding-sequence download failed: {exc}")
            cds_path = ""

    gene_table = _resolve.build_crosswalk(cds_frame, gff_frame)
    warnings.extend(gene_table.notes)
    emit(f"        gene table: {len(gene_table.frame)} genes ({gene_table.source}).")

    # ----------------------------------------------------------- evidence E1
    emit("[4/7] Harvesting curated virulence-factor evidence (VFDB, VICTORS, PATRIC_VF) ...")
    overrides = dict(config.scoring_overrides or {})
    # An explicit scoring override wins over the top-level convenience field,
    # so a caller that sets both gets the one it spelled out.
    params = _classify.params_for_stringency(config.stringency, **overrides)
    genome_ids, bv_notes = _evidence.bvbrc_genome_ids(target, web)
    for note in bv_notes:
        emit(f"        {note}")
    if not genome_ids:
        warnings.append(
            "No BV-BRC genome record matched this assembly, so VFDB/VICTORS/PATRIC_VF "
            "evidence is missing. Results rest on UniProtKB and annotation evidence only."
        )
    vf_raw = _evidence.fetch_bvbrc_virulence(genome_ids, web) if genome_ids else pd.DataFrame()
    vf_frame, vf_notes = _evidence.map_bvbrc_to_genes(
        vf_raw,
        gene_table,
        genome_ids,
        web,
        min_identity=params.min_identity,
        min_coverage=params.min_coverage,
    )
    warnings.extend(vf_notes)
    for note in vf_notes:
        emit(f"        {note}")
    emit(
        f"        {len(vf_frame)} usable records on "
        f"{vf_frame['locus_tag'].nunique() if not vf_frame.empty else 0} genes."
    )

    # ----------------------------------------------------------- evidence E2
    emit("[5/7] Harvesting UniProtKB curation and Gene Ontology closures ...")
    uniprot_frame, up_notes = _evidence.map_uniprot_to_genes(
        target, gene_table, web, include_species_transfer=config.include_species_transfer
    )
    for note in up_notes:
        emit(f"        {note}")
    go_closures, go_warnings = _evidence.fetch_go_closures(go_axes, web)
    warnings.extend(go_warnings)
    for item in go_warnings:
        emit(f"        [warning] {item}")
    emit(
        "        Gene Ontology closures: "
        + ", ".join(f"{k} = {len(v)} terms" for k, v in sorted(go_closures.items()))
    )

    # ------------------------------------------------------- evidence E4/E5
    emit("[6/7] Applying the controlled vocabularies ...")
    vocabulary_hits = _evidence.apply_vocabulary(gene_table, vocabulary)
    core_hits = _evidence.apply_core_vocabulary(gene_table, core_vocabulary)
    emit(
        f"        {vocabulary_hits['locus_tag'].nunique() if not vocabulary_hits.empty else 0} genes "
        f"matched the virulence lexicon; "
        f"{core_hits['locus_tag'].nunique() if not core_hits.empty else 0} matched the core-process lexicon."
    )

    # ------------------------------------------------------------ evidence E6
    emit("        Applying species-specific acquired virulence-locus definitions ...")
    species_loci = _species_loci.detect_species_loci(
        gene_table, target, resource_dir=resource_dir, definitions=locus_definitions
    )
    for note in species_loci.notes:
        emit(f"        {note}")
    if not config.include_species_loci:
        emit("        the species-locus layer was switched off for this run.")

    # An experimental gene-requirement table. A table the user points at wins
    # over a bundled one; a bundled one is offered only when its recorded
    # taxon or assembly matches the genome that was actually resolved, never
    # on a name resemblance.
    essentiality = None
    essentiality_provenance = ""
    if config.essentiality_table:
        path = os.path.abspath(os.path.expanduser(config.essentiality_table))
        try:
            # Two accepted shapes. A workbook carrying a "Screen definitions"
            # sheet is the USER TEMPLATE: raw screen columns plus a declaration
            # of what each one means, which the reference builder interprets
            # into the controlled vocabularies. Anything else is read as a
            # table already in canonical form. Detecting by sheet rather than
            # asking the user to say which is which removes a step they would
            # otherwise get wrong on first use.
            if _is_template_workbook(path):
                built = _reference.build_reference_from_template(path)
                essentiality = built["frame"]
                for note in built.get("notes", []):
                    warnings.append(note)
                emit(f"        experimental gene-requirement table: built from the "
                     f"template {os.path.basename(path)} "
                     f"({built['n_columns_declared']} declared screen column(s), "
                     f"{built['n_attenuated']} attenuated, "
                     f"{built['n_essential']} essential in vitro).")
                essentiality_provenance = (
                    f"user template {path} "
                    f"({built['n_columns_declared']} declared screen columns)"
                )
            else:
                essentiality = _reference.validate_reference(_read_table(path), source=path)
                emit(f"        experimental gene-requirement table: {len(essentiality)} rows "
                     f"from {os.path.basename(path)}.")
                essentiality_provenance = f"user-supplied table {path}"
        except Exception as exc:
            warnings.append(f"Could not read the essentiality table {path}: {exc}")
            emit(f"        [warning] essentiality table unreadable: {exc}")
    elif config.use_bundled_reference:
        entry = _reference.find_reference(
            taxon_id=str(target.tax_id or ""),
            species_taxon_id=str(target.species_tax_id or ""),
            assembly=str(target.assembly_accession or ""),
            organism=str(target.organism_label or ""),
        ) or _reference.find_reference(assembly=str(target.paired_accession or ""))
        if entry:
            try:
                essentiality = _reference.load_reference(entry)
                cov = _reference.reference_coverage(essentiality)
                essentiality_provenance = (
                    f"bundled reference '{entry['key']}' ({os.path.basename(entry['path'])})"
                )
                emit(f"        bundled experimental reference for {entry['organism']}: "
                     f"{cov['rows']} genes, {cov['essential']} essential in vitro, "
                     f"{cov['attenuated']} attenuated in vivo, "
                     f"{cov['no_phenotype']} assayed without a phenotype.")
            except Exception as exc:
                warnings.append(f"Bundled reference {entry['key']} unreadable: {exc}")
                emit(f"        [warning] bundled reference unreadable: {exc}")

    # ------------------------------------------------------------- classify
    emit("[7/7] Scoring and assigning the six-rung ladder ...")
    genes, summary = _classify.classify_genes(
        gene_table,
        vf_frame,
        uniprot_frame,
        vocabulary_hits,
        core_hits,
        go_closures,
        keyword_axes,
        vf_category_axes,
        params,
        essentiality=essentiality,
        species_locus_hits=species_loci.hits,
    )
    summary["experimental_reference"] = essentiality_provenance or "none"
    columns = _outputs.build_cluster_columns(
        genes,
        tiers=_requested_tiers(config),
        include_subcategories=config.include_subcategories,
        min_subcategory_size=config.min_subcategory_size,
        subcategory_tier=params.stringency,
    )
    for n, rung in enumerate(_classify.LADDER, start=1):
        count = int(genes[f"tier_{rung}"].sum())
        note = ""
        if rung in ("in_vivo_contributing", "in_vivo_required", "in_vitro_neutral") and count == 0:
            if int(genes["tier_in_vivo_contributing"].sum()) == 0:
                note = "   [needs experimental mutant-screen data for this organism]"
            else:
                note = "   [needs per-gene effect sizes; the table reports calls without magnitudes]"
        elif rung == "in_vitro_neutral" and not summary.get("in_vitro_fitness_data_present", False) and count:
            note = "   [identical to rung 5: no graded in vitro fitness data]"
        emit(
            f"        {n}. {rung:<18s} {count:5d} genes{note}"
        )
        emit(f"              {_classify.LADDER_DEFINITIONS[rung]}")
    emit(f"        ladder nested (asserted): {summary.get('tier_ladder_nested')}")
    # Rungs 4 to 6 inherit the calling threshold of whichever screen supplied
    # the organism's table, and those thresholds differ by more than an order
    # of magnitude between published studies. Their counts therefore describe
    # the screen as much as the pathogen, and must not be read as a measure of
    # how specialised one organism is relative to another. Rung 3 carries no
    # such dependence and is the axis to compare across genomes.
    call_rate = summary.get("in_vivo_screen_call_rate")
    if call_rate is not None:
        warnings.append(
            f"This organism's screens called {call_rate}% of the genes they assayed "
            f"attenuated ({summary.get('in_vivo_attenuated')} of "
            f"{summary.get('in_vivo_assayed')}). That rate is set by the original "
            "authors' thresholds and ranges from under 2% to nearly 30% across the "
            "bundled tables, so rung 4, 5 and 6 COUNTS ARE NOT COMPARABLE BETWEEN "
            "ORGANISMS. Compare rung 3, which is computed identically for every "
            "genome, or compare rungs 4 to 6 only within one organism."
        )
        emit(
            f"        in vivo screen call rate: {call_rate}% of assayed genes "
            "- rungs 4-6 are within-organism measures only"
        )
    for name, members in columns.items():
        emit(f"        column '{name}': {len(members)} genes")

    layer_counts = {
        "Virulence-factor database records retained": int(len(vf_frame)),
        "Genes with virulence-factor database evidence": int(
            vf_frame["locus_tag"].nunique() if not vf_frame.empty else 0
        ),
        "Literature-curated database records": int(
            vf_frame["is_literature"].sum() if "is_literature" in vf_frame.columns else 0
        ),
        "UniProtKB entries mapped": int(len(uniprot_frame)),
        "Genes with the UniProtKB Virulence keyword": int(
            summary.get("n_uniprot_virulence_keyword", 0)
        ),
        "Genes matching the virulence vocabulary": int(
            vocabulary_hits["locus_tag"].nunique() if not vocabulary_hits.empty else 0
        ),
        "Genes matching the core-process vocabulary": int(
            core_hits["locus_tag"].nunique() if not core_hits.empty else 0
        ),
        "Genes added by genomic context": int(summary.get("n_context_added", 0)),
        "Genes in a species-specific acquired virulence locus": int(
            summary.get("n_species_locus_genes", 0)
        ),
        "- of which matched a published locus definition directly": int(
            summary.get("n_species_locus_marker_genes", 0)
        ),
    }
    summary.update(species_loci.summary_dict())

    evidence_tables = {
        "Virulence factor records": vf_frame,
        "UniProt annotation": uniprot_frame,
        "Vocabulary matches": vocabulary_hits,
        "Core process matches": core_hits,
        "Species virulence loci": species_loci.as_frame(),
        "Species locus genes": species_loci.hits,
    }

    result = DiscoveryResult(
        config=config,
        target=target,
        gene_table=gene_table,
        genes=genes,
        summary=summary,
        params=params,
        columns=columns,
        evidence_tables=evidence_tables,
        warnings=warnings,
        layer_counts=layer_counts,
        species_loci=species_loci,
    )

    if write_outputs:
        os.makedirs(out_dir, exist_ok=True)
        prefix = config.output_prefix or _outputs.safe_name(
            target.organism_label or target.scientific_name or config.target
        )
        stem = f"{prefix} virulence"
        clusters_path = os.path.join(out_dir, f"{stem} clusters.xlsx")
        evidence_path = os.path.join(out_dir, f"{stem} evidence.xlsx")
        report_path = os.path.join(out_dir, f"{stem} discovery report.txt")
        manifest_path = os.path.join(out_dir, f"{stem} run.json")

        # ---- provenance-report section for the species-locus layer -------
        locus_section: List[str] = []
        if not species_loci.in_scope:
            locus_section.append(
                "No species-specific acquired-virulence-locus definition is scoped to "
                "this genome, so this layer contributed nothing: no gene score, no tier"
            )
            locus_section.append(
                "assignment and no reported virulence score were affected by it."
            )
            locus_section.append("")
            locus_section.extend(f"  {n}" for n in species_loci.notes)
        else:
            locus_section.append(
                "Published locus definitions were applied to the RefSeq product names "
                "through the identifier crosswalk and the operon/island machinery."
            )
            locus_section.append(f"  Scope matched           : {species_loci.scope_description}")
            locus_section.append(
                f"  Resource table          : {_species_loci.SPECIES_LOCI_TABLE}"
            )
            locus_section.append("")
            locus_section.append(
                f"  {'locus':<10s}{'present':<9s}{'markers':<9s}{'recruited':<11s}"
                f"{'complete':<10s}{'scores':<8s}machinery"
            )
            for call in species_loci.calls:
                locus_section.append(
                    f"  {call.locus:<10s}{str(call.present):<9s}"
                    f"{str(call.n_markers) + '/' + str(call.expected_gene_count):<9s}"
                    f"{len(call.recruited_genes):<11d}{call.completeness:<10.2f}"
                    f"{str(bool(call.counts_toward_score)):<8s}{call.machinery}"
                )
            locus_section.append("")
            for call in species_loci.calls:
                if not call.present:
                    continue
                locus_section.append(f"  {call.locus} - {call.locus_label}")
                locus_section.append(f"    replicon              : {call.seqid}")
                locus_section.append(
                    f"    genes matching the published definition ({len(call.marker_genes)}):"
                )
                for tag, product in zip(call.marker_genes, call.marker_products):
                    locus_section.append(f"      {tag:<16s} {product[:92]}")
                if call.recruited_genes:
                    locus_section.append(
                        f"    recruited within the locus block ({len(call.recruited_genes)}):"
                    )
                    for tag in call.recruited_genes:
                        locus_section.append(f"      {tag}")
                locus_section.append(f"    reference             : {call.reference}")
                locus_section.append("")
            for call in species_loci.calls:
                if call.present or call.n_markers == 0:
                    continue
                locus_section.append(
                    f"  {call.locus} NOT called: {call.reason}"
                )
                for tag, product in zip(call.marker_genes, call.marker_products):
                    locus_section.append(f"      {tag:<16s} {product[:92]}")
                locus_section.append("")
            if species_loci.score is not None:
                locus_section.append(
                    f"  SCORE ({species_loci.score_scheme}) : {species_loci.score}"
                )
                locus_section.append(
                    f"    contributing loci     : "
                    f"{', '.join(species_loci.score_tokens) or 'none'}"
                )
                locus_section.append("")
                locus_section.extend(
                    f"  {line}"
                    for line in _species_loci.SCORE_SCHEME_RULES.get(
                        species_loci.score_scheme, ()
                    )
                )
                locus_section.append("")
            locus_section.append(
                "  WHAT THIS LAYER IS NOT: the Klebsiella definitions reproduce the five"
            )
            locus_section.append(
                "  acquired loci and the 0-5 score rule published for Kleborate (Lam et al."
            )
            locus_section.append(
                "  2021, Nat Commun 12:4188), applied by ANNOTATION MATCHING ONLY. This is"
            )
            locus_section.append(
                "  not the Kleborate program: no MLST, no Kaptive K-locus or O-locus typing,"
            )
            locus_section.append(
                "  no sequence-identity or coverage scoring against the Kleborate allele"
            )
            locus_section.append(
                "  databases and therefore no ybST/cbST/AbST/SmST lineage assignment, and no"
            )
            locus_section.append(
                "  external binary (minimap2, mash, BLAST, Kaptive) of any kind. A locus is"
            )
            locus_section.append(
                "  called from the presence of ANNOTATED genes, so a sparsely annotated,"
            )
            locus_section.append(
                "  pseudogenised or truncated locus may be called differently from how"
            )
            locus_section.append(
                "  Kleborate, which reads the nucleotide sequence, would call it. Report any"
            )
            locus_section.append(
                "  value from this layer as a \"Kleborate-style\" score, never as Kleborate"
            )
            locus_section.append("  output.")

        summary_rows = pd.DataFrame(
            [{"metric": k, "value": v} for k, v in layer_counts.items()]
            + [{"metric": f"cluster: {k}", "value": len(v)} for k, v in columns.items()]
            + [
                {"metric": f"tier: {t}", "value": int(genes[f"tier_{t}"].sum())}
                for t in _classify.TIERS
            ]
            + [
                {"metric": f"rung {i}: {name}", "value": int(genes[f"tier_{name}"].sum())}
                for i, name in enumerate(_classify.LADDER, start=1)
            ]
            # Caveats the classifier raised about this run -- notably whether
            # rung 6 could separate anything at all. A reader of the workbook
            # must see these next to the counts, not only in the text report.
            + [
                {"metric": f"caveat {i}", "value": note}
                for i, note in enumerate(summary.get("notes", []) or [], start=1)
            ]
            + (
                [
                    {
                        "metric": "in vivo discordant (attenuated call vs negative effect)",
                        "value": int(genes["in_vivo_discordant"].sum()),
                    },
                    {
                        "metric": "in vivo discordant, inside rung 2",
                        "value": int(
                            (
                                genes["in_vivo_discordant"].fillna(False).astype(bool)
                                & genes["tier_standard"].fillna(False).astype(bool)
                            ).sum()
                        ),
                    },
                ]
                if "in_vivo_discordant" in genes.columns
                else []
            )
            + [
                {"metric": "headline tier", "value": params.stringency},
                {"metric": "threshold permissive", "value": params.threshold_permissive},
                {"metric": "threshold standard", "value": params.threshold_standard},
                # The three in vivo cut-offs are conventions the operator can
                # change, so each run records the ones it actually used. A
                # reader who does not have the interface in front of them can
                # still tell what line was drawn.
                {"metric": "in vivo effect scale", "value": params.effect_scale},
                {
                    "metric": "in vivo effect cut-off, contributing",
                    "value": summary.get("effect_contributing_log2"),
                },
                {
                    "metric": "in vivo effect cut-off, required",
                    "value": summary.get("effect_required_log2"),
                },
                {
                    "metric": "in vivo effect cut-off basis",
                    "value": summary.get("effect_scale_note"),
                },
                {
                    "metric": "in vivo significance cut-off",
                    "value": params.max_in_vivo_significance,
                },
                {
                    "metric": "rung 4 genes also in rung 3",
                    "value": summary.get("n_rung4_also_in_rung3"),
                },
                {
                    "metric": "attenuated in vivo but outside rung 3",
                    "value": summary.get("n_in_vivo_outside_rung3"),
                },
                {
                    "metric": "rung 3 genes without in vivo support",
                    "value": summary.get("n_rung3_without_in_vivo_support"),
                },
                {"metric": "in vivo genes assayed", "value": summary.get("in_vivo_assayed")},
                {
                    "metric": "in vivo screen call rate (%)",
                    "value": summary.get("in_vivo_screen_call_rate"),
                },
                {"metric": "tier ladder nested", "value": summary.get("tier_ladder_nested")},
                {"metric": "assembly", "value": target.assembly_accession},
                {"metric": "species virulence loci in scope", "value": species_loci.in_scope},
                {
                    "metric": "species virulence loci present",
                    "value": ", ".join(species_loci.loci_present) or "none",
                },
                {
                    "metric": f"virulence score ({species_loci.score_scheme or 'n/a'})",
                    "value": "not applicable" if species_loci.score is None else species_loci.score,
                },
                {"metric": "module version", "value": MODULE_VERSION},
            ]
        )
        _outputs.write_cluster_workbook(clusters_path, columns)
        _outputs.write_evidence_workbook(evidence_path, genes, evidence_tables, summary_rows)
        _outputs.write_report(
            report_path,
            target,
            params,
            summary,
            columns,
            gene_table.source,
            layer_counts,
            warnings,
            web.summary(),
            web.accesses,
            extra_sections={
                "Species-specific acquired virulence loci (evidence layer E6)": locus_section
            },
        )
        outputs = {
            "clusters_xlsx": clusters_path,
            "evidence_xlsx": evidence_path,
            "report_txt": report_path,
            "manifest_json": manifest_path,
            "cds_fasta": cds_path,
            "cache_dir": cache_dir,
        }
        _outputs.write_manifest(
            manifest_path,
            target,
            params,
            summary,
            columns,
            resource_dir,
            web.accesses,
            outputs,
            warnings,
            MODULE_VERSION,
            extra={
                "species_virulence_loci": {
                    "in_scope": bool(species_loci.in_scope),
                    "scope": species_loci.scope_description,
                    "resource_table": _species_loci.SPECIES_LOCI_TABLE,
                    "score_scheme": species_loci.score_scheme,
                    "score": species_loci.score,
                    "score_rule": list(
                        _species_loci.SCORE_SCHEME_RULES.get(species_loci.score_scheme, ())
                    ),
                    "score_tokens": list(species_loci.score_tokens),
                    "loci_present": list(species_loci.loci_present),
                    "calls": species_loci.as_frame().to_dict(orient="records"),
                    "is_kleborate_program": False,
                    "disclaimer": (
                        "Published Kleborate locus definitions (Lam MMC et al. 2021, "
                        "Nat Commun 12:4188) applied by annotation matching. NOT the "
                        "Kleborate program: no MLST, no Kaptive K/O typing, no "
                        "sequence-identity scoring, no external binaries."
                    ),
                    "notes": list(species_loci.notes),
                }
            },
        )
        result.outputs = outputs
        emit(f"        outputs written to {out_dir}")

    emit(f"Completed in {time.time() - started:.1f} s.")
    return result
