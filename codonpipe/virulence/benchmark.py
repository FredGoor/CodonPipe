#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Concordance benchmark against manually curated virulence-gene sets.

What this measures, and what it does not
----------------------------------------
There is no gold-standard virulence gene set for any bacterium, so this is
**not** a precision/recall evaluation against truth.  What it measures is
*concordance* with a previously published manual curation, decomposed so that
every disagreement is attributable:

``recall_of_reference``
    the fraction of the manually curated genes the automated pipeline
    recovers.  High recall means the automated set is not losing genes the
    curators considered important.

``concordance``
    the fraction of automated calls that appear in the manual reference.  This
    is deliberately *not* called precision: a call outside the reference may be
    a false positive, but it may equally be a gene the curators placed in a
    different column of their workbook (motility, iron uptake, T6SS), or one
    they simply did not include.  The breakdown below separates these.

``in_other_curated_column`` / ``in_curated_housekeeping``
    where the extra calls land.  Some of the curated workbooks carry the
    authors' own ``Housekeeping / stress`` exclusion lists, and for those the
    second number is a genuine false-positive rate against the curators' own
    negative set -- the sharpest available test of the core-process filter.
    Only three of the six default cases declare such a column (EHEC O157:H7
    Sakai, UPEC UTI89 and *S. flexneri* 2a); for the others the field is
    reported as ``None``, meaning *not available*, which is not the same as
    zero.

In-sample caveat
----------------
The *Salmonella* SL1344 case is **not** held out: the default score weights
and the three stringency thresholds were calibrated on exactly that curated
set.  Its concordance is therefore in-sample and will be optimistic relative
to the other five, which were never used for tuning.  Read the SL1344 row as
an upper bound and the remaining rows as the out-of-sample estimate.

Running it
----------
>>> from codonpipe.virulence.benchmark import default_cases, run_benchmark
>>> summary, details = run_benchmark(default_cases("Preloaded genomes"))

The default cases are the five Enterobacterales pathogens and
*Staphylococcus aureus* analysed in the CodonPipe manuscripts, read from the
``Preloaded genomes`` folder shipped with CodonPipe.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Set, Tuple

import pandas as pd

from .classify import TIERS, TIER_DEFINITIONS
from .pipeline import DiscoveryConfig, run_discovery

__all__ = ["BenchmarkCase", "DEFAULT_CASES", "default_cases", "run_benchmark",
           "write_benchmark_workbook"]


@dataclass
class BenchmarkCase:
    """One species: a genome, an assembly, and a curated reference set."""

    name: str
    assembly: str
    fasta_name: str
    workbook_name: str
    reference_column: str
    housekeeping_column: str = ""
    other_columns: Tuple[str, ...] = ()
    fasta_path: str = ""
    workbook_path: str = ""


#: The species analysed in the CodonPipe manuscripts, with the curated column
#: that defines the manual virulence set in each workbook.
DEFAULT_CASES: Tuple[BenchmarkCase, ...] = (
    BenchmarkCase(
        name="Salmonella enterica SL1344",
        assembly="GCF_000210855.2",
        fasta_name="Salmonella enterica SL1344.fna",
        workbook_name="Salmonella enterica SL1344 clusters.xlsx",
        reference_column="Virulence (full)",
        other_columns=("SPI-1", "SPI-2", "O-antigen", "Adhesins / fimbriae / curli",
                       "Toxins", "T6SS", "Flagellum / motility / chemotaxis",
                       "Siderophores / iron uptake"),
    ),
    BenchmarkCase(
        name="Escherichia coli O157:H7 Sakai",
        assembly="GCF_000008865.2",
        fasta_name="Escherichia coli O157-H7 Sakai.fna",
        workbook_name="Escherichia coli O157-H7 Sakai clusters.xlsx",
        reference_column="Virulence (all)",
        housekeeping_column="Housekeeping / stress",
        other_columns=("T3SS", "Toxins", "Adhesin / fimbriae / pilus / curli", "T6SS",
                       "LPS / O-antigen / capsule", "Flagellum / motility / chemotaxis",
                       "Siderophores / iron uptake", "T2SS"),
    ),
    BenchmarkCase(
        name="Escherichia coli UPEC UTI89",
        assembly="GCF_000013265.1",
        fasta_name="Escherichia coli UPEC UTI89.fna",
        workbook_name="Escherichia coli UPEC UTI89 clusters.xlsx",
        reference_column="Virulence (all)",
        housekeeping_column="Housekeeping / stress",
        other_columns=("Toxins", "Adhesin / fimbriae / pilus / curli",
                       "LPS / O-antigen / capsule", "T6SS",
                       "Flagellum / motility / chemotaxis", "Siderophores / iron uptake"),
    ),
    BenchmarkCase(
        name="Shigella flexneri 2a",
        assembly="GCF_000006925.2",
        fasta_name="Shigella flexneri 2a.fna",
        workbook_name="Shigella flexneri 2a clusters.xlsx",
        reference_column="Virulence (all)",
        housekeeping_column="Housekeeping / stress",
        other_columns=("T3SS", "Adhesin / fimbriae / pilus / curli",
                       "LPS / O-antigen / capsule", "Flagellum / motility / chemotaxis",
                       "Siderophores / iron uptake"),
    ),
    BenchmarkCase(
        name="Klebsiella pneumoniae HS11286",
        assembly="GCF_000240185.1",
        fasta_name="Klebsiella pneumoniae.fna",
        workbook_name="Klebsiella pneumoniae clusters.xlsx",
        reference_column="Virulence (all)",
        other_columns=("O-antigen & Capsule", "Toxins - Yersiniabactin", "T6SS",
                       "Adhesion / fimbriae / pili", "Enterobactin / ferrienterobactin",
                       "Flagella / motility / chemotaxis", "Iron acquisition / siderophore"),
    ),
    BenchmarkCase(
        name="Staphylococcus aureus NCTC 8325",
        assembly="GCF_000013425.1",
        fasta_name="S. aureus NCTC8325.fna",
        workbook_name="S. aureus NCTC8325 clusters.xlsx",
        reference_column="Virulence",
        other_columns=("Secretion", "Adhesion"),
    ),
)


def default_cases(genomes_dir: str) -> List[BenchmarkCase]:
    """Return the default cases whose genome and workbook are both present."""
    folder = os.path.abspath(os.path.expanduser(genomes_dir))
    out: List[BenchmarkCase] = []
    for case in DEFAULT_CASES:
        fasta = os.path.join(folder, case.fasta_name)
        workbook = os.path.join(folder, case.workbook_name)
        if os.path.isfile(fasta) and os.path.isfile(workbook):
            out.append(
                BenchmarkCase(
                    **{**case.__dict__, "fasta_path": fasta, "workbook_path": workbook}
                )
            )
    return out


def _column_set(workbook: str, column: str) -> Set[str]:
    try:
        frame = pd.ExcelFile(workbook).parse(0)
    except Exception:
        return set()
    if column not in frame.columns:
        return set()
    return {str(v).strip() for v in frame[column].dropna() if str(v).strip()}


def _metrics(called: Set[str], reference: Set[str]) -> Dict[str, float]:
    overlap = called & reference
    recall = len(overlap) / len(reference) if reference else float("nan")
    concordance = len(overlap) / len(called) if called else float("nan")
    if recall and concordance and (recall + concordance):
        harmonic = 2 * recall * concordance / (recall + concordance)
    else:
        harmonic = float("nan")
    return {
        "n_called": len(called),
        "n_reference": len(reference),
        "n_overlap": len(overlap),
        "recall_of_reference": round(recall, 4),
        "concordance": round(concordance, 4),
        "harmonic_mean": round(harmonic, 4),
    }


def run_benchmark(
    cases: Sequence[BenchmarkCase],
    tiers: Sequence[str] = TIERS,
    cache_dir: str = "",
    log: Optional[Callable[[str], None]] = None,
    **config_kwargs,
) -> Tuple[pd.DataFrame, Dict[str, pd.DataFrame]]:
    """Score the tier ladder against each curated reference set.

    Because all three tiers come out of one scoring pass, each species is
    analysed once.  Returns ``(summary, details)``: one summary row per species
    and tier, and a per-species frame listing every disagreement with the
    reason the classifier gave.
    """
    emit = log if callable(log) else (lambda _m: None)
    rows: List[dict] = []
    details: Dict[str, pd.DataFrame] = {}

    for case in cases:
        reference = _column_set(case.workbook_path, case.reference_column)
        housekeeping = (
            _column_set(case.workbook_path, case.housekeeping_column)
            if case.housekeeping_column
            else set()
        )
        other: Dict[str, Set[str]] = {
            col: _column_set(case.workbook_path, col) for col in case.other_columns
        }
        other = {k: v for k, v in other.items() if v}
        hk_note = (
            f"{len(housekeeping)} curator-excluded housekeeping genes"
            if case.housekeeping_column
            else "no curator housekeeping list in this workbook"
        )
        in_sample = (
            " [in-sample: thresholds were calibrated on this set]"
            if case.name.startswith("Salmonella enterica SL1344")
            else ""
        )
        emit(
            f"[{case.name}] reference '{case.reference_column}' = {len(reference)} genes; "
            f"{len(other)} further curated columns; {hk_note}.{in_sample}"
        )
        if not reference:
            emit(f"[{case.name}] reference column not found; case skipped.")
            continue

        # The benchmark compares genomes with each other, so the strict tier
        # must mean the same thing in every row. Rung 3 now has a single,
        # annotation-based definition for every genome, so nothing has to be
        # pinned here; the experimentally gated rungs 4 to 6 are reported
        # separately where a reference table exists.
        case_kwargs = dict(config_kwargs)
        config = DiscoveryConfig(
            target=case.assembly,
            assembly=case.assembly,
            cds_fasta=case.fasta_path,
            cache_dir=cache_dir,
            **case_kwargs,
        )
        result = run_discovery(config, log=None, write_outputs=False)
        genes = result.genes

        for tier in TIERS:
            if tier not in tiers:
                continue
            called = set(genes.loc[genes[f"tier_{tier}"].fillna(False).astype(bool), "locus_tag"])
            metrics = _metrics(called, reference)
            extra = called - reference
            placed = {t for t in extra if any(t in members for members in other.values())}
            rows.append(
                {
                    "species": case.name,
                    "assembly": case.assembly,
                    "tier": tier,
                    "tier_definition": TIER_DEFINITIONS[tier],
                    **metrics,
                    "n_extra": len(extra),
                    "n_extra_in_other_curated_column": len(placed),
                    "n_extra_unplaced": len(extra) - len(placed),
                    "n_in_curated_housekeeping": len(called & housekeeping) if housekeeping else None,
                    "n_curated_housekeeping": len(housekeeping) or None,
                    "housekeeping_list_available": bool(case.housekeeping_column),
                    "in_sample_for_calibration": case.name.startswith("Salmonella enterica SL1344"),
                    "strict_definition": "annotation (ungated, comparable across genomes)",
                    "n_genes_in_genome": int(len(genes)),
                }
            )
            emit(
                f"[{case.name}] {tier:<11s} {metrics['n_called']:5d} genes, "
                f"recall {metrics['recall_of_reference']:.1%}, "
                f"concordance {metrics['concordance']:.1%}"
            )

        # Disagreement detail against the headline (standard) tier.
        called = set(genes.loc[genes["tier_standard"].fillna(False).astype(bool), "locus_tag"])
        status = []
        for row in genes.itertuples(index=False):
            tag = row.locus_tag
            in_ref, in_call = tag in reference, tag in called
            if in_ref == in_call:
                continue
            where = [name for name, members in other.items() if tag in members]
            status.append(
                {
                    "species": case.name,
                    "locus_tag": tag,
                    "gene": row.gene,
                    "product": row.product,
                    "disagreement": "called, not in reference" if in_call else "in reference, not called",
                    "tier": row.tier,
                    "evidence_score": row.evidence_score,
                    "in_permissive_tier": bool(row.tier_permissive),
                    "in_strict_tier": bool(row.tier_host_interaction),
                    "other_curated_columns": "; ".join(where),
                    "curator_housekeeping": tag in housekeeping,
                    "classification_reason": row.classification_reason,
                    "database_evidence": row.database_evidence,
                }
            )
        details[case.name] = pd.DataFrame(status)

    return pd.DataFrame(rows), details


def write_benchmark_workbook(
    path: str, summary: pd.DataFrame, details: Dict[str, pd.DataFrame]
) -> str:
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with pd.ExcelWriter(path, engine="xlsxwriter") as writer:
        summary.to_excel(writer, sheet_name="Concordance summary", index=False)
        for name, frame in details.items():
            if frame is None or frame.empty:
                continue
            sheet = "".join(c for c in name if c.isalnum() or c in " _-")[:31]
            frame.to_excel(writer, sheet_name=sheet, index=False)
    return path
