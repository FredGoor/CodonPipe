#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Automated identification of bacterial virulence genes for CodonPipe.

This subpackage answers a question that any codon-usage analysis of virulence
genes has to answer first: *which genes are the virulence genes?*  Until now
that set came from a database screen followed by manual curation, which is
slow, hard to report, and -- because the curation rules live in the curator's
head -- not reproducible by a reader.

The module replaces that step with a reproducible pipeline driven by a species
name, an NCBI taxonomy identifier or a genome assembly accession.  It harvests
curated virulence-factor evidence (VFDB, VICTORS and PATRIC_VF via BV-BRC),
UniProtKB curation, Gene Ontology closures, RefSeq annotation and genomic
context, scores each gene with a declared additive rule, and assigns a strictly nested
ladder of three tiers grounded in Falkow's molecular Koch's postulates:

``permissive``
    any virulence evidence at all, including annotation-only candidates --
    the deliberately over-inclusive set;

``standard``
    evidence anchored in database curation, or in operon context around a
    curated gene;

``strict``
    the host-interaction set: genes acting directly on the host with no
    evidence of a core cellular function, so that postulate 2 -- loss of
    virulence upon specific inactivation -- is interpretable rather than
    confounded by a general growth defect.

One run yields all three, so the sensitivity of a downstream result to the
choice of virulence-gene definition can be read off directly.

Every call carries its evidence trail, its per-postulate flags and the reason
it was included or excluded, and every run writes a provenance report and a
machine-readable manifest.

Quick start
-----------
>>> from codonpipe.virulence import DiscoveryConfig, run_discovery
>>> result = run_discovery(DiscoveryConfig(
...     target="Salmonella Typhimurium SL1344",
...     stringency="standard",
...     output_dir="virulence_out",
... ))
>>> len(result.permissive), len(result.standard), len(result.strict)  # doctest: +SKIP

The decision rules live in editable tab-separated tables under
``codonpipe/virulence/resources``; nothing is hidden in the code.
"""

from __future__ import annotations

from .classify import (
    CATEGORY_GENERAL,
    CATEGORY_HOST,
    STRINGENCY_PRESETS,
    TIER_DEFINITIONS,
    TIER_LABELS,
    TIERS,
    ScoringParams,
    classify_genes,
    default_params,
    params_for_stringency,
)
from .context import extend_by_context
from .outputs import CITATIONS
from .pipeline import (
    MODULE_VERSION,
    DiscoveryConfig,
    DiscoveryResult,
    run_discovery,
)
from .resolve import TargetInfo, GeneTable, resolve_target
from .species_loci import (
    KLEBORATE_SCORE_RULE,
    SPECIES_LOCI_TABLE,
    LocusCall,
    LocusDefinition,
    SpeciesLociResult,
    definitions_in_scope,
    detect_species_loci,
    kleborate_virulence_score,
    load_species_loci,
)
from .webcache import NetworkUnavailable, OfflineCacheMiss, WebCache

__version__ = MODULE_VERSION

__all__ = [
    "MODULE_VERSION",
    "__version__",
    "DiscoveryConfig",
    "DiscoveryResult",
    "run_discovery",
    "ScoringParams",
    "STRINGENCY_PRESETS",
    "TIERS",
    "TIER_LABELS",
    "TIER_DEFINITIONS",
    "default_params",
    "params_for_stringency",
    "classify_genes",
    "extend_by_context",
    "detect_species_loci",
    "load_species_loci",
    "definitions_in_scope",
    "kleborate_virulence_score",
    "KLEBORATE_SCORE_RULE",
    "SPECIES_LOCI_TABLE",
    "LocusDefinition",
    "LocusCall",
    "SpeciesLociResult",
    "CATEGORY_GENERAL",
    "CATEGORY_HOST",
    "CITATIONS",
    "TargetInfo",
    "GeneTable",
    "resolve_target",
    "WebCache",
    "OfflineCacheMiss",
    "NetworkUnavailable",
]
