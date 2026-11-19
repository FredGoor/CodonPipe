#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Evidence scoring and Falkow-anchored classification of virulence genes.

Why two categories
------------------
Codon-usage analyses of "virulence genes" are sensitive to what the gene set
contains.  A set assembled from curated virulence-factor databases mixes two
biologically different things: proteins that act on the host, and proteins
that a pathogen needs in order to function at all but that happen to have been
reported in an infection phenotype -- chaperones, proteases, central
metabolism, global regulators.  The two behave differently, and conflating
them weakens any statement about virulence-associated codon usage.

Stanley Falkow's molecular Koch's postulates (Falkow 1988; revisited in
Falkow 2004) give the distinction a principled basis:

1. the phenotype under investigation should be associated with pathogenic
   members of a genus, or with pathogenic strains of a species;
2. specific inactivation of the gene should cause a measurable loss of
   pathogenicity or virulence;
3. reversion or allelic replacement of the mutated gene should restore
   pathogenicity.

Postulate 2 is the discriminating one.  If disrupting a gene also impairs
growth in vitro, the loss of virulence cannot be attributed to the
host-interaction function of that gene: the attenuation is confounded by a
general growth defect.  Such a gene can still be *virulence-associated*, but
it cannot satisfy postulate 2 cleanly.  That is exactly the distinction this
module implements:

The three tiers
---------------
A single scoring pass assigns every gene to a **strictly nested ladder** of
three tiers, each of which is a *definition* rather than a tuning knob:

``permissive`` -- postulate 1, read inclusively
    Every gene with any virulence evidence at all: a curated virulence-factor
    record, a UniProtKB virulence annotation, operon context around one of
    those, or a supporting annotation match on its own.  Deliberately
    over-inclusive -- "all virulence genes, not just those that fit particular
    criteria" -- and therefore the right set for asking whether a result
    survives the most generous possible definition.

``standard`` -- postulate 1, anchored in curation
    The subset whose evidence is anchored in database curation, or in operon
    context around a curated gene, and that reaches the evidence threshold.
    Annotation-only candidates are excluded.  This is the default.

``strict`` -- postulate 2 interpretable (the host-interaction set)
    The subset of ``standard`` that additionally (i) carries positive evidence
    of acting on the host -- translocation into host cells, adhesion,
    invasion, toxicity, immune interference -- and (ii) shows no evidence of a
    core cellular function, so that a loss of virulence upon specific
    inactivation cannot be explained by a generic growth defect.

``strict`` is a subset of ``standard`` is a subset of ``permissive`` **by
construction**: each tier adds a predicate to the previous one, so editing a
weight or a rule table cannot break the nesting.  :func:`classify_genes`
asserts it regardless.

Because all three tiers come out of one scoring pass, a single run yields the
whole ladder, and the sensitivity of any downstream result to the choice of
virulence-gene definition can be read off directly.

Every gene carries explicit per-postulate flags and the evidence behind them,
so the basis of its assignment can be inspected rather than trusted.

A note on what cannot be computed
---------------------------------
Postulate 3 (complementation) is a wet-lab result and is not recorded in a
machine-readable form by any of the sources used here.  It is reported as
``reported_in_source_literature`` when the underlying database entry is
literature-curated, and ``not_assessed`` otherwise; it never drives the
classification.  Likewise, in vitro dispensability is inferred from functional
annotation, not measured: where an experimental essentiality or fitness table
is supplied it overrides the annotation-based inference.

References
----------
Falkow S (1988) Molecular Koch's postulates applied to microbial
pathogenicity. *Rev Infect Dis* 10(Suppl 2):S274-S276.
Falkow S (2004) Molecular Koch's postulates applied to bacterial
pathogenicity -- a personal recollection 15 years later.
*Nat Rev Microbiol* 2:67-72.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Sequence, Set, Tuple

import pandas as pd

from .context import extend_by_context
from .resolve import GeneTable

__all__ = [
    "ScoringParams",
    "STRINGENCY_PRESETS",
    "TIERS",
    "TIER_LABELS",
    "LADDER",
    "LADDER_LABELS",
    "PUBLISHED_TIERS",
    "GATES",
    "TIER_NUMBERS",
    "LADDER_DEFINITIONS",
    "IN_VIVO_SPECIFIC_LABEL",
    "TIER_ALIASES",
    "canonical_tier",
    "TIER_DEFINITIONS",
    "tier_definitions",
    "params_for_stringency",
    "default_params",
    "classify_genes",
    "CATEGORY_GENERAL",
    "CATEGORY_HOST",
]

#: The full ladder, from most inclusive to most restrictive. Every rung is
#: computed in every run. Containment holds along two chains, asserted at run
#: time (see LADDER_BRANCHES for why it is two and not one):
#:
#:   host_interaction < standard < permissive
#:   in_vitro_neutral < in_vivo_required < in_vivo_contributing < standard
#:
#: Rungs 4 to 6 need experimental mutant-screen data and are empty for a genome
#: that has none; rung 6 additionally needs *graded* in vitro fitness data, not
#: just essentiality, and collapses onto rung 5 until such data are supplied
#: (see reference.py, which supplies it for Salmonella from three screens).
LADDER = (
    "permissive",
    "standard",
    "host_interaction",
    "in_vivo_contributing",
    "in_vivo_required",
    "in_vitro_neutral",
)

#: The ladder is a nested chain, with the experimental rungs refining rung 2.
#:
#: An earlier revision let rungs 4-6 stand outside the annotation rungs, on the
#: reasoning that a measured in vivo phenotype outranks an annotation-based
#: inference. True in principle, unusable in practice: published screens
#: disagree by fifteen-fold about what counts as attenuated, so admitting every
#: screen hit regardless of annotation gave 1,242 rung-4 genes in UPEC UTI89 --
#: a quarter of that genome -- against 337 in Salmonella, inverting the known
#: biology and measuring library depth rather than the pathogen. Requiring
#: curation-anchored membership first makes the experimental rungs a refinement
#: of a set computed identically for every genome, which is what makes them
#: interpretable.
#:
#: The gate is rung 2, not rung 3. Rung 3 adds a host-DIRECTEDNESS test on top
#: of curation, and that test legitimately excludes pleiotropic regulators:
#: phoP and phoQ carry only the generic UniProt Virulence keyword on the host
#: axis and also regulate magnesium transport and other core functions, so they
#: stop at rung 2. Gating the experimental rungs on rung 3 would therefore have
#: dropped phoP, phoQ, ompR, ssrA and ssrB from the required set, which no
#: reading of the Salmonella literature supports. Rung 2 is also the faithful
#: statement of the requirement: the gene must have been found by the database
#: and annotation analysis AND by the in vivo screen.
#:
#: Rungs 4-6 are each nested in the one before. Genes attenuated in vivo but
#: outside rung 2 are recorded in the `in_vivo_outside_rung2` column rather
#: than discarded.
LADDER_BRANCHES = {
    "annotation": ("permissive", "standard", "host_interaction"),
    "experimental": (
        "standard", "in_vivo_contributing", "in_vivo_required", "in_vitro_neutral",
    ),
}

#: Column names written to the cluster workbook, numbered so that the ladder
#: order survives alphabetical sorting in a spreadsheet.
#: The six rungs are the computation; the published scheme has FOUR tiers.
#: Rungs 2 and 4 are gates, not tiers: "standard" is the curation-anchored gate
#: that every experimental rung is computed inside, and "in vivo contributing"
#: is the 0.0-log2 cut that the required tier sharpens to 1.0 log2. They are
#: still written, under gate names, because the sets are needed to audit the
#: tiers - but they are not numbered, so a reader cannot mistake them for
#: tiers of the published ladder.
PUBLISHED_TIERS = (
    "permissive", "host_interaction", "in_vivo_required", "in_vitro_neutral",
)
GATES = ("standard", "in_vivo_contributing")

LADDER_LABELS = {
    "permissive": "Virulence tier 1",
    "standard": "Gate - curation anchored",
    "host_interaction": "Virulence tier 2",
    "in_vivo_contributing": "Gate - in vivo contributing",
    "in_vivo_required": "Virulence tier 3",
    "in_vitro_neutral": "Virulence tier 4",
}

#: The published tier number of each rung, for reports. Gates have none.
TIER_NUMBERS = {
    "permissive": 1, "host_interaction": 2,
    "in_vivo_required": 3, "in_vitro_neutral": 4,
}

#: What each rung adds to the one below it.
LADDER_DEFINITIONS = {
    "permissive": (
        "postulate 1 read inclusively: any virulence evidence, including a "
        "product-name annotation match standing on its own"
    ),
    "standard": (
        "+ the evidence must be anchored in curation: a curated virulence "
        "record, a UniProtKB virulence keyword, a published acquired-locus "
        "definition, or operon context around one of those. Annotation-only "
        "candidates drop out"
    ),
    "host_interaction": (
        "+ postulate 2 interpretable: positive host-directed evidence, and no "
        "core cellular role that would confound a loss-of-virulence phenotype. "
        "Experimentally essential genes are excluded here; experimentally "
        "dispensable ones are admitted even against core-process annotation"
    ),
    "in_vivo_contributing": (
        "postulate 2 satisfied, at contributing strength: the mutant shows a "
        "statistically significant fitness defect in an infection model with "
        "at least the contributing effect size. Experimental, and NOT filtered "
        "by rung 3 -- a measured phenotype outranks an annotation-based "
        "inference -- but still excluding genes that cannot be inactivated at "
        "all, since those yield no interpretable phenotype"
    ),
    "in_vivo_required": (
        "+ at required strength: the defect reaches the larger effect size, "
        "the difference between a gene that contributes to fitness in the host "
        "and one the pathogen cannot do without"
    ),
    "in_vitro_neutral": (
        "+ the in vivo defect is SPECIFIC to the host interaction: the gene "
        "also passes the rung-3 test (host-directed evidence, no confounding "
        "core cellular role), so its attenuation is unlikely to be a general "
        "fitness effect. This is the only rung that takes rung 3 as well as "
        "rung 2; rungs 4 and 5 deliberately do not, because gating them "
        "would delete the O-antigen block and the global regulators from the "
        "required set. The genes rung 6 removes from rung 5 stay in rung 5 - "
        "they are not demoted, they are marked as not demonstrably specific. "
        "AND no measured in vitro fitness cost. READ THAT SECOND CLAUSE: it is "
        "only as strong as the in vitro column of the reference table. Where "
        "that column is BINARY essentiality, 'no measured cost' means no more "
        "than 'a viable mutant was obtained', which does not exclude a growth "
        "or stress-survival defect short of lethality; the clause then "
        "contributes almost nothing and the rung is carried by the "
        "host-directedness clause, and the run summary says so. It becomes "
        "load-bearing only when the table supplies GRADED in vitro fitness, "
        "i.e. an in_vitro_requirement of 'fitness_cost' for "
        "viable-but-impaired mutants. The bundled Salmonella table does "
        "supply that, from three screens -- one transferred from S. Typhi and "
        "two measured natively in Typhimurium -- so for SL1344 the clause "
        "fires and removes 4 genes that pass the host-directedness test. "
        "It removed 7 before the D23580 screen was added; phoP, phoQ and "
        "dsbA were returned to the rung when a native screen measured their "
        "mutants at or above neutral abundance in rich medium and so "
        "positively contradicted the transferred Typhi call. rpoS remains "
        "the case worth knowing: viable, non-essential, plainly impaired in "
        "stationary phase, and passed by the in vitro clause -- it is "
        "excluded by the host-directedness clause instead (host axis 0.00, "
        "core axis 1.10)"
    ),
}

#: The legacy three-tier view, retained so that saved sessions, scripts and
#: the concordance benchmark keep working. ``strict`` aliases whichever rung
#: one definition; nothing selects between variants.
#:
#: The third tier has two names in common use and they mean the same thing.
#: "Strict" names its position on the ladder; "host-interaction" names what it
#: contains. Internally the key is always ``"strict"``; every user-facing label
#: says "host-interaction" with "(strict tier)" attached, so that the panel,
#: the workbook column and the report cannot be read as three different sets.
TIERS = ("permissive", "standard", "host_interaction")

#: Accepted spellings for the third tier on the command line and in config.
TIER_ALIASES = {
    "strict": "strict",
    "host-interaction": "strict",
    "host_interaction": "strict",
    "hostinteraction": "strict",
    "host interaction": "strict",
}


def canonical_tier(name: str) -> str:
    """Map any accepted spelling of a tier name onto its internal key."""

    key = str(name or "").strip().lower()
    return TIER_ALIASES.get(key, key)

#: A set defined purely by experiment, standing outside the tier ladder.
#:
#: It holds every gene whose mutant is attenuated in an animal infection model
#: and which is *not* required for growth in vitro. No annotation enters it:
#: it is the two experimental axes of :mod:`codonpipe.virulence.reference` and
#: nothing else.
#:
#: It is deliberately NOT called a host-interaction set. "Required in vivo and
#: dispensable in vitro" is a weaker statement than "acts on the host": an
#: auxotroph is attenuated because the host is nutrient-poor, not because its
#: product engages host tissue, and in rich medium it carries no in vitro
#: fitness cost at all. In S. Typhimurium this set contains aroA, aroC, aroD,
#: carA, carB, guaB, purA, pyrC and pyrD -- the classic attenuating
#: auxotrophs. Reporting it as host-interaction would smuggle them in.
#:
#: The gated host-interaction tier is a subset of it (both require a measured
#: in vivo phenotype, and the host-interaction tier adds the annotation-based
#: tests); that containment is asserted at run time.
#: A non-ladder column, deliberately UNGATED: every gene the screens called
#: attenuated in a host and that carries no measured in vitro cost, whatever
#: its annotation says and without the significance or effect-size filters the
#: ladder applies. It is larger than rung 4 for that reason -- 411 against 138
#: in SL1344 -- and the old name "In vivo required, in vitro dispensable"
#: invited that difference to be read as an inconsistency. The name now states
#: the two things that distinguish it from rung 5: no annotation requirement,
#: and no effect-size requirement.
IN_VIVO_SPECIFIC_LABEL = (
    "Attenuated in vivo, no in vitro cost (ungated, any effect size)"
)

#: Column names used in the output cluster workbook.
TIER_LABELS = {
    "permissive": "Virulence permissive",
    "standard": "Virulence standard",
    "strict": "Host-interaction (strict tier)",
}

#: One-line definitions, reproduced into the provenance report.
TIER_DEFINITIONS = {
    "permissive": (
        "postulate 1 read inclusively: any virulence evidence, including "
        "annotation-only candidates"
    ),
    "standard": (
        "postulate 1 anchored in curation: database-curated evidence, or operon "
        "context around a curated gene, above the evidence threshold"
    ),
    "strict": (
        "postulate 2 interpretable: the standard set restricted to genes with "
        "positive host-directed evidence and no core-cellular-function evidence"
    ),
}



def tier_definitions(strict_requires: str = "") -> Dict[str, str]:
    """The tier definitions. The argument is accepted and ignored.

    It used to specialise the strict line to whichever of four rules was in
    force; rung 3 now has one definition, so there is nothing to specialise.
    The parameter is kept so existing callers do not break.
    """
    return dict(TIER_DEFINITIONS)

# Retained for continuity with CodonPipe 1.2.0 outputs, where the two
# categories were reported instead of the tier ladder. CATEGORY_GENERAL is
# the standard tier and CATEGORY_HOST is the strict tier.
CATEGORY_GENERAL = "Virulence (general)"
CATEGORY_HOST = "Host-interaction"

_UNINFORMATIVE = {"uninformative", "negative"}


@dataclass
class ScoringParams:
    """Every weight and threshold used by the classifier, in one place.

    The whole object is serialised into the run manifest, so a published
    result states its own decision rule.  Score weights are additive; the
    default values and the three stringency thresholds were calibrated on a
    single genome -- the manually curated *Salmonella* Typhimurium SL1344
    virulence-gene set of the accompanying study (see CITATION.cff) -- by
    sweeping the score
    threshold and reading off the concordance curve.  Whether that calibration
    transfers to other pathogens is an empirical question, answered by
    :mod:`codonpipe.virulence.benchmark`, which recomputes concordance against
    independently curated sets for several species.  Treat the defaults as a
    documented starting point rather than a tuned optimum, and re-examine them
    for a new clade.
    """

    # --- virulence-factor database evidence -------------------------------
    w_vfdb_categorised: float = 2.0
    w_vfdb_uncategorised: float = 1.0
    w_other_db_categorised: float = 1.2
    w_other_db_uncategorised: float = 0.4
    w_literature_bonus: float = 0.8
    w_negative_record: float = -1.5

    # --- UniProtKB curation -----------------------------------------------
    w_uniprot_virulence_strain: float = 2.0
    w_uniprot_virulence_transfer: float = 1.4
    w_host_axis_scale: float = 0.4
    w_core_axis_scale: float = -0.8
    cap_host_axis: float = 2.0
    cap_core_axis: float = 3.0
    species_transfer_scale: float = 0.7

    # --- controlled vocabulary --------------------------------------------
    w_vocabulary_scale: float = 0.6
    w_core_vocabulary: float = -1.2

    # --- species-specific acquired virulence loci (E6) --------------------
    #: Multiplier applied to the per-locus ``weight`` column of
    #: ``species_virulence_loci.tsv`` for a *marker* gene, i.e. one whose own
    #: gene symbol or RefSeq product name matched the published locus
    #: definition. The shipped weights are 2.0, so a marker gene clears the
    #: standard-tier threshold on the strength of the locus call alone -- which
    #: is the intent: the locus has been established by several co-located
    #: annotation matches, and the published locus definition is the curation.
    w_species_locus_scale: float = 1.0
    #: Flat score for a gene recruited into an already-accepted locus block by
    #: the operon/island machinery rather than by its own annotation. Set below
    #: the shipped marker weight, because the evidence for such a gene is its
    #: membership of the locus rather than anything said about the gene itself.
    w_species_locus_context: float = 1.5
    #: Host-interaction axis credit for a member of a locus whose ``axis_hint``
    #: is ``host_interaction``. Acquired siderophore, genotoxin and capsule
    #: up-regulator loci act on the host, so they are allowed to satisfy the
    #: strict tier's host-directedness test; the siderophore and capsule loci
    #: are nevertheless classified ``dual_use``, so a run with
    #: ``strict_requires_dedicated=True`` still excludes them.
    w_species_locus_host_axis: float = 1.0

    # --- genomic context ---------------------------------------------------
    use_context: bool = True
    w_context_operon: float = 0.9
    w_context_island: float = 0.6
    context_seed_threshold: float = 1.5
    max_intergenic_gap: int = 50
    require_same_strand: bool = True
    island_bridge_genes: int = 2
    min_island_seeds: int = 3

    # --- database alignment gate ------------------------------------------
    min_identity: float = 80.0
    min_coverage: float = 60.0

    # --- tier thresholds ---------------------------------------------------
    #: Minimum score for the permissive tier. Low by design: the permissive
    #: tier is defined by the *presence* of evidence, not by its strength.
    threshold_permissive: float = 0.5
    #: Minimum score for the standard tier, which additionally requires the
    #: evidence to be anchored in curation.
    threshold_standard: float = 1.5
    #: Host-directedness and core-process tolerances for the strict tier.
    min_host_axis: float = 0.5
    max_core_axis: float = 1.0
    require_host_axis_exceeds_core: bool = True

    #: When true, the strict tier additionally requires *dedicated* host-attack
    #: machinery, excluding genes whose only host-directed evidence comes from
    #: dual-use structures (fimbriae and curli, flagella and chemotaxis,
    #: siderophores, O-antigen, capsule). These serve biofilm, environmental
    #: persistence and abiotic-surface attachment as well as host attachment,
    #: so whether they satisfy Falkow's second postulate is a biological
    #: judgement rather than a computable fact. Reported both ways by the
    #: benchmark; left permissive by default so that nothing is silently
    #: discarded.
    strict_requires_dedicated: bool = False

    #: Rung 3 has a single definition and no selector. The removed
    #: ``strict_requires`` parameter chose between four readings of the
    #: tier and defaulted to ``"auto"``, which switched between them on
    #: reference-table coverage -- so the tier meant different things on
    #: different genomes without the operator choosing. Experimental
    #: evidence now has rungs 4 to 6 of its own.


    #: Effect-size and significance cuts that separate rung 4 (contributing)
    #: from rung 5 (required). Both are conventions, not facts derivable from
    #: the data, and both are reported in the provenance report of every run so
    #: that a result can be re-read under a different cut.
    #:
    #: ``in_vivo_effect_log2`` in a reference table is the magnitude of the
    #: mutant's fitness defect in log2 units, positive meaning depleted in the
    #: host. 1.0 is a two-fold depletion, the usual line below which a
    #: transposon-sequencing defect is not called. 3.32 is log2(10), the
    #: classical ten-fold competitive-index threshold for calling a mutant
    #: attenuated, and is what separates "contributes to fitness in the host"
    #: from "the pathogen cannot do without it".
    #:
    #: Where a table reports an attenuation call but no effect size, the gene
    #: enters rung 4 on the call alone and cannot reach rung 5; that is
    #: recorded per gene rather than guessed at.
    #: How the two cut-offs below are read. "absolute" treats them as log2
    #: depletion values, which is the right choice when effect sizes come from
    #: competitive indices. "percentile" treats them as percentiles of the
    #: measured in vivo effect sizes in THIS organism, which is the right
    #: choice when comparing organisms whose screens are on different scales:
    #: a pooled transposon library is bottlenecked, so its log fold changes
    #: are systematically larger than a competitive index of the same
    #: biological severity, and one absolute line cannot serve both.
    #:
    #: The defaults are the thresholds the bundled Salmonella compilation uses
    #: on its own data: `COMPETITIVE_INDEX_ATTENUATED_AT` is 1.0, a two-fold
    #: under-representation of the mutant. A ten-fold line was tried and is
    #: wrong for this biology -- the SPI-1 competitive indices in that
    #: compilation run from 0.4 to 2.1 log2, because its screens are systemic
    #: infection models in which SPI-1 mutants are only mildly attenuated, so
    #: a ten-fold bar removed the entire SPI-1 apparatus from rung 5 while
    #: keeping SPI-2. Rung 4 takes any positive measured defect.
    #: Whether a gene needs host-axis evidence from a source that says
    #: something SPECIFIC about host interaction, rather than only the generic
    #: UniProt "Virulence" keyword or a generic virulence-factor category.
    #:
    #: Both settings are self-consistent; the previous behaviour was not,
    #: because it vetoed generic-only evidence when a virulence-factor
    #: category was present and admitted it when none was. Requiring specific
    #: evidence is the stricter reading but costs 10.5 points of rung-3 recall
    #: in S. aureus NCTC 8325, where 31 curated genes carry the generic keyword
    #: and nothing else, so it is off by default and the generic keyword is
    #: accepted uniformly. The `host_evidence` column always names the source,
    #: so a reader can see which genes rest on it.
    require_specific_host_evidence: bool = False
    effect_scale: str = "absolute"
    effect_contributing_log2: float = 0.0
    effect_required_log2: float = 1.0
    #: Maximum significance value accepted from a screen. Genes whose source
    #: reports no significance value are admitted on the attenuation call.
    max_in_vivo_significance: float = 0.05

    #: Which infection models may supply an in vivo call for rungs 4 and 5.
    #: Empty (the default) accepts every model and simply records which one
    #: each call came from; a gene's model is already resolved by precedence
    #: when the reference table is built, so an animal result outranks a
    #: cell-infection one for the same gene. Set this to
    #: ``("animal_infection",)`` to restrict the ladder to whole-animal
    #: evidence, at the cost of every gene known only from cell infection,
    #: organoid or ex vivo work.
    in_vivo_accepted_models: tuple = ()

    #: Which tier the run reports as its primary set. Does not change the
    #: scoring -- all three tiers are always computed.
    stringency: str = "standard"

    def to_dict(self) -> dict:
        return asdict(self)


#: Retained so that saved sessions and scripts naming a stringency keep
#: working.  Tiers are no longer selected by swapping weight sets -- all three
#: come out of one scoring pass -- so these presets only record which tier a
#: run reports as its primary set, plus the historical tolerances.
STRINGENCY_PRESETS: Dict[str, dict] = {
    "permissive": {"stringency": "permissive", "min_host_axis": 0.3, "max_core_axis": 1.5},
    "standard": {"stringency": "standard"},
    "strict": {"stringency": "strict"},
}


def default_params(**overrides) -> ScoringParams:
    """Build the scoring parameters, with optional explicit overrides.

    There is one weight set.  The tier ladder is produced by applying three
    nested predicates to the resulting score, so changing a weight shifts all
    three tiers coherently instead of redefining them independently.
    """
    params = ScoringParams()
    for name, value in overrides.items():
        if value is None:
            continue
        if not hasattr(params, name):
            raise ValueError(f"Unknown scoring parameter {name!r}.")
        setattr(params, name, value)
    return params


def params_for_stringency(stringency: str = "standard", **overrides) -> ScoringParams:
    """Scoring parameters whose *primary reported tier* is ``stringency``.

    All three tiers are computed in every run; this only records which one the
    run treats as its headline set, and applies the historical host/core
    tolerances associated with that name.
    """
    key = canonical_tier(stringency or "standard")
    if key not in STRINGENCY_PRESETS:
        raise ValueError(
            f"Unknown tier {stringency!r}; choose one of {', '.join(TIERS)}."
        )
    merged = dict(STRINGENCY_PRESETS[key])
    merged.update({k: v for k, v in overrides.items() if v is not None})
    return default_params(**merged)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _split_categories(text) -> List[str]:
    """VF classification fields mix ``;`` and ``,`` separators."""
    out: List[str] = []
    for chunk in str(text or "").replace(";", ",").split(","):
        token = chunk.strip()
        if token and token.lower() not in ("nan", "none"):
            out.append(token)
    return out


def _axes_for_categories(categories: Sequence[str], axis_map: Dict[str, str]) -> Set[str]:
    return {axis_map[c.lower()] for c in categories if c.lower() in axis_map}


class _Accumulator:
    """Per-gene additive score with a human-readable component trail."""

    def __init__(self) -> None:
        self.score: Dict[str, float] = {}
        self.components: Dict[str, List[str]] = {}

    def add(self, tag: str, value: float, reason: str) -> None:
        if not tag or not value:
            return
        self.score[tag] = self.score.get(tag, 0.0) + float(value)
        self.components.setdefault(tag, []).append(f"{reason} ({value:+.2f})")


# ---------------------------------------------------------------------------
# main entry point
# ---------------------------------------------------------------------------
def classify_genes(
    gene_table: GeneTable,
    vf_frame: pd.DataFrame,
    uniprot_frame: pd.DataFrame,
    vocabulary_hits: pd.DataFrame,
    core_vocabulary_hits: pd.DataFrame,
    go_closures: Dict[str, Set[str]],
    keyword_axes: pd.DataFrame,
    vf_category_axes: Dict[str, str],
    params: ScoringParams,
    essentiality: Optional[pd.DataFrame] = None,
    species_locus_hits: Optional[pd.DataFrame] = None,
) -> Tuple[pd.DataFrame, Dict[str, object]]:
    """Score every gene and assign the two virulence categories.

    Parameters
    ----------
    gene_table:
        The gene universe (see :func:`codonpipe.virulence.resolve.build_crosswalk`).
    vf_frame:
        Mapped virulence-factor records from BV-BRC (VFDB / VICTORS / PATRIC_VF).
    uniprot_frame:
        Mapped UniProtKB entries.
    vocabulary_hits, core_vocabulary_hits:
        Controlled-vocabulary matches for virulence and core-process lexicons.
    go_closures:
        ``{"host_interaction": {...GO ids...}, "core_process": {...}}``.
    keyword_axes:
        UniProt keyword-to-axis table.
    vf_category_axes:
        Virulence-factor category-to-axis map.
    essentiality:
        Optional experimental table with columns ``locus_tag`` and
        ``in_vitro_requirement`` (``essential``, ``fitness_cost`` or
        ``dispensable``).  When supplied it overrides the annotation-based
        in vitro inference for the genes it covers.
    species_locus_hits:
        Optional per-gene membership table for the species-specific acquired
        virulence loci (evidence layer E6), as produced by
        :func:`codonpipe.virulence.species_loci.detect_species_loci`.  Columns
        used: ``locus_tag``, ``locus``, ``role`` (``marker`` or
        ``locus_context``), ``machinery``, ``axis_hint``, ``weight``.  Passing
        ``None``, or an empty frame, leaves every score and every tier
        unchanged, which is what happens for a genome outside the configured
        taxonomic scopes.

    Returns
    -------
    (frame, summary)
        ``frame`` has one row per gene of the universe; ``summary`` carries
        counts and the notes generated during classification.
    """
    tags = gene_table.frame["locus_tag"].astype(str).tolist()
    tag_set = set(tags)
    notes: List[str] = []

    acc = _Accumulator()
    host_axis: Dict[str, float] = {}
    core_axis: Dict[str, float] = {}
    # UniProt-derived portion of the two axes, tracked separately because it is
    # the only portion that contributes to the additive evidence score.
    uniprot_host: Dict[str, float] = {}
    uniprot_core: Dict[str, float] = {}
    host_ev: Dict[str, List[str]] = {}
    core_ev: Dict[str, List[str]] = {}
    db_ev: Dict[str, List[str]] = {}
    subcats: Dict[str, Set[str]] = {}
    cluster_subcats: Dict[str, Set[str]] = {}
    machinery_classes: Dict[str, Set[str]] = {}
    vf_axes: Dict[str, Set[str]] = {}
    #: Tags that received host-axis credit from a SOURCE THAT SAYS SOMETHING
    #: SPECIFIC about host interaction -- a host-interaction VF category, a
    #: host GO term, a named vocabulary subcategory, an acquired locus, or a
    #: UniProt keyword other than the generic KW-0843 "Virulence".
    #:
    #: Without this, host-directedness depended on whether BV-BRC happened to
    #: return a record: the `general_only` veto is computed from VF categories,
    #: so a gene whose only host evidence was the generic Virulence keyword
    #: was vetoed when it had a VF category and ADMITTED when it had none.
    #: 16 of 245 rung-3 genes in SL1344 entered that way.
    host_specific: Set[str] = set()
    has_literature: Set[str] = set()
    has_database: Set[str] = set()
    has_uniprot_virulence: Set[str] = set()
    negative_record: Set[str] = set()
    # E6: species-specific acquired virulence loci
    species_locus_member: Dict[str, Set[str]] = {}
    species_locus_role: Dict[str, str] = {}

    def bump(store: Dict[str, float], tag: str, value: float) -> None:
        store[tag] = store.get(tag, 0.0) + float(value)

    # ---------------------------------------------- E1 virulence databases
    if vf_frame is not None and not vf_frame.empty:
        # BV-BRC can return the same specialty-gene record more than once for
        # one feature, and the loop below credits each row, so the same
        # database accession was counted twice: 133 exactly duplicated
        # (locus_tag, source, source_id) rows in SL1344, 249 in S. flexneri,
        # 125 in Sakai, 102 in UTI89. pefC showed "VFDB record with
        # virulence-factor category (+2.00)" twice for accession VFG000434.
        # One accession is one fact.
        before = len(vf_frame)
        if {"locus_tag", "source", "source_id"} <= set(vf_frame.columns):
            vf_frame = vf_frame.drop_duplicates(["locus_tag", "source", "source_id"])
            if len(vf_frame) < before:
                notes.append(
                    f"{before - len(vf_frame)} duplicated virulence-factor record(s) "
                    "were collapsed; each database accession is credited once per gene."
                )
        for row in vf_frame.itertuples(index=False):
            tag = str(getattr(row, "locus_tag", "") or "")
            if tag not in tag_set:
                continue
            source = str(getattr(row, "source", "") or "unknown")
            categories = _split_categories(getattr(row, "classification", ""))
            axes = _axes_for_categories(categories, vf_category_axes)
            informative = bool(axes - _UNINFORMATIVE)
            is_vfdb = source.upper() == "VFDB"
            if informative:
                weight = params.w_vfdb_categorised if is_vfdb else params.w_other_db_categorised
                label = f"{source} record with virulence-factor category"
            else:
                weight = params.w_vfdb_uncategorised if is_vfdb else params.w_other_db_uncategorised
                label = f"{source} record without functional category"
            acc.add(tag, weight, label)
            has_database.add(tag)
            source_id = str(getattr(row, "source_id", "") or "")
            identity = getattr(row, "identity", None)
            detail = f"{source}:{source_id}" if source_id else source
            if categories:
                detail += f" [{'; '.join(categories)}]"
            try:
                if identity is not None and identity == identity:
                    detail += f" {float(identity):.0f}% id"
            except Exception:
                pass
            db_ev.setdefault(tag, []).append(detail)
            if bool(getattr(row, "is_literature", False)):
                acc.add(tag, params.w_literature_bonus, f"{source} literature curation")
                has_literature.add(tag)
            if "negative" in axes:
                acc.add(tag, params.w_negative_record, f"{source} states no virulence evidence")
                negative_record.add(tag)
            vf_axes.setdefault(tag, set()).update(axes)
            if "host_interaction" in axes:
                bump(host_axis, tag, 1.0)
                host_specific.add(tag)
                host_ev.setdefault(tag, []).append(
                    f"virulence-factor category: {'; '.join(categories)}"
                )
            if "virulence_general" in axes:
                core_ev.setdefault(tag, []).append(
                    f"virulence-factor category engaged in vitro: {'; '.join(categories)}"
                )

    # ------------------------------------------------------ E2 UniProtKB
    kw_by_axis: Dict[str, Dict[str, float]] = {}
    for row in keyword_axes.itertuples(index=False):
        kw_by_axis.setdefault(str(row.axis), {})[str(row.keyword_id).strip()] = float(row.weight)
    host_terms = set(go_closures.get("host_interaction") or ())
    core_terms = set(go_closures.get("core_process") or ())

    if uniprot_frame is not None and not uniprot_frame.empty:
        # UniProtKB evidence is aggregated PER GENE before it is scored, not
        # per entry.
        #
        # A gene symbol routinely maps to several UniProt entries -- one for
        # the strain plus one for each reviewed entry of the parent species --
        # and phoP in SL1344 maps to seven. Scoring row by row credited the
        # same curatorial fact once per entry, so the host and core axes grew
        # with how heavily the parent species had been sequenced rather than
        # with the gene's biology: phoP reached a host axis of 5.20 from a
        # single keyword (1.0 for the strain entry plus 6 x 0.7 for the
        # transferred ones), and genes with six or more entries averaged a
        # core axis of 8.01 against 0.22 for single-entry genes, which pushed
        # well-curated genes out of the host-interaction rung for no
        # biological reason.
        #
        # Each distinct keyword and GO term is therefore counted once, at the
        # better of the scopes it was seen in: full weight if any strain-scope
        # entry carries it, transfer weight if it is only in the parent
        # species. No evidence is lost -- the union of terms is used -- and the
        # accession list is still recorded in full.
        strain_kw: Dict[str, Set[str]] = {}
        transfer_kw: Dict[str, Set[str]] = {}
        strain_go: Dict[str, Set[str]] = {}
        transfer_go: Dict[str, Set[str]] = {}
        accessions: Dict[str, List[str]] = {}
        for row in uniprot_frame.itertuples(index=False):
            tag = str(getattr(row, "locus_tag", "") or "")
            if tag not in tag_set:
                continue
            is_strain = str(getattr(row, "scope", "strain")) == "strain"
            keyword_ids = {
                k.strip() for k in str(getattr(row, "keyword_ids", "")).split(";") if k.strip()
            }
            go_ids = {
                token.strip()
                for token in str(getattr(row, "go_ids", "")).replace(";", " ").split()
                if token.strip().startswith("GO:")
            }
            (strain_kw if is_strain else transfer_kw).setdefault(tag, set()).update(keyword_ids)
            (strain_go if is_strain else transfer_go).setdefault(tag, set()).update(go_ids)
            accession = str(getattr(row, "uniprot", "") or "")
            if accession:
                accessions.setdefault(tag, []).append(
                    accession + ("" if is_strain else " (transferred)")
                )

        for tag in set(strain_kw) | set(transfer_kw) | set(strain_go) | set(transfer_go):
            kw_strain = strain_kw.get(tag, set())
            # A term already seen at strain scope is not also counted as
            # transferred.
            kw_transfer = transfer_kw.get(tag, set()) - kw_strain
            go_strain = strain_go.get(tag, set())
            go_transfer = transfer_go.get(tag, set()) - go_strain

            if "KW-0843" in kw_strain or "KW-0843" in kw_transfer:
                at_strain = "KW-0843" in kw_strain
                weight = (
                    params.w_uniprot_virulence_strain
                    if at_strain
                    else params.w_uniprot_virulence_transfer
                )
                scope_label = "strain" if at_strain else "species_reviewed"
                acc.add(tag, weight, f"UniProtKB Virulence keyword ({scope_label})")
                has_uniprot_virulence.add(tag)
                shown = accessions.get(tag, [])
                db_ev.setdefault(tag, []).append(
                    "UniProtKB:"
                    + ", ".join(shown[:4])
                    + (f" and {len(shown) - 4} more" if len(shown) > 4 else "")
                    + f" KW-0843 ({scope_label})"
                )

            for axis, table in kw_by_axis.items():
                matched = [
                    (k, w, 1.0) for k, w in table.items() if k in kw_strain
                ] + [
                    (k, w, params.species_transfer_scale)
                    for k, w in table.items()
                    if k in kw_transfer
                ]
                if not matched:
                    continue
                total = sum(w * s for _, w, s in matched)
                is_host = axis == "host_interaction"
                if is_host and any(k != "KW-0843" for k, _, _ in matched):
                    host_specific.add(tag)
                bump(host_axis if is_host else core_axis, tag, total)
                bump(uniprot_host if is_host else uniprot_core, tag, total)
                (host_ev if is_host else core_ev).setdefault(tag, []).append(
                    "UniProt keywords " + ", ".join(sorted({k for k, _, _ in matched}))
                )

            for terms, label, target, ev in (
                (host_terms, "host-interaction", host_axis, host_ev),
                (core_terms, "core-process", core_axis, core_ev),
            ):
                n_strain = len(go_strain & terms)
                n_transfer = len((go_transfer & terms) - (go_strain & terms))
                if not (n_strain or n_transfer):
                    continue
                # The cap applies to the DISTINCT term count, as before, so a
                # richly annotated protein cannot dominate.
                value = (
                    min(n_strain, 3) * 0.5
                    + min(n_transfer, 3) * 0.5 * params.species_transfer_scale
                )
                bump(target, tag, value)
                bump(uniprot_host if target is host_axis else uniprot_core, tag, value)
                if target is host_axis:
                    host_specific.add(tag)
                ev.setdefault(tag, []).append(
                    f"{n_strain + n_transfer} {label} GO term(s)"
                )

    # The UniProt axes feed the evidence score as well as the category gates:
    # host-directed curation raises the score, core-process curation lowers it.
    # The cap applies to THIS score contribution only. The host_axis and
    # core_axis values that the rung-3 gates read are deliberately uncapped
    # sums, so their magnitude still tracks how deeply a genome is curated
    # (observed maxima 4.6 in KPPR1 to 9.7 in S. aureus). Capping the gate
    # inputs is not a safe fix: the two caps differ (2.0 and 3.0), so clamping
    # would invert require_host_axis_exceeds_core for genes whose uncapped
    # host axis exceeds their core axis. The limitation is declared instead,
    # and rung 3 retains 70-93% of rung 2 depending on the organism.
    for tag, value in uniprot_host.items():
        acc.add(
            tag,
            min(value, params.cap_host_axis) * params.w_host_axis_scale,
            "UniProt host-interaction annotation",
        )
    for tag, value in uniprot_core.items():
        acc.add(
            tag,
            min(value, params.cap_core_axis) * params.w_core_axis_scale,
            "UniProt core-process annotation",
        )

    # --------------------------------------------- E4 controlled vocabulary
    if vocabulary_hits is not None and not vocabulary_hits.empty:
        for row in vocabulary_hits.itertuples(index=False):
            tag = str(row.locus_tag)
            if tag not in tag_set:
                continue
            acc.add(
                tag,
                float(row.weight) * params.w_vocabulary_scale,
                f"annotation matches '{row.subcategory}' vocabulary",
            )
            subcats.setdefault(tag, set()).add(str(row.subcategory))
            if bool(getattr(row, "as_cluster", True)):
                cluster_subcats.setdefault(tag, set()).add(str(row.subcategory))
            machinery_classes.setdefault(tag, set()).add(
                str(getattr(row, "machinery", "unclassified"))
            )
            if str(row.axis_hint) == "host_interaction":
                bump(host_axis, tag, float(row.weight) * 0.5)
                host_specific.add(tag)
                host_ev.setdefault(tag, []).append(f"vocabulary: {row.subcategory}")
            elif str(row.axis_hint) == "virulence_general":
                core_ev.setdefault(tag, []).append(f"vocabulary: {row.subcategory}")

    if core_vocabulary_hits is not None and not core_vocabulary_hits.empty:
        for row in core_vocabulary_hits.itertuples(index=False):
            tag = str(row.locus_tag)
            if tag not in tag_set:
                continue
            acc.add(tag, params.w_core_vocabulary, "annotation indicates a core cellular function")
            bump(core_axis, tag, 1.0)
            core_ev.setdefault(tag, []).append(f"core-process annotation: {row.note or row.pattern}")

    # ------------------------------- E6 species-specific acquired loci
    # Published, species-scoped locus definitions applied to the RefSeq
    # annotation (see codonpipe.virulence.species_loci). The locus has already
    # been established there by several co-located annotation matches, so a
    # member gene enters the score as curation-anchored evidence rather than as
    # a weak annotation hint -- the curation being the published locus
    # definition itself. Empty or None leaves everything below unchanged.
    if species_locus_hits is not None and not species_locus_hits.empty:
        for row in species_locus_hits.itertuples(index=False):
            tag = str(getattr(row, "locus_tag", "") or "")
            if tag not in tag_set:
                continue
            locus = str(getattr(row, "locus", "") or "")
            role = str(getattr(row, "role", "marker") or "marker")
            if role == "marker":
                value = float(getattr(row, "weight", 2.0)) * params.w_species_locus_scale
                label = f"species virulence locus {locus} (annotation matches the published definition)"
            else:
                value = params.w_species_locus_context
                label = f"species virulence locus {locus} (member of the detected locus block)"
            acc.add(tag, value, label)
            species_locus_member.setdefault(tag, set()).add(locus)
            # A marker role always wins over a recruited role in the report.
            if role == "marker" or tag not in species_locus_role:
                species_locus_role[tag] = role
            machinery_classes.setdefault(tag, set()).add(
                str(getattr(row, "machinery", "unclassified") or "unclassified")
            )
            db_ev.setdefault(tag, []).append(
                f"species locus {locus} [{role}] - {getattr(row, 'reference', '')}"
            )
            if str(getattr(row, "axis_hint", "")) == "host_interaction":
                bump(host_axis, tag, params.w_species_locus_host_axis)
                host_specific.add(tag)
                host_ev.setdefault(tag, []).append(
                    f"member of the {locus} acquired virulence locus"
                )

    # ------------------------------------------------- E5 genomic context
    context_layer: Dict[str, str] = {}
    context_detail: Dict[str, str] = {}
    if params.use_context:
        seeds = [t for t, v in acc.score.items() if v >= params.context_seed_threshold]
        ctx = extend_by_context(
            gene_table,
            seeds,
            operon_extension=True,
            max_intergenic_gap=params.max_intergenic_gap,
            require_same_strand=params.require_same_strand,
            island_bridging=True,
            island_bridge_genes=params.island_bridge_genes,
            min_island_seeds=params.min_island_seeds,
        )
        notes.extend(ctx.notes)
        for tag, info in ctx.added.items():
            if tag not in tag_set:
                continue
            layer = str(info.get("layer"))
            weight = (
                params.w_context_operon if layer == "operon" else params.w_context_island
            )
            acc.add(tag, weight, f"genomic context ({layer})")
            context_layer[tag] = layer
            context_detail[tag] = str(info.get("detail", ""))
    else:
        notes.append("Genomic-context extension was disabled at this stringency.")

    # ------------------------------------------------- experimental override
    #
    # A measured gene requirement beats an inferred one, in both directions:
    # an essential gene cannot yield an interpretable loss-of-virulence
    # phenotype however host-directed its annotation looks, and a gene whose
    # viable mutant has been built and assayed is dispensable in vitro
    # however much core-process annotation it carries.
    essential_map: Dict[str, str] = {}
    invivo_map: Dict[str, str] = {}
    invivo_basis: Dict[str, str] = {}
    invivo_model: Dict[str, str] = {}
    invivo_effect: Dict[str, float] = {}
    invivo_sig: Dict[str, float] = {}
    invivo_rank: Dict[str, str] = {}
    if essentiality is not None and not essentiality.empty:
        cols = {c.lower(): c for c in essentiality.columns}
        tag_col = cols.get("locus_tag")
        req_col = cols.get("in_vitro_requirement")
        vivo_col = cols.get("in_vivo_fitness")
        basis_col = cols.get("in_vivo_basis")
        model_col = cols.get("in_vivo_model")
        effect_col = cols.get("in_vivo_effect_log2")
        sig_col = cols.get("in_vivo_significance")
        rank_col = cols.get("in_vivo_effect_rank")
        if tag_col and (req_col or vivo_col):
            for row in essentiality.itertuples(index=False):
                tag = gene_table.lookup(getattr(row, tag_col))
                if not tag:
                    continue
                if req_col:
                    value = str(getattr(row, req_col)).strip().lower()
                    if value and value != "unknown":
                        essential_map[tag] = value
                if vivo_col:
                    value = str(getattr(row, vivo_col)).strip().lower()
                    if value and value != "not_assayed":
                        invivo_map[tag] = value
                if basis_col:
                    invivo_basis[tag] = str(getattr(row, basis_col)).strip()
                if effect_col:
                    try:
                        value = float(getattr(row, effect_col))
                        if value == value:
                            invivo_effect[tag] = value
                    except (TypeError, ValueError):
                        pass
                if sig_col:
                    try:
                        value = float(getattr(row, sig_col))
                        if value == value:
                            invivo_sig[tag] = value
                    except (TypeError, ValueError):
                        pass
                if rank_col:
                    value = str(getattr(row, rank_col)).strip().lower()
                    if value and value not in ("nan", "none"):
                        invivo_rank[tag] = value
                if model_col:
                    value = str(getattr(row, model_col)).strip().lower()
                    if value and value != "not_applicable":
                        invivo_model[tag] = value
            notes.append(
                f"An experimental gene-requirement table was applied: in vitro "
                f"requirement for {len(essential_map)} gene(s), in vivo fitness "
                f"for {len(invivo_map)} gene(s). Measured values override the "
                f"annotation-based inference."
            )
        else:
            notes.append(
                "The supplied experimental table lacks a 'locus_tag' column "
                "together with 'in_vitro_requirement' or 'in_vivo_fitness', "
                "and was ignored."
            )

    # Resolve what the strict tier will demand of the third postulate.
    coverage = (len(set(essential_map) | set(invivo_map)) / len(gene_table.frame)) if len(gene_table.frame) else 0.0
    # Rung 3 has a single definition; see the note at the tier assignment.
    # `coverage` is still reported, because a reader needs to know how much of
    # the genome the experimental rungs could speak for.
    strict_rule = "evidence"
    if essential_map and not any(v == "fitness_cost" for v in essential_map.values()):
        notes.append(
            "RUNG 6'S IN VITRO CLAUSE DID NOT FIRE IN THIS RUN. The experimental "
            "table reports in vitro requirement as BINARY essentiality, not as "
            "graded in vitro fitness, so no gene is marked 'fitness_cost' and the "
            "clause removes nothing; 'dispensable' here means only that a viable "
            "mutant was obtained, which does not exclude a growth or stress-"
            "survival defect short of lethality. Rung 6 is therefore carried "
            "entirely by its host-directedness clause in this run. Do not read it "
            "as evidence that these genes have no in vitro cost. A quantitative "
            "in vitro fitness score, ideally in a medium resembling the host "
            "niche, loaded as in_vitro_requirement='fitness_cost', is what would "
            "make the in vitro clause mean anything."
        )

    # ---- resolve the effect-size cut-offs -----------------------------------
    # In percentile mode the two cut-offs are read against the distribution of
    # measured effect sizes among genes this organism's screens actually called
    # attenuated. That keeps the rung 4 / rung 5 split comparable between a
    # competitive-index study and a pooled transposon screen, which sit on
    # different numeric scales. The resolved absolute values are reported, so a
    # percentile run can still be stated in log2 units.
    cut_contributing = float(params.effect_contributing_log2)
    cut_required = float(params.effect_required_log2)
    effect_scale_note = "absolute log2 depletion"
    if str(params.effect_scale).lower().startswith("perc"):
        # The two cut-offs are reinterpreted as percentiles, so leaving them at
        # their ABSOLUTE defaults (0.0 and 1.0) would mean P0 and P1 -- the
        # minimum of the observed distribution -- and rung 5 would collapse
        # into rung 4. Refuse rather than silently produce that.
        if (
            params.effect_contributing_log2 == ScoringParams.effect_contributing_log2
            and params.effect_required_log2 == ScoringParams.effect_required_log2
        ):
            raise ValueError(
                "effect_scale='percentile' reinterprets the two cut-offs as "
                "percentiles, but they are still at their absolute-scale defaults "
                f"({params.effect_contributing_log2} and {params.effect_required_log2}), "
                "which would mean the 0th and 1st percentiles. Set them to "
                "percentiles, for example 0 and 90."
            )
        measured = sorted(
            value
            for tag, value in invivo_effect.items()
            if value is not None and invivo_map.get(tag) == "attenuated"
        )
        if measured:
            def _pct(p: float) -> float:
                p = min(max(float(p), 0.0), 100.0)
                if len(measured) == 1:
                    return measured[0]
                pos = (p / 100.0) * (len(measured) - 1)
                low = int(pos)
                high = min(low + 1, len(measured) - 1)
                return measured[low] + (measured[high] - measured[low]) * (pos - low)

            cut_contributing = _pct(params.effect_contributing_log2)
            cut_required = _pct(params.effect_required_log2)
            effect_scale_note = (
                f"percentile of {len(measured)} measured effect sizes; "
                f"P{params.effect_contributing_log2:g} = {cut_contributing:.2f} log2, "
                f"P{params.effect_required_log2:g} = {cut_required:.2f} log2"
            )
        else:
            effect_scale_note = (
                "percentile requested, but this organism has no measured effect "
                "sizes; the cut-offs were left as absolute log2 values"
            )

    # ------------------------------------------------------------ assembly
    rows: List[dict] = []
    for row in gene_table.frame.itertuples(index=False):
        tag = str(row.locus_tag)
        score = round(acc.score.get(tag, 0.0), 3)
        h_axis = round(host_axis.get(tag, 0.0), 3)
        c_axis = round(core_axis.get(tag, 0.0), 3)
        axes = vf_axes.get(tag, set())

        # Postulate 1: association with pathogenic strains, from curation.
        if tag in has_literature or tag in has_uniprot_virulence:
            p1 = "literature_curated"
        elif tag in has_database:
            p1 = "database_homology"
        elif tag in species_locus_member:
            # Membership of a published, species-scoped acquired virulence
            # locus. The curation is the published locus definition rather than
            # a per-protein database record, so it gets its own value instead
            # of being reported as database homology.
            p1 = "species_locus_definition"
        elif tag in context_layer:
            p1 = "genomic_context_only"
        elif tag in subcats:
            p1 = "annotation_vocabulary_only"
        else:
            p1 = "none"

        # Postulate 2: is a loss-of-virulence phenotype interpretable?
        override = essential_map.get(tag)
        if override == "essential":
            # No viable mutant, so no loss-of-virulence phenotype is
            # obtainable and postulate 2 cannot be read at all.
            p2 = "confounded"
            p2_basis = "experimental in vitro requirement: essential"
        elif override == "fitness_cost":
            # A viable but growth-impaired mutant. Postulate 2 IS readable --
            # the mutant exists and can be assayed -- so this gene belongs in
            # rung 3. What its in vitro cost threatens is the interpretation
            # of an in vivo defect, which is precisely what rung 5 removes.
            # Treating it as confounded here would collapse rung 5 into rung 4
            # and make the fifth rung redundant by construction.
            p2 = "interpretable"
            p2_basis = (
                "experimental in vitro requirement: fitness_cost - viable but "
                "impaired in vitro, so a mutant phenotype is obtainable; "
                "excluded at rung 5, where the in vivo defect must not be "
                "attributable to a general growth defect"
            )
        elif override == "dispensable":
            p2 = "interpretable"
            p2_basis = "experimental in vitro requirement: dispensable"
        elif c_axis > params.max_core_axis:
            p2 = "confounded"
            p2_basis = "; ".join(core_ev.get(tag, [])[:4]) or "core-process annotation"
        elif h_axis >= params.min_host_axis or tag in has_database:
            p2 = "interpretable"
            if core_ev.get(tag):
                # Some core-process annotation exists but stays within the
                # tolerance, so it does not block a host-specific reading.
                p2_basis = (
                    f"core-process evidence present but below the tolerance "
                    f"({c_axis:.2f} <= {params.max_core_axis:.2f}): "
                    + "; ".join(dict.fromkeys(core_ev[tag]))[:240]
                )
            elif core_axis.get(tag, 0.0) > 0:
                # Core-process evidence exists but sits exactly at the
                # tolerance rather than above it. Saying "no evidence" here
                # would contradict the gene's own core_evidence column.
                p2_basis = (
                    f"core-process evidence present at the tolerance "
                    f"({c_axis:.2f} = {params.max_core_axis:.2f}), which does "
                    f"not block a host-specific reading"
                )
            else:
                p2_basis = "no evidence of a core cellular function"
        else:
            p2 = "not_assessed"
            p2_basis = "insufficient functional annotation"

        # Postulate 2 SATISFIED: a measured loss of virulence on specific
        # inactivation. Falkow's second postulate is the one that asks for an
        # attenuation phenotype; his third asks for restoration by reversion
        # or allelic replacement, which no database records per gene and
        # which this pipeline therefore reports as not_assessed throughout.
        # Rung 3 tests whether postulate 2 is *interpretable*; this field
        # records whether it is *satisfied*.
        # Never computable from sequence data, but readable from a published
        # mutant screen where one covers the gene. The three experimental
        # states are kept distinct: a gene that was assayed and showed no
        # phenotype is a different thing from a gene nobody assayed, and
        # effector redundancy makes the difference matter.
        vivo = invivo_map.get(tag, "")
        model = invivo_model.get(tag, "")
        # An in vivo call can be restricted to particular infection models.
        # Left empty, every model counts and the model is merely recorded.
        model_filtered = bool(
            vivo
            and params.in_vivo_accepted_models
            and model
            and model not in params.in_vivo_accepted_models
        )
        if model_filtered:
            # Record WHY the call was dropped. Reporting it as "not_assayed"
            # would conflate "assayed in a model you excluded" with "nobody
            # assayed it", which is the exact distinction the reference tables
            # exist to preserve.
            vivo = ""
        if vivo == "attenuated":
            p3 = "demonstrated_in_vivo"
            p3_basis = invivo_basis.get(tag, "") or "mutant attenuated in an animal infection screen"
        elif vivo == "no_phenotype":
            p3 = "refuted_in_vivo"
            p3_basis = "mutant assayed in an animal infection screen with no fitness defect"
        elif tag in has_literature:
            p3 = "reported_in_source_literature"
            p3_basis = "asserted by a literature-curated virulence record"
        else:
            p3 = "not_assessed"
            p3_basis = "no mutant screen covers this gene"

        # ---------------- the tier ladder, three nested predicates -----------
        # Each tier adds one predicate to the previous one, so
        # strict subset standard subset permissive holds structurally.

        # permissive: postulate 1 read inclusively. Any evidence source counts,
        # including an annotation match standing on its own.
        t_permissive = bool(p1 != "none" and score >= params.threshold_permissive)

        # standard: postulate 1 anchored in curation. Annotation-only
        # candidates drop out; operon context around a curated gene is kept,
        # because the anchoring curation is the neighbouring gene's.
        curation_anchored = p1 in (
            "literature_curated",
            "database_homology",
            "species_locus_definition",
            "genomic_context_only",
        )
        t_standard = bool(
            t_permissive and curation_anchored and score >= params.threshold_standard
        )

        # strict: postulate 2 interpretable. Host-directed function, and no
        # core cellular role that would confound a loss-of-virulence phenotype.
        # Retained for the classification_reason text only. It is NOT part of
        # the gate any more: it was computed from the virulence-factor axis set
        # alone, and vf_category_axes.tsv maps exactly one term to
        # "uninformative" -- BV-BRC's generic classification "virulence",
        # carried by 213 of 1003 SL1344 records. A gene with an uncategorised
        # curated record was therefore vetoed while a gene with NO curated
        # record at all was never tested, making an uncategorised record
        # strictly worse than none. The host_specific set below tests the same
        # idea from the right input: the provenance of the host-axis credit.
        general_only = bool(axes) and not (axes - {"virulence_general", "uninformative", "negative"})
        host_directed = (
            h_axis >= params.min_host_axis
            and (tag in host_specific or not params.require_specific_host_evidence)
            and not (params.require_host_axis_exceeds_core and h_axis <= c_axis)
        )
        mach = machinery_classes.get(tag, set())
        machinery = (
            "dedicated" if "dedicated" in mach
            else ("dual_use" if "dual_use" in mach
                  else ("regulatory" if "regulatory" in mach
                        else ("unclassified" if mach else "")))
        )
        # The dual-use question: fimbriae, flagella, siderophores, O-antigen and
        # capsule are host-facing but also serve biofilm and environmental
        # persistence, so whether they satisfy postulate 2 is a judgement call.
        # Genes with no machinery classification are never excluded here.
        dedicated_ok = not (
            params.strict_requires_dedicated and machinery == "dual_use"
        )
        # The annotation-based host-interaction set: curation-anchored,
        # host-directed, and free of a confounding core cellular role. This
        # column is always written and never gated, because it is the one
        # that stays comparable between a genome with mutant-screen data and
        # a genome without.
        t_strict_annotation = bool(
            t_standard and host_directed and p2 == "interpretable" and dedicated_ok
        )

        # ---- the experimental rungs 4 to 6 --------------------------------
        # These sit inside rung 2, not rung 3: they refine the
        # curation-anchored set with experimental evidence. See
        # LADDER_BRANCHES for why the gate is rung 2. A gene must
        # also be inactivatable, since an essential gene yields no
        # interpretable phenotype and its apparent in vivo defect is a growth
        # defect.
        effect = invivo_effect.get(tag)
        signif = invivo_sig.get(tag)
        significant = (
            signif is None or signif <= params.max_in_vivo_significance
        )
        inactivatable = essential_map.get(tag) != "essential"
        attenuated = p3 == "demonstrated_in_vivo"

        # Rung 4: contributing. With an effect size, it must clear the
        # contributing cut; without one, the source's own attenuation call is
        # all there is and the gene enters on that.
        t_in_vivo_contributing = bool(
            t_standard
            and attenuated
            and inactivatable
            and significant
            and (effect is None or effect >= cut_contributing)
        )
        # Rung 5: required. Two routes in, because the screens behind these
        # tables do not all report the same kind of number.
        #   - a magnitude that clears the larger cut, or
        #   - the source's own TOP attenuation category, for screens that
        #     grade attenuation ordinally and publish no ratio.
        # Without the second route this rung would silently mean "has a large
        # competitive index" rather than "is required in vivo", and would drop
        # genes such as aroC, purA and guaB that no competitive index covers.
        # Silence is still not promotion: a gene with neither a magnitude nor
        # a strength call stays at rung 4.
        rank = invivo_rank.get(tag, "")
        t_in_vivo_required = bool(
            t_in_vivo_contributing
            and (
                (effect is not None and effect >= cut_required)
                or rank == "strong"
            )
        )
        # Rung 5: and no measured fitness cost in vitro, so the in vivo
        # phenotype cannot be a general growth defect. Rung 3 already excluded
        # genes called 'essential', so the only value this rung removes is
        # 'fitness_cost' -- which needs GRADED in vitro fitness data. Where a
        # reference table reports essentiality only, this rung is identical to
        # rung 4, and the run summary says so rather than implying a filter
        # that did not actually run.
        # An ABSENT in vitro measurement is not a measurement of no cost.
        # Requiring a positive call keeps "no measured fitness cost" honest;
        # where a table carries no in vitro column at all the rung is reported
        # as degenerate rather than silently filled.
        # ...AND host-directed with no confounding core cellular role, i.e.
        # the rung-3 test. This is the ONLY rung that takes rung 3 as well as
        # rung 2, and it is deliberate. Rungs 4 and 5 must not, because a
        # measured in vivo phenotype outranks an annotation-based inference
        # and gating them would delete the O-antigen block and the global
        # regulators from the required set. Rung 6 asks a different question:
        # is the in vivo defect SPECIFIC to the host interaction, or could it
        # be a general fitness effect? rpoS is the case that forces this.
        # With binary essentiality it is 'dispensable', so the in vitro test
        # passes it; and there is no threshold on the TraDIS in vitro log P
        # that would catch it, because rpoS sits at -0.0000137 while every
        # one of the 28 SPI-2 genes in rung 5 reads more depleted than that
        # (least depleted: ssaP, -0.54). The ordering is inverted relative to
        # the biology, so any cut that removes rpoS removes the whole SPI-2
        # apparatus. The annotation axes separate them cleanly instead:
        # rpoS carries host 0.00 / core 1.10, the SPI-1 and SPI-2 genes a
        # median host axis of 2.55. Keeps 28/28 SPI-2 and 22/22 SPI-1.
        t_in_vitro_neutral = bool(
            t_in_vivo_required
            and t_strict_annotation
            and tag in essential_map
            and essential_map.get(tag) not in ("essential", "fitness_cost")
        )

        # Rung 3 has ONE definition. It used to have four, selected by a
        # `strict_requires` parameter whose `auto` default switched between
        # them depending on how much of the genome a reference table covered --
        # so the same tier name meant different things on different genomes
        # without the operator choosing. It also emitted a second column,
        # `tier_strict`, that disagreed with `tier_host_interaction` in the
        # same file (88 against 255 for SL1344). Experimental evidence now has
        # its own rungs, 4 to 6, so the parameter was redundant as well as
        # opaque and has been removed.
        tier = "strict" if t_strict_annotation else ("standard" if t_standard else ("permissive" if t_permissive else ""))
        in_general, in_host = t_standard, t_strict_annotation

        if t_strict_annotation:
            reason = "curation-anchored virulence evidence, host-directed function, no core-process evidence"
        elif t_standard and general_only:
            reason = "virulence-associated, but the annotated role is also engaged during growth in vitro"
        elif t_standard and p2 == "confounded":
            reason = "virulence-associated, but core cellular function prevents a host-specific interpretation"
        elif t_standard and host_directed and p2 == "interpretable" and not dedicated_ok:
            reason = (
                "host-directed, but only through dual-use machinery that also serves "
                "biofilm or environmental persistence"
            )
        elif t_standard and h_axis < params.min_host_axis:
            reason = "virulence-associated, but no positive evidence of direct host engagement"
        elif t_standard:
            reason = "virulence-associated"
        elif t_permissive and not curation_anchored:
            reason = "annotation evidence only, with no curated database record to anchor it"
        elif t_permissive:
            reason = (
                f"evidence score {score:.2f} below the standard-tier threshold "
                f"{params.threshold_standard:.2f}"
            )
        elif p1 == "none":
            reason = "no virulence evidence"
        else:
            reason = (
                f"evidence score {score:.2f} below the permissive-tier threshold "
                f"{params.threshold_permissive:.2f}"
            )

        rows.append(
            {
                "order": getattr(row, "order", None),
                "locus_tag": tag,
                "gene": row.gene,
                "product": row.product,
                "old_locus_tag": row.old_locus_tag,
                "protein_id": row.protein_id,
                "seqid": row.seqid,
                "start": row.start,
                "end": row.end,
                "strand": row.strand,
                "evidence_score": score,
                "host_interaction_axis": h_axis,
                "core_process_axis": c_axis,
                "tier": tier,
                "tier_permissive": t_permissive,
                "tier_standard": t_standard,
                # Retained for continuity with CodonPipe 1.2.0 column names.
                "virulence_general": in_general,
                "host_interaction": in_host,
                "koch_postulate_1": p1,
                "koch_postulate_2": p2,
                "koch_postulate_2_basis": p2_basis,
                "koch_postulate_2_satisfied": p3,
                "koch_postulate_3": "not_assessed",
                "koch_postulate_2_satisfied_basis": p3_basis,
                "koch_postulate_3_basis": "complementation evidence is not machine-readable",
                # The six-rung ladder.
                "tier_host_interaction": t_strict_annotation,
                "tier_in_vivo_contributing": t_in_vivo_contributing,
                "tier_in_vivo_required": t_in_vivo_required,
                # Attenuated in vivo but outside rung 3, so not admitted to rungs
            # 4-6. Kept as its own column because these are exactly the genes
            # the annotation gate costs us, and a reader is entitled to see
            # them rather than be told they do not exist.
            "in_vivo_outside_rung2": bool(
                attenuated and inactivatable and significant and not t_standard
            ),
            # A gene the reference table calls attenuated while the effect
            # size it carries is NEGATIVE -- i.e. one screen reports a defect
            # and another reports the mutant was ENRICHED in the host. The
            # effect gate then silently keeps the gene out of rungs 4-6. That
            # may be right (a direct competitive index contradicts a
            # categorical call) or wrong (the two screens are different
            # models: 28-day persistence and acute competition can genuinely
            # disagree), so the conflict is reported rather than hidden.
            # 36 SL1344 genes, 25 of them in rung 2.
            "in_vivo_discordant": bool(
                attenuated and effect is not None and effect < 0
            ),
            "in_vivo_effect_log2": effect,
            "in_vivo_effect_rank": rank,
                "in_vivo_significance": signif,
                "tier_in_vitro_neutral": t_in_vitro_neutral,
                "in_vivo_fitness": (
                    vivo
                    or ("excluded_by_model_filter" if model_filtered else "not_assayed")
                ),
                "in_vivo_model": model if vivo else "not_applicable",
                # Purely experimental, outside the ladder. See
                # IN_VIVO_SPECIFIC_LABEL for why it is not a host-interaction
                # set.
                "in_vivo_specific": bool(
                    vivo == "attenuated"
                    and essential_map.get(tag) not in ("essential", "fitness_cost")
                ),
                "machinery": machinery,
                "species_locus": "; ".join(sorted(species_locus_member.get(tag, set()))),
                "species_locus_role": species_locus_role.get(tag, ""),
                "vf_axes": "; ".join(sorted(axes)) if axes else "",
                "subcategories": "; ".join(sorted(subcats.get(tag, set()))),
                "cluster_subcategories": "; ".join(sorted(cluster_subcats.get(tag, set()))),
                "context_layer": context_layer.get(tag, ""),
                "context_detail": context_detail.get(tag, ""),
                "database_evidence": " | ".join(db_ev.get(tag, [])[:12]),
                "host_evidence": " | ".join(dict.fromkeys(host_ev.get(tag, []))),
                "core_evidence": " | ".join(dict.fromkeys(core_ev.get(tag, []))),
                "score_components": " | ".join(acc.components.get(tag, [])[:20]),
                "classification_reason": reason,
                "in_vitro_requirement_source": (
                    "experimental" if tag in essential_map else "annotation_inferred"
                ),
            }
        )

    frame = pd.DataFrame(rows)

    # The ladder is nested by construction; verify it so that a future edit to
    # the rule tables or the weights cannot quietly produce three unrelated
    # gene sets while the documentation still promises a ladder.
    n_perm = int(frame["tier_permissive"].sum())
    n_std = int(frame["tier_standard"].sum())
    # Every rung must sit inside the one below it. Checked on the frame, not
    # argued from the predicates, so that a future edit cannot quietly produce
    # five unrelated gene sets while the documentation still promises a ladder.
    # Each BRANCH nests internally. Rung 4 is deliberately not checked
    # against rung 3: the experimental branch supersedes the annotation one,
    # so a gene can be in rung 4 and not rung 3, and that is the point.
    chains = [
        ["tier_" + r for r in LADDER_BRANCHES["annotation"]],
        ["tier_" + r for r in LADDER_BRANCHES["experimental"]],
    ]
    broken = []
    for chain in chains:
        for i in range(len(chain) - 1):
            n_bad = int((frame[chain[i + 1]] & ~frame[chain[i]]).sum())
            if n_bad:
                broken.append((chain[i + 1], chain[i], n_bad))
    if broken:
        raise AssertionError(
            "The virulence ladder is not nested: "
            + "; ".join(f"{n} gene(s) in {a} but not in {b}" for a, b, n in broken)
            + ". This indicates a bug in the rung predicates, not a data problem."
        )

    violations = (
        int((frame["tier_standard"] & ~frame["tier_permissive"]).sum()),
        int((frame["tier_host_interaction"] & ~frame["tier_standard"]).sum()),
    )
    # Rung 5 must sit inside the purely experimental set: both require a
    # measured in vivo phenotype AND in vitro neutrality, and rung 5 adds the
    # annotation-based tests on top. Rung 4 is deliberately NOT checked this
    # way -- it may contain a gene with an in vitro fitness cost, which is
    # exactly what rung 5 and the experimental set both remove.
    escaped = int((frame["tier_in_vitro_neutral"] & ~frame["in_vivo_specific"]).sum())
    if escaped:
        raise AssertionError(
            f"{escaped} gene(s) are in rung 5 but not in the in-vivo-specific "
            f"set, which cannot happen: both require a measured in vivo "
            f"phenotype and no measured in vitro fitness cost."
        )
    if any(violations):
        raise AssertionError(
            "The virulence tier ladder is not nested: "
            f"{violations[0]} gene(s) are standard but not permissive, "
            f"{violations[1]} are strict but not standard, "
            f"{violations[2]} are annotation-strict but not standard and "
            f"{violations[3]} are strict but not annotation-strict. This "
            "indicates a bug in the tier predicates rather than a data problem."
        )

    general = frame[frame["virulence_general"]]
    host = frame[frame["host_interaction"]]
    summary = {
        "n_genes": int(len(frame)),
        "n_tier_permissive": n_perm,
        "n_tier_standard": n_std,
                "tier_ladder_nested": True,
        "n_annotation_only_candidates": int(n_perm - n_std),
        "strict_requires_dedicated": bool(params.strict_requires_dedicated),
        "n_tier_host_interaction": int(frame["tier_host_interaction"].sum()),
        "experimental_table_coverage": round(coverage, 4),
        "n_postulate2_demonstrated": int(
            (frame["koch_postulate_2_satisfied"] == "demonstrated_in_vivo").sum()
        ),
        "ladder": {
            name: int(frame[f"tier_{name}"].sum()) for name in LADDER
        },
        "ladder_branches": {k: list(v) for k, v in LADDER_BRANCHES.items()},
        # How much of the experimental branch the annotation branch had
        # already found, and how much it had missed. This is the number that
        # says what rung 3 is worth.
        "n_rung4_also_in_rung3": int(
            (frame["tier_in_vivo_contributing"] & frame["tier_host_interaction"]).sum()
        ),
        # Attenuated in vivo but excluded from rung 4 by the rung-3 gate. The
        # price of the gate, reported so it is never silent.
        "n_in_vivo_outside_rung2": int(frame["in_vivo_outside_rung2"].sum()),
        # How much of the genome this organism's screens actually called.
        # Rungs 4 to 6 inherit the ORIGINAL authors' calling thresholds, which
        # differ enormously between studies, so this rate is the first thing
        # to read before comparing those rungs with another organism's.
        "in_vivo_assayed": int(sum(1 for v in invivo_map.values() if v in ("attenuated", "no_phenotype"))),
        "in_vivo_attenuated": int(sum(1 for v in invivo_map.values() if v == "attenuated")),
        "in_vivo_screen_call_rate": (
            round(
                100.0
                * sum(1 for v in invivo_map.values() if v == "attenuated")
                / max(sum(1 for v in invivo_map.values() if v in ("attenuated", "no_phenotype")), 1),
                1,
            )
            if invivo_map
            else None
        ),
        "n_rung3_without_in_vivo_support": int(
            (frame["tier_host_interaction"] & ~frame["tier_in_vivo_contributing"]).sum()
        ),
        "effect_scale": params.effect_scale,
        "effect_scale_note": effect_scale_note,
        "effect_contributing_log2": cut_contributing,
        "effect_required_log2": cut_required,
        "effect_contributing_setting": params.effect_contributing_log2,
        "effect_required_setting": params.effect_required_log2,
        "max_in_vivo_significance": params.max_in_vivo_significance,
        "ladder_definitions": dict(LADDER_DEFINITIONS),
        "n_tier_host_interaction": int(frame["tier_host_interaction"].sum()),
        "n_tier_in_vivo_contributing": int(frame["tier_in_vivo_contributing"].sum()),
        "n_tier_in_vivo_required": int(frame["tier_in_vivo_required"].sum()),
        "n_tier_in_vitro_neutral": int(frame["tier_in_vitro_neutral"].sum()),
        # Rung 5 can only differ from rung 4 if the experimental table reports
        # graded in vitro fitness. Say whether it actually did.
        "in_vitro_fitness_data_present": bool(
            any(v == "fitness_cost" for v in essential_map.values())
        ),
        # Genes where the categorical in vivo call and the effect size point
        # in opposite directions. Reported so that a reader can see how many
        # rung-4 exclusions rest on one screen overruling another.
        "n_in_vivo_discordant": int(frame["in_vivo_discordant"].sum()),
        "n_in_vivo_discordant_in_rung2": int(
            (frame["in_vivo_discordant"] & frame["tier_standard"]).sum()
        ),
        "in_vivo_models": {
            m: int(((frame["koch_postulate_2_satisfied"] == "demonstrated_in_vivo")
                    & (frame["in_vivo_model"] == m)).sum())
            for m in sorted(set(frame["in_vivo_model"]) - {"not_applicable", ""})
        },
        "in_vivo_accepted_models": list(params.in_vivo_accepted_models or ()),
        "n_in_vivo_specific": int(frame["in_vivo_specific"].sum()),
        # Named for the common case; the denominator also contains genes
        # excluded for a measured in vitro fitness cost, not only essentials.
        "n_in_vivo_specific_excluded_as_not_inactivatable": int(
            ((frame["koch_postulate_2_satisfied"] == "demonstrated_in_vivo")
             & ~frame["in_vivo_specific"]).sum()
        ),
        "n_postulate2_refuted": int(
            (frame["koch_postulate_2_satisfied"] == "refuted_in_vivo").sum()
        ),
        "n_host_interaction_dropped_by_postulate2": int(
            (frame["tier_host_interaction"] & ~frame["tier_host_interaction"]).sum()
        ),
        "n_host_interaction_dedicated": int(
            ((frame["tier_host_interaction"]) & (frame["machinery"] == "dedicated")).sum()
        ),
        "n_host_interaction_dual_use": int(
            ((frame["tier_host_interaction"]) & (frame["machinery"] == "dual_use")).sum()
        ),
        "n_virulence_general": int(len(general)),
        "n_host_interaction": int(len(host)),
        "n_general_only": int(len(general) - len(host)),
        "n_with_database_evidence": int(len(has_database & tag_set)),
        "n_literature_curated": int(len(has_literature & tag_set)),
        "n_uniprot_virulence_keyword": int(len(has_uniprot_virulence & tag_set)),
        "n_context_added": int(len(context_layer)),
        "n_species_locus_genes": int(len(set(species_locus_member) & tag_set)),
        "n_species_locus_marker_genes": int(
            sum(1 for t, r in species_locus_role.items() if r == "marker" and t in tag_set)
        ),
        "n_negative_records": int(len(negative_record & tag_set)),
        "stringency": params.stringency,
        "threshold_permissive": params.threshold_permissive,
        "threshold_standard": params.threshold_standard,
        "tier_definitions": tier_definitions(strict_rule),
        "notes": notes,
    }
    return frame, summary
