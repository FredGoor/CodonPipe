# codonpipe/codon_sensitivity.py
"""Single source of truth for which codons count as starvation-sensitive or
starvation-insensitive, and a switch between the available definitions.

Background
----------
Elf, Nilsson & Ehrenberg (PMID 12805541) measured, for individual codons, how
strongly ribosome occupancy rises when the cognate amino acid runs out. Codons
read by isoacceptors that stay charged during starvation barely respond
(low sensitivity); codons read by isoacceptors that are deacylated respond
strongly (high sensitivity). Amino-acid families showing this split are
Ala, Gly, Ile, Leu, Pro, Arg, Ser and Thr.

Every percentage computed from these lists is a *synonymous* quantity: the
denominator is all codons of the amino-acid families the list spans, never all
61 sense codons. A gene cannot score high merely by being rich in one amino
acid.

Definitions
-----------
``legacy``
    What the pipeline used before 2026-11-12. Preserved verbatim so that
    re-running an old analysis reproduces its numbers exactly.

    Note that "legacy" is not one list. Two parts of the pipeline carried
    different twelve-codon sets, and this module keeps both rather than
    silently reconciling them:

      * per-gene metrics (``perc_insensitive``) used
        AGA AGG CGG GGA GGG ATA TTA TTG AGC AGT TCG ACG
      * the per-cluster bars (PLOT 4) used
        TTG TTA ATA ACG AGT GCA GCT CCT AGA AGG GGA GGG

    They share nine codons (ACG AGA AGG AGT ATA GGA GGG TTA TTG) and differ in
    three on each side: the gene metric adds AGC, CGG and TCG, the bar set adds
    CCT, GCA and GCT. The bar set came from an exhaustive subset search against
    curated virulence genes; the gene-metric set did not. Results computed with one are not interchangeable with the
    other, which is the reason this module exists.

``custom``
    The refined definition supplied by the user on 2026-11-12, taken from
    "Sensitive and insensitive codons table.xlsx". Each codon carries its
    published Elf sensitivity (see ``ELF_SENSITIVITY``) and is classified by
    that number alone:

      * insensitive: sensitivity < 4   (ten codons, range 0.0 - 3.7)
      * sensitive:   sensitivity > 12  (fourteen codons, range 12.5 - 35.5)

    Codons whose sensitivity falls between the two cut-offs, and every codon of
    the Ala family, are left unclassified: they contribute to the denominator
    but to neither numerator. Alanine is absent from the table entirely, so the
    custom definition spans seven families, not eight, and Ala codons are
    excluded from the denominator as well.

    Unlike ``legacy``, one list serves both the per-gene metric and PLOT 4.

Three further sets are kept for reference and remain selectable by name:
``published`` (all sixteen of PMID 12805541), ``optimum`` (the eight-codon
subset search optimum) and ``acgfree`` (that optimum without ACG).

Selecting a definition
----------------------
``Plotting_Pipeline.CODON_SENSITIVITY_DEFINITION`` or the environment variable
``CODONPIPE_CODON_DEFINITION`` chooses the active definition by key. Code that
needs the lists calls :func:`get_definition`; code that needs the active one
calls :func:`active_definition`.
"""

from __future__ import annotations

import os
from typing import Dict, List, Optional, Sequence


# --------------------------------------------------------------------------
# Amino-acid families with selective isoacceptor charging (PMID 12805541)
# --------------------------------------------------------------------------

SELECTIVE_FAMILIES: Dict[str, List[str]] = {
    "Ala": ["GCT", "GCC", "GCA", "GCG"],
    "Gly": ["GGT", "GGC", "GGA", "GGG"],
    "Ile": ["ATT", "ATC", "ATA"],
    "Leu": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "Pro": ["CCT", "CCC", "CCA", "CCG"],
    "Arg": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "Ser": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "Thr": ["ACT", "ACC", "ACA", "ACG"],
}

#: Family membership for every codon of the eight families.
_CODON_FAMILY: Dict[str, str] = {
    c: fam for fam, codons in SELECTIVE_FAMILIES.items() for c in codons
}


def to_dna(codon: str) -> str:
    """Normalise a codon label to upper-case DNA."""
    return str(codon).strip().upper().replace("U", "T")


def family_codons(families: Sequence[str]) -> List[str]:
    """All codons of the named amino-acid families, in family order."""
    return [c for fam in families for c in SELECTIVE_FAMILIES[fam]]


# --------------------------------------------------------------------------
# Published per-codon sensitivity (Elf et al. 2003, PMID 12805541)
# --------------------------------------------------------------------------

#: Measured starvation sensitivity per codon. Lower means the cognate tRNA
#: stays charged when its amino acid runs out. Values as tabulated by the user
#: in "Sensitive and insensitive codons table.xlsx"; codons absent from that
#: table have no published value here and are left unclassified by the
#: threshold rule.
ELF_SENSITIVITY: Dict[str, float] = {
    # insensitive end
    "TTG": 0.59, "TTA": 3.00, "ATA": 2.00, "ACG": 2.40, "AGT": 3.40,
    "CCT": 3.70, "AGA": 0.50, "AGG": 0.00, "GGA": 0.60, "GGG": 0.30,
    # sensitive end
    "CTA": 26.90, "CTT": 24.80, "CTC": 24.80, "ATT": 15.20, "ATC": 15.20,
    "ACC": 20.90, "TCC": 35.50, "CCA": 32.00, "CCG": 13.30, "CGA": 18.00,
    "CGT": 18.00, "CGC": 18.00, "GGC": 12.50, "GGT": 12.50,
}

#: Cut-offs applied to ELF_SENSITIVITY by the ``custom`` definition.
INSENSITIVE_MAX_SENSITIVITY = 4.0
SENSITIVE_MIN_SENSITIVITY = 12.0


def classify_by_sensitivity(
    insensitive_max: float = INSENSITIVE_MAX_SENSITIVITY,
    sensitive_min: float = SENSITIVE_MIN_SENSITIVITY,
):
    """Split ``ELF_SENSITIVITY`` into (insensitive, sensitive) by threshold.

    Returns two codon lists sorted by sensitivity. Codons whose value lies
    between the cut-offs belong to neither list.
    """
    ins = sorted((c for c, v in ELF_SENSITIVITY.items() if v < insensitive_max),
                 key=lambda c: ELF_SENSITIVITY[c])
    sen = sorted((c for c, v in ELF_SENSITIVITY.items() if v > sensitive_min),
                 key=lambda c: ELF_SENSITIVITY[c])
    return ins, sen


# --------------------------------------------------------------------------
# Codon lists as they stood before 2026-11-12
# --------------------------------------------------------------------------

#: All sixteen codons read by isoacceptors that stay charged (PMID 12805541).
PUBLISHED_INSENSITIVE_CODONS: List[str] = [
    "TTG", "TTA", "AGG", "AGA", "CGG", "ATA", "ACG", "AGT",
    "AGC", "GCG", "GCA", "GCT", "GGG", "GGA", "CCC", "CCT",
]

#: The GC-decoupled subset returned by the exhaustive search (eight codons).
SEARCH_OPTIMUM_INSENSITIVE_CODONS: List[str] = [
    "TTG", "TTA", "ATA", "ACG", "AGT", "GCA", "GCT", "CCT",
]

#: The eight-codon optimum plus the four rare Arg and Gly codons. This drove
#: PLOT 4 before 2026-11-12.
LEGACY_BAR_INSENSITIVE_CODONS: List[str] = [
    "TTG", "TTA", "ATA", "ACG", "AGT", "GCA", "GCT", "CCT",
    "AGA", "AGG", "GGA", "GGG",
]

#: The same set with ACG removed.
ACG_FREE_INSENSITIVE_CODONS: List[str] = [
    "TTG", "TTA", "ATA", "AGT", "GCA", "GCT", "CCT",
]

#: The insensitive set behind the per-gene ``perc_insensitive`` metric before
#: 2026-11-12. Differs from LEGACY_BAR_INSENSITIVE_CODONS in four codons.
LEGACY_METRIC_INSENSITIVE_CODONS: List[str] = [
    "AGA", "AGG", "CGG", "GGA", "GGG", "ATA", "TTA", "TTG",
    "AGC", "AGT", "TCG", "ACG",
]

#: The sensitive set behind the per-gene ``perc_sensitive`` metric before
#: 2026-11-12.
LEGACY_SENSITIVE_CODONS: List[str] = [
    "ACC", "ATC", "CCG", "ATT", "CTC", "GGC", "CTA", "CCA",
    "CGC", "CGA", "CGT", "GCC", "GTC", "CTT", "TCC", "CAA", "GTG", "GGT",
]


# --------------------------------------------------------------------------
# Definition objects
# --------------------------------------------------------------------------

class CodonDefinition:
    """One named way of deciding which codons are sensitive or insensitive.

    Attributes
    ----------
    key, label, description
        Identification, for config files and the GUI.
    metric_insensitive, metric_sensitive
        Lists behind the per-gene ``perc_insensitive`` / ``perc_sensitive``
        metrics in :mod:`codonpipe.fasta_metrics`.
    bar_insensitive
        List behind the per-cluster bars of PLOT 4.
    families
        Amino-acid families spanned by this definition. PLOT 4 uses their
        codons as its denominator.
    """

    __slots__ = ("key", "label", "description", "metric_insensitive",
                 "metric_sensitive", "bar_insensitive", "families")

    def __init__(self, key, label, description, metric_insensitive,
                 metric_sensitive, bar_insensitive, families):
        self.key = str(key)
        self.label = str(label)
        self.description = str(description)
        self.metric_insensitive = [to_dna(c) for c in metric_insensitive]
        self.metric_sensitive = [to_dna(c) for c in metric_sensitive]
        self.bar_insensitive = [to_dna(c) for c in bar_insensitive]
        self.families = list(families)

    @property
    def denominator_codons(self) -> List[str]:
        """Codons forming the PLOT 4 denominator for this definition."""
        return family_codons(self.families)

    def sensitivity_of(self, codon: str) -> Optional[float]:
        """Published Elf sensitivity for a codon, or None when untabulated."""
        return ELF_SENSITIVITY.get(to_dna(codon))

    def __repr__(self):
        return (f"CodonDefinition({self.key!r}, "
                f"{len(self.bar_insensitive)} insensitive, "
                f"{len(self.metric_sensitive)} sensitive, "
                f"{len(self.families)} families)")


_CUSTOM_INSENSITIVE, _CUSTOM_SENSITIVE = classify_by_sensitivity()

#: Families spanned by the custom definition. Ala is absent from the user's
#: table, so it is excluded from both numerator and denominator.
_CUSTOM_FAMILIES = sorted(
    {_CODON_FAMILY[c] for c in (_CUSTOM_INSENSITIVE + _CUSTOM_SENSITIVE)},
    key=list(SELECTIVE_FAMILIES).index,
)


DEFINITIONS: Dict[str, CodonDefinition] = {
    "legacy": CodonDefinition(
        key="legacy",
        label="Legacy (as used before 2026-11-12)",
        description=(
            "The sets the pipeline shipped with. Per-gene metrics and PLOT 4 "
            "use different twelve-codon insensitive lists; both are preserved "
            "so earlier analyses reproduce exactly."
        ),
        metric_insensitive=LEGACY_METRIC_INSENSITIVE_CODONS,
        metric_sensitive=LEGACY_SENSITIVE_CODONS,
        bar_insensitive=LEGACY_BAR_INSENSITIVE_CODONS,
        families=list(SELECTIVE_FAMILIES),
    ),
    "custom": CodonDefinition(
        key="custom",
        label="Custom (Elf sensitivity < 4 insensitive, > 12 sensitive)",
        description=(
            "Refined definition of 2026-11-12. Ten insensitive and fourteen "
            "sensitive codons classified purely by published Elf sensitivity. "
            "Alanine is excluded entirely; intermediate codons are "
            "unclassified. One list serves both the per-gene metric and PLOT 4."
        ),
        metric_insensitive=_CUSTOM_INSENSITIVE,
        metric_sensitive=_CUSTOM_SENSITIVE,
        bar_insensitive=_CUSTOM_INSENSITIVE,
        families=_CUSTOM_FAMILIES,
    ),
    "published": CodonDefinition(
        key="published",
        label="Published 16 (PMID 12805541, unrestricted)",
        description="All sixteen codons read by isoacceptors that stay charged.",
        metric_insensitive=PUBLISHED_INSENSITIVE_CODONS,
        metric_sensitive=LEGACY_SENSITIVE_CODONS,
        bar_insensitive=PUBLISHED_INSENSITIVE_CODONS,
        families=list(SELECTIVE_FAMILIES),
    ),
    "optimum": CodonDefinition(
        key="optimum",
        label="Subset-search optimum (eight codons)",
        description="The GC-decoupled eight-codon optimum of the exhaustive search.",
        metric_insensitive=SEARCH_OPTIMUM_INSENSITIVE_CODONS,
        metric_sensitive=LEGACY_SENSITIVE_CODONS,
        bar_insensitive=SEARCH_OPTIMUM_INSENSITIVE_CODONS,
        families=list(SELECTIVE_FAMILIES),
    ),
    "acgfree": CodonDefinition(
        key="acgfree",
        label="Subset-search optimum without ACG (seven codons)",
        description="The eight-codon optimum with ACG dropped.",
        metric_insensitive=ACG_FREE_INSENSITIVE_CODONS,
        metric_sensitive=LEGACY_SENSITIVE_CODONS,
        bar_insensitive=ACG_FREE_INSENSITIVE_CODONS,
        families=list(SELECTIVE_FAMILIES),
    ),
}

#: Keys offered in the GUI, in display order.
SELECTABLE_KEYS: List[str] = ["legacy", "custom", "published", "optimum", "acgfree"]

DEFAULT_DEFINITION_KEY = "legacy"

#: Environment variable consulted by :func:`active_definition_key`.
ENV_VAR = "CODONPIPE_CODON_DEFINITION"

_ALIASES = {
    "": DEFAULT_DEFINITION_KEY,
    "default": DEFAULT_DEFINITION_KEY,
    "before": "legacy",
    "old": "legacy",
    "previous": "legacy",
    "lt4": "legacy",
    "current": "custom",
    "refined": "custom",
    "new": "custom",
    "elf": "published",
    "elf16": "published",
    "search": "optimum",
    "acg_free": "acgfree",
    "acg-free": "acgfree",
}


def normalise_key(name: Optional[str]) -> Optional[str]:
    """Map a user-written definition name to a registry key, or None."""
    k = str(name or "").strip().lower().replace(" ", "")
    k = _ALIASES.get(k, k)
    return k if k in DEFINITIONS else None


def get_definition(name: Optional[str] = None) -> CodonDefinition:
    """Look up a definition by key or alias. Falls back to the default."""
    key = normalise_key(name)
    if key is None:
        if name:
            print(f"[WARN] Unknown codon definition {name!r}; "
                  f"using {DEFAULT_DEFINITION_KEY!r}. "
                  f"Known: {', '.join(SELECTABLE_KEYS)}.")
        key = DEFAULT_DEFINITION_KEY
    return DEFINITIONS[key]


def active_definition_key() -> str:
    """Definition key currently selected, from the environment."""
    return normalise_key(os.environ.get(ENV_VAR, "")) or DEFAULT_DEFINITION_KEY


def active_definition() -> CodonDefinition:
    """The currently selected definition."""
    return DEFINITIONS[active_definition_key()]


def set_active_definition(name: Optional[str]) -> CodonDefinition:
    """Select a definition for this process and its children.

    Also refreshes :data:`codonpipe.fasta_metrics.CODON_GROUPS` in place so
    per-gene metrics follow the selection.
    """
    definition = get_definition(name)
    os.environ[ENV_VAR] = definition.key
    try:
        from codonpipe import fasta_metrics
    except Exception:
        pass
    else:
        fasta_metrics.apply_codon_definition(definition)
    return definition


def describe(name: Optional[str] = None) -> str:
    """One-paragraph summary of a definition, for logs and figure footers."""
    d = get_definition(name)
    ins = " ".join(d.bar_insensitive)
    sen = " ".join(d.metric_sensitive)
    return (f"{d.label}\n  insensitive ({len(d.bar_insensitive)}): {ins}"
            f"\n  sensitive ({len(d.metric_sensitive)}): {sen}"
            f"\n  families ({len(d.families)}): {', '.join(d.families)}")
