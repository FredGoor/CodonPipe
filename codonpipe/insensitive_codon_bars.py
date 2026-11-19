# codonpipe/insensitive_codon_bars.py
"""Starvation-insensitive codon use per cluster (PLOT 4).

One horizontal bar per cluster: the cluster's use of starvation-insensitive
codons minus the whole-genome value, in percentage points.

Definition. Elf, Nilsson & Ehrenberg (PMID 12805541) identified eight amino-acid
families -- Ala, Gly, Ile, Leu, Pro, Arg, Ser, Thr -- whose tRNA isoacceptors
are charged selectively during amino-acid starvation, and 16 codons within them
read by the isoacceptors that stay charged. The quantity plotted here is

    (codons from the insensitive set) / (all codons of those eight families)

computed per cluster from absolute codon usage, minus the same ratio for the
whole genome. The denominator is the eight families' total, NOT all 61 sense
codons, so the measure is insensitive to amino-acid composition and reports
synonymous choice only.

Codon set. Which codons count is decided by the definition selected in
:mod:`codonpipe.codon_sensitivity`, not here; this module only draws what that
registry hands it. Pass ``definition="legacy"`` or ``definition="custom"`` to a
drawing call to override the active selection for one figure.

Under ``legacy`` the set is DEFAULT_INSENSITIVE_CODONS, the twelve-codon
reported set: an eight-codon search optimum plus AGA, AGG, GGA and GGG. An
exhaustive screen of all 65,535 non-empty subsets against curated virulence
genes in S. Typhimurium SL1344 and E. coli O157:H7 Sakai showed that a codon's
individual contribution is predicted by its correlation with gene GC content:
the four lowest-scoring codons (GCG, CGG, AGC, ACG) are the four most
positively GC-correlated, and the two that help most (TTA, ATA) are the most
negatively correlated. Restricting the set removes the route by which gene base
composition, rather than translational physiology, enters the measurement. The
twelve-codon set excludes three of those four but RETAINS ACG, which the
exhaustive search keeps because it raises the score of the particular subset it
sits in even though its marginal contribution averaged over all subsets is
slightly negative. ACG_FREE_INSENSITIVE_CODONS drops it: that costs about 8% of
the match to curated virulence genes and slightly improves their enrichment
against GC-matched controls, so the choice changes no conclusion. The four
added Arg and Gly codons are not GC-confounded, only uninformative for the
virulence discrimination, so including them broadens the claim without changing
any conclusion. SEARCH_OPTIMUM_INSENSITIVE_CODONS is the eight-codon optimum
and PUBLISHED_INSENSITIVE_CODONS the unrestricted 16.

Under ``custom`` the set is the ten codons whose published Elf sensitivity is
below 4, and the denominator narrows from eight amino-acid families to seven:
alanine carries no tabulated sensitivity in that definition and is excluded
from numerator and denominator alike. Percentages are therefore NOT comparable
between definitions, only within one.

the eight-codon optimum and PUBLISHED_INSENSITIVE_CODONS the unrestricted 16.
"""

from __future__ import annotations

import os
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


from codonpipe.codon_sensitivity import (  # noqa: F401  (re-exported)
    SELECTIVE_FAMILIES,
    PUBLISHED_INSENSITIVE_CODONS,
    SEARCH_OPTIMUM_INSENSITIVE_CODONS,
    ACG_FREE_INSENSITIVE_CODONS,
    LEGACY_BAR_INSENSITIVE_CODONS,
    active_definition,
    get_definition,
    family_codons,
)

#: Codons of all eight selective families. Kept for callers that pass an
#: explicit codon list and expect the historical eight-family denominator.
FAMILY_CODONS: List[str] = family_codons(list(SELECTIVE_FAMILIES))

#: Backwards-compatible alias for the twelve-codon legacy bar set.
DEFAULT_INSENSITIVE_CODONS: List[str] = LEGACY_BAR_INSENSITIVE_CODONS


DEFAULT_BAR_COLOR = "0.62"
DEFAULT_HIGHLIGHT_COLOR = "#b2182b"


def _to_dna(codon: str) -> str:
    return str(codon).strip().upper().replace("U", "T")


def _codon_key(label) -> str:
    """Pull a bare DNA codon out of whatever the workbook row label looks like.

    Codon-usage sheets label rows variously 'GCT', 'Ala_GCU', 'A (GCU)'. Take
    the last run of three A/C/G/T/U characters and normalise it.
    """
    s = str(label).strip().upper().replace("U", "T")
    best = ""
    run = ""
    for ch in s:
        if ch in "ACGT":
            run += ch
            if len(run) >= 3:
                best = run[-3:]
        else:
            run = ""
    return best


def _resolve(insensitive=None, definition=None):
    """Return (codon list, denominator codon list) for a drawing call.

    ``definition`` names a registry entry; ``insensitive`` overrides its codon
    list. When only an explicit codon list is given the denominator stays the
    historical eight-family total, which is what callers passing a hand-written
    list have always meant.
    """
    if definition is not None:
        d = definition if hasattr(definition, "bar_insensitive") else get_definition(definition)
    elif insensitive is None:
        d = active_definition()
    else:
        d = None
    if d is None:
        return [_to_dna(c) for c in insensitive], list(FAMILY_CODONS)
    codons = [_to_dna(c) for c in (insensitive if insensitive is not None else d.bar_insensitive)]
    return codons, list(d.denominator_codons)


def insensitive_fraction(counts: pd.Series,
                         insensitive: Sequence[str],
                         denominator: Optional[Sequence[str]] = None) -> float:
    """Insensitive-codon share of the selective families, as a percentage.

    ``denominator`` defaults to all codons of the eight selective families.
    Definitions spanning fewer families (``custom`` excludes alanine) pass
    their own narrower list.
    """
    ins = {_to_dna(c) for c in insensitive}
    den_codons = [_to_dna(c) for c in (denominator if denominator is not None else FAMILY_CODONS)]
    num = float(sum(counts.get(c, 0.0) for c in ins))
    den = float(sum(counts.get(c, 0.0) for c in den_codons))
    if den <= 0:
        return float("nan")
    return 100.0 * num / den


def cluster_values(
    usage_df: pd.DataFrame,
    genome_columns: Sequence[str],
    insensitive: Optional[Sequence[str]] = None,
    definition=None,
) -> pd.Series:
    """Per-column insensitive-codon percentage from an absolute-usage table.

    `usage_df` is a codon-usage sheet: first column codon labels, remaining
    columns one per cluster plus one or more whole-genome columns.

    With neither `insensitive` nor `definition` given, the definition selected
    in codonpipe.codon_sensitivity applies.
    """
    insensitive, denominator = _resolve(insensitive, definition)
    frame = usage_df.copy()
    frame.index = [_codon_key(v) for v in frame.iloc[:, 0].tolist()]
    frame = frame.iloc[:, 1:]
    frame = frame[~frame.index.duplicated(keep="first")]
    frame = frame.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    out = {}
    for col in frame.columns:
        out[str(col)] = insensitive_fraction(frame[col], insensitive, denominator)
    return pd.Series(out, dtype=float)


def plot_insensitive_codon_bars(
    usage_df: pd.DataFrame,
    output_file: str,
    genome_columns: Sequence[str] = (),
    insensitive: Optional[Sequence[str]] = None,
    definition=None,
    include_columns: Optional[Sequence[str]] = None,
    highlight_columns: Optional[Sequence[str]] = None,
    cluster_sizes: Optional[Dict[str, int]] = None,
    fig_w_in: float = 7.0,
    height_per_cluster: float = 0.26,
    label_fontsize: int = 10,
    dpi: int = 300,
    show: bool = False,
    title: str = "Starvation-insensitive codon use by cluster",
) -> Optional[str]:
    """Draw PLOT 4 and return the path written, or None when there is nothing to draw.

    `definition` selects a codon definition by key ("legacy", "custom", ...);
    `insensitive` overrides its codon list. Omit both to follow the active
    selection.
    """
    insensitive, denominator = _resolve(insensitive, definition)
    vals = cluster_values(usage_df, genome_columns, insensitive, definition)
    gnames = {str(g).strip().lower() for g in genome_columns}
    gnames |= {"whole genome", "genome", "all genes", "all"}

    genome_vals = [v for k, v in vals.items() if k.strip().lower() in gnames and np.isfinite(v)]
    if not genome_vals:
        print("[WARN] PLOT 4: no whole-genome column found in the codon-usage table; "
              "using the mean over clusters as the reference.")
        base = float(np.nanmean(vals.values))
    else:
        base = float(np.mean(genome_vals))

    clusters = {k: v for k, v in vals.items()
                if k.strip().lower() not in gnames and np.isfinite(v)}
    if include_columns:
        want = {str(c).strip().lower() for c in include_columns}
        clusters = {k: v for k, v in clusters.items() if k.strip().lower() in want}
    if not clusters:
        print("[WARN] PLOT 4: no cluster columns left to plot; skipping.")
        return None

    data = pd.Series(clusters).sort_values()
    delta = data - base
    hi = {str(c).strip().lower() for c in (highlight_columns or [])}
    colors = [DEFAULT_HIGHLIGHT_COLOR if k.strip().lower() in hi else DEFAULT_BAR_COLOR
              for k in delta.index]

    labels = []
    for k in delta.index:
        n = (cluster_sizes or {}).get(k)
        labels.append(f"{k}  ({n})" if n else str(k))

    fig_h = max(2.4, 1.15 + height_per_cluster * len(delta))
    fig = plt.figure(figsize=(fig_w_in, fig_h))
    ax = fig.add_subplot(111)
    y = np.arange(len(delta))
    ax.barh(y, delta.values, color=colors, height=0.74, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=label_fontsize)
    ax.axvline(0.0, color="0.25", lw=0.9, zorder=4)
    ax.grid(axis="x", alpha=0.25, zorder=0)
    ax.tick_params(axis="x", labelsize=label_fontsize - 1)
    ax.set_xlabel("insensitive-codon use above genome mean (percentage points)",
                  fontsize=label_fontsize)
    ax.set_title(title, loc="left", fontsize=label_fontsize + 2)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    span = float(np.nanmax(np.abs(delta.values))) or 1.0
    for yy, v in zip(y, delta.values):
        ax.text(v + (0.03 * span if v >= 0 else -0.03 * span), yy, f"{v:+.1f}",
                va="center", ha="left" if v >= 0 else "right",
                fontsize=max(5.0, label_fontsize - 3), color="0.3")
    ax.set_xlim(min(delta.min() * 1.35, -0.5), max(delta.max() * 1.35, 0.5))
    ax.margins(y=0.01)

    codon_note = " ".join(_to_dna(c) for c in insensitive)
    den_set = {_to_dna(c) for c in denominator}
    fam_note = ", ".join(f for f, cs in SELECTIVE_FAMILIES.items()
                         if den_set.issuperset(_to_dna(c) for c in cs))
    fig.text(0.0, 0.004,
             f"Insensitive codons ({len(insensitive)}): {codon_note}.  "
             f"Denominator: all codons of {fam_note}.  Genome value {base:.1f}%.",
             fontsize=max(5.0, label_fontsize - 3.5), color="0.35", ha="left", va="bottom")

    folder = os.path.dirname(os.path.abspath(output_file))
    if folder:
        os.makedirs(folder, exist_ok=True)
    fig.savefig(output_file, dpi=int(dpi), bbox_inches="tight")
    print(f"[INFO] Insensitive-codon bar plot written: {output_file}")
    if show:
        plt.show()
    else:
        plt.close(fig)
    return output_file
