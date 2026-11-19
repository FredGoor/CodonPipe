#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Genomic-context extension of virulence evidence.

Curated virulence-factor databases are built from individual characterised
proteins, so a database screen typically hits the conspicuous members of a
locus -- the secretion-system ATPase, the named effector, the major fimbrial
subunit -- and misses the chaperones, minor subunits and uncharacterised open
reading frames transcribed with them.  Manual curation closes that gap by
walking the annotation around each hit, which is slow and hard to report
reproducibly.

Two automated, parameterised layers replace that step:

**Operon extension** propagates evidence to co-directional neighbours
separated by short intergenic distances.  Bacterial operons are defined by
exactly this signature, and the default 50 bp threshold is the conventional
cut-off for co-transcription in gamma-proteobacteria.

**Island bridging** merges runs of evidence-positive genes that are separated
by a small number of evidence-negative genes, recovering complete
pathogenicity islands such as SPI-1 and SPI-2 whose internal genes are not all
individually characterised.

Both layers are reported separately from direct database evidence, so a
reviewer can see exactly which genes entered a set by context rather than by
curation, and both can be switched off.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import pandas as pd

from .resolve import GeneTable

__all__ = ["ContextResult", "extend_by_context"]


@dataclass
class ContextResult:
    """Genes added by context, with the reason for each addition."""

    added: Dict[str, dict] = field(default_factory=dict)
    notes: List[str] = field(default_factory=list)

    def as_frame(self) -> pd.DataFrame:
        if not self.added:
            return pd.DataFrame(
                columns=["locus_tag", "layer", "seed_locus_tags", "block_size", "detail"]
            )
        rows = []
        for tag, info in self.added.items():
            rows.append(
                {
                    "locus_tag": tag,
                    "layer": info.get("layer"),
                    "seed_locus_tags": "; ".join(info.get("seeds", [])[:8]),
                    "block_size": info.get("block_size"),
                    "detail": info.get("detail", ""),
                }
            )
        return pd.DataFrame(rows).sort_values("locus_tag").reset_index(drop=True)


def _ordered_by_replicon(gene_table: GeneTable) -> Dict[str, List[Tuple[str, float, float, str]]]:
    frame = gene_table.frame
    if frame.empty or frame["start"].isna().all():
        return {}
    usable = frame.dropna(subset=["start", "end"]).copy()
    usable["seqid"] = usable["seqid"].fillna("unknown").astype(str)
    out: Dict[str, List[Tuple[str, float, float, str]]] = {}
    for seqid, group in usable.sort_values(["seqid", "start"]).groupby("seqid"):
        out[str(seqid)] = [
            (str(r.locus_tag), float(r.start), float(r.end), str(r.strand or "?"))
            for r in group.itertuples(index=False)
        ]
    return out


def extend_by_context(
    gene_table: GeneTable,
    seed_tags: Iterable[str],
    operon_extension: bool = True,
    max_intergenic_gap: int = 50,
    require_same_strand: bool = True,
    island_bridging: bool = True,
    island_bridge_genes: int = 2,
    min_island_seeds: int = 3,
) -> ContextResult:
    """Propagate virulence evidence into operons and pathogenicity islands.

    Parameters
    ----------
    seed_tags:
        Locus tags carrying direct database evidence.
    operon_extension, max_intergenic_gap, require_same_strand:
        Co-directional neighbours whose intergenic distance does not exceed
        ``max_intergenic_gap`` base pairs are treated as co-transcribed with a
        seed.  Overlapping genes (negative distance) always qualify.
    island_bridging, island_bridge_genes, min_island_seeds:
        Seed genes separated by at most ``island_bridge_genes``
        evidence-negative genes are merged into one block; a block is accepted
        as an island only if it contains at least ``min_island_seeds`` seeds.

    Returns
    -------
    ContextResult
        Only genes that were *not* already seeds appear in ``added``.
    """
    result = ContextResult()
    seeds: Set[str] = {str(t) for t in seed_tags if t}
    if not seeds:
        result.notes.append("No direct evidence was available, so no context extension was attempted.")
        return result

    per_replicon = _ordered_by_replicon(gene_table)
    if not per_replicon:
        result.notes.append(
            "Genomic coordinates are unavailable for this genome, so operon and island "
            "extension were skipped."
        )
        return result

    # ---------------------------------------------------------------- operons
    if operon_extension:
        n_before = len(result.added)
        for seqid, genes in per_replicon.items():
            n = len(genes)
            blocks: List[List[int]] = []
            current = [0]
            for i in range(1, n):
                prev_tag, prev_start, prev_end, prev_strand = genes[i - 1]
                tag, start, end, strand = genes[i]
                gap = start - prev_end - 1
                same_strand = (strand == prev_strand) or not require_same_strand
                if same_strand and gap <= float(max_intergenic_gap):
                    current.append(i)
                else:
                    blocks.append(current)
                    current = [i]
            blocks.append(current)
            for block in blocks:
                tags = [genes[i][0] for i in block]
                block_seeds = [t for t in tags if t in seeds]
                if not block_seeds:
                    continue
                for tag in tags:
                    if tag in seeds or tag in result.added:
                        continue
                    result.added[tag] = {
                        "layer": "operon",
                        "seeds": block_seeds,
                        "block_size": len(tags),
                        "detail": (
                            f"co-directional with {len(block_seeds)} evidence-positive "
                            f"gene(s) within {max_intergenic_gap} bp on {seqid}"
                        ),
                    }
        result.notes.append(
            f"Operon extension added {len(result.added) - n_before} gene(s) "
            f"(co-directional, intergenic gap <= {max_intergenic_gap} bp)."
        )

    # ---------------------------------------------------------------- islands
    if island_bridging:
        n_before = len(result.added)
        for seqid, genes in per_replicon.items():
            tags = [g[0] for g in genes]
            seed_idx = [i for i, t in enumerate(tags) if t in seeds]
            if len(seed_idx) < max(1, int(min_island_seeds)):
                continue
            groups: List[List[int]] = []
            current = [seed_idx[0]]
            for prev, cur in zip(seed_idx, seed_idx[1:]):
                if (cur - prev - 1) <= int(island_bridge_genes):
                    current.append(cur)
                else:
                    groups.append(current)
                    current = [cur]
            groups.append(current)
            for group in groups:
                if len(group) < int(min_island_seeds):
                    continue
                lo, hi = group[0], group[-1]
                block_seeds = [tags[i] for i in group]
                for i in range(lo, hi + 1):
                    tag = tags[i]
                    if tag in seeds or tag in result.added:
                        continue
                    result.added[tag] = {
                        "layer": "island",
                        "seeds": block_seeds,
                        "block_size": hi - lo + 1,
                        "detail": (
                            f"inside a run of {len(group)} evidence-positive genes "
                            f"spanning {hi - lo + 1} consecutive genes on {seqid}"
                        ),
                    }
        result.notes.append(
            f"Island bridging added {len(result.added) - n_before} gene(s) "
            f"(runs of >= {min_island_seeds} evidence-positive genes, bridging up to "
            f"{island_bridge_genes} negative gene(s))."
        )

    return result
