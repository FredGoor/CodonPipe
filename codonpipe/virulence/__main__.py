#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Command-line entry point for automated virulence-gene discovery.

Examples
--------
Identify virulence genes for a named strain::

    python -m codonpipe.virulence "Salmonella Typhimurium SL1344" -o out

Use an assembly accession and a coding-sequence FASTA you already have, so the
output locus tags match the file CodonPipe will analyse::

    python -m codonpipe.virulence GCF_000210855.2 \\
        --cds "Preloaded genomes/Salmonella enterica SL1344.fna" \\
        --stringency standard -o out

Reproduce a previous run with no network access::

    python -m codonpipe.virulence GCF_000210855.2 -o out --offline

Run the concordance benchmark against curated reference sets::

    python -m codonpipe.virulence --benchmark "Preloaded genomes" -o out
"""

from __future__ import annotations

import argparse
import os
import sys

from .classify import canonical_tier


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m codonpipe.virulence",
        description="Identify bacterial virulence genes and split them into "
        "virulence-associated and host-interaction sets.",
    )
    parser.add_argument(
        "target",
        nargs="?",
        default="",
        help="Species or strain name, NCBI taxid (216597 / txid:216597), or "
        "assembly accession (GCF_000210855.2).",
    )
    parser.add_argument("-o", "--output-dir", default="CodonPipe virulence",
                        help="Folder for the four output artefacts.")
    parser.add_argument("--cds", default="",
                        help="Coding-sequence FASTA defining the gene universe. "
                             "Omit to download the RefSeq CDS FASTA.")
    parser.add_argument("--assembly", default="",
                        help="Pin the assembly explicitly.")
    parser.add_argument("--stringency", default="standard",
                        choices=["permissive", "standard", "strict", "host-interaction"],
                        help="Headline tier. All three are always computed. "
                             "'strict' and 'host-interaction' are two names for "
                             "the same tier.")
    parser.add_argument("--tiers", default="all",
                        help="Tiers to write as cluster columns: 'all' (default) or a "
                             "comma-separated subset of permissive,standard,strict.")
    parser.add_argument("--cache-dir", default="",
                        help="Cache folder; defaults to <output-dir>/virulence_cache.")
    parser.add_argument("--offline", action="store_true",
                        help="Replay from the cache without any network access.")
    parser.add_argument("--refresh", action="store_true",
                        help="Re-download every record even if cached.")
    parser.add_argument("--no-species-transfer", action="store_true",
                        help="Do not transfer reviewed UniProtKB curation from the "
                             "parent species by gene symbol.")
    parser.add_argument("--essentiality", default="",
                        help="Experimental gene-requirement table (locus_tag, "
                             "in_vitro_requirement, in_vivo_fitness). Overrides "
                             "the bundled reference for the organism.")

    parser.add_argument("--effect-contributing", type=float, default=None,
                        help="log2 depletion at or above which an in vivo defect counts as "
                             "contributing (rung 4). Default 0.0: any positive measured "
                             "defect. On the percentile scale, a percentile.")
    parser.add_argument("--effect-required", type=float, default=None,
                        help="log2 depletion at or above which the gene counts as required "
                             "in vivo (rung 5). Default 1.0, a two-fold defect, which is the "
                             "threshold the bundled Salmonella compilation applies to its own "
                             "competitive indices. On the percentile scale, a percentile.")
    parser.add_argument("--max-in-vivo-significance", type=float, default=None,
                        help="largest p or q value accepted where a screen reports one. "
                             "Default 0.05.")
    parser.add_argument("--effect-scale", choices=["absolute", "percentile"], default=None,
                        help="read the two effect cut-offs as absolute log2 depletion "
                             "(default) or as percentiles of this organism's measured "
                             "in vivo effect sizes.")
    parser.add_argument("--no-bundled-reference", action="store_true",
                        help="Do not load the bundled experimental reference table "
                             "even when one matches the resolved organism.")
    parser.add_argument("--ncbi-email", default="",
                        help="Contact address sent to NCBI Entrez, which asks "
                             "automated clients to identify themselves.")
    parser.add_argument("--resource-dir", default="",
                        help="Alternative folder of classification rule tables.")
    parser.add_argument("--prefix", default="",
                        help="Prefix for the output file names.")
    parser.add_argument("--benchmark", metavar="GENOMES_DIR", default="",
                        help="Run the concordance benchmark using the curated "
                             "workbooks in this 'Preloaded genomes' folder.")
    parser.add_argument("--quiet", action="store_true")
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    log = (lambda _m: None) if args.quiet else (lambda m: print(m, flush=True))

    if args.benchmark:
        from .benchmark import default_cases, run_benchmark, write_benchmark_workbook

        cases = default_cases(args.benchmark)
        if not cases:
            print(
                f"No benchmark case found in {args.benchmark!r}: expected the CodonPipe "
                "genome FASTA files and their curated cluster workbooks side by side.",
                file=sys.stderr,
            )
            return 2
        summary, details = run_benchmark(cases, cache_dir=args.cache_dir, log=log)
        os.makedirs(args.output_dir, exist_ok=True)
        path = write_benchmark_workbook(
            os.path.join(args.output_dir, "CodonPipe virulence discovery benchmark.xlsx"),
            summary,
            details,
        )
        log(f"benchmark written to {path}")
        return 0

    if not args.target:
        build_parser().print_help()
        return 2

    from .pipeline import DiscoveryConfig, run_discovery

    # Only cut-offs the caller actually typed are passed through, so an
    # unspecified flag leaves the published default in place rather than
    # overwriting it with None.
    overrides = {
        key: value
        for key, value in (
            ("effect_contributing_log2", args.effect_contributing),
            ("effect_required_log2", args.effect_required),
            ("max_in_vivo_significance", args.max_in_vivo_significance),
            ("effect_scale", args.effect_scale),
        )
        if value is not None
    }

    config = DiscoveryConfig(
        target=args.target,
        assembly=args.assembly,
        cds_fasta=args.cds,
        stringency=canonical_tier(args.stringency),
        tiers=args.tiers,
        output_dir=args.output_dir,
        cache_dir=args.cache_dir,
        offline=args.offline,
        refresh=args.refresh,
        include_species_transfer=not args.no_species_transfer,
        essentiality_table=args.essentiality,
        use_bundled_reference=not args.no_bundled_reference,
        ncbi_email=args.ncbi_email,
        resource_dir=args.resource_dir,
        output_prefix=args.prefix,
        scoring_overrides=overrides,
    )
    result = run_discovery(config, log=log, write_outputs=True)
    for name, path in result.outputs.items():
        log(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
