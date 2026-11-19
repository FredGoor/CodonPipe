"""Experimentally measured gene-requirement tables.

Every other evidence layer in this package infers a gene's role from curation
or from annotation.  None of them can decide Falkow's third postulate, which
asks for a *measured* loss of virulence when the gene is inactivated, nor can
they decide whether such a measurement is even possible -- a gene that cannot
be deleted without killing the bacterium in broth can never produce an
interpretable loss-of-virulence phenotype.

This module supplies that missing layer for organisms where published
mutant screens exist.  A reference table is a per-gene record of two
independent experimental facts:

``in_vitro_requirement``
    ``essential``    the gene cannot be inactivated without a growth defect
                     in rich medium, so postulate 2 is not interpretable;
    ``dispensable``  a viable insertion or deletion mutant has been made and
                     assayed, which settles the question in the other
                     direction regardless of what the annotation suggests;
    ``unknown``      not determined.

``in_vivo_fitness``
    ``attenuated``   the mutant is under-represented in at least one animal
                     infection screen -- direct evidence for postulate 3;
    ``no_phenotype`` the mutant was assayed in vivo and was *not*
                     under-represented;
    ``not_assayed``  the gene was never covered by a screen.

The distinction between ``no_phenotype`` and ``not_assayed`` is the point of
the three-state encoding.  Screens cover part of a genome, and effector
families are functionally redundant, so a gene without an attenuation
phenotype is not thereby shown to be irrelevant.  Collapsing the two states
would silently convert "nobody looked" into "we looked and found nothing".

A reference table *overrides* the inferred layers wherever the two disagree,
because a measurement beats an annotation.  Tables are organism-scoped and
are matched on taxon id or assembly accession, never applied across species.

Users can supply their own table for any organism through the GUI panel or
``--essentiality`` on the command line; the file format is the same as the
bundled tables, and a user-supplied table always wins over a bundled one.
"""

from __future__ import annotations

import hashlib
import os
import warnings
from typing import Dict, List, Optional

import numpy as np

import pandas as pd

RESOURCE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources")

#: Vocabularies written into the tables.  Kept here so that a typo in a
#: hand-edited table is caught on load rather than silently ignored.
IN_VITRO_VALUES = ("essential", "fitness_cost", "dispensable", "unknown")
IN_VIVO_VALUES = ("attenuated", "no_phenotype", "not_assayed")

#: The infection model an in vivo call came from, with the rank used to settle
#: disagreements between screens.
#:
#: "In vivo" is not a synonym for "in an animal". A macrophage-survival screen,
#: an enteroid infection and a mouse systemic infection all report on a gene's
#: behaviour inside a host, and they do not report the same thing -- SPI-2 is
#: the obvious case, where the magnitude of the requirement differs sharply
#: between a cultured macrophage and a mouse. Recording which model produced a
#: call keeps that difference visible instead of flattening it.
#:
#: Higher rank wins. Where a gene carries calls from more than one model class,
#: the highest-ranked class decides, and the lower-ranked results are kept in
#: ``in_vivo_basis`` rather than discarded, so an override can be audited. The
#: ordering is a judgement about proximity to natural infection, and it is
#: stated here rather than buried in the builder so that it can be argued with.
#: ``animal_infection_enteric`` is a PEER of ``animal_infection``, not a
#: lower class: an oral-challenge colonisation screen is as natural an
#: infection as an intravenous one, it simply interrogates a different route.
#: Equal rank means the two cannot outrank each other, so a gene attenuated
#: by either is attenuated, and a conflict between them is resolved by the
#: within-class rule (attenuated beats no_phenotype) and reported in
#: ``in_vivo_basis``. Keeping the route in the model name is what lets a
#: reader see that SPI-1's requirement comes from the enteric screens and
#: SPI-2's from the systemic ones, which is the real biology.
IN_VIVO_MODEL_RANK = {
    "animal_infection": 4,
    "animal_infection_enteric": 4,
    "organoid": 3,
    "ex_vivo": 3,
    "cell_infection": 2,
    "non_mammalian_host": 1,
    "not_applicable": 0,
    "": 0,
}
IN_VIVO_MODELS = tuple(k for k in IN_VIVO_MODEL_RANK if k)

#: Bundled reference tables.  ``match`` lists the identifiers a resolved
#: genome may carry; a table is used only when one of them matches.
REFERENCE_INDEX: List[dict] = [
    {
        "key": "salmonella_typhimurium_sl1344",
        "organism": "Salmonella enterica subsp. enterica serovar Typhimurium SL1344",
        "table": "reference_salmonella_sl1344.tsv",
        "id_column": "locus_tag",
        # Strain-scoped, like every other entry. The species and genus ids
        # were previously listed here, so any S. enterica serovar -- Typhi
        # CT18, Enteritidis P125109, Typhimurium LT2 -- inherited a
        # compilation keyed on SL1344 locus tags, which would silently match
        # almost nothing and report the resulting emptiness as biology. The
        # other four tables already carry an explicit warning against exactly
        # this; this one did not follow its own rule.
        "match": {
            "taxon_id": ["216597"],
            "assembly": ["GCF_000210855.2", "GCA_000210855.2"],
        },
        "description": (
            "Curated compilation of published Salmonella Typhimurium mutant "
            "screens, keyed on SL1344 locus tags. SL1344 only: other serovars "
            "and other S. Typhimurium strains do not share these locus tags."
        ),
    },
    {
        "key": "shigella_flexneri_2a_301",
        "organism": "Shigella flexneri 2a str. 301",
        "table": "reference_shigella_flexneri_2a_301.tsv",
        "id_column": "locus_tag",
        "match": {
            "taxon_id": ["198214"],
            "species_taxon_id": ["623"],
            "assembly": ["GCF_000006925.2", "GCA_000006925.2"],
        },
        "description": (
            "In vitro requirement from the Keio/PEC orthologue calls and Keio "
            "growth yields carried by the S. flexneri 2a 2457T Tn-seq "
            "compilation; in vivo fitness from the M90T human gut organoid "
            "colonisation screen. Both axes are orthologue-transferred, which "
            "the per-gene basis text states."
        ),
    },
    {
        "key": "ecoli_o157h7_sakai",
        "organism": "Escherichia coli O157:H7 str. Sakai",
        "table": "reference_ecoli_o157h7_sakai.tsv",
        "id_column": "locus_tag",
        # Deliberately not matched on the E. coli species taxon: the EHEC
        # screens behind this table do not describe UPEC or K-12, and a
        # species-level match would silently apply them to both.
        "match": {
            "taxon_id": ["386585"],
            "assembly": ["GCF_000008865.2", "GCA_000008865.2"],
        },
        "description": (
            "Transposon-insertion sequencing of EHEC O157:H7 EDL933 in vitro "
            "and in the infant rabbit colon, transferred to Sakai by 1:1 "
            "orthologue family, plus Keio growth yields for graded in vitro "
            "fitness."
        ),
    },
    {
        "key": "ecoli_upec_uti89",
        "organism": "Escherichia coli UTI89 (UPEC)",
        "table": "reference_ecoli_upec_uti89.tsv",
        "id_column": "locus_tag",
        "match": {
            "taxon_id": ["364106"],
            "assembly": ["GCF_000013265.1", "GCA_000013265.1"],
        },
        "description": (
            "TraDIS measured in UTI89 itself: essentiality, graded fitness in "
            "LB and in minimal medium, growth in human urine, three murine "
            "urinary tract infections, and a bladder cell infection cycle "
            "screen. The only one of the four tables whose every call is a "
            "measurement in the target strain."
        ),
    },
    {
        "key": "klebsiella_pneumoniae_kppr1",
        "organism": "Klebsiella pneumoniae subsp. pneumoniae ATCC 43816 KPPR1",
        "table": "reference_klebsiella_pneumoniae_kppr1.tsv",
        "id_column": "locus_tag",
        # Not matched on the K. pneumoniae species taxon (573): HS11286 and
        # NTUH-K2044 are different strains and must not inherit these calls.
        "match": {
            "taxon_id": ["72407"],
            "assembly": ["GCF_000742755.1", "GCA_000742755.1"],
        },
        "description": (
            "Five in vivo screens keyed on VK055 locus tags -- murine "
            "pneumonia, bacteremia and gut colonisation, Galleria mellonella, "
            "and ECL8 urine/serum transferred through the publishing authors' "
            "own homology table -- with TraDIS essentiality and graded LB "
            "fitness measured in ATCC 43816."
        ),
    },
]


def available_references() -> List[dict]:
    """Return the bundled reference descriptors that have a readable table."""

    out = []
    for entry in REFERENCE_INDEX:
        path = os.path.join(RESOURCE_DIR, entry["table"])
        if os.path.exists(path):
            out.append(dict(entry, path=path))
    return out


def find_reference(
    taxon_id: str = "",
    species_taxon_id: str = "",
    assembly: str = "",
    organism: str = "",
) -> Optional[dict]:
    """Find a bundled reference table for a resolved genome.

    Matching is by identifier, not by name similarity: a reference measured in
    one strain is only offered for that strain's own taxon or assembly, or for
    the species it belongs to.  ``organism`` is used only as a last resort and
    only on an exact, case-folded match of the recorded organism string.
    """

    wanted = {
        "taxon_id": str(taxon_id or "").strip(),
        "species_taxon_id": str(species_taxon_id or "").strip(),
        "assembly": str(assembly or "").strip(),
    }
    for entry in available_references():
        rules = entry.get("match", {})
        for field, value in wanted.items():
            if value and value in rules.get(field, []):
                return entry
        if organism and organism.strip().lower() == entry["organism"].lower():
            return entry
    return None


def load_reference(entry_or_key) -> pd.DataFrame:
    """Load a bundled reference table and validate its vocabularies."""

    if isinstance(entry_or_key, str):
        matches = [e for e in available_references() if e["key"] == entry_or_key]
        if not matches:
            raise KeyError(f"No bundled reference table named {entry_or_key!r}.")
        entry = matches[0]
    else:
        entry = entry_or_key
    frame = pd.read_csv(entry["path"], sep="\t", dtype=str).fillna("")
    return validate_reference(frame, source=entry["path"])


def validate_reference(frame: pd.DataFrame, source: str = "") -> pd.DataFrame:
    """Check a reference table's required columns and controlled vocabularies.

    Raises on an unrecognised category rather than dropping it, so that a
    mistyped value cannot quietly remove genes from the host-interaction tier.
    """

    cols = {c.lower(): c for c in frame.columns}
    if "locus_tag" not in cols:
        raise ValueError(f"Reference table {source or '<frame>'} has no 'locus_tag' column.")
    for name, allowed in (("in_vitro_requirement", IN_VITRO_VALUES),
                          ("in_vivo_fitness", IN_VIVO_VALUES),
                          ("in_vivo_model", IN_VIVO_MODELS)):
        col = cols.get(name)
        if not col:
            continue
        bad = sorted({
            str(v).strip().lower()
            for v in frame[col].unique()
            if str(v).strip() and str(v).strip().lower() not in allowed
        })
        if bad:
            raise ValueError(
                f"Reference table {source or '<frame>'} column {name!r} contains "
                f"unrecognised value(s) {bad}; allowed: {list(allowed)}."
            )
    return frame


def resolve_in_vivo_call(records: "list[dict]") -> Dict[str, str]:
    """Settle an in vivo call from several screens, by model precedence.

    ``records`` is a list of ``{"fitness": ..., "model": ..., "basis": ...}``
    for one gene. The highest-ranked model class decides; within that class an
    ``attenuated`` result wins over ``no_phenotype``, because a screen that
    detected a defect is positive evidence while one that did not may simply
    have lacked power. Every contributing screen, including overridden ones,
    is listed in the returned basis.
    """

    usable = [r for r in records if str(r.get("fitness") or "") in ("attenuated", "no_phenotype")]
    if not usable:
        return {"in_vivo_fitness": "not_assayed", "in_vivo_model": "not_applicable",
                "in_vivo_basis": ""}
    best = max(IN_VIVO_MODEL_RANK.get(str(r.get("model") or ""), 0) for r in usable)
    top = [r for r in usable if IN_VIVO_MODEL_RANK.get(str(r.get("model") or ""), 0) == best]
    fitness = "attenuated" if any(r["fitness"] == "attenuated" for r in top) else "no_phenotype"
    model = next(r.get("model") for r in top if r["fitness"] == fitness)
    parts = []
    for r in sorted(usable, key=lambda r: -IN_VIVO_MODEL_RANK.get(str(r.get("model") or ""), 0)):
        tag = "" if IN_VIVO_MODEL_RANK.get(str(r.get("model") or ""), 0) == best else " [outranked]"
        parts.append(f"{r.get('basis') or r['fitness']} ({r.get('model')}){tag}")
    return {"in_vivo_fitness": fitness, "in_vivo_model": str(model),
            "in_vivo_basis": "; ".join(parts)}


def reference_coverage(frame: pd.DataFrame) -> Dict[str, int]:
    """Summarise how much of a reference table carries usable calls."""

    cols = {c.lower(): c for c in frame.columns}
    vitro = frame[cols["in_vitro_requirement"]].str.lower() if "in_vitro_requirement" in cols else pd.Series(dtype=str)
    vivo = frame[cols["in_vivo_fitness"]].str.lower() if "in_vivo_fitness" in cols else pd.Series(dtype=str)
    return {
        "rows": int(len(frame)),
        "essential": int((vitro == "essential").sum()) if len(vitro) else 0,
        "dispensable": int((vitro == "dispensable").sum()) if len(vitro) else 0,
        "attenuated": int((vivo == "attenuated").sum()) if len(vivo) else 0,
        "no_phenotype": int((vivo == "no_phenotype").sum()) if len(vivo) else 0,
        "not_assayed": int((vivo == "not_assayed").sum()) if len(vivo) else 0,
    }


def file_checksum(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------
# Builder for the bundled Salmonella table
# --------------------------------------------------------------------------
#
# The source is a compilation workbook maintained by the author, one row per
# S. Typhimurium gene, with one column per published screen.  The mapping from
# those columns to the two controlled vocabularies is recorded here rather
# than in a notebook, so that the derivation can be re-read and re-run.

#: Provenance for every source column that the builder reads.
SALMONELLA_SOURCES = [
    {
        "column": "screen_in_vitro_khatiwara2012.tsv (separate resource, not a "
                  "column of the compilation)",
        "field": "in_vitro_requirement / khatiwara_fi_*",
        "assay": "Tn-seq in S. Typhimurium ATCC 14028; fitness index = output "
                 "reads / input reads, neutral 1.0, the authors call a gene "
                 "required at <= 0.2. Four arms: LB 37 C (optimal growth) is "
                 "the only one used to call in vitro fitness cost; diluted LB, "
                 "LB + bile acid and LB 42 C are deliberate proxies for host "
                 "stressors and are carried for audit only. Added because the "
                 "compilation's graded in vitro column is transferred from "
                 "S. Typhi, leaving rung 6 resting on one cross-serovar "
                 "screen; this one is measured in Typhimurium. 110 of the 113 "
                 "genes in Table S3 transfer to SL1344 one-to-one. The screen "
                 "analysed 1,931 of ~4,500 genes and published only those "
                 "crossing threshold, so absence carries no information",
        "citation": "Khatiwara A et al. (2012) Appl Environ Microbiol 78:3098",
        "pmid": "22367088",
    },
    {
        "column": "screen_in_vitro_canals2019.tsv (separate resource, not a "
                  "column of the compilation)",
        "field": "in_vitro_requirement / canals_fi_*",
        "assay": "TIS in S. Typhimurium ST313 D23580; index = insertion-site "
                 "density in the output over the input, median-normalised so "
                 "neutral is 1.0 and a cost is called at <= 0.2, the same "
                 "five-fold depletion as the ATCC 14028 screen. Three arms: LB "
                 "37 C (optimal growth) is the only one used to call in vitro "
                 "fitness cost; InSPI2 is a SPI-2-inducing minimal medium and "
                 "the macrophage arm is a host condition, both carried for "
                 "audit only -- aroA reads 1.18 in LB and 0.05 in InSPI2, so "
                 "gating on InSPI2 would discard the best-validated attenuated "
                 "mutants in the literature. The published three-level "
                 "requirement call is NOT read: it marks 17 of 41 SPI-2 genes "
                 "required in LB because the input library covers SPI-2 "
                 "sparsely, whereas on the normalised index 0 of 36 SPI-2 "
                 "genes reach <= 0.2 in LB and 20 of 36 fall below 0.5 in "
                 "macrophages. Every gene is reported, not only those crossing "
                 "threshold, so this is the one in vitro resource whose "
                 "coverage is genome-wide: 3,998 SL1344 genes receive a value. "
                 "Transferred from ST313 to SL1344 (ST19) by protein identity "
                 "(3,459), >=90% identity (400) or unique gene symbol (139)",
        "citation": "Canals R et al. (2019) PLoS Pathog 15:e1007948",
        "pmid": "31560731",
    },
    {
        "column": "screen_oral_tradis_chaudhuri2013.tsv (separate resource, not a "
                  "column of the compilation)",
        "field": "in_vivo_fitness / in_vivo_effect_rank / in_vivo_model",
        "assay": "Transposon-directed insertion-site sequencing after ORAL challenge; "
                 "the authors' core set of genes significantly attenuated in all "
                 "three host species assayed (chicken, pig, cattle). Supplies an "
                 "attenuation call and the source's own top ordinal category; the "
                 "paper publishes no per-gene fold change for the core set, so no "
                 "magnitude is recorded. 591 of the published 611 genes resolve to "
                 "an SL1344 coding sequence. Added because every other screen here "
                 "is a systemic mouse infection, in which SPI-1 is largely bypassed",
        "citation": "Chaudhuri RR et al. (2013) PLoS Genet 9:e1003456",
        "pmid": "23637626",
    },
    {
        "column": "screen_in_vivo_typhi_31447308.tsv (separate resource, not a "
                  "column of the compilation)",
        "field": "in_vivo_fitness / in_vivo_effect_log2 / in_vivo_significance",
        "assay": "Genome-wide TraDIS of S. Typhi Ty2 in HUMANIZED MICE, the only "
                 "animal screen here from another serovar and the only one "
                 "reporting a per-gene FDR. effect_log2 = -logFC, positive = "
                 "depleted in the animal; attenuation is called at the authors' "
                 "own q <= 0.05, not retuned. Controls pass cleanly: 0 of the 44 "
                 "assayed ribosomal proteins are called, and the complete 10-gene "
                 "Vi antigen locus is recovered. 74 of 4,589 genes are significant "
                 "(1.6%). Transferred to SL1344 by identical protein accession, "
                 "unique gene symbol, or a single >=90% identity match: 3,034 of "
                 "4,868 genes (62.3%) receive a value and 47 of the 74 significant "
                 "genes transfer. The 27 that do not are mostly structures "
                 "Typhimurium does not share -- the Vi locus, 5 O-antigen rfb and "
                 "4 waa outer-core genes. Absence of a value is NOT evidence of no "
                 "phenotype. Not used to license an in vitro dispensable call, "
                 "because the viable mutant existed in Typhi",
        "citation": "Nickerson KP et al. (2019) Cell Host Microbe 26:245-258",
        "pmid": "31447308",
    },
    {
        "column": "Maskell TMDH (importance -1 to 2)",
        "field": "tmdh_importance",
        "assay": "Transposon-mediated differential hybridisation, acute systemic "
                 "infection of BALB/c mice, attenuation score on a -1 to 2 scale",
        "citation": "Chaudhuri RR et al. (2009) PLoS Pathog 5:e1000529",
        "pmid": "19649318",
    },
    {
        "column": "Santiviago -log2(CI)",
        "field": "competitive_index_log2",
        "assay": "Pooled targeted deletion mutants, intraperitoneal infection of "
                 "BALB/c mice, spleen recovery; -log2 competitive index",
        "citation": "Santiviago CA et al. (2009) PLoS Pathog 5:e1000477",
        "doi": "10.1371/journal.ppat.1000477",
    },
    {
        "column": "Monack 16518469 NRAMP",
        "field": "persistence_screen",
        "assay": "Genome-wide screen for genes required for long-term systemic "
                 "infection of Nramp1-positive mice",
        "pmid": "16518469",
    },
    {
        "column": "Falkow 16113260",
        "field": "falkow_screen",
        "assay": "Microarray-based negative selection of a signature-tagged "
                 "transposon mutant library in RAW 264.7 macrophages and in "
                 "the spleens of BALB/cJ mice. The column name is a "
                 "misnomer carried from the source workbook: the senior "
                 "author is Falkow but the first author is Chan, and the "
                 "assay is a genome-wide negative-selection screen, not a "
                 "targeted one. It is the THINNEST input in this table at 18 "
                 "genes (0.4% of the genome), so it is carried for its "
                 "direction rather than its coverage",
        "citation": "Chan K, Kim CC, Falkow S (2005) Infect Immun 73:5438",
        "pmid": "16113260",
    },
    {
        "column": "STM 15009898",
        "field": "stm_screen -> in_vitro_requirement (essential)",
        "assay": "Insertion-duplication mutagenesis trapping lethal insertions; identifies "
                 "genes essential for Salmonella Typhimurium UNDER LABORATORY CONDITIONS. "
                 "This is an IN VITRO essentiality screen, not an in vivo one, and is read "
                 "as essentiality measured in Salmonella itself rather than transferred "
                 "from an E. coli orthologue",
        "pmid": "15009898",
    },
    {
        "column": "essential in E. coli (Keio)",
        "field": "keio_essential",
        "assay": "Essential in the Escherichia coli K-12 single-gene knockout "
                 "collection (Keio), transferred to the Salmonella orthologue",
        "citation": "Baba T et al. (2006) Mol Syst Biol 2:2006.0008",
        "pmid": "16738554",
    },
    {
        "column": "essentialiaty PEC",
        "field": "pec_essential",
        "assay": "Essential in the Profiling of Escherichia coli Chromosome "
                 "database, transferred to the Salmonella orthologue",
    },
    {
        "column": "essential in STm 23470992 log(P-value)",
        "field": "tradis_stm_logp",
        "assay": "Dense transposon insertion library in S. Typhimurium, log "
                 "P-value; RECORDED BUT NOT USED FOR GATING (see note)",
        "pmid": "23470992",
    },
    {
        "column": "essential in vitro in STY 19826075",
        "field": "tradis_sty_logp",
        "assay": "Dense transposon insertion library in S. Typhi, log P-value; "
                 "RECORDED BUT NOT USED FOR GATING (see note)",
        "pmid": "19826075",
    },
    {
        "column": "Proteom Milz",
        "field": "spleen_proteome",
        "assay": "Protein abundance in the spleen proteome",
    },
    {
        "column": "detected in vivo (Neil proteomes)",
        "field": "in_vivo_proteome",
        "assay": "Protein abundance in in vivo proteomes",
    },
]

#: Why the two transposon-insertion columns are carried but not used to call
#: essentiality.  They assign strongly negative log P-values to ssaD, ssaJ,
#: ssaT, spiA, sseJ and sseK2 -- intact SPI-2 genes that the same workbook
#: scores as maximally attenuating in vivo.  A gene cannot simultaneously be
#: required for growth in broth and be a dispensable mutant that was built and
#: shown to be attenuated in mice, so using these columns as an essentiality
#: gate would delete the core of SPI-2 from the host-interaction tier.  The
#: values are retained in the table for audit and for anyone who wants to
#: apply their own threshold.
TRADIS_NOTE = (
    "Transposon-insertion log P-values from the S. Typhimurium and S. Typhi "
    "libraries are carried in the table but are not used to call in vitro "
    "essentiality, because they flag intact SPI-2 genes that the same source "
    "scores as attenuated in vivo."
)

#: Thresholds on the two quantitative in vivo screens.  Both are the natural
#: cut points of their own scales: 1 on the TMDH -1..2 importance scale is the
#: lowest score the source calls important, and -log2(CI) >= 1 is a two-fold
#: under-representation of the mutant in the output pool.
TMDH_ATTENUATED_AT = 1.0
COMPETITIVE_INDEX_ATTENUATED_AT = 1.0


#: Oral-route colonisation core set, shipped as its own resource.
ORAL_SCREEN_FILE = os.path.join(
    RESOURCE_DIR, "screen_oral_tradis_chaudhuri2013.tsv"
)

#: Second in vitro fitness screen, shipped as its own resource: Tn-seq fitness
#: indices in S. Typhimurium ATCC 14028 (Khatiwara et al. 2012, PMID 22367088).
#: Its value is that it is measured in Typhimurium itself, where the graded in
#: vitro column of the compilation is transferred from S. Typhi, so the two are
#: independent on the one axis rung 6 depends on.
KHATIWARA_SCREEN_FILE = os.path.join(
    RESOURCE_DIR, "screen_in_vitro_khatiwara2012.tsv"
)

#: The paper's own cut point: a fitness index of 0.2 is a five-fold depletion
#: of the gene's insertion mutants during selection. The index is a ratio whose
#: neutral value is 1.0, so SMALLER is worse.
KHATIWARA_FITNESS_COST_AT = 0.2

#: Where the native screen OVERRIDES the transferred Typhi call. A value that
#: merely failed the paper's own 0.2 cutoff is not evidence of no cost -- a
#: fitness index of 0.33 is still a three-fold depletion, and it agrees in
#: direction with a Typhi "required" call rather than contradicting it. Only a
#: gene whose mutants were at or above neutral abundance (index >= 1.0) in the
#: target serovar genuinely contradicts "required for growth in rich medium",
#: and only those override. Of the 19 genes where the two sources differ on the
#: published cutoff, 15 are sub-threshold depletions that agree in direction and
#: 4 are at or above neutral.
KHATIWARA_NEUTRAL_AT = 1.0

#: Only the LB-37C arm is used to call in vitro fitness cost. The other three
#: arms -- diluted LB, LB + bile acid and LB at 42 C -- are in vitro in the
#: trivial sense that no animal is involved, but each is a deliberate proxy for
#: a host stressor (the nutrient-poor macrophage vacuole, intestinal bile, avian
#: body temperature), and the bile arm is dominated by LPS core and O-antigen
#: genes. Calling those a "fitness cost" and using them to demote genes out of
#: rung 6 would subtract the host's own selective pressures from the definition
#: of host specificity, which is backwards. They are carried for audit.
KHATIWARA_GATING_ARM = "fitness_index_lb"


def _load_khatiwara_in_vitro_screen(locus_tags: pd.Series) -> pd.DataFrame:
    """Per-gene in vitro fitness indices from the Typhimurium Tn-seq screen.

    Takes the locus-tag COLUMN of the table being built and returns a frame on
    that same index with one column per growth condition, NaN where the screen
    has no value. Matches either locus-tag convention, because the compilation
    uses old SL1344_#### tags and a live run uses RefSeq SL1344_RS##### ones.

    Absence is NOT neutrality. The screen analysed 1,931 of about 4,500 genes --
    only those with at least three recovered insertion sites -- and the
    published table lists only the genes that crossed the threshold somewhere.
    A gene missing from the resource may have been above threshold or may never
    have been assessed, and the published table cannot tell the two apart, so it
    contributes no in vitro measurement at all rather than a reassuring one.
    """
    return _load_screen_arms(
        KHATIWARA_SCREEN_FILE,
        ["fitness_index_lb", "fitness_index_dlb",
         "fitness_index_bile", "fitness_index_42c"],
        locus_tags,
        "screen_in_vitro_khatiwara2012.tsv is missing; the Salmonella "
        "reference table will lose its ATCC 14028 in vitro cross-check.",
        "The Typhimurium ATCC 14028 Tn-seq in vitro screen matched no gene in "
        "this table; check that its locus tags use the same convention.",
    )


#: Third in vitro fitness screen: transposon-insertion-sequencing indices in
#: S. Typhimurium ST313 strain D23580 (Canals et al. 2019, PMID 31560731). It
#: is the best-covered in vitro resource in the package -- every gene in the
#: genome is reported, not only those crossing a threshold -- so unlike the
#: Khatiwara table, absence from it IS informative for the genes it omits only
#: because they had fewer than three recovered insertions.
CANALS_SCREEN_FILE = os.path.join(
    RESOURCE_DIR, "screen_in_vitro_canals2019.tsv"
)

#: As for Khatiwara: only the optimal-growth arm gates. The InSPI2 arm is a
#: SPI-2-inducing minimal medium and the macrophage arm is a host condition;
#: both are carried for audit. The file's own header shows why -- aroA scores
#: 1.18 in LB and 0.05 in InSPI2, and its cost in minimal medium IS its
#: mechanism of attenuation in a host poor in aromatic amino acids.
CANALS_GATING_ARM = "fitness_index_lb"


def _load_canals_in_vitro_screen(locus_tags: pd.Series) -> pd.DataFrame:
    """Per-gene in vitro fitness indices from the D23580 TIS screen.

    Same contract as :func:`_load_khatiwara_in_vitro_screen`. The published
    three-level requirement call is deliberately NOT read: it marks 17 of 41
    SPI-2 genes "required" in LB because the input library covers the SPI-2
    region sparsely, and the input-normalised index in this resource does not.
    """
    return _load_screen_arms(
        CANALS_SCREEN_FILE,
        ["fitness_index_lb", "fitness_index_inspi2",
         "fitness_index_macrophage"],
        locus_tags,
        "screen_in_vitro_canals2019.tsv is missing; the Salmonella reference "
        "table will lose its best-covered in vitro fitness column.",
        "The D23580 TIS in vitro screen matched no gene in this table; check "
        "that its locus tags use the same convention.",
    )


def _load_screen_arms(
    path: str,
    arms: List[str],
    locus_tags: pd.Series,
    missing_message: str,
    nomatch_message: str,
) -> pd.DataFrame:
    """Shared reader for the per-arm in vitro fitness screen resources.

    Takes the locus-tag COLUMN of the table being built and returns a frame on
    that same index with one column per arm, NaN where the screen has no value.
    Matches either locus-tag convention, because a bundled compilation uses old
    SL1344_#### tags and a live run uses RefSeq SL1344_RS##### ones.

    Absence is NOT neutrality in general: a screen that published only the
    genes crossing its threshold cannot distinguish "above threshold" from
    "never assessed", so a missing gene contributes no measurement at all
    rather than a reassuring one. Each resource's header states which case it
    is.
    """
    empty = pd.DataFrame(
        {a: pd.Series(np.nan, index=locus_tags.index) for a in arms}
    )
    if not os.path.exists(path):
        warnings.warn(missing_message, RuntimeWarning)
        return empty
    table = pd.read_csv(path, sep="\t", comment="#")
    lookup: Dict[str, Dict[str, float]] = {}
    for _, row in table.iterrows():
        values = {a: pd.to_numeric(row.get(a), errors="coerce") for a in arms}
        for key in (row.get("locus_tag"), row.get("old_locus_tag")):
            text = "" if key is None else str(key).strip()
            if text and text.lower() != "nan":
                lookup[text] = values
    found = empty.copy()
    n = 0
    for i, tag in zip(locus_tags.index, locus_tags):
        values = lookup.get(str(tag).strip())
        if values is None:
            continue
        n += 1
        for a in arms:
            found.at[i, a] = values[a]
    if not n:
        warnings.warn(nomatch_message, RuntimeWarning)
    return found


#: In vivo TraDIS of S. Typhi Ty2 in humanized mice (Nickerson et al. 2019,
#: PMID 31447308). The only genome-wide IN VIVO screen available for Typhi, and
#: the only animal screen in this table from a different serovar. Transferred to
#: SL1344 by protein identity, unique gene symbol, or >= 90% protein identity.
TYPHI_INVIVO_SCREEN_FILE = os.path.join(
    RESOURCE_DIR, "screen_in_vivo_typhi_31447308.tsv"
)

#: The authors' own FDR cutoff. Not retuned: a gene is attenuated here only if
#: the source called it significant.
TYPHI_INVIVO_MAX_Q = 0.05


def _load_typhi_in_vivo_screen(locus_tags: pd.Series) -> pd.DataFrame:
    """Transferred in vivo effect sizes from the Typhi humanized-mouse screen.

    Returns a frame indexed like ``locus_tags`` with ``effect_log2`` (positive =
    mutant depleted in the animal) and ``q_value``. An all-empty frame is
    returned when the resource is absent, so a stripped-down install still
    builds; the caller reports it.

    Coverage is 62.3% of SL1344 genes. A gene with no value was not measured
    after transfer, which is NOT the same as measured and unaffected, so the
    absent rows must never be read as evidence against attenuation.
    """
    empty = pd.DataFrame(
        {"effect_log2": np.nan, "q_value": np.nan}, index=locus_tags.index
    )
    if not os.path.exists(TYPHI_INVIVO_SCREEN_FILE):
        warnings.warn(
            "screen_in_vivo_typhi_31447308.tsv is missing; the Salmonella "
            "reference table will carry one fewer in vivo source.",
            RuntimeWarning,
        )
        return empty
    table = pd.read_csv(TYPHI_INVIVO_SCREEN_FILE, sep="\t", comment="#")
    table = table.dropna(subset=["locus_tag"])
    # The compilation is keyed on legacy SL1344_#### tags while the screen was
    # transferred onto current RefSeq tags, so both are accepted. Legacy first:
    # it is what this table actually uses.
    frames = []
    for key in ("old_locus_tag", "locus_tag"):
        if key not in table.columns:
            continue
        part = table.dropna(subset=[key]).copy()
        part[key] = part[key].astype(str).str.strip()
        part = part[part[key] != ""].drop_duplicates(subset=[key], keep="first")
        frames.append(part.set_index(key)[["effect_log2", "q_value"]])
    if not frames:
        return empty
    lut = pd.concat(frames)
    lut = lut[~lut.index.duplicated(keep="first")]
    keys = pd.Index([str(t).strip() for t in locus_tags])
    out = pd.DataFrame(
        {
            "effect_log2": lut["effect_log2"].reindex(keys).to_numpy(),
            "q_value": lut["q_value"].reindex(keys).to_numpy(),
        },
        index=locus_tags.index,
    )
    if out["effect_log2"].notna().sum() == 0:
        warnings.warn(
            "The Typhi humanized-mouse screen matched no gene in this table; "
            "check that its locus tags use the same convention.",
            RuntimeWarning,
        )
    return out


def _load_oral_colonisation_screen(locus_tags: pd.Series) -> pd.Series:
    """Membership of the published oral-colonisation core set.

    ``locus_tags`` is the locus-tag COLUMN of the table being built, not its
    index: the builder indexes rows by the source workbook's row number, so
    matching on the index would silently match nothing. The returned series
    keeps that same index so it can be combined with the other screen flags.

    Returns an all-False series when the resource is absent, so that a
    stripped-down install still builds a table -- but says so, because a
    silently missing screen would make SPI-1 look dispensable again.
    """
    if not os.path.exists(ORAL_SCREEN_FILE):
        warnings.warn(
            "screen_oral_tradis_chaudhuri2013.tsv is missing; the Salmonella "
            "reference table will carry systemic-infection evidence only and "
            "will under-report SPI-1.",
            RuntimeWarning,
        )
        return pd.Series(False, index=locus_tags.index)
    table = pd.read_csv(ORAL_SCREEN_FILE, sep="\t", dtype=str, comment="#")
    hits = set(str(x).strip() for x in table["locus_tag"].dropna())
    hits |= set(
        str(x).strip()
        for x in table.get("old_locus_tag", pd.Series(dtype=str)).dropna()
    )
    found = pd.Series(
        [str(t).strip() in hits for t in locus_tags], index=locus_tags.index
    )
    if not found.any():
        warnings.warn(
            "The oral-colonisation core set matched no gene in this table; "
            "check that its locus tags use the same convention.",
            RuntimeWarning,
        )
    return found


def build_salmonella_sl1344_reference(
    workbook_path: str,
    out_path: str = "",
    sheet_name: str = "STM genome",
) -> Dict[str, object]:
    """Derive the bundled Salmonella reference table from the source workbook.

    Returns a summary dict; writes the table and a provenance sidecar.
    """

    raw = pd.read_excel(workbook_path, sheet_name=sheet_name)

    def col(name):
        if name not in raw.columns:
            raise KeyError(f"Source workbook has no column {name!r}.")
        return raw[name]

    locus = col("SL1344 locus     22538806").astype(str).str.strip()
    keep = locus.str.match(r"^SL\d{4}$", na=False)
    out = pd.DataFrame(index=raw.index[keep])
    # The workbook writes SL1344 locus tags without the strain prefix.
    out["locus_tag"] = "SL1344_" + locus[keep].str[2:]
    out["stm_number"] = col("Gene number")[keep].astype(str).str.strip()
    out["gene_symbol"] = col("Gene symbol")[keep].fillna("").astype(str).str.strip()

    num = lambda name: pd.to_numeric(col(name)[keep], errors="coerce")
    flag = lambda name: col(name)[keep].notna()

    tmdh = num("Maskell TMDH (importance -1 to 2)")
    cindex = num("Santiviago -log2(CI)")
    persistence = flag("Monack 16518469 NRAMP")
    falkow = flag("Falkow 16113260")
    stm = flag("STM 15009898")
    keio = flag("essential in E. coli (Keio)")
    pec = flag("essentialiaty PEC")

    # ---- in vitro requirement ------------------------------------------
    # Oral-route colonisation screen, loaded from its own resource file
    # because it is not a column of the compilation. Every other screen here
    # is a SYSTEMIC mouse infection, in which SPI-1 is largely bypassed: its
    # competitive indices run 0.4-2.1 and no SPI-1 gene reaches TMDH
    # importance 2, so a ladder built on those screens alone under-reports
    # the intestinal invasion apparatus. This set is the authors' own core
    # set -- significant after oral challenge in chicken, pig AND cattle.
    oral = _load_oral_colonisation_screen(out["locus_tag"])
    # Fourth in vivo source, and the only one from another serovar: TraDIS of
    # S. Typhi Ty2 in humanized mice. It supplies a graded effect size with the
    # authors' own FDR, so it can reach the required rung on a number rather
    # than an ordinal call. It is deliberately NOT added to `assayed`: the
    # viable mutant existed in Typhi, and that does not license a Typhimurium
    # in vitro "dispensable" call.
    typhi = _load_typhi_in_vivo_screen(out["locus_tag"])
    typhi_eff = typhi["effect_log2"]
    typhi_q = typhi["q_value"]
    typhi_atten = (
        typhi_eff.notna()
        & (typhi_q <= TYPHI_INVIVO_MAX_Q)
        & (typhi_eff > 0)
    )
    # Second, Typhimurium-native in vitro fitness screen. See
    # KHATIWARA_GATING_ARM for why only the optimal-growth arm gates.
    khat = _load_khatiwara_in_vitro_screen(out["locus_tag"])
    khat_lb = khat[KHATIWARA_GATING_ARM]
    khat_cost = khat_lb.notna() & (khat_lb <= KHATIWARA_FITNESS_COST_AT)
    # A measurement in THIS serovar outranks one transferred from another, the
    # same precedence the in vivo side already applies to infection models. Where
    # the Typhimurium screen measured a gene and found its mutants undepleted in
    # rich medium, that overrides a Typhi "required in rich medium" call -- and
    # the override is recorded per gene rather than applied silently, because it
    # is the only place in the table where two sources are allowed to disagree
    # on the in vitro axis.
    khat_neutral = khat_lb.notna() & (khat_lb >= KHATIWARA_NEUTRAL_AT)
    # Third screen, also Typhimurium-native: the D23580 TIS landscape. It is
    # read on the same scale and with the same gating rule, and it matters
    # because it covers 3,998 genes against Khatiwara's 110, so it supplies an
    # in vitro observation for most of the genome rather than for a list.
    can = _load_canals_in_vitro_screen(out["locus_tag"])
    can_lb = can[CANALS_GATING_ARM]
    can_cost = can_lb.notna() & (can_lb <= KHATIWARA_FITNESS_COST_AT)
    can_neutral = can_lb.notna() & (can_lb >= KHATIWARA_NEUTRAL_AT)
    # The two native screens have equal standing -- same serovar, same assay
    # family, neither transferred -- so neither outranks the other on
    # provenance, and the rule that resolves them is symmetric: a cost call
    # stands unless an equally native screen puts the same gene at or above
    # neutral abundance, which is a positive contradiction rather than a mere
    # failure to cross a cutoff. On the 106 genes both measure this fires
    # twice, for rfaH (0.20 against 1.55, Khatiwara's value sitting exactly on
    # its own cutoff) and SL1344_RS09150 (0.19 against 1.19). The two screens
    # otherwise agree in rank, at Spearman rho 0.518 (P = 1.3e-8), with
    # Khatiwara's cost genes at a median 0.554 in the D23580 index against
    # 0.995 for its neutral ones (Mann-Whitney P = 6.9e-6).
    native_neutral = khat_neutral | can_neutral
    native_cost = (khat_cost | can_cost) & ~native_neutral
    native_measured = khat_lb.notna() | can_lb.notna()
    # Essentiality is taken from the two orthologue-based binary calls only.
    # Graded in vitro fitness, transferred from the Salmonella Typhi TraDIS
    # screen (Langridge et al. 2009, PMID 19826075; 560 citations, the most
    # cited Salmonella transposon screen). The transfer is gated on orthology
    # quality: the compilation's own STY columns must show a PRESENT homologue
    # (>92% DNA homology) AND a reciprocal best hit, which holds for 3,748 of
    # 4,271 genes. A gene the Typhi authors called required in rich medium,
    # that is not already essential, is a viable-but-impaired mutant --
    # exactly the `fitness_cost` class, and the only thing rung 6's in vitro
    # clause can act on.
    #
    # Why this screen and not the Typhimurium one (PMID 23470992), which is
    # the same assay in the same species: the Typhimurium log P flags 16 of the
    # 36 annotated SPI-2 genes at log P <= -10, and 30 of 36 at the threshold
    # used here, which no reading of SPI-2 biology supports. The Typhi screen
    # flags ZERO SPI-1 and ZERO SPI-2 genes, and it is not for want of the
    # genes being there -- 32 of the 36 have a Typhi orthologue, at a median
    # 97% DNA and 98.8% amino acid identity, with a reciprocal best hit, so the
    # screen could have called them and did not. The SPI-2 signal is specific to that one library. Against
    # Keio/PEC the Typhi screen reads 93% sensitivity and 8.3% false-positive
    # rate (odds ratio 154, P = 7e-245) and recovers 52 of 54 ribosomal
    # proteins, so it is a sound essentiality/fitness measure in its own right.
    sty_logp = num("essential in vitro in STY 19826075")
    sty_present = (
        col("STY homologue status (present = >92 DNA homology)")[keep]
        .astype(str).str.strip().eq("P")
        & col("STY reciprocal best hit (amino acid)")[keep].notna()
    )
    sty_required = sty_logp.notna() & sty_present

    # PMID 15009898 is an IN VITRO essentiality screen -- insertion-duplication
    # mutagenesis identifying 'essential genes of Salmonella enterica serovar
    # Typhimurium under laboratory conditions'. Earlier revisions of this
    # builder read it as a signature-tagged mutagenesis IN VIVO screen and
    # credited its 15 genes as attenuated in a host, which is the wrong
    # direction entirely: the paper says they cannot be inactivated in the
    # lab. It now contributes essentiality, measured in Salmonella itself
    # rather than transferred from an E. coli orthologue.
    essential = keio | pec | stm
    # Required in Typhi, or depleted five-fold in optimal growth in
    # Typhimurium, but not essential = viable but impaired.
    in_vitro_discordant = sty_required & native_neutral
    fitness_cost = ((sty_required & ~native_neutral) | native_cost) & ~essential
    # A gene that was assayed as a viable insertion or deletion mutant is
    # dispensable in vitro by construction: the mutant exists.  This is the
    # strongest statement in the table and it overrides annotation-based
    # inference in the opposite direction.
    # A gene assayed in vivo had a viable mutant, or there would have been
    # nothing to assay, so the oral screen's coverage licenses `dispensable`
    # for the 93 genes it reaches that no systemic screen covered -- which is
    # why `oral` belongs in this union and not only in the in vivo one.
    assayed = tmdh.notna() | cindex.notna() | persistence | falkow | oral
    # A gene the Typhimurium screen MEASURED and did not deplete is a positive
    # in vitro observation of a viable, unimpaired mutant -- which is what
    # rung 6 asks for, and the only evidence of that kind measured in this
    # serovar. It licenses `dispensable` on its own.
    in_vitro_measured = native_measured
    requirement = pd.Series("unknown", index=out.index)
    requirement[(assayed | in_vitro_measured) & ~essential] = "dispensable"
    requirement[fitness_cost] = "fitness_cost"
    requirement[essential] = "essential"
    out["in_vitro_requirement"] = requirement
    basis = pd.Series("", index=out.index)
    basis[keio] = "essential in the Keio collection (orthologue)"
    basis[pec & ~keio] = "essential in PEC (orthologue)"
    basis[keio & pec] = "essential in the Keio collection and in PEC (orthologue)"
    basis[(requirement == "dispensable")] = "viable mutant constructed and assayed in vivo"
    measured_only = in_vitro_measured & ~assayed & (requirement == "dispensable")
    basis[measured_only] = (
        "insertion mutants recovered and not depleted in optimal growth in "
        "S. Typhimurium (Tn-seq fitness index "
        + khat_lb[measured_only].round(2).astype(str) + ")"
    )
    basis[stm & essential] = "essential in vitro (insertion-duplication mutagenesis, Salmonella)"
    basis[fitness_cost & sty_required] = (
        "required for growth in rich medium in S. Typhi (TraDIS log P "
        + sty_logp[fitness_cost & sty_required].round(1).astype(str)
        + "), orthologue present with a reciprocal best hit"
    )
    both_cost = fitness_cost & sty_required & khat_cost
    basis[both_cost] = (
        "required for growth in rich medium in S. Typhi (TraDIS log P "
        + sty_logp[both_cost].round(1).astype(str)
        + ") and depleted in optimal growth in S. Typhimurium (Tn-seq fitness "
        "index " + khat_lb[both_cost].round(2).astype(str) + ")"
    )
    basis[in_vitro_discordant & ~essential] = (
        "called required in rich medium in S. Typhi (TraDIS log P "
        + sty_logp[in_vitro_discordant & ~essential].round(1).astype(str)
        + ") but measured undepleted in optimal growth in S. Typhimurium "
        "(Tn-seq fitness index "
        + khat_lb[in_vitro_discordant & ~essential].round(2).astype(str)
        + "); the serovar-native measurement takes precedence"
    )
    khat_only_cost = fitness_cost & khat_cost & ~sty_required
    basis[khat_only_cost] = (
        "depleted five-fold in optimal growth in S. Typhimurium (Tn-seq "
        "fitness index " + khat_lb[khat_only_cost].round(2).astype(str) + ")"
    )
    out["in_vitro_basis"] = basis

    # ---- in vivo fitness -------------------------------------------------
    attenuated = (
        (tmdh >= TMDH_ATTENUATED_AT)
        | (cindex >= COMPETITIVE_INDEX_ATTENUATED_AT)
        | persistence
        | falkow
        | oral
        | typhi_atten
    )
    fitness = pd.Series("not_assayed", index=out.index)
    fitness[assayed] = "no_phenotype"
    fitness[attenuated] = "attenuated"
    out["in_vivo_fitness"] = fitness

    srcs = []
    for idx in out.index:
        hits = []
        if pd.notna(tmdh.get(idx)) and tmdh.get(idx) >= TMDH_ATTENUATED_AT:
            hits.append(f"TMDH importance {tmdh.get(idx):.0f}")
        if pd.notna(cindex.get(idx)) and cindex.get(idx) >= COMPETITIVE_INDEX_ATTENUATED_AT:
            hits.append(f"-log2(CI) {cindex.get(idx):.1f}")
        if persistence.get(idx):
            hits.append("persistent infection screen")
        if falkow.get(idx):
            hits.append("targeted in vivo screen")
        if oral.get(idx):
            hits.append("oral colonisation core set (chicken, pig, cattle)")
        if typhi_atten.get(idx):
            hits.append(
                f"S. Typhi in humanized mice, log2 depletion "
                f"{typhi_eff.get(idx):.1f} (q {typhi_q.get(idx):.2g}), orthologue-transferred"
            )
        srcs.append("; ".join(hits))
    out["in_vivo_basis"] = srcs
    # Every screen in this compilation is a whole-animal infection: BALB/c
    # mice for the TMDH and competitive-index screens, Nramp1-positive mice
    # for the persistence screen. There is no cell-infection, organoid or
    # ex vivo evidence here, so the model column is constant, and the
    # precedence machinery in resolve_in_vivo_call has nothing to arbitrate.
    # Route matters and is recorded. A gene called by the oral screen and by
    # no systemic screen is labelled enteric; one called by both keeps the
    # systemic label with the oral hit visible in in_vivo_basis. The two
    # classes rank equally, so neither silently overrides the other.
    systemic_call = (
        (tmdh >= TMDH_ATTENUATED_AT)
        | (cindex >= COMPETITIVE_INDEX_ATTENUATED_AT)
        | persistence
        | falkow
        | typhi_atten
    )
    out["in_vivo_model"] = np.where(
        fitness == "not_assayed",
        "not_applicable",
        np.where(oral & ~systemic_call, "animal_infection_enteric", "animal_infection"),
    )

    # Normalised effect size: positive = mutant depleted in vivo. The
    # competitive-index column is already reported as -log2(CI), so a mutant
    # recovered ten-fold below the wild type sits at 3.32 with no conversion
    # needed. The TMDH column is an ordinal importance score on a -1..2 scale,
    # not a fold change, so it supplies an attenuation CALL but no magnitude;
    # genes attenuated only by TMDH, the persistence screen or the targeted
    # screens therefore carry no effect size and stop at the contributing
    # rung rather than being promoted on a number that does not exist.
    # The competitive index is measured in Typhimurium and takes precedence.
    # Where the compilation supplies none, the transferred Typhi depletion is
    # used instead, so a gene is never promoted on a cross-serovar number when
    # a serovar-native one exists.
    effect = cindex.copy()
    use_typhi = effect.isna() & typhi_atten
    effect[use_typhi] = typhi_eff[use_typhi]
    out["in_vivo_effect_log2"] = effect
    out.loc[fitness == "not_assayed", "in_vivo_effect_log2"] = np.nan
    # The source reports no per-gene significance value for the competitive
    # index, so the field is left empty and the attenuation call stands on
    # the source's own thresholding.
    # Ordinal strength, for sources that grade attenuation without reporting a
    # magnitude. The TMDH screen scores every gene -1..2 and 2 is its own top
    # attenuation category; refusing to read that, and admitting only genes
    # carrying a competitive index, would make the required rung mean "has a
    # ten-fold competitive index" rather than "is required in vivo", and would
    # silently drop the classic attenuating auxotrophs (aroC, aroD, purA,
    # guaB) that no competitive index in this compilation covers. The score is
    # NOT converted to log2 units -- it is not a ratio -- it is carried on its
    # own ordinal axis and read alongside the numeric one.
    rank = pd.Series("", index=out.index, dtype=object)
    rank[tmdh >= 2] = "strong"
    rank[(tmdh >= 1) & (tmdh < 2)] = "moderate"
    # The oral core set is that paper's own top category: required in all
    # three host species, not in one. It is read exactly as TMDH importance 2
    # is read -- as an ordinal top call, never converted to a magnitude.
    rank[oral] = "strong"
    rank[fitness != "attenuated"] = ""
    out["in_vivo_effect_rank"] = rank
    # Only the Typhi screen reports a per-gene significance value; the
    # competitive index does not, and its attenuation call stands on the
    # source's own thresholding.
    signif = pd.Series(np.nan, index=out.index)
    signif[typhi_atten] = typhi_q[typhi_atten]
    out["in_vivo_significance"] = signif
    out["in_vivo_effect_source"] = np.where(
        out["in_vivo_effect_log2"].notna(),
        "Santiviago competitive index, -log2(CI), positive = depleted in vivo",
        "",
    )
    out.loc[use_typhi, "in_vivo_effect_source"] = (
        "S. Typhi humanized-mouse TraDIS, -log2 fold change, "
        "positive = depleted in vivo (orthologue-transferred)"
    )
    out.loc[
        (out["in_vivo_effect_source"] == "") & (out["in_vivo_effect_rank"] != ""),
        "in_vivo_effect_source",
    ] = "Chaudhuri TMDH importance, ordinal -1..2, 2 = strongest attenuation"

    # ---- raw values, carried for audit ----------------------------------
    out["tmdh_importance"] = tmdh
    out["competitive_index_log2"] = cindex
    out["persistence_screen"] = persistence.map({True: "1", False: ""})
    out["falkow_screen"] = falkow.map({True: "1", False: ""})
    out["stm_screen"] = stm.map({True: "1", False: ""})
    out["keio_essential"] = keio.map({True: "1", False: ""})
    out["pec_essential"] = pec.map({True: "1", False: ""})
    out["tradis_stm_logp"] = num("essential in STm 23470992 log(P-value)")
    out["tradis_sty_logp"] = num("essential in vitro in STY 19826075")
    out["in_vitro_discordant"] = np.where(in_vitro_discordant & ~essential, "1", "")
    out["khatiwara_fi_lb"] = khat["fitness_index_lb"]
    out["khatiwara_fi_dlb"] = khat["fitness_index_dlb"]
    out["khatiwara_fi_bile"] = khat["fitness_index_bile"]
    out["khatiwara_fi_42c"] = khat["fitness_index_42c"]
    out["canals_fi_lb"] = can["fitness_index_lb"]
    out["canals_fi_inspi2"] = can["fitness_index_inspi2"]
    out["canals_fi_macrophage"] = can["fitness_index_macrophage"]
    out["typhi_invivo_log2"] = typhi_eff
    out["typhi_invivo_q"] = typhi_q
    out["spleen_proteome"] = num("Proteom Milz")
    out["in_vivo_proteome"] = num("detected in vivo (Neil proteomes)")

    out = out.drop_duplicates(subset=["locus_tag"], keep="first").reset_index(drop=True)
    validate_reference(out, source=workbook_path)

    out_path = out_path or os.path.join(RESOURCE_DIR, "reference_salmonella_sl1344.tsv")
    out.to_csv(out_path, sep="\t", index=False)

    sidecar = os.path.splitext(out_path)[0] + "_sources.tsv"
    pd.DataFrame(SALMONELLA_SOURCES).to_csv(sidecar, sep="\t", index=False)

    summary = reference_coverage(out)
    summary.update({
        "table": out_path,
        "sources": sidecar,
        "source_workbook": os.path.basename(workbook_path),
        "source_checksum_sha256": file_checksum(workbook_path),
        "note": TRADIS_NOTE,
    })
    return summary



# ---------------------------------------------------------------------------
# Building a reference table from the user template
# ---------------------------------------------------------------------------
#: Roles a column of a user template may declare. The point of the template is
#: that the USER supplies raw screen data -- the numbers and flags their source
#: actually publishes -- and declares what each column means, while the module
#: does the interpretation into the controlled vocabularies. Asking a user to
#: hand-convert their screens into `in_vitro_requirement` / `in_vivo_fitness`
#: puts the scientific judgement in the least auditable place: a spreadsheet
#: edit nobody records. A declaration is recorded, travels into the provenance
#: sidecar and can be re-read.
TEMPLATE_ROLES = {
    "in_vitro_essential_flag": "any non-empty value means no viable mutant was obtained",
    "in_vitro_fitness_score": "numeric in vitro fitness; beyond `attenuated_at` in the "
                              "declared direction the gene is called fitness_cost "
                              "(viable but impaired) -- this is what rung 6 needs. "
                              "`attenuated_at` is read in the column's own units, as "
                              "printed in the paper (0.2 for a fitness index where "
                              "smaller is worse, -2.0 for a log P)",
    "in_vivo_attenuation_flag": "any non-empty value means the mutant was attenuated in a host",
    "in_vivo_effect_log2": "numeric in vivo effect size; the sign convention is declared",
    "in_vivo_ordinal": "ordinal attenuation grade; `top_category` names its strongest value",
    "in_vivo_significance": "per-gene significance value for the in vivo call",
    "annotation_only": "carried into the table for audit, never used for gating",
}

#: Columns the Screen definitions sheet may carry. Only `column` and `role` are
#: required; the rest default.
TEMPLATE_DEFINITION_COLUMNS = (
    "column", "role", "model", "attenuated_at", "positive_means",
    "top_category", "assay", "citation", "pmid", "neutral_at", "overrides",
)


def _template_text(value) -> str:
    """A declaration cell as text, with pandas' blanks normalised away.

    A missing cell arrives as float NaN, which is truthy, so ``str(x or "")``
    yields the string "nan" and an absent declaration reads as a present one.
    """
    if value is None:
        return ""
    if isinstance(value, float) and np.isnan(value):
        return ""
    text = str(value).strip()
    return "" if text.lower() in ("nan", "none") else text


def _template_direction(value, role: str = "", column: str = "") -> int:
    """+1 when a LARGER value means a worse mutant, -1 when smaller does.

    Required for the numeric roles and ignored for the rest. Blank is NOT
    silently defaulted on a numeric column: getting the sign wrong inverts the
    result, turning enriched mutants into attenuated ones, so a missing
    declaration raises with the two options spelled out.
    """
    text = "" if value is None or (isinstance(value, float) and np.isnan(value)) else str(value).strip().lower()
    if text in ("nan", "none"):
        text = ""
    numeric_role = role in ("in_vitro_fitness_score", "in_vivo_effect_log2")
    if text == "":
        if numeric_role:
            raise ValueError(
                f"Column {column!r} has role {role!r}, which is numeric, but declares no "
                "positive_means. Set it to 'depleted_in_host' if a LARGER number means the "
                "mutant did worse (e.g. -log2(CI)), or 'enriched_in_host' if a SMALLER or "
                "more negative number does (e.g. a raw log2 fold change). It is not "
                "defaulted, because the wrong sign silently inverts the result."
            )
        return 1
    if text in ("larger_is_worse", "depleted_in_host", "positive_means_depleted",
                "higher_is_depleted", "depleted"):
        return 1
    if text in ("smaller_is_worse", "enriched_in_host", "negative_means_depleted",
                "lower_is_depleted", "log2_fold_change", "competitive_index"):
        return -1
    raise ValueError(
        f"positive_means={value!r} is not recognised. Use 'depleted_in_host' when a "
        "LARGER number means the mutant did worse (e.g. -log2(CI)), or "
        "'enriched_in_host' when a SMALLER/more negative number does (e.g. a raw "
        "log2 fold change)."
    )


def build_reference_from_template(
    workbook_path: str,
    out_path: str = "",
    data_sheet: str = "Gene data",
    definition_sheet: str = "Screen definitions",
) -> Dict[str, object]:
    """Derive a reference table from a filled-in user template workbook.

    The template carries raw screen columns plus a sheet declaring what each
    column is. Returns a summary dict; writes the table and a provenance
    sidecar next to ``out_path``.

    Minimum that works: a ``locus_tag`` column and ONE declared screen column.
    Everything else is optional -- with no in vitro column every gene is
    ``unknown`` in vitro, which leaves rungs 4 and 5 computable and rung 6
    empty, and the run log says so.
    """

    data = pd.read_excel(workbook_path, sheet_name=data_sheet)
    defs = pd.read_excel(workbook_path, sheet_name=definition_sheet)
    defs = defs.dropna(subset=["column"]) if "column" in defs.columns else defs
    if "column" not in defs.columns or "role" not in defs.columns:
        raise ValueError(
            f"{definition_sheet!r} must have at least 'column' and 'role' columns; "
            f"found {list(defs.columns)}."
        )

    cols = {str(c).strip().lower(): c for c in data.columns}
    if "locus_tag" not in cols:
        raise ValueError(
            f"{data_sheet!r} has no 'locus_tag' column. It must carry the locus tags of "
            "the assembly the run resolves to (old locus tags are also accepted, in an "
            "'old_locus_tag' column)."
        )

    out = pd.DataFrame(index=data.index)
    out["locus_tag"] = data[cols["locus_tag"]].astype(str).str.strip()
    for optional in ("old_locus_tag", "gene_symbol", "stm_number"):
        if optional in cols:
            out[optional] = data[cols[optional]].fillna("").astype(str).str.strip()

    declared, notes, sources = [], [], []
    for _, row in defs.iterrows():
        name = str(row["column"]).strip()
        role = str(row["role"]).strip().lower()
        if role not in TEMPLATE_ROLES:
            raise ValueError(
                f"Column {name!r} declares role {role!r}, which is not recognised. "
                f"Allowed roles: {sorted(TEMPLATE_ROLES)}."
            )
        if name not in data.columns:
            notes.append(f"Declared column {name!r} is absent from {data_sheet!r}; skipped.")
            continue
        model = str(row.get("model") or "animal_infection").strip() or "animal_infection"
        if role.startswith("in_vivo") and model not in IN_VIVO_MODEL_RANK:
            raise ValueError(
                f"Column {name!r} declares model {model!r}; allowed: {list(IN_VIVO_MODELS)}."
            )
        declared.append({
            "name": name, "role": role, "model": model,
            "attenuated_at": row.get("attenuated_at"),
            "direction": _template_direction(row.get("positive_means"), role, name),
            "top_category": row.get("top_category"),
            # Optional, in vitro fitness scores only: `overrides` names another
            # declared column whose fitness-cost call this one outranks, and
            # `neutral_at` is the value at which this column positively
            # contradicts it. Used when one screen is measured in the target
            # organism and the other is transferred from a relative: a direct
            # measurement outranks a transfer, the same precedence the in vivo
            # side applies to infection models. Both must be given together.
            "neutral_at": row.get("neutral_at"),
            "overrides": _template_text(row.get("overrides")),
        })
        sources.append({
            "column": name, "field": role, "assay": str(row.get("assay") or ""),
            "citation": str(row.get("citation") or ""), "pmid": str(row.get("pmid") or ""),
        })
    if not declared:
        raise ValueError(
            f"{definition_sheet!r} declares no column that is present in {data_sheet!r}. "
            "At least one screen column is required."
        )

    n = len(out)
    false_series = lambda: pd.Series(False, index=out.index)
    essential, fitness_cost, in_vitro_assayed = false_series(), false_series(), false_series()
    attenuated, in_vivo_assayed = false_series(), false_series()
    basis = [[] for _ in range(n)]
    fit_hits: Dict[str, pd.Series] = {}
    fit_neutral: Dict[str, pd.Series] = {}
    effect = pd.Series(np.nan, index=out.index)
    effect_model = pd.Series("", index=out.index, dtype=object)
    effect_source = pd.Series("", index=out.index, dtype=object)
    rank = pd.Series("", index=out.index, dtype=object)
    signif = pd.Series(np.nan, index=out.index)
    models_seen = pd.Series("", index=out.index, dtype=object)

    for d in declared:
        raw = data[d["name"]]
        role, label = d["role"], d["name"]
        if role == "annotation_only":
            out[label] = raw
            continue
        if role == "in_vitro_essential_flag":
            hit = raw.notna() & (raw.astype(str).str.strip() != "")
            essential |= hit
            in_vitro_assayed |= hit
            for i in out.index[hit]:
                basis[out.index.get_loc(i)].append(f"{label}: essential")
            continue
        if role == "in_vitro_fitness_score":
            # Compared in the column's OWN units, on the side the declaration
            # names -- not in a sign-flipped space. An in vitro fitness score
            # need not be centred on zero: a Tn-seq fitness index and a
            # competitive index are ratios whose neutral value is 1.0, so
            # multiplying the value by the direction and testing `>= cut`
            # against a positive cut flags nothing at all, silently. The
            # declaration therefore sets the threshold where the user read it
            # in the paper (0.2 for a five-fold depleted fitness index, -2.0
            # for a log P) and the direction says which side is worse.
            num = pd.to_numeric(raw, errors="coerce")
            cut = pd.to_numeric(pd.Series([d["attenuated_at"]]), errors="coerce").iloc[0]
            in_vitro_assayed |= num.notna()
            if pd.isna(cut):
                notes.append(
                    f"Column {label!r} is a graded in vitro fitness score but declares no "
                    "attenuated_at threshold, so no fitness_cost call was made from it and "
                    "rung 6 cannot use it."
                )
            else:
                hit = num.notna() & (num <= cut if d["direction"] < 0 else num >= cut)
                fit_hits[label] = hit
            neutral_cut = pd.to_numeric(
                pd.Series([d.get("neutral_at")]), errors="coerce").iloc[0]
            if d.get("overrides"):
                if pd.isna(neutral_cut):
                    raise ValueError(
                        f"Column {label!r} declares overrides={d['overrides']!r} but no "
                        "neutral_at. A column may only override another where it "
                        "positively contradicts it, so the value at which it reads as "
                        "unimpaired has to be stated: a measurement that merely failed "
                        "the source's own cut-off is a weaker depletion, not a "
                        "contradiction."
                    )
                fit_neutral[label] = num.notna() & (
                    num >= neutral_cut if d["direction"] < 0 else num <= neutral_cut)
            continue
        if role == "in_vivo_attenuation_flag":
            hit = raw.notna() & (raw.astype(str).str.strip() != "")
            attenuated |= hit
            in_vivo_assayed |= hit
            for i in out.index[hit]:
                basis[out.index.get_loc(i)].append(f"{label} ({d['model']})")
                models_seen[i] = models_seen[i] or d["model"]
            continue
        if role == "in_vivo_effect_log2":
            num = pd.to_numeric(raw, errors="coerce") * d["direction"]
            in_vivo_assayed |= num.notna()
            cut = pd.to_numeric(pd.Series([d["attenuated_at"]]), errors="coerce").iloc[0]
            cut = 1.0 if pd.isna(cut) else cut
            hit = num.notna() & (num >= cut)
            attenuated |= hit
            # Model precedence decides which magnitude is reported when two
            # screens disagree; the loser stays visible in in_vivo_basis.
            better = num.notna() & (
                effect.isna()
                | (effect_model.map(lambda m: IN_VIVO_MODEL_RANK.get(m, 0))
                   < IN_VIVO_MODEL_RANK.get(d["model"], 0))
            )
            effect[better] = num[better]
            effect_model[better] = d["model"]
            effect_source[better] = f"{label} (positive = depleted in host)"
            for i in out.index[hit]:
                basis[out.index.get_loc(i)].append(f"{label} {num[i]:.2f} ({d['model']})")
                models_seen[i] = models_seen[i] or d["model"]
            continue
        if role == "in_vivo_ordinal":
            top = str(d["top_category"] or "").strip()
            text = raw.fillna("").astype(str).str.strip()
            in_vivo_assayed |= text != ""
            if not top:
                notes.append(
                    f"Column {label!r} is an ordinal grade but declares no top_category, "
                    "so it supplies no attenuation call."
                )
                continue
            # Compare numerically when BOTH the cell and top_category parse as
            # numbers: a spreadsheet will happily store the grade 2 as "2.0",
            # and a string comparison against "2" would then silently demote
            # every top-category gene to the intermediate rank.
            top_num = pd.to_numeric(pd.Series([top]), errors="coerce").iloc[0]
            cell_num = pd.to_numeric(text, errors="coerce")
            if pd.notna(top_num):
                strong = cell_num.notna() & (cell_num == top_num)
                strong |= (cell_num.isna()) & (text.str.lower() == top.lower())
            else:
                strong = text.str.lower() == top.lower()
            # An ordinal scale may have intermediate attenuated grades below its
            # top one (the TMDH scale runs -1..2 and calls 1 attenuated as well
            # as 2). `attenuated_at` sets where the attenuation call starts on a
            # NUMERIC scale; without it only the top category counts, which is
            # the right default for a scale whose other values are unordered
            # labels.
            cut = pd.to_numeric(pd.Series([d["attenuated_at"]]), errors="coerce").iloc[0]
            num = pd.to_numeric(text, errors="coerce")
            if pd.notna(cut) and num.notna().any():
                hit = num.notna() & (num >= cut)
                moderate = hit & ~strong
                rank[moderate] = "moderate"
            else:
                hit = strong
            attenuated |= hit
            rank[strong] = "strong"
            for i in out.index[hit]:
                basis[out.index.get_loc(i)].append(
                    f"{label} = {text[i]} ({d['model']})"
                )
                models_seen[i] = models_seen[i] or d["model"]
            continue
        if role == "in_vivo_significance":
            num = pd.to_numeric(raw, errors="coerce")
            signif = signif.where(signif.notna(), num)

    # Combine the graded in vitro columns, applying any declared precedence.
    # Done after the loop rather than inline, because an override can only be
    # resolved once both columns have been read.
    overridden = pd.Series(False, index=out.index)
    for d in declared:
        label = d["name"]
        if d["role"] != "in_vitro_fitness_score":
            continue
        # `overrides` may name SEVERAL columns, separated by semicolons. Two
        # screens of equal standing -- both measured natively, neither
        # transferred -- each name the other as well as the transferred column,
        # which makes the precedence symmetric: a cost call stands unless some
        # other declared column positively contradicts it.
        targets = [t.strip() for t in (d.get("overrides") or "").split(";") if t.strip()]
        for target in targets:
            if target not in fit_hits:
                notes.append(
                    f"Column {label!r} declares overrides={target!r}, which is not a "
                    "declared in vitro fitness column with a threshold; no precedence "
                    "was applied.")
                continue
            beats = fit_neutral.get(label, pd.Series(False, index=out.index)) & fit_hits[target]
            fit_hits[target] = fit_hits[target] & ~beats
            overridden |= beats
            for i in out.index[beats]:
                basis[out.index.get_loc(i)].append(
                    f"{target}: in vitro fitness cost, overridden by {label}")
    for label, hit in fit_hits.items():
        fitness_cost |= hit
        for i in out.index[hit]:
            basis[out.index.get_loc(i)].append(f"{label}: in vitro fitness cost")
    if overridden.any():
        notes.append(
            f"{int(overridden.sum())} gene(s) had a fitness-cost call overridden by a "
            "column declared to outrank it.")

    # A gene that was ASSAYED IN VIVO had a viable mutant, or there would have
    # been nothing to assay. That is a measurement of viability, not an
    # inference from annotation, so it licenses `dispensable` even with no in
    # vitro column at all -- and it is why a template carrying only in vivo
    # screens still produces a usable in vitro call for the genes it covers.
    requirement = pd.Series("unknown", index=out.index, dtype=object)
    requirement[in_vivo_assayed] = "dispensable"
    requirement[in_vitro_assayed] = "dispensable"
    requirement[fitness_cost] = "fitness_cost"
    requirement[essential] = "essential"
    out["in_vitro_requirement"] = requirement
    out["in_vitro_basis"] = ["; ".join(b) or "" for b in basis]

    fitness = pd.Series("not_assayed", index=out.index, dtype=object)
    fitness[in_vivo_assayed] = "no_phenotype"
    fitness[attenuated] = "attenuated"
    out["in_vivo_fitness"] = fitness
    out["in_vivo_basis"] = ["; ".join(b) or "" for b in basis]
    out["in_vivo_model"] = np.where(
        fitness == "not_assayed", "not_applicable",
        models_seen.replace("", "animal_infection"),
    )
    out["in_vivo_effect_log2"] = effect.where(fitness == "attenuated")
    out["in_vivo_effect_rank"] = rank.where(fitness == "attenuated", "")
    out["in_vivo_significance"] = signif
    out["in_vivo_effect_source"] = effect_source

    if all(d["role"] == "annotation_only" for d in declared):
        notes.append(
            "Every declared column has role 'annotation_only', so this table gates "
            "nothing: every gene reads 'unknown' in vitro and 'not_assayed' in vivo, "
            "and rungs 4-6 will be empty. Declare at least one in vitro or in vivo "
            "role to populate them."
        )
    if not any(d["role"].startswith("in_vivo") for d in declared):
        notes.append(
            "No in vivo column was declared, so rungs 4-6 will be empty: they require "
            "a measured phenotype in a host and are never filled by inference."
        )
    if not any(d["role"] == "in_vitro_fitness_score" for d in declared):
        notes.append(
            "No graded in vitro fitness column was declared, so rung 6's in vitro "
            "clause cannot fire and that rung rests entirely on its host-directedness "
            "clause. Declaring a quantitative in vitro fitness score is what makes the "
            "in vitro half of rung 6 mean anything."
        )
    out = validate_reference(out, source=os.path.basename(workbook_path))
    summary: Dict[str, object] = {
        "workbook": workbook_path,
        "sha256": file_checksum(workbook_path),
        "n_rows": int(len(out)),
        "n_columns_declared": len(declared),
        "n_essential": int((out["in_vitro_requirement"] == "essential").sum()),
        "n_fitness_cost": int((out["in_vitro_requirement"] == "fitness_cost").sum()),
        "n_dispensable": int((out["in_vitro_requirement"] == "dispensable").sum()),
        "n_attenuated": int((out["in_vivo_fitness"] == "attenuated").sum()),
        "n_no_phenotype": int((out["in_vivo_fitness"] == "no_phenotype").sum()),
        "n_not_assayed": int((out["in_vivo_fitness"] == "not_assayed").sum()),
        "in_vivo_models": {m: int((out["in_vivo_model"] == m).sum())
                           for m in sorted(set(out["in_vivo_model"]))},
        "graded_in_vitro_fitness_present": bool((out["in_vitro_requirement"] == "fitness_cost").any()),
        "n_in_vitro_overridden": int(overridden.sum()),
        "notes": notes,
    }
    if out_path:
        os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
        out.to_csv(out_path, sep="\t", index=False)
        pd.DataFrame(sources).to_csv(
            os.path.splitext(out_path)[0] + "_sources.tsv", sep="\t", index=False
        )
        summary["out_path"] = out_path
    summary["frame"] = out
    return summary
