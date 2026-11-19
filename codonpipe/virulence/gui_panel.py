#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The 'Automated virulence gene discovery' panel for the CodonPipe interface.

The panel is deliberately self-contained and independent of the rest of the
interface: it does not read or write any other CodonPipe setting, and the
discovery run neither needs nor touches a clustering result.  The only point of
contact is optional and explicit -- a button that copies the resulting cluster
workbook into the *Gene clusters* field once the user has looked at it.

``attach_virulence_panel`` is written to be impossible to break the host
interface with.  It is called inside a ``try`` block from
``CodonPipe_GUI._build_input_tab``, catches every exception while building, and
runs the discovery in a worker thread so the window stays responsive.  If this
module is missing or raises, CodonPipe behaves exactly as it did before.
"""

from __future__ import annotations

import os
import queue
import threading
import traceback
import webbrowser
from typing import Callable, Optional

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

__all__ = ["attach_virulence_panel", "VirulencePanel"]

#: Shown in the panel. The third entry names the tier both ways on purpose:
#: the dropdown below it is labelled "Host-interaction (strict) tier requires",
#: and a reader must be able to see at a glance that the two refer to the same
#: set of genes.
STRINGENCY_CHOICES = [
    "2 standard",
    "1 permissive",
    "3 host-interaction",
]
_STRINGENCY_TO_INTERNAL = {
    "2 standard": "standard",
    "1 permissive": "permissive",
    "3 host-interaction": "strict",
}
CATEGORY_CHOICES = [
    "Published tiers 1-4",
    "All ladder rungs (tiers + both gates)",
    "curation gate + tier 2",
    "tier 2 only (host-interaction)",
    "curation gate only",
    "in vivo evidence only (gate 4, tiers 3 and 4)",
    "tier 3 only (in vivo required)",
]
_CATEGORY_TO_INTERNAL = {
    "Published tiers 1-4": "published",
    "All ladder rungs (tiers + both gates)": "all",
    "curation gate + tier 2": "standard,host_interaction",
    "tier 2 only (host-interaction)": "host_interaction",
    "curation gate only": "standard",
    "in vivo evidence only (gate 4, tiers 3 and 4)": (
        "in_vivo_contributing,in_vivo_required,in_vitro_neutral"
    ),
    "tier 3 only (in vivo required)": "in_vivo_required",
}

#: Defaults for the two effect-size cuts and the significance cut, exposed in
#: the panel because they are reporting CONVENTIONS rather than measurements.
#: 1.0 on a log2 scale is a two-fold depletion, the line below which a
#: transposon-sequencing defect is usually not called; 3.32 is log2(10), the
#: ten-fold competitive index that the Salmonella literature has long treated
#: as the boundary of attenuation. Both travel into every output file so a run
#: can be read without knowing what was typed here.
DEFAULT_EFFECT_CONTRIBUTING = "0.0"
DEFAULT_EFFECT_REQUIRED = "1.0"
DEFAULT_MAX_SIGNIFICANCE = "0.05"

#: Absolute log2 cut-offs are not comparable between assay types: a pooled
#: transposon library is bottlenecked, so its log fold changes run far larger
#: than a competitive index of the same biological severity. Switching the
#: scale to percentile reads the same two boxes against this organism's own
#: distribution of measured effects, which is what makes a rung 5 count mean
#: the same thing in two genomes.
EFFECT_SCALE_CHOICES = ["absolute log2 depletion", "percentile within organism"]
_EFFECT_SCALE_TO_INTERNAL = {
    "absolute log2 depletion": "absolute",
    "percentile within organism": "percentile",
}
_EFFECT_SCALE_DEFAULTS = {
    "absolute": (DEFAULT_EFFECT_CONTRIBUTING, DEFAULT_EFFECT_REQUIRED),
    "percentile": ("0", "90"),
}

_INTRO = (
    "Identify virulence genes automatically from a species name or NCBI identifier. "
    "This tool is independent of the rest of the interface: it needs no clustering "
    "result and changes no other setting. One run returns a strictly nested ladder "
    "of three gene sets defined by Falkow's molecular Koch's postulates - "
    "permissive (any virulence evidence), standard (evidence anchored in database "
    "curation) and strict (the host-interaction set, where loss of virulence on "
    "specific inactivation cannot be explained by a generic growth defect)."
)

_TIP_TARGET = (
    "Accepts a species or strain name (Salmonella Typhimurium SL1344, "
    "Klebsiella pneumoniae), an NCBI taxonomy identifier (216597), or a genome "
    "assembly accession (GCF_000210855.2). Assembly accessions are the most "
    "precise because they pin the exact annotation used."
)
_TIP_STRINGENCY = (
    "Which tier the run highlights and which tier supplies the functional "
    "subcategory columns. All rungs are computed either way. standard is the "
    "default: evidence anchored in database curation. permissive adds "
    "annotation-only candidates and is the set to quote when asked whether a "
    "result survives the most inclusive possible definition.\n\n"
    "Every run writes the whole ladder, which is a TREE, not a single chain. "
    "Rungs 1 to 3 are built from annotation and curation; rungs 4 to 6 are "
    "built from published mutant screens and are read independently, because "
    "a measured phenotype outranks an inference drawn from annotation.\n\n"
    "  annotation branch\n"
    "  1 permissive          any evidence, annotation alone included\n"
    "  2 standard            + anchored in curation\n"
    "  3 host-interaction    + host-directed, not confounded by core processes\n\n"
    "  experimental branch\n"
    "  4 in vivo contributing  a significant, non-trivial defect in a host,\n"
    "                          and not essential in vitro\n"
    "  5 in vivo required      + the defect is large (ten-fold by default)\n"
    "  6 in vivo only          + no measured fitness cost in vitro\n\n"
    "Each branch nests internally. Rung 4 is deliberately NOT inside rung 3: a "
    "gene with a measured in vivo defect belongs in the experimental branch "
    "whether or not annotation had flagged it as host-directed. The report "
    "states how many genes the two branches share.\n\n"
    "Rungs 4 to 6 need published mutant-screen data for the organism and are "
    "absent when none exists. Rung 5 additionally needs per-gene effect sizes, "
    "and rung 6 graded in vitro fitness data; where those are missing the rung "
    "is reported as empty rather than filled by assumption.\n\n"
    "'Host-interaction' and 'strict' are two names for rung 3."
)
_TIP_THRESHOLDS = (
    "Where to draw the line between a gene that CONTRIBUTES to virulence and "
    "one that is REQUIRED for it. Both are log2 depletion of the mutant in the "
    "host: 1.0 is two-fold, 2.0 four-fold, 3.32 ten-fold.\n\n"
    "The defaults are the thresholds the bundled Salmonella compilation "
    "applies to its own data: any positive measured defect for contributing, "
    "and two-fold for required.\n\n"
    "A ten-fold line was tried and is wrong for this biology. The SPI-1 "
    "competitive indices in that compilation run from 0.4 to 2.1 log2, because "
    "its screens are systemic infection models in which SPI-1 mutants are only "
    "mildly attenuated, so a ten-fold bar removed the whole SPI-1 apparatus "
    "from the required rung while keeping SPI-2. At two-fold the required set "
    "contains the SPI-2 apparatus, the core of SPI-1, phoP, phoQ, ompR, envZ, "
    "crp, fur and the O-antigen cluster.\n\n"
    "Lowering the required cut-off moves genes from rung 4 into rung 5; "
    "raising it does the reverse. Nothing is discarded by changing them. The "
    "values used are written into the summary sheet of every run.\n\n"
    "The contributing cut-off must not exceed the required one."
)
_TIP_EFFECT_SCALE = (
    "How the two boxes to the left are read.\n\n"
    "absolute log2 depletion -- the cut-offs are fold changes: 1.0 is "
    "two-fold, 3.32 ten-fold. Use this when effect sizes come from "
    "competitive indices, and when you want a threshold with a literature "
    "meaning.\n\n"
    "percentile within organism -- the cut-offs are percentiles of the "
    "measured in vivo effect sizes in THIS organism, so 90 means the "
    "strongest tenth. Use this when comparing organisms.\n\n"
    "The distinction is not cosmetic. A pooled transposon library is "
    "bottlenecked, so its log fold changes run much larger than a competitive "
    "index of the same biological severity. A ten-fold line takes the "
    "strongest 7% of the Salmonella competitive indices but roughly 75% of "
    "the transposon-screen values, which would make a rung 5 count mean "
    "different things in different genomes. On the percentile scale the same "
    "setting takes 5 to 10% in every organism tested.\n\n"
    "Switching the scale resets the two boxes to that scale's defaults. The "
    "scale, the settings and the resolved log2 values are all written to the "
    "run summary."
)
_TIP_SIGNIFICANCE = (
    "The largest p or q value accepted as a real in vivo defect. Applied only "
    "where the underlying screen reports a per-gene value; where it does not, "
    "the screen's own published attenuation call stands, rather than the gene "
    "being dropped for lacking a number its authors never printed.\n\n"
    "A gene that clears significance but has no reported effect size stops at "
    "rung 4: it is never promoted to 'required' on the strength of silence."
)
_TIP_OFFLINE = (
    "Replay a previous run entirely from the cache folder, without any network "
    "access. Use this to reproduce a published result, or to re-run after the "
    "first successful run."
)
_TIP_ESSENTIALITY = (
    "Optional table of experimentally measured gene requirements (transposon "
    "insertion sequencing, essential-gene screens, animal infection screens). "
    "Needs a locus_tag column, plus in_vitro_requirement (essential, "
    "fitness_cost, dispensable) and/or in_vivo_fitness (attenuated, "
    "no_phenotype, not_assayed). Where supplied it overrides the "
    "annotation-based inference for Falkow's second and third postulates, and "
    "it takes precedence over the reference table bundled for this organism."
)


class VirulencePanel:
    """Build and drive the discovery panel inside an existing CodonPipe GUI."""

    def __init__(self, app, parent):
        self.app = app
        self.parent = parent
        self.root = getattr(app, "root", None)
        self._queue: "queue.Queue[str]" = queue.Queue()
        self._thread: Optional[threading.Thread] = None
        self._result = None
        self._last_outputs: dict = {}

        self.target = tk.StringVar(value="")
        self.stringency = tk.StringVar(value=STRINGENCY_CHOICES[0])
        self.categories = tk.StringVar(value=CATEGORY_CHOICES[0])
        self.output_dir = tk.StringVar(value="")
        self.offline = tk.BooleanVar(value=False)
        self.refresh = tk.BooleanVar(value=False)
        self.use_selected_genome = tk.BooleanVar(value=True)
        self.species_transfer = tk.BooleanVar(value=True)
        self.essentiality_table = tk.StringVar(value="")
        self.effect_contributing = tk.StringVar(value=DEFAULT_EFFECT_CONTRIBUTING)
        self.effect_required = tk.StringVar(value=DEFAULT_EFFECT_REQUIRED)
        self.max_significance = tk.StringVar(value=DEFAULT_MAX_SIGNIFICANCE)
        self.effect_scale = tk.StringVar(value=EFFECT_SCALE_CHOICES[0])
        self.ncbi_email = tk.StringVar(value="")
        self.status = tk.StringVar(value="Not run yet.")

        self._build()
        self._poll()

    # ----------------------------------------------------------------- build
    def _tip(self, widget, text):
        helper = getattr(self.app, "_tip", None)
        if callable(helper):
            try:
                helper(widget, text)
            except Exception:
                pass

    def _build(self) -> None:
        box = ttk.LabelFrame(
            self.parent, text="Automated virulence gene discovery (independent)", padding=10
        )
        box.pack(fill="x", pady=6)
        box.columnconfigure(1, weight=1)
        self.box = box

        ttk.Label(box, text=_INTRO, foreground="#444444", wraplength=760, justify="left").grid(
            row=0, column=0, columnspan=4, sticky="w", pady=(0, 8)
        )

        label = ttk.Label(box, text="Species or NCBI identifier")
        label.grid(row=1, column=0, sticky="w", padx=(0, 8), pady=3)
        self._tip(label, _TIP_TARGET)
        entry = ttk.Entry(box, textvariable=self.target)
        entry.grid(row=1, column=1, columnspan=2, sticky="ew", pady=3)
        self._tip(entry, _TIP_TARGET)
        button = ttk.Button(box, text="Use selected genome", command=self._fill_from_selected_genome)
        button.grid(row=1, column=3, sticky="w", padx=(8, 0), pady=3)
        self._tip(
            button,
            "Read the assembly accession out of the genome currently chosen in "
            "'Upload genome' above, by looking up its replicon accession at NCBI.",
        )

        check = ttk.Checkbutton(
            box,
            text="Analyse the coding sequences of the genome selected above",
            variable=self.use_selected_genome,
        )
        check.grid(row=2, column=0, columnspan=4, sticky="w", pady=(2, 4))
        self._tip(
            check,
            "When enabled, the gene set is reported with the locus tags of the FASTA "
            "already loaded in CodonPipe, so the cluster workbook drops straight into "
            "the Gene clusters field. When disabled, the RefSeq coding sequences of "
            "the resolved assembly are downloaded instead.",
        )

        label = ttk.Label(box, text="Headline tier")
        label.grid(row=3, column=0, sticky="w", padx=(0, 8), pady=3)
        self._tip(label, _TIP_STRINGENCY)
        combo = ttk.Combobox(
            box, textvariable=self.stringency, state="readonly",
            values=STRINGENCY_CHOICES, width=18,
        )
        combo.grid(row=3, column=1, sticky="w", pady=3)
        self._tip(combo, _TIP_STRINGENCY)

        label = ttk.Label(box, text="Tiers to write out")
        label.grid(row=3, column=2, sticky="e", padx=(8, 8), pady=3)
        combo = ttk.Combobox(
            box, textvariable=self.categories, state="readonly",
            values=CATEGORY_CHOICES, width=26,
        )
        combo.grid(row=3, column=3, sticky="w", pady=3)
        self._tip(
            combo,
            "All three tiers are always computed; this chooses which become columns "
            "in the cluster workbook. The default writes all three so that the "
            "sensitivity of a downstream analysis to the virulence-gene definition "
            "can be read off directly.",
        )

        label = ttk.Label(box, text="Output folder")
        label.grid(row=4, column=0, sticky="w", padx=(0, 8), pady=3)
        entry = ttk.Entry(box, textvariable=self.output_dir)
        entry.grid(row=4, column=1, columnspan=2, sticky="ew", pady=3)
        self._tip(
            entry,
            "Where the cluster workbook, the evidence audit workbook, the provenance "
            "report and the run manifest are written. Defaults to a "
            "'Virulence discovery' folder inside the default export root.",
        )
        ttk.Button(box, text="Browse", command=self._browse_output_dir).grid(
            row=4, column=3, sticky="w", padx=(8, 0), pady=3
        )

        label = ttk.Label(box, text="In vitro fitness table (optional)")
        label.grid(row=5, column=0, sticky="w", padx=(0, 8), pady=3)
        self._tip(label, _TIP_ESSENTIALITY)
        entry = ttk.Entry(box, textvariable=self.essentiality_table)
        entry.grid(row=5, column=1, columnspan=2, sticky="ew", pady=3)
        self._tip(entry, _TIP_ESSENTIALITY)
        ttk.Button(box, text="Browse", command=self._browse_essentiality).grid(
            row=5, column=3, sticky="w", padx=(8, 0), pady=3
        )

        label = ttk.Label(box, text="In vivo effect thresholds:")
        label.grid(row=7, column=0, sticky="w", padx=(0, 8), pady=3)
        self._tip(label, _TIP_THRESHOLDS)
        cuts = ttk.Frame(box)
        cuts.grid(row=7, column=1, columnspan=3, sticky="w", pady=3)
        for text, var, width, tip in (
            ("contributing \u2265", self.effect_contributing, 6, _TIP_THRESHOLDS),
            ("required \u2265", self.effect_required, 6, _TIP_THRESHOLDS),
            ("p/q \u2264", self.max_significance, 6, _TIP_SIGNIFICANCE),
        ):
            sub = ttk.Label(cuts, text=text)
            sub.pack(side="left", padx=(0, 3))
            self._tip(sub, tip)
            entry = ttk.Entry(cuts, textvariable=var, width=width, justify="right")
            entry.pack(side="left", padx=(0, 12))
            self._tip(entry, tip)
        scale = ttk.Combobox(
            cuts, textvariable=self.effect_scale, state="readonly",
            values=EFFECT_SCALE_CHOICES, width=24,
        )
        scale.pack(side="left", padx=(0, 10))
        scale.bind("<<ComboboxSelected>>", lambda _e: self._reset_thresholds())
        self._tip(scale, _TIP_EFFECT_SCALE)

        reset = ttk.Button(cuts, text="Defaults", width=9, command=self._reset_thresholds)
        reset.pack(side="left")
        self._tip(
            reset,
            "Restore the published conventions: two-fold (1.0) for a contributing "
            "gene, ten-fold (3.32) for a required one, 0.05 for significance.",
        )

        options = ttk.Frame(box)
        options.grid(row=8, column=0, columnspan=4, sticky="w", pady=(4, 2))
        check = ttk.Checkbutton(options, text="Offline (use cache only)", variable=self.offline)
        check.pack(side="left", padx=(0, 14))
        self._tip(check, _TIP_OFFLINE)
        check = ttk.Checkbutton(options, text="Force refresh of cached records", variable=self.refresh)
        check.pack(side="left", padx=(0, 14))
        self._tip(
            check,
            "Re-download every record even if it is cached. Use when a database has "
            "been updated since the last run.",
        )
        check = ttk.Checkbutton(
            options, text="Transfer curation from the parent species", variable=self.species_transfer
        )
        check.pack(side="left")
        self._tip(
            check,
            "UniProtKB curation is concentrated in a few reference strains. When "
            "enabled, reviewed entries of the parent species are matched by gene "
            "symbol and marked as a cross-strain transfer in the audit table.",
        )

        label = ttk.Label(box, text="Contact e-mail for NCBI (optional)")
        label.grid(row=8, column=0, sticky="w", padx=(0, 8), pady=3)
        entry = ttk.Entry(box, textvariable=self.ncbi_email)
        entry.grid(row=8, column=1, sticky="ew", pady=3)
        self._tip(
            entry,
            "NCBI asks automated clients to identify themselves. Supplying an address "
            "is optional; it is sent only to NCBI Entrez and is stored nowhere else.",
        )

        actions = ttk.Frame(box)
        actions.grid(row=9, column=0, columnspan=4, sticky="w", pady=(10, 2))
        self.run_button = ttk.Button(
            actions, text="Identify virulence genes", command=self.run
        )
        self.run_button.pack(side="left", padx=(0, 8))
        self._tip(
            self.run_button,
            "Run the discovery pipeline. The first run for a genome fetches records "
            "from NCBI, BV-BRC, UniProtKB and QuickGO and takes roughly half a "
            "minute; later runs replay from the cache.",
        )
        self.use_button = ttk.Button(
            actions, text="Use as cluster file", command=self._use_as_cluster_file, state="disabled"
        )
        self.use_button.pack(side="left", padx=(0, 8))
        self._tip(
            self.use_button,
            "Copy the resulting cluster workbook into the 'Gene clusters' field above "
            "and switch the cluster source to the user-provided workbook. Nothing is "
            "changed until you press this.",
        )
        self.merge_button = ttk.Button(
            actions, text="Merge into cluster file", command=self._merge_into_cluster_file,
            state="disabled",
        )
        self.merge_button.pack(side="left", padx=(0, 8))
        self._tip(
            self.merge_button,
            "Write a new workbook holding your currently selected gene clusters "
            "PLUS the virulence columns from this run, and select it. Both "
            "originals are left untouched, and a column name present in both "
            "workbooks is kept separate rather than merged, because the same "
            "heading in two workbooks is not necessarily the same gene set.",
        )
        self.open_button = ttk.Button(
            actions, text="Open output folder", command=self._open_output_folder, state="disabled"
        )
        self.open_button.pack(side="left", padx=(0, 8))
        self.report_button = ttk.Button(
            actions, text="Open report", command=self._open_report, state="disabled"
        )
        self.report_button.pack(side="left")

        ttk.Label(box, textvariable=self.status, foreground="#444444", wraplength=760).grid(
            row=10, column=0, columnspan=4, sticky="w", pady=(6, 0)
        )

    # --------------------------------------------------------------- helpers
    def _log(self, message: str) -> None:
        self._queue.put(str(message).rstrip() + "\n")

    def _poll(self) -> None:
        try:
            while True:
                text = self._queue.get_nowait()
                appender = getattr(self.app, "_append_log", None)
                if callable(appender):
                    appender(text)
        except queue.Empty:
            pass
        except Exception:
            pass
        if self.root is not None:
            try:
                self.root.after(150, self._poll)
            except Exception:
                pass

    def _default_output_dir(self) -> str:
        root_var = getattr(self.app, "default_root", None)
        base = ""
        try:
            base = (root_var.get() or "").strip() if root_var is not None else ""
        except Exception:
            base = ""
        if not base:
            base = os.getcwd()
        return os.path.join(base, "Virulence discovery")

    def _selected_genome_fasta(self) -> str:
        getter = getattr(self.app, "_current_genome_fasta_path", None)
        if callable(getter):
            try:
                path = getter()
                if path and os.path.isfile(path):
                    return path
            except Exception:
                pass
        try:
            path = (self.app.fasta_path.get() or "").strip()
            if path and os.path.isfile(path):
                return path
        except Exception:
            pass
        return ""

    def _browse_output_dir(self) -> None:
        initial = self.output_dir.get().strip() or self._default_output_dir()
        chosen = filedialog.askdirectory(title="Virulence discovery output folder", initialdir=initial)
        if chosen:
            self.output_dir.set(chosen)

    def _browse_essentiality(self) -> None:
        chosen = filedialog.askopenfilename(
            title="Experimental in vitro fitness table",
            filetypes=[("Tables", "*.xlsx *.xls *.csv *.tsv *.txt"), ("All files", "*.*")],
        )
        if chosen:
            self.essentiality_table.set(chosen)

    def _fill_from_selected_genome(self) -> None:
        fasta = self._selected_genome_fasta()
        if not fasta:
            messagebox.showinfo(
                "Virulence discovery",
                "No genome FASTA is currently selected in the 'Upload genome' section.",
            )
            return
        self.status.set("Looking up the assembly for the selected genome ...")
        try:
            from .resolve import assembly_from_cds_fasta
            from .webcache import WebCache

            web = WebCache(
                self._cache_dir(), offline=self.offline.get(), log=self._log
            )
            accession, notes = assembly_from_cds_fasta(
                fasta, web, ncbi_email=self.ncbi_email.get().strip()
            )
            for note in notes:
                self._log(f"[virulence] {note}")
            if accession:
                self.target.set(accession)
                self.status.set(f"Assembly {accession} taken from the selected genome.")
            else:
                self.status.set(
                    "Could not identify the assembly from the FASTA headers; "
                    "type the species name instead."
                )
        except Exception as exc:
            self._log("[virulence] assembly lookup failed:\n" + traceback.format_exc())
            self.status.set(f"Assembly lookup failed: {exc}")

    def _cache_dir(self) -> str:
        out = self.output_dir.get().strip() or self._default_output_dir()
        return os.path.join(out, "virulence_cache")

    # ------------------------------------------------- in vivo effect cut-offs
    def _current_effect_scale(self) -> str:
        return _EFFECT_SCALE_TO_INTERNAL.get(self.effect_scale.get(), "absolute")

    def _reset_thresholds(self) -> None:
        """Put the cut-offs back to the defaults for the selected scale."""
        contributing, required = _EFFECT_SCALE_DEFAULTS[self._current_effect_scale()]
        self.effect_contributing.set(contributing)
        self.effect_required.set(required)
        self.max_significance.set(DEFAULT_MAX_SIGNIFICANCE)

    def _read_thresholds(self) -> dict:
        """Validate the three boxes and return them as scoring overrides.

        Raises ``ValueError`` with a message meant for the user. The ordering
        rule matters: if the contributing cut were allowed above the required
        one, rung 5 would stop being a subset of rung 4 and the ladder would
        no longer nest, so it is refused here rather than caught downstream.
        """
        fields = (
            ("contributing effect size", self.effect_contributing, "effect_contributing_log2"),
            ("required effect size", self.effect_required, "effect_required_log2"),
            ("significance", self.max_significance, "max_in_vivo_significance"),
        )
        values = {}
        for label, var, key in fields:
            raw = (var.get() or "").strip()
            try:
                value = float(raw)
            except ValueError:
                raise ValueError(f"The {label} box must contain a number, not {raw!r}.")
            if value < 0:
                raise ValueError(f"The {label} cut-off cannot be negative.")
            values[key] = value
        values["effect_scale"] = self._current_effect_scale()
        if values["effect_scale"] == "percentile" and values["effect_required_log2"] > 100:
            raise ValueError(
                "On the percentile scale the cut-offs must lie between 0 and 100."
            )
        if values["max_in_vivo_significance"] > 1:
            raise ValueError(
                "The significance cut-off is a p or q value and must lie between 0 and 1."
            )
        if values["effect_contributing_log2"] > values["effect_required_log2"]:
            raise ValueError(
                "The contributing effect size must be smaller than or equal to the "
                "required one: rung 5 is the subset of rung 4 that clears the higher "
                "bar, so reversing them would break the ladder."
            )
        return values

    # ------------------------------------------------------------------- run
    def run(self) -> None:
        if self._thread is not None and self._thread.is_alive():
            messagebox.showinfo("Virulence discovery", "A discovery run is already in progress.")
            return
        target = self.target.get().strip()
        if not target:
            messagebox.showerror(
                "Virulence discovery",
                "Enter a species name, an NCBI taxonomy identifier or a genome "
                "assembly accession, or press 'Use selected genome'.",
            )
            return

        try:
            thresholds = self._read_thresholds()
        except ValueError as exc:
            messagebox.showerror("Virulence discovery", str(exc))
            return

        out_dir = self.output_dir.get().strip() or self._default_output_dir()
        self.output_dir.set(out_dir)
        cds = self._selected_genome_fasta() if self.use_selected_genome.get() else ""
        prefix = ""
        try:
            label = (self.app.preloaded_genome_choice.get() or "").strip()
            if cds and label:
                prefix = label
        except Exception:
            prefix = ""

        self.run_button.configure(state="disabled")
        self.use_button.configure(state="disabled")
        self.open_button.configure(state="disabled")
        self.report_button.configure(state="disabled")
        self.status.set("Running ...")
        self._log("\n" + "=" * 84)
        self._log(f"[virulence] automated virulence gene discovery - target: {target}")

        def worker() -> None:
            try:
                from .outputs import safe_name
                from .pipeline import DiscoveryConfig, run_discovery

                config = DiscoveryConfig(
                    target=target,
                    cds_fasta=cds,
                    stringency=_STRINGENCY_TO_INTERNAL.get(
                        self.stringency.get().strip(), "standard"
                    ),
                    tiers=_CATEGORY_TO_INTERNAL.get(self.categories.get(), "all"),
                    output_dir=out_dir,
                    offline=bool(self.offline.get()),
                    refresh=bool(self.refresh.get()),
                    include_species_transfer=bool(self.species_transfer.get()),
                    essentiality_table=self.essentiality_table.get().strip(),
                    scoring_overrides=thresholds,
                    ncbi_email=self.ncbi_email.get().strip(),
                    output_prefix=safe_name(prefix) if prefix else "",
                )
                result = run_discovery(config, log=self._log, write_outputs=True)
                self._result = result
                self._last_outputs = dict(result.outputs)
                counts = ", ".join(f"{k}: {len(v)}" for k, v in result.columns.items())
                self._finish(
                    True,
                    f"{result.target.organism_label or result.target.scientific_name} - {counts}",
                )
            except Exception as exc:
                self._log("[virulence] run failed:\n" + traceback.format_exc())
                self._finish(False, str(exc))

        self._thread = threading.Thread(target=worker, name="virulence-discovery", daemon=True)
        self._thread.start()

    def _finish(self, ok: bool, message: str) -> None:
        # Tk must only be touched from the main thread, so the worker hands the
        # update back through root.after().
        def _apply() -> None:
            try:
                self.run_button.configure(state="normal")
                if ok:
                    self.status.set(message)
                    self.use_button.configure(state="normal")
                    self.merge_button.configure(state="normal")
                    self.open_button.configure(state="normal")
                    self.report_button.configure(state="normal")
                    # Nothing above this panel has been touched. Say so, and
                    # say what the workbook does and does not contain, because
                    # a file called "clusters" that holds only virulence
                    # columns is exactly the kind of thing that gets handed to
                    # the clustering step by mistake.
                    messagebox.showinfo(
                        "Virulence discovery",
                        "Discovery completed.\n\n"
                        + message
                        + "\n\nNothing in the interface above has been changed. The "
                        "cluster workbook just written contains ONLY the virulence "
                        "columns from this run; it does not contain any of your own "
                        "gene clusters.\n\n"
                        "To use these genes in the codon-usage pipeline, merge them "
                        "into your own cluster workbook first -- press 'Merge into "
                        "cluster file' to write a combined workbook, which leaves "
                        "both originals untouched. 'Use as cluster file' REPLACES "
                        "the current selection with the virulence-only workbook and "
                        "should be used only when you want the virulence columns "
                        "alone.",
                    )
                else:
                    self.status.set(f"Failed: {message}")
                    messagebox.showerror(
                        "Virulence discovery",
                        "The discovery run failed. See the console log for details.\n\n" + message,
                    )
            except Exception:
                pass

        if self.root is not None:
            try:
                self.root.after(0, _apply)
                return
            except Exception:
                pass
        _apply()

    # --------------------------------------------------------------- results
    def _merge_into_cluster_file(self) -> None:
        """Combine the user's own cluster workbook with this run's columns."""
        path = (self._last_outputs or {}).get("clusters_xlsx", "")
        if not path or not os.path.isfile(path):
            messagebox.showinfo("Virulence discovery", "No cluster workbook is available yet.")
            return
        current = self._current_cluster_file()
        if not current:
            messagebox.showinfo(
                "Virulence discovery",
                "No cluster workbook is selected in the interface above, so "
                "there is nothing to merge with.\n\nSelect your own cluster "
                "workbook first, then press this again. If the virulence "
                "columns are the only clusters you want, use 'Use as cluster "
                "file' instead.",
            )
            return
        try:
            from . import outputs as _outputs

            result = _outputs.merge_cluster_workbooks(current, path)
        except Exception as exc:
            messagebox.showerror(
                "Virulence discovery", f"Could not merge the cluster workbooks: {exc}"
            )
            return

        renamed = result.get("renamed") or {}
        note = ""
        if renamed:
            shown = "; ".join(f"{k} -> {v}" for k, v in list(renamed.items())[:6])
            note = (
                f"\n\n{len(renamed)} column name(s) appear in both workbooks and were "
                f"kept separate rather than combined, because the same heading in two "
                f"workbooks is not necessarily the same gene set:\n    {shown}"
            )
        try:
            setter = getattr(self.app, "_set_refined_cluster_file", None)
            if callable(setter):
                setter(result["path"])
            else:
                self.app.refined_cluster_file.set(result["path"])
            try:
                self.app.user_cluster_mode.set(
                    "Provided by user, xlsx file with 1 column per cluster"
                )
                handler = getattr(self.app, "_on_cluster_mode_changed", None)
                if callable(handler):
                    handler()
            except Exception:
                pass
        except Exception as exc:
            messagebox.showerror(
                "Virulence discovery",
                f"The merged workbook was written to\n\n    {result['path']}\n\n"
                f"but it could not be selected automatically ({exc}). Point the "
                f"'Gene clusters' field at it by hand.",
            )
            return

        self.status.set(
            f"Merged cluster file set to {os.path.basename(result['path'])} "
            f"({result['n_columns_total']} columns)."
        )
        self._log(f"[virulence] merged cluster workbook written to {result['path']}")
        messagebox.showinfo(
            "Virulence discovery",
            "A merged cluster workbook has been written and selected:\n\n"
            f"    {os.path.basename(result['path'])}\n\n"
            f"{result['n_columns_base']} column(s) from your own workbook plus "
            f"{result['n_columns_added']} from this discovery run, "
            f"{result['n_columns_total']} in total.\n\nNeither original file was "
            "modified." + note,
        )

    def _current_cluster_file(self) -> str:
        """Whatever cluster workbook the interface above is pointing at."""
        try:
            value = self.app.refined_cluster_file.get()
        except Exception:
            return ""
        value = str(value or "").strip()
        return value if value and os.path.isfile(value) else ""

    def _use_as_cluster_file(self) -> None:
        path = (self._last_outputs or {}).get("clusters_xlsx", "")
        if not path or not os.path.isfile(path):
            messagebox.showinfo("Virulence discovery", "No cluster workbook is available yet.")
            return
        # Replacing a selection the user made themselves is destructive to
        # their intent even though no file is deleted, so it is confirmed
        # explicitly and the alternative is named.
        current = self._current_cluster_file()
        if current and os.path.abspath(current) != os.path.abspath(path):
            proceed = messagebox.askokcancel(
                "Replace the selected cluster file?",
                "A cluster workbook is already selected:\n\n"
                f"    {os.path.basename(current)}\n\n"
                "Pressing OK REPLACES it with the virulence-only workbook\n\n"
                f"    {os.path.basename(path)}\n\n"
                "The two are NOT merged, and every cluster defined in your own "
                "workbook will be absent from the run.\n\n"
                "If you meant to keep your own clusters and add the virulence "
                "columns to them, press Cancel and use 'Merge into cluster "
                "file' instead.\n\n"
                "Neither file on disk is modified either way.",
                icon="warning",
            )
            if not proceed:
                return
        try:
            setter = getattr(self.app, "_set_refined_cluster_file", None)
            if callable(setter):
                setter(path)
            else:
                self.app.refined_cluster_file.set(path)
            try:
                self.app.user_cluster_mode.set(
                    "Provided by user, xlsx file with 1 column per cluster"
                )
                handler = getattr(self.app, "_on_cluster_mode_changed", None)
                if callable(handler):
                    handler()
            except Exception:
                pass
            self.status.set(f"Cluster file set to {os.path.basename(path)}.")
            self._log(f"[virulence] cluster file set to {path}")
        except Exception as exc:
            messagebox.showerror("Virulence discovery", f"Could not set the cluster file: {exc}")

    def _open_output_folder(self) -> None:
        folder = self.output_dir.get().strip()
        if folder and os.path.isdir(folder):
            try:
                os.startfile(folder)  # noqa: S606 - Windows shell open
            except Exception:
                webbrowser.open(folder)

    def _open_report(self) -> None:
        path = (self._last_outputs or {}).get("report_txt", "")
        if path and os.path.isfile(path):
            try:
                os.startfile(path)  # noqa: S606 - Windows shell open
            except Exception:
                webbrowser.open(path)


def attach_virulence_panel(app, parent) -> Optional[VirulencePanel]:
    """Attach the panel to ``parent``; never raise into the host interface."""
    try:
        return VirulencePanel(app, parent)
    except Exception:
        try:
            import sys

            print(
                "[virulence] the discovery panel could not be built:\n"
                + traceback.format_exc(),
                file=sys.stderr,
                flush=True,
            )
        except Exception:
            pass
        return None
