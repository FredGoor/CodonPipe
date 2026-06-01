"""codonpipe.density_bridge

Utility to execute Plotting_Pipeline.py programmatically.
"""

import os
import sys
import uuid
import importlib.util
import traceback


def run_density_plot_script(excel_path, pipeline_cfg):
    """
    Executes the plotting pipeline script (e.g. Plotting_Pipeline.py) as a standalone module.

    The bridge wires the workbook/sheet names and can optionally override selected
    plotting globals from ``pipeline_cfg`` so a GUI can drive the figure settings
    without editing Plotting_Pipeline.py directly.
    """
    script_path = pipeline_cfg.get("plotting_pipeline_script_path", "")
    if not script_path:
        print("[WARN] No density_plot_script_path specified; skipping plotting pipeline.")
        return

    script_path = os.path.abspath(os.path.expanduser(script_path))
    if not os.path.exists(script_path):
        print(f"[WARN] Plotting pipeline script not found at:\n  {script_path}\nSkipping.")
        return

    script_name = os.path.basename(script_path)
    print(f"[INFO] Running plotting pipeline: {script_name}")

    module_name = f"_codonpipe_plot_{uuid.uuid4().hex}"
    try:
        spec = importlib.util.spec_from_file_location(module_name, script_path)
        if spec is None or spec.loader is None:
            raise RuntimeError("Could not create import spec for plotting script.")

        mod = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = mod
        spec.loader.exec_module(mod)

        # -------- Required wiring --------
        if hasattr(mod, "EXCEL_PATH"):
            mod.EXCEL_PATH = excel_path
        if hasattr(mod, "USE_FILE_DIALOG"):
            mod.USE_FILE_DIALOG = False
        if hasattr(mod, "UMAP_SHEET") and "sheet_coordinates" in pipeline_cfg:
            mod.UMAP_SHEET = pipeline_cfg["sheet_coordinates"]
        if hasattr(mod, "DATA_SHEET") and "sheet_binary" in pipeline_cfg:
            mod.DATA_SHEET = pipeline_cfg["sheet_binary"]
        if hasattr(mod, "INCLUDE_COLUMNS") and ("plot_include_columns" in pipeline_cfg):
            mod.INCLUDE_COLUMNS = pipeline_cfg.get("plot_include_columns")
        if hasattr(mod, "HIGHLIGHT_COLUMNS") and ("plot_highlight_columns" in pipeline_cfg):
            mod.HIGHLIGHT_COLUMNS = pipeline_cfg.get("plot_highlight_columns")
        elif hasattr(mod, "HIGHLIGHT_COLUMNS") and ("plot_include_columns" in pipeline_cfg):
            mod.HIGHLIGHT_COLUMNS = pipeline_cfg.get("plot_include_columns")
        if hasattr(mod, "GCHM_INCLUDE_COLUMNS") and ("gchm_include_columns" in pipeline_cfg):
            mod.GCHM_INCLUDE_COLUMNS = pipeline_cfg.get("gchm_include_columns")
        if hasattr(mod, "MAX_NROWS") and ("plot_max_nrows" in pipeline_cfg):
            try:
                mod.MAX_NROWS = max(1, int(pipeline_cfg.get("plot_max_nrows", getattr(mod, "MAX_NROWS", 2))))
            except Exception:
                pass

        # -------- Optional global overrides from GUI / caller --------
        overridable_globals = [
            "RUN_2D_DENSITY_PLOTS",
            "INCLUDE_GENOMIC_DENSITY_PANEL",
            "PLOT_MODE",
            "PANEL_W_IN",
            "PANEL_H_IN",
            "PNG_DPI",
            "FIGURE_FORMAT",
            "OUTPUT_EXT",
            "FIGURE_SUPTITLE",
            "SHOW_FIG",
            "SHOW_COLORBAR",
            "SUBPLOT_WSPACE",
            "SUBPLOT_HSPACE",
            "COLOR_MODE",
            "DENSITY_CMAP_NAME",
            "ENRICHMENT_CMAP_NAME",
            "ENRICHMENT_SCALE",
            "ENRICHMENT_STYLE",
            "ENRICHMENT_VMAX",
            "ENRICHMENT_PERCENTILE",
            "ENRICHMENT_SYMMETRIC",
            "DENSITY_OUTPUT_BASENAME",
            "RUN_GENE_CLUSTER_LOCALIZATION_HEATMAP",
            "GCHM_COLORMAP",
            "GCHM_SIGMA",
            "GCHM_SPREAD_FACTOR",
            "GCHM_HEIGHT_PER_CLUSTER",
            "GCHM_LABEL_FONTSIZE",
            "GCHM_DPI",
            "GCHM_CMAP_MIN_REL",
            "GCHM_CMAP_MAX_REL",
            "GCHM_OUTPUT_FILENAME",
            "GCHM_OUTPUT_BASENAME",
            "GCHM_SHOW_FIG",
            "GCHM_INCLUDE_COLUMNS",
            "CODON_USAGE_PLOT_MODE",
            "CODON_USAGE_WORKBOOK",
            "CODON_USAGE_PNG_DPI",
            "CODON_USAGE_PANEL_W_IN",
            "CODON_USAGE_PANEL_H_IN",
            "CODON_USAGE_SHOW_FIG",
            "CODON_USAGE_OUTPUT_FORMAT",
            "CODON_USAGE_OUTPUT_BASENAME",
            "DENSITY_XMIN",
            "DENSITY_XMAX",
            "DENSITY_YMIN",
            "DENSITY_YMAX",
            "GCHM_XMIN",
            "GCHM_XMAX",
            "GCHM_YMIN",
            "GCHM_YMAX",
            "CODON_USAGE_XMIN",
            "CODON_USAGE_XMAX",
            "CODON_USAGE_YMIN",
            "CODON_USAGE_YMAX",
        ]
        for key in overridable_globals:
            if key in pipeline_cfg and hasattr(mod, key):
                setattr(mod, key, pipeline_cfg[key])

        # -------- Debug prints --------
        print(f"[INFO]   EXCEL_PATH : {getattr(mod, 'EXCEL_PATH', '(missing)')}")
        print(f"[INFO]   UMAP_SHEET : {getattr(mod, 'UMAP_SHEET', '(missing)')}")
        print(f"[INFO]   DATA_SHEET : {getattr(mod, 'DATA_SHEET', '(missing)')}")
        print(f"[INFO]   MAX_NROWS  : {getattr(mod, 'MAX_NROWS', '(missing)')}")
        if hasattr(mod, "COLOR_MODE"):
            print(f"[INFO]   {script_name}.COLOR_MODE = {getattr(mod, 'COLOR_MODE')}")
        if hasattr(mod, "FIGURE_FORMAT"):
            print(f"[INFO]   {script_name}.FIGURE_FORMAT = {getattr(mod, 'FIGURE_FORMAT')}")

        if hasattr(mod, "main") and callable(mod.main):
            mod.main()
        else:
            raise RuntimeError(f"{script_name} has no callable main(); nothing executed.")

    except Exception as e:
        print(f"[WARN] Failed to run plotting pipeline script: {e}")
        print(traceback.format_exc())
