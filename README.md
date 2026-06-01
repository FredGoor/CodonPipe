[README_CodonPipe.md](https://github.com/user-attachments/files/28458020/README_CodonPipe.md)
# CodonPipe

**CodonPipe** is a Python GUI pipeline for genome-wide analysis of bacterial codon usage, gene clustering, and tRNA decoding strategies. It was designed to explore how codon usage varies across bacterial genomes and how these patterns relate to functional gene groups, wobble decoding, and tRNA modification-dependent decoding.

## Main features

- Import coding-sequence FASTA files and compute per-gene codon usage.
- Analyze codon usage using absolute codon usage (ACU), relative codon usage (RCU), or amino-acid usage.
- Perform dimensionality reduction using UMAP, t-SNE, or PCA.
- Cluster genes using k-means, k-medoids, hierarchical clustering, DBSCAN, or spectral clustering.
- Generate publication-ready figures, including:
  - genome-wide codon-usage heatmaps,
  - 2D embedded gene-density plots,
  - cluster-specific codon-usage plots,
  - gene-cluster localization heatmaps,
  - tRNA decoding and wobble-decoding enrichment heatmaps.
- Export structured Excel workbooks containing reordered genes, cluster assignments, codon-usage tables, and tRNA-decoding summaries.
- Optional DAVID functional enrichment integration for sliding-window analyses.

## Repository structure

```text
CodonPipe_GUI.py          Main graphical interface
Clustering_Pipeline.py    Main analysis pipeline
Plotting_Pipeline.py      Standalone plotting pipeline
codonpipe/                Core modules used by the pipeline
```

Core modules include codon counting, clustering, Excel output generation, density plotting, 2D KS testing, gene-cluster heatmaps, tRNA decoding analyses, and DAVID enrichment utilities.

## Installation

Create and activate a conda environment:

```bash
conda create -n codonpipe python=3.10
conda activate codonpipe
```

Install the required packages:

```bash
pip install numpy pandas scipy matplotlib scikit-learn umap-learn openpyxl xlsxwriter
```

Optional packages:

```bash
pip install scikit-learn-extra suds-py3
```

`scikit-learn-extra` enables k-medoids clustering. `suds-py3` is only required for DAVID web-service analyses.

## Running CodonPipe

From the repository folder, launch the graphical interface:

```bash
python CodonPipe_GUI.py
```

The GUI allows the user to select input files, analysis parameters, clustering options, tRNA-decoding settings, and figure-export options.

## Input files

The minimal input is a coding-sequence FASTA file containing bacterial CDS sequences.

Optional inputs include:

- a user-defined gene-cluster file in Excel, CSV, TSV, or TXT format;
- a tRNA decoding table defining codon-to-tRNA decoding relationships;
- a DAVID gene-to-term file for functional enrichment analyses;
- a tRNA abundance table for abundance-correlation analyses.

## Main outputs

CodonPipe writes a structured output folder containing:

- `Clustering analysis results.xlsx` — main workbook with reordered genes, coordinates, binary cluster membership, and metadata;
- `Gene lists per cluster.xlsx` — gene lists and annotations for each cluster;
- `Codon usage tables per cluster.xlsx` — per-cluster codon-usage summaries;
- `Gene IDs.xlsx` — parsed gene identifiers and annotations;
- `Figures/` — exported heatmaps, density plots, and cluster-level figures;
- `Methods/` — automatically generated methods and figure-legend text files;
- optional DAVID enrichment outputs.

## Typical workflow

1. Launch `CodonPipe_GUI.py`.
2. Select a CDS FASTA file.
3. Choose the codon-usage metric: ACU, RCU, or amino-acid usage.
4. Select a dimensionality-reduction method and clustering method.
5. Optionally provide gene clusters or tRNA decoding tables.
6. Run the analysis.
7. Inspect the exported Excel files and figures.

## Notes

- CodonPipe was developed for bacterial CDS datasets.
- Relative codon usage is calculated within synonymous codon families.
- Missing relative-usage values are preserved in exported tables but are imputed with genome-average values only for dimensionality reduction and clustering.
- tRNA-decoding analyses depend on the supplied decoding table and should be interpreted according to the assumptions encoded in that table.

## Citation

If you use CodonPipe in published work, please cite the associated manuscript or repository release.

## License

Please specify the license for this repository before public release.
