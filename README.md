# CodonPipe

CodonPipe measures synonymous codon usage and tRNA decoding across a bacterial
genome, groups genes by the codon-usage landscape they occupy, and tests whether
a gene set — a functional category, a genomic island, a set of virulence genes —
uses the genetic code differently from the rest of the genome.

It is a desktop application with a Tk interface, plus two scripts that run the
same pipeline without one.

## Installation

Python 3.10 or newer.

```
git clone <repository-url>
cd CodonPipe
pip install -r requirements.txt
python "Preloaded genomes/fetch_genomes.py"   # retrieves the reference genomes
```

## Running

```
python CodonPipe_GUI.py            # interface: load a genome, cluster, plot
python Clustering_Pipeline.py      # clustering and workbook export, no interface
python Plotting_Pipeline.py        # figures from an existing clustering workbook
```

A run takes a coding-sequence FASTA and, optionally, a cluster workbook listing
gene sets as columns of locus tags. It writes one Excel workbook holding the
codon table, the gene order, the cluster membership, the per-gene metrics and a
`Meta` sheet recording the parameters the run used, together with the figures.

## What it computes

- **Codon usage per gene**, as relative synonymous codon usage or as a
  percentage within selected amino-acid families. Group percentages are
  normalised by the synonymous families the selected codons span, never by all
  codons, so amino-acid composition cannot enter the measurement.
- **A codon-usage landscape**: a 2-D embedding (UMAP, t-SNE or PCA) of the
  standardised per-gene codon usage, clustered, and projected onto a genome axis.
- **Cluster statistics**: enrichment of each cluster along the genome axis
  against a permutation null, and comparisons of codon or tRNA-decoding values
  between clusters.
- **Starvation-sensitive and -insensitive codon usage**, under a choice of
  published definitions (`codonpipe/codon_sensitivity.py`). Percentages computed
  under different definitions are not comparable; the active definition is
  recorded in the `Meta` sheet of every workbook.
- **tRNA decoding**: wobble pairing, isoacceptor shift and modification-dependent
  codon groups, from a user-supplied decoding table.
- **Virulence-gene discovery** (`codonpipe/virulence/`): an optional add-on that
  assembles a virulence-gene set for a species from curated databases and
  published genome-wide screens, and reports it as four tiers of decreasing
  permissiveness.

## Statistics

Comparisons between two gene sets use the Mann-Whitney U test by default; the
distributions are heavy-tailed and the groups differ in size, so a rank test is
the appropriate default. Student's and Welch's t-tests remain selectable.
Genome-axis enrichment uses a permutation null with the add-one estimator,
`(#{null >= observed} + 1) / (permutations + 1)`, so a p-value is never reported
as zero.

## Reference genomes

`Preloaded genomes/genomes.tsv` lists each genome with its INSDC accessions and
CDS count; `fetch_genomes.py` downloads them. The files are not kept in the
repository because they total roughly 230 MB. NCBI re-annotates genomes over
time, so a freshly downloaded file may contain slightly more or fewer coding
sequences than the manifest records; the downloader reports the difference.
One genome, *Chrysiobacterium nematophagum*, is an in-house assembly with no
public accession and is shipped with the repository.

## Dependencies

Required: numpy, pandas, scipy, matplotlib, scikit-learn, openpyxl.

Optional, imported defensively so the rest of the application runs without them:

| package | enables |
| --- | --- |
| umap-learn | UMAP embedding (t-SNE and PCA are built in) |
| scikit-learn-extra | k-medoids clustering |
| astroML | the two-dimensional Kolmogorov-Smirnov test |
| suds-community, certifi | DAVID functional-annotation queries |
| requests | the virulence add-on's database lookups |

`scikit-learn-extra` is no longer maintained and does not build against recent
scikit-learn releases. Leave it out unless you need k-medoids; every other
clustering method works without it.

## Repository layout

```
CodonPipe_GUI.py          interface
Clustering_Pipeline.py    clustering and workbook export
Plotting_Pipeline.py      figures from an existing workbook
codonpipe/                library: codon metrics, clustering, plots, exports
codonpipe/virulence/      virulence-gene discovery add-on
Preloaded genomes/        genome manifest, downloader, cluster workbooks
```

## Citation

See `CITATION.cff`.

## Licence

MIT. See `LICENSE`.
