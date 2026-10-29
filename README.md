# Next-Generation Sequencing — Genomics & DNA Analysis

**Elena Dalboni da Rocha · [@enterlina](https://github.com/enterlina)**  
[Design & Product Portfolio](https://buildnxt.app) · [Machine Learning](https://github.com/enterlina/MachineLearning) · [AI Engineering](https://github.com/enterlina/ai-engineering)

A computational genomics portfolio documenting **9 Python and Jupyter implementations** from my bioinformatics studies. This work focuses on analyzing DNA sequencing reads, understanding errors introduced during sequencing and alignment, evaluating correction methods, and constructing graphs for genome assembly.

**Core areas:** FASTQ quality control · DNA read alignment · genomic coverage · insertion/deletion analysis · long-read sequencing · error correction · de Bruijn graphs

## The workflow

`FASTQ reads` → **base quality & GC analysis** → `SAM/BAM alignment` → **coverage & error profiles** → **correction assessment** → **de Bruijn graph assembly experiments**

## Algorithms & experiments

| Topic | What the code explores | Source |
| --- | --- | --- |
| 01 · FASTQ quality | FASTQ read QC: GC-content distributions under quality filters and per-base Phred quality boxplots. | [`01-fastq-gc-content-and-phred-quality.py`](algorithms/01-fastq-gc-content-and-phred-quality.py) |
| 01 · FASTQ quality | Biopython FASTQ parsing, quality-aware GC histograms, per-position quality and FastQC comparison. | [`01-fastq-quality-control-notebook.ipynb`](algorithms/01-fastq-quality-control-notebook.ipynb) |
| 02 · Alignment analysis | SAM/BAM alignment depth via difference-array prefix sums; insert-size distribution and reference substitution counts. | [`02-alignment-coverage-insert-size-substitutions.py`](algorithms/02-alignment-coverage-insert-size-substitutions.py) |
| 02 · Error profiling | Aligned-read substitutions and indels, coverage, error percentages, homopolymer length discrepancies and base-quality comparisons. | [`02-sequencing-error-profiles-indels-homopolymers.py`](algorithms/02-sequencing-error-profiles-indels-homopolymers.py) |
| 03 · Long-read sequencing | Comparison of PacBio and Oxford Nanopore alignments using genome coverage, substitution tables and indel length profiles. | [`03-long-read-platform-coverage-error-comparison.py`](algorithms/03-long-read-platform-coverage-error-comparison.py) |
| 04 · Error correction | Base-level comparison of raw versus corrected alignments, partitioning true corrections, missed errors, false corrections and removed bases. | [`04-read-error-correction-evaluation.py`](algorithms/04-read-error-correction-evaluation.py) |
| 05 · Genome assembly | k-mer / reverse-complement De Bruijn graph construction, coverage-aware edge compaction, tip simplification and graph exports. | [`05-debruijn-graph-genome-assembly-notebook.ipynb`](algorithms/05-debruijn-graph-genome-assembly-notebook.ipynb) |
| 05 · Genome assembly | Graph vertex/edge construction from sequence reads, reverse complements, compaction and DOT/PNG visualization. | [`05-debruijn-graph-construction-visualization.py`](algorithms/05-debruijn-graph-construction-visualization.py) |
| 05 · Genome assembly | Extended genome assembly graph operations, path compaction, low-coverage/short-tip removal and sequence export. | [`05-debruijn-graph-tip-removal-assembly.py`](algorithms/05-debruijn-graph-tip-removal-assembly.py) |

## Engineering methods

- **Sequence data processing:** Parsing FASTQ/FASTA formats, Phred scores, nucleotides and reverse complements with Biopython and pysam.
- **Alignment analytics:** Per-position coverage using a difference array and prefix sums; insert-size histograms; substitution and indel matrices.
- **Sequencing technology comparison:** Code for evaluating read-alignment error distributions from PacBio and Oxford Nanopore datasets.
- **Error-correction audit:** Classifying bases as corrected, uncorrected, falsely modified or removed by comparing aligned reads.
- **Graph algorithms for genome assembly:** Building a k-mer De Bruijn graph, maintaining directed edges and read coverage, merging non-branching paths, pruning tips and exporting graph visualizations.

## Tools

Python · Biopython · pysam · NumPy · pandas · Matplotlib · seaborn · Jupyter · PyGraphviz

## Reproducibility

These are **historical algorithm exercises, not a clinical diagnostic pipeline or production-validated genome assembler**. The original FASTQ/FASTA/SAM/BAM datasets and evaluation outputs are not included in the repository. Several scripts use absolute local paths and legacy dependencies; not every implementation has been tested for current Python versions. Code is preserved as authored; the descriptive filenames and this index do not imply that all analyses run end-to-end or have validated accuracy.

The repository also retains historical IDE and virtual-environment artifacts for provenance; they are not recommended as a reproducible environment.

## Author

Elena Dalboni da Rocha — AI Product Designer, Design Engineer and bioinformatics graduate, combining computational biology foundations with modern AI product development.

[Portfolio](https://buildnxt.app) · [GitHub](https://github.com/enterlina)
