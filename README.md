# Next-Generation Sequencing Analysis

Python and Jupyter coursework in sequencing-data analysis, completed during my bioinformatics studies. The exercises examine read quality, alignment errors, genome coverage, and graph-based assembly.

## What is included

| Area | Code |
| --- | --- |
| FASTQ quality scores, GC-content distributions, per-position quality plots | [HW1](HW1/) |
| Alignment statistics, substitutions, insertions/deletions, and quality-related errors | [NGS_HW2_2.py](NGS_HW2_2.py), [HW2](HW2/) |
| Coverage plots and substitution/indel summaries | [NGS_HW3.py](NGS_HW3.py) |
| Comparing errors before and after read correction | [NGS_HW_4.py](NGS_HW_4.py) |
| De Bruijn graph construction and visualization | [HW5_1](HW5_1/) |

## Tools

Python, Biopython, pysam, NumPy, pandas, Matplotlib, seaborn, Jupyter, and PyGraphviz for graph exercises.

## Reproducibility

The scripts reference external FASTQ, FASTA, and alignment files through local paths. Those input datasets are not included here; paths and dependencies need to be adapted before running. Historical virtual-environment files are present, but there is no maintained reproducible environment. This repository documents educational analyses and algorithm experiments.

## Author

Elena Dalboni da Rocha · [Portfolio and current work](https://buildnxt.app) · [GitHub](https://github.com/enterlina)
