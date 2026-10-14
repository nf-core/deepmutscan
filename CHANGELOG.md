# nf-core/deepmutscan: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0 - [unreleased<!-- TODO nf-core: replace with date on release -->]

Initial release of nf-core/deepmutscan, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Samplesheet input of paired-end FASTQ files per library, annotated with `sample`, `type` (`input`, `output`, `wildtype`) and `replicate`
- Raw read QC with FastQC and MultiQC
- Read alignment to a gene-sized reference ORF with BWA-MEM (`--fasta`, `--reading_frame`)
- Filtering of unmapped, secondary, low mapping-quality (MAPQ < 30) and indel-containing alignments with samtools; wildtype reads are retained for error correction
- Read-pair merging with `vsearch --fastq_mergepairs`, re-alignment, coordinate sorting and indexing
- Light-weight variant counter built on the `pysam` and `polars` libraries in Python, replacing `GATK AnalyzeSaturationMutagenesis` with a column-compatible output, including a minimum base quality (`--base_qual`, default Q40) and read-edge exclusion window (`--min_flank`, default 2 bp)
- Annotation and filtering of variant counts against the programmed mutagenesis library (`--mutagenesis_type` `nnk`, `nns`, `nnh`, `nnn`, `nnk_nns`, `nnk_nns_nnh` or `custom` with `--custom_codon_library`)
- Single-nucleotide variant sequencing-error correction (`--error_correction`):
  - `false_doubles` (default): error rates estimated from read-linked false double mutant codons, by maximum likelihood (`--false_doubles_method mle`, default) or empirical Bayes (`eb`), within a configurable codon window (`--false_doubles_codon_window`)
  - `wildtype`: error profile from additional deep sequencing of the unmutated template (`type: wildtype` samplesheet rows)
  - `none`
- Interactive, run-level HTML error-correction report
- DMS library QC per library: count and count-per-coverage heatmaps, sorted count distributions, sliding-window coverage and count profiles (`--sliding_window_size`, `--aimed_cov`), and sequencing-depth rarefaction on a geometric depth grid (`--run_seqdepth`)
- Optional fitness estimation (`--fitness`): merged count tables, experimental design, synonymous wildtype proxy selection, default log-ratio fitness with replicate rescaling and summary statistics, replicate correlation plots and fitness heatmap
- Optional fitness estimation with DiMSum (`--dimsum`) and mutscan edgeR / limma (`--mutscan`)
- Minimum input read count per replicate for fitness estimation (`--min_counts`, default 10), applied to the default estimator, DiMSum and mutscan
- Optional output pseudocount for dropout variants (`--output_pseudocount`, default 1; 0 disables), applied to the default estimator and DiMSum
- Optional interactive 3D variant effect inspection tool built from a user-supplied wildtype structure (`--pdb`)
- Self-contained, all-in-one run report (`deepmutscan_report.html`) embedding QC, error-correction, fitness, MultiQC and run statistics
- `test` profile on a 50,000 read-pair subsample of a GID1A nicking-mutagenesis GluePCA experiment, with nf-test snapshot testing, and `test_full` profile on the complete dataset (ENA PRJEB110196); test data and the PDB structure are hosted on [nf-core/test-datasets](https://github.com/nf-core/test-datasets/tree/deepmutscan)

### `Fixed`

### `Dependencies`

| Dependency    | Version               |
| ------------- | --------------------- |
| `fastqc`      | 0.12.1                |
| `multiqc`     | 1.35                  |
| `bwa`         | 0.7.19                |
| `samtools`    | 1.21, 1.22.1, 1.23.1  |
| `vsearch`     | 2.30.0                |
| `python`      | 3.12                  |
| `pysam`       | 0.24.0                |
| `polars`      | 1.33.1 (lts-cpu)      |
| `pyarrow`     | 24.0.0                |
| `biopython`   | 1.87                  |
| `numpy`       | 1.26.4, 2.5.1         |
| `pandas`      | 2.2.1                 |
| `r-base`      | 4.4.2 (DiMSum), 4.5.1 |
| `biostrings`  | 2.74.0, 2.78.0        |
| `r-tidyverse` | 2.0.0                 |
| `r-ggplot2`   | 3.5.1, 4.0.2          |
| `r-zoo`       | 1.8_15                |
| `r-dimsum`    | 1.4                   |
| `mutscan`     | 1.0.0                 |

### `Deprecated`
