# nf-core/deepmutscan: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0dev - [unreleased]

### `Changed`

- Nextflow >= 26.10.0 is required (was 25.10.4)
- Raw data QC: `RAW_DATA_QC` checks the two FASTQ files of a library with `fq lint` and, for libraries that pass, writes the FastQC report of each file with [Falco](https://github.com/smithlabcode/falco) 2.0.3, in one task per library, replacing `FASTQ_CHECK` and nf-core's `FASTQC`. Falco computes FastQC's modules and writes FastQC's `fastqc_data.txt`, which MultiQC (and with it the all-in-one report) shows as a section named Falco with the same plots. On the test reads the numbers match FastQC's (sequence counts, duplication, overrepresented sequences, quality and N content, up to rounding) and the module grades agree except overrepresented sequences and adapter content, which Falco 2.0.3 grades against cutoffs in the wrong units (PASS where FastQC says WARN or FAIL; reported as [smithlabcode/falco#145](https://github.com/smithlabcode/falco/issues/145), fix proposed in [#146](https://github.com/smithlabcode/falco/pull/146); the values are correct). Per 1.2 M reads Falco takes 1.4 CPU s and 14 MB instead of FastQC's 14.5 CPU s and 0.6 GB; on the 50 000-pair test libraries the QC of a library takes 0.4-0.7 s instead of 6-9 s. The reports are published to `raw_data_qc/` (`*_fastqc_data.txt`, `*_falco_report.html`, `*_falco_summary.txt`) instead of `fastqc/`
- Aligner: bwa-mem3 (`--fast` preset) replaces BWA 0.7.19 in `MERGE_ALIGN`, and `BWAMEM3_INDEX` (nf-core module) replaces `BWA_INDEX`. For a gene-sized reference bwa-mem3 produces alignments identical to BWA with the same scoring in about a third of the CPU time; `--fast` halves that again and changes a single read in one of the four test libraries. Set `ext.args3 = '-L 100 -O 12,12 -B 3'` (no `--fast`) in a custom config for alignments byte-identical to BWA.
- Alignment scoring in `MERGE_ALIGN` (`ext.args3`, default `--fast -L 100 -O 12,12 -B 3`): substitutions are aligned end-to-end instead of locally. With BWA-MEM's defaults, 2-base codon changes within 3 bp and 3-base changes within 7 bp of a read end were soft-clipped away together with their coverage, under-counting full-codon substitutions by about 13% relative to single-base changes on 150 bp merged reads and making the recovery of read-linked double mutants depend on their codon distance. Count tables gain the previously clipped variants; the nf-test snapshot is updated accordingly.
- Variant counter: the Python template (`pysam`/`polars`) is replaced by `dms_count`, a single-pass C/htslib counter. Count tables are byte-identical to v1.0.0; `VARIANTCOUNTING` needs about 2% of the CPU time and 14 MB of memory, and its resource label drops from `process_high` to `process_single`.

### `Added`

- Carbon-footprint estimate of every run with the [nf-co2footprint](https://nextflow-io.github.io/nf-co2footprint/) plugin: `pipeline_info/co2footprint_{report,summary,trace}_<timestamp>.*`, with `--co2_location` for the compute site's grid zone

### `Dependencies`

| Dependency        | Version          |
| ----------------- | ---------------- |
| `nf-co2footprint` | 1.4.1 (new)      |
| `falco`           | 2.0.3 (new)      |
| `fastqc`          | 0.13.0 (removed) |
| `bwa-mem3`        | 0.17.0 (new)     |
| `bwa`             | 0.7.19 (removed) |
| `samtools`        | 1.25 (was 1.21)  |
| `dms-count`       | 0.1.0 (new)      |
| `htslib`          | 1.25 (new)       |
| `pysam`           | 0.24.0 (removed) |
| `polars`          | 1.33.1 (removed) |
| `pyarrow`         | 24.0.0 (removed) |
| `biopython`       | 1.87 (removed)   |

## v1.0.0 - Puyol - [2026-11-14]

Initial release of nf-core/deepmutscan, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Samplesheet input of paired-end FASTQ files per library, annotated with `sample`, `type` (`input`, `output`, `wildtype`) and `replicate`, validated at start-up (one sample per run, unique `sample`/`type`/`replicate`, matching `input` and `output` replicate numbers, one `wildtype` library for `--error_correction wildtype`), together with `--reading_frame` against the reference
- Input FASTQ checks before any reads are processed: quick checks of the first read of every file at start-up (empty, not gzip or not FASTQ, same file or same reads twice, mates with different names, file shared between libraries), then `fq lint` on the complete files (truncated or corrupted gzip, malformed records, unequal read numbers, mates out of order); all faulty libraries are reported together
- Raw read QC with FastQC and MultiQC
- Read-pair merging into consensus reads with `vsearch --fastq_mergepairs`
- Alignment of the merged reads to a gene-sized reference ORF with BWA-MEM (`--fasta`, `--reading_frame`)
- Filtering of unmapped, secondary, supplementary, low mapping-quality (MAPQ < 30) and indel-containing alignments with samtools, then coordinate sorting and indexing; wildtype reads are retained for error correction
- Light-weight variant counter built on the `pysam` and `polars` libraries in Python, replacing `GATK AnalyzeSaturationMutagenesis` with a column-compatible output, including a minimum base quality (`--base_qual`, default Q40) and read-edge exclusion window (`--min_flank`, default 2 bp)
- Annotation and filtering of variant counts against the programmed mutagenesis library (`--mutagenesis_type` `nnk`, `nns`, `nnh`, `nnn`, `nnk_nns`, `nnk_nns_nnh` or `custom` with `--custom_codon_library`)
- Single-nucleotide variant sequencing-error correction (`--error_correction`):
  - `false_doubles` (default): error rates estimated from read-linked false double mutant codons, by maximum likelihood (`--false_doubles_method mle`, default) or empirical Bayes (`eb`), within a configurable codon window (`--false_doubles_codon_window`)
  - `wildtype`: error profile from additional deep sequencing of the unmutated template (`type: wildtype` samplesheet rows)
  - `none`
- Interactive, run-level HTML error-correction report
- DMS library QC per library: count and count-per-coverage heatmaps, sorted count distributions with LogDiff (variants without counts drawn at a floor derived from the data and counted separately; LogDiff not reported when more than 10% of variants have no counts), sliding-window coverage and count profiles (`--sliding_window_size`, `--aimed_cov`), and sequencing-depth rarefaction on a geometric depth grid (`--run_seqdepth`)
- Optional fitness estimation (`--fitness`): merged count tables, experimental design, synonymous wildtype proxy selection, default log-ratio fitness with replicate rescaling and summary statistics, replicate correlation plots and fitness heatmap
- Optional fitness estimation with DiMSum (`--dimsum`) and mutscan edgeR / limma (`--mutscan`)
- Minimum number of reads a variant needs in every `input` replicate to receive a fitness estimate (`--min_counts`, default 10), applied in the same way by the default estimator, DiMSum and mutscan; variant counting, the sequencing-error correction and the library QC always use all reads
- Pseudocount for variants that drop out of an `output` replicate (`--output_pseudocount`, default 1; 0 gives dropouts no fitness value)
- Optional interactive 3D variant effect inspection tool built from a user-supplied wildtype structure (`--pdb`)
- Self-contained, all-in-one run report (`deepmutscan_report.html`) embedding QC, error-correction, fitness and MultiQC results, with links to Nextflow's execution report and timeline
- Docker and Singularity images (amd64 and arm64) and conda lock files for all modules, built from each module's `environment.yml`; DiMSum and mutscan are amd64 only for now
- `test` profile on a 50,000 read-pair subsample of a GID1A nicking-mutagenesis GluePCA experiment, with three nf-test pipeline tests (all options; without fitness, error correction, SeqDepth and structure; stub run), and `test_full` profile on the complete dataset (ENA PRJEB110196); test data and the PDB structure are hosted on [nf-core/test-datasets](https://github.com/nf-core/test-datasets/tree/deepmutscan)

### `Dependencies`

| Dependency    | Version               |
| ------------- | --------------------- |
| `fastqc`      | 0.13.0                |
| `fq`          | 0.12.0                |
| `multiqc`     | 1.35                  |
| `bwa`         | 0.7.19                |
| `samtools`    | 1.21, 1.22.1          |
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
