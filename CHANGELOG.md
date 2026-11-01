# nf-core/deepmutscan: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0dev - [unreleased]

### `Changed`

- Nextflow >= 26.10.0 is required (was 25.10.4)
- Raw data QC: `RAW_DATA_QC` checks the two FASTQ files of a library with `fq lint` and, for libraries that pass, writes the FastQC report of each file with [Falco](https://github.com/smithlabcode/falco) 2.0.3, in one task per library, replacing `FASTQ_CHECK` and nf-core's `FASTQC`. Falco computes FastQC's modules and writes FastQC's `fastqc_data.txt`, which MultiQC (and with it the all-in-one report) shows as a section named Falco with the same plots. On the test reads the numbers match FastQC's (sequence counts, duplication, overrepresented sequences, quality and N content, up to rounding) and the module grades agree except overrepresented sequences and adapter content, which Falco 2.0.3 grades against cutoffs in the wrong units (PASS where FastQC says WARN or FAIL; reported as [smithlabcode/falco#145](https://github.com/smithlabcode/falco/issues/145), fix proposed in [#146](https://github.com/smithlabcode/falco/pull/146); the values are correct). Per 1.2 M reads Falco takes 1.4 CPU s and 14 MB instead of FastQC's 14.5 CPU s and 0.6 GB; on the 50 000-pair test libraries the QC of a library takes 0.4-0.7 s instead of 6-9 s. The reports are published to `raw_data_qc/` (`*_fastqc_data.txt`, `*_falco_report.html`, `*_falco_summary.txt`) instead of `fastqc/`
- Aligner: bwa-mem3 (`--fast` preset) replaces BWA 0.7.19 in `MERGE_ALIGN`, and its index is built by `PREPARE_REFERENCE` instead of `BWA_INDEX`. For a gene-sized reference bwa-mem3 produces alignments identical to BWA with the same scoring in about a third of the CPU time; `--fast` halves that again and changes a single read in one of the four test libraries. Set `ext.args3 = '-L 100 -O 12,12 -B 3'` (no `--fast`) in a custom config for alignments byte-identical to BWA.
- Alignment scoring in `MERGE_ALIGN` (`ext.args3`, default `--fast -L 100 -O 12,12 -B 3`): substitutions are aligned end-to-end instead of locally. With BWA-MEM's defaults, 2-base codon changes within 3 bp and 3-base changes within 7 bp of a read end were soft-clipped away together with their coverage, under-counting full-codon substitutions by about 13% relative to single-base changes on 150 bp merged reads and making the recovery of read-linked double mutants depend on their codon distance. Count tables gain the previously clipped variants; the nf-test snapshot is updated accordingly.
- Variant counter: the Python template (`pysam`/`polars`) is replaced by `dms_count`, a single-pass C/htslib counter. Count tables are byte-identical to v1.0.0; `VARIANTCOUNTING` needs about 2% of the CPU time and 14 MB of memory, and its resource label drops from `process_high` to `process_single`.
- `MERGE_ALIGN_COUNT`: merging, alignment, filtering and variant counting run as one stream per library (`vsearch | bwa-mem3 | samtools view -u | dms_count --bam -`), replacing `MERGE_ALIGN` and `VARIANTCOUNTING`. `dms_count` reads the alignments as they arrive and never needed them sorted or indexed, so the `samtools sort`, the BAM (~2 GB per full-size library) and one task per library are gone; reads with an insertion, deletion or reference skip are skipped by `dms_count --skip_indel_reads` instead of an `awk` filter. Read merging uses vsearch 2.32.0 (from 2.30.0): its directly addressed k-mer index halves the CPU time of the merge (10.2 -> 5.0 CPU s per 1.2 M pairs) with the same merged reads, discarded pairs and statistics, and `--no_progress` keeps the merging log as it was. Count tables and merging logs are byte-identical on the four test libraries and a 1.2 M read-pair benchmark library; on two CPUs that library takes 6.9 s and 0.43 GB instead of 9.7 s and 0.84 GB for the two separate tasks
- `--save_bam` (default off): alignment files are written only on request. The filtered alignments are then also sorted (compression level 1) and indexed by a background `samtools sort` fed through `tee`, and published to `intermediate_files/bam_files/sorted/`; they now include the reads with insertions or deletions, which the counter skips, so the BAM is about 25% larger than before
- False-doubles error correction rebuilt as three processes - evidence per library (`DMSANALYSIS_FDM_EVIDENCE`), estimation per pool (`DMSANALYSIS_FDM_ESTIMATE`), correction per library - with evidence pooled across the libraries of a sample (`--false_doubles_pooling all|replicates|none`), empirical Bayes as the default estimator (`--false_doubles_method eb`) with a trinucleotide-context prior (`--false_doubles_prior`), and the error-correction report setting per-library estimates against the pooled rates. Corrected counts change relative to v1.0.0, in a depth-dependent way for shallow libraries
- `--false_doubles_codon_window` renamed to `--max_FDdist` (default 40), the companion of `--min_FDdist`: anchors are used when `min_FDdist` < distance <= `max_FDdist`, and the upper bound is now inclusive. Swept from 15 to 120 codons against matched wildtype-template sequencing on five libraries, 40 is the best setting; beyond it the evidence carries an additive background of order 1e-5 that the genuine-double subtraction does not remove. Passing the old name is an error with a message pointing at the new one
- `--min_FDdist` (default 13): anchors closer than this many codons to a variant are excluded from the false-doubles evidence, where the primer's own recurrent errors dominate; replaces the earlier depth-dependent filter on the highest false-double counts
- Default error-signature candidates curated from 23 to 14: dropped SBS54 (germline-variant contamination), SBS30 (endogenous base-excision-repair deficiency, a shape proxy for the deamination signature already in the set), SBS95 (single study, no replication, 0.90 to SBS18) and the six single-channel shapes SBS48, SBS49, SBS52, SBS56, SBS59, SBS60. The full set is bundled as `candidate_signatures_GRCh37_extended.csv`
- False-doubles error correction from the fdmcorrect R package (`bioconda::dms-fdmcorrect` 0.1.2, [github.com/MaximilianStammnitz/dms-fdmcorrect](https://github.com/MaximilianStammnitz/dms-fdmcorrect)): `DMSANALYSIS_FDM_EVIDENCE`, `DMSANALYSIS_FDM_ESTIMATE` and `DMSANALYSIS_ERROR_CORRECTION_FALSE_DOUBLES` call the package, which carries the model and also runs on its own, outside the pipeline. Evidence extraction takes seconds per library instead of minutes
- Error-correction report: the spectrum section now prints what each selected signature stands for, from the bundled `signature_descriptions.csv`
- Error-correction report: the plotted false-double spectrum is now built on the genuine-double-adjusted reads (`F_adj`), the same evidence the correction applies, rather than on raw `F` - the panel, its tile and the rates used downstream are one quantity. The six substitution classes follow the COSMIC/SigProfiler SBS-96 colour convention throughout the report (spectrum, selection panels, ORF panel, correction-magnitude histogram and every legend), with T>A and T>G darkened for legibility on white; the spectrum carries the COSMIC class band, whole-percent axis ticks and the per-channel error rate on hover, and the signature reconstruction is drawn in grey tones. New per-spectrum fields `rate_per_mb`, `burden_per_mb` and the per-channel genuine-double range in `error_signatures.json`
- Error-correction report: sections named and ordered ERROR REPRODUCIBILITY, ERROR EVIDENCE, ERROR SPECTRUM, ERROR IMPACT, SINGLE FILE; every error-frequency axis fixed at 10^-7 to 10^-2; the spectrum section shows only the error model the run corrects with, and the anchors section only under `--error_correction false_doubles`; raw-versus-corrected is now one panel per sequencing file; the ORF panel shows the error rate by substitution class; the genuine- and false-double distance profiles share one rate axis, with the 2+2/3+2/3+3 versus 2+1/3+1 distinction spelled out; new `gd_anchor_exposure` column in `fdm_distance_profile.csv`
- `--subtract_genuine_doubles` (default true): genuine double integrations are subtracted from the false-double evidence before the error rate is estimated. kappa(d) is fitted per pool on reads carrying two programmed multi-nucleotide variants and applied per anchor distance; _F_ is split between error and genuine doubles by expectation-maximisation. New per-pool output `fdm_kappa_profile.csv` and new `fdm_rates.csv` columns `kappa_exposure`, `G`, `F_error`
- Mutational-signature fit of the error spectrum (false-doubles correction; `--skip_error_signatures`): the false-double, wildtype-template and excess spectra of every pool are fitted to COSMIC artefact and in-vitro DNA-damage signatures (`--error_signatures_candidates`) by two-stage non-negative least squares with a site bootstrap (`--error_signatures_bootstraps`, `--error_signatures_max_cosine`), after removing genuine double mutants measured on reads with two programmed multi-nucleotide variants (`--error_signatures_genuine_doubles`); shown in the error-correction report and written to `error_signatures/`
- Library QC: the six per-library plot processes (`VISUALIZATION_COUNTS_HEATMAP`, `VISUALIZATION_COUNTS_PER_COV`, `VISUALIZATION_GLOBAL_POS_BIASES_COUNTS`, `VISUALIZATION_GLOBAL_POS_BIASES_COV`, `VISUALIZATION_LOGDIFF`, `VISUALIZATION_SEQDEPTH`) are one process, `VISUALIZATION_LIBRARY_QC`, which runs the same R code in one task per library; the same files are published to `library_QC/<library>/`. Each plot takes seconds, so scheduling six tasks per library for them cost more than drawing them
- Reference tables: `DMSANALYSIS_REFERENCE_TABLES` writes the wildtype amino-acid sequence and the table of programmed codon variants in one task, replacing `DMSANALYSIS_AASEQ` and `DMSANALYSIS_POSSIBLE_MUTATIONS` (same R code, same files in `intermediate_files/`)
- Fitness: `FITNESS_ESTIMATION` runs the default fitness estimation of a sample in one task - merging the per-library counts, deriving the experimental design, picking the synonymous wildtype proxy, estimating fitness and drawing the correlation plots and heatmaps - replacing the six-process chain `MERGE_COUNTS` → `EXPDESIGN_FITNESS` → `FIND_SYNONYMOUS_MUTATION` → `FITNESS_CALCULATION` → `FITNESS_QC` → `FITNESS_HEATMAP` (same R code, same files under `fitness/`). DiMSum and mutscan stay separate processes and take the bundle's tables
- `--mutscan`: the step no longer uses the mutscan package. It calls edgeR and limma directly, exactly as mutscan's `calculateRelativeFC()` does (wildtype counts as per-sample offsets, quasi-likelihood F-test, voom), builds the per-mutant count matrix the way `summarizeExperiment()` aggregates it, and draws the same pairs and volcano plots; the two result tables are identical to mutscan's, except that a TTG or CTG variant codon at the first position of the ORF is now named L, as the genetic code reads it (mutscan's names inherit Biostrings' initiation-codon rule and call it M). The result tables now carry the mutant name as their first column (`mutant`), which the previous template did not write. The variant naming is vectorised (15 min → 10 s on a full-depth sample), the environment shrinks to edgeR, limma, ggplot2 and GGally, and the step runs natively on arm64
- `--dimsum`: DiMSum's STEAM stages (4 and 5, the only ones the pipeline runs) are run from `bin/dimsum_fitness.R`, DiMSum 1.4's own R code for those stages (MIT, credited) with a driver in place of `dimsum()`, instead of through the `r-dimsum` Bioconda package, whose recipe pins R 4.0 and depends on the WRAP-stage tools (fastqc, cutadapt, vsearch, starcode) and Bioconductor packages the pipeline never calls. Biostrings and Cairo are replaced by base R, so the environment is data.table, ggplot2, GGally, hexbin, reshape2, plyr, seqinr, cowplot, gridExtra, optparse, rmarkdown and pandoc - all with arm64 builds. Outputs are identical (every table, the MaveDB export and the `.RData` objects, with and without the error model); the per-variant Hamming distances and the count-histogram jitter are vectorised and the preliminary report render is skipped, so the step takes about two thirds of the time and half the memory
- `--dimsum`: figures that DiMSum writes as both PNG and PDF are laid out once and drawn into both files, R's just-in-time compiler is off for the script, and with the error model on one CPU (the default) the 100 bootstraps run in the main R process - from the random-number stream DiMSum's one-worker cluster used - with the per-row variances and the model columns computed without per-call overhead. Every table, `.RData` object and PNG is identical (the PDFs show the same figures, laid out with the PNG device's font metrics); a full-depth sample takes 28 s instead of 39 s, and 51 s instead of 79 s with `--fitnessErrorModel T`
- `--dimsum`: DiMSum's error model is fitted by default (`--dimsum_error_model`, default `true`; previously `--fitnessErrorModel F` in `ext.args`): replicate-specific error terms and replicate normalisation, as in DiMSum's own default. The `test` profile switches it off, because its subsampled libraries have too few variants covered in every replicate; a library that is too sparse makes `RUN_DIMSUM` stop with DiMSum's message and a pointer to `--dimsum_error_model false`

### `Added`

- Carbon-footprint estimate of every run with the [nf-co2footprint](https://nextflow-io.github.io/nf-co2footprint/) plugin: `pipeline_info/co2footprint_{report,summary,trace}_<timestamp>.*`, with `--co2_location` for the compute site's grid zone
- `--dimsum_error_model` (default `true`): whether DiMSum fits its error model (`--fitnessErrorModel`)

### `Fixed`

- Translation of the first codon of the mutagenised window: `Biostrings::translate()` reads TTG and CTG at position 1 as methionine (its initiation-codon rule). The default fitness estimation and the reference tables now translate every codon by the genetic code (`no.init.codon = TRUE`), so a TTG/CTG variant at the first position is reported as L; `--mutscan` does the same (see Changed)
- False-doubles correction: single-nucleotide variants without an observed false double (_F_ = 0) are no longer dropped from the error-rate fit and left uncorrected when false doubles were seen for them that did not pair with a programmed variant; this made the fitted error rate depend on sequencing depth
- False-doubles correction: the maximum possible joint coverage of a false double is the numeric, not the alphabetical, minimum of the two codons' coverages
- False-doubles correction: wildtype-template libraries are no longer corrected with rates estimated from their handful of chance false doubles
- False-doubles correction: the 99.5 % count filter on false doubles is removed. It discarded reads concentrated inside the mutagenesis primer footprint (79-99 % of them, at a median separation of 5-7 codons) rather than the genuine double mutants it was meant to remove, and the amount it removed depended on sequencing depth; `--min_FDdist` targets the same artefact directly
- False-doubles correction: a library without double mutants or without usable anchors, or with an empty count table, no longer stops the evidence, estimation or correction step; its pool is estimated from the other libraries, or left uncorrected with a message in the log
- False-doubles correction: the libraries of a pool are combined in the order of their names rather than the order in which their tasks finished, so the pooled rates and priors no longer change in the last digits from one run to the next
- False-doubles correction: `--false_doubles_method mle` needs at least 1,000 false-double reads in a pool (`eb`: 100); measured on GID1A and OsTIR1 pools thinned to 0.02-30 % of their reads, the per-variant _F_/_M_ is worse than no correction below that, while `eb` beats no correction at every depth
- Library QC: a plot that cannot be drawn for a very sparse library (for example the sliding-window or sequencing-depth curves of a library with only a handful of covered positions) is replaced by a one-page PDF giving the reason, with a warning in the log, instead of stopping the run
- A reference FASTA with lower-case (soft-masked) bases no longer stops the false-doubles correction with a reference/counts coordinate mismatch, nor gives its positions substitution classes outside the 6 types and 96 contexts: the reference is converted to upper case before any step reads it. `PREPARE_REFERENCE` does this with GNU Awk and builds the bwa-mem3 index of the upper-case reference, in one task per run; the index files are named after the reference (`intermediate_files/bwamem3/<name>.uppercase.fasta.*`)

### `Dependencies`

| Dependency        | Version             |
| ----------------- | ------------------- |
| `nf-co2footprint` | 1.4.1 (new)         |
| `falco`           | 2.0.3 (new)         |
| `fastqc`          | 0.13.0 (removed)    |
| `bwa-mem3`        | 0.17.0 (new)        |
| `bwa`             | 0.7.19 (removed)    |
| `samtools`        | 1.25 (was 1.21)     |
| `dms-count`       | 0.1.0 (new)         |
| `dms-fdmcorrect`  | 0.1.2 (new)         |
| `htslib`          | 1.25 (new)          |
| `pysam`           | 0.24.0 (removed)    |
| `polars`          | 1.33.1 (removed)    |
| `pyarrow`         | 24.0.0 (removed)    |
| `biopython`       | 1.87 (removed)      |
| `vsearch`         | 2.32.0 (was 2.30.0) |
| `mutscan`         | 1.0.0 (removed)     |
| `edgeR`           | 4.8.2 (direct)      |
| `limma`           | 3.66.0 (direct)     |
| `r-dimsum`        | 1.4 (removed)       |
| `gawk`            | 5.3.1 (new)         |

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
