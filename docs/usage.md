# nf-core/deepmutscan: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/deepmutscan/usage](https://nf-co.re/deepmutscan/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

**nf-core/deepmutscan** is a workflow for the analysis of deep mutational scanning (DMS) data. It takes short-read next-generation DNA sequencing data of a site-saturation mutagenesis library – typically before (`input`) and after (`output`) a variant selection experiment – and turns them into quality-controlled, error-corrected variant frequencies and, optionally, relative variant fitness estimates.

The pipeline was mainly designed for shotgun sequencing of long open reading frames (ORFs), where the mutated gene is randomly fragmented before sequencing and where resulting variant count signals are sparse, but it processes amplicon (tile) sequencing data in the same way.

This page explains how to prepare the inputs, how to run the pipeline, and what each processing step does and why. It is intended both for new DMS users and for developers who want to understand the exact design choices.

## Samplesheet input

You will need to create a samplesheet with information about the sequencing libraries you would like to analyse before running the pipeline. Use the `--input` parameter to specify its location. It has to be a comma-separated file with 5 columns and a header row, as shown in the example below.

```bash
--input '[path to samplesheet file]'
```

```csv title="samplesheet.csv"
sample,type,replicate,file1,file2
ORF1,input,1,/reads/input1_R1.fastq.gz,/reads/input1_R2.fastq.gz
ORF1,input,2,/reads/input2_R1.fastq.gz,/reads/input2_R2.fastq.gz
ORF1,output,1,/reads/output1_R1.fastq.gz,/reads/output1_R2.fastq.gz
ORF1,output,2,/reads/output2_R1.fastq.gz,/reads/output2_R2.fastq.gz
```

| Column      | Description                                                                                                                                                                                                              |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `sample`    | Name of the biological sample or experiment (e.g. the gene or condition). Must not contain spaces or `/`. As of v1.0.0, all rows of a samplesheet must have the same `sample`; analyse several samples in separate runs. |
| `type`      | Role of the library: `input` (before selection), `output` (after selection) or `wildtype` (deep sequencing of the unmutated template, used by `--error_correction wildtype`).                                            |
| `replicate` | Replicate number (positive integer). For fitness estimation, `input` and `output` libraries are paired by replicate number, so both need the same replicate numbers (e.g. 1, 2, 3).                                      |
| `file1`     | Full path to the FASTQ file for read 1, with extension `.fastq.gz`, `.fq.gz`, `.fastq` or `.fq`.                                                                                                                         |
| `file2`     | Full path to the FASTQ file for read 2, with extension `.fastq.gz`, `.fq.gz`, `.fastq` or `.fq`. As of v1.0.0, paired-end reads are required.                                                                            |

Each combination of `sample`, `type` and `replicate` may appear only once. Each library gets an internal identifier `<sample>_<type>_<replicate>_pe`, which is used to name its output files and folders.

The samplesheet is checked when the pipeline starts, and the run stops with an error message naming the problem if, for example, a file is missing, a row is duplicated or the replicate numbers of `input` and `output` libraries do not match.

The FASTQ files themselves are checked in two steps, before any reads are processed. Problems are collected across all libraries and reported together, so that all faulty files can be replaced at once.

1. **At start-up**, from the first read of each file: a file is empty, is named `.gz` but is not gzip-compressed, or is not in FASTQ format; `file1` and `file2` are the same file, start with the same read, or start with reads of different names; a file is used by more than one library.
2. **As the first step of the workflow** ([**fq lint**](https://github.com/stjude-rust-labs/fq)), on the complete files: a compressed file is truncated or corrupted (e.g. by an interrupted download or copy), a FASTQ record is malformed, `file1` and `file2` have different numbers of reads, or the read names of `file1` and `file2` do not match record by record (mates not in the same order). All other steps wait until every library has passed.

Mates must be in the same order in `file1` and `file2`, because vsearch pairs reads by position: out-of-order mates would not merge, or merge into chimeric reads that look like double mutants. If both files contain the same reads in a different order (for example, lanes concatenated in a different order), re-pair them before running the pipeline with [`seqkit pair`](https://github.com/shenwei356/seqkit) or BBTools [`repair.sh`](https://github.com/bbushnell/BBTools).

## Reference and reading frame

As of v1.0.0, all libraries in one run of `nf-core/deepmutscan` must stem from the same gene. Provide its wildtype sequence as a FASTA file via `--fasta`. We recommend that the reference contains N- and C-terminal flanking sequences of approximately one read length (e.g. from surrounding vector regions amplified during the sequencing library preparation), which helps to align reads more uniformly at the ends of the ORF.

Use `--reading_frame` to give the position of the mutagenised ORF within that sequence, as 1-based, inclusive nucleotide coordinates in the format `start-stop`: `start` is the first nucleotide of the first codon and `stop` the last nucleotide of the last codon to be analysed, so the length of the range must be divisible by three. For example, `--reading_frame 151-450` analyses exactly 100 codons of the reference FASTA.

Use `--mutagenesis_type` to define which codons were programmed at each position (default `nnk`). This determines which observed variants count as library variants for filtering, QC and fitness estimation:

- `nnk`, `nns`, `nnh`, `nnn`: degenerate codons, as used e.g. in nicking mutagenesis
- `nnk_nns`, `nnk_nns_nnh`: position-dependent degenerate codons, chosen by the third base of the wildtype codon: NNS where it is T, NNH where it is G (`nnk_nns_nnh` only), and NNK otherwise
- `custom`: a user-defined codon library, provided with `--custom_codon_library` as a `.csv` file. Either give one global comma-separated list of codons without a header (e.g. `AAA,AAC,AAG`), or a position-wise list: a header line containing the word `Position`, followed by one row per codon position with the position and its allowed codons (e.g. `1,ACG,AAA,ACA`).

> [!IMPORTANT]
> Shotgun sequencing cannot distinguish reads from unmutated plasmids from wildtype-matching fragments of mutated plasmids, so exact wildtype reads are not used as the fitness reference. Instead, the pipeline uses a synonymous variant as a wildtype proxy. We therefore highly recommend that your library design choice includes at least some synonymous wildtype amino acid codons – degenerate NNK/NNS codons do so at almost every position.

## Running the pipeline

The typical command for running the pipeline (here on a hypothetical gene of exactly 100 codons) is as follows:

```bash title="example_run.sh"
nextflow run nf-core/deepmutscan \
   -profile <docker/singularity/.../institute> \
   --input ./samplesheet.csv \
   --fasta ./ref.fa \
   --reading_frame 151-450 \
   --outdir ./results
```

This will launch the pipeline with the `docker` configuration profile (or the profile of your choice, see [below](#-profile)). It performs the default read QC, read merging, alignment, filtering, variant counting, sequencing-error correction and library QC for every library in the samplesheet. Add `--fitness` to also estimate variant fitness from matched input and output libraries, and `--dimsum` and/or `--mutscan` for the additional fitness estimators.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow.log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/deepmutscan -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
fasta: './ref.fa'
reading_frame: '151-450'
fitness: true
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

## Pipeline steps in detail

### 1. Read merging

All paired-end raw reads are first quality-checked with [**FastQC**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

Even the highest-quality sequencing platforms do not have perfect base accuracy. To minimise the effect of base errors, which would otherwise be counted as false mutations, `nf-core/deepmutscan` uses the overlap of each read pair. As base errors on the forward and reverse read are largely independent, [**vsearch --fastq_mergepairs**](https://github.com/torognes/vsearch) converts each raw read pair into a single consensus read with adjusted base qualities. Read pairs that cannot be merged (overlap shorter than 10 bp) are discarded. The merging statistics of every library are saved in `intermediate_files/merged_reads/`.

### 2. Alignment

The merged reads are aligned to the reference with [**BWA-MEM**](https://github.com/lh3/bwa). BWA-MEM is a highly efficient mapping algorithm for reads ≥ 70 bp.

### 3. Filtering

Alignments are filtered with [**samtools view**](https://www.htslib.org/doc/samtools-view.html): only primary alignments (no unmapped, secondary or supplementary records) with a mapping quality ≥ 30 and without insertions, deletions or splice gaps are kept, as indel-containing reads are most likely artefacts in substitution libraries. The filtered alignments are coordinate-sorted and indexed with samtools, and used for variant counting.

For long ORF site-saturation mutagenesis libraries, most aligned shotgun sequencing reads match the reference exactly. These wildtype reads are retained: they are ignored by the variant counter but contribute to the overall sequencing coverage and sequencing-error estimation.

Merging, alignment, filtering, sorting and indexing run as one streamed, multi-threaded task per library (`MERGE_ALIGN`, using its `task.cpus`). The vsearch options (`ext.args`) and the samtools filter flags (`ext.args2`) can be changed in a custom config.

> [!TIP]
> The merging benefit is largest if the average DNA fragment size matches the read length. For example, libraries sequenced with 150 bp paired-end reads should ideally also be sheared or tagmented to a mean size of about 150 bp.

### 4. Variant counting

Merged reads are screened for base-level mismatches against the reference. `nf-core/deepmutscan` uses a light-weight Python counter (built on [**pysam**](https://pysam.readthedocs.io) and [**polars**](https://pola.rs)) that counts all single, double, triple and higher-order nucleotide changes per read, together with the sequencing coverage of every position. This replaces the GATK `AnalyzeSaturationMutagenesis` tool used in earlier versions; a description of every column is written next to each count table (`variant_counts_columns.tsv`).

Two parameters tune the counter:

- `--base_qual` (default `40`): minimum Phred base quality for a mismatch to be counted, and for a base to contribute to coverage.
  Note what this threshold implies: `vsearch` recalculates qualities only where the two mates overlap, so a base in the overlap
  that both mates agree on can exceed the instrument's own maximum, while bases outside the overlap keep their instrument quality.
  On quality-binned instrument data (NovaSeq bins at Q37 or below) the default of `40` therefore restricts counting to the
  overlapping, double-observed part of each merged read. Lower it (e.g. `--base_qual 30`) to count the full read length.
- `--min_flank` (default `2`): minimum distance (in nucleotides) of a mismatch or covered base from the read end. Set to `0` to disable.

### 5. Library annotation and filtering

Variant counts are annotated with their codon and amino acid changes within `--reading_frame`, and filtered for the variants that were programmed by the chosen `--mutagenesis_type`. In addition, a library-completed table lists every programmed variant, including those that were not observed (with a near-zero placeholder count), for the count distribution plots.

### 6. Single-nucleotide variant error correction

Sequencing errors occur at typically low, but position- and substitution-specific rates. With very deep sequencing, they inflate the counts of variants that differ from the wildtype by only a single nucleotide. Variants with two or three nucleotide changes in one codon, however, are hardly affected. Select a correction strategy with `--error_correction`:

- `false_doubles` (**default**): the library only contains single-codon changes, so reads carrying a programmed multi-nucleotide codon variant plus an unprogrammed single-nucleotide change in a nearby codon ("false double mutants") reveal the sequencing error rate of that single-nucleotide change. The error rate is estimated for every single-nucleotide variant from these read-linked observations within `--false_doubles_codon_window` codons (default `40`), and the expected number of erroneous counts is subtracted. No extra data are required. If a library contains too few false double mutants to model their coverage along the read (e.g. very small or shallow libraries), no correction is applied and a warning is written to the log. Two estimators are available via `--false_doubles_method`:
  - `mle` (**default**): a per-variant maximum-likelihood error rate.
  - `eb`: an empirical-Bayes estimate that shrinks the per-variant rate towards the mean of its substitution class. It is steadier for sparsely observed variants, e.g. in shallowly sequenced libraries.
- `wildtype`: subtracts the variant-specific background measured by **additional deep sequencing of the unmutated template** from every library variant. Add these libraries to the samplesheet with `type` set to `wildtype` and the same `sample` name as the libraries they should correct. Exactly one `wildtype` library is supported; the pipeline stops with an error if there is none or more than one. In this mode, the `wildtype` libraries themselves are not included in the library QC and fitness outputs.
- `none`: no correction; counts are used as observed.

Corrected counts are used for every count-dependent step downstream (count heatmaps, positional QC and fitness). The uncorrected tables are kept next to the corrected ones, and an interactive HTML report shows the effect of the correction across all libraries of the run.

### 7. DMS library quality control

Based on the library-filtered (and error-corrected) counts, `nf-core/deepmutscan` produces per-library visualisations to assess mutation efficiency, coverage and saturation along the ORF:

- heatmaps of variant counts and counts per coverage, per position and amino acid
- sorted, log-scale distributions of counts per coverage of all programmed variants, also stratified by the number of changed nucleotides, summarised by the log10 ratio of their 90th and 10th percentiles
- sliding-window profiles of coverage and counts along the ORF (window size set by `--sliding_window_size`, default `10` codons), with a reference line for the coverage needed to observe every amino acid variant `--aimed_cov` times (default `100`)
- a sequencing-depth rarefaction curve that estimates the fraction of programmed variants recovered at lower sequencing depths (`--run_seqdepth`, on by default)

### 8. Fitness estimation (optional)

With `--fitness`, libraries of each `sample` are grouped into a merged count table of all `input` and `output` replicates. Fitness is then estimated per replicate as the natural logarithm of each variant's output-to-input count ratio, relative to that of a synonymous wildtype proxy variant. As the proxy, the pipeline picks the synonymous variant with two nucleotide changes in one codon (or, if there is none, one nucleotide change) with the highest mean input count. Non-synonymous variants encoding the same amino acid sequence are aggregated. Per-replicate fitness values are linearly rescaled so that the median of synonymous variants is `0` and the median of stop codon variants is `-1` (if one of these anchors is missing, fitness values are only centred), and summarised as the mean and, with more than one replicate, the standard deviation across replicates.

Two parameters control which variants are scored, and how:

- `--min_counts` (default `10`): minimum number of reads a variant needs in every `input` replicate to receive a fitness estimate. Variants with fewer input reads are poorly measured and are excluded. The threshold applies to each codon variant before synonymous codon variants are merged, in the same way for the default estimator, DiMSum and mutscan.
- `--output_pseudocount` (default `1`): pseudocount added to `output` replicates in which a variant observed in the input was not detected ("dropouts"), so that strongly depleted variants still receive a finite fitness value. Set it to `0` to disable; dropout variants then receive no fitness estimate in that replicate.

Two established statistical frameworks can be run on the same merged counts:

- `--dimsum`: [DiMSum](https://github.com/lehner-lab/DiMSum) fitness estimates (stages 4–5 of DiMSum, without its own read processing), on the same merged count table as the default estimator. The pipeline passes `--min_counts` as `--fitnessMinInputCountAll` and `--output_pseudocount` as `--fitnessDropoutPseudocount`, so DiMSum scores the same variants. By default, DiMSum also runs with:
  - `--mixedSubstitutions T`: the synonymous wildtype proxy can differ from the reference by one codon. With DiMSum's own default (`F`), single mutants relative to the proxy would be discarded as "mixed substitutions". With `T`, DiMSum also adds up synonymous codon variants before estimating fitness, as the default estimator does.
  - `--fitnessErrorModel F`: DiMSum's error model needs enough variants that are well covered in every replicate, which sparse shotgun libraries often lack. Without it, DiMSum uses count-based errors only and no replicate normalisation. For deeply sequenced libraries, the error model can be switched on.
  - `--retainIntermediateFiles T`.

  These options are set with `ext.args` and can be changed in a custom config, for example:

  ```groovy
  process {
      withName: 'RUN_DIMSUM' { ext.args = '--mixedSubstitutions T --fitnessErrorModel T --retainIntermediateFiles T' }
  }
  ```

  DiMSum's container is currently only available for x86-64 (amd64), so on ARM machines (e.g. Apple Silicon) `--dimsum` runs under amd64 emulation, and the `arm64` profile has no native image for it.

- `--mutscan`: [mutscan](https://github.com/fmicompbio/mutscan) enrichment estimates with edgeR and limma, on the same merged count table. mutscan is currently only packaged for x86-64 (amd64), so on ARM machines (e.g. Apple Silicon) `--mutscan` runs under amd64 emulation, and the `arm64` profile has no native image for it.

### Interactive variant effect inspection tool (`--pdb`)

When `--fitness` is set and a wildtype 3D structure is supplied via `--pdb <structure.pdb>` (an experimental structure, an AlphaFold DB model, or any PDB file matching the wildtype ORF), the pipeline additionally builds a portable, interactive HTML tool. It projects per-residue fitness, coverage, counts, counts per coverage and – when error correction is on – the positional error bias onto the structure; clicking a residue shows the effect of each substitution. Structure prediction is not part of the v1.0.0 release, so a user-specified structure file is required to enable the tool.

### 9. Reporting

`nf-core/deepmutscan` writes a single, self-contained `deepmutscan_report.html` to the top of the results directory. It embeds the library QC plots, the error-correction report, the fitness results (and DiMSum/mutscan output, if run), the [MultiQC](https://multiqc.info) report, run statistics and the citations of all tools that were used, so it can be shared as a single file.

## Updating the pipeline

When you run the original command above, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use this cached version if available - even if the pipeline has been updated since. To make sure that you are running the latest version of the pipeline, make sure that you regularly update the cached version:

```bash
nextflow pull nf-core/deepmutscan
```

## Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used consistently. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/deepmutscan releases page](https://github.com/nf-core/deepmutscan/releases) and find the latest pipeline version - numeric only (eg. `1.0.0`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.0.0`.

This version number will also be logged in reports when you run the pipeline, for example at the bottom of the MultiQC reports.

To further assist in reproducibility, you can share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such a profile (e.g. providing it as supplementary material for academic publications), make sure to _not_ include your cluster specific file paths or institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it may lead to varying or even irreproducible results across users' different computer environments.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results (from within the `/work` directory) from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~/.bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
