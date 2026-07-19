# nf-core/deepmutscan: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/deepmutscan/usage](https://nf-co.re/deepmutscan/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

**nf-core/deepmutscan** is a workflow designed for the analysis of deep mutational scanning (DMS) data. DMS enables researchers to experimentally measure the fitness effects of thousands of genes or gene variants simultaneously, helping to classify disease causing mutants in human and animal populations, to learn the fundamental rules of virus evolution, protein architecture, splicing, small-molecule interactions and many other phenotypes.

This page provides in-depth descriptions of the data processing modules implemented in `nf-core/deepmutscan`. It is similarly intended for new deep mutational scanning aficionados, for advanced users and developers who want to understand the rationale behind certain design choices, to explore implementation details and consider potential future extensions.

## Running the pipeline

The typical command for running the pipeline (on an example protein-coding gene with 100 amino acids) is as follows:

```bash title="example_run.sh"
nextflow run nf-core/deepmutscan \
   -profile <docker/singularity/.../institute> \
   --input ./samplesheet.csv \
   --fasta ./ref.fa \
   --reading_frame 1-300 \
   --outdir ./results
```

The `-profile <PROFILE>` specification is mandatory and should reflect either your own institutional profile or any pipeline profile specified in the [profile section](##-profile).

This will launch the pipeline by performing sequencing read alignments, various raw data QC analyses, optional mutant count error corrections, fitness and fitness error estimations.

Note that the pipeline will create the following files in your working directory:

```console title="working directory"
work               # Directory containing the nextflow working files
results            # Finished results in specified location (defined with --outdir); needs full writing access
.nextflow_log      # Log file from Nextflow
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash title="example_run.sh"
nextflow run nf-core/deepmutscan -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
gene reference: 'ref.fa'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

## Inputs

Users need to first prepare a samplesheet with your input/output data in which each row represents a pair of matched fastq files (paired end). This should look as follows:

```csv title="samplesheet.csv"
sample,type,replicate,file1,file2
ORF1,input,1,/reads/forward1.fastq.gz,/reads/reverse1.fastq.gz
ORF1,input,2,/reads/forward2.fastq.gz,/reads/reverse2.fastq.gz
ORF1,output,1,/reads/forward3.fastq.gz,/reads/reverse3.fastq.gz
ORF1,output,2,/reads/forward4.fastq.gz,/reads/reverse4.fastq.gz
```

Secondly, users need to specify the gene or gene region of interest using a reference FASTA file via `--fasta`. Provide the exact codon coordinates using `--reading_frame`.

## Optional parameters

Several optional parameters are available for `nf-core/deepmutscan`, some of which are currently _(in development)_.

| Parameter            | Default         | Description                                                                    |
| -------------------- | --------------- | ------------------------------------------------------------------------------ |
| `--run_seqdepth`     | `true`          | Estimate sequencing-depth saturation by closed-form hypergeometric rarefaction |
| `--fitness`          | `false`         | Default fitness inference module                                               |
| `--dimsum`           | `false`         | Optional fitness inference module _(AMD/x86_64 systems only)_                  |
| `--mutagenesis`      | `nnk`           | Deep mutational scanning strategy used                                         |
| `--error-estimation` | `wt_sequencing` | Error model used to correct 1nt counts _(in development)_                      |
| `--read-align`       | `bwa-mem`       | Customised read aligner _(in development)_                                     |

## Pipeline output

After execution, the pipeline creates the following directory structure:

```tree title="nf-core/deepmutscan results"
results/
├── fastqc/              # Individual HTML reports for specified fastq files, raw sequencing QC
├── fitness/             # Merged variant count tables, fitness and error estimates, replicate correlations and heatmaps
├── intermediate_files/  # Raw alignments, raw and pre-filtered variant count tables, QC reports
├── library_QC/          # Sample-specific PDF visualizations: position-wise sequencing coverage, count heatmaps, etc.
├── multiqc/             # Shared HTML reports for all fastq files, raw sequencing QC
├── pipelineinfo/        # Nextflow helper files for timeline and summary report generation
├── timeline.html        # Nextflow timeline for all tasks
└── report.html          # Nextflow summary report incl. detailed CPU and memory usage per for all tasks
```

## Detailed steps

### 1. Alignment

All paired-end raw reads are first aligned to the provided reference ORF using [**bwa-mem**](http://bio-bwa.sourceforge.net/). This is a highly efficient mapping algorithm for reads ≥100 bp, with its multi-threading support automatically handled by nf-core.

In future versions of `nf-core/deepmutscan`, we consider the use of [**bwa-mem2**](https://github.com/bwa-mem2/bwa-mem2), which provides similar alignment quality at moderate speed increase ([Vasimuddin et al., _IPDPS_ 2019](https://ieeexplore.ieee.org/document/8820962)). With the increasing diversity of sequencing platforms for DMS new read length, throughput and error profiles may require further alignment options to be implemented.

### 2. Filtering

For long ORF site-saturation mutagenesis libraries, most aligned shotgun sequencing reads contain exact matches against the reference. Reads with likely artefactual indel-containing alignments are removed, while exact-match (wildtype) reads are retained — they are ignored by the variant counter and kept available for downstream sequencing-error correction.

To this end, we use [**samtools view**](https://www.htslib.org/doc/samtools.html).

### 3. Read Merging

Even the highest-quality next-generation sequencing platforms do not feature perfect base accuracy. To minimise the effect of base errors (which would otherwise be counted as "false mutations"), `nf-core/deepmutscan` uses the overlap of each aligned read pair. With base errors on the forward and reverse read being independent, the pipeline applies the [**vsearch fastq_mergepairs**](https://github.com/torognes/vsearch) function to convert each read pair into a single consensus molecule with adjusted base error scores.

> [!TIP]
> Optimal merging benefit is usually obtained if the average DNA fragment size matches the read size. For example, libraries sequenced with 150 bp paired-end reads should ideally also be sheared/tagmented to a mean size of 150 bp.

Future versions may offer additional options depending on sequencing type and error profiles.

### 4. Variant Counting

Aligned consensus reads are screened for exact, base-level mismatches. `nf-core/deepmutscan` uses a lightweight, dependency-free Python counter (built on [**pysam**](https://pysam.readthedocs.io) and [**polars**](https://pola.rs)) to count occurrences of all single, double, triple, and higher-order nucleotide changes between each read and the reference ORF, from a coordinate-sorted, indexed BAM. It produces a variant count table that is column-compatible with the previously used GATK `AnalyzeSaturationMutagenesis` output.

Two options tune the counter: `--base_qual` sets the minimum base quality for a mutation to be counted (default `30`; not available in GATK), and `--min_flank` sets the minimum distance (in nt) a mutation or covered base must be from a read edge (default `2`, matching GATK's `--min-flanking-length`; set to `0` to disable). The `--min_flank` window is applied consistently to both variant calls and the coverage calculation.

### 5. DMS Library Quality Control

By integrating the reference ORF coordinates and the chosen DMS library type (default: NNK degenerate codons), `nf-core/deepmutscan` calculates a number of mutation count summary statistics.

Custom visualisations allow for inspection of (1) mutation efficiency along the ORF, (2) position-specific recovery of mutant amino acid diversity, and (3) overall sequencing coverage evenness and library saturation.

### 6. Data Summarisation for Fitness Estimation

Steps 1-5 are run in parallel across all individual samples defined in the `.csv` spreadsheet. Once read alignment, filtering, merging, variant counting, and DMS library QC have been completed for the full list of samples – if input/output sample pairs are available – users can opt to proceed towards fitness estimation. To this end, the pipeline generates all the necessary preparatory files by generating a merged mutation count table across samples.

### 7. Single Nucleotide Variant Error Correction

Very deep sequencing introduces position-dependent false counts. Select a strategy with `--error_correction`:

- `false_doubles` (**default**): the library only contains single-codon changes, so observed multi-codon variants ("false doubles") are sequencing errors. Their counts, together with a distance-dependent coverage model, estimate the per-nucleotide error rate, which is subtracted from the single-codon counts. No extra data required. Two estimators are available via `--false_doubles_method`:

  - `mle` (**default**): a per-variant maximum-likelihood error rate.
  - `eb`: an empirical-Bayes estimate that pools across variants and shrinks per substitution class, which is steadier for sparsely observed variants.

  The matching window (in codons) is set by `--false_doubles_codon_window` (default 40).

- `none`: disables error correction (counts pass through unchanged).
- `wildtype`: subtracts a position-specific error profile measured from an **additional deep wildtype-only sequencing** sample. Provide it in the samplesheet with `type: wildtype`, sharing the same `sample` name as the input/output samples it should correct. This replaces the false-doubles correction.

The corrected counts are used for every count-dependent step (fitness, count heatmaps, positional-bias and coverage QC). Whenever correction is applied, a self-contained `*_error_correction_report.html` is written per sample so its effect can be inspected interactively.

### Interactive variant effect inspection tool (`--pdb`)

When `--fitness` is set and a wildtype 3D structure is supplied via `--pdb <structure.pdb>` (a crystal structure, an AlphaFold DB model, or any PDB matching the wildtype ORF), the pipeline additionally builds a single, portable, interactive HTML tool per sample. It projects per-residue fitness, coverage, counts, counts/coverage and — when error correction is on — the positional error bias onto the structure; clicking a residue shows each substitution's effect. Structure prediction (nf-core/proteinfold) is not wired in this release, so `--pdb` is required to enable the tool.

### 8. Fitness Estimation

The final step of the pipeline will perform fitness estimation based on mutation counts. By default, we calculate fitness scores as the logarithm of variants' output to input ratio, normalised to that of the provided wildtype nucleotide sequence.

Future expansions may include:

- Integration of other popular fitness inference tools, including [DiMSum](https://github.com/lehner-lab/DiMSum), [Enrich2](https://github.com/FowlerLab/Enrich2), [rosace](https://github.com/pimentellab/rosace/) and [mutscan](https://github.com/fmicompbio/mutscan)
- Standardised output formats for downstream analyses and comparison

> [!IMPORTANT]
> We note that exact wildtype sequence reads are filtered out in stage 2. Including synonymous wildtype codons in the original mutagenesis design is therefore essential when it comes to calibrating the fitness calculations.

## Notes for Developers

- Custom R and Python scripts used in variant counting, filtering, error correction and QC visualisation live in each module's own `templates/` directory (e.g. `modules/local/dmsanalysis/*/templates/`).
- Modules are implemented in Nextflow DSL2 and follow the nf-core community guidelines.
- Contributions, optimisations, and additional analysis modules are welcome: please open a Github [issue](https://github.com/nf-core/deepmutscan/issues/new) or [pull request](https://github.com/nf-core/deepmutscan/compare) to discuss or suggest ideas.

_This document is meant as a living reference. As the pipeline evolves, the descriptions of steps 7 and 8 will be further expanded with concrete implementation details._

## Updating the pipeline

When you run the original command above, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use this cached version if available - even if the pipeline has been updated since. To make sure that you are running the latest version of the pipeline, make sure that you regularly update the cached version:

```bash
nextflow pull nf-core/deepmutscan
```

## Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/deepmutscan releases page](https://github.com/nf-core/deepmutscan/releases) and find the latest pipeline version - numeric only (eg. `1.0.0`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.0.0`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

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
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
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
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
