# nf-core/deepmutscan: Detailed Pipeline Steps

This page provides in-depth descriptions of the data processing modules implemented in **nf-core/deepmutscan**. It is mainly intended for advanced users and developers who want to understand the rationale behind design choices, explore implementation details, and consider potential future extensions.

---

## Overview

The pipeline processes deep mutational scanning (DMS) sequencing data in several stages:

1. Alignment of reads to the reference open reading frame (ORF)
2. Filtering of wildtype and erroneous reads
3. Read merging for base error reduction
4. Mutation counting
5. DMS library quality control
6. Data summarisation across samples
7. Single nucleotide variant error correction
8. Fitness estimation _(in development)_

![pipeline](/docs/pipeline.png)

Each step is explained below. Links are provided to the primary tools and libraries used, where applicable.

---

## 1. Alignment

All paired-end raw reads are first aligned to the provided reference ORF using [**bwa-mem**](http://bio-bwa.sourceforge.net/). This is a highly efficient mapping algorithm for reads ≥100 bp, with its multi-threading support automatically handled by nf-core.

In future versions of nf-core/deepmutscan, we consider the use of [**bwa-mem2**](https://github.com/bwa-mem2/bwa-mem2), which provides similar alignment rates with a moderate speed increase ([Vasimuddin et al., _IPDPS_ 2019](https://ieeexplore.ieee.org/document/8820962)). With the increasing diversity of sequencing platforms for DMS, new throughput, read length, and error profiles may require further alignment options to be implemented.

---

## 2. Filtering

For long ORF site-saturation mutagenesis libraries, most aligned shotgun sequencing reads contain exact matches against the reference. Erroneous reads with unexpected indels are removed, while exact-match (wildtype) reads are retained — they are ignored by the variant counter and kept available for downstream sequencing-error correction.

To this end, we use [**samtools view**](https://www.htslib.org/doc/samtools.html).

---

## 3. Read Merging

Even the highest-accuracy next-generation sequencing platforms do not have perfect base accuracy. To minimise the effect of base errors (which would otherwise be counted as "false mutations"), nf-core/deepmutscan uses the overlap of each aligned read pair. With base errors on the forward and reverse read being independent, the pipeline applies the [**vsearch fastq_mergepairs**](https://github.com/torognes/vsearch) function to convert each read pair into a single consensus molecule with adjusted base error scores.

> [!TIP]
> Optimal merging performance is usually obtained if the average DNA fragment size matches the read size. For example, libraries sequenced with 150 bp paired-end reads should ideally also be sheared/tagmented to a mean size of 150 bp.

Future versions may offer additional options depending on sequencing type and error profiles.

---

## 4. Variant Counting

Aligned consensus reads are screened for exact, base-level mismatches. nf-core/deepmutscan uses a lightweight, dependency-free Python counter (built on [**pysam**](https://pysam.readthedocs.io) and [**polars**](https://pola.rs)) to count occurrences of all single, double, triple, and higher-order nucleotide changes between each read and the reference ORF, from a coordinate-sorted, indexed BAM. The counter replaces GATK `AnalyzeSaturationMutagenesis` (keeping the container light) while producing a column-compatible count table.

Two options tune the counter: `--base_qual` sets the minimum base quality for a mutation to be counted (default Q30; not available in GATK), and `--min_flank` sets the minimum distance (nt) a mutation or covered base must lie from a read edge (default 2, matching GATK's `--min-flanking-length`; 0 disables). The `--min_flank` window is applied to both variant calls and coverage.

---

## 5. DMS Library Quality Control

By integrating the reference ORF coordinates and the chosen DMS library type (default: NNK/NNs degenerate codon-based nicking), nf-core/deepmutscan calculates a number of mutation count summary statistics.

Custom visualisations allow for inspection of (1) mutation efficiency along the ORF, (2) position-specific recovery of amino acid diversity, and (3) overall sequencing coverage evenness and saturation.

---

## 6. Data Summarisation for Fitness Estimation

Steps 1-5 are iteratively run across all samples defined in the `.csv` spreadsheet. Once read alignment, merging, mutation counting, and library QC have been completed for the full list of samples, users can opt to proceed with fitness estimation. To this end, the pipeline generates all the necessary input files by merging mutation counts across samples.

---

## 7. Single Nucleotide Variant Error Correction

With very deep sequencing, single-codon counts show position-dependent bias from sequencing errors. `nf-core/deepmutscan` distinguishes true single-nucleotide variants from artefacts via the `--error_correction` parameter, with three options:

- `false_doubles` (**default**): because the library only contains single-codon changes, any observed multi-codon variant ("false double") must arise from sequencing error. Their counts, with a distance-dependent coverage model, estimate the per-nucleotide error rate that is subtracted from the single-codon counts. Needs no extra data. `--false_doubles_method` picks the estimator: `mle` (default, per-variant maximum likelihood) or `eb` (empirical Bayes, pooled and shrunk per substitution class — steadier for sparse variants). `--false_doubles_codon_window` (default 40) sets the codon matching window.
- `none`: no correction.
- `wildtype`: empirical error-rate modelling from **additional deep wildtype-only sequencing**. Add the wildtype sample(s) to the samplesheet with `type: wildtype`; their position-specific error profile is subtracted from each library sample of the same `sample`. This replaces the false-doubles correction.

The corrected counts feed the fitness estimation and the mutation-count heatmaps; the corrected count tables are written next to the uncorrected ones under `intermediate_files/processed_variant_counts/`.

---

## 8. Fitness Estimation _(in development)_

The final step of the pipeline will perform fitness estimation based on mutation counts. By default, we calculate fitness scores as the logarithm of variants' output to input ratio, normalised to that of the provided wildtype sequence. Future expansions may include:

- Integration of other popular fitness inference tools, including [DiMSum](https://github.com/lehner-lab/DiMSum), [Enrich2](https://github.com/FowlerLab/Enrich2), [rosace](https://github.com/pimentellab/rosace/) and [mutscan](https://github.com/fmicompbio/mutscan)
- Standardised output formats for downstream analyses and comparison

> [!IMPORTANT]
> We note that exact wildtype sequence reads are filtered out in stage 2. Including synonymous wildtype codons in the original mutagenesis design is therefore essential when it comes to calibrating the fitness calculations.

---

## Notes for Developers

- Custom scripts used in filtering and mutation counting are available in the `bin/` directory of the repository.
- Modules are implemented in Nextflow DSL2 and follow the nf-core community guidelines.
- Contributions, optimisations, and additional analysis modules are welcome - please open a pull request or GitHub issue to discuss ideas.

_This document is meant as a living reference. As the pipeline evolves, the descriptions of steps 7 and 8 will be expanded with concrete implementation details._

---

## Contact

For detailled scientific or technical questions, feedback and experimental discussions, feel free to contact us directly:

- Benjamin Wehnert — wehnertbenjamin@gmail.com
- Maximilian Stammnitz — maximilian.stammnitz@crg.eu
