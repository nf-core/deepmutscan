# nf-core/dms: Detailed Pipeline Steps

This page provides an in-depth description of the data processing modules implemented in **nf-core/dms**. It is intended for advanced users and developers who want to understand the rationale behind design choices, explore implementation details, and consider potential future extensions.

------------------------------------------------------------------------

## Overview

The pipeline processes deep mutational scanning (DMS) sequencing data in
several stages:\
1. Alignment of reads to the reference ORF\
2. Filtering of wildtype and erroneous reads\
3. Read merging for error reduction\
4. Mutation counting\
5. DMS library quality control\
6. Data summarisation across samples\
7. Single nucleotide variant error correction *(in development)*\
8. Fitness estimation *(in development)*

Each step is explained below. Links are provided to the primary tools
and libraries used, where applicable.

------------------------------------------------------------------------

## 1. Alignment

All paired-end raw reads are aligned to the reference ORF using
[**bwa-mem**](http://bio-bwa.sourceforge.net/).\
This is a highly efficient mapping algorithm for reads ≥ 50 bp, with its
multi-threading support automatically handled by nf-core.

In future versions of this pipeline, we consider the use of
[**bwa-mem2**](https://github.com/bwa-mem2/bwa-mem2), which provides
similar alignment rates with a moderate speed increase (preprint).\
With the increasing diversity of sequencing platforms for DMS, new
throughput, read length, and error profiles may require further
alignment options to be implemented.

------------------------------------------------------------------------

## 2. Filtering

For long ORF site-saturation mutagenesis libraries, most aligned shotgun
sequencing reads contain exact matches against the reference.\
It is not possible to infer which of these stem from mutant versus
wildtype DNA molecules prior to fragmentation, hence they are filtered
out.

Similarly, erroneous reads with unexpected indels are also removed.\
To this end, we use a custom Python script that relies on the
[**pysam**](https://pysam.readthedocs.io/en/latest/) library.

------------------------------------------------------------------------

## 3. Read Merging

Even the highest-accuracy next-generation sequencing platforms do not
have perfect base accuracy.\
To minimise the effect of base errors (which would otherwise be counted
as "false mutations"), nf-core/dms uses the overlap of each aligned read
pair.

With base errors on the forward and reverse read being independent, the
pipeline applies the [**vsearch
mergereads**](https://github.com/torognes/vsearch) function to convert
each read pair into a single consensus molecule with adjusted base error
scores.

> **Tip**: Optimal merging performance is usually obtained if the
> average DNA fragment size matches the read size. For example,
> libraries sequenced with 150 bp paired-end reads should ideally also
> be sheared/tagmented to a mean size of 150 bp.

Future versions may offer additional options depending on sequencing
type and error profiles.

------------------------------------------------------------------------

## 4. Mutation Counting

Aligned, non-wildtype consensus reads are screened for exact, base-level
mismatches.\
nf-core/dms uses a custom Python implementation to count the occurrences
of all high-quality single, double, triple, and higher-order nucleotide
changes between each read and the reference ORF.

Unlike the popular [**GATK
AnalyzeSaturationMutagenesis**](https://gatk.broadinstitute.org/hc/en-us/articles/360037057952-AnalyzeSaturationMutagenesis),
nf-core/dms allows users to specify a minimum base quality cutoff for
mutations to be included in the final count table (default: Q30).

------------------------------------------------------------------------

## 5. DMS Library Quality Control

By integrating the reference ORF coordinates and the chosen DMS library
type (default: NNK/NNs degenerate codon-based nicking), nf-core/dms
calculates a number of mutation count summary statistics.

Custom visualisations allow for inspection of: 1. Mutation efficiency
along the ORF\
2. Position-specific recovery of amino acid diversity\
3. Overall sequencing coverage evenness and saturation

------------------------------------------------------------------------

## 6. Data Summarisation for Fitness Estimation

Steps 1--5 are iteratively run across all samples defined in the `.csv`
spreadsheet.\
Once read alignment, merging, mutation counting, and library QC have
been completed for the full list of samples, users can opt to proceed
with fitness estimation using a number of tools (default:
[**DiMSum**](https://github.com/lehner-lab/DiMSum)).

The pipeline generates the tool-specific input files by merging mutation
counts across samples.

------------------------------------------------------------------------

## 7. Single Nucleotide Variant Error Correction *(in development)*

This module will implement strategies to distinguish true low-frequency
variants from sequencing artefacts.\
Methods under consideration include:\
- Empirical error rate modelling based on wildtype controls\
- Position-specific background correction\
- Integration of UMIs (unique molecular identifiers), where available

------------------------------------------------------------------------

## 8. Fitness Estimation *(in development)*

The final step of the pipeline will perform fitness estimation based on
mutation counts. Planned features include:\
- Integration of multiple fitness inference tools beyond DiMSum\
- Options for replicates and condition-specific comparisons\
- Standardised output formats for downstream analysis and visualisation

------------------------------------------------------------------------

## Notes for Developers

-   Custom scripts used in filtering and mutation counting are available
    in the `bin/` directory of the repository.\
-   Modules are implemented in Nextflow DSL2 and follow the nf-core
    community guidelines.\
-   Contributions, optimisations, and additional analysis modules are
    welcome --- please open a pull request or GitHub issue to discuss
    ideas.

------------------------------------------------------------------------

*This document is meant as a living reference. As the pipeline evolves,
the descriptions of steps 7 and 8 will be expanded with concrete
implementation details.*
