process MERGE_COUNTS {
  tag "${sample.sample}"
  label 'process_single'

  conda "${moduleDir}/environment.yml"

  container "${ workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a6/a66af652f7b1c6d9dab9f081f5676bd9452d653eed43034f20a2cf172921a4cf/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-mutscan_pysam_biopython_pruned:fb9a2095922ddd59' }"

  input:
    tuple val(sample), path(input_counts), path(output_counts)

  output:
    tuple val(sample), path("counts_merged.tsv"), emit: merged_counts
    path "versions.yml", emit: versions

  script:
  template 'merge_counts.R'
}
