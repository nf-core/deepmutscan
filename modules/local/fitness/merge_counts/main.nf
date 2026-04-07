process MERGE_COUNTS {
  tag "${sample.sample}"
  label 'process_single'

  conda "${moduleDir}/environment.yml"

  container "${ workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/73/73a72ec77725aeb67678a74228938fdd6827b669d01a8c96951b1a8ef96eeb0f/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-mutscan_r-base_r-biocmanager_pruned:c65036d76406f342' }"

  input:
    tuple val(sample), val(metas), path(input_counts), path(output_counts)

  output:
    tuple val(sample), path("counts_merged.tsv"), emit: merged_counts
    path "versions.yml", emit: versions

  script:
  template 'merge_counts.R'
}
