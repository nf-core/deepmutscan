process MERGE_COUNTS {
  tag "${sample.sample}"
  label 'process_single'

  conda "${moduleDir}/environment.yml"

  container "${ workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/78/7815f74d3680b875aa8c7be8f78dc77921ec970b99b45a539c5d2225749c0b90/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-dplyr_r-ggplot2_pruned:a65d6fd383d85714' }"

  input:
    tuple val(sample), path(input_counts), path(output_counts)

  output:
    tuple val(sample), path("counts_merged.tsv"), emit: merged_counts
    path "versions.yml", emit: versions

  script:
  template 'merge_counts.R'

  stub:
  """
  touch counts_merged.tsv
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      r-base: "0.0.0"
  END_VERSIONS
  """
}
