process MERGE_COUNTS {
  tag "${sample.sample}"
  label 'process_single'
  label 'r_env'

  conda "${moduleDir}/environment.yml"

  input:
    tuple val(sample), val(metas), path(input_counts), path(output_counts)

  output:
    tuple val(sample), path("counts_merged.tsv"), emit: merged_counts
    path "versions.yml", emit: versions

  script:
  template 'merge_counts.R'
}
