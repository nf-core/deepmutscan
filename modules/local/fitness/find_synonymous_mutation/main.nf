process FIND_SYNONYMOUS_MUTATION {
  tag { sample.sample }
  label 'process_single'

  conda "${moduleDir}/environment.yml"

  container "${ workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/78/7815f74d3680b875aa8c7be8f78dc77921ec970b99b45a539c5d2225749c0b90/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-dplyr_r-ggplot2_pruned:a65d6fd383d85714' }"

  input:
    tuple val(sample), path(counts_merged)   // from MERGE_COUNTS.out.merged_counts
    path wt_fasta                            // broadcast singleton
    val pos_range                            // "start-end", broadcast singleton

  output:
    tuple val(sample), path("synonymous_wt.txt"), emit: synonymous_wt
    path "versions.yml", emit: versions

  script:
  template 'find_syn_mutation.R'

  stub:
  """
  touch synonymous_wt.txt
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      r-base: "0.0.0"
  END_VERSIONS
  """
}
