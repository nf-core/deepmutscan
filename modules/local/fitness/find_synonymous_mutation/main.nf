process FIND_SYNONYMOUS_MUTATION {
  tag { sample.sample }
  label 'process_single'

  conda "${moduleDir}/environment.yml"

  container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
      ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/7988a2e369e98cec897d48bab3e55daddd8a973d9b7491ad27b6b15d550a4c5b/data'
      : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-dplyr_r-ggplot2_pruned:32f1af0a83a018e7'}"

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
