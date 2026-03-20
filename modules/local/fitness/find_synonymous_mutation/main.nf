process FIND_SYNONYMOUS_MUTATION {
  tag { sample.sample }
  label 'process_single'
  label 'r_env'

  conda "${moduleDir}/environment.yml"

  input:
    tuple val(sample), path(counts_merged)   // from MERGE_COUNTS.out.merged_counts
    path wt_fasta                            // broadcast singleton
    val pos_range                            // "start-end", broadcast singleton

  output:
    tuple val(sample), path("synonymous_wt.txt"), emit: synonymous_wt
    path "versions.yml", emit: versions

  script:
  template 'find_syn_mutation.R'
}
