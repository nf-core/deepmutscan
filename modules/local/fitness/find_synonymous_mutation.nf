process FIND_SYNONYMOUS_MUTATION {
  tag { sample.sample }
  label 'process_single'

  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity'
      ? 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:ce2ba7ad7f6e7f2c'
      : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:0fd2e39a5bf2ecaa' }"

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
