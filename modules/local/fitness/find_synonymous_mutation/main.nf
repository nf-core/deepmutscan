process FIND_SYNONYMOUS_MUTATION {
  tag { sample.sample }
  label 'process_single'

  conda "${moduleDir}/environment.yml"

  container "${ workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a6/a66af652f7b1c6d9dab9f081f5676bd9452d653eed43034f20a2cf172921a4cf/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-mutscan_pysam_biopython_pruned:fb9a2095922ddd59' }"

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
