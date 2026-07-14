process RUN_MUTSCAN {
  tag "${sample.sample}"
  label 'process_medium'

  conda "${moduleDir}/environment.yml"

  container "${ workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a6/a66af652f7b1c6d9dab9f081f5676bd9452d653eed43034f20a2cf172921a4cf/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-mutscan_pysam_biopython_pruned:fb9a2095922ddd59' }"

  input:
    tuple val(sample), path(counts_merged)
    path(syn_wt_txt)
    path(exp_design)

  output:
    tuple val(sample), path("fitness_estimation_mutscan_edgeR.tsv"), emit: fitness_mutscan_edgeR
    tuple val(sample), path("fitness_estimation_mutscan_limma.tsv"), emit: fitness_mutscan_limma
    tuple val(sample), path("*.pdf"), emit: qc_plots
    path "versions.yml", emit: versions

  script:
  template 'fitness_calculation_mutscan.R'
}
