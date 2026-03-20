process RUN_MUTSCAN {
  tag "${sample.sample}"
  label 'process_medium'

  conda "${moduleDir}/environment.yml"

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
