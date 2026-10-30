process RUN_MUTSCAN {
  tag "${sample.sample}"
  label 'process_single'

  conda "${moduleDir}/environment.yml"

  container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
      ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0a/0a32a5f637634c064a3caecc77d9e99f4367efd459bc28d1d4971ed5d05b73a7/data'
      : 'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-mutscan_r-base:c5cefc95c9f1002e'}"

  input:
    tuple val(sample), path(counts_merged)
    path(syn_wt_txt)
    path(exp_design)
    val(min_counts)

  output:
    tuple val(sample), path("fitness_estimation_mutscan_edgeR.tsv"), emit: fitness_mutscan_edgeR
    tuple val(sample), path("fitness_estimation_mutscan_limma.tsv"), emit: fitness_mutscan_limma
    tuple val(sample), path("*.pdf"), emit: qc_plots
    path "versions.yml", emit: versions

  script:
  template 'fitness_calculation_mutscan.R'

  stub:
  """
  touch fitness_estimation_mutscan_edgeR.tsv fitness_estimation_mutscan_limma.tsv mutscan_stub.pdf
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      r-base: "0.0.0"
  END_VERSIONS
  """
}
