process RUN_MUTSCAN {
  tag "${sample.sample}"
  label 'process_single'

  conda "${moduleDir}/environment.yml"

  container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
      ? 'REPLACE_WITH_SINGULARITY_HTTPS_URL'
      : 'REPLACE_WITH_DOCKER_IMAGE'}"

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
  template 'fitness_mutscan.R'

  stub:
  """
  touch fitness_estimation_mutscan_edgeR.tsv fitness_estimation_mutscan_limma.tsv mutscan_stub.pdf
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      r-base: "0.0.0"
      r-edgeR: "0.0.0"
      r-limma: "0.0.0"
  END_VERSIONS
  """
}
