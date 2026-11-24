process RUN_MUTSCAN {
  tag "${sample.sample}"
  label 'process_single'

  conda "${moduleDir}/environment.yml"

  container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
      ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e2/e2d657c330e5ac4dcc3b40ef513d377367e60c293e01003906877e40f6865d75/data'
      : 'community.wave.seqera.io/library/bioconductor-edger_bioconductor-impute_bioconductor-limma_dms-fdmcorrect_pruned:be9da4b29f56a3d0'}"

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
