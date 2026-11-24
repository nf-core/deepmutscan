process RUN_MUTSCAN {
  tag "${sample.sample}"
  label 'process_single'

  conda "${moduleDir}/environment.yml"

  container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
      ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4a/4a4c1441760e6f27f01d8930305d8a9770e864b951b46cc4cf946e4509da6140/data'
      : 'community.wave.seqera.io/library/bioconductor-edger_bioconductor-limma_r-base_r-ggally_r-ggplot2:826c68a10b1231ea'}"

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
