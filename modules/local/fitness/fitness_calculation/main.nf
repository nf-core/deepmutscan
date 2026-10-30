process FITNESS_CALCULATION {
  tag { sample.sample }
  label 'process_single'

  conda "${moduleDir}/environment.yml"

  container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
      ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/7988a2e369e98cec897d48bab3e55daddd8a973d9b7491ad27b6b15d550a4c5b/data'
      : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-dplyr_r-ggplot2_pruned:32f1af0a83a018e7'}"

  input:
    tuple val(sample), path(counts_merged)
    path(exp_design)
    path(syn_wt_txt)
    val(min_counts)
    val(output_pseudocount)

  output:
    tuple val(sample), path("fitness_estimation.tsv"), emit: fitness_estimation
    path "versions.yml", emit: versions

  script:
  template 'fitness_calculation.R'

  stub:
  """
    touch fitness_estimation.tsv
    echo "FITNESS_CALCULATION:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
  """
}
