process FITNESS_CALCULATION {
  tag { sample.sample }
  label 'process_single'

  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity'
      ? 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:ce2ba7ad7f6e7f2c'
      : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:0fd2e39a5bf2ecaa' }"

  input:
    tuple val(sample), path(counts_merged)
    path(exp_design)
    path(syn_wt_txt)

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
