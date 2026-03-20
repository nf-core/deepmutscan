process FITNESS_CALCULATION {
  tag { sample.sample }
  label 'process_single'
  label 'r_env'

  conda "${moduleDir}/environment.yml"

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
