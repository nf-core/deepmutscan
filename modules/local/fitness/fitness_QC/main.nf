process FITNESS_QC {
  tag { sample.sample }
  label 'process_single'

  conda "${moduleDir}/environment.yml"

  container "${ workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/78/7815f74d3680b875aa8c7be8f78dc77921ec970b99b45a539c5d2225749c0b90/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-dplyr_r-ggplot2_pruned:a65d6fd383d85714' }"

  input:
    tuple val(sample), path(fitness_estimation_tsv)   // from FITNESS_CALCULATION

  output:
    tuple val(sample), path("fitness_estimation_count_correlation.pdf"), emit: counts_corr_pdf
    tuple val(sample), path("fitness_estimation_fitness_correlation.pdf"), emit: fitness_corr_pdf
    path "versions.yml", emit: versions

  script:
  template 'fitness_QC.R'

  stub:
  """
    touch fitness_estimation_count_correlation.pdf
    touch fitness_estimation_fitness_correlation.pdf
    cat > versions.yml <<'EOF'
    "${task.process}":
      stub-version: "0.0.0"
    EOF
  """
}
