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
    path script  // fitness_calculation.R

  output:
    tuple val(sample), path("fitness_estimation.tsv"), emit: fitness_estimation
    path "versions.yml", emit: versions

  script:
  """
    set -euo pipefail

    R_version=\$(R --version | head -n 1 | sed 's/^R version //')
    
    Rscript -e "source('$script'); run_fitness_estimation('$counts_merged', '$exp_design', '$syn_wt_txt', 'fitness_estimation.tsv')"

    cat > versions.yml <<EOF
    FITNESS_CALCULATION:
      r: "\$R_version"
    EOF
  """

  stub:
  """
    touch fitness_estimation.tsv
    echo "FITNESS_CALCULATION:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
  """
}



process FITNESS_QC {
  tag { sample.sample }
  label 'process_single'

  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity' 
      ? 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:ce2ba7ad7f6e7f2c' 
      : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:0fd2e39a5bf2ecaa' }"

  input:
    tuple val(sample), path(fitness_estimation_tsv)   // from FITNESS_CALCULATION
    path script                                       // fitness_plots.R

  output:
    tuple val(sample), path("fitness_estimation_count_correlation.pdf"), emit: counts_corr_pdf
    tuple val(sample), path("fitness_estimation_fitness_correlation.pdf"), emit: fitness_corr_pdf
    path "versions.yml", emit: versions

  script:
  """
    set -euo pipefail

    R_version=\$(R --version | head -n 1 | sed 's/^R version //')
    
    Rscript -e "source('$script'); run_fitness_plots(
    '$fitness_estimation_tsv',
    'fitness_estimation_count_correlation.pdf',
    'fitness_estimation_fitness_correlation.pdf')"

    cat > versions.yml <<EOF
    FITNESS_CALCULATION:
      r: "\$R_version"
    EOF
  """

  stub:
  """
    touch fitness_estimation_count_correlation.pdf
    touch fitness_estimation_fitness_correlation.pdf
    cat > versions.yml <<'EOF'
    FITNESS_PLOTS:
      stub-version: "0.0.0"
    EOF
  """
}
