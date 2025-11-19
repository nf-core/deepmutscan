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



process FITNESS_QC {
  tag { sample.sample }
  label 'process_single'

  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity' 
      ? 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:ce2ba7ad7f6e7f2c' 
      : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:0fd2e39a5bf2ecaa' }"

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
    FITNESS_PLOTS:
      stub-version: "0.0.0"
    EOF
  """
}



process FITNESS_HEATMAP {
  tag { sample.sample }
  label 'process_single'

  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity' 
      ? 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:ce2ba7ad7f6e7f2c' 
      : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:0fd2e39a5bf2ecaa' }"

  input:
    tuple val(sample), path(fitness_estimation_tsv)   // from FITNESS_CALCULATION
    tuple val(sample), path(wt_seq)   		      // WT sequence

  output:
    tuple val(sample), path("fitness_heatmap.pdf"), emit: fitness_heatmap
    path "versions.yml", emit: versions

  script:
    template 'fitness_heatmap.R'

  stub:
  """
    touch fitness_heatmap.pdf
    cat > versions.yml <<'EOF'
    FITNESS_HEATMAP:
      stub-version: "0.0.0"
    EOF
  """
}
