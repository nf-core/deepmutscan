process FITNESS_HEATMAP {
  tag { sample.sample }
  label 'process_single'

  conda "${moduleDir}/environment.yml"

  container "${ workflow.containerEngine == 'singularity'
        ? 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:ce2ba7ad7f6e7f2c'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:ce2ba7ad7f6e7f2c' }"

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
