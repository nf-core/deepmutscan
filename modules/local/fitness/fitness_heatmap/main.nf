process FITNESS_HEATMAP {
  tag { sample.sample }
  label 'process_single'
  label 'r_env'

  conda "${moduleDir}/environment.yml"

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
