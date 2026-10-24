process FITNESS_HEATMAP {
  tag { sample.sample }
  label 'process_single'

  conda "${moduleDir}/environment.yml"

  container "${ workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/78/7815f74d3680b875aa8c7be8f78dc77921ec970b99b45a539c5d2225749c0b90/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-dplyr_r-ggplot2_pruned:a65d6fd383d85714' }"

  input:
    tuple val(sample), path(fitness_estimation_tsv)   // from FITNESS_CALCULATION
    tuple val(sample2), path(wt_seq)   		      // WT sequence

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
