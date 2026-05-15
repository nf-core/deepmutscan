process FITNESS_HEATMAP {
  tag { sample.sample }
  label 'process_single'

  conda "${moduleDir}/environment.yml"

  container "${ workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/73/73a72ec77725aeb67678a74228938fdd6827b669d01a8c96951b1a8ef96eeb0f/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-mutscan_r-base_r-biocmanager_pruned:c65036d76406f342' }"

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
