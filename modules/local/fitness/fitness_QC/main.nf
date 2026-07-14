process FITNESS_QC {
  tag { sample.sample }
  label 'process_single'

  conda "${moduleDir}/environment.yml"

  container "${ workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a6/a66af652f7b1c6d9dab9f081f5676bd9452d653eed43034f20a2cf172921a4cf/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-mutscan_pysam_biopython_pruned:fb9a2095922ddd59' }"

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
