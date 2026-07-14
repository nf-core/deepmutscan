process FITNESS_CALCULATION {
  tag { sample.sample }
  label 'process_single'

  conda "${moduleDir}/environment.yml"

  container "${ workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a6/a66af652f7b1c6d9dab9f081f5676bd9452d653eed43034f20a2cf172921a4cf/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-mutscan_pysam_biopython_pruned:fb9a2095922ddd59' }"

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
