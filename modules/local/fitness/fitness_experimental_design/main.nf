process EXPDESIGN_FITNESS {
  tag "experimentalDesign"
  label 'process_single'

  conda "${moduleDir}/environment.yml"

  container "${ workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/73/73a72ec77725aeb67678a74228938fdd6827b669d01a8c96951b1a8ef96eeb0f/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-mutscan_r-base_r-biocmanager_pruned:c65036d76406f342' }"

  input:
    path samplesheet_csv

  output:
    path "experimentalDesign.tsv", emit: experimental_design
    path "versions.yml",          emit: versions

  script:
  template 'dimsum_experimentalDesign.R'

  stub:
  """
  touch experimentalDesign.tsv
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      r-base: "0.0.0"
  END_VERSIONS
  """
}
