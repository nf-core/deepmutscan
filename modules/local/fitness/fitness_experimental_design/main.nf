process EXPDESIGN_FITNESS {
  tag "experimentalDesign"
  label 'process_single'

  conda "${moduleDir}/environment.yml"

  container "${ workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a6/a66af652f7b1c6d9dab9f081f5676bd9452d653eed43034f20a2cf172921a4cf/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-mutscan_pysam_biopython_pruned:fb9a2095922ddd59' }"

  input:
    path samplesheet_csv

  output:
    path "experimentalDesign.tsv", emit: experimental_design
    path "versions.yml",          emit: versions

  script:
  template 'dimsum_experimentalDesign.R'
}
