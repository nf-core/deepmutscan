process EXPDESIGN_FITNESS {
  tag "experimentalDesign"
  label 'process_single'
  label 'r_env'

  conda "${moduleDir}/environment.yml"

  input:
    path samplesheet_csv

  output:
    path "experimentalDesign.tsv", emit: experimental_design
    path "versions.yml",          emit: versions

  script:
  template 'dimsum_experimentalDesign.R'
}
