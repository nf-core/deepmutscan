process EXPDESIGN_FITNESS {
  tag "experimentalDesign"
  label 'process_single'

  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity'
      ? 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:ce2ba7ad7f6e7f2c'
      : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:0fd2e39a5bf2ecaa' }"

  input:
    path samplesheet_csv
    path script   // R file that defines: make_dimsum_experimental_design(input_csv, out_path)

  output:
    path "experimentalDesign.tsv", emit: experimental_design
    path "versions.yml",          emit: versions

  script:
  """
  set -euo pipefail

  Rscript -e "source('$script'); make_dimsum_experimental_design('$samplesheet_csv', 'experimentalDesign.tsv')"

  R_VERSION=\$(R --version | head -n 1 | sed -E 's/^R version ([0-9.]+).*/\\1/')
  cat <<-END_VERSIONS > versions.yml
  DIMSUM_EXPDESIGN:
    r-base: \$R_VERSION
  END_VERSIONS
  """
}