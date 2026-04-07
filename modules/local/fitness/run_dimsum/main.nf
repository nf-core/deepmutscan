process RUN_DIMSUM {
  tag { sample.sample }
  label 'process_single'

  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity'
      ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/9298c391f285d7f2c1155f6a519a96c4f7591971780ba4f10569558282f40b6f/data'
      : 'community.wave.seqera.io/library/bioconductor-biostrings_r-dimsum_r-base_r-data.table_pruned:cca13eed371a9d84' }"

  input:
    tuple val(sample), path(counts_merged)
    path(wt_txt)
    path(exp_design)

  output:
    path "dimsum_results**", emit: results_dir
    path "versions.yml", emit: versions

  script:
  """
  set -euo pipefail

  # DiMSum expects the sequence string, not a file path
  WT=\$(tr -d ' \r\\n\\t' < "$wt_txt")

  DiMSum \
      --experimentDesignPath "$exp_design" \
      --wildtypeSequence "\$WT" \
      --countPath "$counts_merged" \
      --startStage 4 \
      --stopStage 5 \
      --fitnessErrorModel F \
      --retainIntermediateFiles T \
      --projectName "dimsum_results" \
      --fastqFileDir . \

  R_VERSION=\$(R --version | head -n 1 | sed -E 's/^R version ([0-9.]+).*/\\1/')
  cat <<-END_VERSIONS > versions.yml
  DIMSUM_RUN:
    r-base: \$R_VERSION
  END_VERSIONS
  """
}
