process RUN_DIMSUM {
  tag { sample.sample }
  label 'process_single'

  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity'
      ? 'oras://community.wave.seqera.io/library/r-dimsum:1.4--4357734d345c8ccc'
      : 'docker.io/bwehnert1008/dms_qc_dimsum_environment@sha256:08f3bd8441df7b4a7e05aadeca178862153cf723e64097a48a2744b2698b15dd' }"

  input:
    tuple val(sample), path(counts_merged)
    path(wt_txt)
    path(exp_design)

  output:
    path "reports/**", emit: reports_folder
    path "dimsum_results**", emit: results_dir
    path "report.html", emit: report
    path "tmp/**", emit: tmp_folder
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
