process RUN_DIMSUM {
  tag { sample.sample }
  label 'process_single'

  conda "${moduleDir}/environment.yml"

  container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
      ? 'REPLACE_WITH_SINGULARITY_HTTPS_URL'
      : 'REPLACE_WITH_DOCKER_IMAGE'}"

  input:
    tuple val(sample), path(counts_merged)
    path(wt_txt)
    path(exp_design)
    val(min_counts)
    val(output_pseudocount)
    path(report_rmd)      // assets/dimsum/report.Rmd
    path(logo_png)        // assets/dimsum/Dumpling.png

  output:
    path "dimsum_results**", emit: results_dir
    path "versions.yml", emit: versions

  script:
  def args = task.ext.args ?: ''
  """
  set -euo pipefail

  # DiMSum expects the sequence string, not a file path
  WT=\$(tr -d ' \r\\n\\t' < "$wt_txt")

  # DiMSum's STEAM stages (4: process variants, 5: fitness) from bin/dimsum_fitness.R - DiMSum 1.4's own
  # code, run without the package; the same options as DiMSum's command line. Its error messages are
  # kept so that a library too sparse for the error model gets a pointer to --dimsum_error_model.
  status=0
  dimsum_fitness.R \\
      --experimentDesignPath "$exp_design" \\
      --wildtypeSequence "\$WT" \\
      --countPath "$counts_merged" \\
      --startStage 4 \\
      --stopStage 5 \\
      --fitnessMinInputCountAll ${min_counts} \\
      --fitnessDropoutPseudocount ${output_pseudocount} \\
      --numCores $task.cpus \\
      --projectName "dimsum_results" \\
      --fastqFileDir . \\
      --reportRmd "$report_rmd" \\
      --logoPng "$logo_png" \\
      $args \\
      2> dimsum_fitness.err || status=\$?
  cat dimsum_fitness.err >&2
  if [ \$status -ne 0 ]; then
      if grep -q "insufficent number of variants satisfying full fitness range" dimsum_fitness.err; then
          echo "DiMSum's error model needs at least 30 variants per replicate (90 for three replicates) that pass --min_counts and are observed in every replicate. Rerun with --dimsum_error_model false (count-based errors only)." >&2
      fi
      exit \$status
  fi

  R_VERSION=\$(R --version | head -n 1 | sed -E 's/^R version ([0-9.]+).*/\\1/')
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      r-base: \$R_VERSION
      DiMSum: 1.4 (stages 4-5, bin/dimsum_fitness.R)
END_VERSIONS
  """

  stub:
  """
  mkdir -p dimsum_results
  touch dimsum_results/report.html
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      r-base: "0.0.0"
      DiMSum: "0.0.0"
  END_VERSIONS
  """
}
