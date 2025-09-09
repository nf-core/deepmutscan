process MERGE_COUNTS {
  tag "${sample.sample}"
  label 'process_single'
  publishDir "${params.outdir}/dimsum/${sample.sample}", mode: 'copy'

  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity' 
      ? 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:ce2ba7ad7f6e7f2c' 
      : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:0fd2e39a5bf2ecaa' }"

  input:
    tuple val(sample), val(metas), path(input_counts), path(output_counts)
    path merge_script

  output:
    path "counts_merged.tsv", emit: merged_counts
    path "versions.yml", emit: versions

  script:
  def in_list  = input_counts .collect { it.getName() }.join(' ')
  def out_list = output_counts.collect { it.getName() }.join(' ')
  """
  set -euo pipefail

  Rscript "${merge_script}" \\
      --inputs  ${in_list} \\
      --outputs ${out_list} \\
      --out counts_merged.tsv

  R_VERSION=\$(R --version | head -n 1 | sed -E 's/^R version ([0-9.]+).*/\\1/')
    cat <<-END_VERSIONS > versions.yml
    GATK_GATKTODIMSUM:
      r-base: \$R_VERSION
    END_VERSIONS
  """
}
