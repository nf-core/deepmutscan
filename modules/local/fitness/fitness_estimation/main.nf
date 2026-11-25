process FITNESS_ESTIMATION {
  tag "${sample.sample}"
  label 'process_single'

  conda "${moduleDir}/environment.yml"

  container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
      ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e2/e2d657c330e5ac4dcc3b40ef513d377367e60c293e01003906877e40f6865d75/data'
      : 'community.wave.seqera.io/library/bioconductor-edger_bioconductor-impute_bioconductor-limma_dms-fdmcorrect_pruned:be9da4b29f56a3d0'}"

  input:
    tuple val(sample), val(library_ids), val(library_types), path(counts, stageAs: 'counts_??/*')
                                                                 // the libraries of one biological sample (ids and types,
                                                                 // space-separated: inputs and outputs by replicate, then any
                                                                 // others) and their variant-count tables, in that order
    path samplesheet_csv
    path wt_fasta
    val pos_range                                                // "start-end", 1-based inclusive
    path wt_seq                                                  // wildtype amino-acid sequence (aa_seq.txt)
    val min_counts
    val output_pseudocount

  output:
    tuple val(sample), path("*_fitness_input.tsv"),                        emit: fitness_input
    tuple val(sample), path("counts_merged.tsv"),                          emit: merged_counts
    tuple val(sample), path("experimentalDesign.tsv"),                     emit: experimental_design
    tuple val(sample), path("synonymous_wt.txt"),                          emit: synonymous_wt
    tuple val(sample), path("fitness_estimation.tsv"),                     emit: fitness_estimation
    tuple val(sample), path("fitness_estimation_count_correlation.pdf"),   emit: counts_corr_pdf
    tuple val(sample), path("fitness_estimation_fitness_correlation.pdf"), emit: fitness_corr_pdf
    tuple val(sample), path("fitness_heatmap.pdf"),                        emit: fitness_heatmap
    path "versions.yml", emit: versions

  script:
  template 'fitness_estimation.R'

  stub:
  def stub_inputs = library_ids.tokenize(' ').collect { id -> "${id}_fitness_input.tsv" }.join(' ')
  """
  touch ${stub_inputs}
  touch counts_merged.tsv experimentalDesign.tsv synonymous_wt.txt fitness_estimation.tsv \\
        fitness_estimation_count_correlation.pdf fitness_estimation_fitness_correlation.pdf fitness_heatmap.pdf
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      r-base: "0.0.0"
      r-ggplot2: "0.0.0"
  END_VERSIONS
  """
}
