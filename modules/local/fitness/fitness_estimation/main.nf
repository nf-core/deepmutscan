process FITNESS_ESTIMATION {
  tag "${sample.sample}"
  label 'process_single'

  conda "${moduleDir}/environment.yml"

  container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
      ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/7988a2e369e98cec897d48bab3e55daddd8a973d9b7491ad27b6b15d550a4c5b/data'
      : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-dplyr_r-ggplot2_pruned:32f1af0a83a018e7'}"

  input:
    tuple val(sample), path(input_counts), path(output_counts)   // per-library DiMSum-format counts of one biological sample
    path samplesheet_csv
    path wt_fasta
    val pos_range                                                // "start-end", 1-based inclusive
    path wt_seq                                                  // wildtype amino-acid sequence (aa_seq.txt)
    val min_counts
    val output_pseudocount

  output:
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
  """
  touch counts_merged.tsv experimentalDesign.tsv synonymous_wt.txt fitness_estimation.tsv \\
        fitness_estimation_count_correlation.pdf fitness_estimation_fitness_correlation.pdf fitness_heatmap.pdf
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      r-base: "0.0.0"
      r-Biostrings: "0.0.0"
      r-dplyr: "0.0.0"
      r-ggplot2: "0.0.0"
  END_VERSIONS
  """
}
