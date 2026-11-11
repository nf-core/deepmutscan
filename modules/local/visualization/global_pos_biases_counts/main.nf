process VISUALIZATION_GLOBAL_POS_BIASES_COUNTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/7988a2e369e98cec897d48bab3e55daddd8a973d9b7491ad27b6b15d550a4c5b/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-dplyr_r-ggplot2_pruned:32f1af0a83a018e7'}"

    input:
    tuple val(meta), path(variantCounts_filtered_by_library)
    path aa_seq
    val sliding_window_size

    output:
    tuple val(meta), path("rolling_counts.pdf"), emit: rolling_counts
    tuple val(meta), path("rolling_counts_per_cov.pdf"), emit: rolling_counts_per_cov
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'global_position_biases_counts_and_counts_per_cov.R'

    stub:
    """
    touch rolling_counts.pdf
    touch rolling_counts_per_cov.pdf
    echo "VISUALIZATION_GLOBAL_POS_BIASES_COUNTS:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
    """
}
