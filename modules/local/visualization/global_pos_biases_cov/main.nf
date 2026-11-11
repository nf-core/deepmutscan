process VISUALIZATION_GLOBAL_POS_BIASES_COV {
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
    val aimed_cov

    output:
    tuple val(meta), path("rolling_coverage.pdf"), emit: rolling_coverage
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'global_position_biases_cov.R'

    stub:
    """
    touch rolling_coverage.pdf
    echo "VISUALIZATION_GLOBAL_POS_BIASES_COV:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
    """
}
