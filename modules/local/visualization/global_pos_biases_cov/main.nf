process VISUALIZATION_GLOBAL_POS_BIASES_COV {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/78/7815f74d3680b875aa8c7be8f78dc77921ec970b99b45a539c5d2225749c0b90/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-dplyr_r-ggplot2_pruned:a65d6fd383d85714' }"

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
