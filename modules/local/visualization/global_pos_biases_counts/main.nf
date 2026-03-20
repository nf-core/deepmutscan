process VISUALIZATION_GLOBAL_POS_BIASES_COUNTS {
    tag "$meta.id"
    label 'process_single'
    label 'r_env'

    conda "${moduleDir}/environment.yml"

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
    echo "VISUALIZATION_COUNTS_HEATMAP:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
    """
}
