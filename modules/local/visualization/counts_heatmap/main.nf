process VISUALIZATION_COUNTS_HEATMAP {
    tag "$meta.id"
    label 'process_single'
    label 'r_env'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(variantCounts_for_heatmaps)
    val min_counts

    output:
    tuple val(meta), path("counts_heatmap.pdf"), emit: counts_heatmap
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'counts_heatmap.R'

    stub:
    """
    touch counts_heatmap.pdf
    echo "VISUALIZATION_COUNTS_HEATMAP:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
    """
}
