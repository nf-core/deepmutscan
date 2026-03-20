process VISUALIZATION_COUNTS_PER_COV {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity'
        ? 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:ce2ba7ad7f6e7f2c'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:0fd2e39a5bf2ecaa' }"

    input:
    tuple val(meta), path(variantCounts_for_heatmaps)
    val min_counts

    output:
    tuple val(meta), path("counts_per_cov_heatmap.pdf"), emit: counts_per_cov_heatmap
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'counts_per_cov_heatmap.R'

    stub:
    """
    touch counts_per_cov_heatmap.pdf
    echo "VISUALIZATION_COUNTS_PER_COV:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
    """
}
