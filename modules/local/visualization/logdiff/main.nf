process VISUALIZATION_LOGDIFF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity'
        ? 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:ce2ba7ad7f6e7f2c'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:0fd2e39a5bf2ecaa' }"

    input:
    tuple val(meta), path(library_completed_variantCounts)

    output:
    tuple val(meta), path("logdiff_plot.pdf"), emit: logdiff_plot
    tuple val(meta), path("logdiff_varying_bases.pdf"), emit: logdiff_varying_bases
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'logdiff.R'

    stub:
    """
    touch logdiff_plot.pdf
    touch logdiff_varying_bases.pdf
    echo "VISUALIZATION_COUNTS_HEATMAP:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
    """
}
