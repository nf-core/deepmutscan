process VISUALIZATION_SEQDEPTH {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity'
        ? 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:ce2ba7ad7f6e7f2c'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:0fd2e39a5bf2ecaa' }"

    input:
    tuple val(meta), path(variantCounts_filtered_by_library)
    path possible_mutations
    val min_counts

    output:
    tuple val(meta), path("SeqDepth.pdf"), emit: SeqDepth
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'SeqDepth_simulation.R'

    stub:
    """
    touch SeqDepth.pdf
    """
}
