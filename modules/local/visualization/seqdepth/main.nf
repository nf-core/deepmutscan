process VISUALIZATION_SEQDEPTH {
    tag "$meta.id"
    label 'process_high'
    label 'r_env'

    conda "${moduleDir}/environment.yml"

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
