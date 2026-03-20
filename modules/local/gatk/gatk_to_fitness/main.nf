process GATK_GATKTOFITNESS {
    tag "$meta.id"
    label 'process_single'
    label 'r_env'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(variantCounts_filtered_by_library)
    path wt_seq
    val pos_range

    output:
    tuple val(meta), path("${meta.id}_fitness_input.tsv"), emit: fitness_input
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'gatk_to_fitness.R'

    stub:
    """
    touch ${meta.id}_fitness_input.tsv
    echo "GATK_GATKTOFITNESS:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
    """
}
