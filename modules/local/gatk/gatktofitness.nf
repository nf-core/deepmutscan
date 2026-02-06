process GATK_GATKTOFITNESS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity'
        ? 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:ce2ba7ad7f6e7f2c'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:0fd2e39a5bf2ecaa' }"

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
