process DMSANALYSIS_PROCESS_GATK {
    tag "$meta.id"
    label 'process_single'
    label 'r_env'

    conda "${moduleDir}/environment.yml"

    publishDir "${params.outdir}/intermediate_files", mode: 'copy'

    input:
    tuple val(meta), path(variantCounts)
    path possible_mutations
    path aa_seq
    val min_counts

    output:
    tuple val(meta),
        path("annotated_variantCounts.csv"),
        path("variantCounts_filtered_by_library.csv"),
        path("library_completed_variantCounts.csv"),
        path("variantCounts_for_heatmaps.csv"),
        emit: processed_variantCounts

    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'process_gatk.R'

    stub:
    """
    touch annotated_variantCounts.csv variantCounts_filtered_by_library.csv library_completed_variantCounts.csv variantCounts_for_heatmaps.csv
    echo "DMSANALYSIS_PROCESSGATK:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
    """
}
