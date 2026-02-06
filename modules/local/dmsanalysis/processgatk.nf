process DMSANALYSIS_PROCESS_GATK {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity'
        ? 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:ce2ba7ad7f6e7f2c'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:0fd2e39a5bf2ecaa' }"

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
