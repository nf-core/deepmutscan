process DMSANALYSIS_PROCESS_VARIANT_COUNTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/73/73a72ec77725aeb67678a74228938fdd6827b669d01a8c96951b1a8ef96eeb0f/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-mutscan_r-base_r-biocmanager_pruned:c65036d76406f342' }"

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
    template 'process_variant_counts.R'

    stub:
    """
    touch annotated_variantCounts.csv variantCounts_filtered_by_library.csv library_completed_variantCounts.csv variantCounts_for_heatmaps.csv
    echo "DMSANALYSIS_PROCESS_VARIANT_COUNTS:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
    """
}
