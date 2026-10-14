process DMSANALYSIS_ERROR_CORRECTION_FALSE_DOUBLES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/73/73a72ec77725aeb67678a74228938fdd6827b669d01a8c96951b1a8ef96eeb0f/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-mutscan_r-base_r-biocmanager_pruned:c65036d76406f342' }"

    input:
    tuple val(meta), path(raw_counts), path(filtered), path(completed)
    path wt_fasta
    path aa_seq
    val min_counts
    val method
    val codon_window

    output:
    tuple val(meta),
        path("variantCounts_filtered_by_library_error_corrected.csv"),
        path("variantCounts_for_heatmaps_error_corrected.csv"),
        path("library_completed_variantCounts_error_corrected.csv"),
        emit: corrected
    tuple val(meta), path("seq_error_rate.csv"), emit: seq_error_rate
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'error_correction_false_doubles.R'

    stub:
    """
    touch variantCounts_filtered_by_library_error_corrected.csv
    touch variantCounts_for_heatmaps_error_corrected.csv
    touch library_completed_variantCounts_error_corrected.csv
    touch seq_error_rate.csv
    echo "DMSANALYSIS_ERROR_CORRECTION_FALSE_DOUBLES:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
    """
}
