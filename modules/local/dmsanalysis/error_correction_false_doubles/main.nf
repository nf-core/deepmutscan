process DMSANALYSIS_ERROR_CORRECTION_FALSE_DOUBLES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bc/bcd9a556c2ee2e9ce6d21749b6a15954d2578c6d4738e1ab72ea5de148ab2ab9/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-stringr:b64d48973b58dc0c' }"

    input:
    tuple val(meta), path(raw_counts), path(filtered), path(completed)
    path wt_fasta
    path aa_seq
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
