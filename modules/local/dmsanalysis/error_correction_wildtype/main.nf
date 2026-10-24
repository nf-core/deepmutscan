process DMSANALYSIS_ERROR_CORRECTION_WILDTYPE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/78/7815f74d3680b875aa8c7be8f78dc77921ec970b99b45a539c5d2225749c0b90/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-dplyr_r-ggplot2_pruned:a65d6fd383d85714' }"

    input:
    tuple val(meta), path(filtered, stageAs: 'input_counts.csv'), path(wt_counts, stageAs: 'wildtype_counts.csv'), path(completed)
    path aa_seq

    output:
    tuple val(meta),
        path("variantCounts_filtered_by_library_error_corrected.csv"),
        path("variantCounts_for_heatmaps_error_corrected.csv"),
        path("library_completed_variantCounts_error_corrected.csv"),
        emit: corrected
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'error_correction_wildtype.R'

    stub:
    """
    touch variantCounts_filtered_by_library_error_corrected.csv
    touch variantCounts_for_heatmaps_error_corrected.csv
    touch library_completed_variantCounts_error_corrected.csv
    echo "DMSANALYSIS_ERROR_CORRECTION_WILDTYPE:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
    """
}
