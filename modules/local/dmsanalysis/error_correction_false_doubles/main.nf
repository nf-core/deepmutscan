process DMSANALYSIS_ERROR_CORRECTION_FALSE_DOUBLES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e5/e5a2e5419dff47f9dbcea8823692f60f9a1b1c1a1ee32ae56eef0ce2989f7b72/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-stringr:801bac9c074a1f9b'}"

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
