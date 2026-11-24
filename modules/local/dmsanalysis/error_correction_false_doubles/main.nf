// Applies the FDM error rates of a library's pool (DMSANALYSIS_FDM_ESTIMATE) to that library's counts.
process DMSANALYSIS_ERROR_CORRECTION_FALSE_DOUBLES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e2/e2d657c330e5ac4dcc3b40ef513d377367e60c293e01003906877e40f6865d75/data'
        : 'community.wave.seqera.io/library/bioconductor-edger_bioconductor-impute_bioconductor-limma_dms-fdmcorrect_pruned:be9da4b29f56a3d0'}"

    input:
    tuple val(meta), path(filtered), path(completed), path(evidence), path(rates)
    path wt_fasta
    path aa_seq

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
