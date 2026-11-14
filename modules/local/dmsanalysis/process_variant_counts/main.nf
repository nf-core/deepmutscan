process DMSANALYSIS_PROCESS_VARIANT_COUNTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/7988a2e369e98cec897d48bab3e55daddd8a973d9b7491ad27b6b15d550a4c5b/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-dplyr_r-ggplot2_pruned:32f1af0a83a018e7'}"

    input:
    tuple val(meta), path(variantCounts)
    path possible_mutations
    path aa_seq

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
