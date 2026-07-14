process VISUALIZATION_GLOBAL_POS_BIASES_COV {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a6/a66af652f7b1c6d9dab9f081f5676bd9452d653eed43034f20a2cf172921a4cf/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-mutscan_pysam_biopython_pruned:fb9a2095922ddd59' }"

    input:
    tuple val(meta), path(variantCounts_filtered_by_library)
    path aa_seq
    val sliding_window_size
    val aimed_cov

    output:
    tuple val(meta), path("rolling_coverage.pdf"), emit: rolling_coverage
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'global_position_biases_cov.R'

    stub:
    """
    touch rolling_coverage.pdf
    echo "VISUALIZATION_COUNTS_HEATMAP:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
    """
}
