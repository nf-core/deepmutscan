process VISUALIZATION_COUNTS_PER_COV {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a6/a66af652f7b1c6d9dab9f081f5676bd9452d653eed43034f20a2cf172921a4cf/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-mutscan_pysam_biopython_pruned:fb9a2095922ddd59' }"

    input:
    tuple val(meta), path(variantCounts_for_heatmaps)
    val min_counts

    output:
    tuple val(meta), path("counts_per_cov_heatmap.pdf"), emit: counts_per_cov_heatmap
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'counts_per_cov_heatmap.R'

    stub:
    """
    touch counts_per_cov_heatmap.pdf
    echo "VISUALIZATION_COUNTS_PER_COV:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
    """
}
