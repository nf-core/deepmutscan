process VARIANTCOUNTING {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a6/a66af652f7b1c6d9dab9f081f5676bd9452d653eed43034f20a2cf172921a4cf/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-mutscan_pysam_biopython_pruned:fb9a2095922ddd59' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    val reading_frame
    val min_counts
    val base_qual
    val min_flank

    output:
    tuple val(meta), path("*.variant_counts.tsv"), emit: variant_counts
    path "variant_counts_columns.tsv"            , emit: columns
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'count_variants.py'

    stub:
    """
    touch ${meta.id}.variant_counts.tsv
    touch variant_counts_columns.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        variantcounting: stub
    END_VERSIONS
    """
}
