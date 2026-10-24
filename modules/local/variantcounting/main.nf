process VARIANTCOUNTING {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/19/199686bb01fc18539d6ca35597ce574329a29ccae5f0b4bb1a6982f2bd3ce19a/data'
        : 'community.wave.seqera.io/library/pysam_biopython_numpy_polars-lts-cpu_pruned:ac0353fd8d903112' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    val reading_frame
    val base_qual
    val min_flank

    output:
    tuple val(meta), path("*.variant_counts.tsv")          , emit: variant_counts
    path "variant_counts_columns.tsv"                      , emit: columns
    path "versions.yml"                                    , emit: versions

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
