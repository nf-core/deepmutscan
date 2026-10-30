process VARIANTCOUNTING {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/62/6247bcc77159c0adea94e0fdc439006f94d499e9d32004470fa7e402ede08538/data'
        : 'community.wave.seqera.io/library/pysam_biopython_numpy_polars-lts-cpu_pruned:f223db6b9aa476d1'}"

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
