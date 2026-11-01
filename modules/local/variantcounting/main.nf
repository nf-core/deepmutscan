process VARIANTCOUNTING {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'REPLACE_WITH_SINGULARITY_HTTPS_URL'
        : 'REPLACE_WITH_DOCKER_IMAGE'}"

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
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # One pass over the sorted BAM: per-base coverage of the reading frame (quality and read-edge filters),
    # per-read nucleotide variants from CIGAR + MD, codon and amino-acid annotation, aggregation. The table
    # is byte-identical to the one the Python counter of v1.0.0 produced.
    dms_count \\
        --bam ${bam} \\
        --fasta ${fasta} \\
        --reading_frame ${reading_frame} \\
        --base_qual ${base_qual} \\
        --min_flank ${min_flank} \\
        --threads ${task.cpus} \\
        --id ${meta.id} \\
        --columns variant_counts_columns.tsv \\
        --out ${prefix}.variant_counts.tsv \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dms_count: \$(dms_count --version | sed 's/^dms_count //; s/ (.*//')
        htslib: \$(dms_count --version | sed 's/.*htslib //; s/)//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.variant_counts.tsv
    touch variant_counts_columns.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dms_count: stub
        htslib: stub
    END_VERSIONS
    """
}
