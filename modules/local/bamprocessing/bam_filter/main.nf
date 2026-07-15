
process BAMFILTER_DMS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h96c455f_1':
        'biocontainers/samtools:1.21--h96c455f_1' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Keep primary, mapped, MAPQ>=30 reads and drop indel/splice alignments.
    # Wildtype (NM:i:0) reads are intentionally retained (needed for downstream
    # sequencing-error correction); the variant counter ignores them for calling.
    samtools view -h -F 4 -F 256 -q 30 $bam | \
    samtools view -h | \
    awk '{if(\$6 !~ /I/ && \$6 !~ /D/ && \$6 !~ /N/) print \$0}' | \
    samtools view -bS > ${meta.id}_filtered.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version |& sed '1!d ; s/samtools //')
END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamfilteringdms: \$(samtools --version |& sed '1!d ; s/samtools //')
END_VERSIONS
    """
}
