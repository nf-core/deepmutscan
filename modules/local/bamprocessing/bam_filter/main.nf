
process BAMFILTER_DMS {
    tag "$meta.id"
    // samtools | awk | samtools is a three-process pipe: measured ~199% cpu on the full-size
    // GID1A run, where it was the longest single task. process_single declared 1 cpu for it.
    label 'process_low'

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
    samtools view -h $args $bam | \
    awk '/^@/ || (\$6 !~ /I/ && \$6 !~ /D/ && \$6 !~ /N/)' | \
    samtools view -bS --threads $task.cpus -o ${prefix}_filtered.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version |& sed '1!d ; s/samtools //')
END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Must match the name the script block writes: the input bam is staged under
    # \${meta.id}.bam, so touching that name here would just re-touch the input and
    # leave the process with no output of its own.
    touch ${prefix}_filtered.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version |& sed '1!d ; s/samtools //')
END_VERSIONS
    """
}
