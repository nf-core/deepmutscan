
process PREMERGE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/bwa_samtools_vsearch:28e8640725d3d8e9"

    input:
    tuple val(meta), path(bam)
    path wt_seq

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Convert BAM to paired FASTQ files
    samtools fastq -1 forward_reads.fastq -2 reverse_reads.fastq -0 /dev/null -s /dev/null -n $bam

    # Merge paired reads
    vsearch --fastq_mergepairs forward_reads.fastq --reverse reverse_reads.fastq --fastqout merged_reads.fastq --fastq_minovlen 10 --fastq_allowmergestagger

    # Re-align merged reads
    bwa index $wt_seq
    bwa mem $wt_seq merged_reads.fastq | samtools view -Sb - > ${prefix}_merged.bam

    # Save version information
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        premerge: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch merged_reads.fastq
    touch merged_reads.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        premerge: dummy_version
    END_VERSIONS
    """
}
