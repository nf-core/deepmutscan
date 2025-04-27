// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process PREMERGE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/bwa_samtools_vsearch:28e8640725d3d8e9"

    input:
    tuple val(meta), path(bam)
    path wt_seq

    output:
    tuple val(meta), path("merged_reads.bam"), emit: bam
    tuple val(meta), path("merged_reads.fastq"), emit: fastq
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
    bwa mem $wt_seq merged_reads.fastq | samtools view -Sb - > merged_reads.bam

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
