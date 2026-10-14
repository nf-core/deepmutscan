
process PREMERGE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/bwa_samtools_vsearch:28e8640725d3d8e9"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*_merged.bam"), emit: bam
    tuple val("${task.process}"), val('bwa'), eval('bwa 2>&1 | sed -n "s/^Version: //p"'), topic: versions, emit: versions_bwa
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools
    tuple val("${task.process}"), val('vsearch'), eval("vsearch --version 2>&1 | head -n1 | sed 's/vsearch v//; s/_.*//'"), topic: versions, emit: versions_vsearch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    # Convert the filtered paired alignments back to FASTQ. Written gzipped: these are pure
    # intermediates, and at full scale (~80M read pairs per library) each one is ~27 GB
    # uncompressed, which is what made the work directory grow to ~600 GB on the GID1A test_full.
    samtools fastq -@ $task.cpus -1 forward_reads.fastq.gz -2 reverse_reads.fastq.gz -0 /dev/null -s /dev/null -n $bam

    # Merge overlapping mates into one consensus read, so the overlap is counted once and its
    # bases carry a recalculated quality, then re-align against the same index BWA_MEM used and
    # re-apply the filter BAMFILTER_DMS applies to the paired alignments. Without that filter the
    # counter sees the re-aligned BAM unfiltered: supplementary alignments of chimeric merged
    # reads get counted twice, low-MAPQ hits get counted at all, and a merged read carrying a
    # deletion is reported as whatever substitution it also happens to carry (gatk_formating
    # keeps only X events).
    #
    # The merged reads are streamed into bwa rather than written out: vsearch cannot gzip its
    # output (given a .gz name it writes plain text under it), and the merged FASTQ is another
    # ~27 GB per library. Streaming keeps it off disk entirely. Verified to produce a
    # byte-identical alignment set to the file-based form.
    vsearch --fastq_mergepairs forward_reads.fastq.gz \\
        --reverse reverse_reads.fastq.gz \\
        --fastqout - \\
        --threads $task.cpus \\
        $args \\
        | bwa mem -t $task.cpus \$INDEX - \\
        | samtools view -h $args2 - \\
        | awk '/^@/ || \$6 !~ /[IDN]/' \\
        | samtools view -b --threads $task.cpus -o ${prefix}_merged.bam -
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_merged.bam
    """
}
