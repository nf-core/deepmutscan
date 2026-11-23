process MERGE_ALIGN_COUNT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/67/6749b447ff36bedbae41e1cd651673e714ec69d8f67facd631d6fd373e098f44/data'
        : 'community.wave.seqera.io/library/bwa-mem3_dms-count_samtools_vsearch:0a05801b076ff393'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    path fasta
    val reading_frame
    val base_qual
    val min_flank
    val save_bam

    output:
    tuple val(meta), path("*.variant_counts.tsv")                       , emit: variant_counts
    path "variant_counts_columns.tsv"                                   , emit: columns
    tuple val(meta), path("*.vsearch_mergepairs.log")                   , emit: log
    tuple val(meta), path("*.sorted.bam"), path("*.sorted.bam.bai")     , emit: bam, optional: true
    tuple val("${task.process}"), val('vsearch'), eval("vsearch --version 2>&1 | head -n1 | sed 's/vsearch v//; s/_.*//'"), topic: versions, emit: versions_vsearch
    tuple val("${task.process}"), val('bwa-mem3'), eval("bwa-mem3 version | sed -nE '1 s/^([0-9]+(\\.[0-9]+)+).*/\\1/p'"), topic: versions, emit: versions_bwamem3
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools
    tuple val("${task.process}"), val('dms_count'), eval("dms_count --version | sed 's/^dms_count //; s/ (.*//'"), topic: versions, emit: versions_dms_count
    tuple val("${task.process}"), val('htslib'), eval("dms_count --version | sed 's/.*htslib //; s/)//'"), topic: versions, emit: versions_htslib

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args  ?: ''
    def args2  = task.ext.args2 ?: ''
    def args3  = task.ext.args3 ?: ''
    def args4  = task.ext.args4 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // With --save_bam the filtered alignments are also copied off the stream (tee) into a background
    // samtools sort; the counter does not need them sorted, so without it nothing is sorted or written.
    def bam_start = save_bam ? """
    mkfifo alignments.fifo
    samtools sort -l 1 --threads ${task.cpus} -o ${prefix}.sorted.bam##idx##${prefix}.sorted.bam.bai --write-index - < alignments.fifo &
    sort_pid=\$!
    """ : ''
    def bam_tee  = save_bam ? '| tee alignments.fifo \\' : '\\'
    def bam_wait = save_bam ? 'wait \$sort_pid' : ''
    """
    # The bwa-mem3 index prefix (same files as a bwa-mem2 index: .0123 .amb .ann .bwt.2bit.64 .pac),
    # resolved with a glob and bash parameter expansion instead of `find`.
    # Seqera Containers builds the Docker image of a conda spec on an Ubuntu base but the
    # Singularity image on an Amazon Linux 2023 one, and that base ships no findutils, so the usual
    # `find -L ./ -name "*.amb"` fails with exit 127 under -profile singularity. Globbing needs no
    # external binary. The leading ./ is kept so the prefix, which bwa records in the BAM's @PG
    # header, stays exactly what find returned.
    INDEX=''
    for amb in ./${index}/*.amb; do
        if [ -e "\$amb" ]; then INDEX="\${amb%.amb}"; fi
    done
    if [ -z "\$INDEX" ]; then
        echo "No bwa-mem3 index (*.amb) found in './${index}'" >&2
        exit 1
    fi

    # One stream per library: merge -> align -> filter -> count.
    #
    # Merge the overlapping mates of each raw read pair into one consensus read, so the overlap is
    # counted once and its bases carry a recalculated quality. Merging comes first, before any
    # alignment: the reads are aligned once (as merged reads), and the filter below acts on exactly
    # the alignments the counter reads. Pairs that cannot be merged are discarded. From 2.32.0 on,
    # vsearch looks up the seeding 5-mers of --fastq_mergepairs in a directly addressed index rather
    # than a hash table: the same merged reads and statistics as 2.30.0 at half the CPU time.
    # --no_progress: since 2.31.0 vsearch writes one progress line per percent when stderr is not a
    # terminal (as here); without them the log is the statistics block 2.30.0 wrote.
    #
    # bwa-mem3 mem (a maintained bwa-mem2 derivative; alignments identical to bwa 0.7.19 for the same
    # scoring, 3-4x faster on a gene-sized reference, ~7x with --fast) runs with the options in ext.args3
    # (default: --fast, no soft-clipping, gap open 12, mismatch 3). With bwa's default scoring a local
    # alignment clips a read end rather than paying for two or three adjacent mismatches, so a 2-base
    # codon change within 3 bp, or a 3-base change within 7 bp, of a read end disappears from the
    # alignment (and from the coverage): full-codon substitutions were under-counted by ~13% relative
    # to single-base changes on 150 bp merged reads. A clipping penalty of 100 makes the alignment
    # end-to-end for substitutions; the gap-open penalty is raised from 6 to 12 so that three
    # mismatches are not realigned as an indel, and the mismatch penalty lowered from 4 to 3 for the
    # same reason.
    #
    # The alignments are filtered on the fly (ext.args2): no unmapped, secondary or supplementary
    # records and MAPQ >= 30. Without that filter, supplementary alignments of chimeric merged reads
    # get counted twice and low-MAPQ hits get counted at all. Wildtype reads are kept: the error
    # correction needs them. samtools view hands uncompressed BAM (-u) straight to the counter.
    #
    # dms_count reads the alignments once, in the order they arrive (it never needs them sorted or
    # indexed): per-base coverage of the reading frame (quality and read-edge filters), per-read
    # nucleotide variants from CIGAR + MD, codon and amino-acid annotation, aggregation. Reads whose
    # alignment contains an insertion, deletion or reference skip are ignored (--skip_indel_reads):
    # a merged read carrying a deletion would otherwise be reported as whatever substitution it also
    # happens to carry. Counting from the stream gives the same table as counting from the sorted BAM
    # and saves the sort, the BAM (~2 GB per full-size library) and a task per library.
    #
    # Threads: vsearch and bwa-mem3 each get task.cpus and run concurrently, so the two compute stages
    # share the task's CPUs; samtools view and dms_count need well under a CPU together. With a
    # single CPU the stream still works, but bwa-mem3 then reads and aligns in turns and the
    # stages run one after the other.
    ${bam_start}
    vsearch --fastq_mergepairs ${reads[0]} \\
        --reverse ${reads[1]} \\
        --fastqout - \\
        --threads $task.cpus \\
        --no_progress \\
        $args \\
        2> ${prefix}.vsearch_mergepairs.log \\
        | bwa-mem3 mem -t $task.cpus $args3 \$INDEX - \\
        | samtools view -u $args2 - \\
        ${bam_tee}
        | dms_count \\
            --bam - \\
            --fasta ${fasta} \\
            --reading_frame ${reading_frame} \\
            --base_qual ${base_qual} \\
            --min_flank ${min_flank} \\
            --skip_indel_reads \\
            --threads 1 \\
            --id ${meta.id} \\
            --columns variant_counts_columns.tsv \\
            --out ${prefix}.variant_counts.tsv \\
            $args4
    ${bam_wait}

    # vsearch's first log line reports the host's RAM and core count, so the file differs
    # between machines and the nf-test snapshot can never match across them. Drop that line and
    # keep the statistics, which are reproducible; the version is reported through the versions
    # topic anyway.
    sed '1d' ${prefix}.vsearch_mergepairs.log > mergepairs.stats && mv mergepairs.stats ${prefix}.vsearch_mergepairs.log

    # vsearch's merging statistics, also in the task log
    cat ${prefix}.vsearch_mergepairs.log >&2
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.variant_counts.tsv
    touch variant_counts_columns.tsv
    touch ${prefix}.vsearch_mergepairs.log
    ${save_bam ? "touch ${prefix}.sorted.bam ${prefix}.sorted.bam.bai" : ''}
    """
}
