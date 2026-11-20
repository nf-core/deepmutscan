
process MERGE_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/40/4016ac329eacb2871c416a0324088d7905f22804c99e97c6badc6fe0d432bf5f/data'
        : 'community.wave.seqera.io/library/bwa-mem3_samtools_vsearch:a4bfb8df646b8637'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: bam
    tuple val(meta), path("*.vsearch_mergepairs.log"), emit: log
    tuple val("${task.process}"), val('bwa-mem3'), eval("bwa-mem3 version | sed -nE '1 s/^([0-9]+(\\.[0-9]+)+).*/\\1/p'"), topic: versions, emit: versions_bwamem3
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools
    tuple val("${task.process}"), val('vsearch'), eval("vsearch --version 2>&1 | head -n1 | sed 's/vsearch v//; s/_.*//'"), topic: versions, emit: versions_vsearch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
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

    # Merge the overlapping mates of each raw read pair into one consensus read, so the overlap is
    # counted once and its bases carry a recalculated quality. Merging comes first, before any
    # alignment: the reads are aligned once (as merged reads), and the filter below acts on exactly
    # the alignments the counter reads. Pairs that cannot be merged are discarded.
    #
    # bwa-mem3 mem (a maintained bwa-mem2 derivative; alignments identical to bwa 0.7.19 for the same
    # scoring, 3-4x faster on a gene-sized reference, ~7x with --fast) runs with the options in ext.args3
    # (default: --fast, no soft-clipping, gap open 12, mismatch 3). With bwa's default scoring a local
    # alignment clips a read end rather than paying for two or
    # three adjacent mismatches, so a 2-base codon change within 3 bp, or a 3-base change within 7 bp,
    # of a read end disappears from the alignment (and from the coverage): full-codon substitutions
    # were under-counted by ~13% relative to single-base changes on 150 bp merged reads. A clipping
    # penalty of 100 makes the alignment end-to-end for substitutions; the gap-open penalty is raised
    # from 6 to 12 so that three mismatches are not realigned as an indel (those reads are discarded
    # below), and the mismatch penalty lowered from 4 to 3 for the same reason.
    #
    # The merged reads are streamed into bwa-mem3 rather than written out (vsearch cannot gzip its
    # output, and the merged FASTQ is ~27 GB per full-size library), and the alignments are filtered
    # on the fly: no unmapped, secondary or supplementary records and MAPQ >= 30 (ext.args2), and no
    # insertions, deletions or splice gaps. Without that filter, supplementary alignments of
    # chimeric merged reads get counted twice, low-MAPQ hits get counted at all, and a merged read
    # carrying a deletion is reported as whatever substitution it also happens to carry
    # (gatk_formating keeps only X events). Wildtype reads are kept: the error correction needs them.
    # The filtered alignments are coordinate-sorted and indexed in the same pipe: an unsorted BAM is
    # about 3.5 times larger, and writing it out for a separate sort and index task only costs time.
    vsearch --fastq_mergepairs ${reads[0]} \\
        --reverse ${reads[1]} \\
        --fastqout - \\
        --threads $task.cpus \\
        $args \\
        2> ${prefix}.vsearch_mergepairs.log \\
        | bwa-mem3 mem -t $task.cpus $args3 \$INDEX - \\
        | samtools view -h $args2 - \\
        | awk '/^@/ || \$6 !~ /[IDN]/' \\
        | samtools sort --threads $task.cpus -o ${prefix}.sorted.bam##idx##${prefix}.sorted.bam.bai --write-index -

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
    touch ${prefix}.sorted.bam
    touch ${prefix}.sorted.bam.bai
    touch ${prefix}.vsearch_mergepairs.log
    """
}
