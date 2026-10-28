
process MERGE_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/85/856173edaf176a724c7e5fe23987e42e5ab0bc476d20e458bf2a68d738e74293/data'
        : 'community.wave.seqera.io/library/bwa_samtools_vsearch:593a959d25ab8798' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: bam
    tuple val(meta), path("*.vsearch_mergepairs.log"), emit: log
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
    # The BWA index prefix, resolved with a glob and bash parameter expansion instead of `find`.
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
        echo "No BWA index (*.amb) found in './${index}'" >&2
        exit 1
    fi

    # Merge the overlapping mates of each raw read pair into one consensus read, so the overlap is
    # counted once and its bases carry a recalculated quality. Merging comes first, before any
    # alignment: the reads are aligned once (as merged reads), and the filter below acts on exactly
    # the alignments the counter reads. Pairs that cannot be merged are discarded.
    #
    # The merged reads are streamed into bwa mem rather than written out (vsearch cannot gzip its
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
        | bwa mem -t $task.cpus \$INDEX - \\
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
