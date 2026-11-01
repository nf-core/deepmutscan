process RAW_DATA_QC {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7a/7aca4e9958243151b7f31fe9df0639e20bc73d97c3ae22aa06a44d23c135e6ea/data'
        : 'community.wave.seqera.io/library/falco_fq:db74719b1b1ce638'}"

    input:
    tuple val(meta), path(reads, arity: '2')

    output:
    tuple val(meta), path("*.fastq_check.tsv")                     , emit: report
    tuple val(meta), path("*_fastqc_data.txt")                     , emit: data, optional: true
    tuple val(meta), path("*_falco_report.html")                   , emit: html, optional: true
    tuple val(meta), path("*_falco_summary.txt")                   , emit: summary, optional: true
    tuple val("${task.process}"), val('fq'), eval("fq lint --version | sed 's/fq-lint //; s/ .*//'"), emit: versions_fq, topic: versions
    tuple val("${task.process}"), val('falco'), eval("falco --version 2>&1 | head -n1 | sed 's/^[^0-9]*//'"), emit: versions_falco, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args  ?: ''
    def args2  = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def r1     = reads[0]
    def r2     = reads[1]
    // Falco reports carry the file name; the reads are linked as <prefix>_1/_2 (as nf-core's FastQC
    // module renames them), so MultiQC shows one sample per library and mate
    def ext1   = r1.name.endsWith('.gz') ? '.fastq.gz' : '.fastq'
    def ext2   = r2.name.endsWith('.gz') ? '.fastq.gz' : '.fastq'
    """
    # Raw data QC of one paired-end library in one task: fq lint (completeness and consistency of the
    # pair; the run stops if any library fails) and Falco (the FastQC report, in C++).
    #
    # fq lint reads both files completely and checks gzip integrity, the format of every FASTQ
    # record, equal read numbers, and matching read names in the same order. Empty files, and the
    # same file given twice, are caught before this, at start-up (validateInputFiles). The task never
    # fails: it writes PASS or FAIL with a plain-language reason to ${prefix}.fastq_check.tsv, and the
    # workflow stops with one message that lists every library that failed.
    status=0
    fq lint ${args} ${r1} ${r2} > /dev/null 2> fq_lint.log || status=\$?
    sed 's/\\x1b\\[[0-9;]*m//g' fq_lint.log > fq_lint.clean.log
    cat fq_lint.clean.log >&2

    # fq's own message: the first lines that are not log lines or stack trace
    detail=\$(grep -v -E ' INFO |^ *[0-9]+: |^ +at |^Stack backtrace|^ *\$' fq_lint.clean.log | head -n 3 | tr '\\n' ' ' | sed 's/  */ /g; s/ *\$//' || true)

    if [ "\$status" -ne 0 ]; then
        if grep -q -E 'deflate|gzip|I/O error|unexpected end' fq_lint.clean.log; then
            # Find out which of the two files is broken
            broken=""
            for fastq in ${r1} ${r2}; do
                fq lint "\$fastq" > /dev/null 2>&1 || broken="\$broken \$fastq"
            done
            result="FAIL\tcorrupted, truncated or not gzip-compressed file:\${broken:- file1 or file2} (\$detail)"
        elif grep -q 'r1-src unexpectedly ended' fq_lint.clean.log; then
            result="FAIL\tfile1 (${r1}) has fewer reads than file2 (${r2}): one of them is incomplete, or they are not mates"
        elif grep -q 'r2-src unexpectedly ended' fq_lint.clean.log; then
            result="FAIL\tfile2 (${r2}) has fewer reads than file1 (${r1}): one of them is incomplete, or they are not mates"
        elif grep -q '\\[P001\\]' fq_lint.clean.log; then
            result="FAIL\tread names of file1 and file2 do not match (\$detail). vsearch pairs reads by position, so the mates must be in the same order. If both files contain the same reads in a different order (e.g. lanes concatenated in a different order), re-pair them with 'seqkit pair' (https://github.com/shenwei356/seqkit) or BBTools 'repair.sh' (https://github.com/bbushnell/BBTools); otherwise check that the two files belong to the same library"
        elif grep -q -E '\\[S[0-9]+\\]' fq_lint.clean.log; then
            result="FAIL\tmalformed FASTQ record (\$detail)"
        else
            result="FAIL\tfq lint failed (\$detail)"
        fi
    else
        result="PASS\tcomplete and consistent"
    fi
    printf '%b\\n' "\$result" > ${prefix}.fastq_check.tsv

    # Falco (de Sena Brandine & Smith 2019) computes FastQC's modules and writes FastQC's
    # fastqc_data.txt, which MultiQC reads as it reads FastQC's (its section is named Falco), at a
    # fraction of FastQC's cost: 1.4 instead of 14.5 CPU s per 1.2 M reads, 14 MB instead of 0.6 GB.
    # Only files that passed fq lint are analysed. The two mates run side by side, one thread each:
    # with more threads per file, the per-tile section of fastqc_data.txt varies from run to run.
    # --groups (ext.args2) bins positions after the first 9 bp, as FastQC does by default.
    #
    # Falco 2.0.3 grades two modules against cutoffs in the wrong units (smithlabcode/falco#145, fix in #146):
    # overrepresented sequences and adapter content can read PASS where FastQC's cutoffs give WARN or
    # FAIL. The values themselves, and the plots, are right.
    if [ "\$status" -eq 0 ]; then
        mkdir falco_input
        ln -s ../${r1} falco_input/${prefix}_1${ext1}
        ln -s ../${r2} falco_input/${prefix}_2${ext2}
        falco -t 1 $args2 -o falco_output falco_input/${prefix}_1${ext1} &
        falco_1=\$!
        falco -t 1 $args2 -o falco_output falco_input/${prefix}_2${ext2} &
        falco_2=\$!
        wait \$falco_1
        wait \$falco_2

        # Falco 2.0.3 writes its version line as "##Falco 2.0.3"; MultiQC expects a tab there to show
        # the section as Falco, with its version
        for mate in 1 2; do
            sed '1s/^##Falco /##Falco\t/' falco_output/${prefix}_\${mate}/fastqc_data.txt > ${prefix}_\${mate}_fastqc_data.txt
            mv falco_output/${prefix}_\${mate}/fastqc_report.html ${prefix}_\${mate}_falco_report.html
            mv falco_output/${prefix}_\${mate}/summary.txt ${prefix}_\${mate}_falco_summary.txt
        done
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'PASS\\tstub\\n' > ${prefix}.fastq_check.tsv
    for mate in 1 2; do
        touch ${prefix}_\${mate}_fastqc_data.txt ${prefix}_\${mate}_falco_report.html ${prefix}_\${mate}_falco_summary.txt
    done
    """
}
