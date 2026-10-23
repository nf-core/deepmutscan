process FASTQ_CHECK {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fq:0.12.0--h9ee0642_0':
        'quay.io/biocontainers/fq:0.12.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(reads, arity: '2')

    output:
    tuple val(meta), path("*.fastq_check.tsv"), emit: report
    tuple val("${task.process}"), val('fq'), eval("fq lint --version | sed 's/fq-lint //; s/ .*//'"), emit: versions_fq, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def r1     = reads[0]
    def r2     = reads[1]
    """
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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'PASS\\tstub\\n' > ${prefix}.fastq_check.tsv
    """
}
