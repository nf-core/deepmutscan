process PREPARE_REFERENCE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2c/2cf682587a97332817d9ebb3ea7fdf5e0e73bbb47e084982871f35a0c87a7043/data'
        : 'community.wave.seqera.io/library/bwa-mem3_gawk:d06b0c4e4fd28ed3'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.uppercase.fasta")                          , emit: fasta
    tuple val(meta), path("bwamem3")                                    , emit: index
    tuple val("${task.process}"), val('gawk'), eval("gawk --version | sed '1!d; s/.*Awk //; s/,.*//'"), topic: versions, emit: versions_gawk
    tuple val("${task.process}"), val('bwa-mem3'), eval("bwa-mem3 version | sed -nE '1 s/^([0-9]+(\\.[0-9]+)+).*/\\1/p'"), topic: versions, emit: versions_bwamem3

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''   // bwa-mem3 index
    def prefix = task.ext.prefix ?: fasta.name.replaceFirst(/\.gz$/, '').replaceFirst(/\.(fasta|fas|fna|fa)$/, '')
    def uppercase = "${prefix}.uppercase.fasta"
    if (fasta.name == uppercase) {
        error "PREPARE_REFERENCE: the reference is already named ${uppercase}; rename it, or set ext.prefix"
    }
    // The reference in upper case, before any step reads it. Lower-case (soft-masked) bases are the
    // same bases, but the counter, the error correction and the fitness steps compare reference and
    // read bases as characters. Header lines are left as they are; a Windows line end is dropped.
    // gzip -cdf reads a gzipped reference and passes a plain one through.
    """
    gzip -cdf ${fasta} \\
        | gawk '{ sub(/\\r\$/, "") } /^>/ { print; next } { print toupper(\$0) }' \\
        > ${uppercase}

    mkdir bwamem3
    bwa-mem3 \\
        index \\
        ${args} \\
        -t ${task.cpus} \\
        -p bwamem3/${uppercase} \\
        ${uppercase}
    """

    stub:
    def prefix = task.ext.prefix ?: fasta.name.replaceFirst(/\.gz$/, '').replaceFirst(/\.(fasta|fas|fna|fa)$/, '')
    def uppercase = "${prefix}.uppercase.fasta"
    """
    touch ${uppercase}
    mkdir bwamem3
    touch bwamem3/${uppercase}.{amb,ann,bwt.2bit.64,pac}
    """
}
