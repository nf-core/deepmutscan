process DMSANALYSIS_AASEQ {
    tag "amino_acid_sequence"
    label 'process_single'
    label 'r_env'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(wt_seq)
    val pos_range

    output:
    tuple val(meta), path("aa_seq.txt"), emit: aa_seq
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'aa_seq.R'

    stub:
    """
    touch aa_seq.txt
    echo "DMSANALYSIS_AASEQ:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
    """
}
