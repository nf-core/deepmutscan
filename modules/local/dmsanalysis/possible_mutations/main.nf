process DMSANALYSIS_POSSIBLE_MUTATIONS {
    tag "table /w all possible variants"
    label 'process_single'
    label 'r_env'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(wt_seq)
    val pos_range
    val mutagenesis_type
    path custom_codon_library

    output:
    tuple val(meta), path("possible_mutations.csv"), emit: possible_mutations
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'possible_mutations.R'

    stub:
    """
    touch possible_mutations.csv
    echo "DMSANALYSIS_POSSIBLE_MUTATIONS:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
    """
}
