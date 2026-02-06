process DMSANALYSIS_POSSIBLE_MUTATIONS {
    tag "table /w all possible variants"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity'
        ? 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:ce2ba7ad7f6e7f2c'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:0fd2e39a5bf2ecaa' }"

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
