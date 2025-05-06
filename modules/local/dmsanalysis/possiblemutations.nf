process DMSANALYSIS_POSSIBLE_MUTATIONS {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:0fd2e39a5bf2ecaa"

    input:
    path wt_seq
    val pos_range
    val mutagenesis_type
    path custom_codon_library
    path script // possible_mutations.R

    output:
    path "possible_mutations.csv", emit: possible_mutations
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    start_stop_codon="$pos_range"

    if [[ "$custom_codon_library" == "/NULL" ]]; then
        Rscript -e "source('$script'); generate_possible_variants('$wt_seq', '\$start_stop_codon', '$mutagenesis_type', NULL, 'possible_mutations.csv')"
    else
        Rscript -e "source('$script'); generate_possible_variants('$wt_seq', '\$start_stop_codon', '$mutagenesis_type', '$custom_codon_library', 'possible_mutations.csv')"
    fi

    echo "." > versions.yml
    """

    stub:
    """
    touch possible_mutations.csv
    echo "stub" > versions.yml
    """
}
