process DMSANALYSIS_POSSIBLE_MUTATIONS {
    tag "table /w all possible variants"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:0fd2e39a5bf2ecaa"

    publishDir "${params.outdir}/intermediate_files", mode: 'copy'

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

    # Capture R base and packages versions
    echo "DMSANALYSIS_POSSIBLE_MUTATIONS:" > versions.yml
    echo "  r-base: \$(R --version | head -n 1 | sed 's/^R version //')" >> versions.yml
    echo "  biostrings: \$(Rscript -e 'packageVersion(\"Biostrings\")' | tail -n 1 | sed 's/\\[1\\] //')" >> versions.yml    
    """

    stub:
    """
    touch possible_mutations.csv
    echo "DMSANALYSIS_POSSIBLE_MUTATIONS:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
    """
}
