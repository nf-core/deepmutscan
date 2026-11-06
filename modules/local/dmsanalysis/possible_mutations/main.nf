process DMSANALYSIS_POSSIBLE_MUTATIONS {
    tag "table /w all possible variants"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/7988a2e369e98cec897d48bab3e55daddd8a973d9b7491ad27b6b15d550a4c5b/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-dplyr_r-ggplot2_pruned:32f1af0a83a018e7'}"

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
