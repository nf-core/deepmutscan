process DMSANALYSIS_AASEQ {
    tag "amino_acid_sequence"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/78/7815f74d3680b875aa8c7be8f78dc77921ec970b99b45a539c5d2225749c0b90/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-dplyr_r-ggplot2_pruned:a65d6fd383d85714' }"

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
