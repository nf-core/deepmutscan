process VISUALIZATION_SEQDEPTH {
    tag "$meta.id"
    // Closed-form hypergeometric rarefaction - fast, so no longer a heavy process.
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/78/7815f74d3680b875aa8c7be8f78dc77921ec970b99b45a539c5d2225749c0b90/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-dplyr_r-ggplot2_pruned:a65d6fd383d85714' }"

    input:
    tuple val(meta), path(variantCounts_filtered_by_library)
    path possible_mutations

    output:
    tuple val(meta), path("SeqDepth.pdf"), emit: SeqDepth
    tuple val(meta), path("seqdepth_curve.csv"), emit: curve
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'SeqDepth_simulation.R'

    stub:
    """
    touch SeqDepth.pdf
    touch seqdepth_curve.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: "0.0.0"
        r-ggplot2: "0.0.0"
    END_VERSIONS
    """
}
