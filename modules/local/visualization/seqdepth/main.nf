process VISUALIZATION_SEQDEPTH {
    tag "$meta.id"
    // Closed-form hypergeometric rarefaction - fast, so no longer a heavy process.
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/7988a2e369e98cec897d48bab3e55daddd8a973d9b7491ad27b6b15d550a4c5b/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-dplyr_r-ggplot2_pruned:32f1af0a83a018e7'}"

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
