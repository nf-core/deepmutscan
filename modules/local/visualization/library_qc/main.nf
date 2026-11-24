process VISUALIZATION_LIBRARY_QC {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/7988a2e369e98cec897d48bab3e55daddd8a973d9b7491ad27b6b15d550a4c5b/data'
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-dplyr_r-ggplot2_pruned:32f1af0a83a018e7'}"

    input:
    tuple val(meta), path(variantCounts_filtered_by_library), path(variantCounts_for_heatmaps), path(library_completed_variantCounts)
    path aa_seq
    path possible_mutations
    val sliding_window_size
    val aimed_cov
    val run_seqdepth

    output:
    tuple val(meta), path("counts_heatmap.pdf"),         emit: counts_heatmap
    tuple val(meta), path("counts_per_cov_heatmap.pdf"), emit: counts_per_cov_heatmap
    tuple val(meta), path("rolling_counts.pdf"),         emit: rolling_counts
    tuple val(meta), path("rolling_counts_per_cov.pdf"), emit: rolling_counts_per_cov
    tuple val(meta), path("rolling_coverage.pdf"),       emit: rolling_coverage
    tuple val(meta), path("logdiff_plot.pdf"),           emit: logdiff_plot
    tuple val(meta), path("logdiff_varying_bases.pdf"),  emit: logdiff_varying_bases
    tuple val(meta), path("SeqDepth.pdf"),               emit: SeqDepth, optional: true
    tuple val(meta), path("seqdepth_curve.csv"),         emit: curve,    optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'library_qc.R'

    stub:
    """
    touch counts_heatmap.pdf counts_per_cov_heatmap.pdf rolling_counts.pdf rolling_counts_per_cov.pdf \\
          rolling_coverage.pdf logdiff_plot.pdf logdiff_varying_bases.pdf
    if [ "${run_seqdepth}" = "true" ]; then
        touch SeqDepth.pdf seqdepth_curve.csv
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: "0.0.0"
        r-dplyr: "0.0.0"
        r-ggplot2: "0.0.0"
        r-scales: "0.0.0"
        r-zoo: "0.0.0"
    END_VERSIONS
    """
}
