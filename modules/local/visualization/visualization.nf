process VISUALIZATION_COUNTS_PER_COV {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' 
        ? 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:ce2ba7ad7f6e7f2c' 
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:0fd2e39a5bf2ecaa' }"

    input:
    tuple val(meta), path(variantCounts_for_heatmaps)
    val min_counts
    path script // counts_per_cov_heatmap.R

    output:
    tuple val(meta), path("counts_per_cov_heatmap.pdf"), emit: counts_per_cov_heatmap
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    Rscript -e "source('$script'); counts_per_cov_heatmap('$variantCounts_for_heatmaps', $min_counts, 'counts_per_cov_heatmap.pdf')"

    R_VERSION=\$(R --version | head -n 1 | sed -E 's/^R version ([0-9.]+).*/\\1/')
    cat <<-END_VERSIONS > versions.yml
    VISUALIZATION_COUNTS_PER_COV:
      r-base: \$R_VERSION
    END_VERSIONS
    """

    stub:
    """
    touch counts_per_cov_heatmap.pdf
    echo "VISUALIZATION_COUNTS_PER_COV:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
    """
}

process VISUALIZATION_COUNTS_HEATMAP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' 
        ? 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:ce2ba7ad7f6e7f2c' 
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:0fd2e39a5bf2ecaa' }"

    input:
    tuple val(meta), path(variantCounts_for_heatmaps)
    val min_counts
    path script // counts_heatmap.R

    output:
    tuple val(meta), path("counts_heatmap.pdf"), emit: counts_heatmap
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    Rscript -e "source('$script'); counts_heatmap('$variantCounts_for_heatmaps', $min_counts, 'counts_heatmap.pdf')"

    R_VERSION=\$(R --version | head -n 1 | sed -E 's/^R version ([0-9.]+).*/\\1/')
    cat <<-END_VERSIONS > versions.yml
    VISUALIZATION_COUNTS_PER_COV:
      r-base: \$R_VERSION
    END_VERSIONS
    """

    stub:
    """
    touch counts_heatmap.pdf
    echo "VISUALIZATION_COUNTS_HEATMAP:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
    """
}

process VISUALIZATION_GLOBAL_POS_BIASES_COUNTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' 
        ? 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:ce2ba7ad7f6e7f2c' 
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:0fd2e39a5bf2ecaa' }"

    input:
    tuple val(meta), path(variantCounts_filtered_by_library)
    path aa_seq
    val sliding_window_size
    path script // global_position_biases_counts_and_counts_per_cov.R

    output:
    tuple val(meta), path("rolling_counts.pdf"), emit: rolling_counts
    tuple val(meta), path("rolling_counts_per_cov.pdf"), emit: rolling_counts_per_cov
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    Rscript -e "source('$script'); position_biases('$variantCounts_filtered_by_library', '$aa_seq', $sliding_window_size)"

    R_VERSION=\$(R --version | head -n 1 | sed -E 's/^R version ([0-9.]+).*/\\1/')
    cat <<-END_VERSIONS > versions.yml
    VISUALIZATION_COUNTS_PER_COV:
      r-base: \$R_VERSION
    END_VERSIONS
    """

    stub:
    """
    touch rolling_counts.pdf
    touch rolling_counts_per_cov.pdf
    echo "VISUALIZATION_COUNTS_HEATMAP:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
    """
}

process VISUALIZATION_GLOBAL_POS_BIASES_COV {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' 
        ? 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:ce2ba7ad7f6e7f2c' 
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:0fd2e39a5bf2ecaa' }"

    input:
    tuple val(meta), path(variantCounts_filtered_by_library)
    path aa_seq
    val sliding_window_size
    val aimed_cov
    path script // R script

    output:
    tuple val(meta), path("rolling_coverage.pdf"), emit: rolling_coverage
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    Rscript -e "source('$script'); position_biases('$variantCounts_filtered_by_library', '$aa_seq', $sliding_window_size, 'rolling_coverage.pdf', $aimed_cov)"

    R_VERSION=\$(R --version | head -n 1 | sed -E 's/^R version ([0-9.]+).*/\\1/')
    cat <<-END_VERSIONS > versions.yml
    VISUALIZATION_COUNTS_PER_COV:
      r-base: \$R_VERSION
    END_VERSIONS
    """

    stub:
    """
    touch rolling_coverage.pdf
    echo "VISUALIZATION_COUNTS_HEATMAP:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
    """
}

process VISUALIZATION_LOGDIFF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' 
        ? 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:ce2ba7ad7f6e7f2c' 
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:0fd2e39a5bf2ecaa' }"

    input:
    tuple val(meta), path(library_completed_variantCounts)
    path script // R script

    output:
    tuple val(meta), path("logdiff_plot.pdf"), emit: logdiff_plot
    tuple val(meta), path("logdiff_varying_bases.pdf"), emit: logdiff_varying_bases
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    Rscript -e "source('$script'); logdiff_plot('$library_completed_variantCounts')"

    R_VERSION=\$(R --version | head -n 1 | sed -E 's/^R version ([0-9.]+).*/\\1/')
    cat <<-END_VERSIONS > versions.yml
    VISUALIZATION_COUNTS_PER_COV:
      r-base: \$R_VERSION
    END_VERSIONS
    """

    stub:
    """
    touch logdiff_plot.pdf
    touch logdiff_varying_bases.pdf
    echo "VISUALIZATION_COUNTS_HEATMAP:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
    """
}

process VISUALIZATION_SEQDEPTH {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' 
        ? 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:ce2ba7ad7f6e7f2c' 
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:0fd2e39a5bf2ecaa' }"

    input:
    tuple val(meta), path(variantCounts_filtered_by_library)
    path possible_mutations
    val min_counts
    path script // R script

    output:
    tuple val(meta), path("SeqDepth.pdf"), emit: SeqDepth
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    Rscript -e "source('$script'); SeqDepth_simulation_plot('$variantCounts_filtered_by_library', '$possible_mutations', 'SeqDepth.pdf', 0.01, $min_counts)"

    R_VERSION=\$(R --version | head -n 1 | sed -E 's/^R version ([0-9.]+).*/\\1/')
    cat <<-END_VERSIONS > versions.yml
    VISUALIZATION_COUNTS_PER_COV:
      r-base: \$R_VERSION
    END_VERSIONS
    """

    stub:
    """
    touch SeqDepth.pdf
    echo "VISUALIZATION_COUNTS_HEATMAP:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
    """
}
