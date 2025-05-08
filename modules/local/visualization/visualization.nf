process DMSANALYSIS_VISUALIZATION_COUNTS_PER_COV {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:0fd2e39a5bf2ecaa"

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
    DMSANALYSIS_VISUALIZATION_COUNTS_PER_COV:
      r-base: \$R_VERSION
    END_VERSIONS
    """

    stub:
    """
    touch counts_per_cov_heatmap.pdf
    echo "DMSANALYSIS_VISUALIZATION_COUNTS_PER_COV:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
    """
}

