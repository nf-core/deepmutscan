process DMSANALYSIS_PROCESSGATK {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:0fd2e39a5bf2ecaa"

    publishDir "${params.outdir}/intermediate_files", mode: 'copy'

    input:
    tuple val(meta), path(variantCounts)
    path possible_mutations
    path aa_seq
    val min_counts
    path process_raw_gatk_script
    path filter_by_library_script
    path complete_gatk_script
    path prepare_counts_heatmap_script

    output:
    tuple val(meta),
        path("annotated_variantCounts.csv"),
        path("variantCounts_filtered_by_library.csv"),
        path("library_completed_variantCounts.csv"),
        path("variantCounts_for_heatmaps.csv"),
        emit: processed_variantCounts

    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    Rscript -e "source('$process_raw_gatk_script'); process_raw_gatk('$variantCounts', 'annotated_variantCounts.csv')"
    Rscript -e "source('$filter_by_library_script'); filter_gatk_by_codon_library('annotated_variantCounts.csv', '$possible_mutations', 'variantCounts_filtered_by_library.csv')"
    Rscript -e "source('$complete_gatk_script'); complete_prefiltered_gatk('$possible_mutations', 'variantCounts_filtered_by_library.csv', 'library_completed_variantCounts.csv')"
    Rscript -e "source('$prepare_counts_heatmap_script'); prepare_gatk_data_for_counts_heatmaps('variantCounts_filtered_by_library.csv', '$aa_seq', 'variantCounts_for_heatmaps.csv', $min_counts)"

    R_VERSION=\$(R --version | head -n 1 | sed -E 's/^R version ([0-9.]+).*/\\1/')
    cat <<-END_VERSIONS > versions.yml
    DMSANALYSIS_PROCESSGATK:
      r-base: \$R_VERSION
    END_VERSIONS
    """

    stub:
    """
    touch annotated_variantCounts.csv variantCounts_filtered_by_library.csv library_completed_variantCounts.csv variantCounts_for_heatmaps.csv
    echo "DMSANALYSIS_PROCESSGATK:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
    """
}
