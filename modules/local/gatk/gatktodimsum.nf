process GATK_GATKTODIMSUM {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' 
        ? 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:ce2ba7ad7f6e7f2c' 
        : 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:0fd2e39a5bf2ecaa' }"

    input:
    tuple val(meta), path(variantCounts_filtered_by_library)
    path wt_seq
    val pos_range
    path script // R script

    output:
    tuple val(meta), path("dimsum_input.tsv"), emit: dimsum_input
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    start_stop_codon="$pos_range"

    Rscript -e "source('$script'); generate_dimsum_input('$wt_seq', '$variantCounts_filtered_by_library', '\$start_stop_codon', 'dimsum_input.tsv')"

    R_VERSION=\$(R --version | head -n 1 | sed -E 's/^R version ([0-9.]+).*/\\1/')
    cat <<-END_VERSIONS > versions.yml
    GATK_GATKTODIMSUM:
      r-base: \$R_VERSION
    END_VERSIONS
    """

    stub:
    """
    touch dimsum_input.tsv
    echo "GATK_GATKTODIMSUM:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
    """
}
