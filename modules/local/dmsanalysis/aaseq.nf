process DMSANALYSIS_AASEQ {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:0fd2e39a5bf2ecaa"

    input:
    tuple val(meta), path(wt_seq)
    val pos_range
    path script  // aa_seq.R

    output:
    tuple val(meta), path("aa_seq.txt"), emit: aa_seq
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    start_stop_codon="$pos_range"

    Rscript -e "source('$script'); aa_seq('$wt_seq', '\$start_stop_codon', 'aa_seq.txt')"

    cat <<-END_VERSIONS > versions.yml
    "DMSANALYSIS_AASEQ":
        r: $(R --version | head -n 1 | sed 's/^R version //')
        aa_seq_script: custom_R_script_bin
    END_VERSIONS

    """

    stub:
    """
    touch aa_seq.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dmsanalysis_aaseq: dummy_version
    END_VERSIONS
    """
}
