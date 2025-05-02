process DMSANALYSIS_AASEQ {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-biocmanager_r-dplyr_pruned:8fa94107068d5af9"

    input:
    tuple val(meta), path(wt_seq)
    val pos_range

    output:
    tuple val(meta), path("aa_seq.txt"), emit: aa_seq
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    start_stop_codon="$pos_range"

    Rscript -e "source('../bin/aa_seq.R'); aa_seq('$wt_seq', '\$start_stop_codon', 'aa_seq.txt')"

    cat <<-END_VERSIONS > versions.yml
    "\${task.process}":
        r: \$(R --version | head -n 1 | sed 's/^R version //')
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
