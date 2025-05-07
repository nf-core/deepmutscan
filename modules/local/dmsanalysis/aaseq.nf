process DMSANALYSIS_AASEQ {
    tag "amino_acid_sequence"
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
    
    R_version=\$(R --version | head -n 1 | sed 's/^R version //')
    
    Rscript -e "source('$script'); aa_seq('$wt_seq', '\$start_stop_codon', 'aa_seq.txt')"
    
    # Extract R base and packages versions
    R_VERSION=\$(R --version | head -n 1 | sed -E 's/^R version ([0-9.]+).*/\\1/')
    BIOSTRINGS_VERSION=\$(Rscript -e "packageVersion('Biostrings')" | grep -Eo '[0-9]+(\\.[0-9]+)+')
    cat <<-END_VERSIONS > versions.yml
    DMSANALYSIS_AASEQ:
      r-base: \$R_VERSION
      biostrings: \$BIOSTRINGS_VERSION
    END_VERSIONS   
    """

    stub:
    """
    touch aa_seq.txt
    echo "DMSANALYSIS_AASEQ:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
    """
}
