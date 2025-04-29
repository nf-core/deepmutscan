process GATK_SATURATIONMUTAGENESIS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/gatk_samtools:YOUR_TAG_HERE"

    input:
    tuple val(meta), path(premerged_reads)
    path wt_seq
    val pos_range
    val min_counts

    output:
    tuple val(meta), path("gatk_output.*"), emit: gatk_output
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Index reference
    samtools faidx $wt_seq
    gatk CreateSequenceDictionary -R $wt_seq

    # Read start and stop codon from input
    start_stop_codon=\$(cat $pos_range)

    # Run GATK AnalyzeSaturationMutagenesis
    gatk AnalyzeSaturationMutagenesis \
        -I $premerged_reads \
        -R $wt_seq \
        --orf \$start_stop_codon \
        --paired-mode false \
        --min-q 30 \
        --min-variant-obs $min_counts \
        -O gatk_output

    # Save versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version |& sed '1!d ; s/samtools //')
        gatk: \$(gatk --version |& sed 's/^.*GATK/\1/')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch gatk_output.variantCounts
    touch versions.yml
    """
}
