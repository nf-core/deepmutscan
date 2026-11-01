// Mutational-signature fit of the error spectrum of one FDM pool: the false-double error rates, the
// wildtype-template error rates (if the sample has a wildtype library) and their difference, each fitted to a
// candidate signature set by two-stage NNLS (see templates/error_signatures.py).
process DMSANALYSIS_ERROR_SIGNATURES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pandas:2.2.1'
        : 'quay.io/biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(evidence, stageAs: 'ev??/*'), path(summaries, stageAs: 'sm??/*'), path(wt_evidence, stageAs: 'wt??/*'), path(rates, stageAs: 'rt/*')
    path candidates
    path trinucleotides
    path descriptions
    val max_cosine
    val bootstraps
    val genuine_doubles

    output:
    tuple val(meta), path("error_signatures.json"), emit: json
    tuple val(meta), path("error_signatures.csv") , emit: fit
    tuple val(meta), path("error_spectra.csv")    , emit: spectra
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'error_signatures.py'

    stub:
    """
    echo '{}' > error_signatures.json
    echo 'pool,spectrum,signature,exposure,ci_low,ci_high,cosine,reads_accounted,n_signatures' > error_signatures.csv
    echo 'channel' > error_spectra.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "0.0.0"
        pandas: "0.0.0"
        numpy: "0.0.0"
    END_VERSIONS
    """
}
