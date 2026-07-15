// Build the self-contained per-sample error-correction HTML report.
process VISUALIZATION_ERROR_CORRECTION_REPORT {
    tag "${meta.sample ?: meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pandas:2.2.1'
        : 'quay.io/biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(corrected)
    val method
    path template_html

    output:
    tuple val(meta), path("*_error_correction_report.html"), emit: report
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'build_ec_report.py'

    stub:
    """
    touch ${meta.id}_error_correction_report.html
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "0.0.0"
        pandas: "0.0.0"
    END_VERSIONS
    """
}
