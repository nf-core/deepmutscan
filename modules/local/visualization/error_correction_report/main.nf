// Build the self-contained error-correction HTML report for the whole run.
//
// One report, not one per file: every sequencing file's numbers travel in the page and each panel
// recomputes in the browser from whichever files the reader selects, so a single-file view is still
// one click away while cross-replicate spread (the error bars) becomes possible at all.
process VISUALIZATION_ERROR_CORRECTION_REPORT {
    tag "${sample}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pandas:2.2.1'
        : 'quay.io/biocontainers/pandas:2.2.1' }"

    input:
    val sample                               // run-level name, for the heading
    val ids_b64                              // base64(JSON) [file id, ...] aligned with `corrected`
    path corrected, stageAs: 'c*/*'          // per-file tables share a basename -> numbered dirs
    val method
    path template_html

    output:
    path "error_correction_report.html", emit: report
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'build_ec_report.py'

    stub:
    """
    touch error_correction_report.html
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "0.0.0"
        pandas: "0.0.0"
    END_VERSIONS
    """
}
