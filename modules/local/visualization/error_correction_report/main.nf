// Build the self-contained error-correction HTML report for the whole run.
//
// One report, not one per file: every sequencing file's numbers travel in the page and each panel
// recomputes in the browser from whichever files the reader selects, so a single-file view is still
// one click away while cross-replicate spread (the error bars) becomes possible at all.
process VISUALIZATION_ERROR_CORRECTION_REPORT {
    tag "${sample}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f7/f727259d394cff047a26f3e3fa3a43e1a7a1c9dc51a35e6def7aa90a2085bb7b/data'
        : 'community.wave.seqera.io/library/numpy_pandas_python:cfb38dbf0574e16d'}"

    input:
    val sample                               // run-level name, for the heading
    val ids_b64                              // base64(JSON) [file id, ...] aligned with `corrected`
    path corrected, stageAs: 'c*/*'          // per-file tables share a basename -> numbered dirs
    val method
    val false_doubles_method                 // 'mle' | 'eb' - names the estimator in the report heading
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
