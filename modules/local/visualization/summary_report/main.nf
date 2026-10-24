// Build the single, self-contained all-in-one HTML report for a run.
//
// Every artefact the report shows is passed through one generic pair of inputs: a base64-encoded
// JSON manifest and the matching staged files, in the same order. That keeps the process signature
// stable no matter which optional parts of the pipeline ran - a section simply disappears when its
// entries are absent from the manifest, so no channel juggling is needed for --fitness/--pdb/EC.
process VISUALIZATION_SUMMARY_REPORT {
    tag "${meta.sample ?: meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/97/977734114fdf9ecc0fbe96d3b2f09901eb1c4490de64d863a18679d926a028c4/data'
        : 'community.wave.seqera.io/library/numpy_pandas_python:43a80707feb1fac3' }"

    input:
    tuple val(meta), val(run_b64)            // base64(JSON) run metadata (params, samples, versions)
    val manifest_b64                         // base64(JSON) [{kind,group,name}, ...] aligned with `assets`
    path assets, stageAs: 'a*/*'             // pdf/csv/tsv/html artefacts; numbered dirs avoid basename clashes
    path template_html
    path threedmol_js                        // 3Dmol.js library, inlined only when the run has a structure

    output:
    tuple val(meta), path("deepmutscan_report.html"), emit: report
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'build_summary_report.py'

    stub:
    """
    touch deepmutscan_report.html
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "0.0.0"
        pandas: "0.0.0"
    END_VERSIONS
    """
}
