// Build the single, self-contained interactive variant effect inspection tool (HTML).
// Aggregates per-residue metrics (fitness, pLDDT accuracy, coverage, counts, counts/cov, error bias),
// aligns the wildtype sequence onto the structure, and inlines 3Dmol.js + PDB + data into
// one portable HTML file (plus a structure_data.json sidecar).
process VARIANT_EFFECT_INSPECTION_TOOL {
    tag "${meta.sample ?: meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pandas:2.2.1'
        : 'quay.io/biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(pdb), path(fitness_tsv)
    path(variant_counts, stageAs: 'vc*/*')   // per-sample files share a basename -> stage in numbered dirs
    path aa_seq
    path threedmol_js
    path template_html

    output:
    tuple val(meta), path("variant_effect_inspection_tool.html"), emit: viewer
    tuple val(meta), path("structure_data.json")                , emit: data
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'build_variant_effect_inspection_tool.py'

    stub:
    """
    touch variant_effect_inspection_tool.html structure_data.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "0.0.0"
        pandas: "0.0.0"
    END_VERSIONS
    """
}
