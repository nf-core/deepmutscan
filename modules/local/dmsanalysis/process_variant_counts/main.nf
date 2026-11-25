// One library's count tables: annotated, filtered by the library design, completed with the design's
// unobserved variants, and summarised for the count heatmaps; with --error_correction false_doubles also
// its false-double-mutant (FDM) error evidence (see templates/process_variant_counts.R).
process DMSANALYSIS_PROCESS_VARIANT_COUNTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e2/e2d657c330e5ac4dcc3b40ef513d377367e60c293e01003906877e40f6865d75/data'
        : 'community.wave.seqera.io/library/bioconductor-edger_bioconductor-impute_bioconductor-limma_dms-fdmcorrect_pruned:be9da4b29f56a3d0'}"

    input:
    tuple val(meta), path(variantCounts)
    path possible_mutations
    path aa_seq
    path wt_fasta              // the reference the reads were aligned to, for the FDM evidence
    val fdm_evidence           // true: also the FDM evidence (--error_correction false_doubles)
    val max_fd_dist
    val min_fd_dist

    output:
    tuple val(meta),
        path("annotated_variantCounts.csv"),
        path("variantCounts_filtered_by_library.csv"),
        path("library_completed_variantCounts.csv"),
        path("variantCounts_for_heatmaps.csv"),
        emit: processed_variantCounts
    tuple val(meta), path("fdm_evidence.csv")        , emit: evidence, optional: true
    tuple val(meta), path("fdm_evidence_summary.csv"), emit: summary , optional: true
    tuple val(meta), path("fdm_distance_profile.csv"), emit: profile , optional: true
    tuple val(meta), path("fdm_exposure_by_distance.csv"), emit: exposure, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'process_variant_counts.R'

    stub:
    def fdm = fdm_evidence.toString() == 'true'
    """
    touch annotated_variantCounts.csv variantCounts_filtered_by_library.csv library_completed_variantCounts.csv variantCounts_for_heatmaps.csv
    if ${fdm}; then
        echo "base_mut,codon_mut,codon,counts,cov,type,context,F,M,n_anchor_codons,n_fd_codons,n_fd_candidates" > fdm_evidence.csv
        touch fdm_evidence_summary.csv
        echo "dist,fd_reads,fd_exposure,gd_reads,gd_exposure,library" > fdm_distance_profile.csv
        echo "base_mut,dist,m,library" > fdm_exposure_by_distance.csv
    fi
    echo "DMSANALYSIS_PROCESS_VARIANT_COUNTS:" > versions.yml
    echo "  stub-version: 0.0.0" >> versions.yml
    """
}
