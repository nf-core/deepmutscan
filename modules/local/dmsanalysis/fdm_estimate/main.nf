// False-double-mutant (FDM) error rates for one pool of libraries (`--false_doubles_pooling`), and the
// correction of each of its libraries: sums F and M per variant across the pool, fits the estimator once
// and applies the rates to every library of the pool (see templates/fdm_estimate.R).
process DMSANALYSIS_FDM_ESTIMATE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e2/e2d657c330e5ac4dcc3b40ef513d377367e60c293e01003906877e40f6865d75/data'
        : 'community.wave.seqera.io/library/bioconductor-edger_bioconductor-impute_bioconductor-limma_dms-fdmcorrect_pruned:be9da4b29f56a3d0'}"

    input:
    tuple val(meta), val(libraries), val(library_ids),        // the pool's library metas and their ids (space-separated), by id
        path(evidence, stageAs: 'evidence_??/*'), path(profiles, stageAs: 'profile_??/*'), path(exposures, stageAs: 'exposure_??/*'),
        path(filtered, stageAs: 'filtered_??/*'), path(completed, stageAs: 'completed_??/*')
    path wt_fasta
    path aa_seq
    val method
    val prior
    val max_fd_dist
    val subtract_genuine_doubles

    output:
    tuple val(meta), path("fdm_rates.csv")    , emit: rates
    tuple val(meta), path("fdm_prior_fit.csv"), emit: prior_fit
    tuple val(meta), path("fdm_rates_per_library.csv"), emit: rates_per_library
    tuple val(meta), path("fdm_kappa_profile.csv"), emit: kappa_profile
    // per library of the pool, each in a folder named after the library (meta.id)
    tuple val(meta), val(libraries),
        path("*/variantCounts_filtered_by_library_error_corrected.csv"),
        path("*/variantCounts_for_heatmaps_error_corrected.csv"),
        path("*/library_completed_variantCounts_error_corrected.csv"),
        emit: corrected
    tuple val(meta), val(libraries), path("*/seq_error_rate.csv"), emit: seq_error_rate
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'fdm_estimate.R'

    stub:
    def stub_libs = library_ids.tokenize(' ')
    def stub_corrected = stub_libs.collect { id ->
        "mkdir -p ${id} && touch ${id}/variantCounts_filtered_by_library_error_corrected.csv ${id}/variantCounts_for_heatmaps_error_corrected.csv ${id}/library_completed_variantCounts_error_corrected.csv ${id}/seq_error_rate.csv"
    }.join('\n    ')
    """
    ${stub_corrected}
    echo "base_mut,type,context,F,M,kappa_exposure,G,F_error,prior_mean,prior_nu,shrinkage_weight,e" > fdm_rates.csv
    touch fdm_prior_fit.csv
    echo "base_mut,F,M,e,shrinkage_weight,library,pool" > fdm_rates_per_library.csv
    echo "dist,gd_reads,gd_exposure,kappa_raw,kappa,pool" > fdm_kappa_profile.csv
    echo '"${task.process}":' > versions.yml
    echo '    r-base: 0.0.0' >> versions.yml
    echo '    dms-fdmcorrect: 0.0.0' >> versions.yml
    """
}
