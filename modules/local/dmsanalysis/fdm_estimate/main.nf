// False-double-mutant (FDM) error rates for one pool of libraries (`--false_doubles_pooling`):
// sums F and M per variant across the pool and fits the estimator once (see templates/fdm_estimate.R).
process DMSANALYSIS_FDM_ESTIMATE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'REPLACE_WITH_SINGULARITY_HTTPS_URL'
        : 'REPLACE_WITH_DOCKER_IMAGE'}"

    input:
    tuple val(meta), path(evidence, stageAs: 'evidence_??/*'), path(profiles, stageAs: 'profile_??/*'), path(exposures, stageAs: 'exposure_??/*')
    val method
    val prior
    val max_fd_dist
    val subtract_genuine_doubles

    output:
    tuple val(meta), path("fdm_rates.csv")    , emit: rates
    tuple val(meta), path("fdm_prior_fit.csv"), emit: prior_fit
    tuple val(meta), path("fdm_rates_per_library.csv"), emit: rates_per_library
    tuple val(meta), path("fdm_kappa_profile.csv"), emit: kappa_profile
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'fdm_estimate.R'

    stub:
    """
    echo "base_mut,type,context,F,M,kappa_exposure,G,F_error,prior_mean,prior_nu,shrinkage_weight,e" > fdm_rates.csv
    touch fdm_prior_fit.csv
    echo "base_mut,F,M,e,shrinkage_weight,library,pool" > fdm_rates_per_library.csv
    echo "dist,gd_reads,gd_exposure,kappa_raw,kappa,pool" > fdm_kappa_profile.csv
    echo '"${task.process}":' > versions.yml
    echo '    r-base: 0.0.0' >> versions.yml
    echo '    dms-fdmcorrect: 0.0.0' >> versions.yml
    """
}
