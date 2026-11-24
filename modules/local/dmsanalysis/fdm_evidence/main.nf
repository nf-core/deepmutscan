// False-double-mutant (FDM) error evidence for one library: per single-nucleotide variant, the
// false-double reads F and the anchor-read exposure M (see templates/fdm_evidence.R).
process DMSANALYSIS_FDM_EVIDENCE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/76/76821c826dabd8d88eb42b9075e8648d371a14b560ac3dc755dafd5dfc46fb07/data'
        : 'community.wave.seqera.io/library/dms-fdmcorrect_r-base:373aad1776732a61'}"

    input:
    tuple val(meta), path(raw_counts), path(filtered)
    path wt_fasta
    val max_fd_dist
    val min_fd_dist

    output:
    tuple val(meta), path("fdm_evidence.csv")        , emit: evidence
    tuple val(meta), path("fdm_evidence_summary.csv"), emit: summary
    tuple val(meta), path("fdm_distance_profile.csv"), emit: profile
    tuple val(meta), path("fdm_exposure_by_distance.csv"), emit: exposure
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'fdm_evidence.R'

    stub:
    """
    echo "base_mut,codon_mut,codon,counts,cov,type,context,F,M,n_anchor_codons,n_fd_codons,n_fd_candidates" > fdm_evidence.csv
    touch fdm_evidence_summary.csv
    echo "dist,fd_reads,fd_exposure,gd_reads,gd_exposure,library" > fdm_distance_profile.csv
    echo "base_mut,dist,m,library" > fdm_exposure_by_distance.csv
    echo '"${task.process}":' > versions.yml
    echo '    r-base: 0.0.0' >> versions.yml
    echo '    dms-fdmcorrect: 0.0.0' >> versions.yml
    """
}
