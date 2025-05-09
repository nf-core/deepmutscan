/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { BWA_INDEX              } from '../modules/nf-core/bwa/index/main'
include { BWA_MEM                } from '../modules/nf-core/bwa/mem/main'
include { BAMFILTER_DMS          } from '../modules/local/bamprocessing/bamfilteringdms'
include { PREMERGE               } from '../modules/local/bamprocessing/premerge'
include { GATK_SATURATIONMUTAGENESIS          } from '../modules/local/gatk/saturationmutagenesis'
include { DMSANALYSIS_AASEQ      } from '../modules/local/dmsanalysis/aaseq'
include { DMSANALYSIS_POSSIBLE_MUTATIONS      } from '../modules/local/dmsanalysis/possiblemutations'
include { DMSANALYSIS_PROCESS_GATK      } from '../modules/local/dmsanalysis/processgatk'
include { VISUALIZATION_COUNTS_PER_COV      } from '../modules/local/visualization/visualization'
include { VISUALIZATION_COUNTS_HEATMAP      } from '../modules/local/visualization/visualization'
include { VISUALIZATION_GLOBAL_POS_BIASES_COUNTS      } from '../modules/local/visualization/visualization'
include { VISUALIZATION_GLOBAL_POS_BIASES_COV      } from '../modules/local/visualization/visualization'
include { VISUALIZATION_LOGDIFF      } from '../modules/local/visualization/visualization'
include { VISUALIZATION_SEQDEPTH      } from '../modules/local/visualization/visualization'
include { GATK_GATKTODIMSUM          } from '../modules/local/gatk/gatktodimsum'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_dmscore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Params defaults
params.min_counts = params.min_counts ?: 3                                      // minimum counts for variant to be recognized. All variants<min_counts will be set to 0
params.mutagenesis_type = params.mutagenesis_type ?: 'nnk'                      // default library is set to nnk
params.custom_codon_library = params.custom_codon_library ?: '/NULL'            // when mutagenesis_type is set to >>custom<< this variable has to be path to .txt with custom library
params.sliding_window_size = params.sliding_window_size ?: 10                   // sliding window size to flatten graphs in plots (e.g. GLOBAL_POS_BIASES_COUNTS function)
params.aimed_cov = params.aimed_cov ?: 100                                      // aimed coverage (assuming equal spread) to visualize threshold in plots
params.run_seqdepth = params.run_seqdepth ?: false                              // creating seqdepth simulation plot, is computationally quite heavy. per default disabled.


// Define fasta file as channel (e.g. for BWA index)
Channel
    .fromPath(params.fasta, checkIfExists: true)
    .map { fasta -> tuple( { id: "ref" }, fasta ) }
    .set { ch_fasta }

// Define reading_frame as channel (e.g. for gatk function)
Channel
    .value(params.reading_frame)
    .set { reading_frame_ch }

// Define min_counts as channel (e.g. for gatk function) -> if not set - default: 3
Channel
    .value(params.min_counts)
    .set { min_counts_ch }

// Define custom library as channel
Channel
    .value(params.custom_codon_library)
    .set { custom_codon_library_ch }

// Define mutagenesis type as channel
Channel
    .value(params.mutagenesis_type)
    .set { mutagenesis_type_ch }

// Define sliding_window_size as channel (e.g. for GLOBAL_POS_BIASES_COUNTS function) -> if not set - default: 10
Channel
    .value(params.sliding_window_size)
    .set { sliding_window_size_ch }

// Define aimed_cov as channel (e.g. for GLOBAL_POS_BIASES_COVERAGE function) -> if not set - default: 100
Channel
    .value(params.aimed_cov)
    .set { aimed_cov_ch }

// Define if seqdepth plot should be executed
Channel
    .value(params.run_seqdepth)
    .set { run_seqdepth_ch }


// Define R scripts as channels
Channel.fromPath("modules/local/dmsanalysis/bin/SeqDepth_simulation.R", checkIfExists: true).set { seqdepth_simulation_script_ch }
Channel.fromPath("modules/local/dmsanalysis/bin/aa_seq.R", checkIfExists: true).set { aa_seq_script_ch }
Channel.fromPath("modules/local/dmsanalysis/bin/complete_prefiltered_gatk.R", checkIfExists: true).set { complete_gatk_script_ch }
Channel.fromPath("modules/local/dmsanalysis/bin/counts_heatmap.R", checkIfExists: true).set { counts_heatmap_script_ch }
Channel.fromPath("modules/local/dmsanalysis/bin/counts_per_cov_heatmap.R", checkIfExists: true).set { counts_per_cov_heatmap_script_ch }
Channel.fromPath("modules/local/dmsanalysis/bin/detect_codons.R", checkIfExists: true).set { detect_codons_script_ch }
Channel.fromPath("modules/local/dmsanalysis/bin/filter_gatk_by_codon_library.R", checkIfExists: true).set { filter_by_library_script_ch }
Channel.fromPath("modules/local/dmsanalysis/bin/fitness_heatmap.R", checkIfExists: true).set { fitness_heatmap_script_ch }
Channel.fromPath("modules/local/dmsanalysis/bin/gatk_to_dimsum.R", checkIfExists: true).set { gatk_to_dimsum_script_ch }
Channel.fromPath("modules/local/dmsanalysis/bin/global_position_biases_counts_and_counts_per_cov.R", checkIfExists: true).set { global_bias_counts_cov_script_ch }
Channel.fromPath("modules/local/dmsanalysis/bin/global_position_biases_cov.R", checkIfExists: true).set { global_bias_cov_script_ch }
Channel.fromPath("modules/local/dmsanalysis/bin/install_packages.R", checkIfExists: true).set { install_packages_script_ch }
Channel.fromPath("modules/local/dmsanalysis/bin/logdiff.R", checkIfExists: true).set { logdiff_script_ch }
Channel.fromPath("modules/local/dmsanalysis/bin/low_count_variants.R", checkIfExists: true).set { low_count_variants_script_ch }
Channel.fromPath("modules/local/dmsanalysis/bin/possible_mutations.R", checkIfExists: true).set { possible_mutations_script_ch }
Channel.fromPath("modules/local/dmsanalysis/bin/prepare_gatk_data_for_count_heatmaps.R", checkIfExists: true).set { prepare_counts_heatmap_script_ch }
Channel.fromPath("modules/local/dmsanalysis/bin/prepare_gatk_data_for_fitness_heatmap.R", checkIfExists: true).set { prepare_fitness_heatmap_script_ch }
Channel.fromPath("modules/local/dmsanalysis/bin/process_raw_gatk.R", checkIfExists: true).set { process_raw_gatk_script_ch }


workflow DMSCORE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'dmscore_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    //
    // MODULE: BWA Index
    //
    BWA_INDEX (
        ch_fasta
    )

    // Broadcast index to all samples
    ch_bwa_index = BWA_INDEX.out.index

    BWA_MEM (
    ch_samplesheet,           // each sample's reads: tuple [meta, reads]
    ch_bwa_index,             // broadcasted index files
    ch_fasta,                 // broadcasted fasta file
    false                     // sort_bam
    )

    BAMFILTER_DMS (
    BWA_MEM.out.bam
    )

    PREMERGE (
    BAMFILTER_DMS.out.bam,
    ch_fasta.map{ it[1] }     // extract fasta path from Tuple
    )

    GATK_SATURATIONMUTAGENESIS (
    PREMERGE.out.bam,             // merged reads (tuple(meta, merged_reads.fastq))
    ch_fasta.map{ it[1] },        // fasta path (path)
    reading_frame_ch,             // codon range (path)
    min_counts_ch                 // min_counts (val)
    )

    DMSANALYSIS_AASEQ (
    ch_fasta,
    reading_frame_ch,
    aa_seq_script_ch // path to aa_seq.R (defined at the top)
    )
    ch_versions = ch_versions.mix(DMSANALYSIS_AASEQ.out.versions)

    DMSANALYSIS_POSSIBLE_MUTATIONS(
    ch_fasta,
    reading_frame_ch,          // pos_range (as val)
    mutagenesis_type_ch,       // mutagenesis_type (as val)
    custom_codon_library_ch,   // custom_codon_library (as path)
    possible_mutations_script_ch  // path to R script
    )
    ch_versions = ch_versions.mix(DMSANALYSIS_POSSIBLE_MUTATIONS.out.versions)

    DMSANALYSIS_PROCESS_GATK(
    GATK_SATURATIONMUTAGENESIS.out.variantCounts,
    DMSANALYSIS_POSSIBLE_MUTATIONS.out.possible_mutations.map{ it[1] },
    DMSANALYSIS_AASEQ.out.aa_seq.map{ it[1] },
    min_counts_ch,
    process_raw_gatk_script_ch,
    filter_by_library_script_ch,
    complete_gatk_script_ch,
    prepare_counts_heatmap_script_ch
    )

    annotated_variantCounts_ch = DMSANALYSIS_PROCESS_GATK.out.processed_variantCounts.map { meta, a, b, c, d -> tuple(meta, a) }
    variantCounts_filtered_by_library_ch = DMSANALYSIS_PROCESS_GATK.out.processed_variantCounts.map { meta, a, b, c, d -> tuple(meta, b) }
    library_completed_variantCounts_ch = DMSANALYSIS_PROCESS_GATK.out.processed_variantCounts.map { meta, a, b, c, d -> tuple(meta, c) }
    variantCounts_for_heatmaps_ch = DMSANALYSIS_PROCESS_GATK.out.processed_variantCounts.map { meta, a, b, c, d -> tuple(meta, d) }

    VISUALIZATION_COUNTS_PER_COV(
    variantCounts_for_heatmaps_ch,
    min_counts_ch,
    counts_per_cov_heatmap_script_ch
    )

    VISUALIZATION_COUNTS_HEATMAP(
    variantCounts_for_heatmaps_ch,
    min_counts_ch,
    counts_heatmap_script_ch
    )

    VISUALIZATION_GLOBAL_POS_BIASES_COUNTS(
    variantCounts_filtered_by_library_ch,
    DMSANALYSIS_AASEQ.out.aa_seq.map{ it[1] },
    sliding_window_size_ch,
    global_bias_counts_cov_script_ch
    )

    VISUALIZATION_GLOBAL_POS_BIASES_COV(
    variantCounts_filtered_by_library_ch,
    DMSANALYSIS_AASEQ.out.aa_seq.map{ it[1] },
    sliding_window_size_ch,
    aimed_cov_ch,
    global_bias_cov_script_ch
    )

    VISUALIZATION_LOGDIFF(
    library_completed_variantCounts_ch,
    logdiff_script_ch
    )

    if (params.run_seqdepth) {
    VISUALIZATION_SEQDEPTH(
        variantCounts_filtered_by_library_ch,
        DMSANALYSIS_POSSIBLE_MUTATIONS.out.possible_mutations.map{ it[1] },
        min_counts_ch,
        seqdepth_simulation_script_ch
        )
    }

    GATK_GATKTODIMSUM(
    variantCounts_filtered_by_library_ch,
    ch_fasta.map{ it[1] },
    reading_frame_ch,
    gatk_to_dimsum_script_ch
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
    bwa_index      = BWA_INDEX.out.index
    aligned_bam    = BWA_MEM.out.bam
    filtered_bam   = BAMFILTER_DMS.out.bam
    premerged_bam  = PREMERGE.out.bam

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
