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
include { GATK_GATKTOFITNESS          } from '../modules/local/gatk/gatktofitness'

include { CALCULATEFITNESS } from '../subworkflows/local/calculatefitness'

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_deepmutscan_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Define input variables as channels
Channel
    .fromPath(params.fasta, checkIfExists: true)
    .map { fasta -> tuple( [id: 'ref'], fasta ) }
    .set { ch_fasta }
Channel
    .value(params.reading_frame)
    .set { reading_frame_ch }
Channel
    .value(params.min_counts)
    .set { min_counts_ch }
Channel
    .value(params.custom_codon_library)
    .set { custom_codon_library_ch }
Channel
    .value(params.mutagenesis_type)
    .set { mutagenesis_type_ch }
Channel
    .value(params.sliding_window_size)
    .set { sliding_window_size_ch }
Channel
    .value(params.aimed_cov)
    .set { aimed_cov_ch }
Channel
    .value(params.run_seqdepth)
    .set { run_seqdepth_ch }
Channel
  .fromPath(params.input, checkIfExists: true)
  .set { ch_samplesheet_csv }


workflow DEEPMUTSCAN {

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
            name: 'nf_core_'  +  'deepmutscan_software_'  + 'mqc_'  + 'versions.yml',
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

    // Broadcast the index to all samples
    ch_bwa_index_broadcast = ch_samplesheet
      .combine(ch_bwa_index)
      .map { [it[2], it[3]] }

    // Broadcast the fasta to all samples
    ch_fasta_broadcast = ch_fasta
      .combine(ch_samplesheet)
      .map { [it[0], it[1]] }

    // Broadcast the sort flag to all samples
    ch_sort_bam = ch_samplesheet.map { false }

    // Run BWA_MEM with all four inputs aligned
    BWA_MEM(
      ch_samplesheet,
      ch_bwa_index_broadcast,
      ch_fasta_broadcast,
      ch_sort_bam
    )

    BAMFILTER_DMS (
    BWA_MEM.out.bam
    )

    // Broadcast the FASTA path to every BAM emitted by BAMFILTER_DMS
    ch_fasta_path_broadcast = ch_fasta
      .combine(BAMFILTER_DMS.out.bam)   // flattened item: [meta3, fasta, meta, bam]
      .map { it[1] }                    // keep only the fasta path (N emissions)

    PREMERGE(
      BAMFILTER_DMS.out.bam,     // tuple(val(meta), path(bam))
      ch_fasta_path_broadcast    // path(fasta)
    )

    // FASTA path for GATK: broadcast to N
    ch_fasta_for_gatk  = ch_fasta.combine(PREMERGE.out.bam).map { it[1] }		// path -- N
    // Reading frame for GATK: broadcast to N (it's a val string)
    ch_rf_for_gatk     = reading_frame_ch.combine(PREMERGE.out.bam).map { it[0] }	// val  -- N
    // min_counts for GATK: broadcast to N (also a val)
    ch_min_for_gatk    = min_counts_ch.combine(PREMERGE.out.bam).map { it[0] }		// val  -- N

    GATK_SATURATIONMUTAGENESIS(
      PREMERGE.out.bam,   // merged reads - tuple(val(meta), path(bam))
      ch_fasta_for_gatk,  // path(fasta)
      ch_rf_for_gatk,     // val(reading_frame string)
      ch_min_for_gatk     // val(min_counts)
    )

    DMSANALYSIS_AASEQ (
    ch_fasta,
    reading_frame_ch
    )
    ch_versions = ch_versions.mix(DMSANALYSIS_AASEQ.out.versions)

    DMSANALYSIS_POSSIBLE_MUTATIONS(
    ch_fasta,
    reading_frame_ch,			// pos_range (as val)
    mutagenesis_type_ch,		// mutagenesis_type (as val)
    custom_codon_library_ch		// custom_codon_library (as path)
    )
    ch_versions = ch_versions.mix(DMSANALYSIS_POSSIBLE_MUTATIONS.out.versions)

    // Anchor (N items; one per sample)
    def ch_vc = GATK_SATURATIONMUTAGENESIS.out.variantCounts   // tuple(val(meta), path)

    // Fan-out helpers (broadcast singleton → N)
    def fanout = { ch_singleton -> ch_singleton.combine(ch_vc).map { it[0] } }

    // Build per-sample inputs
    ch_possible_mut_for_proc = fanout( DMSANALYSIS_POSSIBLE_MUTATIONS.out.possible_mutations.map { it[1] } )
    ch_aa_seq_for_proc       = fanout( DMSANALYSIS_AASEQ.out.aa_seq.map { it[1] } )
    ch_min_counts_for_proc   = fanout( min_counts_ch )

    // Call with all inputs aligned (each has N items now)
    DMSANALYSIS_PROCESS_GATK(
      ch_vc,                        // tuple(val(meta), path(variantCounts))  -- N
      ch_possible_mut_for_proc,     // path(possible_mutations)               -- N
      ch_aa_seq_for_proc,           // path(aa_seq)                           -- N
      ch_min_counts_for_proc        // val(min_counts)                        -- N
    )    

    annotated_variantCounts_ch = DMSANALYSIS_PROCESS_GATK.out.processed_variantCounts.map { meta, a, b, c, d -> tuple(meta, a) }
    variantCounts_filtered_by_library_ch = DMSANALYSIS_PROCESS_GATK.out.processed_variantCounts.map { meta, a, b, c, d -> tuple(meta, b) }
    library_completed_variantCounts_ch = DMSANALYSIS_PROCESS_GATK.out.processed_variantCounts.map { meta, a, b, c, d -> tuple(meta, c) }
    variantCounts_for_heatmaps_ch = DMSANALYSIS_PROCESS_GATK.out.processed_variantCounts.map { meta, a, b, c, d -> tuple(meta, d) }

    // Broadcast `singleton` so it emits once per item in `anchorN`
    def fanoutTo = { anchorN, singleton -> singleton.combine(anchorN).map { it[0] } }

    // --- For VISUALIZATION_COUNTS_PER_COV (anchor: variantCounts_for_heatmaps_ch)
    min_counts_for_cov_ch          = fanoutTo(variantCounts_for_heatmaps_ch, min_counts_ch)

    // --- For VISUALIZATION_COUNTS_HEATMAP (anchor: variantCounts_for_heatmaps_ch)
    min_counts_for_heatmap_ch      = fanoutTo(variantCounts_for_heatmaps_ch, min_counts_ch)

    // --- For VISUALIZATION_GLOBAL_POS_BIASES_* (anchor: variantCounts_filtered_by_library_ch)
    aa_seq_for_bias_ch             = fanoutTo(variantCounts_filtered_by_library_ch, DMSANALYSIS_AASEQ.out.aa_seq.map { it[1] })
    sliding_window_size_N          = fanoutTo(variantCounts_filtered_by_library_ch, sliding_window_size_ch)
    aimed_cov_N                    = fanoutTo(variantCounts_filtered_by_library_ch, aimed_cov_ch)

    // --- For VISUALIZATION_SEQDEPTH (anchor: variantCounts_filtered_by_library_ch)
    possible_mutations_N           = fanoutTo(variantCounts_filtered_by_library_ch, DMSANALYSIS_POSSIBLE_MUTATIONS.out.possible_mutations.map { it[1] })
    min_counts_for_seqdepth_ch     = fanoutTo(variantCounts_filtered_by_library_ch, min_counts_ch)

    VISUALIZATION_COUNTS_PER_COV(
      variantCounts_for_heatmaps_ch,
      min_counts_for_cov_ch
    )

    VISUALIZATION_COUNTS_HEATMAP(
      variantCounts_for_heatmaps_ch,
      min_counts_for_heatmap_ch
    )

    VISUALIZATION_GLOBAL_POS_BIASES_COUNTS(
      variantCounts_filtered_by_library_ch,
      aa_seq_for_bias_ch,
      sliding_window_size_N
    )

    VISUALIZATION_GLOBAL_POS_BIASES_COV(
      variantCounts_filtered_by_library_ch,
      aa_seq_for_bias_ch,
      sliding_window_size_N,
      aimed_cov_N
    )

    VISUALIZATION_LOGDIFF(
      library_completed_variantCounts_ch
    )

    if (params.run_seqdepth) {
      VISUALIZATION_SEQDEPTH(
        variantCounts_filtered_by_library_ch,
        possible_mutations_N,
        min_counts_for_seqdepth_ch
      )
    }

    // Broadcast singletons to N (one per sample), anchored on variantCounts_filtered_by_library_ch
    ch_fasta_for_fitness    = ch_fasta.combine(variantCounts_filtered_by_library_ch).map { it[1] }		// path(fasta) -- N
    ch_rf_for_fitness       = reading_frame_ch.combine(variantCounts_filtered_by_library_ch).map { it[0] }	// val(range) -- N

    // Call with aligned inputs
    GATK_GATKTOFITNESS(
      variantCounts_filtered_by_library_ch,	// tuple(val(meta), path)
      ch_fasta_for_fitness,			// path(fasta)
      ch_rf_for_fitness 			// val(reading_frame)
    )

    // Execution of fitness subworkflow, if --fitness true
    if (params.fitness) {
        
        CALCULATEFITNESS (
            GATK_GATKTOFITNESS.out.fitness_input, // Input vom vorherigen Schritt
            ch_samplesheet_csv,                   // Pfad zum Samplesheet
            ch_fasta,                             // Das originale Fasta Tuple
            reading_frame_ch,                     // Reading frame value channel
            DMSANALYSIS_AASEQ.out.aa_seq          // Amino Acid Sequence (für Heatmap)
        )

        // Versionen einsammeln
        ch_versions = ch_versions.mix(CALCULATEFITNESS.out.versions)
    }



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
