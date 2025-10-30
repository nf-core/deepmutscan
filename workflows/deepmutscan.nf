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
include { MERGE_COUNTS               } from '../modules/local/fitness/merge_counts'
include { EXPDESIGN_FITNESS               } from '../modules/local/fitness/fitness_experimental_design'
include { FIND_SYNONYMOUS_MUTATION } from '../modules/local/fitness/find_synonymous_mutation'
include { RUN_DIMSUM } from '../modules/local/dimsum/run_dimsum'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_deepmutscan_pipeline'

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
params.fitness       = params.fitness       ?: false				// run basic fitness calculation from selection input & output samples
params.dimsum       = params.dimsum       ?: false				// run DiMSum for fitness/functionality scores from selection input & output samples


// Define fasta file as channel (e.g. for BWA index)
Channel
    .fromPath(params.fasta, checkIfExists: true)
    .map { fasta -> tuple( [id: 'ref'], fasta ) }
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

// The deepmutscan samplesheet file provided via --input
Channel
  .fromPath(params.input, checkIfExists: true)
  .set { ch_samplesheet_csv }


// Define R scripts as channels
// Helper to anchor to the repo root
def R = { rel -> Channel.fromPath("${workflow.projectDir}/${rel}", checkIfExists: true) }

// Define R scripts as channels (repo-anchored, works from any launch dir)
R("modules/local/dmsanalysis/bin/SeqDepth_simulation.R").set { seqdepth_simulation_script_ch }
R("modules/local/dmsanalysis/bin/aa_seq.R").set { aa_seq_script_ch }
R("modules/local/dmsanalysis/bin/complete_prefiltered_gatk.R").set { complete_gatk_script_ch }
R("modules/local/dmsanalysis/bin/counts_heatmap.R").set { counts_heatmap_script_ch }
R("modules/local/dmsanalysis/bin/counts_per_cov_heatmap.R").set { counts_per_cov_heatmap_script_ch }
R("modules/local/dmsanalysis/bin/detect_codons.R").set { detect_codons_script_ch }
R("modules/local/dmsanalysis/bin/filter_gatk_by_codon_library.R").set { filter_by_library_script_ch }
R("modules/local/dmsanalysis/bin/fitness_heatmap.R").set { fitness_heatmap_script_ch }
R("modules/local/dmsanalysis/bin/gatk_to_fitness.R").set { gatk_to_fitness_script_ch }
R("modules/local/dmsanalysis/bin/global_position_biases_counts_and_counts_per_cov.R").set { global_bias_counts_cov_script_ch }
R("modules/local/dmsanalysis/bin/global_position_biases_cov.R").set { global_bias_cov_script_ch }
R("modules/local/dmsanalysis/bin/install_packages.R").set { install_packages_script_ch }
R("modules/local/dmsanalysis/bin/logdiff.R").set { logdiff_script_ch }
R("modules/local/dmsanalysis/bin/low_count_variants.R").set { low_count_variants_script_ch }
R("modules/local/dmsanalysis/bin/possible_mutations.R").set { possible_mutations_script_ch }
R("modules/local/dmsanalysis/bin/prepare_gatk_data_for_count_heatmaps.R").set { prepare_counts_heatmap_script_ch }
R("modules/local/dmsanalysis/bin/prepare_gatk_data_for_fitness_heatmap.R").set { prepare_fitness_heatmap_script_ch }
R("modules/local/dmsanalysis/bin/process_raw_gatk.R").set { process_raw_gatk_script_ch }
R("modules/local/dmsanalysis/bin/merge_counts.R").set { merge_counts_script_ch }
R("modules/local/dmsanalysis/bin/dimsum_experimentalDesign.R").set { exp_design_ch }
R("modules/local/dmsanalysis/bin/find_syn_mutation.R").set { syn_mut_ch }



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
    ch_fasta_for_gatk  = ch_fasta.combine(PREMERGE.out.bam).map { it[1] }     // path -- N
    // Reading frame for GATK: broadcast to N (it's a val string)
    ch_rf_for_gatk     = reading_frame_ch.combine(PREMERGE.out.bam).map { it[0] }  // val  -- N
    // min_counts for GATK: broadcast to N (also a val)
    ch_min_for_gatk    = min_counts_ch.combine(PREMERGE.out.bam).map { it[0] }     // val  -- N

    GATK_SATURATIONMUTAGENESIS(
      PREMERGE.out.bam,   // merged reads - tuple(val(meta), path(bam))
      ch_fasta_for_gatk,  // path(fasta)
      ch_rf_for_gatk,     // val(reading_frame string)
      ch_min_for_gatk     // val(min_counts)
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

    // Anchor (N items; one per sample)
    def ch_vc = GATK_SATURATIONMUTAGENESIS.out.variantCounts   // tuple(val(meta), path)

    // Fan-out helpers (broadcast singleton → N)
    def fanout = { ch_singleton -> ch_singleton.combine(ch_vc).map { it[0] } }

    // Build per-sample inputs
    ch_possible_mut_for_proc = fanout( DMSANALYSIS_POSSIBLE_MUTATIONS.out.possible_mutations.map { it[1] } )
    ch_aa_seq_for_proc       = fanout( DMSANALYSIS_AASEQ.out.aa_seq.map { it[1] } )
    ch_min_counts_for_proc   = fanout( min_counts_ch )
    ch_proc_raw_script       = fanout( process_raw_gatk_script_ch )
    ch_filter_lib_script     = fanout( filter_by_library_script_ch )
    ch_complete_script       = fanout( complete_gatk_script_ch )
    ch_prepare_heatmap_script= fanout( prepare_counts_heatmap_script_ch )

    // Call with all inputs aligned (each has N items now)
    DMSANALYSIS_PROCESS_GATK(
      ch_vc,                        // tuple(val(meta), path(variantCounts))  -- N
      ch_possible_mut_for_proc,     // path(possible_mutations)               -- N
      ch_aa_seq_for_proc,           // path(aa_seq)                           -- N
      ch_min_counts_for_proc,       // val(min_counts)                        -- N
      ch_proc_raw_script,           // path(R script)                         -- N
      ch_filter_lib_script,         // path(R script)                         -- N
      ch_complete_script,           // path(R script)                         -- N
      ch_prepare_heatmap_script     // path(R script)                         -- N
    )    

    annotated_variantCounts_ch = DMSANALYSIS_PROCESS_GATK.out.processed_variantCounts.map { meta, a, b, c, d -> tuple(meta, a) }
    variantCounts_filtered_by_library_ch = DMSANALYSIS_PROCESS_GATK.out.processed_variantCounts.map { meta, a, b, c, d -> tuple(meta, b) }
    library_completed_variantCounts_ch = DMSANALYSIS_PROCESS_GATK.out.processed_variantCounts.map { meta, a, b, c, d -> tuple(meta, c) }
    variantCounts_for_heatmaps_ch = DMSANALYSIS_PROCESS_GATK.out.processed_variantCounts.map { meta, a, b, c, d -> tuple(meta, d) }

    // Broadcast `singleton` so it emits once per item in `anchorN`
    def fanoutTo = { anchorN, singleton -> singleton.combine(anchorN).map { it[0] } }

    // --- For VISUALIZATION_COUNTS_PER_COV (anchor: variantCounts_for_heatmaps_ch)
    min_counts_for_cov_ch          = fanoutTo(variantCounts_for_heatmaps_ch, min_counts_ch)
    counts_per_cov_heatmap_scriptN = fanoutTo(variantCounts_for_heatmaps_ch, counts_per_cov_heatmap_script_ch)

    // --- For VISUALIZATION_COUNTS_HEATMAP (anchor: variantCounts_for_heatmaps_ch)
    min_counts_for_heatmap_ch      = fanoutTo(variantCounts_for_heatmaps_ch, min_counts_ch)
    counts_heatmap_scriptN         = fanoutTo(variantCounts_for_heatmaps_ch, counts_heatmap_script_ch)

    // --- For VISUALIZATION_GLOBAL_POS_BIASES_* (anchor: variantCounts_filtered_by_library_ch)
    aa_seq_for_bias_ch             = fanoutTo(variantCounts_filtered_by_library_ch, DMSANALYSIS_AASEQ.out.aa_seq.map { it[1] })
    sliding_window_size_N          = fanoutTo(variantCounts_filtered_by_library_ch, sliding_window_size_ch)
    aimed_cov_N                    = fanoutTo(variantCounts_filtered_by_library_ch, aimed_cov_ch)
    global_bias_counts_cov_scriptN = fanoutTo(variantCounts_filtered_by_library_ch, global_bias_counts_cov_script_ch)
    global_bias_cov_scriptN        = fanoutTo(variantCounts_filtered_by_library_ch, global_bias_cov_script_ch)

    // --- For VISUALIZATION_LOGDIFF (anchor: library_completed_variantCounts_ch)
logdiff_scriptN                = fanoutTo(library_completed_variantCounts_ch, logdiff_script_ch)

    // --- For VISUALIZATION_SEQDEPTH (anchor: variantCounts_filtered_by_library_ch)
    possible_mutations_N           = fanoutTo(variantCounts_filtered_by_library_ch, DMSANALYSIS_POSSIBLE_MUTATIONS.out.possible_mutations.map { it[1] })
    min_counts_for_seqdepth_ch     = fanoutTo(variantCounts_filtered_by_library_ch, min_counts_ch)
    seqdepth_simulation_scriptN    = fanoutTo(variantCounts_filtered_by_library_ch, seqdepth_simulation_script_ch)

    VISUALIZATION_COUNTS_PER_COV(
      variantCounts_for_heatmaps_ch,
      min_counts_for_cov_ch,
      counts_per_cov_heatmap_scriptN
    )

    VISUALIZATION_COUNTS_HEATMAP(
      variantCounts_for_heatmaps_ch,
      min_counts_for_heatmap_ch,
      counts_heatmap_scriptN
    )

    VISUALIZATION_GLOBAL_POS_BIASES_COUNTS(
      variantCounts_filtered_by_library_ch,
      aa_seq_for_bias_ch,
      sliding_window_size_N,
      global_bias_counts_cov_scriptN
    )

    VISUALIZATION_GLOBAL_POS_BIASES_COV(
      variantCounts_filtered_by_library_ch,
      aa_seq_for_bias_ch,
      sliding_window_size_N,
      aimed_cov_N,
      global_bias_cov_scriptN
    )

    VISUALIZATION_LOGDIFF(
      library_completed_variantCounts_ch,
      logdiff_scriptN
    )

    if (params.run_seqdepth) {
      VISUALIZATION_SEQDEPTH(
        variantCounts_filtered_by_library_ch,
        possible_mutations_N,
        min_counts_for_seqdepth_ch,
        seqdepth_simulation_scriptN
      )
    }

    // Broadcast singletons to N (one per sample), anchored on variantCounts_filtered_by_library_ch
    ch_fasta_for_fitness    = ch_fasta.combine(variantCounts_filtered_by_library_ch).map { it[1] }   // path(fasta) -- N
    ch_rf_for_fitness       = reading_frame_ch.combine(variantCounts_filtered_by_library_ch).map { it[0] }  // val(range) -- N
    ch_script_for_fitness   = gatk_to_fitness_script_ch.combine(variantCounts_filtered_by_library_ch).map { it[0] } // path(script) -- N

    // Call with aligned inputs
    GATK_GATKTOFITNESS(
      variantCounts_filtered_by_library_ch,  // tuple(val(meta), path)
      ch_fasta_for_fitness,                  // path(fasta)
      ch_rf_for_fitness,                     // val(reading_frame)
      ch_script_for_fitness                  // path(R script)
    )



    // ----- DiMSum: group per biological sample (from samplesheet) and merge counts to use for DiMSum input -----

    GATK_GATKTOFITNESS.out.fitness_input
      .map { meta, tsv ->
          def s  = meta.sample as String
          def id = meta.id     as String
          def base = s ? (s.replaceFirst(/_(input|output|quality)\d+$/, ''))
                   : (id?.tokenize('_')?.first())
          tuple(base as String, tuple(meta, tsv))
      }
      .groupTuple()
      .map { base, pairs ->
          def metas   = pairs.collect { it[0] }
          def inputs  = pairs.findAll { it[0].type == 'input'  }.sort { it[0].replicate }.collect { it[1] }
          def outputs = pairs.findAll { it[0].type == 'output' }.sort { it[0].replicate }.collect { it[1] }
          tuple([sample: base], metas, inputs, outputs)
      }
      .filter { smeta, metas, ins, outs -> ins && outs }
      .set { ch_fitness_bundled }

    // Broadcast the singleton script path to match each bundle
    def ch_merge_script_for_each = merge_counts_script_ch
      .combine(ch_fitness_bundled)
      .map { it[0] }   // keep the script path, one per bundle

    // Launch the merge of counts in DiMSum input format
    if (params.fitness) {
      MERGE_COUNTS(
        ch_fitness_bundled,         // tuple val(sample), val(metas), path(input_counts), path(output_counts)
        ch_merge_script_for_each   // path merge_script (broadcast)
      )
    }
    
    // Create experimental design file to use for DiMSum
    if (params.fitness) {
      EXPDESIGN_FITNESS(
        ch_samplesheet_csv,   // path to CSV
        exp_design_ch         // path to R script
      )
    }


// --- Synonymous WT selection (runs only when --dimsum) ---
// Strip meta once: keep only the fasta path
ch_fasta.map { it[1] }.set { ch_fasta_path }   // path(/…/GID1A.fasta)

if (params.fitness) {
  // MERGE_COUNTS.out.merged_counts shape: tuple( val([sample:'GID1A']), path("counts_merged.tsv") )
  FIND_SYNONYMOUS_MUTATION(
    MERGE_COUNTS.out.merged_counts,                                      // tuple(val(sample), path counts_merged.tsv)
    ch_fasta_path.combine(MERGE_COUNTS.out.merged_counts).map { it[0] }, // path wt_fasta (broadcast to N)
    reading_frame_ch.combine(MERGE_COUNTS.out.merged_counts).map { it[0] }, // val pos_range (broadcast)
    syn_mut_ch.combine(MERGE_COUNTS.out.merged_counts).map { it[0] }     // path R script (broadcast)
  )
}


// --- fitness calculation (only when --fitness true) ---
if (params.fitness) {
  // Shapes:
  // MERGE_COUNTS.out.merged_counts              -> tuple(val([sample:'GID1A']), path('counts_merged.tsv'))
  // FIND_SYNONYMOUS_MUTATION.out.synonymous_wt -> tuple(val([sample:'GID1A']), path('synonymous_wt.txt'))
  // EXPDESIGN_FITNESS.out.experimental_design   -> path('experimentalDesign.tsv')   (singleton)

  // 1) Key counts and WT by biological sample name to align them robustly
  def ch_counts_keyed_d = MERGE_COUNTS.out.merged_counts
      .map { smp, counts -> tuple(smp.sample as String, smp, counts) }

  def ch_wt_keyed_d = FIND_SYNONYMOUS_MUTATION.out.synonymous_wt
      .map { smp, wt -> tuple(smp.sample as String, wt) }

  // 2) Join by key -> (val(sample), path(counts), path(wt))
  def ch_counts_wt_d = ch_counts_keyed_d.join(ch_wt_keyed_d)
      .map { key, smp, counts, wt -> tuple(smp, counts, wt) }

  // 3) Broadcast experimental design (singleton) to each sample triple
  def ch_exp_for_each_d = EXPDESIGN_FITNESS.out.experimental_design
      .combine(ch_counts_wt_d)
      .map { it[0] }

  // 4) Final aligned channels for FITNESS and RUN_DIMSUM
  def ch_run_counts_d = ch_counts_wt_d.map { smp, counts, wt -> tuple(smp, counts) }  // matches: tuple val(sample), path(counts_merged)
  def ch_run_wt_d     = ch_counts_wt_d.map { smp, counts, wt -> wt }                  // matches: path(wt_txt)
  def ch_run_exp_d    = ch_exp_for_each_d                                             // matches: path(exp_design)
}

// --- DiMSum execution (only when --dimsum) ---
if (params.dimsum) {

  RUN_DIMSUM(
    ch_run_counts_d,   // tuple val(sample), path(counts_merged)
    ch_run_wt_d,       // path wt_txt
    ch_run_exp_d       // path experimentalDesign.tsv
  )
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
