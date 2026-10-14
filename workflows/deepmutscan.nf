/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { BWA_INDEX              } from '../modules/nf-core/bwa/index/main'
include { BWA_MEM                } from '../modules/nf-core/bwa/mem/main'
include { BAMFILTER_DMS          } from '../modules/local/bamprocessing/bam_filter/main'
include { PREMERGE               } from '../modules/local/bamprocessing/premerge/main'
include { SAMTOOLS_SORT              } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX             } from '../modules/nf-core/samtools/index/main'
include { VARIANTCOUNTING            } from '../modules/local/variantcounting/main'
include { DMSANALYSIS_AASEQ      } from '../modules/local/dmsanalysis/aa_seq/main'
include { DMSANALYSIS_POSSIBLE_MUTATIONS      } from '../modules/local/dmsanalysis/possible_mutations/main'
include { DMSANALYSIS_PROCESS_VARIANT_COUNTS      } from '../modules/local/dmsanalysis/process_variant_counts/main'
include { DMSANALYSIS_ERROR_CORRECTION_FALSE_DOUBLES } from '../modules/local/dmsanalysis/error_correction_false_doubles/main'
include { DMSANALYSIS_ERROR_CORRECTION_WILDTYPE      } from '../modules/local/dmsanalysis/error_correction_wildtype/main'
include { VISUALIZATION_COUNTS_PER_COV      } from '../modules/local/visualization/counts_per_cov/main'
include { VISUALIZATION_COUNTS_HEATMAP      } from '../modules/local/visualization/counts_heatmap/main'
include { VISUALIZATION_GLOBAL_POS_BIASES_COUNTS      } from '../modules/local/visualization/global_pos_biases_counts/main'
include { VISUALIZATION_GLOBAL_POS_BIASES_COV      } from '../modules/local/visualization/global_pos_biases_cov/main'
include { VISUALIZATION_LOGDIFF      } from '../modules/local/visualization/logdiff/main'
include { VISUALIZATION_SEQDEPTH      } from '../modules/local/visualization/seqdepth/main'
include { VISUALIZATION_ERROR_CORRECTION_REPORT } from '../modules/local/visualization/error_correction_report/main'
include { VISUALIZATION_SUMMARY_REPORT } from '../modules/local/visualization/summary_report/main'
include { VARIANT_EFFECT_INSPECTION_TOOL } from '../modules/local/structure/variant_effect_inspection_tool/main'
include { COUNTS_TO_FITNESS          } from '../modules/local/fitness/counts_to_fitness/main'

include { CALCULATE_FITNESS } from '../subworkflows/local/calculate_fitness/main'

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_deepmutscan_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DEEPMUTSCAN {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = Channel.empty()
    def ch_multiqc_files = Channel.empty()

    // Define input channels from parameters
    def ch_fasta = Channel
        .fromPath(params.fasta, checkIfExists: true)
        .map { fasta -> tuple( [id: 'ref'], fasta ) }

    def reading_frame_ch        = Channel.value(params.reading_frame)
    def min_counts_ch           = Channel.value(params.min_counts)
    def custom_codon_library_ch = Channel.value(params.custom_codon_library)
    def mutagenesis_type_ch     = Channel.value(params.mutagenesis_type)
    def sliding_window_size_ch  = Channel.value(params.sliding_window_size)
    def aimed_cov_ch            = Channel.value(params.aimed_cov)
    def run_seqdepth_ch         = Channel.value(params.run_seqdepth)
    def base_qual_ch            = Channel.value(params.base_qual)
    def min_flank_ch            = Channel.value(params.min_flank)

    // Raw samplesheet path channel for downstream subworkflows
    def ch_samplesheet_csv      = Channel.fromPath(params.input, checkIfExists: true)

    // '--dimsum' only has an effect together with '--fitness'
    if (params.dimsum && !params.fitness) {
        log.warn("'--dimsum true' only works together with '--fitness true'. DiMSum will be skipped.")
    }

    // '--error_correction wildtype' needs a 'wildtype' sample for every biological sample
    if (params.error_correction == 'wildtype') {
        ch_samplesheet
            .map { meta, _reads -> tuple(meta.sample as String, meta.type as String) }
            .toList()
            .map { pairs ->
                pairs.groupBy { it[0] }.each { sample, rows ->
                    def types = rows.collect { it[1] }
                    if (types.any { it != 'wildtype' } && !types.contains('wildtype')) {
                        error(
                            "[--error_correction wildtype] No 'wildtype' sample found for '${sample}' in the samplesheet.\n" +
                            "  Add a samplesheet row with type 'wildtype' and sample '${sample}' (deep wildtype-only sequencing),\n" +
                            "  or use the internal correction instead (default: --error_correction false_doubles),\n" +
                            "  or disable error correction with --error_correction none."
                        )
                    }
                }
                pairs
            }
            .subscribe { }
    }

    // Nudge fitness users who did not supply a structure towards the interactive tool
    if (params.fitness && !params.pdb) {
        log.info(
            "\n" +
            "  ┌─ Optional: interactive variant effect inspection tool ───────────────────┐\n" +
            "     Supply a wildtype structure with '--pdb <structure.pdb>' to also build an\n" +
            "     interactive HTML tool that projects fitness, counts and error-correction\n" +
            "     biases onto the 3D structure to support data interpretation.\n" +
            "  └───────────────────────────────────────────────────────────────────────────┘\n"
        )
    }

    //
    // MODULE: Run FastQC
    //
    FASTQC(ch_samplesheet)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map{ _meta, file -> file })

    //
    // MODULE: BWA Index
    //
    BWA_INDEX (
        ch_fasta
    )

    // Broadcast index to all samples
    def ch_bwa_index = BWA_INDEX.out.index

    // Broadcast the index to all samples
    def ch_bwa_index_broadcast = ch_samplesheet
      .combine(ch_bwa_index)
      .map { [it[2], it[3]] }

    // Broadcast the fasta to all samples
    def ch_fasta_broadcast = ch_fasta
      .combine(ch_samplesheet)
      .map { [it[0], it[1]] }

    // Broadcast the sort flag to all samples
    def ch_sort_bam = ch_samplesheet.map { false }

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
    def ch_fasta_path_broadcast = ch_fasta
      .combine(BAMFILTER_DMS.out.bam)   // flattened item: [meta3, fasta, meta, bam]
      .map { it[1] }                    // keep only the fasta path (N emissions)

    PREMERGE(
      BAMFILTER_DMS.out.bam,     // tuple(val(meta), path(bam))
      ch_fasta_path_broadcast    // path(fasta)
    )

    // Coordinate-sort + index the premerged BAM (required by the read-based counter)
    SAMTOOLS_SORT(
      PREMERGE.out.bam,                                        // tuple(val(meta), path(bam))
      PREMERGE.out.bam.map { _meta, _bam -> [ [:], [], [] ] }, // empty fasta tuple (BAM sort)
      Channel.value('')                                        // index_format '' -> index separately
    )
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)

    // Join sorted BAM with its index -> tuple(val(meta), path(bam), path(bai))
    def ch_sorted_indexed = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.index)

    // Broadcast singletons to N (one per sample), anchored on the sorted+indexed channel
    def ch_fasta_for_counts = ch_fasta.combine(ch_sorted_indexed).map { it[1] }          // path -- N
    def ch_rf_for_counts    = reading_frame_ch.combine(ch_sorted_indexed).map { it[0] }  // val  -- N
    def ch_min_for_counts   = min_counts_ch.combine(ch_sorted_indexed).map { it[0] }     // val  -- N
    def ch_bq_for_counts    = base_qual_ch.combine(ch_sorted_indexed).map { it[0] }      // val  -- N
    def ch_flank_for_counts = min_flank_ch.combine(ch_sorted_indexed).map { it[0] }      // val  -- N

    VARIANTCOUNTING(
      ch_sorted_indexed,   // tuple(val(meta), path(bam), path(bai))
      ch_fasta_for_counts, // path(fasta)
      ch_rf_for_counts,    // val(reading_frame string, 1-based inclusive)
      ch_min_for_counts,   // val(min_counts)
      ch_bq_for_counts,    // val(base_qual)
      ch_flank_for_counts  // val(min_flank)
    )
    ch_versions = ch_versions.mix(VARIANTCOUNTING.out.versions)

    DMSANALYSIS_AASEQ (
      ch_fasta,
      reading_frame_ch
    )
    ch_versions = ch_versions.mix(DMSANALYSIS_AASEQ.out.versions)

    DMSANALYSIS_POSSIBLE_MUTATIONS(
      ch_fasta,
      reading_frame_ch,           // pos_range (as val)
      mutagenesis_type_ch,        // mutagenesis_type (as val)
      custom_codon_library_ch     // custom_codon_library (as path)
    )
    ch_versions = ch_versions.mix(DMSANALYSIS_POSSIBLE_MUTATIONS.out.versions)

    // Anchor (N items; one per sample)
    def ch_vc = VARIANTCOUNTING.out.variant_counts   // tuple(val(meta), path)

    // Build per-sample inputs using inline combinations (replaces fanout)
    def ch_possible_mut_for_proc = DMSANALYSIS_POSSIBLE_MUTATIONS.out.possible_mutations.map { it[1] }.combine(ch_vc).map { it[0] }
    def ch_aa_seq_for_proc       = DMSANALYSIS_AASEQ.out.aa_seq.map { it[1] }.combine(ch_vc).map { it[0] }
    def ch_min_counts_for_proc   = min_counts_ch.combine(ch_vc).map { it[0] }

    // Call with all inputs aligned (each has N items now)
    DMSANALYSIS_PROCESS_VARIANT_COUNTS(
      ch_vc,                        // tuple(val(meta), path(variantCounts))  -- N
      ch_possible_mut_for_proc,     // path(possible_mutations)               -- N
      ch_aa_seq_for_proc,           // path(aa_seq)                           -- N
      ch_min_counts_for_proc        // val(min_counts)                        -- N
    )

    def annotated_variantCounts_ch           = DMSANALYSIS_PROCESS_VARIANT_COUNTS.out.processed_variantCounts.map { meta, a, b, c, d -> tuple(meta, a) }
    def variantCounts_filtered_by_library_ch = DMSANALYSIS_PROCESS_VARIANT_COUNTS.out.processed_variantCounts.map { meta, a, b, c, d -> tuple(meta, b) }
    def library_completed_variantCounts_ch   = DMSANALYSIS_PROCESS_VARIANT_COUNTS.out.processed_variantCounts.map { meta, a, b, c, d -> tuple(meta, c) }
    def variantCounts_for_heatmaps_ch        = DMSANALYSIS_PROCESS_VARIANT_COUNTS.out.processed_variantCounts.map { meta, a, b, c, d -> tuple(meta, d) }

    //
    // Sequencing-error correction of single-codon counts (default: false_doubles).
    // Replaces the filtered + heatmap channels with corrected versions; downstream fitness
    // picks up correction via the `counts_corrected` column, count heatmaps via canonical names.
    //
    if (params.error_correction == 'false_doubles') {
      // per sample: raw (pre-filter) counts + library-filtered counts + completed library table
      def ch_fd_in       = ch_vc.join(variantCounts_filtered_by_library_ch).join(library_completed_variantCounts_ch)  // (meta, raw, filtered, completed)
      def aa_seq_for_ec  = DMSANALYSIS_AASEQ.out.aa_seq.map { it[1] }.combine(ch_fd_in).map { it[0] }
      def min_for_ec     = min_counts_ch.combine(ch_fd_in).map { it[0] }
      // the reference the reads were aligned to (base_mut is in this fasta's coordinate frame),
      // broadcast one-per-sample; the estimator and the mode/window are run-level constants.
      // ch_fasta is tuple(meta, fasta), so the path is it[1] - it[0] is the [id:ref] meta map.
      def fasta_for_ec   = ch_fasta.combine(ch_fd_in).map { it[1] }
      DMSANALYSIS_ERROR_CORRECTION_FALSE_DOUBLES(
        ch_fd_in, fasta_for_ec, aa_seq_for_ec, min_for_ec,
        params.false_doubles_method, params.false_doubles_codon_window
      )
      ch_versions = ch_versions.mix(DMSANALYSIS_ERROR_CORRECTION_FALSE_DOUBLES.out.versions)
      variantCounts_filtered_by_library_ch = DMSANALYSIS_ERROR_CORRECTION_FALSE_DOUBLES.out.corrected.map { meta, filt, heat, comp -> tuple(meta, filt) }
      variantCounts_for_heatmaps_ch        = DMSANALYSIS_ERROR_CORRECTION_FALSE_DOUBLES.out.corrected.map { meta, filt, heat, comp -> tuple(meta, heat) }
      library_completed_variantCounts_ch   = DMSANALYSIS_ERROR_CORRECTION_FALSE_DOUBLES.out.corrected.map { meta, filt, heat, comp -> tuple(meta, comp) }
    }
    else if (params.error_correction == 'wildtype') {
      // pair each non-wildtype sample with the wildtype sample of the same `sample` base
      def filt_comp = variantCounts_filtered_by_library_ch.join(library_completed_variantCounts_ch)  // (meta, filtered, completed)
      def wt_ch     = filt_comp.filter { meta, f, c -> meta.type == 'wildtype' }.map { meta, f, c -> tuple(meta.sample, f) }
      def tgt_ch    = filt_comp.filter { meta, f, c -> meta.type != 'wildtype' }.map { meta, f, c -> tuple(meta.sample, meta, f, c) }
      def ch_wt_in  = tgt_ch.combine(wt_ch, by: 0).map { _base, meta, filt, comp, wtf -> tuple(meta, filt, wtf, comp) }  // (meta, filtered, wt_filtered, completed)
      def aa_seq_for_ec = DMSANALYSIS_AASEQ.out.aa_seq.map { it[1] }.combine(ch_wt_in).map { it[0] }
      def min_for_ec    = min_counts_ch.combine(ch_wt_in).map { it[0] }
      DMSANALYSIS_ERROR_CORRECTION_WILDTYPE( ch_wt_in, aa_seq_for_ec, min_for_ec )
      ch_versions = ch_versions.mix(DMSANALYSIS_ERROR_CORRECTION_WILDTYPE.out.versions)
      variantCounts_filtered_by_library_ch = DMSANALYSIS_ERROR_CORRECTION_WILDTYPE.out.corrected.map { meta, filt, heat, comp -> tuple(meta, filt) }
      variantCounts_for_heatmaps_ch        = DMSANALYSIS_ERROR_CORRECTION_WILDTYPE.out.corrected.map { meta, filt, heat, comp -> tuple(meta, heat) }
      library_completed_variantCounts_ch   = DMSANALYSIS_ERROR_CORRECTION_WILDTYPE.out.corrected.map { meta, filt, heat, comp -> tuple(meta, comp) }
    }
    // else 'none': keep the uncorrected channels

    // Artefacts the all-in-one report embeds, collected as [[kind, group, name], file]. Optional
    // parts simply never mix anything in, and the report drops the corresponding section.
    def ch_report_assets = Channel.empty()

    //
    // MODULE: per-sample error-correction report (only when correction was applied)
    //
    if (params.error_correction != 'none') {
      // One run-level report over every file. Sorting first keeps the id list and the staged files in
      // the same order, which is the only thing tying a column of numbers to the file it came from.
      def ch_ec = variantCounts_filtered_by_library_ch.toSortedList { a, b -> a[0].id <=> b[0].id }

      VISUALIZATION_ERROR_CORRECTION_REPORT(
        ch_ec.map { rows -> rows[0][0].sample ?: 'run' },
        ch_ec.map { rows -> groovy.json.JsonOutput.toJson(rows.collect { [id: it[0].id, type: it[0].type, replicate: it[0].replicate] }).bytes.encodeBase64().toString() },
        ch_ec.map { rows -> rows.collect { it[1] } },
        params.error_correction,
        params.false_doubles_method,
        file("${projectDir}/assets/error_correction_report/ec_report_template.html", checkIfExists: true)
      )
      ch_versions = ch_versions.mix(VISUALIZATION_ERROR_CORRECTION_REPORT.out.versions)
      ch_report_assets = ch_report_assets.mix(
        VISUALIZATION_ERROR_CORRECTION_REPORT.out.report.map { f -> [['report', 'ec', 'Error correction'], f] }
      )
    }

    // --- For VISUALIZATION_COUNTS_PER_COV & HEATMAP (replaces fanoutTo)
    def min_counts_for_cov_ch          = min_counts_ch.combine(variantCounts_for_heatmaps_ch).map { it[0] }
    def min_counts_for_heatmap_ch      = min_counts_ch.combine(variantCounts_for_heatmaps_ch).map { it[0] }

    // --- For VISUALIZATION_GLOBAL_POS_BIASES_*
    def aa_seq_for_bias_ch             = DMSANALYSIS_AASEQ.out.aa_seq.map { it[1] }.combine(variantCounts_filtered_by_library_ch).map { it[0] }
    def sliding_window_size_N          = sliding_window_size_ch.combine(variantCounts_filtered_by_library_ch).map { it[0] }
    def aimed_cov_N                    = aimed_cov_ch.combine(variantCounts_filtered_by_library_ch).map { it[0] }

    // --- For VISUALIZATION_SEQDEPTH
    def possible_mutations_N           = DMSANALYSIS_POSSIBLE_MUTATIONS.out.possible_mutations.map { it[1] }.combine(variantCounts_filtered_by_library_ch).map { it[0] }
    def min_counts_for_seqdepth_ch     = min_counts_ch.combine(variantCounts_filtered_by_library_ch).map { it[0] }

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

    // library QC plots (original R PDFs) + the tables the report re-plots interactively
    ch_report_assets = ch_report_assets
      .mix(VISUALIZATION_GLOBAL_POS_BIASES_COV.out.rolling_coverage.map          { meta, f -> [['qc', meta.id, 'rolling_coverage.pdf'], f] })
      .mix(VISUALIZATION_GLOBAL_POS_BIASES_COUNTS.out.rolling_counts.map         { meta, f -> [['qc', meta.id, 'rolling_counts.pdf'], f] })
      .mix(VISUALIZATION_GLOBAL_POS_BIASES_COUNTS.out.rolling_counts_per_cov.map { meta, f -> [['qc', meta.id, 'rolling_counts_per_cov.pdf'], f] })
      .mix(VISUALIZATION_COUNTS_HEATMAP.out.counts_heatmap.map                   { meta, f -> [['qc', meta.id, 'counts_heatmap.pdf'], f] })
      .mix(VISUALIZATION_COUNTS_PER_COV.out.counts_per_cov_heatmap.map           { meta, f -> [['qc', meta.id, 'counts_per_cov_heatmap.pdf'], f] })
      .mix(VISUALIZATION_LOGDIFF.out.logdiff_plot.map                            { meta, f -> [['qc', meta.id, 'logdiff_plot.pdf'], f] })
      .mix(VISUALIZATION_LOGDIFF.out.logdiff_varying_bases.map                   { meta, f -> [['qc', meta.id, 'logdiff_varying_bases.pdf'], f] })
      .mix(variantCounts_filtered_by_library_ch.map                              { meta, f -> [['filtered', meta.id, 'counts.csv'], f] })
      .mix(variantCounts_for_heatmaps_ch.map                                     { meta, f -> [['heatmap', meta.id, 'heatmap.csv'], f] })
      // Fixes the ORF axis at the wildtype length: positions with no library coverage still have to
      // occupy a slot, or every rolling window downstream is computed over the wrong neighbourhood.
      .mix(DMSANALYSIS_AASEQ.out.aa_seq.map                                      { _meta, f -> [['aa_seq', 'run', 'aa_seq.txt'], f] })
      // Needed to name the exact nucleotide change behind each heatmap cell: the fitness table gives
      // the variant ORFs, and the wildtype to diff them against only exists in the reference.
      .mix(ch_fasta.map                                                          { _meta, f -> [['fasta', 'run', 'reference.fasta'], f] })

    if (params.run_seqdepth) {
      VISUALIZATION_SEQDEPTH(
        variantCounts_filtered_by_library_ch,
        possible_mutations_N,
        min_counts_for_seqdepth_ch
      )
      ch_versions = ch_versions.mix(VISUALIZATION_SEQDEPTH.out.versions)
      // The rarefaction curve travels to the summary report, which re-plots it per file under
      // Sequencing QC (the PDF is still written to the results folder for the record).
      ch_report_assets = ch_report_assets.mix(
        VISUALIZATION_SEQDEPTH.out.curve.map { meta, f -> [['seqdepth', meta.id, 'seqdepth.csv'], f] }
      )
    }

    // Broadcast singletons to N (one per sample), anchored on variantCounts_filtered_by_library_ch
    def ch_fasta_for_fitness    = ch_fasta.combine(variantCounts_filtered_by_library_ch).map { it[1] }      // path(fasta) -- N
    def ch_rf_for_fitness       = reading_frame_ch.combine(variantCounts_filtered_by_library_ch).map { it[0] }  // val(range) -- N

    // Call with aligned inputs
    COUNTS_TO_FITNESS(
      variantCounts_filtered_by_library_ch, // tuple(val(meta), path)
      ch_fasta_for_fitness,         // path(fasta)
      ch_rf_for_fitness             // val(reading_frame)
    )

    // Execution of fitness subworkflow, if --fitness true
    if (params.fitness) {

        CALCULATE_FITNESS (
            COUNTS_TO_FITNESS.out.fitness_input, // Input from previous step
            ch_samplesheet_csv,                   // Path to samplesheet
            ch_fasta,                             // The original Fasta tuple
            reading_frame_ch,                     // Reading frame value channel
            DMSANALYSIS_AASEQ.out.aa_seq          // Amino Acid Sequence (for Heatmap)
        )

        // Collect versions
        ch_versions = ch_versions.mix(CALCULATE_FITNESS.out.versions)

        ch_report_assets = ch_report_assets
          .mix(CALCULATE_FITNESS.out.fitness_estimation.first().map { _m, f -> [['fitness_table', 'default', 'fitness_estimation.tsv'], f] })
          .mix(CALCULATE_FITNESS.out.fitness_plots.map { _m, f -> [['fitness_pdf', 'default', f.name], f] })
          // mutscan emits all of its plots as one list, so flatten to a file per entry
          .mix(CALCULATE_FITNESS.out.mutscan_plots.transpose().map { _m, f -> [['mutscan_pdf', 'mutscan', f.name], f] })
          .mix(CALCULATE_FITNESS.out.dimsum_dir.flatten().filter { it.name == 'report.html' }.map { f -> [['report', 'dimsum', 'DiMSum'], f] })

        //
        // MODULE: interactive variant effect inspection tool (one self-contained HTML per sample).
        // Runs only when a wildtype structure is supplied via --pdb. Projects fitness, counts and
        // error-correction biases onto the 3D structure. (Structure prediction is not wired here.)
        //
        if (params.pdb) {
            def ch_structure = Channel
                .fromPath(params.pdb, checkIfExists: true)
                .map { pdb -> tuple([id: 'wildtype'], pdb) }
            def ch_viewer_input = CALCULATE_FITNESS.out.fitness_estimation
                .combine(ch_structure)
                .map { smeta, fitness_tsv, _stmeta, pdb -> tuple(smeta, pdb, fitness_tsv) }
            VARIANT_EFFECT_INSPECTION_TOOL(
                ch_viewer_input,
                variantCounts_filtered_by_library_ch.map { _meta, file -> file }.collect(),
                DMSANALYSIS_AASEQ.out.aa_seq.map { _meta, file -> file }.first(),
                file("${projectDir}/assets/variant_effect_inspection_tool/3Dmol-min.js",         checkIfExists: true),
                file("${projectDir}/assets/variant_effect_inspection_tool/viewer_template.html", checkIfExists: true),
            )
            ch_versions = ch_versions.mix(VARIANT_EFFECT_INSPECTION_TOOL.out.versions)
            ch_report_assets = ch_report_assets.mix(
                VARIANT_EFFECT_INSPECTION_TOOL.out.viewer.map { _meta, f -> [['report', 'viewer', 'Variant effect inspection tool'], f] }
            )
            // The raw structure also travels to the summary report, which builds its own slow-spinning
            // fitness-coloured viewer on the Variant effects page from it.
            ch_report_assets = ch_report_assets.mix(
                ch_structure.map { _meta, pdb -> [['pdb', 'structure', 'structure.pdb'], pdb] }
            )
        }
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'deepmutscan_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'deepmutscan'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    //
    // MODULE: the single all-in-one report. Runs last, because it embeds every other artefact -
    // the MultiQC report included - as a base64 data URI, so the file stands alone when shared.
    //
    ch_report_assets = ch_report_assets.mix(
        MULTIQC.out.report.map { _meta, f -> [['report', 'multiqc', 'MultiQC'], f] }
    )
    // The collated software versions drive the Citation page: it cites only the tools this run used.
    ch_report_assets = ch_report_assets.mix(
        ch_collated_versions.map { f -> [['versions', 'run', 'versions.yml'], f] }
    )
    // Nextflow's own execution trace drives the Overview run-statistics table. It is still being
    // appended to as the report is built, so it captures every compute-heavy task (all but the report
    // task itself) - the page labels it as of report generation. checkIfExists is off because the file
    // is engine-managed, not a process output.
    ch_report_assets = ch_report_assets.mix(
        Channel.fromPath("${outdir}/pipeline_info/execution_trace_${params.trace_report_suffix}.txt", checkIfExists: false)
            .map { f -> [['trace', 'run', 'execution_trace.txt'], f] }
    )

    // Sorting keeps the manifest and the staged files in one deterministic, matching order.
    def ch_summary_input = ch_report_assets
        .toSortedList { a, b -> a[0].join('|') <=> b[0].join('|') }
        .map { rows ->
            def manifest = rows.collect { [kind: it[0][0], group: it[0][1], name: it[0][2]] }
            tuple(
                groovy.json.JsonOutput.toJson(manifest).bytes.encodeBase64().toString(),
                rows.collect { it[1] },
            )
        }

    def ch_run_meta = variantCounts_filtered_by_library_ch
        .first()
        .map { meta, _f ->
            tuple(
                [id: 'run', sample: meta.sample],
                groovy.json.JsonOutput.toJson([
                    sample: meta.sample,
                    reading_frame: params.reading_frame,
                    mutagenesis_type: params.mutagenesis_type,
                    error_correction: params.error_correction,
                    version: workflow.manifest.version,
                    // The report smooths client-side with the same window the R plots use, and draws
                    // the same required-coverage reference line, so both must travel with the run.
                    sliding_window_size: params.sliding_window_size,
                    aimed_cov: params.aimed_cov,
                    fitness: params.fitness,
                    dimsum: params.dimsum,
                    mutscan: params.mutscan,
                    outdir: params.outdir,
                    // Lets the Overview link straight to Nextflow's own execution report + timeline,
                    // which sit next to this file as execution_{report,timeline}_<suffix>.html.
                    trace_report_suffix: params.trace_report_suffix,
                    false_doubles_method: params.false_doubles_method,
                ]).bytes.encodeBase64().toString(),
            )
        }

    VISUALIZATION_SUMMARY_REPORT(
        ch_run_meta,
        ch_summary_input.map { manifest, _files -> manifest },
        ch_summary_input.map { _manifest, files -> files },
        file("${projectDir}/assets/summary_report/summary_report_template.html", checkIfExists: true),
        file("${projectDir}/assets/variant_effect_inspection_tool/3Dmol-min.js", checkIfExists: true),
    )
    ch_versions = ch_versions.mix(VISUALIZATION_SUMMARY_REPORT.out.versions)

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
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
