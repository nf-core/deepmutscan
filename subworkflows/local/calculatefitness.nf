/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MERGE_COUNTS }             from '../../modules/local/fitness/merge_counts/main'
include { EXPDESIGN_FITNESS }        from '../../modules/local/fitness/fitness_experimental_design/main'
include { FIND_SYNONYMOUS_MUTATION } from '../../modules/local/fitness/find_synonymous_mutation/main'
include { FITNESS_CALCULATION }      from '../../modules/local/fitness/fitness_calculation/main'
include { FITNESS_QC }               from '../../modules/local/fitness/fitness_QC/main'
include { FITNESS_HEATMAP }          from '../../modules/local/fitness/fitness_heatmap/main'
include { RUN_DIMSUM }               from '../../modules/local/fitness/run_dimsum/main'
include { RUN_MUTSCAN }              from '../../modules/local/fitness/run_mutscan/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW DEFINITION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CALCULATEFITNESS {

    take:
    ch_fitness_input    // channel: output from GATK_GATKTOFITNESS (grouped by sample/replicate logic)
    ch_samplesheet_csv  // channel: path to samplesheet csv
    ch_fasta            // channel: tuple([meta], path(fasta))
    val_reading_frame   // value: reading frame param
    ch_aa_seq           // channel: output from DMSANALYSIS_AASEQ (for heatmap)

    main:
    ch_versions = Channel.empty()

    // 1. Group per biological sample and merge counts
    ch_fitness_input
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

    // 2. Merge Counts
    MERGE_COUNTS( ch_fitness_bundled )
    ch_versions = ch_versions.mix(MERGE_COUNTS.out.versions)

    // 3. Experimental Design
    EXPDESIGN_FITNESS( ch_samplesheet_csv )
    ch_versions = ch_versions.mix(EXPDESIGN_FITNESS.out.versions)

    // 4. Find Synonymous Mutations
    // Prepare inputs (broadcast logic)
    def ch_fasta_path = ch_fasta.map { it[1] } // strip meta

    FIND_SYNONYMOUS_MUTATION(
        MERGE_COUNTS.out.merged_counts,
        ch_fasta_path.combine(MERGE_COUNTS.out.merged_counts).map { it[0] },
        val_reading_frame.combine(MERGE_COUNTS.out.merged_counts).map { it[0] }
    )
    ch_versions = ch_versions.mix(FIND_SYNONYMOUS_MUTATION.out.versions)


    // 5. Align Channels for Fitness Calculation

    // Key counts and WT by biological sample name
    def ch_counts_keyed_d = MERGE_COUNTS.out.merged_counts
        .map { smp, counts -> tuple(smp.sample as String, smp, counts) }

    def ch_wt_keyed_d = FIND_SYNONYMOUS_MUTATION.out.synonymous_wt
        .map { smp, wt -> tuple(smp.sample as String, wt) }

    // Join by key
    def ch_counts_wt_d = ch_counts_keyed_d.join(ch_wt_keyed_d)
        .map { key, smp, counts, wt -> tuple(smp, counts, wt) }

    // Broadcast experimental design
    def ch_exp_for_each_d = EXPDESIGN_FITNESS.out.experimental_design
        .combine(ch_counts_wt_d)
        .map { it[0] }

    // Final aligned channels
    def ch_run_counts_d = ch_counts_wt_d.map { smp, counts, wt -> tuple(smp, counts) }
    def ch_run_wt_d     = ch_counts_wt_d.map { smp, counts, wt -> wt }
    def ch_run_exp_d    = ch_exp_for_each_d

    // 6. Run Fitness Calculation & QC
    FITNESS_CALCULATION(
        ch_run_counts_d,
        ch_run_exp_d,
        ch_run_wt_d
    )
    ch_versions = ch_versions.mix(FITNESS_CALCULATION.out.versions)

    FITNESS_QC( FITNESS_CALCULATION.out.fitness_estimation )
    ch_versions = ch_versions.mix(FITNESS_QC.out.versions)

    FITNESS_HEATMAP(
        FITNESS_CALCULATION.out.fitness_estimation,
        ch_aa_seq
    )
    ch_versions = ch_versions.mix(FITNESS_HEATMAP.out.versions)

    // 7. Run DiMSum (optional based on params inside subworkflow or handled by control logic)
    // Note: Logic checking for 'params.dimsum' needs to be handled.
    // Since subworkflows inherit params, we can check params.dimsum here.

    if (params.dimsum) {
        log.warn("""
        '--dimsum true' only works together with '--fitness true'
        and is currently (30 Oct 2025) NOT supported on ARM processors.
        Use AMD/x86_64 systems for DiMSum execution.
        """)

        RUN_DIMSUM(
            ch_run_counts_d,
            ch_run_wt_d,
            ch_run_exp_d
        )
        ch_versions = ch_versions.mix(RUN_DIMSUM.out.versions)
    }

    // 8. Run Mutscan
    if (params.mutscan) {
        RUN_MUTSCAN(
            ch_run_counts_d,
            ch_run_wt_d,
            ch_run_exp_d
        )
        ch_versions = ch_versions.mix(RUN_MUTSCAN.out.versions)
    }

    emit:
    fitness_estimation = FITNESS_CALCULATION.out.fitness_estimation
    versions           = ch_versions
}
