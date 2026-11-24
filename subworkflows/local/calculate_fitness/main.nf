/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FITNESS_ESTIMATION } from '../../../modules/local/fitness/fitness_estimation/main'
include { RUN_DIMSUM }         from '../../../modules/local/fitness/run_dimsum/main'
include { RUN_MUTSCAN }        from '../../../modules/local/fitness/run_mutscan/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW DEFINITION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CALCULATE_FITNESS {

    take:
    ch_fitness_input    // channel: output from COUNTS_TO_FITNESS (grouped by sample/replicate logic)
    ch_samplesheet_csv  // value channel: path to samplesheet csv
    ch_fasta            // value channel: tuple([meta], path(fasta))
    val_reading_frame   // value channel: reading frame param
    ch_aa_seq           // value channel: tuple([meta], path(aa_seq)) from DMSANALYSIS_REFERENCE_TABLES (for heatmap)
    val_run_dimsum      // value: boolean, run DiMSum
    val_run_mutscan     // value: boolean, run mutscan

    main:
    ch_versions = Channel.empty()

    // 1. Group per biological sample
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
            // sort by replicate so the staged file order is deterministic (cacheable)
            def inputs  = pairs.findAll { it[0].type == 'input'  }.sort { it[0].replicate }.collect { it[1] }
            def outputs = pairs.findAll { it[0].type == 'output' }.sort { it[0].replicate }.collect { it[1] }
            tuple([sample: base], inputs, outputs)
        }
        .filter { smeta, ins, outs -> ins && outs }
        .set { ch_fitness_bundled }

    // 2. Default fitness estimation: one task per sample merges the counts, derives the experimental
    //    design, picks the synonymous wildtype proxy, estimates fitness and draws the QC plots and
    //    heatmaps. The run-level inputs are value channels, so every sample gets them.
    FITNESS_ESTIMATION(
        ch_fitness_bundled,                       // tuple(val(sample), path(inputs), path(outputs))
        ch_samplesheet_csv,                       // path(samplesheet)
        ch_fasta.map { _meta, fasta -> fasta },   // path(fasta)
        val_reading_frame,                        // val(pos_range)
        ch_aa_seq.map { _meta, aa -> aa },        // path(aa_seq.txt)
        params.min_counts,
        params.output_pseudocount
    )
    ch_versions = ch_versions.mix(FITNESS_ESTIMATION.out.versions)

    // 3. The tables DiMSum and mutscan take, aligned per sample: the merged counts with their own
    //    synonymous wildtype proxy and experimental design (all three are outputs of the same task)
    def ch_tables = FITNESS_ESTIMATION.out.merged_counts
        .join(FITNESS_ESTIMATION.out.synonymous_wt)
        .join(FITNESS_ESTIMATION.out.experimental_design)   // (sample, counts, wt, design)

    def ch_run_counts_d = ch_tables.map { smp, counts, _wt, _design -> tuple(smp, counts) }
    def ch_run_wt_d     = ch_tables.map { _smp, _counts, wt, _design -> wt }
    def ch_run_exp_d    = ch_tables.map { _smp, _counts, _wt, design -> design }

    // DiMSum and mutscan are optional, so their artefact channels start empty and are only filled
    // when the tool actually ran - the summary report then simply omits the section.
    def ch_dimsum_dir    = Channel.empty()
    def ch_mutscan_plots = Channel.empty()

    if (val_run_dimsum) {
        RUN_DIMSUM(
            ch_run_counts_d,
            ch_run_wt_d,
            ch_run_exp_d,
            params.min_counts,
            params.output_pseudocount,
            file("${projectDir}/assets/dimsum/report.Rmd", checkIfExists: true),
            file("${projectDir}/assets/dimsum/Dumpling.png", checkIfExists: true)
        )
        ch_versions = ch_versions.mix(RUN_DIMSUM.out.versions)
        ch_dimsum_dir = RUN_DIMSUM.out.results_dir
    }

    // 4. Run Mutscan
    if (val_run_mutscan) {
        RUN_MUTSCAN(
            ch_run_counts_d,
            ch_run_wt_d,
            ch_run_exp_d,
            params.min_counts
        )
        ch_versions = ch_versions.mix(RUN_MUTSCAN.out.versions)
        ch_mutscan_plots = RUN_MUTSCAN.out.qc_plots
    }

    emit:
    fitness_estimation = FITNESS_ESTIMATION.out.fitness_estimation
    // artefacts the all-in-one report embeds
    fitness_plots      = FITNESS_ESTIMATION.out.counts_corr_pdf
                            .mix(FITNESS_ESTIMATION.out.fitness_corr_pdf)
                            .mix(FITNESS_ESTIMATION.out.fitness_heatmap)
    mutscan_plots      = ch_mutscan_plots
    dimsum_dir         = ch_dimsum_dir
    versions           = ch_versions
}
