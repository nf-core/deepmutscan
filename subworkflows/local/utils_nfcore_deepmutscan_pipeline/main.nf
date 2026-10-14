//
// Subworkflow with functionality specific to the nf-core/deepmutscan pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList        } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/deepmutscan ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/deepmutscan/blob/master/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,                                 // input_workflow
        validate_params,                          // validate_params
        "${projectDir}/nextflow_schema.json",     // parameters_schema
        params.help,                              // help
        params.help_full,                         // help_full
        params.show_hidden,                       // show_hidden
        before_text,                              // before_text
        after_text,                               // after_text
        command,                                  // command
        true                                      // cli_typecast: cast CLI strings for validation
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //
    // samplesheetToList validates against assets/schema_input.json and hands back proper
    // file() objects (so relative paths in the samplesheet work), the meta mapping, and
    // duplicate detection via the schema's uniqueEntries.
    //
    def samplesheet_rows = samplesheetToList(params.input, "${projectDir}/assets/schema_input.json")

    validateInputSamplesheet(samplesheet_rows)

    def ch_samplesheet = Channel
        .fromList(samplesheet_rows)
        .map { meta, file1, file2 ->
            // Every library is paired-end (the schema requires file2). The _pe suffix is kept
            // because it is part of every published output name.
            def row_meta = meta + [
                id        : "${meta.sample}_${meta.type}_${meta.replicate}_pe".toString(),
                replicate : meta.replicate as int,
            ]
            tuple(row_meta, [file1, file2])
        }

    // Emit the samplesheet channel and an empty version channel for use in the workflow
    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Params given on the command line arrive as Strings, so `--fitness false` yields the String
// "false", which is truthy in Groovy - `if (params.fitness)` would switch fitness ON. Values set
// in a config or a -params-file keep their real type. These two helpers accept either shape.
// nf-schema's cast_cli_params only casts for validation; it does not change params itself.
//
def asBool(value) {
    if (value == null) {
        return false
    }
    if (value instanceof Boolean) {
        return value
    }
    return !(value.toString().trim().toLowerCase() in ['false', '0', 'no', 'null', ''])
}

def asInt(value) {
    if (value == null) {
        return null
    }
    if (value instanceof Number) {
        return value as int
    }
    return value.toString().trim() as int
}

//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    validateReadingFrame()
    validateCustomCodonLibrary()
}

//
// --reading_frame must describe whole codons and lie inside the reference. A partial last codon
// otherwise surfaces as a KeyError inside a counter worker process, which is hard to read.
//
def validateReadingFrame() {
    def matcher = (params.reading_frame as String) =~ /^(\d+)-(\d+)$/
    if (!matcher.matches()) {
        error("--reading_frame must be given as '<start>-<end>', 1-based and inclusive (e.g. '352-1383'); got '${params.reading_frame}'.")
    }
    def start = matcher[0][1] as int
    def end   = matcher[0][2] as int

    if (start < 1) {
        error("--reading_frame start must be >= 1 (1-based, inclusive); got '${params.reading_frame}'.")
    }
    if (end < start) {
        error("--reading_frame end must not be smaller than its start; got '${params.reading_frame}'.")
    }
    if ((end - start + 1) % 3 != 0) {
        error(
            "--reading_frame '${params.reading_frame}' spans ${end - start + 1} nt, which is not a whole " +
            "number of codons. Give a range whose length is divisible by 3."
        )
    }

    // Length check needs the sequence; plain FASTA only, so a gzipped reference simply skips it.
    if (params.fasta && !(params.fasta as String).endsWith('.gz')) {
        def ref = file(params.fasta)
        if (ref.exists()) {
            def n_seqs = 0
            def seq_len = 0
            ref.eachLine { line ->
                if (line.startsWith('>')) {
                    n_seqs = n_seqs + 1
                }
                else {
                    seq_len += line.trim().length()
                }
            }
            if (n_seqs > 1) {
                error("--fasta must contain exactly one sequence (found ${n_seqs}); deepmutscan counts variants against a single ORF.")
            }
            if (seq_len > 0 && end > seq_len) {
                error("--reading_frame '${params.reading_frame}' ends past the end of --fasta (${seq_len} nt).")
            }
        }
    }
}

//
// A custom library is the only way to describe a synthesised (non-degenerate) oligo library,
// so asking for one without supplying it is always a mistake.
//
def validateCustomCodonLibrary() {
    if (params.mutagenesis_type == 'custom' && !params.custom_codon_library) {
        error("--mutagenesis_type custom requires --custom_codon_library <file>.")
    }
    if (params.custom_codon_library && params.mutagenesis_type != 'custom') {
        log.warn("--custom_codon_library is ignored unless --mutagenesis_type custom is set.")
    }
}

//
// Validate the parsed samplesheet as a whole
//
def validateInputSamplesheet(rows) {
    // One biological sample per run. The DiMSum experimental design, the fitness output names
    // and the summary report are all run-level, so a second sample silently overwrites the
    // first rather than producing a second set of results.
    def samples = rows.collect { row -> row[0].sample }.unique()
    if (samples.size() > 1) {
        error(
            "nf-core/deepmutscan ${workflow.manifest.version} supports one biological sample per run " +
            "(found: ${samples.join(', ')}).\n" +
            "  Please split the samplesheet and run the pipeline once per sample."
        )
    }

    if (params.error_correction == 'wildtype') {
        def wt_rows = rows.findAll { row -> row[0].type == 'wildtype' }
        if (wt_rows.isEmpty()) {
            error(
                "[--error_correction wildtype] No 'wildtype' row found in the samplesheet.\n" +
                "  Add a row with type 'wildtype' (deep sequencing of the unmutated template),\n" +
                "  or use the internal correction instead (--error_correction false_doubles),\n" +
                "  or disable error correction with --error_correction none."
            )
        }
        if (wt_rows.size() > 1) {
            error(
                "[--error_correction wildtype] Found ${wt_rows.size()} 'wildtype' rows, but the correction " +
                "pairs each library with exactly one wildtype profile.\n" +
                "  Please supply a single wildtype library per run."
            )
        }
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // Only the tools this run actually used, so the methods text matches the command line.
    def citations = [
        "Tools used in the workflow included:",
        "FastQC (Andrews 2010),",
        "BWA-MEM (Li 2013),",
        "SAMtools (Danecek et al. 2021),",
        "VSEARCH (Rognes et al. 2016),",
        "pysam (Heger et al.),",
        "Polars (Vink et al.),",
        "Biopython (Cock et al. 2009),",
        "Biostrings (Pages et al.),",
        asBool(params.dimsum) ? "DiMSum (Faure et al. 2020)," : "",
        asBool(params.mutscan) ? "mutscan (Soneson et al. 2023)," : "",
        asBool(params.mutscan) ? "edgeR (Robinson et al. 2010)," : "",
        asBool(params.mutscan) ? "limma (Ritchie et al. 2015)," : "",
        params.pdb ? "3Dmol.js (Rego and Koes 2015)," : "",
        "MultiQC (Ewels et al. 2016)",
        "."
    ]

    return citations.findAll { entry -> entry != "" }.join(' ').trim()
}

def toolBibliographyText() {
    def references = [
        "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/.</li>",
        "<li>Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997. doi: 10.48550/arXiv.1303.3997</li>",
        "<li>Danecek P, Bonfield JK, Liddle J, et al. (2021) Twelve years of SAMtools and BCFtools. Gigascience, 10(2), giab008. doi: 10.1093/gigascience/giab008</li>",
        "<li>Rognes T, Flouri T, Nichols B, Quince C, Mahe F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ, 4, e2584. doi: 10.7717/peerj.2584</li>",
        "<li>Heger A, Jacobs K, et al. pysam: a Python module for reading, manipulating and writing genomic data sets. URL: https://github.com/pysam-developers/pysam</li>",
        "<li>Vink R, et al. Polars: a lightning-fast DataFrame library for Rust and Python. URL: https://github.com/pola-rs/polars</li>",
        "<li>Cock PJ, Antao T, Chang JT, et al. (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422-1423. doi: 10.1093/bioinformatics/btp163</li>",
        "<li>Pages H, Aboyoun P, Gentleman R, DebRoy S. Biostrings: efficient manipulation of biological strings. R package. doi: 10.18129/B9.bioc.Biostrings</li>",
        asBool(params.dimsum) ? "<li>Faure AJ, Schmiedel JM, Baeza-Centurion P, Lehner B. (2020) DiMSum: an error model and pipeline for analyzing deep mutational scanning data and diagnosing common experimental pathologies. Genome Biology, 21(1), 207. doi: 10.1186/s13059-020-02091-3</li>" : "",
        asBool(params.mutscan) ? "<li>Soneson C, Bendel AM, Diss G, Stadler MB. (2023) mutscan - a flexible R package for efficient end-to-end analysis of multiplexed assays of variant effect data. Genome Biology, 24(1), 132. doi: 10.1186/s13059-023-02967-0</li>" : "",
        asBool(params.mutscan) ? "<li>Robinson MD, McCarthy DJ, Smyth GK. (2010) edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26(1), 139-140. doi: 10.1093/bioinformatics/btp616</li>" : "",
        asBool(params.mutscan) ? "<li>Ritchie ME, Phipson B, Wu D, et al. (2015) limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research, 43(7), e47. doi: 10.1093/nar/gkv007</li>" : "",
        params.pdb ? "<li>Rego N, Koes D. (2015) 3Dmol.js: molecular visualization with WebGL. Bioinformatics, 31(8), 1322-1324. doi: 10.1093/bioinformatics/btu829</li>" : "",
        "<li>Ewels, P., Magnusson, M., Lundin, S., & Kaller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047-3048. doi: 10.1093/bioinformatics/btw354</li>"
    ]

    return references.findAll { entry -> entry != "" }.join(' ').trim()
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
