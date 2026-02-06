#!/usr/bin/env Rscript

# input: wildtype-seq (string or fasta file), start&stopp pos., output_folder_path, mutagenesis_type (choose from nnk, nns, max_diff_to_wt, custom), if you choose custom, add: custom_codon_library -> comma-separated .txt (-> "AAA, AAC, AAG, AAT, ...")
# output: .csv with all possible programmed codons for each position

# ----------------------------------------------------------------------
# 1. SETUP AND NEXTFLOW INPUT INJECTION
# ----------------------------------------------------------------------

suppressMessages(library(Biostrings))
suppressMessages(library(methods))

generate_possible_variants <- function(wt_seq_input, start_stop_pos, mutagenesis_type,
                                       custom_codon_library_path, output_file) {

    # Parse the start and stop positions from the input format "start-stop"
    positions <- unlist(strsplit(start_stop_pos, "-"))
    start_pos <- as.numeric(positions[1])
    stop_pos <- as.numeric(positions[2])

    # Check if the input is a file or a string
    if (file.exists(wt_seq_input)) {
        seq_data <- Biostrings::readDNAStringSet(filepath = wt_seq_input)
        wt_seq <- seq_data[[1]]  # Extract the sequence
    } else {
        wt_seq <- Biostrings::DNAString(wt_seq_input)
    }

    # Extract the sequence between start and stop codons
    coding_seq <- Biostrings::subseq(wt_seq, start = start_pos, end = stop_pos)
    coding_seq <- as.character(coding_seq)

    # List of predefined NNK & NNS codons (logic kept the same)
    nnk_codons <- c('AAG', 'AAT', 'ATG', 'ATT', 'AGG', 'AGT', 'ACG', 'ACT', 'TAG', 'TAT', 'TTG', 'TTT', 'TGG', 'TGT', 'TCG', 'TCT', 'GAG', 'GAT', 'GTG', 'GTT', 'GGG', 'GGT', 'GCG', 'GCT', 'CAG', 'CAT', 'CTG', 'CTT', 'CGG', 'CGT', 'CCG', 'CCT')
    nns_codons <- c('AAG', 'AAC', 'ATG', 'ATC', 'AGG', 'AGC', 'ACG', 'ACC', 'TAG', 'TAC', 'TTG', 'TTC', 'TGG', 'TGC', 'TCG', 'TCC', 'GAG', 'GAC', 'GTG', 'GTC', 'GGG', 'GGC', 'GCG', 'GCC', 'CAG', 'CAC', 'CTG', 'CTC', 'CGG', 'CGC', 'CCG', 'CCC')

    # Read custom codons if mode is 'custom' (LOGIC ADAPTED)
    if (mutagenesis_type == "custom") {
        if (!file.exists(custom_codon_library_path) || is.null(custom_codon_library_path)) {
            # This should not happen if the Nextflow filter works, but R needs a check.
            stop("Custom codons file must be provided and valid when using 'custom' mutagenesis_type.")
        }
        # Read and parse the custom codons from the file
        custom_codons <- unlist(strsplit(readLines(custom_codon_library_path), ","))
        custom_codons <- trimws(custom_codons)  # Remove any whitespace
    } else {
        custom_codons <- NULL # Ensure it's defined
    }

    # Function to split a DNA sequence into codons (triplets)
    split_into_codons <- function(seq) {
        # Important: double escaping for perl regex
        return(strsplit(seq, "(?<=.{3})", perl = TRUE)[[1]])
    }

    # Split wild-type sequence into codons
    wt_codons <- split_into_codons(coding_seq)

    # Initialize DataFrame to store mutated variants
    # WICHTIG: \$ escaping
    result <- data.frame(Codon_Number = integer(), wt_codon = character(), Variant = character(), stringsAsFactors = FALSE)

    # Determine the codon list based on the mutagenesis_type
    get_codon_list <- function(wt_codon) {
        if (mutagenesis_type == "nnk") {
            return(nnk_codons)
        } else if (mutagenesis_type == "nns") {
            return(nns_codons)
        } else if (mutagenesis_type == "max_diff_to_wt") {
            if (substr(wt_codon, 3, 3) == "T") {
                return(nns_codons)
            } else {
                return(nnk_codons)
            }
        } else if (mutagenesis_type == "custom") {
            if (is.null(custom_codons)) stop("Custom codons list not defined.")
            return(custom_codons)
        } else {
            stop("Invalid mutagenesis_type Choose from 'nnk', 'nns', 'max_diff_to_wt', or 'custom'.")
        }
    }

    # Iterate over each codon in the wild-type sequence
    for (i in seq_along(wt_codons)) {
        wt_codon <- wt_codons[i]
        codon_list <- get_codon_list(wt_codon)

        # Filter codons that are different from the wild-type codon
        possible_variants <- codon_list[codon_list != wt_codon]

        # Add all variants to the result list, including the wild-type codon
        for (variant in possible_variants) {
            # WICHTIG: \$ escaping
            result <- rbind(result, data.frame(Codon_Number = i, wt_codon = wt_codon, Variant = variant, stringsAsFactors = FALSE))
        }
    }

    # Save the variants into a CSV file
    write.csv(result, output_file, row.names = FALSE)
}

# ----------------------------------------------------------------------
# 2. MAIN EXECUTION BLOCK (Nextflow substitution occurs here)
# ----------------------------------------------------------------------

# The Bash 'if/else' is replaced by R's conditional logic on the string "/NULL"
# If the file input was optional and missing, Nextflow passes the string "/NULL".
# If it's "/NULL", we pass the R NULL object to the function.

custom_lib_arg <- if ("$custom_codon_library" == "/NULL") {
    NULL
} else {
    "$custom_codon_library"
}

#####
# run function
#####
generate_possible_variants(
    wt_seq_input          = "$wt_seq",
    start_stop_pos        = "$pos_range",
    mutagenesis_type      = "$mutagenesis_type",
    custom_codon_library_path = "$custom_codon_library",	# Will be '/NULL' or a file path
    output_file           = "possible_mutations.csv"
)

# ----------------------------------------------------------------------
# 3. VERSIONING
# ----------------------------------------------------------------------

# Extract R base and Biostrings versions
r_version <- strsplit(version[['version.string']], ' ')[[1]][3]
biostrings_version <- as.character(packageVersion("Biostrings"))

f <- file("versions.yml", "w")
writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r_version),
        paste('    biostrings:', biostrings_version)
    ),
    f
)
close(f)

