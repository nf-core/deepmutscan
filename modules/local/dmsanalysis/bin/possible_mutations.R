# input: wildtype-seq (string or fasta file), start&stopp pos., output_folder_path, mutagenesis_type (choose from nnk, nns, max_diff_to_wt, custom), if you choose custom, add: custom_codon_library -> comma-separated .txt (-> "AAA, AAC, AAG, AAT, ...")
# output: .csv with all possible programmed codons for each position

suppressMessages(library(Biostrings))

# Define the function
generate_possible_variants <- function(wt_seq_input, start_stop_pos, mutagenesis_type = "nnk",
                                       custom_codon_library = NULL, output_file) {
  # Parse the start and stop positions from the input format "start-stop"
  positions <- unlist(strsplit(start_stop_pos, "-"))
  start_pos <- as.numeric(positions[1])
  stop_pos <- as.numeric(positions[2])

  # Check if the input is a file or a string
  if (file.exists(wt_seq_input)) {
    # If it's a file, read the sequence from the fasta file
    seq_data <- readDNAStringSet(filepath = wt_seq_input)
    wt_seq <- seq_data[[1]]  # Extract the sequence
  } else {
    # Otherwise, treat the input as a sequence string
    wt_seq <- DNAString(wt_seq_input)
  }

  # Extract the sequence between start and stop codons
  coding_seq <- subseq(wt_seq, start = start_pos, end = stop_pos)
  coding_seq <- as.character(coding_seq)

  # List of predefined NNK & NNS codons
  nnk_codons <- c('AAG', 'AAT', 'ATG', 'ATT', 'AGG', 'AGT', 'ACG', 'ACT',
                  'TAG', 'TAT', 'TTG', 'TTT', 'TGG', 'TGT', 'TCG', 'TCT',
                  'GAG', 'GAT', 'GTG', 'GTT', 'GGG', 'GGT', 'GCG', 'GCT',
                  'CAG', 'CAT', 'CTG', 'CTT', 'CGG', 'CGT', 'CCG', 'CCT')

  nns_codons <- c('AAG', 'AAC', 'ATG', 'ATC', 'AGG', 'AGC', 'ACG', 'ACC',
                  'TAG', 'TAC', 'TTG', 'TTC', 'TGG', 'TGC', 'TCG', 'TCC',
                  'GAG', 'GAC', 'GTG', 'GTC', 'GGG', 'GGC', 'GCG', 'GCC',
                  'CAG', 'CAC', 'CTG', 'CTC', 'CGG', 'CGC', 'CCG', 'CCC')

  # Function to split a DNA sequence into codons (triplets)
  split_into_codons <- function(seq) {
    return(strsplit(seq, "(?<=.{3})", perl = TRUE)[[1]])
  }

  # Read custom codons if mode is 'custom'
  if (mutagenesis_type == "custom") {
    if (is.null(custom_codon_library) || !file.exists(custom_codon_library)) {
      stop("Custom codons file must be provided and valid when using 'custom' mutagenesis_type")
    }
    # Read and parse the custom codons from the file
    custom_codons <- unlist(strsplit(readLines(custom_codon_library), ","))
    custom_codons <- trimws(custom_codons)  # Remove any whitespace
  }

  # Split wild-type sequence into codons
  wt_codons <- split_into_codons(coding_seq)

  # Initialize DataFrame to store mutated variants
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
      result <- rbind(result, data.frame(Codon_Number = i, wt_codon = wt_codon, Variant = variant, stringsAsFactors = FALSE))
    }
  }

  # Save the variants into a CSV file
  write.csv(result, output_file, row.names = FALSE)
}

# Example usage
# Possibly generate a custom codons file: "AAA, AAC, AAG, AAT, ..."
# generate_possible_variants("/Users/benjaminwehnert/CRG/DMS_QC/testing_data/MORtn5_reference.fa", "23-1225", "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/possible_NNK_mutations.csv", mutagenesis_type = "nnk")    ### this one's correct for the dataset
# generate_possible_variants("/Users/benjaminwehnert/CRG/DMS_QC/testing_data/MORtn5_reference.fa", "23-1225", "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/testing_outputs/possible_NNK_mutations_taylors_nnk_and_nns.csv", mutagenesis_type = "max_diff_to_wt")





