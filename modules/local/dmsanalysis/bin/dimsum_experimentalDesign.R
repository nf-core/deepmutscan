# Make a DiMSum experimental design from a deepmutscan samplesheet.
# - samplesheet_csv: path to CSV with columns sample,type,replicate,file1,file2
# - out_path: where to write the TSV (default "experimentalDesign.tsv")
# Returns: the experimental design as a data.frame
make_dimsum_experimental_design <- function(samplesheet_csv, out_path = "experimentalDesign.tsv") {
  # ---- read & normalize ----
  ss <- read.csv(samplesheet_csv, stringsAsFactors = FALSE, check.names = FALSE)
  names(ss) <- tolower(names(ss))
  
  # tolerate missing file2 column (single-end)
  if (!"file2" %in% names(ss)) ss$file2 <- ""
  
  required <- c("sample", "type", "replicate", "file1", "file2")
  missing  <- setdiff(required, names(ss))
  if (length(missing) > 0) stop("Samplesheet missing columns: ", paste(missing, collapse = ", "))
  
  # coerce types
  ss$replicate <- as.integer(ss$replicate)
  
  # ---- derive sample_name strategy ----
  # If only one biological sample present (e.g. one protein), use "input1", "output2", ...
  # If multiple biological samples present, prefix with 'sample' to avoid collisions:
  # "GID1A_input1", "GID1B_output2", ...
  multi_base <- length(unique(ss$sample)) > 1
  if (multi_base) {
    sample_name <- paste(ss$sample, ss$type, ss$replicate, sep = "_")
  } else {
    sample_name <- paste0(ss$type, ss$replicate)
  }
  
  # ---- build DiMSum columns ----
  experiment_replicate <- ss$replicate
  selection_id <- ifelse(ss$type == "input", 0L,
                         ifelse(ss$type == "output", 1L, NA_integer_))
  # assume one selection batch
  selection_replicate <- ifelse(ss$type == "output", 1L, NA_integer_)
  # assume one technical batch
  technical_replicate <- rep(1L, nrow(ss))
  
  pair1 <- basename(ss$file1)
  # keep empty string for single-end / missing file2
  pair2 <- ifelse(is.na(ss$file2) | ss$file2 == "", "", basename(ss$file2))
  
  ed <- data.frame(
    sample_name          = sample_name,
    experiment_replicate = experiment_replicate,
    selection_id         = selection_id,
    selection_replicate  = selection_replicate,
    technical_replicate  = technical_replicate,
    pair1                = pair1,
    pair2                = pair2,
    stringsAsFactors     = FALSE
  )
  
  # ---- order rows: by sample (if multiple), type (input, output, quality), then replicate ----
  type_rank <- match(ss$type, c("input", "output", "quality"))
  ord <- if (multi_base) {
    order(ss$sample, type_rank, ss$replicate, na.last = TRUE)
  } else {
    order(type_rank, ss$replicate, na.last = TRUE)
  }
  ed <- ed[ord, , drop = FALSE]
  rownames(ed) <- NULL
  
  # ---- write & return ----
  write.table(ed, file = out_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "")
  return(ed)
}

