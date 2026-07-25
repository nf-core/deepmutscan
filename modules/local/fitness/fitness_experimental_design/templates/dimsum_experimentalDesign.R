#!/usr/bin/env Rscript

# Make a DiMSum experimental design from a deepmutscan samplesheet.
# - samplesheet_csv: path to CSV with columns sample,type,replicate,file1,file2
# - out_path: where to write the TSV (default "experimentalDesign.tsv")
# Returns: the experimental design as a data.frame
make_dimsum_experimental_design <- function(samplesheet_csv, out_path = "experimentalDesign.tsv") {
  # ---- read & normalize ----
  ss <- read.csv(samplesheet_csv, stringsAsFactors = FALSE, check.names = FALSE)
  names(ss) <- tolower(names(ss))

  # tolerate missing file2 column (single-end)
  if (!"file2" %in% names(ss)) ss\$file2 <- ""

  required <- c("sample", "type", "replicate", "file1", "file2")
  missing  <- setdiff(required, names(ss))
  if (length(missing) > 0) stop("Samplesheet missing columns: ", paste(missing, collapse = ", "))

  # coerce types
  ss\$replicate <- as.integer(ss\$replicate)

  # ---- keep only the samples DiMSum scores ----
  # A samplesheet legitimately carries other row types: `wildtype` (required by
  # --error_correction wildtype) and `quality`. They are real samples, but they are not part of the
  # selection design - they feed error correction and QC - and MERGE_COUNTS drops them too. Leaving
  # them in gave them an empty `selection_id`, which DiMSum reads back as NA and then dies in
  # dimsum__check_experiment_design on `sum(exp_design[,"selection_id"] < 0)` with the opaque
  # "missing value where TRUE/FALSE needed".
  ss <- ss[ss\$type %in% c("input", "output"), , drop = FALSE]
  if (nrow(ss) == 0) {
    stop("Samplesheet has no 'input'/'output' rows - DiMSum needs both to compute fitness.")
  }

  # ---- order rows: by sample (if multiple), input before output, then replicate ----
  # Sorting happens before naming, because the names are derived from the row order below.
  multi_base <- length(unique(ss\$sample)) > 1
  type_rank  <- match(ss\$type, c("input", "output"))
  ord <- if (multi_base) {
    order(ss\$sample, type_rank, ss\$replicate)
  } else {
    order(type_rank, ss\$replicate)
  }
  ss <- ss[ord, , drop = FALSE]

  # ---- sample_name: must equal the count column MERGE_COUNTS writes for this row ----
  # merge_counts.R names its columns positionally - paste0("input", seq_len(n_in)) over the files
  # sorted by replicate - so the name has to come from each row's POSITION within its type, not from
  # the samplesheet's `replicate` value. The two only coincide when the replicates happen to be a
  # contiguous 1..N; replicates 1,2,4 would otherwise produce "input4" for a column named "input3"
  # and DiMSum could not match them. DiMSum also requires sample_name to be alphanumeric and to
  # start with a letter, so the parts are joined without a separator.
  idx <- stats::ave(seq_len(nrow(ss)), ss\$sample, ss\$type, FUN = seq_along)
  # NOTE: for multiple biological samples these names still will not match the counts file, because
  # one design is built from the whole samplesheet and broadcast to every per-sample DiMSum run
  # (see subworkflows/local/calculate_fitness/main.nf). Multi-sample + DiMSum is a separate,
  # known-broken case; single-sample runs are correct.
  sample_name <- if (multi_base) paste0(ss\$sample, ss\$type, idx) else paste0(ss\$type, idx)

  # ---- build DiMSum columns ----
  experiment_replicate <- ss\$replicate           # pairs an input with an output; need not be 1..N
  selection_id <- ifelse(ss\$type == "input", 0L, 1L)
  # assume one selection batch; DiMSum requires this to be blank for input samples
  selection_replicate <- ifelse(ss\$type == "output", 1L, NA_integer_)
  # assume one technical batch
  technical_replicate <- rep(1L, nrow(ss))

  pair1 <- basename(ss\$file1)
  # keep empty string for single-end / missing file2
  pair2 <- ifelse(is.na(ss\$file2) | ss\$file2 == "", "", basename(ss\$file2))

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
  rownames(ed) <- NULL

  # ---- write & return ----
  write.table(ed, file = out_path, sep = "\\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "")
  return(ed)
}

#####
# run function
#####
make_dimsum_experimental_design(
  samplesheet_csv = "$samplesheet_csv",
  out_path = "experimentalDesign.tsv"
)

####
# create versions.yml
####
r_version <- strsplit(version[['version.string']], ' ')[[1]][3]

if (is.null(r_version)) r_version <- "unknown"

f <- file("versions.yml", "w")
writeLines(
  c(
    '"${task.process}":',
    paste('    r-base:', r_version)
  ),
  f
)
close(f)
