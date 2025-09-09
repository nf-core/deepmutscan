# merge_counts.R
# Merges DiMSum-ready count tables into a single matrix with columns:
# nt_seq, input1..inputN, output1..outputM

suppressMessages({
  # keep base R; if you prefer data.table, we can swap in later
})

# ---- Core function ----
merge_dimsum_counts <- function(input_paths, output_paths, out_path = "counts.tsv") {
  # input_paths, output_paths: character vectors of file paths
  # out_path: output TSV path
  
  # Helper to read a 2-col TSV without header -> data.frame(nt_seq, count)
  read_counts <- function(fp) {
    df <- utils::read.table(
      fp, header = FALSE, sep = "\t", quote = "",
      col.names = c("nt_seq", "count"),
      colClasses = c("character", "numeric"),
      comment.char = "", check.names = FALSE
    )
    df
  }
  
  # Read all inputs / outputs
  input_list  <- lapply(input_paths,  read_counts)
  output_list <- lapply(output_paths, read_counts)
  
  # Collect universe of sequences
  all_seqs <- unique(c(
    unlist(lapply(input_list,  function(x) x$nt_seq)),
    unlist(lapply(output_list, function(x) x$nt_seq))
  ))
  
  # Pre-allocate output frame
  n_in  <- length(input_list)
  n_out <- length(output_list)
  col_names <- c(
    "nt_seq",
    if (n_in  > 0) paste0("input",  seq_len(n_in))  else character(0),
    if (n_out > 0) paste0("output", seq_len(n_out)) else character(0)
  )
  # initialize numeric with 0 counts
  out <- data.frame(
    nt_seq = all_seqs,
    matrix(0, nrow = length(all_seqs), ncol = n_in + n_out),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  names(out) <- col_names
  
  # Fill inputs
  if (n_in > 0) {
    for (i in seq_len(n_in)) {
      df <- input_list[[i]]
      idx <- match(df$nt_seq, out$nt_seq)
      out[idx, paste0("input", i)] <- df$count
    }
  }
  
  # Fill outputs
  if (n_out > 0) {
    for (j in seq_len(n_out)) {
      df <- output_list[[j]]
      idx <- match(df$nt_seq, out$nt_seq)
      out[idx, paste0("output", j)] <- df$count
    }
  }
  
  # Write TSV with header, no row names, no quotes
  utils::write.table(out, file = out_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  invisible(out)
}

# ---- Lightweight CLI wrapper ----
# Usage:
# Rscript merge_counts.R --inputs <i1> <i2> ... --outputs <o1> <o2> ... --out counts.tsv
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  
  # simple flag parser that supports space-separated lists
  get_vals <- function(flag) {
    if (!(flag %in% args)) return(character(0))
    start <- which(args == flag)
    stop  <- which(args %in% c("--inputs","--outputs","--out"))
    stop  <- min(c(stop[stop > start], length(args) + 1)) - 1
    if (stop <= start) return(character(0))
    args[(start + 1):stop]
  }
  
  input_paths  <- get_vals("--inputs")
  output_paths <- get_vals("--outputs")
  out_path     <- get_vals("--out")
  out_path     <- if (length(out_path)) out_path[1] else "counts.tsv"
  
  if (!length(input_paths) && !length(output_paths)) {
    stop("No inputs/outputs provided. Use --inputs <files> and/or --outputs <files>.")
  }
  
  merge_dimsum_counts(input_paths, output_paths, out_path)
}
