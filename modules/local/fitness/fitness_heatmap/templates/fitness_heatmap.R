#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(methods)
  library(dplyr)
  library(ggplot2)
  library(grid)        # for unit()
})

# ---------- helper functions ----------
find_col <- function(df, candidates) {
  norm <- function(x) gsub("[^a-z0-9]+", "_", tolower(x))
  nms <- colnames(df); nn <- norm(nms)
  for (cand in candidates) {
    hit <- which(nn == norm(cand))
    if (length(hit) == 1) return(nms[hit])
  }
  stop(sprintf("Could not find any of columns: %s", paste(candidates, collapse = ", ")))
}

get_rescaled_cols <- function(df) {
  nms  <- colnames(df)
  hits <- grep("^rescaled[_ ]?fitness", nms, ignore.case = TRUE, value = TRUE)
  if (!length(hits)) stop("No 'rescaled_fitness' columns found.")
  idx <- suppressWarnings(as.integer(gsub(".*?([0-9]+)\$", "\\\\1", hits))) 	# additional backslashs to make it groovy readable
  hits[order(is.na(idx), idx, hits)]
}

# Find "mean fitness" column (if present)
find_mean_col <- function(df) {
  nms <- colnames(df)
  key <- tolower(gsub("[^a-z0-9]+", "_", nms))
  hit <- which(key == "mean_fitness")
  if (length(hit) == 1) nms[hit] else NULL
}

# NEW: read WT amino acid sequence from .txt file (single line)
read_wt_seq_aa_txt <- function(path) {
  if (is.null(path)) stop("wt_seq_aa_txt_path must be provided.")
  x <- readLines(path, warn = FALSE)
  x <- x[nzchar(x)]
  if (!length(x)) stop("WT AA TXT is empty.")
  aa <- toupper(gsub("\\\\s+", "", x[which.max(nchar(x))]))
  # Keep only valid amino acid letters (including stop '*')
  aa <- gsub("[^ACDEFGHIKLMNPQRSTVWY*]", "", aa)
  if (!nchar(aa)) stop("WT AA TXT contains no valid AA letters.")
  aa
}

# Build full AA×position grid for positions 1..wt_len (from WT sequence),
# join fitness values; any missing combos stay NA (-> grey).
# Positions > wt_len are considered padded and will be white.
build_heatmap_long <- function(df,
                               wt_aa_col,
                               pos_col,
                               mut_aa_col,
                               fitness_col,
                               positions_per_row = 75,
                               wt_seq_aa,
                               fill_missing_as_zero = FALSE) {

  # authoritative WT length from provided sequence
  letters <- strsplit(wt_seq_aa, "", fixed = TRUE)[[1]]
  wt_len  <- length(letters)

  # normalize data
  df0 <- df %>%
    transmute(
      position = suppressWarnings(as.numeric(.data[[pos_col]])),
      wt_aa_in = .data[[wt_aa_col]],
      mut_aa   = .data[[mut_aa_col]],
      fitness  = suppressWarnings(as.numeric(.data[[fitness_col]]))
    ) %>%
    filter(is.finite(position))

  # drop any rows that claim positions beyond WT length
  if (nrow(df0) && any(df0\$position > wt_len, na.rm = TRUE)) {
    dropped <- sum(df0\$position > wt_len, na.rm = TRUE)
    warning(sprintf("Dropping %d row(s) with position > WT length (%d).", dropped, wt_len))
    df0 <- df0 %>% filter(position <= wt_len)
  }

  # pad to next multiple of 75 (by rows)
  rem       <- wt_len %% positions_per_row
  pad_need  <- if (rem == 0) 0 else positions_per_row - rem
  max_paded <- wt_len + pad_need

  # full grid: positions 1..max_paded (so the tail exists), AA set of 21
  all_positions   <- seq_len(max_paded)
  all_amino_acids <- c("A","C","D","E","F","G","H","I","K","L",
                       "M","N","P","Q","R","S","T","V","W","Y","*")

  grid_df <- expand.grid(position = all_positions,
                         mut_aa   = all_amino_acids,
                         KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) %>%
    mutate(is_padded = position > wt_len)

  # join fitness only for real positions (<= wt_len)
  fit_df <- df0 %>% select(position, mut_aa, fitness)
  d <- grid_df %>%
    left_join(fit_df, by = c("position","mut_aa"))

  if (fill_missing_as_zero) {
    d\$fitness[is.na(d\$fitness) & d\$position <= wt_len] <- 0
  }

  # authoritative WT AA per real position; tail gets placeholder 'Y'
  wt_map <- tibble(position = seq_len(wt_len), wt_aa = letters)
  d <- d %>%
    left_join(wt_map, by = "position") %>%
    mutate(wt_aa = ifelse(is.na(wt_aa) & position > wt_len, "Y", wt_aa))

  # layout fields
  d <- d %>%
    mutate(
      row_group = ((position - 1) %/% positions_per_row) + 1,
      wt_aa_pos = paste0(wt_aa, position),
      wt_aa_pos = factor(wt_aa_pos, levels = unique(wt_aa_pos)),
      synonymous = mut_aa == wt_aa
    )

  # IMPORTANT: use WT length as the true end of the protein
  d\$max_pos <- wt_len
  d
}

syn_segments <- function(d, positions_per_row = 75) {
  amino_order <- rev(c("G", "A", "V", "L", "M", "I", "F",
                         "Y", "W", "K", "R", "H", "D", "E",
                         "S", "T", "C", "N", "Q", "P", "*"))
  d %>%
    mutate(
      mut_aa = factor(mut_aa, levels = amino_order),
      x = as.numeric(factor(wt_aa_pos, levels = levels(wt_aa_pos))) -
          ((row_group - 1) * positions_per_row),
      y = as.numeric(factor(mut_aa, levels = amino_order))
    ) %>%
    filter(synonymous, position <= max_pos)
}

# Draw one solid white rectangle per row group covering the padded tail region
white_tail_rects <- function(d, positions_per_row = 75) {
  wt_len <- unique(d\$max_pos)[1]
  if (!is.finite(wt_len)) return(dplyr::tibble()[0,])

  # if perfectly divisible by 75, there is no tail to cover
  if (wt_len %% positions_per_row == 0) return(dplyr::tibble()[0,])

  # which facet (row group) contains the last real position?
  last_group     <- ((wt_len - 1) %/% positions_per_row) + 1
  last_local_idx <- ((wt_len - 1) %% positions_per_row) + 1  # 1..75 within the facet

  tibble::tibble(
    row_group = last_group,
    xmin = last_local_idx + 0.5 - 0.025,     # tiny epsilon to avoid hairlines
    xmax = positions_per_row + 0.5 + 0.025,
    ymin = 0.5  - 0.025,
    ymax = 21.5 + 0.025
  )
}

plot_heatmap <- function(d, title_text, positions_per_row = 75) {
  amino_order <- rev(c("G", "A", "V", "L", "M", "I", "F",
                         "Y", "W", "K", "R", "H", "D", "E",
                         "S", "T", "C", "N", "Q", "P", "*"))
  d <- d %>% mutate(mut_aa = factor(mut_aa, levels = amino_order))

  min_f <- suppressWarnings(min(d\$fitness, na.rm = TRUE)); if (!is.finite(min_f)) min_f <- 0
  max_f <- suppressWarnings(max(d\$fitness, na.rm = TRUE)); if (!is.finite(max_f)) max_f <- 0
  max_orig_pos <- unique(d\$max_pos)[1]

  syn  <- syn_segments(d, positions_per_row)
  rect <- white_tail_rects(d, positions_per_row)

  ggplot(d, aes(x = wt_aa_pos, y = mut_aa, fill = fitness)) +
    scale_fill_gradientn(
      colours = c("#D73027", "#F0F0F0", "#4575B4"),
      values  = if ((abs(min_f) + max_f) > 0) c(0, abs(min_f)/(abs(min_f)+max_f), 1) else c(0, 0.5, 1),
      na.value = "grey35",
      limits   = c(min_f, max_f)
    ) +
    scale_x_discrete(
      labels = function(x) {
        num <- suppressWarnings(as.numeric(gsub("[^0-9]", "", x)))
        ifelse(num > max_orig_pos, " ", x)
      },
      expand = expansion(mult = c(0, 0))  # no extra margin area
    ) +
    geom_tile() +
    # Solid white block covering the tail (no pattern / no seams)
    { if (nrow(rect)) geom_rect(data = rect, inherit.aes = FALSE,
                                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                                fill = "white", color = NA) } +
    geom_segment(
      data = syn,
      aes(x = x - 0.485, xend = x + 0.485, y = y - 0.485, yend = y + 0.485),
      linewidth = 0.2, inherit.aes = FALSE, color = "grey10"
    ) +
    theme_minimal() +
    labs(title = title_text, x = "Wild-type amino acid", y = "Mutant amino acid", fill = "Fitness") +
    theme(
      plot.title   = element_text(size = 16, face = "bold"),
      axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
      axis.text.y  = element_text(size = 10),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
      panel.spacing = grid::unit(0.2, "lines")
    ) +
    facet_wrap(~ row_group, scales = "free_x", ncol = 1)
}

# ---------- main callable ----------
# fitness_table_path : path to fitness_estimation.tsv
# wt_seq_aa_txt_path : path to TXT file containing WT AA sequence (one line)
# output_pdf_path    : output PDF (default "fitness_heatmap.pdf")
# positions_per_row  : default 75
run_fitness_rescaled_heatmaps <- function(fitness_table_path,
                                          wt_seq_aa_txt_path,
                                          output_pdf_path = "fitness_heatmap.pdf",
                                          positions_per_row = 75) {

  df <- read.table(
    fitness_table_path, sep = "\t", header = TRUE,
    check.names = FALSE, quote = "", comment.char = ""
  )

  wt_aa_col     <- find_col(df, c("wt aa", "wt_aa", "wt"))
  pos_col       <- find_col(df, c("pos", "position"))
  mut_aa_col    <- find_col(df, c("mut aa", "mut_aa", "aa"))
  rescaled_cols <- get_rescaled_cols(df)

  wt_seq_aa <- read_wt_seq_aa_txt(wt_seq_aa_txt_path)

  plots <- list()

  ## 1) Mean first – use existing "mean fitness" column if available
  mean_col <- find_mean_col(df)
  if (is.null(mean_col)) {
    df\$`rescaled_fitness_mean` <- if (length(rescaled_cols) == 1) df[[rescaled_cols[1]]] else rowMeans(df[, rescaled_cols], na.rm = TRUE)
    mean_col <- "rescaled_fitness_mean"
  }
  long_df_mean <- build_heatmap_long(df, wt_aa_col, pos_col, mut_aa_col, mean_col,
                                     positions_per_row, wt_seq_aa)
  plots[[length(plots) + 1]] <- list(
    title = sprintf("Fitness — mean of %d replicate(s)", length(rescaled_cols)),
    data  = long_df_mean
  )

  ## 2) Then individual replicates
  for (i in seq_along(rescaled_cols)) {
    col <- rescaled_cols[i]
    long_df <- build_heatmap_long(df, wt_aa_col, pos_col, mut_aa_col, col,
                                  positions_per_row, wt_seq_aa)
    plots[[length(plots) + 1]] <- list(
      title = sprintf("Fitness — rep%d", i),
      data  = long_df
    )
  }

  # Device height: (#row groups × 4)
  page_heights <- vapply(plots, function(p) max(p\$data\$row_group, na.rm = TRUE), numeric(1))
  device_height <- max(4, as.numeric(page_heights) * 4, na.rm = TRUE)

  grDevices::pdf(output_pdf_path, width = 16, height = device_height)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  for (p in plots) print(plot_heatmap(p\$data, p\$title, positions_per_row))
  invisible(TRUE)
}

####
# run function
####
run_fitness_rescaled_heatmaps(
  fitness_table_path = "$fitness_estimation_tsv",
  wt_seq_aa_txt_path = "$wt_seq",
  output_pdf_path    = "fitness_heatmap.pdf"
)

####
# create versions.yml
####
r_version <- strsplit(version[['version.string']], ' ')[[1]][3]
dplyr_version <- as.character(packageVersion("dplyr"))
ggplot2_version <- as.character(packageVersion("ggplot2"))
methods_version <- as.character(packageVersion("methods"))
grid_version <- as.character(packageVersion("grid"))

if (is.null(r_version)) r_version <- "unknown"
if (length(dplyr_version) == 0) dplyr_version <- "unknown"
if (length(ggplot2_version) == 0) ggplot2_version <- "unknown"
if (length(methods_version) == 0) methods_version <- "unknown"
if (length(grid_version) == 0) grid_version <- "unknown"

f <- file("versions.yml", "w")
writeLines(
  c(
    '"${task.process}":',
    paste('    r-base:', r_version),
    paste('    r-dplyr:', dplyr_version),
    paste('    r-ggplot2:', ggplot2_version),
    paste('    r-methods:', methods_version),
    paste('    r-grid:', grid_version)
  ),
  f
)
close(f)
