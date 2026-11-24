#!/usr/bin/env Rscript

## mutscan-style variant enrichment (edgeR and limma-voom) for nf-core/deepmutscan
##
## Follows mutscan (Soneson C, Bendel AM, Diss G, Stadler MB. Genome Biology 2023, 24:132;
## https://github.com/fmicompbio/mutscan, MIT licence, (c) Charlotte Soneson and Michael Stadler, Friedrich
## Miescher Institute) step by step, without the mutscan package:
##
##   1. the per-variant table mutscan's summarizeExperiment() aggregates - counts summed over the
##      nucleotide variants that share an amino-acid level mutant name, rows in byte order of the name,
##      as the C++ mergeValues() leaves them;
##   2. calculateRelativeFC() with WTrows = "WT" and its default normMethod = "sum" - the wildtype
##      row's counts as per-sample offsets (edgeR::scaleOffset), then the edgeR quasi-likelihood fit
##      (estimateDisp, glmQLFit, glmQLFTest, predFC for the shrunken fold change) and the limma-voom
##      fit (voom with the offsets as library sizes, lmFit, eBayes, topTable with confidence intervals);
##   3. plotPairs() and plotVolcano() as drawn in mutscan (GGally::ggpairs with mutscan's correlation,
##      scatter and histogram panels; the volcano scatter centred on logFC = 0).
##
## The statistics are therefore edgeR's and limma's, called exactly as mutscan calls them, and the two
## result tables are identical to mutscan's for the same input - with one deliberate difference: a
## TTG or CTG variant codon at the first position is named L, as the genetic code reads it, where
## Biostrings' initiation-codon rule (which mutscan's names inherit) called it M; the functions in sections 2 and 3 are
## transcriptions of mutscan's R code (calculateRelativeFC.R, plotPairs.R, plotResults.R; mutscan
## 1.3.1 and 1.0.0 agree on them). What is not used is mutscan's read-processing core, which the
## pipeline never called: counts come from the pipeline's own tables.
##
## Inputs:  counts_merged (nt_seq, input1..N, output1..M), exp_design (DiMSum experimental design:
##          sample_name, experiment_replicate, selection_id), syn_wt_txt (the synonymous wildtype
##          proxy, one sequence), min_counts
## Needs:   edgeR, limma, ggplot2, GGally (translation is done with a codon table here)
## Outputs: fitness_estimation_mutscan_edgeR.tsv, fitness_estimation_mutscan_limma.tsv,
##          mutscan_counts_corr.pdf, mutscan_edgeR_volcano.pdf, mutscan_limma_volcano.pdf

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(GGally)
})

# The standard genetic code, for the translation of the variant ORFs (every sequence the pipeline
# hands to this step is built from the reference and the programmed codons, so it is plain ACGT).
GENETIC_CODE <- c(
  TTT = "F", TTC = "F", TTA = "L", TTG = "L", CTT = "L", CTC = "L", CTA = "L", CTG = "L",
  ATT = "I", ATC = "I", ATA = "I", ATG = "M", GTT = "V", GTC = "V", GTA = "V", GTG = "V",
  TCT = "S", TCC = "S", TCA = "S", TCG = "S", CCT = "P", CCC = "P", CCA = "P", CCG = "P",
  ACT = "T", ACC = "T", ACA = "T", ACG = "T", GCT = "A", GCC = "A", GCA = "A", GCG = "A",
  TAT = "Y", TAC = "Y", TAA = "*", TAG = "*", CAT = "H", CAC = "H", CAA = "Q", CAG = "Q",
  AAT = "N", AAC = "N", AAA = "K", AAG = "K", GAT = "D", GAC = "D", GAA = "E", GAG = "E",
  TGT = "C", TGC = "C", TGA = "*", TGG = "W", CGT = "R", CGC = "R", CGA = "R", CGG = "R",
  AGT = "S", AGC = "S", AGA = "R", AGG = "R", GGT = "G", GGC = "G", GGA = "G", GGG = "G"
)

# Translate equally long in-frame ORFs: one character matrix of amino acids, a row per sequence.
# Every codon is read by the standard genetic code, the first one included: a TTG or CTG in the first
# position is L, not the methionine that Biostrings::translate() returns by default for the
# alternative initiation codons (and that mutscan's mutant names inherited). In a deep mutational scan
# a variant codon at position 1 is a designed substitution, not a start codon.
translate_matrix <- function(seqs, orf_length) {
  starts <- seq(1, orf_length, by = 3)
  codons <- vapply(starts, function(i) substr(seqs, i, i + 2), character(length(seqs)))
  codons <- matrix(toupper(codons), nrow = length(seqs))
  aa <- GENETIC_CODE[codons]
  if (anyNA(aa)) stop("sequence with a codon outside ACGT cannot be translated")
  matrix(aa, nrow = length(seqs))
}

####################################################################################################
# 1. Variant table: amino-acid level mutant names and the aggregated count matrix
####################################################################################################

# One name per nucleotide variant: "WT" for synonymous variants of the wildtype proxy, otherwise the
# amino-acid changes as <wt><position><mutant>, joined with "_" if there are several. Vectorised over
# the table: one character matrix comparison instead of a strsplit() per row.
mutant_names <- function(nt_seq, wt_seq) {
  if (any(nchar(nt_seq) != nchar(wt_seq))) stop("variant sequences and the wildtype differ in length")
  wt_aa  <- translate_matrix(wt_seq, nchar(wt_seq))[1, ]
  aa_mat <- translate_matrix(nt_seq, nchar(wt_seq))
  differs <- aa_mat != matrix(wt_aa, nrow = length(nt_seq), ncol = length(wt_aa), byrow = TRUE)
  hits <- which(differs, arr.ind = TRUE)
  labels <- paste0(wt_aa[hits[, "col"]], hits[, "col"], aa_mat[hits])
  # one label per variant, the positions in increasing order
  o <- order(hits[, "row"], hits[, "col"])
  per_row <- split(labels[o], factor(hits[o, "row"], levels = seq_along(nt_seq)))
  name <- vapply(per_row, paste, character(1), collapse = "_")
  name[!nzchar(name)] <- "WT"
  unname(name)
}

# Counts per mutant name and sample: the nucleotide variants of one name summed (mutscan sums the
# duplicated (row, column) entries of its sparse matrix), rows in byte order of the name (the order
# of mutscan's std::map), columns in the experimental design's order.
aggregate_counts <- function(counts, mutant_name, sample_names) {
  rows <- sort(unique(mutant_name), method = "radix")        # C-locale byte order, like std::map
  idx  <- match(mutant_name, rows)
  mat  <- rowsum(as.matrix(counts[, sample_names, drop = FALSE]), idx, reorder = TRUE)
  mat  <- mat[match(seq_along(rows), as.integer(rownames(mat))), , drop = FALSE]
  dimnames(mat) <- list(rows, sample_names)
  storage.mode(mat) <- "double"
  mat
}

####################################################################################################
# 2. Relative fold changes (mutscan::calculateRelativeFC, WTrows = "WT", normMethod = "sum")
####################################################################################################

relative_fc <- function(counts, coldata, design, coef, WTrows, pseudocount = 1, method = "edgeR") {
  dge <- DGEList(counts = counts, samples = coldata)
  # size factors from the sum of the WT rows
  tmp0 <- dge\$counts[WTrows, , drop = FALSE]
  dge  <- scaleOffset(dge, log(colSums(tmp0)))
  if (method == "edgeR") {
    dge <- estimateDisp(dge, design = design)
    fit <- glmQLFit(dge, design = design)
    qlf <- glmQLFTest(fit, coef = coef)
    tt  <- topTags(qlf, n = Inf, sort.by = "none")\$table
    # shrunken fold changes
    predfc <- predFC(dge, design = design, prior.count = pseudocount)
    tt\$logFC_shrunk <- predfc[, coef]
    tt\$df.total <- qlf\$df.total
    tt\$df.prior <- qlf\$df.prior
    tt\$df.test  <- qlf\$df.test
  } else if (method == "limma") {
    # mutscan hands voom the DGEList with its offset and lib.size = exp(offset); whether voom reads
    # the offset or the lib.size argument depends on the limma version, so pass the library sizes
    # explicitly and nothing else - the same numbers either way
    vm  <- voom(modifyList(dge, list(offset = NULL)), design = design, lib.size = exp(dge\$offset))
    fit <- lmFit(vm, design = design)
    fit <- eBayes(fit)
    tt  <- topTable(fit, coef = coef, confint = TRUE, number = Inf, sort.by = "none")
    tt\$se.logFC <- sqrt(fit\$s2.post) * fit\$stdev.unscaled[, coef]
    tt\$df.total <- fit\$df.total
    tt\$df.prior <- fit\$df.prior
  }
  tt
}

####################################################################################################
# 3. Plots (mutscan::plotPairs with addIdentityLine = TRUE; mutscan::plotVolcano with pointSize = "large")
####################################################################################################

plot_pairs <- function(counts, pseudocount = 1, corMethod = "pearson", histBreaks = 40,
                       corSizeMult = 5, corSizeAdd = 2, pointSize = 0.1, pointAlpha = 0.3) {
  ggtheme <- list(theme_bw(), theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
  mat <- as.data.frame(log10(counts + pseudocount))
  title <- paste0("log10(counts + ", pseudocount, "), ",
                  toupper(substr(corMethod, 1, 1)), substr(corMethod, 2, nchar(corMethod)), " correlation")
  allCors  <- cor(mat, method = corMethod, use = "pairwise.complete.obs")
  corRange <- range(abs(allCors[upper.tri(allCors)]))
  if (length(unique(corRange)) == 1) corRange <- c(0, 1)
  cor_fcn <- function(data, mapping, ...) {
    xData <- eval_data_col(data, mapping\$x)
    yData <- eval_data_col(data, mapping\$y)
    mainCor <- cor(xData, yData, method = corMethod, use = "pairwise.complete.obs")
    transfCor <- (abs(mainCor) - min(corRange)) / (max(corRange) - min(corRange))
    transfCor <- max(0.0, min(1.0, transfCor))
    cols <- if (mainCor >= 0) hcl.colors(n = 11, palette = "RdBu")[6:2] else hcl.colors(n = 11, palette = "RdBu")[6:10]
    col  <- rgb(colorRamp(cols)(transfCor), maxColorValue = 255)
    ggplot(data = data, mapping = mapping) +
      annotate(x = 0.5, y = 0.5, label = round(mainCor, digits = 3), geom = "text",
               size = abs(mainCor) * corSizeMult + corSizeAdd) +
      ggtheme + ylim(c(0, 1)) +
      theme(panel.background = element_rect(fill = col))
  }
  plotpoints <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
      geom_point(alpha = pointAlpha, size = pointSize) +
      ggtheme
  }
  diaghist <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
      geom_histogram(fill = "#F5C710", color = "grey50", bins = histBreaks) +
      ggtheme
  }
  ggpairs(data = mat, mapping = NULL, columns = seq_len(ncol(mat)), title = title, xlab = NULL, ylab = NULL,
          upper = list(continuous = cor_fcn), lower = list(continuous = plotpoints),
          diag = list(continuous = diaghist), progress = FALSE, axisLabels = "show")
}

plot_volcano <- function(res, pvalCol, padjCol, padjThreshold = 0.05) {
  res\$negLog10P <- -log10(res[[pvalCol]])
  res\$feature   <- rownames(res)
  xr <- range(res\$logFC, na.rm = TRUE)
  xr <- c(-max(abs(xr), na.rm = TRUE), max(abs(xr), na.rm = TRUE))
  yr <- range(res\$negLog10P, na.rm = TRUE)
  gg <- ggplot(res, aes(x = .data[["logFC"]], y = .data[["negLog10P"]], label = .data[["feature"]])) +
    theme_minimal() +
    coord_cartesian(xlim = xr, ylim = yr) +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14), title = element_text(size = 14)) +
    labs(x = "logFC", y = paste0("-log10(", pvalCol, ")")) +
    geom_point(fill = "lightgrey", color = "grey", pch = 21, size = 1.5)
  if (any(res[[padjCol]] < padjThreshold)) {
    gg <- gg +
      labs(subtitle = paste0("Points are indicated in red if ", padjCol, "<", padjThreshold)) +
      geom_point(data = res[res[[padjCol]] < padjThreshold, , drop = FALSE],
                 fill = "red", color = "grey", pch = 21, size = 1.5)
  }
  gg
}

####################################################################################################
# 4. Run
####################################################################################################

run_mutscan_fitness_estimation <- function(counts_path, design_path, wt_seq_path,
                                           output_path_edgeR, output_path_limma, min_counts = 10) {

  ## 1. Import key files
  merged.counts <- read.table(counts_path, sep = "\\t", header = TRUE, check.names = FALSE)
  exp.design    <- read.table(design_path, sep = "\\t", header = TRUE, check.names = FALSE)
  wt.seq        <- toupper(trimws(readLines(wt_seq_path, warn = FALSE)[1]))

  ## keep variants with at least `min_counts` reads in every input replicate (and the wildtype),
  ## so that the same variants are analysed as by the default estimator and DiMSum
  input.cols <- grep("^input[0-9]+\$", colnames(merged.counts), value = TRUE)
  keep <- apply(merged.counts[, input.cols, drop = FALSE] >= min_counts, 1, all)
  keep[is.na(keep)] <- FALSE
  keep[merged.counts[, 1] == wt.seq] <- TRUE
  merged.counts <- merged.counts[keep, , drop = FALSE]

  ## 2. Count matrix per amino-acid mutant (mutscan's summarizeExperiment) and the column data
  sample_names <- as.character(exp.design\$sample_name)
  counts <- aggregate_counts(merged.counts, mutant_names(merged.counts[, 1], wt.seq), sample_names)
  condition <- ifelse(exp.design\$selection_id == 0, "input", "output")
  coldata <- data.frame(Name = sample_names, Condition = condition,
                        Replicate = factor(as.character(exp.design\$experiment_replicate)),
                        stringsAsFactors = FALSE)
  rownames(coldata) <- sample_names

  ## 3. Model: replicate as a blocking factor, the input -> output change as the coefficient
  model.design <- model.matrix(~ Replicate + Condition, data = coldata)

  ## 4. logFC estimates (mutscan's calculateRelativeFC with the WT row as the normalisation)
  logFC.edgeR <- relative_fc(counts, coldata, model.design, coef = "Conditionoutput", WTrows = "WT",
                             pseudocount = 1, method = "edgeR")
  logFC.limma <- relative_fc(counts, coldata, model.design, coef = "Conditionoutput", WTrows = "WT",
                             pseudocount = 1, method = "limma")

  ## 5. QC plots
  pdf("mutscan_counts_corr.pdf", height = 9, width = 14)
  print(plot_pairs(counts))
  dev.off()
  pdf("mutscan_edgeR_volcano.pdf", height = 9, width = 14)
  print(plot_volcano(logFC.edgeR, pvalCol = "PValue", padjCol = "FDR"))
  dev.off()
  pdf("mutscan_limma_volcano.pdf", height = 9, width = 14)
  print(plot_volcano(logFC.limma, pvalCol = "P.Value", padjCol = "adj.P.Val"))
  dev.off()

  ## 6. Data export: the mutant name (edgeR's and limma's row name) as the first column - the
  ##    previous template wrote the tables without it, which left the rows unidentifiable
  write.table(cbind(mutant = rownames(logFC.edgeR), logFC.edgeR), output_path_edgeR,
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\\t", na = "")
  write.table(cbind(mutant = rownames(logFC.limma), logFC.limma), output_path_limma,
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\\t", na = "")

  invisible(list(edgeR = logFC.edgeR, limma = logFC.limma))
}

#####
# run function
#####
run_mutscan_fitness_estimation(
  counts_path = "$counts_merged",
  design_path = "$exp_design",
  wt_seq_path = "$syn_wt_txt",
  output_path_edgeR = "fitness_estimation_mutscan_edgeR.tsv",
  output_path_limma = "fitness_estimation_mutscan_limma.tsv",
  min_counts = as.numeric("$min_counts")
)

####
# create versions.yml
####
r_version <- strsplit(version[['version.string']], ' ')[[1]][3]
edger_version <- as.character(packageVersion("edgeR"))
limma_version <- as.character(packageVersion("limma"))
ggplot2_version <- as.character(packageVersion("ggplot2"))
ggally_version <- as.character(packageVersion("GGally"))

if (is.null(r_version)) r_version <- "unknown"
if (length(edger_version) == 0) edger_version <- "unknown"
if (length(limma_version) == 0) limma_version <- "unknown"
if (length(ggplot2_version) == 0) ggplot2_version <- "unknown"
if (length(ggally_version) == 0) ggally_version <- "unknown"

f <- file("versions.yml", "w")
writeLines(
  c(
    '"${task.process}":',
    paste('    r-base:', r_version),
    paste('    r-edgeR:', edger_version),
    paste('    r-limma:', limma_version),
    paste('    r-ggplot2:', ggplot2_version),
    paste('    r-GGally:', ggally_version)
  ),
  f
)
close(f)
