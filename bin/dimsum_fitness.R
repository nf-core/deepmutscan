#!/usr/bin/env Rscript

# dimsum_fitness.R - DiMSum's STEAM stages as one standalone script
#
# DiMSum (Faure AJ, Schmiedel JM, Baeza-Centurion P, Lehner B. Genome Biology 2020, 21:207;
# https://github.com/lehner-lab/DiMSum) runs in two halves: WRAP (stages 0-3: demultiplexing, FastQC,
# trimming, read merging and tallying with cutadapt, vsearch and starcode) and STEAM (stage 4: variant
# processing and filtering; stage 5: fitness and error model). nf-core/deepmutscan produces the variant
# count table itself and only ever ran STEAM (`--countPath --startStage 4 --stopStage 5`). This file is
# that path, lifted from DiMSum 1.4 (commit a441cdf): every DiMSum function the
# STEAM path reaches is below, verbatim, followed by a driver that does what dimsum() does when a count
# file is given. The outputs - the variant tables, fitness_*.txt, the MaveDB export, the .RData files
# and report.html - are the same files with the same content.
#
# What is different from the package, and why:
#   - no DiMSum package: the version string, the report R markdown and the logo are passed in;
#   - the Biostrings, Cairo, IRanges, ShortRead and stringr dependencies are gone: STEAM used
#     Biostrings for translation, reverse complement and the three-letter amino-acid code (base-R
#     versions below) and Cairo for PNG devices (grDevices::png(type = "cairo")); the other three were
#     only used by WRAP. What remains is data.table, ggplot2, GGally, hexbin, reshape2, plyr, seqinr,
#     cowplot, gridExtra, optparse, mgcv and rmarkdown (with pandoc) - all of them with linux-aarch64
#     builds;
#   - it runs on ggplot2 4 (DiMSum 1.4 was written for 3.5): the four ggplot2 changes that would move
#     marks in its figures are undone below, so that the figures are those of ggplot2 3.5;
#   - the per-variant Hamming distances are a byte comparison instead of one strsplit() per variant,
#     the count-histogram jitter is drawn in one vectorised step, and the report is rendered once, at
#     the end, instead of also before stage 4 (DiMSum renders a preliminary report so that users of a
#     long WRAP run have something to look at);
#   - figures written as both PNG and PDF are laid out once and drawn into both files (the PNGs are
#     unchanged; the PDFs use the PNG's layout), and R's just-in-time compiler is off;
#   - error model: on one core the bootstraps run in this process instead of a one-worker cluster,
#     from the same random-number stream, and the per-row variances call the C routine of var()
#     directly - the same numbers, in a fraction of the time;
#   - the command line is DiMSum's (the same option names and defaults), plus `--reportRmd` and
#     `--logoPng` for the two assets.
# Every edit of DiMSum's code is recorded in the generator (make_dimsum_fitness.py) and applied to the
# exact upstream lines, so that this file can be re-derived from a newer DiMSum release.
#
# DiMSum is MIT-licensed:
#
#   MIT License
#   Copyright (c) 2020 Lehner Lab
#   Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
#   associated documentation files (the "Software"), to deal in the Software without restriction,
#   including without limitation the rights to use, copy, modify, merge, publish, distribute,
#   sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
#   furnished to do so, subject to the following conditions:
#   The above copyright notice and this permission notice shall be included in all copies or
#   substantial portions of the Software.
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
#   NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES
#   OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
#   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

suppressPackageStartupMessages(library(data.table))
# DiMSum calls reshape2::melt() on data.tables and relies on getting data.tables back. Inside the
# package that works because its namespace imports data.table, so S3 dispatch finds melt.data.table;
# from a script it does not (R looks past the search path), so the methods are registered explicitly.
registerS3method("melt", "data.table", data.table::melt.data.table, envir = asNamespace("reshape2"))
registerS3method("dcast", "data.table", data.table::dcast.data.table, envir = asNamespace("reshape2"))
# Every function below runs once or a few times per call of the script; compiling each of them to byte
# code on its first calls (R's default) costs more than it saves. Package code is compiled already.
invisible(compiler::enableJIT(0))

####################################################################################################
# Figures: one layout, two files
####################################################################################################
# DiMSum writes most figures twice, as PNG (dimsum__save_png) and PDF (ggplot2::ggsave), and each call
# lays the figure out from scratch - for the all-vs-all scatter-plot matrices most of the run time. Here
# the layout (the gtable) is built once, while the PNG device is open, as print() inside
# dimsum__save_png builds it, and drawn into both files: the PNGs are byte-identical to DiMSum's, the
# PDFs show the same figure laid out with the PNG device's font metrics (sub-point shifts).
dimsum_fitness__save_png_pdf <- function(
  png_file,
  pdf_file,
  plot_object,
  width = 7,
  height = 7,
  plot_dpi = 300){
  #Open the PNG device (as dimsum__save_png)
  grDevices::png(
    png_file,
    width = width,
    height = height,
    units = "in",
    res = plot_dpi,
    bg = "white",
    type = "cairo")
  if(isTRUE(names(plot_object)[1]=="grobs")){
    #Already a gtable (gridExtra::arrangeGrob): drawn as dimsum__save_png draws it
    g <- plot_object
    plot(g)
  }else{
    #ggplot or GGally ggmatrix: the layout print() would build, built once
    if(inherits(plot_object, "ggmatrix")){
      g <- GGally::ggmatrix_gtable(plot_object)
    }else{
      g <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(plot_object))
    }
    grid::grid.newpage()
    grid::grid.draw(g)
  }
  #Close device
  { sink("/dev/null")
    dev.off()
    sink()
  }
  #The same layout into the PDF
  ggplot2::ggsave(pdf_file, g, width = width, height = height)
}

# Evaluate `expr` with the random-number state that parallel::clusterSetRNGStream(cl, iseed) gives the
# first (here: only) worker of a cluster - L'Ecuyer-CMRG, set.seed(iseed) - and restore the caller's
# state afterwards, as clusterSetRNGStream() restores it in the main process
dimsum_fitness__with_worker_stream <- function(iseed, expr){
  oldseed <- if(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit(
    if(!is.null(oldseed)){
      assign(".Random.seed", oldseed, envir = .GlobalEnv)
    }else if(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)){
      rm(".Random.seed", envir = .GlobalEnv)
    })
  RNGkind("L'Ecuyer-CMRG")
  set.seed(iseed)
  expr
}

DIMSUM_FITNESS_VERSION <- "1.4 (STEAM stages, nf-core/deepmutscan)"

####################################################################################################
# ggplot2 4: DiMSum's figures as ggplot2 3.5 drew them
####################################################################################################

# DiMSum 1.4 was written for ggplot2 3.5; the pipeline's R environment has ggplot2 4. Four of the
# changes move marks in DiMSum's figures, and are undone here:
#   - geom_abline() now ends the line at the panel's y range too, not only at its x range: the line is
#     the same, but a dashed line's dashes start elsewhere. dimsum__geom_abline() draws it from the
#     end points ggplot2 3.5 used;
#   - a log10 scale's minor breaks were taken from the major breaks of the expanded range before they
#     were cut to that range, so the minor break beyond the last major one sat halfway (in log space)
#     to the next major break; ggplot2 4 cuts first and extrapolates. dimsum__scale_x_log10() and
#     dimsum__scale_y_log10() compute them as ggplot2 3.5 did;
#   - a density's legend key is now drawn with the curve's round line joins (3.5: mitred corners):
#     dimsum__draw_key_polygon() mitres them again;
#   - text is now 11 / .pt = 3.87 mm by default (from the theme's font size), not 3.88 mm.
# (ggplot objects are also no longer lists, which the plot-saving functions above account for.)
dimsum__GeomAbline <- ggplot2::ggproto("dimsum__GeomAbline", ggplot2::GeomAbline,
  draw_panel = function(data, panel_params, coord, lineend = "butt"){
    ranges <- coord$backtransform_range(panel_params)
    if(coord$clip == "on" && coord$is_linear()){
      ranges$x <- ranges$x + c(-1, 1) * diff(ranges$x)
    }
    data$x <- ranges$x[1]
    data$xend <- ranges$x[2]
    data$y <- ranges$x[1] * data$slope + data$intercept
    data$yend <- ranges$x[2] * data$slope + data$intercept
    ggplot2::GeomSegment$draw_panel(unique(data), panel_params, coord, lineend = lineend)
  })
dimsum__geom_abline <- function(...){
  layer <- ggplot2::geom_abline(...)
  layer$geom <- dimsum__GeomAbline
  layer
}
dimsum__minor_breaks_log10 <- function(limits, breaks){
  10^scales::regular_minor_breaks()(log10(scales::breaks_log(base = 10)(limits)), log10(limits), 2)
}
dimsum__scale_x_log10 <- function(...){ggplot2::scale_x_log10(..., minor_breaks = dimsum__minor_breaks_log10)}
dimsum__scale_y_log10 <- function(...){ggplot2::scale_y_log10(..., minor_breaks = dimsum__minor_breaks_log10)}
dimsum__draw_key_polygon <- function(data, params, size){
  params$linejoin <- "mitre"
  ggplot2::draw_key_polygon(data, params, size)
}
ggplot2::update_geom_defaults("text", list(size = 3.88))

####################################################################################################
# Base-R replacements for the Biostrings calls
####################################################################################################

# The standard genetic code; stop codons are "*". Codons with a letter outside ACGT translate to "X".
DIMSUM_FITNESS_GENETIC_CODE <- c(
  TTT = "F", TTC = "F", TTA = "L", TTG = "L", CTT = "L", CTC = "L", CTA = "L", CTG = "L",
  ATT = "I", ATC = "I", ATA = "I", ATG = "M", GTT = "V", GTC = "V", GTA = "V", GTG = "V",
  TCT = "S", TCC = "S", TCA = "S", TCG = "S", CCT = "P", CCC = "P", CCA = "P", CCG = "P",
  ACT = "T", ACC = "T", ACA = "T", ACG = "T", GCT = "A", GCC = "A", GCA = "A", GCG = "A",
  TAT = "Y", TAC = "Y", TAA = "*", TAG = "*", CAT = "H", CAC = "H", CAA = "Q", CAG = "Q",
  AAT = "N", AAC = "N", AAA = "K", AAG = "K", GAT = "D", GAC = "D", GAA = "E", GAG = "E",
  TGT = "C", TGC = "C", TGA = "*", TGG = "W", CGT = "R", CGC = "R", CGA = "R", CGG = "R",
  AGT = "S", AGC = "S", AGA = "R", AGG = "R", GGT = "G", GGC = "G", GGA = "G", GGG = "G"
)

# Biostrings::AMINO_ACID_CODE for the twenty standard amino acids
DIMSUM_FITNESS_AMINO_ACID_CODE <- c(
  G = "Gly", A = "Ala", V = "Val", L = "Leu", M = "Met", I = "Ile", F = "Phe", Y = "Tyr", W = "Trp",
  K = "Lys", R = "Arg", H = "His", D = "Asp", E = "Glu", S = "Ser", T = "Thr", C = "Cys", N = "Asn",
  Q = "Gln", P = "Pro"
)

# as.character(Biostrings::translate(DNAStringSet(x), no.init.codon = TRUE)): case-insensitive, the
# trailing one or two bases of a sequence whose length is not a multiple of three are ignored
dimsum_fitness__translate <- function(x) {
  x <- toupper(x)
  n <- nchar(x) %/% 3L
  out <- character(length(x))
  for (len in unique(n)) {
    i <- which(n == len)
    if (len == 0L) { out[i] <- ""; next }
    starts <- seq.int(1L, by = 3L, length.out = len)
    codons <- vapply(starts, function(s) substr(x[i], s, s + 2L), character(length(i)))
    aa <- DIMSUM_FITNESS_GENETIC_CODE[codons]
    aa[is.na(aa)] <- "X"
    out[i] <- apply(matrix(aa, nrow = length(i)), 1, paste0, collapse = "")
  }
  out
}

# as.character(Biostrings::reverseComplement(DNAStringSet(x))), for A/C/G/T in either case
dimsum_fitness__reverse_complement <- function(x) {
  comp <- chartr("ACGTacgt", "TGCAtgca", x)
  vapply(strsplit(comp, "", fixed = TRUE), function(s) paste0(rev(s), collapse = ""), character(1))
}

# mapply(dimsum__hamming_distance, x, ref) for equally long strings, as a byte comparison (15x faster
# than one strsplit() per sequence on a 1 kb ORF; the same integers)
dimsum_fitness__hamming_distances <- function(x, ref) {
  if (length(x) == 0) return(integer(0))
  if (length(ref) != 1) return(mapply(dimsum__hamming_distance, x, ref))
  r <- charToRaw(ref)
  vapply(x, function(s) sum(charToRaw(s) != r), integer(1), USE.NAMES = FALSE)
}

####################################################################################################
# DiMSum functions (verbatim, see the header)
####################################################################################################

# ---- dimsum__add_dropout_pseudocount
dimsum__add_dropout_pseudocount <- function(
  dimsum_meta,
  input_dt,
  all_reps,
  verbose = T
  ){

  #Skip pseudocounts if not required
  if(dimsum_meta[["fitnessDropoutPseudocount"]]==0){
    return(input_dt)
  }

  if(verbose){dimsum__status_message("Adding pseudocount to output samples with dropout...\n")}

  #Add psuedocounts for all output replicates
  for(j in all_reps){
    input_dt[get(paste0("count_e", j, "_s0"))>0 & get(paste0("count_e", j, "_s1"))==0, paste0("count_e", j, "_s1") := .SD[[1]] + dimsum_meta[["fitnessDropoutPseudocount"]],,.SDcols = paste0("count_e", j, "_s1")]
  }

  if(verbose){dimsum__status_message("Done\n")}

  return(input_dt)

}

# ---- dimsum__aggregate_AA_variants
dimsum__aggregate_AA_variants <- function(
  dimsum_meta,
  input_dt,
  all_reps
  ){

  dimsum__status_message("Aggregating counts for identical amino acid variants...\n")

  #Number of input and output replicates
  all_reps_str <- paste0(all_reps, collapse="|")

  #Sample names
  input_samples <- names(input_dt)[grep(paste0("e(", all_reps_str, ")_s0_b.*_count$"), names(input_dt))]
  output_samples <- names(input_dt)[grep(paste0("e(", all_reps_str, ")_s1_b.*_count$"), names(input_dt))]

  dimsum__check_variants(dimsum_meta = dimsum_meta, input_dt = input_dt)

  #Set merge_seq to nt_seq
  input_dt[,merge_seq := aa_seq,nt_seq]

  #Reference AA sequences for synonymous sequences
  synseq_refs <- unique(c(dimsum_meta[["synonym_sequences"]], input_dt[WT==T,aa_seq]))

  #Pick out reference-like variants (without nucleotide indels)
  syn_dt <- input_dt[aa_seq %in% synseq_refs & nchar(nt_seq)==(3*nchar(aa_seq))]
  syn_dt[,merge_seq := nt_seq,nt_seq]
  syn_dt[, error_model := F]
  syn_dt <- syn_dt[,.SD,merge_seq,.SDcols = c(
    "nt_seq","aa_seq","Nham_nt","Nham_aa",
    "Nmut_codons","WT","indel","STOP","STOP_readthrough","error_model",
    c(input_samples, output_samples))]

  #For all count columns
  for(j in c(input_samples, output_samples)){
    #Aggregate counts accross identical AA variants
    input_dt[!is.na(get(j)),paste0(j,"_agg") := sum(.SD, na.rm = T),merge_seq,.SDcols = j]
    #Set all synonyms to the same count
    input_dt[,paste0(j,"_agg") := unique(na.omit(.SD)),merge_seq,.SDcols = paste0(j,"_agg")]
  }

  #Retain only one row per AA variant
  output_dt <- input_dt[!duplicated(merge_seq),.SD,merge_seq,.SDcols = c(
    "nt_seq","aa_seq","Nham_nt","Nham_aa",
    "Nmut_codons","WT","indel","STOP","STOP_readthrough","error_model",
    names(input_dt)[grep(names(input_dt),pattern="_agg$")])]

  #Revert to original names of aggregated count columns
  names(output_dt)[grep(names(output_dt),pattern="_agg$")] <- gsub("_agg$", "", names(output_dt)[grep(names(output_dt),pattern="_agg$")])

  #Nham_nt and Nmut_codons columns are meaningless after aggregation at AA level
  output_dt[, Nham_nt := NA]
  output_dt[, Nmut_codons := NA]

  #Define WT at AA level
  output_dt[Nham_aa==0, WT := T]

  #nt_seq column is meaningless after aggregation
  output_dt[, nt_seq := NA]

  #Add back silent/synonymous variants
  output_dt <- rbind(syn_dt, output_dt)

  dimsum__status_message("Done\n")

  return(output_dt)

}

# ---- dimsum__aggregate_AA_variants_fitness
dimsum__aggregate_AA_variants_fitness <- function(
  dimsum_meta,
  input_dt,
  all_reps
  ){

  dimsum__status_message("Aggregating fitness for identical amino acid variants...\n")

  #Number of input and output replicates
  all_reps_str <- paste0(all_reps, collapse="|")

  #Sample names
  input_samples <- names(input_dt)[grep(paste0("count_e(", all_reps_str, ")_s0$"), names(input_dt))]
  output_samples <- names(input_dt)[grep(paste0("count_e(", all_reps_str, ")_s1$"), names(input_dt))]

  dimsum__check_variants(dimsum_meta = dimsum_meta, input_dt = input_dt)

  #Set merge_seq to nt_seq
  input_dt[,merge_seq := aa_seq,nt_seq]

  #Set merge_seq to nt_seq if WT or silent/synonymous variants
  input_dt[Nham_aa==0,merge_seq := nt_seq,nt_seq]

  #Reference AA sequences for synonymous sequences
  synseq_refs <- unique(c(dimsum_meta[["synonym_sequences"]], input_dt[WT==T,aa_seq]))

  #Pick out reference-like variants
  syn_dt <- input_dt[aa_seq %in% synseq_refs & nchar(nt_seq)==(3*nchar(aa_seq))]
  syn_dt[,merge_seq := nt_seq,nt_seq]
  syn_dt[, error_model := F]
  syn_dt <- syn_dt[,.SD,merge_seq,.SDcols = c(
    "aa_seq","Nham_nt","Nham_aa",
    "Nmut_codons","WT","indel","STOP","STOP_readthrough","error_model",
    "mean_count",
    c(input_samples, output_samples),
    names(syn_dt)[grep(names(syn_dt),pattern="^fitness|^sigma")])]

  #For all count columns
  for(j in c(input_samples, output_samples)){
    #Aggregate counts accross identical AA variants
    input_dt[!is.na(get(j)),paste0(j,"_agg") := sum(.SD, na.rm = T),merge_seq,.SDcols = j]
    input_dt[,paste0(j,"_agg") := unique(na.omit(.SD)),merge_seq,.SDcols = paste0(j,"_agg")]
  }

  #Calculate mean input counts
  input_dt[,mean_count := rowMeans(.SD, na.rm = T),,.SDcols = paste0("count_e", all_reps, "_s0_agg")]

  #For all fitness columns
  idx <- names(input_dt)[grep(names(input_dt),pattern="^fitness")]
  for (i in seq_along(idx)) {
    #Aggregate fitness accross identical AA variants
    input_dt[!is.na(get(idx[i])) & !is.na(get(gsub("fitness", "sigma", idx[i]))),paste0(idx[i],"_agg") := sum(.SD[[1]]/(.SD[[2]]^2), na.rm = T)/sum(1/(.SD[[2]]^2), na.rm = T),merge_seq,.SDcols = c(idx[i], gsub("fitness", "sigma", idx[i]))]
    #Set all synonyms to the same fitness
    input_dt[,paste0(idx[i],"_agg") := unique(na.omit(.SD)),merge_seq,.SDcols = paste0(idx[i],"_agg")]
    #Aggregate fitness error accross identical AA variants
    input_dt[!is.na(get(gsub("fitness", "sigma", idx[i]))),paste0(gsub("fitness", "sigma", idx[i]),"_agg") := sqrt(1/sum(1/(.SD[[1]]^2), na.rm = T)),merge_seq,.SDcols = gsub("fitness", "sigma", idx[i])]
    #Set all synonyms to the same fitness error
    input_dt[,paste0(gsub("fitness", "sigma", idx[i]),"_agg") := unique(na.omit(.SD)),merge_seq,.SDcols = paste0(gsub("fitness", "sigma", idx[i]),"_agg")]
  }

  #Retain only one row per AA variant
  output_dt <- input_dt[!duplicated(merge_seq),.SD,merge_seq,.SDcols = c(
    "aa_seq","Nham_nt","Nham_aa",
    "Nmut_codons","WT","indel","STOP","STOP_readthrough","error_model",
    "mean_count",
    names(input_dt)[grep(names(input_dt),pattern="_agg$")])]

  #Revert to original names of aggregated count columns
  names(output_dt)[grep(names(output_dt),pattern="_agg$")] <- gsub("_agg$", "", names(output_dt)[grep(names(output_dt),pattern="_agg$")])

  #Nham_nt column is meaningless for nonsynonymous variants after aggregation
  output_dt[Nham_aa!=0, Nham_nt := NA]
  output_dt[Nham_aa!=0, Nmut_codons := NA]

  #Add back silent/synonymous variants
  output_dt <- rbind(syn_dt, output_dt)

  dimsum__status_message("Done\n")

  return(output_dt)

}

# ---- dimsum__bayesian_double_fitness [modified, see the header and the generator]
dimsum__bayesian_double_fitness <- function(
  dimsum_meta,
  doubles_dt,
  singles_dt,
  wt_dt,
  all_reps,
  report = TRUE,
  report_outpath = NULL
  ){

  dimsum__status_message("Calculating Bayesian double mutant fitness estimates...\n")

  #Poisson distribution for score likelihood
  lam_d <- dimsum_meta[["bayesianDoubleFitnessLamD"]]

  #Number of input and output replicates
  all_reps_str <- paste0(all_reps, collapse="|")

  #Plot fitness densities for different mean count bins (replicate 1)
  doubles_dt[,fitness_for_bins := .SD[[1]],,.SDcols = paste0("fitness",all_reps[1],"_uncorr")]
  if(report){
    d <- ggplot2::ggplot(doubles_dt[between(bin_count,2,8) & !is.infinite(fitness_for_bins)],ggplot2::aes(fitness_for_bins,..density..,color=factor(bin_count))) +
      ggplot2::geom_density(adjust=1, key_glyph = dimsum__draw_key_polygon)
    dimsum__save_png(file.path(report_outpath, "dimsum_stage_fitness_report_4_doubles_bayesian_framework1.png"), d, width = 7, height = 5)
  }

  #Plot fitness densities for mean counts greater/less than 50 (replicate 1)
  if(report){
    d <- ggplot2::ggplot(doubles_dt[!is.infinite(fitness_for_bins) & !is.na(fitness_for_bins)],ggplot2::aes(fitness_for_bins,..density..,color=counts_for_bins >= dimsum_meta[["fitnessDoubleHighConfidenceCount"]])) +
      ggplot2::geom_density(adjust=1, key_glyph = dimsum__draw_key_polygon)
    dimsum__save_png(file.path(report_outpath, "dimsum_stage_fitness_report_4_doubles_bayesian_framework2.png"), d, width = 7, height = 5)
  }
  #>> try to estimate what fitness scores are for variants with low sequence coverage
  # use double mutants with variants >= dimsum_meta[["fitnessDoubleHighConfidenceCount"]] counts

  ## calculate posterior double mutant fitness based on prior from single mutants
  postpois_conditioned_singleF <- function(i){
    require(data.table)
    count_in <- double_data[i,count_in]
    count_out <- double_data[i,count_out]
    lam_in <- exp(seq(floor(log(count_in+0.1)-max(c(0.5,1/log10(count_in+1.75)))),(log(count_in+0.1)+max(c(0.5,1/log10(count_in+1.75)))),lam_d))
    lam_out <- exp(seq(floor(log(count_out+0.1)-max(c(0.5,1/log10(count_out+1.75)))),(log(count_out+0.1)+max(c(0.5,1/log10(count_out+1.75)))),lam_d))
    lam_low <- range(log(lam_out))[1] - range(log(lam_in))[2]
    lam_high <- range(log(lam_out))[2] - range(log(lam_in))[1]
    idx <- row(matrix(NA,nrow=length(lam_out),ncol=length(lam_in))) - col(matrix(NA,nrow=length(lam_out),ncol=length(lam_in)))
    likelihood <- sapply(split(outer(dpois(count_out,lambda = lam_out),dpois(count_in,lambda = lam_in)),idx),sum)
    score_prior <- density(score_prior_cond[,.(fdist = sqrt((double_data[i,F1]-F1)^2+(double_data[i,F2]-F2)^2),F)][
      order(fdist)][1:Nneighbours,F],
      from = (lam_low-wt_corr),
      to = (lam_high-wt_corr),
      n = as.integer(as.character(round((lam_high-lam_low)/lam_d + 1)))) #super weird bug
    posterior <- score_prior$y*likelihood

    moments <- list()
    moments[1] <- weighted.mean(x = score_prior$x,w = posterior)
    moments[2] <- sqrt(sum(( moments[[1]]-score_prior$x)^2 * posterior)/
                        sum(posterior))
    return(moments)
  }

  # Setup cluster
  clust <- parallel::makeCluster(dimsum_meta[['numCores']]) #This line will take time

  #Calculate conditional fitness and sigma
  for (E in all_reps) {
    #wildtype "correction" to calculate scores
    wt_corr <- wt_dt[,log(unlist(.SD[,1]) / unlist(.SD[,2])),,.SDcols = c(paste0("count_e",E,"_s1"),paste0("count_e",E,"_s0"))]
    #data for prior calculation
    double_data <- doubles_dt[,.(Pos1,Mut1,Pos2,Mut2,count_in = unlist(.SD[,1]),count_out = unlist(.SD[,2]),
                             F = unlist(.SD[,3])),,
                          .SDcols = c(paste0("count_e",E,"_s0"),paste0("count_e",E,"_s1"),paste0("fitness",E,"_uncorr"))]
    # double_data = merge(double_data,singles_dt[,.(Pos,Mut,F1 = .SD),,.SDcols = paste0("fitness",E)],by.x = c("Pos1","Mut1"),by.y = c("Pos","Mut"))
    double_data <- merge(double_data,singles_dt[!is.na(singles_dt[,paste0("fitness",E)]) & mean_count >= dimsum_meta[["fitnessHighConfidenceCount"]],.(Pos,Mut,F1 = .SD),,.SDcols = paste0("fitness",E)],by.x = c("Pos1","Mut1"),by.y = c("Pos","Mut"))
    # double_data = merge(double_data,singles_dt[,.(Pos,Mut,F2 = .SD),,.SDcols = paste0("fitness",E)],by.x = c("Pos2","Mut2"),by.y = c("Pos","Mut"))
    double_data <- merge(double_data,singles_dt[!is.na(singles_dt[,paste0("fitness",E)]) & mean_count >= dimsum_meta[["fitnessHighConfidenceCount"]],.(Pos,Mut,F2 = .SD),,.SDcols = paste0("fitness",E)],by.x = c("Pos2","Mut2"),by.y = c("Pos","Mut"))

    Nneighbours <- 500
    score_prior_cond <- double_data[count_in >= dimsum_meta[["fitnessDoubleHighConfidenceCount"]] & F > -Inf & F1 > -Inf & F2 > -Inf]

    # make variables available to each core's workspace
    parallel::clusterExport(clust, list("double_data","lam_d","wt_corr","score_prior_cond","Nneighbours"), envir = environment())

    #posterior fitness conditioned on single fitness
    t<-proc.time()
    helper <- parallel::parSapply(clust,X = 1:nrow(double_data), postpois_conditioned_singleF)
    print(proc.time()-t)
    helper1 <- matrix(unlist(helper),nrow=2)
    double_data[,paste0("fitness",E,"_cond") := helper1[1,]]
    double_data[,paste0("sigma",E,"_cond") := helper1[2,]]
    doubles_dt <- merge(doubles_dt, double_data[,.SD,,.SDcols = c("Pos1", "Pos2", "Mut1", "Mut2", paste0("fitness",E,"_cond"), paste0("sigma",E,"_cond"))], by = c("Pos1", "Pos2", "Mut1", "Mut2"), all.x = T)
  }
  parallel::stopCluster(clust)

  if(report){
    #Scatterplot matrix - singles
    if(dimsum_meta[["sequenceType"]]=="coding"){
      d <- GGally::ggpairs(singles_dt[Nham_aa==1 & apply(abs(singles_dt[,.SD,,.SDcols = paste0("fitness",all_reps)])==Inf, 1, sum)==0,
        grep(names(singles_dt),pattern="fitness"),with=F],
        upper=list(continuous = "cor"))
    }else{
      d <- GGally::ggpairs(singles_dt[apply(abs(singles_dt[,.SD,,.SDcols = paste0("fitness",all_reps)])==Inf, 1, sum)==0,
        grep(names(singles_dt),pattern="fitness"),with=F],
        upper=list(continuous = "cor"))
    }
    dimsum__save_png(file.path(report_outpath, "dimsum_stage_fitness_report_4_doubles_bayesian_framework_scattermatrix_singles.png"), d, width = 10, height = 10)

    #Scatterplot matrix - doubles, uncorrected
    set.seed(1)
    d <- GGally::ggpairs(doubles_dt[apply(abs(doubles_dt[,.SD,,.SDcols = paste0("fitness",all_reps,"_uncorr")])==Inf, 1, sum)==0
      ][sample(x = .N,min(c(.N,1000))),grep(names(doubles_dt),pattern=paste0("fitness(", all_reps_str, ")_uncorr")),with=F],
      upper=list(continuous = "cor"))
    dimsum__save_png(file.path(report_outpath, "dimsum_stage_fitness_report_4_doubles_bayesian_framework_scattermatrix_doubles_uncorr.png"), d, width = 10, height = 10)

    #Scatterplot matrix - doubles, conditional
    set.seed(1)
    d <- GGally::ggpairs(doubles_dt[apply(abs(doubles_dt[,.SD,,.SDcols = paste0("fitness",all_reps,"_cond")])==Inf, 1, sum)==0
      ][sample(x = .N,min(c(.N,1000))),grep(names(doubles_dt),pattern=paste0("fitness(", all_reps_str, ")_cond")),with=F],
      upper=list(continuous = "cor"))
    dimsum__save_png(file.path(report_outpath, "dimsum_stage_fitness_report_4_doubles_bayesian_framework_scattermatrix_doubles_cond.png"), d, width = 10, height = 10)
  }

  dimsum__status_message("Done\n")

  return(doubles_dt)

}

# ---- dimsum__calculate_fitness
dimsum__calculate_fitness <- function(
  dimsum_meta,
  input_dt,
  all_reps,
  error_model_dt,
  norm_model_dt,
  verbose = T
  ){

  if(verbose){dimsum__status_message("Calculating fitness and error...\n")}

  #Number of input and output replicates
  all_reps_str <- paste0(all_reps, collapse="|")

  #Calculate fitness
  for(j in all_reps){
    wt_corr <- as.numeric(input_dt[WT==T & error_model==T, log(.SD[,2]/.SD[,1]),,
      .SDcols = c(grep(paste0("count_e", j, "_s0"), names(input_dt)), grep(paste0("count_e", j, "_s1"), names(input_dt)))])
    input_dt[, paste0("fitness", j, "_uncorr") := log(.SD[,2]/.SD[,1]) - wt_corr,,
      .SDcols = c(
        grep(paste0("count_e", j, "_s0"), names(input_dt)),
        grep(paste0("count_e", j, "_s1"), names(input_dt)))]
    #Set infinite or undefined fitness values to NA
    input_dt[is.nan(get(paste0("fitness",j, "_uncorr"))) | is.infinite(get(paste0("fitness",j, "_uncorr"))), paste0("fitness",j, "_uncorr") := NA]
  }

  #Normalize fitness if necessary
  if(dimsum_meta[["fitnessNormalise"]]){
    #Wild-type correction such that mean(wild-type) = 0
    wt_corr <- input_dt[WT == T & error_model==T, rowMeans((.SD +
        unlist(norm_model_dt[,.SD,,.SDcols = grep(paste0("shift_(", all_reps_str, ")$"), names(norm_model_dt))])) *
        unlist(norm_model_dt[,.SD,,.SDcols = grep(paste0("scale_(", all_reps_str, ")$"), names(norm_model_dt))])),,
      .SDcols = grep(paste0("fitness(", all_reps_str, ")_uncorr$"), names(input_dt))]
    #Shift and scale
    for(j in all_reps){
      input_dt[, paste0("fitness", j, "_uncorr") := (.SD +
          unlist(norm_model_dt[,.SD,,.SDcols = paste0("shift_", j)])) *
          unlist(norm_model_dt[,.SD,,.SDcols = paste0("scale_", j)]) - wt_corr,,
        .SDcols = paste0("fitness", j, "_uncorr")]
    }
  }

  #Use error model parameters (or count-based error only) to calculate replicate-specific errors per variant
  for(j in all_reps){
    if(dimsum_meta[["fitnessNormalise"]]){
      Corr <- matrix(unlist(norm_model_dt[,.SD,.SDcols = paste0("scale_", j)]), ncol = 1, nrow = nrow(input_dt))
    }else{
      Corr <- matrix(1, ncol = 1, nrow = nrow(input_dt))
    }

    #Error model was fit
    if(dimsum_meta[["fitnessErrorModel"]]){
      input_dt[, paste0("sigma", j, "_uncorr") := sqrt(Corr * rowSums(matrix(unlist(error_model_dt[parameter %in% c("input", "output") & rep == j, mean_value]), nrow = .N, ncol = 2, byrow = T)/.SD) +
        matrix(error_model_dt[parameter %in% c("reperror") & rep == j, mean_value], nrow = .N, ncol = 1, byrow = T)),,
        .SDcols = c(
          grep(paste0("count_e", j, "_s0"), names(input_dt)),
          grep(paste0("count_e", j, "_s1"), names(input_dt)))]
    }else{
      #Calculate count-based error per replicate
      wt_corr = as.numeric(input_dt[WT==T & error_model==T,1/.SD[,2] + 1/.SD[,1],,
        .SDcols = c(grep(paste0("count_e", j, "_s0"),names(input_dt)),grep(paste0("count_e", j, "_s1"),names(input_dt)))])

      input_dt[,paste0("sigma", j, "_uncorr") := sqrt(1/.SD[,2] + 1/.SD[,1] + wt_corr),,
        .SDcols = c(grep(paste0("count_e", j, "_s0"),names(input_dt)),grep(paste0("count_e", j, "_s1"),names(input_dt)))]
    }

    #Set error of NA fitness values to NA
    input_dt[is.na(get(paste0("fitness", j, "_uncorr"))),paste0("sigma", j, "_uncorr") := NA]
  }

  #Remove unnecessary columns
  output_dt <- input_dt[,.SD,,.SDcols = c(
    "merge_seq","nt_seq","aa_seq","Nham_nt","Nham_aa",
    "Nmut_codons","WT","indel","STOP","STOP_readthrough","error_model",names(input_dt)[grep(names(input_dt),pattern="^count|^fitness|^sigma")])]

  #Remove variants without fitness estimates in any replicates
  fitness_cols <- names(output_dt)[grep(paste0("fitness(", all_reps_str, ")_uncorr"), names(output_dt))]
  output_dt <- output_dt[!is.nan(rowMeans(output_dt[,fitness_cols,with=F], na.rm = T))]

  #Calculate mean input counts
  output_dt[,mean_count := rowMeans(.SD, na.rm = T),,.SDcols = paste0("count_e", all_reps, "_s0")]

  #Save fitness estimates
  all_variants <- output_dt
  save(all_variants,
    file = file.path(dimsum_meta[["fitness_path"]], paste0(dimsum_meta[["projectName"]], '_fitness_intermediate.RData')),
    version = 2)

  if(verbose){dimsum__status_message("Done\n")}

  return(output_dt)

}

# ---- dimsum__check_barcode_design
dimsum__check_barcode_design <- function(
  barcode_design
  ){

  ### Column existence checks
  #Check if mandatory columns present
  mandatory_cols <- c("pair1", "pair2", "barcode1", "barcode2", "new_pair_prefix")
  if(sum(unlist(lapply(mandatory_cols, "%in%", colnames(barcode_design)))==FALSE)!=0){
    stop(paste0("One or more mandatory columns missing from barcodeDesign file ('pair1', 'pair2', 'barcode1', 'barcode2', 'new_pair_prefix')"), call. = FALSE)
  }

  ### FASTQ file checks (pair1/pair2 columns)
  #Check file name columns are of type character
  if(typeof(barcode_design[,"pair1"])!="character" | typeof(barcode_design[,"pair2"])!="character"){
    stop(paste0("One or more invalid FASTQ file name values in barcodeDesign file (only characters allowed)"), call. = FALSE)
  }
  #Check for incorrect pairing of FASTQ files
  unique_pairs <- unique(barcode_design[,c("pair1", "pair2")])
  if(sum(duplicated(unique_pairs[,"pair1"]))!=0 | sum(duplicated(unique_pairs[,"pair2"]))!=0){
    stop(paste0("FASTQ files not correctly paired in barcodeDesign file columns: 'pair1' and 'pair2'"), call. = FALSE)
  }
  #Check for FASTQ files with same name in different directories - pair1
  if(sum(duplicated(basename(unique(file.path(barcode_design[,"pair_directory"], barcode_design[,"pair1"])))))!=0){
    stop(paste0("FASTQ files with same name in different directories not allowed in barcodeDesign file columns: 'pair1' and 'pair2'"), call. = FALSE)
  }
  #Check for FASTQ files with same name in different directories - pair2
  if(sum(duplicated(basename(unique(file.path(barcode_design[,"pair_directory"], barcode_design[,"pair2"])))))!=0){
    stop(paste0("FASTQ files with same name in different directories not allowed in barcodeDesign file columns: 'pair1' and 'pair2'"), call. = FALSE)
  }

  ### FASTQ file directory checks (pair_directory column)
  #Check pair_directory column is of type character
  if(typeof(barcode_design[,"pair_directory"])!="character"){
    stop(paste0("One or more invalid 'fastqFileDir' values in barcodeDesign file (only characters allowed)"), call. = FALSE)
  }

  ### Barcode checks
  #Check barcode column is of type character
  if(typeof(barcode_design[,"barcode1"])!="character" | typeof(barcode_design[,"barcode2"])!="character"){
    stop(paste0("One or more invalid barcode values in barcodeDesign file (only characters allowed)"), call. = FALSE)
  }
  #Check barcode sequences are valid (ACGT characters only)
  all_characters <- unique(unlist(strsplit(as.character(unlist(barcode_design[,c("barcode1", "barcode2")])), "")))
  if(sum(!all_characters %in% c("A", "C", "G", "T", NA))!=0){
    stop("Invalid barcode sequences. Only valid nucleotide sequences allowed (A/C/T/G).", call. = FALSE)
  }

  ### New pair prefix name checks (new_pair_prefix column)
  #Check new_pair_prefix column is of type character
  if(typeof(barcode_design[,"new_pair_prefix"])!="character"){
    stop(paste0("One or more invalid 'new_pair_prefix' values in barcodeDesign file (only alphanumeric and underscore characters allowed)"), call. = FALSE)
  }
  #Check for duplicated prefices
  if(sum(duplicated(barcode_design[,"new_pair_prefix"]))!=0){
    stop(paste0("Duplicate 'new_pair_prefix' values not allowed in barcodeDesign file"), call. = FALSE)
  }
  #Check only alphanumeric characters (and underscore)
  if(sum(!grepl("^[A-Za-z0-9_]+$", barcode_design[,"new_pair_prefix"], perl = T))!=0){
    stop(paste0("One or more invalid 'new_pair_prefix' values in barcodeDesign file (only alphanumeric and underscore characters allowed)"), call. = FALSE)
  }

  ### Duplicate matrix row checks
  #Check for duplicated rows in the following column pairs: "pair1":"barcode" and "pair2":"barcode"
  if(sum(duplicated(barcode_design[,c("pair1", "barcode1", "barcode2")]))!=0 | sum(duplicated(barcode_design[,c("pair2", "barcode1", "barcode2")]))!=0){
    stop(paste0("One or more duplicated mappings in barcodeDesign file matrix (FASTQ file to barcode mapping should be unique)"), call. = FALSE)
  }

}

# ---- dimsum__check_barcodeidentityfile
dimsum__check_barcodeidentityfile <- function(
  dimsum_meta,
  return_data = FALSE
  ){

  ### Abort if no count file supplied
  if(is.null(dimsum_meta[["barcodeIdentityPath"]])){return(NULL)}

  #Load file
  barcode_dt <- fread(dimsum_meta[["barcodeIdentityPath"]])
  #Check if mandatory columns present
  mandatory_cols <- c("barcode", "variant")
  if(sum(unlist(lapply(mandatory_cols, "%in%", colnames(barcode_dt)))==FALSE)!=0){
    stop(paste0("One or more mandatory columns missing from barcodeIdentity file ('barcode', 'variant')"), call. = FALSE)
  }
  #Check all columns are of type character
  if(typeof(barcode_dt[,barcode])!="character" | typeof(barcode_dt[,variant])!="character"){
    stop(paste0("One or more invalid values in barcodeIdentity file (only A/C/G/T characters allowed)"), call. = FALSE)
  }

  ### Return data if required
  if(return_data){
    return(barcode_dt)
  }
}

# ---- dimsum__check_countfile
dimsum__check_countfile <- function(
  dimsum_meta,
  return_data = FALSE
  ){

  ### Abort if no count file supplied
  if(is.null(dimsum_meta[["countPath"]])){return(NULL)}

  ### Load variant data
  input_dt <- data.table::fread(dimsum_meta[["countPath"]])

  ### Nucleotide sequence checks (nt_seq column)
  #Check if mandatory columns present
  mandatory_cols <- c("nt_seq")
  if(sum(unlist(lapply(mandatory_cols, "%in%", colnames(input_dt)))==FALSE)!=0){
    stop(paste0("One or more mandatory columns missing from file specified by countPath ('nt_seq')"), call. = FALSE)
  }
  #Check nucleotide sequence column is of type character
  if(sapply(input_dt, typeof)["nt_seq"]!="character"){
    stop("One or more invalid 'nt_seq' values in variant count file specified by 'countPath'. Only valid nucleotide sequences allowed (A/C/T/G).", call. = FALSE)
  }
  #Set nucleotide sequence to lower case
  input_dt[, nt_seq := tolower(nt_seq)]
  #Check nucleotide sequences are valid (ACGT characters only)
  if(sum(input_dt[,grepl("[^acgt]", nt_seq)])!=0){
    stop("One or more invalid 'nt_seq' values in variant count file specified by 'countPath'. Only valid nucleotide sequences allowed (A/C/T/G).", call. = FALSE)
  }

  ### Count column checks
  #Check if sample name columns present
  mandatory_cols <- unique(dimsum_meta[["exp_design"]][,"sample_name"])
  if(sum(unlist(lapply(mandatory_cols, "%in%", colnames(input_dt)))==FALSE)!=0){
    stop(paste0("One or more sample names in experimentDesign file missing from column names in variant count file specified by 'countPath'"), call. = FALSE)
  }
  #Check all count columns are of type integer
  typeof_cols <- sapply(input_dt[,.SD,,.SDcols = names(input_dt)[names(input_dt)!="nt_seq"]], typeof)
  if(sum(typeof_cols!="integer")!=0){
    stop(paste0("Invalid type of sample count column in variant count file specified by 'countPath'. Only positive integers allowed (zero inclusive)."), call. = FALSE)
  }
  #Check all count columns positive integer zero inclusive
  if(input_dt[,min(.SD, na.rm = T),,.SDcols = names(input_dt)[names(input_dt)!="nt_seq"]]<0){
    stop(paste0("Invalid type of sample count column in variant count file specified by 'countPath'. Only positive integers allowed (zero inclusive)."), call. = FALSE)
  }

  ### Duplicated variants check
  if(input_dt[,sum(duplicated(nt_seq))]!=0){
    stop(paste0("Duplicated 'nt_seq' values not allowed in variant count file specified by 'countPath'."), call. = FALSE)
  }

  ### Return data if required
  if(return_data){
    #Sample names (ignore 'technical_replicate' column)
    sample_names <- as.list(paste0(
      dimsum_meta[["exp_design"]][,"sample_name"], '_e',
      dimsum_meta[["exp_design"]][,"experiment"], '_s',
      dimsum_meta[["exp_design"]][,"selection_id"], '_b',
      dimsum_meta[["exp_design"]][,"biological_replicate"], '_tNA_count', sep = ""))
    names(sample_names) <- dimsum_meta[["exp_design"]][,"sample_name"]

    #Reformat count column names
    names(input_dt)[names(input_dt)!="nt_seq"] <- unlist(sample_names[names(input_dt)[names(input_dt)!="nt_seq"]])

    return(input_dt)
  }
}

# ---- dimsum__check_experiment_design
dimsum__check_experiment_design <- function(
  dimsum_meta,
  exp_design
  ){

  ### Column existence checks
  #Check if mandatory columns present
  mandatory_cols <- c("sample_name", "experiment_replicate", "selection_id", "selection_replicate", "technical_replicate", "pair1", "pair2")
  if(sum(unlist(lapply(mandatory_cols, "%in%", colnames(exp_design)))==FALSE)!=0){
    stop(paste0("One or more mandatory columns missing from experimentDesign file ('sample_name', 'experiment_replicate', 'selection_id', 'selection_replicate', 'technical_replicate', 'pair1', 'pair2')"), call. = FALSE)
  }

  ### FASTQ file checks (pair1/pair2 columns)
  #Check file name columns are of type character
  if(typeof(exp_design[,"pair1"])!="character" | typeof(exp_design[,"pair2"])!="character"){
    stop(paste0("One or more invalid FASTQ file name values in experimentDesign file (only characters allowed)"), call. = FALSE)
  }
  #Check for duplicated FASTQ files
  if(sum(duplicated(exp_design[,c("pair1", "pair2")]))!=0 & !dimsum_meta[["experimentDesignPairDuplicates"]]){
    stop(paste0("Duplicate FASTQ files not allowed in experimentDesign file columns: 'pair1' and 'pair2'"), call. = FALSE)
  }
  #Check for FASTQ files with same name in different directories - pair1
  if(sum(duplicated(basename(unique(file.path(exp_design[,"pair_directory"], exp_design[,"pair1"])))))!=0){
    stop(paste0("FASTQ files with same name in different directories not allowed in experimentDesign file columns: 'pair1' and 'pair2'"), call. = FALSE)
  }
  #Check for FASTQ files with same name in different directories - pair2
  if(sum(duplicated(basename(unique(file.path(exp_design[,"pair_directory"], exp_design[,"pair2"])))))!=0){
    stop(paste0("FASTQ files with same name in different directories not allowed in experimentDesign file columns: 'pair1' and 'pair2'"), call. = FALSE)
  }

  ### FASTQ file directory checks (pair_directory column)
  #Check pair_directory column is of type character
  if(typeof(exp_design[,"pair_directory"])!="character"){
    stop(paste0("One or more invalid 'fastqFileDir' values in experimentDesign file (only characters allowed)"), call. = FALSE)
  }

  ### Sample name checks (sample_name column)
  #Check sample name column is of type character
  if(typeof(exp_design[,"sample_name"])!="character"){
    stop(paste0("One or more invalid 'sample_name' values in experimentDesign file (only characters allowed)"), call. = FALSE)
  }
  #Check for duplicated sample names (after collapsing on technical replicate id)
  sample_names_collapsed <- unique(exp_design[,c("sample_name", "experiment_replicate", "selection_id", "selection_replicate")])[,"sample_name"]
  if(sum(duplicated(sample_names_collapsed))!=0){
    stop(paste0("Duplicate 'sample_name' values not allowed in experimentDesign file (if not technical replicates)"), call. = FALSE)
  }
  #Check only alphanumeric characters
  if(sum(!grepl("^[A-Za-z0-9]+$", exp_design[,"sample_name"], perl = T))!=0){
    stop(paste0("One or more invalid 'sample_name' values in experimentDesign file (only alphanumeric characters allowed)"), call. = FALSE)
  }
  #Check first character alphabetic
  if(sum(!grepl("^[A-Za-z]+$", substr(exp_design[,"sample_name"], 1, 1), perl = T))!=0){
    stop(paste0("One or more invalid 'sample_name' values in experimentDesign file (first character must be alphabetic)"), call. = FALSE)
  }
  #Check that all input sample names for the same experiment_replicate are identical
  experiment_input <- unique(exp_design[exp_design[,"selection_id"]==0,c("sample_name", "experiment_replicate")])
  if(sum(duplicated(experiment_input[,"experiment_replicate"]))!=0){
    stop(paste0("One or more invalid input 'sample_name' values in experimentDesign file (must be identical for the same experiment replicate)"), call. = FALSE)
  }
  #Check that all output sample names for the same experiment_replicate and selection_replicate are identical
  experiment_output <- unique(exp_design[exp_design[,"selection_id"]!=0,c("sample_name", "experiment_replicate", "selection_replicate")])
  if(sum(duplicated(experiment_output[,c("experiment_replicate", "selection_replicate")]))!=0){
    stop(paste0("One or more invalid output 'sample_name' values in experimentDesign file (must be identical for the same experiment replicate and selection replicate)"), call. = FALSE)
  }

  ### Experiment id checks (experiment_replicate column)
  #Check experiment_replicate strictly positive integer
  if(typeof(exp_design[,"experiment_replicate"])!="integer"){
    stop(paste0("One or more invalid 'experiment_replicate' values in experimentDesign file (only positive integers allowed)"), call. = FALSE)
  }
  if(sum(exp_design[,"experiment_replicate"]<=0)!=0){
    stop(paste0("One or more invalid 'experiment_replicate' values in experimentDesign file (only positive integers allowed)"), call. = FALSE)
  }
  # #Check that at least two experimental replicates exist (requirement to run error model)
  # if(length(unique(exp_design[,"experiment_replicate"]))<2){
  #   stop(paste0("One or more invalid 'experiment_replicate' values in experimentDesign file (at least two replicates required)"), call. = FALSE)
  # }
  # #Check that less than 10 experimental replicates exist (more than 9 not supported)
  # if(length(unique(exp_design[,"experiment_replicate"]))>9){
  #   stop(paste0("One or more invalid 'experiment_replicate' values in experimentDesign file (more than 9 not supported)"), call. = FALSE)
  # }

  ### Selection id checks (selection_id column)
  #Check selection_id strictly positive integer (zero inclusive)
  if(typeof(exp_design[,"selection_id"])!="integer"){
    stop(paste0("One or more invalid 'selection_id' values in experimentDesign file (only 0 or positive integers allowed)"), call. = FALSE)
  }
  if(sum(exp_design[,"selection_id"]<0)!=0){
    stop(paste0("One or more invalid 'selection_id' values in experimentDesign file (only 0 or positive integers allowed)"), call. = FALSE)
  }

  ### Minimum required matrix data checks
  #Check if at least one input and output sample present for each experiment_replicate
  experiment_input <- unique(exp_design[exp_design[,"selection_id"]==0,"experiment_replicate"])
  experiment_output <- unique(exp_design[exp_design[,"selection_id"]!=0,"experiment_replicate"])
  if(length(experiment_input)!=length(experiment_output) | sum(!experiment_input %in% experiment_output)!=0){
    stop(paste0("All experiment replicates should have at least one sample before and after selection in experimentDesign file"), call. = FALSE)
  }

  ### Biological replicate id checks (selection_replicate column)
  #Check that selection_replicate column blank (NA) for input samples
  if(sum(!is.na(exp_design[exp_design[,"selection_id"]==0,"selection_replicate"]))!=0){
    stop(paste0("One or more invalid 'selection_replicate' values in experimentDesign file (leave blank for input samples)"), call. = FALSE)
  }
  #Check selection_replicate strictly positive integer
  if(typeof(exp_design[,"selection_replicate"])!="integer"){
    stop(paste0("One or more invalid 'selection_replicate' values in experimentDesign file (only positive integers allowed)"), call. = FALSE)
  }
  if(sum(exp_design[,"selection_replicate"]<=0, na.rm = T)!=0){
    stop(paste0("One or more invalid 'selection_replicate' values in experimentDesign file (only positive integers allowed)"), call. = FALSE)
  }

  ### Technical replicate id checks (technical_replicate column)
  #Check technical_replicate strictly positive integer (or all missing/NA)
  if(typeof(exp_design[,"technical_replicate"])!="integer" & sum(!is.na(exp_design[,"technical_replicate"]))!=0){
    stop(paste0("One or more invalid 'technical_replicate' values in experimentDesign file (only positive integers allowed)"), call. = FALSE)
  }
  if(sum(exp_design[,"technical_replicate"]<=0, na.rm = T)!=0){
    stop(paste0("One or more invalid 'technical_replicate' values in experimentDesign file (only positive integers allowed)"), call. = FALSE)
  }

  ### Duplicate matrix row checks
  #Check for duplicated rows in the following columns: "experiment_replicate", "selection_id", "selection_replicate", "technical_replicate"
  if(sum(duplicated(exp_design[,c("experiment_replicate", "selection_id", "selection_replicate", "technical_replicate")]))!=0){
    stop(paste0("One or more duplicated rows in experimentDesign file matrix (sample rows should be unique)"), call. = FALSE)
  }

  # ### 5'/3' adapter/cut existence checks
  # #Check that at least one cutadapt argument specified for all samples
  # if(sum(is.na(exp_design[,"cutadapt5First"]))!=0 & sum(is.na(exp_design[,"cutadapt5Second"]))!=0 & sum(is.na(exp_design[,"cutadapt3First"]))!=0 & sum(is.na(exp_design[,"cutadapt3Second"]))!=0){
  #   if(sum(is.na(exp_design[,"cutadaptCut5First"]))!=0 & sum(is.na(exp_design[,"cutadaptCut5Second"]))!=0 & sum(is.na(exp_design[,"cutadaptCut3First"]))!=0 & sum(is.na(exp_design[,"cutadaptCut3Second"]))!=0){
  #     stop("At least one of cutadapt5First, cutadapt5Second, cutadapt3First, cutadapt3Second, cutadaptCut5First, cutadaptCut5Second, cutadaptCut3First, cutadaptCut3Second must be specified for all samples. Please check that the corresponding experiment design file columns are correct.", call. = FALSE)
  #   }
  # }

  ### Constant region checks
  #Check constant region columns are of type character (or logical i.e. all empty/NA)
  if(sum(!unlist(lapply(exp_design[,c("cutadapt5First", "cutadapt5Second", "cutadapt3First", "cutadapt3Second")], typeof)) %in% c("character", "logical"))!=0){
    stop(paste0("One or more invalid constant region sequences. Only valid nucleotide sequences allowed (A/C/T/G)."), call. = FALSE)
  }
  #Check constant region sequences are valid (ACGT characters only)
  all_constant_regions <- as.character(unlist(exp_design[,c("cutadapt5First", "cutadapt5Second", "cutadapt3First", "cutadapt3Second")]))
  all_constant_regions <- gsub(";required\\.\\.\\.|;optional\\.\\.\\.", "\\.\\.\\.", all_constant_regions)
  all_constant_regions <- gsub(";required$|;optional$", "", all_constant_regions)
  all_characters <- unique(unlist(strsplit(gsub("[ACGT]\\.\\.\\.[ACGT]", "", all_constant_regions), "")))
  if(sum(!all_characters %in% c("A", "C", "G", "T", NA))!=0){
    stop("Invalid constant region sequences. Only valid nucleotide sequences allowed (A/C/T/G).", call. = FALSE)
  }

  ### Misc cutadapt argument checks
  #Check cutadaptCut... columns are of type integer (or logical i.e. all empty/NA)
  if(sum(!unlist(lapply(exp_design[,c("cutadaptCut5First", "cutadaptCut5Second", "cutadaptCut3First", "cutadaptCut3Second")], typeof)) %in% c("integer", "logical"))!=0){
    stop(paste0("One or more invalid 'cutadaptCut...' arguments. Only positive integers allowed (zero exclusive)."), call. = FALSE)
  }
  #Check strictly positive integer cutadaptCut... arguments (if not NA)
  if(sum(unlist(exp_design[,c("cutadaptCut5First", "cutadaptCut5Second", "cutadaptCut3First", "cutadaptCut3Second")])<=0, na.rm = T)!=0){
    stop("Invalid 'cutadaptCut...' arguments. Only positive integers allowed (zero exclusive).", call. = FALSE)
  }
  #Check cutadaptMinLength strictly positive integer
  if(typeof(exp_design[,"cutadaptMinLength"])!="integer"){
    stop("Invalid 'cutadaptMinLength' argument. Only positive integers allowed (zero exclusive).", call. = FALSE)
  }
  #Check cutadaptMinLength argument
  if(sum(exp_design[,c("cutadaptMinLength")]<1)!=0){
    stop("Invalid 'cutadaptMinLength' argument. Only positive integers allowed (zero exclusive).", call. = FALSE)
  }
  #Check cutadaptErrorRate double
  if(typeof(exp_design[,"cutadaptErrorRate"])!="double"){
    stop("Invalid 'cutadaptErrorRate' argument. Only positive doubles less than 1 allowed (zero inclusive).", call. = FALSE)
  }
  #Check cutadaptErrorRate positive double less than 1
  if(sum(exp_design[,c("cutadaptErrorRate")]<0)!=0 | sum(exp_design[,c("cutadaptErrorRate")]>=1)!=0){
    stop("Invalid 'cutadaptErrorRate' argument. Only positive doubles less than 1 allowed (zero inclusive).", call. = FALSE)
  }
  #Check cutadaptOverlap integer
  if(typeof(exp_design[,"cutadaptOverlap"])!="integer"){
    stop("Invalid 'cutadaptOverlap' argument. Only positive integers allowed (zero inclusive).", call. = FALSE)
  }
  #Check cutadaptOverlap positive integer zero inclusive
  if(sum(exp_design[,c("cutadaptOverlap")]<0)!=0){
    stop("Invalid 'cutadaptOverlap' argument. Only positive integers allowed (zero inclusive).", call. = FALSE)
  }

  ### Generations checks
  #Check generations column is of type double (or logical i.e. all empty/NA)
  if(!typeof(exp_design[,"generations"]) %in% c("double", "logical")){
    stop(paste0("One or more invalid generations values. Only positive doubles allowed (zero exclusive)."), call. = FALSE)
  }
  #Check generations column is strictly positive
  if(sum(exp_design[,"generations"]<=0, na.rm=T)!=0){
    stop(paste0("One or more invalid generations values. Only positive doubles allowed (zero exclusive)."), call. = FALSE)
  }
  #Check that all output samples have a generations value (if one or more specified)
  if(typeof(exp_design[,"generations"])=="double" & sum(is.na(exp_design[exp_design[,"selection_id"]==1,"generations"]))!=0){
    stop(paste0("One or more missing generations values. Generations values must be specified for all output samples (or none)."), call. = FALSE)
  }
  #Check that generations values are identical for technical output replicates (if one or more specified)
  temp_n_output_replicates <- length(unique(exp_design[exp_design[,"selection_id"]!=0,c("sample_name")]))
  if(typeof(exp_design[,"generations"])=="double" & dim(unique(exp_design[exp_design[,"selection_id"]!=0,c("sample_name", "generations")]))[1]!=temp_n_output_replicates){
    stop(paste0("Generations values not identical for technical output replicates. Generations values must be specified for all output samples (or none)."), call. = FALSE)
  }

  ### Selection time checks
  #Check selection_time column is of type double (or logical i.e. all empty/NA)
  if(!typeof(exp_design[,"selection_time"]) %in% c("double", "logical")){
    stop(paste0("One or more invalid selection_time values. Only positive doubles allowed (zero exclusive)."), call. = FALSE)
  }
  #Check selection_time column is strictly positive
  if(sum(exp_design[,"selection_time"]<=0, na.rm=T)!=0){
    stop(paste0("One or more invalid selection_time values. Only positive doubles allowed (zero exclusive)."), call. = FALSE)
  }
  #Check that all output samples have a selection_time value (if one or more specified)
  if(typeof(exp_design[,"selection_time"])=="double" & sum(is.na(exp_design[exp_design[,"selection_id"]==1,"selection_time"]))!=0){
    stop(paste0("One or more missing selection_time values. Selection time values must be specified for all output samples (or none)."), call. = FALSE)
  }
  #Check that selection_time values are identical for technical output replicates (if one or more specified)
  temp_n_output_replicates <- length(unique(exp_design[exp_design[,"selection_id"]!=0,c("sample_name")]))
  if(typeof(exp_design[,"selection_time"])=="double" & dim(unique(exp_design[exp_design[,"selection_id"]!=0,c("sample_name", "selection_time")]))[1]!=temp_n_output_replicates){
    stop(paste0("Selection time values not identical for technical output replicates. Selection time values must be specified for all output samples (or none)."), call. = FALSE)
  }

  ### Cell density checks
  #Check cell_density column is of type double (or logical i.e. all empty/NA)
  if(!typeof(exp_design[,"cell_density"]) %in% c("double", "logical")){
    stop(paste0("One or more invalid cell_density values. Only positive doubles allowed (zero exclusive)."), call. = FALSE)
  }
  #Check cell_density column is strictly positive
  if(sum(exp_design[,"cell_density"]<=0, na.rm=T)!=0){
    stop(paste0("One or more invalid cell_density values. Only positive doubles allowed (zero exclusive)."), call. = FALSE)
  }
  #Check that all samples have a cell_density value (if one or more specified)
  if(typeof(exp_design[,"cell_density"])=="double" & sum(is.na(exp_design[,"cell_density"]))!=0){
    stop(paste0("One or more missing cell_density values. Cell density values must be specified for all samples (or none)."), call. = FALSE)
  }
  #Check that cell_density values are identical for all input technical replicates (if one or more specified)
  temp_n_input_replicates <- length(unique(exp_design[exp_design[,"selection_id"]==0,c("sample_name")]))
  if(typeof(exp_design[,"cell_density"])=="double" & dim(unique(exp_design[exp_design[,"selection_id"]==0,c("sample_name", "cell_density")]))[1]!=temp_n_input_replicates){
    stop(paste0("Cell density values not identical for technical input replicates. Cell density values must be specified for all samples (or none)."), call. = FALSE)
  }
  #Check that cell_density values are identical for all output technical replicates (if one or more specified)
  temp_n_output_replicates <- length(unique(exp_design[exp_design[,"selection_id"]!=0,c("sample_name")]))
  if(typeof(exp_design[,"cell_density"])=="double" & dim(unique(exp_design[exp_design[,"selection_id"]!=0,c("sample_name", "cell_density")]))[1]!=temp_n_output_replicates){
    stop(paste0("Cell density values not identical for technical output replicates. Cell density values must be specified for all samples (or none)."), call. = FALSE)
  }

}

# ---- dimsum__check_experiment_design_countfile
dimsum__check_experiment_design_countfile <- function(
  dimsum_meta,
  exp_design
  ){

  ### Column existence checks
  #Check if mandatory columns present
  mandatory_cols <- c("sample_name", "experiment_replicate", "selection_id", "selection_replicate", "technical_replicate", "pair1", "pair2")
  if(sum(unlist(lapply(mandatory_cols, "%in%", colnames(exp_design)))==FALSE)!=0){
    stop(paste0("One or more mandatory columns missing from experimentDesign file ('sample_name', 'experiment_replicate', 'selection_id', 'selection_replicate', 'technical_replicate', 'pair1', 'pair2')"), call. = FALSE)
  }

  # ### FASTQ file checks (pair1/pair2 columns)
  # #Check file name columns are of type character
  # if(typeof(exp_design[,"pair1"])!="character" | typeof(exp_design[,"pair2"])!="character"){
  #   stop(paste0("One or more invalid FASTQ file name values in experimentDesign file (only characters allowed)"), call. = FALSE)
  # }
  # #Check for duplicated FASTQ files
  # if(sum(duplicated(exp_design[,c("pair1", "pair2")]))!=0 & !dimsum_meta[["experimentDesignPairDuplicates"]]){
  #   stop(paste0("Duplicate FASTQ files not allowed in experimentDesign file columns: 'pair1' and 'pair2'"), call. = FALSE)
  # }

  ### Sample name checks (sample_name column)
  #Check sample name column is of type character
  if(typeof(exp_design[,"sample_name"])!="character"){
    stop(paste0("One or more invalid 'sample_name' values in experimentDesign file (only characters allowed)"), call. = FALSE)
  }
  #Check for duplicated sample names (after collapsing on technical replicate id)
  sample_names_collapsed <- unique(exp_design[,c("sample_name", "experiment_replicate", "selection_id", "selection_replicate")])[,"sample_name"]
  if(sum(duplicated(sample_names_collapsed))!=0){
    stop(paste0("Duplicate 'sample_name' values not allowed in experimentDesign file (if not technical replicates)"), call. = FALSE)
  }
  #Check only alphanumeric characters
  if(sum(!grepl("^[A-Za-z0-9]+$", exp_design[,"sample_name"], perl = T))!=0){
    stop(paste0("One or more invalid 'sample_name' values in experimentDesign file (only alphanumeric characters allowed)"), call. = FALSE)
  }
  #Check first character alphabetic
  if(sum(!grepl("^[A-Za-z]+$", substr(exp_design[,"sample_name"], 1, 1), perl = T))!=0){
    stop(paste0("One or more invalid 'sample_name' values in experimentDesign file (first character must be alphabetic)"), call. = FALSE)
  }
  #Check that all input sample names for the same experiment_replicate are identical
  experiment_input <- unique(exp_design[exp_design[,"selection_id"]==0,c("sample_name", "experiment_replicate")])
  if(sum(duplicated(experiment_input[,"experiment_replicate"]))!=0){
    stop(paste0("One or more invalid input 'sample_name' values in experimentDesign file (must be identical for the same experiment replicate)"), call. = FALSE)
  }
  #Check that all output sample names for the same experiment_replicate and selection_replicate are identical
  experiment_output <- unique(exp_design[exp_design[,"selection_id"]!=0,c("sample_name", "experiment_replicate", "selection_replicate")])
  if(sum(duplicated(experiment_output[,c("experiment_replicate", "selection_replicate")]))!=0){
    stop(paste0("One or more invalid output 'sample_name' values in experimentDesign file (must be identical for the same experiment replicate and selection replicate)"), call. = FALSE)
  }

  ### Experiment id checks (experiment_replicate column)
  #Check experiment_replicate strictly positive integer
  if(typeof(exp_design[,"experiment_replicate"])!="integer"){
    stop(paste0("One or more invalid 'experiment_replicate' values in experimentDesign file (only positive integers allowed)"), call. = FALSE)
  }
  if(sum(exp_design[,"experiment_replicate"]<=0)!=0){
    stop(paste0("One or more invalid 'experiment_replicate' values in experimentDesign file (only positive integers allowed)"), call. = FALSE)
  }
  # #Check that at least two experimental replicates exist (requirement to run error model)
  # if(length(unique(exp_design[,"experiment_replicate"]))<2){
  #   stop(paste0("One or more invalid 'experiment_replicate' values in experimentDesign file (at least two replicates required)"), call. = FALSE)
  # }
  #Check that less than 10 experimental replicates exist (more than 9 not supported)
  if(length(unique(exp_design[,"experiment_replicate"]))>9){
    stop(paste0("One or more invalid 'experiment_replicate' values in experimentDesign file (more than 9 not supported)"), call. = FALSE)
  }

  ### Selection id checks (selection_id column)
  #Check selection_id strictly positive integer (zero inclusive)
  if(typeof(exp_design[,"selection_id"])!="integer"){
    stop(paste0("One or more invalid 'selection_id' values in experimentDesign file (only 0 or positive integers allowed)"), call. = FALSE)
  }
  if(sum(exp_design[,"selection_id"]<0)!=0){
    stop(paste0("One or more invalid 'selection_id' values in experimentDesign file (only 0 or positive integers allowed)"), call. = FALSE)
  }

  ### Minimum required matrix data checks
  #Check if at least one input and output sample present for each experiment_replicate
  experiment_input <- unique(exp_design[exp_design[,"selection_id"]==0,"experiment_replicate"])
  experiment_output <- unique(exp_design[exp_design[,"selection_id"]!=0,"experiment_replicate"])
  if(length(experiment_input)!=length(experiment_output) | sum(!experiment_input %in% experiment_output)!=0){
    stop(paste0("All experiment replicates should have at least one sample before and after selection in experimentDesign file"), call. = FALSE)
  }

  ### Biological replicate id checks (selection_replicate column)
  #Check that selection_replicate column blank (NA) for input samples
  if(sum(!is.na(exp_design[exp_design[,"selection_id"]==0,"selection_replicate"]))!=0){
    stop(paste0("One or more invalid 'selection_replicate' values in experimentDesign file (leave blank for input samples)"), call. = FALSE)
  }
  #Check selection_replicate strictly positive integer
  if(typeof(exp_design[,"selection_replicate"])!="integer"){
    stop(paste0("One or more invalid 'selection_replicate' values in experimentDesign file (only positive integers allowed)"), call. = FALSE)
  }
  if(sum(exp_design[,"selection_replicate"]<=0, na.rm = T)!=0){
    stop(paste0("One or more invalid 'selection_replicate' values in experimentDesign file (only positive integers allowed)"), call. = FALSE)
  }

  # ### Technical replicate id checks (technical_replicate column)
  # #Check technical_replicate strictly positive integer (or all missing/NA)
  # if(typeof(exp_design[,"technical_replicate"])!="integer" & sum(!is.na(exp_design[,"technical_replicate"]))!=0){
  #   stop(paste0("One or more invalid 'technical_replicate' values in experimentDesign file (only positive integers allowed)"), call. = FALSE)
  # }
  # if(sum(exp_design[,"technical_replicate"]<=0, na.rm = T)!=0){
  #   stop(paste0("One or more invalid 'technical_replicate' values in experimentDesign file (only positive integers allowed)"), call. = FALSE)
  # }

  ### Duplicate matrix row checks
  #Check for duplicated rows in the following columns: "experiment_replicate", "selection_id", "selection_replicate"
  if(sum(duplicated(exp_design[,c("experiment_replicate", "selection_id", "selection_replicate")]))!=0){
    stop(paste0("One or more duplicated rows in experimentDesign file matrix (sample rows should be unique)"), call. = FALSE)
  }

  # ### 5'/3' adapter/cut existence checks
  # #Check that at least one cutadapt argument specified for all samples
  # if(sum(is.na(exp_design[,"cutadapt5First"]))!=0 & sum(is.na(exp_design[,"cutadapt5Second"]))!=0 & sum(is.na(exp_design[,"cutadapt3First"]))!=0 & sum(is.na(exp_design[,"cutadapt3Second"]))!=0){
  #   if(sum(is.na(exp_design[,"cutadaptCut5First"]))!=0 & sum(is.na(exp_design[,"cutadaptCut5Second"]))!=0 & sum(is.na(exp_design[,"cutadaptCut3First"]))!=0 & sum(is.na(exp_design[,"cutadaptCut3Second"]))!=0){
  #     stop("At least one of cutadapt5First, cutadapt5Second, cutadapt3First, cutadapt3Second, cutadaptCut5First, cutadaptCut5Second, cutadaptCut3First, cutadaptCut3Second must be specified for all samples. Please check that the corresponding experiment design file columns are correct.", call. = FALSE)
  #   }
  # }

  # ### Constant region checks
  # #Check constant region columns are of type character (or logical i.e. all empty/NA)
  # if(sum(!unlist(lapply(exp_design[,c("cutadapt5First", "cutadapt5Second", "cutadapt3First", "cutadapt3Second")], typeof)) %in% c("character", "logical"))!=0){
  #   stop(paste0("One or more invalid constant region sequences. Only valid nucleotide sequences allowed (A/C/T/G)."), call. = FALSE)
  # }
  # #Check constant region sequences are valid (ACGT characters only)
  # all_characters <- unique(unlist(strsplit(gsub("[ACGT]\\.\\.\\.[ACGT]", "", as.character(unlist(exp_design[,c("cutadapt5First", "cutadapt5Second", "cutadapt3First", "cutadapt3Second")]))), "")))
  # if(sum(!all_characters %in% c("A", "C", "G", "T", NA))!=0){
  #   stop("Invalid constant region sequences. Only valid nucleotide sequences allowed (A/C/T/G).", call. = FALSE)
  # }

  # ### Misc cutadapt argument checks
  # #Check cutadaptCut... columns are of type integer (or logical i.e. all empty/NA)
  # if(sum(!unlist(lapply(exp_design[,c("cutadaptCut5First", "cutadaptCut5Second", "cutadaptCut3First", "cutadaptCut3Second")], typeof)) %in% c("integer", "logical"))!=0){
  #   stop(paste0("One or more invalid 'cutadaptCut...' arguments. Only positive integers allowed (zero exclusive)."), call. = FALSE)
  # }
  # #Check strictly positive integer cutadaptCut... arguments (if not NA)
  # if(sum(unlist(exp_design[,c("cutadaptCut5First", "cutadaptCut5Second", "cutadaptCut3First", "cutadaptCut3Second")])<=0, na.rm = T)!=0){
  #   stop("Invalid 'cutadaptCut...' arguments. Only positive integers allowed (zero exclusive).", call. = FALSE)
  # }
  # #Check cutadaptMinLength strictly positive integer
  # if(typeof(exp_design[,"cutadaptMinLength"])!="integer"){
  #   stop("Invalid 'cutadaptMinLength' argument. Only positive integers allowed (zero exclusive).", call. = FALSE)
  # }
  # #Check cutadaptMinLength argument
  # if(sum(exp_design[,c("cutadaptMinLength")]<1)!=0){
  #   stop("Invalid 'cutadaptMinLength' argument. Only positive integers allowed (zero exclusive).", call. = FALSE)
  # }
  # #Check cutadaptErrorRate double
  # if(typeof(exp_design[,"cutadaptErrorRate"])!="double"){
  #   stop("Invalid 'cutadaptErrorRate' argument. Only positive doubles less than 1 allowed (zero inclusive).", call. = FALSE)
  # }
  # #Check cutadaptErrorRate positive double less than 1
  # if(sum(exp_design[,c("cutadaptErrorRate")]<0)!=0 | sum(exp_design[,c("cutadaptErrorRate")]>=1)!=0){
  #   stop("Invalid 'cutadaptErrorRate' argument. Only positive doubles less than 1 allowed (zero inclusive).", call. = FALSE)
  # }
  # #Check cutadaptOverlap integer
  # if(typeof(exp_design[,"cutadaptOverlap"])!="integer"){
  #   stop("Invalid 'cutadaptOverlap' argument. Only positive integers allowed (zero inclusive).", call. = FALSE)
  # }
  # #Check cutadaptOverlap positive integer zero inclusive
  # if(sum(exp_design[,c("cutadaptOverlap")]<0)!=0){
  #   stop("Invalid 'cutadaptOverlap' argument. Only positive integers allowed (zero inclusive).", call. = FALSE)
  # }

  ### Generations checks
  #Check generations column is of type double (or logical i.e. all empty/NA)
  if(!typeof(exp_design[,"generations"]) %in% c("double", "logical")){
    stop(paste0("One or more invalid generations values. Only positive doubles allowed (zero exclusive)."), call. = FALSE)
  }
  #Check generations column is strictly positive
  if(sum(exp_design[,"generations"]<=0, na.rm=T)!=0){
    stop(paste0("One or more invalid generations values. Only positive doubles allowed (zero exclusive)."), call. = FALSE)
  }
  #Check that all output samples have a generations value (if one or more specified)
  if(typeof(exp_design[,"generations"])=="double" & sum(is.na(exp_design[exp_design[,"selection_id"]==1,"generations"]))!=0){
    stop(paste0("One or more missing generations values. Generations values must be specified for all output samples (or none)."), call. = FALSE)
  }
  #Check that generations values are identical for technical output replicates (if one or more specified)
  temp_n_output_replicates <- length(unique(exp_design[exp_design[,"selection_id"]!=0,c("sample_name")]))
  if(typeof(exp_design[,"generations"])=="double" & dim(unique(exp_design[exp_design[,"selection_id"]!=0,c("sample_name", "generations")]))[1]!=temp_n_output_replicates){
    stop(paste0("Generations values not identical for technical output replicates. Generations values must be specified for all output samples (or none)."), call. = FALSE)
  }

}

# ---- dimsum__check_files_exist
dimsum__check_files_exist <- function(
  required_files,
  stage_number = NULL,
  execute = TRUE,
  exit = TRUE
  ){
	missing_files <- required_files[!file.exists(required_files)]
  if(length(missing_files)!=0 & execute){
    if(is.null(stage_number)){
      #No stage number supplied
      dimsum__status_message(paste0("The following files do not exist.\n"))
      dimsum__status_message(paste0(missing_files, "\n"))
      if(exit){
        stop(paste0("Cannot proceed with DiMSum. Required files do not exist."), call. = FALSE)
      }
    }else{
      #Stage number supplied
      dimsum__status_message(paste0("The following files required for DiMSum STAGE ", stage_number, " do not exist.\n"))
      dimsum__status_message(paste0(missing_files, "\n"))
      dimsum__status_message(paste0("Ensure that the 'startStage' argument has been correctly specified and that all previous stages have successfully completed.\n"))
      dimsum__status_message(paste0("If you previously ran DiMSum with 'retainIntermediateFiles' set to FALSE, you will need set 'startStage' to 0 in order to regenerate the files required for this stage.\n"))
      if(exit){
        stop(paste0("Cannot proceed with DiMSum STAGE ", stage_number, ". Required files do not exist."), call. = FALSE)
      }
    }
  }
}

# ---- dimsum__check_variants
dimsum__check_variants <- function(
  dimsum_meta,
  input_dt
  ){

  #Determine replicates to retain
  all_reps <- unique(dimsum_meta[["exp_design"]][,"experiment"])
  if(dimsum_meta[["retainedReplicates"]]!="all"){
    all_reps <- unique(as.integer(strsplit(dimsum_meta[["retainedReplicates"]], ",")[[1]]))
  }

  #Number of input and output replicates
  all_reps_str <- paste0(all_reps, collapse="|")

  #Check if WT sequence and at least 2 substitutions variants exist
  if(input_dt[WT==T,.N]==0){

    if( length(grep(paste0("count_e(", all_reps_str, ")_s[01]"), names(input_dt))) != 0 ){

	    input_dt[, all_reads := rowSums(.SD > 0) == (2*length(all_reps)),,.SDcols = grep(paste0("count_e(", all_reps_str, ")_s[01]"), names(input_dt))]
	    input_dt[, mean_count := rowMeans(.SD),,.SDcols = grep(paste0("count_e(", all_reps_str, ")_s0"), names(input_dt))]
	    dimsum__status_message(paste0("WT variant has zero count in at least one input/output replicate. Did you mean to specify one of the following?\n"))
	    if(dimsum_meta[["sequenceType"]]=="coding" & dimsum_meta[["mixedSubstitutions"]]){
	      print(input_dt[all_reads == T,][order(mean_count, decreasing = T)[1:5],.(aa_seq, all_reads, mean_count)])
	    }else if(dimsum_meta[["sequenceType"]]=="coding"){
	      print(input_dt[all_reads == T,][order(mean_count, decreasing = T)[1:5],.(nt_seq = toupper(nt_seq), aa_seq, all_reads, mean_count)])
	    }else{
	      print(input_dt[all_reads == T,][order(mean_count, decreasing = T)[1:5],.(nt_seq = toupper(nt_seq), all_reads, mean_count)])
	    }
	    stop(paste0("Cannot proceed with fitness estimation: WT variant not found"), call. = FALSE)

    }else{

	    input_dt[, all_reads := rowSums(.SD > 0) == length(grep(paste0("e(", all_reps_str, ")_s[01].*_count$"), names(input_dt))),,.SDcols = grep(paste0("e(", all_reps_str, ")_s[01].*_count$"), names(input_dt))]
	    input_dt[, mean_count := rowMeans(.SD),,.SDcols = grep(paste0("e(", all_reps_str, ")_s[0].*_count$"), names(input_dt))]
	    dimsum__status_message(paste0("WT variant not found. Did you mean to specify one of the following?\n"))
	    print(input_dt[all_reads == T,][order(mean_count, decreasing = T)[1:5],.(nt_seq = toupper(nt_seq), all_reads, mean_count)])
	    stop(paste0("Cannot proceed with fitness estimation: WT variant not found"), call. = FALSE)

    }

  }else if(input_dt[is.na(WT),.N]<2){
    stop(paste0("Cannot proceed with fitness estimation: insufficient substitution variants"), call. = FALSE)
  }

}

# ---- dimsum__create_dir
dimsum__create_dir <- function(
  dimsum_dir,
  execute,
  message = NULL,
  overwrite_dir = TRUE){
  if(!execute){
    if(!is.null(message)){
      dimsum__status_message(paste("\n\n\n*******", message, "(not executed)", "*******\n\n\n"))
    }
    return()
  }
  if(!is.null(message)){
    dimsum__status_message(paste("\n\n\n*******", message, "*******\n\n\n"))
  }
  if(dir.exists(dimsum_dir) & overwrite_dir){
    unlink(dimsum_dir, recursive = TRUE)
  }
  dir.create(dimsum_dir, showWarnings = FALSE)
}

# ---- dimsum__debarcode_variants
dimsum__debarcode_variants <- function(
  dimsum_meta,
  input_dt
  ){

  #Load barcode to variant mapping
  barcode_dt <- dimsum__check_barcodeidentityfile(dimsum_meta = dimsum_meta, return_data = TRUE)
  barcode_list <- as.list(tolower(unlist(barcode_dt[,variant])))
  names(barcode_list) <- tolower(unlist(barcode_dt[,barcode]))

  #Variants with valid barcodes
  input_dt[!nt_seq %in% names(barcode_list), barcode_valid := FALSE]

  #Check if no valid barcodes
  if(input_dt[barcode_valid==T,.N]==0){
    stop(paste0("Cannot proceed with variant processing: No valid barcodes found"), call. = FALSE)
  }

  #Debarcode variants
  debarcode_dt <- input_dt[barcode_valid==T,]
  debarcode_dt[, nt_seq := unlist(barcode_list[nt_seq])]
  #Aggregate counts for identical variants
  idx <- names(debarcode_dt)[grep(names(debarcode_dt),pattern="_count$")]
  for (i in seq_along(idx)) {
    #Aggregate counts accross identical AA variants
    debarcode_dt[,paste0(idx[i],"_agg") := sum(.SD),nt_seq,.SDcols = idx[i]]
  }
  #Retain only one row per AA variant
  debarcode_dt <- debarcode_dt[!duplicated(nt_seq),.SD,nt_seq,.SDcols = c(
    names(debarcode_dt)[grep(names(debarcode_dt),pattern="_agg$")],
    "barcode_valid")]
  #Rename count columns
  names(debarcode_dt) <- gsub("_agg$", "", names(debarcode_dt))

  #Merge variant data
  output_dt <- rbind(debarcode_dt, input_dt[barcode_valid==F,])

  return(output_dt)

}

# ---- dimsum__detect_sequence_type
dimsum__detect_sequence_type <- function(
  input_sequence
  ){

  #Sequence length a multiple of 3?
  if(nchar(input_sequence) %% 3 == 0){
    aa_seq_split <- seqinr::translate(strsplit(input_sequence,split="")[[1]])
    #Translated sequence contains internal STOP?
    if(!"*" %in% aa_seq_split[-length(aa_seq_split)]){
      return("coding")
    }
  }

  #Sequence length not a multiple of 3 or contains internal STOP
  return("noncoding")

}

# ---- dimsum__diagnostics_report
dimsum__diagnostics_report <- function(
  dimsum_meta,
  report_outpath,
  input_samples_pattern = "_s0_",
  output_samples_pattern = "_s1_"
  ){
  #Create report directory (if doesn't already exist)
  report_outpath <- gsub("/$", "", report_outpath)
  suppressWarnings(dir.create(report_outpath))

  #Check if all input files exist
  dimsum__check_files_exist(
    required_files = dimsum_meta[["variant_data_merge_path"]],
    execute = TRUE,
    exit = FALSE)

  #load variant data from RData file
  load(dimsum_meta[["variant_data_merge_path"]])

  #Sample names
  input_samples <- colnames(variant_data_merge)[grep(input_samples_pattern, colnames(variant_data_merge))]
  output_samples <- colnames(variant_data_merge)[grep(output_samples_pattern, colnames(variant_data_merge))]

  #Max number of nucleotide substitutions to plot
  max_nsubs <- dimsum_meta[["maxSubstitutions"]]
  if(dimsum_meta[["sequenceType"]]=="coding"){
    max_nsubs <- max_nsubs*3
  }
  #Set maximum limit of 12 substitutions to display
  max_nsubs <- min(c(max_nsubs, 12))

  #Error rate
  error_rate <- 10^min(-dimsum_meta[["vsearchMinQual"]]/10, log10(dimsum_meta[["vsearchMaxee"]] / nchar(dimsum_meta[["wildtypeSequence"]])))

  #Plot 1: Histogram of input counts split by number of nucleotide mutations
  if(length(input_samples)!=0){
    dimsum__sample_count_distributions(
      input_dt = variant_data_merge[data.table::between(Nham_nt,0,max_nsubs),.SD,.SDcols=c(input_samples, "Nham_nt", "WT")],
      output_file_prefix = file.path(report_outpath, "dimsum__diagnostics_report_count_hist_input_nt"),
      #title = paste0("Nucleotide substitution variant count distributions (input samples)"),
      error_rate = error_rate,
      seq_length = nchar(dimsum_meta[["wildtypeSequence"]]))
  }

  #Plot 2: Histogram of input counts split by number of nucleotide mutations and number of amino acid mutations
  max_nsubs <- min(dimsum_meta[["maxSubstitutions"]], 2)
  if(dimsum_meta[["sequenceType"]]=="coding" & max_nsubs<=2){
    max_nsubs <- max_nsubs*3
    if(length(input_samples)!=0){
      dimsum__sample_count_distributions(
        input_dt = variant_data_merge[data.table::between(Nham_nt,0,max_nsubs),.SD,.SDcols=c(input_samples, "Nham_nt", "Nham_aa", "WT")],
        output_file_prefix = file.path(report_outpath, "dimsum__diagnostics_report_count_hist_input_aa"),
        #title = paste0("Nucleotide substitution variant count distributions (input samples)"),
        error_rate = error_rate,
        seq_length = nchar(dimsum_meta[["wildtypeSequence"]]), height = 12)
    }
  }

  #Threshold lists
  threshold_list_dotted <- lapply(dimsum_meta[["fitnessMinInputCountAny"]], function(x){log10(x + 1)})
  names(threshold_list_dotted) <- rep("dotted", length(threshold_list_dotted))
  threshold_list_dashed <- lapply(dimsum_meta[["fitnessMinInputCountAll"]], function(x){log10(x + 1)})
  names(threshold_list_dashed) <- rep("dashed", length(threshold_list_dashed))

  #Plot 3a: all-vs-all sample count correlations - all variants
  if(length(input_samples)!=0 | length(output_samples)!=0){
    temp_dt <- variant_data_merge[,.SD,.SDcols = c(input_samples, output_samples, "Nham_nt")]
    names(temp_dt)[grep("_count", names(temp_dt))] <- dimsum__plot_samplename(names(temp_dt)[grep("_count", names(temp_dt))])
    dimsum__ggpairs_binhex(
      input_dt = log10(temp_dt[,.SD,.SDcols = dimsum__plot_samplename(c(input_samples, output_samples))]+1),
      output_file_prefix = file.path(report_outpath, "dimsum__diagnostics_report_scatterplotmatrix_all"),
      xlab = "log10(variant count + 1)",
      ylab = "log10(variant count + 1)",
      #title = "Substitution variant inter-sample count correlations (all samples)",
      size = 0.1,
      thresholds = c(threshold_list_dotted, threshold_list_dashed))
  }

  #Plot 3b: all-vs-all sample count correlations - only singles and doubles
  if(length(input_samples)!=0 | length(output_samples)!=0){
    temp_dt <- variant_data_merge[,.SD,.SDcols = c(input_samples, output_samples, "Nham_nt")]
    names(temp_dt)[grep("_count", names(temp_dt))] <- dimsum__plot_samplename(names(temp_dt)[grep("_count", names(temp_dt))])
    dimsum__ggpairs_binhex(
      input_dt = log10(temp_dt[data.table::between(Nham_nt,1,2),.SD,.SDcols = dimsum__plot_samplename(c(input_samples, output_samples))]+1),
      output_file_prefix = file.path(report_outpath, "dimsum__diagnostics_report_scatterplotmatrix_singles_doubles"),
      xlab = "log10(variant count + 1)",
      ylab = "log10(variant count + 1)",
      #title = "Single and double substitution variant inter-sample count correlations (all samples)",
      cut = as.factor(temp_dt[data.table::between(Nham_nt,1,2),Nham_nt]),
      size = 0.1,
      thresholds = c(threshold_list_dotted, threshold_list_dashed))
  }

  #Render report
  dimsum__render_report(dimsum_meta = dimsum_meta)

}

# ---- dimsum__error_model [modified, see the header and the generator]
dimsum__error_model <- function(
  dimsum_meta,
  input_dt,
  all_reps,
  report = TRUE,
  report_outpath = NULL
  ){

  #Return NULL model if no fit required
  if(!dimsum_meta[["fitnessErrorModel"]]){
    return(list(
      "error_model" = NULL,
      "norm_model" = NULL))
  }

  dimsum__status_message("Fit error model...\n")

  #Number of input and output replicates
  all_reps_str <- paste0(all_reps, collapse="|")

  work_data <- input_dt[,.SD,.SDcols = c("Nham_nt","Nham_aa","WT",names(input_dt)[grep(names(input_dt),pattern="^count")])]

  ### Calculate fitness
  ###########################

  for(j in all_reps){
    work_data[, paste0("fitness",j) := log(.SD[,2]/.SD[,1]),,
      .SDcols = c(
        grep(paste0("count_e", j, "_s0"), names(work_data)),
        grep(paste0("count_e", j, "_s1"), names(work_data)))]
    #Set infinite or undefined fitness values to NA
    work_data[is.nan(get(paste0("fitness",j))) | is.infinite(get(paste0("fitness",j))), paste0("fitness",j) := NA]
  }

  #Flag variants that don't have reads in all input/output replicates
  work_data[, all_reads := rowSums(.SD > 0) == (2*length(all_reps)),,.SDcols = grep(paste0("count_e(", all_reps_str, ")_s[01]"), names(work_data))]

  #Check if WT data present in all input/output replicates
  if(work_data[WT == T & all_reads == T,.N]==0){
    input_dt[, all_reads := rowSums(.SD > 0) == (2*length(all_reps)),,.SDcols = grep(paste0("count_e(", all_reps_str, ")_s[01]"), names(input_dt))]
    input_dt[, mean_count := rowMeans(.SD),,.SDcols = grep(paste0("count_e(", all_reps_str, ")_s0"), names(input_dt))]
    dimsum__status_message(paste0("WT variant has zero count in at least one input/output replicate.\n"))
    all_reads_n = input_dt[all_reads == T,.N]
    if(all_reads_n==0){
      dimsum__status_message(paste0("No variants without zero counts in all input/output replicates.\n"))
    }else{
      if(dimsum_meta[["sequenceType"]]=="coding" & dimsum_meta[["mixedSubstitutions"]]){
        dimsum__status_message(paste0("Did you mean to specify one of the following?\n"))
        print(input_dt[all_reads == T,][order(mean_count, decreasing = T)[1:all_reads_n],.(aa_seq, all_reads, mean_count)])
      }else if(dimsum_meta[["sequenceType"]]=="coding"){
        dimsum__status_message(paste0("Did you mean to specify one of the following?\n"))
        print(input_dt[all_reads == T,][order(mean_count, decreasing = T)[1:all_reads_n],.(nt_seq = toupper(nt_seq), aa_seq, all_reads, mean_count)])
      }else{
        dimsum__status_message(paste0("Did you mean to specify one of the following?\n"))
        print(input_dt[all_reads == T,][order(mean_count, decreasing = T)[1:all_reads_n],.(nt_seq = toupper(nt_seq), all_reads, mean_count)])
      }
    }
    stop(paste0("Cannot proceed with error modelling: WT variant has zero count in at least one input/output replicate"), call. = FALSE)
  }

  ### Find input read threshold for full fitness range
  ###########################

  input_count_threshold <- work_data[all_reads == T,exp(-quantile(.SD,probs = 0.01,na.rm=T)),,.SDcols = grep("fitness",names(work_data))]
  #Define variants above threshold for later use
  work_data[,input_above_threshold := rowSums(.SD > input_count_threshold) == length(all_reps),,.SDcols = grep(paste0("count_e(", all_reps_str, ")_s0"),names(work_data))]

  #Correct for WT fitness
  for(j in all_reps){
    wt_corr <- as.numeric(work_data[WT==T, .SD,,.SDcols = paste0("fitness",j)])
    work_data[, paste0("fitness",j) := .SD - wt_corr,,.SDcols = paste0("fitness",j)]
  }

  #Plot input counts versus fitness with threshold
  if(report){
    X <- data.table::data.table()
    for(j in all_reps){
      if(length(X) == 0) {
        X = work_data[,.SD,,.SDcols = c(paste0("count_e", j, "_s0"),paste0("fitness",j))]
        names(X) = c("input","fitness")
        X[,rep := paste0("replicate", j)]
      }else{
        Y = work_data[,.SD,,.SDcols = c(paste0("count_e", j, "_s0"),paste0("fitness",j))]
        names(Y) = c("input","fitness")
        Y[,rep := paste0("replicate", j)]
        X = rbind(X,Y)
      }
    }
    d <- ggplot2::ggplot(X[input>0 & !is.na(fitness),],ggplot2::aes(input,fitness)) +
      ggplot2::geom_hex(size = 0.2, color = "lightgrey") +
      # ggplot2::scale_fill_viridis_c(trans="log10") +
      ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
      dimsum__scale_x_log10() +
      ggplot2::geom_vline(xintercept = input_count_threshold, lty = 2) +
      ggplot2::geom_hline(yintercept = 0, lty = 2, color = "darkgrey") +
      ggplot2::facet_wrap(rep ~ .) + ggplot2::theme_bw() +
      ggplot2::labs(x = "Input variant count (log scale)", y = "Fitness")#, title = "Replicate fitness versus input variant counts")
    dimsum_fitness__save_png_pdf(file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_inputcounts.png"), file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_inputcounts.pdf"), d, width = 6, height = 6)
  }

  #Check if sufficient data remains above input read threshold (and data present in all input/output replicates)
  #Define minimum #variants for fitting = 10 x number of fitted parameters (3 x #reps)
  min_n_variants <- 10 * 3 * length(all_reps)
  if(work_data[input_above_threshold == T & all_reads == T & is.na(WT),.N] < min_n_variants){
    stop(paste0("Cannot proceed with error modelling: insufficent number of variants satisfying full fitness range"), call. = FALSE)
  }

  ### QC plots for fitness between replicates (before normalisation)
  ###########################

  if(report){
    #Fitness distributions (all replicates)
    X <- reshape2::melt(work_data[input_above_threshold == T & all_reads == T,.SD,.SDcols=grep("fitness",names(work_data))],measure.vars=grep("fitness",names(work_data),value=T))
    X[,replicate := gsub("fitness", "", as.character(variable)),variable]
    d <- ggplot2::ggplot(X,ggplot2::aes(value,color=replicate)) +
      ggplot2::geom_density(key_glyph = dimsum__draw_key_polygon) +
      ggplot2::labs(x = "Fitness", y = "Density", color = "replicate") +#, title = "Before inter-replicate normalisation") +
      ggplot2::geom_vline(xintercept = 0, lty = 2, color = "darkgrey") +
      ggplot2::theme_bw()
    dimsum_fitness__save_png_pdf(file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_density.png"), file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_density.pdf"), d, width = 6, height = 4)

    #Fitness correlations (all replicates)
    dimsum__ggpairs_binhex(
      input_dt = work_data[all_reads == T,.SD,,.SDcols = paste0("fitness", all_reps)],
      output_file_prefix = file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_scatter"),
      xlab = "Fitness",
      ylab = "Fitness")#, title = "Before inter-replicate normalisation")
  }

  ### Calculate replicate normalisation parameters (scale and center/shift) to minimise inter-replicate differences
  ###########################

  #Fitness data
  F_data <- work_data[input_above_threshold == T & all_reads == T, as.matrix(.SD),.SDcols = grep(paste0("fitness(", all_reps_str, ")$"),names(work_data))]

  #Calculate replicate normalisation parameters using non-linear minimization
  set.seed(1603)
  p <- nlm(f = dimsum__replicate_fitness_deviation, p = rep(c(1,0), each = length(all_reps)), fitness_mat = F_data, all_reps = all_reps)[["estimate"]]

  #Normalise to first replicate (set scaling factor of first replicate to unity)
  p[1:length(all_reps)] <- p[1:length(all_reps)]/p[1]

  #Check for negative scaling factors
  if(sum(p[1:length(all_reps)]<0)!=0){
    warning("dimsum__error_model.R: Some scaling factors from replicate normalisation are less than zero. Check that input/output samples are not mistakenly switched/misspecified in experimentDesign file!", call. = FALSE, immediate. = TRUE, noBreaks. = TRUE)
  }

  #Save normalisation model
  fitness_norm_model <- data.table(t(p))
  names(fitness_norm_model) <- c(paste0("scale_", all_reps), paste0("shift_", all_reps))
  write.table(round(fitness_norm_model, digits = 4),
    file = file.path(dimsum_meta[["tmp_path"]], "normalisationmodel.txt"), row.names = F, quote = F)

  #Wild-type correction such that mean(wild-type) = 0
  wt_corr <- work_data[WT == T, rowMeans((.SD +
      unlist(fitness_norm_model[,.SD,,.SDcols = grep(paste0("shift_(", all_reps_str, ")$"), names(fitness_norm_model))])) *
      unlist(fitness_norm_model[,.SD,,.SDcols = grep(paste0("scale_(", all_reps_str, ")$"), names(fitness_norm_model))])),,
    .SDcols = grep(paste0("fitness(", all_reps_str, ")$"), names(work_data))]

  #Normalize fitness
  for (j in all_reps){
    work_data[all_reads == T, paste0("fitness", j, "_norm") := (.SD +
        unlist(fitness_norm_model[,.SD,,.SDcols = paste0("shift_", j)])) *
        unlist(fitness_norm_model[,.SD,,.SDcols = paste0("scale_", j)]) - wt_corr,,
      .SDcols = paste0("fitness", j)]
  }

  if(report){
    #Fitness distributions (all replicates)
    X <- reshape2::melt(work_data[input_above_threshold == T & all_reads==T,.SD,.SDcols=grep("_norm",names(work_data))],measure.vars=grep("_norm",names(work_data),value=T))
    X[,replicate := gsub("_norm", "", gsub("fitness", "", as.character(variable))),variable]
    d <- ggplot2::ggplot(X,ggplot2::aes(value,color=replicate)) +
      ggplot2::geom_density(key_glyph = dimsum__draw_key_polygon) +
      ggplot2::labs(x = "Fitness", y = "Density", color = "replicate") +#, title = "After inter-replicate normalisation") +
      ggplot2::geom_vline(xintercept = 0, lty = 2, color = "darkgrey") +
      ggplot2::theme_bw()
    dimsum_fitness__save_png_pdf(file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_density_norm.png"), file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_density_norm.pdf"), d, width = 6, height = 4)

    #Fitness correlations (all replicates)
    dimsum__ggpairs_binhex(
      input_dt = work_data[all_reads == T,.SD,,.SDcols = paste0("fitness", all_reps, "_norm")],
      output_file_prefix = file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_scatter_norm"),
      xlab = "Fitness",
      ylab = "Fitness")#, title = "After inter-replicate normalisation")

    #Fitness replicate deviations (all replicates)
    X <- work_data[all_reads == T & input_above_threshold == T & is.na(WT), cbind(.SD - rowMeans(.SD), M = rowMeans(.SD)),
      .SDcols = grep(paste0("fitness(", all_reps_str, ")$"), names(work_data))]
    Xnorm <- work_data[all_reads == T & input_above_threshold == T & is.na(WT), cbind(.SD - rowMeans(.SD), M = rowMeans(.SD)),
      .SDcols = grep(paste0("fitness(", all_reps_str, ")_norm$"), names(work_data))]
    Y <- rbind(reshape2::melt(X, id.vars = "M"), reshape2::melt(Xnorm, id.vars = "M"))
    Y[, replicate := strsplit(as.character(variable), "_")[[1]][1], variable]
    Y[, normalised := grepl('norm', variable), variable]
    d <- ggplot2::ggplot(Y,ggplot2::aes(M, value, color = normalised)) +
      ggplot2::geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs")) +
      ggplot2::geom_hline(yintercept = 0, lty = 2, color = "darkgrey") +
      ggplot2::labs(x = "Mean fitness", y = "Deviation from mean fitness") +#, title = "Fitness replicate deviations (all replicates)") +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(. ~ replicate)
    dimsum_fitness__save_png_pdf(file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_replicate_deviation_scatter.png"), file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_replicate_deviation_scatter.pdf"), d, width = 6, height = 6)
  }

  #Replace fitness with normalised fitness
  if(dimsum_meta[["fitnessNormalise"]]){
    for(j in all_reps){
      work_data[all_reads ==T,paste0("fitness", j) := .SD,
                ,.SDcols = paste0("fitness", j, "_norm")]
      # work_data[all_reads == T, paste0("fitness", j, "_norm") := NULL]
    }
  }

  ### Error model
  ###########################

  #First calculate count-based error per replicate
  for (j in all_reps) {
    wt_corr = as.numeric(work_data[WT==T,1/.SD[,2] + 1/.SD[,1],,
      .SDcols = c(grep(paste0("count_e", j, "_s0"),names(work_data)),grep(paste0("count_e", j, "_s1"),names(work_data)))])

    if(dimsum_meta[["fitnessNormalise"]]){
      work_data[,paste0("cbe",j) := sqrt(abs(unlist(fitness_norm_model[,.SD,.SDcols = paste0("scale_", j)]))) * sqrt(1/.SD[,2] + 1/.SD[,1] + wt_corr),,
        .SDcols = c(grep(paste0("count_e", j, "_s0"),names(work_data)),grep(paste0("count_e", j, "_s1"),names(work_data)))]
    }else{
      work_data[,paste0("cbe",j) := sqrt(1/.SD[,2] + 1/.SD[,1] + wt_corr),,
        .SDcols = c(grep(paste0("count_e", j, "_s0"),names(work_data)),grep(paste0("count_e", j, "_s1"),names(work_data)))]
    }
    #Set error of NA fitness values to NA
    work_data[is.na(get(paste0("fitness",j))),paste0("cbe",j) := NA]
  }

  #Calculate density of data along mean count based error
  bins <- 50
  work_data[,mean_cbe := rowMeans(.SD),.SDcols = grep(paste0("^cbe(",all_reps_str,")$"), names(work_data))]
  error_range <- seq(work_data[input_above_threshold == T & all_reads == T, log10(quantile(mean_cbe^2, probs=0.001))], 0, length.out = bins)
  work_data[,bin_error := findInterval(mean_cbe^2,vec = 10^error_range)]
  #weight variants according to how many other variants with same # of mutations are around
  if(dimsum_meta[["sequenceType"]]=="coding" & dimsum_meta[["mixedSubstitutions"]]){
    #variant counts already aggregated at AA level
    work_data[, error_model_weighting := sqrt(max(.N, sqrt(nrow(work_data)))), Nham_aa]
  }else{
    work_data[, error_model_weighting := sqrt(max(.N, sqrt(nrow(work_data)))), Nham_nt]
  }

  #Fit and write error model to file
  error_model <- dimsum__fit_error_model(
    dimsum_meta = dimsum_meta,
    input_dt = data.table::copy(work_data),
    all_reps = all_reps,
    norm_dt = fitness_norm_model)
  write.table(error_model, file = file.path(dimsum_meta[["tmp_path"]], "errormodel.txt"), row.names = F, quote = F)

  #Perform leave one out cross validation on replicates and generate QQ plot
  if(report){
    suppressWarnings(dimsum__error_model_qqplot(
      dimsum_meta = dimsum_meta,
      input_dt = work_data,
      all_reps = all_reps,
      norm_dt = fitness_norm_model,
      error_dt = error_model,
      report_outpath = report_outpath))
  }

  #Use error model parameters to calculate replicate-specific errors per variant
  for(j in all_reps){
    if(dimsum_meta[["fitnessNormalise"]]){
      Corr <- matrix(unlist(fitness_norm_model[,.SD,.SDcols = paste0("scale_", j)]), ncol = 1, nrow = nrow(work_data))
    }else{
      Corr <- matrix(1, ncol = 1, nrow = nrow(work_data))
    }

    #Full model
    work_data[,paste0("error", j) := sqrt(Corr * rowSums(matrix(unlist(error_model[parameter %in% c("input", "output") & rep == j, mean_value]), nrow = .N, ncol = 2, byrow = T)/.SD) +
      matrix(error_model[parameter %in% c("reperror") & rep == j, mean_value], nrow = .N, ncol = 1, byrow = T)),,
      .SDcols = c(
        grep(paste0("count_e", j, "_s0"), names(work_data)),
        grep(paste0("count_e", j, "_s1"), names(work_data)))]

    #Just replicate error
    work_data[,paste0("error_reponly", j) := sqrt(matrix(error_model[parameter %in% c("reperror") & rep == j, mean_value], nrow = .N, ncol = 1, byrow = T))]

    #Just input over-sequencing factor
    work_data[,paste0("error_inputosf", j) := sqrt(Corr * rowSums(matrix(unlist(error_model[parameter %in% c("input") & rep == j, mean_value]), nrow = .N, ncol = 1, byrow = T)/.SD)),,
      .SDcols = c(
        grep(paste0("count_e", j, "_s0"), names(work_data)))]

    #Just input (over-sequencing factor==1)
    work_data[,paste0("error_inputonly", j) := sqrt(Corr * rowSums(matrix(1, nrow = .N, ncol = 1, byrow = T)/.SD)),,
      .SDcols = c(
        grep(paste0("count_e", j, "_s0"), names(work_data)))]

    #Just output over-sequencing factor
    work_data[,paste0("error_outputosf", j) := sqrt(Corr * rowSums(matrix(unlist(error_model[parameter %in% c("output") & rep == j, mean_value]), nrow = .N, ncol = 1, byrow = T)/.SD)),,
      .SDcols = c(
        grep(paste0("count_e", j, "_s1"), names(work_data)))]

    #Just output (over-sequencing factor==1)
    work_data[,paste0("error_outputonly", j) := sqrt(Corr * rowSums(matrix(1, nrow = .N, ncol = 1, byrow = T)/.SD)),,
      .SDcols = c(
        grep(paste0("count_e", j, "_s1"), names(work_data)))]

  }

  ### Plot model parameters and fit
  ###########################

  if(report){
    #Add columns for upper/lower bounds of parameter estimates (for plotting)
    plot_error_model <- data.table::copy(error_model)
    # print(plot_error_model)

    #Error model limits
    rep_error_intercept <- mean(plot_error_model[parameter=="reperror",mean_value])
    mult_error_slope <- median(plot_error_model[parameter!="reperror",mean_value])

    #For plot, calculate average variance of fitness values per bin
    bs_data <- work_data[is.na(WT) & input_above_threshold == T & all_reads ==T]
    bs_data[,v := apply(.SD,1,var),.SDcols = grep(paste0("^fitness(", all_reps_str, ")$"), names(bs_data))]
    bs_data[,bin_mean_var := mean(v), bin_error]
    bs_data[,bin_N := .N, bin_error]
    bs_data[,bin_mean_error := mean(mean_cbe^2), bin_error]
    #Average merged error per bin
    bs_data[,bin_pred_error := mean(rowMeans(.SD^2, na.rm = T)), bin_error, .SDcols = paste0("error", all_reps)]
    bs_data[,bin_pred_error_reponly := mean(rowMeans(.SD^2, na.rm = T)), bin_error, .SDcols = paste0("error_reponly", all_reps)]
    bs_data[,bin_pred_error_inputosf := mean(rowMeans(.SD^2, na.rm = T)), bin_error, .SDcols = paste0("error_inputosf", all_reps)]
    bs_data[,bin_pred_error_inputonly := mean(rowMeans(.SD^2, na.rm = T)), bin_error, .SDcols = paste0("error_inputonly", all_reps)]
    bs_data[,bin_pred_error_outputosf := mean(rowMeans(.SD^2, na.rm = T)), bin_error, .SDcols = paste0("error_outputosf", all_reps)]
    bs_data[,bin_pred_error_outputonly := mean(rowMeans(.SD^2, na.rm = T)), bin_error, .SDcols = paste0("error_outputonly", all_reps)]
    melt_bs_data <- unique(data.table::melt(bs_data[,.SD,.SDcols = c(
      "bin_pred_error",
      "bin_pred_error_reponly",
      "bin_pred_error_inputosf",
      "bin_pred_error_inputonly",
      "bin_pred_error_outputosf",
      "bin_pred_error_outputonly",
      "bin_error",
      "bin_mean_error")], id.vars = c("bin_error","bin_mean_error")))

    #Plot1: input and output over-sequencing factor parameters +- sd
    a <- ggplot2::ggplot(plot_error_model[parameter %in% c("input","output")],
      ggplot2::aes(x=interaction(parameter, rep), mean_value, ymin = CI90_lower, ymax = CI90_upper, color = factor(rep))) +
      ggplot2::geom_pointrange(show.legend = F) +
      dimsum__scale_y_log10(limits = c(min(1, plot_error_model[parameter %in% c("input","output"), mean_value], plot_error_model[parameter %in% c("input", "output"), CI90_lower]),
        max(3, plot_error_model[parameter %in% c("input","output"), mean_value], plot_error_model[parameter %in% c("input", "output"), CI90_upper]))) +
      # ggplot2::scale_y_log10() +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
      ggplot2::labs(y = "Multiplicative\nerror terms", x = "Replicate (input or output)")

    #Plot2: replicate error parameters +- sd
    b <- ggplot2::ggplot(plot_error_model[parameter == "reperror"], ggplot2::aes(x = factor(rep), y = sqrt(mean_value), ymin = sqrt(CI90_lower), ymax = sqrt(CI90_upper), color = factor(rep))) +
      ggplot2::geom_pointrange(show.legend = F) +
      dimsum__scale_y_log10(limits = c(min(0.01, plot_error_model[parameter == "reperror", sqrt(CI90_lower)]),
        max(c(0.1, plot_error_model[parameter == "reperror", sqrt(CI90_upper)])))) +
      # ggplot2::scale_y_log10() +
      ggplot2::theme_bw() +
      ggplot2::labs(y = "Additive\nerror terms", x = "Replicate")

    #Plot3: count based error against variance of fitness + fit
    plot_cols <- dimsum__gg_color_hue(6)
    temp_mean_error <- unique(bs_data[,.(bin_mean_error, bin_mean_var, bin_N)])
    temp_mean_error[, ymin := bin_mean_var * (1-2/bin_N)]
    temp_mean_error[, ymax := bin_mean_var * (1+2/bin_N)]
    temp_mean_error[ymin < bs_data[,quantile(v, 0.001)], ymin := bs_data[,quantile(v, 0.001)]]
    temp_mean_error[ymax > bs_data[,quantile(v, 0.999)], ymax := bs_data[,quantile(v, 0.999)]]
    c <- ggplot2::ggplot(bs_data[v>=bs_data[,quantile(v, 0.001)] & v<=bs_data[,quantile(v, 0.999)],]) +
      ggplot2::stat_binhex(ggplot2::aes(mean_cbe^2, v), bins = 100, size = 0.2, color = "lightgrey") +
      dimsum__geom_abline(linetype = 2, size = 1) +
      # ggplot2::geom_abline(intercept = log10(mult_error_slope), linetype = 2, size = 1, color = plot_cols[2]) +
      # ggplot2::geom_hline(yintercept = rep_error_intercept, linetype = 2, size = 1, color = plot_cols[2]) +
      ggplot2::geom_pointrange(inherit.aes = F, data = temp_mean_error,
        ggplot2::aes(x = bin_mean_error, y = bin_mean_var, ymin = ymin, ymax = ymax), color = plot_cols[5]) +
      ggplot2::geom_line(inherit.aes = F, data = melt_bs_data[variable=="bin_pred_error",],
        ggplot2::aes(bin_mean_error, value), size = 1, color = plot_cols[1]) +
      dimsum__scale_y_log10(limits = c(bs_data[,quantile(v, 0.001)], bs_data[,quantile(v, 0.999)])) +
      dimsum__scale_x_log10() +
      ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
      ggplot2::coord_cartesian(xlim = c(min((bs_data[,mean_cbe])^2), 10^max(error_range))) +
      ggplot2::theme_bw() + ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::labs(x = "Count-based error estimate", y = "Variance of replicate fitness scores")

    #Plot4: error model decomposition
    d <- ggplot2::ggplot(melt_bs_data[variable=="bin_pred_error",]) +
      ggplot2::geom_line(inherit.aes = F,
        ggplot2::aes(bin_mean_error, value), size = 1, color = plot_cols[1]) +
      ggplot2::geom_line(inherit.aes = F, data = melt_bs_data[variable=="bin_pred_error_reponly",],
        ggplot2::aes(bin_mean_error, value), size = 1, color = plot_cols[3]) +
      ggplot2::geom_line(inherit.aes = F, data = melt_bs_data[variable=="bin_pred_error_inputosf",],
        ggplot2::aes(bin_mean_error, value), size = 1, color = plot_cols[4]) +
      ggplot2::geom_line(inherit.aes = F, data = melt_bs_data[variable=="bin_pred_error_inputonly",],
        ggplot2::aes(bin_mean_error, value), size = 0.3, color = plot_cols[4], linetype = 2) +
      ggplot2::geom_line(inherit.aes = F, data = melt_bs_data[variable=="bin_pred_error_outputosf",],
        ggplot2::aes(bin_mean_error, value), size = 1, color = plot_cols[6]) +
      ggplot2::geom_line(inherit.aes = F, data = melt_bs_data[variable=="bin_pred_error_outputonly",],
        ggplot2::aes(bin_mean_error, value), size = 0.3, color = plot_cols[6], linetype = 2) +
      dimsum__scale_y_log10(limits = c(bs_data[,quantile(v, 0.001)], bs_data[,quantile(v, 0.999)])) +
      dimsum__scale_x_log10() +
      ggplot2::coord_cartesian(xlim = c(min((bs_data[,mean_cbe])^2), 10^max(error_range))) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_text(size = 8),
        axis.title.y = ggplot2::element_text(size = 8),
        axis.text.x = ggplot2::element_text(size = 7),
        axis.text.y = ggplot2::element_text(size = 7)) +
      ggplot2::labs(x = "Count-based error estimate", y = "Variance of replicate fitness scores")

    #Plot5: legend
    legend_names <- c(
        "Empirical variance (average per bin)",
        "Poisson variance (null expectation)",
        "Full model variance",
        "Input multiplicative terms only\n(and null expectation)",
        "Output multiplicative terms only\n(and null expectation)",
        "Additive terms only")
    legend_data <- data.frame(
      x = rep(1:10, times = 6),
      y = rep(1:10, times = 6),
      variance_type = rep(legend_names, each = 10))
    legend_data[, "variance_type"] <- factor(legend_data[, "variance_type"], levels = legend_names)
    legend_cols <- c(
      plot_cols[5],
      "black",
      plot_cols[1],
      plot_cols[4],
      plot_cols[6],
      plot_cols[3])
    names(legend_cols) <- legend_names
    e <- ggplot2::ggplot(data = legend_data, ggplot2::aes(x,y)) +
      ggplot2::geom_point(ggplot2::aes(color=variance_type)) +
      ggplot2::geom_line(ggplot2::aes(color=variance_type)) +
      ggplot2::scale_colour_manual(name=c("Legend"), values=legend_cols, guide='legend') +
      ggplot2::guides(
        colour = ggplot2::guide_legend(override.aes = list(
          linetype=c(NA,2,1,1,1,1),
          size=c(2,0.75,0.75,0.75,0.75,0.75),
          shape=c(19, NA, NA, NA, NA, NA)))) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(
          fill = "white",
          colour = "white",
          color = "white"),
        plot.background=ggplot2::element_rect(
          fill = "white",
          colour = "white",
          color = "white"),
        rect = ggplot2::element_rect(
          fill = "white",
          colour = "white",
          color = "white"),
        legend.text = ggplot2::element_text(size = 9))
    e <- cowplot::get_legend(e)

    #Combine plots
    p <- gridExtra::arrangeGrob(a, b, c, d, e, nrow = 3,
      layout_matrix = rbind(c(1, 2), c(3, 4), c(3, 5)), heights = c(1, 1, 1), widths = c(2, 0.8))
    #Save
    dimsum_fitness__save_png_pdf(file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_repspec.png"), file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_repspec.pdf"), p, width = 9, height = 7)
  }

  dimsum__status_message("Done\n")

  #Render report
  dimsum__render_report(dimsum_meta = dimsum_meta)

  return(list(
    "error_model" = error_model,
    "norm_model" = fitness_norm_model))
}

# ---- dimsum__error_model_qqplot [modified, see the header and the generator]
dimsum__error_model_qqplot <- function(
  dimsum_meta,
  input_dt,
  all_reps,
  norm_dt,
  error_dt,
  report_outpath = NULL
  ){

  #Calculate error and fitness on training replicates and compare to fitness in leftout replicate
  for (i in seq_along(all_reps)) {

    training_reps <- paste0(all_reps[-i],collapse="|")
    training_reps_num <- all_reps[-i]
    NTreps <- length(training_reps_num)
    test_rep <- as.character(all_reps[i])
    training_data <- data.table::copy(input_dt)

    #Use error model parameters to calculate replicate-specific errors per variant
    for (j in seq_along(training_reps_num)) {
      if(dimsum_meta[["fitnessNormalise"]]){
        Corr <- matrix(unlist(norm_dt[,.SD,.SDcols = paste0("scale_", training_reps_num[j])]), ncol = 1, nrow = nrow(input_dt))
      }else{
        Corr <- matrix(1, ncol = 1, nrow = nrow(input_dt))
      }

      training_data[,paste0("error", training_reps_num[j]) := sqrt(Corr * rowSums(matrix(unlist(error_dt[parameter %in% c("input", "output") & rep == training_reps_num[j], mean_value]), nrow = .N, ncol = 2, byrow = T)/.SD) +
        matrix(error_dt[parameter %in% c("reperror") & rep == training_reps_num[j], mean_value], nrow = .N, ncol = 1, byrow = T)),,
        .SDcols = c(
          grep(paste0("count_e", training_reps_num[j], "_s0"), names(training_data)),
          grep(paste0("count_e", training_reps_num[j], "_s1"), names(training_data)))]
    }

    #Merged fitness values
    training_data[,fitness := rowSums(.SD[,1:NTreps]/(.SD[,(NTreps+1):(2*NTreps)]^2),na.rm=T) /
                    rowSums(1/(.SD[,(NTreps+1):(2*NTreps)]^2),na.rm=T),
                  ,.SDcols = c(grep(paste0("^fitness(", training_reps, ")$"),names(training_data)),
                               grep(paste0("^error(", training_reps, ")$"),names(training_data)))]
    #Merged error
    training_data[,error := sqrt(1/rowSums(1/.SD^2)),,
                  .SDcols = grep(paste0("^error(", training_reps, ")$"),names(training_data))]

    #Predict test replicate ###
    #Only use variants that have non-NA fitness and error values
    training_data[,test_rep_ok := !is.na(.SD) & is.finite(unlist(.SD)) & all_reads == T & input_above_threshold == T,,.SDcols = paste0("fitness",test_rep)]

    Mi_full <- error_dt[parameter == "input" & rep == test_rep,mean(mean_value)]
    Mo_full <- error_dt[parameter == "output" & rep == test_rep,mean(mean_value)]
    A_full <- error_dt[parameter == "reperror" & rep == test_rep, mean_value]

    #Scaling factor from fitness normalization
    if(dimsum_meta[["fitnessNormalise"]]){
      fitness_scale <- norm_dt[,unlist(.SD),.SDcols = paste0("scale_",test_rep)]
    }else{
      fitness_scale <- 1
    }

    #Set sequence hamming distance
    if(dimsum_meta[["sequenceType"]]=="coding" & dimsum_meta[["mixedSubstitutions"]]){
      training_data[, Nham_dist := Nham_aa]
    }else{
      training_data[, Nham_dist := Nham_nt]
    }

    dt_zscore <- training_data[test_rep_ok==T,
        .(test_rep, Nham_dist,
            MioA_full = (fitness-unlist(.SD[,1]))/sqrt(error^2 + fitness_scale*(Mi_full/unlist(.SD[,3]) + Mo_full/unlist(.SD[,4])) + A_full)),
       ,.SDcols = c(paste0("fitness",test_rep),
                    paste0("cbe",test_rep),
                    names(training_data)[grep(paste0("count_e", test_rep, "_s0"), names(training_data))],
                    names(training_data)[grep(paste0("count_e", test_rep, "_s1"), names(training_data))])]

    if (i == 1) {
      dt_zscore_all = dt_zscore
    } else {
      dt_zscore_all = rbind(dt_zscore_all, dt_zscore)
    }
  }

  #Plot results for dataset
  dt_zscore_melt <- data.table::melt(dt_zscore_all[sample(.N, min(1e4, .N))], id.vars = c("test_rep", "Nham_dist"))
  #Repeat data for all replicates combined
  dt_zscore_melt_all <- data.table::copy(dt_zscore_melt)
  dt_zscore_melt_all[, test_rep := "All"]
  dt_zscore_melt <- rbind(dt_zscore_melt, dt_zscore_melt_all)
  #Plot colours
  plot_cols <- c('darkgrey', dimsum__gg_color_hue(length(all_reps)))
  names(plot_cols) <- c("All", all_reps)
  zrange = 3.5
  p <- ggplot2::ggplot(dt_zscore_melt, ggplot2::aes(sample = value, color = test_rep)) +
    dimsum__geom_abline(lty = 2) +
    ggplot2::geom_qq(geom = "line") +
    ggplot2::scale_colour_manual(name=c("Replicate"), values=plot_cols, guide='legend') +
    ggplot2::coord_cartesian(xlim = c(-zrange, zrange), ylim = c(-zrange, zrange)) +
    ggplot2::labs(x = 'Theoretical quantiles (Normal)', y = 'Observed z-score quantiles') +
    ggplot2::theme_bw()
  dimsum_fitness__save_png_pdf(file.path(report_outpath,"dimsum_stage_fitness_report_1_errormodel_leaveoneout_qqplot.png"), file.path(report_outpath,"dimsum_stage_fitness_report_1_errormodel_leaveoneout_qqplot.pdf"), p, width = 6, height = 5)
}

# ---- dimsum__filter_nuc_variants
dimsum__filter_nuc_variants <- function(
  dimsum_meta,
  input_dt,
  all_reps
  ){

  dimsum__status_message("Filtering out low count nucleotide variants...\n")

  #Number of input and output replicates
  all_reps_str <- paste0(all_reps, collapse="|")

  #Sample names
  input_samples <- names(input_dt)[grep(paste0("e(", all_reps_str, ")_s0_b.*_count$"), names(input_dt))]
  output_samples <- names(input_dt)[grep(paste0("e(", all_reps_str, ")_s1_b.*_count$"), names(input_dt))]

  #### Retain variants with input/output read counts > "fitnessMinInputCountAny" in ANY biological replicates
  #### Retain variants with input/output read counts > "fitnessMinInputCountAll" in ALL biological replicates
  output_dt <- copy(input_dt)
  #fitnessMinInputCountAny
  if(is.list(dimsum_meta[["fitnessMinInputCountAny"]])){
    for(i in names(dimsum_meta[["fitnessMinInputCountAny"]])){
      output_dt <- output_dt[(Nham_nt == as.integer(i) & rowSums(output_dt[,input_samples,with=F]>=dimsum_meta[["fitnessMinInputCountAny"]][[i]]) != 0) | Nham_nt != as.integer(i)]
    }
  }else{
    output_dt <- output_dt[rowSums(output_dt[,input_samples,with=F]>=dimsum_meta[["fitnessMinInputCountAny"]]) != 0]
  }
  #fitnessMinOutputCountAny
  if(is.list(dimsum_meta[["fitnessMinOutputCountAny"]])){
    for(i in names(dimsum_meta[["fitnessMinOutputCountAny"]])){
      output_dt <- output_dt[(Nham_nt == as.integer(i) & rowSums(output_dt[,output_samples,with=F]>=dimsum_meta[["fitnessMinOutputCountAny"]][[i]]) != 0) | Nham_nt != as.integer(i)]
    }
  }else{
    output_dt <- output_dt[rowSums(output_dt[,output_samples,with=F]>=dimsum_meta[["fitnessMinOutputCountAny"]]) != 0]
  }
  #fitnessMinInputCountAll
  if(is.list(dimsum_meta[["fitnessMinInputCountAll"]])){
    for(i in names(dimsum_meta[["fitnessMinInputCountAll"]])){
      output_dt <- output_dt[(Nham_nt == as.integer(i) & rowSums(output_dt[,input_samples,with=F]<dimsum_meta[["fitnessMinInputCountAll"]][[i]]) == 0) | Nham_nt != as.integer(i)]
    }
  }else{
    output_dt <- output_dt[rowSums(output_dt[,input_samples,with=F]<dimsum_meta[["fitnessMinInputCountAll"]]) == 0]
  }
  #fitnessMinOutputCountAll
  if(is.list(dimsum_meta[["fitnessMinOutputCountAll"]])){
    for(i in names(dimsum_meta[["fitnessMinOutputCountAll"]])){
      output_dt <- output_dt[(Nham_nt == as.integer(i) & rowSums(output_dt[,output_samples,with=F]<dimsum_meta[["fitnessMinOutputCountAll"]][[i]]) == 0) | Nham_nt != as.integer(i)]
    }
  }else{
    output_dt <- output_dt[rowSums(output_dt[,output_samples,with=F]<dimsum_meta[["fitnessMinOutputCountAll"]]) == 0]
  }

  #### Set input read counts < "fitnessMinInputCountAny" to NA
  for(i in input_samples){
    if(is.list(dimsum_meta[["fitnessMinInputCountAny"]])){
      for(j in names(dimsum_meta[["fitnessMinInputCountAny"]])){
        output_dt[Nham_nt == as.integer(j) & get(i) < dimsum_meta[["fitnessMinInputCountAny"]][[j]], as.character(i) := NA]
      }
    }else{
      output_dt[get(i) < dimsum_meta[["fitnessMinInputCountAny"]], as.character(i) := NA]
    }
  }

  #### Set output read counts < "fitnessMinOutputCountAny" to NA
  for(i in output_samples){
    if(is.list(dimsum_meta[["fitnessMinOutputCountAny"]])){
      for(j in names(dimsum_meta[["fitnessMinOutputCountAny"]])){
        output_dt[Nham_nt == as.integer(j) & get(i) < dimsum_meta[["fitnessMinOutputCountAny"]][[j]], as.character(i) := NA]
      }
    }else{
      output_dt[get(i) < dimsum_meta[["fitnessMinOutputCountAny"]], as.character(i) := NA]
    }
  }

  dimsum__check_variants(dimsum_meta = dimsum_meta, input_dt = output_dt)

  dimsum__status_message("Done\n")

  return(output_dt)

}

# ---- dimsum__fit_error_model [modified, see the header and the generator]
dimsum__fit_error_model <- function(
  dimsum_meta,
  input_dt,
  all_reps,
  norm_dt,
  Nbootstraps = 100,
  maxN = 10000,
  max_tries_per_fit = 20,
  lower_rep = 1e-4
  ){

  #Number of input and output replicates
  all_reps_str <- paste0(all_reps, collapse="|")

  #Fitness correction for normalisation (scaling) factors
  Fcorr <- NULL
  if(dimsum_meta[["fitnessNormalise"]]){
    Fcorr <- unlist(norm_dt[,.SD,.SDcols = grep(paste0("scale_(", all_reps_str, ")$"), names(norm_dt))])
  }

  #Set up fitting parameters
  Nreps <- length(all_reps)
  #Calculate all combinations of replicates of length >= 2
  idx <- list()
  for(i in (length(all_reps)):2){
    idx <- c(idx, combn(length(all_reps), i, function(x) list(x)))
  }

  #Subset to maximum 500 combinations (preferably highest order combinations)
  if(length(idx)>500){
    idx <- idx[1:500]
  }

  if(dimsum_meta[['numCores']] == 1){
    #One core: the bootstraps run in this process, from the random-number stream that
    #clusterSetRNGStream(cl, 1234567) hands the only worker of a one-core cluster - the same draws in
    #the same order, without starting a second R process and sending it the data
    parameters <- t(dimsum_fitness__with_worker_stream(1234567, sapply(
      X = 1:Nbootstraps,
      FUN = dimsum__fit_error_model_bootstrap,
      input_dt = input_dt,
      all_reps = all_reps,
      idx = idx,
      maxN = maxN,
      max_tries_per_fit = max_tries_per_fit,
      lower_rep = lower_rep,
      Fcorr = Fcorr)))
  }else{
    #Setup cluster
    clust <- parallel::makeCluster(dimsum_meta[['numCores']])
    #Set cluster seed
    parallel::clusterSetRNGStream(cl = clust, 1234567)
    #Fit
    parameters <- t(parallel::parSapply(
      cl = clust,
      X = 1:Nbootstraps,
      FUN = dimsum__fit_error_model_bootstrap,
      input_dt = input_dt,
      all_reps = all_reps,
      idx = idx,
      maxN = maxN,
      max_tries_per_fit = max_tries_per_fit,
      lower_rep = lower_rep,
      Fcorr = Fcorr))
    #Stop cluster
    parallel::stopCluster(clust)
  }

  #Return table with error model parameters
  error_model <- data.table(
    parameter = rep(c("input", "output", "reperror"), each = Nreps),
    rep = rep(all_reps, 3),
    mean_value = colMeans(parameters, na.rm = T),
    CI90_lower = apply(parameters, 2, quantile, na.rm = T, probs = 0.1),
    CI90_upper = apply(parameters, 2, quantile, na.rm = T, probs = 0.9),
    ensemble = sum(!is.na(parameters[,1])))
  return(error_model)
}

# ---- dimsum__fit_error_model_bootstrap [modified, see the header and the generator]
dimsum__fit_error_model_bootstrap <- function(
  m,
  input_dt,
  all_reps,
  idx,
  maxN,
  max_tries_per_fit,
  lower_rep,
  Fcorr = NULL
  ){
  require(data.table)

  Nreps <- length(all_reps)
  y <- NA
  counter <- 0
  while(!is.list(y) & counter < max_tries_per_fit){ #this is a routine to catch failed fitting attempts
    tryCatch({

      #Create data concatenation for all combinations of replicates (can be improved?!)
      bs_data <- input_dt[input_above_threshold == T & all_reads == T & is.na(WT)][sample(.N,min(.N,maxN),replace = T)]
      F_data_list <- list() #fitness data
      E_data_list <- list() #count based error data (for weighting of datapoints)
      C_data_list <- list() #inverse count data
      NR <- list() # simple counter for how many replicates are in each data row
      for(i in 1:length(idx)){
        F_data_list[[i]] <- bs_data[,.SD,.SDcols = grep(paste0("^fitness(", paste0(all_reps[idx[[i]]], collapse = "|"), ")$"), names(input_dt))]
        E_data_list[[i]] <- bs_data[,.SD,.SDcols = grep(paste0("^cbe(", paste0(all_reps[idx[[i]]], collapse = "|"), ")$"), names(input_dt))]
        C_data_list[[i]] <- bs_data[,1/.SD,.SDcols = grep(paste0("count_e(", paste0(all_reps[idx[[i]]], collapse = "|"), ")_s"), names(input_dt))]
        NR[[i]] <- data.table(rep(length(idx[[i]]),nrow(bs_data)))
      }

      #make matrices with #replicates = #columns
      F_data <- as.matrix(rbindlist(F_data_list, fill = T))
      #Calculate variance from fitness data
      #(row by row through the C routine that var(x, na.rm = T) calls - .Call(C_cov, x, NULL,
      #"na.or.complete", FALSE) - so the same arithmetic and the same result as
      #apply(F_data, 1, var, na.rm = T), without var()'s argument handling for each of the rows)
      C_cov <- getFromNamespace("C_cov", "stats")
      tF_data <- t(F_data)
      V_data <- vapply(seq_len(ncol(tF_data)), function(k){.Call(C_cov, tF_data[,k], NULL, 5L, FALSE)}, numeric(1))

      #count based error weighting of data points (needed because deviation is proportional to expcectation)
      E_data <- as.matrix(rbindlist(E_data_list, fill = T))
      Ew <- rowMeans(E_data, na.rm = T)^2

      C_data <- as.matrix(rbindlist(C_data_list, fill = T))
      #Avoid NA for fitting
      C_data[is.na(C_data)] <- 0
      #Correct count terms for fitness normalization factors
      if(!is.null(Fcorr)){
        C_data = C_data * matrix(rep(Fcorr, 2), nrow = nrow(C_data), ncol = Nreps*2, byrow = T)
      }

      #weighting according to #variants with same mutation
      Dw <- bs_data$error_model_weighting
      #how many replicates are in the replicate combinations
      NRT <- as.matrix(rbindlist(NR))

      #Binary variable to avoid replicate error fitting for replicates that are not present in a data row
      BV <- F_data
      BV[!is.na(BV)] <- 1
      BV[is.na(BV)] <- 0

      #fit formula (the columns as plain vectors BV1.., C1.. rather than BV[,i] and C_data[,i], which
      #were extracted anew at every evaluation of the model; the same values and the same arithmetic)
      for(k in 1:Nreps){assign(paste0("BV", k), BV[,k])}
      for(k in 1:(2*Nreps)){assign(paste0("C", k), C_data[,k])}
      fit_formula <- as.formula(paste0(
        "V_data ~ (",
        paste0("BV",
          1:Nreps,
          " * (p[",
          1:Nreps,
          "] * C",
          1:Nreps,
          " +  p[",
          Nreps + 1:Nreps,
          "] * C",
          Nreps + 1:Nreps,
          " + p[",
          2*Nreps + 1:Nreps,
          "])",collapse = " + "),
        ") / NRT"))

      #Fitting
      y <- nls(formula = fit_formula,
            start = list(p = c(
              rep(1, each = Nreps*2) + 10^rnorm(2*Nreps),
              rep(0, each = Nreps) + 0.01*10^rnorm(1*Nreps))),
            lower = c(rep(1, each = Nreps*2), rep(lower_rep, each = Nreps)),
            weights = 1 / Ew / Dw,
            algorithm = "port")
    }, error = function(cond){})
    counter <- counter + 1}

  #If fit was succesful after less than max_tries_per_fit attempts, transfer parameters to p
  if(counter < max_tries_per_fit){
    p <- summary(y)[["parameters"]][,1]
  }else{ #otherwise move on to next fit
    p <- rep(NA, 3*Nreps)
  }
  return(p)
}

# ---- dimsum__get_barcode_design
dimsum__get_barcode_design <- function(
  dimsum_meta
  ){
  #Check if barcodeDesignPath specified
  if(is.null(dimsum_meta[["barcodeDesignPath"]])){
    return(NULL)
  }

  #Load barcode design
  if(!file.exists(dimsum_meta[["barcodeDesignPath"]])){
    stop(paste0("Invalid '", "barcodeDesignPath", "' argument (file not found)"), call. = FALSE)
  }
  barcode_design <- read.table(dimsum_meta[["barcodeDesignPath"]], header = T, stringsAsFactors = F, sep="\t")

  #Set pair2 column equal to pair1 column if single-end library (contents of existing pair2 column will be ignored)
  if(!dimsum_meta[["paired"]]){
    if(!"pair1" %in% colnames(barcode_design)){
      stop(paste0("Mandatory column missing from barcodeDesign file ('pair1')"), call. = FALSE)
    }else{
      barcode_design[,"pair2"] <- barcode_design[,"pair1"]
    }
  }

  #Add original FASTQ directory
  if("fastqFileDir" %in% names(barcode_design)){
    names(barcode_design)[names(barcode_design)=="fastqFileDir"] <- "pair_directory"
  }else{
    barcode_design[,"pair_directory"] <- dimsum_meta[["fastqFileDir"]]
  }

  #Set barcode1 column equal to barcode column (backwards compatibility)
  if("barcode" %in% colnames(barcode_design)){
    barcode_design[,"barcode1"] <- barcode_design[,"barcode"]
  }

  #Set barcode1 column equal to barcode column (backwards compatibility)
  if(!"barcode2" %in% colnames(barcode_design)){
    barcode_design[,"barcode2"] <- barcode_design[,"barcode1"]
  }

  #Check whether barcode design is valid
  dimsum__check_barcode_design(barcode_design)

  #Check FASTQ files exist
  #Pair1 files
  for(i in unlist(file.path(barcode_design[,"pair_directory"], barcode_design[,"pair1"]))){
    if(!file.exists(i)){
      stop(paste0("Invalid FASTQ file name '", i, "' in barcodeDesign file (file not found)"), call. = FALSE)
    }
  }
  #Pair2 files
  for(i in unlist(file.path(barcode_design[,"pair_directory"], barcode_design[,"pair2"]))){
    if(!file.exists(i)){
      stop(paste0("Invalid FASTQ file name '", i, "' in barcodeDesign file (file not found)"), call. = FALSE)
    }
  }

  return(barcode_design)
}

# ---- dimsum__get_experiment_design [modified, see the header and the generator]
dimsum__get_experiment_design <- function(
  dimsum_meta
  ){
  #Load experimental design
  if(!file.exists(dimsum_meta[["experimentDesignPath"]])){
    stop(paste0("Invalid '", "experimentDesignPath", "' argument (file not found)"), call. = FALSE)
  }
  exp_design <- read.table(dimsum_meta[["experimentDesignPath"]], header = T, stringsAsFactors = F, sep="\t")
  #Remove rows with missing sample_name
  exp_design <- exp_design[exp_design[,"sample_name"]!="",]

  #Clear pair1, pair2 and technical_replicate columns if count file provided
  if(!is.null(dimsum_meta[["countPath"]])){
    exp_design[,"technical_replicate"] <- NA
    exp_design[,"pair1"] <- NA
    exp_design[,"pair2"] <- NA
  }

  #Set pair2 column equal to pair1 column if single-end library (contents of existing pair2 column will be ignored)
  if(!dimsum_meta[["paired"]]){
    if(!"pair1" %in% colnames(exp_design)){
      stop(paste0("Mandatory column missing from experimentDesign file ('pair1')"), call. = FALSE)
    }else{
      exp_design[,"pair2"] <- exp_design[,"pair1"]
    }
  }

  #Add original FASTQ directory
  if("fastqFileDir" %in% names(exp_design)){
    names(exp_design)[names(exp_design)=="fastqFileDir"] <- "pair_directory"
  }else{
    exp_design[,"pair_directory"] <- dimsum_meta[["fastqFileDir"]]
  }

  #Add sample-specific cutadapt options
  # if((!"cutadapt5First" %in% colnames(exp_design) & is.null(dimsum_meta[["cutadapt5First"]])) | (!"cutadapt5Second" %in% colnames(exp_design) & is.null(dimsum_meta[["cutadapt5Second"]]))){
  #   stop("Sequence of 5' constant region not found. Please check that the 'cutadapt5First' and 'cutadapt5Second' arguments are correctly set.", call. = FALSE)
  # }
  if(!"cutadaptCut5First" %in% colnames(exp_design)){exp_design[,"cutadaptCut5First"] <- ifelse(is.null(dimsum_meta[["cutadaptCut5First"]]), NA, dimsum_meta[["cutadaptCut5First"]])}
  if(!"cutadaptCut5Second" %in% colnames(exp_design)){exp_design[,"cutadaptCut5Second"] <- ifelse(is.null(dimsum_meta[["cutadaptCut5Second"]]), NA, dimsum_meta[["cutadaptCut5Second"]])}
  if(!"cutadaptCut3First" %in% colnames(exp_design)){exp_design[,"cutadaptCut3First"] <- ifelse(is.null(dimsum_meta[["cutadaptCut3First"]]), NA, dimsum_meta[["cutadaptCut3First"]])}
  if(!"cutadaptCut3Second" %in% colnames(exp_design)){exp_design[,"cutadaptCut3Second"] <- ifelse(is.null(dimsum_meta[["cutadaptCut3Second"]]), NA, dimsum_meta[["cutadaptCut3Second"]])}
  if(!"cutadapt5First" %in% colnames(exp_design)){exp_design[,"cutadapt5First"] <- ifelse(is.null(dimsum_meta[["cutadapt5First"]]), NA, dimsum_meta[["cutadapt5First"]])}
  if(!"cutadapt5Second" %in% colnames(exp_design)){exp_design[,"cutadapt5Second"] <- ifelse(is.null(dimsum_meta[["cutadapt5Second"]]), NA, dimsum_meta[["cutadapt5Second"]])}
  if(!"cutadapt3First" %in% colnames(exp_design)){exp_design[,"cutadapt3First"] <- ifelse(is.null(dimsum_meta[["cutadapt3First"]]), NA, dimsum_meta[["cutadapt3First"]])}
  if(!"cutadapt3Second" %in% colnames(exp_design)){exp_design[,"cutadapt3Second"] <- ifelse(is.null(dimsum_meta[["cutadapt3Second"]]), NA, dimsum_meta[["cutadapt3Second"]])}
  if(!"cutadaptMinLength" %in% colnames(exp_design)){exp_design[,"cutadaptMinLength"] <- ifelse(is.null(dimsum_meta[["cutadaptMinLength"]]), NA, dimsum_meta[["cutadaptMinLength"]])}
  if(!"cutadaptErrorRate" %in% colnames(exp_design)){exp_design[,"cutadaptErrorRate"] <- ifelse(is.null(dimsum_meta[["cutadaptErrorRate"]]), NA, dimsum_meta[["cutadaptErrorRate"]])}
  if(!"cutadaptOverlap" %in% colnames(exp_design)){exp_design[,"cutadaptOverlap"] <- ifelse(is.null(dimsum_meta[["cutadaptOverlap"]]), NA, dimsum_meta[["cutadaptOverlap"]])}
  if(!"generations" %in% colnames(exp_design)){exp_design[,"generations"] <- NA}
  if(!"cell_density" %in% colnames(exp_design)){exp_design[,"cell_density"] <- NA}
  if(!"selection_time" %in% colnames(exp_design)){exp_design[,"selection_time"] <- NA}
  #Convert empty string constant region sequences to NA
  exp_design[which(exp_design[,"cutadapt5First"]==""),"cutadapt5First"] <- NA
  exp_design[which(exp_design[,"cutadapt5Second"]==""),"cutadapt5Second"] <- NA
  exp_design[which(exp_design[,"cutadapt3First"]==""),"cutadapt3First"] <- NA
  exp_design[which(exp_design[,"cutadapt3Second"]==""),"cutadapt3Second"] <- NA

  ### Backwards compatibility
  #Copy "experiment" and "biological_replicate" columns to "experiment_replicate" and "selection_replicate" respectively
  if("experiment" %in% colnames(exp_design)){exp_design[,"experiment_replicate"] <- exp_design[,"experiment"]}
  if("biological_replicate" %in% colnames(exp_design)){exp_design[,"selection_replicate"] <- exp_design[,"biological_replicate"]}

  #Check whether experiment design is valid
  if(!is.null(dimsum_meta[["countPath"]])){
    dimsum__check_experiment_design_countfile(dimsum_meta, exp_design)
  }else{
    dimsum__check_experiment_design(dimsum_meta, exp_design)
  }

  #Linked adapters override non-linked adapters (and 5' linked adapter takes preference over 3' linked adapter if both specified)
  exp_design[grepl("\\.\\.\\.", exp_design[,"cutadapt5First"]),"cutadapt3First"] <- exp_design[grepl("\\.\\.\\.", exp_design[,"cutadapt5First"]),"cutadapt5First"]
  exp_design[grepl("\\.\\.\\.", exp_design[,"cutadapt3First"]),"cutadapt5First"] <- NA
  exp_design[grepl("\\.\\.\\.", exp_design[,"cutadapt5Second"]),"cutadapt3Second"] <- exp_design[grepl("\\.\\.\\.", exp_design[,"cutadapt5Second"]),"cutadapt5Second"]
  exp_design[grepl("\\.\\.\\.", exp_design[,"cutadapt3Second"]),"cutadapt5Second"] <- NA

  #Linked adapters are both required by default - read1
  lFirst <- grepl("\\.\\.\\.", exp_design[,"cutadapt3First"])
  l5First_required <- grepl(";required\\.\\.\\.", exp_design[,"cutadapt3First"])
  l5First_optional <- grepl(";optional\\.\\.\\.", exp_design[,"cutadapt3First"])
  l3First_required <- grepl(";required$", exp_design[,"cutadapt3First"])
  l3First_optional <- grepl(";optional$", exp_design[,"cutadapt3First"])
  exp_design[lFirst & !l5First_required & !l5First_optional,"cutadapt3First"] <- gsub("\\.\\.\\.", ";required\\.\\.\\.", exp_design[lFirst & !l5First_required & !l5First_optional,"cutadapt3First"])
  exp_design[lFirst & !l3First_required & !l3First_optional,"cutadapt3First"] <- paste0(exp_design[lFirst & !l3First_required & !l3First_optional,"cutadapt3First"], ";required")

  #Linked adapters are both required by default - read2
  lSecond <- grepl("\\.\\.\\.", exp_design[,"cutadapt3Second"])
  l5Second_required <- grepl(";required\\.\\.\\.", exp_design[,"cutadapt3Second"])
  l5Second_optional <- grepl(";optional\\.\\.\\.", exp_design[,"cutadapt3Second"])
  l3Second_required <- grepl(";required$", exp_design[,"cutadapt3Second"])
  l3Second_optional <- grepl(";optional$", exp_design[,"cutadapt3Second"])
  exp_design[lSecond & !l5Second_required & !l5Second_optional,"cutadapt3Second"] <- gsub("\\.\\.\\.", ";required\\.\\.\\.", exp_design[lSecond & !l5Second_required & !l5Second_optional,"cutadapt3Second"])
  exp_design[lSecond & !l3Second_required & !l3Second_optional,"cutadapt3Second"] <- paste0(exp_design[lSecond & !l3Second_required & !l3Second_optional,"cutadapt3Second"], ";required")

  #Indicate for which samples cutadapt should be run
  num_cutadapt_options <- apply(!is.na(exp_design[,c("cutadapt5First", "cutadapt5Second", "cutadapt3First", "cutadapt3Second")]), 1, sum)
  num_cutadaptCut_options <- apply(!is.na(exp_design[,c("cutadaptCut5First", "cutadaptCut5Second", "cutadaptCut3First", "cutadaptCut3Second")]), 1, sum)
  exp_design[,"run_cutadapt"] <- num_cutadapt_options | num_cutadaptCut_options
  #Indicate for which samples cutadapt should be run in "Cut" mode only
  exp_design[,"run_cutadapt_cutonly"] <- !num_cutadapt_options & num_cutadaptCut_options

  ### Temporary fix for stages relying on "experiment" and "biological_replicate" columns
  #Copy "experiment_replicate" and "selection_replicate" to "experiment" and "biological_replicate" colums respectively
  exp_design[,"experiment"] <- exp_design[,"experiment_replicate"]
  exp_design[,"biological_replicate"] <- exp_design[,"selection_replicate"]

  #Check that each sample has 5' adapters (constant regions) specified (if unstranded library)
  if(!dimsum_meta[["stranded"]]){
    if(sum(is.na(exp_design[,"cutadapt5First"]) & is.na(exp_design[,"cutadapt3First"]))!=0 | sum(is.na(exp_design[,"cutadapt5Second"]) & is.na(exp_design[,"cutadapt3Second"]))!=0){
      stop("Unstranded library but constant regions not found for one or more FASTQ files. Please check that the corresponding experiment design file columns are correct.", call. = FALSE)
    }
  }

  #Check FASTQ files exist (if no count file provided and demultiplexed FASTQ files supplied i.e. no barcodeDesignPath supplied)
  if(is.null(dimsum_meta[["barcodeDesignPath"]]) & is.null(dimsum_meta[["countPath"]])){
    #Pair1 files
    for(i in unlist(file.path(exp_design[,"pair_directory"], exp_design[,"pair1"]))){
      if(!file.exists(i)){
        stop(paste0("Invalid FASTQ file name '", i, "' in experimentDesign file (file not found)"), call. = FALSE)
      }
    }
    #Pair2 files
    for(i in unlist(file.path(exp_design[,"pair_directory"], exp_design[,"pair2"]))){
      if(!file.exists(i)){
        stop(paste0("Invalid FASTQ file name '", i, "' in experimentDesign file (file not found)"), call. = FALSE)
      }
    }
  }

  #Add .gz extension to demultiplexed FASTQ file names for backwards compatibility
  if(!is.null(dimsum_meta[["barcodeDesignPath"]])){
    exp_design[,"pair1"] = paste0(exp_design[,"pair1"], c(".gz", "")[as.numeric(grepl(".gz$", exp_design[,"pair1"]))+1])
    exp_design[,"pair2"] = paste0(exp_design[,"pair2"], c(".gz", "")[as.numeric(grepl(".gz$", exp_design[,"pair2"]))+1])
  }

  #Check that all FASTQ file prefices exist in barcodeDesign file (if barcodeDesignPath supplied)
  if(!is.null(dimsum_meta[["barcodeDesignPath"]])){
    all_prefices <- unique(gsub("1.fastq.gz|2.fastq.gz", "", unlist(exp_design[,c("pair1", "pair2")])))
    if(sum(!all_prefices %in% dimsum_meta[["barcode_design"]][,"new_pair_prefix"])!=0){
      stop(paste0("One or more FASTQ file names in experimentDesign file didn't match a 'new_pair_prefix' value in barcodeDesign file"), call. = FALSE)
    }
  }

  #If not trans library: reverse complement cutadapt 5' constant regions to obtain 3' constant regions (if not already supplied)
  if(!dimsum_meta[["transLibrary"]]){
    if(sum(is.na(exp_design[,"cutadapt3First"]) & !is.na(exp_design[,"cutadapt5Second"]))!=0){
      exp_design[is.na(exp_design[,"cutadapt3First"]) & !is.na(exp_design[,"cutadapt5Second"]),"cutadapt3First"] <- dimsum_fitness__reverse_complement(exp_design[is.na(exp_design[,"cutadapt3First"]) & !is.na(exp_design[,"cutadapt5Second"]),"cutadapt5Second"])
    }
    if(sum(is.na(exp_design[,"cutadapt3Second"]) & !is.na(exp_design[,"cutadapt5First"]))!=0){
      exp_design[is.na(exp_design[,"cutadapt3Second"]) & !is.na(exp_design[,"cutadapt5First"]),"cutadapt3Second"] <- dimsum_fitness__reverse_complement(exp_design[is.na(exp_design[,"cutadapt3Second"]) & !is.na(exp_design[,"cutadapt5First"]),"cutadapt5First"])
    }
  }

  #Check retainedReplicates exist in experiment design file
  if(dimsum_meta[["retainedReplicates"]]!="all"){
    all_reps <- as.integer(strsplit(dimsum_meta[["retainedReplicates"]], ",")[[1]])
    if(sum(!all_reps %in% exp_design[,"experiment"])!=0){
      stop("Invalid 'retainedReplicates' argument. One or more replicates not present in experiment design file.", call. = FALSE)
    }
  }

  return(exp_design)
}

# ---- dimsum__get_synonymsequences [modified, see the header and the generator]
dimsum__get_synonymsequences <- function(
  dimsum_meta
  ){

  ### Abort if no count file supplied
  if(is.null(dimsum_meta[["synonymSequencePath"]])){return(NULL)}

  input_dt <- data.table::fread(dimsum_meta[["synonymSequencePath"]], header = F)

  ### Nucleotide sequence checks (single column)
  #Check if just a single columns present
  if(nrow(input_dt)==0 | ncol(input_dt)!=1){
    stop(paste0("Incorrect number of columns in file specified by synonymSequencePath (only one allowed)"), call. = FALSE)
  }
  #Check nucleotide sequence column is of type character
  if(sapply(input_dt, typeof)["V1"]!="character"){
    stop("One or more invalid sequences in synonym sequences file specified by 'synonymSequencePath'. Only valid coding nucleotide sequences allowed (A/C/T/G).", call. = FALSE)
  }
  #Set nucleotide sequence to lower case
  input_dt[, V1 := tolower(V1)]
  #Check nucleotide sequences are valid (ACGT characters only)
  if(sum(!input_dt[,unique(unlist(strsplit(V1, "")))] %in% c('a', 'c', 'g', 't'))!=0){
    stop("One or more invalid sequences in synonym sequences file specified by 'synonymSequencePath'. Only valid coding nucleotide sequences allowed (A/C/T/G).", call. = FALSE)
  }

  #Check if sequences are coding
  if(sum(unlist(lapply(as.list(input_dt[,V1]), dimsum__detect_sequence_type))!='coding')!=0){
    stop("One or more non-coding sequences in synonym sequences file specified by 'synonymSequencePath'. Only valid coding nucleotide sequences allowed (A/C/T/G).", call. = FALSE)
  }

  #Translate sequences
  input_dt[, aa_seq := dimsum_fitness__translate(V1)]

  return(input_dt[,aa_seq])
}

# ---- dimsum__gg_color_hue
dimsum__gg_color_hue <- function(
  n
  ){
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# ---- dimsum__ggpairs_binhex [modified, see the header and the generator]
dimsum__ggpairs_binhex <- function(
  input_dt,
  output_file_prefix,
  width = 12,
  height = 12,
  bins = 50,
  xlab = "x",
  ylab = "y",
  title = "",
  cut = NULL,
  size = 0.5,
  thresholds = NULL){
  #Check if something to plot
  if(dim(input_dt)[1] == 0){
    warning("dimsum__ggpairs_binhex.R: No data to plot (empty data.table 'input_dt').", call. = FALSE, immediate. = TRUE, noBreaks. = TRUE)
    return(NULL)
  }
  d <- GGally::ggpairs(input_dt,
    columns = 1:dim(input_dt)[2],
    upper=list(continuous = "cor"),
    lower="blank", xlab = xlab, ylab = ylab, title = title)
  if(is.null(cut)){
    plot_dt <- input_dt
  }else{
    plot_dt <- input_dt[, cut := cut]
  }
  for (x in 1:dim(input_dt)[2]){
    for (y in 1:dim(input_dt)[2]){
      if (y>x) {
        dpiece <- ggplot2::ggplot(plot_dt, ggplot2::aes_string(x = colnames(input_dt)[x], y = colnames(input_dt)[y]))
        if(!is.null(thresholds)){
          dpiece <- dpiece +
            ggplot2::geom_hline(yintercept = unlist(thresholds), linetype = unlist(names(thresholds)), color = "darkgrey") +
            ggplot2::geom_vline(xintercept = unlist(thresholds), linetype = unlist(names(thresholds)), color = "darkgrey")
        }
        if(is.null(cut)){
          dpiece <- dpiece + ggplot2::stat_binhex(bins = bins, size = size, color = "lightgrey")
        }else{
          dpiece <- dpiece + ggplot2::stat_binhex(ggplot2::aes(color = cut), bins = bins, size = size)
        }
        dpiece <- dpiece +
          ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
          dimsum__geom_abline(color = "darkgrey", lty = 2)
      d <- GGally::putPlot(d, dpiece, y, x)
      }
    }
  }
  d <- d + ggplot2::theme_bw() + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  suppressWarnings(suppressMessages(dimsum_fitness__save_png_pdf(paste0(output_file_prefix, ".png"), paste0(output_file_prefix, ".pdf"), d, width = width, height = height)))
}

# ---- dimsum__hamming_distance
dimsum__hamming_distance <- function(x, y){
  return(sum(strsplit(x, "")[[1]]!=strsplit(y, "")[[1]]))
}

# ---- dimsum__identify_STOP_mutations
dimsum__identify_STOP_mutations <- function(
  input_dt
  ){

  #WT AA sequence
  wt_AAseq <- input_dt[WT==T,aa_seq]

  #If WT has terminating STOP remove last residue in all sequences
  if(substr(wt_AAseq, nchar(wt_AAseq), nchar(wt_AAseq))=="*"){
    input_dt[, aa_seq_ns := substr(aa_seq, 1, nchar(aa_seq)-1)]
  }else{
    input_dt[, aa_seq_ns := aa_seq]
  }
  #Indicate premature STOP (regardless of presence in WT)
  input_dt[,STOP := ifelse(length(grep(aa_seq_ns,pattern="\\*"))==1,TRUE,FALSE),aa_seq_ns]

  #Indicate readthrough STOP mutations (if terminating STOP exists in WT)
  input_dt[, STOP_readthrough := FALSE]
  if(substr(wt_AAseq, nchar(wt_AAseq), nchar(wt_AAseq))=="*"){
    input_dt[, STOP_readthrough := (substr(aa_seq, nchar(aa_seq), nchar(aa_seq))!="*")]
  }

  return(input_dt[,.SD,,.SDcols = names(input_dt)[names(input_dt)!="aa_seq_ns"]])

}

# ---- dimsum__identify_double_aa_mutations
dimsum__identify_double_aa_mutations <- function(
  input_dt,
  singles_dt,
  wt_AAseq
  ){

  #WT AA sequences
  wt_AAseq_split <- strsplit(wt_AAseq,"")[[1]]

  ### Identify position and identity of double AA mutations
  ###########################

  #Double AA mutants
  doubles <- input_dt[Nham_aa==2]

  #If no double AA mutants, return empty data.table
  if(nrow(doubles)==0){
    return(data.table())
  }

  #Add position, mutant AA, WT AA and mean input count
  doubles[,Pos1 := which(strsplit(aa_seq,"")[[1]] !=wt_AAseq_split)[1],aa_seq]
  doubles[,Pos2 := which(strsplit(aa_seq,"")[[1]] !=wt_AAseq_split)[2],aa_seq]
  doubles[,Mut1 := strsplit(aa_seq,"")[[1]][Pos1],aa_seq]
  doubles[,Mut2 := strsplit(aa_seq,"")[[1]][Pos2],aa_seq]
  doubles[,WT_AA1 := wt_AAseq_split[Pos1],aa_seq]
  doubles[,WT_AA2 := wt_AAseq_split[Pos2],aa_seq]

  #Mean counts
  if(nrow(singles_dt)!=0){
    doubles <- merge(doubles,singles_dt[,.(Pos,Mut,s1_mean_count = mean_count)],by.x = c("Pos1","Mut1"),by.y = c("Pos","Mut"), all.x = T)
    doubles <- merge(doubles,singles_dt[,.(Pos,Mut,s2_mean_count = mean_count)],by.x = c("Pos2","Mut2"),by.y = c("Pos","Mut"), all.x = T)
  }else{
    doubles[, s1_mean_count := NA]
    doubles[, s2_mean_count := NA]
  }

  return(doubles)

}

# ---- dimsum__identify_double_nt_mutations
dimsum__identify_double_nt_mutations <- function(
  input_dt,
  singles_dt,
  wt_ntseq
  ){

  #WT nucleotide sequences
  wt_ntseq_split <- strsplit(wt_ntseq,"")[[1]]

  ### Identify position and identity of double nucleotide mutations
  ###########################

  #Double nucleotide mutants
  doubles <- input_dt[Nham_nt==2]

  #If no double nucleotide mutants, return empty data.table
  if(nrow(doubles)==0){
    return(data.table())
  }

  #Add positions, mutant nucleotides, WT nucleotides (WT_AA1, WT_AA2) and mean input count
  doubles[,Pos1 := which(strsplit(merge_seq,"")[[1]] !=wt_ntseq_split)[1],merge_seq]
  doubles[,Pos2 := which(strsplit(merge_seq,"")[[1]] !=wt_ntseq_split)[2],merge_seq]
  doubles[,Mut1 := strsplit(merge_seq,"")[[1]][Pos1],merge_seq]
  doubles[,Mut2 := strsplit(merge_seq,"")[[1]][Pos2],merge_seq]
  doubles[,WT_AA1 := wt_ntseq_split[Pos1],merge_seq]
  doubles[,WT_AA2 := wt_ntseq_split[Pos2],merge_seq]

  #Mean counts
  if(nrow(singles_dt)!=0){
    doubles <- merge(doubles,singles_dt[,.(Pos,Mut,s1_mean_count = mean_count)],by.x = c("Pos1","Mut1"),by.y = c("Pos","Mut"), all.x = T)
    doubles <- merge(doubles,singles_dt[,.(Pos,Mut,s2_mean_count = mean_count)],by.x = c("Pos2","Mut2"),by.y = c("Pos","Mut"), all.x = T)
  }else{
    doubles[, s1_mean_count := NA]
    doubles[, s2_mean_count := NA]
  }

  return(doubles)

}

# ---- dimsum__identify_permitted_mutations
dimsum__identify_permitted_mutations <- function(
  dimsum_meta,
  input_dt
  ){

  nuc_codes <- list(
    "A" = "A",
    "C" = "C",
    "G" = "G",
    "T" = "T",
    "R" = c("A", "G"),
    "Y" = c("C", "T"),
    "S" = c("C", "G"),
    "W" = c("A", "T"),
    "K" = c("G", "T"),
    "M" = c("A", "C"),
    "B" = c("C", "G", "T"),
    "D" = c("A", "G", "T"),
    "H" = c("A", "C", "T"),
    "V" = c("A", "C", "G"),
    "N" = c("A", "C", "G", "T")
    )

  #Loop over all nucleotide sequence positions and subset to permitted mutations
  input_dt[indel==F, permitted := T]
  for(i in 1:nchar(dimsum_meta[["permittedSequences"]])){
    input_dt[indel==F & !toupper(substr(nt_seq, i, i)) %in% nuc_codes[[substr(dimsum_meta[["permittedSequences"]], i, i)]], permitted := F]
  }

  return(input_dt)

}

# ---- dimsum__identify_single_aa_mutations
dimsum__identify_single_aa_mutations <- function(
  input_dt,
  wt_AAseq
  ){

  #WT AA sequences
  wt_AAseq_split <- strsplit(wt_AAseq,"")[[1]]

  ### Identify position and identity of single AA mutations
  ###########################

  #Single AA mutants
  singles <- input_dt[Nham_aa==1,]

  #If no single AA mutants, return empty data.table
  if(nrow(singles)==0){
    return(data.table())
  }

  #Add position, mutant AA, WT AA and mean input count
  singles[,Pos := which(strsplit(aa_seq,"")[[1]] !=wt_AAseq_split),aa_seq]
  singles[,Mut := strsplit(aa_seq,"")[[1]][Pos],aa_seq]
  singles[,WT_AA := wt_AAseq_split[Pos],aa_seq]

  #Remove unnecessary columns and rename fitness and sigma columns
  singles <- singles[,cbind(Pos,WT_AA,Mut,merge_seq,aa_seq,Nham_nt,Nham_aa,Nmut_codons,STOP,STOP_readthrough,error_model,mean_count,.SD),,.SDcols = grep("_uncorr", names(singles))]
  names(singles) <- gsub("_uncorr", "", names(singles))

  return(singles)

}

# ---- dimsum__identify_single_nt_mutations
dimsum__identify_single_nt_mutations <- function(
  input_dt,
  wt_ntseq
  ){

  #WT nucleotide sequences
  wt_ntseq_split <- strsplit(wt_ntseq,"")[[1]]

  ### Identify position and identity of single nucleotide mutations
  ###########################

  #Single nucleotide mutants
  singles_silent <- input_dt[Nham_nt==1,]

  #If no single nucleotide mutants, return empty data.table
  if(nrow(singles_silent)==0){
    return(data.table())
  }

  #Add position, mutant nucleotide, WT nucleotide (WT_AA)
  singles_silent[,Pos := which(strsplit(merge_seq,"")[[1]] !=wt_ntseq_split),merge_seq]
  singles_silent[,Mut := strsplit(merge_seq,"")[[1]][Pos],merge_seq]
  singles_silent[,WT_AA := wt_ntseq_split[Pos],merge_seq]

  #Remove unnecessary columns and rename fitness and sigma columns
  singles_silent <- singles_silent[,cbind(Pos,WT_AA,Mut,merge_seq,Nham_nt,Nham_aa,Nmut_codons,STOP,STOP_readthrough,error_model,mean_count,.SD),,.SDcols = grep("_uncorr", names(singles_silent))]
  names(singles_silent) <- gsub("_uncorr", "", names(singles_silent))

  return(singles_silent)

}

# ---- dimsum__infer_growth_rates
dimsum__infer_growth_rates <- function(
  dimsum_meta,
  input_dt,
  all_reps
  ){

  #If input data.table empty, return empty data.table
  if(nrow(input_dt)==0){
    return(data.table())
  }

  #Get selection_time data
  selection_time_unique <- unique(dimsum_meta[["exp_design"]][dimsum_meta[["exp_design"]][,"selection_id"]==1,c("sample_name", "experiment", "selection_time")])
  #Mean selection_time per biological replicate
  selection_time_mean <- tapply(selection_time_unique[,"selection_time"], paste0("e", selection_time_unique[,"experiment"]), mean)

  #Get cell_density data for output
  cell_density_output_unique <- unique(dimsum_meta[["exp_design"]][dimsum_meta[["exp_design"]][,"selection_id"]==1,c("sample_name", "experiment", "cell_density")])
  #Mean cell_density per output biological replicate
  cell_density_output_mean <- tapply(cell_density_output_unique[,"cell_density"], paste0("e", cell_density_output_unique[,"experiment"]), mean)

  #Get cell_density data for input
  cell_density_input_unique <- unique(dimsum_meta[["exp_design"]][dimsum_meta[["exp_design"]][,"selection_id"]==0,c("sample_name", "experiment", "cell_density")])
  #Mean cell_density per input biological replicate
  cell_density_input_mean <- tapply(cell_density_input_unique[,"cell_density"], paste0("e", cell_density_input_unique[,"experiment"]), mean)

  #Infer growth rates from fitness
  for (E in all_reps) {
    input_dt[error_model==T, paste0("growthrate", E) := log((.SD[[2]]/sum(.SD[[2]], na.rm = T)*cell_density_output_mean[paste0("e", E)]) / (.SD[[1]]/sum(.SD[[1]], na.rm = T)*cell_density_input_mean[paste0("e", E)])) / selection_time_mean[paste0("e", E)],,.SDcols = paste0("count_e", E, "_s", c(0, 1))]

    gr_cor <- input_dt[!is.na(get(paste0("growthrate", E))) & !is.infinite(get(paste0("growthrate", E))) & error_model==T,cor(.SD[[1]], .SD[[2]], use = "pairwise.complete"),,.SDcols = c(paste0("growthrate", E), paste0("fitness", E, "_uncorr"))]
    dimsum__status_message(paste0("Fitness vs. growth rate correlation for replicate ", E, ": ", round(gr_cor, 2), "\n"))

    gr_lm <- lm(y~x, input_dt[!is.na(get(paste0("growthrate", E))) & !is.infinite(get(paste0("growthrate", E))) & error_model==T,.(y = .SD[[1]], x = .SD[[2]]),,.SDcols = c(paste0("growthrate", E), paste0("fitness", E, "_uncorr"))])$coefficients
    input_dt[, paste0("growthrate", E) := .SD[[1]]*gr_lm[2]+gr_lm[1],,.SDcols = paste0("fitness", E, "_uncorr")]
    input_dt[, paste0("growthrate", E, "_sigma") := .SD[[1]]*gr_lm[2],,.SDcols = paste0("sigma", E, "_uncorr")]
  }

  #Merge replicate growth rates and errors
  all_reps_str <- paste0(all_reps, collapse="|")
  fitness_rx <- input_dt[,.SD,.SDcols = grep(paste0("growthrate(", all_reps_str, ")$"),colnames(input_dt))]
  sigma_rx <- input_dt[,.SD,.SDcols = grep(paste0("growthrate(", all_reps_str, ")_sigma"),colnames(input_dt))]
  input_dt[,growthrate := rowSums(fitness_rx/(sigma_rx^2),na.rm=T)/rowSums(1/(sigma_rx^2),na.rm=T)]
  gr_lm <- lm(growthrate~fitness, input_dt[error_model==T])$coefficients
  input_dt[, growthrate := .SD[[1]]*gr_lm[2]+gr_lm[1],,.SDcols = "fitness"]
  input_dt[, growthrate_sigma := .SD[[1]]*gr_lm[2],,.SDcols = "sigma"]

  return(input_dt)

}

# ---- dimsum__merge_fitness [modified, see the header and the generator]
dimsum__merge_fitness <- function(
  dimsum_meta,
  input_dt,
  doubles_dt,
  singles_dt,
  all_reps,
  fitness_outpath,
  report = TRUE,
  report_outpath = NULL
  ){

  dimsum__status_message("Merging fitness estimates from biological replicates...\n")

  aa_list <- DIMSUM_FITNESS_AMINO_ACID_CODE[strsplit("GAVLMIFYWKRHDESTCNQP", NULL)[[1]]]
  aa_list["*"] <- "*"

  #Number of input and output replicates
  all_reps_str <- paste0(all_reps, collapse="|")

  #### all variants
  fitness_rx <- input_dt[,.SD,.SDcols = grep(paste0("fitness(", all_reps_str, ")_uncorr"),colnames(input_dt))]
  sigma_rx <- input_dt[,.SD,.SDcols = grep(paste0("sigma(", all_reps_str, ")_uncorr"),colnames(input_dt))]
  input_dt[,fitness := rowSums(fitness_rx/(sigma_rx^2),na.rm=T)/rowSums(1/(sigma_rx^2),na.rm=T)]
  input_dt[,sigma := sqrt(1/rowSums(1/(sigma_rx^2),na.rm=T))]

  #### singles
  if(nrow(singles_dt)!=0){
    fitness_rx <- singles_dt[,.SD,.SDcols = grep(paste0("fitness(", all_reps_str, ")$"),colnames(singles_dt))]
    sigma_rx <- singles_dt[,.SD,.SDcols = grep(paste0("sigma(", all_reps_str, ")$"),colnames(singles_dt))]
    singles_dt[,fitness := rowSums(fitness_rx/(sigma_rx^2),na.rm=T)/rowSums(1/(sigma_rx^2),na.rm=T)]
    singles_dt[,sigma := sqrt(1/rowSums(1/(sigma_rx^2),na.rm=T))]
  }

  #### doubles
  #uncorrected fitness
  if(nrow(doubles_dt)!=0){
    fitness_rx <- doubles_dt[,.SD,.SDcols = grep(paste0("fitness(", all_reps_str, ")_uncorr"),colnames(doubles_dt))]
    sigma_rx <- doubles_dt[,.SD,.SDcols = grep(paste0("sigma(", all_reps_str, ")_uncorr"),colnames(doubles_dt))]
    doubles_dt[,fitness_uncorr := rowSums(fitness_rx/(sigma_rx^2),na.rm=T)/rowSums(1/(sigma_rx^2),na.rm=T)]
    doubles_dt[,sigma_uncorr := sqrt(1/rowSums(1/(sigma_rx^2),na.rm=T))]
  }

  #conditioned fitness
  if(dimsum_meta[["bayesianDoubleFitness"]] & nrow(doubles_dt)!=0){
    fitness_rx <- doubles_dt[,.SD,.SDcols = grep(paste0("fitness(", all_reps_str, ")_cond"),colnames(doubles_dt))]
    sigma_rx <- doubles_dt[,.SD,.SDcols = grep(paste0("sigma(", all_reps_str, ")_cond"),colnames(doubles_dt))]
    doubles_dt[,fitness_cond := rowSums(fitness_rx/(sigma_rx^2),na.rm=T)/rowSums(1/(sigma_rx^2),na.rm=T)]
    doubles_dt[,sigma_cond := sqrt(1/rowSums(1/(sigma_rx^2),na.rm=T))]
  }

  dimsum__status_message("Done\n")

  ### Add growth rates if cell_density and selection_time (in hours) supplied
  ###########################

  #Check that selection_time supplied for all output samples
  if(sum(is.na(dimsum_meta[["exp_design"]][dimsum_meta[["exp_design"]][,"selection_id"]==1,"selection_time"]))==0){
    #Check that cell_density supplied for all samples
    if(sum(is.na(dimsum_meta[["exp_design"]][,"cell_density"]))==0){

      dimsum__status_message("Inferring growth rates...\n")

      input_dt <- dimsum__infer_growth_rates(
        dimsum_meta = dimsum_meta,
        input_dt = input_dt,
        all_reps = all_reps
        )

      dimsum__status_message("Done\n")

    }
  }

  ### Output replicate data files
  ###########################

  dimsum__status_message("Saving fitness estimates...\n")

  #Rename columns
  names(input_dt)[names(input_dt)=="merge_seq"] <- "nt_seq"
  names(singles_dt)[names(singles_dt)=="merge_seq"] <- "nt_seq"
  names(doubles_dt)[names(doubles_dt)=="merge_seq"] <- "nt_seq"

  #Reformat columns
  if(dimsum_meta[["sequenceType"]]=="coding"){
    #Set nucleotide sequence to NA when fitness and error aggregated at AA level
    if(dimsum_meta[["mixedSubstitutions"]]){
      input_dt[error_model==T, nt_seq := NA]
      if(nrow(singles_dt)!=0){singles_dt[, nt_seq := NA]}
    }else{
      input_dt[(Nham_aa>0 | indel==T) & error_model==T, nt_seq := NA]
      if(nrow(singles_dt)!=0){singles_dt[Nham_aa>0, nt_seq := NA]}
    }
    if(nrow(doubles_dt)!=0){doubles_dt[, nt_seq := NA]}
    #Remove unnecessary columns
    unnecessary_cols <- c(
      "var_fitness",
      "avg_sigma",
      "isNA",
      "avg_sigma_fit",
      "merge_seq",
      "counts_for_bins",
      "bin_count")
    input_dt <- input_dt[,.SD,,.SDcols = names(input_dt)[!names(input_dt) %in% unnecessary_cols]]
    singles_dt <- singles_dt[,.SD,,.SDcols = names(singles_dt)[!names(singles_dt) %in% unnecessary_cols]]
    doubles_dt <- doubles_dt[,.SD,,.SDcols = names(doubles_dt)[!names(doubles_dt) %in% unnecessary_cols]]
  }else{
    #Remove unnecessary columns
    unnecessary_cols <- c(
      "aa_seq",
      "Nham_aa",
      "Nmut_codons",
      "STOP",
      "STOP_readthrough",
      "var_fitness",
      "avg_sigma",
      "isNA",
      "avg_sigma_fit",
      "merge_seq",
      "counts_for_bins",
      "bin_count")
    input_dt <- input_dt[,.SD,,.SDcols = names(input_dt)[!names(input_dt) %in% unnecessary_cols]]
    singles_dt <- singles_dt[,.SD,,.SDcols = names(singles_dt)[!names(singles_dt) %in% unnecessary_cols]]
    doubles_dt <- doubles_dt[,.SD,,.SDcols = names(doubles_dt)[!names(doubles_dt) %in% unnecessary_cols]]
    #Rename WT_AA columns to WT_nt
    names(input_dt)[grep("^WT_AA", names(input_dt))] <- gsub("WT_AA", "WT_nt", names(input_dt)[grep("^WT_AA", names(input_dt))])
    names(singles_dt)[grep("^WT_AA", names(singles_dt))] <- gsub("WT_AA", "WT_nt", names(singles_dt)[grep("^WT_AA", names(singles_dt))])
    names(doubles_dt)[grep("^WT_AA", names(doubles_dt))] <- gsub("WT_AA", "WT_nt", names(doubles_dt)[grep("^WT_AA", names(doubles_dt))])
  }

  #Separate synonymous variants
  synonymous <- input_dt[error_model==F,.SD,,.SDcols = names(input_dt)[names(input_dt)!="error_model"]]
  input_dt <- input_dt[error_model==T,.SD,,.SDcols = names(input_dt)[names(input_dt)!="error_model"]]

  #Rename objects
  all_variants <- input_dt
  wildtype <- all_variants[WT==T]
  doubles <- doubles_dt

  #Remove synonymous variants (in WT) from all_variants for coding sequences (with mixedSubstitutions==F, these are included in error modelling)
  if(dimsum_meta[["sequenceType"]]=="coding"){
    all_variants <- all_variants[WT==T | Nham_aa>0 | indel==T]
  }

  ### Output plain text files
  ###########################

  ##### finalize data.tables
  singles <- data.table()
  singles_mavedb <- data.table()
  if(nrow(singles_dt)!=0){
    if(dimsum_meta[["sequenceType"]]=="coding"){
      singles <- singles_dt[Nham_aa==1,.(Pos,WT_AA,Mut,nt_seq,aa_seq,Nham_nt,Nham_aa,Nmut_codons,STOP,STOP_readthrough,mean_count,fitness,sigma)]
      singles_mavedb <- singles[,.(hgvs_pro = NA, score = fitness, se = sigma)]
      singles_mavedb[, hgvs_pro := paste0("p.", aa_list[singles[,WT_AA]], singles[,Pos], aa_list[singles[,Mut]])]
    }else{
      singles <- singles_dt[,.(Pos,WT_nt,Mut,nt_seq,Nham_nt,mean_count,fitness,sigma)]
      singles_mavedb <- singles[,.(hgvs_nt = NA, score = fitness, se = sigma)]
      singles_mavedb[, hgvs_nt := paste0("n.", singles[,Pos], toupper(singles[,WT_nt]), ">", toupper(singles[,Mut]))]
    }
  }

  #for doubles #add single mutant fitness/sigma values to double mutant table
  if(nrow(doubles_dt)!=0 & nrow(singles)!=0){
    doubles_dt[,fitness1 := singles[Pos == Pos1 & Mut == Mut1,fitness],.(Pos1,Mut1)]
    doubles_dt[,sigma1 := singles[Pos == Pos1 & Mut == Mut1,sigma],.(Pos1,Mut1)]
    doubles_dt[,fitness2 := singles[Pos == Pos2 & Mut == Mut2,fitness],.(Pos2,Mut2)]
    doubles_dt[,sigma2 := singles[Pos == Pos2 & Mut == Mut2,sigma],.(Pos2,Mut2)]
  }

  if(dimsum_meta[["sequenceType"]]=="coding"){
    retained_cols <- c("Pos1","Pos2","WT_AA1","WT_AA2","Mut1","Mut2","nt_seq","aa_seq","Nham_nt","Nham_aa","Nmut_codons","STOP","STOP_readthrough","mean_count",
      "fitness1","sigma1","fitness2","sigma2",
      "fitness_uncorr","sigma_uncorr",
      "fitness_cond","sigma_cond")
    doubles <- doubles_dt[,.SD,,.SDcols = retained_cols[retained_cols %in% names(doubles_dt)]]

    #write data to files
    write.table(x = wildtype, file = file.path(fitness_outpath, "fitness_wildtype.txt"),
                quote = F,row.names = F, col.names = T)
    write.table(x = synonymous, file = file.path(fitness_outpath, "fitness_synonymous.txt"),
                quote = F,row.names = F, col.names = T)
    write.table(x = singles, file = file.path(fitness_outpath, "fitness_singles.txt"),
                quote = F,row.names = F, col.names = T)
    write.table(x = doubles, file = file.path(fitness_outpath, "fitness_doubles.txt"),
                quote = F,row.names = F, col.names = T)

  }else{
    retained_cols <- c("Pos1","Pos2","WT_nt1","WT_nt2","Mut1","Mut2","nt_seq","Nham_nt","mean_count",
      "fitness1","sigma1","fitness2","sigma2",
      "fitness_uncorr","sigma_uncorr",
      "fitness_cond","sigma_cond")
    doubles <- doubles_dt[,.SD,,.SDcols = retained_cols[retained_cols %in% names(doubles_dt)]]

    #write data to files
    write.table(x = wildtype, file = file.path(fitness_outpath, "fitness_wildtype.txt"),
                quote = F,row.names = F, col.names = T)
    write.table(x = singles, file = file.path(fitness_outpath, "fitness_singles.txt"),
                quote = F,row.names = F, col.names = T)
    write.table(x = doubles, file = file.path(fitness_outpath, "fitness_doubles.txt"),
                quote = F,row.names = F, col.names = T)
  }

  #Write MaveDB formatted singles
  write.table(x = singles_mavedb, file = file.path(fitness_outpath, "fitness_singles_MaveDB.csv"),
              quote = F,row.names = F, col.names = T, sep = ",")

  ### Output RData file
  ###########################

  if(dimsum_meta[["sequenceType"]]=="coding"){
    save(all_variants, wildtype, synonymous, singles, doubles,
      file = file.path(fitness_outpath, paste0(dimsum_meta[["projectName"]], '_fitness_replicates.RData')),
      version = 2)
  }else{
    save(all_variants, wildtype, singles, doubles,
      file = file.path(fitness_outpath, paste0(dimsum_meta[["projectName"]], '_fitness_replicates.RData')),
      version = 2)
  }

  dimsum__status_message("Done\n")

}

# ---- dimsum__merge_report
dimsum__merge_report <- function(
  dimsum_meta,
  report_outpath
  ){

  #Create report directory (if doesn't already exist)
  report_outpath <- gsub("/$", "", report_outpath)
  suppressWarnings(dir.create(report_outpath))

  #Add 'aligned_pair' column if user-specified count file supplied
  if(!is.null(dimsum_meta[["countPath"]])){
    dimsum_meta[['exp_design']][,'aligned_pair'] <- paste0(
      dimsum_meta[["exp_design"]][,"sample_name"], '_e',
      dimsum_meta[["exp_design"]][,"experiment"], '_s',
      dimsum_meta[["exp_design"]][,"selection_id"], '_b',
      dimsum_meta[["exp_design"]][,"biological_replicate"], '_tNA', sep = "")
  }

  #Variant processing statistics - nucleotide substitutions
  merge_df <- dimsum_meta[['exp_design']]
  merge_df[,'pairname'] <- sapply(strsplit(merge_df[,'aligned_pair'], '_t'), '[', 1)
  nuc_subst_df <- data.frame(
    'nuc_subst_0'=sapply(dimsum_meta[['nuc_subst_counts']], '[', '0'),
    'nuc_subst_1'=sapply(dimsum_meta[['nuc_subst_counts']], '[', '1'),
    'nuc_subst_2'=sapply(dimsum_meta[['nuc_subst_counts']], '[', '2'),
    'nuc_subst_sum'=sapply(dimsum_meta[['nuc_subst_counts']], sum),
    'nuc_indrt_sum'=sapply(dimsum_meta[['nuc_indrt_counts']], sum),
    'nuc_mxsub_sum'=sapply(dimsum_meta[['nuc_mxsub_counts']], sum),
    'nuc_tmsub_sum'=sapply(dimsum_meta[['nuc_tmsub_counts']], sum),
    'nuc_frbdn_sum'=sapply(dimsum_meta[['nuc_frbdn_counts']], sum),
    'nuc_const_sum'=sapply(dimsum_meta[['nuc_const_counts']], sum),
    'nuc_indel_sum'=sapply(dimsum_meta[['nuc_indel_counts']], sum),
    'nuc_nbarc_sum'=sapply(dimsum_meta[['nuc_nbarc_counts']], sum))
  nuc_subst_df[is.na(nuc_subst_df)] <- 0
  # nuc_subst_df[,'pairname'] <- dimsum__plot_samplename(unique(sapply(strsplit(merge_df[,'aligned_pair'], '_t'), '[', 1)))
  nuc_subst_df[,'pairname'] <- dimsum__plot_samplename(unique(sapply(strsplit(rownames(nuc_subst_df), '_t'), '[', 1)))
  nuc_subst_df_collapse <- plyr::ddply(nuc_subst_df, "pairname", plyr::summarise,
    nuc_subst_0 = sum(nuc_subst_0),
    nuc_subst_1 = sum(nuc_subst_1),
    nuc_subst_2 = sum(nuc_subst_2),
    nuc_subst_sum = sum(nuc_subst_sum),
    nuc_indrt_sum = sum(nuc_indrt_sum),
    nuc_mxsub_sum = sum(nuc_mxsub_sum),
    nuc_tmsub_sum = sum(nuc_tmsub_sum),
    nuc_frbdn_sum = sum(nuc_frbdn_sum),
    nuc_const_sum = sum(nuc_const_sum),
    nuc_indel_sum = sum(nuc_indel_sum),
    nuc_nbarc_sum = sum(nuc_nbarc_sum))
  nuc_subst_df_collapse[,'nuc_subst_3plus'] <- nuc_subst_df_collapse[,'nuc_subst_sum']-nuc_subst_df_collapse[,'nuc_subst_0']-nuc_subst_df_collapse[,'nuc_subst_1']-nuc_subst_df_collapse[,'nuc_subst_2']
  nuc_subst_df_collapse_perc <- nuc_subst_df_collapse
  nuc_subst_df_collapse_perc[,colnames(nuc_subst_df_collapse_perc)[-1]] <- nuc_subst_df_collapse_perc[,colnames(nuc_subst_df_collapse_perc)[-1]]/(
    nuc_subst_df_collapse_perc[,'nuc_subst_sum'] +
    nuc_subst_df_collapse_perc[,'nuc_indrt_sum'] +
    nuc_subst_df_collapse_perc[,'nuc_mxsub_sum'] +
    nuc_subst_df_collapse_perc[,'nuc_tmsub_sum'] +
    nuc_subst_df_collapse_perc[,'nuc_frbdn_sum'] +
    nuc_subst_df_collapse_perc[,'nuc_const_sum'] +
    nuc_subst_df_collapse_perc[,'nuc_indel_sum'] +
    nuc_subst_df_collapse_perc[,'nuc_nbarc_sum']) * 100
  nuc_subst_df_collapse_perc <- nuc_subst_df_collapse_perc[,-which(colnames(nuc_subst_df_collapse_perc) == "nuc_subst_sum")]
  temp_colnames <- c(
    "0 hamming dist.",
    "1 hamming dist.",
    "2 hamming dist.",
    "indel",
    "mixed",
    "too many",
    "not permitted",
    "internal constant region",
    "indel discarded",
    "invalid barcode",
    "3+ hamming dist.")
  colnames(nuc_subst_df_collapse_perc)[-1] <- temp_colnames
  nuc_subst_df_collapse <- nuc_subst_df_collapse[,-which(colnames(nuc_subst_df_collapse) == "nuc_subst_sum")]
  colnames(nuc_subst_df_collapse)[-1] <- temp_colnames

  #Plot 1: Nucleotide mutation counts
  plot_df <- reshape2::melt(nuc_subst_df_collapse, id="pairname")
  plot_df[,'Mutation_type'] <- factor(plot_df[,'variable'], levels=c(
    "0 hamming dist.",
    "1 hamming dist.",
    "2 hamming dist.",
    "3+ hamming dist.",
    "indel",
    "mixed",
    "too many",
    "not permitted",
    "internal constant region",
    "indel discarded",
    "invalid barcode"))
  d <- ggplot2::ggplot(plot_df, ggplot2::aes(pairname, value)) +
    ggplot2::geom_col(ggplot2::aes(fill = Mutation_type)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::labs(x = "Sample names", y = "Total reads with mutation")#, title = paste0("Read nucleotide mutation statistics"))
  dimsum__save_png(file.path(report_outpath, paste0('dimsum__merge_report_nucmutationcounts.png')), d, width=12, height=8)

  #Plot 2: Nucleotide mutation percentages
  plot_df <- reshape2::melt(nuc_subst_df_collapse_perc, id="pairname")
  plot_df[,'Mutation_type'] <- factor(plot_df[,'variable'], levels=c(
    "0 hamming dist.",
    "1 hamming dist.",
    "2 hamming dist.",
    "3+ hamming dist.",
    "indel",
    "mixed",
    "too many",
    "not permitted",
    "internal constant region",
    "indel discarded",
    "invalid barcode"))
  d <- ggplot2::ggplot(plot_df, ggplot2::aes(pairname, value)) +
    ggplot2::geom_col(ggplot2::aes(fill = Mutation_type)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::labs(x = "Sample names", y = "Percentage of reads with mutation")#, title = paste0("Read nucleotide mutation statistics (percentage)"))
  dimsum__save_png(file.path(report_outpath, paste0('dimsum__merge_report_nucmutationpercentages.png')), d, width=12, height=8)

  #Variant processing statistics - amino acid substitutions
  merge_df <- dimsum_meta[['exp_design']]
  merge_df[,'pairname'] <- sapply(strsplit(merge_df[,'aligned_pair'], '_t'), '[', 1)
  aa_subst_df <- data.frame(
    'aa_subst_0'=sapply(dimsum_meta[['aa_subst_counts']], '[', '0'),
    'aa_subst_1'=sapply(dimsum_meta[['aa_subst_counts']], '[', '1'),
    'aa_subst_2'=sapply(dimsum_meta[['aa_subst_counts']], '[', '2'),
    'aa_subst_sum'=sapply(dimsum_meta[['aa_subst_counts']], sum),
    'nuc_indrt_sum'=sapply(dimsum_meta[['nuc_indrt_counts']], sum),
    'nuc_mxsub_sum'=sapply(dimsum_meta[['nuc_mxsub_counts']], sum),
    'nuc_tmsub_sum'=sapply(dimsum_meta[['nuc_tmsub_counts']], sum),
    'nuc_frbdn_sum'=sapply(dimsum_meta[['nuc_frbdn_counts']], sum),
    'nuc_const_sum'=sapply(dimsum_meta[['nuc_const_counts']], sum),
    'nuc_indel_sum'=sapply(dimsum_meta[['nuc_indel_counts']], sum),
    'nuc_nbarc_sum'=sapply(dimsum_meta[['nuc_nbarc_counts']], sum))
  aa_subst_df[is.na(aa_subst_df)] <- 0
  # aa_subst_df[,'pairname'] <- dimsum__plot_samplename(unique(sapply(strsplit(merge_df[,'aligned_pair'], '_t'), '[', 1)))
  aa_subst_df[,'pairname'] <- dimsum__plot_samplename(unique(sapply(strsplit(rownames(aa_subst_df), '_t'), '[', 1)))
  aa_subst_df_collapse <- plyr::ddply(aa_subst_df, "pairname", plyr::summarise,
    aa_subst_0 = sum(aa_subst_0),
    aa_subst_1 = sum(aa_subst_1),
    aa_subst_2 = sum(aa_subst_2),
    aa_subst_sum = sum(aa_subst_sum),
    nuc_indrt_sum = sum(nuc_indrt_sum),
    nuc_mxsub_sum = sum(nuc_mxsub_sum),
    nuc_tmsub_sum = sum(nuc_tmsub_sum),
    nuc_frbdn_sum = sum(nuc_frbdn_sum),
    nuc_const_sum = sum(nuc_const_sum),
    nuc_indel_sum = sum(nuc_indel_sum),
    nuc_nbarc_sum = sum(nuc_nbarc_sum))
  aa_subst_df_collapse[,'aa_subst_3plus'] <- aa_subst_df_collapse[,'aa_subst_sum']-aa_subst_df_collapse[,'aa_subst_0']-aa_subst_df_collapse[,'aa_subst_1']-aa_subst_df_collapse[,'aa_subst_2']
  aa_subst_df_collapse_perc <- aa_subst_df_collapse
  aa_subst_df_collapse_perc[,colnames(aa_subst_df_collapse_perc)[-1]] <- aa_subst_df_collapse_perc[,colnames(aa_subst_df_collapse_perc)[-1]]/(
    aa_subst_df_collapse_perc[,'aa_subst_sum'] +
    aa_subst_df_collapse_perc[,'nuc_indrt_sum'] +
    aa_subst_df_collapse_perc[,'nuc_mxsub_sum'] +
    aa_subst_df_collapse_perc[,'nuc_tmsub_sum'] +
    aa_subst_df_collapse_perc[,'nuc_frbdn_sum'] +
    aa_subst_df_collapse_perc[,'nuc_const_sum'] +
    aa_subst_df_collapse_perc[,'nuc_indel_sum'] +
    aa_subst_df_collapse_perc[,'nuc_nbarc_sum']) * 100
  aa_subst_df_collapse_perc <- aa_subst_df_collapse_perc[,-which(colnames(aa_subst_df_collapse_perc) == "aa_subst_sum")]
  temp_colnames <- c(
    "0 hamming dist.",
    "1 hamming dist.",
    "2 hamming dist.",
    "indel",
    "mixed",
    "too many",
    "not permitted",
    "internal constant region",
    "indel discarded",
    "invalid barcode",
    "3+ hamming dist.")
  colnames(aa_subst_df_collapse_perc)[-1] <- temp_colnames
  aa_subst_df_collapse <- aa_subst_df_collapse[,-which(colnames(aa_subst_df_collapse) == "aa_subst_sum")]
  colnames(aa_subst_df_collapse)[-1] <- temp_colnames

  #Plot 3: Nucleotide mutation counts
  plot_df <- reshape2::melt(aa_subst_df_collapse, id="pairname")
  plot_df[,'Mutation_type'] <- factor(plot_df[,'variable'], levels=c(
    "0 hamming dist.",
    "1 hamming dist.",
    "2 hamming dist.",
    "3+ hamming dist.",
    "indel",
    "mixed",
    "too many",
    "not permitted",
    "internal constant region",
    "indel discarded",
    "invalid barcode"))
  d <- ggplot2::ggplot(plot_df, ggplot2::aes(pairname, value)) +
    ggplot2::geom_col(ggplot2::aes(fill = Mutation_type)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::labs(x = "Sample names", y = "Total reads with mutation")#, title = paste0("Read amino acid mutation statistics"))
  dimsum__save_png(file.path(report_outpath, paste0('dimsum__merge_report_aamutationcounts.png')), d, width=12, height=8)

  #Plot 4: Nucleotide mutation percentages
  plot_df <- reshape2::melt(aa_subst_df_collapse_perc, id="pairname")
  plot_df[,'Mutation_type'] <- factor(plot_df[,'variable'], levels=c(
    "0 hamming dist.",
    "1 hamming dist.",
    "2 hamming dist.",
    "3+ hamming dist.",
    "indel",
    "mixed",
    "too many",
    "not permitted",
    "internal constant region",
    "indel discarded",
    "invalid barcode"))
  d <- ggplot2::ggplot(plot_df, ggplot2::aes(pairname, value)) +
    ggplot2::geom_col(ggplot2::aes(fill = Mutation_type)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::labs(x = "Sample names", y = "Percentage of reads with mutation")#, title = paste0("Read amino acid mutation statistics (percentage)"))
  dimsum__save_png(file.path(report_outpath, paste0('dimsum__merge_report_aamutationpercentages.png')), d, width=12, height=8)

  #Render report
  dimsum__render_report(dimsum_meta = dimsum_meta)

  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  return(dimsum_meta_new)
}

# ---- dimsum__normalise_fitness
dimsum__normalise_fitness <- function(
  dimsum_meta,
  input_dt,
  all_reps,
  fitness_suffix=""
  ){

  #If input data.table empty, return empty data.table
  if(nrow(input_dt)==0){
    return(data.table())
  }

  #Get generations data
  generations_unique <- unique(dimsum_meta[["exp_design"]][dimsum_meta[["exp_design"]][,"selection_id"]==1,c("sample_name", "experiment", "generations")])
  #Mean generations per biological replicate
  generations_mean <- tapply(generations_unique[,"generations"], paste0("e", generations_unique[,"experiment"]), mean)

  #Normalise fitness and sigma by mean generations
  for (E in all_reps) {
    for(i in c("fitness", "sigma")){
      input_dt[,paste0(i, E, fitness_suffix) := log2(exp(.SD/generations_mean[paste0("e", E)])),,.SDcols = paste0(i, E, fitness_suffix)]
    }
  }

  return(input_dt)

}

# ---- dimsum__parse_minimum_read_count_arguments
dimsum__parse_minimum_read_count_arguments <- function(
  input_arg
  ){
  #A single integer supplied
  if(!is.na(as.integer(input_arg))){
    return(as.integer(input_arg))
  }

  #Multiple integers supplied
  temp_list <- sapply(as.list(unlist(strsplit(input_arg, ","))), strsplit, ":")
  names(temp_list) <- sapply(temp_list, '[', 1)
  #List of list of 2 integers
  if(sum(is.na(as.integer(unlist(temp_list))))!=0 | sum(sapply(temp_list, length)!=2)!=0 | sum(duplicated(names(temp_list)))!=0){
    stop("Invalid 'fitness...Count...' arguments. Multiple thresholds should be specified in the form 'a:b,c:d' where a,b,c and d are positive integers.", call. = FALSE)
  }
  #Set names
  names(temp_list) <- sapply(temp_list, '[', 1)
  #Return list
  return(lapply(lapply(temp_list, '[', 2), as.integer))
}

# ---- dimsum__plot_samplename
dimsum__plot_samplename <- function(
  input_names){
  if(sum(grepl("_t", input_names)) != 0){
    output_names <- sapply(lapply(lapply(strsplit(input_names, "_"), unlist), '[', c(1,5)), paste, collapse = "_")
    output_names <- gsub("_tNA", "", output_names)
  }else{
    output_names <- sapply(lapply(strsplit(input_names, "_"), unlist), '[', 1)
  }
  return(output_names)
}

# ---- dimsum__process_merged_variants [modified, see the header and the generator]
dimsum__process_merged_variants <- function(
  dimsum_meta,
  input_dt
  ){

  dimsum__status_message("Processing merged variants...\n")

  variant_dt <- input_dt

  #WT nucleotide sequence
  wt_ntseq <- tolower(dimsum_meta[['wildtypeSequence']])
  #WT AA sequence
  wt_AAseq <- paste0(seqinr::translate(strsplit(wt_ntseq,split="")[[1]]),collapse="")

  #Nucleotide mutation count dictionaries
  nuc_nbarc_dict <- list()
  nuc_indel_dict <- list()
  nuc_const_dict <- list()
  nuc_frbdn_dict <- list()
  nuc_tmsub_dict <- list()
  nuc_mxsub_dict <- list()
  nuc_indrt_dict <- list()
  nuc_subst_dict <- list()
  aa_subst_dict <- list()

  ### Reverse complement (if necessary)
  ###########################

  if(dimsum_meta[["reverseComplement"]]){
    variant_dt[,nt_seq := tolower(dimsum_fitness__reverse_complement(nt_seq))]
  }

  ### Debarcode variants (if barcode identity file supplied)
  ###########################

  variant_dt[, barcode_valid := TRUE]

  if(!is.null(dimsum_meta[["barcodeIdentityPath"]])){
    variant_dt <- dimsum__debarcode_variants(
      dimsum_meta = dimsum_meta,
      input_dt = variant_dt)
  }

  #No barcode variant count statistics
  nuc_nbarc_dict <- as.list(apply(
    variant_dt[barcode_valid==F,.SD,,.SDcols = names(variant_dt)[grepl("_count$", names(variant_dt))]],
    2, sum))

  #No barcode variants
  nobarcode_variant_dt <- data.table::copy(variant_dt[barcode_valid==F,.SD,,.SDcols = c(
    "nt_seq",
    names(variant_dt)[grepl("_count$", names(variant_dt))],
    "barcode_valid")])

  ### Update remaining variants
  ###########################

  #Remaining variants with valid barcodes
  variant_dt <- variant_dt[barcode_valid==T,.SD,,.SDcols = c(
    "nt_seq",
    names(variant_dt)[grepl("_count$", names(variant_dt))],
    "barcode_valid")]

  #Check if variants remaining
  if(variant_dt[,.N]==0){
    stop(paste0("Cannot proceed with variant processing: No valid variant barcodes found"), call. = FALSE)
  }

  #Define indel variants as those with length different from WT
  variant_dt[, indel := F]
  variant_dt[nchar(nt_seq)!=nchar(wt_ntseq), indel := T]

  #Check if substitution variants exist (if indels discarded)
  if(variant_dt[indel==F,.N]==0 & !dimsum_meta[["indels"]]){
    stop(paste0("Cannot proceed with variant processing: No substitution variants found"), call. = FALSE)
  }

  #Calculate nucleotide hamming distance
  variant_dt[indel==F, Nham_nt := dimsum_fitness__hamming_distances(nt_seq, wt_ntseq)]

  #Determine amino acid sequence
  variant_dt[, aa_seq := dimsum_fitness__translate(nt_seq)]
  #Calculate amino acid hamming distance
  variant_dt[indel==F, Nham_aa := dimsum_fitness__hamming_distances(aa_seq, wt_AAseq)]
  #Reorder
  variant_dt <- variant_dt[,.SD,,.SDcols = c(
    "nt_seq", "aa_seq", "Nham_nt", "Nham_aa",
    names(variant_dt)[grepl("_count$", names(variant_dt))],
    "barcode_valid", "indel")]

  ### Indel variants
  ###########################

  if(!dimsum_meta[["indels"]]){
    variant_dt[indel==T, indel_discarded := T]
  }else{
    if(is.na(dimsum_meta[["indelLengths"]][1])){
      variant_dt[indel==T, indel_discarded := F]
    }else{
      variant_dt[indel==T, indel_discarded := !nchar(nt_seq) %in% dimsum_meta[["indelLengths"]]]
    }
  }

  #Dicarded indel variant count statistics
  nuc_indel_dict <- as.list(apply(
    variant_dt[indel_discarded==T,.SD,,.SDcols = names(variant_dt)[grepl("_count$", names(variant_dt))]],
    2, sum))

  #Retained indel variant count statistics
  nuc_indrt_dict <- as.list(apply(
    variant_dt[indel_discarded==F,.SD,,.SDcols = names(variant_dt)[grepl("_count$", names(variant_dt))]],
    2, sum))

  #Indel variants
  indel_variant_dt <- data.table::copy(variant_dt[indel_discarded==T,.SD,,.SDcols = c(
    "nt_seq", "aa_seq",
    names(variant_dt)[grepl("_count$", names(variant_dt))],
    "barcode_valid", "indel_discarded")])

  ### Update remaining variants
  ###########################

  #Remaining variants without discarded indels
  variant_dt <- variant_dt[indel==F | indel_discarded==F,.SD,,.SDcols = c(
    "nt_seq", "aa_seq", "Nham_nt", "Nham_aa",
    names(variant_dt)[grepl("_count$", names(variant_dt))],
    "barcode_valid", "indel")]
  #Indicate WT sequence
  variant_dt[nt_seq == wt_ntseq,WT := TRUE]

  #Check if variants remaining
  if(variant_dt[,.N]==0){
    stop(paste0("Cannot proceed with variant processing: No variants found after indel variant filtering"), call. = FALSE)
  }

  #Check if WT sequence exists
  if(variant_dt[WT==T,.N]==0){
    variant_dt[, all_reads := rowSums(.SD > 0) == length(grep("_count$", names(variant_dt))),,.SDcols = grep("_count$", names(variant_dt))]
    variant_dt[, mean_count := rowMeans(.SD),,.SDcols = grep(paste0("_s[0].*_count$"), names(variant_dt))]
    dimsum__status_message(paste0("WT variant not found. Did you mean to specify one of the following?\n"))
    print(variant_dt[all_reads == T,][order(mean_count, decreasing = T)[1:5],.(nt_seq = toupper(nt_seq), all_reads, mean_count)])
    stop(paste0("Cannot proceed with variant processing: WT variant not found"), call. = FALSE)
  }

  #Indicate STOPs
  variant_dt <- dimsum__identify_STOP_mutations(variant_dt)

  ### Rejected variants (mutated constant region)
  ###########################

  if(sum(!strsplit(dimsum_meta[["wildtypeSequenceCoded"]], "")[[1]] %in% c("A", "C", "G", "T"))!=0){
    variant_dt <- dimsum__remove_internal_constant_region(
      input_dt = variant_dt,
      wt_ccntseq = dimsum_meta[["wildtypeSequenceCoded"]])
  }else{
    variant_dt[, constant_region := T]
  }

  #Mutated constant region variant count statistics
  nuc_const_dict <- as.list(apply(
    variant_dt[constant_region==F,.SD,,.SDcols = names(variant_dt)[grepl("_count$", names(variant_dt))]],
    2, sum))

  #Constant region variants
  rejected_variant_dt <- data.table::copy(variant_dt[constant_region==F,.SD,,.SDcols = c(
    "nt_seq", "aa_seq",
    names(variant_dt)[grepl("_count$", names(variant_dt))],
    "barcode_valid", "indel", "constant_region")])

  ### Update remaining variants
  ###########################

  #Remaining variants without mutated constant region (or with indels if retained)
  variant_dt <- variant_dt[constant_region==T | (indel==T & dimsum_meta[["indels"]]),]
  #Indicate STOPs
  variant_dt <- dimsum__identify_STOP_mutations(variant_dt)
  #Recalculate nucleotide hamming distance
  variant_dt[indel==F, Nham_nt := dimsum_fitness__hamming_distances(nt_seq, variant_dt[WT==T,nt_seq])]
  #Recalculate amino acid hamming distance
  variant_dt[indel==F, Nham_aa := dimsum_fitness__hamming_distances(aa_seq, variant_dt[WT==T,aa_seq])]

  #Check if variants remaining
  if(variant_dt[,.N]==0){
    stop(paste0("Cannot proceed with variant processing: No variants found after mutated constant region filtering"), call. = FALSE)
  }

  ### Rejected variants (forbidden mutations)
  ###########################

  variant_dt <- dimsum__identify_permitted_mutations(
    dimsum_meta = dimsum_meta,
    input_dt = variant_dt)

  #Forbidden variant count statistics
  nuc_frbdn_dict <- as.list(apply(
    variant_dt[permitted==F,.SD,,.SDcols = names(variant_dt)[grepl("_count$", names(variant_dt))]],
    2, sum))

  #Forbidden variants
  rejected_variant_dt <- rbind(rejected_variant_dt,
    data.table::copy(variant_dt[permitted==F,.SD,,.SDcols = c(
    "nt_seq", "aa_seq",
    names(variant_dt)[grepl("_count$", names(variant_dt))],
    "barcode_valid", "indel", "constant_region", "permitted")]), fill = T)

  ### Update remaining variants
  ###########################

  #Remaining variants with permitted mutations (or with indels if retained)
  variant_dt <- variant_dt[permitted==T | (indel==T & dimsum_meta[["indels"]]),]

  #Check if variants remaining
  if(variant_dt[,.N]==0){
    stop(paste0("Cannot proceed with variant processing: No variants found after forbidden mutation filtering"), call. = FALSE)
  }

  ### Rejected variants (too many substitutions)
  ###########################

  variant_dt[indel==F, too_many_substitutions := F]
  if(dimsum_meta[["sequenceType"]]=="coding"){
    variant_dt[indel==F & Nham_aa>dimsum_meta[["maxSubstitutions"]], too_many_substitutions := T]
  }else{
    variant_dt[indel==F & Nham_nt>dimsum_meta[["maxSubstitutions"]], too_many_substitutions := T]
  }

  #Too many substitution variant count statistics
  nuc_tmsub_dict <- as.list(apply(
    variant_dt[too_many_substitutions==T,.SD,,.SDcols = names(variant_dt)[grepl("_count$", names(variant_dt))]],
    2, sum))

  #Too many substitution variants
  rejected_variant_dt <- rbind(rejected_variant_dt,
    data.table::copy(variant_dt[too_many_substitutions==T,.SD,,.SDcols = c(
    "nt_seq", "aa_seq",
    names(variant_dt)[grepl("_count$", names(variant_dt))],
    "barcode_valid", "indel", "constant_region", "permitted", "too_many_substitutions")]), fill = T)

  ### Update remaining variants
  ###########################

  #Remaining variants without too many substitutions (or with indels if retained)
  variant_dt <- variant_dt[too_many_substitutions==F | (indel==T & dimsum_meta[["indels"]]),]

  #Check if variants remaining
  if(variant_dt[,.N]==0){
    stop(paste0("Cannot proceed with variant processing: No variants found after too many substitutions filtering"), call. = FALSE)
  }

  ### Rejected variants (mixed substitutions)
  ###########################

  #Check if WT sequence exists
  if(variant_dt[WT==T,.N]==0){
    variant_dt[, all_reads := rowSums(.SD > 0) == length(grep("_count$", names(variant_dt))),,.SDcols = grep("_count$", names(variant_dt))]
    variant_dt[, mean_count := rowMeans(.SD),,.SDcols = grep(paste0("_s[0].*_count$"), names(variant_dt))]
    dimsum__status_message(paste0("WT variant not found. Did you mean to specify one of the following?\n"))
    print(variant_dt[all_reads == T,][order(mean_count, decreasing = T)[1:5],.(nt_seq = toupper(nt_seq), all_reads, mean_count)])
    stop(paste0("Cannot proceed with variant processing: WT variant not found"), call. = FALSE)
  }

  #Add number of codons affected by mutations
  wt_ntseq <- variant_dt[WT==T,nt_seq]
  wt_ntseq_split <- strsplit(wt_ntseq,"")[[1]]
  variant_dt[indel==F, Nmut_codons := length(unique(ceiling(which(strsplit(nt_seq,"")[[1]] != wt_ntseq_split)/3))),nt_seq]

  #Identify mixed substitutions
  variant_dt[indel==F, mixed_substitutions := FALSE]
  if(dimsum_meta[["sequenceType"]]=="coding" & !dimsum_meta[["mixedSubstitutions"]]){
    variant_dt[indel==F & (Nmut_codons-Nham_aa) != 0 & Nham_aa != 0, mixed_substitutions := TRUE]
  }

  #Mixed substitution variant count statistics
  nuc_mxsub_dict <- as.list(apply(
    variant_dt[mixed_substitutions==T,.SD,,.SDcols = names(variant_dt)[grepl("_count$", names(variant_dt))]],
    2, sum))

  #Mixed substitution variants
  rejected_variant_dt <- rbind(rejected_variant_dt,
    data.table::copy(variant_dt[mixed_substitutions==T,.SD,,.SDcols = c(
    "nt_seq", "aa_seq",
    names(variant_dt)[grepl("_count$", names(variant_dt))],
    "barcode_valid", "indel", "constant_region", "permitted", "too_many_substitutions", "mixed_substitutions")]), fill = T)

  ### Remaining variants
  ###########################

  #Remaining variants without mixed substitutions (or with indels if retained)
  variant_dt <- variant_dt[mixed_substitutions==F | (indel==T & dimsum_meta[["indels"]]),.SD,,.SDcols = c(
    "nt_seq", "aa_seq", "WT", "STOP", "STOP_readthrough", "Nham_nt", "Nham_aa", "Nmut_codons",
    names(variant_dt)[grepl("_count$", names(variant_dt))],
    "barcode_valid", "indel", "constant_region", "permitted", "too_many_substitutions", "mixed_substitutions")]

  #Check if variants remaining
  if(variant_dt[,.N]==0){
    stop(paste0("Cannot proceed with variant processing: No variants found after mixed substitutions filtering"), call. = FALSE)
  }

  #Filtered variant nucleotide distributions
  for(count_column in names(variant_dt)[grepl("_count$", names(variant_dt))]){
    nuc_subst_dict[[count_column]] <- tapply(
      unlist(variant_dt[,.SD,,.SDcols = count_column]),
      variant_dt[, Nham_nt],
      sum, na.rm = T)
  }

  #Filtered variant amino acid distributions
  for(count_column in names(variant_dt)[grepl("_count$", names(variant_dt))]){
    aa_subst_dict[[count_column]] <- tapply(
      unlist(variant_dt[,.SD,,.SDcols = count_column]),
      variant_dt[, Nham_aa],
      sum, na.rm = T)
  }

  ### Save variant statistics
  ###########################

  #Save mutation statistics dictionaries
  mutation_stats_dicts <- list(
    "aa_subst_dict" = aa_subst_dict,
    "nuc_subst_dict" = nuc_subst_dict,
    "nuc_indrt_dict" = nuc_indrt_dict,
    "nuc_mxsub_dict" = nuc_mxsub_dict,
    "nuc_tmsub_dict" = nuc_tmsub_dict,
    "nuc_frbdn_dict" = nuc_frbdn_dict,
    "nuc_const_dict" = nuc_const_dict,
    "nuc_indel_dict" = nuc_indel_dict,
    "nuc_nbarc_dict" = nuc_nbarc_dict)
  save(mutation_stats_dicts,
    file = file.path(dimsum_meta[["tmp_path"]], "mutation_stats_dicts.RData"),
    version = 2)

  dimsum__status_message("Done\n")

  return(list(
    "nobarcode_variants" = nobarcode_variant_dt,
    "indel_variants" = indel_variant_dt,
    "rejected_variants" = rejected_variant_dt,
    "retained_variants" = variant_dt))

}

# ---- dimsum__remove_internal_constant_region [modified, see the header and the generator]
dimsum__remove_internal_constant_region <- function(
  input_dt,
  wt_ccntseq
  ){

  #Separate nucleotide sequences into pieces with and without constant region
  input_dt[indel==F, nt_seq_cnst := ""]
  input_dt[indel==F, nt_seq_wcnst := ""]
  idx <- 0
  for(b in strsplit(wt_ccntseq, "")[[1]]){
    idx <- idx + 1
    if(b %in% c("A", "C", "G", "T")){
      input_dt[indel==F, nt_seq_wcnst := paste0(nt_seq_wcnst, substr(nt_seq, idx, idx))]
    }else{
      input_dt[indel==F, nt_seq_cnst := paste0(nt_seq_cnst, substr(nt_seq, idx, idx))]
    }
  }

  #Indicate variants with WT internal constant region
  input_dt[indel==F, constant_region := F]
  input_dt[indel==F & nt_seq_cnst==input_dt[WT==T,nt_seq_cnst], constant_region := T]

  #Update nucleotide and amino acid sequences
  input_dt[indel==F & constant_region==T, nt_seq := nt_seq_wcnst]
  input_dt[indel==F & constant_region==T, aa_seq := dimsum_fitness__translate(nt_seq)]

  #Remove unnecessary columns
  output_dt <- input_dt[,.SD,,.SDcols = names(input_dt)[!names(input_dt) %in% c("nt_seq_cnst", "nt_seq_wcnst")]]

  return(output_dt)

}

# ---- dimsum__render_report [modified, see the header and the generator]
dimsum__render_report <- function(
  dimsum_meta,
  initialise=FALSE,
  finalise=FALSE
  ){

  #Initialise plots
  if(initialise){
    all_files <- list(
      "dimsum__fastqc_report_pair1_fastqc.png",
      "dimsum__fastqc_report_pair2_fastqc.png",
      "dimsum__cutadapt_report_pair1.png",
      "dimsum__cutadapt_report_pair2.png",
      "dimsum__vsearch_report_paircounts.png",
      "dimsum__vsearch_report_mergedlength.png",
      "dimsum__merge_report_nucmutationpercentages.png",
      "dimsum__merge_report_nucmutationcounts.png",
      "dimsum__merge_report_aamutationpercentages.png",
      "dimsum__merge_report_aamutationcounts.png",
      "dimsum__diagnostics_report_count_hist_input_nt.png",
      "dimsum__diagnostics_report_count_hist_input_aa.png",
      "dimsum__diagnostics_report_scatterplotmatrix_all.png",
      "dimsum_stage_fitness_report_1_errormodel_fitness_inputcounts.png",
      "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_density.png",
      "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_density_norm.png",
      "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_scatter.png",
      "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_scatter_norm.png",
      "dimsum_stage_fitness_report_1_errormodel_repspec.png",
      "dimsum_stage_fitness_report_1_errormodel_leaveoneout_qqplot.png")
    lapply(all_files, function(x){write(NULL, file = file.path(dimsum_meta[["project_path"]], "reports", x))})
  }

  #Set end time
  if(finalise){
    #Set run end time
    dimsum_meta[["end_time"]] <- Sys.time()
  }

  #Document settings
  doc_settings <- list()
  doc_settings[["dimsum_version"]] <- DIMSUM_FITNESS_VERSION
  doc_settings[["project_name"]] <- dimsum_meta[["projectName"]]
  doc_settings[["start_time"]] <- dimsum_meta[["start_time"]]
  doc_settings[["end_time"]] <- dimsum_meta[["end_time"]]
  doc_settings[["arg_list"]] <- dimsum_meta[["arg_list"]]
  doc_settings[["show_qualitycontrol1"]] <- TRUE
  doc_settings[["show_qualitycontrol2"]] <- TRUE
  doc_settings[["show_trim1"]] <- TRUE
  doc_settings[["show_trim2"]] <- TRUE
  doc_settings[["show_align1"]] <- TRUE
  doc_settings[["show_process1"]] <- TRUE
  doc_settings[["show_processA"]] <- TRUE
  doc_settings[["show_processC"]] <- TRUE
  doc_settings[["show_analyse1"]] <- TRUE
  doc_settings[["show_analyseN"]] <- TRUE

  #Remove unnecessary items
  #Noncoding sequence
  if(dimsum_meta[["sequenceType"]]!="coding"){
    doc_settings[["show_processA"]] <- FALSE
    doc_settings[["show_processC"]] <- FALSE
  }
  #Random mutagenesis
  if(dimsum_meta[["mutagenesisType"]]=="random"){
    doc_settings[["show_processC"]] <- FALSE
  }
  #Single-end sequencing
  if(!dimsum_meta[["paired"]]){
    doc_settings[["show_qualitycontrol2"]] <- FALSE
    doc_settings[["show_trim2"]] <- FALSE
  }
  #No fitness normalisation
  if(!dimsum_meta[["fitnessNormalise"]]){
    doc_settings[["show_analyseN"]] <- FALSE
  }
  #No error model or No biological replicates (requirement for normalisation and error model fit)
  if(length(unique(dimsum_meta[["exp_design"]][,"experiment_replicate"]))<2 | !dimsum_meta[["fitnessErrorModel"]]){
    doc_settings[["show_analyse1"]] <- FALSE
    doc_settings[["show_analyseN"]] <- FALSE
  }
  #STEAM reports only
  if(!is.null(dimsum_meta[["countPath"]])){
    doc_settings[["show_qualitycontrol1"]] <- FALSE
    doc_settings[["show_qualitycontrol2"]] <- FALSE
    doc_settings[["show_trim1"]] <- FALSE
    doc_settings[["show_trim2"]] <- FALSE
    doc_settings[["show_align1"]] <- FALSE
  }

  #Save document settings
  save(doc_settings, file = file.path(dimsum_meta[["project_path"]], "reports", "report_settings.RData"))

  #Copy report R markdown
  file.copy(
    from = DIMSUM_FITNESS_REPORT_RMD,
    to = file.path(dimsum_meta[["project_path"]], "reports", "report.Rmd"))

  #Copy DiMSum image
  file.copy(
    from = DIMSUM_FITNESS_LOGO_PNG,
    to = file.path(dimsum_meta[["project_path"]], "reports", "Dumpling.png"))

  #Render R markdown
  markdown_file <- file.path(dimsum_meta[["project_path"]], "reports", "report.Rmd")
  { sink("/dev/null")
    suppressMessages(rmarkdown::render(
      input = markdown_file,
      output_file = "report.html",
      output_dir = dimsum_meta[["project_path"]]))
    sink()
  }

}

# ---- dimsum__replicate_fitness_deviation
dimsum__replicate_fitness_deviation <- function(
  p,
  fitness_mat,
  all_reps
  ){
  #Normalised fitness (scale and center/shift)
  F_norm <- (fitness_mat +
    matrix(p[(length(all_reps) + 1):(2*length(all_reps))],
      nrow = nrow(fitness_mat), ncol = length(all_reps), byrow = T)) *
    matrix(p[1:length(all_reps)],
      nrow = nrow(fitness_mat), ncol = length(all_reps), byrow = T)
  #Average fitness (center/shift only)
  F_avg <- rowMeans(fitness_mat + matrix(p[(length(all_reps) + 1):(2*length(all_reps))],
    nrow = nrow(fitness_mat), ncol = length(all_reps), byrow = T))
  diffF <- sqrt(rowSums((F_norm - F_avg)^2))
  return(sum(diffF))
}

# ---- dimsum__sample_count_distributions [modified, see the header and the generator]
dimsum__sample_count_distributions <- function(
  input_dt,
  output_file_prefix,
  width = 12,
  height = 8,
  title = "",
  error_rate = 1e-3,
  seq_length){

  #Check if something to plot
  if(input_dt[is.na(WT),.N]==0){
    warning("dimsum__sample_count_distributions.R: No data to plot (empty data.table 'input_dt').", call. = FALSE, immediate. = TRUE, noBreaks. = TRUE)
    return(NULL)
  }

  #Column used for facetting
  fcols <- names(input_dt)[!grepl("_count|WT", names(input_dt))]

  #Add 1 + jitter
  set.seed(1)
  count_cols <- names(input_dt)[grepl("_count", names(input_dt))]
  for(i in count_cols){
    input_dt[, paste0(i) := as.numeric(.SD[[1]]+1),,.SDcols = i]
    #one uniform draw per group of rows (row index modulo 1000, in order of first appearance), as the
    #original by-group evaluation drew them; the arithmetic is unchanged
    x <- input_dt[[i]]
    grp <- seq_len(nrow(input_dt)) %% 1000
    u <- runif(length(unique(grp)))[match(grp, unique(grp))]
    data.table::set(input_dt, j = i, value = 10^(log10(x) + log10((x+1)/x) * u))
    # input_dt[, paste0(i) := .SD + runif(1),1:nrow(input_dt) %% 1000,.SDcols = i]
  }

  #Rename count columns
  names(input_dt)[grep("_count", names(input_dt))] <- dimsum__plot_samplename(names(input_dt)[grep("_count", names(input_dt))])

  #Real variants
  plot_dt <- reshape2::melt(input_dt, id = c(fcols, "WT"))
  plot_dt[, Hamming_distance := factor(.SD[[1]]),,.SDcols = fcols]
  # plot_dt[,value := value+1]

  #WT frequency
  plot_dt_wt <- plot_dt[WT==T,.SD,,.SDcols = c("variable", "value")]

  # #Most abundant variant frequency
  # plot_dt_ma <- plot_dt[,.(value=max(value)),variable]

  #Fake single expected frequency: fake_single_counts = WT_count*error_rate/4
  plot_dt_fsingles <- data.table::copy(plot_dt_wt)
  plot_dt_fsingles[, value := ((value - 1)*error_rate/4 + 1)]
  plot_dt_fsingles[, Hamming_distance := factor(1)]

  #Fake double expected frequency distribution: fake_double_counts = real_single_counts*error_rate*2/4
  plot_dt_fdoubles <- plot_dt[get(fcols[1]) == 1]
  for(i in unique(plot_dt_fdoubles[,variable])){
    #Adjust observed single counts for outflow to doubles
    plot_dt_fdoubles[variable==i,real_value := value - 1 + (value - 1 - plot_dt_fsingles[variable==i,value - 1])*seq_length*error_rate]
    #Adjust observed single counts for inflow from WT
    plot_dt_fdoubles[variable==i,real_value := real_value - plot_dt_fsingles[variable==i,value - 1]]
    #Fake doubles
    plot_dt_fdoubles[variable==i,value := real_value*error_rate*2/4 + 1]
  }
  plot_dt_fdoubles <- plot_dt_fdoubles[,.(value = median(value)),variable]
  plot_dt_fdoubles[, Hamming_distance := factor(2)]

  #Remove real variants with hamming distance ==0
  plot_dt <- plot_dt[get(fcols[1])>0,]

  #Plot
  if(length(fcols)==1){
    #Plot without subplot stratification (Nham_nt)
    d <- ggplot2::ggplot(plot_dt, ggplot2::aes(value, color = Hamming_distance)) +
      ggplot2::geom_density(key_glyph = dimsum__draw_key_polygon) + dimsum__scale_x_log10() + ggplot2::theme_bw() +
      ggplot2::labs(x = "variant count + 1 (log scale)", y = "Density", title = title, color='Nucleotide\nHamming dist.') +
      # ggplot2::geom_vline(data = plot_dt_ma, ggplot2::aes(xintercept = value), linetype = 2) +
      ggplot2::geom_vline(data = plot_dt_wt, ggplot2::aes(xintercept = value), linetype = 2) +
      ggplot2::geom_vline(data = plot_dt_fsingles, ggplot2::aes(xintercept = value, color = Hamming_distance), linetype = 2) +
      ggplot2::geom_vline(data = plot_dt_fdoubles, ggplot2::aes(xintercept = value, color = Hamming_distance), linetype = 2) +
      ggplot2::facet_grid(Hamming_distance ~ variable, scales = "free")
    dimsum_fitness__save_png_pdf(paste0(output_file_prefix, ".png"), paste0(output_file_prefix, ".pdf"), d, width = width, height = height)
  }else{
    #Plot with subplot stratification (using second facet column, Nham_aa)
    plot_list <- list()
    plot_indices <- unlist(plot_dt[order(get(fcols[2]), decreasing = F),unique(.SD),,.SDcols = fcols[2]])
    for(i in plot_indices){
      this_plot_dt <- plot_dt[get(fcols[2])==i & get(fcols[1])<=max((i*3), 3) & get(fcols[1])>=i,]
      this_plot_dt[, Hamming_distance := factor(.SD[[1]], levels = 1:max((i*3), 3)),,.SDcols = fcols]
      plot_list[[as.character(i)]] <- ggplot2::ggplot(this_plot_dt, ggplot2::aes(value, color = Hamming_distance)) +
        ggplot2::geom_density(key_glyph = dimsum__draw_key_polygon) + dimsum__scale_x_log10() + ggplot2::theme_bw() +
        ggplot2::labs(x = "variant count + 1 (log scale)", y = "Density", title = paste0("Amino acid Hamming distance: ", i), color='Nucleotide\nHamming dist.') +
        # ggplot2::geom_vline(data = plot_dt_ma, ggplot2::aes(xintercept = value), linetype = 2) +
        ggplot2::geom_vline(data = plot_dt_wt, ggplot2::aes(xintercept = value), linetype = 2) +
        ggplot2::geom_vline(data = plot_dt_fsingles, ggplot2::aes(xintercept = value, color = Hamming_distance), linetype = 2) +
        ggplot2::geom_vline(data = plot_dt_fdoubles, ggplot2::aes(xintercept = value, color = Hamming_distance), linetype = 2) +
        ggplot2::facet_grid(Hamming_distance ~ variable, scales = "free")
    }
    #Combine plots
    temp_grid.arrange <- gridExtra::arrangeGrob
    p <- do.call("temp_grid.arrange", c(plot_list, nrow = length(plot_indices)))
    dimsum_fitness__save_png_pdf(paste0(output_file_prefix, ".png"), paste0(output_file_prefix, ".pdf"), p, width = width, height = height)
  }
}

# ---- dimsum__save_metadata
dimsum__save_metadata <- function(
  dimsum_meta,
  n
  ){
  save(
    list = c(ls(parent.frame(n = n))),
    file = file.path(dimsum_meta[["project_path"]], paste0(dimsum_meta[["projectName"]], '_workspace.RData')),
    envir = parent.frame(n = n),
    version = 2)
  session_info <- sessionInfo()
  save(session_info,
  	file = file.path(dimsum_meta[["project_path"]], paste0(dimsum_meta[["projectName"]], '_sessionInfo.RData')),
  	version = 2)
}

# ---- dimsum__save_png [modified, see the header and the generator]
dimsum__save_png <- function(
  output_file,
  plot_object,
  width = 7,
  height = 7,
  plot_dpi = 300){

  #Open device
  grDevices::png(
    output_file,
    width = width,
    height = height,
    units = "in",
    res = plot_dpi,
    bg = "white",
    type = "cairo")
  #Plot
  if(isTRUE(names(plot_object)[1]=="grobs")){
    plot(plot_object)
  }else{
    print(plot_object)
  }
  #Close device
  { sink("/dev/null")
    dev.off()
    sink()
  }
}

# ---- dimsum__status_message
dimsum__status_message <- function(
	...,
	newline=F
	){
	cat(sprintf(...), sep='', file=stdout())
	if(newline){
		cat("\n")
	}
}

# ---- dimsum__sum_datatable_columns
dimsum__sum_datatable_columns <- function(
  dt,
  column_patterns,
  suffix=""){
  for(this_column_pattern in column_patterns){
    #Sum columns with given pattern
    temp_data <- apply(dt[,grep(this_column_pattern, colnames(dt)), with=FALSE], 1, sum)
    #Remove summed columns
    dt <- dt[,-grep(this_column_pattern, colnames(dt)), with=FALSE]
    #Append result column
    dt[,paste0(this_column_pattern, suffix)] <- temp_data
  }
  #Return data.table
  return(dt)
}

# ---- dimsum__validate_input
dimsum__validate_input <- function(
  input_list
  ){

  input_list_original <- input_list

  #Check all input parameter types correct
  for(input_arg in names(input_list)){
    #Integer arguments
    if("integer" %in% input_list[[input_arg]][[2]]){
      #Check if numeric
      if(!typeof(input_list[[input_arg]][[1]]) %in% c(input_list[[input_arg]][[2]], "double")){
        stop(paste0("Invalid type of argument '", input_arg, "' (", input_list[[input_arg]][[2]], ")"), call. = FALSE)
      }else{
        #Check if whole number (if not NULL)
        if(!is.null(input_list[[input_arg]][[1]])){
          if(input_list[[input_arg]][[1]]%%1!=0){
            stop(paste0("Invalid type of argument '", input_arg, "' (", input_list[[input_arg]][[2]], ")"), call. = FALSE)
          }
          #Convert to integer
          input_list[[input_arg]][[1]] <- as.integer(input_list[[input_arg]][[1]])
        }
      }
    }else{
      #Other arguments
      if(!typeof(input_list[[input_arg]][[1]]) %in% input_list[[input_arg]][[2]]){
        stop(paste0("Invalid type of argument '", input_arg, "' (", input_list[[input_arg]][[2]], ")"), call. = FALSE)
      }
    }
  }
  #Metadata object
  dimsum_meta <- sapply(input_list, '[', 1)

  #Save original command-line arguments
  dimsum_meta[["arg_list"]] <- sapply(input_list_original, '[', 1)

  #Either fastqFileDir or countPath need to be specified
  if(!is.null(dimsum_meta[["fastqFileDir"]])){
    #Reformat fastqFileDir and check if exists (remove trailing "/" if present)
    dimsum_meta[["fastqFileDir"]] <- gsub("/$", "", dimsum_meta[["fastqFileDir"]])
    if(!file.exists(dimsum_meta[["fastqFileDir"]])){
      stop(paste0("Invalid 'fastqFileDir' argument (directory not found)"), call. = FALSE)
    }
  }else if(!is.null(dimsum_meta[["countPath"]])){
    #Check if exists
    if(!file.exists(dimsum_meta[["countPath"]])){
      stop(paste0("Invalid 'countPath' argument (file not found)"), call. = FALSE)
    }
  }else{
    stop(paste0("Either 'fastqFileDir' or 'countPath' arguments must be specified"), call. = FALSE)
  }

  #Check if synonymSequencePath specified
  if(!is.null(dimsum_meta[["synonymSequencePath"]])){
    #Check if exists
    if(!file.exists(dimsum_meta[["synonymSequencePath"]])){
      stop(paste0("Invalid 'synonymSequencePath' argument (file not found)"), call. = FALSE)
    }
  }

  #Reformat outputPath and check if exists (remove trailing "/" if present)
  dimsum_meta[["outputPath"]] <- gsub("/$", "", dimsum_meta[["outputPath"]])
  if(!file.exists(dimsum_meta[["outputPath"]])){
    stop(paste0("Invalid 'outputPath' argument (directory not found)"), call. = FALSE)
  }

  #FASTQ file extension
  #Append "." if not first character
  dimsum_meta[["fastqFileExtension"]] <- ifelse(
    substr(dimsum_meta[["fastqFileExtension"]], 1, 1)!=".",
    paste0(".", dimsum_meta[["fastqFileExtension"]]),
    dimsum_meta[["fastqFileExtension"]])
  #Check if alphanumeric
  if(!grepl("^[A-Za-z0-9]+$", substr(dimsum_meta[["fastqFileExtension"]], 2, nchar(dimsum_meta[["fastqFileExtension"]])), perl = T)){
    stop(paste0("Invalid '", "fastqFileExtension", "' argument (only alphanumeric characters allowed)"), call. = FALSE)
  }

  #WT sequence
  #Check WT sequence is valid (ACGT characters only) and save case-coded sequence
  all_characters <- unique(unlist(strsplit(dimsum_meta[["wildtypeSequence"]], "")))
  if(sum(!all_characters %in% c("A", "a", "C", "c", "G", "g", "T", "t"))!=0){
    stop("Invalid wild-type nucleotide sequence. Only valid nucleotide sequences allowed (A/C/G/T). Use lower-case letters (a/c/g/t) to indicate internal constant regions (removed before fitness calculations).", call. = FALSE)
  }
  #Check WT sequence contains at least one mutated position
  if(sum(all_characters %in% c("A", "C", "G", "T"))==0){
    stop("Invalid wild-type nucleotide sequence. Must contain at least one mutated position (A/C/G/T).", call. = FALSE)
  }
  #Case-coded WT nucleotide sequence
  dimsum_meta[["wildtypeSequenceCoded"]] <- dimsum_meta[["wildtypeSequence"]]
  dimsum_meta[["wildtypeSequence"]] <- toupper(dimsum_meta[["wildtypeSequence"]])

  #permittedSequences
  all_characters <- unlist(strsplit(dimsum_meta[["wildtypeSequenceCoded"]], ""))
  #Set permittedSequences argument if not specified
  if(is.null(dimsum_meta[["permittedSequences"]])){
    dimsum_meta[["permittedSequences"]] <- paste0(rep("N", sum(all_characters %in% c("A", "C", "G", "T"))), collapse = "")
  }
  #Check permittedSequences argument has same length as mutated sequence
  if(nchar(dimsum_meta[["permittedSequences"]]) != sum(all_characters %in% c("A", "C", "G", "T"))){
    stop("Invalid 'permittedSequences' argument length. Length should match the number of mutated positions in the wild-type nucleotide sequence (A/C/G/T).", call. = FALSE)
  }
  #Check permittedSequences argument sequence is valid
  all_characters <- unique(unlist(strsplit(dimsum_meta[["permittedSequences"]], "")))
  if(sum(!all_characters %in% c("A", "C", "G", "T", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"))!=0){
    stop("Invalid 'permittedSequences' argument. Only valid nucleotide codes allowed (A/C/G/T/R/Y/S/W/K/M/B/D/H/V/N).", call. = FALSE)
  }

  #Check strictly positive integer vsearchMinQual arguments
  if(dimsum_meta[["vsearchMinQual"]]<=0){
    stop("Invalid 'vsearchMinQual' argument. Only positive integers allowed (zero exclusive).", call. = FALSE)
  }

  #Check strictly positive integer vsearchMinovlen arguments
  if(dimsum_meta[["vsearchMinovlen"]]<5){
    stop("Invalid 'vsearchMinovlen' argument. Only integers greater than 4 allowed.", call. = FALSE)
  }

  #Check positive integer fitnessHighConfidenceCount and fitnessDoubleHighConfidenceCount arguments
  if(sum(unlist(dimsum_meta[c("fitnessHighConfidenceCount", "fitnessDoubleHighConfidenceCount")])<0)!=0){
    stop("Invalid 'fitnessHighConfidenceCount' or 'fitnessDoubleHighConfidenceCount' argument. Only positive integers allowed (zero inclusive).", call. = FALSE)
  }

  #Check positive integer minimum read count arguments
  all_counts <- unlist(dimsum_meta[c("fitnessMinInputCountAll", "fitnessMinOutputCountAll", "fitnessMinInputCountAny", "fitnessMinOutputCountAny")])
  all_counts <- unlist(strsplit(all_counts, ",|:"))
  if(sum(!unique(unlist(strsplit(all_counts, ""))) %in% as.character(0:9))!=0){
    stop("Invalid 'fitness...Count...' arguments. Only positive integers allowed (zero inclusive).", call. = FALSE)
  }else{
    #Parse minimum read count arguments
    for(i in c("fitnessMinInputCountAll", "fitnessMinOutputCountAll", "fitnessMinInputCountAny", "fitnessMinOutputCountAny")){
      dimsum_meta[[i]] <- suppressWarnings(dimsum__parse_minimum_read_count_arguments(dimsum_meta[[i]]))
    }
  }

  #Check fitnessDropoutPseudocount argument
  if(dimsum_meta[["fitnessDropoutPseudocount"]]<0){
    stop("Invalid 'fitnessDropoutPseudocount' argument. Only positive integers allowed (zero inclusive).", call. = FALSE)
  }

  #Check strictly positive double splitChunkSize argument
  if(dimsum_meta[["splitChunkSize"]]<=0){
    stop("Invalid 'splitChunkSize' argument. Only positive integers allowed (zero exclusive).", call. = FALSE)
  }

  #Check maxSubstitutions argument greater than 1
  if(dimsum_meta[["maxSubstitutions"]]<2){
    stop("Invalid 'maxSubstitutions' argument. Only integers greater than 1 allowed.", call. = FALSE)
  }

  #Check indels argument valid and set indelLengths
  dimsum_meta[["indelLengths"]] <- NA
  if(dimsum_meta[["indels"]]=="none"){
    #All indels discarded
    dimsum_meta[["indels"]] <- FALSE
  }else if(dimsum_meta[["indels"]]=="all"){
    #All indels retained (regardless of length)
    dimsum_meta[["indels"]] <- TRUE
  }else{
    #List of indel lengths specified
    dimsum_meta[["indelLengths"]] <- strtoi(unlist(strsplit(dimsum_meta[["indels"]], ",")))
    dimsum_meta[["indelLengths"]] <- dimsum_meta[["indelLengths"]][!is.na(dimsum_meta[["indelLengths"]])]
    if(length(dimsum_meta[["indelLengths"]])==0){
      stop("Invalid 'indels' argument. Only integers greater than 1 allowed.", call. = FALSE)
    }
    dimsum_meta[["indels"]] <- TRUE
  }

  #Check barcodeErrorRate argument
  if(dimsum_meta[["barcodeErrorRate"]]<0 | dimsum_meta[["barcodeErrorRate"]]>=1){
    stop("Invalid 'barcodeErrorRate' argument. Only positive doubles less than 1 allowed (zero inclusive).", call. = FALSE)
  }

  #Check vsearchMaxee argument
  if(dimsum_meta[["vsearchMaxee"]]<=0){
    stop("Invalid 'vsearchMaxee' argument. Only positive doubles allowed (zero exclusive).", call. = FALSE)
  }

  #Check vsearchMaxQual argument
  if(dimsum_meta[["vsearchMaxQual"]]>93){
    stop("Invalid 'vsearchMaxQual' argument. Cannot be greater than 93.", call. = FALSE)
  }

  #Disable bayesianDoubleFitness option (TEMPORARY FIX: still in development)
  dimsum_meta[["bayesianDoubleFitness"]] <- FALSE

  #Check bayesianDoubleFitnessLamD argument
  if(dimsum_meta[["bayesianDoubleFitnessLamD"]]<=0){
    stop("Invalid 'bayesianDoubleFitnessLamD' argument. Only positive doubles allowed (zero exclusive).", call. = FALSE)
  }

  #Check startStage argument
  if(dimsum_meta[["startStage"]]<0){
    stop("Invalid 'startStage' argument. Only positive integers allowed (zero inclusive).", call. = FALSE)
  }

  #Check stopStage argument
  if(dimsum_meta[["stopStage"]]<=0){
    stop("Invalid 'stopStage' argument. Only positive integers allowed (zero exclusive).", call. = FALSE)
  }

  #Check numCores argument
  if(dimsum_meta[["numCores"]]<=0){
    stop("Invalid 'numCores' argument. Only positive integers allowed (zero exclusive).", call. = FALSE)
  }

  #Check library design paired if trans library specified
  if(dimsum_meta[["transLibrary"]] & !dimsum_meta[["paired"]]){
    stop("Invalid 'paired' argument. Only paired-end trans libraries allowed.", call. = FALSE)
  }

  #Check library design paired if unstranded library specified
  if(!dimsum_meta[["stranded"]] & !dimsum_meta[["paired"]]){
    stop("Invalid 'stranded' argument. Only unstranded paired-end libraries allowed.", call. = FALSE)
  }

  #Check sequenceType one of noncoding/coding/auto
  if(!dimsum_meta[["sequenceType"]] %in% c("noncoding", "coding", "auto")){
    stop("Invalid 'sequenceType' argument. Only noncoding/coding/auto allowed.", call. = FALSE)
  }else if(dimsum_meta[["sequenceType"]]=="auto"){
    dimsum_meta[["sequenceType"]] <- dimsum__detect_sequence_type(dimsum_meta[["wildtypeSequence"]])
  }

  #Check mutagenesisType one of noncoding/coding/auto
  if(!dimsum_meta[["mutagenesisType"]] %in% c("random", "codon")){
    stop("Invalid 'mutagenesisType' argument. Only random/codon allowed.", call. = FALSE)
  }

  #Check retainedReplicates comma-separated list of integers or "all"
  if(dimsum_meta[["retainedReplicates"]] != "all"){
    if(sum(!strsplit(dimsum_meta[["retainedReplicates"]], ",")[[1]] %in% strsplit("0123456789", "")[[1]])!=0){
      stop("Invalid 'retainedReplicates' argument. Only comma-separated list of integers or 'all' allowed.", call. = FALSE)
    }
  }

  #Check if barcodeIdentityPath specified
  if(!is.null(dimsum_meta[["barcodeIdentityPath"]])){
    #Load barcode identity table
    if(!file.exists(dimsum_meta[["barcodeIdentityPath"]])){
      stop(paste0("Invalid '", "barcodeIdentityPath", "' argument (file not found)"), call. = FALSE)
    }
  }

  #Set fitnessNormalise to FALSE if fitnessErrorModel is FALSE
  if(!dimsum_meta[["fitnessErrorModel"]]){
    dimsum_meta[["fitnessNormalise"]] <- FALSE
  }

  #Set run start/end time
  dimsum_meta[["start_time"]] <- Sys.time()
  dimsum_meta[["end_time"]] <- "In progress"

  #Return
  return(dimsum_meta)
}

# ---- dimsum_stage_counts_to_fitness
dimsum_stage_counts_to_fitness <- function(
  dimsum_meta,
  fitness_outpath,
  report = TRUE,
  report_outpath = NULL,
  save_workspace = TRUE
  ){

  #Whether or not to execute the system command
  this_stage <- 5
  execute <- (dimsum_meta[["startStage"]] <= this_stage & dimsum_meta[["stopStage"]] >= this_stage)

  #This stage after stopStage
  if(dimsum_meta[["stopStage"]] < this_stage){
    return(dimsum_meta)
  }

  #Save current workspace for debugging purposes
  if(save_workspace){dimsum__save_metadata(dimsum_meta = dimsum_meta, n = 2)}

  #Fitness path
  dimsum_meta[["fitness_path"]] <- file.path(fitness_outpath)

  #Do nothing if analysis not executed
  if(!execute){
    return(dimsum_meta)
  }

  #Create/overwrite fitness directory (if executed)
  fitness_outpath <- gsub("/$", "", fitness_outpath)
  dimsum__create_dir(fitness_outpath, execute = execute, message = "DiMSum STAGE 5 (STEAM): ANALYSE VARIANT COUNTS", overwrite_dir = FALSE)

  #Check if all input files exist
  dimsum__check_files_exist(
    required_files = dimsum_meta[["variant_data_merge_path"]],
    stage_number = this_stage,
    execute = execute)

  ### Load variant data
  ###########################

  #load variant data from RData file
  e1 <- new.env()
  load(dimsum_meta[["variant_data_merge_path"]])
  all_data <- get('variant_data_merge', e1)
  dimsum__check_variants(dimsum_meta = dimsum_meta, input_dt = all_data)
  rm(e1)

  #WT nucleotide sequences
  wt_ntseq <- all_data[WT==T,nt_seq]

  #WT AA sequences
  wt_AAseq <- all_data[WT==T,aa_seq]

  #Determine replicates to retain
  all_reps <- unique(dimsum_meta[["exp_design"]][,"experiment"])
  if(dimsum_meta[["retainedReplicates"]]!="all"){
    all_reps <- unique(as.integer(strsplit(dimsum_meta[["retainedReplicates"]], ",")[[1]]))
  }

  #Check if only one experimental replicate exists (requirement for normalisation and error model fit)
  if(length(unique(dimsum_meta[["exp_design"]][,"experiment_replicate"]))<2){
    dimsum_meta[["fitnessNormalise"]] <- FALSE
    dimsum_meta[["fitnessErrorModel"]] <- FALSE
  }

  ### Filter out low count nucleotide variants
  ###########################

  nf_data <- dimsum__filter_nuc_variants(
    dimsum_meta = dimsum_meta,
    input_dt = all_data,
    all_reps = all_reps)

  ### Identify variants to use for error modelling
  ###########################

  nf_data[, error_model := T]

  ### Aggregate counts at the AA level (if coding sequence and "mixedSubstitutions"==T)
  ###########################

  #For coding sequences aggregate variant counts at the AA level
  if(dimsum_meta[["sequenceType"]]=="coding" & dimsum_meta[["mixedSubstitutions"]]){
    nf_data <- dimsum__aggregate_AA_variants(
      dimsum_meta = dimsum_meta,
      input_dt = nf_data,
      all_reps = all_reps)
  }else{
    nf_data[,merge_seq := nt_seq]
  }

  ### Aggregate counts for biological output replicates
  ###########################

  dimsum__status_message("Aggregating counts for biological output replicates...\n")

  #Add up counts for biological output reps
  for (E in all_reps) {
    idx <- names(nf_data)[grep(names(nf_data),pattern = paste0("e",E,"_s1_b"))]
    nf_data[!rowSums(is.na(nf_data[,.SD,.SDcols = idx]))==length(idx),paste0("count_e",E,"_s1") := rowSums(.SD, na.rm = T),,.SDcols = idx]
    names(nf_data)[grep(names(nf_data),pattern = paste0("e",E,"_s0_b"))] <- paste0("count_e",E,"_s0")
  }
  nf_data <- nf_data[,.SD,,.SDcols = c(
    "merge_seq","nt_seq","aa_seq","Nham_nt","Nham_aa",
    "Nmut_codons","WT","indel","STOP","STOP_readthrough","error_model",names(nf_data)[grep(names(nf_data),pattern="^count")])]

  dimsum__status_message("Done\n")

  ### Fit error model
  ###########################

  #Fit error model
  model_result <- dimsum__error_model(
    dimsum_meta = dimsum_meta,
    input_dt = data.table::copy(nf_data[error_model==T]),
    all_reps = all_reps,
    report_outpath = report_outpath)

  ### Add psuedocounts to output samples with dropout
  ###########################

  nf_data_pseudo <- dimsum__add_dropout_pseudocount(
    dimsum_meta = dimsum_meta,
    input_dt = nf_data,
    all_reps = all_reps)

  ### Calculate fitness and count-based error (and remove variants without fitness estimates in any replicate)
  ###########################

  nff_data <- dimsum__calculate_fitness(
    dimsum_meta = dimsum_meta,
    input_dt = nf_data_pseudo,
    all_reps = all_reps,
    error_model_dt = model_result[["error_model"]],
    norm_model_dt = model_result[["norm_model"]])

  ### Merge fitness and error at the AA level (if coding sequence and "mixedSubstitutions"==F)
  ###########################

  #Aggregate nonsynonymous variant fitness and error at the AA level
  if(dimsum_meta[["sequenceType"]]=="coding" & !dimsum_meta[["mixedSubstitutions"]]){
    nff_data <- dimsum__aggregate_AA_variants_fitness(
      dimsum_meta = dimsum_meta,
      input_dt = nff_data,
      all_reps = all_reps)
  }else{
    if(dimsum_meta[["sequenceType"]]!="coding"){nff_data[,merge_seq := nt_seq,nt_seq]}
    nff_data <- nff_data[,.SD,merge_seq,.SDcols = c(
      "aa_seq","Nham_nt","Nham_aa",
      "Nmut_codons","WT","indel","STOP","STOP_readthrough", "error_model",
      names(nff_data)[grep(names(nff_data),pattern="^count")],
      "mean_count",
      names(nff_data)[grep(names(nff_data),pattern="^fitness|sigma")])]
  }

  ### Wild type
  ###########################

  wildtype <- nff_data[WT==TRUE & error_model==TRUE,]

  ### Identify position and identity of single AA/NT mutations (and all silent mutants in the case of AA mutations)
  ###########################

  if(dimsum_meta[["sequenceType"]]=="coding"){
    singles <- dimsum__identify_single_aa_mutations(
      input_dt = nff_data[error_model==T],
      wt_AAseq = wt_AAseq)
  }else{
    singles <- dimsum__identify_single_nt_mutations(
      input_dt = nff_data[error_model==T],
      wt_ntseq = wt_ntseq)
  }

  ### Identify position and identity of double AA/NT mutations
  ###########################

  if(dimsum_meta[["sequenceType"]]=="coding"){
    doubles <- dimsum__identify_double_aa_mutations(
      input_dt = nff_data[error_model==T],
      singles_dt = singles,
      wt_AAseq = wt_AAseq)
  }else{
    doubles <- dimsum__identify_double_nt_mutations(
      input_dt = nff_data[error_model==T],
      singles_dt = singles,
      wt_ntseq = wt_ntseq)
  }

  ### Bayesian framework for fitness estimation for double mutants
  ###########################

  #Bin mean counts for replicate 1
  # doubles[,counts_for_bins := .SD[[1]],,.SDcols = paste0("count_e",all_reps[1],"_s0")]
  # doubles[,bin_count := findInterval(log10(counts_for_bins),seq(0.5,4,0.25))]
  # doubles[,.(.N,mean(counts_for_bins)),bin_count][order(bin_count)]

  if(!dimsum_meta[["bayesianDoubleFitness"]]){
    dimsum__status_message("Skipping Bayesian double mutant fitness estimation\n")
  }else{
    doubles <- dimsum__bayesian_double_fitness(
      dimsum_meta = dimsum_meta,
      doubles_dt = doubles,
      singles_dt = singles,
      wt_dt = wildtype,
      all_reps = all_reps,
      report_outpath = report_outpath)
  }

  ### Normalise fitness and sigma for differences in number of generations (between biological replicates)
  ###########################

  #Check that number of generations supplied for all output samples
  if(sum(is.na(dimsum_meta[["exp_design"]][dimsum_meta[["exp_design"]][,"selection_id"]==1,"generations"]))!=0){

    dimsum__status_message("Skipping normalisation by number of generations\n")

  }else{

    dimsum__status_message("Normalising fitness and error for differences in number of generations...\n")

    nff_data <- dimsum__normalise_fitness(
      dimsum_meta = dimsum_meta,
      input_dt = nff_data,
      all_reps = all_reps,
      fitness_suffix="_uncorr"
      )
    singles <- dimsum__normalise_fitness(
      dimsum_meta = dimsum_meta,
      input_dt = singles,
      all_reps = all_reps
      )
    doubles <- dimsum__normalise_fitness(
      dimsum_meta = dimsum_meta,
      input_dt = doubles,
      all_reps = all_reps,
      fitness_suffix="_uncorr"
      )
    if(dimsum_meta[["bayesianDoubleFitness"]]){
      doubles <- dimsum__normalise_fitness(
        dimsum_meta = dimsum_meta,
        input_dt = doubles,
        all_reps = all_reps,
        fitness_suffix="_cond"
        )
    }

    dimsum__status_message("Done\n")

  }

  ### Merge fitness values and save fitness data
  ###########################

  dimsum__merge_fitness(
    dimsum_meta = dimsum_meta,
    input_dt = nff_data,
    doubles_dt = doubles,
    singles_dt = singles,
    all_reps = all_reps,
    fitness_outpath = fitness_outpath,
    report = TRUE,
    report_outpath = report_outpath)

  #Delete files when last stage complete
  if(!dimsum_meta[["retainIntermediateFiles"]]){
    if(dimsum_meta[["stopStage"]]==this_stage){
      if(!is.null(dimsum_meta[["deleteIntermediateFiles"]])){suppressWarnings(temp_out <- file.remove(dimsum_meta[["deleteIntermediateFiles"]]))}
      if(!is.null(dimsum_meta[["touchIntermediateFiles"]])){suppressWarnings(temp_out <- file.create(dimsum_meta[["touchIntermediateFiles"]]))}
    }
  }

  return(dimsum_meta)
}

# ---- dimsum_stage_merge
dimsum_stage_merge <- function(
  dimsum_meta,
  merge_outpath,
  report = TRUE,
  report_outpath = NULL,
  save_workspace = TRUE
  ){

  #Whether or not to execute the system command
  this_stage <- 4
  execute <- (dimsum_meta[["startStage"]] <= this_stage & dimsum_meta[["stopStage"]] >= this_stage)

  #This stage after stopStage
  if(dimsum_meta[["stopStage"]] < this_stage){
    return(dimsum_meta)
  }

  #Save current workspace for debugging purposes
  if(save_workspace){dimsum__save_metadata(dimsum_meta = dimsum_meta, n = 2)}

  #Create merge directory (if doesn't already exist)
  merge_outpath <- gsub("/$", "", merge_outpath)
  dimsum__create_dir(merge_outpath, execute = execute, message = "DiMSum STAGE 4 (STEAM): PROCESS VARIANT SEQUENCES", overwrite_dir = FALSE)

  #Input files
  if(!is.null(dimsum_meta[["countPath"]])){
    all_count <- dimsum_meta[["countPath"]]
  }else{
    all_count <- unique(file.path(dimsum_meta[['exp_design']][,"aligned_pair_unique_directory"], dimsum_meta[['exp_design']][,"aligned_pair_unique"]))
  }

  #Check if all input files exist
  dimsum__check_files_exist(
    required_files = all_count,
    stage_number = this_stage,
    execute = execute)

  #WT nucleotide sequence
  wt_ntseq <- tolower(dimsum_meta[['wildtypeSequence']])
  #WT AA sequence
  wt_AAseq <- paste0(seqinr::translate(strsplit(wt_ntseq,split="")[[1]]),collapse="")
  #Initialise merged variant data object
  variant_data <- NULL

  #Load variant count files
  dimsum__status_message("Loading variant count files:\n")
  dimsum__status_message(paste0(all_count, "\n"))
  dimsum__status_message("Processing...\n")
  for(count_file in all_count){
    dimsum__status_message(paste0("\t", basename(count_file), "\n"))
    #Check if this code should be executed
    if(execute){

      #load variant data from user-specified count file
      if(!is.null(dimsum_meta[["countPath"]])){
        variant_data <- dimsum__check_countfile(dimsum_meta = dimsum_meta, return_data = TRUE)
        break
      }

      #load variant data from dimsum variant count files
      file_id <- gsub('.vsearch.unique', '', basename(count_file))
      #Initialise count table
      count_dt <- data.table(
        nt_seq = character(),
        count = numeric())
      #Load fasta file (if exists)
      if(!file.exists(count_file)){
        warning("dimsum_stage_merge.R: File does not exist.", call. = FALSE, immediate. = TRUE, noBreaks. = TRUE)
      }else{
        #Create variant count table with nucleotide and amino acid sequence
        temp_dt <- fread(count_file)
        if(dim(temp_dt)[1]!=0){count_dt <- temp_dt}
        names(count_dt) <- c("nt_seq", "count")
        count_dt[, nt_seq := tolower(nt_seq)]
      }
      #First file loaded
      if(is.null(variant_data)){
        variant_data <- count_dt
      #Not first file loaded (merge with previous data)
      }else{
        variant_data <- merge(variant_data,count_dt,by=names(variant_data)[-grep("count", names(variant_data))],all = T)
      }
      names(variant_data)[names(variant_data)=="count"] = paste0(file_id, "_count")
    }
  }
  #Check if this code should be executed
  if(execute){
    #Replace NA counts with zeros
    variant_data[is.na(variant_data)] <- 0

    #Merge counts from technical replicates (if count file not specified)
    if(!is.null(dimsum_meta[["countPath"]])){
      variant_data_merge <- variant_data
    }else{
      split_base <- unique(sapply(strsplit(colnames(variant_data)[grep('_count', colnames(variant_data))], "_t"), '[', 1))
      variant_data_merge <- dimsum__sum_datatable_columns(dt=variant_data, column_patterns=split_base, suffix="_count")
    }

    #Process variants to separate into indel, rejected and retained variants data.tables
    variant_list <- dimsum__process_merged_variants(dimsum_meta = dimsum_meta, input_dt = variant_data_merge)

    #Save merged variant data
    dimsum__status_message("Saving merged variant data...\n")
    write.table(variant_list[["nobarcode_variants"]], file=file.path(merge_outpath, paste0(dimsum_meta[["projectName"]], '_nobarcode_variant_data_merge.tsv')), sep = "\t", quote = F, row.names = F)
    write.table(variant_list[["indel_variants"]], file=file.path(merge_outpath, paste0(dimsum_meta[["projectName"]], '_indel_variant_data_merge.tsv')), sep = "\t", quote = F, row.names = F)
    write.table(variant_list[["rejected_variants"]], file=file.path(merge_outpath, paste0(dimsum_meta[["projectName"]], '_rejected_variant_data_merge.tsv')), sep = "\t", quote = F, row.names = F)
    variant_data_merge <- variant_list[["retained_variants"]]
    save(variant_data_merge,
      file=file.path(merge_outpath, paste0(dimsum_meta[["projectName"]], '_variant_data_merge.RData')),
      version = 2)
    write.table(variant_list[["retained_variants"]], file=file.path(merge_outpath, paste0(dimsum_meta[["projectName"]], '_variant_data_merge.tsv')), sep = "\t", quote = F, row.names = F)
    dimsum__status_message("Done\n")
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Merged variant data path
  dimsum_meta_new[["variant_data_merge_path"]] <- file.path(merge_outpath, paste0(dimsum_meta[["projectName"]], '_variant_data_merge.RData'))

  #Check if mutation stats file exists
  dimsum__check_files_exist(
    required_files = file.path(dimsum_meta[["tmp_path"]], "mutation_stats_dicts.RData"),
    stage_number = this_stage,
    execute = execute)

  #Load previously saved mutation statistics
  load(file.path(dimsum_meta[["tmp_path"]], "mutation_stats_dicts.RData"))

  #Add to metadata
  dimsum_meta_new[["aa_subst_counts"]] <- mutation_stats_dicts[["aa_subst_dict"]]
  dimsum_meta_new[["nuc_subst_counts"]] <- mutation_stats_dicts[["nuc_subst_dict"]]
  dimsum_meta_new[["nuc_indrt_counts"]] <- mutation_stats_dicts[["nuc_indrt_dict"]]
  dimsum_meta_new[["nuc_mxsub_counts"]] <- mutation_stats_dicts[["nuc_mxsub_dict"]]
  dimsum_meta_new[["nuc_tmsub_counts"]] <- mutation_stats_dicts[["nuc_tmsub_dict"]]
  dimsum_meta_new[["nuc_frbdn_counts"]] <- mutation_stats_dicts[["nuc_frbdn_dict"]]
  dimsum_meta_new[["nuc_const_counts"]] <- mutation_stats_dicts[["nuc_const_dict"]]
  dimsum_meta_new[["nuc_indel_counts"]] <- mutation_stats_dicts[["nuc_indel_dict"]]
  dimsum_meta_new[["nuc_nbarc_counts"]] <- mutation_stats_dicts[["nuc_nbarc_dict"]]
  #Delete files when last stage complete
  if(!dimsum_meta_new[["retainIntermediateFiles"]]){
    if(dimsum_meta_new[["stopStage"]]==this_stage){
      if(!is.null(dimsum_meta_new[["deleteIntermediateFiles"]])){suppressWarnings(temp_out <- file.remove(dimsum_meta_new[["deleteIntermediateFiles"]]))}
      if(!is.null(dimsum_meta_new[["touchIntermediateFiles"]])){suppressWarnings(temp_out <- file.create(dimsum_meta_new[["touchIntermediateFiles"]]))}
    }
  }
  #Generate merge report
  if(report){
    tryCatch({
      dimsum_meta_new_report <- dimsum__merge_report(dimsum_meta = dimsum_meta_new, report_outpath = report_outpath)
      }, error=function(e){
        dimsum__status_message("There were problems while running 'dimsum__merge_report'\n")
        dimsum_meta_new_report <- dimsum_meta_new
        })
    tryCatch({
      suppressWarnings(dimsum__diagnostics_report(dimsum_meta = dimsum_meta_new_report, report_outpath = report_outpath))
      }, error=function(e){
        dimsum__status_message("There were problems while running 'dimsum__diagnostics_report'\n")
        })
    return(dimsum_meta_new_report)
  }
  return(dimsum_meta_new)
}

####################################################################################################
# Driver: what dimsum() does when a count file is given, then STEAM stages 4 and 5
####################################################################################################

dimsum_fitness <- function(arg_list) {

  #Display welcome
  dimsum__status_message(paste("\n\n\n*******", "Running DiMSum STEAM stages", "*******\n\n\n"))
  dimsum__status_message(paste(formatDL(unlist(list(
    "DiMSum version" = DIMSUM_FITNESS_VERSION,
    "R version" = version$version.string))), collapse = "\n"), newline = T)

  ### Basic checks
  ###########################

  #Required binaries
  required_binaries <- c(
    "pandoc")
  which_binaries <- Sys.which(required_binaries)
  missing_binaries <- names(which_binaries)[which_binaries==""]
  if(length(missing_binaries)!=0){
    stop(paste0("Required executables not installed. Please install the following software: ", paste(missing_binaries, sep = ", ")), call. = FALSE)
  }

  #Binary versions
  suppressMessages(suppressWarnings(binary_versions <- list(
    Pandoc = system("pandoc -v", intern = TRUE)[1])))

  #Display binary versions
  dimsum__status_message(paste("\n\n\n*******", "Binary dependency versions", "*******\n\n\n"))
  dimsum__status_message(paste(formatDL(unlist(binary_versions)), collapse = "\n"), newline = T)

  ### Setup
  ###########################

  dimsum_arg_list <- list(
    "runDemo" = list(arg_list[["runDemo"]], c("logical")), #logical -- checked in dimsum__validate_input
    "fastqFileDir" = list(arg_list[["fastqFileDir"]], c("character", "NULL")), #directory exists if running WRAP stages -- checked in dimsum__validate_input
    "fastqFileExtension" = list(arg_list[["fastqFileExtension"]], c("character")), #alphanumeric character string starting with '.' -- checked in dimsum__validate_input
    "gzipped" = list(arg_list[["gzipped"]], c("logical")), #logical -- checked in dimsum__validate_input
    "stranded" = list(arg_list[["stranded"]], c("logical")), #logical -- checked in dimsum__validate_input
    "paired" = list(arg_list[["paired"]], c("logical")), #logical -- checked in dimsum__validate_input
    "barcodeDesignPath" = list(arg_list[["barcodeDesignPath"]], c("character", "NULL")), #file exists (if not NULL)
    "barcodeErrorRate" = list(arg_list[["barcodeErrorRate"]], c("double")), #positive double less than 1 (zero inclusive) -- checked in dimsum__validate_input
    "experimentDesignPath" = list(arg_list[["experimentDesignPath"]], c("character")), #file exists -- checked in dimsum__get_experiment_design
    "experimentDesignPairDuplicates" = list(arg_list[["experimentDesignPairDuplicates"]], c("logical")), #logical -- checked in dimsum__validate_input
    "barcodeIdentityPath" = list(arg_list[["barcodeIdentityPath"]], c("character", "NULL")), #file exists (if not NULL) -- checked in dimsum__validate_input
    "countPath" = list(arg_list[["countPath"]], c("character", "NULL")), #file exists (if not NULL) -- checked in dimsum__validate_input
    "synonymSequencePath" = list(arg_list[["synonymSequencePath"]], c("character", "NULL")), #file exists (if not NULL) -- checked in dimsum__validate_input
    "cutadaptCut5First" = list(arg_list[["cutadaptCut5First"]], c("integer", "NULL")), #strictly positive integer (if not NULL) -- checked in dimsum__check_experiment_design
    "cutadaptCut5Second" = list(arg_list[["cutadaptCut5Second"]], c("integer", "NULL")), #strictly positive integer (if not NULL) -- checked in dimsum__check_experiment_design
    "cutadaptCut3First" = list(arg_list[["cutadaptCut3First"]], c("integer", "NULL")), #strictly positive integer (if not NULL) -- checked in dimsum__check_experiment_design
    "cutadaptCut3Second" = list(arg_list[["cutadaptCut3Second"]], c("integer", "NULL")), #strictly positive integer (if not NULL) -- checked in dimsum__check_experiment_design
    "cutadapt5First" = list(arg_list[["cutadapt5First"]], c("character", "NULL")), #AGCT character string (if not NULL) -- checked in dimsum__check_experiment_design
    "cutadapt5Second" = list(arg_list[["cutadapt5Second"]], c("character", "NULL")), #AGCT character string (if not NULL) -- checked in dimsum__check_experiment_design
    "cutadapt3First" = list(arg_list[["cutadapt3First"]], c("character", "NULL")), #AGCT character string (if not NULL) -- checked in dimsum__check_experiment_design
    "cutadapt3Second" = list(arg_list[["cutadapt3Second"]], c("character", "NULL")), #AGCT character string (if not NULL) -- checked in dimsum__check_experiment_design
    "cutadaptMinLength" = list(arg_list[["cutadaptMinLength"]], c("integer")), #strictly positive integer -- checked in dimsum__check_experiment_design
    "cutadaptErrorRate" = list(arg_list[["cutadaptErrorRate"]], c("double")), #positive double less than 1 (zero inclusive) -- checked in dimsum__check_experiment_design
    "cutadaptOverlap" = list(arg_list[["cutadaptOverlap"]], c("integer")), #positive integer (zero inclusive) -- checked in dimsum__check_experiment_design
    "vsearchMinQual" = list(arg_list[["vsearchMinQual"]], c("integer")), #strictly positive integer -- checked in dimsum__validate_input
    "vsearchMaxQual" = list(arg_list[["vsearchMaxQual"]], c("integer")), #strictly positive integer -- checked in dimsum__validate_input
    "vsearchMaxee" = list(arg_list[["vsearchMaxee"]], c("double")), #strictly positive double -- checked in dimsum__validate_input
    "vsearchMinovlen" = list(arg_list[["vsearchMinovlen"]], c("integer")), #strictly positive integer -- checked in dimsum__validate_input
    "outputPath" = list(arg_list[["outputPath"]], c("character")), #directory exists -- checked in dimsum__validate_input
    "projectName" = list(arg_list[["projectName"]], c("character")), #character string -- checked in dimsum__validate_input
    "wildtypeSequence" = list(arg_list[["wildtypeSequence"]], c("character")), #AaGgCcTt character string -- checked in dimsum__validate_input
    "permittedSequences" = list(arg_list[["permittedSequences"]], c("character", "NULL")), #ACGTRYSWKMBDHVN character string -- checked in dimsum__validate_input
    "reverseComplement" = list(arg_list[["reverseComplement"]], c("logical")), #logical -- checked in dimsum__validate_input
    "sequenceType" = list(arg_list[["sequenceType"]], c("character")), #character string; either noncoding/coding/auto -- checked in dimsum__validate_input
    "mutagenesisType" = list(arg_list[["mutagenesisType"]], c("character")), #character string; either random/codon -- checked in dimsum__validate_input
    "transLibrary" = list(arg_list[["transLibrary"]], c("logical")), #logical -- checked in dimsum__validate_input
    "transLibraryReverseComplement" = list(arg_list[["transLibraryReverseComplement"]], c("logical")), #logical -- checked in dimsum__validate_input
    "bayesianDoubleFitness" = list(arg_list[["bayesianDoubleFitness"]], c("logical")), #logical -- checked in dimsum__validate_input
    "bayesianDoubleFitnessLamD" = list(arg_list[["bayesianDoubleFitnessLamD"]], c("double")), #strictly positive double -- checked in dimsum__validate_input
    "fitnessMinInputCountAll" = list(arg_list[["fitnessMinInputCountAll"]], c("character")), #character string; positive integer (zero inclusive) -- checked in dimsum__validate_input
    "fitnessMinInputCountAny" = list(arg_list[["fitnessMinInputCountAny"]], c("character")), #character string; positive integer (zero inclusive) -- checked in dimsum__validate_input
    "fitnessMinOutputCountAll" = list(arg_list[["fitnessMinOutputCountAll"]], c("character")), #character string; positive integer (zero inclusive) -- checked in dimsum__validate_input
    "fitnessMinOutputCountAny" = list(arg_list[["fitnessMinOutputCountAny"]], c("character")), #character string; positive integer (zero inclusive) -- checked in dimsum__validate_input
    "fitnessHighConfidenceCount" = list(arg_list[["fitnessHighConfidenceCount"]], c("integer")), #positive integer (zero inclusive) -- checked in dimsum__validate_input
    "fitnessDoubleHighConfidenceCount" = list(arg_list[["fitnessDoubleHighConfidenceCount"]], c("integer")), #positive integer (zero inclusive) -- checked in dimsum__validate_input
    "fitnessNormalise" = list(arg_list[["fitnessNormalise"]], c("logical")), #logical -- checked in dimsum__validate_input
    "fitnessErrorModel" = list(arg_list[["fitnessErrorModel"]], c("logical")), #logical -- checked in dimsum__validate_input
    "fitnessDropoutPseudocount" = list(arg_list[["fitnessDropoutPseudocount"]], c("integer")), #positive integer (zero inclusive) -- checked in dimsum__validate_input
    "indels" = list(arg_list[["indels"]], c("character")), #character string; either all/none or a comma-separated list of sequence lengths -- checked in dimsum__validate_input
    "maxSubstitutions" = list(arg_list[["maxSubstitutions"]], c("integer")), #positive integer (greater than 1) -- checked in dimsum__validate_input
    "mixedSubstitutions" = list(arg_list[["mixedSubstitutions"]], c("logical")), #logical -- checked in dimsum__validate_input
    "retainIntermediateFiles" = list(arg_list[["retainIntermediateFiles"]], c("logical")), #logical -- checked in dimsum__validate_input
    "splitChunkSize" = list(arg_list[["splitChunkSize"]], c("double")), #strictly positive double -- checked in dimsum__validate_input
    "retainedReplicates" = list(arg_list[["retainedReplicates"]], c("character")), #comma-separated list of integers or "all" -- checked in dimsum__validate_input
    "startStage" = list(arg_list[["startStage"]], c("integer")), #strictly positive integer -- checked in dimsum__validate_input
    "stopStage" = list(arg_list[["stopStage"]], c("integer")), #positive integer (zero inclusive) -- checked in dimsum__validate_input
    "numCores" = list(arg_list[["numCores"]], c("integer")) #strictly positive integer -- checked in dimsum__validate_input
    )

  #Validate input
  exp_metadata <- dimsum__validate_input(dimsum_arg_list)

  #Get and check barcode design (if provided)
  exp_metadata[["barcode_design"]] <- dimsum__get_barcode_design(exp_metadata)

  #Get and check experiment design
  exp_metadata[["exp_design"]] <- dimsum__get_experiment_design(exp_metadata)

  #Check count file
  dimsum__check_countfile(exp_metadata)

  #Check barcode identity file
  dimsum__check_barcodeidentityfile(exp_metadata)

  #Get and check synonym sequences
  exp_metadata[["synonym_sequences"]] <- dimsum__get_synonymsequences(exp_metadata)

  #Create project directory (if doesn't already exist)
  exp_metadata[["project_path"]] <- file.path(exp_metadata[["outputPath"]], exp_metadata[["projectName"]])
  suppressWarnings(dir.create(exp_metadata[["project_path"]]))
  #Create temp directory (if doesn't already exist)
  exp_metadata[["tmp_path"]] <- file.path(exp_metadata[["project_path"]], "tmp")
  suppressWarnings(dir.create(exp_metadata[["tmp_path"]]))

  #Initialise reports directory (and erase if already exists)
  dimsum__create_dir(file.path(exp_metadata[["project_path"]], "reports"), execute = T)
  #Placeholder plot files (DiMSum writes them when it renders its preliminary report; the report
  #itself is rendered once, at the end)
  dimsum__render_report_placeholders(dimsum_meta = exp_metadata)

  ### Pipeline stages
  ###########################

  #Start pipeline tracking
  pipeline <- list()
  pipeline[['initial']] <- exp_metadata

  #Stages 0 to 3 (WRAP) return the metadata untouched when a count file is given
  pipeline[['3_tally']] <- pipeline[['initial']]

  ### Stage 4 (STEAM): Merge count files, process variant sequence and filter for desired substitution variants
  pipeline[['4_process']] <- dimsum_stage_merge(dimsum_meta = pipeline[['3_tally']], merge_outpath = pipeline[['3_tally']][["project_path"]],
    report_outpath = file.path(pipeline[['3_tally']][["project_path"]], "reports"))

  ### Stage 5 (STEAM): Calculate fitness and model error
  pipeline[['5_analyse']] <- dimsum_stage_counts_to_fitness(dimsum_meta = pipeline[['4_process']], fitness_outpath = pipeline[['4_process']][["project_path"]],
    report_outpath = file.path(pipeline[['4_process']][["project_path"]], "reports"))

  ### Save workspace
  ###########################

  dimsum__status_message("\n\n\nSaving workspace image...\n")
  dimsum__save_metadata(dimsum_meta = pipeline[['5_analyse']], n = 1)
  dimsum__render_report(dimsum_meta = pipeline[['5_analyse']], finalise = TRUE)
  dimsum__status_message("Done\n")
}

# The empty plot files dimsum__render_report(initialise = TRUE) creates, without rendering the report
dimsum__render_report_placeholders <- function(dimsum_meta) {
  all_files <- list(
    "dimsum__fastqc_report_pair1_fastqc.png",
    "dimsum__fastqc_report_pair2_fastqc.png",
    "dimsum__cutadapt_report_pair1.png",
    "dimsum__cutadapt_report_pair2.png",
    "dimsum__vsearch_report_paircounts.png",
    "dimsum__vsearch_report_mergedlength.png",
    "dimsum__merge_report_nucmutationpercentages.png",
    "dimsum__merge_report_nucmutationcounts.png",
    "dimsum__merge_report_aamutationpercentages.png",
    "dimsum__merge_report_aamutationcounts.png",
    "dimsum__diagnostics_report_count_hist_input_nt.png",
    "dimsum__diagnostics_report_count_hist_input_aa.png",
    "dimsum__diagnostics_report_scatterplotmatrix_all.png",
    "dimsum_stage_fitness_report_1_errormodel_fitness_inputcounts.png",
    "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_density.png",
    "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_density_norm.png",
    "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_scatter.png",
    "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_scatter_norm.png",
    "dimsum_stage_fitness_report_1_errormodel_repspec.png",
    "dimsum_stage_fitness_report_1_errormodel_leaveoneout_qqplot.png")
  lapply(all_files, function(x){write(NULL, file = file.path(dimsum_meta[["project_path"]], "reports", x))})
  invisible(NULL)
}

####################################################################################################
# Command line (DiMSum's options, plus the two assets)
####################################################################################################

status_message <- function(...){cat(sprintf(...), sep='', file=stdout())}

option_list <- list(
  optparse::make_option(opt_str=c("--runDemo"), type="logical", default=F, help = "Run the DiMSum demo (default:F)"),
  optparse::make_option(opt_str=c("--fastqFileDir", "-i"), default='./', help = "Path to directory containing input FASTQ files (default:'./' i.e. current working directory)"),
  optparse::make_option(opt_str=c("--fastqFileExtension", "-l"), default='.fastq', help = "FASTQ file extension (default:'.fastq')"),
  optparse::make_option(opt_str=c("--gzipped", "-g"), type="logical", default=T, help = "Are FASTQ files are gzipped? (default:T)"),
  optparse::make_option(opt_str=c("--stranded"), type="logical", default=T, help = "Is the library design stranded? (default:T)"),
  optparse::make_option(opt_str=c("--paired"), type="logical", default=T, help = "Is the library design paired-end? (default:T)"),
  optparse::make_option(opt_str=c("--barcodeDesignPath", "-b"), help = "Path to barcode design file (tab-separated plain text file with barcode design)"),
  optparse::make_option(opt_str=c("--barcodeErrorRate"), type="double", default=0.25, help = "Maximum allowed error rate for barcode to be matched (default:0.25)"),
  optparse::make_option(opt_str=c("--experimentDesignPath", "-e"), help = "Path to Experimental Design File (required if '--runDemo'=F)"),
  optparse::make_option(opt_str=c("--experimentDesignPairDuplicates"), type="logical", default=F, help = "Are multiple instances of FASTQ files in the Experimental Design File permitted? (default:F)"),
  optparse::make_option(opt_str=c("--barcodeIdentityPath"), help = "Path to Variant Identity File (tab-separated plain text file mapping barcodes to variants)"),
  optparse::make_option(opt_str=c("--countPath"), help = "Path to Variant Count File for analysis with STEAM only (tab-separated plain text file with sample counts for all variants)"),
  optparse::make_option(opt_str=c("--synonymSequencePath"), help = "Path to Synonym Sequences File with coding sequences for which synonymous variant fitness should be quantified (default: plain text file with one coding nucleotide sequence per line)"),
  optparse::make_option(opt_str=c("--cutadaptCut5First"), type="integer", help = "Remove fixed number of bases from start (5') of first (or only) read before constant region trimming (optional)"),
  optparse::make_option(opt_str=c("--cutadaptCut5Second"), type="integer", help = "Remove fixed number of bases from start (5') of second read in pair before constant region trimming (optional)"),
  optparse::make_option(opt_str=c("--cutadaptCut3First"), type="integer", help = "Remove fixed number of bases from end (3') of first (or only) read before constant region trimming (optional)"),
  optparse::make_option(opt_str=c("--cutadaptCut3Second"), type="integer", help = "Remove fixed number of bases from end (3') of second read in pair before constant region trimming (optional)"),
  optparse::make_option(opt_str=c("--cutadapt5First"), help = "Sequence of 5' constant region to be trimmed from first (or only) read (optional)"),
  optparse::make_option(opt_str=c("--cutadapt5Second"), help = "Sequence of 5' constant region to be trimmed from second read in pair (optional)"),
  optparse::make_option(opt_str=c("--cutadapt3First"), help = "Sequence of 3' constant region to be trimmed from first (or only) read (default: reverse complement of '--cutadapt5Second')"),
  optparse::make_option(opt_str=c("--cutadapt3Second"), help = "Sequence of 3' constant region to be trimmed from second read in pair (default: reverse complement of '--cutadapt5First')"),
  optparse::make_option(opt_str=c("--cutadaptMinLength", "-n"), type="integer", default=50, help = "Discard reads shorter than LENGTH after trimming (default:50)"),
  optparse::make_option(opt_str=c("--cutadaptErrorRate", "-a"), type="double", default=0.2, help = "Maximum allowed error rate for trimming constant regions (default:0.2)"),
  optparse::make_option(opt_str=c("--cutadaptOverlap"), type="integer", default=3, help = "Minimum overlap between read and constant region for trimming (default:3)"),
  optparse::make_option(opt_str=c("--vsearchMinQual", "-q"), type="integer", default=30, help = "Minimum Phred base quality score required to retain read or read pair (default:30)"),
  optparse::make_option(opt_str=c("--vsearchMaxQual"), type="integer", default=41, help = "Maximum Phred base quality score accepted when reading (and used when writing) FASTQ files; cannot be greater than 93 (default:41)"),
  optparse::make_option(opt_str=c("--vsearchMaxee", "-m"), type="double", default=0.5, help = "Maximum number of expected errors tolerated to retain read or read pair (default:0.5)"),
  optparse::make_option(opt_str=c("--vsearchMinovlen"), type="integer", default=10, help = "Discard read pair if the alignment length is shorter than this (default:10)"),
  optparse::make_option(opt_str=c("--outputPath", "-o"), default='./', help = "Path to directory to use for output files (default:'./' i.e. current working directory)"),
  optparse::make_option(opt_str=c("--projectName", "-p"), default='DiMSum_Project', help = "Project name and directory where results are to be saved (default:'DiMSum_Project')"),
  optparse::make_option(opt_str=c("--wildtypeSequence", "-w"), help = "Wild-type nucleotide sequence (A/C/G/T). Lower-case bases (a/c/g/t) indicate internal constant regions to be removed (required if '--runDemo'=F)"),
  optparse::make_option(opt_str=c("--permittedSequences"), help = "Nucleotide sequence of IUPAC ambiguity codes (A/C/G/T/R/Y/S/W/K/M/B/D/H/V/N) with length matching the number of mutated positions (i.e upper-case letters) in '--wildtypeSequence' (default:N i.e. any substitution mutation allowed)"),
  optparse::make_option(opt_str=c("--reverseComplement"), type="logical", default=F, help = "Reverse complement sequence (default:F)"),
  optparse::make_option(opt_str=c("--sequenceType", "-u"), default="auto", help = "Coding potential of sequence: either 'noncoding', 'coding' or 'auto'. If the specified wild-type nucleotide sequence ('--wildtypeSequence') has a valid translation without a premature STOP codon, it is assumed to be 'coding' (default:'auto')"),
  optparse::make_option(opt_str=c("--mutagenesisType"), default="random", help = "Whether mutagenesis was performed at the nucleotide or codon/amino acid level; either 'random' or 'codon' (default:'random')"),
  optparse::make_option(opt_str=c("--transLibrary", "-r"), type="logical", default=F, help = "Paired-end reads correspond to distinct molecules? (default:F)"),
  optparse::make_option(opt_str=c("--transLibraryReverseComplement"), type="logical", default=F, help = "Reverse complement second read in pair (default:F)"),
  optparse::make_option(opt_str=c("--bayesianDoubleFitness", "-y"), type="logical", default=F, help = "In development: improve double mutant fitness estimates using Bayesian framework (DISABLED: still in development)"),
  optparse::make_option(opt_str=c("--bayesianDoubleFitnessLamD", "-d"), type="double", default=0.025, help = "In development: Poisson distribution for score likelihood (default:0.025)"),
  optparse::make_option(opt_str=c("--fitnessMinInputCountAll"), default="0", help = "Minimum input read count (in all replicates) to be retained during fitness calculations (default:0)"),
  optparse::make_option(opt_str=c("--fitnessMinInputCountAny"), default="0", help = "Minimum input read count (in any replicate) to be retained during fitness calculations (default:0)"),
  optparse::make_option(opt_str=c("--fitnessMinOutputCountAll"), default="0", help = "Minimum output read count (in all replicates) to be retained during fitness calculations (default:0)"),
  optparse::make_option(opt_str=c("--fitnessMinOutputCountAny"), default="0", help = "Minimum output read count (in any replicates) to be retained during fitness calculations (default:0)"),
  optparse::make_option(opt_str=c("--fitnessHighConfidenceCount"), type="integer", default=10, help = "In development: minimum mean input read count for high confidence variants (default:10)"),
  optparse::make_option(opt_str=c("--fitnessDoubleHighConfidenceCount"), type="integer", default=50, help = "In development: minimum input replicate read count for doubles used to derive prior for Bayesian doubles correction (default:50)"),
  optparse::make_option(opt_str=c("--fitnessNormalise"), type="logical", default=T, help = "Normalise fitness values to minimise inter-replicate differences (default:T)"),
  optparse::make_option(opt_str=c("--fitnessErrorModel"), type="logical", default=T, help = "Fit fitness error model (default:T)"),
  optparse::make_option(opt_str=c("--fitnessDropoutPseudocount"), type="integer", default=0, help = "Pseudocount added to output replicates with dropout i.e. variants present in input but absent from output (default:0)"),
  optparse::make_option(opt_str=c("--indels"), default="none", help = "Indel variants to be retained: either 'all', 'none' or a comma-separated list of sequence lengths (default:'none')"),
  optparse::make_option(opt_str=c("--maxSubstitutions"), type="integer", default=2, help = "Maximum number of nucleotide or amino acid substitutions for coding or non-coding sequences respectively (default:2)"),
  optparse::make_option(opt_str=c("--mixedSubstitutions"), type="logical", default=F, help = "For coding sequences, are nonsynonymous variants with silent/synonymous substitutions in other codons allowed? (default:F)"),
  optparse::make_option(opt_str=c("--retainIntermediateFiles", "-f"), type="logical", default=F, help = "Should intermediate files be retained? Intermediate files can be many gigabytes, but are required to rerun DiMSum starting at intermediate pipeline stages (default:F)"),
  optparse::make_option(opt_str=c("--splitChunkSize", "-k"), type="double", default=3758096384, help = "Internal: FASTQ file split chunk size in bytes (default:3758096384)"),
  optparse::make_option(opt_str=c("--retainedReplicates"), default='all', help = "Comma-separated list of (integer) experiment replicates to retain or 'all' (default:'all')"),
  optparse::make_option(opt_str=c("--startStage", "-s"), type="integer", default=0, help = "(Re-)Start DiMSum at a specific pipeline stage (default:0)"),
  optparse::make_option(opt_str=c("--stopStage", "-t"), type="integer", default=5, help = "Stop DiMSum at a specific pipeline stage (default:5)"),
  optparse::make_option(opt_str=c("--numCores", "-c"), type="integer", default=1, help = "Number of available CPU cores. All pipeline stages make use of parallel computing to decrease runtime if multiple cores are available (default:1)"),
  optparse::make_option(opt_str=c("--reportRmd"), help = "Path to DiMSum's report.Rmd (required)"),
  optparse::make_option(opt_str=c("--logoPng"), help = "Path to DiMSum's Dumpling.png (required)")
)

arg_list <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
#Display arguments
status_message(paste("\n\n\n*******", "DiMSum STEAM command-line arguments", "*******\n\n\n"))
status_message(paste(formatDL(unlist(arg_list)), collapse = "\n"))

if(is.null(arg_list[["countPath"]])){
  stop("dimsum_fitness.R runs DiMSum's STEAM stages on a variant count table: '--countPath' is required (the WRAP stages are not included).", call. = FALSE)
}
if(is.null(arg_list[["reportRmd"]]) | is.null(arg_list[["logoPng"]])){
  stop("'--reportRmd' and '--logoPng' are required.", call. = FALSE)
}
DIMSUM_FITNESS_REPORT_RMD <- normalizePath(arg_list[["reportRmd"]], mustWork = TRUE)
DIMSUM_FITNESS_LOGO_PNG   <- normalizePath(arg_list[["logoPng"]], mustWork = TRUE)
arg_list[["reportRmd"]] <- NULL
arg_list[["logoPng"]] <- NULL
#DiMSum's STEAM entry point: stages 4 and 5
if(arg_list[["startStage"]] < 4){arg_list[["startStage"]] <- 4L}

dimsum_fitness(arg_list)
