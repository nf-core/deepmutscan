#!/usr/bin/env Rscript

## mutscan (Soneson et al., Genome Biology 2023) fitness estimation for nf-core/deepmutscan
## 19.03.2026

## 0. Libraries ##
##################

suppressPackageStartupMessages({
  library(mutscan)
  library(Biostrings)
})

## --- Helper functions ---

# calculate nt hamming distances from the specified WT
nbrMutBases <- function(x, wt.seq) {
  x <- cbind(x, "nbrMutBases" = rep(NA, nrow(x)))
  for (i in 1:nrow(x)){
    tmp.wt <- strsplit(as.character(wt.seq), "")[[1]]
    tmp.mut <- strsplit(as.character(x\$sequence[i]), "")[[1]]
    if(length(which(tmp.mut != tmp.wt)) == 0){
      x\$nbrMutBases[i] <- 0
      rm(tmp.mut, tmp.wt)
      next
    }else{
      x\$nbrMutBases[i] <- length(which(tmp.mut != tmp.wt))
      rm(tmp.mut, tmp.wt)
      next
    }
  }
  x\$nbrMutBases <- as.character(x\$nbrMutBases)
  return(x)
}

# calculate codon hamming distances from the specified WT
nbrMutCodons <- function(x, wt.seq) {
  x <- cbind(x, "nbrMutCodons" = rep(NA, nrow(x)))
  for (i in 1:nrow(x)){
    tmp.wt <- strsplit(as.character(wt.seq), "")[[1]]
    tmp.mut <- strsplit(as.character(x\$sequence[i]), "")[[1]]
    if(length(which(tmp.mut != tmp.wt)) == 0){
      x\$nbrMutCodons[i] <- 0
      rm(tmp.mut, tmp.wt)
      next
    }else if(length(which(tmp.mut != tmp.wt)) == 1){
      x\$nbrMutCodons[i] <- 1
      rm(tmp.mut, tmp.wt)
    }else if(length(which(tmp.mut != tmp.wt)) > 1){
      mut.pos <- which(tmp.mut != tmp.wt)
      x\$nbrMutCodons[i] <- length(unique(ceiling(mut.pos / 3)))
      rm(tmp.mut, tmp.wt)
      next
    }
  }
  x\$nbrMutCodons <- as.character(x\$nbrMutCodons)
  return(x)
}

# translate sequences and add aa_seq
sequenceAA <- function(x) {
  x <- cbind(x, "sequenceAA" = as.character(translate(DNAStringSet(x\$sequence))))
  return(x)
}

# calculate AA hamming distances from the WT
nbrMutAAs <- function(x, wt.seq) {
  x <- cbind(x, "nbrMutAAs" = rep(NA, nrow(x)))
  for (i in 1:nrow(x)){
    tmp.wt <- strsplit(as.character(translate(wt.seq)), "")[[1]]
    tmp.mut <- strsplit(as.character(x\$sequenceAA[i]), "")[[1]]
    if(length(which(tmp.mut != tmp.wt)) == 0){
      x\$nbrMutAAs[i] <- 0
      rm(tmp.mut, tmp.wt)
      next
    }else{
      x\$nbrMutAAs[i] <- length(which(tmp.mut != tmp.wt))
      rm(tmp.mut, tmp.wt)
      next
    }
  }
  x\$nbrMutAAs <- as.character(x\$nbrMutAAs)
  return(x)
}

# mutant base labels
mutantNameBase <- function(x, wt.seq){
  x <- cbind(x, "mutantNameBase" = rep(NA, nrow(x)))
  for(i in 1:nrow(x)){
    tmp.wt <- strsplit(as.character(wt.seq), "")[[1]]
    tmp.mut <- strsplit(as.character(x\$sequence[i]), "")[[1]]
    if(length(which(tmp.mut != tmp.wt)) == 0){
      x\$mutantNameBase[i] <- "f.0.WT"
      rm(tmp.mut, tmp.wt)
      next
    }else if(length(which(tmp.mut != tmp.wt)) > 0){
      mut.pos <- which(tmp.mut != tmp.wt)
      x\$mutantNameBase[i] <- paste0("f.", mut.pos, ".", tmp.mut[mut.pos], collapse = "_")
      rm(tmp.mut, tmp.wt, mut.pos)
      next
    }
  }
  return(x)
}

# mutant codon labels
mutantNameCodon <- function(x, wt.seq){
  x <- cbind(x, "mutantNameCodon" = rep(NA, nrow(x)))
  for(i in 1:nrow(x)){
    tmp.wt <- strsplit(as.character(wt.seq), "(?<=.{3})", perl = TRUE)[[1]]
    tmp.mut <- strsplit(as.character(x\$sequence[i]), "(?<=.{3})", perl = TRUE)[[1]]
    if(length(which(tmp.mut != tmp.wt)) == 0){
      x\$mutantNameCodon[i] <- "f.0.WT"
      rm(tmp.mut, tmp.wt)
      next
    }else if(length(which(tmp.mut != tmp.wt)) > 0){
      mut.pos <- which(tmp.mut != tmp.wt)
      x\$mutantNameCodon[i] <- paste0("f.", mut.pos, ".", tmp.mut[mut.pos], collapse = "_")
      rm(tmp.mut, tmp.wt, mut.pos)
      next
    }
  }
  return(x)
}

# name the mutations (AA level)
mutantNameAA <- function(x, wt.seq) {
  x <- cbind(x, "mutantNameAA" = rep(NA, nrow(x)))
  for (i in 1:nrow(x)){
    if(x\$nbrMutAAs[i] == 0){
      x\$mutantNameAA[i] <- "f.0.WT"
    }else{
      tmp.wt <- strsplit(as.character(translate(wt.seq)), "")[[1]]
      tmp.mut <- strsplit(as.character(x\$sequenceAA[i]), "")[[1]]
      tmp.pos <- which(tmp.mut != tmp.wt)
      x\$mutantNameAA[i] <- paste0("f.", tmp.pos, ".", tmp.mut[tmp.pos])
      rm(tmp.mut, tmp.wt,tmp.pos)
    }
  }
  return(x)
}

# categorise the mutation types
mutationTypes <- function(x){
  x <- cbind(x, "mutationTypes" = rep(NA, nrow(x)))
  x[grep("WT", x\$mutantNameAA), "mutationTypes"] <- "silent" ## silent/synonymous
  x[grep("[*]", x\$mutantNameAA), "mutationTypes"] <- "stop" ## stop
  x[-grep("[*]|WT", x\$mutantNameAA), "mutationTypes"] <- "nonsynonymous" ## nonsynonymous
  x
}

# mutant name
mutantName <- function(x, wt.seq){
  x <- cbind("mutantName" = rep(NA, nrow(x)), x)
  for (i in 1:nrow(x)){
    if(x\$nbrMutAAs[i] == 0){
      x\$mutantName[i] <- "WT"
    }else{
      tmp.wt <- strsplit(as.character(translate(wt.seq)), "")[[1]]
      tmp.mut <- strsplit(as.character(x\$sequenceAA[i]), "")[[1]]
      tmp.pos <- which(tmp.mut != tmp.wt)
      x\$mutantName[i] <- paste0(tmp.wt[tmp.pos], tmp.pos,tmp.mut[tmp.pos])
      rm(tmp.mut, tmp.wt,tmp.pos)
    }
  }
  return(x)
}

# mutscan 'summaryTable' building (from merged counts TSV file)
mutscan.summaryTable.from.counts <- function(x, sample, wt.seq){

  ### pre-process
  x <- x[,c(1, grep(sample, colnames(x)))]
  colnames(x) <- c("sequence", "nbrReads")
  x <- cbind(x, "maxNbrReads" = x[,"nbrReads"], "nbrUmis" = x[,"nbrReads"])

  ### nbrMutBases
  x <- nbrMutBases(x, wt.seq)

  ### nbrMutCodons
  x <- nbrMutCodons(x, wt.seq)

  ### sequenceAA
  x <- sequenceAA(x)

  ### nbrMutAAs
  x <- nbrMutAAs(x, wt.seq)

  ### varLengths
  x <- cbind(x, "varLengths" = as.character(nchar(x\$sequence)))

  ### mutantNameBase
  x <- mutantNameBase(x, wt.seq)

  ### mutantNameCodon
  x <- mutantNameCodon(x, wt.seq)

  ### mutantNameBaseHGVS (would ideally look like this: "f:c", "f:c.32_33delinsAC", etc.)
  x <- cbind(x, "mutantNameBaseHGVS" = x\$mutantNameCodon)

  ### mutantNameAA
  x <- mutantNameAA(x, wt.seq)

  ### mutantNameAAHGVS (would ideally look like this: "f:p" or "f:p.(Leu11His)")
  x <- cbind(x, "mutantNameAAHGVS" = x\$mutantNameAA)

  ### mutationTypes: silent, stop, nonsynonymous
  x <- mutationTypes(x)

  ### mutantName
  x <- mutantName(x, wt.seq)

  ### re-order columns
  x <- x[,c(1:7,9:16,8)]

  ### output
  return(x)
}

## --- Main function ---

#' Run mutscan fitness estimation with configurable I/O paths
#'
#' @param counts_path Path to counts_merged.tsv
#' @param design_path Path to experimentalDesign.tsv
#' @param wt_seq_path Path to synonymous_wt.txt (single line DNA sequence)
#' @param output_path Path to write fitness_estimation.tsv
#'
#' @return Invisibly returns the final data.frame; writes the output to output_path.
run_mutscan_fitness_estimation <- function(counts_path,
                                           design_path,
                                           wt_seq_path,
                                           output_path_edgeR,
                                           output_path_limma){

  ## 1. Import key files ##
  #########################

  merged.counts <- read.table(counts_path, sep = "\t", header = T, check.names = F)
  exp.design <- read.table(design_path, sep = "\t", header = T, check.names = F)
  wt.seq <- DNAString(as.character(read.table(wt_seq_path)))

  ## 2. Variant count matrix reformatting ##
  ##########################################

  var.tables <- vector(mode = "list", length = nrow(exp.design))
  names(var.tables) <- exp.design[,"sample_name"]
  var.tables <- lapply(var.tables, function(x){x <- vector(mode = "list", length = 4);
  names(x) <- c("summaryTable", "filterSummary", "errorStatistics", "parameters");
  return(x)})

  # mutscan 'summaryTable'
  for(i in 1:length(var.tables)){
    print(i)
    var.tables[[i]]\$summaryTable <- mutscan.summaryTable.from.counts(merged.counts, names(var.tables)[i], wt.seq)
  }

  # mutscan 'filterSummary' (fill with minimal decoy)
  var.tables <- lapply(var.tables, function(x){x\$filterSummary <- data.frame(NA); return(x)})

  # mutscan 'errorStatistics' (fill with minimal decoy)
  var.tables <- lapply(var.tables, function(x){x\$errorStatistics <- NA; return(x)})

  # mutscan 'parameters' (fill with minimal decoy)
  var.tables <- lapply(var.tables, function(x){x\$parameters\$mutNameDelimiter <- '.'; return(x)})

  ## 3. summarizeExperiment object ##
  ###################################

  # mutscan 'coldata' object (from experimental design TSV file)
  condition <- exp.design\$selection_id
  condition[which(condition == '0')] <- "input"
  condition[which(condition == '1')] <- "output"
  coldata <- as.data.frame(cbind("Name" = exp.design\$sample_name,
                                 "Condition" = condition,
                                 "Replicate" = exp.design\$experiment_replicate))
  class(coldata\$Replicate) <- "integer"

  # mutscan 'summarizeExperiment' object
  se <- summarizeExperiment(x = var.tables,
                            coldata = coldata,
                            countType = "reads")

  ## 4. logFC calculations ##
  ###########################

  # mutscan 'model.matrix' object to calculate logFC values
  model.design <- model.matrix(~ Replicate + Condition, data = se@colData)

  # edgeR
  logFC.edgeR <- calculateRelativeFC(se = se,
                                     design = model.design,
                                     coef = "Conditionoutput",
                                     WTrows = "WT",
                                     selAssay = "counts",
                                     pseudocount = 1,
                                     method = "edgeR")

  # limma
  logFC.limma <- calculateRelativeFC(se = se,
                                     design = model.design,
                                     coef = "Conditionoutput",
                                     WTrows = "WT",
                                     selAssay = "counts",
                                     pseudocount = 1,
                                     method = "limma")

  ## 5. QC plots ##
  #################

  # raw counts comparison
  pdf('mutscan_counts_corr.pdf', height = 9, width = 14)
  print(plotPairs(se, selAssay = "counts", addIdentityLine = TRUE))
  dev.off()

  # edgeR volcano plot
  pdf('mutscan_edgeR_volcano.pdf', height = 9, width = 14)
  print(plotVolcano(logFC.edgeR, pointSize = "large"))
  dev.off()

  # limma volcano plot
  pdf('mutscan_limma_volcano.pdf', height = 9, width = 14)
  print(plotVolcano(logFC.limma, pointSize = "large"))
  dev.off()

  ## 6. Data export ##
  ####################

  write.table(logFC.edgeR, output_path_edgeR,
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t", na = "")
  write.table(logFC.limma, output_path_limma,
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t", na = "")

  invisible(logFC.edgeR)
  invisible(logFC.limma)

}

#####
# run function
#####
run_mutscan_fitness_estimation(
  counts_path = "$counts_merged",
  design_path = "$exp_design",
  wt_seq_path = "$syn_wt_txt",
  output_path_edgeR = "fitness_estimation_mutscan_edgeR.tsv",
  output_path_limma = "fitness_estimation_mutscan_limma.tsv"
)

####
# create versions.yml
####
r_version <- strsplit(version[['version.string']], ' ')[[1]][3]
mutscan_version <- as.character(packageVersion("mutscan"))
Biostrings_version <- as.character(packageVersion("Biostrings"))

if (is.null(r_version)) r_version <- "unknown"
if (length(mutscan_version) == 0) mutscan_version <- "unknown"
if (length(Biostrings_version) == 0) Biostrings_version <- "unknown"

f <- file("versions.yml", "w")
writeLines(
  c(
    '"\${task.process}":',
    paste('    r-base:', r_version),
    paste('    r-mutscan:', mutscan_version),
    paste('    r-Biostrings:', Biostrings_version)
  ),
  f
)
close(f)
