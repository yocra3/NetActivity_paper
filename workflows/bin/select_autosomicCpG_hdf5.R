#! /usr/local/bin/Rscript
#'#################################################################################
#'#################################################################################
#' Select autosomic probes
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
h5n <- args[1]
prefix <- args[2]
setPrefix <- gsub("assays.h5", "", h5n)

## Load libraries
library(DelayedMatrixStats)
library(HDF5Array)
library(SummarizedExperiment)

SE <- loadHDF5SummarizedExperiment(dir = "./", prefix = setPrefix)

SE <- SE[seqnames(rowRanges(SE)) %in% c(1:22, paste0("chr", 1:22)), ]

saveHDF5SummarizedExperiment(SE, ".", prefix = paste0(prefix, "autosomicProbes_"))
