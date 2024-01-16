#! /usr/local/bin/Rscript
#'#################################################################################
#'#################################################################################
#' Remove probes with methylation range < 0.1
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

quant <- rowQuantiles(assay(SE), probs = c(0.99, 0.01), na.rm = TRUE)
ranges <- apply(quant, 1, function(x) x[1] - x[2])
SE <- SE[ranges > 0.1, ]

saveHDF5SummarizedExperiment(SE, ".", prefix = paste0(prefix, "variantProbes_"))
