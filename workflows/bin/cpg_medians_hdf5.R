#! /usr/local/bin/Rscript
#'#################################################################################
#'#################################################################################
#' Remove probes with methylation range < 0.1
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
h5n <- args[1]
setPrefix <- gsub("assays.h5", "", h5n)

## Load libraries
library(DelayedMatrixStats)
library(HDF5Array)
library(SummarizedExperiment)

SE <- loadHDF5SummarizedExperiment(dir = "./", prefix = setPrefix)

medians <- DelayedMatrixStats::rowMedians(assay(SE), na.rm = TRUE)
names(medians) <- rownames(SE)
save(medians, file = "cpg_medians.Rdata")