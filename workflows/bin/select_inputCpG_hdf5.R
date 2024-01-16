#! /usr/local/bin/Rscript
#'#################################################################################
#'#################################################################################
#' Select probes from input model
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
h5n <- args[1]
probes_path <- args[2]

setPrefix <- gsub("assays.h5", "", h5n)

## Load libraries
library(DelayedMatrixStats)
library(HDF5Array)
library(SummarizedExperiment)

SE <- loadHDF5SummarizedExperiment(dir = "./", prefix = setPrefix)
load(probes_path)

out_probes <- setdiff(names(medians), rownames(SE))
out <- DelayedArray(matrix(rep(medians[out_probes], ncol(SE)), 
                           nrow = length(out_probes), ncol = ncol(SE),
              dimnames = list(out_probes, colnames(SE))))
beta <- rbind(assay(SE), out)
beta <- beta[names(medians), ]

SE <-  SummarizedExperiment(beta, colData = colData(SE))
saveHDF5SummarizedExperiment(SE, "./", prefix = paste0(setPrefix, "inputProbes_"))
