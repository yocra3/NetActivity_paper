#! /usr/local/bin/Rscript
#'#################################################################################
#'#################################################################################
#' Remove probes with call rate < 90%
#' Remove non-CpG probes
#' Sort object by genomic coordinates
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

## Remove non-cg positions
SE <- SE[grep("cg", rownames(SE)), ]

## Filter CpGs with all missings
pNA <- rowMeans(is.na(assay(SE)))
SE <- SE[pNA < 0.8, ]

## sort SE by GenomicCoordinates
SE <- sort(SE)
saveHDF5SummarizedExperiment(SE, ".", prefix = paste0(prefix, "probesFilt_"))

