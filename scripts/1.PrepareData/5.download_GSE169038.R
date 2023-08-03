#'#################################################################################
#'#################################################################################
#' Process GSE169038 gene expression for deep learning
#'#################################################################################
#'#################################################################################


## Load libraries ####
library(GEOquery)
library(SummarizedExperiment)
library(tidyverse)
library(HDF5Array)

## Download GEO
geo <- getGEO("GSE169038")

se <- SummarizedExperiment(exprs(geo[[1]]), rowData = fData(geo[[1]]), colData = pData(geo[[1]]))
rowData(se)$gene <- as.character(sapply(rowData(se)$mrna_assignment, function(x) unique(str_extract( x, "ENSG\\d*"))))

save(se, file = "results/GSE169038/allGenes.se.RData")


## Subset to genes present in TCGA
genes <- as.character(read.table("./results/TCGA_gexp_combat_coding/input_genes.txt")$V1)

se.filt <- se[rowData(se)$gene %in% genes, ]
rownames(se.filt) <- rowData(se.filt)$gene

out_probes <- setdiff(genes, rownames(se.filt))
out <- matrix(0, ncol(se.filt),
              nrow = length(out_probes), ncol = ncol(se.filt),
              dimnames = list(out_probes, colnames(se.filt)))

new_assay <- rbind(assay(se.filt), out)
se.tcga_genes <- SummarizedExperiment(new_assay, colData = colData(se))
se.tcga_genes <- se.tcga_genes[genes, ]
saveHDF5SummarizedExperiment(se.tcga_genes, "results/GSE169038/", prefix = "network_genes")
