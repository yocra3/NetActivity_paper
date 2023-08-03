#'#################################################################################
#'#################################################################################
#' Process GSE57945 gene expression for deep learning
#'#################################################################################
#'#################################################################################


## Load libraries ####
library(recount)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)
library(HDF5Array)
library(BiocParallel)

## Download the RangedSummarizedExperiment object at the gene level for study SRP049593
url <- download_study('SRP042228', outdir = 'data/SRP042228')
load(file.path('data/SRP042228', 'rse_gene.Rdata'))

proj <- "results/SRP042228/"
## Add phenotypes
rse_gene$sex <- gsub("Sex: ", "", sapply(rse_gene$characteristics, function(x) x[2]))
rse_gene$age <- as.numeric(gsub("age at diagnosis:", "", sapply(rse_gene$characteristics, function(x) x[3])))
rse_gene$age_group <- ifelse(sapply(rse_gene$characteristics, function(x) x[4]) == "paris age: A1a", "Children", "Adolescent" )
rse_gene$diagnosis <- gsub("diagnosis: ", "", sapply(rse_gene$characteristics, function(x) x[5]))
rse_gene$diagnosis <- ifelse(rse_gene$diagnosis %in% c("Not IBD", "not IBD"), "Control", rse_gene$diagnosis)

rse_gene$diagnosis2 <- gsub("l2 type: ", "", sapply(rse_gene$characteristics, function(x) x[6]))
rse_gene$diagnosis2 <- ifelse(rse_gene$diagnosis2 %in% c("Not IBD", "not IBD"), "Control", rse_gene$diagnosis2)

rse_gene$histology <- gsub("histopathology: ", "", sapply(rse_gene$characteristics, function(x) x[7]))
rse_gene$ulcer <- sapply(rse_gene$characteristics, function(x) x[8])

save(rse_gene, file = paste0(proj, "RSE_phenotypes.Rdata"))


## Compute vst
dds <- DESeqDataSetFromMatrix(countData = assay(rse_gene),
                              colData = colData(rse_gene),
                              design = ~ 1)
rownames(dds) <- gsub("\\.[0-9]*", "", rownames(dds), perl = TRUE)
vsd <- vst(dds, blind=FALSE)
saveHDF5SummarizedExperiment(vsd, proj, prefix = "vsd_norm")


## Filter TCGA probes
sel.genes <- read.table("results/TCGA_gexp/input_genes.txt")

vsd.filt <- vsd[as.character(sel.genes$V1), ]
saveHDF5SummarizedExperiment(vsd.filt, proj, prefix = "vsd_norm_TCGAgenes_")


## Filter TCGA probes in go
go.genes <- read.table("results/TCGA_gexp_go/input_genes.txt")

vsd.filt <- vsd[as.character(go.genes$V1), ]
saveHDF5SummarizedExperiment(vsd.filt, proj, prefix = "vsd_norm_GOgenes_")

## Filter TCGA probes in go an autosomics
auto.go.genes <- read.table("results/TCGA_gexp_go/input_genes_autosomics.txt")
vsd.filt2 <- vsd.filt[as.character(auto.go.genes$V1), ]
saveHDF5SummarizedExperiment(vsd.filt2, proj, prefix = "vsd_norm_autosom_GOgenes_")


## Filter TCGA coding probes
sel.genes2 <- read.table("results/TCGA_gexp_combat_coding/input_genes.txt")

vsd.filt <- vsd[as.character(sel.genes2$V1), ]
saveHDF5SummarizedExperiment(vsd.filt, proj, prefix = "vsd_norm_TCGA_codingGenes_")
