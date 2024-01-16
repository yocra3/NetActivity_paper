#'#################################################################################
#'#################################################################################
#' Prepare GSE179848 dataset
#'#################################################################################
#'#################################################################################

## Load libraries
library(GEOquery)
library(tidyverse)
library(SummarizedExperiment)
library(org.Hs.eg.db)
library(DESeq2)
library(NetActivity)

## Download gene expression
geo <- getGEO("GSE179848")[[1]]
colData <- pData(geo)
rownames(colData) <- colData$description

counts <- read_csv("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179848/suppl/GSE179848%5Fraw%5Fcounts%5Fcell%5Flifespan%5FRNAseq%5Fdata.csv.gz")

## Merge phenotypes
phenotypes <- read_csv("data/Lifespan_Study_selected_data.csv")
colData$SampleID <- colData$description
phenotypes$SampleID <- paste0("Sample_", phenotypes$RNAseq_sampleID)

comb_phenotypes <- left_join(colData, phenotypes, by = "SampleID")
rownames(comb_phenotypes) <- comb_phenotypes$SampleID

## Prepare SE
count_matrix <- data.matrix(counts[, -1])
geneName <-  unlist(counts[, 1])
names(geneName) <- NULL
rownames(count_matrix) <- geneName
gse179848 <- SummarizedExperiment(count_matrix, colData = comb_phenotypes)
save(gse179848, file = "results/GSE179848/rawCounts_GSE179848.Rdata")

## Compute vst
counts <- assay(gse179848)
counts[] <- as.integer(counts)
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData(gse179848),
                              design = ~ 1)
gse179848_vst <- vst(dds, blind=FALSE)
save(gse179848_vst, file = "results/GSE179848/VST_GSE179848.Rdata")


## Convert gene names
rownames(gse179848_vst) <- mapIds(org.Hs.eg.db,
    keys = rownames(gse179848_vst),
    column = 'ENSEMBL',
    keytype = 'REFSEQ')

gse179848_ensembl <- gse179848_vst[!is.na(rownames(gse179848_vst)), ]
gse179848_ensembl <- gse179848_ensembl[!duplicated(rownames(gse179848_ensembl)), ]
save(gse179848_ensembl, file = "results/GSE179848/GSE179848_ENSEMBL.Rdata")

gse179848_prep <- prepareSummarizedExperiment(gse179848_ensembl, "gtex_gokegg")
gse179848_scores <- computeGeneSetScores(gse179848_prep, "gtex_gokegg")

save(gse179848_scores, file = "results/GSE179848/GSE179848_scores.Rdata")
