#'#################################################################################
#'#################################################################################
# Download TCGA gene expression data
#'#################################################################################
#'#################################################################################

# Load libraries ####
library(GEOquery)
library(TCGAbiolinks)
library(HDF5Array)
library(SummarizedExperiment)

# TCGA ####
## Expression data ####
querye <- GDCquery(project = projects,
                  data.category = "Transcriptome Profiling",
                  legacy = FALSE,
                  workflow.type = "HTSeq - Counts"
)
GDCdownload(querye, method = "api", files.per.chunk = 10)
gexp_tcga <- GDCprepare(querye)
save(gexp_tcga, file = "data/tcga_gexp.Rdata")
