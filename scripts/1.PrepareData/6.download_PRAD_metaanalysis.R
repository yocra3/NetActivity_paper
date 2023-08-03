#'#################################################################################
#'#################################################################################
#' Download data for a meta-analysis in PRAD
#'#################################################################################
#'#################################################################################

## Load libraries
library(limma)
library(oligo)
library(SummarizedExperiment)
library(DESeq2)
library(HDF5Array)
library(org.Hs.eg.db)
library(parallel)
library(GEOquery)
library(tidyverse)
library(huex10sttranscriptcluster.db)
library(NetActivity)
library(NetActivityData)
library(recount3)
options(timeout=10000)

data(gtex_gokegg)
sel_genes <- colnames(gtex_gokegg)

data(tcga_gokegg)
sel_genes_tcga <- colnames(tcga_gokegg)

## GSE46691 - array
### Download CEL files from web

### Preprocess CEL files
gse46691_celFiles <- list.celfiles("data/GSE46691", full.names=TRUE, listGzipped = TRUE)
gse46691_rawData <- read.celfiles(gse46691_celFiles)
gse46691_es <- rma(gse46691_rawData, target="core")

## Add annotation
annot_ensembl <- contents(huex10sttranscriptclusterENSEMBL)
annot_ensembl_filt <- annot_ensembl[rownames(gse46691_es)]
fData(gse46691_es) <- data.frame(id = names(annot_ensembl_filt))
fData(gse46691_es)$ENSEMBL <- annot_ensembl_filt
rownames(fData(gse46691_es)) <- fData(gse46691_es)$id

## Add phenotypes
colnames(gse46691_es) <- sapply(strsplit(colnames(gse46691_es), "_"), `[`, 1)
geo46691 <- getGEO("GSE46691")

pData(gse46691_es) <- pData(geo46691[[1]])[colnames(gse46691_es), ]
create_dir("results/preprocess/GSE46691")
save(gse46691_es, file = "results/preprocess/GSE46691/GSE46691_normalized.Rdata")

## Create NetActivity scores
gse46691_se <- SummarizedExperiment(exprs(gse46691_es), rowData = fData(gse46691_es), colData = pData(gse46691_es))
rowData(gse46691_se)$sel_gene <- sapply(rowData(gse46691_se)$ENSEMBL, function(x) {
  if (length(x) == 1 && is.na(x)){
    return(NA)
  } else if (all(!x %in% sel_genes)){
    return(NA)
  } else{
    x[x %in% sel_genes]
  }
})
gse46691_ensemb <- gse46691_se[!is.na(rowData(gse46691_se)$sel_gene), ]
gse46691_ensemb <- gse46691_ensemb[!duplicated(rowData(gse46691_ensemb)$sel_gene), ]
rownames(gse46691_ensemb) <- rowData(gse46691_ensemb)$sel_gene

save(gse46691_ensemb, file = "results/preprocess/GSE46691/GSE46691_ENSEMBL.Rdata")


gse46691_prep <- prepareSummarizedExperiment(gse46691_ensemb, "gtex_gokegg")
gse46691_scores <- computeGeneSetScores(gse46691_prep, "gtex_gokegg")
save(gse46691_scores, file = "results/preprocess/GSE46691/GSE46691_scores.Rdata")

### TCGA
gse46691_prep_tcga <- prepareSummarizedExperiment(gse46691_ensemb, "tcga_gokegg")
gse46691_scores_tcga <- computeGeneSetScores(gse46691_prep_tcga, "tcga_gokegg")
save(gse46691_scores_tcga, file = "results/preprocess/GSE46691/GSE46691_scores_tcga.Rdata")

## GSE141551 - array
gse141551 <- getGEO("GSE141551")[[1]]
gse141551 <- gse141551[[1]]
gse141551_se <- SummarizedExperiment(exprs(gse141551), rowData = fData(gse141551), colData = pData(gse141551))

rownames(gse141551_se) <- rowData(gse141551_se)$Symbol
gse141551_se <- gse141551_se[!is.na(rownames(gse141551_se)), ]
gse141551_se <- gse141551_se[!duplicated(rownames(gse141551_se)), ]

rownames(gse141551_se) <- mapIds(org.Hs.eg.db,
    keys = rownames(gse141551_se),
    column = 'ENSEMBL',
    keytype = 'SYMBOL')

gse141551_se <- gse141551_se[!is.na(rownames(gse141551_se)), ]
gse141551_se <- gse141551_se[!duplicated(rownames(gse141551_se)), ]

dir.create("results/preprocess/GSE141551")
save(gse141551_se, file = "results/preprocess/GSE141551/GSE141551_ENSEMBL.Rdata")

gse141551_prep <- prepareSummarizedExperiment(gse141551_se, "gtex_gokegg")
gse141551_scores <- computeGeneSetScores(gse141551_prep, "gtex_gokegg")

save(gse141551_scores, file = "results/preprocess/GSE141551/GSE141551_scores.Rdata")

### TCGA
gse141551_prep_tcga <- prepareSummarizedExperiment(gse141551_se, "tcga_gokegg")
gse141551_scores_tcga <- computeGeneSetScores(gse141551_prep_tcga, "tcga_gokegg")
save(gse141551_scores_tcga, file = "results/preprocess/GSE141551/GSE141551_scores_tcga.Rdata")


## GSE21034 - array
### Download CEL files from web

### Preprocess CEL files
gse21034_celFiles <- list.celfiles("data/GSE21034", full.names=TRUE, listGzipped = TRUE)
gse21034_rawData <- read.celfiles(gse21034_celFiles)
gse21034_es <- rma(gse21034_rawData, target="core")

## Add annotation
annot_ensembl <- contents(huex10sttranscriptclusterENSEMBL)
annot_ensembl_filt <- annot_ensembl[rownames(gse21034_es)]
fData(gse21034_es) <- data.frame(id = names(annot_ensembl_filt))
fData(gse21034_es)$ENSEMBL <- annot_ensembl_filt
rownames(fData(gse21034_es)) <- fData(gse21034_es)$id

## Add phenotypes
geo21034 <- getGEO("GSE21034")
colnames(gse21034_es) <- sapply(strsplit(colnames(gse21034_es), "_"), `[`, 1)

pData(gse21034_es) <- pData(geo21034)[colnames(gse21034_es), ]

## Select prostate cancer samples
gse21034_es <- gse21034_es[, !is.na(gse21034_es$`disease status:ch1`) & gse21034_es$`disease status:ch1` == "prostate cancer"]

create_dir("results/preprocess/GSE21034")
save(gse21034_es, file = "results/preprocess/GSE21034/GSE21034_normalized.Rdata")

## Create NetActivity scores
gse21034_se <- SummarizedExperiment(exprs(gse21034_es), rowData = fData(gse21034_es), colData = pData(gse21034_es))
rowData(gse21034_se)$sel_gene <- sapply(rowData(gse21034_se)$ENSEMBL, function(x) {
  if (length(x) == 1 && is.na(x)){
    return(NA)
  } else if (all(!x %in% sel_genes)){
    return(NA)
  } else{
    x[x %in% sel_genes]
  }
})
gse21034_ensemb <- gse21034_se[!is.na(rowData(gse21034_se)$sel_gene), ]
gse21034_ensemb <- gse21034_ensemb[!duplicated(rowData(gse21034_ensemb)$sel_gene), ]
rownames(gse21034_ensemb) <- rowData(gse21034_ensemb)$sel_gene

save(gse21034_ensemb, file = "results/preprocess/GSE21034/GSE21034_ENSEMBL.Rdata")


gse21034_prep <- prepareSummarizedExperiment(gse21034_ensemb, "gtex_gokegg")
gse21034_scores <- computeGeneSetScores(gse21034_prep, "gtex_gokegg")
save(gse21034_scores, file = "results/preprocess/GSE21034/GSE21034_scores.Rdata")

### TCGA
gse21034_prep_tcga <- prepareSummarizedExperiment(gse21034_ensemb, "tcga_gokegg")
gse21034_scores_tcga <- computeGeneSetScores(gse21034_prep_tcga, "tcga_gokegg")
save(gse21034_scores_tcga, file = "results/preprocess/GSE21034/GSE21034_scores_tcga.Rdata")

## GSE70768 - array
gse70768 <- getGEO("GSE70768")[[1]]
gse70768_se <- SummarizedExperiment(exprs(gse70768), rowData = fData(gse70768), colData = pData(gse70768))

rownames(gse70768_se) <- rowData(gse70768_se)$Symbol
gse70768_se <- gse70768_se[!is.na(rownames(gse70768_se)), ]
gse70768_se <- gse70768_se[!duplicated(rownames(gse70768_se)), ]

rownames(gse70768_se) <- mapIds(org.Hs.eg.db,
    keys = rownames(gse70768_se),
    column = 'ENSEMBL',
    keytype = 'SYMBOL')

gse70768_se <- gse70768_se[!is.na(rownames(gse70768_se)), ]
gse70768_se <- gse70768_se[!duplicated(rownames(gse70768_se)), ]

dir.create("results/preprocess/GSE70768")

## Filter non-tumor samples or samples without gleason
gse70768_se_filt <- gse70768_se[, !is.na(gse70768_se$`tumour gleason:ch1`) & gse70768_se$`tumour gleason:ch1` != "N/A" ]
save(gse70768_se_filt, file = "results/preprocess/GSE70768/GSE70768_ENSEMBL.Rdata")


gse70768_prep <- prepareSummarizedExperiment(gse70768_se_filt[rowMeans(is.na(assay(gse70768_se_filt))) == 0, ], "gtex_gokegg")
gse70768_scores <- computeGeneSetScores(gse70768_prep, "gtex_gokegg")

save(gse70768_scores, file = "results/preprocess/GSE70768/GSE70768_scores.Rdata")

### TCGA
gse70768_prep_tcga <- prepareSummarizedExperiment(gse70768_se_filt, "tcga_gokegg")
gse70768_scores_tcga <- computeGeneSetScores(gse70768_prep_tcga, "tcga_gokegg")
save(gse70768_scores_tcga, file = "results/preprocess/GSE70768/GSE70768_scores_tcga.Rdata")

## GSE70769 - array
gse70769 <- getGEO("GSE70769")[[1]]
gse70769_se <- SummarizedExperiment(exprs(gse70769), rowData = fData(gse70769), colData = pData(gse70769))

rownames(gse70769_se) <- rowData(gse70769_se)$Symbol
gse70769_se <- gse70769_se[!is.na(rownames(gse70769_se)), ]
gse70769_se <- gse70769_se[!duplicated(rownames(gse70769_se)), ]

rownames(gse70769_se) <- mapIds(org.Hs.eg.db,
    keys = rownames(gse70769_se),
    column = 'ENSEMBL',
    keytype = 'SYMBOL')

gse70769_se <- gse70769_se[!is.na(rownames(gse70769_se)), ]
gse70769_se <- gse70769_se[!duplicated(rownames(gse70769_se)), ]

dir.create("results/preprocess/GSE70769")

## Filter non-tumor samples or samples without gleason
gse70769_se_filt <- gse70769_se[, !is.na(gse70769_se$`tumour gleason:ch1`) & gse70769_se$`tumour gleason:ch1` != "unknown" ]
save(gse70769_se_filt, file = "results/preprocess/GSE70769/GSE70769_ENSEMBL.Rdata")

gse70769_prep <- prepareSummarizedExperiment(gse70769_se_filt, "gtex_gokegg")
gse70769_scores <- computeGeneSetScores(gse70769_prep, "gtex_gokegg")

save(gse70769_scores, file = "results/preprocess/GSE70769/GSE70769_scores.Rdata")

### TCGA
gse70769_prep_tcga <- prepareSummarizedExperiment(gse70769_se_filt, "tcga_gokegg")
gse70769_scores_tcga <- computeGeneSetScores(gse70769_prep_tcga, "tcga_gokegg")
save(gse70769_scores_tcga, file = "results/preprocess/GSE70769/GSE70769_scores_tcga.Rdata")

## GSE183019 - RNAseq
gse183019  <- getGEO("GSE183019")[[1]]
gse183019$id <- sapply(strsplit(gse183019$source_name_ch1, " "), `[`, 7)
colnames(gse183019) <- gse183019$id

gse183019_counts <- read_delim("data/GSE183019/GSE183019_processed_counts.txt.gz", delim = "\t")
gse183019_counts_mat <- data.matrix(gse183019_counts[, -c(1:2)])
mode(gse183019_counts_mat) <- "integer"
rownames(gse183019_counts_mat) <- gse183019_counts$GENE

gse183019_se <- SummarizedExperiment(gse183019_counts_mat, colData = pData(gse183019)[colnames(gse183019_counts_mat), ])
dir.create("results/preprocess/GSE183019")
save(gse183019_se, file = "results/preprocess/GSE183019/GSE183019_counts.Rdata")

gse183019_dds <- DESeqDataSet(gse183019_se, design = ~ 1)
gse183019_vst <- varianceStabilizingTransformation(gse183019_dds)


rownames(gse183019_vst) <- mapIds(org.Hs.eg.db,
    keys = rownames(gse183019_vst),
    column = 'ENSEMBL',
    keytype = 'SYMBOL')

gse183019_vst <- gse183019_vst[!is.na(rownames(gse183019_vst)), ]
gse183019_vst <- gse183019_vst[!duplicated(rownames(gse183019_vst)), ]

gse183019_prep <- prepareSummarizedExperiment(gse183019_vst, "gtex_gokegg")
gse183019_scores <- computeGeneSetScores(gse183019_prep, "gtex_gokegg")

save(gse183019_scores, file = "results/preprocess/GSE183019/GSE183019_scores.Rdata")

### TCGA
gse183019_prep_tcga <- prepareSummarizedExperiment(gse183019_vst, "tcga_gokegg")
gse183019_scores_tcga <- computeGeneSetScores(gse183019_prep_tcga, "tcga_gokegg")
save(gse183019_scores_tcga, file = "results/preprocess/GSE183019/GSE183019_scores_tcga.Rdata")



## GSE201284 - RNAseq
gse201284  <- getGEO("GSE201284")[[1]]
colnames(gse201284) <- gse201284$description.3

gse201284_counts <- read_delim("data/GSE201284/GSE201284_processed_counts.txt.gz", delim = "\t")
gse201284_counts_mat <- data.matrix(gse201284_counts[, -1])
mode(gse201284_counts_mat) <- "integer"
rownames(gse201284_counts_mat) <- gse201284_counts$gene_id

gse201284_se <- SummarizedExperiment(gse201284_counts_mat, colData = pData(gse201284)[colnames(gse201284_counts_mat), ])
dir.create("results/preprocess/GSE201284")
save(gse201284_se, file = "results/preprocess/GSE201284/GSE201284_counts.Rdata")

gse201284_dds <- DESeqDataSet(gse201284_se, design = ~ 1)
gse201284_vst <- varianceStabilizingTransformation(gse201284_dds)


rownames(gse201284_vst) <- mapIds(org.Hs.eg.db,
    keys = rownames(gse201284_vst),
    column = 'ENSEMBL',
    keytype = 'SYMBOL')

gse201284_vst <- gse201284_vst[!is.na(rownames(gse201284_vst)), ]
gse201284_vst <- gse201284_vst[!duplicated(rownames(gse201284_vst)), ]

gse201284_vst_filt <- gse201284_vst[, gse201284_vst$`sample type:ch1` == "TUMOR" ]
gse201284_prep <- prepareSummarizedExperiment(gse201284_vst_filt, "gtex_gokegg")
gse201284_scores <- computeGeneSetScores(gse201284_prep, "gtex_gokegg")

save(gse201284_scores, file = "results/preprocess/GSE201284/GSE201284_scores.Rdata")

### TCGA
gse201284_prep_tcga <- prepareSummarizedExperiment(gse201284_vst, "tcga_gokegg")
gse201284_scores_tcga <- computeGeneSetScores(gse201284_prep_tcga, "tcga_gokegg")
save(gse201284_scores_tcga, file = "results/preprocess/GSE201284/GSE201284_scores_tcga.Rdata")
