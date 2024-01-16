#'#################################################################################
#'#################################################################################
#' Process TCGA gene expression for deep learning
#'#################################################################################
#'#################################################################################


## Load libraries ####
library(TCGAbiolinks)
library(limma)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)
library(HDF5Array)
library(BiocParallel)
library(sva)
library(rtracklayer)

register(MulticoreParam(10))

## Load data
load("data/tcga_gexp.Rdata")
gexp_fold <- "results/TCGA_gexp_combat/"

### Correct sequencing center with ComBat
phenos <- get_IDs(gexp_tcga)


gexp_tcga$sample_type2 <- ifelse(gexp_tcga$sample_type == "Solid Tissue Normal", "normal", "tumor")

adj.mod <- model.matrix(~ project_id, colData(gexp_tcga))[, -14]## Avoid confounding
adj_counts <- ComBat_seq(assay(gexp_tcga), batch = phenos$center , group = gexp_tcga$sample_type2, covar_mod = adj.mod, full_mod=TRUE)

gexp_tcga_combat <- gexp_tcga
assay(gexp_tcga_combat) <- adj_counts
save(gexp_tcga_combat, file = "data/tcga_gexp_combat.Rdata")

## Load gencode annotation
annot <- readGFF("/home/SHARED/DATA/REFERENCES/GRCh37/GenesAnnotation/gencode.v33lift37.annotation.gtf.gz")
annot <- filter(annot, type == "gene") %>%
  as_tibble()

annot_cod <- subset(annot, gene_type == "protein_coding")   %>%
 mutate(ensembl = gsub("\\.[0-9]*_[0-9]*", "", gene_id , perl = TRUE))

adj_counts2 <- adj_counts
mode(adj_counts2) <- "integer"
adj_counts2[is.na(adj_counts2)] <- max(adj_counts2, na.rm = TRUE)
dds <- DESeqDataSetFromMatrix(countData = adj_counts2,
                              colData = colData(gexp_tcga),
                              design = ~ project_id)
## Select genes expressed in at least 20 samples
keep <- rowSums(counts(dds) > 0) >= 20
vst <- vst(dds[keep, ], blind=FALSE)

## Output vsd values
saveHDF5SummarizedExperiment(vst, "results/TCGA_gexp_combat/", prefix = "vsd_norm")

## Write gene names
genes <- rownames(vst)
write.table(genes, file =  paste0(gexp_fold, "input_genes.txt"), quote = FALSE,
            row.names = FALSE, col.names = FALSE)

## Create labels file
project <- gexp_tcga$project_id
project[gexp_tcga$sample_type == "Solid Tissue Normal"] <- "Normal"


write.table(project, file = paste0(gexp_fold, "individuals_labels.txt"), quote = FALSE,
            row.names = FALSE, col.names = FALSE)

## Select protein coding genes
annot <- readGFF("/home/SHARED/DATA/REFERENCES/GRCh37/GenesAnnotation/gencode.v33lift37.annotation.gtf.gz")
annot <- filter(annot, type == "gene") %>%
  as_tibble() %>%
     mutate(Symbol = gsub("_[0-9]*", "", gene_id , perl = TRUE),
          Symbol = gsub("\\.[0-9]*", "", Symbol , perl = TRUE))

cod_genes <- subset(annot, gene_type == "protein_coding")$Symbol
vst_coding <- vst[rownames(vst) %in% cod_genes, ]

gexp_fold_cod <- "results/TCGA_gexp_combat_coding/"

saveHDF5SummarizedExperiment(vst_coding, "results/TCGA_gexp_combat_coding/", prefix = "vsd_norm")

vst_coding$Group <- as.character(vst_coding$project_id)
saveHDF5SummarizedExperiment(vst_coding, "results/TCGA_gexp_combat_coding/", prefix = "vsd_norm_group")


## Write gene names
genes <- rownames(vst_coding)
write.table(genes, file =  paste0(gexp_fold_cod, "input_genes.txt"), quote = FALSE,
            row.names = FALSE, col.names = FALSE)


## Create labels file
project <- as.character(vst_coding$project_id)
project[vst_coding$sample_type == "Solid Tissue Normal"] <- "Normal"
write.table(project, file = paste0(gexp_fold_cod, "individuals_labels.txt"), quote = FALSE,
            row.names = FALSE, col.names = FALSE)

gexp_fold_cod2 <- "results/TCGA_gexp_coding_noPRAD/"

vst_prad <- vst_coding[, vst_coding$project_id == "TCGA-PRAD"]
saveHDF5SummarizedExperiment(vst_prad, "results/TCGA_gexp_coding_noPRAD/", prefix = "vsd_norm_prad")

vst_prad_tum <- vst_prad[, !is.na(vst_prad$paper_Reviewed_Gleason_category)]
saveHDF5SummarizedExperiment(vst_prad_tum, "results/TCGA_gexp_coding_noPRAD/", prefix = "vsd_norm_prad_tumor")

vst_prad_control <- vst_prad[, vst_prad$sample_type == "Solid Tissue Normal"]
saveHDF5SummarizedExperiment(vst_prad_control, "results/TCGA_gexp_coding_noPRAD/", prefix = "vsd_norm_prad_control")



genes <- rownames(vst_train)
write.table(genes, file =  paste0(gexp_fold_cod2, "input_genes.txt"), quote = FALSE,
            row.names = FALSE, col.names = FALSE)

# Get labels
project <- as.character(vst_train$project_id)
project[vst_train$sample_type == "Solid Tissue Normal"] <- "Normal"

write.table(project, file = paste0(gexp_fold_cod2, "individuals_labels.txt"), quote = FALSE,
            row.names = FALSE, col.names = FALSE)


## Create object with only control samples
vst_control <- vst_coding[, vst_coding$sample_type == "Solid Tissue Normal"]
project <- as.character(vst_control$project_id)
project[project %in% names(table(project)[table(project) < 5])] <- "Other"


saveHDF5SummarizedExperiment(vst_control, "results/TCGA_gexp_coding_control/", prefix = "vsd_norm_control")
write.table(project, file = "results/TCGA_gexp_coding_control/individuals_labels.txt", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
