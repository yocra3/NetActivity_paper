#'#################################################################################
#'#################################################################################
#' Analyze BRCA
#'#################################################################################
#'#################################################################################

## Load libraries
library(limma)
library(DESeq2)
library(tidyverse)
library(cowplot)
library(GSVA)
library(HDF5Array)
library(hipathia)
library(org.Hs.eg.db)
library(parallel)
library(NetActivity)
library(NetActivityData)
library(GSVA)
library(hipathia)


load("data/tcga_gexp_combat.Rdata")

## Prepare BRCA data
brca.all <- gexp_tcga_combat[, gexp_tcga_combat$project_id == "TCGA-BRCA"]
brca <- brca.all[, !is.na(brca.all$paper_pathologic_stage) & brca.all$paper_pathologic_stage != "NA"]

brca$stage <- ifelse(brca$paper_pathologic_stage %in% c("Stage_I", "Stage_II"), "Stage I-II", "Stage III-IV")
brca_dds <- DESeqDataSet(brca, design = ~ paper_BRCA_Pathology  + age_at_index + race + stage )
vst_brca <- vst(brca_dds, blind=FALSE)
save(vst_brca, file = "results/TCGA_BRCA/vst_SE.Rdata")

## Raw DE analysis
dds <- DESeq(brca_dds)
res_brca <- results(dds)
res_brca$p.adj.bf <- p.adjust(res_brca$pvalue )
save(res_brca, file = "results/TCGA_BRCA/genes_results.Rdata")


## NetActivity
preproc_brca <- prepareSummarizedExperiment(vst_brca, "gtex_gokegg")
scores_brca <- computeGeneSetScores(preproc_brca, "gtex_gokegg")
save(scores_brca, file = "results/TCGA_BRCA/NetActivity_scores.Rdata")

mod <- model.matrix(~ stage + paper_BRCA_Pathology + age_at_index + race, colData(scores_brca))
lm.path <- lmFit(assay(scores_brca), mod) %>% eBayes()
tab.path_brca <- topTable(lm.path, coef = 2, n = Inf)
tab.path_brca$category <- rownames(tab.path_brca)
save(tab.path_brca,  file = "results/TCGA_BRCA/pathways_results.Rdata")

## GSVA
path.map <- read.table("results/GTEx_coding/go_kegg_filt2_gene_map.tsv", header = TRUE)
paths.vec <- rownames(scores_brca)
path_genes <- mclapply(paths.vec, function(x) subset(path.map, PathwayID == x & !is.na(Symbol))$Symbol, mc.cores = 10)
names(path_genes) <- paths.vec
brca_gsva <- gsva(brca, path_genes, min.sz=5, max.sz=500, kcdf = "Poisson")

lm.gsva <- lmFit(assay(brca_gsva), mod) %>% eBayes()
tab.gsva_brca <- topTable(lm.gsva, coef = 2, n = Inf)
tab.gsva_brca$category <- rownames(tab.gsva_brca)
save(brca_gsva, tab.gsva_brca, file = "results/TCGA_BRCA/GSVA_results.Rdata")

## hipathia
trans_data <- translate_data(vst_brca, "hsa")
exp_data <- normalize_data(trans_data)
hip_pathways <- load_pathways(species = "hsa")

hip.res_brca <- hipathia(exp_data, hip_pathways, decompose = FALSE, verbose = TRUE)
hip.brca_vals <- get_paths_data(hip.res_brca )
hip.comp_brca <- do_wilcoxon(hip.brca_vals, hip.brca_vals$stage, g1 = "Stage I-II", g2 = "Stage III-IV")

save(hip.comp_brca, hip.brca_vals, file = "results/TCGA_BRCA/hipathia.res.Rdata")
