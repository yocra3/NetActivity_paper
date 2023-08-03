#'#################################################################################
#'#################################################################################
#' Analyze PRAD
#'#################################################################################
#'#################################################################################

## Load libraries
library(limma)
library(DESeq2)
library(tidyverse)
library(cowplot)
library(GSVA)
library(e1071)
library(HDF5Array)
library(hipathia)
library(org.Hs.eg.db)
library(parallel)
library(rjson)
library(NetActivity)


load("results/GSE169038/allGenes.se.RData")
se$race <- ifelse(grepl( "White", se$characteristics_ch1.4), "EUR", "AFR")
se$decipher <- factor(gsub("decipher risk group: ", "", se$characteristics_ch1.3), levels = c("Lower", "Average", "Higher"))
se$primary <- gsub("primary gleason: ", "", se$characteristics_ch1.1)
se$secondary <- gsub("secondary gleason: ", "", se$characteristics_ch1.2)
se$primary <- as.numeric(ifelse(se$primary == "--", 1, se$primary))
se$secondary <- as.numeric(ifelse(se$secondary == "--", 1, se$secondary))
se$gleason_cat <- paste(se$primary, se$secondary, sep = "-")
se$gleason <- ifelse(se$primary == 5 |  se$secondary == 5 | se$gleason_cat == "4+4", "High", "Low")

## Subset samples with gleason < 3
se.filt <- se[, !(se$primary == 1 |  se$secondary == 1)]
save(se.filt, file = "results/GSE169038/SE_filt.Rdata")

## DE genes
mod <- model.matrix(~  gleason + race + decipher, colData(se.filt))
lm.genes <- lmFit(assay(se.filt), mod) %>% eBayes()
tab.genes_geo <- topTable(lm.genes, coef = 2, n = Inf)

tab.genes_geo$gene <- rowData(se)[as.character(rownames(tab.genes_geo )), "gene"]

save(tab.genes_geo, file = "results/GSE169038/de_genes_results.Rdata")


## NetActivity
gse169038_se <- loadHDF5SummarizedExperiment("results/GSE169038/", prefix = "network_genes")
assay(gse169038_se) <- data.matrix(assay(gse169038_se))
gse169038_se$race <- ifelse(grepl( "White", gse169038_se$characteristics_ch1.4), "EUR", "AFR")
gse169038_se$decipher <- factor(gsub("decipher risk group: ", "", gse169038_se$characteristics_ch1.3), levels = c("Lower", "Average", "Higher"))
gse169038_se$primary <- gsub("primary gleason: ", "", gse169038_se$characteristics_ch1.1)
gse169038_se$secondary <- gsub("secondary gleason: ", "", gse169038_se$characteristics_ch1.2)
gse169038_se$primary <- as.numeric(ifelse(gse169038_se$primary == "--", 1, gse169038_se$primary))
gse169038_se$secondary <- as.numeric(ifelse(gse169038_se$secondary == "--", 1, gse169038_se$secondary))
gse169038_se$gleason_cat <- paste(gse169038_se$primary, gse169038_se$secondary, sep = "-")
gse169038_se$gleason <- factor(ifelse(gse169038_se$primary == 5 |  gse169038_se$secondary == 5 | gse169038_se$gleason_cat == "4+4", "High", "Low"), levels = c("Low", "High"))

gse169038_se_filt <- gse169038_se[, !(gse169038_se$primary == 1 |  gse169038_se$secondary == 1)]
gse169038_prep <- prepareSummarizedExperiment(gse169038_se_filt, "gtex_gokegg")
gse169038_scores <- computeGeneSetScores(gse169038_prep, "gtex_gokegg")
save(gse169038_scores, file = "results/GSE169038/NetActivity_scores.Rdata")

mod_gse169038 <- model.matrix(~ gleason + race + decipher, colData(gse169038_scores))
lm.paths <- lmFit(assay(gse169038_scores), mod_gse169038) %>% eBayes()
tab.paths_geo <- topTable(lm.paths, coef = 2, n = Inf)
tab.paths_geo$category <- rownames(tab.paths_geo)

tab.paths_geo$DE_prop <- sapply( tab.paths_geo$category, function(cat) {
  genes <- subset(path.map, PathwayID  == cat )$Symbol
  mini_tab <- subset(tab.genes_geo, gene %in% genes)
  mean(mini_tab$adj.P.Val < 0.05, na.rm = TRUE)
})

cor(tab.paths_geo$DE_prop, -log10(tab.paths_geo$P.Value))
# [1] 0.3426017

cor(abs(tab.paths_geo$logFC), tab.paths_geo$DE_prop)
# [1] 0.3162817

cor(tab.paths_geo$DE_prop[tab.paths_geo$DE_prop > 0 ], abs(tab.paths_geo$logFC)[tab.paths_geo$DE_prop > 0 ], use = "complete")
# 0.3053439
save(tab.paths_geo, file = "results/GSE169038/pathways_results.Rdata")

## GSVA
path_genes <- mclapply(paths.vec, function(x) subset(path.map, PathwayID == x & !is.na(Symbol))$Symbol, mc.cores = 10)
names(path_genes) <- paths.vec
geo_gsva <- gsva(data.matrix(assay(se.tcga_genes_filt)), path_genes, min.sz=5, max.sz=500)

lm.gsva <- lmFit(geo_gsva, mod) %>% eBayes()
tab.gsva_geo <- topTable(lm.gsva, coef = 2, n = Inf)
tab.gsva_geo$category <- rownames(tab.gsva_geo)

save(tab.gsva_geo, geo_gsva, file = "results/GSE169038/GSVA_results.Rdata")

tab.gsva_geo_race <- topTable(lm.gsva, coef = 3, n = Inf)

tab.paths_geo_race <- topTable(lm.paths, coef = 3, n = Inf)


## hipathia
rownames(se) <- rowData(se)$gene
trans_data <- translate_data(se, "hsa")
exp_data <- normalize_data(trans_data)
hip_pathways <- load_pathways(species = "hsa")

hip.res_geo <- hipathia(exp_data, hip_pathways, decompose = FALSE, verbose = TRUE)
hip.geo_path <- get_paths_data(hip.res_geo )
hip.comp_geo <- do_wilcoxon(hip.geo_path, hip.geo_path$gleason, g1 = "High", g2 = "Low")

save(hip.comp_geo, hip.geo_path, file = "results/GSE169038/hipathia.res.Rdata")