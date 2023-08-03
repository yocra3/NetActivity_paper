#'#################################################################################
#'#################################################################################
#' Run a meta-analysis in PRAD
#'#################################################################################
#'#################################################################################

## Load libraries
library(limma)
library(cowplot)
library(tidyverse)
library(SummarizedExperiment)
library(HDF5Array)
library(NetActivity)
library(cowplot)

createMETAL_tab <- function(fit){

  SE <- sqrt(fit$s2.post) * fit$stdev.unscaled %>% data.frame()
  res <- topTable(fit, coef = 2, n = Inf)
  res$SE <- SE[rownames(res), 2]
  res$Marker <- rownames(res)
  res$TestedAllele <- "C"
  res$OtherAllele <- "T"
  res$N <- nrow(fit$qr$qr)
  res
}
dir.create("results/PRAD_metanalysis/")

## GSE46691
load("results/preprocess/GSE46691/GSE46691_scores.Rdata")

gse46691_scores$gleason <- factor(ifelse(gse46691_scores$`gleason score:ch1` >= 8, "High", "Low"), levels = c("Low", "High"))
gse46691_scores$metastasis <- gse46691_scores$`metastatic event:ch1`

mod_gse46691 <- model.matrix(~ gleason + metastasis, colData(gse46691_scores))
lm.gse46691 <- lmFit(assay(gse46691_scores), mod_gse46691) %>% eBayes()
tab.gse46691 <- topTable(lm.gse46691, coef = 2, n = Inf)

metal.gse46691 <- createMETAL_tab(lm.gse46691)
write.table(metal.gse46691, file = "results/PRAD_metanalysis/GSE46691_metal.txt",
  sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

gse46691_pc <- prcomp(t(assay(gse46691_scores)), rank. = 4)
gse46691_pc_plot <- gse46691_pc$x %>%
  data.frame() %>%
  mutate(Gleason = gse46691_scores$gleason) %>%
  ggplot(aes(x = PC1, y = PC2, color = Gleason)) +
  geom_point() +
  theme_bw() +
  ggtitle("GSE46691") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


## GSE141551
load("results/preprocess/GSE141551/GSE141551_scores.Rdata")
gse141551_scores$gleason <- factor(ifelse(gse141551_scores$`gleason_group:ch1` == "8-10", "High", "Low"), levels = c("Low", "High"))
gse141551_scores$race <- gse141551_scores$`race:ch1`
gse141551_scores$age <- gse141551_scores$`age_group:ch1`
gse141551_scores$stage <- gse141551_scores$`stage_disease:ch1`

mod_gse141551 <- model.matrix(~ gleason + race + age + stage, colData(gse141551_scores))
lm.gse141551 <- lmFit(assay(gse141551_scores), mod_gse141551) %>% eBayes()
tab.gse141551 <- topTable(lm.gse141551, coef = 2, n = Inf)

metal.gse141551 <- createMETAL_tab(lm.gse141551)
write.table(metal.gse141551, file = "results/PRAD_metanalysis/GSE141551_metal.txt",
  sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### PC
gse141551_pc <- prcomp(t(assay(gse141551_scores)), rank. = 4)
gse141551_pc_plot <- gse141551_pc$x %>%
  data.frame() %>%
  mutate(Gleason = gse141551_scores$gleason) %>%
  ggplot(aes(x = PC1, y = PC2, color = Gleason)) +
  geom_point() +
  theme_bw() +
  ggtitle("GSE141551") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

## GSE21034 - array
load("results/preprocess/GSE21034/GSE21034_scores.Rdata")
gse21034_scores$gleason <- factor(ifelse(gse21034_scores$`biopsy_gleason_grade:ch1` >= 8, "High", "Low"), levels = c("Low", "High"))
gse21034_scores$metastasis <- gse21034_scores$`tumor type:ch1`

mod_gse21034 <- model.matrix(~ gleason + metastasis, colData(gse21034_scores))
lm.gse21034 <- lmFit(assay(gse21034_scores), mod_gse21034) %>% eBayes()
tab.gse21034 <- topTable(lm.gse21034, coef = 2, n = Inf)

metal.gse21034 <- createMETAL_tab(lm.gse21034)
write.table(metal.gse21034, file = "results/PRAD_metanalysis/GSE21034_metal.txt",
  sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### PC
gse21034_pc <- prcomp(t(assay(gse21034_scores)), rank. = 4)
gse21034_pc_plot <- gse21034_pc$x %>%
  data.frame() %>%
  mutate(Gleason = gse21034_scores$gleason) %>%
  ggplot(aes(x = PC1, y = PC2, color = Gleason)) +
  geom_point() +
  theme_bw() +
  ggtitle("GSE21034") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

## GSE70768 - array
load("results/preprocess/GSE70768/GSE70768_scores.Rdata")
gse70768_scores$gleason <- factor(ifelse(sapply(strsplit(gse70768_scores$`tumour gleason:ch1`, "="), `[`, 1) %>% as.numeric() >= 8, "High", "Low"), levels = c("Low", "High"))
gse70768_scores$age <- as.numeric(gse70768_scores$`age at diag:ch1`)
gse70768_scores$tum_prop <- as.numeric(gsub("%", "", gse70768_scores$`tumour %:ch1`))

mod_gse70768 <- model.matrix(~ gleason + age + tum_prop, colData(gse70768_scores))
lm.gse70768 <- lmFit(assay(gse70768_scores)[, rownames(mod_gse70768)], mod_gse70768) %>% eBayes()
tab.gse70768 <- topTable(lm.gse70768, coef = 2, n = Inf)

metal.gse70768 <- createMETAL_tab(lm.gse70768)
write.table(metal.gse70768, file = "results/PRAD_metanalysis/GSE70768_metal.txt",
  sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### PC
gse70768_pc <- prcomp(t(assay(gse70768_scores)), rank. = 4)
gse70768_pc_plot <- gse70768_pc$x %>%
  data.frame() %>%
  mutate(Gleason = gse70768_scores$gleason) %>%
  ggplot(aes(x = PC1, y = PC2, color = Gleason)) +
  geom_point() +
  theme_bw() +
  ggtitle("GSE70768") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

## GSE70769 - array
load("results/preprocess/GSE70769/GSE70769_scores.Rdata")
gse70769_scores$gleason <- factor(ifelse(sapply(strsplit(gse70769_scores$`tumour gleason:ch1`, "="), `[`, 1) %>% as.numeric() >= 8, "High", "Low"), levels = c("Low", "High"))

mod_gse70769 <- model.matrix(~ gleason, colData(gse70769_scores))
lm.gse70769 <- lmFit(assay(gse70769_scores), mod_gse70769) %>% eBayes()
tab.gse70769 <- topTable(lm.gse70769, coef = 2, n = Inf)

metal.gse70769 <- createMETAL_tab(lm.gse70769)
write.table(metal.gse70769, file = "results/PRAD_metanalysis/GSE70769_metal.txt",
  sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### PC
gse70769_pc <- prcomp(t(assay(gse70769_scores)), rank. = 4)
gse70769_pc_plot <- gse70769_pc$x %>%
  data.frame() %>%
  mutate(Gleason = gse70769_scores$gleason) %>%
  ggplot(aes(x = PC1, y = PC2, color = Gleason)) +
  geom_point() +
  theme_bw() +
  ggtitle("GSE70769") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))



## GSE183019 - RNAseq
load("results/preprocess/GSE183019/GSE183019_scores.Rdata")
gse183019_scores$gleason <- factor(ifelse(gse183019_scores$`gleason score:ch1` %in% c("4+4", "4+5", "5+4", "5+5 T4"), "High", "Low"), levels = c("Low", "High"))
gse183019_scores$age <- as.numeric(gse183019_scores$`age:ch1`)


mod_gse183019 <- model.matrix(~ gleason + age, colData(gse183019_scores))
lm.gse183019 <- lmFit(assay(gse183019_scores), mod_gse183019) %>% eBayes()
tab.gse183019 <- topTable(lm.gse183019, coef = 2, n = Inf)

metal.gse183019 <- createMETAL_tab(lm.gse183019)
write.table(metal.gse183019, file = "results/PRAD_metanalysis/GSE183019_metal.txt",
  sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

n.sv <- num.sv(assay(gse183019_scores), mod_gse183019, method="leek")
mod0_gse183019 <- model.matrix(~ age, colData(gse183019_scores))
svobj <- sva(assay(gse183019_scores), mod_gse183019, mod0_gse183019, n.sv=n.sv)
modSv <- cbind(mod_gse183019, svobj$sv)
lm.gse183019_sva <- lmFit(assay(gse183019_scores), modSv) %>% eBayes()
tab.gse183019_sva <- topTable(lm.gse183019_sva, coef = 2, n = Inf)


metal.gse183019_sva <- createMETAL_tab(lm.gse183019_sva)
write.table(metal.gse183019_sva, file = "results/PRAD_metanalysis/GSE183019_metal.txt",
  sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### PC
gse183019_pc <- prcomp(t(assay(gse183019_scores)), rank. = 4)
gse183019_pc_plot <- gse183019_pc$x %>%
  data.frame() %>%
  mutate(Gleason = gse183019_scores$gleason) %>%
  ggplot(aes(x = PC1, y = PC2, color = Gleason)) +
  geom_point() +
  theme_bw() +
  ggtitle("GSE183019") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


## PRAD - RNAseq
prad <- loadHDF5SummarizedExperiment("results/TCGA_gexp_coding_noPRAD/", prefix = "vsd_norm_prad_tumor")
prad$gleason <- factor(ifelse(prad$paper_Reviewed_Gleason_category == ">=8", "High", "Low"), levels = c("Low", "High"))
assay(prad) <- data.matrix(assay(prad))

prad_prep <- prepareSummarizedExperiment(prad, "gtex_gokegg")
prad_scores <- computeGeneSetScores(prad_prep, "gtex_gokegg")

mod_prad <- model.matrix(~ gleason + paper_Subtype + age_at_index + race, colData(prad))
lm.prad <- lmFit(assay(prad_scores), mod_prad) %>% eBayes()
tab.prad <- topTable(lm.prad, coef = 2, n = Inf)

metal.prad <- createMETAL_tab(lm.prad)
write.table(metal.prad, file = "results/PRAD_metanalysis/PRAD_metal.txt",
  sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### PC
prad_pc <- prcomp(t(assay(prad_scores)), rank. = 4)
prad_pc_plot <- prad_pc$x %>%
  data.frame() %>%
  mutate(Gleason = prad_scores$gleason) %>%
  ggplot(aes(x = PC1, y = PC2, color = Gleason)) +
  geom_point() +
  theme_bw() +
  ggtitle("PRAD") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


## GSE169038
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

mod_gse169038 <- model.matrix(~ gleason + race + decipher, colData(gse169038_scores))
lm.gse169038 <- lmFit(assay(gse169038_scores), mod_gse169038) %>% eBayes()
tab.gse169038 <- topTable(lm.gse169038, coef = 2, n = Inf)

metal.gse169038 <- createMETAL_tab(lm.gse169038)
write.table(metal.gse169038, file = "results/PRAD_metanalysis/GSE169038_metal.txt",
  sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### PC
gse169038_pc <- prcomp(t(assay(gse169038_scores)), rank. = 4)
gse169038_pc_plot <- gse169038_pc$x %>%
  data.frame() %>%
  mutate(Gleason = gse169038_scores$gleason) %>%
  ggplot(aes(x = PC1, y = PC2, color = Gleason)) +
  geom_point() +
  theme_bw() +
  ggtitle("GSE169038") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


## GSE201284 - RNAseq
load("results/preprocess/GSE201284/GSE201284_scores.Rdata")
gse201284_scores$gleason <- factor(ifelse(gse201284_scores$`patient gleason score:ch1` %in% c("3+5", "3+5 T4", "4+4", "4+5", "5+4", "5+5", "4+4 T5", "4+5 T3", "5+4 T3"), "High", "Low"), levels = c("Low", "High"))
gse201284_scores$age <- as.numeric(gse201284_scores$`age:ch1`)
gse201284_scores$tumor_prop <- as.numeric(gse201284_scores$`% tumor:ch1`)


mod_gse201284 <- model.matrix(~ gleason + age + tumor_prop, colData(gse201284_scores))
lm.gse201284 <- lmFit(assay(gse201284_scores)[, rownames(mod_gse201284)], mod_gse201284) %>% eBayes()
tab.gse201284 <- topTable(lm.gse201284, coef = 2, n = Inf)

metal.gse201284 <- createMETAL_tab(lm.gse201284)
write.table(metal.gse201284, file = "results/PRAD_metanalysis/GSE201284_metal.txt",
  sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### PC
gse201284_pc <- prcomp(t(assay(gse201284_scores)), rank. = 4)
gse201284_pc_plot <- gse201284_pc$x %>%
  data.frame() %>%
  mutate(Gleason = gse201284_scores$gleason) %>%
  ggplot(aes(x = PC1, y = PC2, color = Gleason)) +
  geom_point() +
  theme_bw() +
  ggtitle("GSE201284") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


## PCs of scores
png("figures/PRAD_metanalysis_pcs.png", width = 1000, height = 1000)
plot_grid(prad_pc_plot, gse201284_pc_plot, gse183019_pc_plot,
          gse70769_pc_plot, gse70768_pc_plot, gse141551_pc_plot,
          gse46691_pc_plot, gse21034_pc_plot, gse169038_pc_plot,
        nrow = 3, labels = LETTERS[1:9])
dev.off()

## Evaluation
prad_meta <- read_table("results/PRAD_metanalysis/PRAD_metanalysis1.tbl")
prad_meta$p_adj <- p.adjust(prad_meta$`P-value`)
dir_df <- strsplit(prad_meta$Direction, "") %>% Reduce(f = rbind) %>%
  data.frame()
colnames(dir_df) <- paste0("Dir_", c("GSE183019", "GSE70769", "GSE70768", "GSE46691",
  "GSE141551", "GSE21034", "PRAD", "GSE169038", "GSE201284"))

prad_meta <- cbind(prad_meta, dir_df) %>% as_tibble()
prad_meta <- mutate(prad_meta,
    pos = rowSums(across(starts_with("Dir_")) == "+"),
    concor = ifelse(pos < 5, 9 - pos, pos),
    com_dir = ifelse(sign(Effect) == -1, "-", "+"))



prad_meta_txt <- dplyr::select(prad_meta, MarkerName, Direction, Effect, StdErr,  `P-value`, p_adj)
library(NetActivityData)
data(gtex_gokegg_annot)
colnames(prad_meta_txt) <- c("GeneSet", "Direction", "Effect", "StdErr", "pvalue", "adj_pvalue")
prad_meta_txt <- right_join(dplyr::select(gtex_gokegg_annot, GeneSet, Term), prad_meta_txt, by = "GeneSet") %>%
  as_tibble() %>% arrange(pvalue)

write.table(prad_meta_txt, file = "results/PRAD_metanalysis/PRAD_metanalysis.tsv",
          sep = "\t", quote = FALSE, row.names = FALSE)


ggplot(prad_meta, aes(x = factor(concor), y = -log10(`P-value`))) +
 geom_boxplot() +
 theme_bw()


meta_sig <- subset(prad_meta, p_adj < 0.05)

dataset_meta <- data.frame(Dataset = c("GSE141551", "GSE169038", "GSE183019", "GSE201284", "GSE21034", "GSE46691", "GSE70768", "GSE70769", "TCGA-PRAD"),
                          Technology = c("HumanHT-12", "HuEx1.0", "RNAseq", "RNAseq", "HuEx1.0", "HuEx1.0", "HumanHT-12", "HumanHT-12", "RNAseq"),
                          N = c(503, 1131, 226, 123, 150, 545, 122, 91, 334))

meta <- select(meta_sig, MarkerName, starts_with("Dir_"), com_dir) %>%
  gather(Dataset, Direction, 2:10) %>%
  group_by(Dataset) %>%
  summarize(gene_set_cor = mean(com_dir == Direction)) %>%
  mutate(Dataset = gsub("Dir_", "", Dataset)) %>%
  left_join(dataset_meta, by = "Dataset")

ggplot(meta, aes(x = Dataset, y = gene_set_cor, fill = Technology)) +
  geom_bar(stat = "identity") +
  facet_grid(~ Technology, scale = "free_x") +
  theme_bw()


ggplot(meta, aes(x = N, y = gene_set_cor, color = Technology)) +
  geom_point() +
  theme_bw() +
  scale_x_log10()


datasets <- dir("results/PRAD_metanalysis/", pattern = "_metal.txt", full.names = TRUE)
names(datasets) <- gsub("results/PRAD_metanalysis//", "", gsub("_metal.txt", "", datasets))
metal_dfs <- lapply(datasets, read_table)

names(datasets)[9] <- "TCGA-PRAD"

#
library(metafor)

plotForest <- function(geneset, xlim = c(-1, 1.5)){
  res <- rma(sapply(metal_dfs, function(x) subset(x, Marker == geneset)$logFC),
            sei = sapply(metal_dfs, function(x) subset(x, Marker == geneset)$SE),
            weights = sapply(metal_dfs, function(x) subset(x, Marker == geneset)$N),
          slab = names(datasets))

  forest(res, xlim=xlim, rows = c(9, 3, 15, 16, 4, 5, 10, 11, 17), ylim=c(-1, 22),
    main = paste0(geneset, "\n", gtex_gokegg_annot[geneset, "Term"]), header = TRUE)
  par(font=4)
  text(min(xlim), c(18,12,6), pos = 4, c("RNAseq",
                                 "HumanHT-12",
                                 "HuEx1.0"))
}
plotWeights <- function(geneset){
  weights <- gtex_gokegg_annot[geneset, ]$Weights_SYMBOL[[1]]
  df <- data.frame(weight = weights, gene = names(weights)) %>%
      mutate(Direction = ifelse(weight > 0, "Positive", "Negative")) %>%
      arrange(desc(abs(weight)))
  df$gene <- factor(df$gene, levels = df$gene)
  df   %>%
      ggplot(aes(x = gene, y = abs(weight), fill = Direction)) +
      geom_bar(stat = "identity") +
      theme_bw() +
      ylab("Weight") +
      xlab("Gene") +
      theme(axis.text.x  = element_text(angle=45, vjust=0.5),
            text = element_text(size = 16))

}

p1 <- function()  plotForest("GO:0048012", xlim = c(-1.5, 1))
p2 <- function()  plotForest("GO:0051984", xlim = c(-1.5, 1))


plotForest("GO:0051984", xlim = c(-1.5, 1))
plotWeights("GO:0051984")


png("figures/PRAD_metaanalysis_GO1.png")
forest(res, xlim=c(-1, 1.5), rows = c(9, 3, 15, 16, 4, 5, 10, 11, 17), ylim=c(-1, 22),
  main = "GO:0051984\nPositive regulation of chromosome segregration", header = TRUE)
par(font=4)
text(-1, c(18,12,6), pos = 4, c("RNAseq",
                               "HumanHT-12",
                               "HuEx1.0"))
dev.off()

res <- rma(sapply(metal_dfs, function(x) subset(x, Marker == "GO:1904424")$logFC),
          sei = sapply(metal_dfs, function(x) subset(x, Marker == "GO:1904424")$SE),
          weights = sapply(metal_dfs, function(x) subset(x, Marker == "GO:1904424")$N),
        slab = names(datasets))

png("figures/PRAD_metaanalysis_GO2.png")
forest(res, xlim=c(-1, 1.5), rows = c(9, 3, 15, 16, 4, 5, 10, 11, 17), ylim=c(-1, 22),
  main = "GO:1904424\nregulation of GTP binding", header = TRUE)
par(font=4)
text(-1, c(18,12,6), pos = 4, c("RNAseq",
                               "HumanHT-12",
                               "HuEx1.0"))
dev.off()

## Check genes in original data
##### Compute normalized gene expression per gene as in NetActivity

sel_genes <- c(ESM1 = "ENSG00000164283", PAK1 = "ENSG00000149269", MET = "ENSG00000105976", SMC6 = "ENSG00000163029", CDT1 = "ENSG00000167513", RAD18 = "ENSG00000070950" )

### GSE46691
load("results/preprocess/GSE46691/GSE46691_ENSEMBL.Rdata")
gse46691_prep <- prepareSummarizedExperiment(gse46691_ensemb, "gtex_gokegg")

gse46691_prep$gleason <- factor(ifelse(gse46691_prep$`gleason score:ch1` >= 8, "High", "Low"), levels = c("Low", "High"))
gse46691_prep$metastasis <- gse46691_prep$`metastatic event:ch1`

mod_gse46691 <- model.matrix(~ gleason + metastasis, colData(gse46691_prep))
lm.gse46691 <- lmFit(assay(gse46691_prep), mod_gse46691) %>% eBayes()
tab.gse46691_gene <- topTable(lm.gse46691, coef = 2, n = Inf)

### GSE141551
load("results/preprocess/GSE141551/GSE141551_ENSEMBL.Rdata")
gse141551_prep <- prepareSummarizedExperiment(gse141551_se, "gtex_gokegg")
gse141551_prep$gleason <- factor(ifelse(gse141551_prep$`gleason_group:ch1` == "8-10", "High", "Low"), levels = c("Low", "High"))
gse141551_prep$race <- gse141551_prep$`race:ch1`
gse141551_prep$age <- gse141551_prep$`age_group:ch1`
gse141551_prep$stage <- gse141551_prep$`stage_disease:ch1`

mod_gse141551 <- model.matrix(~ gleason + race + age + stage, colData(gse141551_prep))
lm.gse141551 <- lmFit(assay(gse141551_prep), mod_gse141551) %>% eBayes()
tab.gse141551_gene <- topTable(lm.gse141551, coef = 2, n = Inf)

## GSE21034
load("results/preprocess/GSE21034/GSE21034_ENSEMBL.Rdata")
gse21034_prep <- prepareSummarizedExperiment(gse21034_ensemb, "gtex_gokegg")
gse21034_prep$gleason <- factor(ifelse(gse21034_prep$`biopsy_gleason_grade:ch1` >= 8, "High", "Low"), levels = c("Low", "High"))
gse21034_prep$metastasis <- gse21034_prep$`tumor type:ch1`

mod_gse21034 <- model.matrix(~ gleason + metastasis, colData(gse21034_prep))
lm.gse21034 <- lmFit(assay(gse21034_prep), mod_gse21034) %>% eBayes()
tab.gse21034_gene <- topTable(lm.gse21034, coef = 2, n = Inf)

## GSE70768
load("results/preprocess/GSE70768/GSE70768_ENSEMBL.Rdata")
gse70768_prep <- prepareSummarizedExperiment(gse70768_se_filt[rowMeans(is.na(assay(gse70768_se_filt))) == 0, ], "gtex_gokegg")
gse70768_prep$gleason <- factor(ifelse(sapply(strsplit(gse70768_prep$`tumour gleason:ch1`, "="), `[`, 1) %>% as.numeric() >= 8, "High", "Low"), levels = c("Low", "High"))
gse70768_prep$age <- as.numeric(gse70768_prep$`age at diag:ch1`)
gse70768_prep$tum_prop <- as.numeric(gsub("%", "", gse70768_prep$`tumour %:ch1`))

mod_gse70768 <- model.matrix(~ gleason + age + tum_prop, colData(gse70768_prep))
lm.gse70768 <- lmFit(assay(gse70768_prep)[, rownames(mod_gse70768)], mod_gse70768) %>% eBayes()
tab.gse70768_gene <- topTable(lm.gse70768, coef = 2, n = Inf)

## GSE70769
load("results/preprocess/GSE70769/GSE70769_ENSEMBL.Rdata")
gse70769_prep <- prepareSummarizedExperiment(gse70769_se_filt, "gtex_gokegg")
gse70769_prep$gleason <- factor(ifelse(sapply(strsplit(gse70769_prep$`tumour gleason:ch1`, "="), `[`, 1) %>% as.numeric() >= 8, "High", "Low"), levels = c("Low", "High"))

mod_gse70769 <- model.matrix(~ gleason, colData(gse70769_prep))
lm.gse70769 <- lmFit(assay(gse70769_prep), mod_gse70769) %>% eBayes()
tab.gse70769_gene <- topTable(lm.gse70769, coef = 2, n = Inf)

## GSE183019
load("results/preprocess/GSE183019/GSE183019_counts.Rdata")
gse183019_dds <- DESeqDataSet(gse183019_se, design = ~ 1)
gse183019_vst <- varianceStabilizingTransformation(gse183019_dds)

rownames(gse183019_vst) <- mapIds(org.Hs.eg.db,
    keys = rownames(gse183019_vst),
    column = 'ENSEMBL',
    keytype = 'SYMBOL')

gse183019_vst <- gse183019_vst[!is.na(rownames(gse183019_vst)), ]
gse183019_vst <- gse183019_vst[!duplicated(rownames(gse183019_vst)), ]

gse183019_prep <- prepareSummarizedExperiment(gse183019_vst, "gtex_gokegg")
gse183019_prep$gleason <- factor(ifelse(gse183019_prep$`gleason score:ch1` %in% c("4+4", "4+5", "5+4", "5+5 T4"), "High", "Low"), levels = c("Low", "High"))
gse183019_prep$age <- as.numeric(gse183019_prep$`age:ch1`)


mod_gse183019 <- model.matrix(~ gleason + age, colData(gse183019_prep))
lm.gse183019 <- lmFit(assay(gse183019_prep), mod_gse183019) %>% eBayes()
tab.gse183019_gene <- topTable(lm.gse183019, coef = 2, n = Inf)

## GSE201284
load("results/preprocess/GSE201284/GSE201284_counts.Rdata")

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

gse201284_prep$gleason <- factor(ifelse(gse201284_prep$`patient gleason score:ch1` %in% c("3+5", "3+5 T4", "4+4", "4+5", "5+4", "5+5", "4+4 T5", "4+5 T3", "5+4 T3"), "High", "Low"), levels = c("Low", "High"))
gse201284_prep$age <- as.numeric(gse201284_prep$`age:ch1`)
gse201284_prep$tumor_prop <- as.numeric(gse201284_prep$`% tumor:ch1`)


mod_gse201284 <- model.matrix(~ gleason + age + tumor_prop, colData(gse201284_prep))
lm.gse201284 <- lmFit(assay(gse201284_prep)[, rownames(mod_gse201284)], mod_gse201284) %>% eBayes()
tab.gse201284_gene <- topTable(lm.gse201284, coef = 2, n = Inf)

## TCGA
prad <- loadHDF5SummarizedExperiment("results/TCGA_gexp_coding_PRAD/", prefix = "vsd_norm_prad_tumor")
prad$gleason <- factor(ifelse(prad$paper_Reviewed_Gleason_category == ">=8", "High", "Low"), levels = c("Low", "High"))
assay(prad) <- data.matrix(assay(prad))

prad_prep <- prepareSummarizedExperiment(prad, "gtex_gokegg")

mod_prad <- model.matrix(~ gleason + paper_Subtype + age_at_index + race, colData(prad_prep))
lm.prad <- lmFit(assay(prad_prep), mod_prad) %>% eBayes()
tab.prad_gene <- topTable(lm.prad, coef = 2, n = Inf)


## GSE169038
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

mod_gse169038 <- model.matrix(~ gleason + race + decipher, colData(gse169038_prep))
lm.gse169038 <- lmFit(assay(gse169038_prep), mod_gse169038) %>% eBayes()
tab.gse169038_gene <- topTable(lm.gse169038, coef = 2, n = Inf)

sel_gos <- c("GO:0048012", "GO:0051984")
sel_tab_cols <- function(df1, df2, dataset_name){
  rbind(df1[sel_genes, ] %>%
    mutate(gene = names(sel_genes),
           dataset =  dataset_name),
       df2[sel_gos, ] %>%
    mutate(gene = sel_gos,
           dataset =  dataset_name)    
           ) %>%
    dplyr::select(logFC, P.Value, gene, dataset)
}

merge_tabs <- function(df1, df2){

  df1_sel <- dplyr::select(df1, -dataset)
  if (unique(df1$dataset) != ""){
    colnames(df1_sel)[1:2] <- paste0(colnames(df1_sel)[1:2], "_", unique(df1$dataset))
  }
  df2_sel <- dplyr::select(df2, -dataset)
  colnames(df2_sel)[1:2] <- paste0(colnames(df2_sel)[1:2], "_", unique(df2$dataset))

  left_join(df1_sel, df2_sel, by = "gene") %>%
    mutate(dataset = "")
}

gene_tabs <- Map(sel_tab_cols,  list(tab.gse141551_gene, tab.gse169038_gene, tab.gse183019_gene, tab.gse201284_gene, tab.gse21034_gene, tab.gse46691_gene, tab.gse70768_gene, tab.gse70769_gene, tab.prad_gene),
                  list(tab.gse141551, tab.gse169038, tab.gse183019, tab.gse201284, tab.gse21034, tab.gse46691, tab.gse70768, tab.gse70769, tab.prad), 
                  c("GSE141551", "GSE169038", "GSE183019", "GSE201284", "GSE21034", "GSE46691", "GSE70768", "GSE70769", "TCGA-PRAD"))
gene_tab <- Reduce(merge_tabs, gene_tabs) %>%
  gather(Name, Value, c(1:2, 4:19)) %>%
  mutate(Dataset = sapply(strsplit(Name, "_"), `[`, 2), 
         Measure = sapply(strsplit(Name, "_"), `[`, 1)) %>%
  dplyr::select(-Name) %>%
  spread(Measure, Value) %>%
  mutate(Signif = ifelse(P.Value < 0.05, "Significant", "Non-Significant")) %>%
  as_tibble()


go1_plot <- gene_tab %>%
  filter(gene %in% c("GO:0048012", "ESM1", "PAK1", "MET")) %>%
  mutate(gene = factor(gene, levels = c("GO:0048012", "ESM1", "PAK1", "MET")),
        Dataset = factor(Dataset, levels = c("TCGA-PRAD", "GSE201284", "GSE183019", "GSE70769", "GSE70768", "GSE141551", "GSE46691", "GSE21034", "GSE169038" ))) %>%
  ggplot(aes(x = gene, y = logFC, fill = Dataset, color = Signif )) +
  ylab("Standardized logFC") +
  xlab("") +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("green1", "green3", "green4", "lightblue", "blue", "blue4", "pink", "red", "red4")) +
  scale_color_manual(name = "", values = c("white", "black")) +
  theme_bw()



go2_plot <- gene_tab %>%
  filter(gene %in% c("GO:0051984", "SMC6", "CDT1", "RAD18")) %>%
  mutate(gene = factor(gene, levels = c("GO:0051984", "SMC6", "CDT1", "RAD18")),
        Dataset = factor(Dataset, levels = c("TCGA-PRAD", "GSE201284", "GSE183019", "GSE70769", "GSE70768", "GSE141551", "GSE46691", "GSE21034", "GSE169038" ))) %>%
  ggplot(aes(x = gene, y = logFC, fill = Dataset, color = Signif )) +
  geom_bar(stat = "identity", position = position_dodge()) +
  ylab("Standardized logFC") +
  xlab("") +
  scale_fill_manual(values = c("green1", "green3", "green4", "lightblue", "blue", "blue4", "pink", "red", "red4")) +
  scale_color_manual(name = "", values = c("white", "black")) +
  theme_bw()

## Figure 3
png("figures/PRAD_metaanalysis.png", width = 1500, height = 800)
plot_grid(
  plot_grid(ggdraw(p1), plotWeights("GO:0048012"), go1_plot, ncol = 3, labels = LETTERS[1:3]),
  plot_grid(ggdraw(p2), plotWeights("GO:0051984"), go2_plot, ncol = 3, labels = LETTERS[4:6]),
  nrow = 2
)
dev.off()
