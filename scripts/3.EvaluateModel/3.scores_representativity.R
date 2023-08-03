#'#################################################################################
#'#################################################################################
#' Representativity of GSAS
#'#################################################################################
#'#################################################################################


## Load libraries
library(SummarizedExperiment)
library(HDF5Array)
library(DESeq2)
library(NetActivity)
library(NetActivityData)
library(tidyverse)
library(cowplot)
library(parallel)
library(vegan)
library(ggcorrplot)

data(gtex_gokegg)
input_genes <- read.table("results/GTEx/input_genes.txt")

gtex.vst <- loadHDF5SummarizedExperiment("results/GTEx/", prefix = "vst_all_")
assay(gtex.vst) <- data.matrix(assay(gtex.vst))
preproc <- prepareSummarizedExperiment(gtex.vst, "gtex_gokegg")
gtex_scores <- computeGeneSetScores(preproc, "gtex_gokegg")
rda_gtex <- rda(t(assay(gtex.vst)) ~ t(assay(gtex_scores)))
save(rda_gtex, file = "results/manuscript/mainModel_GTEx_explainability.Rdata")


# TCGA
## Load  data
load("data/tcga_gexp_combat.Rdata")
## Restrict analysis to input genes
tcga_input <- gexp_tcga_combat[as.character(input_genes$V1), ]

## Overall variability
counts <- assay(tcga_input)
mode(counts) <- "integer"
counts[is.na(counts)] <- max(counts, na.rm = TRUE)
deseq <- DESeqDataSetFromMatrix(countData = counts,
              colData = colData(gexp_tcga_combat),
              design = ~ 1)
vst <-  vst(deseq, blind=FALSE)

preproc <- prepareSummarizedExperiment(vst, "gtex_gokegg")
all_scores <- computeGeneSetScores(preproc, "gtex_gokegg")
rda_tcga <- rda(t(assay(vst)) ~ t(assay(all_scores)))
save(rda_tcga, file = "results/manuscript/mainModel_TCGA_explainability.Rdata")

## Select tumors
tumor_tab <- table(gexp_tcga_combat$project_id)
sel_tumors <- names(tumor_tab)
names(sel_tumors) <- sel_tumors

## Compute PCs
tumor_pcs <- mclapply(sel_tumors, function(tum){
  message(tum)
  set <- tcga_input[, tcga_input$project_id == tum]

  counts <- assay(set)
  mode(counts) <- "integer"
  counts[is.na(counts)] <- max(counts, na.rm = TRUE)
  deseq <- DESeqDataSetFromMatrix(countData = counts,
                      colData = colData(set),
                      design = ~ 1)
  vst <-  vst(deseq, blind=FALSE)
  pc_vst <- prcomp(t(assay(vst)), rank. = 10)

  preproc <- prepareSummarizedExperiment(vst, "gtex_gokegg")
  scores <- computeGeneSetScores(preproc, "gtex_gokegg")
  pc_scores <- prcomp(t(assay(scores)), rank. = 10)

  list(
    raw = vst,
    vst = pc_vst,
    scores = pc_scores)
}, mc.cores = 5)
rda_pc10 <- mclapply(tumor_pcs, function(y) rda(t(assay(y$raw)) ~ y$scores$x), mc.cores = 10)
var_prop10 <- sapply(rda_pc10, function(i) i$CCA$tot.chi/i$tot.chi)
pc_prop10 <- sapply(tumor_pcs, function(x) sum(x$vst$sdev[1:10]**2)/sum(x$vst$sdev**2))


cor_prad <- cor(tumor_pcs[["TCGA-PRAD"]]$vst$x, tumor_pcs[["TCGA-PRAD"]]$scores$x)

## Figure 1I
corplot_prad <- ggcorrplot(abs(cor_prad), method = "circle", hc.order = FALSE,
      title = "TCGA-PRAD") +
  scale_x_discrete("Original gene expression") +
  scale_y_discrete(name = "NetActivity GSAS") +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
   axis.title.x = element_text(angle = 0, color = 'grey20', size = 14),
   axis.title.y = element_text(angle = 90, color = 'grey20', size = 14))


png("figures/corPlot_PRAD.png", width = 2000, height = 1000, res = 300)
corplot_prad
dev.off()


## IBD cohort
vst_ibd <- loadHDF5SummarizedExperiment("results/SRP042228/", prefix = "vsd_norm")
vst_ibd_input <- vst_ibd[as.character(input_genes$V1), ]
pc_vst_ibd <- prcomp(t(assay(vst_ibd_input)), rank. = 10)

assay(vst_ibd_input) <- data.matrix(assay(vst_ibd_input))
preproc_ibd <- prepareSummarizedExperiment(vst_ibd_input, "gtex_gokegg")
scores_ibd <- computeGeneSetScores(preproc_ibd, "gtex_gokegg")
pc_scores_ibd <- prcomp(t(assay(scores_ibd)), rank. = 10)
rda_ibd <- rda(t(assay(vst_ibd_input)) ~ pc_scores_ibd$x)

pc_summary <- data.frame(prop10 = c(var_prop10, rda_ibd$CCA$tot.chi/rda_ibd$tot.chi),
          dataset = c(names(var_prop10), "IBD"),
          pc10 = c(pc_prop10, sum(pc_vst_ibd$sdev[1:10]**2)/sum(pc_vst_ibd$sdev**2)))


summary(lm(prop10 ~ pc10, pc_summary))


# Figure 1H
corplot_ibd <- ggcorrplot(abs(cor(pc_vst_ibd$x, pc_scores_ibd$x)), method = "circle", hc.order = FALSE,
      title = "IBD") +
  scale_x_discrete("Original gene expression") +
  scale_y_discrete(name = "NetActivity GSAS") +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
   axis.title.x = element_text(angle = 0, color = 'grey20', size = 14),
   axis.title.y = element_text(angle = 90, color = 'grey20', size = 14))
png("figures/corPlot_IBD.png", width = 2000, height = 1000, res = 300)
corplot_ibd
dev.off()

load("results/manuscript/mainModel_GTEx_explainability.Rdata")
load("results/manuscript/mainModel_TCGA_explainability.Rdata")

# Figure 1G
plot_represent <- data.frame(Perc = c(rda_gtex$CCA$tot.chi/rda_gtex$tot.chi, rda_tcga$CCA$tot.chi/rda_tcga$tot.chi), 
                             Dataset = c("GTEx", "TCGA")) %>%
  ggplot(aes(x = Dataset, y = Perc*100)) +
  geom_bar(stat = "identity") +
  xlab("Dataset") +
  ylab("Representativity") +
  theme_bw() +
  ylim(c(0, 100)) +
  theme(text = element_text(size = 20))

png("figures/mainModel_representativity.png", width = 2000, height = 1000, res = 300)
plot_represent
dev.off()
