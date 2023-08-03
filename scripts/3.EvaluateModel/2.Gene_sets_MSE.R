#'#################################################################################
#'#################################################################################
#' Explore GSE57945 features from different models
#'#################################################################################
#'#################################################################################


## Load libraries ####
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)
library(HDF5Array)
library(BiocParallel)
library(rjson)
library(cowplot)
library(ggcorrplot)

makePCdf <- function(seobj, vars){
  pc <- prcomp(t(assay(seobj)))
  pcdf <- data.frame(pc$x[, 1:10])
  pcdf <- cbind(pcdf, colData(seobj)[, vars])
}
makePCplot <- function(pcdf, var){
  ggplot(pcdf, aes_string(x = "PC1", y = "PC2", col = var)) +
  geom_point() +
  theme_bw()
}
readFeatures <- function(path, seobj){
  tab <- read.table(path, header = TRUE)
  se <- SummarizedExperiment(assay = t(data.matrix(tab)),
                              colData = colData(seobj))
}

read_training <- function(path, name){

  mat <- read.table(path, row.names = 1) %>% t()
  df <- mat[, c(1, 3)] %>%
      data.frame() %>%
      mutate(epoch = seq_len(nrow(mat))) %>%
      gather(measure, mse, 1:2) %>%
      mutate(model = name)
}

## Compare trainings
all.paths <- read_training("results/GTEx_coding/paths_all_full_v3.11/model_trained/GTEx_coding_training_evaluation.txt", "Gene Sets (Initial)")
base <- read_training("results/GTEx_coding/paths_filt2_full_v3.11/model_trained/GTEx_coding_training_evaluation.txt", "Gene Set")
base2 <- read_training("results/GTEx_coding/paths_filt2_full_v3.11/model_trained/GTEx_coding_training_evaluation.txt", "Whole training")

drop <- read_training("results/GTEx_coding/paths_filt2_full_drop_noprime_v3.7/model_trained/GTEx_coding_training_evaluation.txt", "Step 1 + step 3 + dropout")
drop.full <- read_training("results/GTEx_coding/paths_filt2_full_drop_prime_v3.9/model_trained/GTEx_coding_training_evaluation.txt", "Whole training + dropout")

pre <- read_training("results/GTEx_coding/paths_filt2_full_predense_v6.2/model_trained/GTEx_coding_training_evaluation.txt", "Dense + Gene Set")
post <- read_training("results/GTEx_coding/paths_filt2_full_postdense_v4.3/model_trained/GTEx_coding_training_evaluation.txt", "Gene Set + Dense")
pre_post <- read_training("results/GTEx_coding/paths_filt2_full_prepostdense_v5.3/model_trained/GTEx_coding_training_evaluation.txt", "Dense + Gene Set + Dense")


df.models <- Reduce(rbind, list(base, post, pre, pre_post)) %>%
  mutate(model = factor(model, levels = c("Gene Set", "Gene Set + Dense", "Dense + Gene Set", "Dense + Gene Set + Dense" )),
          dataset = ifelse(measure == "loss", "Training", "Validation"))


# Sup Figure 3
png("figures/TCGA.pathways.trainingeval_models.png", width = 2700, height = 1000, res = 300)
df.models %>%
filter(epoch > 1) %>%
  ggplot(aes(x = epoch, y = mse, color = dataset, group = dataset)) +
  geom_point() +
  geom_line() +
  facet_grid(~ model) +
  theme_bw() +
  scale_color_discrete(name = "") +
  scale_y_continuous(name = "MSE")
dev.off()

df.train <- Reduce(rbind, list(base2, drop, drop.full)) %>%
  mutate(model = factor(model, levels = c("Whole training", "Whole training + dropout", "Step 1 + step 3 + dropout")),
          dataset = ifelse(measure == "loss", "Training", "Validation"))

# Sup Figure 5
png("figures/TCGA.pathways.trainingeval_training.png", width = 2100, height = 900, res = 300)
df.train  %>%
filter(epoch > 1) %>%
ggplot(aes(x = epoch, y = mse, color = dataset, group = dataset)) +
  geom_point() +
  geom_line() +
  facet_grid( ~ model) +
  theme_bw() +
  scale_color_discrete(name = "") +
  scale_y_continuous(name = "MSE")
dev.off()
