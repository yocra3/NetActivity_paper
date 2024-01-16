#'#################################################################################
#'#################################################################################
#' Explore replicability in pathways activations from different models
#'#################################################################################
#'#################################################################################


## Load libraries ####
# library(topGO)
library(SummarizedExperiment)
library(tidyverse)
library(rjson)
library(rhdf5)
library(HDF5Array)
library(cowplot)



readPathways <- function(model, sufix, path_name){
  lapply(sufix, function(i){
    path <- paste0("results/GTEx_coding/", model, i, "/model_features/prune_low_magnitude_dense.tsv")
    tab <- read_table(path, col_names = TRUE)
    tab <- data.matrix(tab)
    colnames(tab) <- path_name
    tab
  })
}


readPathways2 <- function(model, sufix, path_name){
  lapply(sufix, function(i){
    path <- paste0("results/GTEx_coding/", model, i, "/model_features/prune_low_magnitude_dense_1.tsv")
    tab <- read_table(path, col_names = TRUE)
    tab <- data.matrix(tab)
    colnames(tab) <- path_name
    tab
  })
}


pathwayCorr <- function(path_list, col){
  path_mat <- sapply(path_list, function(x) x[, col])
  cors <- cor(path_mat, method = "spearman")
  cors[upper.tri(cors)]
}

makeDFsum <- function(cors, mod_name){
  df <- data.frame(path = colnames(cors), minCor = colMins(abs(cors)), medCor = colMedians(abs(cors))) %>%
    mutate(class = ifelse(minCor > 0.8, "High", ifelse(minCor < 0.3, "low", "intermediate")),
            model = mod_name) %>%
            as_tibble()

}

kegg.map <- read.table("results/preprocess/go_kegg_gene_map.tsv", header = TRUE)
kegg.N <- table(kegg.map$PathwayID)

kegg.annot <- fromJSON(file = "data/kegg_pathways.json")
kegg.df <- lapply(kegg.annot$children, function(x) {
  top_cat <- x$name
  paths.df <- lapply(x$children, function(y){
    cat2 <- y$name
    paths <- sapply(y$children, function(z) z$name)
    data.frame(top_cat = top_cat, category = cat2, path = paths)
  })
  Reduce(rbind, paths.df) %>%
    as_tibble()
}) %>%
 Reduce(rbind, .)
kegg.df <- kegg.df %>%
  mutate(pathID = paste0("path:hsa", substring(path, 0, 5)),
          pathName = gsub("^[0-9]*  ", "", path))

paths <- read.table("results/GTEx_coding/paths_filt2_full_v3.11/model_trained/pathways_names.txt", header = TRUE)
paths.ini <- read.table("results/GTEx_coding/paths_all_full_v3.11/model_trained/pathways_names.txt", header = TRUE)

paths.vec <- as.character(paths[, 1])
paths.ini <- as.character(paths.ini[, 1])

kegg.df.com <- subset(kegg.df, pathID %in% paths.vec)
kegg.genes.N <- kegg.map %>%
  group_by(Symbol) %>%
  summarize(N = n())

## Load correlations
readCors <- function(base, models, paths.name, model.name){
  base <- readPathways(base, models, paths.name)
  base.cors <- sapply(paths.name, pathwayCorr, path_list = base)
  colnames(base.cors) <- paths.name
  # rownames(base.cors) <- c("1-2", "1-3", "2-3", "1-4", "2-4", "3-4", "1-5", "2-5", "3-5", "4-5")
  df.base <- makeDFsum(base.cors, model.name) %>%
    left_join(data.frame(kegg.N) %>% mutate(path = Var1) %>% select(-Var1), by = "path")

}

all_train <- readCors("paths_all_full_v3.11", c("", letters[1:5]), paths.ini, "All gene sets") %>%
  mutate(training = "Step 1 + step 2 + step 3")

all_nofrozen <- readCors("paths_all_pretrain_v3.10", c("", letters[1:5]), paths.ini, "All gene sets") %>%
  mutate(training = "Step 1 + step 3")

all_init <- readCors("paths_all_pretrain_v3.8", c("", letters[1:5]), paths.ini, "All gene sets") %>%
  mutate(training = "Step 1")


main <- readCors("paths_filt2_full_v3.11", c("", letters[1:5]), paths.vec, "Selected gene sets") %>%
  mutate(training = "Step 1 + step 2 + step 3")
drop <- readCors("paths_filt2_full_drop_noprime_v3.7", c("", letters[1:5]), paths.vec, "Selected gene sets") %>%
  mutate(training = "Dropout")
pre <- readCors("paths_filt2_pre_v3.8", c("", letters[1:5]), paths.vec, "Selected gene sets") %>%
  mutate(training = "Step 1")
drop_full <- readCors("paths_filt2_full_drop_prime_v3.9", c("", letters[1:5]), paths.vec, "Selected gene sets") %>%
  mutate(training = "Step 2 + dropout")

unfrozen <- readCors("paths_filt2_unfrozen_v3.10", c("", letters[1:5]), paths.vec, "Selected gene sets") %>%
  mutate(training = "Step 1 + step 3")

no_prime <- readCors("paths_filt2_full_noprime_v3.12", c("", letters[1:5]), paths.vec, "Selected gene sets") %>%
  mutate(training = "Step 3")

post2 <- readCors("paths_filt2_full_postdense_v4.3", c("", letters[1:5]), paths.vec, "Selected gene sets")  %>%
  mutate(training = "Gene Set + Dense")


readCors2 <- function(base, models, paths.name, model.name){
  base <- readPathways2(base, models, paths.name)
  base.cors <- sapply(paths.name, pathwayCorr, path_list = base)
  colnames(base.cors) <- paths.name
  # rownames(base.cors) <- c("1-2", "1-3", "2-3", "1-4", "2-4", "3-4", "1-5", "2-5", "3-5", "4-5")
  df.base <- makeDFsum(base.cors, model.name) %>%
    left_join(data.frame(kegg.N) %>% mutate(path = Var1) %>% select(-Var1), by = "path")

}
pre2 <- readCors2("paths_filt2_full_predense_v6.2",c("", letters[1:5]), paths.vec, "Selected gene sets") %>%
  mutate(training = "Dense + Gene Set")
pre2.post <- readCors2("paths_filt2_full_prepostdense_v5.3", c("", letters[1:5]), paths.vec, "Selected gene sets") %>%
  mutate(training = "Dense + Gene Set + Dense")

df.path_sel <- Reduce(rbind, list(main, pre, unfrozen, all_train, all_nofrozen, all_init)) %>%
  left_join(mutate(kegg.df.com, path = pathID) %>% select(path, top_cat, category), by = "path") %>%
  mutate(group = ifelse(model == "Selected gene sets", "Selected GOs and KEGGs", "All GOs + KEGGs"),
         training = recode(training, "pretrained" = "Step 1", "pretrained only" = "Step 1",
                            "whole training" = "Step 1 + step 2 + step 3", "primed + pretrained" = "Step 1 + step 2 + step 3",
                            "unfrozen" = "Step 1 + step 3", "Whole training, without adaptation" = "Step 1 + step 3"),
         training = factor(training, levels = c("Step 1", "Step 1 + step 3",  "Step 1 + step 2 + step 3")))


### Sup Figure 1
png("figures/robustness_goSelection.png", height = 900, width = 2500, res = 300)
df.path_sel %>%
  filter(training == "Step 1 + step 2 + step 3") %>%
  ggplot(aes(x = group, y = minCor)) +
 geom_boxplot() +
 scale_x_discrete(name = "") +
 scale_y_continuous(name = "Robustness") +
 theme_bw() +
 theme(plot.title = element_text(hjust = 0.5),
                  text = element_text(size = 20))
dev.off()

## Figure 1D
plot_rep2 <- Reduce(rbind, list(main, unfrozen, no_prime)) %>%
  left_join(mutate(kegg.df.com, path = pathID) %>% select(path, top_cat, category), by = "path") %>%
  mutate(group = ifelse(model == "Selected gene sets", "Selected GOs and KEGGs", "All GOs + KEGGs"),
         training = recode(training, "Step 1 + step 3" = "Step 1 + Step 3", 
                                    "Step 1 + step 2 + step 3" = "Step 1 + Step 2 + Step 3"),
         training = factor(training, levels = c("Step 3", "Step 1 + Step 3",  "Step 1 + Step 2 + Step 3"))) %>%
  ggplot(aes(x = training, y = minCor)) +
  geom_boxplot() +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Robustness") +
  theme_bw() +
 theme(text = element_text(size = 20))

png("figures/mainModel_robustness.png", height = 900, width = 2500, res = 300)
plot_rep2
dev.off()

## Sup Figure 1
plot_genes <-  df.path_sel %>%
  filter(group == "All GOs + KEGGs") %>%
  ggplot(aes(x = Freq, y = minCor)) +
    geom_point() +
    scale_x_log10(name = "Genes per gene set") +
    scale_y_continuous(name = "Robustness") +
    theme_bw() +
    facet_wrap(~ training) +
    geom_vline(xintercept = 30, linetype = "dashed", color = "grey")

png("figures/minCor_pretraning_Ngenes.png", height = 600, width = 1500, res = 300)
plot_genes
dev.off()

table(ifelse(grepl("GO", main$path), "GO", "KEGG"))
df.path_sel %>% group_by(training, group) %>%
summarize(p = mean(minCor > 0.7))

## Models comparison
df.mod <- Reduce(rbind, list(main,post2, pre2, pre2.post)) %>%
  left_join(mutate(kegg.df.com, path = pathID) %>% select(path, top_cat, category), by = "path") %>%
  mutate( training = recode(training, `Step 1 + step 2 + step 3` = "Gene Set"),
          training = factor(training , levels = c("Gene Set", "Gene Set + Dense", "Dense + Gene Set", "Dense + Gene Set + Dense")))

## Sup Figure 4
png("figures/minCor_models.png", width = 2400, height = 900, res = 300)
df.mod %>%
  ggplot(aes(x = training, y = minCor)) +
  geom_boxplot() +
  theme_bw() +
  scale_y_continuous(name = "Robustness") +
  xlab("Network structure")
dev.off()

## Training comparison
df.train <- Reduce(rbind, list(main, drop, drop_full)) %>%
  left_join(mutate(kegg.df.com, path = pathID) %>% select(path, top_cat, category), by = "path") %>%
  mutate(Training = recode(training, "Step 1 + step 2 + step 3" = "Whole training",
        "Step 2 + dropout" = "Whole training + dropout",
        Dropout = "Step 1 + step 3 + dropout"),
        Training = factor(Training , levels = c("Whole training", "Whole training + dropout", "Step 1 + step 3 + dropout")),)

## Save correlations in data.frames to speed up computations
save(df.path_sel, df.mod, df.train, file = "results/manuscript/df_gtex_training_replicability.Rdata")


## Sup Figure 6
png("figures/minCor_training.png", width = 1800, height = 900, res = 300)
df.train %>%
  ggplot(aes(x = Training, y = minCor)) +
  geom_boxplot() +
  theme_bw() +
  scale_y_continuous(name = "Robustness")
dev.off()

## Plot correlation of worse path
path_vals <- readPathways("paths_filt3_full_v3.6", sufix = c("", letters[1:5]), path_name = paths.vec)
path_mat <- sapply(path_vals, function(m) m[, 427 ])
path_mat2 <- path_mat
path_mat2[, 5] <-  - path_mat2[, 5]

## Evaluate training robustness
## Random gene sets
library(NetActivity)
## Load data ####
gtex <- loadHDF5SummarizedExperiment("results/GTEx/", prefix = "vst_all_group_")
norm <- prepareSummarizedExperiment(gtex, "gtex_gokegg")

pathwayCorr2 <- function(path_list, row){
  path_mat <- sapply(path_list, function(x) x[row, ])
  cors <- cor(path_mat, method = "spearman")
  cors[upper.tri(cors)]
}

randomGSAS <- lapply(letters[1:6], function(i){
    load(paste0("results/randomGeneSets/", i, "/NetActivity_weights/Dataset_NetActivity_weights.Rdata"))
    scores <- computeGeneSetScores(norm, weights)
    tab <- data.matrix(assay(scores))
    tab
  })
paths.name2 <- rownames(randomGSAS[[1]]) 
rand.cors <- sapply(paths.name2, pathwayCorr2, path_list = randomGSAS)
colnames(rand.cors) <- paths.name2
df.random <- makeDFsum(rand.cors, "Random")

df.gene_sets <- rbind(main %>% select(-Freq, -training), df.random) %>%
  rbind(all_train %>% select(-Freq, -training))


### Sup Figure 3
png("figures/robustness_randomGeneSets.png", height = 900, width = 2500, res = 300)
df.gene_sets %>%
  mutate(model = factor(model, levels = c("All gene sets", "Selected gene sets", "Random"))) %>%
  ggplot(aes(x = model, y = minCor)) +
   geom_boxplot() +
 scale_x_discrete(name = "") +
 scale_y_continuous(name = "Robustness") +
 theme_bw() +
 theme(plot.title = element_text(hjust = 0.5),
                  text = element_text(size = 20))
dev.off()
