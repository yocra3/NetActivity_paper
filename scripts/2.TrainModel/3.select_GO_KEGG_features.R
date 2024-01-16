#'#################################################################################
#'#################################################################################
#' Select GO and KEGG features for model
#'#################################################################################
#'#################################################################################


## Load libraries ####
library(tidyverse)
library(parallel)
library(matrixStats)
library(ggrepel)
library(ggcorrplot)
library(cowplot)

## Load functions
readPathways <- function(model, sufix, path_name){
  lapply(sufix, function(i){
    path <- paste0("results/GTEx_coding/", model, i, "/model_features/prune_low_magnitude_dense.tsv")
    tab <- read.table(path, header = TRUE)
    tab <- data.matrix(tab)
    colnames(tab) <- path_name
    tab
  })
}

pathwayCorr <- function(path_list, col){
  path_mat <- sapply(path_list, function(x) x[, col])
  cors <- cor(path_mat)
  cors[upper.tri(cors)]
}


makeDFsum <- function(cors, mod_name){
  df <- data.frame(path = colnames(cors), minCor = colMins(abs(cors)), medCor = colMedians(abs(cors))) %>%
    mutate(class = ifelse(minCor > 0.8, "High", ifelse(minCor < 0.3, "low", "intermediate")),
            model = mod_name) %>%
            as_tibble()

}

readCors <- function(base, models, paths.name, model.name){
  base <- readPathways(base, models, paths.name)
  base.cors <- sapply(paths.name, pathwayCorr, path_list = base)
  colnames(base.cors) <- paths.name
  # rownames(base.cors) <- c("1-2", "1-3", "2-3", "1-4", "2-4", "3-4", "1-5", "2-5", "3-5", "4-5")
  df.base <- makeDFsum(base.cors, model.name) %>%
    left_join(data.frame(kegg.N) %>% mutate(path = Var1) %>% select(-Var1), by = "path")

}
## Load pathways annotations
kegg.map <- read.table("results/preprocess/go_kegg_gene_map.tsv", header = TRUE)
paths <- read.table("results/GTEx_coding/paths_all_full_v3.11/model_trained/pathways_names.txt", header = TRUE)
paths.vec <- as.character(paths[, 1])
input_genes <- read.table("results/GTEx_coding/input_genes.txt", header = FALSE)
paths.ini <- read.table("results/GTEx_coding/paths_all_full_v3.11/model_trained/pathways_names.txt", header = TRUE)
paths.ini <- as.character(paths.ini[, 1])



kegg.map.com <- subset(kegg.map, PathwayID %in% paths.vec & Symbol %in% input_genes$V1)
kegg.N <- table(kegg.map.com$PathwayID) %>%
  data.frame()

### Define similarity between pathways
path_genes <- mclapply(paths.vec, function(x) subset(kegg.map.com, PathwayID == x & !is.na(Symbol))$Symbol, mc.cores = 10)
names(path_genes) <- paths.vec

gene_mat <- matrix(0, length(paths.vec), length(unique( kegg.map.com$Symbol)), dimnames = list(paths.vec, unique(kegg.map.com$Symbol)))
 for (i in paths.vec) gene_mat[i, path_genes[[i]]] <- 1
gene_d <- dist(gene_mat, "binary")
save(gene_d, file = "results/GTEx_coding/go_kegg_pathways_distance.Rdata")
gene_dmat <- as.matrix(gene_d)


## Compare similarity and replicability
all_train <- readCors("paths_all_full_v3.11", c("", letters[1:5]), paths.ini, "All gene sets") %>%
  mutate(training = "Step 1 + step 2 + step 3")

gene_sim <- gene_dmat 
diag(gene_sim) <- 1
maxSim <- rowMins(gene_sim)

all_train$Sim <- 1 - maxSim[all_train$path]

png("figures/geneset_similary_robustness.png", width = 1800, height = 900, res = 300)
all_train %>% 
  mutate(Selection = ifelse(Freq > 30, "Large", "Small"),
  Select = factor(Selection, levels = c("Small", "Large"))) %>%
ggplot(aes(x = Sim, y = minCor)) +
  geom_point(alpha = 0.2) +
  geom_smooth() +  
  ylab("Robustness") +
  xlab("Similarity") +
  theme_bw() +
  facet_grid(~ Selection)
dev.off()

## Select gene sets by number of genes and similarity
good_paths <- as.character(subset(kegg.N, Freq < 30 & Freq >= 10)$Var)
gene_dmat_sel <- gene_dmat[good_paths, good_paths]
gene_dmini <- as.dist(gene_dmat_sel)

sel_paths <- good_paths
gene_dloop <- gene_dmini
gene_dmat_loop <- gene_dmat_sel
path_cls <- cutree(hclust(gene_dmini), h = 0.5)
while(length(sel_paths) != length(unique(path_cls))){
  sel_paths <- sapply(unique(path_cls), function(cl){
    paths <- rownames( gene_dmat_loop)[path_cls == cl]
    df.sub <- subset(kegg.N, Var1 %in% paths)
    as.character(df.sub$Var1[which.max(df.sub$Freq)])
  })
  gene_dmat_loop <- gene_dmat_loop[sel_paths, sel_paths]
  gene_dloop <- as.dist(gene_dmat_loop)
  path_cls <- cutree(hclust(gene_dloop), h = 0.5)
}
gene_dmat_filt <- gene_dmat[sel_paths, sel_paths]
diag(gene_dmat_filt) <- 1
summary(rowMins(gene_dmat_filt))

kegg.map.filt <- subset(kegg.map, PathwayID %in% sel_paths)
write.table(kegg.map.filt, file = "results/GTEx_coding/go_kegg_filt_gene_map.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")


## Select gene sets with a replicability > 0.7
paths2 <- read.table("results/GTEx_coding/paths_filt1_full_v3.11/model_trained/pathways_names.txt", header = TRUE)
paths.vec2 <- as.character(paths2[, 1])

path_vals_full2b <- readPathways("paths_filt1_full_v3.11", sufix = c("", letters[1:5]), path_name = paths.vec2)

full.cors2b <- sapply(paths.vec2, pathwayCorr, path_list = path_vals_full2b)
colnames(full.cors2b) <- paths.vec2
df.full2b <- makeDFsum(full.cors2b,"full") %>%
  left_join(data.frame(kegg.N) %>% mutate(path = Var1) %>% dplyr::select(-Var1), by = "path") %>%
  mutate(Database = ifelse(substring(path, 1, 2) == "GO", "GO", "KEGG"))


sel_paths2 <- subset(df.full2b, minCor > 0.7)$path
kegg.map.filt2 <- subset(kegg.map, PathwayID %in% sel_paths2 & Symbol %in% input_genes$V1)
write.table(kegg.map.filt2, file = "results/GTEx_coding/go_kegg_filt2_gene_map.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")


## Evaluate final gene sets 
paths3b <- read.table("results/GTEx_coding/paths_filt2_full_v3.11/model_trained/pathways_names.txt", header = TRUE)
paths.vec3b <- as.character(paths3b[, 1])

path_vals_full3b <- readPathways("paths_filt2_full_v3.11", sufix = c("", letters[1:5]), path_name = paths.vec3b)

full.cors3b <- sapply(paths.vec3b, pathwayCorr, path_list = path_vals_full3b)
colnames(full.cors3b) <- paths.vec3b
df.full3b <- makeDFsum(full.cors3b,"full") %>%
  left_join(data.frame(kegg.N) %>% mutate(path = Var1) %>% dplyr::select(-Var1), by = "path") %>%
  mutate(Database = ifelse(substring(path, 1, 2) == "GO", "GO", "KEGG"))

full.cors3b_main <- sapply(paths.vec3b, function(x) median(abs(cor(sapply(path_vals_full3b, function(y) y[, x]))[1,-1])))

df.3b_worse <- df.full3b %>%
  mutate(main_cor = full.cors3b_main) %>%
  filter(minCor < 0.7)