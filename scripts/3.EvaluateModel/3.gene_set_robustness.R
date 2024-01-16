#'#################################################################################
#'#################################################################################
#' Explore replicability in gene set activity scores from different subsets and random gene sets
#'#################################################################################
#'#################################################################################


## Load libraries ####
library(SummarizedExperiment)
library(tidyverse)
library(HDF5Array)
library(NetActivity)

## Load data ####
gtex <- loadHDF5SummarizedExperiment("results/GTEx/", prefix = "vst_all_group_")
norm <- prepareSummarizedExperiment(gtex, "gtex_gokegg")



computeGSAS <- function(size, model, sufix){
  lapply(sufix, function(i){
    load(paste0("results/subsetGeneSets/models/size_", size, "/", model, "/", i, "/NetActivity_weights/Dataset_NetActivity_weights.Rdata"))
    scores <- computeGeneSetScores(norm, weights)
    tab <- data.matrix(assay(scores))
    tab
  })
}

pathwayCorr <- function(path_list, row){
  path_mat <- sapply(path_list, function(x) x[row, ])
  cors <- cor(path_mat, method = "spearman")
  cors[upper.tri(cors)]
}

readCors <- function(size, model){
  base <- computeGSAS(size, model, letters[1:6])
  paths.name <- rownames(base[[1]]) 
  base.cors <- sapply(paths.name, pathwayCorr, path_list = base)
  colnames(base.cors) <- paths.name
   df.base <- makeDFsum(base.cors, size, model)
}

makeDFsum <- function(cors, size, model){
  df <- data.frame(path = colnames(cors), minCor = colMins(abs(cors)), medCor = colMedians(abs(cors))) %>%
    mutate(class = ifelse(minCor > 0.8, "High", ifelse(minCor < 0.3, "low", "intermediate")),
            Size = size,
            Replicate = model) %>%
            as_tibble()

}

## Subset gene sets
sizes <- c(30, 100, 300, 1000)
full_df <- lapply(sizes, function(size){
    mod_l <- lapply(1:10, function(model){
        message("Size: ", size, " Model: ", model)
        corDF <- readCors(as.character(size), as.character(model))
    })
    mod_l
})
geneset_df <- Reduce(rbind, unlist(full_df, recursive = FALSE)) %>%
  mutate(Size = factor(Size, levels = c("30", "100", "300", "1000")))
save(geneset_df, file = "results/subsetGeneSets/replicability_df.Rdata")

png("figures/geneset_subsets_replicability.png", width = 1800, height = 900, res = 300)
ggplot(geneset_df, aes(x = Replicate, y = minCor)) +
  geom_boxplot() +
  ylab("Robustness") +
  xlab("Sample Size") +
  guides(color = FALSE) +
  theme_bw() +
  facet_grid(~ Size) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
dev.off()

