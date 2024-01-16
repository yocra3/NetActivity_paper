#'#################################################################################
#'#################################################################################
#' Prepare files to generate models with less or random gene sets 
#'#################################################################################
#'#################################################################################

## Load libraries
library(tidyverse)

original <- read_table("results/GTEx_coding/go_kegg_filt2_gene_map.tsv")
colnames(original) <- c("GeneSetID", "Gene")

## Subset gene sets
rep <- 10
sizes <- c(30, 100, 300, 1000)

gene_sets <- unique(original$GeneSetID)

for (size in sizes){
    for (i in seq_len(rep)){
        sel_gene_sets <- sample(gene_sets, size)
        sub <- subset(original, GeneSetID %in% sel_gene_sets)
        write_delim(sub, file = paste0("results/subsetGeneSets/subset_", i, "_N", size, "_geneMap.tsv"),
        delim = "\t")
    }
}


## Random gene set
sizes <- as.vector(table(original$GeneSetID))
size_freq <- table(sizes)

### Get gene set sizes
set.seed(27)
rand_sizes <- sample(as.numeric(names(size_freq)), size = 1518, replace = TRUE, prob = size_freq/sum(size_freq))

un_genes <- unique(original$Gene)
genes <- lapply(rand_sizes, function(size){
    genes <- sample(un_genes, size = size, replace = FALSE)
})
rand_gene_sets <- data.frame(GeneSetID = paste0("GeneSet", rep(seq_len(length(genes)), lengths(genes))),
                    Gene = unlist(genes))
write_delim(rand_gene_sets, file = "results/randomGeneSets/random_geneMap.tsv", delim = "\t")
