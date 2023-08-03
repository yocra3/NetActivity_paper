#'#################################################################################
#'#################################################################################
#' Define GO and KEGG terms to include in the model
#'#################################################################################
#'#################################################################################

library(GOfuncR)
library(tidyverse)
library(parallel)


## Load new GO graph
godir <-  "data/GO_terms/"
term <- read.table(paste0(godir, "/term.txt"), sep = "\t", quote = "", comment.char = "", as.is = TRUE)
graph <- read.table(paste0(godir, "/graph_path.txt"), sep = "\t", quote = "", comment.char = "", as.is = TRUE)

genes <- as.character(read.table("./results/TCGA_gexp_combat_coding/input_genes.txt")$V1)


## Get all GO terms
all_gos <- get_child_nodes("GO:0008150", term, graph)

## Remove obsolete terms
term_filt <- subset(term, V5 == 0)
all_gos <- subset(all_gos, child_go_id %in% term_filt$V4)
genes_pairs <- get_anno_genes(go_ids = all_gos$child_go_id, term_df = term, graph_path_df = graph)

tab <- left_join(mutate(all_gos, go_id = child_go_id), genes_pairs, by = "go_id") %>%
  as_tibble() %>%
  filter(!is.na(gene))

tab$Symbol <- mapIds(org.Hs.eg.db, tab$gene , column= "ENSEMBL", keytype="SYMBOL")
tab$PathwayID <- tab$go_id
tab_filt <- subset(tab, Symbol %in% genes)


## Discard GOs with too few or too many genes
gos_gene_tab <- tab_filt %>%
  group_by(go_id) %>%
  summarize(n = n())

bad_gos <- subset(gos_gene_tab,  n < 10 )$go_id
tab_final <- filter(tab_filt, !go_id %in% bad_gos)


write.table(tab[, c("PathwayID", "Symbol")], file = "results/preprocess/go_gene_map.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

tab2 <- read.table(file = "results/preprocess/kegg_filt_manual_gene_map.tsv", header = TRUE)
tab_com <- rbind(tab_final[, c("PathwayID", "Symbol")], tab2[, c("PathwayID", "Symbol")])

write.table(tab_com, file = "results/preprocess/go_kegg_gene_map.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

