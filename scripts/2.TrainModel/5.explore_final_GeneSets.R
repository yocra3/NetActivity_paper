#'#################################################################################
#'#################################################################################
#' Explore gene sets included in the model
#'#################################################################################
#'#################################################################################

library(GOfuncR)
library(tidyverse)
library(parallel)
library(rjson)
library(NetActivityData)

## Load final pathways
ori_paths <- read.table("results/GTEx_coding/paths_all_full_v3.11/model_trained/pathways_names.txt", header = TRUE)$X0
final_paths <- read.table("results/GTEx_coding/paths_filt2_full_v3.11/model_trained/pathways_names.txt", header = TRUE)$X0
genes <- as.character(read.table("./results/TCGA_gexp_combat_coding/input_genes.txt")$V1)

# GO terms
## Load new GO graph
godir <-  "data/GO_terms/"
term <- read.table(paste0(godir, "/term.txt"), sep = "\t", quote = "", comment.char = "", as.is = TRUE)
graph <- read.table(paste0(godir, "/graph_path.txt"), sep = "\t", quote = "", comment.char = "", as.is = TRUE)

## Get all GO terms
all_gos <- get_child_nodes("GO:0008150", term, graph)

## Remove obsolete terms
term_filt <- subset(term, V5 == 0)
all_gos <- subset(all_gos, child_go_id %in% term_filt$V4)
all_gos$Model <- ifelse(all_gos$child_go_id %in% final_paths, "Included", "Excluded")
all_gos$Original <- ifelse(all_gos$child_go_id %in% ori_paths, "Initial", "Filtered")

genes_pairs <- get_anno_genes(go_ids = all_gos$child_go_id, term_df = term, graph_path_df = graph)

tab <- left_join(mutate(all_gos, go_id = child_go_id), genes_pairs, by = "go_id") %>%
  as_tibble() %>%
  filter(!is.na(gene))

tab$Symbol <- mapIds(org.Hs.eg.db, tab$gene , column= "ENSEMBL", keytype="SYMBOL")
tab$PathwayID <- tab$go_id
tab_filt <- subset(tab, Symbol %in% genes)

gos_filt <- subset(all_gos, child_go_id %in% unique(tab_filt$child_go_id))

top_terms <- subset(gos_filt, distance == 1)$child_go_id
top_terms <- top_terms[!top_terms %in% c("GO:0048518", "GO:0048519")]
top_terms_name <- subset(all_gos, child_go_id %in% top_terms)$child_name
names(top_terms_name) <- top_terms_name

main_gos <- subset(gos_filt, distance > 1)
get_top_parent_go <- function(go){
  print(go)
  parent <- tryCatch(get_parent_nodes(go), error = function(e) NULL)
  if (is.null(parent)){
    return(NA)
  }
  top_parent <- subset(parent, parent_go_id %in% top_terms)
  top_parent <- subset(top_parent, distance == min(top_parent$distance))
  return(top_parent$parent_name)
}
main_gos$top_parent  <- lapply(main_gos$child_go_id, get_top_parent_go)

top_included <- lapply(top_terms_name, function(x) sapply(main_gos$top_parent, function(y) x %in% y)) %>%
      data.frame()



main_gos_df <-  cbind(main_gos, top_included) %>%
  gather(Top_path, Included, 8:30) %>%
  as_tibble() %>%
  filter(Included)
#
top_summary <- group_by(main_gos_df, Top_path ) %>%
  summarize(Total = n(), Initial = sum(Original == "Initial"), Final = sum(Model == "Included")) %>%
  mutate(Prop_Total = Final/Total,
        Prop_Initial = Final/Initial,
        Prop_Initial_Total = Initial/Total) %>%
        arrange(desc(Total)) %>%
        mutate(Top_path = gsub(".", " ", Top_path, fixed = TRUE))
top_summary <- top_summary[, c(1:3, 7, 4:6)]
write.table(top_summary, file = "figures/sel_pathways_GO_topCats.txt", quote = FALSE,
  col.names = TRUE, row.names = FALSE, sep = "\t")

#
# top_paths <- tail(names(sort(table(main_gos_df$Top_path))), 10)
#
# a <- subset(top_summary, Final > 0)
# cor(a$Prop_Initial_Total, a$Prop_Initial, use = "complete")
#
#
# png("figures/cot.png")
# ggplot(top_summary, aes(x = Prop_Initial, y = Prop_Initial_Total, size = log10(Total))) +
#   geom_point() +
#   theme_bw()
# dev.off()
#
# png("figures/sel_pathways_GO_top_distr.png", width = 2000)
# filter(main_gos_df, Top_path %in% top_paths) %>%
#   group_by(Top_path ) %>%
#   summarize(Total = n(), Initial = sum(Original == "Initial"), Final = sum(Model == "Included")) %>%
#   gather(Dataset, N, 2:4) %>%
#   mutate(Dataset = factor(Dataset, levels = c("Total", "Initial", "Final"))) %>%
#   ggplot(aes(x = Dataset, y = N, fill = Dataset)) +
#   geom_bar(stat = "identity") +
#   theme_bw() +
#   facet_wrap(~ Top_path)
# dev.off()
#
# prop.table(data.matrix(subset(top_summary, Top_path %in% top_paths)[, 4:5]), margin = 1)
# chisq.test(subset(top_summary, Top_path %in% top_paths)[-c(4, 8), 4:5])
#
# #
# group_by(main_gos_df, distance ) %>%
#   summarize(Total = n(), Initial = sum(Original == "Initial"), Final = sum(Model == "Included")) %>%
#   data.frame()
#

## KEGG
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

#
kegg_top_summary <- kegg.df %>%
  mutate(Model = ifelse(pathID %in% final_paths, "Included", "Excluded"),
        Original = ifelse(pathID %in% ori_paths, "Initial", "Filtered")) %>%
  group_by(top_cat  ) %>%
  summarize(Total = n(), Initial = sum(Original == "Initial"), Final = sum(Model == "Included")) %>%
  mutate(Prop_Total = Final/Total,
        Prop_Initial = Final/Initial,
        Prop_Initial_Total = Initial/Total) %>%
        arrange(desc(Final)) %>%
        filter(Initial > 0)
#
kegg_top_summary <- kegg_top_summary[, c(1:3, 7, 4:6)]
write.table(kegg_top_summary, file = "figures/sel_pathways_KEGG_topCats.txt", quote = FALSE,
  col.names = TRUE, row.names = FALSE, sep = "\t")

## Create Sup Table
data(gtex_gokegg_annot)

main_gos$category <- apply(top_included, 1, function(i) paste(sort(colnames(top_included)[i]), collapse = ";"))
main_gos$GeneSet <- main_gos$child_go_id
model_gos <- subset(main_gos, GeneSet %in% rownames(gtex_gokegg_annot))

kegg.df$GeneSet <- kegg.df$pathID 
model_kegg <- subset(kegg.df, GeneSet %in% rownames(gtex_gokegg_annot))

model_feats <- rbind(dplyr::select(model_gos, GeneSet, category), dplyr::select(model_kegg, GeneSet, category))

gtex_gokegg_sup <- left_join(gtex_gokegg_annot, model_feats, by = "GeneSet") %>%
  mutate(N = lengths(Weights)) %>%
  dplyr::select(-starts_with("Weights")) %>%
  mutate(Ontology = ifelse(grepl("GO", GeneSet), "GO", "KEGG"),
        Category = gsub(".", " ", category, fixed = TRUE)) %>%
  dplyr::select(GeneSet, Term, Ontology, Category, N)
write.table(gtex_gokegg_sup, file = "figures/sel_pathways_annot.txt", quote = FALSE,
  col.names = TRUE, row.names = FALSE, sep = "\t")

