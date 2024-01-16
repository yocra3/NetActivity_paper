#'#################################################################################
#'#################################################################################
#' Compare stability of scores
#'#################################################################################
#'#################################################################################

##########################################################################

## Load libraries
library(limma)
library(DESeq2)
library(tidyverse)
library(cowplot)
library(GSVA)
library(HDF5Array)
library(hipathia)
library(org.Hs.eg.db)
library(parallel)
library(rjson)
library(ggVennDiagram)
library(NetActivity)
library(NetActivityData)
library(ggExtra)

## Load data 
load("data/tcga_gexp_combat.Rdata")
genes <- read.table("./results/GTEx_coding/input_genes.txt")

path.map <- read.table("results/GTEx_coding/go_kegg_filt2_gene_map.tsv", header = TRUE)
path_N <- group_by(path.map, PathwayID) %>% summarize(N = n()) %>% mutate(category = PathwayID)

data(gtex_gokegg)
gset_genes <- lapply(seq_len(nrow(gtex_gokegg)), function(x) colnames(gtex_gokegg)[which(gtex_gokegg[x, ] != 0)])



## NetActivity
prad.all <- gexp_tcga_combat[, gexp_tcga_combat$project_id == "TCGA-PRAD"]
ddsSE <- DESeqDataSet(prad.all, design = ~  age_at_index + race  )
vst.prad <- vst(ddsSE, blind=FALSE)

preproc_prad_ctrl <- prepareSummarizedExperiment(vst.prad[, vst.prad$sample_type == "Solid Tissue Normal"], "gtex_gokegg")
scores_prad_ctrl <- computeGeneSetScores(preproc_prad_ctrl, "gtex_gokegg")

preproc_prad_all <- prepareSummarizedExperiment(vst.prad, "gtex_gokegg")
scores_prad_all <- computeGeneSetScores(preproc_prad_all, "gtex_gokegg")
path.cor <- sapply(seq_len(nrow(scores_prad_ctrl)), function(i) 
  cor(t(assay(scores_prad_all[i, scores_prad_all$sample_type == "Solid Tissue Normal"])), t(assay(scores_prad_ctrl[i, ] ))))

vst.prad_control <- vst.prad[, vst.prad$sample_type == "Solid Tissue Normal"]
scores_prad_all.ctrl <- scores_prad_all[, scores_prad_all$sample_type == "Solid Tissue Normal"]


cor_all_netactivity <- lapply(seq_len(length(gset_genes)), function(i) cor(t(assay(vst.prad_control[gset_genes[[i]], ])), t(assay(scores_prad_all.ctrl [i, ]))))
cor_ctrl_netactivity <- lapply(seq_len(length(gset_genes)), function(i)
  cor(t(assay(vst.prad_control[gset_genes[[i]], ])), t(assay(scores_prad_ctrl[i, ]))))
cor_netActivity_internal <- unlist(Map(cor, cor_all_netactivity, cor_ctrl_netactivity, use = "complete.obs"))



## GSVA
path_genes <- mclapply(paths.vec, function(x) subset(path.map, PathwayID == x & !is.na(Symbol))$Symbol, mc.cores = 10)
names(path_genes) <- paths.vec
gsva.all <- gsva(prad.all, path_genes, min.sz=5, max.sz=500, kcdf = "Poisson")
gsva.ctrl <- gsva(prad.all[, prad.all$sample_type == "Solid Tissue Normal"], path_genes, min.sz=5, max.sz=500, kcdf = "Poisson")
gsva.all.ctrl <- gsva.all[, gsva.all$sample_type == "Solid Tissue Normal"]

save(gsva.all, gsva.ctrl, file = "results/TCGA_PRAD/GSVA_allPRAD_values.Rdata")
gsva.cor <- sapply(seq_len(nrow(gsva.all.ctrl)), function(i) cor(t(assay(gsva.all.ctrl[i, ])), t(assay(gsva.ctrl[i,])) ))

cor_all_gsva <- lapply(seq_len(length(gset_genes)), function(i) cor(t(assay(vst.prad_control[gset_genes[[i]], ])), t(assay(gsva.all.ctrl[i, ]))))
cor_ctrl_gsva <- lapply(seq_len(length(gset_genes)), function(i)
  cor(t(assay(vst.prad_control[gset_genes[[i]], ])), t(assay(gsva.ctrl[i, ]))))
cor_gsva_internal <- unlist(Map(cor, cor_all_gsva, cor_ctrl_gsva, use = "complete.obs"))


## hipathia
hip_pathways <- load_pathways(species = "hsa")

trans_data.all <- translate_data(prad.all, "hsa")
exp_data.all <- normalize_data(trans_data.all)
hip.all <- hipathia(exp_data.all, hip_pathways, decompose = FALSE, verbose = TRUE)
hip.all_vals <- get_paths_data(hip.all )

trans_data.ctrl <- translate_data(prad.all[, prad.all$sample_type == "Solid Tissue Normal" ], "hsa")
exp_data.ctrl <- normalize_data(trans_data.ctrl)
hip.ctrl <- hipathia(exp_data.ctrl, hip_pathways, decompose = FALSE, verbose = TRUE)
hip.ctrl_vals <- get_paths_data(hip.ctrl )

save(hip.all_vals, hip.ctrl_vals, file = "results/TCGA_PRAD/hipathia_allPRAD_values.Rdata")

hip.all.ctrl <- hip.all_vals[, hip.all_vals$sample_type == "Solid Tissue Normal"]
hip.cor <- sapply(seq_len(nrow(hip.all.ctrl)), function(i) cor(t(assay(hip.all.ctrl[i, ])), t(assay(hip.ctrl_vals[i,])) ))

## Compute correlations
hip_pathways <- load_pathways(species = "hsa")

whole_genes <- lapply(hip_pathways$pathigraphs, function(x){
  ig <- x$graph
  genes_list <- V(ig)$genesList
  names(genes_list) <- V(ig)$name
  genes_list
})
names(whole_genes) <- NULL
whole_genes <- unlist(whole_genes, recursive = FALSE)
whole_genes_df <- data.frame(node = rep(names(whole_genes), lengths(whole_genes)),
            GeneID = unlist(whole_genes)) %>%
            filter(!is.na(GeneID) & GeneID != "/") %>%
            as_tibble()

all_effector <- lapply(hip_pathways$pathigraphs, function(x){
  eff <- lapply(x$effector.subgraphs, function(x){
    names(V(x))
  })
  data.frame(PathwayID = rep(names(x$effector.subgraphs), lengths(eff)),
              node = unlist(eff))
})
all_effector <- Reduce(rbind, all_effector)

hipathia_map <- right_join(all_effector, whole_genes_df, by = "node") %>%
  as_tibble()

hipathia_genes <- lapply(unique(hipathia_map$PathwayID), function(path) subset(hipathia_map, PathwayID == path)$GeneID)

trans_data.all_ctrl <- trans_data.all[,  trans_data.all$sample_type == "Solid Tissue Normal"]
cor_all_hip <- lapply(seq_len(length(hipathia_genes)),
 function(i) cor(t(assay(trans_data.all_ctrl[rownames(trans_data.all_ctrl) %in% hipathia_genes[[i]], ])),  t(assay(hip.all.ctrl[i, ]))))

cor_ctrl_hip <- lapply(seq_len(length(hipathia_genes)),
 function(i) cor(t(assay(trans_data.ctrl[rownames(trans_data.ctrl) %in% hipathia_genes[[i]], ])),  t(assay(hip.ctrl_vals[i, ]))))


cor_mod <- function(x, y){
  if (sum(!is.na(x)) < 2 | sum(!is.na(y)) < 2){
    return(NA)
  }
  if (length(x) != length(y)){
    return(NA)
  }
  cor(x, y, use = "complete.obs")
}

cor_hip_internal <- unlist(Map(cor_mod, cor_all_hip, cor_ctrl_hip))


df.cor <- tibble(replicability = c(path.cor, gsva.cor, hip.cor),
  consistency = c(cor_netActivity_internal, cor_gsva_internal, cor_hip_internal),
  Method = rep(c("NetActivity", "GSVA", "Hipathia"), lengths(list(path.cor, gsva.cor, hip.cor)))) %>%
  mutate(Method = factor(Method, levels = c("GSVA", "Hipathia", "NetActivity"))) 

#


plot_stab <- df.cor %>%
  mutate(Replicability = cut(replicability, c(0, 0.7, 0.9, 1.1), c("Low/Inter.", "High", "Equivalent")),
         Consistency = cut(consistency, c(-1.1, 0, 0.7, 0.9, 1.1), c("Opposite", "Low/Inter.", "High", "Equivalent"))) %>%
  gather(Measure, Value, 4:5) %>%
  filter(!is.na(Value)) %>%
  mutate(Measure = recode(Measure, Replicability = "Score replicability", Consistency = "Definition consistency"),
         Measure = factor(Measure, levels = c("Score replicability", "Definition consistency"))) %>%
  group_by(Method, Measure, Value) %>% 
  summarize(n = n()) %>%
  group_by(Method, Measure) %>%
  mutate(prop = n/sum(n)) %>%
  ungroup() %>%
  complete(Method, Measure, Value, fill = list(n = 0, prop = 0)) %>%
  filter(!(Measure == "Score replicability" & Value == "Opposite")) %>%
  ggplot(aes(x = Value, fill = Method, y = prop*100)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  xlab("Gene Sets") +
  ylab("Proportion (%)") +
  facet_grid(. ~ Measure, scales = "free") +
  ggtitle("Same dataset") +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(size = 16),
        legend.position = "none")
  
png("figures/PRAD_score_stability.png", height = 1000, width = 2000, res = 300)
plot_stab
dev.off()

summary(lm(cor ~ Method, df.cor))

tapply(df.cor$cor, df.cor$Method, summary)
df.cor %>%
  group_by(Method) %>%
  summarize(m = mean(cor > 0.95, na.rm = TRUE))




## Make panel
load("results/GSE169038/de_genes_results.Rdata")
load("results/GSE169038/GSVA_results.Rdata")
load("results/GSE169038/pathways_results.Rdata")
load("results/GSE169038/hipathia.res.Rdata")



# Comparison between PRAD and GEO
## Pathways
load("results/TCGA_PRAD/pathways_results.Rdata")
comb_paths <- left_join(tab.path_prad, tab.paths_geo, by = "category", suffix = c(".TCGA", ".GEO")) %>%
  as_tibble() %>%
  mutate(Signif = ifelse(adj.P.Val.TCGA < 0.05, ifelse(adj.P.Val.GEO < 0.05, "Both", "TCGA"),
                              ifelse(adj.P.Val.GEO < 0.05, "GEO", "None")))

path.plot <- ggplot(comb_paths, aes(x = logFC.TCGA, y = logFC.GEO, col = Signif)) +
  geom_point() +
  theme_bw()  +
  scale_color_manual(values = c("#004D40", "#1E88E5", "#9E9E9E", "#FFC107")) +
  ggtitle("NetActivity") +
  xlab("logFC in TCGA-PRAD") +
  ylab("logFC in GEO-PRAD") +
  geom_text(data =  data.frame(label = sprintf("N = %d \n r = %.2f", nrow(comb_paths), cor(comb_paths$logFC.TCGA, comb_paths$logFC.GEO)),
   x = -Inf, y = Inf, hjust = -0.3, vjust = 1.5), aes(label = label, x = x, y = y, hjust = hjust, vjust = vjust), col = "black", size = 6) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
    text = element_text(size = 16))



png("figures/TCGAvsGEO_logFC.png")
path.plot
dev.off()

## Add GO names
term <- read.table("data/GO_terms/term.txt", sep = "\t", quote = "", comment.char = "", as.is = TRUE)
term <- term[, c(2, 4)]
colnames(term) <- c("Name", "category")



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
  mutate(category = paste0("path:hsa", substring(path, 0, 5)),
          Name = gsub("^[0-9]*  ", "", path))

path_map <- rbind(term, dplyr::select(kegg.df, Name, category))
comb_paths_annot <- left_join(comb_paths, path_map, by = "category") %>%
  dplyr::select(category, Name, ends_with("TCGA"), ends_with("GEO"), Signif) %>%
  dplyr::select(-starts_with("AveExpr"), -starts_with("t"), -starts_with("B"))
write.table(comb_paths_annot, file = "results/GSE169038/paths_results_comb.txt",
  sep = "\t", quote = FALSE, row.names = FALSE )

#

summary(lm(logFC.TCGA ~ logFC.GEO, comb_paths ))
summary(lm(logFC.TCGA ~ logFC.GEO, comb_paths, subset = Signif != "None" ))

cor(comb_paths$logFC.TCGA, comb_paths$logFC.GEO)
# [1] 0.5047529
cor(comb_paths[comb_paths$Signif != "None", ]$logFC.TCGA, comb_paths[comb_paths$Signif != "None", ]$logFC.GEO)
# [1] 0.7034897
cor(comb_paths[comb_paths$Signif == "Both", ]$logFC.TCGA, comb_paths[comb_paths$Signif == "Both", ]$logFC.GEO)
# [1] 0.9121616

table(sign(comb_paths$logFC.TCGA) == sign(comb_paths$logFC.GEO), comb_paths$Signif)
prop.table(table(sign(comb_paths$logFC.TCGA) == sign(comb_paths$logFC.GEO), comb_paths$Signif), margin = 2)
prop.table(table(sign(comb_paths$logFC.TCGA) == sign(comb_paths$logFC.GEO), comb_paths$Signif != "None"), margin = 2)

png("figures/GSE169038_TCGA_propDE_vs_logFCPath.png", height = 300)
rbind(mutate(tab.path_prad, Dataset = "TCGA"),
    mutate(tab.paths_geo, Dataset = "GEO")) %>%
    mutate(Dataset = factor(Dataset, levels = c("TCGA", "GEO"))) %>%
    ggplot(aes(x = DE_prop , y = abs(logFC))) +
    geom_point() +
    scale_x_continuous(name = "Proportion of genes DE") +
    scale_y_continuous(name = "logFC Gene Set (absolute value)") +
    facet_wrap(~ Dataset) +
    theme_bw()
dev.off()



## Compare Gene DE
load("results/TCGA_PRAD/genes_results.Rdata")
res_prad$gene <- rownames(res_prad)
comb.genes <- left_join(data.frame(res_prad), tab.genes_geo , by = "gene")   %>%
  mutate(Signif = ifelse(!is.na(padj) & padj  < 0.05, ifelse(adj.P.Val < 0.05, "Both", "TCGA"),
                              ifelse(adj.P.Val < 0.05, "GEO", "None"))) %>%
  filter(!is.na(pvalue ) & !is.na(P.Value  ))

gene.plot <- ggplot(comb.genes, aes(x = log2FoldChange, y = logFC, col = Signif)) +
  geom_point() +
  theme_bw()  +
  scale_color_manual(values = c("#004D40", "#1E88E5", "#9E9E9E", "#FFC107")) +
  ggtitle("Genes") +
  xlab("log2FC in TCGA-PRAD") +
  ylab("logFC in GEO-PRAD") +
  geom_text(data =  data.frame(label = sprintf("N = %d \n r = %.2f", nrow(comb.genes), cor(comb.genes$log2FoldChange, comb.genes$logFC)),
   x = -Inf, y = Inf, hjust = -0.3, vjust = 1.5), aes(label = label, x = x, y = y, hjust = hjust, vjust = vjust), col = "black", size = 6) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
    text = element_text(size = 16))

png("figures/TCGAvsGEO_genes_logFC.png")
gene.plot
dev.off()


cor(comb.genes$log2FoldChange, comb.genes$logFC)
# [1] 0.2488165
cor(comb.genes[comb.genes$Signif != "None", ]$log2FoldChange, comb.genes[comb.genes$Signif != "None", ]$logFC, use = "complete")
# [1] 0.3250043

table(sign(comb.genes$log2FoldChange) == sign(comb.genes$logFC), comb.genes$Signif)
prop.table(table(sign(comb.genes$log2FoldChange) == sign(comb.genes$logFC), comb.genes$Signif), margin = 2)
prop.table(table(sign(comb.genes$log2FoldChange) == sign(comb.genes$logFC), comb.genes$Signif != "None"), margin = 2)


## GSEA
load("results/TCGA_PRAD/GSVA_results.Rdata")
comb.gsva <- inner_join(tab.gsva_prad, tab.gsva_geo, by = "category", suffix = c(".TCGA", ".GEO"))   %>%
  mutate(Signif = ifelse(!is.na(adj.P.Val.TCGA) & adj.P.Val.TCGA  < 0.05, ifelse(adj.P.Val.GEO < 0.05, "Both", "TCGA"),
                              ifelse(adj.P.Val.GEO < 0.05, "GEO", "None")))


#
gsva.plot <- ggplot(comb.gsva, aes(x = logFC.TCGA, y = logFC.GEO, col = Signif)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = c("#004D40", "#1E88E5", "#9E9E9E", "#FFC107")) +
  ggtitle("GSVA") +
  xlab("logFC in TCGA-PRAD") +
  ylab("logFC in GEO-PRAD") +
  geom_text(data =  data.frame(label = sprintf("N = %d \n r = %.2f", nrow(comb.gsva), cor(comb.gsva$logFC.TCGA, comb.gsva$logFC.GEO)),
   x = -Inf, y = Inf, hjust = -0.3, vjust = 1.5), aes(label = label, x = x, y = y, hjust = hjust, vjust = vjust), col = "black", size = 6) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
    text = element_text(size = 16))

png("figures/TCGAvsGEO_GSVA_logFC.png")
gsva.plot
dev.off()

cor(comb.gsva$logFC.TCGA, comb.gsva$logFC.GEO)
# [1] 0.6390369
cor(comb.gsva[comb.gsva$Signif != "None", ]$logFC.TCGA, comb.gsva[comb.gsva$Signif != "None", ]$logFC.GEO, use = "complete")
# [1] 0.809432
cor(comb.gsva[comb.gsva$Signif == "Both", ]$logFC.TCGA, comb.gsva[comb.gsva$Signif == "Both", ]$logFC.GEO, use = "complete")
# [1] 0.9415681


table(sign(comb.gsva$logFC.TCGA) == sign(comb.gsva$logFC.GEO), comb.gsva$Signif)
prop.table(table(sign(comb.gsva$logFC.TCGA) == sign(comb.gsva$logFC.GEO), comb.gsva$Signif), margin = 2)
prop.table(table(sign(comb.gsva$logFC.TCGA) == sign(comb.gsva$logFC.GEO), comb.gsva$Signif != "None"), margin = 2)

## hipathia
load("results/TCGA_PRAD/hipathia.res.Rdata")
comb.hipathia <- inner_join(hip.comp_prad, hip.comp_geo, by = "name", suffix = c(".TCGA", ".GEO"))   %>%
  mutate(Signif = ifelse(!is.na(FDRp.value.TCGA) & FDRp.value.TCGA  < 0.05, ifelse(FDRp.value.GEO < 0.05, "Both", "TCGA"),
                              ifelse(FDRp.value.GEO < 0.05, "GEO", "None")))


#
hip.plot <- ggplot(comb.hipathia, aes(x = statistic.TCGA, y = statistic.GEO, col = Signif)) +
  geom_point() +
  theme_bw()  +
  scale_color_manual(name = "Significance", values = c("#004D40", "#1E88E5", "#9E9E9E", "#FFC107")) +
  ggtitle("Hipathia") +
  xlab("U statistic in TCGA-PRAD") +
  ylab("U statistic in GEO-PRAD") +
  geom_text(data =  data.frame(label = sprintf("N = %d \n r = %.2f", nrow(comb.hipathia), cor(comb.hipathia$statistic.TCGA, comb.hipathia$statistic.GEO)),
   x = -Inf, y = Inf, hjust = -0.3, vjust = 1.5), aes(label = label, x = x, y = y, hjust = hjust, vjust = vjust), col = "black", size = 6) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
    text = element_text(size = 16))

png("figures/TCGAvsGEO_hipathia_stat.png")
hip.plot
dev.off()

cor(comb.hipathia$statistic.TCGA, comb.hipathia$statistic.GEO)
# [1] 0.234694
cor(comb.hipathia[comb.hipathia$Signif != "None", ]$statistic.TCGA, comb.hipathia[comb.hipathia$Signif != "None", ]$statistic.GEO, use = "complete")
# [1] 0.4226619

table(sign(comb.hipathia$statistic.TCGA) == sign(comb.hipathia$statistic.GEO), comb.hipathia$Signif)
prop.table(table(sign(comb.hipathia$statistic.TCGA) == sign(comb.hipathia$statistic.GEO), comb.hipathia$Signif), margin = 2)
prop.table(table(sign(comb.hipathia$statistic.TCGA) == sign(comb.hipathia$statistic.GEO), comb.hipathia$Signif != "None"), margin = 2)

legend <- get_legend(
  # create some space to the left of the legend
  hip.plot + theme(legend.box.margin = margin(0, 0, 0, 12),
                    text = element_text(size = 16))+
                    guides(color = guide_legend(override.aes = list(size = 8)))
)



## Compare GSVA vs GO+kegg
# paths.vec2 <- gsub(":", "_",paths.vec)
# cor_measures <- sapply(paths.vec2, function(i) cor( geo.feat.filt[, i], geo_gsva[i, ]))
# names(cor_measures) <- paths.vec

gsva_path_df <- left_join(dplyr::select(comb_paths, category, Signif, starts_with("log")),
                          dplyr::select(comb.gsva, category, Signif, starts_with("log")), by = "category", suffix = c(".path", ".GSVA")) %>%
                          mutate(Signif = ifelse(Signif.path == "Both",
                                                    ifelse(Signif.GSVA == "Both", "Both",
                                                            ifelse(Signif.GSVA == "None", "Both paths - None GSVA", "Both paths - 1 GSVA")),
                                                    ifelse(Signif.path == "None",
                                                            ifelse(Signif.GSVA == "Both", "None paths - Both GSVA",
                                                                    ifelse(Signif.GSVA == "None", "None", "None paths - 1 GSVA")),
                                                              ifelse(Signif.GSVA == "Both", "1 paths - Both GSVA",
                                                                    ifelse(Signif.GSVA == "None", "1 paths - None GSVA", "1 paths - 1 GSVA")))),
                                  Signif.GSVA = factor(Signif.GSVA, levels = c("Both", "GEO", "TCGA", "None")),
                                  Signif.path = factor(Signif.path, levels = c("Both", "GEO", "TCGA", "None"))) %>%
                  left_join(path_N, by = "category")


gsva_path_df$GEO_prop <- sapply(gsva_path_df$category, function(path) {
  sel <- subset(comb.genes, gene %in% subset(path.map, PathwayID == path)$Symbol)
  p <- mean(sign(sel$logFC) == 1)
})
gsva_path_df$TCGA_prop <- sapply(gsva_path_df$category, function(path) {
  sel <- subset(comb.genes, gene %in% subset(path.map, PathwayID == path)$Symbol)
  p <- mean(sign(sel$log2FoldChange) == 1)
})

## Sup Figure 12
png("figures/TCGAvsGEO_GSVA_NetActivity_comp_prop.png", width = 2000, height = 1500, res = 300)
gsva_path_df %>%
    gather(Dataset, Proportion, 11:12) %>%
    mutate(GSVA.Sig = ifelse(Dataset == "GEO_prop",
                        ifelse(Signif.GSVA %in% c("Both", "GEO"), "Significant", "Non-significant"),
                        ifelse(Signif.GSVA %in% c("Both", "TCGA"), "Significant", "Non-significant")),
                      path.Sig = ifelse(Dataset == "GEO_prop",
                        ifelse(Signif.path %in% c("Both", "GEO"), "Significant", "Non-significant"),
                        ifelse(Signif.path %in% c("Both", "TCGA"), "Significant", "Non-significant"))) %>%
    gather(Method, Significance, 13:14) %>%
    mutate(Method = recode(Method, GSVA.Sig = "GSVA", path.Sig = "NetActivity"),
            Dataset = recode(Dataset, GEO_prop = "GEO-PRAD", TCGA_prop = "TCGA-PRAD" ),
            Dataset = factor(Dataset, levels = c("TCGA-PRAD", "GEO-PRAD"))) %>%
    ggplot(aes(color = Significance, x = Proportion)) +
    geom_density() +
    facet_grid(Method ~ Dataset) +
    theme_bw() +
    ylab("Density") +
    xlab("Prop. of genes with higher expression in Gleason high") +
    geom_vline(xintercept = 0.5)
dev.off()


venn_geo <- ggVennDiagram(list(GSVA = subset(gsva_path_df, Signif.GSVA %in% c("Both", "GEO"))$category,
                              NetActivity = subset(gsva_path_df, Signif.path %in% c("Both", "GEO"))$category),
                            set_size = 7, label_size = 7) +
            scale_fill_gradient(low = "#FFFFFF", high = "#FFFFFF") +
            ggtitle("GEO-PRAD") +
            theme(plot.title = element_text(hjust = 0.5, size = 25),
                  legend.position = "none")
#
venn_tcga <- ggVennDiagram(list(GSVA = subset(gsva_path_df, Signif.GSVA %in% c("Both", "TCGA"))$category,
                              NetActivity = subset(gsva_path_df, Signif.path %in% c("Both", "TCGA"))$category),
                            set_size = 7, label_size = 7) +
            scale_fill_gradient(low = "#FFFFFF", high = "#FFFFFF") +
            ggtitle("TCGA-PRAD") +
            theme(plot.title = element_text(hjust = 0.5, size = 25),
                  legend.position = "none")

#
venn_both <- ggVennDiagram(list(GSVA = subset(gsva_path_df, Signif.GSVA %in% c("Both"))$category,
                              NetActivity = subset(gsva_path_df, Signif.path %in% c("Both"))$category),
                            set_size = 7, label_size = 7) +
            scale_fill_gradient(low = "#FFFFFF", high = "#FFFFFF") +
            ggtitle("Both") +
            theme(plot.title = element_text(hjust = 0.5, size = 25),
                  legend.position = "none")

## Sup Figure 11
png("figures/TCGAvsGEO_GOKEGG_overlap.png", width = 1200, height = 300)
plot_grid(venn_tcga, venn_geo, venn_both, nrow = 1, labels = LETTERS[1:3])
dev.off()



## Compare concordance in gene set score definition
load("results/TCGA_PRAD/vst_SE.Rdata")
load("results/TCGA_PRAD/NetActivity_scores.Rdata")
load("results/GSE169038/NetActivity_scores.Rdata")
load("results/GSE169038/GSVA_results.Rdata")

load("results/TCGA_BRCA/vst_SE.Rdata")
load("results/TCGA_BRCA/NetActivity_scores.Rdata")
load("results/TCGA_PRAD/GSVA_results.Rdata")
load("results/TCGA_BRCA/GSVA_results.Rdata")

gse169038_se <- loadHDF5SummarizedExperiment("results/GSE169038/", prefix = "network_genes")
assay(gse169038_se) <- data.matrix(assay(gse169038_se))
gse169038_se$primary <- gsub("primary gleason: ", "", gse169038_se$characteristics_ch1.1)
gse169038_se$secondary <- gsub("secondary gleason: ", "", gse169038_se$characteristics_ch1.2)
gse169038_se$primary <- as.numeric(ifelse(gse169038_se$primary == "--", 1, gse169038_se$primary))
gse169038_se$secondary <- as.numeric(ifelse(gse169038_se$secondary == "--", 1, gse169038_se$secondary))
gse169038_se_filt <- gse169038_se[, !(gse169038_se$primary == 1 |  gse169038_se$secondary == 1)]

data(gtex_gokegg)
gset_genes <- lapply(seq_len(nrow(gtex_gokegg)), function(x) colnames(gtex_gokegg)[which(gtex_gokegg[x, ] != 0)])

prad_input <- vst.prad[colnames(gtex_gokegg), ]
brca_input <- vst_brca[colnames(gtex_gokegg), ]

cor_prad_netactivity <- lapply(seq_len(length(gset_genes)), function(i) cor(t(assay(prad_input[gset_genes[[i]], ])), t(assay(scores_prad[i, ]))))
cor_geo_netactivity <- lapply(seq_len(length(gset_genes)), function(i) cor(t(assay(gse169038_se_filt[gset_genes[[i]], ])), t(assay(gse169038_scores[i, ]))))
cor_brca_netactivity <- lapply(seq_len(length(gset_genes)), function(i) cor(t(assay(brca_input[gset_genes[[i]], ])), t(assay(scores_brca[i, ]))))
cor_netActivity_rnaseq <- unlist(Map(cor, cor_prad_netactivity, cor_brca_netactivity, use = "complete.obs"))
cor_netActivity_prad <- unlist(Map(cor, cor_prad_netactivity, cor_geo_netactivity, use = "complete.obs"))


cor_prad_gsva <- lapply(seq_len(length(gset_genes)), function(i) cor(t(assay(prad_input[gset_genes[[i]], ])),  t(assay(prad_gsva[i, ]))))
cor_geo_gsva <- lapply(seq_len(length(gset_genes)), function(i) cor(t(assay(gse169038_se_filt[gset_genes[[i]], ])), geo_gsva[i, ]))
cor_brca_gsva <- lapply(seq_len(length(gset_genes)), function(i) cor(t(assay(brca_input[gset_genes[[i]], ])), t(assay(brca_gsva[i, ]))))
cor_gsva_rnaseq <- unlist(Map(cor, cor_prad_gsva, cor_brca_gsva, use = "complete.obs"))
cor_gsva_prad <- unlist(Map(cor, cor_prad_gsva, cor_geo_gsva, use = "complete.obs"))

### Hipathia
load("results/TCGA_PRAD/hipathia.res.Rdata")
load("results/TCGA_BRCA/hipathia.res.Rdata")
load("results/GSE169038/hipathia.res.Rdata")

prad_trans <- translate_data(vst.prad, "hsa")
brca_trans <- translate_data(vst_brca, "hsa")
geo_trans <- translate_data(gse169038_se_filt, "hsa")


hip.geo <- hip.geo_path[, colnames(geo_trans) ]

cor_prad_hip <- lapply(seq_len(length(hipathia_genes)),
 function(i) cor(t(assay(prad_trans[rownames(prad_trans) %in% hipathia_genes[[i]], ])),  t(assay(hip.prad_vals[i, ]))))
cor_geo_hip <- lapply(seq_len(length(hipathia_genes)), 
  function(i) cor(t(assay(geo_trans[rownames(geo_trans) %in% hipathia_genes[[i]], ])), t(assay(hip.geo[i, ]))))
cor_brca_hip <- lapply(seq_len(length(hipathia_genes)), function(i) cor(t(assay(brca_trans[rownames(brca_trans) %in% hipathia_genes[[i]], ])), t(assay(hip.brca_vals[i, ]))))

cor_hip_rnaseq <- unlist(Map(cor_mod, cor_prad_hip, cor_brca_hip))
cor_hip_prad <- unlist(Map(cor_mod, cor_prad_hip, cor_geo_hip))

df_reprod <- tibble(Technology = c(cor_netActivity_prad, cor_gsva_prad, cor_hip_prad),
                    Tumor = c(cor_netActivity_rnaseq, cor_gsva_rnaseq, cor_hip_rnaseq),
                    Method = rep(c("NetActivity", "GSVA", "Hipathia"), lengths(list(cor_netActivity_prad, cor_gsva_prad, cor_hip_prad))),
                    GeneSet = c(rownames(scores_prad), rownames(prad_gsva), rownames(hip.prad_vals))) %>%
              gather(Dataset, Reproducibility, 1:2) %>%
              filter(!is.na(Reproducibility)) %>%
              mutate(Replicability = cut(Reproducibility, c(-1.1, 0, 0.7, 0.9, 1.1), c("Opposite", "Low/Inter.", "High", "Equivalent")),
                     Replicability = factor(Replicability, levels = c("Equivalent", "High", "Low/Inter.", "Opposite"))) %>%
              group_by(Method, Dataset, Replicability) %>% 
              summarize(n = n()) %>%
              group_by(Method, Dataset) %>%
              mutate(prop = n/sum(n)) %>%
              ungroup() %>%
              complete(Method, Dataset, Replicability, fill = list(n = 0, prop = 0))

plot_reprod <- ggplot(df_reprod, aes(x = Replicability, fill = Method, y = prop*100)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  xlab("") +
  ylab("Proportion (%)") +
  facet_grid(. ~ Dataset, scales = "free") +
  ggtitle("Different datasets") +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(size = 16))
  
## Figure 2
fig2_panel <- plot_grid(
  plot_grid(plot_stab, plot_reprod, nrow = 1, labels = c("A", "B"), label_size = 16, rel_widths = c(1, 1.3)),
  plot_grid(
    plot_grid(gene.plot, gsva.plot, hip.plot + theme(legend.position = "none"), path.plot, ncol = 2, labels = LETTERS[3:6], label_size = 16),
    legend, ncol = 2, rel_widths = c(5, 1)
  )
, ncol = 1, rel_heights = c(1, 3))
png("figures/TCGAvsGEO_panel.png", width = 4000, height = 3200, res = 300)
fig2_panel
dev.off()
ggsave("figures/Figure2.eps", plot = fig2_panel, device = "eps", width = 4000, height = 3200, units = "px")



df_reprod %>%
  group_by(Method, Dataset) %>%
  summarize(m = mean(Reproducibility > 0, na.rm = TRUE),
            m30 = mean(Reproducibility > 0.3, na.rm = TRUE),
            m70 = mean(Reproducibility > 0.7, na.rm = TRUE))

tapply(cor_netActivity_prad[comb_paths$category], comb_paths$Signif, summary)
tapply(cor_gsva_prad[comb.gsva$category], comb.gsva$Signif, summary)

comb_paths$repr <- cor_netActivity_prad[comb_paths$category]
