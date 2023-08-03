#'#################################################################################
#'#################################################################################
#' Evaluate biological relevance of NetActivity weights
#'#################################################################################
#'#################################################################################


# Load libraries
library(tidyverse)
library(cowplot)
library(NetActivityData)
library(AnnotationDbi)
library(org.Hs.eg.db)


## Load gnomad constraints ####
constraint <- read_table("data/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz")
constraint$OMIM <- mapIds(org.Hs.eg.db,
                             keys = constraint$transcript_type,
                             column = 'OMIM',
                             keytype = 'ENSEMBL')
## Load disgenet
disgenet <- read_delim("data/gene_associations.tsv.gz") %>%
  mutate(Ndisease = as.numeric(NofDiseases))
disgenet$gene <- mapIds(org.Hs.eg.db,
                          keys = disgenet$geneSymbol,
                          column = 'ENSEMBL',
                          keytype = 'SYMBOL')


## Process weights ####
data("gtex_gokegg")

weigth_trans <- abs(gtex_gokegg)/rowSums(abs(gtex_gokegg))

weight_df <- tibble(gene = colnames(gtex_gokegg), 
                        Npaths = colSums(gtex_gokegg != 0),
                        relevance = colSums(weigth_trans)) %>%
  left_join(constraint, by = join_by(gene == transcript_type)) %>%
  left_join(disgenet, by = "gene") %>%
  mutate(pLIcat = cut(pLI, c(-1, 0.1, 0.9, 2), c("Low", "Medium", "High")),
         pLIcat = factor(pLIcat, levels = c("Low", "Medium", "High")),
         Npathscat = cut(Npaths, c(0, 1, 2, 5, 60), c("1", "2", "3-5", "5+")),
         Ndisease = ifelse(is.na(Ndisease), 0, Ndisease))


## pLI ####
### Figure 1F
pli_plot <- filter(weight_df, !is.na(pLI)) %>%
  ggplot(aes(x = Npathscat, y = relevance, color = pLIcat)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Gene Sets per gene") +
  ylab("NetActivity relevance") +
  scale_color_discrete(name = "Gene constraint")

png("figures/weights_pLI.png", height = 1200, width = 2000, res = 300)
pli_plot
dev.off()

weight_df %>%
  filter(Npaths > 4) %>%
ggplot(aes(x = relevance, y = pLI, color = pLIcat)) +
  geom_point() +
  theme_bw()


logit <- function(x) log10(x/(1 - x))
weight_df$pLI_logit <- ifelse(weight_df$pLI < 1e-3, 1e-3, weight_df$pLI )
weight_df$pLI_logit <- ifelse(weight_df$pLI_logit > 1 - 1e-3, 1 - 1e-3, weight_df$pLI_logit )
weight_df$pLI_logit <- logit(weight_df$pLI_logit)

summary(lm(pLI_logit ~ Npaths + relevance, weight_df))
summary(lm(pLI_logit ~ relevance, weight_df, subset = Npaths == 1))
summary(lm(pLI_logit ~ relevance, weight_df, subset = Npaths == 2))
summary(lm(pLI_logit ~ Npaths + relevance, weight_df, subset = Npaths > 2 & Npaths < 6))
summary(lm(pLI_logit ~ Npaths + relevance, weight_df, subset = Npaths > 6))

## Disgenet ####
### Figure 1E
disgenet_plot <- ggplot(weight_df, aes(x = Ndisease, y = relevance)) +
  geom_point() +
  scale_x_log10() +
  theme_bw() +
  xlab("Diseases per gene") +
  ylab("NetActivity relevance") +
  facet_grid(~ Npathscat)
png("figures/weights_disgenet.png", height = 1200, width = 2000, res = 300)
disgenet_plot
dev.off()

summary(lm(log10(Ndisease + 0.1) ~ Npaths + relevance, weight_df))
summary(lm(log10(Ndisease  + 0.1) ~ relevance, weight_df, subset = Npaths == 1))
summary(lm(log10(Ndisease  + 0.1) ~ relevance, weight_df, subset = Npaths == 2))
summary(lm(log10(Ndisease + 0.1) ~ Npaths + relevance, weight_df, subset = Npaths > 2 & Npaths < 6))
summary(lm(log10(Ndisease  + 0.1) ~ Npaths + relevance, weight_df, subset = Npaths > 6))

