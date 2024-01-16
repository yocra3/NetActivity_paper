#'#################################################################################
#'#################################################################################
#' Evaluate correlation between GSE179848 measurements and scores
#'#################################################################################
#'#################################################################################

## Load libraries
library(tidyverse)
library(SummarizedExperiment)
library(NetActivity)
library(cowplot)

## Load data
load("results/GSE179848/GSE179848_scores.Rdata")

## Correlate scores with measurements
phenos <- colData(gse179848_scores)
### Select phenotypes with more than 120 values
phenos <- phenos[, colSums(!is.na(colData(gse179848_scores))) > 120]

### Select numeric phenotypes 
phenos <- phenos[, sapply(phenos, class) == "numeric"]

### Remove not interesting phenotypes
bad_phenos <- c("Sample", "Donor_age", "Study_part", "Replicate_line", "Percent_oxygen", "Passage", "sizeFactor", "Days_grown_Udays", 
    "Days_after_treatment_Udays", "Seahorse_run_number", "pre_study_passages_Udivisions")
bad_phenos <- c(bad_phenos, colnames(phenos)[grep("RNAseq", colnames(phenos))], 
    colnames(phenos)[grep("DNAm", colnames(phenos))],  colnames(phenos)[grep("clock", colnames(phenos))],
    colnames(phenos)[grep("Days_Grown", colnames(phenos))])

phenos <- phenos[, !colnames(phenos) %in% bad_phenos]


cor_doub <- cor(t(assay(gse179848_scores)), data.matrix(phenos[, "Doubling_rate_Udays_per_division"]), use = "pairwise.complete.obs")

doup_tib <- tibble(cors = cor_doub, GeneSet = rownames(gse179848_scores), 
    GeneSetName = rowData(gse179848_scores)$Term)
arrange(doup_tib, desc(abs(cors))) %>% data.frame() %>% head(20)

doub_df <- data.frame(GSAS1 = as.vector(assay(gse179848_scores["GO:1902969", ])),
            GSAS2 = as.vector(assay(gse179848_scores["GO:0060623", ])),
            Division_Rate = gse179848_scores$Doubling_rate_Udays_per_division, 
            Treatment = gse179848_scores$Treatments,
            Cell = gse179848_scores$Cell_line)  %>%
            filter(!is.na(Division_Rate)) %>%
            filter(!is.na(Cell)) %>%
            mutate(Cell_type = ifelse(Cell %in% paste0("HC", 1:4), "Control", "Mutated"))
cor.test(doub_df$GSAS1, doub_df$Division_Rate)
## r = 0.56 p-value < 2.2e-16

p1 <- ggplot(doub_df, aes(x = GSAS1, y = Division_Rate)) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_bw() +
    xlab("mitotic DNA replication GSAS") +
    ylab("Doubling Rate")

png("figures/validation_score_doubling.png", height = 1500, width = 1500, res = 300)
p1
dev.off()


## Glycolisis
glyc <- colnames(pheno_cors)[grep("gly", colnames(pheno_cors))]
pheno_cors[glyc, glyc]

cor_glyc <- cor(t(assay(gse179848_scores)), data.matrix(phenos[, "ATPglyc_UpmolATP_per_min_per_20kcells"]), use = "complete.obs")
metab_glyc <- tibble(cors = cor_glyc[, 1], GeneSet = rownames(gse179848_scores), 
    GeneSetName = rowData(gse179848_scores)$Term)
arrange(metab_glyc, desc(abs(cors))) %>% data.frame() %>% head(10)

glyco_df <- data.frame(GSAS1 = as.vector(assay(gse179848_scores["GO:0051156", ])),
            GSAS2 = as.vector(assay(gse179848_scores["GO:0006007", ])),
            GSAS3 =  as.vector(assay(gse179848_scores["path:hsa00020", ])),
            Glycolisis = gse179848_scores$ATPglyc_UpmolATP_per_min_per_20kcells, 
            TotalATP = gse179848_scores$ATPtotal_UpmolATP_per_min_per_20kcells, 
            Treatment = factor(gse179848_scores$Treatments),
            Cell = gse179848_scores$Cell_line)  %>%
            filter(!is.na(Glycolisis)) %>%
            mutate(Cell_type = ifelse(Cell %in% paste0("HC", 1:4), "Control", "Mutated"),
                Treatment = fct_relevel(Treatment, "Control"),
                State = ifelse(Treatment == "Control", "Basal", 
                    ifelse(grepl("^Oligo|mito", Treatment), "OXPHOS inhibition", 
                        ifelse(Treatment %in% c("2DG", "Galactose"), "Glycolisis inhibition", "Other"))),
                State = ifelse(Cell_type == "Mutated", "OXPHOS mutation", State),
                GlycProp = Glycolisis/TotalATP)

glyco_df$logGlyc <- log(glyco_df$Glycolisis)

states <- unique(glyco_df$State)
names(states) <- states
lapply(states, function(state){
    summary(lm(GSAS1 ~ Glycolisis, glyco_df,  subset = State == state))
})
lapply(states, function(state){
    sub <- subset(glyco_df, State == state)
    cor.test(sub$GSAS1, sub$Glycolisis)
})


lapply(states, function(state){
    summary(lm(GSAS2 ~ Glycolisis, glyco_df,  subset = State == state))
})
lapply(states, function(state){
    sub <- subset(glyco_df, State == state)
    cor.test(sub$GSAS2, sub$Glycolisis)
})


lapply(states, function(state){
    summary(lm(GSAS1 ~ logGlyc, glyco_df,  subset = State == state))
})
lapply(states, function(state){
    sub <- subset(glyco_df, State == state)
    cor.test(sub$GSAS1, sub$logGlyc)
})


cor.test(glyco_df$GSAS1, glyco_df$Glycolisis)
cor.test(glyco_df$GSAS2, glyco_df$Glycolisis)
cor.test(glyco_df$GSAS3, glyco_df$Glycolisis)

plot(glyco_df$GSAS2, log(glyco_df$Glycolisis))


p3 <- glyco_df %>%
    filter(State != "Other") %>%
ggplot(aes(x = GSAS1, y = Glycolisis)) +
            geom_point() +
            geom_smooth(method = "lm") +
            xlab("glucose 6-phosphate metabolic process GSAS") +
            ylab("Glycolisis Activity") +
            theme_bw() +
            facet_grid(~ State)

png("figures/validation_score_glycolisis.png", height = 1500, width = 2100, res = 300)
p3
dev.off()

p4 <- glyco_df %>%
    ggplot(aes(x = GSAS2, y = Glycolisis)) +
        geom_point(aes(color = State)) +
        geom_smooth(method = "lm") +
        xlab("glucose catabolic process GSAS") +
        ylab("Glycolisis Activity") +
        theme_bw() 

png("figures/validation_score_glycolisis2.png", height = 1500, width = 1500, res = 300)
p4
dev.off()


summary(lm(GSAS2 ~ Glycolisis + State, glyco_df))
cor.test(glyco_df$GSAS2, glyco_df$Glycolisis)



summary(lm(GSAS1 ~ Glycolisis, glyco_df,  subset = Cell_type == "Mutated" & Treatment == "Control"))
summary(lm(GSAS1 ~ Glycolisis, glyco_df,  subset = Cell_type == "Mutated" & Treatment == "DEX"))
summary(lm(GSAS1 ~ Glycolisis, glyco_df,  subset = Cell_type == "Control" & Treatment == "DEX"))

summary(lm(GSAS1 ~ Glycolisis, glyco_df,  subset = Cell_type == "Control" & Treatment == "Control"))
summary(lm(GSAS1 ~ Glycolisis, glyco_df,  subset = Cell_type == "Control" & Treatment == "Oligomycin"))

summary(lm(Glycolisis ~ GSAS1 + Treatment*Cell_type, glyco_df))


summary(lm(Glycolisis ~ Treatment + Cell, glyco_df, subset = Cell_type == "Control"))



## Oxidative phosphorylation
ox_df <- data.frame(GSAS = as.vector(assay(gse179848_scores["GO:0002082", ])),
            GSAS2 = as.vector(assay(gse179848_scores["path:hsa00020", ])),
            oxPhos = gse179848_scores$ATPox_UpmolATP_per_min_per_20kcells, 
            TotalATP = gse179848_scores$ATPtotal_UpmolATP_per_min_per_20kcells, 
            Treatment = factor(gse179848_scores$Treatments),
            Cell = gse179848_scores$Cell_line)  %>%
            filter(!is.na(oxPhos)) %>%
            mutate(Cell_type = ifelse(Cell %in% paste0("HC", 1:4), "Control", "Mutated"),
                Treatment = fct_relevel(Treatment, "Control"),
                State = ifelse(Treatment == "Control", "Basal", 
                    ifelse(grepl("^Oligo|mito", Treatment), "OXPHOS inhibition", 
                        ifelse(Treatment %in% c("2DG", "Galactose"), "Glycolisis inhibition", "Other"))),
                State = ifelse(Cell_type == "Mutated", "OXPHOS mutation", State),
                OxProp = oxPhos/TotalATP)

ox_df %>%
    ggplot(aes(x = GSAS, y = oxPhos)) +
    geom_point() +
        geom_smooth(method = "lm") +
        xlab("OxPhos GSAS") +
        ylab("OxPhos") +
        facet_grid(~ State) +
        theme_bw() 
ox_df %>%
    ggplot(aes(x = GSAS2, y = oxPhos)) +
    geom_point() +
        geom_smooth(method = "lm") +
        xlab("OxPhos GSAS") +
        ylab("OxPhos") +
        facet_grid(~ State) +
        theme_bw() 

states <- unique(ox_df$State)
names(states) <- states
lapply(states, function(state){
    summary(lm(GSAS ~ oxPhos, ox_df,  subset = State == state))
})
lapply(states, function(state){
    sub <- subset(ox_df, State == state)
    cor.test(sub$GSAS, sub$oxPhos)
})

sub <- subset(ox_df, State %in% c("Basal", "Other", "Glycolisis inhibition"))
cor.test(sub$GSAS2, sub$oxPhos)

### GSVA
library(GSVA)
library(parallel)
library(org.Hs.eg.db)
paths.vec <- rownames(gse179848_scores)
path.map <- read.table("results/GTEx_coding/go_kegg_filt2_gene_map.tsv", header = TRUE)
path_genes <- mclapply(paths.vec, function(x) subset(path.map, PathwayID == x & !is.na(Symbol))$Symbol, mc.cores = 10)
names(path_genes) <- paths.vec

load("results/GSE179848/rawCounts_GSE179848.Rdata")
rownames(gse179848) <- mapIds(org.Hs.eg.db,
    keys = rownames(gse179848),
    column = 'ENSEMBL',
    keytype = 'REFSEQ')

gse179848 <- gse179848[!is.na(rownames(gse179848)), ]
gse179848 <- gse179848[!duplicated(rownames(gse179848)), ]
names(assays(gse179848)) <- "count"
gsva_res <- gsva(gse179848, path_genes,  kcdf = "Poisson")
save(gsva_res, file = "results/GSE179848/gsva_GSE179848.Rdata")

doub_df <- data.frame(GSAS1 = as.vector(assay(gse179848_scores["GO:1902969", ])),
            GSAS2 = as.vector(assay(gse179848_scores["GO:0060623", ])),
            GSAS1g = as.vector(assay(gsva_res["GO:1902969", ])),
            Division_Rate = gse179848_scores$Doubling_rate_Udays_per_division, 
            Treatment = gse179848_scores$Treatments,
            Cell = gse179848_scores$Cell_line)  %>%
            filter(!is.na(Division_Rate)) %>%
            filter(!is.na(Cell)) %>%
            mutate(Cell_type = ifelse(Cell %in% paste0("HC", 1:4), "Control", "Mutated"))
cor.test(doub_df$GSAS1, doub_df$Division_Rate)
cor.test(doub_df$GSAS1g, doub_df$Division_Rate)



glyco_df <- data.frame(GSAS1 = as.vector(assay(gse179848_scores["GO:0051156", ])),
            GSAS2 = as.vector(assay(gse179848_scores["GO:0006007", ])),
            GSAS1g = as.vector(assay(gsva_res["GO:0051156", ])),
            GSAS2g = as.vector(assay(gsva_res["GO:0006007", ])),
            Glycolisis = gse179848_scores$ATPglyc_UpmolATP_per_min_per_20kcells, 
            TotalATP = gse179848_scores$ATPtotal_UpmolATP_per_min_per_20kcells, 
            Treatment = factor(gse179848_scores$Treatments),
            Cell = gse179848_scores$Cell_line)  %>%
            filter(!is.na(Glycolisis)) %>%
            mutate(Cell_type = ifelse(Cell %in% paste0("HC", 1:4), "Control", "Mutated"),
                Treatment = fct_relevel(Treatment, "Control"),
                State = ifelse(Treatment == "Control", "Basal", 
                    ifelse(grepl("^Oligo|mito", Treatment), "OXPHOS inhibition", 
                        ifelse(Treatment %in% c("2DG", "Galactose"), "Glycolisis inhibition", "Other"))),
                State = ifelse(Cell_type == "Mutated", "OXPHOS mutation", State),
                GlycProp = Glycolisis/TotalATP)
cor.test(glyco_df$GSAS2, glyco_df$Glycolisis)
cor.test(glyco_df$GSAS2g, glyco_df$Glycolisis)

cor.test(glyco_df$GSAS1, glyco_df$Glycolisis)
cor.test(glyco_df$GSAS1g, glyco_df$Glycolisis)


ox_df <- data.frame(GSAS = as.vector(assay(gse179848_scores["GO:0002082", ])),
            GSASg = as.vector(assay(gsva_res["GO:0002082", ])),
            oxPhos = gse179848_scores$ATPox_UpmolATP_per_min_per_20kcells, 
            TotalATP = gse179848_scores$ATPtotal_UpmolATP_per_min_per_20kcells, 
            Treatment = factor(gse179848_scores$Treatments),
            Cell = gse179848_scores$Cell_line)  %>%
            filter(!is.na(oxPhos)) %>%
            filter(oxPhos < 500) %>%
            mutate(Cell_type = ifelse(Cell %in% paste0("HC", 1:4), "Control", "Mutated"),
                Treatment = fct_relevel(Treatment, "Control"),
                State = ifelse(Treatment == "Control", "Basal", 
                    ifelse(grepl("^Oligo|mito", Treatment), "OXPHOS inhibition", 
                        ifelse(Treatment %in% c("2DG", "Galactose"), "Glycolisis inhibition", "Other"))),
                State = ifelse(Cell_type == "Mutated", "OXPHOS mutation", State),
                OxProp = oxPhos/TotalATP)


cor.test(ox_df$GSAS, ox_df$oxPhos)
cor.test(ox_df$GSASg, ox_df$oxPhos)

sub <- subset(ox_df, State %in% c("Basal", "Other", "Glycolisis inhibition"))
cor.test(sub$GSAS, sub$oxPhos)
cor.test(sub$GSASg, sub$oxPhos)


ox_df %>%
    ggplot(aes(x = GSAS, y = oxPhos)) +
    geom_point() +
        geom_smooth(method = "lm") +
        xlab("OxPhos GSAS") +
        ylab("OxPhos") +
        facet_grid(~ State) +
        theme_bw() 


ox_df %>%
    ggplot(aes(x = GSASg, y = oxPhos)) +
    geom_point() +
        geom_smooth(method = "lm") +
        xlab("OxPhos GSVA") +
        ylab("OxPhos") +
        facet_grid(~ State) +
        theme_bw() 

par(mfrow = c(1, 2))
plot(ox_df$GSAS, ox_df$oxPhos, main = "NetActivity")
plot(ox_df$GSASg, ox_df$oxPhos, main = "GSVA")

## hipathia
hip_pathways <- load_pathways(species = "hsa")
translated_names <- get_path_names(hip_pathways, names(hip_pathways$path.norm))

