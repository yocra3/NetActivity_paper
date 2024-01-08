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
pheno_cors <- cor(data.matrix(phenos), use = "pairwise.complete.obs")
## Correlate with scores
score_cors <- cor(t(assay(gse179848_scores)), data.matrix(phenos), use = "pairwise.complete.obs")
top_cors <- apply(score_cors, 1, function(x) c(x[which.max(abs(x))], colnames(phenos)[which.max(abs(x))]))

top_cors_tib <- tibble(cors = as.numeric(top_cors[1, ]), GeneSet = colnames(top_cors), Var = top_cors[2, ]) %>%
    mutate(GeneSetName = rowData(gse179848_scores)[GeneSet, "Term"])

arrange(top_cors_tib, desc(abs(cors))) %>% data.frame() %>% head(10)

arrange(top_cors_tib, desc(abs(cors))) %>% data.frame() %>% head(40)
plot(assay(gse179848_scores["GO:1902969", ]), gse179848_scores$Doubling_rate_Udays_per_division)


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
doub_sub <- subset(doub_df, Cell_type == "Control")
cor.test(doub_df$GSAS1, doub_df$Division_Rate)
## r = 0.56 p-value < 2.2e-16


p1 <- ggplot(doub_df, aes(x = GSAS1, y = Division_Rate)) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_bw() +
    xlab("mitotic DNA replication GSAS") +
    ylab("Doubling Rate")

png("figures/validation_score_doubling.png", height = 1500, width = 900, res = 300)
p1
dev.off()


lapply(cells, function(cell){
    summary(lm(GSAS2 ~ Division_Rate, doub_df,  subset = Cell == cell))
})

lapply(cells, function(cell){
    sub <- subset(doub_df, Cell == cell)
    cor.test(sub$GSAS2, sub$Division_Rate)
})
cor.test(doub_df$GSAS2, doub_df$Division_Rate)
## r = 0.59

p2 <- ggplot(doub_df, aes(x = GSAS2, y = Division_Rate)) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_bw() +
    xlab("regulation of chromosome condensation GSAS") +
    ylab("Doubling Rate") +
    facet_wrap(~ Cell_type)

## Glycolisis
glyc <- colnames(pheno_cors)[grep("gly", colnames(pheno_cors))]
pheno_cors[glyc, glyc]

cor_glyc <- cor(t(assay(gse179848_scores)), data.matrix(phenos[, "ATPglyc_UpmolATP_per_min_per_20kcells"]), use = "complete.obs")
metab_glyc <- tibble(cors = cor_glyc[, 1], GeneSet = rownames(gse179848_scores), 
    GeneSetName = rowData(gse179848_scores)$Term)
arrange(metab_glyc, desc(abs(cors))) %>% data.frame() %>% head(10)
filter(metab_glyc, GeneSet %in% c("GO:0051156","GO:0006007"))

plot(assay(gse179848_scores["GO:0051156", ]), gse179848_scores$ATPglyc_UpmolATP_per_min_per_20kcells)
plot(assay(gse179848_scores["GO:0051156", ]), gse179848_scores$ATPglyc_max_UpmolATP_per_min_per_20kcells)
plot(assay(gse179848_scores["GO:0051156", ]), gse179848_scores$ATPglyc_spare_UpmolATP_per_min_per_20kcells)
dev.off()

plot(assay(gse179848_scores["GO:0006007", ]), gse179848_scores$ATPglyc_UpmolATP_per_min_per_20kcells)
plot(assay(gse179848_scores["GO:0006007", ]), gse179848_scores$ATPglyc_max_UpmolATP_per_min_per_20kcells)
plot(assay(gse179848_scores["GO:0006007", ]), gse179848_scores$ATPglyc_spare_UpmolATP_per_min_per_20kcells)
dev.off()

score_cors[c("GO:0051156","GO:0006007"), glyc]


glyco_df <- data.frame(GSAS1 = as.vector(assay(gse179848_scores["GO:0051156", ])),
            GSAS2 = as.vector(assay(gse179848_scores["GO:0006007", ])),
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

cells <- unique(glyco_df$Cell)
names(cells) <- cells
lapply(cells, function(cell){
    summary(lm(GSAS1 ~ Glycolisis, glyco_df,  subset = Cell == cell))
})

lapply(cells, function(cell){
    sub <- subset(glyco_df, Cell == cell)
    cor.test(sub$GSAS1, sub$Glycolisis)
})
cor.test(glyco_df$GSAS1, glyco_df$Glycolisis)
# r = 0.28

glyco_ctrl <- subset(glyco_df, Cell_type == "Control")
cor.test(glyco_ctrl$GSAS1, glyco_ctrl$Glycolisis)


states <- unique(glyco_df$State)
names(states) <- states
lapply(states, function(state){
    summary(lm(GSAS1 ~ Glycolisis, glyco_df,  subset = State == state))
})
lapply(states, function(state){
    sub <- subset(glyco_df, State == state)
    cor.test(sub$GSAS1, sub$Glycolisis)
})

p3 <- ggplot(glyco_df, aes(x = GSAS1, y = Glycolisis, col = Cell)) +
            geom_point() +
            geom_smooth(method = "lm") +
            xlab("glucose 6-phosphate metabolic process GSAS") +
            ylab("Glycolisis Activity") +
            theme_bw() +
            facet_grid(Treatment ~ Cell_type, scales = "free")

glyco_df %>%
    filter(State != "Other") %>%
ggplot(aes(x = GSAS1, y = Glycolisis, col = State)) +
            geom_point() +
            geom_smooth(method = "lm") +
            xlab("glucose 6-phosphate metabolic process GSAS") +
            ylab("Glycolisis Activity") +
            theme_bw() 
dev.off()

summary(lm(GSAS1 ~ Glycolisis, glyco_df,  subset = Cell_type == "Mutated" & Treatment == "Control"))
summary(lm(GSAS1 ~ Glycolisis, glyco_df,  subset = Cell_type == "Mutated" & Treatment == "DEX"))
summary(lm(GSAS1 ~ Glycolisis, glyco_df,  subset = Cell_type == "Control" & Treatment == "DEX"))

summary(lm(GSAS1 ~ Glycolisis, glyco_df,  subset = Cell_type == "Control" & Treatment == "Control"))
summary(lm(GSAS1 ~ Glycolisis, glyco_df,  subset = Cell_type == "Control" & Treatment == "Oligomycin"))

summary(lm(Glycolisis ~ GSAS1 + Treatment*Cell_type, glyco_df))


summary(lm(Glycolisis ~ Treatment + Cell, glyco_df, subset = Cell_type == "Control"))


lapply(cells, function(cell){
    summary(lm(GSAS2 ~ Glycolisis, glyco_df,  subset = Cell == cell))
})

lapply(cells, function(cell){
    sub <- subset(glyco_df, Cell == cell)
    cor.test(sub$GSAS2, sub$Glycolisis)
})
cor.test(glyco_df$GSAS2, glyco_df$Glycolisis)
# r = 0.21
cor.test(glyco_ctrl$GSAS2, glyco_ctrl$Glycolisis)

lapply(states, function(state){
    summary(lm(GSAS2 ~ Glycolisis, glyco_df,  subset = State == state))
})
lapply(states, function(state){
    sub <- subset(glyco_df, State == state)
    cor.test(sub$GSAS2, sub$Glycolisis)
})

glyco_df %>%
    filter(State != "Other") %>%
ggplot(aes(x = GSAS2, y = Glycolisis, col = State)) +
            geom_point() +
            geom_smooth(method = "lm") +
            xlab("glucose 6-phosphate metabolic process GSAS") +
            ylab("Glycolisis Activity") +
            theme_bw() 
dev.off()


p4 <- ggplot(glyco_df, aes(x = GSAS2, y = Glycolisis)) +
            geom_point() +
            geom_smooth(method = "lm") +
            xlab("glucose catabolic process GSAS") +
            ylab("Glycolisis Activity") +
            theme_bw() +
            facet_wrap(~ Cell_type)
png("figures/validation_scores.png", height = 500, width = 700)
plot_grid(p1, p2, p3, p4, ncol = 2, labels = "auto")
dev.off()

