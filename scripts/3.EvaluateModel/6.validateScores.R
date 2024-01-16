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

## Figure 1G
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


cor.test(glyco_df$GSAS1, glyco_df$Glycolisis)
cor.test(glyco_df$GSAS2, glyco_df$Glycolisis)


## Figure 1H
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

## Sup Figure 15
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


