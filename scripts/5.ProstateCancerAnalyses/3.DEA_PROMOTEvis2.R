#'#################################################################################
#'#################################################################################
#'        PROMOTE vis 2 Differential gene expression analysis and enrichment
#'                Used in Figure 5; Supp. Fig 13 and Supp. table 7,8 and 9
#'                       
#' 1. Load data
#' 2. Filter bone samples 
#' 3. Differential expression analysis with DESeq2
#' 4. Save results
#' 5. Overrepresentation analysis
#' 6. Plot results
#' 
#'#################################################################################
#'#################################################################################

# Load libraries ####
library(ggplot2)
library(viridis)
library(dplyr)
library(DESeq2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(SummarizedExperiment)
library(NetActivity)
library(NetActivityData)
data(gtex_gokegg)

# Load data ####
message("Loading SE")

ddsSE <- readRDS('data/promoteV2_SE.rds')

# remove non-bone samples Vis2 ####
message("Filtering bone samples")

ddsSE <- subset(ddsSE, select =  Biopsy_site_Visit_2  == "bone")
design(ddsSE) # TTC.adjust_log

# DEA with DESeq2 ####
dds <- DESeq(ddsSE)

res_promote <- results(dds, pAdjustMethod = "BH") # Default is "BH" == FDR (Less conservative than Bonferroni)
dds_DGE_results <- res_promote[order(res_promote$padj),] # order by FDR/BH
dds_DGE_results_df <- as.data.frame(dds_DGE_results)
dds_DGE_results_df$Gene <- rownames(dds_DGE_results_df)
dds_DGE_results_df <- dds_DGE_results_df[,c(7,2,4,5,6)]
colnames(dds_DGE_results_df) <- c("Gene","log2FC","stat", "pvalue" , "p.adj.BH")

# Number of Significant genes: 1064-BH
n_sig <- sum(dds_DGE_results_df$p.adj.BH < 0.05, na.rm=T)

# Number of Significant genes abs(log2FC)>1 592-BH
n_sig_log2FC <- sum(dds_DGE_results_df$p.adj.BH < 0.05 & abs(dds_DGE_results_df$log2FC) > 1, na.rm = T)

message(paste(n_sig,"differentially expressed (FDR < 0.05) of which",n_sig_log2FC,"present |log2 Fold Change| > 1"))
head(dds_DGE_results_df)

# Select genes that have log2 fold change higher than 1 and adjusted pvalue (BH) lower than 0.05
DEgenes_DESeq <- dds_DGE_results_df[which(abs(dds_DGE_results_df$log2FC) > log2(2) & dds_DGE_results_df$p.adj.BH  < 0.05),]
dim(DEgenes_DESeq)

# Save results ####

dir.create("results/mCRPC_analyses/")
save(dds_DGE_results, file = 'results/mCRPC_analyses/promoteV2bone_DEallgenes_results_DESeq.Rdata')
write.csv(DEgenes_DESeq,'results/mCRPC_analyses/SupTable7_promoteV2bone_DEgenes_TTC&2LFC_DESeq.csv', row.names = F)

# Overrepresentation analysis of differentially expressed genes in GO terms ####
# Update symbol genes of DE genes (does not match the order)
# DEgenes_DESeq <- read.csv('./NetActivity_results/manuscript_tables/SupTable7_promoteV2bone_DEgenes_TTC&2LFC_DESeq.csv',header = T,row.names = 1)

DEgenes_DESeq2 <- as.character(rowData(ddsSE)$new_symbol[which(names(rowData(ddsSE)$new_symbol)%in%DEgenes_DESeq$Gene)])
# DEgenes_DESeq2 <- as.character(rowData(ddsSE)$new_symbol[which(names(rowData(ddsSE)$new_symbol)%in%rownames(dds_DGE_results_df)[dds_DGE_results_df$p.adj.BH < 0.05])])

# TRADITIONAL DEA (all GO terms) # Universe promote
enrich_all <- enrichGO(DEgenes_DESeq2,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       universe = rowData(ddsSE)$new_symbol, 
                       ont = "BP",
                       minGSSize = 10,
                       maxGSSize = 5000)

enrich_all_df <- enrich_all@result[,-9]
colnames(enrich_all_df) <- c("GeneSet","Term","GeneRatio","BgRatio","pvalue","p.adj.BH","qvalue", "geneID")
# write.csv(enrich_all_df, file = './NetActivity_results/manuscript_tables/SupTable8_enrichment_universe_PROMOTE_DEG.csv',quote = T,row.names = F)
write.csv(enrich_all_df, file = 'results/mCRPC_analyses/SupTable8_enrichment_universe_PROMOTE_DEG.csv',quote = T,row.names = F)

#PLOT
ORA_dotplot <- clusterProfiler::dotplot(clusterProfiler::simplify(enrich_all))+scale_colour_viridis_b()+
  ggtitle("Over-representation Analysis on differentially expressed genes (DESeq2)",
          subtitle = "All biological process GOterms, all sample genes (traditional analysis)")

# png(filename="./NetActivity_plots/manuscript_plots/ORAdotplot_DEgenes_traditional_analysis.png",res = 300,height = 25,width = 25,units = "cm")
png(filename="figures/SupFig13_ORAdotplot_DEgenes_traditional_analysis.png",res = 300,height = 25,width = 25,units = "cm")
ORA_dotplot
dev.off()


# Universe model genes; test intersect DEA and model genes (within universe)
load(file = 'results/mCRPC_analyses/NetActivity_GTEX_PROMOTEv2.Rdata') # input_SE, scores
model_genes <- unique(names(unlist(rowData(scores)$Weights_SYMBOL)))

enrich_all2 <- enrichGO(intersect(DEgenes_DESeq2, model_genes),
                        OrgDb = org.Hs.eg.db,
                        keyType = "SYMBOL",
                        universe = model_genes, 
                        ont = "BP",
                        minGSSize = 10,
                        maxGSSize = 5000)

enrich_all2_df <- enrich_all2@result[,-9]
colnames(enrich_all2_df) <- c("GeneSet","Term","GeneRatio","BgRatio","pvalue","p.adj.BH","qvalue", "geneID")

write.csv(enrich_all2_df, file='results/mCRPC_analyses/SupTable9_enrichment_universe_NetActivity_DEG.csv',quote = T,row.names = F)

message("Done!")
# EOF