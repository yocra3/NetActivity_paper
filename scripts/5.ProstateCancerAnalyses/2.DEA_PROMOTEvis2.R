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
library(ggpattern)
library(ggrepel)
library(dplyr)
library(readxl)
library(DESeq2)
library(AnnotationDbi)
library(xlsx)
library(org.Hs.eg.db)
library(clusterProfiler)

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
# Manuscript
dir.create("results/mCRPC_analyses/")
save(dds_DGE_results,'results/mCRPC_analyses/promoteV2bone_allgenes_results_DESeq.Rdata')
write.csv(DEgenes_DESeq,'results/mCRPC_analyses/SupTable7_promoteV2bone_DEgenes_TTC&2LFC_DESeq.csv',row.names = F,col.names = T)

# Overrepresentation analysis of differentially expressed genes in GO terms ####
# Ipdate symbol genes of DE genes (does not match the order)
DEgenes_DESeq2 <- rowData(ddsSE)$new_symbol[which(names(rowData(ddsSE)$new_symbol)%in%DEgenes_DESeq)]

# TRADITIONAL DEA (all GO terms) # Universe promote
enrich_all <- enrichGO(DEgenes_DESeq2,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       universe = rowData(ddsSE)$new_symbol, 
                       ont="BP",
                       minGSSize = 10,
                       maxGSSize = 5000)
write.csv(enrich_all@result,file='./NetActivity_results/SupTable8_enrichment_universe_PROMOTE_DEG.csv',sep = ",",row.names = T,quote = T,col.names = T)

# Universe model genes; test intersect DEA and model genes (within universe)
model_genes <- unique(names(unlist(rowData(scores)$Weights_SYMBOL)))
enrich_all2 <- enrichGO(intersect(DEgenes_DESeq2, model_genes),
                        OrgDb = org.Hs.eg.db,
                        keyType = "SYMBOL",
                        universe = model_genes, 
                        ont="BP",
                        minGSSize = 10,
                        maxGSSize = 5000)
write.csv(enrich_all2@result,file='./NetActivity_results/SupTable9_enrichment_universe_NetActivity_DEG.csv',sep = ",",row.names = T,quote = T,col.names = T)



dotp3 <- clusterProfiler::dotplot(clusterProfiler::simplify(enrich_all2),showCategory=30)+scale_colour_viridis_b()+
  ggtitle("Differential Expression Analysis",subtitle = "All terms (traditional analysis)")

enrich_all_simplify <- clusterProfiler::simplify(enrich_all)
enrich_all_subset <- enrich_all_simplify@result %>% filter(p.adjust<0.05)
myGeneRatio <- sapply(enrich_all_subset$GeneRatio , function(x){
  ratio <- as.numeric(stringr::str_split(x,"/")[[1]][1])/as.numeric(stringr::str_split(x,"/")[[1]][2])
})

