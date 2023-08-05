#'#################################################################################
#'#################################################################################
#'                    PROMOTE vis 2 results with NetActivity
#'            Used in Figure 5; Supp. Fig 14 and Supp. table 6
#'                       
#' 1. Pre-Process PROMOTE gene expression data for NetActivity R package
#' 1.1 Load data and filter bone samples
#' 1.2 Transform counts (vst)
#' 1.3 Change gene SYMBOL to ENSEMBL ID 
#'     Removed genes from expression matrix that are not present in NetActivity 
#'     model and duplicated expression of 1SYMBOL=2ENSEMBLE
#' 
#' 2. Use prepareSummarizedExperiment (Standardized and arranging matrix)
#' 
#' 3. Use computeGeneSetScores (Calculate activation scores with NetActivity)
#' 
#' 4. Differential gene set scores analysis
#' 
#' 5. Save results
#' 
#'#################################################################################
#'#################################################################################


# Load libraries ####

library(SummarizedExperiment)
library(dplyr)
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(HDF5Array)
library(limma)
library(NetActivity)
library(NetActivityData)

# Load data ####
message("Loading SE")

ddsSE <- readRDS('data/promoteV2_SE.rds')

## remove non-bone samples Vis2
message("Filtering bone samples")

ddsSE <- subset(ddsSE, select =  Biopsy_site_Visit_2  == "bone")

## VST
vst_out <- varianceStabilizingTransformation(ddsSE)


## UPDATE SYMBOL NAMES
message("Change GENE ANNOTATIONS")

# Matrix with the NetActivity computed Weights (GO pathways X  genes)
data(gtex_gokegg)
str(gtex_gokegg)
model_genes_ENSG <- data.frame(ModelGenes=colnames(gtex_gokegg))

# Get symbol IDs for MODEL

# THIS IS FROM THE LATEST RELEASE v40 ### CAREFUL IF RE-RUNING THIS FUNCTION
model_genes_ENSG$GeneSymbol <- mapIds(x = org.Hs.eg.db,
                                      keys = model_genes_ENSG$ModelGenes,
                                      column = "SYMBOL",
                                      keytype = "ENSEMBL")
# dim(model_genes_ENSG) # [1] 8758    2


# Although prepareSummarizedExperiment will remove them we trim the input 
# matrix with the genes that are included in the model (using HGNC updated IDs)
keep_genes <- intersect(rowData(vst_out)$new_symbol,model_genes_ENSG$GeneSymbol)
# length(keep_genes)  # 8632


## FIRST change rownames and subset the matrix!
rownames(vst_out) <- rowData(vst_out)$new_symbol
trimmed_expr_matrix <- vst_out[keep_genes,]
# dim(trimmed_expr_matrix)

# Check possible errors/mismatches because of duplicated genes
sum(duplicated(model_genes_ENSG$GeneSymbol))

sub_model_genes_ENSG <- model_genes_ENSG[model_genes_ENSG$GeneSymbol%in%keep_genes,]
dim(sub_model_genes_ENSG) # diff of 1

# There is 1 SYMBOL id repeated, matching 2 different ENSEMBL ids -> 
# Will duplicate this row in expression_mat with the two different ENSEMBLE ids
sub_model_genes_ENSG[duplicated(sub_model_genes_ENSG$GeneSymbol),]# SMN1 will raise error in future scripts
rep_gene <- sub_model_genes_ENSG[duplicated(sub_model_genes_ENSG$GeneSymbol),2]
rep_gene_idx <- which(rownames(trimmed_expr_matrix)==rep_gene)
rep_e_gene <- model_genes_ENSG[which(model_genes_ENSG$GeneSymbol==rep_gene),]$ModelGenes

# Change rownames in expression matrix (be careful with order)
idx <- match(rownames(trimmed_expr_matrix), sub_model_genes_ENSG$GeneSymbol)
rownames(trimmed_expr_matrix) <- sub_model_genes_ENSG$ModelGenes[idx]

# Add rep gene
rep_expr <- trimmed_expr_matrix[rep_gene_idx,]
rownames(rep_expr) <- rep_e_gene[!rep_e_gene%in%rownames(trimmed_expr_matrix)]
input_SE <- rbind(trimmed_expr_matrix,rep_expr)


# prepareSummarizedExperiment Netactivity ####
out <- prepareSummarizedExperiment(input_SE, "gtex_gokegg")

message("Computing gene set scores")
# computeGeneSetScores Netactivity ####
scores <- computeGeneSetScores(out, "gtex_gokegg")

# Save data
save(input_SE, scores, file='results/mCRPC_analyses/NetActivity_GTEX_PROMOTEv2.Rdata')

# write.csv(x=assay(scores), file = "results/mCRPC_analyses/promoteV2_scores_GTEx.csv", sep = ",", row.names = T, quote = T, col.names = T)

# Differential gene set scores analysis ####
# LINEAR MODEL WITH LIMMA
mod <- model.matrix(~TTC.adjust_log, colData(scores))
fit <- lmFit(assay(scores), mod) %>% eBayes()
results <- decideTests(fit)

# Test LM
topTabraw <- topTable(fit, coef=2, n = Inf)

# Arrange data.frame
topTab <- topTabraw
topTab$Term <- rowData(scores)[rownames(topTab), "Term"]
topTab$GeneSet <- rownames(topTab)
i <- max(which(topTab$adj.P.Val<0.05))

topTab <- topTab[c(8,7,1,3,4,5)]
colnames(topTab) <- c( "GeneSet", "Term", "log2FC", "t","pvalue", "p.adj.BH")

# Save results ####

message("Saving DGSA results")

save(topTabraw, "results/mCRPC_analyses/promoteV2_DGSA_TTClog_GTEx.Rdata")
write.csv(head(topTab,i),file='results/mCRPC_analyses/SupTable6_promoteV2_diff_paths_TTCajdustlog_GTEx.csv')

message("Done!")



