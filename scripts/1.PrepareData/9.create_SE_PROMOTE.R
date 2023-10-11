#'#################################################################################
#'#################################################################################
#'                    
#'         Prepare mCRPC PROMOTE vis2 SummarizedExperiment object 
#'                      with updated gene SYMBOL names                    
#' 1. Load files: RAW COUNTS and CLINICAL
#' 2. Convert clinical into appropriated factor variables
#' 3. Add updated symbol names
#' 4. Create SE and save
# 
#'#################################################################################
#'#################################################################################

# Load libraries ####
rm(list=ls())
library(readxl)
library(dplyr)
library(SummarizedExperiment)
library(DESeq2)

# Helper Function ####

get_old_symbol <- function(gene=NULL,symbol_db_modx=NULL){
  # Gene to search
  pattern <- paste0("^",gene,"$")
  # Loop over prev_symbols of symbol_db (rows) string_splitted by "|" and check 
  # for matches with the promote gene (pattern for exact match)
  out <- unlist(apply(symbol_db_modx,1,function(x){
    
    if (grepl("\\|",x[3])){
      genes <- unlist(strsplit(x[3],"\\|"))
    }else{
      genes <- x[3]
    }
    
    if (sum(grepl(pattern=pattern,genes))>=1){
      out <- x[1] #output new_symbol gene name
    }else{
      NULL # can't find promote gene in any previous symbol
    }
  }))
  
  if (is.null(unlist(out))){
    new_promote_genes <- gene # when you can't find it then you keep the same symbol
  }else if (length(unlist(out))==1){
    new_promote_genes <- as.character(out) # this line could be deleted
  }else{
    new_promote_genes <- as.character(out)[1]} 
  return(new_promote_genes)
}

# Load files ####
message("Loading files")
# Expression matrix PROMOTE VIS2 RAW COUNTS
PROMOTE_gcounts_Visit2 <- read_excel("data/PROMOTE_gcounts_visit2.xls",sheet = "gcrawV2")
prostate_raw_counts <- as.data.frame(PROMOTE_gcounts_Visit2[,2:length(PROMOTE_gcounts_Visit2)])
rownames(prostate_raw_counts) <-  PROMOTE_gcounts_Visit2[,1, drop=T]

# Clinical outcome (TTC)
clinical_promote <- read_excel("data/PROMOTE_gcounts_visit2.xls",sheet = "promote_clinical")
coldata <- as.data.frame(clinical_promote[,c(2,5,6)])
rownames(coldata) <- clinical_promote$EXN
coldata <-na.omit(coldata)

# Tidy data.frame ####
coldata$condition <- factor(clinical_promote$TTC.Quartile,levels=c("High","Middle", "Low"))
coldata$Biopsy_site_Visit_2 <- factor(coldata$Biopsy_site_Visit_2,
                                      levels=c("bone","lymph node","liver", "prostate bed","lung","Soft Tissue"),
                                      labels=c("bone","lymph node","liver", "prostate bed","lung","Soft Tissue"))

coldata$TTC.adjust_log <- coldata$TTC.adjust

coldata$TTC.adjust_log <- log10(coldata$TTC.adjust_log + (1-min(coldata$TTC.adjust_log)))

# Create SE ####
ddsSE <- DESeqDataSetFromMatrix(countData = prostate_raw_counts[,rownames(coldata)],
                                colData = coldata,
                                design = ~ TTC.adjust_log)

promote_genes <- rownames(ddsSE)

# Update gene names Slow step! ####
message("Matching gene symbol IDs for update")

new_promote_genes <- sapply(promote_genes, get_old_symbol, symbol_db_modx)

rowData(ddsSE)$new_symbol <- new_promote_genes

# Save ####
saveRDS(ddsSE,'data/promoteV2_SE.rds')
message("Done!")

# Checks ####
# identical(rownames(ddsSE),names(new_promote_genes))
# # [1] TRUE
# identical(rownames(ddsSE),new_promote_genes)
# # [1] FALSE

# Save in case
# save_promote_symbol_genes <- data.frame(promote_genes=promote_genes,new_symbol=new_promote_genes)
# write.table(save_promote_symbol_genes, 'data/save_promote_symbol_genes.tsv', sep='\t',row.names = F,col=T)




