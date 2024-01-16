#'#################################################################################
#'#################################################################################
#'                    
#' Download and filter newer version of gene SYMBOL ID to update PROMOTE gene names
#'                      
#'#################################################################################
#'#################################################################################

# Load libraries ####
library(dplyr)
library(readxl)
library(SummarizedExperiment)
library(dplyr)
library(DESeq2)

# HGNC ####
download.file("https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/hgnc_complete_set_2022-10-01.txt",destfile = "data/hgnc_complete_set_2022-10-01.tsv")

symbol_db <- read.csv('data/hgnc_complete_set_2022-10-01.tsv', sep="\t")
symbol_db <- symbol_db %>%
  dplyr::filter(locus_group =="protein-coding gene") %>%
  dplyr::select(ensembl_gene_id, symbol,entrez_id, name, prev_symbol, date_symbol_changed,
                locus_group, location, gene_group, gene_group_id, uniprot_ids, pubmed_id, enzyme_id )

# Select only lines with previous symbol (that could be in PROMOTE annotations)
symbol_db_modx <- symbol_db[symbol_db$prev_symbol!="",]

# Remove lines with duplicated prev_symbol! Here we lose some new names that could be in GTEx model..
symbol_db_modx <- symbol_db_modx[!duplicated(symbol_db_modx$prev_symbol),]
head(symbol_db_modx,10)

# Save ####
write.table(symbol_db_modx,"data/hgnc_filter_set_2022-10-01.tsv", sep"\t", row.names = F, col=T)

