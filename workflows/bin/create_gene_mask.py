#! /usr/local/bin/python

#'#################################################################################
#'#################################################################################
#'  Train a neural network defined in Keras
#'#################################################################################
#'#################################################################################

import pickle
import csv
import numpy as np
import sys
import scipy
import functools
import operator
import pandas as pd
import h5py
from scipy.sparse import coo_matrix

sys.path.append('./')
import params

with open('inputcpgs.txt','r') as file:
    cpgs = file.read()
cpgs = cpgs.split('\n')[0:-1]

## Select map rows with cpgs in input
with open('cpgs_map.txt','r') as file:
    cpg_map = pd.read_csv(file, delimiter = ' ', index_col = False)


cpg_map_filt = cpg_map[cpg_map['cpgs'].isin(cpgs)]
gene_count = cpg_map_filt['gene'].value_counts()
small_genes = gene_count.index[gene_count < 3]
mask = cpg_map_filt['gene'].isin(small_genes.to_list())

cpg_map_filt.loc[cpg_map_filt.index[mask], 'gene'] = "Intergenic"
gene_count2 = cpg_map_filt['gene'].value_counts()

ff = dict(zip(cpgs,range(len(cpgs))))
cpg_map_filt['idx'] = [ff[cpg] for cpg in cpg_map_filt['cpgs']]

## Create list of inputs - gene
gene_dict = {}
for g,idx in zip(cpg_map_filt['gene'].values,cpg_map_filt['idx'].values):
  if g not in gene_dict:
    gene_dict[g] = [idx]
  else:
    gene_dict[g].append(idx)
    

## Create list of inputs - chr
if params.chr_genes:
  cpg_map_filt = cpg_map_filt[cpg_map_filt['gene'] == "Intergenic"]
  

cpg_map_filt2 = cpg_map_filt[['cpgs', 'chromosome', 'idx']].drop_duplicates()

chr_dict = {}
for ch,idx in zip(cpg_map_filt2['chromosome'].values, cpg_map_filt2['idx'].values):
  if ch not in chr_dict:
    chr_dict[ch] = [idx]
  else:
    chr_dict[ch].append(idx)

## Define gene mask #### 
ngenes = len(gene_dict.keys())*params.gene_neurons
nchr = len(chr_dict.keys())*params.chr_neurons

mask_d  = np.zeros((len(cpgs), ngenes + nchr), bool)

gene_names = []

gcol = 0
for gene in gene_dict.keys():
  i = gcol*params.gene_neurons
  mask_d[gene_dict[gene], i:i + params.gene_neurons] = True
  gcol += 1
  gene_names.append(gene)
  
ccol = 0
chr_names = [] 
for chrom in chr_dict.keys():
  i = ngenes + ccol*params.chr_neurons
  mask_d[chr_dict[chrom], i:i + params.chr_neurons] = True
  ccol += 1
  chr_names.append(chrom)
  
gene_mask =  coo_matrix(mask_d)
pickle.dump( gene_mask, open("gene_mask.pb", "wb" ), protocol = 4 )

