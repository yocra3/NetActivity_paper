#! /usr/local/bin/python

#'#################################################################################
#'#################################################################################
#' Transform hdf5 data to list with input for biodnn model
#'#################################################################################
#'#################################################################################

import h5py
import pickle
import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelEncoder

# Load and reshape data
f = h5py.File('assays.h5', 'r')
meth_matrix = f['assay001']
input_matrix = meth_matrix[...]
f.close()

A = open('input_cpgs.pb', 'rb')
input_cpgs =  pickle.load(A)
A.close()

with open('input_cpgs.txt','r') as file:
    cpgs = file.read()
cpgs = cpgs.split('\n')[0:-1]

### Create input
tab = pd.DataFrame({'cpgs': cpgs, 'idx': range(len(cpgs))})
tab = tab.set_index('cpgs')

input_list = []
for cpgs_l in input_cpgs:
  subset = tab.loc[cpgs_l]
  selcpgs = subset['idx'].to_list()
  input_list.append(input_matrix[:, selcpgs])

pickle.dump( input_list, open( "input_list.pb", "wb" ), protocol = 4 )
