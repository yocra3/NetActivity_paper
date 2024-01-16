#'#################################################################################
#'#################################################################################
#' Standardize Gene expression data
#'#################################################################################
#'#################################################################################

import h5py
import numpy as np
from sklearn.preprocessing import LabelEncoder
import pickle

## GTEx
f = h5py.File('results/GTEx/vst_all_assays.h5', 'r')
meth_matrix = f['assay001']
gtex = meth_matrix[...]
f.close()

with open('results/GTEx/individuals_labels.txt','r') as file:
    tissue = file.read()
tissue = tissue.split('\n')[0:-1]

## Convert labels to integers labels
### integer encode
label_encoder = LabelEncoder()
label_int = label_encoder.fit_transform(tissue)

means_gtex = np.mean(gtex, axis = 0)
stds_gtex = np.std(gtex, axis = 0)
x_gtex = (gtex - means_gtex)/stds_gtex

### Save reshaped training data
f = h5py.File('results/GTEx/all_reshaped_standardized.h5', 'w')
dataset_input = f.create_dataset('methy', (x_gtex.shape[0], x_gtex.shape[1]))
dataset_input[...] = x_gtex
dataset_label = f.create_dataset('label', (len(label_int),))
dataset_label[...] = label_int
dataset_label.attrs['labels'] = label_encoder.classes_.tolist()
f.close()


## Transform PRAD
f = h5py.File('results/TCGA_gexp_coding_noPRAD/vsd_norm_pradassays.h5', 'r')
meth_matrix = f['assay001']
x_prad = meth_matrix[...]
f.close()

means_prad = np.mean(x_prad, axis = 0)
stds_prad = np.std(x_prad, axis = 0)
x_prad_std = (x_prad - means_prad)/stds_prad

### Save reshaped training data
f = h5py.File('results/TCGA_gexp_coding_noPRAD/prad_assay_reshaped_standardized.h5', 'w')
dataset_input = f.create_dataset('methy', (x_prad_std.shape[0], x_prad_std.shape[1]))
dataset_input[...] = x_prad_std
f.close()

## Transform GSE169038
f = h5py.File('results/GSE169038/network_genesassays.h5', 'r')
meth_matrix = f['assay001']
x_array = meth_matrix[...]
f.close()

means_array = np.mean(x_array, axis = 0)
stds_array = np.std(x_array, axis = 0)
x_prad_array = (x_array - means_array)/stds_array
x_prad_array[:, stds_array == 0] = 0

### Save reshaped training data
f = h5py.File('results/GSE169038/prad_array_reshaped_standardized.h5', 'w')
dataset_input = f.create_dataset('methy', (x_prad_array.shape[0], x_prad_array.shape[1]))
dataset_input[...] = x_prad_array
f.close()

## Transform SRP042228
f = h5py.File('results/SRP042228/vsd_norm_TCGA_codingGenes_assays.h5', 'r')
gexp = f['assay001']
x_train = gexp[...]
f.close()

### Standardized based on GSE values
x_train_std = (x_train - np.mean(x_train, axis = 0))/np.std(x_train, axis = 0)

f = h5py.File('results/SRP042228/assay_reshaped_coding_std_gse.h5', 'w')
dataset_input = f.create_dataset('methy', (x_train.shape[0], x_train.shape[1]))
dataset_input[...] = x_train_std
f.close()
